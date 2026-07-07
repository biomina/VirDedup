"""
Step 2: Cross-database matching (GenBank ↔ GISAID).
Uses fixed-rule matching (no scoring) to identify duplicate records.

A cross-database MATCH requires ALL of the following:
  [1] Sequence hash identical (exact nucleotide match)
  [2] Subtype matches (A ↔ A or B ↔ B)
  [3] Country matches (skip if missing on either side)
  [4] Collection date matches at best available granularity
  [5] Isolate name is similar enough:
      - Substring match (partial_ratio == 100), OR
      - Shared sample-identifying token AND token_sort_ratio >= 30

      Shared token detection (has_shared_id) uses five rules:
        a) Shared token containing digits (len>=2, not a year)
        b) Two or more shared alphabetic tokens (len>=2, not years)
        c) Shared numeric substring >=3 digits after letter-digit split
        d) Token substring containment (2022JECVB / JECVB)
        e) Single shared alphabetic token >=5 chars (YXLSLF type)

Records meeting [1]+[2] but failing any of [3]-[5] are flagged as EDGE_CASE
for manual review, EXCEPT:
  - Country mismatch (both valid but differ) → not an edge case
  - Collection date same granularity but differs → not an edge case

Edge case categories:
  - "Isolate: numeric code"    — GB isolate is pure 7+ digit code (UMC lab ID)
  - "Isolate: structured name" — alphanumeric isolate, no shared ID token
  - "Date: mismatch" / "Length: mismatch" — other metadata conflicts

Outputs:
  cross_database_matches.csv  -- all confirmed match pairs
  edge_cases.csv             -- records needing manual review
"""
import argparse
import csv
import re
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

import config
from harmonize_metadata import (
    dates_match_at_best_granularity,
    normalize_date_genbank,
    normalize_date_gisaid,
    parse_subtype_genbank,
    parse_subtype_gisaid,
    parse_location_genbank,
    parse_location_gisaid,
    normalize_country,
    normalize_host,
    parse_isolate_gisaid,
)

try:
    from rapidfuzz import fuzz
except ImportError:
    raise ImportError("Please install rapidfuzz: pip install rapidfuzz")


# ═══════════════════════════════════════════════════════════════════════════
#  Isolate name similarity
# ═══════════════════════════════════════════════════════════════════════════

def has_shared_id(gb_isolate, gi_isolate):
    """
    Check if two isolate names share a sample-identifying token.
    Avoids false matches from coincidental tokens (years, lab prefixes).

    Five rules (any one can trigger a match):
      1. Shared token containing digits (len>=2, not a year, not a lab/visit code)
      2. Two or more shared alphabetic tokens (len>=2, not years)
      3. Shared numeric substring >=3 digits after letter-digit split
      4. Token substring containment (e.g. 2022JECVB / JECVB)
      5. Single shared alphabetic token >=5 chars (e.g. YXLSLF)
    """
    def is_year(t):
        return len(t) == 4 and t.isdigit() and 1900 <= int(t) <= 2030

    gb_tok = [t for t in re.split(r'[^a-zA-Z0-9]+', gb_isolate.lower()) if t]
    gi_tok = [t for t in re.split(r'[^a-zA-Z0-9]+', gi_isolate.lower()) if t]

    # If both names have excluded lab/visit tokens but those tokens differ,
    # they are distinct samples from a structured-lab naming scheme
    # (e.g. SE02-0021-D04 vs SE02-0021-D02 — same study, different visits).
    gb_ex = {t for t in gb_tok if t in config.EXCLUDED_SHARED_TOKENS}
    gi_ex = {t for t in gi_tok if t in config.EXCLUDED_SHARED_TOKENS}
    if gb_ex and gi_ex and gb_ex != gi_ex:
        return False

    common = set(gb_tok) & set(gi_tok)

    # Rule 1: Shared token containing digits, len>=2, not a year
    # Exclude known lab/institution/visit codes (e.g. ox02, v01) that
    # create false cross-sample matches.
    digit_tokens = {
        t for t in common
        if len(t) >= 2 and not is_year(t) and re.search(r'\d', t)
        and t not in config.EXCLUDED_SHARED_TOKENS
    }
    if digit_tokens:
        return True

    # Rule 2: Multiple (>=2) shared alphabetic tokens, len>=2, not years
    alpha_tokens = {t for t in common if len(t) >= 2 and not is_year(t) and not re.search(r'\d', t)}
    if len(alpha_tokens) >= 2:
        return True

    # Rule 3: Shared numeric substring >=3 digits after letter-digit split
    # (catches cases like PV058/058 where tokenisation differs)
    gb_parts = [p for t in gb_tok for p in re.findall(r'[a-zA-Z]+|\d+', t)]
    gi_parts = [p for t in gi_tok for p in re.findall(r'[a-zA-Z]+|\d+', t)]
    gb_nums = [t for t in gb_parts if t.isdigit() and len(t) >= 3 and not is_year(t)]
    gi_nums = [t for t in gi_parts if t.isdigit() and len(t) >= 3 and not is_year(t)]
    for gn in gb_nums:
        for gin in gi_nums:
            if gn == gin or gn in gin or gin in gn:
                return True

    # Rule 4: Token substring containment (catches year+ID like 2022JECVB / JECVB)
    for gt in gb_tok:
        if len(gt) >= 3 and not is_year(gt):
            for git in gi_tok:
                if git != gt and (gt in git or git in gt):
                    return True

    # Rule 5: Single shared alphabetic token >= 5 chars (likely a sample ID
    # like YXLSLF, not a lab prefix like OHSU)
    long_alpha = {t for t in common if len(t) >= 5 and not is_year(t) and not re.search(r'\d', t)}
    if long_alpha:
        return True

    return False


# ═══════════════════════════════════════════════════════════════════════════
#  Matching logic
# ═══════════════════════════════════════════════════════════════════════════

def check_match(gb_row, gi_row):
    """
    Check if a GenBank record and a GISAID record are duplicates.
    Returns (is_match: bool, is_edge: bool, reasons: list).

    is_match = True if ALL 5 rules pass.
    is_edge  = True if sequence+subtype match but metadata rules fail.
    """
    reasons = []

    # [1] Sequence hash match (always checked before calling this function)
    # (Handled by the caller -- we only get called when hashes are equal)
    reasons.append("seq_hash_match")

    # [2] Subtype match
    _EMPTY = frozenset(("", "nan", "na", "n/a", "none"))
    gb_subtype_raw = str(gb_row.get("Subtype", "")).strip().upper()
    gi_subtype_raw = str(gi_row.get("Subtype", "")).strip().upper()
    gb_subtype = "" if gb_subtype_raw.lower() in _EMPTY else gb_subtype_raw
    gi_subtype = "" if gi_subtype_raw.lower() in _EMPTY else gi_subtype_raw
    both_present = bool(gb_subtype) and bool(gi_subtype)
    subtype_match = both_present and (gb_subtype == gi_subtype)
    if subtype_match:
        reasons.append(f"subtype_match:{gb_subtype}")
    elif not both_present:
        reasons.append("subtype_skipped_missing")
    else:
        reasons.append(f"subtype_mismatch:{gb_subtype}vs{gi_subtype}")
        return False, False, reasons

    # [3] Country match (skip if missing on either side)
    gb_country = str(gb_row.get("Country", "")).strip()
    gi_country = str(gi_row.get("Country", "")).strip()
    gb_country_missing = not gb_country or gb_country.lower() in ("nan", "na", "n/a", "none", "")
    gi_country_missing = not gi_country or gi_country.lower() in ("nan", "na", "n/a", "none", "")
    if not gb_country_missing and not gi_country_missing:
        country_match = (gb_country.lower() == gi_country.lower())
        if country_match:
            reasons.append(f"country_match:{gb_country}")
        else:
            reasons.append(f"country_mismatch:{gb_country}vs{gi_country}")
            return False, False, reasons  # different countries → different samples
    else:
        reasons.append("country_skipped_missing")

    # [4] Collection date match at best granularity
    gb_date = str(gb_row.get("Collection_Date_Norm", "")).strip()
    gi_date = str(gi_row.get("Collection_Date_Norm", "")).strip()
    date_match = dates_match_at_best_granularity(gb_date, gi_date)
    if date_match:
        reasons.append(f"date_match:{gb_date}|{gi_date}")
    else:
        # Determine if both dates have the same precision.
        # Same granularity + different values = clearly different samples, not edge case.
        # Different granularity = ambiguous due to differing precision, flag as edge case.
        def _gran(d):
            if re.match(r'^\d{4}$', d): return 3        # year (YYYY)
            if re.match(r'^\d{4}-\d{2}$', d): return 2  # month (YYYY-MM)
            if re.match(r'^\d{4}-\d{2}-\d{2}$', d): return 1  # day (YYYY-MM-DD)
            return 0
        same_gran = _gran(gb_date) == _gran(gi_date) and _gran(gb_date) > 0
        gb_iso = str(gb_row.get("Isolate", "")).strip()
        gi_iso = str(gi_row.get("Isolate", "")).strip()
        isolates_alike = (
            gb_iso and gi_iso
            and gb_iso.lower() not in ("nan", "na", "n/a", "none", "")
            and gi_iso.lower() not in ("nan", "na", "n/a", "none", "")
            and (fuzz.partial_ratio(gb_iso.lower(), gi_iso.lower()) == 100
                 or (has_shared_id(gb_iso, gi_iso)
                     and fuzz.token_sort_ratio(gb_iso.lower(), gi_iso.lower()) >= config.ISOLATE_SIMILARITY_THRESHOLD))
        )
        if isolates_alike:
            reasons.append(f"date_mismatch:{gb_date}vs{gi_date}")
            return False, True, reasons  # edge case — isolates share sample ID but dates differ
        if same_gran:
            reasons.append(f"date_mismatch:{gb_date}vs{gi_date}")
            return False, False, reasons  # clearly different dates
        else:
            reasons.append(f"date_mismatch:{gb_date}vs{gi_date}")
            return False, True, reasons  # edge case

    # [5] Isolate name similarity check (length is redundant — same hash = same length).
    gb_isolate = str(gb_row.get("Isolate", "")).strip()
    gi_isolate = str(gi_row.get("Isolate", "")).strip()
    # Treat missing/NaN-like values as missing data, skip the check
    if gb_isolate and gi_isolate and gb_isolate.lower() not in ("nan", "na", "n/a", "none", "") and gi_isolate.lower() not in ("nan", "na", "n/a", "none", ""):
        # GB isolate is a pure numeric code (7+ digits) — lab-assigned ID with
        # no semantic relationship to GISAID names. Skip the comparison and rely
        # on sequence + core metadata (subtype, country, date, length) instead.
        if re.match(r'^\d{7,}$', gb_isolate):
            reasons.append(f"isolate_mismatch:numeric_code|{gb_isolate}|{gi_isolate}")
            return False, True, reasons  # edge case — pure numeric code has no semantic relationship
        else:
            pr = fuzz.partial_ratio(gb_isolate.lower(), gi_isolate.lower())
            ts = fuzz.token_sort_ratio(gb_isolate.lower(), gi_isolate.lower())
            shared = has_shared_id(gb_isolate, gi_isolate)
            if pr == 100 or (shared and ts >= config.ISOLATE_SIMILARITY_THRESHOLD):
                reasons.append(f"isolate_match:{max(pr, ts):.0f}|{gb_isolate}|{gi_isolate}")
            else:
                reasons.append(f"isolate_mismatch:{max(pr, ts):.0f}|{gb_isolate}|{gi_isolate}")
                return False, True, reasons  # edge case
    else:
        reasons.append("isolate_skipped_missing")

    # All checks passed
    return True, False, reasons


# ═══════════════════════════════════════════════════════════════════════════
#  Main matching
# ═══════════════════════════════════════════════════════════════════════════

def match_cross_database():
    """Run cross-database matching between deduped GenBank and GISAID records."""
    print("=" * 60)
    print("CROSS-DATABASE MATCHING")
    print("=" * 60)

    # ── 1. Load deduped metadata ─────────────────────────────────────────
    print("\n[1/4] Loading deduped metadata ...")
    gb_meta = pd.read_csv(config.DEDUPED_GENBANK_METADATA, low_memory=False)
    gi_meta = pd.read_csv(config.DEDUPED_GISAID_METADATA, low_memory=False)
    print(f"  GenBank records: {len(gb_meta):,}")
    print(f"  GISAID records:  {len(gi_meta):,}")

    # Ensure SeqHash is string
    gb_meta["SeqHash"] = gb_meta["SeqHash"].astype(str)
    gi_meta["SeqHash"] = gi_meta["SeqHash"].astype(str)

    # ── 2. Build hash -> records index ────────────────────────────────────
    print("\n[2/4] Building sequence hash index ...")
    gb_by_hash = defaultdict(list)
    for idx, row in gb_meta.iterrows():
        gb_by_hash[row["SeqHash"]].append(row)

    gi_by_hash = defaultdict(list)
    for idx, row in gi_meta.iterrows():
        gi_by_hash[row["SeqHash"]].append(row)

    gb_hashes = set(gb_by_hash.keys())
    gi_hashes = set(gi_by_hash.keys())
    shared_hashes = gb_hashes & gi_hashes
    gb_only_hashes = gb_hashes - gi_hashes
    gi_only_hashes = gi_hashes - gb_hashes

    print(f"  Shared sequence hashes: {len(shared_hashes):,}")
    print(f"  GenBank-only hashes:   {len(gb_only_hashes):,}")
    print(f"  GISAID-only hashes:    {len(gi_only_hashes):,}")

    # ── 3. Match ─────────────────────────────────────────────────────────
    print("\n[3/4] Matching records (fixed rules) ...")

    matches = []              # confirmed duplicates
    edge_cases = []           # same sequence+subtype, ambiguous metadata conflict
    considered_rejected = []  # same sequence+subtype, clearly different metadata
    gb_only_records = []
    gi_only_records = []

    for seq_hash in tqdm(shared_hashes, desc="  Processing shared hashes"):
        gb_rows = gb_by_hash[seq_hash]
        gi_rows = gi_by_hash[seq_hash]

        matched_gi = set()  # track which GI records already matched

        def _make_record(gb_row, gi_row, reasons, cat, seq_hash, gi_acc):
            """Build a standard dict for pair-level output."""
            return {
                "SeqHash": seq_hash,
                "GenBank_Accession": gb_row.get("Accession", ""),
                "GISAID_Accession": gi_acc,
                "Category": cat,
                "Reason": "; ".join(reasons),
                "GB_Subtype": gb_row.get("Subtype", ""),
                "GI_Subtype": gi_row.get("Subtype", ""),
                "GB_Country": gb_row.get("Country", ""),
                "GI_Country": gi_row.get("Country", ""),
                "GB_Date": gb_row.get("Collection_Date_Norm", ""),
                "GI_Date": gi_row.get("Collection_Date_Norm", ""),
                "GB_Length": gb_row.get("Length", ""),
                "GI_Length": gi_row.get("Length", ""),
                "GB_Isolate": gb_row.get("Isolate", ""),
                "GI_Isolate": gi_row.get("Isolate", ""),
            }

        for gb_row in gb_rows:
            best_match = None
            best_match_tuple = None

            for gi_row in gi_rows:
                gi_acc = gi_row["Accession"]
                if gi_acc in matched_gi:
                    continue  # already matched to another GB record

                is_match, is_edge, reasons = check_match(gb_row, gi_row)

                if is_match:
                    best_match = (gb_row, gi_row, reasons)
                    best_match_tuple = gi_acc
                    break  # take first match

                # Track all non-match pairs (edge cases + clean rejections)
                reasons_str = " ".join(r for r in reasons if isinstance(r, str))
                if "subtype_match" in reasons_str or "subtype_skipped_missing" in reasons_str:
                    gb_iso = str(gb_row.get("Isolate", "")).strip().lower()
                    has_subtype_skip = any("subtype_skipped_missing" in r for r in reasons if isinstance(r, str))
                    has_iso_issue = any("isolate_mismatch" in r for r in reasons if isinstance(r, str))
                    has_date_issue = any("date_mismatch" in r for r in reasons if isinstance(r, str))
                    has_country_issue = any("country_mismatch" in r for r in reasons if isinstance(r, str))
                    if has_subtype_skip:
                        cat = "Subtype: missing"
                    elif has_iso_issue and re.match(r'^\d{7,}$', gb_iso):
                        cat = "Isolate: numeric code"
                    elif has_iso_issue:
                        cat = "Isolate: structured name"
                    elif has_date_issue:
                        cat = "Date: mismatch"
                    elif has_country_issue:
                        cat = "Country: mismatch"
                    else:
                        cat = "Other"
                    rec = _make_record(gb_row, gi_row, reasons, cat, seq_hash, gi_acc)
                    if is_edge:
                        edge_cases.append(rec)
                    else:
                        considered_rejected.append(rec)

            if best_match is not None:
                gb_row, gi_row, reasons = best_match
                matches.append({
                    "SeqHash": seq_hash,
                    "GenBank_Accession": gb_row.get("Accession", ""),
                    "GISAID_Accession": gi_row.get("Accession", ""),
                    "Match_Reasons": "; ".join(reasons),
                    "GB_Subtype": gb_row.get("Subtype", ""),
                    "GI_Subtype": gi_row.get("Subtype", ""),
                    "GB_Country": gb_row.get("Country", ""),
                    "GI_Country": gi_row.get("Country", ""),
                    "GB_Date": gb_row.get("Collection_Date_Norm", ""),
                    "GI_Date": gi_row.get("Collection_Date_Norm", ""),
                    "GB_Length": gb_row.get("Length", ""),
                    "GI_Length": gi_row.get("Length", ""),
                    "GB_Isolate": gb_row.get("Isolate", ""),
                    "GI_Isolate": gi_row.get("Isolate", ""),
                })
                matched_gi.add(best_match_tuple)

        # Remaining unmatched GI records for this hash are GISAID-only
        for gi_row in gi_rows:
            if gi_row["Accession"] not in matched_gi:
                gi_only_records.append(gi_row)

    # GB-only hashes
    for seq_hash in tqdm(gb_only_hashes, desc="  Processing GB-only hashes"):
        for gb_row in gb_by_hash[seq_hash]:
            gb_only_records.append(gb_row)

    # GI-only hashes (already partially handled above, but also from GI-only hashes)
    for seq_hash in tqdm(gi_only_hashes, desc="  Processing GI-only hashes"):
        for gi_row in gi_by_hash[seq_hash]:
            gi_only_records.append(gi_row)

    # Remove any edge case / considered-rejected record whose GenBank or
    # GISAID accession already appears in a confirmed match pair.
    matched_gb = {m["GenBank_Accession"] for m in matches}
    matched_gi = {m["GISAID_Accession"] for m in matches}
    edge_cases = [
        e for e in edge_cases
        if e["GenBank_Accession"] not in matched_gb
        and e["GISAID_Accession"] not in matched_gi
    ]
    considered_rejected = [
        c for c in considered_rejected
        if c["GenBank_Accession"] not in matched_gb
        and c["GISAID_Accession"] not in matched_gi
    ]

    gb_only_seq_hashes = {str(r.get("SeqHash", "")) for r in gb_only_records}
    gb_only_seq_hashes.discard("")
    gi_only_seq_hashes = {str(r.get("SeqHash", "")) for r in gi_only_records}
    gi_only_seq_hashes.discard("")

    print(f"\n  Matches (confirmed duplicates): {len(matches):,}")
    print(f"  Edge cases (needs review):    {len(edge_cases):,}")
    print(f"  Considered but rejected:      {len(considered_rejected):,}")
    print(f"  GenBank-only unique sequences: {len(gb_only_seq_hashes):,}")
    print(f"  GenBank-only records:         {len(gb_only_records):,}")
    print(f"  GISAID-only unique sequences: {len(gi_only_seq_hashes):,}")
    print(f"  GISAID-only records:          {len(gi_only_records):,}")

    # ── 4. Write outputs ─────────────────────────────────────────────────
    print("\n[4/4] Writing output files ...")

    # Cross-database matches
    if matches:
        matches_df = pd.DataFrame(matches)
        matches_df.to_csv(config.CROSS_MATCHES_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.CROSS_MATCHES_CSV} ({len(matches_df):,} match pairs)")
    else:
        pd.DataFrame(columns=[
            "SeqHash", "GenBank_Accession", "GISAID_Accession", "Match_Reasons",
            "GB_Subtype", "GI_Subtype", "GB_Country", "GI_Country",
            "GB_Date", "GI_Date", "GB_Length", "GI_Length",
            "GB_Isolate", "GI_Isolate",
        ]).to_csv(config.CROSS_MATCHES_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.CROSS_MATCHES_CSV} (empty)")

    # Edge cases
    if edge_cases:
        edge_df = pd.DataFrame(edge_cases)
        edge_df.to_csv(config.EDGE_CASES_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.EDGE_CASES_CSV} ({len(edge_df):,} edge cases)")
    else:
        pd.DataFrame(columns=[
            "SeqHash", "GenBank_Accession", "GISAID_Accession", "Category", "Reason",
            "GB_Subtype", "GI_Subtype", "GB_Country", "GI_Country",
            "GB_Date", "GI_Date", "GB_Length", "GI_Length",
            "GB_Isolate", "GI_Isolate",
        ]).to_csv(config.EDGE_CASES_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.EDGE_CASES_CSV} (empty)")

    # Considered but rejected pairs
    if considered_rejected:
        cr_df = pd.DataFrame(considered_rejected)
        cr_df.to_csv(config.CONSIDERED_REJECTED_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.CONSIDERED_REJECTED_CSV} ({len(cr_df):,} considered-rejected pairs)")
    else:
        pd.DataFrame(columns=[
            "SeqHash", "GenBank_Accession", "GISAID_Accession", "Category", "Reason",
            "GB_Subtype", "GI_Subtype", "GB_Country", "GI_Country",
            "GB_Date", "GI_Date", "GB_Length", "GI_Length",
            "GB_Isolate", "GI_Isolate",
        ]).to_csv(config.CONSIDERED_REJECTED_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.CONSIDERED_REJECTED_CSV} (empty)")

    return matches, edge_cases, considered_rejected, gb_only_records, gi_only_records


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cross-database matching (Step 2).")
    parser.add_argument("--deduped-gb-meta", default=None, help="Deduped GenBank metadata CSV")
    parser.add_argument("--deduped-gi-meta", default=None, help="Deduped GISAID metadata CSV")
    parser.add_argument("--output-dir", default=None, help="Output directory for match files")
    args = parser.parse_args()

    if args.deduped_gb_meta is not None:
        config.DEDUPED_GENBANK_METADATA = args.deduped_gb_meta
    if args.deduped_gi_meta is not None:
        config.DEDUPED_GISAID_METADATA = args.deduped_gi_meta
    if args.output_dir is not None:
        config.set_output_dir(args.output_dir)

    matches, edge_cases, considered_rejected, gb_only, gi_only = match_cross_database()
    gb_only_seq_hashes = {str(r.get("SeqHash", "")) for r in gb_only}
    gb_only_seq_hashes.discard("")
    gi_only_seq_hashes = {str(r.get("SeqHash", "")) for r in gi_only}
    gi_only_seq_hashes.discard("")
    print("\n" + "=" * 60)
    print("CROSS-DATABASE MATCH SUMMARY")
    print("=" * 60)
    print(f"  Confirmed duplicates (to remove GISAID copy): {len(matches):,}")
    print(f"  Edge cases (manual review required):          {len(edge_cases):,}")
    print(f"  Considered but rejected:                      {len(considered_rejected):,}")
    print(f"  GenBank-only unique sequences:                {len(gb_only_seq_hashes):,}")
    print(f"  GenBank-only records:                         {len(gb_only):,}")
    print(f"  GISAID-only unique sequences:                 {len(gi_only_seq_hashes):,}")
    print(f"  GISAID-only records:                          {len(gi_only):,}")
    print("\nDone with cross-database matching.")
