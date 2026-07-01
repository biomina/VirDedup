"""
Step 3: Generate final deduplicated output files.

Reads intermediate outputs from steps 1 and 2, then produces:
  - deduplicated_sequences.fasta   -- clean FASTA (one sequence per unique entry)
  - deduplicated_metadata.csv      -- harmonized metadata
  - removed_sequences.fasta        -- all removed sequences with original headers
  - removed_sequences.csv          -- side-by-side removed + kept metadata
  - deduplication_report.txt       -- summary statistics
"""
import argparse
import csv
import os
from collections import defaultdict
from datetime import datetime

import re
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

import config
from harmonize_metadata import (
    extract_genbank_accession,
    extract_gisaid_accession,
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


# ═══════════════════════════════════════════════════════════════════════════
#  Output schema helpers
# ═══════════════════════════════════════════════════════════════════════════

# Harmonized metadata columns for the final deduplicated output
HARMONIZED_COLS = [
    "Primary_Accession",
    "GISAID_Accession",
    "Source",
    "Integration_Status",
    "Subtype",
    "Isolate",
    "Country",
    "Region",
    "Collection_Date",
    "Length",
    "Host",
    "Genotype",
    "Lineage",
    "Release_Date",
    "Submitters",
]


def _clean_column_name(name) -> str:
    """Make source metadata column names stable for CSV consumers."""
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_")
    return re.sub(r"_+", "_", cleaned) or "Column"


def _prefixed(prefix: str, field: str) -> str:
    """Return a harmonized field name, with no leading underscore for final rows."""
    return f"{prefix}_{field}" if prefix else field


def add_source_metadata(data: dict, row, source_prefix: str) -> None:
    """Append all raw/intermediate metadata columns with a source-specific prefix."""
    if row is None:
        return
    for col, val in row.items():
        data[f"{source_prefix}_{_clean_column_name(col)}"] = _v(val)

# Side-by-side columns for removed_sequences.csv
REMOVED_COLS = []
for side in ("Removed", "Kept"):
    for field in ("Accession", "Source", "Subtype", "Isolate", "Country", "Region",
                  "Collection_Date", "Length", "Host", "Genotype", "Release_Date", "Submitter"):
        REMOVED_COLS.append(f"{side}_{field}")
REMOVED_COLS.append("Removal_Reason")


def _v(val, default=""):
    """Return val if it's a non-empty string, else default."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return default
    s = str(val).strip()
    return s if s and s.lower() not in ("nan", "nat", "none", "") else default


def extract_metadata_row(row, source_db, prefix=""):
    """Extract standardized metadata fields from a row (Series or dict)."""
    data = {}
    acc = _v(row.get("Accession", "")) or _v(row.get("GenBank_Accession", "")) or _v(row.get("GISAID_Accession", ""))
    if prefix:
        data[f"{prefix}_Accession"] = acc
    else:
        data["Primary_Accession"] = acc
    data[_prefixed(prefix, "Source")] = source_db
    data[_prefixed(prefix, "Subtype")] = _v(row.get("Subtype", ""))
    data[_prefixed(prefix, "Isolate")] = _v(row.get("Isolate", ""))
    data[_prefixed(prefix, "Country")] = _v(row.get("Country", ""))
    data[_prefixed(prefix, "Region")] = _v(row.get("Region", ""))
    dt = _v(row.get("Collection_Date_Norm", "")) or _v(row.get("Collection_Date", ""))
    data[_prefixed(prefix, "Collection_Date")] = dt
    data[_prefixed(prefix, "Length")] = _v(row.get("Length", ""))
    data[_prefixed(prefix, "Host")] = _v(row.get("Host", ""))
    data[_prefixed(prefix, "Genotype")] = _v(row.get("Genotype", ""))
    lineage = _v(row.get("Lineage", "")) or _v(row.get("GISAID_Lineage", ""))
    data[_prefixed(prefix, "Lineage")] = lineage
    rdt = _v(row.get("Release_Date_Norm", "")) or _v(row.get("Release_Date", ""))
    data[_prefixed(prefix, "Release_Date")] = rdt
    submitters = _v(row.get("Submitters", "")) or _v(row.get("Submitter", ""))
    submitter_key = f"{prefix}_Submitter" if prefix else "Submitters"
    data[submitter_key] = submitters
    return data


# ═══════════════════════════════════════════════════════════════════════════
#  Main assembly
# ═══════════════════════════════════════════════════════════════════════════

def generate_output(remove_edge_copies=False):
    """Assemble final deduplicated output files."""
    print("=" * 60)
    print("FINAL OUTPUT GENERATION")
    print("=" * 60)

    # ── 1. Load all intermediate files ───────────────────────────────────
    print("\n[1/5] Loading intermediate files ...")

    # Deduped metadata
    gb_deduped = pd.read_csv(config.DEDUPED_GENBANK_METADATA, low_memory=False)
    gi_deduped = pd.read_csv(config.DEDUPED_GISAID_METADATA, low_memory=False)
    print(f"  GenBank deduped: {len(gb_deduped):,} records")
    print(f"  GISAID deduped:  {len(gi_deduped):,} records")

    # Cross-database matches
    cross_matches = pd.read_csv(config.CROSS_MATCHES_CSV, low_memory=False)
    print(f"  Cross-database matches: {len(cross_matches):,}")

    # Edge cases (loaded for including in dedup output, keeping one copy per pair)
    edge_cases = pd.read_csv(config.EDGE_CASES_CSV, low_memory=False)
    print(f"  Edge cases (review):    {len(edge_cases):,}")

    # Build indices for quick lookup
    gb_by_acc = {}
    for idx, row in gb_deduped.iterrows():
        gb_by_acc[row.get("Accession", "")] = row

    gi_by_acc = {}
    for idx, row in gi_deduped.iterrows():
        gi_by_acc[row.get("Accession", "")] = row

    gb_by_hash = defaultdict(list)
    for idx, row in gb_deduped.iterrows():
        gb_by_hash[str(row.get("SeqHash", ""))].append(row)

    gi_by_hash = defaultdict(list)
    for idx, row in gi_deduped.iterrows():
        gi_by_hash[str(row.get("SeqHash", ""))].append(row)

    # ── 2. Build harmonized deduplicated metadata ────────────────────────
    print("\n[2/5] Building deduplicated metadata ...")

    dedup_records = []
    matched_gisaid_accessions = set()
    matched_genbank_accessions = set()

    # Cross-database matches: keep GenBank record, note GISAID accession
    for idx, match_row in cross_matches.iterrows():
        gb_acc = match_row.get("GenBank_Accession", "")
        gi_acc = match_row.get("GISAID_Accession", "")
        if not gb_acc or not gi_acc:
            continue

        gb_row = gb_by_acc.get(gb_acc)
        if gb_row is None:
            # GenBank record not in deduped metadata (shouldn't happen)
            continue
        gi_row = gi_by_acc.get(gi_acc)

        # Build harmonized record
        rec = extract_metadata_row(gb_row, "GenBank", prefix="")
        rec["GISAID_Accession"] = gi_acc
        rec["Source"] = "GenBank"
        rec["Integration_Status"] = "MATCH"
        # Replace Primary_Accession with GenBank accession
        rec["Primary_Accession"] = gb_acc
        if gi_row is not None and not rec.get("Lineage"):
            rec["Lineage"] = _v(gi_row.get("GISAID_Lineage", "")) or _v(gi_row.get("Lineage", ""))
        add_source_metadata(rec, gb_row, "GenBank")
        add_source_metadata(rec, gi_row, "GISAID")
        add_source_metadata(rec, match_row, "Match")
        dedup_records.append(rec)

        matched_gisaid_accessions.add(gi_acc)
        matched_genbank_accessions.add(gb_acc)

    # Edge cases: optionally keep GenBank record, note GISAID accession (one copy per pair)
    if remove_edge_copies:
        for idx, edge_row in edge_cases.iterrows():
            gb_acc = edge_row.get("GenBank_Accession", "")
            gi_acc = edge_row.get("GISAID_Accession", "")
            if not gb_acc or not gi_acc:
                continue

            gb_row = gb_by_acc.get(gb_acc)
            if gb_row is None:
                continue
            gi_row = gi_by_acc.get(gi_acc)

            rec = extract_metadata_row(gb_row, "GenBank", prefix="")
            rec["GISAID_Accession"] = gi_acc
            rec["Source"] = "GenBank"
            rec["Integration_Status"] = "EDGE"
            rec["Primary_Accession"] = gb_acc
            if gi_row is not None and not rec.get("Lineage"):
                rec["Lineage"] = _v(gi_row.get("GISAID_Lineage", "")) or _v(gi_row.get("Lineage", ""))
            add_source_metadata(rec, gb_row, "GenBank")
            add_source_metadata(rec, gi_row, "GISAID")
            add_source_metadata(rec, edge_row, "Edge")
            dedup_records.append(rec)

            matched_gisaid_accessions.add(gi_acc)
            matched_genbank_accessions.add(gb_acc)

    # GenBank-only records (not matched)
    for idx, row in gb_deduped.iterrows():
        acc = row.get("Accession", "")
        if acc in matched_genbank_accessions:
            continue
        rec = extract_metadata_row(row, "GenBank", prefix="")
        rec["GISAID_Accession"] = ""
        rec["Source"] = "GenBank"
        rec["Integration_Status"] = "GENBANK_ONLY"
        rec["Primary_Accession"] = acc
        add_source_metadata(rec, row, "GenBank")
        dedup_records.append(rec)

    # GISAID-only records (not matched)
    for idx, row in gi_deduped.iterrows():
        acc = row.get("Accession", "")
        if acc in matched_gisaid_accessions:
            continue
        rec = extract_metadata_row(row, "GISAID", prefix="")
        rec["GISAID_Accession"] = acc  # Primary is also GISAID for GI-only
        rec["Source"] = "GISAID"
        rec["Integration_Status"] = "GISAID_ONLY"
        rec["Primary_Accession"] = acc
        add_source_metadata(rec, row, "GISAID")
        dedup_records.append(rec)

    dedup_meta = pd.DataFrame(dedup_records)
    # Keep the compact harmonized schema first, then preserve every available
    # source/intermediate metadata column.
    first_cols = [c for c in HARMONIZED_COLS if c in dedup_meta.columns]
    extra_cols = [c for c in dedup_meta.columns if c not in first_cols]
    dedup_meta = dedup_meta[first_cols + extra_cols]

    # Also add any extra columns we want to preserve
    # Ensure consistent types
    dedup_meta = dedup_meta.fillna("")
    print(f"  Deduplicated metadata: {len(dedup_meta):,} records")

    # ── 3. Build deduplicated FASTA ──────────────────────────────────────
    print("\n[3/5] Building deduplicated FASTA ...")

    # For cross-database matches, use GenBank sequence
    # For GB-only, use GenBank sequence
    # For GI-only, use GISAID sequences

    # Determine which accessions to include in the FASTA
    fasta_entries = {}  # accession -> (header_line, sequence)

    # Matched and edge: use GenBank sequence
    for gb_acc in matched_genbank_accessions:
        fasta_entries[gb_acc] = ("gb_match", None)  # placeholder

    # GB-only
    for idx, row in gb_deduped.iterrows():
        acc = row.get("Accession", "")
        if acc not in matched_genbank_accessions:
            fasta_entries[acc] = ("gb", None)

    # GI-only
    for idx, row in gi_deduped.iterrows():
        acc = row.get("Accession", "")
        if acc not in matched_gisaid_accessions:
            fasta_entries[acc] = ("gi", None)

    # Read and write sequences
    written_count = 0
    with open(config.FINAL_DEDUP_FASTA, "w") as out_fa:
        # Process GenBank FASTA
        for record in tqdm(SeqIO.parse(config.GENBANK_FASTA, "fasta"),
                           desc="  GenBank FASTA", unit=" records"):
            acc = extract_genbank_accession(record.id)
            entry = fasta_entries.get(acc)
            if entry is not None and entry[0] in ("gb", "gb_match"):
                seq_str = str(record.seq)
                if len(seq_str) >= config.MIN_SEQUENCE_LENGTH:
                    out_fa.write(f">{acc}\n{seq_str}\n")
                    fasta_entries[acc] = (entry[0], seq_str)
                    written_count += 1

        # Process GISAID FASTA
        for record in tqdm(SeqIO.parse(config.GISAID_FASTA, "fasta"),
                           desc="  GISAID FASTA", unit=" records"):
            acc = extract_gisaid_accession(record.id)
            entry = fasta_entries.get(acc)
            if entry is not None and entry[0] == "gi":
                seq_str = str(record.seq)
                if len(seq_str) >= config.MIN_SEQUENCE_LENGTH:
                    out_fa.write(f">{acc}\n{seq_str}\n")
                    fasta_entries[acc] = (entry[0], seq_str)
                    written_count += 1

    print(f"  Deduplicated FASTA: {written_count:,} sequences")
    print(f"  -> {config.FINAL_DEDUP_FASTA}")

    # ── 4. Build removed_sequences.csv (side-by-side) ────────────────────
    print("\n[4/5] Building removed sequences output ...")

    removed_records = []

    # 4a. Intra-GenBank removed records
    if os.path.exists(config.REMOVED_INTRA_GENBANK_CSV):
        intra_gb_removed = pd.read_csv(config.REMOVED_INTRA_GENBANK_CSV, low_memory=False)
        for idx, rem_row in intra_gb_removed.iterrows():
            rem_hash = str(rem_row.get("SeqHash", ""))
            # Find the kept record from the same hash group
            kept_rows = gb_by_hash.get(rem_hash, [])
            if len(kept_rows) > 0:
                kept_row = kept_rows[0]  # first kept record with same hash
                rec = {}
                rec.update(extract_metadata_row(rem_row, "GenBank", prefix="Removed"))
                rec.update(extract_metadata_row(kept_row, "GenBank", prefix="Kept"))
                rec["Removal_Reason"] = "intra-GenBank duplicate"
                removed_records.append(rec)

    # 4b. Intra-GISAID removed records
    if os.path.exists(config.REMOVED_INTRA_GISAID_CSV):
        intra_gi_removed = pd.read_csv(config.REMOVED_INTRA_GISAID_CSV, low_memory=False)
        for idx, rem_row in intra_gi_removed.iterrows():
            rem_hash = str(rem_row.get("SeqHash", ""))
            kept_rows = gi_by_hash.get(rem_hash, [])
            if len(kept_rows) > 0:
                kept_row = kept_rows[0]
                rec = {}
                rec.update(extract_metadata_row(rem_row, "GISAID", prefix="Removed"))
                rec.update(extract_metadata_row(kept_row, "GISAID", prefix="Kept"))
                rec["Removal_Reason"] = "intra-GISAID duplicate"
                removed_records.append(rec)

    # 4c. Cross-database removed records (GISAID side removed, GenBank kept)
    for idx, match_row in cross_matches.iterrows():
        gb_acc = match_row.get("GenBank_Accession", "")
        gi_acc = match_row.get("GISAID_Accession", "")
        gb_row = gb_by_acc.get(gb_acc)
        gi_row = gi_by_acc.get(gi_acc)
        if gb_row is None or gi_row is None:
            continue
        rec = {}
        rec.update(extract_metadata_row(gi_row, "GISAID", prefix="Removed"))
        rec.update(extract_metadata_row(gb_row, "GenBank", prefix="Kept"))
        rec["Removal_Reason"] = "cross-database GISAID removed for GenBank match"
        removed_records.append(rec)

    # 4d. Edge case removed records (GISAID side removed, GenBank kept)
    if remove_edge_copies:
        for idx, edge_row in edge_cases.iterrows():
            gb_acc = edge_row.get("GenBank_Accession", "")
            gi_acc = edge_row.get("GISAID_Accession", "")
            gb_row = gb_by_acc.get(gb_acc)
            gi_row = gi_by_acc.get(gi_acc)
            if gb_row is None or gi_row is None:
                continue
            rec = {}
            rec.update(extract_metadata_row(gi_row, "GISAID", prefix="Removed"))
            rec.update(extract_metadata_row(gb_row, "GenBank", prefix="Kept"))
            rec["Removal_Reason"] = "edge case GISAID removed for GenBank match"
            removed_records.append(rec)

    if removed_records:
        removed_df = pd.DataFrame(removed_records)
        # Ensure all columns exist
        for col in REMOVED_COLS:
            if col not in removed_df.columns:
                removed_df[col] = ""
        removed_df = removed_df[[c for c in REMOVED_COLS if c in removed_df.columns]]
        removed_df.to_csv(config.FINAL_REMOVED_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.FINAL_REMOVED_CSV} ({len(removed_df):,} records)")
    else:
        pd.DataFrame(columns=REMOVED_COLS).to_csv(
            config.FINAL_REMOVED_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.FINAL_REMOVED_CSV} (empty)")

    # ── 4d. Build removed_sequences.fasta ────────────────────────────────
    print("  Building removed sequences FASTA ...")

    # Combine all removed accessions
    removed_acc_to_source = {}  # accession -> (source_db, original_header, kept_text)

    # Intra-GenBank
    if os.path.exists(config.REMOVED_INTRA_GENBANK_FASTA):
        for record in SeqIO.parse(config.REMOVED_INTRA_GENBANK_FASTA, "fasta"):
            acc = extract_genbank_accession(record.id)
            removed_acc_to_source[acc] = ("GenBank", record.description, str(record.seq))

    # Intra-GISAID
    if os.path.exists(config.REMOVED_INTRA_GISAID_FASTA):
        for record in SeqIO.parse(config.REMOVED_INTRA_GISAID_FASTA, "fasta"):
            acc = extract_gisaid_accession(record.id)
            removed_acc_to_source[acc] = ("GISAID", record.description, str(record.seq))

    # Cross-database removed GISAID records (matches + edge cases)
    gi_removed_cross_accs = set()
    for idx, match_row in cross_matches.iterrows():
        gi_acc = match_row.get("GISAID_Accession", "")
        if gi_acc:
            gi_removed_cross_accs.add(gi_acc)
    if remove_edge_copies:
        for idx, edge_row in edge_cases.iterrows():
            gi_acc = edge_row.get("GISAID_Accession", "")
            if gi_acc:
                gi_removed_cross_accs.add(gi_acc)

    if gi_removed_cross_accs:
        for record in SeqIO.parse(config.GISAID_FASTA, "fasta"):
            acc = extract_gisaid_accession(record.id)
            if acc in gi_removed_cross_accs and acc not in removed_acc_to_source:
                removed_acc_to_source[acc] = ("GISAID", record.description, str(record.seq))

    removed_fasta_count = 0
    with open(config.FINAL_REMOVED_FASTA, "w") as out_fa:
        for acc, (source, header, seq) in removed_acc_to_source.items():
            out_fa.write(f">{header}\n{seq}\n")
            removed_fasta_count += 1
    print(f"  -> {config.FINAL_REMOVED_FASTA} ({removed_fasta_count:,} sequences)")

    # ── 5. Save deduplicated metadata & write report ─────────────────────
    print("\n[5/5] Saving final metadata and report ...")

    dedup_meta.to_csv(config.FINAL_DEDUP_METADATA, index=False, quoting=csv.QUOTE_ALL)
    print(f"  -> {config.FINAL_DEDUP_METADATA}")

    # Report
    total_removed = len(removed_records)
    edge_count = len(edge_cases) if not edge_cases.empty else 0

    # Count by removal reason
    reason_counts = defaultdict(int)
    for r in removed_records:
        reason_counts[r.get("Removal_Reason", "unknown")] += 1

    # ── Edge case grouping for report ────────────────────────────────────
    def categorize_edge_case(row):
        reason = str(row.get("Reason", ""))
        cat = []
        if "isolate_mismatch" in reason:
            gb_iso = str(row.get("GB_Isolate", "")).strip().lower()
            if gb_iso in ("nan", "na", "n/a", "none", ""):
                cat.append("GB isolate missing")
            elif re.match(r'^\d{7,}$', gb_iso):
                cat.append("GB isolate is numeric code (7+ digits)")
            else:
                cat.append("GB isolate is structured name")
        if "date_mismatch" in reason:
            try:
                gd = datetime.strptime(str(row.get("GB_Date", "")).strip(), "%Y-%m-%d")
                gid = datetime.strptime(str(row.get("GI_Date", "")).strip(), "%Y-%m-%d")
                d = abs((gd - gid).days)
                if d <= 7:
                    cat.append("date differs by 1-7 days")
                elif d <= 30:
                    cat.append("date differs by 8-30 days")
                elif d <= 365:
                    cat.append("date differs by 31-365 days")
                else:
                    cat.append("date differs by >365 days")
            except:
                cat.append("date parse error")
        if "country_mismatch" in reason:
            cat.append("country mismatch")
        return "; ".join(cat) if cat else "other"

    edge_group_counts = defaultdict(int)
    for _, r in edge_cases.iterrows():
        edge_group_counts[categorize_edge_case(r)] += 1

    try:
        gb_input_count = len(pd.read_csv(config.GENBANK_METADATA))
    except Exception:
        gb_input_count = 0
    try:
        gi_input_count = len(pd.read_excel(config.GISAID_METADATA))
    except Exception:
        gi_input_count = 0

    gb_removed_count = 0
    if os.path.exists(config.REMOVED_INTRA_GENBANK_CSV):
        with open(config.REMOVED_INTRA_GENBANK_CSV, encoding="utf-8") as f:
            gb_removed_count = sum(1 for _ in f) - 1  # minus header
    gi_removed_count = 0
    if os.path.exists(config.REMOVED_INTRA_GISAID_CSV):
        with open(config.REMOVED_INTRA_GISAID_CSV, encoding="utf-8") as f:
            gi_removed_count = sum(1 for _ in f) - 1

    report_lines = [
        "=" * 60,
        "DEDUPLICATION REPORT",
        "=" * 60,
        "",
        f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "--- Input -----------------------------------------------------------",
        f"  GenBank records:          {gb_input_count:,}",
        f"  GISAID records:           {gi_input_count:,}",
        "",
        "--- Intra-database dedup --------------------------------------------",
        f"  GenBank kept:             {len(gb_deduped):,}",
        f"  GenBank removed:          {gb_removed_count:,}",
        f"  GISAID kept:              {len(gi_deduped):,}",
        f"  GISAID removed:           {gi_removed_count:,}",
        "",
        "--- Cross-database matching -----------------------------------------",
        f"  Confirmed duplicates:     {len(cross_matches):,}",
        f"  Edge cases (review):      {edge_count:,}",
        "",
        "--- Edge case breakdown ----------------------------------------------",
    ]
    for group_label, count in sorted(edge_group_counts.items(), key=lambda x: -x[1]):
        report_lines.append(f"  {count:>5,}  {group_label}")
    report_lines.extend([
        "",
        "--- Final deduplicated ----------------------------------------------",
        f"  Total unique sequences:   {written_count:,}",
        f"  Total metadata records:   {len(dedup_meta):,}",
        f"  Total removed sequences:  {removed_fasta_count:,}",
        "",
        "--- Removal breakdown -----------------------------------------------",
    ])
    for reason, count in sorted(reason_counts.items(), key=lambda x: -x[1]):
        report_lines.append(f"  {reason}: {count:,}")
    report_lines.extend([
        "",
        "--- Output files --------------------------------------------------",
        f"  {config.FINAL_DEDUP_FASTA}",
        f"  {config.FINAL_DEDUP_METADATA}",
        f"  {config.FINAL_REMOVED_FASTA}",
        f"  {config.FINAL_REMOVED_CSV}",
        f"  {config.CROSS_MATCHES_CSV}",
        f"  {config.EDGE_CASES_CSV}",
        "",
        "Done.",
    ])

    with open(config.FINAL_REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    print(f"  -> {config.FINAL_REPORT}")

    # Print summary to console (use 'replace' error handler for terminal)
    summary_text = "\n".join(report_lines)
    try:
        print("\n" + summary_text)
    except UnicodeEncodeError:
        print("\n" + summary_text.encode("ascii", errors="replace").decode("ascii"))

    return dedup_meta


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate final deduplicated output files.")
    parser.add_argument("--remove-edge-copies", action="store_true", default=False,
                       help="Retain only the GenBank copy for edge-case pairs (remove GISAID side)")
    args = parser.parse_args()
    generate_output(remove_edge_copies=args.remove_edge_copies)
