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

def generate_output(remove_categories=frozenset()):
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
        # If GB subtype is missing, inherit from matched GI record
        if gi_row is not None and not rec.get("Subtype"):
            rec["Subtype"] = _v(gi_row.get("Subtype", ""))
        add_source_metadata(rec, gb_row, "GenBank")
        add_source_metadata(rec, gi_row, "GISAID")
        add_source_metadata(rec, match_row, "Match")
        dedup_records.append(rec)

        matched_gisaid_accessions.add(gi_acc)
        matched_genbank_accessions.add(gb_acc)

    # Edge cases: optionally remove GI copies per category
    if remove_categories:
        edges_to_remove = edge_cases[edge_cases["Category"].isin(remove_categories)].copy()
        edges_to_keep   = edge_cases[~edge_cases["Category"].isin(remove_categories)].copy()
    else:
        edges_to_remove = pd.DataFrame()
        edges_to_keep   = edge_cases

    if not edges_to_remove.empty:
        existing_gb = {r["Primary_Accession"] for r in dedup_records}
        for idx, edge_row in edges_to_remove.iterrows():
            gb_acc = edge_row.get("GenBank_Accession", "")
            gi_acc = edge_row.get("GISAID_Accession", "")
            if not gb_acc or not gi_acc:
                continue

            gb_row = gb_by_acc.get(gb_acc)
            if gb_row is None:
                continue
            gi_row = gi_by_acc.get(gi_acc)

            if gb_acc in existing_gb:
                matched_gisaid_accessions.add(gi_acc)
                continue

            rec = extract_metadata_row(gb_row, "GenBank", prefix="")
            rec["GISAID_Accession"] = gi_acc
            rec["Source"] = "GenBank"
            rec["Integration_Status"] = "EDGE"
            rec["Primary_Accession"] = gb_acc
            if gi_row is not None and not rec.get("Lineage"):
                rec["Lineage"] = _v(gi_row.get("GISAID_Lineage", "")) or _v(gi_row.get("Lineage", ""))
            if gi_row is not None and not rec.get("Subtype"):
                rec["Subtype"] = _v(gi_row.get("Subtype", ""))
            add_source_metadata(rec, gb_row, "GenBank")
            add_source_metadata(rec, gi_row, "GISAID")
            add_source_metadata(rec, edge_row, "Edge")
            dedup_records.append(rec)

            matched_gisaid_accessions.add(gi_acc)
            matched_genbank_accessions.add(gb_acc)
            existing_gb.add(gb_acc)

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
    if remove_categories and not edges_to_remove.empty:
        for idx, edge_row in edges_to_remove.iterrows():
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

    # Build lookup: GI accession -> (original_header, sequence)
    gi_seq_lookup = {}
    for record in SeqIO.parse(config.GISAID_FASTA, "fasta"):
        acc = extract_gisaid_accession(record.id)
        if acc and acc not in gi_seq_lookup:
            gi_seq_lookup[acc] = (record.description, str(record.seq))

    # Build lookup: GB accession -> (original_header, sequence)
    gb_seq_lookup = {}
    for record in SeqIO.parse(config.GENBANK_FASTA, "fasta"):
        acc = extract_genbank_accession(record.id)
        if acc and acc not in gb_seq_lookup:
            gb_seq_lookup[acc] = (record.description, str(record.seq))

    removed_fasta_count = 0
    with open(config.FINAL_REMOVED_FASTA, "w") as out_fa:
        for rec in removed_records:
            gi_acc = rec.get("Removed_Accession", "")
            source = rec.get("Removed_Source", "")
            header, seq = None, None
            if source == "GISAID":
                entry = gi_seq_lookup.get(gi_acc)
                if entry:
                    header, seq = entry
            else:
                entry = gb_seq_lookup.get(gi_acc)
                if entry:
                    header, seq = entry
            if header and seq:
                out_fa.write(f">{header}\n{seq}\n")
                removed_fasta_count += 1
            else:
                # Fallback: try the original header from the CSV
                orig_h = rec.get("Removed_Original_Header", "")
                if orig_h:
                    out_fa.write(f">{orig_h}\n{{missing_sequence}}\n")
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
    edge_group_counts = defaultdict(int)
    for _, r in edge_cases.iterrows():
        cat = str(r.get("Category", ""))
        edge_group_counts[cat] += 1

    edge_removed_counts = defaultdict(int)
    if remove_categories:
        for _, r in edges_to_remove.iterrows():
            cat = str(r.get("Category", ""))
            if cat:
                edge_removed_counts[cat] += 1

    gb_input_count = 0
    gb_len_filtered = 0
    gi_input_count = 0
    gi_len_filtered = 0
    try:
        gb_input = pd.read_csv(config.GENBANK_METADATA)
        gb_input_count = len(gb_input)
        gb_len_filtered = int((gb_input["Length"] < config.MIN_SEQUENCE_LENGTH).sum())
    except Exception:
        pass
    try:
        gi_input = pd.read_excel(config.GISAID_METADATA)
        gi_input_count = len(gi_input)
        gi_len_filtered = int((gi_input["Sequence Length"] < config.MIN_SEQUENCE_LENGTH).sum())
    except Exception:
        pass

    gb_removed_count = 0
    if os.path.exists(config.REMOVED_INTRA_GENBANK_CSV):
        with open(config.REMOVED_INTRA_GENBANK_CSV, encoding="utf-8") as f:
            gb_removed_count = sum(1 for _ in f) - 1  # minus header
    gi_removed_count = 0
    if os.path.exists(config.REMOVED_INTRA_GISAID_CSV):
        with open(config.REMOVED_INTRA_GISAID_CSV, encoding="utf-8") as f:
            gi_removed_count = sum(1 for _ in f) - 1

    gb_no_fasta = gb_input_count - len(gb_deduped) - gb_removed_count - gb_len_filtered
    gi_no_fasta = gi_input_count - len(gi_deduped) - gi_removed_count - gi_len_filtered

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
        f"  GenBank length-filtered:  {gb_len_filtered:,}",
        f"  GenBank no FASTA entry:   {gb_no_fasta:,}",
        f"  GISAID kept:              {len(gi_deduped):,}",
        f"  GISAID removed:           {gi_removed_count:,}",
        f"  GISAID length-filtered:   {gi_len_filtered:,}",
        f"  GISAID no FASTA entry:    {gi_no_fasta:,}",
        "",
        "--- Cross-database matching -----------------------------------------",
        f"  Confirmed duplicates:     {len(cross_matches):,}",
        f"  Edge cases (review):      {edge_count:,}",
        "",
        "--- Edge case breakdown ----------------------------------------------",
    ]
    for group_label, count in sorted(edge_group_counts.items(), key=lambda x: -x[1]):
        report_lines.append(f"  {count:>5,}  {group_label}")
    if edge_removed_counts:
        report_lines.append("")
        for group_label, count in sorted(edge_removed_counts.items(), key=lambda x: -x[1]):
            report_lines.append(f"        ({count:,} removed)")
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
    parser.add_argument("--deduped-gb-meta", default=None, help="Deduped GenBank metadata CSV")
    parser.add_argument("--deduped-gi-meta", default=None, help="Deduped GISAID metadata CSV")
    parser.add_argument("--cross-matches", default=None, help="Cross-database matches CSV")
    parser.add_argument("--edge-cases", default=None, help="Edge cases CSV")
    parser.add_argument("--genbank-fasta", default=None, help="Original GenBank FASTA")
    parser.add_argument("--gisaid-fasta", default=None, help="Original GISAID FASTA")
    parser.add_argument("--output-dir", default=None, help="Output directory for final files")
    parser.add_argument("--remove-edge-copies", action="store_true", default=False,
                       help="Shorthand for all --remove-edge-* flags (retain only GenBank copy for all edge cases)")
    parser.add_argument("--remove-edge-isolate", action="store_true", default=False,
                       help="Remove GISAID copy for isolate-mismatch edge cases")
    parser.add_argument("--remove-edge-date", action="store_true", default=False,
                       help="Remove GISAID copy for date-mismatch edge cases")
    parser.add_argument("--remove-edge-country", action="store_true", default=False,
                       help="Remove GISAID copy for country-mismatch edge cases")
    parser.add_argument("--remove-edge-other", action="store_true", default=False,
                       help="Remove GISAID copy for other edge cases")
    parser.add_argument("--min-seq-length", type=int, default=7000,
                        help="Minimum sequence length used for filtering (default: 7000; 0 = no filter)")
    args = parser.parse_args()

    config.MIN_SEQUENCE_LENGTH = args.min_seq_length

    remove_categories = set()
    if args.remove_edge_copies or args.remove_edge_isolate:
        remove_categories.add("Isolate: numeric code")
        remove_categories.add("Isolate: structured name")
    if args.remove_edge_copies or args.remove_edge_date:
        remove_categories.add("Date: mismatch")
    if args.remove_edge_copies or args.remove_edge_country:
        remove_categories.add("Country: mismatch")
    if args.remove_edge_copies or args.remove_edge_other:
        remove_categories.add("Other")
        remove_categories.add("Subtype: missing")

    if args.deduped_gb_meta is not None:
        config.DEDUPED_GENBANK_METADATA = args.deduped_gb_meta
    if args.deduped_gi_meta is not None:
        config.DEDUPED_GISAID_METADATA = args.deduped_gi_meta
    if args.cross_matches is not None:
        config.CROSS_MATCHES_CSV = args.cross_matches
    if args.edge_cases is not None:
        config.EDGE_CASES_CSV = args.edge_cases
    if args.genbank_fasta is not None:
        config.GENBANK_FASTA = args.genbank_fasta
    if args.gisaid_fasta is not None:
        config.GISAID_FASTA = args.gisaid_fasta
    if args.output_dir is not None:
        config.set_output_dir(args.output_dir)

    generate_output(remove_categories=frozenset(remove_categories))
