"""
Step 1: Intra-database deduplication.
Removes duplicate records within GenBank and within GISAID separately.

Logic:
  - Hash each nucleotide sequence (SHA256) for quick comparison.
  - Group records by sequence hash.
  - Within each group, compare core metadata fields:
      Subtype, Country, Collection_Date (at best granularity), Length
  - If all core fields match -> keep the record with the latest release/submission date.
  - If any core field differs -> keep ALL records (different isolates).

Outputs (per database):
  - deduped_*.csv          -- metadata with duplicates removed
  - deduped_*.fasta        -- corresponding sequences
  - removed_intra_*.csv    -- metadata of removed records
  - removed_intra_*.fasta  -- sequences of removed records
"""
import argparse
import csv
import re
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

import config
from harmonize_metadata import (
    hash_sequence,
    parse_subtype_genbank,
    parse_subtype_gisaid,
    parse_location_genbank,
    parse_location_gisaid,
    normalize_host,
    normalize_date_genbank,
    normalize_date_gisaid,
    dates_match_at_best_granularity,
    parse_isolate_gisaid,
    extract_genbank_accession,
    extract_gisaid_accession,
)


# ═══════════════════════════════════════════════════════════════════════════
#  Shared helpers
# ═══════════════════════════════════════════════════════════════════════════

def build_accession_to_hash(fasta_path: str, extract_acc_fn, min_length: int = 0):
    """
    Read a FASTA file and return:
      acc_to_hash: dict {accession: sequence_hash}
      hash_to_headers: dict {hash: [(accession, original_header), ...]}
    Only includes sequences >= min_length.
    """
    acc_to_hash = {}
    hash_to_headers = defaultdict(list)
    for record in tqdm(SeqIO.parse(fasta_path, "fasta"), desc=f"  Hashing {fasta_path}"):
        seq_str = str(record.seq)
        if min_length > 0 and len(seq_str) < min_length:
            continue
        h = hash_sequence(seq_str)
        acc = extract_acc_fn(record.id)
        original_header = record.description  # includes the '>'
        acc_to_hash[acc] = h
        hash_to_headers[h].append((acc, original_header))
    return acc_to_hash, hash_to_headers


def core_fields_match(group_df, subtype_col, country_col, date_col, length_col):
    """
    Check if all records in a group have the same core metadata fields.
    Returns True if all core fields are effectively identical.
    """
    if len(group_df) <= 1:
        return True
    subtypes = group_df[subtype_col].dropna().unique()
    if len([s for s in subtypes if s]) > 1:
        return False
    countries = group_df[country_col].dropna().unique()
    if len([c for c in countries if c]) > 1:
        return False
    # Dates: compare at best granularity for each pair
    dates = group_df[date_col].dropna().tolist()
    non_empty_dates = [d for d in dates if d]
    if len(non_empty_dates) > 1:
        ref = non_empty_dates[0]
        for d in non_empty_dates[1:]:
            if not dates_match_at_best_granularity(ref, d):
                return False
    # Lengths
    lengths = group_df[length_col].dropna().unique()
    if len(lengths) > 1:
        return False
    return True


# ═══════════════════════════════════════════════════════════════════════════
#  GenBank intra-database dedup
# ═══════════════════════════════════════════════════════════════════════════

def dedup_genbank():
    """Deduplicate GenBank records: same sequence + same metadata -> keep newest."""
    print("=" * 60)
    print("GENBANK INTRA-DATABASE DEDUPLICATION")
    print("=" * 60)

    # ── 1. Load metadata ──────────────────────────────────────────────────
    print("\n[1/5] Loading GenBank metadata ...")
    meta = pd.read_csv(config.GENBANK_METADATA, sep=",", low_memory=False)
    print(f"  Records loaded: {len(meta):,}")
    print(f"  Columns: {list(meta.columns)}")

    # Filter by minimum length
    before_len = len(meta)
    meta = meta[meta["Length"] >= config.MIN_SEQUENCE_LENGTH].copy()
    print(f"  After length filter (>={config.MIN_SEQUENCE_LENGTH}bp): {len(meta):,} ({before_len - len(meta)} removed)")

    # ── 2. Compute sequence hashes ───────────────────────────────────────
    print("\n[2/5] Computing GenBank sequence hashes ...")
    acc_to_hash, hash_to_headers = build_accession_to_hash(
        config.GENBANK_FASTA, extract_genbank_accession,
        min_length=config.MIN_SEQUENCE_LENGTH,
    )
    print(f"  Accessions with sequence data: {len(acc_to_hash):,}")
    print(f"  Unique sequences (by hash): {len(hash_to_headers):,}")

    # ── 3. Merge hash into metadata ──────────────────────────────────────
    print("\n[3/5] Merging hashes into metadata ...")
    meta["SeqHash"] = meta["Accession"].map(acc_to_hash)
    before_merge = len(meta)
    meta = meta.dropna(subset=["SeqHash"]).copy()
    print(f"  Records with sequence data: {len(meta):,} ({before_merge - len(meta)} without FASTA entry dropped)")

    # Parse harmonized fields
    print("  Parsing subtype, location, date, host ...")
    meta["Subtype"] = meta["Organism_Name"].apply(parse_subtype_genbank)
    # Fall back to Isolate column for RSV[A/B] and HRSV[A/B] patterns anywhere in string
    mask_empty = meta["Subtype"] == ""
    iso_sub = meta.loc[mask_empty, "Isolate"].str.extract(
        r"(?:H?RSV)([ABab])(?:[^A-Za-z]|$)", expand=False)
    iso_filled = iso_sub.dropna().str.upper()
    meta.loc[iso_filled.index, "Subtype"] = iso_filled
    # Isolate patterns: Ken/X/A/ → A, Ken/X/B/ → B
    mask_empty = meta["Subtype"] == ""
    iso_sub = meta.loc[mask_empty, "Isolate"].str.extract(r"Ken/\d+/([ABab])/", expand=False)
    iso_filled = iso_sub.dropna().str.upper()
    meta.loc[iso_filled.index, "Subtype"] = iso_filled
    # Isolate pattern: /BA/ → B (BA is an RSV-B genotype)
    mask_empty = meta["Subtype"] == ""
    iso_ba = meta.loc[mask_empty, "Isolate"].str.contains(r"(?<=/)BA(?=/)", na=False, regex=True)
    meta.loc[iso_ba[iso_ba].index, "Subtype"] = "B"
    # Isolate pattern: non-BA → non-BA (Kenyan samples)
    mask_empty = meta["Subtype"] == ""
    iso_nonba = meta.loc[mask_empty, "Isolate"].str.contains(r"non-BA", na=False, regex=True)
    meta.loc[iso_nonba[iso_nonba].index, "Subtype"] = "non-BA"
    # Isolate chunk-based subtype extraction: split on separators, check each chunk
    mask_empty = meta["Subtype"] == ""
    iso_series = meta.loc[mask_empty, "Isolate"].fillna("")
    # Standalone A/B chunk or exact RSVA/RSVB chunk (case-insensitive)
    chunk_a = iso_series.str.contains(
        r"(?:^|[/_\-.:,; ])(?:RSV)?A(?:$|[/_\-.:,; ])", na=False, regex=True,
        flags=re.IGNORECASE)
    chunk_b = iso_series.str.contains(
        r"(?:^|[/_\-.:,; ])(?:RSV)?B(?:$|[/_\-.:,; ])", na=False, regex=True,
        flags=re.IGNORECASE)
    # Prefer B over A if both appear in same isolate (takes last assignment)
    meta.loc[mask_empty & chunk_b, "Subtype"] = "B"
    meta.loc[mask_empty & chunk_a & (meta["Subtype"] == ""), "Subtype"] = "A"
    # Flag chunks that start with RSV/HRSV/HRV followed by a non-A/B letter
    flagged = meta.loc[mask_empty, "Isolate"].str.contains(
        r"(?:^|[/_\-.:,; ])(?:H?RSV|HRV)[A-Za-z](?![ABab]|\s|[/_\-.:,;]|$)",
        na=False, regex=True, flags=re.IGNORECASE)
    if flagged.any():
        print(f"  Warning: {flagged.sum()} records have RSV/HRSV/HRV followed by non-A/B "
              "letter in Isolate (subtype not assigned)")
    # Fall back to Genotype column when organism name and isolate yield no subtype
    mask_empty = meta["Subtype"] == ""
    meta.loc[mask_empty, "Subtype"] = meta.loc[mask_empty, "Genotype"].str.strip().str.upper()
    # Standardise non-standard subtype values
    def map_subtype(val):
        if not isinstance(val, str) or not val.strip():
            return val
        s = val.strip().upper()
        exact = {
            "9320": "B", "CB1": "B", "OA1": "A", "ON1": "A", "SAA1": "A",
            "THB": "B", "HRSV-A": "A", "HRSV-B": "B", "RSV A": "A", "RSV B": "B",
            "RSVA": "A", "RSVB": "B", "RSV-B": "B", "RSVA/NA": "A", "RSVA/ON1": "A",
            "RSVB/BA": "B",
        }
        if s in exact:
            return exact[s]
        if s.startswith("BA"):  return "B"
        if s.startswith("GA"):  return "A"
        if s.startswith("GB"):  return "B"
        if s.startswith("NA"):  return "A"
        if s.startswith("SAB"): return "B"
        if s.startswith("ON"):  return "A"
        return val
    meta["Subtype"] = meta["Subtype"].apply(map_subtype)
    loc_data = meta["Geo_Location"].apply(parse_location_genbank)
    meta["Country"] = loc_data.apply(lambda x: x[0])
    meta["Region"] = loc_data.apply(lambda x: x[1])
    meta["Host"] = meta["Host"].apply(normalize_host)
    meta["Collection_Date_Norm"] = meta["Collection_Date"].apply(normalize_date_genbank)

    # Store original Release_Date for "keep newest" comparison
    meta["Release_Date_Norm"] = meta["Release_Date"].apply(normalize_date_genbank)

    # ── 4. Deduplicate by hash ───────────────────────────────────────────
    print("\n[4/5] Deduplicating by sequence hash ...")

    kept_indices = []
    removed_records = []
    removed_sequences = []  # (accession, original_header, seq_str)

    # Build a lookup: accession -> sequence string (for removed sequences)
    accession_to_seq = {}
    for record in SeqIO.parse(config.GENBANK_FASTA, "fasta"):
        acc = extract_genbank_accession(record.id)
        if len(str(record.seq)) >= config.MIN_SEQUENCE_LENGTH:
            accession_to_seq[acc] = str(record.seq)

    # Group by hash
    groups = meta.groupby("SeqHash")
    for seq_hash, group in tqdm(groups, desc="  Processing hash groups", total=len(hash_to_headers)):
        if len(group) == 1:
            kept_indices.append(group.index[0])
        else:
            if core_fields_match(group, "Subtype", "Country", "Collection_Date_Norm", "Length"):
                # Same metadata -> keep newest by Release_Date
                group = group.copy()
                group["_sort_date"] = group["Release_Date_Norm"].fillna("")
                group = group.sort_values("_sort_date", ascending=False)
                keep_idx = group.index[0]
                kept_indices.append(keep_idx)
                # All others are removed
                for idx in group.index[1:]:
                    removed_records.append(group.loc[idx].to_dict())
                    acc = group.loc[idx, "Accession"]
                    seq = accession_to_seq.get(acc, "")
                    # original_header = seq_hash  # fallback; we'll fetch real headers
                    for hdr_acc, hdr in hash_to_headers.get(seq_hash, []):
                        if hdr_acc == acc:
                            removed_sequences.append((acc, hdr, seq))
                            break
            else:
                # Different metadata -> keep all
                for idx in group.index:
                    kept_indices.append(idx)

    deduped = meta.loc[kept_indices].copy()
    print(f"  Records kept: {len(deduped):,}")
    print(f"  Records removed: {len(removed_records):,}")

    # ── 5. Write outputs ─────────────────────────────────────────────────
    print("\n[5/5] Writing output files ...")

    # Deduplicated metadata
    deduped.to_csv(config.DEDUPED_GENBANK_METADATA, index=False, quoting=csv.QUOTE_ALL)
    print(f"  -> {config.DEDUPED_GENBANK_METADATA}")

    # Deduplicated FASTA
    kept_accessions = set(deduped["Accession"])
    with open(config.DEDUPED_GENBANK_FASTA, "w") as out_fa:
        written = 0
        for record in SeqIO.parse(config.GENBANK_FASTA, "fasta"):
            acc = extract_genbank_accession(record.id)
            if acc in kept_accessions and len(str(record.seq)) >= config.MIN_SEQUENCE_LENGTH:
                out_fa.write(f">{record.description}\n{record.seq}\n")
                written += 1
    print(f"  -> {config.DEDUPED_GENBANK_FASTA} ({written:,} sequences)")

    # Removed metadata CSV
    if removed_records:
        removed_df = pd.DataFrame(removed_records)
        removed_df.to_csv(config.REMOVED_INTRA_GENBANK_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.REMOVED_INTRA_GENBANK_CSV} ({len(removed_df):,} records)")

    # Removed FASTA
    if removed_sequences:
        with open(config.REMOVED_INTRA_GENBANK_FASTA, "w") as out_fa:
            for acc, hdr, seq in removed_sequences:
                out_fa.write(f">{hdr}\n{seq}\n")
        print(f"  -> {config.REMOVED_INTRA_GENBANK_FASTA} ({len(removed_sequences):,} sequences)")
    else:
        # Write empty file
        open(config.REMOVED_INTRA_GENBANK_FASTA, "w").close()
        print(f"  -> {config.REMOVED_INTRA_GENBANK_FASTA} (empty)")

    return deduped, kept_accessions, removed_records


# ═══════════════════════════════════════════════════════════════════════════
#  GISAID intra-database dedup
# ═══════════════════════════════════════════════════════════════════════════

def dedup_gisaid():
    """Deduplicate GISAID records: same sequence + same metadata -> keep newest."""
    print("\n" + "=" * 60)
    print("GISAID INTRA-DATABASE DEDUPLICATION")
    print("=" * 60)

    # ── 1. Load metadata ─────────────────────────────────────────────────
    print("\n[1/5] Loading GISAID metadata ...")
    meta = pd.read_excel(config.GISAID_METADATA)
    print(f"  Records loaded: {len(meta):,}")
    print(f"  Columns: {list(meta.columns)}")

    # Rename columns for consistency
    meta.rename(columns={
        "Accession ID":    "Accession",
        "Virus name":      "Virus_Name",
        "Collection date": "Collection_Date",
        "Sequence Length": "Length",
        "Subtype":         "GISAID_Subtype",
        "Lineage":         "GISAID_Lineage",
    }, inplace=True)

    # Filter by minimum length
    before_len = len(meta)
    meta["Length"] = pd.to_numeric(meta["Length"], errors="coerce")
    meta = meta[meta["Length"] >= config.MIN_SEQUENCE_LENGTH].copy()
    print(f"  After length filter (>={config.MIN_SEQUENCE_LENGTH}bp): {len(meta):,} ({before_len - len(meta)} removed)")

    # ── 2. Compute sequence hashes ───────────────────────────────────────
    print("\n[2/5] Computing GISAID sequence hashes ...")
    acc_to_hash, hash_to_headers = build_accession_to_hash(
        config.GISAID_FASTA, extract_gisaid_accession,
        min_length=config.MIN_SEQUENCE_LENGTH,
    )
    print(f"  Accessions with sequence data: {len(acc_to_hash):,}")
    print(f"  Unique sequences (by hash): {len(hash_to_headers):,}")

    # ── 3. Merge hash into metadata ──────────────────────────────────────
    print("\n[3/5] Merging hashes into metadata ...")
    meta["SeqHash"] = meta["Accession"].map(acc_to_hash)
    before_merge = len(meta)
    meta = meta.dropna(subset=["SeqHash"]).copy()
    print(f"  Records with sequence data: {len(meta):,} ({before_merge - len(meta)} without FASTA entry dropped)")

    # Parse harmonized fields
    print("  Parsing subtype, location, date, host ...")
    meta["Subtype"] = meta["GISAID_Subtype"].apply(parse_subtype_gisaid)
    loc_data = meta["Location"].apply(parse_location_gisaid)
    meta["Country"] = loc_data.apply(lambda x: x[1])
    meta["Region"] = loc_data.apply(lambda x: x[2])
    meta["Host"] = meta["Host"].apply(normalize_host)
    meta["Collection_Date_Norm"] = meta["Collection_Date"].apply(normalize_date_gisaid)
    meta["Isolate"] = meta["Virus_Name"].apply(parse_isolate_gisaid)

    # For GISAID we use the EPI_ISL ID itself as the date proxy (no release date field)
    # but we can use Publishing Embargo until or just preserve the first occurrence.
    # GISAID doesn't have a direct "release date" -- we'll keep the first occurrence.
    meta["_gisaid_order"] = range(len(meta))

    # ── 4. Deduplicate by hash ───────────────────────────────────────────
    print("\n[4/5] Deduplicating by sequence hash ...")

    kept_indices = []
    removed_records = []
    removed_sequences = []

    # Build a lookup: accession -> sequence string
    accession_to_seq = {}
    for record in SeqIO.parse(config.GISAID_FASTA, "fasta"):
        acc = extract_gisaid_accession(record.id)
        if len(str(record.seq)) >= config.MIN_SEQUENCE_LENGTH:
            accession_to_seq[acc] = str(record.seq)

    groups = meta.groupby("SeqHash")
    for seq_hash, group in tqdm(groups, desc="  Processing hash groups", total=len(hash_to_headers)):
        if len(group) == 1:
            kept_indices.append(group.index[0])
        else:
            if core_fields_match(group, "Subtype", "Country", "Collection_Date_Norm", "Length"):
                # Same metadata -> keep the first record (no reliable date field)
                group = group.sort_values("_gisaid_order")
                keep_idx = group.index[0]
                kept_indices.append(keep_idx)
                for idx in group.index[1:]:
                    removed_records.append(group.loc[idx].to_dict())
                    acc = group.loc[idx, "Accession"]
                    seq = accession_to_seq.get(acc, "")
                    for hdr_acc, hdr in hash_to_headers.get(seq_hash, []):
                        if hdr_acc == acc:
                            removed_sequences.append((acc, hdr, seq))
                            break
            else:
                for idx in group.index:
                    kept_indices.append(idx)

    deduped = meta.loc[kept_indices].copy()
    print(f"  Records kept: {len(deduped):,}")
    print(f"  Records removed: {len(removed_records):,}")

    # ── 5. Write outputs ─────────────────────────────────────────────────
    print("\n[5/5] Writing output files ...")

    deduped.to_csv(config.DEDUPED_GISAID_METADATA, index=False, quoting=csv.QUOTE_ALL)
    print(f"  -> {config.DEDUPED_GISAID_METADATA}")

    kept_accessions = set(deduped["Accession"])
    with open(config.DEDUPED_GISAID_FASTA, "w") as out_fa:
        written = 0
        for record in SeqIO.parse(config.GISAID_FASTA, "fasta"):
            acc = extract_gisaid_accession(record.id)
            if acc in kept_accessions and len(str(record.seq)) >= config.MIN_SEQUENCE_LENGTH:
                out_fa.write(f">{record.description}\n{record.seq}\n")
                written += 1
    print(f"  -> {config.DEDUPED_GISAID_FASTA} ({written:,} sequences)")

    if removed_records:
        removed_df = pd.DataFrame(removed_records)
        removed_df.to_csv(config.REMOVED_INTRA_GISAID_CSV, index=False, quoting=csv.QUOTE_ALL)
        print(f"  -> {config.REMOVED_INTRA_GISAID_CSV} ({len(removed_df):,} records)")

    if removed_sequences:
        with open(config.REMOVED_INTRA_GISAID_FASTA, "w") as out_fa:
            for acc, hdr, seq in removed_sequences:
                out_fa.write(f">{hdr}\n{seq}\n")
        print(f"  -> {config.REMOVED_INTRA_GISAID_FASTA} ({len(removed_sequences):,} sequences)")
    else:
        open(config.REMOVED_INTRA_GISAID_FASTA, "w").close()
        print(f"  -> {config.REMOVED_INTRA_GISAID_FASTA} (empty)")

    return deduped, kept_accessions, removed_records


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Intra-database deduplication (Step 1).")
    parser.add_argument("--genbank-meta", default=None, help="GenBank metadata CSV path")
    parser.add_argument("--genbank-fasta", default=None, help="GenBank FASTA path")
    parser.add_argument("--gisaid-meta", default=None, help="GISAID metadata Excel path")
    parser.add_argument("--gisaid-fasta", default=None, help="GISAID FASTA path")
    parser.add_argument("--output-dir", default=None, help="Output directory for deduped files")
    parser.add_argument("--min-seq-length", type=int, default=7000,
                        help="Minimum sequence length to process (default: 7000; 0 = no filter)")
    args = parser.parse_args()

    config.MIN_SEQUENCE_LENGTH = args.min_seq_length

    if args.genbank_meta is not None:
        config.GENBANK_METADATA = args.genbank_meta
    if args.genbank_fasta is not None:
        config.GENBANK_FASTA = args.genbank_fasta
    if args.gisaid_meta is not None:
        config.GISAID_METADATA = args.gisaid_meta
    if args.gisaid_fasta is not None:
        config.GISAID_FASTA = args.gisaid_fasta
    if args.output_dir is not None:
        config.set_output_dir(args.output_dir)

    gb_deduped, gb_kept_accs, gb_removed = dedup_genbank()
    gi_deduped, gi_kept_accs, gi_removed = dedup_gisaid()

    print("\n" + "=" * 60)
    print("INTRA-DATABASE DEDUP SUMMARY")
    print("=" * 60)
    print(f"  GenBank: {len(gb_deduped):,} kept, {len(gb_removed):,} removed "
          f"({len(gb_kept_accs):,} unique accessions in output)")
    print(f"  GISAID:  {len(gi_deduped):,} kept, {len(gi_removed):,} removed "
          f"({len(gi_kept_accs):,} unique accessions in output)")
    print("\nDone with intra-database deduplication.")
