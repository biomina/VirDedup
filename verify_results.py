"""
Verify the deduplication pipeline results independently.
Reads the intermediate and output files and prints a summary.
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config


def count_fasta_records(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return sum(1 for line in f if line.startswith(">"))


def main():
    print("=" * 60)
    print("DEDUPLICATION RESULT VERIFICATION")
    print("=" * 60)

    # ── 1. Input ────────────────────────────────────────────────────────────
    print("\n--- INPUT ---")
    gb_raw = count_fasta_records(config.GENBANK_FASTA)
    gi_raw = count_fasta_records(config.GISAID_FASTA)
    print(f"  GenBank FASTA entries:     {gb_raw:>8,}")
    print(f"  GISAID FASTA entries:      {gi_raw:>8,}")
    print(f"  Total input:               {gb_raw + gi_raw:>8,}")

    # ── 2. Intra-database dedup ─────────────────────────────────────────────
    print("\n--- INTRA-DATABASE DEDUP ---")
    gb_deduped = len(pd.read_csv(config.DEDUPED_GENBANK_METADATA, low_memory=False))
    gi_deduped = len(pd.read_csv(config.DEDUPED_GISAID_METADATA, low_memory=False))
    gb_removed = len(pd.read_csv(config.REMOVED_INTRA_GENBANK_CSV, low_memory=False))
    gi_removed = len(pd.read_csv(config.REMOVED_INTRA_GISAID_CSV, low_memory=False))
    print(f"  GenBank kept:              {gb_deduped:>8,}")
    print(f"  GenBank removed:           {gb_removed:>8,}")
    print(f"  GISAID kept:               {gi_deduped:>8,}")
    print(f"  GISAID removed:            {gi_removed:>8,}")
    print(f"  Total retained (intra):    {gb_deduped + gi_deduped:>8,}")
    print(f"  Total removed (intra):     {gb_removed + gi_removed:>8,}")

    # ── 3. Cross-database matching ──────────────────────────────────────────
    print("\n--- CROSS-DATABASE MATCHING ---")
    # All unique SeqHashes from deduped metadata
    gb_hashes = set(pd.read_csv(config.DEDUPED_GENBANK_METADATA, low_memory=False)["SeqHash"].astype(str))
    gi_hashes = set(pd.read_csv(config.DEDUPED_GISAID_METADATA, low_memory=False)["SeqHash"].astype(str))
    shared_h = gb_hashes & gi_hashes
    gb_only_h = gb_hashes - gi_hashes
    gi_only_h = gi_hashes - gb_hashes

    matches = pd.read_csv(config.CROSS_MATCHES_CSV, low_memory=False)
    edge_cases = pd.read_csv(config.EDGE_CASES_CSV, low_memory=False)
    considered_rej = pd.read_csv(config.CONSIDERED_REJECTED_CSV, low_memory=False)

    unique_match_hashes = set(matches["SeqHash"].astype(str))
    unique_edge_hashes = set(edge_cases["SeqHash"].astype(str))
    unique_cr_hashes = set(considered_rej["SeqHash"].astype(str))
    all_cross_hashes = unique_match_hashes | unique_edge_hashes | unique_cr_hashes

    print(f"  Shared sequence hashes:    {len(shared_h):>8,}")
    print(f"  GenBank-only hashes:       {len(gb_only_h):>8,}")
    print(f"  GISAID-only hashes:        {len(gi_only_h):>8,}")
    print(f"  Total unique hashes:       {len(gb_only_h) + len(gi_only_h) + len(all_cross_hashes):>8,}")
    print()
    print(f"  Confirmed duplicates:      {len(matches):>8,}")
    if not matches.empty:
        print(f"   -> unique hashes:         {len(unique_match_hashes):>8,}")
    print(f"  Edge cases (review):       {len(edge_cases):>8,}")
    if not edge_cases.empty:
        print(f"   -> unique hashes:         {len(unique_edge_hashes):>8,}")
        print(f"   -> by category:")
        for cat, cnt in edge_cases["Category"].value_counts().sort_index().items():
            print(f"       {cat}: {cnt}")
    print(f"  Considered rejected:       {len(considered_rej):>8,}")
    if not considered_rej.empty:
        print(f"   -> unique hashes:         {len(unique_cr_hashes):>8,}")

    # ── 4. GB-only and GI-only records ──────────────────────────────────────
    # Distinguish unique sequence hashes from accession-level records.
    matched_gi_accessions = set(matches["GISAID_Accession"].astype(str)) if not matches.empty else set()
    gi_meta = pd.read_csv(config.DEDUPED_GISAID_METADATA, low_memory=False)
    gi_unmatched = gi_meta[~gi_meta["Accession"].astype(str).isin(matched_gi_accessions)].copy()
    gi_only_seq_count = gi_unmatched["SeqHash"].astype(str).nunique()
    gi_only_record_count = len(gi_unmatched)
    gb_meta = pd.read_csv(config.DEDUPED_GENBANK_METADATA, low_memory=False)
    gb_only_records = gb_meta[gb_meta["SeqHash"].astype(str).isin(gb_only_h)]
    gb_only_seq_count = len(gb_only_h)
    gb_only_record_count = len(gb_only_records)

    print(f"\n  GenBank-only unique seqs:  {gb_only_seq_count:>8,}")
    print(f"  GenBank-only records:      {gb_only_record_count:>8,}")
    print(f"  GISAID-only unique seqs:   {gi_only_seq_count:>8,}")
    print(f"  GISAID-only records:       {gi_only_record_count:>8,}")

    # ── 5. Final output ────────────────────────────────────────────────────
    print("\n--- FINAL OUTPUT ---")
    final_fasta = count_fasta_records(config.FINAL_DEDUP_FASTA)
    final_meta = len(pd.read_csv(config.FINAL_DEDUP_METADATA, low_memory=False))
    removed_fasta = count_fasta_records(config.FINAL_REMOVED_FASTA)
    removed_meta = len(pd.read_csv(config.FINAL_REMOVED_CSV, low_memory=False))
    print(f"  Deduplicated FASTA:        {final_fasta:>8,}")
    print(f"  Deduplicated metadata:     {final_meta:>8,}")
    print(f"  Removed FASTA:             {removed_fasta:>8,}")
    print(f"  Removed metadata:          {removed_meta:>8,}")

    # ── 6. Cross-checks ────────────────────────────────────────────────────
    print("\n--- CROSS-CHECKS ---")
    # 1. Final = GB-only + GI-only + matched GB records + edge case records + considered-rejected records
    #    (edge cases and considered-rejected keep both sides)
    #    deduplicated_metadata.csv should contain one row per kept sequence.
    #    GB copy is kept for matches; for edge cases and considered-rejected, both GB and GI rows exist.
    #    GB-only + GI-only are also kept.
    #
    #    Count: gb_only + gi_only + (match_gb) + (edge_gb + edge_gi) + (cr_gb + cr_gi)
    #    But we can just read the final metadata and FASTAs.
    print(f"  FASTA == metadata:         {'OK' if final_fasta == final_meta else 'MISMATCH!'}")
    print(f"  Removed FASTA == removed:  {'OK' if removed_fasta == removed_meta else 'MISMATCH!'}")

    expected_removed = gb_removed + gi_removed + len(matches)
    print(f"  Expected removed (intra + cross-db matches): {expected_removed:>8,}")
    print(f"  Actual removed in CSV:                       {removed_meta:>8,}")
    print(f"  Removal count check:      {'OK' if expected_removed == removed_meta else 'MISMATCH!'}")

    expected_kept = gb_deduped + gi_deduped - len(matches)
    print(f"  Expected kept (deduped GB+GI minus cross-db matches): {expected_kept:>8,}")
    print(f"  Actual deduplicated:                                  {final_fasta:>8,}")
    print(f"  Kept count check:        {'OK' if expected_kept == final_fasta else 'MISMATCH!'}")

    print("\nDone.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify deduplication pipeline results.")
    parser.add_argument("--input-dir", default=None, help="Input directory (original files)")
    parser.add_argument("--output-dir", default=None, help="Output directory (all intermediate + final files)")
    args = parser.parse_args()

    if args.input_dir is not None:
        config.set_input_dir(args.input_dir)
    if args.output_dir is not None:
        config.set_output_dir(args.output_dir)

    main()
