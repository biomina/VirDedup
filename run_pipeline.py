"""
Orchestrator for the VirDedup cross-database deduplication pipeline.

Runs steps 1-3 sequentially via subprocess, passing file paths from
command-line arguments. Every step can also be run independently.

Usage:
    python run_pipeline.py --input-dir ./data --output-dir ./results [--remove-edge-copies]
"""
import argparse
import os
import subprocess
import sys
import time


def run_step(description: str, cmd: list[str]) -> None:
    """Run a pipeline step and exit on failure."""
    print(f"\n{'=' * 60}")
    print(f"  {description}")
    print(f"{'=' * 60}")
    print(f"  Running: {' '.join(cmd)}")
    t0 = time.time()
    result = subprocess.run(cmd)
    elapsed = time.time() - t0
    if result.returncode != 0:
        print(f"\n  ERROR: '{description}' failed (exit code {result.returncode})")
        sys.exit(result.returncode)
    print(f"  Done ({elapsed:.1f}s)")


def main():
    parser = argparse.ArgumentParser(
        description="Run the full VirDedup cross-database deduplication pipeline."
    )
    parser.add_argument("--input-dir", required=True,
                        help="Directory containing input files (genbank_meta.csv, "
                             "genbank_sequences.fasta, gisaid_metadata.xlsx, "
                             "gisaid_sequences.fasta)")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for all intermediate and final output files")
    parser.add_argument("--remove-edge-copies", action="store_true", default=False,
                        help="Retain only the GenBank copy for edge-case pairs")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    python = sys.executable
    script_dir = os.path.dirname(os.path.abspath(__file__))

    os.makedirs(output_dir, exist_ok=True)

    # Build file paths
    gb_meta = os.path.join(input_dir, "genbank_meta.csv")
    gb_fasta = os.path.join(input_dir, "genbank_sequences.fasta")
    gi_meta = os.path.join(input_dir, "gisaid_metadata.xlsx")
    gi_fasta = os.path.join(input_dir, "gisaid_sequences.fasta")

    deduped_gb_meta = os.path.join(output_dir, "deduped_genbank_metadata.csv")
    deduped_gi_meta = os.path.join(output_dir, "deduped_gisaid_metadata.csv")
    cross_matches = os.path.join(output_dir, "cross_database_matches.csv")
    edge_cases = os.path.join(output_dir, "edge_cases.csv")

    # ── Step 1: Intra-database dedup ──────────────────────────────────────
    run_step(
        "Step 1: Intra-database deduplication",
        [python, os.path.join(script_dir, "01_intra_database_dedup.py"),
         "--genbank-meta", gb_meta,
         "--genbank-fasta", gb_fasta,
         "--gisaid-meta", gi_meta,
         "--gisaid-fasta", gi_fasta,
         "--output-dir", output_dir],
    )

    # ── Step 2: Cross-database matching ───────────────────────────────────
    run_step(
        "Step 2: Cross-database matching",
        [python, os.path.join(script_dir, "02_cross_database_match.py"),
         "--deduped-gb-meta", deduped_gb_meta,
         "--deduped-gi-meta", deduped_gi_meta,
         "--output-dir", output_dir],
    )

    # ── Step 3: Generate output ──────────────────────────────────────────
    cmd3 = [
        python, os.path.join(script_dir, "03_generate_output.py"),
        "--deduped-gb-meta", deduped_gb_meta,
        "--deduped-gi-meta", deduped_gi_meta,
        "--cross-matches", cross_matches,
        "--edge-cases", edge_cases,
        "--genbank-fasta", gb_fasta,
        "--gisaid-fasta", gi_fasta,
        "--output-dir", output_dir,
    ]
    if args.remove_edge_copies:
        cmd3.append("--remove-edge-copies")
    run_step("Step 3: Final output generation", cmd3)

    # ── Verification ──────────────────────────────────────────────────────
    run_step(
        "Verification: Checking results",
        [python, os.path.join(script_dir, "verify_results.py"),
         "--input-dir", input_dir,
         "--output-dir", output_dir],
    )

    print(f"\n{'=' * 60}")
    print("  PIPELINE COMPLETE")
    print(f"{'=' * 60}")
    print(f"  Input:  {input_dir}")
    print(f"  Output: {output_dir}")
    print(f"  Edge copies removed: {args.remove_edge_copies}")


if __name__ == "__main__":
    main()
