"""
Step 4: Generate a clean subset of the deduplicated metadata.
Writes a file with core columns in a fixed order plus genuinely
additive source-prefixed columns that carry data.
"""
import argparse
import os
import pandas as pd
import config


def has_data(series) -> bool:
    return series.notna().any() and not (series.astype(str).str.strip().isin(["", "nan", "na", "n/a", "none"]).all())


# Columns that provide novel information not already in the core set.
_EXTRA_CANDIDATES = [
    # Taxonomy
    "GenBank_Organism_Name",
    "GenBank_Species",
    "GenBank_Genus",
    "GenBank_Family",
    # GenBank specifics
    "GenBank_Nuc_Completeness",
    "GenBank_Geo_Location",
    "GenBank_Tissue_Specimen_Source",
    "GenBank_Organization",
    "GenBank_Org_location",
    "GenBank_Publications",
    "GenBank_Molecule_type",
    # GISAID specifics
    "GISAID_Virus_Name",
    "GISAID_Location",
    "GISAID_Sampling_strategy",
    "GISAID_Gender",
    "GISAID_Patient_age",
    "GISAID_Patient_status",
    "GISAID_Vaccination_history",
    "GISAID_Passage",
    "GISAID_Specimen",
    "GISAID_Additional_host_information",
    "GISAID_AA_Substitutions",
    "GISAID_Sequencing_technology",
    "GISAID_Assembly_method",
    "GISAID_Comment",
    # Diagnostics from match / edge files
    "Match_Match_Reasons",
    "Edge_Category",
    "Edge_Reason",
]


def main():
    parser = argparse.ArgumentParser(description="Generate clean deduplicated output.")
    parser.add_argument("--output-dir", default=None, help="Output directory")
    parser.add_argument("--output-file", default=None,
                        help="Output filename (default: clean_deduplicated_metadata.csv)")
    args = parser.parse_args()

    if args.output_dir is not None:
        config.set_output_dir(args.output_dir)
    output_path = os.path.join(config.OUTPUT_DIR, args.output_file or "clean_deduplicated_metadata.csv")

    print("=" * 60)
    print("CLEAN OUTPUT GENERATION")
    print("=" * 60)

    meta = pd.read_csv(config.FINAL_DEDUP_METADATA, low_memory=False)
    print(f"  Loaded {len(meta):,} records, {len(meta.columns)} columns")

    core = [
        "Primary_Accession",
        "GenBank_Accession",
        "GISAID_Accession",
        "Source",
        "Integration_Status",
        "Subtype",
        "GenBank_Isolate",
        "Country",
        "Region",
        "Collection_Date",
        "Length",
        "Host",
        "GenBank_Genotype",
        "Lineage",
        "Release_Date",
        "Submitters",
    ]

    out = pd.DataFrame()
    for col in core:
        out[col] = meta[col] if col in meta.columns else ""

    extra_cols = [c for c in _EXTRA_CANDIDATES if c in meta.columns and has_data(meta[c])]
    for col in extra_cols:
        out[col] = meta[col]

    out.to_csv(output_path, index=False)
    print(f"  Written: {output_path}")
    print(f"  Columns: {len(core)} core + {len(extra_cols)} extra = {len(out.columns)} total")
    print("Done.")


if __name__ == "__main__":
    main()
