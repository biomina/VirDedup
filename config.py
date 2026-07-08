"""
Configuration for cross-database deduplication pipeline.
All paths, thresholds, and normalization rules are defined here.
"""
import os

# ── Root directory (where this config file lives) ─────────────────────────
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Input files ───────────────────────────────────────────────────────────
GENBANK_METADATA    = os.path.join(ROOT_DIR, "genbank_meta.csv")
GENBANK_FASTA       = os.path.join(ROOT_DIR, "genbank_sequences.fasta")
GISAID_METADATA     = os.path.join(ROOT_DIR, "gisaid_metadata.xlsx")
GISAID_FASTA        = os.path.join(ROOT_DIR, "gisaid_sequences.fasta")

# ── Output directory ──────────────────────────────────────────────────────
OUTPUT_DIR = os.path.join(ROOT_DIR, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Intra-database dedup outputs ──────────────────────────────────────────
DEDUPED_GENBANK_METADATA = os.path.join(OUTPUT_DIR, "deduped_genbank_metadata.csv")
DEDUPED_GENBANK_FASTA    = os.path.join(OUTPUT_DIR, "deduped_genbank_sequences.fasta")
DEDUPED_GISAID_METADATA  = os.path.join(OUTPUT_DIR, "deduped_gisaid_metadata.csv")
DEDUPED_GISAID_FASTA     = os.path.join(OUTPUT_DIR, "deduped_gisaid_sequences.fasta")

REMOVED_INTRA_GENBANK_CSV  = os.path.join(OUTPUT_DIR, "removed_intra_genbank.csv")
REMOVED_INTRA_GENBANK_FASTA = os.path.join(OUTPUT_DIR, "removed_intra_genbank.fasta")
REMOVED_INTRA_GISAID_CSV   = os.path.join(OUTPUT_DIR, "removed_intra_gisaid.csv")
REMOVED_INTRA_GISAID_FASTA = os.path.join(OUTPUT_DIR, "removed_intra_gisaid.fasta")

# ── Cross-database match outputs ──────────────────────────────────────────
CROSS_MATCHES_CSV       = os.path.join(OUTPUT_DIR, "cross_database_matches.csv")
EDGE_CASES_CSV          = os.path.join(OUTPUT_DIR, "edge_cases.csv")
CONSIDERED_REJECTED_CSV = os.path.join(OUTPUT_DIR, "considered_rejected.csv")

# ── Final outputs ─────────────────────────────────────────────────────────
FINAL_DEDUP_FASTA     = os.path.join(OUTPUT_DIR, "deduplicated_sequences.fasta")
FINAL_DEDUP_METADATA  = os.path.join(OUTPUT_DIR, "deduplicated_metadata.csv")
FINAL_REMOVED_FASTA   = os.path.join(OUTPUT_DIR, "removed_sequences.fasta")
FINAL_REMOVED_CSV     = os.path.join(OUTPUT_DIR, "removed_sequences.csv")
FINAL_REPORT          = os.path.join(OUTPUT_DIR, "deduplication_report.txt")


# ── Path helpers for CLI overrides ────────────────────────────────────────
def set_input_dir(path):
    """Recompute all input paths after overriding the input directory."""
    global GENBANK_METADATA, GENBANK_FASTA, GISAID_METADATA, GISAID_FASTA
    GENBANK_METADATA = os.path.join(path, "genbank_meta.csv")
    GENBANK_FASTA    = os.path.join(path, "genbank_sequences.fasta")
    GISAID_METADATA  = os.path.join(path, "gisaid_metadata.xlsx")
    GISAID_FASTA     = os.path.join(path, "gisaid_sequences.fasta")


def set_output_dir(path):
    """Recompute all output paths after overriding the output directory."""
    global OUTPUT_DIR
    global DEDUPED_GENBANK_METADATA, DEDUPED_GENBANK_FASTA
    global DEDUPED_GISAID_METADATA, DEDUPED_GISAID_FASTA
    global REMOVED_INTRA_GENBANK_CSV, REMOVED_INTRA_GENBANK_FASTA
    global REMOVED_INTRA_GISAID_CSV, REMOVED_INTRA_GISAID_FASTA
    global CROSS_MATCHES_CSV, EDGE_CASES_CSV, CONSIDERED_REJECTED_CSV
    global FINAL_DEDUP_FASTA, FINAL_DEDUP_METADATA
    global FINAL_REMOVED_FASTA, FINAL_REMOVED_CSV, FINAL_REPORT
    OUTPUT_DIR = path
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    DEDUPED_GENBANK_METADATA   = os.path.join(OUTPUT_DIR, "deduped_genbank_metadata.csv")
    DEDUPED_GENBANK_FASTA      = os.path.join(OUTPUT_DIR, "deduped_genbank_sequences.fasta")
    DEDUPED_GISAID_METADATA    = os.path.join(OUTPUT_DIR, "deduped_gisaid_metadata.csv")
    DEDUPED_GISAID_FASTA       = os.path.join(OUTPUT_DIR, "deduped_gisaid_sequences.fasta")
    REMOVED_INTRA_GENBANK_CSV  = os.path.join(OUTPUT_DIR, "removed_intra_genbank.csv")
    REMOVED_INTRA_GENBANK_FASTA = os.path.join(OUTPUT_DIR, "removed_intra_genbank.fasta")
    REMOVED_INTRA_GISAID_CSV   = os.path.join(OUTPUT_DIR, "removed_intra_gisaid.csv")
    REMOVED_INTRA_GISAID_FASTA = os.path.join(OUTPUT_DIR, "removed_intra_gisaid.fasta")
    CROSS_MATCHES_CSV          = os.path.join(OUTPUT_DIR, "cross_database_matches.csv")
    EDGE_CASES_CSV             = os.path.join(OUTPUT_DIR, "edge_cases.csv")
    CONSIDERED_REJECTED_CSV    = os.path.join(OUTPUT_DIR, "considered_rejected.csv")
    FINAL_DEDUP_FASTA          = os.path.join(OUTPUT_DIR, "deduplicated_sequences.fasta")
    FINAL_DEDUP_METADATA       = os.path.join(OUTPUT_DIR, "deduplicated_metadata.csv")
    FINAL_REMOVED_FASTA        = os.path.join(OUTPUT_DIR, "removed_sequences.fasta")
    FINAL_REMOVED_CSV          = os.path.join(OUTPUT_DIR, "removed_sequences.csv")
    FINAL_REPORT               = os.path.join(OUTPUT_DIR, "deduplication_report.txt")

# ── Matching rules ────────────────────────────────────────────────────────
MIN_SEQUENCE_LENGTH = 7000  # only process near-complete genomes



# Country name standardization: {raw: standard}
COUNTRY_NORMALIZATION = {
    "Cape Verde": "Cabo Verde",
    "Democratic Republic of the Congo": "Democratic Republic of the Congo",
    "La Reunion": "Reunion",
    "US Virgin Islands": "United States",
    "USA": "United States",
    "U.S.A.": "United States",
    "U.S.": "United States",
    "United States of America": "United States",
    "Viet Nam": "Vietnam",
    "UK": "United Kingdom",
    "The Netherlands": "Netherlands",
    "South Korea": "Korea",
}

# Host name normalization
HOST_NORMALIZATION = {
    "Homo sapiens": "Human",
}

# The core fields used to determine "same metadata" for intra-database dedup
CORE_METADATA_FIELDS = ["Subtype", "Country", "Collection_Date_Norm", "Length"]

# Isolate name similarity threshold (rapidfuzz token_sort_ratio).
# partial_ratio == 100 (one name is a substring of the other) also counts as
# a match regardless of token_sort_ratio score, e.g. "b2" in "PA-WU-b2/2015".
# Used together with has_shared_id (shared sample-identifying token) to avoid
# false matches from coincidental token overlaps.
# Pure numeric GB isolates (7+ digits, UMC lab codes) are always EDGE_CASE.
ISOLATE_SIMILARITY_THRESHOLD = 30

# Tokens excluded from has_shared_id Rule 1 (shared digit-containing tokens).
# These are lab/institutional codes and visit/timepoint codes that create
# false cross-sample matches (e.g. ox02 pairs different Oxford samples).
EXCLUDED_SHARED_TOKENS = {
    # Visit codes (Oxford, Utrecht, etc.)
    "v01", "v02", "v03",
    "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11",
    # Lab/institution codes
    "ox02", "ox03", "uu02", "uu03", "se02", "im02",
    "se01a", "uu01a",
}
