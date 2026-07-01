# Dedup — Cross-Database Sequence Deduplication Pipeline

A deterministic, rule-based pipeline that identifies and removes duplicate viral sequences within and between GenBank and GISAID. Uses exact SHA256 hash matching combined with metadata comparison rules.

## Pipeline Stages

```
01_intra_database_dedup.py        02_cross_database_match.py        03_generate_output.py
┌──────────────────────────┐     ┌──────────────────────────┐     ┌──────────────────────────┐
│ GenBank metadata + FASTA │     │ Deduped GB + GI metadata│      │ Cross-db matches         │
│ GISAID metadata + FASTA  │ ──> │ Shared hash pairs        │ ──> │ Deduped metadata + FASTA │
│                          │     │ checked via 6 rules      │     │ Removed sequences        │
│ Output: deduped per-DB   │     │ Output: match/edge/CR    │     │ Report + verification    │
└──────────────────────────┘     └──────────────────────────┘     └──────────────────────────┘
```

## Quick Start

```bash
pip install pandas numpy biopython rapidfuzz tqdm openpyxl
python 01_intra_database_dedup.py
python 02_cross_database_match.py
python 03_generate_output.py
python verify_results.py          # optional: check consistency
```

## Setup

### Dependencies

- Python 3.8+
- `pandas`, `numpy`, `biopython`, `rapidfuzz`, `tqdm`, `openpyxl`

### Input Files

Place these in the pipeline root directory:

| File | Format | Content |
|------|--------|---------|
| `genbank_meta.csv` | CSV | GenBank metadata with columns: Accession, Organism_Name, Isolate, Length, Geo_Location, Host, Collection_Date, Release_Date, Submitters |
| `genbank_sequences.fasta` | FASTA | Headers matching GenBank accessions |
| `gisaid_metadata.xlsx` | Excel | GISAID metadata with columns: Virus name, Subtype, Accession ID, Collection date, Location, Host, Sequence Length |
| `gisaid_sequences.fasta` | FASTA | Headers matching GISAID accessions (with EPI_ISL IDs) |

### Configuration

All settings are in `config.py`:

- **Paths** — input and output file locations
- **`MIN_SEQUENCE_LENGTH`** — minimum sequence length (default: 7000 bp)
- **`COUNTRY_NORMALIZATION`** — maps raw country names to standard forms
- **`HOST_NORMALIZATION`** — maps raw host names (e.g. `Homo sapiens` -> `Human`)
- **`CORE_METADATA_FIELDS`** — fields used for intra-database dedup comparison
- **`ISOLATE_SIMILARITY_THRESHOLD`** — minimum `token_sort_ratio` for shared-ID isolate matching
- **`EXCLUDED_SHARED_TOKENS`** — lab/visit codes excluded from shared-token detection

## How It Works

### Stage 1: Intra-Database Dedup

Sequences are hashed (SHA256) and grouped by hash. Within each hash group, records with identical core metadata (subtype, country, collection date, length) are considered duplicate submissions. The record with the most recent release date is kept; all others are removed.

Records with identical sequences but different metadata are **kept** at this stage — they may represent different samples with coincidentally identical genomes or independent sequencing by different labs.

### Stage 2: Cross-Database Matching

For sequences whose SHA256 hash appears in both databases, a 6-rule deterministic check classifies each pair:

| Rule | Check | Fail Behavior |
|------|-------|-------------|
| 1. Sequence hash | Identical string | (prerequisite) |
| 2. Subtype | Both present and equal | Not a duplicate |
| 3. Country | Both present and equal | Not a duplicate (rejected) |
| 4. Collection date | Match at coarsest common granularity | Edge case (if granularity differs) or rejected (if same granularity) |
| 5. Sequence length | Both present and equal | Not a duplicate |
| 6. Isolate name | Substring match OR shared ID token + similarity >= threshold | Edge case |

#### Isolate Name Matching (`has_shared_id`)

The shared-token detection uses five independent rules:

1. **Shared digit-containing token** (len >= 2, not a year, not in exclusion list) — e.g. `00249` shared between `MN-MDH-RSVB-00249` and `MN-MDH-00249/2023`
2. **Two or more shared alphabetic tokens** (len >= 2, not years) — shared prefix or institution code
3. **Shared numeric substring >= 3 digits** after splitting letter-digit compounds — catches `PV058`/`058`
4. **Token substring containment** — catches `2022JECVB` / `JECVB`
5. **Single shared alphabetic token >= 5 chars** — likely a lab-assigned sample ID

**Numeric codes:** GenBank isolate names consisting of 7+ digits (e.g. UMC Utrecht lab codes) are always classified as edge cases — they have no semantic relationship to GISAID names.

**Exclusion list:** Tokens like `ox02`, `v01`, `d02`–`d11` are excluded from Rule 1 to prevent false matches between different samples from the same institution.

### Stage 3: Output Generation

Assembles the final deduplicated datasets:

- Cross-database matches -> keep GenBank copy, annotate with GISAID accession
- GenBank-only records -> kept as-is
- GISAID-only records -> kept as-is
- Intra-database duplicates -> removed (tracked with kept counterpart)
- Edge cases and rejected pairs -> tracked in separate files (both sides kept)

## Output Files

All written to `output/`:

| File | Description |
|------|-------------|
| `deduplicated_sequences.fasta` | One sequence per unique entry |
| `deduplicated_metadata.csv` | Harmonized metadata with Integration_Status (MATCH/GENBANK_ONLY/GISAID_ONLY) |
| `removed_sequences.fasta` | All removed sequences with original headers |
| `removed_sequences.csv` | Side-by-side removed + kept metadata with removal reason |
| `cross_database_matches.csv` | Confirmed cross-database duplicates (GB kept, GI removed) |
| `edge_cases.csv` | Same hash + subtype, metadata ambiguous — human review required |
| `considered_rejected.csv` | Same hash + subtype, metadata clearly differs — not duplicates |
| `deduped_genbank_metadata.csv` | GenBank records after intra-database dedup |
| `deduped_gisaid_metadata.csv` | GISAID records after intra-database dedup |
| `deduplication_report.txt` | Summary statistics for the run |

### Verification

```bash
python verify_results.py
```

Reads all intermediate and output files independently and reports:
- Input counts, intra-dedup counts, cross-db match counts
- Edge case breakdown by category
- Cross-checks (FASTA vs metadata counts, expected vs actual kept/removed)

## Adapting to a New Pathogen

1. **Subtype parsing** — override `parse_subtype_genbank` in `harmonize_metadata.py` if the default patterns (`virus <subtype>` or `<subtype> virus`) don't match your organism names
2. **Exclusion list** — add lab/institution/visit codes to `EXCLUDED_SHARED_TOKENS` in `config.py` as needed
3. **Normalization maps** — extend `COUNTRY_NORMALIZATION` and `HOST_NORMALIZATION` for your data

## Output Example

```
--- Input -----------------------------------------------------------
  GenBank records:          60,663
  GISAID records:           92,060
--- Intra-database dedup --------------------------------------------
  GenBank kept:             20,584
  GenBank removed:             113
  GISAID kept:              49,123
  GISAID removed:            2,547
--- Cross-database matching -----------------------------------------
  Confirmed duplicates:     13,480
  Edge cases (review):       1,681
--- Edge case breakdown ----------------------------------------------
  1,384  GB isolate is numeric code (7+ digits)
    173  GB isolate is structured name
     53  date differs by 8-30 days
     42  date differs by 1-7 days
     27  date differs by 31-365 days
      2  date differs by >365 days
--- Final deduplicated ----------------------------------------------
  Total unique sequences:   56,227
  Total metadata records:   56,227
  Total removed sequences:  16,140
```

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
