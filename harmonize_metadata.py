"""
Metadata harmonization utilities for cross-database deduplication.
Normalizes fields from GenBank and GISAID into a common schema.
"""
import re
import hashlib
from datetime import datetime

import config


# ── Sequence hashing ──────────────────────────────────────────────────────

def hash_sequence(seq_str: str) -> str:
    """SHA256 hash of uppercase, whitespace-stripped nucleotide sequence."""
    clean = seq_str.replace(" ", "").replace("\n", "").upper()
    return hashlib.sha256(clean.encode("ascii", errors="replace")).hexdigest()


# ── Subtype parsing ───────────────────────────────────────────────────────

def parse_subtype_genbank(organism_name) -> str:
    """
    Extract subtype from GenBank Organism_Name.
    Tries common patterns: "virus <subtype>" or "<subtype> virus".
    Override this function for pathogen-specific parsing.
    """
    if not isinstance(organism_name, str):
        return ""
    m = re.search(r"virus\s+(?:type|subtype)\s+([A-Za-z0-9]{1,4})\s*$", organism_name, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.search(r"virus\s+([A-Za-z0-9]{1,4})\s*$", organism_name, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.search(r"\b([A-Za-z0-9]{1,4})\s+virus", organism_name, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    return ""


def parse_subtype_gisaid(subtype_val) -> str:
    """GISAID Subtype column contains the subtype directly."""
    if isinstance(subtype_val, str):
        return subtype_val.strip().upper()
    return ""


# ── Country normalization ─────────────────────────────────────────────────

def normalize_country(raw: str) -> str:
    """Standardize country name; returns empty string if missing."""
    if not isinstance(raw, str) or not raw.strip():
        return ""
    raw = raw.strip()
    raw = raw.replace("_", " ").replace("-", " ").replace(".", " ")
    # Apply normalization map
    for pattern, replacement in config.COUNTRY_NORMALIZATION.items():
        if raw.lower() == pattern.lower():
            return replacement
    # Capitalize each word
    return " ".join(w.capitalize() for w in raw.split() if w)


def parse_location_genbank(geo_location: str):
    """
    Parse GenBank Geo_Location: "Country: Region, City" or just "Country".
    Returns (country, region, city).
    """
    if not isinstance(geo_location, str) or not geo_location.strip():
        return ("", "", "")
    parts = [p.strip() for p in geo_location.split(",")]
    main = parts[0] if parts else ""
    if ":" in main:
        country, region = [x.strip() for x in main.split(":", 1)]
    else:
        country = main
        region = ""
    city = parts[1].strip() if len(parts) > 1 else ""
    return (normalize_country(country), region, city)


def parse_location_gisaid(location: str):
    """
    Parse GISAID Location: "Continent / Country / Region / Subregion".
    Returns (continent, country, region, subregion).
    """
    if not isinstance(location, str) or not location.strip():
        return ("", "", "", "")
    parts = [p.strip() for p in location.split("/")]
    continent = parts[0] if len(parts) > 0 else ""
    country = normalize_country(parts[1]) if len(parts) > 1 else ""
    region = parts[2] if len(parts) > 2 else ""
    subregion = parts[3] if len(parts) > 3 else ""
    return (continent, country, region, subregion)


# ── Host normalization ────────────────────────────────────────────────────

def normalize_host(raw: str) -> str:
    """Normalize host names to a common form."""
    if not isinstance(raw, str) or not raw.strip():
        return ""
    raw = raw.strip()
    for pattern, replacement in config.HOST_NORMALIZATION.items():
        if raw.strip().lower() == pattern.lower():
            return replacement
    return raw


# ── Date normalization ────────────────────────────────────────────────────

def _normalize_iso_like_date(date_str: str) -> str:
    """Normalize YYYY-M-D / YYYY-M / YYYY date strings when valid."""
    date_str = date_str.strip()
    timestamp_match = re.match(
        r"^(\d{4}-\d{1,2}-\d{1,2})[ T]\d{1,2}:\d{2}(?::\d{2}(?:\.\d+)?)?$",
        date_str,
    )
    if timestamp_match:
        date_str = timestamp_match.group(1)

    day_match = re.match(r"^(\d{4})-(\d{1,2})-(\d{1,2})$", date_str)
    if day_match:
        year, month, day = map(int, day_match.groups())
        try:
            dt = datetime(year, month, day)
            return dt.strftime("%Y-%m-%d")
        except ValueError:
            return ""

    month_match = re.match(r"^(\d{4})-(\d{1,2})$", date_str)
    if month_match:
        year, month = month_match.groups()
        month = int(month)
        if 1 <= month <= 12:
            return f"{year}-{month:02d}"
        return ""

    if re.match(r"^\d{4}$", date_str):
        return date_str
    return ""


def normalize_date_genbank(date_str) -> str:
    """
    GenBank dates are typically DD/MM/YYYY or just YYYY.
    Returns YYYY-MM-DD, YYYY-MM, or YYYY (or empty if unparseable).
    """
    if not isinstance(date_str, str):
        return ""
    date_str = date_str.strip().replace(".0", "")
    if not date_str or date_str.lower() in ("nan", "unknown", ""):
        return ""
    # Try DD/MM/YYYY
    try:
        dt = datetime.strptime(date_str, "%d/%m/%Y")
        return dt.strftime("%Y-%m-%d")
    except ValueError:
        pass
    # Try DD/MM/YY
    try:
        dt = datetime.strptime(date_str, "%d/%m/%y")
        return dt.strftime("%Y-%m-%d")
    except ValueError:
        pass
    return _normalize_iso_like_date(date_str)


def normalize_date_gisaid(date_val) -> str:
    """
    GISAID dates are given as datetime objects or YYYY-MM-DD strings.
    Returns YYYY-MM-DD, YYYY-MM, or YYYY (or empty if unparseable).
    """
    if date_val is None:
        return ""
    if isinstance(date_val, datetime):
        return date_val.strftime("%Y-%m-%d")
    date_str = str(date_val).strip()
    if not date_str or date_str.lower() in ("nan", "unknown", "none", ""):
        return ""
    return _normalize_iso_like_date(date_str)


def parse_date_granularity(date_str: str):
    """
    Given a normalized date string, return (year, month, day).
    Missing components are None.
    """
    if not date_str:
        return (None, None, None)
    parts = date_str.split("-")
    year = int(parts[0]) if len(parts) >= 1 and parts[0].isdigit() else None
    month = int(parts[1]) if len(parts) >= 2 and parts[1].isdigit() and parts[1] != "XX" else None
    day = int(parts[2]) if len(parts) >= 3 and parts[2].isdigit() and parts[2] != "XX" else None
    return (year, month, day)


def dates_match_at_best_granularity(d1: str, d2: str) -> bool:
    """
    Compare two normalized date strings at the coarsest common granularity.
    If one has YYYY and the other YYYY-MM-DD, compare only the year.
    If one has YYYY-MM and the other YYYY-MM-DD, compare YYYY-MM.
    Returns True if they match at the available level of detail.
    """
    if not d1 or not d2:
        return True  # missing date = no conflict
    y1, m1, day1 = parse_date_granularity(d1)
    y2, m2, day2 = parse_date_granularity(d2)
    # Year must always match if both have it
    if y1 is not None and y2 is not None and y1 != y2:
        return False
    # If both have month, check month
    if m1 is not None and m2 is not None and m1 != m2:
        return False
    # If both have day, check day
    if day1 is not None and day2 is not None and day1 != day2:
        return False
    return True


# ── Isolate name parsing ─────────────────────────────────────────────────

def parse_isolate_gisaid(virus_name) -> str:
    """
    Extract last 2 components of GISAID virus name.
    "hPVI/A/Spain/AS-HUCA-232576852/2025" -> "AS-HUCA-232576852/2025"
    """
    if not isinstance(virus_name, str) or not virus_name.strip():
        return ""
    components = virus_name.strip().split("/")
    if len(components) >= 2:
        return "/".join(components[-2:])
    return virus_name.strip()


# ── Accession extraction from FASTA headers ──────────────────────────────

def extract_genbank_accession(fasta_id: str) -> str:
    """
    GenBank FASTA header format: ">PZ463377.1" or ">PZ463377"
    Strip version number after '.' and leading '>'.
    """
    acc = fasta_id.lstrip(">")
    return acc.split(".")[0]


def extract_gisaid_accession(fasta_id: str) -> str:
    """
    GISAID FASTA header format: ">hPVI/A/USA/TX-79321/2005|EPI_ISL_15752010|2005"
    Extract the EPI_ISL_XXXXX portion.
    """
    fid = fasta_id.lstrip(">")
    if "|" in fid:
        parts = fid.split("|")
        if len(parts) >= 2:
            return parts[1]
    return fid


# ═══════════════════════════════════════════════════════════════════════════
#  RSV-specific subtype extraction (full fallback chain)
# ═══════════════════════════════════════════════════════════════════════════

def extract_rsv_subtype(meta):
    """
    Full RSV/HRSV subtype extraction fallback chain for GenBank records.
    Steps (in order, first non-empty wins):
      1. Organism_Name via parse_subtype_genbank()
      2. Isolate: H?RSV[A/B] anywhere in string
      3. Isolate: Ken/<digits>/[A/B]/
      4. Isolate: /BA/ pattern
      5. Isolate: non-BA pattern
      6. Isolate: chunk-based standalone A/B or RSVA/RSVB at separator boundaries
      7. Genotype column (raw, uppercased)
      8. Conversion map for non-standard values (BA→B, GA→A, ON1→A, etc.)

    Returns the DataFrame with 'Subtype' column populated.
    """
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
    def _map_subtype(val):
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

    meta["Subtype"] = meta["Subtype"].apply(_map_subtype)
    return meta
