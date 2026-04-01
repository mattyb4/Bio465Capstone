"""
build_ptm_score_dict.py

Parses ptm_with_scores.tsv and builds a dictionary keyed by 'GENE_PTMsite'
(e.g. 'SPOP_S119'). Only rows with actual LTP/HTP score values are included
(rows where both scores are empty are skipped). '?' scores are kept as-is.

Each entry stores:
    - uniprot_id  : UniProt accession
    - gene        : gene name
    - ptm_site    : PTM site (e.g. S119)
    - ltp_score   : LTP Score (int, '?', or None if only one score present)
    - htp_score   : HTP Score (int, '?', or None if only one score present)

Usage:
    # Run standalone to preview the dictionary
    python build_ptm_score_dict.py

    # Import into another script
    from build_ptm_score_dict import ptm_scores, build_ptm_dict
"""

import csv
from pathlib import Path


def build_ptm_dict(filepath: str) -> dict:
    """
    Reads the TSV file and returns a dictionary of PTM scores.

    Args:
        filepath: Path to ptm_with_scores.tsv

    Returns:
        dict keyed by 'GENE_PTMsite', e.g.:
        {
            'SPOP_S119': {
                'uniprot_id': 'O43791',
                'gene':       'SPOP',
                'ptm_site':   'S119',
                'ltp_score':  1,
                'htp_score':  0,
            },
            ...
        }
    """
    ptm_dict = {}

    with open(filepath, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            ltp_raw = row["LTP Score"].strip()
            htp_raw = row["HTP Score"].strip()

            # Skip rows where BOTH scores are empty (no data)
            if ltp_raw == "" and htp_raw == "":
                continue

            # Parse score: keep '?' as string, convert digits to int, empty to None
            def parse_score(val):
                if val == "?":
                    return "?"
                if val == "":
                    return None
                return int(val)

            uniprot  = row["UniProt"].strip()
            gene     = row["gene"].strip()
            ptm_site = row["ptm_type"].strip()

            key = f"{gene}_{ptm_site}"

            ptm_dict[key] = {
                "uniprot_id": uniprot,
                "gene":       gene,
                "ptm_site":   ptm_site,
                "ltp_score":  parse_score(ltp_raw),
                "htp_score":  parse_score(htp_raw),
            }

    return ptm_dict


# ---------------------------------------------------------------------------
# Build the dictionary on import so other scripts can do:
#   from build_ptm_score_dict import ptm_scores
# ---------------------------------------------------------------------------
_default_path = Path(__file__).parent / "ptm_with_scores.tsv"
ptm_scores = build_ptm_dict(_default_path)


# ---------------------------------------------------------------------------
# Standalone preview
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print(f"Total PTM entries loaded: {len(ptm_scores)}\n")

    print("--- Sample entries ---")
    for i, (key, val) in enumerate(ptm_scores.items()):
        print(f"{key}: {val}")
        if i >= 9:
            print("...")
            break

    # Example lookup
    example_key = "SPOP_S119"
    if example_key in ptm_scores:
        print(f"\nExample lookup '{example_key}':")
        print(ptm_scores[example_key])
