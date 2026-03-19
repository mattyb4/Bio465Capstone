"""
match_ptm_scores.py

Matches LTP/HTP scores from ptm_with_scores.tsv onto a new database TSV.
Outputs a new TSV with LTP_score and HTP_score columns appended.

SETUP:
    Place this script in the same folder as build_ptm_score_dict.py
    and ptm_with_scores.tsv.

USAGE:
    python match_ptm_scores.py \\
        --input  new_database.tsv \\
        --gene   GENE_COLUMN_NAME \\
        --ptm    PTM_COLUMN_NAME \\
        --output matched_output.tsv

EXAMPLE:
    python match_ptm_scores.py \\
        --input  my_phospho_db.tsv \\
        --gene   gene_name \\
        --ptm    ptm_site \\
        --output my_phospho_db_with_scores.tsv

NOTES:
    - --gene and --ptm should match the exact column headers in your new database.
    - Rows with no match in ptm_with_scores.tsv will have empty LTP/HTP columns.
    - A summary is printed at the end showing how many rows matched.
"""

import csv
import argparse
from pathlib import Path

from build_ptm_score_dict import build_ptm_dict


def match_scores(input_path, gene_col, ptm_col, output_path, scores_tsv):
    # Load the PTM score dictionary
    ptm_scores = build_ptm_dict(scores_tsv)
    print(f"Loaded {len(ptm_scores)} PTM score entries from '{scores_tsv}'")

    matched   = 0
    unmatched = 0
    total     = 0

    with open(input_path, newline="") as infile, \
         open(output_path, "w", newline="") as outfile:

        reader = csv.DictReader(infile, delimiter="\t")

        # Validate that the specified columns exist
        if gene_col not in reader.fieldnames:
            raise ValueError(
                f"Gene column '{gene_col}' not found in input file.\n"
                f"Available columns: {reader.fieldnames}"
            )
        if ptm_col not in reader.fieldnames:
            raise ValueError(
                f"PTM column '{ptm_col}' not found in input file.\n"
                f"Available columns: {reader.fieldnames}"
            )

        # Add LTP/HTP columns to the output
        out_fields = reader.fieldnames + ["LTP_score", "HTP_score", "uniprot_id_matched"]
        writer = csv.DictWriter(outfile, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()

        for row in reader:
            total += 1
            gene     = row[gene_col].strip()
            ptm_site = row[ptm_col].strip()
            key      = f"{gene}_{ptm_site}"

            if key in ptm_scores:
                entry = ptm_scores[key]
                row["LTP_score"]         = entry["ltp_score"] if entry["ltp_score"] is not None else ""
                row["HTP_score"]         = entry["htp_score"] if entry["htp_score"] is not None else ""
                row["uniprot_id_matched"] = entry["uniprot_id"]
                matched += 1
            else:
                row["LTP_score"]         = ""
                row["HTP_score"]         = ""
                row["uniprot_id_matched"] = ""
                unmatched += 1

            writer.writerow(row)

    print(f"\n--- Matching Summary ---")
    print(f"Total rows processed : {total}")
    print(f"Matched              : {matched}")
    print(f"Unmatched            : {unmatched}")
    print(f"Output written to    : '{output_path}'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add LTP/HTP scores to a new PTM database TSV."
    )
    parser.add_argument("--input",   required=True,  help="Path to the new database TSV file")
    parser.add_argument("--gene",    required=True,  help="Column name for gene names in the new database")
    parser.add_argument("--ptm",     required=True,  help="Column name for PTM sites in the new database")
    parser.add_argument("--output",  required=True,  help="Path for the output TSV file")
    parser.add_argument(
        "--scores",
        default=str(Path(__file__).parent / "ptm_with_scores.tsv"),
        help="Path to ptm_with_scores.tsv (default: same folder as this script)"
    )

    args = parser.parse_args()

    match_scores(
        input_path  = args.input,
        gene_col    = args.gene,
        ptm_col     = args.ptm,
        output_path = args.output,
        scores_tsv  = args.scores,
    )
