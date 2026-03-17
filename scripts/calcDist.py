import pandas as pd
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
input_file = PROJECT_ROOT / "Output" / "nearby_mutations_dbv2.0.tsv"
output_file = PROJECT_ROOT / "Output" / "ptm_linear_distances.tsv"

def read_nearby_table(path: Path):
    for encoding in ("utf-16", "utf-8-sig", "utf-8", "latin1"):
        try:
            return pd.read_csv(path, sep="\t", encoding=encoding)
        except UnicodeError:
            continue
    return pd.read_csv(path, sep="\t")


df = read_nearby_table(input_file)

COL_PTM = "ptm_pos"
COL_PTM_TYPE = "ptm_type"
COL_WITHIN = "mutations_within_5_positions"
COL_MORE = "mutations_more_than_5_positions"


def get_ptm_pos(row):
    if COL_PTM in row and pd.notna(row[COL_PTM]):
        return int(row[COL_PTM])

    ptm_type = str(row.get(COL_PTM_TYPE, "")).strip()
    match = re.search(r"\d+", ptm_type)
    if not match:
        raise ValueError(f"Could not parse PTM position from '{ptm_type}'.")
    return int(match.group())


def linear_distances(cell, ptm_pos):
    if pd.isna(cell) or str(cell).strip() == "":
        return ""

    ptm_pos = int(ptm_pos)
    distances = set()

    for item in str(cell).split(","):
        item = item.strip()

        mut_part = item.split("-")[0]  # e.g. A264T
        match = re.search(r"\d+", mut_part)

        if match:
            mut_pos = int(match.group())
            dist = abs(mut_pos - ptm_pos)
            distances.add(dist)

    return ",".join(str(d) for d in sorted(distances))


df["within5_linear_distance"] = df.apply(
    lambda r: linear_distances(r[COL_WITHIN], get_ptm_pos(r)), axis=1
)

df["morethan5_linear_distance"] = df.apply(
    lambda r: linear_distances(r[COL_MORE], get_ptm_pos(r)), axis=1
)

df.to_csv(output_file, sep="\t", index=False)

print("Saved to:", output_file)