import pandas as pd

input_file = "PTM data with scores and mapped distances from ptm site .xlsx"
output_file = "ptm_linear_distances.tsv"

df = pd.read_excel(input_file)

COL_PTM = "ptm_pos"
COL_WITHIN = "mutations_within_5_positions"
COL_MORE = "mutations_more_than_5_positions"


def linear_distances(cell, ptm_pos):
    if pd.isna(cell) or str(cell).strip() == "":
        return ""

    ptm_pos = int(ptm_pos)
    distances = []

    for item in str(cell).split(","):
        item = item.strip()

        mut_pos = int(item.split("-")[0])   # take number before dash

        dist = mut_pos - ptm_pos
        distances.append(str(dist))

    return ",".join(distances)


df["within5_linear_distance"] = df.apply(
    lambda r: linear_distances(r[COL_WITHIN], r[COL_PTM]), axis=1
)

df["morethan5_linear_distance"] = df.apply(
    lambda r: linear_distances(r[COL_MORE], r[COL_PTM]), axis=1
)

# Save EVERYTHING
df.to_csv(output_file, sep="\t", index=False)

print("Saved to:", output_file)
