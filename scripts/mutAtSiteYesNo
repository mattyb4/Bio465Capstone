import pandas as pd 

input_file= "updated LTPHTP.txt"
output_file="ptm_updated.tsv"

df = pd.read_csv(input_file, sep="\t")

df["Mutation at PTM Site?"] = df["within5_linear_distance"].apply(
    lambda x: "yes" if pd.notna(x) and "0" in str(x).split(",") else "no"
)

df.to_csv(output_file, sep="\t", index=False)