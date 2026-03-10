import re
import pandas as pd


def clean_str_list(values):
    cleaned = []
    seen = set()

    for value in values.dropna():
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            cleaned.append(text)

    return "; ".join(cleaned)


def extract_gene_from_protein_change(protein_change):
    if pd.isna(protein_change):
        return None

    text = str(protein_change).strip()
    parts = text.split()

    return parts[0] if parts else None


def extract_aa_change(protein_change):
    if pd.isna(protein_change):
        return None

    text = str(protein_change).strip()
    parts = text.split()

    if len(parts) < 2:
        return None

    return parts[1]


def is_simple_substitution(change):
    if pd.isna(change):
        return False

    text = str(change).strip()
    return bool(re.fullmatch(r"[A-Z]\d+[A-Z]", text))


def build_ptm_site(row):
    residue = str(row["Residue"]).strip() if pd.notna(row["Residue"]) else ""

    position = ""
    if pd.notna(row["Position"]):
        try:
            position = str(int(float(row["Position"])))
        except Exception:
            position = str(row["Position"]).strip()

    return f"{residue}{position}".strip()


def format_mutation_with_count(row):
    return f'{row["mutation"]} ({int(row["affected_cases"])})'


def find_case_count_column(df):
    candidates = [
        "num_cohort_ssm_affected_cases",
        "cohort_ssm_affected_cases",
        "affected_cases",
        "cases",
        "count"
    ]

    for col in candidates:
        if col in df.columns:
            return col

    raise ValueError(
        "Could not find a case-count column in the TCGA file. "
        "Expected one of: num_cohort_ssm_affected_cases, "
        "cohort_ssm_affected_cases, affected_cases, cases, count"
    )


def main():
    ptmd_file = "data/Total_4.tsv"
    tcga_file = "data/frequent-mutations.2026-03-05.tsv"
    idmap_file = "data/idmap.xlsx"

    output_file = "data/PTMD_TCGA_hotspots_by_protein.tsv"

    # -----------------------
    # Load files
    # -----------------------
    ptmd = pd.read_csv(ptmd_file, sep="\t")
    tcga = pd.read_csv(tcga_file, sep="\t")
    idmap = pd.read_excel(idmap_file)

    # -----------------------
    # Filter PTMD disruptions
    # -----------------------
    ptmd = ptmd[ptmd["State"] == "N"].copy()

    # -----------------------
    # Map UniProt -> gene
    # -----------------------
    idmap = idmap[["query", "symbol"]].dropna().drop_duplicates()
    idmap = idmap.rename(columns={"query": "UniProt", "symbol": "gene"})

    ptmd = ptmd.merge(idmap, on="UniProt", how="left")

    if "Gene name" in ptmd.columns:
        ptmd["gene"] = ptmd["gene"].fillna(ptmd["Gene name"])

    ptmd = ptmd[ptmd["gene"].notna()].copy()

    # -----------------------
    # Build PTM site and PTM-disease pair
    # -----------------------
    ptmd["ptm_site"] = ptmd.apply(build_ptm_site, axis=1)
    ptmd["ptm_disease_pair"] = ptmd["ptm_site"] + " | " + ptmd["Disease"].astype(str)

    # -----------------------
    # Prepare TCGA mutations
    # -----------------------
    tcga["gene"] = tcga["protein_change"].apply(extract_gene_from_protein_change)
    tcga["aa_change"] = tcga["protein_change"].apply(extract_aa_change)

    tcga = tcga[tcga["gene"].notna()].copy()
    tcga = tcga[tcga["aa_change"].apply(is_simple_substitution)].copy()

    # Find case count column automatically
    case_count_col = find_case_count_column(tcga)

    tcga["affected_cases"] = pd.to_numeric(tcga[case_count_col], errors="coerce")
    tcga = tcga[tcga["affected_cases"].notna()].copy()

    # Keep only hotspot mutations with affected cases >= 3
    tcga = tcga[tcga["affected_cases"] >= 3].copy()

    # Keep only needed mutation label
    tcga["mutation"] = tcga["aa_change"]
    tcga["mutation_with_count"] = tcga.apply(format_mutation_with_count, axis=1)

    # -----------------------
    # Aggregate PTMs by gene
    # -----------------------
    ptmd_grouped = (
        ptmd.groupby("gene", as_index=False)
        .agg(
            uniprot_id=("UniProt", clean_str_list),
            ptms_on_protein=("ptm_site", clean_str_list),
            ptm_disease_pairs=("ptm_disease_pair", clean_str_list),
        )
    )

    # -----------------------
    # Aggregate hotspot mutations by gene
    # -----------------------
    tcga_grouped = (
        tcga.groupby("gene", as_index=False)
        .agg(
            mutations_on_protein=("mutation_with_count", clean_str_list),
        )
    )

    # -----------------------
    # Merge datasets
    # -----------------------
    merged = ptmd_grouped.merge(tcga_grouped, on="gene", how="left")

    merged = merged[
        [
            "uniprot_id",
            "gene",
            "ptms_on_protein",
            "mutations_on_protein",
            "ptm_disease_pairs",

        ]
    ]

    merged = merged[merged["mutations_on_protein"].notna()].copy()

    # -----------------------
    # Save result
    # -----------------------
    merged.to_csv(output_file, sep="\t", index=False)

    print("Done.")
    print(f"Using case count column: {case_count_col}")
    print(f"PTMD genes kept: {ptmd['gene'].nunique()}")
    print(f"TCGA hotspot genes kept: {tcga['gene'].nunique()}")
    print(f"Proteins with PTMs and mutation hotspots: {len(merged)}")
    print(f"Output saved to: {output_file}")
    print("PTMD disruption genes:", ptmd["gene"].nunique())
    print("TCGA hotspot genes:", tcga["gene"].nunique())
    print("Final merged proteins:", len(merged))


if __name__ == "__main__":
    main()