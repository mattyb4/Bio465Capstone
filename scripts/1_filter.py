import re
import requests
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent


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

    site = f"{residue}{position}".strip()
    ptm_type = str(row["Type"]).strip() if pd.notna(row["Type"]) else ""

    return f"{site}:{ptm_type}" if ptm_type else site


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


def fetch_uniprot_gene_mapping(uniprot_ids, batch_size=100):
    """Fetch UniProt accession -> primary gene symbol via the UniProt REST API."""
    # Strip variant suffixes (e.g. Q16613_VAR_A129T -> Q16613) — AlphaFold models canonical sequences
    ids = list({uid.split("_")[0] for uid in set(uniprot_ids)})
    if not ids:
        return pd.DataFrame(columns=["UniProt", "gene"])

    rows = []
    uniprot_release = None
    for i in range(0, len(ids), batch_size):
        batch = ids[i : i + batch_size]
        query = " OR ".join(f"accession:{uid}" for uid in batch)
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": query, "fields": "accession,gene_names", "format": "tsv", "size": batch_size}

        while url:
            resp = requests.get(url, params=params)
            resp.raise_for_status()
            if uniprot_release is None:
                uniprot_release = resp.headers.get("X-UniProt-Release")
            lines = resp.text.strip().split("\n")
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) >= 2:
                    accession = parts[0].strip()
                    primary_gene = parts[1].strip().split()[0] if parts[1].strip() else None
                    if primary_gene:
                        rows.append({"UniProt": accession, "gene": primary_gene})

            link_header = resp.headers.get("Link", "")
            match = re.search(r'<([^>]+)>; rel="next"', link_header)
            url = match.group(1) if match else None
            params = None

    if uniprot_release:
        print(f"Using UniProt release: {uniprot_release}")
    return pd.DataFrame(rows).drop_duplicates(subset=["UniProt"])


def main():
    ptmd_file = PROJECT_ROOT / "data" / "PTMD_disease_associated_ptms.tsv"
    tcga_file = PROJECT_ROOT / "data" / "TCGA_frequent_mutations.tsv"
    output_file = PROJECT_ROOT / "data" / "steps" / "PTMD_TCGA_hotspots_by_protein.tsv"

    # -----------------------
    # Load files
    # -----------------------
    ptmd = pd.read_csv(ptmd_file, sep="\t")
    tcga = pd.read_csv(tcga_file, sep="\t")

    # -----------------------
    # Filter PTMD disruptions
    # -----------------------
    ptmd = ptmd[ptmd["State"] == "N"].copy()

    # Normalize variant UniProt IDs to canonical accession (e.g. Q16613_VAR_A129T -> Q16613)
    ptmd["UniProt"] = ptmd["UniProt"].str.split("_").str[0]

    # -----------------------
    # Map UniProt -> gene via UniProt REST API
    # -----------------------
    uniprot_ids = ptmd["UniProt"].dropna().unique().tolist()
    idmap = fetch_uniprot_gene_mapping(uniprot_ids)

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
            uniprot_id=("UniProt", "first"),
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