#  Identifying Cancer Driving Mutations that Target Post-Translational Modifications (PTMs)

Created by Matt Banks, Jaden Searle, Tyler Plauche, and Alissa Moulder

Mentored by Dr. Josh Anderson

## Introduction

This is the workflow for our 2026 Senior Bioinformatics Capstone project at Brigham Young University. This project is a collaboration with the Huntsman Cancer Institute. Below are detailed steps for reproducability.

---

## First: Downloading data

We downloaded raw data from two websites, TCGA data from the National Cancer Institute (https://portal.gdc.cancer.gov/) as well as PTM disease associated data (https://ptmd.biocuckoo.cn/download.php).

The TCGA data shows the most frequent mutations across all cancer types. This data is called **TCGA_frequent_mutations.tsv** and is located in the `data/` directory of this repository.

The PTM disrupting data shows PTMs across the genome associated with disease. This data is called **PTMD_disease_associated_ptms.tsv** and is located in the `data/` directory of this repository.

HTP/LTP scores (**htp_ltp_scores.tsv**) were downloaded from [PhosphoSitePlus](https://www.phosphosite.org/). These scores indicate whether a PTM site has been detected by high-throughput (HTP) or low-throughput (LTP) methods, and are used in Step 4 to annotate the final output. This file is also located in the `data/` directory.

Feel free to download these data files directly from this GitHub in the `data/` directory to begin.

## Reproducing the Analysis

### Requirements

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) (handles Python and all dependencies automatically):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Run the pipeline (Steps 1–4)

```bash
uv run main.py
```

This runs all four steps in sequence:

1. **Filter** — merges and filters the PTMD and TCGA datasets
2. **Download structures** — fetches AlphaFold CIF models and PAE files for each protein
3. **Find nearby mutations** — computes 3D distances between PTM sites and nearby cancer mutations
4. **Merge HTP/LTP scores** — annotates the output with PhosphoSitePlus confidence scores

The main output is **`Output/ptm_mutation_proximity_db.tsv`** — a table of PTM sites, their nearby TCGA mutations, 3D distances, and HTP/LTP annotations.

### Generate Figure 3

```bash
uv run scripts/makeFigure3.py
```

Output figures are saved to the **`Output/`** directory as `.png` files.

---

## Notes

**UniProt gene mapping** is fetched live from the UniProt REST API (Step 1). The release version used is printed to the console during Step 1 (`Using UniProt release: ...`).

**AlphaFold structures** are downloaded from AlphaFold DB. The model version is encoded in each downloaded filename (e.g., `AF-P12345-F1-model_v6.cif`). AlphaFold DB v6 covers the full human proteome.

**Input data files** (`PTMD_disease_associated_ptms.tsv`, `TCGA_frequent_mutations.tsv`) are static files downloaded from PTMD and TCGA.

**`ptm_diseases` is pan-cancer:** The `ptm_diseases` column in the output reflects which diseases the PTM site is associated with in PTMD. The nearby TCGA mutations are pan-cancer and were not filtered by cancer type, so a nearby mutation appearing in the output does not imply it co-occurs in the same cancer type as the PTM disease association.
