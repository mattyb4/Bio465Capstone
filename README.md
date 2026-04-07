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

HTP/LTP scores (**htp_ltp_scores.tsv**) were gathered from [PhosphoSitePlus](https://www.phosphosite.org/). These scores indicate whether a PTM site has been detected by high-throughput (HTP) or low-throughput (LTP) methods, and are used in Step 4 to annotate the final output. This file is also located in the `data/` directory.

Feel free to download these data files directly from this GitHub in the `data/` directory to begin.

## Reproducing the Analysis

### Clone the Repository

First, clone this repository to your local machine and navigate into the project directory:

```bash
git clone https://github.com/mattyb4/Bio465Capstone.git
cd Bio465Capstone
```

### Requirements

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) (handles Python and all dependencies automatically):

**macOS/Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

---
### Run the pipeline (Steps 1–4)

```bash
uv run main.py
```

*If you get the error "uv: command not found", see troubleshooting steps below.

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
## Interpereting the Data

### Output Database
The main output of this pipeline is ptm_mutation_proximity_db.tsv, found in the Output folder. This tsv file has the following columns:  

**UniProt** - the UniProt ID  
**gene** - the gene the protein is associated with  
**ptm_site** - position within protein sequence where PTM is  
**ptm_type** - the type of PTM  
**mutations_within_5_positions** - list of all mutation hotspots within 5 residues of PTM site. The formatting is initial amino acid, location, AA it mutates to - distance from PTM site. Then the PAE score* is in parentheses.  
**mutation_count_within_5_positions** - sum of total mutation hotspots in previous column  
**unique_mutation_position_count_within_5_positions** - sum of mutations in unique positions from mutations_within_5_positions  
**mutations_more_than_5_positions** - list of all mutation hotsposts further than 5 residues of PTM site  
**mutation_count_more_than_5_positions** - sum of total mutation hostspots in previous column  
**unique_mutation_position_count_more_than_5_positions** - sum of mutations in unique positions from mutations_more_than_5_positions  
**morethan5_linear_distance** - list of distances on linear amino acid sequence for all mutation hotspots in mutations_more_than_5_positions. This allows for easily seeing entries with mutations that are far on the linear sequence but fold close to PTM site in 3D space  
**mutation_at_ptm_site** - indicates if the PTM site itself is a mutation hotspot  
**ptm_diseases** - lists diseases PTM is associated with according to PTMD 2.0  
**LTP_score** - gives LTP score of that PTM manually retrieved from PhosphoSite  
**HTP_score** - gives HTP score of that PTM manually retrieved from PhosphoSite  



*Predicted Alignment Error (PAE) score is how confident AlphaFold is that those residues are at that position. Lower score = higher confidence

## Error logging
The pipeline also generates logs found in Output/logs to record any issues where the pipeline was unable to download a file for a certain protein from AlphaFold or unable to run calculations for a PTM and why. For more information, see skipped_ptm_summary.md in Output/logs 

## 


---

### Troubleshooting: `uv: command not found`

**macOS/Linux:** After installing, your shell session needs to reload its PATH. Run:

```bash
source "$HOME/.local/bin/env"
```

Then open a new terminal and `uv` should work. If you use conda, ensure `~/.local/bin` is on your PATH by adding this to your shell profile (e.g. `~/.zshrc` or `~/.bash_profile`) and restarting your terminal:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

**Windows:** After installing, close and reopen PowerShell. If `uv` is still not found, add it to your PATH manually:
1. Search for **"Edit the system environment variables"** in the Start menu.
2. Under **User variables**, select `Path` and click **Edit**.
3. Add `%USERPROFILE%\.local\bin`.
4. Click OK and reopen your terminal.

---

## Notes

**UniProt gene mapping** is fetched live from the UniProt REST API (Step 1). The release version used is printed to the console during Step 1 (`Using UniProt release: ...`).

**AlphaFold structures** are downloaded from AlphaFold DB. The model version is encoded in each downloaded filename (e.g., `AF-P12345-F1-model_v6.cif`). AlphaFold DB v6 covers the full human proteome.

**Input data files** (`PTMD_disease_associated_ptms.tsv`, `TCGA_frequent_mutations.tsv`) are static files downloaded from PTMD and TCGA.

**`ptm_diseases` is pan-cancer:** The `ptm_diseases` column in the output reflects which diseases the PTM site is associated with in PTMD. The nearby TCGA mutations are pan-cancer and were not filtered by cancer type, so a nearby mutation appearing in the output does not imply it co-occurs in the same cancer type as the PTM disease association.

**404 / Isoforms Only:** Proteins without available AlphaFold structures or lacking canonical models were excluded from structural analysis.
