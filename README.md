# Bio465Capstone
Repository for BYU Bio 465 Capstone research group

## Notes

**UniProt gene mapping** is fetched live from the UniProt REST API (Step 1). The release version used is printed to the console during Step 1 (`Using UniProt release: ...`).

**AlphaFold structures** are downloaded from AlphaFold DB. The model version is encoded in each downloaded filename (e.g., `AF-P12345-F1-model_v6.cif`). AlphaFold DB v6 covers the full human proteome.

**Input data files** (`PTMD_disease_associated_ptms.tsv`, `TCGA_frequent_mutations.tsv`) are static files downloaded from PTMD and TCGA.

**`ptm_diseases` is pan-cancer:** The `ptm_diseases` column in the output reflects which diseases the PTM site is associated with in PTMD. The nearby TCGA mutations are pan-cancer and were not filtered by cancer type, so a nearby mutation appearing in the output does not imply it co-occurs in the same cancer type as the PTM disease association.

#  Identifying Cancer Driving Mutations that Target Post-Translational Modifications (PTMs)

Created by Matt Banks, Jaden Searle, Tyler Plauche, and Alissa Moulder

Mentored by Dr. Josh Anderson

## Introduction

This is the workflow for our 2026 Senior Bioinformatics Capstone project at Brigham Young University. This project is a collaboration with the Huntsman Cancer Institute. Below are detailed steps for reproducability.

---

## First: Downloading data

We downloaded raw data from two websites, TCGA data from the National Cancer Institute (https://portal.gdc.cancer.gov/) as welll as PTM disrupting mutational data (https://ptmd.biocuckoo.cn/download.php). 

The TCGA data shows the most frequent mutations across all cancer types. This data is called **TCGA_frequent_mutations.tsv** and is located in the data/ directory of this github.

The PTM disrupting data shows all mutations which disrupt a PTM across the genome. This data is called **PTMD_disease_associated_ptms.tsv** and is located in the data/ directory of this github. 

Lastly, make sure to download the idmap (SO I DON"T FORGET WHERE IS THIS LOL). This is an ID conversion  between the uniprot IDs in the TCGA data and the reuglar gene name in the PTMD database. This file will help us merge these two files together.

Feel free to download these two data files directly from this github in the data/ direcytory to begin. 


