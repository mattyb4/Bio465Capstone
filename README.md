# Bio465Capstone
Repository for BYU Bio 465 Capstone research group

## Notes

**UniProt gene mapping** is fetched live from the UniProt REST API (Step 1). The release version used is printed to the console during Step 1 (`Using UniProt release: ...`).

**AlphaFold structures** are downloaded from AlphaFold DB. The model version is encoded in each downloaded filename (e.g., `AF-P12345-F1-model_v6.cif`). AlphaFold DB v6 covers the full human proteome.

**Input data files** (`PTMD_disease_associated_ptms.tsv`, `TCGA_frequent_mutations.tsv`) are static files downloaded from PTMD and TCGA.

**`ptm_diseases` is pan-cancer:** The `ptm_diseases` column in the output reflects which diseases the PTM site is associated with in PTMD. The nearby TCGA mutations are pan-cancer and were not filtered by cancer type, so a nearby mutation appearing in the output does not imply it co-occurs in the same cancer type as the PTM disease association.
