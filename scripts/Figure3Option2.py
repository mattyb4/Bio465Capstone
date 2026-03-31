import pandas as pd
import matplotlib.pyplot as plt
import re


plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11
})
# -----------------------------
# Load data
# -----------------------------
file_path = "ptm_mutation_proximity_db (1).tsv"
df = pd.read_csv(file_path, sep="\t", encoding="utf-16")

# -----------------------------
# Helpers
# -----------------------------
def clean_ptm_type(x):
    if pd.isna(x):
        return "Unknown"

    x = str(x).strip()

    mapping = {
        "serine phosphorylation": "Ser phosphorylation",
        "threonine phosphorylation": "Thr phosphorylation",
        "tyrosine phosphorylation": "Tyr phosphorylation",
        "phosphorylation": "Phosphorylation",
        "ubiquitination": "Ubiquitination",
        "monoubiquitination": "Monoubiquitination",
        "n-linked glycosylation": "N-linked glycosylation",
        "o-linked glycosylation": "O-linked glycosylation",
        "glycosylation": "Glycosylation",
        "methylation": "Methylation",
        "dimethylation": "Dimethylation",
        "acetylation": "Acetylation",
        "sumoylation": "SUMOylation",
        "adp-ribosylation": "ADP-ribosylation",
        "s-nitrosylation": "S-Nitrosylation",
        "deamidation": "Deamidation",
        "citrullination": "Citrullination",
    }

    return mapping.get(x.lower(), x)

def extract_angstroms(cell):
    if pd.isna(cell):
        return []
    text = str(cell)
    matches = re.findall(r'-([0-9]+(?:\.[0-9]+)?)Å', text)
    return [float(x) for x in matches]

def clean_disease(x):
    x = str(x).strip().lower()

    if x in {"", "nan", "none", "na", "n/a"}:
        return None

    if "hereditary" in x or "syndrome" in x or "multiple cancers" in x:
        return None
    if x in {"cancer", "sporadic cancer"}:
        return None
    if "hamartoma tumor syndrome" in x:
        return None
    if "predisposition" in x:
        return None
    if "tumor" in x and all(k not in x for k in [
        "gastrointestinal stromal tumor",
        "brain tumor",
        "testicular tumor",
        "biliary tract tumor",
        "germ cell tumor"
    ]):
        return None

    if "breast" in x:
        return "Breast"
    elif "lung" in x:
        return "Lung"
    elif "colon" in x or "colorectal" in x or "rectum" in x:
        return "Colorectal"
    elif "glioblastoma" in x or "glioma" in x or "brain lower grade glioma" in x or x == "brain tumor":
        return "Brain"
    elif "melanoma" in x:
        return "Melanoma"
    elif "liver" in x or "hepatocellular" in x:
        return "Liver"
    elif "pancreatic" in x:
        return "Pancreatic"
    elif "prostate" in x:
        return "Prostate"
    elif "ovarian" in x:
        return "Ovarian"
    elif "kidney" in x or "renal" in x:
        return "Renal"
    elif "bladder" in x or "urinary bladder" in x or "transitional cell carcinoma" in x:
        return "Bladder"
    elif "stomach" in x or "gastric" in x:
        return "Gastric"
    elif "endometrial" in x or "uterine" in x or "uterus" in x:
        return "Endometrial"
    elif "cervical" in x:
        return "Cervical"
    elif "esophageal" in x:
        return "Esophageal"
    elif "head and neck" in x:
        return "Head & Neck"
    elif "leukemia" in x:
        return "Leukemia"
    elif "lymphoma" in x:
        return "Lymphoma"
    elif "sarcoma" in x or "leiomyosarcoma" in x or "osteosarcoma" in x:
        return "Sarcoma"
    elif "cholangiocarcinoma" in x or "biliary tract tumor" in x:
        return "Biliary"
    elif "thyroid" in x:
        return "Thyroid"
    elif "mesothelioma" in x:
        return "Mesothelioma"
    elif "medulloblastoma" in x:
        return "Medulloblastoma"
    elif "neuroblastoma" in x:
        return "Neuroblastoma"
    elif "gastrointestinal stromal tumor" in x:
        return "GIST"
    elif "pheochromocytoma" in x or "paraganglioma" in x:
        return "Pheo/Paraganglioma"
    elif "testicular" in x or "germ cell" in x:
        return "Testicular/Germ cell"
    elif "adrenocortical" in x:
        return "Adrenocortical"
    elif "skin cancer" in x:
        return "Skin"
    elif "hematologic cancer" in x or "myeloproliferative" in x:
        return "Hematologic"
    else:
        return None

# -----------------------------
# Clean main columns
# -----------------------------
df["ptm_type_clean"] = df["ptm_type"].apply(clean_ptm_type)
df["ptm_residue"] = df["ptm_site"].astype(str).str.extract(r"^([A-Z])")

# -----------------------------
# Build PTM type x cancer type table
# -----------------------------
heatmap_rows = []

for _, row in df.iterrows():
    ptm_type = row["ptm_type_clean"]
    disease_cell = row.get("ptm_diseases", None)

    if pd.isna(disease_cell):
        continue

    disease_parts = re.split(r"[;|,]", str(disease_cell))

    for disease in disease_parts:
        disease_clean = clean_disease(disease)
        if disease_clean is not None:
            heatmap_rows.append({
                "ptm_type": ptm_type,
                "cancer_type": disease_clean
            })

heatmap_df = pd.DataFrame(heatmap_rows)

if heatmap_df.empty:
    raise ValueError("No usable cancer-type data found in ptm_diseases.")

ptm_cancer_counts = pd.crosstab(
    heatmap_df["ptm_type"],
    heatmap_df["cancer_type"]
)

top_cancers = heatmap_df["cancer_type"].value_counts().head(10).index
ptm_cancer_counts = ptm_cancer_counts.loc[:, top_cancers]

ptm_cancer_counts = ptm_cancer_counts.loc[
    :, ptm_cancer_counts.sum(axis=0).sort_values(ascending=False).index
]

ptm_totals = ptm_cancer_counts.sum(axis=1)
ptm_cancer_counts = ptm_cancer_counts.loc[ptm_totals >= 10]

ptm_cancer_counts = ptm_cancer_counts.loc[
    ptm_cancer_counts.sum(axis=1).sort_values(ascending=False).index
]

ptm_cancer_table = ptm_cancer_counts.div(
    ptm_cancer_counts.sum(axis=1), axis=0
).fillna(0)

# -----------------------------
# Other summaries for panels C and D
# -----------------------------
far_distances = []
for val in df["morethan5_linear_distance"].dropna():
    parts = str(val).split(",")
    for p in parts:
        p = p.strip()
        if p:
            try:
                far_distances.append(float(p))
            except ValueError:
                pass

within5_angstroms = []
morethan5_angstroms = []

for val in df["mutations_within_5_positions"].dropna():
    within5_angstroms.extend(extract_angstroms(val))

for val in df["mutations_more_than_5_positions"].dropna():
    morethan5_angstroms.extend(extract_angstroms(val))

# -----------------------------
# Plot
# -----------------------------
fig, axes = plt.subplots(2, 2, figsize=(20, 14))
fig.suptitle("Snapshot of PTM–Mutation Proximity Dataset", fontsize=18, y=0.99)

# ---- Panel A: PTM type distribution ----
ax = axes[0, 0]
ptm_type_counts = df["ptm_type_clean"].value_counts().sort_values(ascending=True)
bars = ax.barh(ptm_type_counts.index, ptm_type_counts.values)
ax.set_title("A. PTM Type Distribution", fontsize=13)
ax.set_xlabel("Count")
ax.set_ylabel("PTM Type")

# FIX 1: extend x-axis so value labels don't get clipped
max_val = ptm_type_counts.values.max()
ax.set_xlim(0, max_val * 1.15)

for i, v in enumerate(ptm_type_counts.values):
    ax.text(v + max_val * 0.01, i, str(v), va="center", fontsize=9)

# ---- Panel B: PTM type x cancer type heatmap ----
ax = axes[0, 1]
im = ax.imshow(ptm_cancer_table.values, aspect="auto")

ax.set_title("B. PTM Type × Cancer Type", fontsize=13)
ax.set_xlabel("Cancer Type", labelpad=60)  # extra pad to make room for rotated labels
ax.set_ylabel("PTM Type")

ax.set_xticks(range(len(ptm_cancer_table.columns)))
# FIX 2: use horizontal labels instead of rotated to avoid cutoff
ax.set_xticklabels(
    ptm_cancer_table.columns,
    rotation=45,
    ha="right",
    rotation_mode="anchor",
    fontsize=9
)
ax.set_yticks(range(len(ptm_cancer_table.index)))
ax.set_yticklabels(ptm_cancer_table.index)

for i in range(ptm_cancer_table.shape[0]):
    for j in range(ptm_cancer_table.shape[1]):
        value = ptm_cancer_table.iloc[i, j]
        if value > 0:
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8)

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Fraction of PTMs")

# ---- Panel C: far-away mutation linear distances ----
ax = axes[1, 0]
if far_distances:
    ax.hist(far_distances, bins=30, edgecolor="black")
    ax.set_title("C. Linear Distance of Far-Away Mutations (>5 aa)", fontsize=13)
    ax.set_xlabel("Linear distance from PTM site (amino acids)")
    ax.set_ylabel("Number of mutations")
else:
    ax.text(0.5, 0.5, "No far-away mutation distances found",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_title("C. Linear Distance of Far-Away Mutations (>5 aa)", fontsize=13)

# ---- Panel D: 3D distances ----
ax = axes[1, 1]
data_to_plot = []
labels = []

if within5_angstroms:
    data_to_plot.append(within5_angstroms)
    labels.append("Within ±5 aa")

if morethan5_angstroms:
    data_to_plot.append(morethan5_angstroms)
    labels.append(">5 aa away")

if data_to_plot:
    ax.boxplot(data_to_plot, tick_labels=labels, showfliers=False)
    ax.set_title("D. 3D Distance from PTM Site", fontsize=13, pad=10)
    ax.set_ylabel("Distance (Å)")
else:
    ax.text(0.5, 0.5, "No 3D distance data found",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_title("D. 3D Distance from PTM Site", fontsize=13, pad=10)

# FIX 3: improved spacing — more bottom margin for Panel B x-labels, more room overall
plt.tight_layout(rect=[0, 0, 1, 0.97])  # rect leaves room for suptitle

plt.savefig("ptm_dataset_snapshot_cancer_focus.png", dpi=300, bbox_inches="tight")
plt.show()