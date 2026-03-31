import pandas as pd
import matplotlib.pyplot as plt
import re

# -----------------------------
# Load data
# -----------------------------
file_path = "ptm_mutation_proximity_db (1).tsv"
df = pd.read_csv(file_path, sep="\t", encoding="utf-16")

# -----------------------------
# Helpers
# -----------------------------
def normalize_ptm_type(x):
    if pd.isna(x):
        return "Unknown"

    x_lower = str(x).strip().lower()

    if "phosphorylation" in x_lower:
        return "Phosphorylation"
    elif "ubiquitination" in x_lower or "monoubiquitination" in x_lower:
        return "Ubiquitination"
    elif "glycosylation" in x_lower:
        return "Glycosylation"
    elif "methylation" in x_lower or "dimethylation" in x_lower:
        return "Methylation"
    elif "acetylation" in x_lower:
        return "Acetylation"
    elif "sumoylation" in x_lower:
        return "SUMOylation"
    elif "adp-ribosylation" in x_lower:
        return "ADP-ribosylation"
    elif "s-nitrosylation" in x_lower:
        return "S-Nitrosylation"
    elif "deamidation" in x_lower:
        return "Deamidation"
    elif "citrullination" in x_lower:
        return "Citrullination"
    else:
        return str(x).strip()

def extract_angstroms(cell):
    """
    Extract all 3D distances in Å from strings like:
    'S516L-0.00Å(PAE:0.0), R263H-7.04Å(PAE:7.0)'
    """
    if pd.isna(cell):
        return []
    text = str(cell)
    matches = re.findall(r'-([0-9]+(?:\.[0-9]+)?)Å', text)
    return [float(x) for x in matches]

# -----------------------------
# Clean / derive useful columns
# -----------------------------
df["ptm_residue"] = df["ptm_site"].astype(str).str.extract(r"^([A-Z])")
df["ptm_type_clean"] = df["ptm_type"].apply(normalize_ptm_type)

# -----------------------------
# Summaries
# -----------------------------
ptm_type_counts = df["ptm_type_clean"].value_counts()
residue_counts = df["ptm_residue"].value_counts()

preferred_order = ["S", "T", "Y", "K", "R", "N", "Q", "C", "H"]
residue_counts = residue_counts.reindex(
    [r for r in preferred_order if r in residue_counts.index] +
    [r for r in residue_counts.index if r not in preferred_order]
)

# Collect far-away linear distances
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

# Collect 3D distances
within5_angstroms = []
morethan5_angstroms = []

for val in df["mutations_within_5_positions"].dropna():
    within5_angstroms.extend(extract_angstroms(val))

for val in df["mutations_more_than_5_positions"].dropna():
    morethan5_angstroms.extend(extract_angstroms(val))

# -----------------------------
# Plot
# -----------------------------
fig, axes = plt.subplots(2, 2, figsize=(16, 11))
fig.suptitle("Overall Snapshot of PTM–Mutation Proximity Dataset", fontsize=18, y=0.98)

# ---- Panel A: PTM type distribution
ax = axes[0, 0]
top_types = ptm_type_counts.sort_values(ascending=True)
ax.barh(top_types.index, top_types.values)
ax.set_title("A. PTM Type Distribution", fontsize=13)
ax.set_xlabel("Count")
ax.set_ylabel("PTM Type")

for i, v in enumerate(top_types.values):
    ax.text(v + 2, i, str(v), va="center", fontsize=9)

# ---- Panel B: Residue distribution
ax = axes[0, 1]
ax.bar(residue_counts.index, residue_counts.values)
ax.set_title("B. Modified Residue Distribution", fontsize=13)
ax.set_xlabel("Residue")
ax.set_ylabel("Count")

for i, v in enumerate(residue_counts.values):
    ax.text(i, v + 5, str(v), ha="center", fontsize=9)

# ---- Panel C: far-away mutation linear distances
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

# ---- Panel D: 3D distances in angstroms
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
    ax.set_title("D. 3D Distance from PTM Site", fontsize=13)
    ax.set_ylabel("Distance (Å)")
else:
    ax.text(0.5, 0.5, "No 3D distance data found",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_title("D. 3D Distance from PTM Site", fontsize=13)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("ptm_dataset_snapshot.png", dpi=300, bbox_inches="tight")
plt.show()