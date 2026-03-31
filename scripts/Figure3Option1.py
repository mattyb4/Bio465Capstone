import pandas as pd
import matplotlib.pyplot as plt
import re


def extract_ptm_position(ptm_site: str) -> int:
    match = re.search(r"(\d+)", str(ptm_site))
    if not match:
        raise ValueError(f"Could not parse PTM site: {ptm_site}")
    return int(match.group(1))


def parse_mutation_entry(entry: str):
    """
    Parse strings like:
        G60V-8.31Å(PAE:2.1)
        S506P-8.20Å(PAE:4.0)

    Returns:
        mutation_name, mutation_position, distance_3d
    """
    if pd.isna(entry) or str(entry).strip() == "":
        return None

    entry = str(entry).strip()

    # matches amino acid mutation + 3D distance
    m = re.match(r"^([A-Z]\d+[A-Z\*])\-([0-9.]+)Å", entry)
    if not m:
        return None

    mutation_name = m.group(1)
    distance_3d = float(m.group(2))

    pos_match = re.search(r"\d+", mutation_name)
    if not pos_match:
        return None

    mutation_position = int(pos_match.group())

    return mutation_name, mutation_position, distance_3d


def split_mutation_column(cell_value: str):
    """
    Split a mutation column that may contain multiple entries separated by commas or semicolons.
    """
    if pd.isna(cell_value) or str(cell_value).strip() == "":
        return []

    text = str(cell_value).strip()
    parts = re.split(r"[;,]\s*", text)
    return [p for p in parts if p.strip()]


def build_ptm_mutation_table(df: pd.DataFrame, gene: str, ptm_site: str) -> pd.DataFrame:
    """
    Build a long-format mutation table for a single gene and PTM site.
    """
    subset = df[(df["gene"] == gene) & (df["ptm_site"] == ptm_site)].copy()

    if subset.empty:
        raise ValueError(f"No rows found for gene={gene}, ptm_site={ptm_site}")

    rows = []

    for _, row in subset.iterrows():
        ptm_pos = extract_ptm_position(row["ptm_site"])

        # within 5 aa
        for entry in split_mutation_column(row.get("mutations_within_5_positions", "")):
            parsed = parse_mutation_entry(entry)
            if parsed is None:
                continue

            mutation_name, mutation_pos, dist_3d = parsed
            rows.append({
                "UniProt": row["UniProt"],
                "gene": row["gene"],
                "ptm_site": row["ptm_site"],
                "ptm_type": row["ptm_type"],
                "mutation": mutation_name,
                "mutation_position": mutation_pos,
                "ptm_position": ptm_pos,
                "linear_distance": abs(mutation_pos - ptm_pos),
                "distance_3d": dist_3d,
                "group": "Within 5 aa"
            })

        # more than 5 aa
        for entry in split_mutation_column(row.get("mutations_more_than_5_positions", "")):
            parsed = parse_mutation_entry(entry)
            if parsed is None:
                continue

            mutation_name, mutation_pos, dist_3d = parsed
            rows.append({
                "UniProt": row["UniProt"],
                "gene": row["gene"],
                "ptm_site": row["ptm_site"],
                "ptm_type": row["ptm_type"],
                "mutation": mutation_name,
                "mutation_position": mutation_pos,
                "ptm_position": ptm_pos,
                "linear_distance": abs(mutation_pos - ptm_pos),
                "distance_3d": dist_3d,
                "group": "More than 5 aa"
            })

    result = pd.DataFrame(rows)

    if result.empty:
        raise ValueError(f"No parseable mutations found for gene={gene}, ptm_site={ptm_site}")

    result = result.sort_values("mutation_position").reset_index(drop=True)
    return result


def plot_two_panel_ptm_figure(
    mutation_df: pd.DataFrame,
    gene: str,
    ptm_site: str,
    save_path: str | None = None,
    local_window: int = 10,
    structural_cutoff: float = 10.0
):
    """
    Create a two-panel figure:
      Panel A: local sequence context around PTM (unlabeled mutation points)
      Panel B: categorical plot of 3D distances with labels
    """
    ptm_pos = mutation_df["ptm_position"].iloc[0]
    mutation_df = mutation_df.sort_values("mutation_position").reset_index(drop=True)

    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(15, 8),
        height_ratios=[1, 1.4],
        constrained_layout=True
    )

    # -------------------------
    # Panel A: local sequence context
    # -------------------------
    local_min = ptm_pos - local_window
    local_max = ptm_pos + local_window

    local_df = mutation_df[
        (mutation_df["mutation_position"] >= local_min) &
        (mutation_df["mutation_position"] <= local_max)
    ].copy().sort_values("mutation_position")

    # protein segment line
    ax1.hlines(0, local_min, local_max, linewidth=2)

    # shade PTM +/- 5 aa window
    ax1.axvspan(ptm_pos - 5, ptm_pos + 5, alpha=0.15)

    # PTM site
    ax1.scatter(
        ptm_pos,
        0,
        s=180,
        zorder=4,
        edgecolor="black"
    )
    ax1.annotate(
        ptm_site,
        (ptm_pos, 0),
        xytext=(0, 16),
        textcoords="offset points",
        ha="center",
        fontsize=13,
        fontweight="bold"
    )

    # plot mutation points
    for _, row in local_df.iterrows():
        marker = "o" if row["group"] == "Within 5 aa" else "s"
        ax1.scatter(
            row["mutation_position"],
            0,
            s=130,
            marker=marker,
            edgecolor="black",
            zorder=3
        )

    # assign labels to alternating rows
    top_y = 0.32
    bottom_y = -0.32

    # tiny x offsets to reduce collisions for very close residues
    x_nudges = [-0.35, 0.35, -0.2, 0.2, -0.45, 0.45, -0.1, 0.1]

    for i, (_, row) in enumerate(local_df.iterrows()):
        x = row["mutation_position"]
        label = row["mutation"]

        nudge = x_nudges[i % len(x_nudges)]

        if i % 2 == 0:
            label_y = top_y
            va = "bottom"
        else:
            label_y = bottom_y
            va = "top"

        ax1.annotate(
            label,
            xy=(x, 0),
            xytext=(x + nudge, label_y),
            textcoords="data",
            ha="center",
            va=va,
            fontsize=10,
            arrowprops=dict(
                arrowstyle="-",
                lw=0.8,
                shrinkA=0,
                shrinkB=4,
                connectionstyle="angle,angleA=90,angleB=0"
            )
        )

    ax1.set_xlim(local_min, local_max)
    ax1.set_ylim(-0.7, 0.7)
    ax1.set_yticks([])
    ax1.set_xlabel("Protein position near PTM")
    ax1.set_title(f"{gene} {ptm_site}: local sequence context and structural proximity")

    ax1.spines["left"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    # -------------------------
    # Panel B: 3D distance plot
    # -------------------------
    plot_df = mutation_df.copy()
    plot_df["x"] = range(len(plot_df))

    within = plot_df[plot_df["group"] == "Within 5 aa"]
    farther = plot_df[plot_df["group"] == "More than 5 aa"]

    if not within.empty:
        ax2.scatter(
            within["x"],
            within["distance_3d"],
            s=100,
            marker="o",
            label="Within 5 aa"
        )

    if not farther.empty:
        ax2.scatter(
            farther["x"],
            farther["distance_3d"],
            s=100,
            marker="s",
            label="More than 5 aa"
        )

    # optional 3D cutoff line
    ax2.axhline(structural_cutoff, linestyle="--", alpha=0.6, linewidth=1)
    ax2.text(
        len(plot_df) - 0.2,
        structural_cutoff + 0.05,
        f"{structural_cutoff:g} Å",
        ha="right",
        va="bottom",
        fontsize=9
    )

    # label all mutations in panel B
    for i, (_, row) in enumerate(plot_df.iterrows()):
        y_offset = 8 if i % 2 == 0 else 14
        ax2.annotate(
            row["mutation"],
            (row["x"], row["distance_3d"]),
            xytext=(0, y_offset),
            textcoords="offset points",
            ha="center",
            fontsize=9
        )

    ax2.set_xticks(plot_df["x"])
    ax2.set_xticklabels(plot_df["mutation"], rotation=45, ha="right")
    ax2.set_ylabel("3D distance from PTM site (Å)")
    ax2.set_xlabel("Mutation")
    ax2.grid(True, axis="y", alpha=0.3)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.show()

# ----------------------------
# Example usage
# ----------------------------
df = pd.read_csv(r"C:\Users\aliss\Downloads\vsCode\ptm_mutation_proximity_db.csv.csv")

gene_of_interest = "PTPN11"
ptm_of_interest = "Y63"

mutation_df = build_ptm_mutation_table(df, gene_of_interest, ptm_of_interest)
print(mutation_df)


plot_two_panel_ptm_figure(
    mutation_df,
    gene=gene_of_interest,
    ptm_site=ptm_of_interest,
    save_path="PTPN11_Y63_two_panel_figure.png",
    local_window=8,
    structural_cutoff=10.0
)
