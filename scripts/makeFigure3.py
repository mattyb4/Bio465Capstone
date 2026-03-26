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
    # Panel A: local sequence context with far-mutation inset
    # -------------------------
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    local_min = ptm_pos - local_window
    local_max = max(ptm_pos + local_window, 72)

    local_df = mutation_df[
        (mutation_df["mutation_position"] >= local_min) &
        (mutation_df["mutation_position"] <= local_max)
    ].copy().sort_values("mutation_position")

    far_df = mutation_df[
        (mutation_df["mutation_position"] >= 500) &
        (mutation_df["mutation_position"] <= 510)
    ].copy().sort_values("mutation_position")

    # main local protein segment
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

    # plot local mutation points
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

    def add_staggered_labels(ax, df, proximity=0.7, fontsize=10):
        top_levels = [0.24, 0.34, 0.44]
        bottom_levels = [-0.24, -0.34, -0.44]
        placed = []

        for _, row in df.iterrows():
            x = row["mutation_position"]
            label = row["mutation"]

            nearby = [p for p in placed if abs(p["x"] - x) < proximity]

            if not nearby:
                side = "top"
                level_idx = 0
            else:
                top_count = sum(1 for p in nearby if p["side"] == "top")
                bottom_count = sum(1 for p in nearby if p["side"] == "bottom")

                if top_count <= bottom_count:
                    side = "top"
                    level_idx = top_count % len(top_levels)
                else:
                    side = "bottom"
                    level_idx = bottom_count % len(bottom_levels)

            if side == "top":
                label_y = top_levels[level_idx]
                va = "bottom"
                line_end_y = label_y - 0.02
            else:
                label_y = bottom_levels[level_idx]
                va = "top"
                line_end_y = label_y + 0.02

            ax.plot(
                [x, x],
                [0, line_end_y],
                color="black",
                lw=0.8,
                zorder=2,
                clip_on=False
            )

            ax.text(
                x,
                label_y,
                label,
                ha="center",
                va=va,
                fontsize=fontsize,
                clip_on=False
            )

            placed.append({"x": x, "side": side})

    # stagger labels in local panel
    add_staggered_labels(ax1, local_df, proximity=0.7, fontsize=10)

    # style main local axis
    ax1.set_xlim(local_min - 0.8, local_max + 0.8)
    ax1.set_ylim(-0.5, 0.5)
    ax1.set_yticks([])
    ax1.set_ylabel("")
    ax1.spines["left"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(False)

    # optional: simplify x ticks
    ax1.set_xticks([60, 62, 64, 66, 68, 70, 72])

    # -------------------------
    # Inset for far-away mutations (506/507 region)
    # -------------------------
    if not far_df.empty:
        ax1_inset = inset_axes(ax1, width="26%", height="42%", loc="upper right", borderpad=1.2)

        ax1_inset.hlines(0, 505.5, 507.5, linewidth=2)

        for _, row in far_df.iterrows():
            marker = "o" if row["group"] == "Within 5 aa" else "s"
            ax1_inset.scatter(
                row["mutation_position"],
                0,
                s=110,
                marker=marker,
                edgecolor="black",
                zorder=3
            )

        add_staggered_labels(ax1_inset, far_df, proximity=0.5, fontsize=9)

        ax1_inset.set_xlim(505.4, 507.6)
        ax1_inset.set_ylim(-0.5, 0.5)
        ax1_inset.set_xticks([506, 507])
        ax1_inset.set_yticks([])
        ax1_inset.set_ylabel("")
        ax1_inset.spines["left"].set_visible(False)
        ax1_inset.spines["right"].set_visible(False)
        ax1_inset.grid(False)
        ax1_inset.set_title("Far-away mutations", fontsize=9)
    
    
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
            label="Within 5 AA in linear sequence"
        )

    if not farther.empty:
        ax2.scatter(
            farther["x"],
            farther["distance_3d"],
            s=100,
            marker="s",
            label="More than 5 AA in linear sequence"
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
        fig.savefig(save_path, dpi=300, bbox_inches="tight", pad_inches=0.25)

    plt.show()

# ----------------------------
# Example usage
# ----------------------------
df = pd.read_csv(r"/Users/alissaallen/gitScratch/Scratch/ptm_mutation_proximity_db.csv.csv")

gene_of_interest = "PTPN11"
ptm_of_interest = "Y62"

mutation_df = build_ptm_mutation_table(df, gene_of_interest, ptm_of_interest)
print(mutation_df)


plot_two_panel_ptm_figure(
    mutation_df,
    gene=gene_of_interest,
    ptm_site=ptm_of_interest,
    save_path="PTPN11_Y62_two_panel_figure.png",
    local_window=8,
    structural_cutoff=10.0
)