import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent


def extract_ptm_position(ptm_site: str) -> int:
    match = re.search(r"(\d+)", str(ptm_site))
    if not match:
        raise ValueError(f"Could not parse PTM site: {ptm_site}")
    return int(match.group(1))


def parse_mutation_entry(entry: str):
    if pd.isna(entry) or str(entry).strip() == "":
        return None
    entry = str(entry).strip()
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
    if pd.isna(cell_value) or str(cell_value).strip() == "":
        return []
    text = str(cell_value).strip()
    parts = re.split(r"[;,]\s*", text)
    return [p for p in parts if p.strip()]


def build_ptm_mutation_table(df: pd.DataFrame, gene: str, ptm_site: str) -> pd.DataFrame:
    subset = df[(df["gene"] == gene) & (df["ptm_site"] == ptm_site)].copy()
    if subset.empty:
        raise ValueError(f"No rows found for gene={gene}, ptm_site={ptm_site}")
    rows = []
    for _, row in subset.iterrows():
        ptm_pos = extract_ptm_position(row["ptm_site"])
        for entry in split_mutation_column(row.get("mutations_within_5_positions", "")):
            parsed = parse_mutation_entry(entry)
            if parsed is None:
                continue
            mutation_name, mutation_pos, dist_3d = parsed
            rows.append({
                "UniProt": row["UniProt"], "gene": row["gene"],
                "ptm_site": row["ptm_site"], "ptm_type": row["ptm_type"],
                "mutation": mutation_name, "mutation_position": mutation_pos,
                "ptm_position": ptm_pos,
                "linear_distance": abs(mutation_pos - ptm_pos),
                "distance_3d": dist_3d, "group": "Within 5 aa"
            })
        for entry in split_mutation_column(row.get("mutations_more_than_5_positions", "")):
            parsed = parse_mutation_entry(entry)
            if parsed is None:
                continue
            mutation_name, mutation_pos, dist_3d = parsed
            rows.append({
                "UniProt": row["UniProt"], "gene": row["gene"],
                "ptm_site": row["ptm_site"], "ptm_type": row["ptm_type"],
                "mutation": mutation_name, "mutation_position": mutation_pos,
                "ptm_position": ptm_pos,
                "linear_distance": abs(mutation_pos - ptm_pos),
                "distance_3d": dist_3d, "group": "More than 5 aa"
            })
    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError(f"No parseable mutations found for gene={gene}, ptm_site={ptm_site}")
    result = result.sort_values("mutation_position").reset_index(drop=True)
    return result


def add_staggered_labels(ax, df, proximity=0.7, fontsize=10):
    """
    Place mutation labels above/below the sequence line.
    When multiple mutations share the same (or very close) position,
    fan them out diagonally so lines and text don't pile on top of each other.

    Strategy:
    - Group mutations by position (rounded to nearest 0.5 to catch near-overlaps).
    - For a lone mutation: straight vertical line, label directly above/below.
    - For a crowded group (2+ at the same position):
        * Split evenly: first half go above, second half below.
        * Within each half, fan diagonally: each successive label is shifted
          both higher (larger y) AND horizontally (larger x_offset), producing
          a staircase / diagonal layout.
    """
    from collections import defaultdict

    # ── 1. bucket mutations into positional groups ──────────────────────────
    def bucket_key(x):
        return round(x / proximity) * proximity  # snap to grid

    groups = defaultdict(list)
    for _, row in df.iterrows():
        key = bucket_key(row["mutation_position"])
        groups[key].append(row)

    # ── 2. draw each group ──────────────────────────────────────────────────
    # Tuning knobs
    y_lone        = 0.22   # height for a solo label
    y_fan_start   = 0.20   # lowest label height in a crowded fan
    y_fan_step    = 0.13   # additional height per label going left→right
    x_fan_step    = 0.55   # horizontal spread between fan labels (data units)

    for key in sorted(groups.keys()):
        members = groups[key]
        n = len(members)

        members_sorted = sorted(members, key=lambda r: r["mutation_position"])
        x_anchor = sum(r["mutation_position"] for r in members_sorted) / n

        if n == 1:
            # ── simple case: straight vertical line ──────────────────────
            row = members_sorted[0]
            x = row["mutation_position"]
            ax.plot([x, x], [0, y_lone - 0.02],
                    color="black", lw=0.8, zorder=2, clip_on=False)
            ax.text(x, y_lone, row["mutation"],
                    ha="center", va="bottom", fontsize=fontsize, clip_on=False)

        elif n == 2:
            # ── two: one straight up, one straight down ───────────────────
            for row, (y_label, y_line_end, va) in zip(
                members_sorted,
                [(y_lone, y_lone - 0.02, "bottom"), (-y_lone, -y_lone + 0.02, "top")]
            ):
                x = row["mutation_position"]
                ax.plot([x, x], [0, y_line_end],
                        color="black", lw=0.8, zorder=2, clip_on=False)
                ax.text(x, y_label, row["mutation"],
                        ha="center", va=va, fontsize=fontsize, clip_on=False)

        else:
            # ── 3+: split top/bottom, fan each half diagonally ───────────
            top_members    = members_sorted[: (n + 1) // 2]
            bottom_members = members_sorted[(n + 1) // 2 :]

            n_top = len(top_members)
            for i, row in enumerate(top_members):
                x_label = x_anchor + (i - (n_top - 1) / 2) * x_fan_step
                y_label = y_fan_start + i * y_fan_step
                ax.plot([row["mutation_position"], x_label], [0, y_label - 0.02],
                        color="black", lw=0.8, zorder=2, clip_on=False)
                ax.text(x_label, y_label, row["mutation"],
                        ha="center", va="bottom", fontsize=fontsize, clip_on=False)

            n_bot = len(bottom_members)
            for i, row in enumerate(bottom_members):
                x_label = x_anchor + (i - (n_bot - 1) / 2) * x_fan_step
                y_label = -(y_fan_start + i * y_fan_step)
                ax.plot([row["mutation_position"], x_label], [0, y_label + 0.02],
                        color="black", lw=0.8, zorder=2, clip_on=False)
                ax.text(x_label, y_label, row["mutation"],
                        ha="center", va="top", fontsize=fontsize, clip_on=False)


def style_sequence_axis(ax, xlim, xticks):
    """Apply common styling to sequence axes."""
    ax.set_xlim(xlim)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_ylabel("")
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.grid(False)
    if xticks is not None:
        ax.set_xticks(xticks)


def plot_two_panel_ptm_figure(
    mutation_df: pd.DataFrame,
    gene: str,
    ptm_site: str,
    save_path: str | None = None,
    local_window: int = 10,
    structural_cutoff: float = 10.0
):
    ptm_pos = mutation_df["ptm_position"].iloc[0]
    mutation_df = mutation_df.sort_values("mutation_position").reset_index(drop=True)

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

    has_far = not far_df.empty

    # ── Layout ──────────────────────────────────────────────────────────────
    # Row 0: sequence panel (local | gap indicator | far-away)
    # Row 1: 3-D distance panel
    # The local region gets ~75% of the top row width, far-away ~25%.

    fig = plt.figure(figsize=(15, 8))

    if has_far:
        gs_top = gridspec.GridSpec(
            1, 3,
            figure=fig,
            width_ratios=[0.72, 0.04, 0.24],   # local | gap | far
            wspace=0.0,
            left=0.06, right=0.97, top=0.93, bottom=0.54
        )
        ax_local = fig.add_subplot(gs_top[0])
        ax_gap   = fig.add_subplot(gs_top[1])
        ax_far   = fig.add_subplot(gs_top[2])
    else:
        gs_top = gridspec.GridSpec(
            1, 1,
            figure=fig,
            left=0.06, right=0.97, top=0.93, bottom=0.54
        )
        ax_local = fig.add_subplot(gs_top[0])

    gs_bot = gridspec.GridSpec(
        1, 1,
        figure=fig,
        left=0.06, right=0.97, top=0.46, bottom=0.10
    )
    ax2 = fig.add_subplot(gs_bot[0])

    # ── Panel A-left: local sequence context ────────────────────────────────
    ax_local.hlines(0, local_min, local_max, linewidth=2, color="#888888")
    ax_local.axvspan(ptm_pos - 5, ptm_pos + 5, alpha=0.12, color="#2C6FAC")

    # Label the ±5 aa window at the top of the shaded region
    ax_local.text(
        ptm_pos, 0.50,
        "±5 aa window",
        ha="center", va="bottom",
        fontsize=9, color="#2C6FAC",
        style="italic"
    )

    # PTM site marker
    ax_local.scatter(ptm_pos, 0, s=200, zorder=4, edgecolor="black",
                     color="#FFE500")
    ax_local.annotate(
        ptm_site, (ptm_pos, 0),
        xytext=(0, 18), textcoords="offset points",
        ha="center", fontsize=14, fontweight="bold"
    )

    for _, row in local_df.iterrows():
        marker = "o" if row["group"] == "Within 5 aa" else "s"
        ax_local.scatter(
            row["mutation_position"], 0,
            s=140, marker=marker,
            color="#888888", edgecolor="black", zorder=3
        )

    add_staggered_labels(ax_local, local_df, proximity=0.7, fontsize=10)

    style_sequence_axis(
        ax_local,
        xlim=(local_min - 0.5, local_max + 0.5),
        xticks=[60, 62, 64, 66, 68, 70, 72]
    )
    # Remove right spine so it blends into the gap column
    ax_local.spines["right"].set_visible(False)

    # ── Panel A-middle: broken-axis gap indicator ────────────────────────────
    if has_far:
        ax_gap.set_xlim(0, 1)
        ax_gap.set_ylim(-0.55, 0.55)
        ax_gap.axis("off")

        # Draw two diagonal slash marks to signal the axis break
        slash_kw = dict(transform=ax_gap.transAxes, color="gray",
                        lw=1.2, clip_on=False)
        d = 0.08
        for x_center in [0.25, 0.75]:
            ax_gap.plot(
                [x_center - d, x_center + d],
                [0.44 - 0.05, 0.44 + 0.05],
                **slash_kw
            )
            ax_gap.plot(
                [x_center - d, x_center + d],
                [0.54 - 0.05, 0.54 + 0.05],
                **slash_kw
            )

        # Dotted connecting line through the gap at y=0.5 (axis mid-point)
        ax_gap.plot([0, 1], [0.5, 0.5], lw=1.5, ls=":",
                    color="#888888", transform=ax_gap.transAxes, clip_on=False)

    # ── Panel A-right: far-away mutations ───────────────────────────────────
    if has_far:
        far_xmin = far_df["mutation_position"].min() - 0.8
        far_xmax = far_df["mutation_position"].max() + 0.8

        ax_far.hlines(0, far_xmin, far_xmax, linewidth=2, color="#888888")

        for _, row in far_df.iterrows():
            marker = "o" if row["group"] == "Within 5 aa" else "s"
            ax_far.scatter(
                row["mutation_position"], 0,
                s=140, marker=marker,
                color="#888888", edgecolor="black", zorder=3
            )

        add_staggered_labels(ax_far, far_df, proximity=0.5, fontsize=10)

        far_xticks = sorted(far_df["mutation_position"].unique().tolist())
        style_sequence_axis(
            ax_far,
            xlim=(far_xmin, far_xmax),
            xticks=far_xticks
        )

        # Box the far panel with a subtle frame to signal it's a separate region
        for spine in ["top", "bottom"]:
            ax_far.spines[spine].set_visible(True)
            ax_far.spines[spine].set_color("lightgray")
            ax_far.spines[spine].set_linewidth(0.8)
        ax_far.spines["left"].set_visible(True)
        ax_far.spines["left"].set_color("lightgray")
        ax_far.spines["left"].set_linewidth(0.8)
        ax_far.spines["right"].set_visible(True)
        ax_far.spines["right"].set_color("lightgray")
        ax_far.spines["right"].set_linewidth(0.8)

        ax_far.set_title("Far-away mutations", fontsize=9, color="gray", pad=4)

    # ── Panel B: 3-D distance ────────────────────────────────────────────────
    plot_df = mutation_df.copy()
    plot_df["x"] = range(len(plot_df))

    within  = plot_df[plot_df["group"] == "Within 5 aa"]
    farther = plot_df[plot_df["group"] == "More than 5 aa"]

    if not within.empty:
        ax2.scatter(within["x"], within["distance_3d"],
                    s=100, marker="o", color="#2C6FAC",
                    label="Within 5 AA in linear sequence", zorder=3)
    if not farther.empty:
        ax2.scatter(farther["x"], farther["distance_3d"],
                    s=100, marker="s", color="#E07B39",
                    label="More than 5 AA in linear sequence", zorder=3)

    ax2.axhline(structural_cutoff, linestyle="--", color="gray",
                alpha=0.6, linewidth=1)
    ax2.text(
        len(plot_df) - 0.2, structural_cutoff + 0.05,
        f"{structural_cutoff:g} Å",
        ha="right", va="bottom", fontsize=9, color="gray"
    )

    for i, (_, row) in enumerate(plot_df.iterrows()):
        y_offset = 8 if i % 2 == 0 else 14
        ax2.annotate(
            row["mutation"],
            (row["x"], row["distance_3d"]),
            xytext=(0, y_offset),
            textcoords="offset points",
            ha="center", fontsize=9
        )

    ax2.set_xticks(plot_df["x"])
    ax2.set_xticklabels(plot_df["mutation"], rotation=45, ha="right")
    ax2.set_ylabel("3D distance from PTM site (Å)")
    ax2.set_xlabel("Mutation")
    ax2.grid(True, axis="y", alpha=0.3)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(loc="upper left")

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight", pad_inches=0.25)

    plt.show()


# ----------------------------
# Example usage
# ----------------------------
df = pd.read_csv(PROJECT_ROOT / "Output" / "ptm_mutation_proximity_db.tsv", sep="\t", encoding="utf-16")

gene_of_interest = "PTPN11"
ptm_of_interest  = "Y62"

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