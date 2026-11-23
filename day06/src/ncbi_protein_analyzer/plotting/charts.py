"""Standalone Matplotlib Figure builders (safe for Tk embedding and headless PNG export)."""

from __future__ import annotations

from collections.abc import Sequence
from statistics import mean, median

from matplotlib.figure import Figure
from matplotlib.patches import Patch

from ncbi_protein_analyzer.models import STANDARD_AAS, ComparisonResult, ProteinStats

_AA_GROUPS: dict[str, tuple[str, str]] = {
    # Hydrophobic
    "A": ("Hydrophobic", "#3b82f6"),
    "V": ("Hydrophobic", "#3b82f6"),
    "I": ("Hydrophobic", "#3b82f6"),
    "L": ("Hydrophobic", "#3b82f6"),
    "M": ("Hydrophobic", "#3b82f6"),
    "F": ("Hydrophobic", "#3b82f6"),
    "W": ("Hydrophobic", "#3b82f6"),
    # Polar uncharged
    "S": ("Polar", "#10b981"),
    "T": ("Polar", "#10b981"),
    "N": ("Polar", "#10b981"),
    "Q": ("Polar", "#10b981"),
    "Y": ("Polar", "#10b981"),
    # Positively charged
    "K": ("Positive", "#f59e0b"),
    "R": ("Positive", "#f59e0b"),
    "H": ("Positive", "#f59e0b"),
    # Negatively charged
    "D": ("Negative", "#ef4444"),
    "E": ("Negative", "#ef4444"),
    # Special
    "C": ("Special", "#8b5cf6"),
    "G": ("Special", "#8b5cf6"),
    "P": ("Special", "#8b5cf6"),
}


def create_length_histogram(stats: Sequence[ProteinStats]) -> Figure:
    """Create a histogram of protein sequence lengths with mean and median markers."""
    fig = Figure(figsize=(7.5, 4.5), dpi=100)
    ax = fig.add_subplot(111)

    lengths = [s.length for s in stats]
    if not lengths:
        ax.text(0.5, 0.5, "No sequence data to plot", ha="center", va="center", fontsize=12)
        ax.set_axis_off()
        return fig

    bins = min(25, max(5, len(set(lengths))))
    ax.hist(lengths, bins=bins, color="#3b82f6", edgecolor="#1e3a8a", alpha=0.82)

    mean_len = mean(lengths)
    median_len = median(lengths)
    ax.axvline(mean_len, color="#dc2626", linestyle="--", linewidth=1.6, label=f"Mean ({mean_len:.1f} aa)")
    ax.axvline(median_len, color="#059669", linestyle=":", linewidth=1.8, label=f"Median ({median_len:.1f} aa)")

    ax.set_title(f"Protein Length Distribution (n = {len(lengths)})", fontsize=12, fontweight="bold")
    ax.set_xlabel("Sequence Length (amino acids)")
    ax.set_ylabel("Number of Proteins")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.legend(frameon=True)
    fig.tight_layout()
    return fig


def create_aa_composition_chart(
    aa_composition: dict[str, float],
    *,
    title: str = "Amino Acid Composition (Dataset Average)",
) -> Figure:
    """Create a bar chart of the 20 standard amino acids colored by chemical group."""
    fig = Figure(figsize=(8.0, 4.5), dpi=100)
    ax = fig.add_subplot(111)

    if not aa_composition or all(v == 0.0 for v in aa_composition.values()):
        ax.text(0.5, 0.5, "No amino acid composition data", ha="center", va="center", fontsize=12)
        ax.set_axis_off()
        return fig

    aas = list(STANDARD_AAS)
    values = [aa_composition.get(aa, 0.0) for aa in aas]
    colors = [_AA_GROUPS[aa][1] for aa in aas]

    bars = ax.bar(aas, values, color=colors, edgecolor="#1f2937", linewidth=0.5)
    for bar, val in zip(bars, values):
        if val >= 0.5:
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                val + 0.15,
                f"{val:.1f}",
                ha="center",
                va="bottom",
                fontsize=7.5,
            )

    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, max(values) * 1.18 + 0.5)
    ax.grid(axis="y", linestyle="--", alpha=0.35)

    legend_handles = [
        Patch(facecolor="#3b82f6", label="Hydrophobic"),
        Patch(facecolor="#10b981", label="Polar"),
        Patch(facecolor="#f59e0b", label="Positive"),
        Patch(facecolor="#ef4444", label="Negative"),
        Patch(facecolor="#8b5cf6", label="Special"),
    ]
    ax.legend(handles=legend_handles, ncol=5, loc="upper right", fontsize=8, frameon=True)
    fig.tight_layout()
    return fig


def create_identity_heatmap(comparison: ComparisonResult) -> Figure:
    """Create a pairwise sequence identity heatmap (0-100%)."""
    fig = Figure(figsize=(7.2, 5.8), dpi=100)
    ax = fig.add_subplot(111)

    matrix = comparison.identity_matrix
    labels = comparison.labels
    n = len(labels)

    if n < 2 or not matrix:
        ax.text(
            0.5,
            0.5,
            "At least 2 unique sequences are required\nto render an identity heatmap.",
            ha="center",
            va="center",
            fontsize=11,
        )
        ax.set_axis_off()
        return fig

    im = ax.imshow(matrix, vmin=0.0, vmax=100.0, cmap="viridis", aspect="auto")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Sequence Identity (%)")

    short_labels = [lbl if len(lbl) <= 16 else f"{lbl[:14]}.." for lbl in labels]
    tick_fontsize = 8 if n <= 20 else 6.5
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(short_labels, rotation=45, ha="right", fontsize=tick_fontsize)
    ax.set_yticklabels(short_labels, fontsize=tick_fontsize)

    if n <= 15:
        for i in range(n):
            for j in range(n):
                val = matrix[i][j]
                text_color = "black" if val > 65 else "white"
                ax.text(j, i, f"{val:.0f}", ha="center", va="center", color=text_color, fontsize=7.5)

    method_label = {
        "clustalo": "Clustal Omega MSA",
        "muscle": "MUSCLE MSA",
        "pairwise": "Biopython Pairwise Fallback",
    }.get(comparison.method_used, comparison.method_used)

    ax.set_title(f"Pairwise Sequence Identity ({method_label})", fontsize=12, fontweight="bold")
    fig.tight_layout()
    return fig
