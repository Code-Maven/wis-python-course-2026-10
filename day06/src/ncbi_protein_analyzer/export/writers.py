"""File exporters for sequence statistics (CSV), sequences/alignments (FASTA), and plots (PNG)."""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.figure import Figure

from ncbi_protein_analyzer.models import STANDARD_AAS, ProteinRecord, ProteinStats


def export_stats_csv(stats: Sequence[ProteinStats], path: str | Path) -> Path:
    """Write per-protein physicochemical statistics and AA composition to CSV."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    base_fields = [
        "accession",
        "description",
        "organism",
        "length",
        "standard_aa_count",
        "ambiguous_aa_count",
        "ambiguous_aa_pct",
        "is_approximate",
        "molecular_weight_da",
        "isoelectric_point",
        "gravy",
        "warning",
    ]
    aa_fields = [f"aa_{aa}_pct" for aa in STANDARD_AAS]
    fieldnames = base_fields + aa_fields

    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for item in stats:
            row: dict[str, object] = {
                "accession": item.accession,
                "description": item.description,
                "organism": item.organism,
                "length": item.length,
                "standard_aa_count": item.standard_aa_count,
                "ambiguous_aa_count": item.ambiguous_aa_count,
                "ambiguous_aa_pct": f"{item.ambiguous_aa_pct:.2f}",
                "is_approximate": item.is_approximate,
                "molecular_weight_da": (
                    f"{item.molecular_weight_da:.2f}" if item.molecular_weight_da is not None else ""
                ),
                "isoelectric_point": (
                    f"{item.isoelectric_point:.3f}" if item.isoelectric_point is not None else ""
                ),
                "gravy": f"{item.gravy:.4f}" if item.gravy is not None else "",
                "warning": item.warning,
            }
            for aa in STANDARD_AAS:
                row[f"aa_{aa}_pct"] = f"{item.aa_composition.get(aa, 0.0):.2f}"
            writer.writerow(row)
    return out_path


def export_fasta(records: Sequence[ProteinRecord], path: str | Path) -> Path:
    """Write protein records (unaligned or aligned) to a FASTA file."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    seq_records = [
        SeqRecord(
            Seq(rec.sequence),
            id=rec.accession,
            description=rec.description,
        )
        for rec in records
    ]
    with out_path.open("w", encoding="utf-8") as fh:
        SeqIO.write(seq_records, fh, "fasta")
    return out_path


def export_figure_png(fig: Figure, path: str | Path, *, dpi: int = 150) -> Path:
    """Save a Matplotlib Figure as a PNG image."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="png", dpi=dpi, bbox_inches="tight")
    return out_path


def export_all_figures_png(
    figures: Mapping[str, Figure],
    directory: str | Path,
    *,
    dpi: int = 150,
) -> list[Path]:
    """Save all named Matplotlib Figures into a directory as PNG files."""
    out_dir = Path(directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
    for name, fig in figures.items():
        safe_name = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in name).strip("_")
        if not safe_name:
            safe_name = "chart"
        target = out_dir / f"{safe_name}.png"
        saved.append(export_figure_png(fig, target, dpi=dpi))
    return saved
