"""Unit tests for CSV, FASTA, and PNG export writers."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
from Bio import SeqIO

matplotlib.use("Agg")

from ncbi_protein_analyzer.analysis import analyze_dataset
from ncbi_protein_analyzer.export import (
    export_all_figures_png,
    export_fasta,
    export_figure_png,
    export_stats_csv,
)
from ncbi_protein_analyzer.models import ProteinRecord
from ncbi_protein_analyzer.plotting import create_aa_composition_chart, create_length_histogram


def test_export_stats_csv(tmp_path: Path, sample_records: list[ProteinRecord]) -> None:
    stats, _summary = analyze_dataset(sample_records)
    csv_path = export_stats_csv(stats, tmp_path / "stats.csv")
    assert csv_path.exists()

    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 3
    assert rows[0]["accession"] == "NP_000198.1"
    assert rows[0]["is_approximate"] == "False"
    assert "aa_A_pct" in rows[0]


def test_export_fasta(tmp_path: Path, sample_records: list[ProteinRecord]) -> None:
    fasta_path = export_fasta(sample_records, tmp_path / "proteins.fasta")
    parsed = list(SeqIO.parse(fasta_path, "fasta"))
    assert len(parsed) == 3
    assert parsed[0].id == "NP_000198.1"
    assert str(parsed[0].seq) == sample_records[0].sequence


def test_export_png_figures(tmp_path: Path, sample_records: list[ProteinRecord]) -> None:
    stats, summary = analyze_dataset(sample_records)
    fig1 = create_length_histogram(stats)
    fig2 = create_aa_composition_chart(summary.aggregate_aa_composition)

    single = export_figure_png(fig1, tmp_path / "hist.png")
    assert single.exists() and single.stat().st_size > 1000

    saved = export_all_figures_png(
        {"length_histogram": fig1, "aa_composition": fig2},
        tmp_path / "charts",
    )
    assert len(saved) == 2
    assert all(p.exists() and p.stat().st_size > 1000 for p in saved)
