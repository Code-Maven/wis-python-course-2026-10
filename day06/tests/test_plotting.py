"""Headless smoke tests for Matplotlib chart builders."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from ncbi_protein_analyzer.analysis import analyze_dataset, compare_records
from ncbi_protein_analyzer.config import AnalysisConfig, MsaConfig
from ncbi_protein_analyzer.models import ProteinRecord
from ncbi_protein_analyzer.plotting import (
    create_aa_composition_chart,
    create_identity_heatmap,
    create_length_histogram,
)


def test_chart_builders_return_figures(sample_records: list[ProteinRecord]) -> None:
    stats, summary = analyze_dataset(sample_records)
    comparison = compare_records(
        sample_records, AnalysisConfig(msa=MsaConfig(preferred_tool="pairwise"))
    )

    fig_hist = create_length_histogram(stats)
    fig_aa = create_aa_composition_chart(summary.aggregate_aa_composition)
    fig_heat = create_identity_heatmap(comparison)

    assert len(fig_hist.axes) == 1
    assert len(fig_aa.axes) == 1
    assert len(fig_heat.axes) >= 1


def test_chart_builders_handle_empty_inputs() -> None:
    fig_hist = create_length_histogram([])
    fig_aa = create_aa_composition_chart({})
    assert len(fig_hist.axes) == 1
    assert len(fig_aa.axes) == 1
