"""Unit tests for physicochemical stats and ambiguous amino acid handling."""

from __future__ import annotations

import pytest

from ncbi_protein_analyzer.analysis.stats import analyze_dataset, analyze_record
from ncbi_protein_analyzer.models import ProteinRecord


def test_analyze_standard_sequence() -> None:
    rec = ProteinRecord(
        accession="P01308",
        description="Insulin",
        organism="Homo sapiens",
        sequence="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT",
    )
    stats = analyze_record(rec)
    assert stats.length == 54
    assert stats.standard_aa_count == 54
    assert stats.ambiguous_aa_count == 0
    assert stats.ambiguous_aa_pct == 0.0
    assert stats.is_approximate is False
    assert stats.molecular_weight_da is not None and stats.molecular_weight_da > 5000
    assert stats.isoelectric_point is not None and 4.0 < stats.isoelectric_point < 10.0
    assert stats.gravy is not None
    assert sum(stats.aa_composition.values()) == pytest.approx(100.0, abs=1e-6)
    assert stats.warning == ""


def test_analyze_sequence_with_ambiguous_residues_flags_approximate() -> None:
    rec = ProteinRecord(
        accession="XP_AMBIG",
        description="Predicted partial protein",
        organism="Danio rerio",
        sequence="ACDEFGHIKLMNPQRSTVWY XBZJUO*",
    )
    stats = analyze_record(rec)
    assert stats.length == 27
    assert stats.standard_aa_count == 20
    assert stats.ambiguous_aa_count == 7
    assert stats.ambiguous_aa_pct == pytest.approx(7 / 27 * 100.0)
    assert stats.is_approximate is True
    assert stats.molecular_weight_da is not None
    assert stats.isoelectric_point is not None
    assert stats.gravy is not None
    assert "Approximate" in stats.warning
    assert all(val == pytest.approx(5.0) for val in stats.aa_composition.values())


def test_analyze_all_ambiguous_sequence_returns_none_metrics() -> None:
    rec = ProteinRecord(
        accession="XP_ALLX",
        description="Unknown fragment",
        organism="",
        sequence="XXXXXBBB*",
    )
    stats = analyze_record(rec)
    assert stats.length == 9
    assert stats.standard_aa_count == 0
    assert stats.ambiguous_aa_count == 9
    assert stats.is_approximate is True
    assert stats.molecular_weight_da is None
    assert stats.isoelectric_point is None
    assert stats.gravy is None
    assert "unavailable" in stats.warning


def test_analyze_dataset_aggregates(sample_records: list[ProteinRecord]) -> None:
    per_protein, summary = analyze_dataset(sample_records)
    assert len(per_protein) == 3
    assert summary.count == 3
    assert summary.min_length == 60
    assert summary.max_length == 60
    assert summary.mean_length == 60.0
    assert summary.approximate_count == 0
    assert summary.mean_molecular_weight_da is not None
    assert sum(summary.aggregate_aa_composition.values()) == pytest.approx(100.0, abs=1e-6)


def test_analyze_empty_dataset() -> None:
    per_protein, summary = analyze_dataset([])
    assert per_protein == []
    assert summary.count == 0
    assert summary.mean_molecular_weight_da is None
