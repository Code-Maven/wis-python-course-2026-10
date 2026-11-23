"""Unit tests for deduplication, external MSA mocking, and pairwise fallback."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import pytest

from ncbi_protein_analyzer.analysis.compare import (
    compare_records,
    compute_consensus,
    compute_msa_identity_matrix,
)
from ncbi_protein_analyzer.analysis.deduplicate import deduplicate_records
from ncbi_protein_analyzer.analysis.msa import DetectedMsaTool, detect_msa_tool, run_msa
from ncbi_protein_analyzer.config import AnalysisConfig, MsaConfig
from ncbi_protein_analyzer.models import ProteinRecord


def test_deduplicate_groups_exact_matches(sample_records: list[ProteinRecord]) -> None:
    result = deduplicate_records(sample_records)
    assert len(result.unique_records) == 2
    assert result.duplicate_count == 1
    assert result.clusters["NP_000198.1"] == ["NP_000198.1", "XP_999999.1"]
    assert result.clusters["NP_032412.3"] == ["NP_032412.3"]


def test_compare_falls_back_to_pairwise_when_no_msa_binary(
    sample_records: list[ProteinRecord], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("shutil.which", lambda _name: None)
    cfg = AnalysisConfig(max_pairwise_sequences=20, msa=MsaConfig(preferred_tool="auto"))
    result = compare_records(sample_records, cfg)
    assert result.method_used == "pairwise"
    assert result.deduplication.duplicate_count == 1
    assert result.labels == ["NP_000198.1", "NP_032412.3"]
    assert result.identity_matrix[0][0] == 100.0
    assert 50.0 < result.identity_matrix[0][1] < 100.0
    assert any("pairwise fallback" in w for w in result.warnings)


def test_pairwise_handles_non_blosum_ambiguous_letters() -> None:
    records = [
        ProteinRecord("A1", "s1", "Org", "ACDEFGHUJO*KLM"),
        ProteinRecord("A2", "s2", "Org", "ACDEFGHIKLMNPQ"),
    ]
    cfg = AnalysisConfig(msa=MsaConfig(preferred_tool="pairwise"))
    result = compare_records(records, cfg)
    assert result.method_used == "pairwise"
    assert 0.0 <= result.identity_matrix[0][1] <= 100.0


def test_run_msa_with_mocked_clustalo(sample_records: list[ProteinRecord]) -> None:
    unique = deduplicate_records(sample_records).unique_records

    def _fake_runner(cmd: list[str], **_kwargs: Any) -> subprocess.CompletedProcess[str]:
        assert "-i" in cmd and "-o" in cmd
        out_path = Path(cmd[cmd.index("-o") + 1])
        out_path.write_text(
            ">seq_0\n"
            f"{unique[0].sequence}\n"
            ">seq_1\n"
            f"{unique[1].sequence}\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    tool = DetectedMsaTool(name="clustalo", executable="/usr/bin/clustalo")
    aligned = run_msa(unique, tool, runner=_fake_runner)
    assert len(aligned) == 2
    assert aligned[0].accession == "NP_000198.1"
    assert aligned[1].accession == "NP_032412.3"


def test_compare_uses_mocked_muscle_v5(
    sample_records: list[ProteinRecord], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "shutil.which", lambda name: "/usr/bin/muscle" if name == "muscle" else None
    )
    seen_cmds: list[list[str]] = []

    def _fake_runner(cmd: list[str], **_kwargs: Any) -> subprocess.CompletedProcess[str]:
        seen_cmds.append(cmd)
        if "-version" in cmd:
            return subprocess.CompletedProcess(cmd, 0, stdout="muscle 5.1.linux64", stderr="")
        assert "-align" in cmd and "-output" in cmd
        in_path = Path(cmd[cmd.index("-align") + 1])
        out_path = Path(cmd[cmd.index("-output") + 1])
        # Echo input FASTA back as aligned
        out_path.write_text(in_path.read_text(encoding="utf-8"), encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    cfg = AnalysisConfig(msa=MsaConfig(preferred_tool="auto"))
    result = compare_records(sample_records, cfg, runner=_fake_runner)
    assert result.method_used == "muscle"
    assert len(result.aligned_records) == 2
    assert len(result.consensus) == 60
    assert any("-align" in c for c in seen_cmds)


def test_msa_failure_falls_back_to_pairwise(
    sample_records: list[ProteinRecord], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "shutil.which", lambda name: "/usr/bin/clustalo" if name == "clustalo" else None
    )

    def _failing_runner(cmd: list[str], **_kwargs: Any) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="segfault")

    cfg = AnalysisConfig(msa=MsaConfig(preferred_tool="auto"))
    result = compare_records(sample_records, cfg, runner=_failing_runner)
    assert result.method_used == "pairwise"
    assert any("falling back to Biopython pairwise alignment" in w for w in result.warnings)


def test_pairwise_cap_enforced() -> None:
    records = [
        ProteinRecord(f"P{i}", f"prot {i}", "Org", f"ACDEFGHIKL{'A' * (i + 1)}")
        for i in range(5)
    ]
    cfg = AnalysisConfig(max_pairwise_sequences=3, msa=MsaConfig(preferred_tool="pairwise"))
    result = compare_records(records, cfg)
    assert len(result.labels) == 3
    assert any("capped at first 3" in w for w in result.warnings)


def test_msa_identity_and_consensus_calculations() -> None:
    aligned = [
        ProteinRecord("A", "", "", "AC-DE"),
        ProteinRecord("B", "", "", "AC-DF"),
        ProteinRecord("C", "", "", "AC-DE"),
    ]
    matrix = compute_msa_identity_matrix(aligned)
    assert matrix[0][2] == pytest.approx(100.0)
    assert matrix[0][1] == pytest.approx(75.0)
    assert compute_consensus(aligned) == "AC-DE"


def test_detect_msa_tool_respects_pairwise_preference(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: f"/usr/bin/{name}")
    assert detect_msa_tool(MsaConfig(preferred_tool="pairwise")) is None
