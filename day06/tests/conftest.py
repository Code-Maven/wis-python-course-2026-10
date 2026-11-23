"""Shared pytest fixtures for the NCBI Protein Analyzer test suite."""

from __future__ import annotations

from pathlib import Path

import pytest

from ncbi_protein_analyzer.models import ProteinRecord, ProteinSummary


@pytest.fixture
def sample_config_yaml(tmp_path: Path) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        """
ncbi:
  email: "researcher@institute.edu"
  api_key: ""
  tool: "test_tool"
  max_results_default: 25
  fetch_batch_size: 10
  timeout_seconds: 15
analysis:
  max_pairwise_sequences: 20
  msa:
    preferred_tool: "auto"
    clustalo_path: ""
    muscle_path: ""
    timeout_seconds: 30
""".strip(),
        encoding="utf-8",
    )
    return cfg


@pytest.fixture
def sample_summaries() -> list[ProteinSummary]:
    return [
        ProteinSummary(
            uid="1001",
            accession="NP_000198.1",
            title="insulin preproprotein [Homo sapiens]",
            organism="Homo sapiens",
            length=110,
            update_date="2025/01/10",
        ),
        ProteinSummary(
            uid="1002",
            accession="NP_032412.3",
            title="insulin-1 preproprotein [Mus musculus]",
            organism="Mus musculus",
            length=108,
            update_date="2025/02/15",
        ),
    ]


@pytest.fixture
def sample_records() -> list[ProteinRecord]:
    return [
        ProteinRecord(
            accession="NP_000198.1",
            description="insulin preproprotein [Homo sapiens]",
            organism="Homo sapiens",
            sequence="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED",
        ),
        ProteinRecord(
            accession="NP_032412.3",
            description="insulin-1 preproprotein [Mus musculus]",
            organism="Mus musculus",
            sequence="MALLVHFLPLLALLALWEPKPTQAFVKQHLCGPHLVEALYLVCGERGFFYTPKSRREVED",
        ),
        ProteinRecord(
            accession="XP_999999.1",
            description="predicted insulin-like protein [Canis lupus]",
            organism="Canis lupus",
            sequence="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED",
        ),
    ]
