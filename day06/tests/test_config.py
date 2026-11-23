"""Unit tests for YAML configuration loading and validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from ncbi_protein_analyzer.config import ConfigError, load_config


def test_load_valid_config(sample_config_yaml: Path) -> None:
    cfg = load_config(sample_config_yaml)
    assert cfg.ncbi.email == "researcher@institute.edu"
    assert cfg.ncbi.max_results_default == 25
    assert cfg.ncbi.fetch_batch_size == 10
    assert cfg.analysis.max_pairwise_sequences == 20
    assert cfg.analysis.msa.preferred_tool == "auto"


def test_missing_config_file_raises(tmp_path: Path) -> None:
    missing = tmp_path / "does_not_exist.yaml"
    with pytest.raises(ConfigError, match="Configuration file not found"):
        load_config(missing)


@pytest.mark.parametrize(
    "email_value",
    ["", "your.email@example.com", "not-an-email"],
)
def test_invalid_or_placeholder_email_raises(tmp_path: Path, email_value: str) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(f'ncbi:\n  email: "{email_value}"\n', encoding="utf-8")
    with pytest.raises(ConfigError, match="ncbi.email"):
        load_config(cfg)


def test_invalid_preferred_msa_tool_raises(tmp_path: Path) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        'ncbi:\n  email: "valid@lab.org"\nanalysis:\n  msa:\n    preferred_tool: "invalid_tool"\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="preferred_tool"):
        load_config(cfg)


def test_invalid_batch_size_raises(tmp_path: Path) -> None:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        'ncbi:\n  email: "valid@lab.org"\n  fetch_batch_size: 0\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="fetch_batch_size"):
        load_config(cfg)
