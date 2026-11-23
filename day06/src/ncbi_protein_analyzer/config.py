"""YAML configuration loading and validation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

DEFAULT_CONFIG_PATH = Path("config.yaml")
VALID_MSA_TOOLS = frozenset({"auto", "clustalo", "muscle", "pairwise"})
PLACEHOLDER_EMAILS = frozenset(
    {
        "your.email@example.com",
        "example@example.com",
        "user@example.com",
        "email@example.com",
    }
)


class ConfigError(Exception):
    """Raised when config.yaml is missing or invalid."""


@dataclass(slots=True)
class NcbiConfig:
    """NCBI E-utilities connection and rate-limit settings."""

    email: str
    api_key: str = ""
    tool: str = "ncbi_protein_analyzer"
    max_results_default: int = 50
    fetch_batch_size: int = 50
    timeout_seconds: int = 30


@dataclass(slots=True)
class MsaConfig:
    """External MSA tool configuration."""

    preferred_tool: str = "auto"
    clustalo_path: str = ""
    muscle_path: str = ""
    timeout_seconds: int = 120


@dataclass(slots=True)
class AnalysisConfig:
    """Analysis guardrails and MSA configuration."""

    max_pairwise_sequences: int = 50
    msa: MsaConfig = field(default_factory=MsaConfig)


@dataclass(slots=True)
class AppConfig:
    """Top-level application configuration."""

    ncbi: NcbiConfig
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    source_path: Path | None = None


def _as_dict(value: Any, section: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError(f"Configuration section '{section}' must be a YAML mapping.")
    return value


def _positive_int(value: Any, field_name: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigError(f"Configuration field '{field_name}' must be an integer >= {minimum}.")
    return value


def load_config(path: str | Path = DEFAULT_CONFIG_PATH) -> AppConfig:
    """Load and validate application configuration from a YAML file."""
    config_path = Path(path)
    if not config_path.exists():
        raise ConfigError(
            f"Configuration file not found: {config_path}\n"
            "Copy 'config.yaml.example' to 'config.yaml' and set your NCBI email address."
        )

    try:
        raw_text = config_path.read_text(encoding="utf-8")
        data = yaml.safe_load(raw_text)
    except OSError as exc:
        raise ConfigError(f"Could not read configuration file '{config_path}': {exc}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in '{config_path}': {exc}") from exc

    if not isinstance(data, dict):
        raise ConfigError(f"Configuration file '{config_path}' must contain a YAML mapping.")

    ncbi_raw = _as_dict(data.get("ncbi"), "ncbi")
    email = str(ncbi_raw.get("email", "") or "").strip()
    if not email or email.lower() in PLACEHOLDER_EMAILS or "@" not in email:
        raise ConfigError(
            "A valid 'ncbi.email' address is required by NCBI E-utilities.\n"
            f"Please edit '{config_path}' and replace the placeholder email with your real email."
        )

    api_key = str(ncbi_raw.get("api_key", "") or "").strip()
    tool = str(ncbi_raw.get("tool", "ncbi_protein_analyzer") or "ncbi_protein_analyzer").strip()
    max_results_default = _positive_int(
        ncbi_raw.get("max_results_default", 50), "ncbi.max_results_default", minimum=1
    )
    fetch_batch_size = _positive_int(
        ncbi_raw.get("fetch_batch_size", 50), "ncbi.fetch_batch_size", minimum=1
    )
    timeout_seconds = _positive_int(
        ncbi_raw.get("timeout_seconds", 30), "ncbi.timeout_seconds", minimum=1
    )

    analysis_raw = _as_dict(data.get("analysis"), "analysis")
    max_pairwise_sequences = _positive_int(
        analysis_raw.get("max_pairwise_sequences", 50),
        "analysis.max_pairwise_sequences",
        minimum=2,
    )

    msa_raw = _as_dict(analysis_raw.get("msa"), "analysis.msa")
    preferred_tool = str(msa_raw.get("preferred_tool", "auto") or "auto").strip().lower()
    if preferred_tool not in VALID_MSA_TOOLS:
        valid_list = ", ".join(sorted(VALID_MSA_TOOLS))
        raise ConfigError(
            f"Invalid 'analysis.msa.preferred_tool': '{preferred_tool}'. "
            f"Must be one of: {valid_list}."
        )

    msa_timeout = _positive_int(
        msa_raw.get("timeout_seconds", 120), "analysis.msa.timeout_seconds", minimum=1
    )

    return AppConfig(
        ncbi=NcbiConfig(
            email=email,
            api_key=api_key,
            tool=tool,
            max_results_default=max_results_default,
            fetch_batch_size=fetch_batch_size,
            timeout_seconds=timeout_seconds,
        ),
        analysis=AnalysisConfig(
            max_pairwise_sequences=max_pairwise_sequences,
            msa=MsaConfig(
                preferred_tool=preferred_tool,
                clustalo_path=str(msa_raw.get("clustalo_path", "") or "").strip(),
                muscle_path=str(msa_raw.get("muscle_path", "") or "").strip(),
                timeout_seconds=msa_timeout,
            ),
        ),
        source_path=config_path,
    )
