"""Shared domain data models for the NCBI Protein Analyzer."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

STANDARD_AAS: tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWY")
STANDARD_AA_SET: frozenset[str] = frozenset(STANDARD_AAS)


@dataclass(slots=True)
class SearchQuery:
    """Parameters for an NCBI Protein Entrez search."""

    term: str
    organism: str = ""
    refseq_only: bool = False
    min_length: int | None = None
    max_length: int | None = None
    max_results: int = 50


@dataclass(slots=True)
class ProteinSummary:
    """Lightweight metadata returned from NCBI esummary (Step 1)."""

    uid: str
    accession: str
    title: str
    organism: str
    length: int
    create_date: str = ""
    update_date: str = ""
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class SearchResult:
    """Result container from an NCBI esearch + esummary query."""

    query_term: str
    total_count: int
    summaries: list[ProteinSummary]


@dataclass(slots=True)
class ProteinRecord:
    """Downloaded protein sequence and metadata (Step 2)."""

    accession: str
    description: str
    organism: str
    sequence: str


@dataclass(slots=True)
class ProteinStats:
    """Physicochemical and composition statistics for a single protein."""

    accession: str
    description: str
    organism: str
    length: int
    standard_aa_count: int
    ambiguous_aa_count: int
    ambiguous_aa_pct: float
    is_approximate: bool
    molecular_weight_da: float | None
    isoelectric_point: float | None
    gravy: float | None
    aa_composition: dict[str, float]
    warning: str = ""


@dataclass(slots=True)
class DatasetSummary:
    """Aggregate summary statistics across an analyzed set of proteins."""

    count: int
    min_length: int
    median_length: float
    mean_length: float
    max_length: int
    approximate_count: int
    total_ambiguous_aas: int
    mean_molecular_weight_da: float | None
    mean_isoelectric_point: float | None
    mean_gravy: float | None
    aggregate_aa_composition: dict[str, float]


@dataclass(slots=True)
class DeduplicationResult:
    """Result of grouping exact duplicate protein sequences."""

    unique_records: list[ProteinRecord]
    clusters: dict[str, list[str]]
    duplicate_count: int


@dataclass(slots=True)
class ComparisonResult:
    """Result of MSA or pairwise sequence comparison."""

    method_used: str
    deduplication: DeduplicationResult
    labels: list[str]
    identity_matrix: list[list[float]]
    aligned_records: list[ProteinRecord] = field(default_factory=list)
    consensus: str = ""
    warnings: list[str] = field(default_factory=list)
