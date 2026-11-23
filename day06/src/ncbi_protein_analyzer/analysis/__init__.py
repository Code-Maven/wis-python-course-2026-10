"""Protein analysis package (stats, deduplication, MSA, pairwise comparison)."""

from .compare import compare_records, compute_consensus, compute_msa_identity_matrix, compute_pairwise_identity_matrix
from .deduplicate import deduplicate_records
from .msa import DetectedMsaTool, MsaError, detect_msa_tool, run_msa
from .stats import analyze_dataset, analyze_record

__all__ = [
    "DetectedMsaTool",
    "MsaError",
    "analyze_dataset",
    "analyze_record",
    "compare_records",
    "compute_consensus",
    "compute_msa_identity_matrix",
    "compute_pairwise_identity_matrix",
    "deduplicate_records",
    "detect_msa_tool",
    "run_msa",
]
