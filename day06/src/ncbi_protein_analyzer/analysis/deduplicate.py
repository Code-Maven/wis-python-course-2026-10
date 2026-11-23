"""Exact sequence deduplication."""

from __future__ import annotations

from collections.abc import Sequence

from ncbi_protein_analyzer.models import DeduplicationResult, ProteinRecord


def deduplicate_records(records: Sequence[ProteinRecord]) -> DeduplicationResult:
    """Group proteins with identical normalized amino acid sequences."""
    seq_to_rep: dict[str, str] = {}
    unique_records: list[ProteinRecord] = []
    clusters: dict[str, list[str]] = {}

    for rec in records:
        norm_seq = "".join(rec.sequence.upper().split())
        if norm_seq not in seq_to_rep:
            seq_to_rep[norm_seq] = rec.accession
            unique_records.append(rec)
            clusters[rec.accession] = [rec.accession]
        else:
            rep = seq_to_rep[norm_seq]
            clusters[rep].append(rec.accession)

    duplicate_count = len(records) - len(unique_records)
    return DeduplicationResult(
        unique_records=unique_records,
        clusters=clusters,
        duplicate_count=duplicate_count,
    )
