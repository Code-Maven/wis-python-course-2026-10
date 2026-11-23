"""Sequence comparison orchestrator (deduplication, external MSA, and pairwise fallback)."""

from __future__ import annotations

import subprocess
import threading
from collections import Counter
from collections.abc import Callable, Sequence

from Bio.Align import PairwiseAligner, substitution_matrices

from ncbi_protein_analyzer.analysis.deduplicate import deduplicate_records
from ncbi_protein_analyzer.analysis.msa import (
    MsaError,
    SubprocessRunner,
    detect_msa_tool,
    run_msa,
)
from ncbi_protein_analyzer.config import AnalysisConfig
from ncbi_protein_analyzer.models import ComparisonResult, ProteinRecord
from ncbi_protein_analyzer.ncbi.client import CancelledError

ProgressCallback = Callable[[str, float | None], None]


def _sanitize_for_matrix(seq: str, valid_alphabet: frozenset[str]) -> str:
    """Replace residues missing from the substitution matrix (e.g., U, O, J, *) with X."""
    fallback = "X" if "X" in valid_alphabet else "A"
    return "".join(ch if ch in valid_alphabet else fallback for ch in seq.upper())


def compute_msa_identity_matrix(aligned_records: Sequence[ProteinRecord]) -> list[list[float]]:
    """Compute an N x N pairwise percent-identity matrix from aligned sequences."""
    n = len(aligned_records)
    matrix: list[list[float]] = [[100.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    seqs = [r.sequence.upper() for r in aligned_records]

    for i in range(n):
        for j in range(i + 1, n):
            s1, s2 = seqs[i], seqs[j]
            compared = 0
            matches = 0
            for a, b in zip(s1, s2):
                if a == "-" and b == "-":
                    continue
                compared += 1
                if a == b:
                    matches += 1
            ident = (matches / compared * 100.0) if compared > 0 else 0.0
            matrix[i][j] = ident
            matrix[j][i] = ident
    return matrix


def compute_consensus(aligned_records: Sequence[ProteinRecord]) -> str:
    """Build a simple majority-rule consensus string from aligned sequences."""
    if not aligned_records:
        return ""
    alignment_len = min(len(r.sequence) for r in aligned_records)
    consensus_chars: list[str] = []
    for col in range(alignment_len):
        column_chars = [r.sequence[col].upper() for r in aligned_records]
        non_gap = [c for c in column_chars if c != "-"]
        if not non_gap or len(non_gap) < len(column_chars) / 2:
            consensus_chars.append("-")
            continue
        most_common_char, _ = Counter(non_gap).most_common(1)[0]
        consensus_chars.append(most_common_char)
    return "".join(consensus_chars)


def compute_pairwise_identity_matrix(
    records: Sequence[ProteinRecord],
    *,
    cancel_event: threading.Event | None = None,
    progress_cb: ProgressCallback | None = None,
) -> list[list[float]]:
    """Compute an N x N pairwise percent-identity matrix using Bio.Align.PairwiseAligner."""
    n = len(records)
    matrix: list[list[float]] = [[100.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    if n < 2:
        return matrix

    blosum62 = substitution_matrices.load("BLOSUM62")
    alphabet = frozenset(blosum62.alphabet)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    sanitized = [_sanitize_for_matrix(r.sequence, alphabet) for r in records]
    total_pairs = n * (n - 1) // 2
    completed = 0

    for i in range(n):
        for j in range(i + 1, n):
            if cancel_event is not None and cancel_event.is_set():
                raise CancelledError("Pairwise comparison cancelled by user.")

            s1, s2 = sanitized[i], sanitized[j]
            if not s1 or not s2:
                ident = 0.0
            else:
                alignment = aligner.align(s1, s2)[0]
                counts = alignment.counts()
                aligned_cols = int(counts.identities + counts.mismatches + counts.gaps)
                ident = (counts.identities / aligned_cols * 100.0) if aligned_cols > 0 else 0.0

            matrix[i][j] = ident
            matrix[j][i] = ident
            completed += 1
            if progress_cb and (completed % 5 == 0 or completed == total_pairs):
                progress_cb(
                    f"Pairwise alignment ({completed}/{total_pairs} pairs)...",
                    completed / max(total_pairs, 1),
                )

    return matrix


def compare_records(
    records: Sequence[ProteinRecord],
    config: AnalysisConfig,
    *,
    mode_override: str | None = None,
    cancel_event: threading.Event | None = None,
    progress_cb: ProgressCallback | None = None,
    runner: SubprocessRunner = subprocess.run,
) -> ComparisonResult:
    """Deduplicate records and run external MSA or pairwise fallback."""
    dedup = deduplicate_records(records)
    unique = dedup.unique_records
    warnings: list[str] = []

    if dedup.duplicate_count > 0:
        warnings.append(
            f"Grouped {dedup.duplicate_count} exact duplicate sequence(s); "
            f"comparing {len(unique)} unique sequence(s)."
        )

    if len(unique) < 2:
        warnings.append("At least 2 unique sequences are required for sequence comparison.")
        labels = [r.accession for r in unique]
        matrix = [[100.0]] if len(unique) == 1 else []
        return ComparisonResult(
            method_used="none",
            deduplication=dedup,
            labels=labels,
            identity_matrix=matrix,
            aligned_records=[],
            consensus="",
            warnings=warnings,
        )

    tool = detect_msa_tool(config.msa, preferred_override=mode_override)
    if tool is not None:
        if progress_cb:
            progress_cb(f"Running {tool.name} on {len(unique)} unique sequences...", 0.2)
        try:
            aligned = run_msa(
                unique,
                tool,
                timeout_seconds=config.msa.timeout_seconds,
                runner=runner,
            )
            if progress_cb:
                progress_cb("Computing identity matrix & consensus from MSA...", 0.85)
            matrix = compute_msa_identity_matrix(aligned)
            consensus = compute_consensus(aligned)
            if progress_cb:
                progress_cb(f"MSA completed with {tool.name}.", 1.0)
            return ComparisonResult(
                method_used=tool.name,
                deduplication=dedup,
                labels=[r.accession for r in aligned],
                identity_matrix=matrix,
                aligned_records=aligned,
                consensus=consensus,
                warnings=warnings,
            )
        except MsaError as exc:
            warnings.append(
                f"External MSA ({tool.name}) failed ({exc}); falling back to Biopython pairwise alignment."
            )
    else:
        requested = (mode_override or config.msa.preferred_tool or "auto").lower()
        if requested != "pairwise":
            warnings.append(
                "No external MSA binary (clustalo/muscle) found on PATH; using Biopython pairwise fallback."
            )

    pairwise_records = list(unique)
    if len(pairwise_records) > config.max_pairwise_sequences:
        warnings.append(
            f"Pairwise fallback capped at first {config.max_pairwise_sequences} of {len(unique)} "
            "unique sequences to keep runtime responsive. Install Clustal Omega or MUSCLE for larger sets."
        )
        pairwise_records = pairwise_records[: config.max_pairwise_sequences]

    matrix = compute_pairwise_identity_matrix(
        pairwise_records,
        cancel_event=cancel_event,
        progress_cb=progress_cb,
    )
    return ComparisonResult(
        method_used="pairwise",
        deduplication=dedup,
        labels=[r.accession for r in pairwise_records],
        identity_matrix=matrix,
        aligned_records=[],
        consensus="",
        warnings=warnings,
    )
