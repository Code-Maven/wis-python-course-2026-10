"""Physicochemical and amino acid composition analysis using Biopython ProtParam."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from statistics import mean, median

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from ncbi_protein_analyzer.models import (
    STANDARD_AA_SET,
    STANDARD_AAS,
    DatasetSummary,
    ProteinRecord,
    ProteinStats,
)


def analyze_record(record: ProteinRecord) -> ProteinStats:
    """Compute physicochemical stats for a single protein, handling ambiguous AAs safely."""
    seq = "".join(record.sequence.upper().split())
    length = len(seq)
    clean_chars = [ch for ch in seq if ch in STANDARD_AA_SET]
    clean_seq = "".join(clean_chars)
    standard_count = len(clean_seq)
    ambiguous_count = length - standard_count
    ambiguous_pct = (ambiguous_count / length * 100.0) if length > 0 else 0.0
    is_approximate = ambiguous_count > 0

    counts = Counter(clean_chars)
    aa_composition = {
        aa: (counts.get(aa, 0) / standard_count * 100.0) if standard_count > 0 else 0.0
        for aa in STANDARD_AAS
    }

    if standard_count == 0:
        warning = (
            "Sequence is empty."
            if length == 0
            else f"Sequence contains only non-standard/ambiguous residues ({ambiguous_count}); MW/pI/GRAVY unavailable."
        )
        return ProteinStats(
            accession=record.accession,
            description=record.description,
            organism=record.organism,
            length=length,
            standard_aa_count=0,
            ambiguous_aa_count=ambiguous_count,
            ambiguous_aa_pct=ambiguous_pct,
            is_approximate=True,
            molecular_weight_da=None,
            isoelectric_point=None,
            gravy=None,
            aa_composition=aa_composition,
            warning=warning,
        )

    analyser = ProteinAnalysis(clean_seq)
    mw = float(analyser.molecular_weight())
    pi = float(analyser.isoelectric_point())
    gravy = float(analyser.gravy())
    warning = (
        f"Approximate: excluded {ambiguous_count} ambiguous/non-standard residue(s) ({ambiguous_pct:.1f}%)."
        if is_approximate
        else ""
    )

    return ProteinStats(
        accession=record.accession,
        description=record.description,
        organism=record.organism,
        length=length,
        standard_aa_count=standard_count,
        ambiguous_aa_count=ambiguous_count,
        ambiguous_aa_pct=ambiguous_pct,
        is_approximate=is_approximate,
        molecular_weight_da=mw,
        isoelectric_point=pi,
        gravy=gravy,
        aa_composition=aa_composition,
        warning=warning,
    )


def analyze_dataset(
    records: Sequence[ProteinRecord],
) -> tuple[list[ProteinStats], DatasetSummary]:
    """Analyze a collection of protein records and compute dataset-level aggregates."""
    per_protein = [analyze_record(rec) for rec in records]
    if not per_protein:
        empty_comp = {aa: 0.0 for aa in STANDARD_AAS}
        return [], DatasetSummary(
            count=0,
            min_length=0,
            median_length=0.0,
            mean_length=0.0,
            max_length=0,
            approximate_count=0,
            total_ambiguous_aas=0,
            mean_molecular_weight_da=None,
            mean_isoelectric_point=None,
            mean_gravy=None,
            aggregate_aa_composition=empty_comp,
        )

    lengths = [s.length for s in per_protein]
    mws = [s.molecular_weight_da for s in per_protein if s.molecular_weight_da is not None]
    pis = [s.isoelectric_point for s in per_protein if s.isoelectric_point is not None]
    gravies = [s.gravy for s in per_protein if s.gravy is not None]

    total_aa_counts: Counter[str] = Counter()
    total_standard = 0
    for rec in records:
        for ch in rec.sequence.upper():
            if ch in STANDARD_AA_SET:
                total_aa_counts[ch] += 1
                total_standard += 1

    aggregate_aa_comp = {
        aa: (total_aa_counts.get(aa, 0) / total_standard * 100.0) if total_standard > 0 else 0.0
        for aa in STANDARD_AAS
    }

    summary = DatasetSummary(
        count=len(per_protein),
        min_length=min(lengths),
        median_length=float(median(lengths)),
        mean_length=float(mean(lengths)),
        max_length=max(lengths),
        approximate_count=sum(1 for s in per_protein if s.is_approximate),
        total_ambiguous_aas=sum(s.ambiguous_aa_count for s in per_protein),
        mean_molecular_weight_da=float(mean(mws)) if mws else None,
        mean_isoelectric_point=float(mean(pis)) if pis else None,
        mean_gravy=float(mean(gravies)) if gravies else None,
        aggregate_aa_composition=aggregate_aa_comp,
    )
    return per_protein, summary
