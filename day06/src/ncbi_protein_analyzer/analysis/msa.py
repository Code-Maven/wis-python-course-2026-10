"""External Multiple Sequence Alignment (Clustal Omega / MUSCLE) integration."""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ncbi_protein_analyzer.config import MsaConfig
from ncbi_protein_analyzer.models import ProteinRecord

SubprocessRunner = Callable[..., subprocess.CompletedProcess[str]]


class MsaError(Exception):
    """Raised when an external MSA binary cannot be executed or produces invalid output."""


@dataclass(slots=True)
class DetectedMsaTool:
    """Resolved external MSA executable."""

    name: str  # "clustalo" | "muscle"
    executable: str


def _resolve_executable(explicit_path: str, binary_name: str) -> str | None:
    if explicit_path:
        path_obj = Path(explicit_path)
        if path_obj.is_file():
            return str(path_obj)
        found = shutil.which(explicit_path)
        return found
    return shutil.which(binary_name)


def detect_msa_tool(
    config: MsaConfig,
    *,
    preferred_override: str | None = None,
) -> DetectedMsaTool | None:
    """Locate an external MSA tool according to config and PATH."""
    preferred = (preferred_override or config.preferred_tool or "auto").strip().lower()
    if preferred == "pairwise":
        return None

    if preferred == "clustalo":
        exe = _resolve_executable(config.clustalo_path, "clustalo")
        return DetectedMsaTool("clustalo", exe) if exe else None

    if preferred == "muscle":
        exe = _resolve_executable(config.muscle_path, "muscle")
        return DetectedMsaTool("muscle", exe) if exe else None

    # "auto": prefer clustalo, then muscle
    clustalo_exe = _resolve_executable(config.clustalo_path, "clustalo")
    if clustalo_exe:
        return DetectedMsaTool("clustalo", clustalo_exe)

    muscle_exe = _resolve_executable(config.muscle_path, "muscle")
    if muscle_exe:
        return DetectedMsaTool("muscle", muscle_exe)

    return None


def detect_muscle_major_version(
    executable: str,
    *,
    runner: SubprocessRunner = subprocess.run,
) -> int:
    """Detect whether a MUSCLE executable is v5+ or v3."""
    try:
        proc = runner(
            [executable, "-version"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return 5

    combined = f"{proc.stdout}\n{proc.stderr}"
    match = re.search(r"muscle\s+v?(\d+)\.", combined, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return 5


def run_msa(
    records: Sequence[ProteinRecord],
    tool: DetectedMsaTool,
    *,
    timeout_seconds: int = 120,
    runner: SubprocessRunner = subprocess.run,
) -> list[ProteinRecord]:
    """Run Clustal Omega or MUSCLE on records and return aligned ProteinRecords."""
    if len(records) < 2:
        raise MsaError("At least two sequences are required for Multiple Sequence Alignment.")

    with tempfile.TemporaryDirectory(prefix="ncbi_msa_") as tmpdir:
        tmp_path = Path(tmpdir)
        input_fasta = tmp_path / "input.fasta"
        output_fasta = tmp_path / "aligned.fasta"

        id_map: dict[str, ProteinRecord] = {}
        seq_records: list[SeqRecord] = []
        for idx, rec in enumerate(records):
            safe_id = f"seq_{idx}"
            id_map[safe_id] = rec
            seq_records.append(SeqRecord(Seq(rec.sequence), id=safe_id, description=""))

        SeqIO.write(seq_records, input_fasta, "fasta")

        if tool.name == "clustalo":
            cmd = [
                tool.executable,
                "-i",
                str(input_fasta),
                "-o",
                str(output_fasta),
                "--outfmt=fasta",
                "--force",
            ]
        elif tool.name == "muscle":
            major = detect_muscle_major_version(tool.executable, runner=runner)
            if major >= 5:
                cmd = [tool.executable, "-align", str(input_fasta), "-output", str(output_fasta)]
            else:
                cmd = [tool.executable, "-in", str(input_fasta), "-out", str(output_fasta)]
        else:
            raise MsaError(f"Unsupported MSA tool: {tool.name}")

        try:
            proc: Any = runner(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise MsaError(f"{tool.name} timed out after {timeout_seconds}s.") from exc
        except OSError as exc:
            raise MsaError(f"Failed to execute {tool.name} ('{tool.executable}'): {exc}") from exc

        if proc.returncode != 0:
            err_msg = (proc.stderr or proc.stdout or "").strip()
            raise MsaError(f"{tool.name} exited with code {proc.returncode}: {err_msg}")

        if not output_fasta.exists() or output_fasta.stat().st_size == 0:
            raise MsaError(f"{tool.name} did not produce an output alignment file.")

        aligned_by_id: dict[str, str] = {}
        for aligned_seq in SeqIO.parse(output_fasta, "fasta"):
            aligned_by_id[aligned_seq.id] = str(aligned_seq.seq).upper()

        if len(aligned_by_id) != len(records):
            raise MsaError(
                f"{tool.name} returned {len(aligned_by_id)} aligned sequences; expected {len(records)}."
            )

        ordered_aligned: list[ProteinRecord] = []
        for idx, rec in enumerate(records):
            safe_id = f"seq_{idx}"
            if safe_id not in aligned_by_id:
                raise MsaError(f"Missing aligned sequence for '{rec.accession}'.")
            ordered_aligned.append(
                ProteinRecord(
                    accession=rec.accession,
                    description=rec.description,
                    organism=rec.organism,
                    sequence=aligned_by_id[safe_id],
                )
            )
        return ordered_aligned
