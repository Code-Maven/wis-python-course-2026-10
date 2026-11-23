"""Biopython Entrez client for searching and downloading NCBI protein records."""

from __future__ import annotations

import io
import json
import re
import socket
import threading
import time
import urllib.error
from collections.abc import Callable, Sequence
from typing import Any

from Bio import Entrez, SeqIO

from ncbi_protein_analyzer.config import NcbiConfig
from ncbi_protein_analyzer.models import ProteinRecord, ProteinSummary, SearchQuery, SearchResult

ProgressCallback = Callable[[str, float | None], None]

_ORGANISM_BRACKET_RE = re.compile(r"\[([^\[\]]+)\]\s*$")
_RETRYABLE_HTTP_CODES = frozenset({429, 500, 502, 503, 504})


class NcbiClientError(Exception):
    """Raised when an NCBI E-utilities request fails."""


class CancelledError(Exception):
    """Raised when an operation is cancelled by the user."""


def build_entrez_term(query: SearchQuery) -> str:
    """Build an Entrez query string from structured search parameters."""
    base = query.term.strip()
    if not base:
        raise ValueError("Search term cannot be empty.")

    parts: list[str] = [f"({base})"]
    organism = query.organism.strip()
    if organism:
        parts.append(f'"{organism}"[Organism]')
    if query.refseq_only:
        parts.append("srcdb_refseq[PROP]")
    if query.min_length is not None or query.max_length is not None:
        low = query.min_length if query.min_length is not None else 1
        high = query.max_length if query.max_length is not None else 1000000
        if low > high:
            raise ValueError("Minimum length cannot be greater than maximum length.")
        parts.append(f"{low}:{high}[SLEN]")

    return " AND ".join(parts)


def _extract_organism_from_description(description: str) -> str:
    match = _ORGANISM_BRACKET_RE.search(description.strip())
    return match.group(1).strip() if match else ""


class NcbiClient:
    """Rate-limited, retrying wrapper around Bio.Entrez for NCBI Protein."""

    def __init__(
        self,
        config: NcbiConfig,
        *,
        entrez_module: Any = Entrez,
        sleep_fn: Callable[[float], None] = time.sleep,
        monotonic_fn: Callable[[], float] = time.monotonic,
        max_retries: int = 3,
    ) -> None:
        self.config = config
        self._entrez = entrez_module
        self._sleep = sleep_fn
        self._monotonic = monotonic_fn
        self._max_retries = max_retries
        self._last_request_ts: float = 0.0
        self._min_interval = 0.11 if config.api_key else 0.35

        self._entrez.email = config.email
        self._entrez.tool = config.tool
        self._entrez.api_key = config.api_key or None

    def _check_cancelled(self, cancel_event: threading.Event | None) -> None:
        if cancel_event is not None and cancel_event.is_set():
            raise CancelledError("Operation cancelled by user.")

    def _throttle(self) -> None:
        now = self._monotonic()
        if self._last_request_ts > 0.0:
            elapsed = now - self._last_request_ts
            remaining = self._min_interval - elapsed
            if remaining > 0:
                self._sleep(remaining)
        self._last_request_ts = self._monotonic()

    def _call_with_retry(
        self,
        fn: Callable[[], Any],
        *,
        cancel_event: threading.Event | None = None,
        action_name: str = "NCBI request",
    ) -> Any:
        previous_timeout = socket.getdefaulttimeout()
        socket.setdefaulttimeout(self.config.timeout_seconds)
        try:
            for attempt in range(1, self._max_retries + 1):
                self._check_cancelled(cancel_event)
                self._throttle()
                try:
                    return fn()
                except urllib.error.HTTPError as exc:
                    if exc.code in _RETRYABLE_HTTP_CODES and attempt < self._max_retries:
                        self._sleep(0.5 * (2 ** (attempt - 1)))
                        continue
                    raise NcbiClientError(
                        f"{action_name} failed with HTTP {exc.code}: {exc.reason}"
                    ) from exc
                except (urllib.error.URLError, TimeoutError, OSError) as exc:
                    if attempt < self._max_retries:
                        self._sleep(0.5 * (2 ** (attempt - 1)))
                        continue
                    raise NcbiClientError(f"{action_name} failed: {exc}") from exc
        finally:
            socket.setdefaulttimeout(previous_timeout)
        raise NcbiClientError(f"{action_name} failed after {self._max_retries} attempts.")

    def search_summaries(
        self,
        query: SearchQuery,
        *,
        cancel_event: threading.Event | None = None,
        progress_cb: ProgressCallback | None = None,
    ) -> SearchResult:
        """Search NCBI protein database and return esummary metadata records."""
        built_term = build_entrez_term(query)
        if progress_cb:
            progress_cb("Searching NCBI Protein database...", 0.15)

        def _do_esearch() -> dict[str, Any]:
            with self._entrez.esearch(
                db="protein",
                term=built_term,
                retmax=query.max_results,
            ) as handle:
                return self._entrez.read(handle)

        search_data = self._call_with_retry(
            _do_esearch,
            cancel_event=cancel_event,
            action_name="NCBI esearch",
        )
        total_count = int(search_data.get("Count", 0))
        id_list: list[str] = [str(uid) for uid in search_data.get("IdList", [])]

        if not id_list:
            if progress_cb:
                progress_cb("No matching proteins found.", 1.0)
            return SearchResult(query_term=built_term, total_count=total_count, summaries=[])

        summaries: list[ProteinSummary] = []
        batch_size = self.config.fetch_batch_size
        total_batches = (len(id_list) + batch_size - 1) // batch_size

        for batch_idx in range(total_batches):
            self._check_cancelled(cancel_event)
            start = batch_idx * batch_size
            chunk = id_list[start : start + batch_size]
            if progress_cb:
                frac = 0.25 + 0.7 * (batch_idx / max(total_batches, 1))
                progress_cb(
                    f"Loading metadata summaries ({start + 1}-{start + len(chunk)} of {len(id_list)})...",
                    frac,
                )

            def _do_esummary(ids: list[str] = chunk) -> dict[str, Any]:
                with self._entrez.esummary(
                    db="protein",
                    id=",".join(ids),
                    retmode="json",
                ) as handle:
                    raw = handle.read()
                    if isinstance(raw, bytes):
                        raw = raw.decode("utf-8", errors="replace")
                    return json.loads(raw)

            summary_payload = self._call_with_retry(
                _do_esummary,
                cancel_event=cancel_event,
                action_name="NCBI esummary",
            )
            result_block = summary_payload.get("result", {})
            for uid in result_block.get("uids", []):
                entry = result_block.get(str(uid), {})
                if not isinstance(entry, dict):
                    continue
                accession = (
                    str(entry.get("accessionversion") or entry.get("caption") or uid).strip()
                )
                title = str(entry.get("title", "")).strip()
                organism = str(entry.get("organism", "")).strip()
                if not organism and title:
                    organism = _extract_organism_from_description(title)
                try:
                    length = int(entry.get("slen", 0) or 0)
                except (TypeError, ValueError):
                    length = 0
                summaries.append(
                    ProteinSummary(
                        uid=str(uid),
                        accession=accession,
                        title=title,
                        organism=organism,
                        length=length,
                        create_date=str(entry.get("createdate", "") or ""),
                        update_date=str(entry.get("updatedate", "") or ""),
                        extra={
                            k: v
                            for k, v in entry.items()
                            if k
                            not in {
                                "uid",
                                "accessionversion",
                                "caption",
                                "title",
                                "organism",
                                "slen",
                                "createdate",
                                "updatedate",
                            }
                        },
                    )
                )

        if progress_cb:
            progress_cb(f"Loaded {len(summaries)} protein summaries.", 1.0)
        return SearchResult(
            query_term=built_term,
            total_count=total_count,
            summaries=summaries,
        )

    def fetch_fasta(
        self,
        summaries: Sequence[ProteinSummary],
        *,
        cancel_event: threading.Event | None = None,
        progress_cb: ProgressCallback | None = None,
    ) -> list[ProteinRecord]:
        """Download FASTA sequences for the selected protein summaries in batches."""
        if not summaries:
            return []

        summary_by_acc = {s.accession: s for s in summaries if s.accession}
        summary_by_acc_no_ver = {
            s.accession.split(".")[0]: s for s in summaries if s.accession
        }

        ids = [s.uid or s.accession for s in summaries]
        batch_size = self.config.fetch_batch_size
        total_batches = (len(ids) + batch_size - 1) // batch_size
        records: list[ProteinRecord] = []

        for batch_idx in range(total_batches):
            self._check_cancelled(cancel_event)
            start = batch_idx * batch_size
            chunk = ids[start : start + batch_size]
            if progress_cb:
                frac = batch_idx / max(total_batches, 1)
                progress_cb(
                    f"Downloading FASTA ({start + 1}-{start + len(chunk)} of {len(ids)})...",
                    frac,
                )

            def _do_efetch(batch_ids: list[str] = chunk) -> str:
                with self._entrez.efetch(
                    db="protein",
                    id=",".join(batch_ids),
                    rettype="fasta",
                    retmode="text",
                ) as handle:
                    raw = handle.read()
                    return raw.decode("utf-8", errors="replace") if isinstance(raw, bytes) else raw

            fasta_text = self._call_with_retry(
                _do_efetch,
                cancel_event=cancel_event,
                action_name="NCBI efetch (FASTA)",
            )
            for seq_record in SeqIO.parse(io.StringIO(fasta_text), "fasta"):
                acc = seq_record.id.strip()
                desc = seq_record.description.strip()
                if desc.startswith(acc):
                    desc_clean = desc[len(acc) :].strip()
                else:
                    desc_clean = desc

                matched_summary = (
                    summary_by_acc.get(acc)
                    or summary_by_acc_no_ver.get(acc.split(".")[0])
                )
                organism = (
                    matched_summary.organism
                    if matched_summary and matched_summary.organism
                    else _extract_organism_from_description(desc_clean)
                )
                seq_str = "".join(str(seq_record.seq).upper().split())
                records.append(
                    ProteinRecord(
                        accession=acc,
                        description=desc_clean or (matched_summary.title if matched_summary else ""),
                        organism=organism,
                        sequence=seq_str,
                    )
                )

        if progress_cb:
            progress_cb(f"Downloaded {len(records)} FASTA sequences.", 1.0)
        return records
