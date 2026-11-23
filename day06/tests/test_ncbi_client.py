"""Unit tests for the mocked Biopython Entrez NCBI client."""

from __future__ import annotations

import io
import json
import threading
import urllib.error
from typing import Any

import pytest

from ncbi_protein_analyzer.config import NcbiConfig
from ncbi_protein_analyzer.models import ProteinSummary, SearchQuery
from ncbi_protein_analyzer.ncbi.client import (
    CancelledError,
    NcbiClient,
    NcbiClientError,
    build_entrez_term,
)


class _FakeEntrez:
    def __init__(
        self,
        *,
        esearch_data: dict[str, Any] | None = None,
        esummary_payload: dict[str, Any] | None = None,
        fasta_batches: list[str] | None = None,
        fail_first_esearch: int = 0,
    ) -> None:
        self.email = ""
        self.tool = ""
        self.api_key: str | None = None
        self._esearch_data = esearch_data or {"Count": "0", "IdList": []}
        self._esummary_payload = esummary_payload or {"result": {"uids": []}}
        self._fasta_batches = list(fasta_batches or [])
        self._fail_first_esearch = fail_first_esearch
        self.esearch_calls: list[dict[str, Any]] = []
        self.esummary_calls: list[dict[str, Any]] = []
        self.efetch_calls: list[dict[str, Any]] = []

    def esearch(self, **kwargs: Any) -> io.StringIO:
        self.esearch_calls.append(kwargs)
        if self._fail_first_esearch > 0:
            self._fail_first_esearch -= 1
            raise urllib.error.HTTPError(
                url="https://eutils.ncbi.nlm.nih.gov",
                code=429,
                msg="Too Many Requests",
                hdrs=None,  # type: ignore[arg-type]
                fp=None,
            )
        return io.StringIO("mock-xml")

    def read(self, _handle: Any) -> dict[str, Any]:
        return self._esearch_data

    def esummary(self, **kwargs: Any) -> io.StringIO:
        self.esummary_calls.append(kwargs)
        return io.StringIO(json.dumps(self._esummary_payload))

    def efetch(self, **kwargs: Any) -> io.StringIO:
        self.efetch_calls.append(kwargs)
        text = self._fasta_batches.pop(0) if self._fasta_batches else ""
        return io.StringIO(text)


def test_build_entrez_term_with_all_filters() -> None:
    q = SearchQuery(
        term="insulin",
        organism="Homo sapiens",
        refseq_only=True,
        min_length=50,
        max_length=200,
    )
    term = build_entrez_term(q)
    assert term == '(insulin) AND "Homo sapiens"[Organism] AND srcdb_refseq[PROP] AND 50:200[SLEN]'


def test_build_entrez_term_rejects_empty_and_bad_bounds() -> None:
    with pytest.raises(ValueError, match="empty"):
        build_entrez_term(SearchQuery(term="   "))
    with pytest.raises(ValueError, match="Minimum length"):
        build_entrez_term(SearchQuery(term="kinase", min_length=300, max_length=100))


def test_search_summaries_parses_esummary_json() -> None:
    fake = _FakeEntrez(
        esearch_data={"Count": "125", "IdList": ["1001", "1002"]},
        esummary_payload={
            "result": {
                "uids": ["1001", "1002"],
                "1001": {
                    "uid": "1001",
                    "accessionversion": "NP_000198.1",
                    "title": "insulin preproprotein [Homo sapiens]",
                    "organism": "Homo sapiens",
                    "slen": 110,
                    "updatedate": "2025/01/10",
                },
                "1002": {
                    "uid": "1002",
                    "caption": "NP_032412",
                    "title": "insulin-1 [Mus musculus]",
                    "slen": 108,
                },
            }
        },
    )
    client = NcbiClient(
        NcbiConfig(email="user@lab.org", fetch_batch_size=10),
        entrez_module=fake,
        sleep_fn=lambda _s: None,
    )
    res = client.search_summaries(SearchQuery(term="insulin", max_results=10))
    assert res.total_count == 125
    assert len(res.summaries) == 2
    assert res.summaries[0].accession == "NP_000198.1"
    assert res.summaries[0].length == 110
    assert res.summaries[1].accession == "NP_032412"
    assert res.summaries[1].organism == "Mus musculus"


def test_search_summaries_empty_results() -> None:
    fake = _FakeEntrez(esearch_data={"Count": "0", "IdList": []})
    client = NcbiClient(
        NcbiConfig(email="user@lab.org"),
        entrez_module=fake,
        sleep_fn=lambda _s: None,
    )
    res = client.search_summaries(SearchQuery(term="nonexistent_protein_xyz"))
    assert res.total_count == 0
    assert res.summaries == []
    assert fake.esummary_calls == []


def test_retry_on_http_429_succeeds() -> None:
    fake = _FakeEntrez(
        esearch_data={"Count": "1", "IdList": ["1001"]},
        esummary_payload={
            "result": {
                "uids": ["1001"],
                "1001": {"accessionversion": "NP_1.1", "title": "test [Homo sapiens]", "slen": 50},
            }
        },
        fail_first_esearch=1,
    )
    sleeps: list[float] = []
    client = NcbiClient(
        NcbiConfig(email="user@lab.org"),
        entrez_module=fake,
        sleep_fn=sleeps.append,
    )
    res = client.search_summaries(SearchQuery(term="test"))
    assert len(res.summaries) == 1
    assert len(fake.esearch_calls) == 2


def test_fetch_fasta_batches_and_maps_organisms(sample_summaries: list[ProteinSummary]) -> None:
    fasta_1 = ">NP_000198.1 insulin preproprotein [Homo sapiens]\nMALWMRLLPL\n"
    fasta_2 = ">NP_032412.3 insulin-1 preproprotein [Mus musculus]\nMALLVHFLPL\n"
    fake = _FakeEntrez(fasta_batches=[fasta_1, fasta_2])
    client = NcbiClient(
        NcbiConfig(email="user@lab.org", fetch_batch_size=1),
        entrez_module=fake,
        sleep_fn=lambda _s: None,
    )
    records = client.fetch_fasta(sample_summaries)
    assert len(fake.efetch_calls) == 2
    assert len(records) == 2
    assert records[0].accession == "NP_000198.1"
    assert records[0].sequence == "MALWMRLLPL"
    assert records[0].organism == "Homo sapiens"
    assert records[1].accession == "NP_032412.3"


def test_cancelled_before_request_raises() -> None:
    fake = _FakeEntrez()
    client = NcbiClient(
        NcbiConfig(email="user@lab.org"),
        entrez_module=fake,
        sleep_fn=lambda _s: None,
    )
    cancel_event = threading.Event()
    cancel_event.set()
    with pytest.raises(CancelledError):
        client.search_summaries(SearchQuery(term="insulin"), cancel_event=cancel_event)


def test_non_retryable_http_error_raises_client_error() -> None:
    class _FailingEntrez(_FakeEntrez):
        def esearch(self, **kwargs: Any) -> io.StringIO:
            raise urllib.error.HTTPError(
                url="https://eutils.ncbi.nlm.nih.gov",
                code=400,
                msg="Bad Request",
                hdrs=None,  # type: ignore[arg-type]
                fp=None,
            )

    client = NcbiClient(
        NcbiConfig(email="user@lab.org"),
        entrez_module=_FailingEntrez(),
        sleep_fn=lambda _s: None,
    )
    with pytest.raises(NcbiClientError, match="HTTP 400"):
        client.search_summaries(SearchQuery(term="bad"))
