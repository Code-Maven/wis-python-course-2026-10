"""Main Tkinter desktop window orchestrating search, download, analysis, charts, and exports."""

from __future__ import annotations

import threading
import tkinter as tk
from tkinter import messagebox, ttk

from ncbi_protein_analyzer.analysis import analyze_dataset, compare_records, detect_msa_tool
from ncbi_protein_analyzer.config import AppConfig
from ncbi_protein_analyzer.gui.charts_tab import ChartsTab
from ncbi_protein_analyzer.gui.compare_tab import CompareTab
from ncbi_protein_analyzer.gui.search_tab import SearchResultsTab
from ncbi_protein_analyzer.gui.stats_tab import StatsTab
from ncbi_protein_analyzer.gui.worker import BackgroundRunner, ProgressFn
from ncbi_protein_analyzer.models import (
    ComparisonResult,
    DatasetSummary,
    ProteinRecord,
    ProteinStats,
    ProteinSummary,
    SearchQuery,
    SearchResult,
)
from ncbi_protein_analyzer.ncbi import NcbiClient


class NcbiProteinAnalyzerApp(tk.Tk):
    """Main Tkinter window for the NCBI Protein Search & Analysis application."""

    def __init__(self, config: AppConfig) -> None:
        super().__init__()
        self.config_data = config
        self.client = NcbiClient(config.ncbi)
        self.runner = BackgroundRunner(self)

        self._downloaded_records: list[ProteinRecord] = []
        self._detected_msa = detect_msa_tool(config.analysis.msa)

        self.title("NCBI Protein Search & Analyzer")
        self.geometry("1160x760")
        self.minsize(920, 620)

        self._term_var = tk.StringVar()
        self._organism_var = tk.StringVar()
        self._refseq_var = tk.BooleanVar(value=False)
        self._min_len_var = tk.StringVar()
        self._max_len_var = tk.StringVar()
        self._max_results_var = tk.IntVar(value=config.ncbi.max_results_default)
        self._status_var = tk.StringVar(value="Ready.")

        self._build_ui()

    def _build_ui(self) -> None:
        # Top search bar
        search_frame = ttk.LabelFrame(self, text="NCBI Protein Search", padding=8)
        search_frame.pack(fill="x", padx=8, pady=(8, 4))

        row1 = ttk.Frame(search_frame)
        row1.pack(fill="x")

        ttk.Label(row1, text="Query:").pack(side="left")
        term_entry = ttk.Entry(row1, textvariable=self._term_var, width=34)
        term_entry.pack(side="left", padx=(4, 12), fill="x", expand=True)
        term_entry.bind("<Return>", lambda _e: self._start_search())
        term_entry.focus_set()

        ttk.Label(row1, text="Organism:").pack(side="left")
        org_entry = ttk.Entry(row1, textvariable=self._organism_var, width=22)
        org_entry.pack(side="left", padx=(4, 12))
        org_entry.bind("<Return>", lambda _e: self._start_search())

        ttk.Label(row1, text="Max results:").pack(side="left")
        max_spin = ttk.Spinbox(row1, from_=1, to=500, textvariable=self._max_results_var, width=6)
        max_spin.pack(side="left", padx=(4, 12))

        self._search_btn = ttk.Button(row1, text="Search NCBI", command=self._start_search)
        self._search_btn.pack(side="left", padx=(0, 4))

        self._cancel_btn = ttk.Button(
            row1, text="Cancel", command=self._cancel_task, state="disabled"
        )
        self._cancel_btn.pack(side="left")

        row2 = ttk.Frame(search_frame)
        row2.pack(fill="x", pady=(6, 0))

        ttk.Checkbutton(row2, text="RefSeq only", variable=self._refseq_var).pack(
            side="left", padx=(0, 14)
        )
        ttk.Label(row2, text="Length filter (aa):").pack(side="left")
        ttk.Entry(row2, textvariable=self._min_len_var, width=8).pack(side="left", padx=(4, 2))
        ttk.Label(row2, text="to").pack(side="left", padx=2)
        ttk.Entry(row2, textvariable=self._max_len_var, width=8).pack(side="left", padx=(2, 12))

        config_info = f"Entrez email: {self.config_data.ncbi.email}"
        if self.config_data.ncbi.api_key:
            config_info += "  (API key active)"
        ttk.Label(row2, text=config_info, foreground="#6b7280").pack(side="right")

        # Main notebook tabs
        self._notebook = ttk.Notebook(self)
        self._notebook.pack(fill="both", expand=True, padx=8, pady=4)

        self._search_tab = SearchResultsTab(
            self._notebook, on_analyze_selected=self._start_download_and_analyze
        )
        self._stats_tab = StatsTab(self._notebook)
        self._charts_tab = ChartsTab(self._notebook)
        self._compare_tab = CompareTab(
            self._notebook,
            detected_tool=self._detected_msa,
            default_mode=self.config_data.analysis.msa.preferred_tool,
            on_rerun_comparison=self._start_rerun_comparison,
        )

        self._notebook.add(self._search_tab, text="1. Search Results")
        self._notebook.add(self._stats_tab, text="2. Sequence Stats")
        self._notebook.add(self._charts_tab, text="3. Charts")
        self._notebook.add(self._compare_tab, text="4. Sequence Comparison")

        # Bottom status & progress bar
        status_frame = ttk.Frame(self, padding=(8, 4))
        status_frame.pack(fill="x")

        ttk.Label(status_frame, textvariable=self._status_var).pack(side="left", fill="x", expand=True)
        self._progress = ttk.Progressbar(status_frame, mode="determinate", length=220, maximum=100)
        self._progress.pack(side="right")

    def _set_busy(self, busy: bool, status_msg: str | None = None) -> None:
        self._search_btn.configure(state="disabled" if busy else "normal")
        self._cancel_btn.configure(state="normal" if busy else "disabled")
        self._search_tab.set_busy(busy)
        self._compare_tab.set_busy(busy)
        if status_msg is not None:
            self._status_var.set(status_msg)
        if not busy:
            self._progress.stop()
            self._progress.configure(mode="determinate", value=0)

    def _on_progress(self, message: str, fraction: float | None) -> None:
        self._status_var.set(message)
        if fraction is None:
            if str(self._progress["mode"]) != "indeterminate":
                self._progress.configure(mode="indeterminate")
                self._progress.start(12)
        else:
            if str(self._progress["mode"]) != "determinate":
                self._progress.stop()
                self._progress.configure(mode="determinate")
            self._progress["value"] = max(0.0, min(100.0, fraction * 100.0))

    def _parse_optional_int(self, raw: str, label: str) -> int | None:
        text = raw.strip()
        if not text:
            return None
        try:
            val = int(text)
        except ValueError as exc:
            raise ValueError(f"{label} must be a positive integer.") from exc
        if val < 1:
            raise ValueError(f"{label} must be >= 1.")
        return val

    def _build_query_from_inputs(self) -> SearchQuery:
        term = self._term_var.get().strip()
        if not term:
            raise ValueError("Please enter a protein search query (e.g., 'insulin' or 'hemoglobin').")
        min_len = self._parse_optional_int(self._min_len_var.get(), "Minimum length")
        max_len = self._parse_optional_int(self._max_len_var.get(), "Maximum length")
        if min_len is not None and max_len is not None and min_len > max_len:
            raise ValueError("Minimum length cannot be greater than maximum length.")
        try:
            max_results = int(self._max_results_var.get())
        except (tk.TclError, ValueError) as exc:
            raise ValueError("Max results must be an integer between 1 and 500.") from exc
        if not (1 <= max_results <= 500):
            raise ValueError("Max results must be between 1 and 500.")

        return SearchQuery(
            term=term,
            organism=self._organism_var.get().strip(),
            refseq_only=self._refseq_var.get(),
            min_length=min_len,
            max_length=max_len,
            max_results=max_results,
        )

    def _start_search(self) -> None:
        if self.runner.is_running:
            return
        try:
            query = self._build_query_from_inputs()
        except ValueError as exc:
            messagebox.showwarning("Invalid Search Input", str(exc))
            return

        self._set_busy(True, "Starting NCBI search...")

        def _task(cancel_event: threading.Event, progress_cb: ProgressFn) -> SearchResult:
            return self.client.search_summaries(
                query, cancel_event=cancel_event, progress_cb=progress_cb
            )

        def _on_done(result: SearchResult) -> None:
            self._set_busy(
                False,
                f"Search complete: loaded {len(result.summaries)} of {result.total_count:,} NCBI hits.",
            )
            self._search_tab.set_results(result.summaries, result.total_count, result.query_term)
            self._notebook.select(self._search_tab)
            if not result.summaries:
                messagebox.showinfo(
                    "No Results",
                    f"No proteins matched your query:\n{result.query_term}",
                )

        self.runner.start(
            _task,
            on_progress=self._on_progress,
            on_success=_on_done,
            on_error=self._on_error,
            on_cancel=self._on_cancel,
        )

    def _start_download_and_analyze(self, selected: list[ProteinSummary]) -> None:
        if self.runner.is_running:
            return
        self._set_busy(True, f"Downloading FASTA for {len(selected)} protein(s)...")

        def _task(
            cancel_event: threading.Event, progress_cb: ProgressFn
        ) -> tuple[list[ProteinRecord], list[ProteinStats], DatasetSummary, ComparisonResult]:
            def _fetch_progress(msg: str, frac: float | None) -> None:
                scaled = (frac * 0.5) if frac is not None else None
                progress_cb(msg, scaled)

            records = self.client.fetch_fasta(
                selected, cancel_event=cancel_event, progress_cb=_fetch_progress
            )
            progress_cb(f"Computing physicochemical stats for {len(records)} proteins...", 0.58)
            stats, summary = analyze_dataset(records)

            def _cmp_progress(msg: str, frac: float | None) -> None:
                scaled = (0.62 + frac * 0.36) if frac is not None else None
                progress_cb(msg, scaled)

            comparison = compare_records(
                records,
                self.config_data.analysis,
                cancel_event=cancel_event,
                progress_cb=_cmp_progress,
            )
            progress_cb("Analysis complete.", 1.0)
            return records, stats, summary, comparison

        def _on_done(
            payload: tuple[
                list[ProteinRecord], list[ProteinStats], DatasetSummary, ComparisonResult
            ],
        ) -> None:
            records, stats, summary, comparison = payload
            self._downloaded_records = records
            self._stats_tab.set_results(
                records, stats, summary, duplicate_count=comparison.deduplication.duplicate_count
            )
            self._charts_tab.set_results(stats, summary, comparison)
            self._compare_tab.set_results(comparison, has_records=bool(records))
            self._set_busy(
                False,
                f"Analyzed {len(records)} protein(s) (comparison: {comparison.method_used}).",
            )
            self._notebook.select(self._stats_tab)

        self.runner.start(
            _task,
            on_progress=self._on_progress,
            on_success=_on_done,
            on_error=self._on_error,
            on_cancel=self._on_cancel,
        )

    def _start_rerun_comparison(self, mode_override: str) -> None:
        if self.runner.is_running or not self._downloaded_records:
            return
        records = list(self._downloaded_records)
        self._set_busy(True, f"Re-running sequence comparison (mode: {mode_override})...")

        def _task(cancel_event: threading.Event, progress_cb: ProgressFn) -> ComparisonResult:
            return compare_records(
                records,
                self.config_data.analysis,
                mode_override=mode_override,
                cancel_event=cancel_event,
                progress_cb=progress_cb,
            )

        def _on_done(comparison: ComparisonResult) -> None:
            self._compare_tab.set_results(comparison, has_records=True)
            self._charts_tab.update_comparison(comparison)
            self._set_busy(False, f"Sequence comparison complete ({comparison.method_used}).")

        self.runner.start(
            _task,
            on_progress=self._on_progress,
            on_success=_on_done,
            on_error=self._on_error,
            on_cancel=self._on_cancel,
        )

    def _cancel_task(self) -> None:
        self.runner.cancel()
        self._status_var.set("Cancelling...")

    def _on_cancel(self, message: str) -> None:
        self._set_busy(False, message or "Operation cancelled.")

    def _on_error(self, exc: Exception) -> None:
        self._set_busy(False, f"Error: {exc}")
        messagebox.showerror("Operation Failed", str(exc))
