"""Step 1 tab: Search results metadata preview and selection table."""

from __future__ import annotations

import tkinter as tk
from collections.abc import Callable
from tkinter import ttk

from ncbi_protein_analyzer.models import ProteinSummary

CHECKED = "[x]"
UNCHECKED = "[ ]"


class SearchResultsTab(ttk.Frame):
    """Displays esummary hits with filtering, sorting, and checkbox-style selection."""

    def __init__(
        self,
        parent: tk.Misc,
        *,
        on_analyze_selected: Callable[[list[ProteinSummary]], None],
    ) -> None:
        super().__init__(parent, padding=8)
        self._on_analyze_selected = on_analyze_selected
        self._all_summaries: list[ProteinSummary] = []
        self._selected_uids: set[str] = set()
        self._visible_by_iid: dict[str, ProteinSummary] = {}
        self._sort_reverse: dict[str, bool] = {}

        self._filter_var = tk.StringVar()
        self._filter_var.trace_add("write", lambda *_: self._refresh_tree())
        self._summary_info_var = tk.StringVar(value="Run a search to preview NCBI protein records.")

        self._build_ui()

    def _build_ui(self) -> None:
        # Top controls bar
        ctrl = ttk.Frame(self)
        ctrl.pack(fill="x", pady=(0, 6))

        ttk.Label(ctrl, text="Filter table:").pack(side="left")
        filter_entry = ttk.Entry(ctrl, textvariable=self._filter_var, width=28)
        filter_entry.pack(side="left", padx=(4, 12))

        ttk.Button(ctrl, text="Select All Visible", command=self._select_all_visible).pack(
            side="left", padx=2
        )
        ttk.Button(ctrl, text="Select Top 20", command=lambda: self._select_top_n(20)).pack(
            side="left", padx=2
        )
        ttk.Button(ctrl, text="Clear Selection", command=self._clear_selection).pack(
            side="left", padx=2
        )

        self._analyze_btn = ttk.Button(
            ctrl,
            text="Download FASTA & Analyze (0 selected)",
            command=self._handle_analyze_click,
            state="disabled",
        )
        self._analyze_btn.pack(side="right")

        info_label = ttk.Label(self, textvariable=self._summary_info_var, foreground="#374151")
        info_label.pack(fill="x", pady=(0, 4))

        # Treeview + scrollbars
        tree_frame = ttk.Frame(self)
        tree_frame.pack(fill="both", expand=True)

        columns = ("selected", "accession", "title", "organism", "length", "updated")
        self._tree = ttk.Treeview(tree_frame, columns=columns, show="headings", selectmode="browse")

        headings = {
            "selected": ("Sel", 50, "center"),
            "accession": ("Accession", 130, "w"),
            "title": ("Title", 420, "w"),
            "organism": ("Organism", 190, "w"),
            "length": ("Length (aa)", 90, "e"),
            "updated": ("Updated", 100, "center"),
        }
        for col, (label, width, anchor) in headings.items():
            self._tree.heading(col, text=label, command=lambda c=col: self._sort_by(c))
            self._tree.column(col, width=width, anchor=anchor, stretch=(col == "title"))

        yscroll = ttk.Scrollbar(tree_frame, orient="vertical", command=self._tree.yview)
        xscroll = ttk.Scrollbar(tree_frame, orient="horizontal", command=self._tree.xview)
        self._tree.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)

        self._tree.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)

        self._tree.bind("<Button-1>", self._on_tree_click)
        self._tree.bind("<space>", self._on_space_toggle)

        hint = ttk.Label(
            self,
            text="Tip: Click any row (or press Space) to toggle selection, then click 'Download FASTA & Analyze'.",
            foreground="#6b7280",
        )
        hint.pack(fill="x", pady=(4, 0))

    def set_results(self, summaries: list[ProteinSummary], total_count: int, query_term: str) -> None:
        """Populate the preview table with new esummary search results."""
        self._all_summaries = list(summaries)
        # Pre-select up to 20 records for convenience
        self._selected_uids = {s.uid for s in summaries[:20]}
        self._summary_info_var.set(
            f"Showing {len(summaries)} of {total_count:,} NCBI hits for: {query_term}"
        )
        self._refresh_tree()

    def set_busy(self, busy: bool) -> None:
        if busy:
            self._analyze_btn.configure(state="disabled")
        else:
            self._update_analyze_button()

    def _filtered_summaries(self) -> list[ProteinSummary]:
        needle = self._filter_var.get().strip().lower()
        if not needle:
            return list(self._all_summaries)
        return [
            s
            for s in self._all_summaries
            if needle in s.accession.lower()
            or needle in s.title.lower()
            or needle in s.organism.lower()
        ]

    def _refresh_tree(self) -> None:
        self._tree.delete(*self._tree.get_children())
        self._visible_by_iid.clear()
        for summary in self._filtered_summaries():
            mark = CHECKED if summary.uid in self._selected_uids else UNCHECKED
            iid = self._tree.insert(
                "",
                "end",
                values=(
                    mark,
                    summary.accession,
                    summary.title,
                    summary.organism,
                    summary.length,
                    summary.update_date,
                ),
            )
            self._visible_by_iid[iid] = summary
        self._update_analyze_button()

    def _on_tree_click(self, event: tk.Event[tk.Misc]) -> None:
        region = self._tree.identify("region", event.x, event.y)
        if region != "cell":
            return
        iid = self._tree.identify_row(event.y)
        if iid:
            self._toggle_iid(iid)

    def _on_space_toggle(self, _event: tk.Event[tk.Misc]) -> None:
        for iid in self._tree.selection():
            self._toggle_iid(iid)

    def _toggle_iid(self, iid: str) -> None:
        summary = self._visible_by_iid.get(iid)
        if summary is None:
            return
        if summary.uid in self._selected_uids:
            self._selected_uids.remove(summary.uid)
            mark = UNCHECKED
        else:
            self._selected_uids.add(summary.uid)
            mark = CHECKED
        self._tree.set(iid, "selected", mark)
        self._update_analyze_button()

    def _select_all_visible(self) -> None:
        for iid, summary in self._visible_by_iid.items():
            self._selected_uids.add(summary.uid)
            self._tree.set(iid, "selected", CHECKED)
        self._update_analyze_button()

    def _select_top_n(self, n: int) -> None:
        self._selected_uids.clear()
        visible = list(self._visible_by_iid.items())
        for idx, (iid, summary) in enumerate(visible):
            if idx < n:
                self._selected_uids.add(summary.uid)
                self._tree.set(iid, "selected", CHECKED)
            else:
                self._tree.set(iid, "selected", UNCHECKED)
        self._update_analyze_button()

    def _clear_selection(self) -> None:
        self._selected_uids.clear()
        for iid in self._visible_by_iid:
            self._tree.set(iid, "selected", UNCHECKED)
        self._update_analyze_button()

    def _sort_by(self, column: str) -> None:
        reverse = self._sort_reverse.get(column, False)
        self._sort_reverse[column] = not reverse

        def _key(s: ProteinSummary) -> object:
            if column == "selected":
                return 0 if s.uid in self._selected_uids else 1
            if column == "length":
                return s.length
            if column == "accession":
                return s.accession.lower()
            if column == "title":
                return s.title.lower()
            if column == "organism":
                return s.organism.lower()
            return s.update_date

        self._all_summaries.sort(key=_key, reverse=reverse)
        self._refresh_tree()

    def _update_analyze_button(self) -> None:
        count = len(self._selected_uids)
        self._analyze_btn.configure(
            text=f"Download FASTA & Analyze ({count} selected)",
            state="normal" if count > 0 else "disabled",
        )

    def _handle_analyze_click(self) -> None:
        selected = [s for s in self._all_summaries if s.uid in self._selected_uids]
        if selected:
            self._on_analyze_selected(selected)
