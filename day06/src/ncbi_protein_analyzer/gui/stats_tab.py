"""Step 2a tab: Physicochemical sequence statistics table and CSV/FASTA export."""

from __future__ import annotations

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from ncbi_protein_analyzer.export import export_fasta, export_stats_csv
from ncbi_protein_analyzer.models import DatasetSummary, ProteinRecord, ProteinStats


class StatsTab(ttk.Frame):
    """Displays per-protein ProtParam statistics, summary KPIs, and export buttons."""

    def __init__(self, parent: tk.Misc) -> None:
        super().__init__(parent, padding=8)
        self._stats: list[ProteinStats] = []
        self._records: list[ProteinRecord] = []
        self._stats_by_iid: dict[str, ProteinStats] = {}
        self._sort_reverse: dict[str, bool] = {}

        self._kpi_var = tk.StringVar(value="No proteins analyzed yet.")
        self._detail_var = tk.StringVar(value="Select a protein row to view notes or ambiguous residue details.")
        self._build_ui()

    def _build_ui(self) -> None:
        # Summary KPI + export controls
        top = ttk.Frame(self)
        top.pack(fill="x", pady=(0, 6))

        ttk.Label(top, textvariable=self._kpi_var, font=("TkDefaultFont", 9, "bold")).pack(
            side="left", fill="x", expand=True
        )
        self._export_csv_btn = ttk.Button(
            top, text="Export CSV...", command=self._export_csv, state="disabled"
        )
        self._export_csv_btn.pack(side="right", padx=(4, 0))

        self._export_fasta_btn = ttk.Button(
            top, text="Export FASTA...", command=self._export_fasta, state="disabled"
        )
        self._export_fasta_btn.pack(side="right")

        # Stats Treeview
        tree_frame = ttk.Frame(self)
        tree_frame.pack(fill="both", expand=True)

        columns = (
            "accession",
            "length",
            "ambiguous",
            "approx",
            "mw_kda",
            "pi",
            "gravy",
            "organism",
            "description",
        )
        self._tree = ttk.Treeview(tree_frame, columns=columns, show="headings", selectmode="browse")
        self._tree.tag_configure("approx", background="#fef3c7")

        col_specs = {
            "accession": ("Accession", 125, "w"),
            "length": ("Length", 70, "e"),
            "ambiguous": ("Ambig. AAs", 80, "e"),
            "approx": ("Approx?", 75, "center"),
            "mw_kda": ("MW (kDa)", 85, "e"),
            "pi": ("pI", 65, "e"),
            "gravy": ("GRAVY", 75, "e"),
            "organism": ("Organism", 165, "w"),
            "description": ("Description", 290, "w"),
        }
        for col, (label, width, anchor) in col_specs.items():
            self._tree.heading(col, text=label, command=lambda c=col: self._sort_by(c))
            self._tree.column(col, width=width, anchor=anchor, stretch=(col == "description"))

        yscroll = ttk.Scrollbar(tree_frame, orient="vertical", command=self._tree.yview)
        xscroll = ttk.Scrollbar(tree_frame, orient="horizontal", command=self._tree.xview)
        self._tree.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)

        self._tree.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)

        self._tree.bind("<<TreeviewSelect>>", self._on_row_select)

        detail_lbl = ttk.Label(self, textvariable=self._detail_var, foreground="#4b5563")
        detail_lbl.pack(fill="x", pady=(6, 0))

    def set_results(
        self,
        records: list[ProteinRecord],
        stats: list[ProteinStats],
        summary: DatasetSummary,
        duplicate_count: int = 0,
    ) -> None:
        self._records = list(records)
        self._stats = list(stats)

        mw_text = (
            f"{summary.mean_molecular_weight_da / 1000.0:.2f} kDa"
            if summary.mean_molecular_weight_da is not None
            else "N/A"
        )
        pi_text = f"{summary.mean_isoelectric_point:.2f}" if summary.mean_isoelectric_point is not None else "N/A"
        gravy_text = f"{summary.mean_gravy:.3f}" if summary.mean_gravy is not None else "N/A"

        self._kpi_var.set(
            f"Proteins: {summary.count}  |  Length: {summary.min_length}–{summary.max_length} aa "
            f"(mean {summary.mean_length:.1f}, median {summary.median_length:.1f})  |  "
            f"Mean MW: {mw_text}  |  Mean pI: {pi_text}  |  Mean GRAVY: {gravy_text}  |  "
            f"Approx: {summary.approximate_count}  |  Exact duplicates: {duplicate_count}"
        )

        has_data = len(stats) > 0
        self._export_csv_btn.configure(state="normal" if has_data else "disabled")
        self._export_fasta_btn.configure(state="normal" if records else "disabled")
        self._refresh_tree()

    def _refresh_tree(self) -> None:
        self._tree.delete(*self._tree.get_children())
        self._stats_by_iid.clear()
        for item in self._stats:
            prefix = "~" if item.is_approximate else ""
            mw_str = (
                f"{prefix}{item.molecular_weight_da / 1000.0:.2f}"
                if item.molecular_weight_da is not None
                else "N/A"
            )
            pi_str = (
                f"{prefix}{item.isoelectric_point:.2f}"
                if item.isoelectric_point is not None
                else "N/A"
            )
            gravy_str = f"{prefix}{item.gravy:.3f}" if item.gravy is not None else "N/A"
            ambig_str = (
                f"{item.ambiguous_aa_count} ({item.ambiguous_aa_pct:.1f}%)"
                if item.ambiguous_aa_count > 0
                else "0"
            )
            tags = ("approx",) if item.is_approximate else ()
            iid = self._tree.insert(
                "",
                "end",
                values=(
                    item.accession,
                    item.length,
                    ambig_str,
                    "Yes (~)" if item.is_approximate else "No",
                    mw_str,
                    pi_str,
                    gravy_str,
                    item.organism,
                    item.description,
                ),
                tags=tags,
            )
            self._stats_by_iid[iid] = item

    def _on_row_select(self, _event: tk.Event[tk.Misc]) -> None:
        sel = self._tree.selection()
        if not sel:
            return
        item = self._stats_by_iid.get(sel[0])
        if item is None:
            return
        if item.warning:
            self._detail_var.set(f"{item.accession}: {item.warning}")
        else:
            self._detail_var.set(
                f"{item.accession} ({item.organism or 'Unknown organism'}): "
                f"{item.length} standard amino acids, no ambiguous residues."
            )

    def _sort_by(self, column: str) -> None:
        reverse = self._sort_reverse.get(column, False)
        self._sort_reverse[column] = not reverse

        def _key(item: ProteinStats) -> object:
            if column == "accession":
                return item.accession.lower()
            if column == "length":
                return item.length
            if column == "ambiguous":
                return item.ambiguous_aa_count
            if column == "approx":
                return 1 if item.is_approximate else 0
            if column == "mw_kda":
                return item.molecular_weight_da if item.molecular_weight_da is not None else -1.0
            if column == "pi":
                return item.isoelectric_point if item.isoelectric_point is not None else -1.0
            if column == "gravy":
                return item.gravy if item.gravy is not None else -999.0
            if column == "organism":
                return item.organism.lower()
            return item.description.lower()

        self._stats.sort(key=_key, reverse=reverse)
        self._refresh_tree()

    def _export_csv(self) -> None:
        if not self._stats:
            return
        path = filedialog.asksaveasfilename(
            title="Export Protein Statistics as CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialfile="protein_stats.csv",
        )
        if not path:
            return
        try:
            saved = export_stats_csv(self._stats, path)
            messagebox.showinfo("Export Complete", f"Saved statistics CSV to:\n{saved}")
        except OSError as exc:
            messagebox.showerror("Export Failed", f"Could not write CSV file:\n{exc}")

    def _export_fasta(self) -> None:
        if not self._records:
            return
        path = filedialog.asksaveasfilename(
            title="Export Downloaded Sequences as FASTA",
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")],
            initialfile="proteins.fasta",
        )
        if not path:
            return
        try:
            saved = export_fasta(self._records, path)
            messagebox.showinfo("Export Complete", f"Saved FASTA sequences to:\n{saved}")
        except OSError as exc:
            messagebox.showerror("Export Failed", f"Could not write FASTA file:\n{exc}")
