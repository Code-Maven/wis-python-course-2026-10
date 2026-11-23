"""Step 2b tab: Embedded Matplotlib visualizations and PNG export."""

from __future__ import annotations

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from ncbi_protein_analyzer.export import export_all_figures_png, export_figure_png
from ncbi_protein_analyzer.models import ComparisonResult, DatasetSummary, ProteinStats
from ncbi_protein_analyzer.plotting import (
    create_aa_composition_chart,
    create_identity_heatmap,
    create_length_histogram,
)


class ChartsTab(ttk.Frame):
    """Hosts embedded Matplotlib plots with selector and PNG export actions."""

    def __init__(self, parent: tk.Misc) -> None:
        super().__init__(parent, padding=8)
        self._stats: list[ProteinStats] = []
        self._summary: DatasetSummary | None = None
        self._comparison: ComparisonResult | None = None
        self._figures: dict[str, Figure] = {}
        self._canvas: FigureCanvasTkAgg | None = None
        self._toolbar: NavigationToolbar2Tk | None = None

        self._chart_choice = tk.StringVar(value="length_histogram")
        self._aa_target_var = tk.StringVar(value="Dataset Average")
        self._build_ui()

    def _build_ui(self) -> None:
        ctrl = ttk.Frame(self)
        ctrl.pack(fill="x", pady=(0, 6))

        ttk.Radiobutton(
            ctrl,
            text="Length Histogram",
            value="length_histogram",
            variable=self._chart_choice,
            command=self._render_active_chart,
        ).pack(side="left", padx=(0, 8))

        ttk.Radiobutton(
            ctrl,
            text="AA Composition",
            value="aa_composition",
            variable=self._chart_choice,
            command=self._render_active_chart,
        ).pack(side="left", padx=(0, 8))

        ttk.Radiobutton(
            ctrl,
            text="Identity Heatmap",
            value="identity_heatmap",
            variable=self._chart_choice,
            command=self._render_active_chart,
        ).pack(side="left", padx=(0, 12))

        ttk.Label(ctrl, text="AA Source:").pack(side="left")
        self._aa_combo = ttk.Combobox(
            ctrl,
            textvariable=self._aa_target_var,
            state="readonly",
            width=20,
            font="TkTextFont",
            values=["Dataset Average"],
        )
        self._aa_combo.pack(side="left", padx=(4, 12))
        self._aa_combo.bind("<<ComboboxSelected>>", self._on_aa_target_changed)

        self._save_all_btn = ttk.Button(
            ctrl, text="Export All PNGs...", command=self._export_all_pngs, state="disabled"
        )
        self._save_all_btn.pack(side="right", padx=(4, 0))

        self._save_one_btn = ttk.Button(
            ctrl, text="Save Chart PNG...", command=self._export_active_png, state="disabled"
        )
        self._save_one_btn.pack(side="right")

        self._plot_host = ttk.Frame(self)
        self._plot_host.pack(fill="both", expand=True)

        self._placeholder = ttk.Label(
            self._plot_host,
            text="Download and analyze proteins to generate charts.",
            foreground="#6b7280",
        )
        self._placeholder.pack(expand=True)

    def set_results(
        self,
        stats: list[ProteinStats],
        summary: DatasetSummary,
        comparison: ComparisonResult,
    ) -> None:
        self._stats = list(stats)
        self._summary = summary
        self._comparison = comparison

        options = ["Dataset Average"] + [s.accession for s in stats]
        self._aa_combo.configure(values=options)
        self._aa_target_var.set("Dataset Average")

        self._rebuild_figures()
        self._render_active_chart()

    def update_comparison(self, comparison: ComparisonResult) -> None:
        """Refresh only the identity heatmap after a re-run of sequence comparison."""
        self._comparison = comparison
        self._figures["identity_heatmap"] = create_identity_heatmap(comparison)
        if self._chart_choice.get() == "identity_heatmap":
            self._render_active_chart()

    def _rebuild_figures(self) -> None:
        if self._summary is None or self._comparison is None:
            return
        self._figures["length_histogram"] = create_length_histogram(self._stats)
        self._figures["aa_composition"] = self._build_aa_figure()
        self._figures["identity_heatmap"] = create_identity_heatmap(self._comparison)
        self._save_one_btn.configure(state="normal")
        self._save_all_btn.configure(state="normal")

    def _build_aa_figure(self) -> Figure:
        target = self._aa_target_var.get()
        if target == "Dataset Average" or not self._stats:
            comp = self._summary.aggregate_aa_composition if self._summary else {}
            return create_aa_composition_chart(comp, title="Amino Acid Composition (Dataset Average)")

        matched = next((s for s in self._stats if s.accession == target), None)
        if matched is None:
            comp = self._summary.aggregate_aa_composition if self._summary else {}
            return create_aa_composition_chart(comp, title="Amino Acid Composition (Dataset Average)")

        suffix = " [Approx]" if matched.is_approximate else ""
        return create_aa_composition_chart(
            matched.aa_composition,
            title=f"Amino Acid Composition — {matched.accession}{suffix}",
        )

    def _on_aa_target_changed(self, _event: tk.Event[tk.Misc]) -> None:
        if self._summary is None:
            return
        self._figures["aa_composition"] = self._build_aa_figure()
        self._chart_choice.set("aa_composition")
        self._render_active_chart()

    def _render_active_chart(self) -> None:
        key = self._chart_choice.get()
        fig = self._figures.get(key)
        if fig is None:
            return

        for child in self._plot_host.winfo_children():
            child.destroy()

        canvas = FigureCanvasTkAgg(fig, master=self._plot_host)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, self._plot_host, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="bottom", fill="x")
        canvas.get_tk_widget().pack(side="top", fill="both", expand=True)

        self._canvas = canvas
        self._toolbar = toolbar

    def _export_active_png(self) -> None:
        key = self._chart_choice.get()
        fig = self._figures.get(key)
        if fig is None:
            return
        path = filedialog.asksaveasfilename(
            title="Save Chart as PNG",
            defaultextension=".png",
            filetypes=[("PNG images", "*.png"), ("All files", "*.*")],
            initialfile=f"{key}.png",
        )
        if not path:
            return
        try:
            saved = export_figure_png(fig, path)
            messagebox.showinfo("Export Complete", f"Saved chart PNG to:\n{saved}")
        except OSError as exc:
            messagebox.showerror("Export Failed", f"Could not save PNG:\n{exc}")

    def _export_all_pngs(self) -> None:
        if not self._figures:
            return
        directory = filedialog.askdirectory(title="Choose Folder to Export All Charts (PNG)")
        if not directory:
            return
        try:
            saved = export_all_figures_png(self._figures, directory)
            files_list = "\n".join(str(p.name) for p in saved)
            messagebox.showinfo("Export Complete", f"Saved {len(saved)} PNG files to {directory}:\n{files_list}")
        except OSError as exc:
            messagebox.showerror("Export Failed", f"Could not export PNG files:\n{exc}")
