"""Step 2c tab: Sequence comparison (MSA / pairwise fallback), deduplication, and alignment viewer."""

from __future__ import annotations

import tkinter as tk
from collections.abc import Callable
from tkinter import filedialog, messagebox, ttk

from ncbi_protein_analyzer.analysis.msa import DetectedMsaTool
from ncbi_protein_analyzer.export import export_fasta
from ncbi_protein_analyzer.models import ComparisonResult


def _format_alignment_blocks(result: ComparisonResult, block_width: int = 60) -> str:
    """Format aligned sequences + consensus into Clustal-style wrapped blocks."""
    aligned = result.aligned_records
    if not aligned:
        return ""
    total_len = len(aligned[0].sequence)
    name_width = max(len(r.accession) for r in aligned)
    name_width = max(name_width, len("Consensus")) + 2

    lines: list[str] = []
    for start in range(0, total_len, block_width):
        end = min(start + block_width, total_len)
        lines.append(f"Columns {start + 1}-{end}:")
        for rec in aligned:
            chunk = rec.sequence[start:end]
            lines.append(f"{rec.accession:<{name_width}}{chunk}")
        if result.consensus:
            cons_chunk = result.consensus[start:end]
            lines.append(f"{'Consensus':<{name_width}}{cons_chunk}")
        lines.append("")
    return "\n".join(lines)


def _format_pairwise_summary(result: ComparisonResult) -> str:
    """Format top pairwise identity pairs when running in pairwise fallback mode."""
    labels = result.labels
    matrix = result.identity_matrix
    n = len(labels)
    if n < 2:
        return "Not enough unique sequences to compare."

    pairs: list[tuple[float, str, str]] = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((matrix[i][j], labels[i], labels[j]))
    pairs.sort(key=lambda item: item[0], reverse=True)

    lines = [
        f"Pairwise Global Alignment Identity Summary ({len(pairs)} unique pairs):",
        f"{'Rank':<6}{'Protein A':<20}{'Protein B':<20}{'Identity (%)':>14}",
        "-" * 62,
    ]
    for rank, (ident, a, b) in enumerate(pairs[:100], start=1):
        lines.append(f"{rank:<6}{a:<20}{b:<20}{ident:>13.2f}%")
    if len(pairs) > 100:
        lines.append(f"... ({len(pairs) - 100} additional pairs omitted; inspect heatmap in Charts tab)")
    return "\n".join(lines)


class CompareTab(ttk.Frame):
    """Displays MSA/pairwise status, deduplication clusters, warnings, and alignment viewer."""

    def __init__(
        self,
        parent: tk.Misc,
        *,
        detected_tool: DetectedMsaTool | None,
        default_mode: str,
        on_rerun_comparison: Callable[[str], None],
    ) -> None:
        super().__init__(parent, padding=8)
        self._on_rerun_comparison = on_rerun_comparison
        self._comparison: ComparisonResult | None = None
        self._has_records = False

        tool_desc = (
            f"Detected MSA tool: {detected_tool.name} ({detected_tool.executable})"
            if detected_tool
            else "Detected MSA tool: None on PATH (will use Biopython pairwise fallback)"
        )
        self._tool_status_var = tk.StringVar(value=tool_desc)
        self._mode_var = tk.StringVar(value=default_mode)
        self._summary_var = tk.StringVar(value="Download and analyze proteins to compare sequences.")
        self._build_ui()

    def _build_ui(self) -> None:
        top = ttk.Frame(self)
        top.pack(fill="x", pady=(0, 6))

        ttk.Label(top, textvariable=self._tool_status_var, foreground="#1f2937").pack(side="left")

        self._export_aln_btn = ttk.Button(
            top,
            text="Export Aligned FASTA...",
            command=self._export_aligned_fasta,
            state="disabled",
        )
        self._export_aln_btn.pack(side="right", padx=(6, 0))

        self._rerun_btn = ttk.Button(
            top,
            text="Re-run Comparison",
            command=self._handle_rerun,
            state="disabled",
        )
        self._rerun_btn.pack(side="right", padx=(6, 0))

        mode_combo = ttk.Combobox(
            top,
            textvariable=self._mode_var,
            state="readonly",
            width=12,
            values=["auto", "clustalo", "muscle", "pairwise"],
        )
        mode_combo.pack(side="right", padx=(4, 0))
        ttk.Label(top, text="Mode:").pack(side="right")

        ttk.Label(self, textvariable=self._summary_var, font=("TkDefaultFont", 9, "bold")).pack(
            fill="x", pady=(0, 4)
        )

        # Monospace scrollable viewer
        text_frame = ttk.Frame(self)
        text_frame.pack(fill="both", expand=True)

        self._text = tk.Text(
            text_frame,
            wrap="none",
            font=("TkFixedFont", 9),
            state="disabled",
            background="#f9fafb",
        )
        yscroll = ttk.Scrollbar(text_frame, orient="vertical", command=self._text.yview)
        xscroll = ttk.Scrollbar(text_frame, orient="horizontal", command=self._text.xview)
        self._text.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)

        self._text.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        text_frame.rowconfigure(0, weight=1)
        text_frame.columnconfigure(0, weight=1)

    def set_busy(self, busy: bool) -> None:
        if busy:
            self._rerun_btn.configure(state="disabled")
        else:
            self._rerun_btn.configure(state="normal" if self._has_records else "disabled")

    def set_results(self, comparison: ComparisonResult, *, has_records: bool = True) -> None:
        self._comparison = comparison
        self._has_records = has_records
        dedup = comparison.deduplication
        total_in = len(dedup.unique_records) + dedup.duplicate_count
        self._summary_var.set(
            f"Method used: {comparison.method_used}  |  "
            f"Input sequences: {total_in}  |  "
            f"Unique: {len(dedup.unique_records)}  |  "
            f"Exact duplicates grouped: {dedup.duplicate_count}"
        )

        self._rerun_btn.configure(state="normal" if has_records else "disabled")
        self._export_aln_btn.configure(
            state="normal" if comparison.aligned_records else "disabled"
        )

        sections: list[str] = []
        if comparison.warnings:
            sections.append("Notes / Warnings:")
            for w in comparison.warnings:
                sections.append(f"  - {w}")
            sections.append("")

        dup_clusters = {rep: members for rep, members in dedup.clusters.items() if len(members) > 1}
        if dup_clusters:
            sections.append("Exact Duplicate Clusters:")
            for rep, members in dup_clusters.items():
                sections.append(f"  {rep} == {', '.join(members[1:])}")
            sections.append("")

        if comparison.aligned_records:
            sections.append(
                f"Multiple Sequence Alignment ({len(comparison.aligned_records)} sequences, "
                f"{len(comparison.aligned_records[0].sequence)} columns):"
            )
            sections.append("")
            sections.append(_format_alignment_blocks(comparison))
        elif comparison.method_used == "pairwise":
            sections.append(_format_pairwise_summary(comparison))
        else:
            sections.append("No sequence comparison available.")

        self._set_text("\n".join(sections))

    def _set_text(self, content: str) -> None:
        self._text.configure(state="normal")
        self._text.delete("1.0", "end")
        self._text.insert("1.0", content)
        self._text.configure(state="disabled")

    def _handle_rerun(self) -> None:
        if self._has_records:
            self._on_rerun_comparison(self._mode_var.get())

    def _export_aligned_fasta(self) -> None:
        if not self._comparison or not self._comparison.aligned_records:
            return
        path = filedialog.asksaveasfilename(
            title="Export Aligned Sequences as FASTA",
            defaultextension=".fasta",
            filetypes=[("Aligned FASTA", "*.fasta *.aln *.fa"), ("All files", "*.*")],
            initialfile="proteins_aligned.fasta",
        )
        if not path:
            return
        try:
            saved = export_fasta(self._comparison.aligned_records, path)
            messagebox.showinfo("Export Complete", f"Saved aligned FASTA to:\n{saved}")
        except OSError as exc:
            messagebox.showerror("Export Failed", f"Could not write aligned FASTA:\n{exc}")
