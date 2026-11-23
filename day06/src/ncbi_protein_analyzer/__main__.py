"""CLI entrypoint for launching the NCBI Protein Search & Analyzer GUI."""

from __future__ import annotations

import argparse
import sys
import tkinter as tk
from pathlib import Path
from tkinter import messagebox

from ncbi_protein_analyzer.config import ConfigError, DEFAULT_CONFIG_PATH, load_config
from ncbi_protein_analyzer.gui import NcbiProteinAnalyzerApp


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="ncbi-protein-analyzer",
        description="Search NCBI Protein, download FASTA sequences, and analyze them in a desktop GUI.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help=f"Path to YAML configuration file (default: {DEFAULT_CONFIG_PATH})",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        config = load_config(args.config)
    except ConfigError as exc:
        print(f"Configuration error:\n{exc}", file=sys.stderr)
        try:
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("Configuration Error", str(exc))
            root.destroy()
        except tk.TclError:
            pass
        return 1

    app = NcbiProteinAnalyzerApp(config)
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
