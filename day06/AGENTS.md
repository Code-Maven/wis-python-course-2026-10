# AGENTS.md

Minimal guide for coding agents working in this repository. See `SPEC.md` for full application details.

## Commands (`uv`)

- Install/sync dependencies: `uv sync`
- Run tests: `uv run pytest`
- Launch GUI: `uv run ncbi-protein-analyzer`

## Architecture Rules

- **Keep core logic GUI-free**: `src/ncbi_protein_analyzer/{config,models,ncbi,analysis,plotting,export}` must never import `tkinter`.
- **Never block the Tk main thread**: Run NCBI requests, MSA subprocesses, and pairwise alignments through `gui/worker.py` (`BackgroundRunner` + `queue.Queue` + `root.after`).
- **Matplotlib**: Build standalone `matplotlib.figure.Figure` instances in `plotting/charts.py`; do not use global `matplotlib.pyplot`.
- **HiDPI Tkinter styling**: Use the scaled Tk named fonts (`TkDefaultFont`, `TkTextFont`, `TkHeadingFont`, `TkFixedFont`) configured in `gui/app.py`—avoid hardcoded small font sizes, and explicitly pass `font="TkTextFont"` to `ttk.Entry`, `ttk.Spinbox`, and `ttk.Combobox`.

## Domain & Testing Conventions

- **Secrets/Config**: Never commit `config.yaml`. Keep `config.yaml.example` in sync with `config.py`.
- **Ambiguous residues**: Always handle non-standard amino acids (`X, B, Z, J, U, O, *`) safely in `ProtParam` (flag `is_approximate = True`) and `PairwiseAligner` (map non-BLOSUM62 chars to `X`).
- **Offline tests**: Unit tests in `tests/` must mock `Bio.Entrez` and external MSA subprocesses (`clustalo`/`muscle`) so `uv run pytest` runs fast and offline.
- **Verify changes**: Run `uv run pytest` after edits.
