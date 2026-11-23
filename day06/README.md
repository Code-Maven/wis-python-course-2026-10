# NCBI Protein Search & Analyzer

A modular Python desktop GUI application (Tkinter + Biopython + Matplotlib) for searching the NCBI Protein database, previewing and selecting matching records, downloading FASTA sequences, computing physicochemical properties, comparing sequences via external MSA or pairwise alignment, visualizing results, and exporting data.

## Features

- **Two-step NCBI E-utilities workflow** (`Bio.Entrez` + `Bio.SeqIO`):
  1. Search NCBI Protein (`esearch` + `esummary` JSON) with optional Organism, RefSeq-only, and sequence length filters, and preview results in a sortable/filterable table.
  2. Select desired proteins to download FASTA sequences in rate-limited batches and run analysis without freezing the UI.
- **Physicochemical statistics** (`Bio.SeqUtils.ProtParam.ProteinAnalysis`):
  - Length, molecular weight (Da/kDa), isoelectric point (pI), GRAVY hydropathy, and 20-amino-acid composition.
  - Safely handles ambiguous/non-standard residues (`X, B, Z, J, U, O, *`) by counting them, computing stats on the 20 standard AAs, and flagging rows as approximate (`is_approximate = True`).
- **Sequence comparison**:
  - Groups exact duplicate sequences before alignment.
  - Automatically detects and uses **Clustal Omega (`clustalo`)** or **MUSCLE (`muscle`, v3 or v5)** if installed on `PATH` (or configured in `config.yaml`).
  - Seamlessly falls back to pure-Biopython **`Bio.Align.PairwiseAligner`** (BLOSUM62) when no external MSA binary is installed.
- **Embedded visualizations & exports**:
  - Length distribution histogram, amino acid composition bar chart (dataset average or individual protein), and pairwise identity heatmap.
  - Export statistics to **CSV**, downloaded or aligned sequences to **FASTA**, and charts to **PNG**.

## Setup (`uv`)

1. Install dependencies:
   ```bash
   uv sync
   ```
2. Create your local `config.yaml` from the example template and set your real email address (required by NCBI E-utilities):
   ```bash
   cp config.yaml.example config.yaml
   # Edit config.yaml and set ncbi.email (and optionally ncbi.api_key)
   ```
3. Launch the desktop application:
   ```bash
   uv run ncbi-protein-analyzer
   ```
   Or specify a custom config path:
   ```bash
   uv run ncbi-protein-analyzer --config /path/to/config.yaml
   ```

## Optional External MSA Tools

If you want Multiple Sequence Alignment (with consensus & aligned FASTA export) instead of the built-in Biopython pairwise fallback, install either Clustal Omega or MUSCLE:

- **Ubuntu / Debian**:
  ```bash
  sudo apt install clustalo
  # or
  sudo apt install muscle
  ```
- **macOS (Homebrew)**:
  ```bash
  brew install clustal-omega
  # or
  brew install brewsci/bio/muscle
  ```

You can also point directly to an executable in `config.yaml` under `analysis.msa.clustalo_path` or `analysis.msa.muscle_path`.

## Running Tests

```bash
uv run pytest
```
