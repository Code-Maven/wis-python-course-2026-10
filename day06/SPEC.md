# NCBI Protein Search & Analyzer — Application Specification

## 1. Purpose

**NCBI Protein Search & Analyzer** (`ncbi-protein-analyzer`) is a Python desktop GUI application that enables users to search the NCBI Protein database, preview and select matching protein records, download their FASTA sequences, analyze physicochemical properties and sequence similarity, visualize results with embedded charts, and export data and figures.

---

## 2. Technology Stack & Environment

- **Language**: Python `>= 3.10`
- **Package & Environment Management**: `uv` (`pyproject.toml` + `uv.lock`)
- **Desktop GUI**: Tkinter + `ttk` (Python standard library)
- **Bioinformatics & NCBI Access**: Biopython (`Bio.Entrez`, `Bio.SeqIO`, `Bio.SeqUtils.ProtParam`, `Bio.Align`)
- **Visualization**: Matplotlib (`matplotlib.figure.Figure` embedded via `FigureCanvasTkAgg` + `NavigationToolbar2Tk`)
- **Configuration**: YAML (`PyYAML`)
- **Testing**: `pytest`
- **Optional External MSA Binaries**: Clustal Omega (`clustalo`) or MUSCLE (`muscle`, v3 or v5)

---

## 3. Configuration (`config.yaml`)

The application reads its configuration from `config.yaml` in the working directory by default, or from a custom path passed via `--config PATH`. A template is provided in `config.yaml.example`, and `config.yaml` is gitignored to avoid committing user credentials.

### Schema

```yaml
ncbi:
  email: "your.email@example.com"        # Required by NCBI E-utilities
  api_key: ""                            # Optional; raises rate limit from 3 to 10 req/sec
  tool: "ncbi_protein_analyzer"          # Tool identifier sent to NCBI
  max_results_default: 50                # Default max search results (>= 1)
  fetch_batch_size: 50                   # Batch size for esummary/efetch requests (>= 1)
  timeout_seconds: 30                    # Network timeout per request (>= 1)

analysis:
  max_pairwise_sequences: 50             # Guardrail cap for O(N^2) pairwise fallback (>= 2)
  msa:
    preferred_tool: "auto"               # "auto" | "clustalo" | "muscle" | "pairwise"
    clustalo_path: ""                    # Optional explicit executable path
    muscle_path: ""                      # Optional explicit executable path
    timeout_seconds: 120                 # Subprocess timeout for external MSA (>= 1)
```

### Validation Rules
- Startup fails fast with a clear error message (printed to `stderr` and shown in a Tk error dialog when a display is available) if:
  - `config.yaml` is missing or invalid YAML.
  - `ncbi.email` is empty, missing `@`, or still set to a placeholder address (`your.email@example.com`, etc.).
  - Numeric bounds or `analysis.msa.preferred_tool` values are invalid.

---

## 4. Functional Specification

### 4.1 Two-Step Search & Selection Workflow

#### Step 1 — Search & Metadata Preview
- Users enter a required **Query** term and optional filters:
  - **Organism** (appends `AND "<organism>"[Organism]`)
  - **RefSeq only** checkbox (appends `AND srcdb_refseq[PROP]`)
  - **Length filter (aa)** min/max bounds (appends `AND <min>:<max>[SLEN]`)
  - **Max results** (`1`–`500`, defaulting to `ncbi.max_results_default`)
- The client queries `Bio.Entrez.esearch(db="protein", ...)` to obtain total hit count and UIDs, then fetches metadata in batches via `Bio.Entrez.esummary(db="protein", retmode="json")`.
- Results populate the **1. Search Results** tab with:
  - Client-side text filtering across accession, title, and organism without losing selection state.
  - Sortable columns: selection status (`[x]` / `[ ]`), Accession, Title, Organism, Length (aa), Updated date.
  - Bulk selection helpers (`Select All Visible`, `Select Top 20`, `Clear Selection`) and row click / `Space` toggling.
  - Up to the first 20 search hits are pre-selected for convenience.

#### Step 2 — FASTA Download & Analysis
- Clicking **Download FASTA & Analyze** fetches FASTA records for selected proteins in batches via `Bio.Entrez.efetch(db="protein", rettype="fasta", retmode="text")` and parses them with `Bio.SeqIO`.
- Organism metadata is preserved from the `esummary` record (or parsed from trailing `[Organism]` brackets in the FASTA header as fallback).
- Once downloaded, the app runs physicochemical statistics and sequence comparison automatically and populates Tabs 2, 3, and 4.

### 4.2 NCBI Client Reliability & Rate Limiting
- Sets `Entrez.email`, `Entrez.tool`, and optional `Entrez.api_key`.
- Enforces minimum inter-request spacing (`0.35s` without API key, `0.11s` with API key) to stay within NCBI E-utilities rate limits.
- Retries transient errors (HTTP `429, 500, 502, 503, 504`, network errors, timeouts) up to 3 attempts with exponential backoff.
- Checks a `threading.Event` cancellation token between batches and retries so long operations can be cancelled cleanly.

### 4.3 Physicochemical Sequence Statistics
Computed per protein and aggregated across the dataset using `Bio.SeqUtils.ProtParam.ProteinAnalysis`:
- **Per-protein metrics**:
  - Total sequence length (aa)
  - Standard AA count (20 canonical amino acids: `A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y`)
  - Ambiguous / non-standard residue count and percentage (`X, B, Z, J, U, O, *`, etc.)
  - `is_approximate` flag (`True` whenever `ambiguous_aa_count > 0`)
  - Molecular weight (Da, displayed in kDa in GUI and Da in CSV)
  - Isoelectric point (`pI`)
  - GRAVY (Grand Average of Hydropathy)
  - Percentage composition across all 20 standard amino acids (summing to `100%` over standard residues)
- **Ambiguous residue handling**:
  - Non-standard residues are stripped before invoking `ProteinAnalysis` so Biopython never crashes.
  - Affected proteins are flagged `is_approximate = True`, highlighted in the UI table with a `~` prefix and amber row highlight, and annotated with a warning string in both UI and CSV.
  - If a sequence contains zero standard amino acids (e.g., all `X`), numeric metrics (`MW`, `pI`, `GRAVY`) are set to `None` (`N/A` in UI, blank in CSV) with an explanatory warning.
- **Dataset summary**:
  - Protein count, min/median/mean/max length, mean MW/pI/GRAVY, approximate protein count, total ambiguous residues, exact duplicate count, and residue-weighted aggregate amino acid composition.

### 4.4 Sequence Comparison (Deduplication, External MSA, Pairwise Fallback)
1. **Exact Sequence Deduplication**:
   - Identical normalized (uppercase, whitespace-stripped) sequences are grouped into clusters before alignment.
   - Cluster representatives are aligned, and duplicate clusters are reported in the **4. Sequence Comparison** tab.
2. **External MSA Integration**:
   - Controlled by `analysis.msa.preferred_tool` (`auto`, `clustalo`, `muscle`, or `pairwise`) and re-runnable directly from Tab 4 without re-downloading FASTA.
   - In `auto` mode, checks `clustalo` first, then `muscle` (using explicit config path if given, otherwise `shutil.which`).
   - Supports both **Clustal Omega** (`clustalo -i ... -o ... --outfmt=fasta --force`) and **MUSCLE** (auto-detecting v5 `-align/-output` vs v3 `-in/-out`).
   - Writes temporary FASTA files using safe internal IDs (`seq_0`, `seq_1`, ...) and maps aligned sequences back to original accessions and metadata.
   - Computes an $N \times N$ percent-identity matrix over non-gap aligned columns and a majority-rule consensus sequence.
3. **Pure-Biopython Pairwise Fallback**:
   - Used automatically when no external MSA binary is found, when `pairwise` mode is selected, or when an external MSA subprocess fails/times out (with an explicit warning).
   - Uses `Bio.Align.PairwiseAligner` in global mode with the `BLOSUM62` substitution matrix and gap scores (`open = -10.0`, `extend = -0.5`).
   - Maps characters outside the `BLOSUM62` alphabet (e.g., `U, O, J, *`) to `X` before scoring.
   - Enforces `analysis.max_pairwise_sequences` (default `50`) to prevent excessive $O(N^2)$ runtime, appending a warning when truncation occurs.

### 4.5 Embedded Visualizations
Built as standalone `matplotlib.figure.Figure` instances (independent of `pyplot` global state):
1. **Length Distribution Histogram**: Distribution of sequence lengths with dashed mean and dotted median reference lines.
2. **Amino Acid Composition Bar Chart**: 20 standard amino acids color-coded by physicochemical group (Hydrophobic, Polar, Positive, Negative, Special), switchable between **Dataset Average** and any individual downloaded protein.
3. **Pairwise Sequence Identity Heatmap**: $N \times N$ heatmap (`0–100%`, `viridis` colormap, annotated cell percentages when $N \le 15$) labeled with the comparison method used.

### 4.6 Data & Plot Exports
- **CSV Statistics (`export_stats_csv`)**: UTF-8 CSV containing `accession`, `description`, `organism`, `length`, `standard_aa_count`, `ambiguous_aa_count`, `ambiguous_aa_pct`, `is_approximate`, `molecular_weight_da`, `isoelectric_point`, `gravy`, `warning`, and `aa_<X>_pct` for all 20 standard amino acids.
- **FASTA Sequences (`export_fasta`)**: Exports downloaded unaligned sequences (from Tab 2) and aligned sequences when an external MSA was run (from Tab 4).
- **PNG Figures (`export_figure_png`, `export_all_figures_png`)**: Exports either the active chart or all generated charts at `150 DPI`.

---

## 5. GUI Architecture, Responsiveness & HiDPI Scaling

### Concurrency Model
- All blocking work (NCBI `esearch`/`esummary`/`efetch`, dataset analysis, MSA subprocess execution, and pairwise alignment) runs on a daemon worker thread via `BackgroundRunner` (`gui/worker.py`).
- Worker-to-GUI communication uses a thread-safe `queue.Queue` polled from the Tk main thread via `root.after(80, ...)`.
- Action buttons disable while busy; a **Cancel** button sets the active `threading.Event`; a bottom status bar displays live progress text and a `ttk.Progressbar`.

### Tab Layout
1. **1. Search Results**: Filterable/sortable `esummary` table + selection controls + Download & Analyze trigger.
2. **2. Sequence Stats**: Dataset KPI bar + sortable ProtParam table + row detail bar + Export CSV / Export FASTA buttons.
3. **3. Charts**: Chart selector (Length Histogram / AA Composition / Identity Heatmap), AA source selector, embedded Matplotlib canvas with navigation toolbar, Save Chart PNG / Export All PNGs buttons.
4. **4. Sequence Comparison**: Detected MSA tool status, mode selector + Re-run Comparison button, deduplication & warnings summary, monospace scrollable Clustal-style alignment viewer (or ranked pairwise identity table in pairwise mode), and Export Aligned FASTA button.

### Responsive Window & HiDPI Scaling
- On startup, `_configure_window_and_theme()` inspects screen dimensions (`winfo_screenheight()`) and configures window size, named Tk fonts (`TkDefaultFont`, `TkTextFont`, `TkHeadingFont`, `TkFixedFont`, etc.), Treeview `rowheight`, column widths, and widget padding:
  - HiDPI (`height >= 1800px`): `2400×1550` window, `20pt` UI/input fonts, `18pt` monospace font, `50px` Treeview rows (`2.0×` column scale).
  - Medium (`height >= 1200px`): `1800×1150` window, `16pt` UI / `15pt` monospace (`1.5×` scale).
  - Standard (`height < 1200px`): `1440×900` window, `14pt` UI / `13pt` monospace (`1.25×` scale).
- All `ttk.Entry`, `ttk.Spinbox`, and `ttk.Combobox` widgets explicitly bind `font="TkTextFont"` for consistent rendering on X11.

---

## 6. Package Structure

```text
day06/
├── pyproject.toml
├── uv.lock
├── config.yaml.example
├── config.yaml                          # User-local (gitignored)
├── README.md
├── SPEC.md
├── src/
│   └── ncbi_protein_analyzer/
│       ├── __init__.py
│       ├── __main__.py                  # CLI entrypoint (`ncbi-protein-analyzer`)
│       ├── config.py                    # YAML config loading & validation
│       ├── models.py                    # Shared domain dataclasses
│       ├── ncbi/
│       │   └── client.py                # Bio.Entrez client (search, summary, FASTA fetch)
│       ├── analysis/
│       │   ├── stats.py                 # ProtParam stats + ambiguous AA handling
│       │   ├── deduplicate.py           # Exact sequence clustering
│       │   ├── msa.py                   # External clustalo / muscle integration
│       │   └── compare.py               # Comparison orchestrator + PairwiseAligner fallback
│       ├── plotting/
│       │   └── charts.py                # Matplotlib Figure builders
│       ├── export/
│       │   └── writers.py               # CSV, FASTA, and PNG exporters
│       └── gui/
│           ├── worker.py                # Thread + Queue bridge to Tk mainloop
│           ├── app.py                   # Main window, search bar, HiDPI theme, status bar
│           ├── search_tab.py            # Tab 1: Search preview & selection
│           ├── stats_tab.py             # Tab 2: Sequence stats & CSV/FASTA export
│           ├── charts_tab.py            # Tab 3: Embedded charts & PNG export
│           └── compare_tab.py           # Tab 4: MSA/pairwise viewer & aligned FASTA export
└── tests/
    ├── conftest.py
    ├── test_config.py
    ├── test_ncbi_client.py
    ├── test_stats.py
    ├── test_compare.py
    ├── test_plotting.py
    └── test_export.py
```

---

## 7. Entrypoints & Verification

- **Install / sync**: `uv sync`
- **Launch GUI**: `uv run ncbi-protein-analyzer [--config PATH]` (or `uv run python -m ncbi_protein_analyzer`)
- **Automated tests**: `uv run pytest`
  - Unit tests cover config validation, mocked `Bio.Entrez` search/summary/fetch/retry/cancellation, ProtParam stats with standard and ambiguous/all-ambiguous sequences, exact deduplication, mocked `clustalo` & `muscle` v5 subprocess execution, MSA-failure fallback to `PairwiseAligner`, pairwise caps, headless Matplotlib chart generation, and CSV/FASTA/PNG export integrity.
