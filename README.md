# decepticonX-py-rust

An R-free, h5ad-first orchestration layer for a focused
[DECEPTICONx](https://github.com/Hao-Zou-lab/DECEPTICONx) workflow.

It keeps all three reference-construction strategies (BayesPrism,
Monocle3, and MuSiC2-style), runs five independently installed accelerated
deconvolution methods (CIBERSORT, CIBERSORT-ABS, EPIC, DeconRNASeq, and
MuSiC), and produces a DECEPTICON-style cross-reference consensus. The core
package is an adapter: it does not vendor algorithm implementations, native
binaries, EPIC reference data, or licensed mRNA-per-cell values.

> **Status: working alpha.** The end-to-end 3 x 5 branch path has been run in
> the Omicos Python 3.11 environment with every optional native component
> available. The supplied `simulation_data_human.zip` is not a valid accuracy
> oracle: its bulk matrix is truncated and its single-cell assay is normalized,
> not raw counts. See [the data audit](docs/simulation-data-audit.md).

中文要点：默认单细胞输入是标准 `h5ad`（cells x genes），默认严格要求原始
counts；三个 ref 构造算法都保留，五个解卷积方法都可用。EPIC 因授权数据
限制是显式启用项，必须由用户提供合法的 `mRNA_cell` 映射。测试 ZIP 目前
只能做接口兼容性冒烟测试，不能据此宣称解卷积精度。

## What is implemented

| Stage | Implementations |
|---|---|
| Reference construction | BayesPrism-style, Monocle3-style, corrected MuSiC2-style |
| Deconvolution | CIBERSORT relative, CIBERSORT `sig.score` absolute, EPIC, DeconRNASeq, MuSiC |
| Integration | Up to 15 method/reference branches and cross-reference correlation consensus |
| Inputs | Backed or in-memory h5ad plus strict genes-by-samples CSV/TSV/TXT bulk matrix |
| Outputs | Every signature, every branch estimate, normalized/raw consensus, timings, provenance, warnings, engines, and fit diagnostics |

CIBERSORT and CIBERSORT-ABS share one fit batch per reference. MuSiC follows
DECEPTICON's custom-signature behavior: each signature becomes five identical
pseudo-donors. The MuSiC2-style reference fixes two defects in the existing
fast port (one scalar library-size multiplier and incorrect gene order) and
does not execute the defective path merely for comparison.

## Requirements

- Python 3.11 or 3.12.
- A Rust toolchain when building native backends from source.
- Enough memory for the source single-cell object; h5ad reads are backed and
  compacted to bulk-overlapping genes before reference construction.
- Human gene symbols. Ensembl IDs require a symbol column selected with
  `gene_symbol_key`.

The exact revisions used for integration testing are in
[`backend-lock.json`](backend-lock.json).

## Install the core package

```bash
python -m venv .venv
source .venv/bin/activate                 # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install .
```

For development:

```bash
python -m pip install -e ".[test]"
python -m pytest
```

The core installation deliberately does not pull the optional algorithm
repositories. Install only software whose license and access terms apply to
your use.

## Install reference builders and deconvolution backends

### Reference builders and CIBERSORT

Both source installs build their Rust extensions when a Rust toolchain is
available:

```bash
python -m pip install \
  "git+https://github.com/LCGaoZzz/decepticon-fast.git@512c9c9d47b3e3a735f17e882aff07460560cb4c"
python -m pip install \
  "git+https://github.com/LCGaoZzz/python-cibersort-rs.git@cd957af06357d2636ba18af1135bf70b488dfbbb"
```

### DeconRNASeq

Install the Python implementation first:

```bash
python -m pip install \
  "git+https://github.com/LCGaoZzz/deconrnaseq-py.git@a2ec1b9ebe341628e945a6d25c6068cad975f6b5"
```

Its NumPy backend is complete. To require the optional Rust solver, build the
separately named `deconrnaseq-rust` wheel from an upstream checkout and install
it alongside `deconrnaseq-py`:

```bash
git clone https://github.com/LCGaoZzz/deconrnaseq-py.git
cd deconrnaseq-py
git checkout a2ec1b9ebe341628e945a6d25c6068cad975f6b5
python -m pip install . maturin
cd rust
python -m maturin build --release
python -m pip install target/wheels/deconrnaseq_rust-*.whl
```

That separate-wheel packaging fix is reviewed independently in
[deconrnaseq-py PR #3](https://github.com/LCGaoZzz/deconrnaseq-py/pull/3).

### MuSiC

```bash
git clone https://github.com/LCGaoZzz/music-py.git
cd music-py
git checkout a4464c84f87fbeb61f0fa7032625bad21add52dd
python -m pip install . maturin
cd rust
python -m maturin build --release
python -m pip install target/wheels/music_py_rust-*.whl
```

### EPIC

`epic-py-rust` is a private, academic/non-commercial licensed repository.
Use only an authorized checkout. Its reference profiles and default
`mRNA_cell` values are not bundled here or in that package.

```bash
git clone https://github.com/LCGaoZzz/epic-py-rust.git
cd epic-py-rust
git checkout 2e661cccb5c9d7c8327bd65253e3495a1d3c1991
python tools/build_native.py
python -m pip install .
```

Verify both Python APIs and native components:

```bash
decepticonx check-backends --require-all
decepticonx check-backends --require-native
```

`--require-all` checks every package plus the reference builders.
`--require-native` additionally fails if any Rust component is unavailable.

## Input contract

### Single-cell h5ad

AnnData orientation is cells x genes. With the default
`counts_layer="auto"`, the loader searches in this order:

1. `adata.layers["counts"]` when present and integer-like;
2. integer-like `adata.X`;
3. integer-like `adata.raw.X`.

Raw, finite, non-negative counts are required by default. A non-integer `X`
is accepted only with `allow_normalized_x=True` (CLI:
`--allow-normalized-x`), and the run is explicitly marked as a compatibility
run. Cell types default to `adata.obs["cell_type"]`; donor/sample metadata is
optional.

### Bulk matrix

Bulk input is genes x samples. The first column contains unique gene symbols:

```text
gene    sample_1    sample_2
GAPDH   120.0       98.0
CD3D    15.0        43.0
```

R-style tables whose header omits the first gene-column label are accepted.
Ragged rows, missing/non-finite values, negative expression, or duplicate
genes/samples are rejected rather than silently repaired.

See [the complete data contract](docs/data-contract.md).

## Run from the CLI

Validate inputs without importing or running an algorithm backend:

```bash
decepticonx validate cells.h5ad bulk.tsv \
  --cell-type-key cell_type \
  --sample-key donor
```

The safe default runs BayesPrism/Monocle3/MuSiC2 references with CIBERSORT,
CIBERSORT-ABS, DeconRNASeq, and MuSiC. EPIC is available but omitted from the
default method list because it cannot run lawfully or reproducibly without an
explicit mRNA-per-cell mapping.

```bash
decepticonx run cells.h5ad bulk.tsv \
  --output-dir results \
  --threads 8
```

To run all five methods, create `config.json` with an authorized EPIC mapping:

```json
{
  "methods": [
    "cibersort",
    "cibersort_abs",
    "epic",
    "deconrnaseq",
    "music"
  ],
  "epic_mrna_cell": {
    "B cell": 1.5,
    "T cell": 1.2,
    "otherCells": 1.0
  },
  "cibersort_engine": "rust",
  "epic_backend": "rust",
  "epic_solver": "auto",
  "deconrnaseq_backend": "rust",
  "music_backend": "rust",
  "threads": 8
}
```

```bash
decepticonx run cells.h5ad bulk.tsv --config config.json -o results-all-five
```

Every custom EPIC cell type plus `otherCells` needs a positive value, unless
the mapping contains an explicit `"default"`. `{"default": 1.0}` is accepted
as a deliberate unscaled approximation for a custom reference; it is **not**
equivalent to EPIC's licensed defaults or to the original R cell fractions.

`epic_solver` defaults to `"auto"`, preserving the accelerated package's
existing behavior. Use `"nm"` (or CLI `--epic-solver nm`) to force the
original R EPIC `constrOptim`/Nelder-Mead path for strict optimizer parity.
`"nmf"` uses the vectorized Nelder-Mead objective, while `"qp"` selects the
exact active-set solution and can differ from the optimizer result produced by
R. The solver choice is recorded in both run provenance and EPIC branch
diagnostics.

EPIC's requested backend and its execution class are recorded separately. The
`nm` and `qp` solvers always execute in Python, including when
`epic_backend="rust"` was requested. `nmf` resolves to Rust only when the
effective backend permits it and the native kernel is available. `auto` first
probes every sample with Python QP and sends only unresolved samples to Rust;
when native dispatch is available it is therefore recorded as `hybrid`, even
if the unresolved subset happens to be empty. Otherwise it resolves to Python.
Branch diagnostics retain
`backend_requested` and add `backend_effective`, `backend_resolved`,
`native_eligible`, and `backend_resolution_reason`. The branch-level `engines`
entry contains `rust`, `python`, or `hybrid`.

Use a new or empty output directory. The writer refuses a non-empty directory
so files from an earlier branch set cannot be mistaken for current results.

## Python API

```python
from decepticonx import PipelineConfig, run_decepticonx

config = PipelineConfig(
    cell_type_key="cell_type",
    sample_key="donor",
    threads=8,
    cibersort_engine="rust",
    deconrnaseq_backend="rust",
    music_backend="rust",
)

result = run_decepticonx(
    "cells.h5ad",
    "bulk.tsv",
    config=config,
    output_dir="results",
)

print(result.consensus)          # samples x cell types, rows sum to one
print(result.signatures.keys())
print(result.estimates.keys())
print(result.diagnostics["backends"]["engines"])
```

`run_decepticonx` also accepts an in-memory `AnnData` and a genes x samples
`pandas.DataFrame`.

## Output tree

```text
results/
├── signatures/
│   ├── bayesprism.tsv
│   ├── monocle3.tsv
│   └── music2.tsv
├── estimates/
│   ├── cibersort__bayesprism.tsv
│   └── ...
├── consensus.tsv
├── consensus_raw.tsv
└── run.json
```

`run.json` records input scale/range, selected count layer, exclusions,
component versions, requested and resolved engines, branch warnings, EPIC
`fit_gof`, timings, consensus pair selection, and skipped branches. EPIC
`mRNA_cell` values and inline command values are redacted; a SHA-256
fingerprint is retained so authorized configurations can be compared without
copying their contents into results.

## Convert the supplied ExpressionSet to h5ad

The one-time converter needs R packages `Biobase`, `Matrix`, and `reticulate`,
plus Python packages `anndata`, `numpy`, `pandas`, and `scipy`:

```bash
RETICULATE_PYTHON=/path/to/python \
Rscript scripts/convert_expressionset_to_h5ad.R \
  single_cell_schelker.rda single_cell_schelker.h5ad \
  --allow-normalized
```

The override is required because that RDA contains normalized non-integer
expression. If reticulate was built against the NumPy 1.x ABI, use a converter
environment with `numpy<2`; the script imports Python before loading the large
RDA so an ABI mismatch fails early.

The archive's original bulk file must fail strict validation at its truncated
row. Removing that partial row creates only a 362-gene orchestration smoke
input, not a scientific benchmark. The label aliases used to inspect the
complete truth table are in
[`examples/schelker_label_map.json`](examples/schelker_label_map.json):

```bash
python scripts/evaluate_against_truth.py \
  results/consensus.tsv raw_data/cell_fractions.txt \
  --label-map examples/schelker_label_map.json
```

Do not interpret those metrics until the complete bulk matrix and raw counts
are available. Mapping `Dendritic cells` to `Plasmacytoid dendritic cell` is
dataset-specific, not a universal biological synonym.

## Omicos environment

The integration run was executed from the Omicos Linux/WSL Python 3.11
environment. A generic local setup is:

```bash
OMICOS_ENV_DIR=/path/to/omicos-env
source "$OMICOS_ENV_DIR/.venv/bin/activate"
python -m pip install -e /path/to/decepticonX-py-rust
decepticonx check-backends
```

No change to `omicos-core` is required; this package runs as a normal Python
workflow inside the existing environment.

## Validation performed

- Unit/contract tests cover strict input handling, sparse/backed h5ad
  compaction, label normalization, three reference paths, five adapter paths,
  EPIC authorization gates and fit diagnostics, consensus masking/fallbacks,
  CLI configuration, provenance redaction, and stale-output prevention.
- A synthetic count dataset (5,000 genes, 1,000 cells, 8 bulk samples, 5 cell
  types) completed all 15 branches. One observed integration run took about
  117 seconds: reference construction about 0.14 seconds and optional
  backends about 116.9 seconds. This is a smoke timing, not a controlled
  cross-implementation benchmark.
- The mechanically cleaned 362-gene portion of the supplied human archive
  completed all 15 branches (100 samples x 11 types) in both the pinned
  original-source selective harness and this package. The reproducible
  [comparison report](benchmarks/original_comparison/report/2026-08-14/README.md)
  records 166.31 seconds for R, 21.85 seconds for accelerated one-thread, and
  10.71 seconds for accelerated eight-thread execution, along with reference,
  branch, consensus, truth, and memory metrics. Its input defects limit those
  numbers to an operational compatibility smoke test; they do not support a
  general speed or scientific accuracy claim.
- CI builds both sdist and wheel, installs the wheel, runs `pip check`, checks
  the CLI and byte-compilation, then runs the core test suite on Python 3.11
  and 3.12. Optional licensed/native integrations remain separate from public
  CI.

## Design notes and limitations

- Human symbols only in the first release.
- EPIC is opt-in and academic/non-commercial under its upstream terms.
- CIBERSORT's current port retains the v1.04 low-range inverse-log heuristic;
  audit the recorded bulk range and choose quantile normalization deliberately.
- `auto` backend selection can use a Python/NumPy fallback. Pass the explicit
  `rust` backend options to require acceleration; the resolved state is saved
  in `run.json`.
- Consensus normalizes heterogeneous branch scales before correlations and
  masks self-pairs and same-reference pairs. Undefined correlations use a
  deterministic median fallback.
- This is research software, not a clinical diagnostic system.

See [architecture and compatibility policy](docs/architecture.md) for the
pipeline details.

## Licensing and citation

This orchestration repository is MIT-licensed. The optional runtime components
are not all MIT-compatible and are not redistributed here: CIBERSORT/MuSiC
ports are GPL-3.0-or-later, DeconRNASeq is GPL-2.0-only, and EPIC has restricted
academic/non-commercial terms. The reference-builder package also contains
BayesPrism-derived annotation tables with separate notices.

Read [`NOTICE.md`](NOTICE.md) before distribution or combined-environment use.
Publications should cite DECEPTICON, DECEPTICONx, and each scientific method
actually used.
