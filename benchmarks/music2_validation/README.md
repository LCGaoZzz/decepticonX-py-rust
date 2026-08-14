# Official MuSiC2 data validation

This directory contains a reproducible harness and small, commit-bound report
snapshots for the official MuSiC2 pancreatic-islet data. It deliberately
contains no source RDS/RData files, converted h5ad or expression matrices,
fitted branch matrices, or `run.json` files.

The three source objects are downloaded from the
[MuSiC repository at commit `0f8e5aa`](https://github.com/xuranw/MuSiC/tree/0f8e5aa71d6fe76d346c913a823b9f1e6d6bb800/data).
The downloader uses commit-pinned raw URLs and checks the three SHA-256 values
recorded in [`report/2026-08-15/audit.json`](report/2026-08-15/audit.json)
before installing any file. Local downloads and all generated work belong
under the ignored `data/` and `work/` directories.

## Requirements

- Python 3.11 or 3.12 with this package and the four public, commit-pinned
  backends installed;
- a Rust toolchain so `decepticon-fast` and CIBERSORT use the native engines
  recorded by the frozen run;
- R with Bioconductor `Biobase` and `SingleCellExperiment`, plus `Matrix`;
- enough disk space for the downloaded objects, the 995-cell by 25,453-gene
  AnnData reference, and 12 fitted branch matrices;
- substantial runtime: the recorded eight-thread run took 37.39 minutes.

From the repository root, install the package and the same public backend
revisions used by the frozen report:

```text
python -m pip install -e .
python -m pip install "git+https://github.com/LCGaoZzz/decepticon-fast.git@a98b66c7da5fad81ba6a1eaddb9c826daac31fad"
python -m pip install "git+https://github.com/LCGaoZzz/python-cibersort-rs.git@cd957af06357d2636ba18af1135bf70b488dfbbb"
python -m pip install "git+https://github.com/LCGaoZzz/deconrnaseq-py.git@a2ec1b9ebe341628e945a6d25c6068cad975f6b5"
python -m pip install "git+https://github.com/LCGaoZzz/music-py.git@a4464c84f87fbeb61f0fa7032625bad21add52dd"
```

The reference run intentionally omits EPIC because no validated, authorized
mRNA-per-cell mapping is available for these six cell types.

## Reproduce the validation

All commands below are run from the repository root and use only relative
paths.

### 1. Download and verify the official objects

```text
python benchmarks/music2_validation/download_music2_data.py --output-dir benchmarks/music2_validation/data
```

Re-running the command verifies existing files instead of downloading them
again. A missing or mismatched hash is a hard failure. Use `--help` to inspect
the downloader's options; `--audit` can point to another frozen report only
when intentionally validating a different snapshot.

### 2. Export the R objects

```text
Rscript benchmarks/music2_validation/export_music2_official.R --data-dir benchmarks/music2_validation/data --output-dir benchmarks/music2_validation/work/prepared
```

This selects the six official cell types, aligns truth by sample ID, and writes
the bulk matrix, truth, sample metadata, sparse reference matrix, and reference
metadata. It never aligns samples by position.

### 3. Construct and validate `reference.h5ad`

```text
python benchmarks/music2_validation/build_music2_h5ad.py --input-dir benchmarks/music2_validation/work/prepared --output benchmarks/music2_validation/work/prepared/reference.h5ad
python -m decepticonx validate benchmarks/music2_validation/work/prepared/reference.h5ad benchmarks/music2_validation/work/prepared/bulk.tsv --cell-type-key cell_type --sample-key donor --counts-layer auto
```

The official `counts` assay contains 9,638 fractional entries. The converter
preserves the exact values in `X` and writes an explicit rounded integer
`layers["counts"]`; the strict run selects that layer. The analysis step below
quantifies the effect of this documented intervention.

### 4. Run the 12 public branches

```text
python -m decepticonx run benchmarks/music2_validation/work/prepared/reference.h5ad benchmarks/music2_validation/work/prepared/bulk.tsv --cell-type-key cell_type --sample-key donor --counts-layer auto --method cibersort cibersort_abs deconrnaseq music --reference bayesprism monocle3 music2 --cibersort-qn --cibersort-seed 0 --cibersort-engine rust --deconrnaseq-backend numpy --music-backend numpy --consensus-mode r_literal --consensus-pairs 2 --threads 8 --output-dir benchmarks/music2_validation/work/music2_r_literal
```

The explicit method, reference, consensus, and thread flags freeze the run
contract rather than relying on changing defaults. The output directory must
be new or empty.

### 5. Reconstruct corrected consensus and score against truth

```text
python benchmarks/music2_validation/analyze_music2_results.py --source-data-dir benchmarks/music2_validation/data --prepared-dir benchmarks/music2_validation/work/prepared --run-dir benchmarks/music2_validation/work/music2_r_literal --output-dir benchmarks/music2_validation/work/analysis
```

The analyzer:

- reconstructs R-literal consensus and verifies it against saved output;
- requires the exact 4-method x 3-reference branch set and a two-pair
  R-literal source run;
- derives corrected consensus from the identical 12 branch estimates;
- closes every native branch before comparing it with fraction truth;
- writes `metrics.tsv`, `cell_type_metrics.tsv`, and a compact `audit.json`;
- hashes all three source objects and the three prepared inputs;
- measures exact-versus-rounded reference sensitivity.

The corrected consensus files are written under
`work/analysis/corrected_consensus/` by default. Override that location with
`--corrected-output-dir` if necessary.

## Versioned snapshot

- [2026-08-15 report](report/2026-08-15/README.md): 100 bulk samples with
  known six-cell-type proportions, evaluated against 12 non-EPIC branches and
  both consensus semantics.

Only the compact report files should be versioned. Do not add anything from
`data/` or `work/`; in particular, do not commit upstream objects, h5ad files,
bulk matrices, fitted estimates, or `run.json`.

These snapshots are scientific validation records, not CI benchmarks. The
network download, R conversion, optional backends, and full deconvolution run
are intentionally not executed in CI.
