# DECEPTICONx original-versus-accelerated comparison

This directory contains the reproducible harness used to compare the common
three-reference/five-method subset of:

- `Hao-Zou-lab/DECEPTICONx` at
  `f93ce0f01a674c31bed5c217bd5072d7d845311c`;
- `Hao-Zou-lab/DECEPTICON` at
  `53342ce0077001b17d100cb17fc711511baf0d90`; and
- this Python/Rust package plus the revisions in `backend-lock.json`.

The high-level R wrapper is not used. Its custom path ignores `light=TRUE`,
runs methods outside the common five, and its consensus assumes the original
ten-method file layout. `run_original.R` calls the original three builders and
five exported custom methods without changing their algorithm code.

## Data limitation

The supplied `simulation_data.txt` is truncated in its final `AFP` row. Base R
reads it as 363 genes by 100 samples, pads 11 missing values, and emits a
warning; the Python loader rejects it. Both sides therefore use the identical
mechanically derived 362-complete-row prefix. The single-cell RDA contains
normalised non-integer expression, not raw counts. Results from this input are
compatibility smoke measurements, not formal speed or accuracy claims.

The exact audited hashes are:

```text
source bulk  77ca6cf40384dbee4fe1114103b42032e5232bfc3829c22c4449090a576ab2c3
362-row file 883bf9ae63507dd16094935b82b1803ac6ed03e28da02b96da9fde5e70d24a43
```

Recreate the prefix without imputation or numeric parsing:

```bash
python make_complete_prefix.py simulation_data.txt simulation_data_complete_rows_only.tsv
```

## CIBERSORT helper

Current DECEPTICONx `main` deleted `CIBERSORT.R`. The official repository's
historical `DECEPTICONx_1.0.0.tar.gz` tag still contains CIBERSORT v1.04. The
harness accepts an external helper directory but never copies that script into
this repository or any committed result. Users remain responsible for its
terms of use.

## Runs

The R harness writes the original references, all available branch outputs,
stage timings, warnings, session information, and a manifest:

```bash
R_LIBS_USER=/path/to/isolated-r-library \
Rscript run_original.R \
  --rda /data/single_cell_schelker.rda \
  --bulk /data/simulation_data_complete_rows_only.tsv \
  --decepticonx-r /src/DECEPTICONx/R/DECEPTICONx.R \
  --decepticon-methods-r /src/DECEPTICON/R/DECEPTICON_custom_methods.R \
  --helper-dir /licensed/or/historical/helper/checkout \
  --output /results/r-original
```

The frozen-reference solver run feeds those exact R signatures to the
accelerated implementations:

```bash
python run_accelerated_from_refs.py \
  --bulk /data/simulation_data_complete_rows_only.tsv \
  --references /results/r-original/custom_signature_matrix \
  --output /results/accelerated-r-refs-fidelity \
  --epic-mrna-json /secure/local/epic-mrna.json \
  --epic-backend python \
  --epic-solver nm \
  --threads 8
```

`nm` is the harness default because it forces EPIC's
`constrOptim`/Nelder-Mead-compatible port and therefore gives the closest
optimizer-fidelity comparison to the original R optimizer. In the pinned EPIC
implementation this path is executed by the Python reference kernel even when
`--epic-backend rust` is requested; the Rust batch kernel applies to `nmf` and
to the unresolved subset after `auto`'s Python QP probe. For a
production-throughput run, use `--epic-solver auto` and keep that result
separate from the optimizer-fidelity run. Both requested values are recorded
in `run.json`, and each successful EPIC branch records its requested solver in
diagnostics.

The original R custom method lets `EPIC::EPIC()` use the package's
`mRNA_cell_default`. An optimizer-fidelity comparison therefore also requires
an authorized local JSON object containing the same *effective* mapping for
these custom cell labels (and `otherCells`, or an explicit `default`). Pass
either that file path or an inline object with `--epic-mrna-json`. The harness
does not contain or redistribute EPIC defaults. Omitting the option
deliberately uses
`{"default": 1.0}` for a compatibility smoke run; that run is not an EPIC
cell-fraction parity claim. `run.json` stores only the mapping SHA-256, entry
count, and a redaction marker; inline values are also redacted from its command
record.

The R and accelerated harnesses hash the full supplied mapping using the same
sorted-key, 17-significant-digit format. `compare_results.py` requires the R
and accelerated hashes and entry counts to match before it accepts the
optimizer-fidelity comparison; mapping values never enter the report.

The frozen-reference run is the legacy-reference/optimizer-fidelity track: it
consumes the three exact, hash-verified R-generated files, including the
original `MuSiC2_base` behavior (column tokens are validated one-to-one and
relabelled to the Monocle3 spelling before dispatch). The full
accelerated product intentionally builds a corrected MuSiC2-style
mean-expression reference. Thus its reference and downstream differences are
product-behavior results, not kernel regressions. A reproducible full
accelerated E2E configuration for this compatibility dataset is:

```json
{
  "allow_normalized_x": true,
  "exclude_cell_types": ["Unknown"],
  "methods": [
    "cibersort", "cibersort_abs", "epic", "deconrnaseq", "music"
  ],
  "references": ["bayesprism", "monocle3", "music2"],
  "consensus_mode": "r_literal",
  "consensus_pairs": 2,
  "epic_mrna_cell": {"default": 1.0},
  "cibersort_seed": 20260814,
  "threads": 8,
  "cibersort_engine": "rust",
  "epic_backend": "rust",
  "epic_solver": "auto",
  "deconrnaseq_backend": "rust",
  "music_backend": "rust"
}
```

```bash
/usr/bin/time -v decepticonx run \
  /data/single_cell_schelker.h5ad \
  /data/simulation_data_complete_rows_only.tsv \
  --config /data/accelerated-e2e.json \
  --output-dir /results/accelerated-e2e
```

Change only `threads` (for example, 1 and 8) when measuring scaling. Keep the
seed, data, mappings, engine selections, and an empty output directory fixed.

Both accelerated tracks use the production `r_literal` selector/positional
weighting kernel with the deterministic C-locale R filename-equivalent branch
order. Branch aggregation uses the documented all-native-column composition
space, including EPIC `otherCells` in its branch denominator; the broken
original 15-strategy `decostand` block is not claimed as a
reproducible R output. Accuracy metrics use the explicitly closed consensus
view, while the mode-primary and unclosed views remain stored separately.

`compare_results.py` refuses a frozen comparison unless the R manifest and
frozen `run.json` agree on bulk SHA-256 and seed, all three recorded reference
hashes equal the actual R files, and EPIC requested the solver named by
`--expected-epic-solver` (`nm` by default). It also rebuilds the frozen
consensus from the 15 branches and checks the primary, unclosed, and closed files before
reporting direct original-versus-frozen consensus metrics.

Use external `/usr/bin/time -v` around both commands for wall time and peak
RSS. Pass the full accelerated E2E wall/RSS separately from the frozen solver
wall/RSS to `summarize_timings.py`; its options use the
`--accelerated-e2e-*` and `--accelerated-frozen-*` prefixes, and the script
never compares the full R wall time with a frozen-reference solver-only wall
time. The historical R
CIBERSORT helper internally forks up to three workers per sample, so the R run
is stock-source-default resource use, not a one-core match to either
`--threads 1` or `--threads 8`. `/usr/bin/time -v` maximum RSS is a
process-level peak whose child-process accounting and peak timing can differ
between R forks and Rust threads; treat it as an operational observation, not
an exact aggregate-memory equivalence. Internal totals intentionally retain
output-writing time as a separate stage.
