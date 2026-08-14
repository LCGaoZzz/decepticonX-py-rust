# Official MuSiC2 data validation (2026-08-15)

## Outcome

`decepticonX-py-rust` commit
[`86ddffaf95da323ded426ce3bc54e5cff7ab6fde`](https://github.com/LCGaoZzz/decepticonX-py-rust/commit/86ddffaf95da323ded426ce3bc54e5cff7ab6fde)
completed all 12 configured non-EPIC branches on the official MuSiC2
pancreatic-islet data. No branch was skipped, and every 100-sample by
six-cell-type result was finite, non-negative, and schema-complete.

This validates the external-data pipeline contract. It does **not** support a
claim that the DECEPTICON correlation consensus improves accuracy: the best
single branch was materially more accurate than either closed consensus mode.
The parameterized download, conversion, run, and analysis commands are in the
[validation harness](../../README.md).

## Frozen data and conversion

The upstream data revision was
[`0f8e5aa71d6fe76d346c913a823b9f1e6d6bb800`](https://github.com/xuranw/MuSiC/tree/0f8e5aa71d6fe76d346c913a823b9f1e6d6bb800/data):

- `EMTABsce_healthy.rds`: healthy single-cell reference;
- `bulk-eset.rds`: healthy and type-2-diabetes bulk mixtures;
- `true_proportion.RData`: known six-cell-type proportions.

The pinned objects contain 100 bulk samples (50 healthy and 50 T2D), rather
than the 200 stated by the current tutorial text. Alignment used stored sample
IDs, never positional order. The selected reference contains 995 cells from
six healthy donors and six cell types (acinar, alpha, beta, delta, ductal, and
gamma). The run used 10,000 shared genes.

The upstream single-cell `counts` assay contains 9,638 fractional entries
(0.038% of the dense matrix; maximum distance from an integer 0.5). The exact
assay was retained in h5ad `X`; the strict primary run used an explicit rounded
integer `layers["counts"]`. Exact-versus-rounded reference sensitivity was
negligible: relative Frobenius differences were `1.44e-5` for BayesPrism,
`1.20e-5` for Monocle3, and `2.58e-6` for the MuSiC2-style template. Source and
converted-input SHA-256 values are in [`audit.json`](audit.json).

No source data, converted matrices, fitted estimates, or `run.json` are
committed in this report.

## Runtime configuration

- Python 3.11.3 on Windows 10 (build 19045).
- Three references: BayesPrism-style, Monocle3-style, corrected MuSiC2-style.
- Four methods: CIBERSORT, CIBERSORT-ABS, DeconRNASeq, and MuSiC.
- Eight configured threads and default CIBERSORT quantile normalization.
- Default `r_literal` consensus; `corrected` was recomputed from the identical
  12 saved estimates to isolate consensus behavior.
- EPIC was omitted because no validated, authorized six-type mRNA-per-cell
  mapping was available; an arbitrary `{"default": 1}` was not used for an
  accuracy run.

Pinned optional components:

| Component | Commit | Engine in this run |
|---|---|---|
| `decepticon-fast` | `a98b66c7da5fad81ba6a1eaddb9c826daac31fad` | Rust |
| `python-cibersort-rs` | `cd957af06357d2636ba18af1135bf70b488dfbbb` | Rust |
| `deconrnaseq-py` | `a2ec1b9ebe341628e945a6d25c6068cad975f6b5` | NumPy |
| `music-py` | `a4464c84f87fbeb61f0fa7032625bad21add52dd` | NumPy/SciPy |

## Structural checks

- Completed branches: 12/12; skipped branches: 0.
- Saved primary equals the saved R-literal unclosed result exactly.
- Recomputed R-literal closed maximum absolute difference: `2.78e-16`.
- Closed row-sum maximum error: `5.55e-16`.
- Corrected mode required no fallback and selected the same pairs and endpoint
  weights as R-literal mode, so their closed results were identical here.

## Accuracy against known proportions

Every branch was closed over its native output columns before scoring. RMSE
below is pooled over all aligned sample-by-cell-type values; it is not the
macro-average of six independently calculated cell-type RMSE values.

| Result | All RMSE | Healthy RMSE | T2D RMSE | All Pearson |
|---|---:|---:|---:|---:|
| CIBERSORT x BayesPrism | **0.054668** | 0.020569 | **0.074526** | 0.947400 |
| MuSiC x BayesPrism | 0.065449 | 0.022236 | 0.089848 | 0.923116 |
| DeconRNASeq x BayesPrism | 0.067308 | 0.040764 | 0.086017 | 0.926306 |
| MuSiC x MuSiC2-style | 0.072774 | **0.001890** | 0.102901 | 0.903955 |
| R-literal closed consensus | 0.094504 | 0.058011 | 0.120403 | 0.836486 |
| Corrected closed consensus | 0.094504 | 0.058011 | 0.120403 | 0.836486 |
| R-literal primary (unclosed) | 0.108639 | 0.076046 | 0.133499 | 0.832979 |

Relative to the best single branch, closed-consensus RMSE was 1.73x higher
overall, 30.7x higher for healthy samples, and 1.62x higher for T2D samples.
Unclosed row sums ranged from 1.025 to 1.369, so `consensus_closed.tsv` is the
appropriate view for comparison with fraction truth.

The near-exact healthy `music__music2` result supports the corrected raw
cell-type-mean template. Its T2D degradation is consistent with applying a
healthy reference across conditions without the full MuSiC2 algorithm's
iterative removal of cell-type-specific DE genes. In this package, `music2`
names corrected signature construction only; it is not the complete
condition-aware MuSiC2 procedure described in the
[official tutorial](https://xuranw.github.io/MuSiC/articles/pages/MuSiC2.html)
and [paper](https://academic.oup.com/bib/article/23/6/bbac430/6751147).

The consensus selected mutually correlated strategies rather than the most
truth-accurate strategy. Agreement is therefore not an accuracy oracle, and
pair deduplication did not alter this non-degenerate case.

## Runtime

| Stage | Seconds |
|---|---:|
| Bulk input | 0.607 |
| Single-cell input | 0.340 |
| Three references | 0.166 |
| Twelve backend branches | 2242.237 |
| Consensus | 0.009 |
| Total | 2243.359 (37.39 min) |

The recorded backend timing is aggregate, so it cannot attribute the 37
minutes to an individual method. Timing is one observed run, not a controlled
performance comparison.

## Machine-readable files

- [`metrics.tsv`](metrics.tsv): every branch and consensus metric for all,
  healthy, and T2D samples.
- [`cell_type_metrics.tsv`](cell_type_metrics.tsv): per-cell-type metrics.
- [`audit.json`](audit.json): input hashes, engines, timings, consensus
  reconstruction tolerances, and count-rounding sensitivity.
- [`provenance.json`](provenance.json): frozen source, package, backend, and
  environment revisions.
