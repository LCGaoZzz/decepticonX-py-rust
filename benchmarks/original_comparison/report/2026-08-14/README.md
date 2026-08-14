# DECEPTICONx comparison report (2026-08-14)

## Outcome

The selective original-source R harness and the accelerated package both
completed the same three-reference/five-method, 100-sample compatibility run.
The accelerated
full h5ad workflow took 21.85 seconds with one thread and 10.71 seconds with
eight threads, compared with 166.31 seconds for the selective original-source
three-reference/five-method harness.
Those are 7.61x and 15.53x wall-clock speedups, respectively.

This is a compatibility smoke benchmark, not a formal biological accuracy
benchmark. The supplied bulk file is truncated and the supplied single-cell
RDA contains normalized, non-integer expression. Both sides used only the 362
complete bulk rows; no value was imputed.

## Comparison tracks

Two tracks prevent reference construction from being confused with solver
fidelity:

1. **Frozen-reference optimizer fidelity:** the accelerated solvers consumed the exact
   three files emitted by R, verified by SHA-256. EPIC used the same effective
   locally authorized mRNA-per-cell mapping, verified by matching redacted
   cross-runtime hashes and entry counts. `solver=nm` resolved to the
   Python `constrOptim`/Nelder-Mead-compatible kernel.
2. **Product throughput:** the full package consumed the h5ad, built its three
   references, and ran all five methods with eight configured threads, native
   kernels enabled, and `EPIC solver=auto`. A second frozen-reference run used
   the same product engines.

On this input, the accelerated frozen-reference `auto` and `nm` runs were
numerically identical across all 15 branch matrices and both consensus
matrices (`max_abs = 0`). Thus the faster `auto` timings did not trade away
any result within the accelerated implementation for this dataset.

## Wall time and peak RSS

All rows are one independent process measured with `/usr/bin/time -v` on WSL2
ext4. The R run set BLAS/OpenMP thread variables to one, but historical
CIBERSORT v1.04 itself forks up to three workers. The frozen solver process is
shown for operational context and is not divided into the full R wall time.

| Run | Threads / EPIC engine | Wall (s) | Speedup vs full R | Max RSS (KB) |
|---|---:|---:|---:|---:|
| Original-source selective 3-ref/5-method harness | stock CIBERSORT forks; R EPIC | 166.31 | 1.00x | 14,252,724 |
| Accelerated, full h5ad | 1; Rust-enabled hybrid `auto` | 21.85 | 7.61x | 302,716 |
| Accelerated, full h5ad | 8; Rust-enabled hybrid `auto` | 10.71 | 15.53x | 353,656 |
| Accelerated, frozen R refs | 8; Python `nm` optimizer fidelity | 36.20 | not comparable E2E | 224,492 |
| Accelerated, frozen R refs | 8; Rust-enabled hybrid `auto` | 7.65 | not comparable E2E | 224,088 |

The full-process RSS ratio was 40.30x for R versus accelerated eight-thread
execution. Treat it as an operational peak, not exact aggregate-memory
equivalence: forked R children and threaded Rust can be accounted differently,
and their peaks need not occur simultaneously.

## Internal method timing

These are stage timers rather than independent processes. The R custom methods
write their results inside each stage, DeconRNASeq uses `fig=TRUE`, and R fits
CIBERSORT relative and ABS separately; accelerated writing is separate and the
CIBERSORT pair shares each fit. Ratios therefore describe the measured product
workflows, not identical wrapper overhead.

| Stage (three references) | Original R (s) | Accelerated (s) | R / accelerated |
|---|---:|---:|---:|
| Reference construction incl. accelerated SC preprocessing | 81.326 | 3.045 | 26.70x |
| CIBERSORT + CIBERSORT-ABS, frozen refs | 62.218 | 1.895 | 32.83x |
| EPIC, frozen refs, Rust-enabled hybrid `auto` | 8.992 | 1.745 | 5.15x |
| DeconRNASeq, frozen refs | 0.715 | 0.075 | 9.48x |
| MuSiC, frozen refs | 2.204 | 0.112 | 19.66x |
| Five-method backend total, Rust-enabled `auto` | 74.129 | 3.827 | 19.37x |
| Five-method backend total, `nm` optimizer fidelity | 74.129 | 32.081 | 2.31x |

The reference-construction ratio includes accelerated `single_cell_input`
(2.932 s) as well as aggregation (0.113 s), so gene-intersection, remainder
preparation, and loading work is not hidden. It also includes a deliberate
product optimization: the accelerated path retains only the 362
bulk-intersection genes plus a synthetic remainder during construction,
whereas original Monocle3 processes all 17,753 genes.

## Reference results

Metrics below compare expression values on common genes and identically
aligned cell types.

| Reference | R genes | Accelerated genes | Common | Flat Pearson | MAE | Max abs |
|---|---:|---:|---:|---:|---:|---:|
| BayesPrism | 334 | 334 | 334 | 1.000000 | 5.83e-18 | 1.04e-16 |
| Monocle3 | 17,753 | 362 | 362 | 0.999857 | 0.002686 | 0.164819 |
| MuSiC2 | 362 | 362 | 362 | 1.000000 | 2.04e-16 | 2.84e-14 |

The corrected MuSiC2 product implementation is not generally identical to the
legacy wrapper, but the two happened to produce machine-precision-equivalent
references on this dataset. BayesPrism was also equivalent. Monocle3 retained
all downstream-relevant genes but differed slightly after its reduced working
matrix and remainder-based normalization.

## Frozen-reference solver results

Each row summarizes row-normalized proportions across the three reference
branches. `Max abs` is the maximum over all samples, cell types, and reference
branches. `Mean MAE` is the mean of the three branch MAEs.

| Method | Max abs | Mean MAE | Mean RMSE | Median flat Pearson |
|---|---:|---:|---:|---:|
| CIBERSORT | 0.060333 | 0.000381 | 0.002210 | 0.999855 |
| CIBERSORT-ABS | 0.060333 | 0.000381 | 0.002210 | 0.999855 |
| DeconRNASeq | 1.64e-13 | 6.08e-15 | 1.47e-14 | 1.000000 |
| EPIC | 0.242588 | 0.010004 | 0.020015 | 0.995056 |
| MuSiC | 1.62e-11 | 8.61e-14 | 4.11e-13 | 1.000000 |

DeconRNASeq and MuSiC were equivalent to numerical precision. Original EPIC
warned that 89/100 Monocle3, 83/100 BayesPrism, and 86/100 MuSiC2 fits did not
fully converge on the 362-gene input. EPIC therefore dominates the largest
branch differences and should not be interpreted as a stable accuracy test on
this damaged input.

## Consensus results

| Comparison | Max abs | MAE | RMSE | Flat Pearson | Median per-cell Pearson |
|---|---:|---:|---:|---:|---:|
| R vs accelerated on frozen R refs | 0.085304 | 0.004911 | 0.011101 | 0.993922 | 0.998404 |
| R vs accelerated full h5ad E2E | 0.084375 | 0.005902 | 0.012069 | 0.992830 | 0.996301 |

The same neutral five-method consensus implementation was applied to both
sides because the original DECEPTICONx custom consensus is hard-coded for its
ten-method file layout and is undefined for this requested five-method subset.

## Compatibility truth metrics

All predictions were row-normalized and aligned to the supplied 100x11 truth
matrix through the committed label map.

| Run | Macro Pearson | Macro Spearman | MAE | Macro RMSE | Pooled RMSE |
|---|---:|---:|---:|---:|---:|
| Original R common-five consensus | 0.760373 | 0.756754 | 0.060883 | 0.079863 | 0.097071 |
| Accelerated on frozen R refs | 0.753929 | 0.751589 | 0.059595 | 0.078485 | 0.095808 |
| Accelerated full h5ad E2E | 0.751420 | 0.746552 | 0.059616 | 0.078617 | 0.095557 |

The original had slightly higher correlation, while accelerated E2E had 2.08%
lower MAE and 1.56% lower pooled RMSE. A 2,000-repeat sample bootstrap gave an
E2E pooled-RMSE delta (accelerated minus R) of -0.001513 with a 95% interval of
[-0.001983, -0.001080]. These values are reported only to detect gross
regressions; the truncated bulk and normalized single-cell input prevent a
scientific accuracy claim.

## Provenance and limitations

- Original DECEPTICONx commit: `f93ce0f01a674c31bed5c217bd5072d7d845311c`.
- Original DECEPTICON commit: `53342ce0077001b17d100cb17fc711511baf0d90`.
- Historical CIBERSORT v1.04 helper SHA-256:
  `8435ae008d5ec373fa8ebdf526c2e948a588f46676ab9d30770c29e143de3eca`
  (the helper is not redistributed here).
- R 4.3.3, Python 3.11.15, WSL2 Linux 6.18.33.2.
- Intel Core i9-13900KF, 24 cores / 32 logical processors; WSL memory 46 GiB.
- Source bulk SHA-256:
  `77ca6cf40384dbee4fe1114103b42032e5232bfc3829c22c4449090a576ab2c3`.
- Complete 362-row prefix SHA-256:
  `883bf9ae63507dd16094935b82b1803ac6ed03e28da02b96da9fde5e70d24a43`.
- RDA SHA-256:
  `c2b484e178cbc49966e26114f641fd3bcde6614e47398ece6f6853a1f83f76fb`.
- h5ad SHA-256:
  `555e2914cf8f95222ab7ccda1a735800fd4a8d54c83ec83ad6a02bd3185d7515`.
- Truth SHA-256:
  `09efe504699bf647215503549b0819784bcb4d3e67af8259a507ffb54828644c`.
- EPIC mapping values are not committed. The optimizer-fidelity run stored
  only a cross-runtime canonical mapping hash, entry count (12), and redaction
  marker; the R and accelerated hashes matched.
- Timings are single runs (`n=1`) and include process startup and native input
  formats (RDA for R, h5ad for the accelerated package).

The neighboring TSV files are the machine-readable values used by this
report, including the reported bootstrap interval. Raw matrices, user data,
and restricted helper/mapping content are not committed.
