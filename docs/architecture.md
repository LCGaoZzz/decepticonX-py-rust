# Architecture and compatibility policy

The pipeline has four explicit stages:

1. Validate a genes-by-samples bulk matrix and an h5ad single-cell object.
2. Materialize only bulk-overlapping single-cell genes plus a synthetic
   remainder feature that preserves each cell's full library total.
3. Build three reference matrices and run five optional deconvolution methods.
4. Normalize branch outputs for comparability and form a DECEPTICON-style
   correlation consensus while preserving every raw branch result.

The synthetic remainder lets the Monocle3-style size-factor calculation and
MuSiC-style library totals see the full-cell total without densifying a large
cells-by-all-genes h5ad matrix. It is never passed to a deconvolution backend.

## Backends are plugins

Optional algorithm packages are imported lazily. The main wheel contains no
algorithm source, compiled extension, signature/reference data, or EPIC mRNA
scaling data. Missing backends cause an actionable error in strict mode or are
recorded as skipped in permissive mode.

EPIC additionally requires an explicit, authorized `mRNA_cell` mapping. The
adapter never attempts to download or silently substitute the restricted
default data. A user-supplied `{"default": 1.0}` mapping is accepted only as
an explicit unscaled custom-reference approximation and is not represented as
equivalent to the original R EPIC cell fractions.

EPIC is therefore supported but omitted from the default method tuple. A run
that explicitly requests EPIC without `epic_mrna_cell` fails during
configuration preflight, before bulk loading or reference construction.
Explicit backend selectors allow callers to require Rust rather than accept an
`auto` fallback; requested/resolved engines and native availability are saved
per branch.

## MuSiC2 correction

The existing accelerated MuSiC2 wrapper uses one scalar mean library size for
all types and restores genes in single-cell rather than bulk order. The public
reference here is the per-cell-type mean expression matrix in filtered bulk
order. The known-defective implementation is not executed as a discarded
diagnostic, avoiding duplicate work. Provenance distinguishes raw-count and
normalized-compatibility inputs.

## MuSiC compatibility

DECEPTICON's custom-signature path did not pass the original single-cell
object into MuSiC. For each signature matrix it created five identical pseudo
donors, one pseudo-cell per cell type per donor. This implementation retains
that behavior so the branch count and consensus semantics remain comparable.

## Consensus

Raw backend outputs have incompatible scales: relative proportions,
CIBERSORT-ABS scores, and EPIC fractions whose fit includes `otherCells`.
The adapter retains the original EPIC `withOtherCells=True` fit but exposes
only signature cell-type columns through the common branch contract. Consensus
views clip numerical negative noise and normalize each sample row. For each cell type, correlations are
computed across samples; self-pairs and branches using the same reference are
masked. The highest unique cross-reference pairs are averaged. A deterministic
median fallback is used when correlation is undefined or too few branches are
available. The final consensus is row-normalized; the pre-final-normalization
matrix is also returned as `consensus_raw`.

## Audit and output safety

Backend warnings and EPIC `fit_gof` tables are retained in `run.json`.
Installed component versions and selected engines are also recorded. EPIC
mRNA-per-cell values and inline CLI arguments are redacted and replaced by a
SHA-256 configuration fingerprint. Output writing refuses non-empty
directories so stale branch files cannot survive into a newer run.
