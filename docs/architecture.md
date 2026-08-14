# Architecture and compatibility policy

The pipeline has four explicit stages:

1. Validate a genes-by-samples bulk matrix and an h5ad single-cell object.
2. Materialize only bulk-overlapping single-cell genes plus a synthetic
   remainder feature that preserves each cell's full library total.
3. Build three reference matrices, verify a common identifiable cell-type
   schema, and run five optional deconvolution methods.
4. Select strategies from native target estimates, aggregate branch
   compositions, and emit mode-primary, unclosed, and closed consensus views
   while preserving every native branch result.

The synthetic remainder lets the Monocle3-style size-factor calculation and
MuSiC-style library totals see the full-cell total without densifying a large
cells-by-all-genes h5ad matrix. It is never passed to a deconvolution backend.

## Backends are plugins

Optional algorithm packages are imported lazily. The main wheel contains no
algorithm source, compiled extension, signature/reference data, or EPIC mRNA
scaling data. Missing backends cause an actionable error in strict mode or are
recorded as skipped in permissive mode. Permissive execution does not by
itself permit a consensus from a reduced strategy set: when any requested
branch was skipped, consensus fails unless `allow_partial_consensus=True`
(`--allow-partial-consensus`) is explicitly selected.

EPIC additionally requires an explicit, authorized `mRNA_cell` mapping. The
adapter never attempts to download or silently substitute the restricted
default data. A user-supplied `{"default": 1.0}` mapping is accepted only as
an explicit unscaled custom-reference approximation and is not represented as
equivalent to the official EPIC/R `mRNA_cell_default`. A uniform default of 1
is not the default used by the original R implementation.

EPIC is therefore supported but omitted from the default method tuple. A run
that explicitly requests EPIC without `epic_mrna_cell` fails during
configuration preflight, before bulk loading or reference construction.
Explicit backend selectors allow callers to require Rust rather than accept an
`auto` fallback; requested/resolved engines and native availability are saved
per branch.

EPIC is run with `withOtherCells=True`. The returned `otherCells` compartment
is preserved in the native branch artifact. It is not a consensus target, but
it remains in that branch's denominator during compositional aggregation.

## Reference schema and identifiability

All constructed signature matrices must expose exactly the same set of target
cell types. Column order may differ and is aligned deterministically, but a
missing or additional type is an error. Likewise, every completed backend
branch must contain every target type; missing estimates are not imputed as
zero. A backend may also expose a declared auxiliary type such as EPIC
`otherCells`.

Before fitting, the common-gene signature must contain at least
`max(2, number_of_cell_types)` genes, every target column must have positive
mass on those genes, and the signature must have full column rank. These are
identifiability checks, not silent data repairs.

## MuSiC2 correction

The legacy accelerated MuSiC2 wrapper uses one scalar mean library size for
all types and restores genes in single-cell rather than bulk order. The
corrected reference is the per-cell-type mean expression matrix on raw counts,
with filtered genes in bulk order. Construction enforces the original
DECEPTICONx 20% overlap gate, using the observed common-gene count divided by
the smaller of the bulk and source single-cell gene counts, and requires at
least `max(2, number_of_cell_types)` usable genes. The known-defective path is
not executed as a discarded diagnostic. Provenance distinguishes raw-count
and explicitly allowed normalized-compatibility inputs.

## MuSiC compatibility

DECEPTICON's custom-signature path did not pass the original single-cell
object into MuSiC. For each signature matrix it created five identical pseudo
donors, one pseudo-cell per cell type per donor. This implementation retains
that behavior so the branch count and consensus semantics remain comparable.

## Consensus

Raw backend outputs have incompatible scales: relative proportions,
CIBERSORT-ABS scores, and EPIC fractions whose fit includes `otherCells`.
Consensus therefore keeps two separate numerical views of each branch:

1. The native target columns are used for per-cell-type Pearson correlation
   and strategy selection. In particular, CIBERSORT relative and absolute
   scores are not normalized into the same vectors before selection.
2. For aggregation, each native branch row is closed over *all* of its output
   columns and is then restricted to the target columns. Thus EPIC
   `otherCells` remains in the denominator rather than converting its targets
   into conditional fractions.

Small negative numerical noise is clipped only after material negatives have
been rejected. Zero-mass branch rows are retained as zero and recorded.

The default `r_literal` mode reproduces the defined R `optimal_id` selection
and positional-weighting kernel. The R filename-equivalent strategy order is
method-major: `deconrnaseq`, `epic`, `cibersort_abs`, `cibersort`, and `music`,
each in configured reference order. Undefined correlations and same-reference
positions become zero, the stable row-major ranking selects positions 1 and 3,
and each selected endpoint receives 0.25. Positions are not deduplicated: a
repeated/self/reversed selection accumulates its literal weight. This mode
fixes `n_pairs=2` and uses the unclosed matrix as its primary result, so
consensus rows need not sum to one.

The `corrected` mode instead considers finite positive, unique unordered
cross-reference pairs and uses up to the requested `n_pairs`. It does not pad
a short candidate list with zero or masked positions. If no eligible pair
exists for a cell type, it uses a reference-balanced median: a median across
methods is computed separately for every reference, followed by a median
across references. Its primary result is the row-closed matrix.

DECEPTICON's R custom-output `decostand` loop has hard-coded/index-precedence
errors and does not define a valid normalization result for the full
15-strategy configuration. Here, `r_literal` consequently means literal
compatibility with the selector and positional weights, not execution of that
broken normalization block. Closing every native branch over all native
columns is the explicit aggregation rule in both modes.

Three files make the output contract unambiguous:

- `consensus.tsv` is the selected mode's primary result;
- `consensus_unclosed.tsv` is the pre-final-closure result;
- `consensus_closed.tsv` closes rows with positive mass to one.

The Python compatibility aliases `consensus_raw` and
`consensus_normalized` refer to the unclosed and closed views respectively.

## Audit and output safety

Backend warnings and EPIC `fit_gof` tables are retained in `run.json`.
Installed component versions and selected engines are also recorded. EPIC
mRNA-per-cell values and inline CLI arguments are redacted and replaced by a
SHA-256 configuration fingerprint. Output writing refuses non-empty
directories so stale branch files cannot survive into a newer run.
