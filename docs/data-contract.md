# Input data contract

## Single-cell h5ad

- Orientation must be AnnData-standard: cells in `obs`, genes in `var`.
- Raw, non-negative, finite, integer-like counts are required by default.
- `counts_layer="auto"` checks `layers["counts"]`, then count-like `X`, then
  count-like `raw.X`. The chosen source is recorded in run provenance.
- A non-integer `X` is rejected unless `allow_normalized_x=True` is supplied
  explicitly. Such a run is marked non-count and is a compatibility smoke run,
  not a reference-equivalent analysis.
- `obs[cell_type_key]` must exist and contain no missing/blank labels after any
  explicitly requested exclusions.
- Human gene symbols are supported in the first release. Ensembl-like
  `var_names` require an explicit `gene_symbol_key`; duplicate symbols are
  aggregated deterministically.
- Negative values, NaN/Inf, duplicate cell IDs, zero-library cells, and an
  empty bulk-gene intersection are rejected.

MuSiC2-style construction uses per-cell-type mean expression on this count
scale and retains usable genes in bulk order. It applies the original 20%
overlap gate: the common-gene count must be at least 20% of the smaller of the
bulk and source single-cell gene counts. After bulk non-zero/marker filtering,
at least `max(2, number_of_cell_types)` genes must remain.

## Bulk expression

- Text input is genes by samples, with the first column containing gene names.
- CSV and tab-separated files are accepted. Every field must be present,
  numeric, and finite; ragged rows that pandas might otherwise pad with NaN are
  rejected.
- Gene and sample names must be non-empty and unique. Duplicate gene symbols
  should be resolved upstream instead of silently changing the reference.
- Bulk values must be non-negative. Integer counts and normalized expression
  are both permitted, but users must choose CIBERSORT quantile normalization
  according to assay scale.

## Reference and branch compatibility

- Every generated signature matrix must contain exactly the same set of
  target cell types. Different column order is aligned; missing or additional
  types are rejected rather than intersected away.
- For each signature/bulk fit, the common-gene signature must contain at least
  `max(2, number_of_cell_types)` genes, have positive expression mass for
  every target type, and have full column rank. A rank-deficient signature
  cannot identify all requested proportions and is rejected.
- Every completed deconvolution branch must return every target type. Missing
  target columns are errors and are not filled with zero.
- Declared auxiliary outputs are preserved. In particular, EPIC runs with
  `withOtherCells=True`, and its `otherCells` column remains in the native
  branch table and in the denominator used to close that branch for consensus
  aggregation.
- EPIC requires an explicit, authorized, finite positive `mRNA_cell` mapping
  for every target plus `otherCells`, or an explicit positive `"default"`.
  `{"default": 1.0}` is a user-chosen unscaled approximation; it is not the
  official EPIC/R default mapping.

## Consensus data semantics

Strategy selection uses native target estimates, not row-normalized target
submatrices. Aggregation closes each branch over all of its native output
columns before selecting the target columns. This separation preserves method
scale information for Pearson selection while preventing EPIC's unmodelled
compartment from disappearing during aggregation.

The default `r_literal` mode implements the R selector/positional-weighting
kernel and returns the unclosed matrix as `consensus.tsv`. The original R
15-strategy `decostand` loop has an indexing/precedence bug and does not define
a complete literal normalization result, so it is not part of this
compatibility claim. The optional `corrected` mode uses only finite positive,
unique unordered cross-reference pairs, a reference-balanced median when none
exist, and returns the closed matrix as its primary output. Both unclosed and
closed views are always written as `consensus_unclosed.tsv` and
`consensus_closed.tsv`.

If permissive backend mode skips any requested branch, consensus remains
disabled by default. Producing a result from that smaller strategy set requires
the explicit `allow_partial_consensus=True`/`--allow-partial-consensus`
opt-in, because the selected pairs and their weights can change.

## Known CIBERSORT scale behavior

The current optional CIBERSORT port follows the v1.04 heuristic that applies
an inverse log transform when the bulk global maximum is below 50. This cannot
currently be disabled in that backend. Low-range linear counts or scaled data
may therefore be misclassified; the run provenance records the observed bulk
range so users can audit this condition.
