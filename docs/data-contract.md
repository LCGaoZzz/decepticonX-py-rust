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

## Known CIBERSORT scale behavior

The current optional CIBERSORT port follows the v1.04 heuristic that applies
an inverse log transform when the bulk global maximum is below 50. This cannot
currently be disabled in that backend. Low-range linear counts or scaled data
may therefore be misclassified; the run provenance records the observed bulk
range so users can audit this condition.
