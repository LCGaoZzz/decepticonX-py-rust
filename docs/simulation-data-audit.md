# `simulation_data_human.zip` audit

Audited archive SHA-256:
`E92EA54022B2EA31DE1D97304956E44753E7540FD39CA23CB16F700BDE695BCD`.

The archive is not committed because it exceeds GitHub's ordinary single-file
limit and has no redistribution license in the archive.

## Single-cell object

`single_cell_schelker.rda` contains one Biobase `ExpressionSet` with 17,753
genes and 11,759 cells. Its only assay is dense double expression; all nonzero
values are non-integer and range up to about 16.85. It is normalized/log-scale
expression, not raw counts. It can be converted to h5ad for I/O testing only by
using the explicit normalized-input override.

Metadata includes `source`, `donor`, `tsneX1`, `tsneX2`, `cell_id`, and
`cell_type`. There are 12 labels including `Unknown`; the truth file has 11
types, so any exclusion and label aliases must be explicit and preserved in
provenance.

## Bulk corruption

`raw_data/simulation_data.txt` is truncated at exactly 614,400 bytes, lacks a
final newline, and stops alphabetically at gene `AFP`. It contains 362 complete
gene rows plus a partial 363rd row; the missing final 11 fields are read as
NaN. It must fail strict input validation.

Dropping the partial `AFP` row yields a 362-gene by 100-sample matrix suitable
only for orchestration smoke tests. It is not a numerical oracle or a valid
performance benchmark. `cell_fractions.txt` itself is complete (100 samples by
11 cell types, rows summing to one), but accuracy claims must wait for the
complete bulk matrix and raw single-cell counts.

The explicit truth-label aliases used by the evaluation helper are stored in
`examples/schelker_label_map.json`. The dendritic-cell alias is specific to
this archive and should not be reused as a general ontology mapping.
