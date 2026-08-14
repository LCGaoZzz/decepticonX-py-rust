"""Strict, memory-conscious input handling for decepticonX.

The public functions in this module deliberately do less guessing than most
expression-matrix readers.  Bulk matrices must be rectangular genes x samples
tables, while single-cell data are read as cells x genes from AnnData.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
import csv
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any
import warnings as python_warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .exceptions import InputValidationError


REMAINDER_GENE = "__DECEPTICONX_REMAINDER__"
_HUMAN_SPECIES = frozenset(
    {"hs", "human", "homo sapiens", "homo_sapiens", "9606"}
)
_ENSEMBL_GENE_RE = re.compile(r"^ENSG\d+(?:\.\d+)?$", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class PreparedSingleCell:
    """The compact, in-memory AnnData passed to reference builders.

    ``adata.X`` contains one column per ``common_genes`` entry followed by one
    synthetic remainder column.  Consequently, every output row has the same
    sum as that cell had across *all* genes in the selected source matrix.
    """

    adata: ad.AnnData
    common_genes: tuple[str, ...]
    source_layer: str
    normalized_input: bool
    cell_type_key: str
    sample_key: str | None
    warnings: tuple[str, ...]
    provenance: dict[str, Any]

    @property
    def remainder_gene(self) -> str:
        """Return the collision-safe synthetic gene name used in ``adata``."""

        return str(self.provenance["remainder_gene"])


def _fail(message: str) -> InputValidationError:
    return InputValidationError(message)


def validate_bulk_expression(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and copy an already loaded genes x samples bulk matrix."""

    if not isinstance(frame, pd.DataFrame):
        raise TypeError("bulk expression must be a pandas DataFrame")
    if frame.shape[0] == 0 or frame.shape[1] == 0:
        raise _fail("bulk expression must contain at least one gene and one sample")

    genes = pd.Index([str(value).strip() for value in frame.index], dtype="object")
    samples = pd.Index([str(value).strip() for value in frame.columns], dtype="object")
    if any(not value for value in genes):
        raise _fail("bulk expression contains an empty gene name")
    if any(not value for value in samples):
        raise _fail("bulk expression contains an empty sample name")
    if genes.has_duplicates:
        duplicates = genes[genes.duplicated(keep=False)].unique().tolist()
        raise _fail(f"bulk expression contains duplicate gene names: {duplicates[:5]}")
    if samples.has_duplicates:
        duplicates = samples[samples.duplicated(keep=False)].unique().tolist()
        raise _fail(f"bulk expression contains duplicate sample names: {duplicates[:5]}")

    try:
        numeric = frame.apply(pd.to_numeric, errors="raise").astype(np.float64)
    except (TypeError, ValueError) as exc:
        raise _fail("bulk expression contains a non-numeric value") from exc
    values = numeric.to_numpy(copy=False)
    if not np.isfinite(values).all():
        raise _fail("bulk expression contains missing, NaN, or infinite values")
    if (values < 0).any():
        raise _fail("bulk expression contains negative values")
    numeric.index = genes
    numeric.columns = samples
    return numeric


def load_bulk_expression(path: str | Path) -> pd.DataFrame:
    """Read a strict CSV/TSV bulk matrix with genes in rows and samples in columns.

    The first column is always interpreted as the gene identifier.  A manual
    rectangularity check is performed before constructing a DataFrame so that
    malformed/ragged input cannot be silently padded with NaNs by pandas.
    """

    source = Path(path)
    suffixes = [suffix.lower() for suffix in source.suffixes]
    if suffixes[-1:] == [".csv"]:
        delimiter = ","
    elif suffixes[-1:] in ([".tsv"], [".txt"]):
        delimiter = "\t"
    else:
        raise _fail("bulk expression path must end in .csv, .tsv, or .txt")
    if not source.is_file():
        raise _fail(f"bulk expression file does not exist: {source}")

    try:
        handle = source.open("r", encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise _fail(f"could not open bulk expression file: {source}") from exc

    with handle:
        reader = csv.reader(handle, delimiter=delimiter, strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise _fail("bulk expression file is empty") from exc
        except csv.Error as exc:
            raise _fail(f"invalid bulk expression table: {exc}") from exc
        if len(header) < 1:
            raise _fail("bulk expression header must contain sample names")

        genes: list[str] = []
        rows: list[list[float]] = []

        try:
            first_row = next(reader)
        except StopIteration as exc:
            raise _fail("bulk expression file contains a header but no gene rows") from exc
        except csv.Error as exc:
            raise _fail(f"invalid bulk expression table: {exc}") from exc

        # Both layouts are common in R-generated expression files:
        #   gene<TAB>s1<TAB>s2    (header includes an index label)
        #   s1<TAB>s2             (first gene column has no header)
        if len(first_row) == len(header):
            if len(header) < 2:
                raise _fail("bulk expression header must contain sample names")
            sample_names = [value.strip() for value in header[1:]]
            expected_width = len(header)
        elif len(first_row) == len(header) + 1:
            sample_names = [value.strip() for value in header]
            expected_width = len(header) + 1
        else:
            raise _fail(
                "ragged bulk expression table at line 2: expected either "
                f"{len(header)} or {len(header) + 1} fields, got {len(first_row)}"
            )

        def parse_row(row: list[str], line_number: int) -> None:
            if len(row) != expected_width:
                raise _fail(
                    "ragged bulk expression table at line "
                    f"{line_number}: expected {expected_width} fields, got {len(row)}"
                )
            gene = row[0].strip()
            if not gene:
                raise _fail(f"empty bulk gene name at line {line_number}")
            numeric_row: list[float] = []
            for column_number, value in enumerate(row[1:], start=2):
                token = value.strip()
                if not token:
                    raise _fail(
                        "missing bulk expression value at line "
                        f"{line_number}, column {column_number}"
                    )
                try:
                    number = float(token)
                except ValueError as exc:
                    raise _fail(
                        "non-numeric bulk expression value at line "
                        f"{line_number}, column {column_number}: {token!r}"
                    ) from exc
                if not np.isfinite(number):
                    raise _fail(
                        "non-finite bulk expression value at line "
                        f"{line_number}, column {column_number}"
                    )
                numeric_row.append(number)
            genes.append(gene)
            rows.append(numeric_row)

        try:
            parse_row(first_row, 2)
            for line_number, row in enumerate(reader, start=3):
                parse_row(row, line_number)
        except csv.Error as exc:
            raise _fail(f"invalid bulk expression table: {exc}") from exc

    frame = pd.DataFrame(rows, index=genes, columns=sample_names, dtype=np.float64)
    return validate_bulk_expression(frame)


# Short aliases are intentionally kept for callers migrating from early drafts.
load_bulk = load_bulk_expression
read_bulk = load_bulk_expression


def _normalise_species(species: str) -> str:
    value = str(species).strip().lower()
    if value not in _HUMAN_SPECIES:
        raise _fail(
            "only human gene symbols are supported (species must be one of "
            "'hs', 'human', 'Homo sapiens', or '9606')"
        )
    return "human"


def _validate_obs_labels(
    adata: ad.AnnData, key: str, *, parameter_name: str
) -> np.ndarray:
    if key not in adata.obs.columns:
        raise _fail(f"{parameter_name} {key!r} is missing from adata.obs")
    labels = adata.obs[key]
    missing = labels.isna().to_numpy()
    as_strings = labels.astype("string").str.strip()
    missing |= as_strings.isna().to_numpy()
    missing |= as_strings.fillna("").eq("").to_numpy()
    if missing.any():
        locations = np.flatnonzero(missing)[:5].tolist()
        raise _fail(
            f"adata.obs[{key!r}] contains missing/empty labels "
            f"at cell positions {locations}"
        )
    return as_strings.astype(str).to_numpy()


def _matrix_rows(matrix: Any, n_obs: int, n_vars: int) -> Iterable[Any]:
    """Yield bounded row chunks without converting a backed matrix wholesale."""

    # Keep dense chunks near 16 MiB. Sparse chunks stay cheap and use at most
    # 1024 rows; accessing ``dtype`` does not materialize backed datasets.
    itemsize = np.dtype(getattr(matrix, "dtype", np.float64)).itemsize
    rows_for_dense_budget = max(1, (16 * 1024 * 1024) // max(1, n_vars * itemsize))
    chunk_size = min(1024, rows_for_dense_budget)
    for start in range(0, n_obs, chunk_size):
        stop = min(n_obs, start + chunk_size)
        yield matrix[start:stop, :]


def _numeric_profile(matrix: Any, n_obs: int, n_vars: int) -> tuple[bool, bool, bool]:
    """Return (finite, nonnegative, integer_like), scanning in bounded chunks."""

    finite = True
    nonnegative = True
    integer_like = True
    for chunk in _matrix_rows(matrix, n_obs, n_vars):
        values = chunk.data if sparse.issparse(chunk) else np.asarray(chunk)
        if values.size == 0:
            continue
        if not np.issubdtype(values.dtype, np.number):
            return False, False, False
        finite = finite and bool(np.isfinite(values).all())
        if not finite:
            return False, False, False
        nonnegative = nonnegative and bool((values >= 0).all())
        if not nonnegative:
            integer_like = False
            continue
        if not np.issubdtype(values.dtype, np.integer):
            integer_like = integer_like and bool(
                np.allclose(values, np.rint(values), rtol=0.0, atol=1e-6)
            )
    return finite, nonnegative, integer_like


def _source_var(adata: ad.AnnData, source_layer: str) -> pd.DataFrame:
    if source_layer == "raw.X":
        assert adata.raw is not None
        return adata.raw.var
    return adata.var


def _resolve_source(
    adata: ad.AnnData,
    counts_layer: str,
    allow_normalized_x: bool,
) -> tuple[Any, str, bool, tuple[bool, bool, bool]]:
    requested = str(counts_layer).strip()
    if not requested:
        raise _fail("counts_layer cannot be empty")

    def profile(matrix: Any, n_vars: int) -> tuple[bool, bool, bool]:
        return _numeric_profile(matrix, adata.n_obs, n_vars)

    if requested != "auto":
        if requested == "X":
            matrix, label, n_vars = adata.X, "X", adata.n_vars
        elif requested == "raw.X":
            if adata.raw is None:
                raise _fail("counts_layer='raw.X' requested but adata.raw is absent")
            matrix, label, n_vars = adata.raw.X, "raw.X", adata.raw.n_vars
        else:
            if requested not in adata.layers:
                raise _fail(f"counts layer {requested!r} is missing from adata.layers")
            matrix, label, n_vars = (
                adata.layers[requested],
                f"layers[{requested!r}]",
                adata.n_vars,
            )
        stats = profile(matrix, n_vars)
        if not stats[0]:
            raise _fail(f"selected single-cell source {label} contains NaN/Inf or non-numeric data")
        if not stats[1]:
            raise _fail(f"selected single-cell source {label} contains negative values")
        if stats[2]:
            return matrix, label, False, stats
        if label == "X" and allow_normalized_x:
            return matrix, label, True, stats
        raise _fail(
            f"selected single-cell source {label} is not integer-like count data; "
            "only non-integer X can be accepted with allow_normalized_x=True"
        )

    if "counts" in adata.layers:
        matrix = adata.layers["counts"]
        stats = profile(matrix, adata.n_vars)
        if not stats[0]:
            raise _fail("adata.layers['counts'] contains NaN/Inf or non-numeric data")
        if not stats[1]:
            raise _fail("adata.layers['counts'] contains negative values")
        if not stats[2]:
            raise _fail("adata.layers['counts'] is not integer-like count data")
        return matrix, "layers['counts']", False, stats

    x_stats = profile(adata.X, adata.n_vars)
    if x_stats == (True, True, True):
        return adata.X, "X", False, x_stats

    if adata.raw is not None:
        raw_stats = profile(adata.raw.X, adata.raw.n_vars)
        if raw_stats == (True, True, True):
            return adata.raw.X, "raw.X", False, raw_stats

    if not x_stats[0]:
        raise _fail("adata.X contains NaN/Inf or non-numeric data and no valid raw counts exist")
    if not x_stats[1]:
        raise _fail("adata.X contains negative values and no valid raw counts exist")
    if allow_normalized_x:
        return adata.X, "X", True, x_stats
    raise _fail(
        "adata.X is not integer-like count data and no count-like adata.raw.X exists; "
        "set allow_normalized_x=True only when normalized X is intentional"
    )


def _gene_symbols(
    adata: ad.AnnData,
    source_layer: str,
    gene_symbol_key: str | None,
) -> np.ndarray:
    var = _source_var(adata, source_layer)
    if gene_symbol_key is None:
        values = pd.Series(var.index.astype(str), index=var.index, dtype="string")
        ensembl = values.str.match(_ENSEMBL_GENE_RE, na=False)
        if bool(ensembl.any()):
            examples = values[ensembl].head(3).tolist()
            raise _fail(
                "single-cell genes appear to be Ensembl IDs; provide gene_symbol_key "
                f"to map them to human symbols (examples: {examples})"
            )
    else:
        if gene_symbol_key not in var.columns:
            raise _fail(
                f"gene_symbol_key {gene_symbol_key!r} is missing from the selected source var"
            )
        values = var[gene_symbol_key].astype("string")

    stripped = values.str.strip()
    invalid = stripped.isna() | stripped.fillna("").eq("")
    # Tabs/newlines would make provenance and tabular outputs ambiguous.
    invalid |= stripped.fillna("").str.contains(r"[\t\r\n]", regex=True)
    if bool(invalid.any()):
        locations = np.flatnonzero(invalid.to_numpy())[:5].tolist()
        raise _fail(f"gene symbols are missing or invalid at source positions {locations}")
    symbols = stripped.astype(str).to_numpy()
    if gene_symbol_key is not None:
        ensembl = np.fromiter(
            (bool(_ENSEMBL_GENE_RE.fullmatch(value)) for value in symbols),
            dtype=bool,
            count=len(symbols),
        )
        if ensembl.any():
            examples = symbols[ensembl][:3].tolist()
            raise _fail(
                f"gene_symbol_key must contain human symbols, not Ensembl IDs: {examples}"
            )
    return symbols


def _bulk_genes(bulk: pd.DataFrame | Sequence[str] | pd.Index) -> tuple[str, ...]:
    if isinstance(bulk, pd.DataFrame):
        return tuple(validate_bulk_expression(bulk).index.astype(str))
    if isinstance(bulk, (str, bytes)):
        raise TypeError("bulk must be a DataFrame or a sequence of gene names, not a string")
    try:
        genes = tuple(str(value).strip() for value in bulk)
    except TypeError as exc:
        raise TypeError("bulk must be a DataFrame or a sequence of gene names") from exc
    if not genes or any(not value for value in genes):
        raise _fail("bulk gene names must be non-empty")
    index = pd.Index(genes)
    if index.has_duplicates:
        duplicates = index[index.duplicated(keep=False)].unique().tolist()
        raise _fail(f"bulk gene names contain duplicates: {duplicates[:5]}")
    return genes


def _unique_remainder_name(symbols: np.ndarray, bulk_genes: Sequence[str]) -> str:
    occupied = set(symbols)
    occupied.update(bulk_genes)
    candidate = REMAINDER_GENE
    suffix = 1
    while candidate in occupied:
        candidate = f"{REMAINDER_GENE}_{suffix}"
        suffix += 1
    return candidate


def _take_rows(matrix: Any, rows: np.ndarray) -> Any:
    """Take sorted rows from dense/sparse/backed matrices with a safe fallback."""

    try:
        return matrix[rows, :]
    except (IndexError, TypeError, ValueError):
        # Some h5py-backed array implementations only support slices.  The
        # caller bounds chunk sizes, so concatenating individual selected rows
        # remains memory-safe.
        pieces = [matrix[int(row) : int(row) + 1, :] for row in rows]
        if not pieces:
            return np.empty((0, matrix.shape[1]), dtype=getattr(matrix, "dtype", float))
        if all(sparse.issparse(piece) for piece in pieces):
            return sparse.vstack(pieces, format="csr")
        return np.concatenate([np.asarray(piece) for piece in pieces], axis=0)


def _materialize_common_with_remainder(
    matrix: Any,
    selected_rows: np.ndarray,
    n_source_genes: int,
    source_symbols: np.ndarray,
    common_genes: tuple[str, ...],
) -> tuple[Any, np.ndarray]:
    """Materialize only common genes plus a library-preserving remainder."""

    common_position = {gene: position for position, gene in enumerate(common_genes)}
    source_positions = np.fromiter(
        (i for i, symbol in enumerate(source_symbols) if symbol in common_position),
        dtype=np.int64,
    )
    mapped_columns = np.fromiter(
        (common_position[source_symbols[i]] for i in source_positions),
        dtype=np.int64,
    )
    mapping = sparse.csr_matrix(
        (
            np.ones(len(source_positions), dtype=np.float64),
            (np.arange(len(source_positions), dtype=np.int64), mapped_columns),
        ),
        shape=(len(source_positions), len(common_genes)),
    )

    output_chunks: list[Any] = []
    library_totals: list[np.ndarray] = []
    itemsize = np.dtype(getattr(matrix, "dtype", np.float64)).itemsize
    rows_for_dense_budget = max(
        1, (16 * 1024 * 1024) // max(1, n_source_genes * itemsize)
    )
    chunk_size = min(1024, rows_for_dense_budget)

    for start in range(0, len(selected_rows), chunk_size):
        rows = selected_rows[start : start + chunk_size]
        raw_chunk = _take_rows(matrix, rows)
        is_sparse = sparse.issparse(raw_chunk)
        if is_sparse:
            raw_chunk = raw_chunk.tocsr()
            totals = np.asarray(raw_chunk.sum(axis=1), dtype=np.float64).ravel()
            selected = raw_chunk[:, source_positions]
            common = (selected @ mapping).tocsr()
            selected_totals = np.asarray(common.sum(axis=1), dtype=np.float64).ravel()
        else:
            raw_chunk = np.asarray(raw_chunk)
            totals = np.asarray(raw_chunk.sum(axis=1, dtype=np.float64)).ravel()
            selected = np.asarray(raw_chunk[:, source_positions], dtype=np.float64)
            common = np.asarray(mapping.T @ selected.T).T
            selected_totals = common.sum(axis=1, dtype=np.float64)

        remainder = totals - selected_totals
        tolerance = np.maximum(1e-8, np.abs(totals) * 1e-10)
        if np.any(remainder < -tolerance):
            raise _fail("internal gene aggregation error produced a negative remainder")
        remainder[(remainder < 0) & (remainder >= -tolerance)] = 0.0
        if is_sparse:
            output = sparse.hstack(
                [common, sparse.csr_matrix(remainder[:, np.newaxis])], format="csr"
            )
        else:
            output = np.column_stack([common, remainder])
        output_chunks.append(output)
        library_totals.append(totals)

    if all(sparse.issparse(chunk) for chunk in output_chunks):
        output_matrix: Any = sparse.vstack(output_chunks, format="csr")
    else:
        output_matrix = np.concatenate(
            [chunk.toarray() if sparse.issparse(chunk) else chunk for chunk in output_chunks],
            axis=0,
        )
    return output_matrix, np.concatenate(library_totals)


def prepare_single_cell(
    single_cell: str | Path | ad.AnnData,
    bulk: pd.DataFrame | Sequence[str] | pd.Index,
    *,
    cell_type_key: str = "cell_type",
    sample_key: str | None = None,
    counts_layer: str = "auto",
    gene_symbol_key: str | None = None,
    exclude_cell_types: Iterable[str] = (),
    allow_normalized_x: bool = False,
    species: str = "hs",
) -> PreparedSingleCell:
    """Validate and compact an h5ad/AnnData input for reference construction.

    Paths are opened in backed read-only mode and closed after the compact
    in-memory AnnData is produced.  At no point is the complete source matrix
    converted to a dense array.
    """

    _normalise_species(species)
    bulk_gene_names = _bulk_genes(bulk)
    opened_here = not isinstance(single_cell, ad.AnnData)
    source_path: str | None = None
    if opened_here:
        source = Path(single_cell)
        if source.suffix.lower() != ".h5ad":
            raise _fail("single-cell input path must end in .h5ad")
        if not source.is_file():
            raise _fail(f"single-cell h5ad file does not exist: {source}")
        source_path = str(source.resolve())
        try:
            adata = ad.read_h5ad(source, backed="r")
        except (OSError, ValueError) as exc:
            raise _fail(f"could not read h5ad input: {source}") from exc
    else:
        adata = single_cell
        filename = getattr(adata, "filename", None)
        source_path = str(filename) if filename is not None else None

    try:
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise _fail("single-cell AnnData must contain cells and genes")
        obs_names = pd.Index(adata.obs_names.astype(str))
        if obs_names.has_duplicates:
            duplicates = obs_names[obs_names.duplicated(keep=False)].unique().tolist()
            raise _fail(f"single-cell input contains duplicate cell IDs: {duplicates[:5]}")
        if any(not value.strip() for value in obs_names):
            raise _fail("single-cell input contains an empty cell ID")
        labels = _validate_obs_labels(
            adata, cell_type_key, parameter_name="cell_type_key"
        )
        sample_labels: np.ndarray | None = None
        if sample_key is not None:
            sample_labels = _validate_obs_labels(
                adata, sample_key, parameter_name="sample_key"
            )
        excluded = frozenset(str(value).strip() for value in exclude_cell_types)
        keep = ~np.isin(labels, tuple(excluded)) if excluded else np.ones(adata.n_obs, bool)
        selected_rows = np.flatnonzero(keep).astype(np.int64, copy=False)
        if selected_rows.size == 0:
            raise _fail("excluding cell types removed every single cell")

        matrix, source_layer, normalized_input, _ = _resolve_source(
            adata, counts_layer, allow_normalized_x
        )
        symbols = _gene_symbols(adata, source_layer, gene_symbol_key)
        source_gene_set = set(symbols)
        common_genes = tuple(gene for gene in bulk_gene_names if gene in source_gene_set)
        if not common_genes:
            raise _fail("bulk and single-cell inputs have no common human gene symbols")
        remainder_gene = _unique_remainder_name(symbols, bulk_gene_names)

        output_matrix, library_totals = _materialize_common_with_remainder(
            matrix,
            selected_rows,
            len(symbols),
            symbols,
            common_genes,
        )
        if (library_totals <= 0).any():
            locations = np.flatnonzero(library_totals <= 0)[:5].tolist()
            raise _fail(f"single-cell input contains zero-library cells at selected positions {locations}")
        selected_obs = adata.obs.iloc[selected_rows].copy()
        # Validation canonicalises leading/trailing whitespace.  Persist the
        # canonical labels in the compact AnnData so reference construction
        # cannot accidentally recreate distinct " T " and "T" groups.
        selected_obs[cell_type_key] = labels[selected_rows]
        if sample_key is not None and sample_labels is not None:
            selected_obs[sample_key] = sample_labels[selected_rows]
        output_var = pd.DataFrame(
            {
                "is_remainder": [False] * len(common_genes) + [True],
            },
            index=pd.Index((*common_genes, remainder_gene), name="gene_symbol"),
        )
        output = ad.AnnData(X=output_matrix, obs=selected_obs, var=output_var)

        duplicate_count = int(pd.Index(symbols).duplicated(keep=False).sum())
        recorded_warnings: list[str] = []
        if normalized_input:
            message = (
                "adata.X is non-integer and was accepted because "
                "allow_normalized_x=True; library totals are totals in normalized space"
            )
            recorded_warnings.append(message)
            python_warnings.warn(message, RuntimeWarning, stacklevel=2)

        provenance: dict[str, Any] = {
            "source_path": source_path,
            "source_layer": source_layer,
            "counts_layer_requested": str(counts_layer),
            "normalized_input": normalized_input,
            "species": "human",
            "gene_symbol_key": gene_symbol_key,
            "cell_type_key": cell_type_key,
            "sample_key": sample_key,
            "excluded_cell_types": sorted(excluded),
            "n_cells_input": int(adata.n_obs),
            "n_cells_output": int(output.n_obs),
            "n_genes_source": int(len(symbols)),
            "n_common_genes": int(len(common_genes)),
            "duplicate_symbol_source_columns": duplicate_count,
            "remainder_gene": remainder_gene,
            "library_total_min": float(library_totals.min()),
            "library_total_max": float(library_totals.max()),
        }
        return PreparedSingleCell(
            adata=output,
            common_genes=common_genes,
            source_layer=source_layer,
            normalized_input=normalized_input,
            cell_type_key=cell_type_key,
            sample_key=sample_key,
            warnings=tuple(recorded_warnings),
            provenance=provenance,
        )
    finally:
        if opened_here:
            file_manager = getattr(adata, "file", None)
            if file_manager is not None:
                file_manager.close()


# Descriptive alias for API discoverability.
prepare_h5ad = prepare_single_cell


__all__ = [
    "PreparedSingleCell",
    "REMAINDER_GENE",
    "load_bulk",
    "load_bulk_expression",
    "prepare_h5ad",
    "prepare_single_cell",
    "read_bulk",
    "validate_bulk_expression",
]
