"""Reference-expression matrix construction.

The heavy DECEPTICONx-compatible builders live in the optional
``decepticon_fast`` package.  It is imported only when a caller asks to
construct references, keeping importing :mod:`decepticonx` lightweight.

``PreparedSingleCell`` normally contains a reduced AnnData matrix with the
bulk/scRNA-seq intersection followed by one synthetic remainder gene.  The
remainder stores all expression outside the intersection and therefore lets
normalisation algorithms retain the original per-cell library sizes without
materialising the full single-cell matrix here.
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from .exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    InputValidationError,
)
from .models import DEFAULT_REFERENCES

if TYPE_CHECKING:  # pragma: no cover - imported only by type checkers
    from .io import PreparedSingleCell


REFERENCE_BUILDERS = ("bayesprism", "monocle3", "music2")


def _load_fast() -> Any:
    try:
        return import_module("decepticon_fast")
    except (ImportError, OSError) as exc:
        raise BackendUnavailableError(
            "Reference construction requires the optional 'decepticon-fast' "
            "package. Install the pinned source or a compatible wheel: "
            "python -m pip install \"git+https://github.com/LCGaoZzz/"
            "decepticon-fast.git@a98b66c7da5fad81ba6a1eaddb9c826daac31fad\""
        ) from exc


def _as_dense(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        result = matrix.toarray()
    else:
        result = np.asarray(matrix)
    return np.asarray(result, dtype=np.float64)


def _take_columns(matrix: Any, positions: np.ndarray) -> np.ndarray:
    """Return selected AnnData variables as a dense cells-by-genes array."""

    return _as_dense(matrix[:, positions])


def _sum_columns(matrix: Any, positions: np.ndarray) -> np.ndarray:
    if positions.size == 0:
        return np.zeros(matrix.shape[0], dtype=np.float64)
    selected = matrix[:, positions]
    summed = selected.sum(axis=1)
    return np.asarray(summed, dtype=np.float64).reshape(-1)


def _normalise_requested(references: Iterable[str]) -> tuple[str, ...]:
    requested = tuple(str(name).lower() for name in references)
    if not requested:
        raise InputValidationError("At least one reference builder is required.")
    unknown = sorted(set(requested).difference(REFERENCE_BUILDERS))
    if unknown:
        raise InputValidationError(
            f"Unknown reference builder(s): {', '.join(unknown)}. "
            f"Supported builders: {', '.join(REFERENCE_BUILDERS)}."
        )
    if len(set(requested)) != len(requested):
        raise InputValidationError("Reference builder names must be unique.")
    return requested


def _prepare_frames(
    prepared: "PreparedSingleCell", bulk: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.Series, str]:
    """Create common-only and common-plus-remainder working matrices.

    ``PreparedSingleCell.adata`` is expected to be reduced already.  The
    fallback remainder calculation also makes this function safe for callers
    that construct the dataclass around a full AnnData object themselves.
    """

    if not isinstance(bulk, pd.DataFrame):
        raise InputValidationError("bulk must be a pandas DataFrame (genes x samples).")
    if bulk.empty:
        raise InputValidationError("bulk must contain at least one gene and sample.")
    if bulk.index.has_duplicates:
        raise InputValidationError("bulk gene identifiers must be unique.")

    bulk_work = bulk.copy(deep=False)
    bulk_work.index = pd.Index([str(value) for value in bulk.index], name=bulk.index.name)
    if bulk_work.index.has_duplicates:
        raise InputValidationError(
            "bulk gene identifiers must remain unique after string conversion."
        )
    try:
        bulk_values = bulk_work.to_numpy(dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise InputValidationError("bulk expression values must be numeric.") from exc
    if not np.isfinite(bulk_values).all():
        raise InputValidationError("bulk expression values must all be finite.")

    adata = prepared.adata
    var_names = pd.Index([str(value) for value in adata.var_names])
    obs_names = pd.Index([str(value) for value in adata.obs_names])
    if var_names.has_duplicates:
        raise InputValidationError("Single-cell gene identifiers must be unique.")
    if obs_names.has_duplicates:
        raise InputValidationError("Single-cell observation identifiers must be unique.")

    remainder_gene = str(prepared.remainder_gene)
    declared_common = {str(value) for value in prepared.common_genes}
    # Always impose bulk order.  This is important for MuSiC2 and also makes
    # every generated artifact deterministic when the input AnnData order
    # differs from the bulk order.
    common = [
        gene
        for gene in bulk_work.index
        if gene != remainder_gene and gene in declared_common and gene in var_names
    ]
    if not common:
        raise InputValidationError(
            "No declared common genes are present in both bulk and single-cell data."
        )

    positions = var_names.get_indexer(common)
    common_values = _take_columns(adata.X, positions)
    if not np.isfinite(common_values).all():
        raise InputValidationError("Single-cell expression values must all be finite.")

    if remainder_gene in var_names:
        remainder_position = var_names.get_loc(remainder_gene)
        remainder = _take_columns(
            adata.X, np.asarray([remainder_position], dtype=np.intp)
        )[:, 0]
    else:
        # Compatibility path for a full, unreduced AnnData.  Sum everything
        # outside the selected common genes into the synthetic remainder.
        common_positions = set(int(value) for value in positions)
        other_positions = np.asarray(
            [i for i in range(len(var_names)) if i not in common_positions],
            dtype=np.intp,
        )
        remainder = _sum_columns(adata.X, other_positions)
    if not np.isfinite(remainder).all():
        raise InputValidationError("Single-cell remainder values must all be finite.")

    sc_common = pd.DataFrame(
        common_values.T,
        index=pd.Index(common, name=bulk_work.index.name),
        columns=obs_names,
    )
    sc_with_remainder = pd.concat(
        [
            sc_common,
            pd.DataFrame(
                remainder.reshape(1, -1),
                index=pd.Index([remainder_gene], name=sc_common.index.name),
                columns=obs_names,
            ),
        ],
        axis=0,
    )

    cell_type_key = str(prepared.cell_type_key)
    if cell_type_key not in adata.obs:
        raise InputValidationError(
            f"Prepared AnnData is missing cell type column {cell_type_key!r}."
        )
    raw_labels = adata.obs[cell_type_key]
    if raw_labels.isna().any():
        raise InputValidationError("Cell type labels may not contain missing values.")
    labels = pd.Series(
        raw_labels.astype(str).to_numpy(),
        index=obs_names,
        name=cell_type_key,
    )
    if (labels.str.len() == 0).any():
        raise InputValidationError("Cell type labels may not be empty strings.")

    bulk_common = bulk_work.loc[common].copy()
    return sc_common, sc_with_remainder, bulk_common, labels, remainder_gene


def _coerce_signature(frame: Any, name: str) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        raise BackendExecutionError(
            f"{name} reference builder returned {type(frame).__name__}, not a DataFrame."
        )
    if frame.empty:
        raise BackendExecutionError(f"{name} reference builder returned an empty matrix.")
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise BackendExecutionError(
            f"{name} reference builder returned duplicate genes or cell types."
        )
    try:
        values = frame.to_numpy(dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise BackendExecutionError(
            f"{name} reference builder returned non-numeric values."
        ) from exc
    if not np.isfinite(values).all():
        raise BackendExecutionError(
            f"{name} reference builder returned non-finite values."
        )
    result = pd.DataFrame(
        values,
        index=pd.Index([str(value) for value in frame.index], name=frame.index.name),
        columns=pd.Index(
            [str(value) for value in frame.columns], name=frame.columns.name
        ),
    )
    return result


def _music2_mean_reference(
    sc_common: pd.DataFrame,
    bulk_common: pd.DataFrame,
    labels: pd.Series,
    *,
    markers: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Return the corrected MuSiC2-style mean-expression reference.

    DECEPTICONx's small ``MuSiC2_base`` wrapper contains two independent
    defects: it multiplies all cell types by one scalar library-size value,
    and restores marker rows in single-cell rather than bulk order.  The
    actual reference needed by downstream deconvolution is simply mean raw
    expression for each gene and cell type.  Filtering and order follow the
    non-zero bulk rows (and, when supplied, the marker set).
    """

    bulk_values = bulk_common.to_numpy(dtype=np.float64)
    nonzero = bulk_values.mean(axis=1) != 0
    marker_set = None if markers is None else {str(value) for value in markers}
    genes = [
        str(gene)
        for gene, keep in zip(bulk_common.index, nonzero, strict=True)
        if keep and (marker_set is None or str(gene) in marker_set)
    ]
    if not genes:
        raise InputValidationError(
            "MuSiC2 has no genes after bulk non-zero and marker filtering."
        )

    label_values = labels.loc[sc_common.columns].to_numpy(dtype=object)
    cell_types = list(pd.unique(label_values))
    minimum = max(2, len(cell_types))
    if len(genes) < minimum:
        raise InputValidationError(
            "MuSiC2 has too few genes after bulk non-zero and marker filtering "
            f"({len(genes)} found; at least {minimum} required for "
            f"{len(cell_types)} cell type(s))."
        )
    values = sc_common.loc[genes].to_numpy(dtype=np.float64)
    result = np.empty((len(genes), len(cell_types)), dtype=np.float64)
    for column, cell_type in enumerate(cell_types):
        mask = label_values == cell_type
        # Every level came from label_values, so this is non-empty by
        # construction.  Keeping the explicit guard produces a useful error
        # if a future categorical implementation changes that invariant.
        if not np.any(mask):  # pragma: no cover - defensive invariant
            raise InputValidationError(f"Cell type {cell_type!r} has no cells.")
        result[:, column] = values[:, mask].mean(axis=1)
    return pd.DataFrame(result, index=genes, columns=cell_types)


def build_references(
    prepared: "PreparedSingleCell",
    bulk: pd.DataFrame,
    *,
    references: Iterable[str] = DEFAULT_REFERENCES,
    species: str = "hs",
    markers: Sequence[str] | None = None,
) -> tuple[dict[str, pd.DataFrame], dict[str, Any]]:
    """Build the requested reference matrices and diagnostics.

    Parameters
    ----------
    prepared:
        Validated single-cell input prepared by :func:`prepare_single_cell`.
    bulk:
        Bulk expression matrix with genes in rows and samples in columns.
    references:
        Any unique subset of ``bayesprism``, ``monocle3``, and ``music2``.
        The default is ``bayesprism`` plus ``music2``; ``monocle3`` remains
        available when requested explicitly.
    species:
        Currently only ``"hs"`` (human) is supported because the bundled
        BayesPrism annotations are human-specific.
    markers:
        Optional MuSiC2 marker genes.  Their order never overrides bulk order.
    """

    requested = _normalise_requested(references)
    if str(species).lower() != "hs":
        raise InputValidationError(
            "Only species='hs' is currently supported for reference construction."
        )

    fast = _load_fast()
    sc_common, sc_remainder, bulk_common, labels, remainder_gene = _prepare_frames(
        prepared, bulk
    )
    diagnostics: dict[str, Any] = {
        "backend": fast.get_backend() if callable(getattr(fast, "get_backend", None)) else "unknown",
        "common_gene_count": int(sc_common.shape[0]),
        "remainder_gene": remainder_gene,
        "references": {},
    }
    output: dict[str, pd.DataFrame] = {}

    for name in requested:
        if name == "bayesprism":
            try:
                built = fast.bayesprism_base(
                    bulk_common, sc_common, labels, species="hs"
                )
            except Exception as exc:
                raise BackendExecutionError(
                    "decepticon-fast BayesPrism reference construction failed."
                ) from exc
            result = _coerce_signature(built, "BayesPrism")

        elif name == "monocle3":
            metadata = labels.to_frame(name=labels.name or "cell_type")
            try:
                built = fast.monocle3_base(sc_remainder, metadata)
            except Exception as exc:
                raise BackendExecutionError(
                    "decepticon-fast Monocle3 reference construction failed."
                ) from exc
            built = built.drop(index=remainder_gene, errors="ignore")
            # Monocle3 does not filter genes, so missing common rows indicate a
            # malformed/incompatible backend result rather than normal output.
            missing = [gene for gene in sc_common.index if gene not in built.index]
            if missing:
                raise BackendExecutionError(
                    "Monocle3 reference omitted common gene(s): "
                    + ", ".join(missing[:5])
                )
            result = _coerce_signature(
                built.reindex(index=sc_common.index), "Monocle3"
            )

        else:  # music2
            common_count = int(sc_common.shape[0])
            provenance = getattr(prepared, "provenance", {}) or {}
            source_gene_count_value = provenance.get("n_genes_source")
            try:
                source_gene_count = int(source_gene_count_value)
            except (TypeError, ValueError):
                source_gene_count = sum(
                    str(gene) != remainder_gene for gene in prepared.adata.var_names
                )
            if source_gene_count < common_count:
                raise InputValidationError(
                    "Prepared single-cell provenance reports fewer source genes "
                    f"({source_gene_count}) than the observed common-gene count "
                    f"({common_count})."
                )
            comparison_gene_count = min(int(bulk.shape[0]), source_gene_count)
            overlap_fraction = common_count / comparison_gene_count
            if common_count < 0.2 * comparison_gene_count:
                raise InputValidationError(
                    "MuSiC2 has insufficient bulk/single-cell gene overlap "
                    f"({common_count}/{comparison_gene_count}, "
                    f"{overlap_fraction:.1%}); the original DECEPTICONx gate "
                    "requires at least 20%."
                )

            cell_type_count = int(labels.nunique())
            minimum = max(2, cell_type_count)
            if common_count < minimum:
                raise InputValidationError(
                    "MuSiC2 common-gene matrix is not identifiable "
                    f"({common_count} genes for {cell_type_count} cell type(s); "
                    f"at least {minimum} genes required)."
                )

            result = _coerce_signature(
                _music2_mean_reference(
                    sc_common, bulk_common, labels, markers=markers
                ),
                "MuSiC2",
            )
            diagnostics["references"][name] = {
                "genes": int(result.shape[0]),
                "cell_types": int(result.shape[1]),
                # Calling decepticon_fast.music2_base here would repeat the
                # expensive work only to discard its known scalar/order-buggy
                # result.  The corrected implementation is intentionally the
                # sole execution path.
                "upstream_legacy_call_skipped": True,
                "corrected_mean_expression": True,
                "source_gene_count": source_gene_count,
                "bulk_gene_count": int(bulk.shape[0]),
                "overlap_denominator_gene_count": comparison_gene_count,
                "overlap_fraction": overlap_fraction,
                "minimum_identifiable_gene_count": minimum,
                "expression_scale": (
                    "normalized_compatibility"
                    if prepared.normalized_input
                    else "raw_counts"
                ),
            }
            output[name] = result
            continue

        diagnostics["references"][name] = {
            "genes": int(result.shape[0]),
            "cell_types": int(result.shape[1]),
        }
        output[name] = result

    return output, diagnostics


__all__ = ["REFERENCE_BUILDERS", "build_references"]
