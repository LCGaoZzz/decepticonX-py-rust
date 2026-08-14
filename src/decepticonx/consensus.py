"""Deterministic correlation-selected consensus deconvolution."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Literal, Mapping, Sequence

import numpy as np
import pandas as pd

from .exceptions import ConsensusError
from .models import BranchKey


_NEGATIVE_TOLERANCE = 1e-8
_CONSENSUS_MODES = {"r_literal", "corrected"}
ConsensusMode = Literal["r_literal", "corrected"]


@dataclass(frozen=True, slots=True)
class ConsensusResult:
    """Unclosed, closed, and mode-primary consensus tables plus diagnostics."""

    unclosed: pd.DataFrame
    closed: pd.DataFrame
    primary: pd.DataFrame
    diagnostics: dict[str, Any]

    @property
    def raw(self) -> pd.DataFrame:
        """Backward-compatible alias for :attr:`unclosed`."""

        return self.unclosed

    @property
    def normalized(self) -> pd.DataFrame:
        """Backward-compatible alias for :attr:`closed`."""

        return self.closed


def _ordered_keys(
    estimates: Mapping[BranchKey, pd.DataFrame],
    strategy_order: Sequence[BranchKey] | None,
) -> list[BranchKey]:
    keys = list(estimates)
    for key in keys:
        if not isinstance(key, BranchKey):
            raise ConsensusError("Every estimate key must be a BranchKey.")
    slugs = [key.slug for key in keys]
    if len(set(slugs)) != len(slugs):
        raise ConsensusError(
            "Branch method/reference names produce duplicate output slugs; "
            "avoid '__' combinations that make branches ambiguous."
        )

    if strategy_order is None:
        # Mapping insertion order is observable in the R-literal ranking when
        # correlations tie.  Keep it instead of silently imposing slug order.
        return keys

    ordered = list(strategy_order)
    for key in ordered:
        if not isinstance(key, BranchKey):
            raise ConsensusError("Every strategy_order entry must be a BranchKey.")
    if len(set(ordered)) != len(ordered):
        raise ConsensusError("strategy_order must not contain duplicate branches.")
    if len(ordered) != len(keys) or set(ordered) != set(keys):
        missing = sorted(key.slug for key in set(keys).difference(ordered))
        extra = sorted(key.slug for key in set(ordered).difference(keys))
        raise ConsensusError(
            "strategy_order must contain every estimate branch exactly once "
            f"(missing={missing}, extra={extra})."
        )
    return ordered


def _prepare_estimates(
    estimates: Mapping[BranchKey, pd.DataFrame],
    cell_types: list[str],
    strategy_order: Sequence[BranchKey] | None,
) -> tuple[
    list[BranchKey],
    pd.Index,
    dict[BranchKey, np.ndarray],
    dict[BranchKey, np.ndarray],
    dict[str, Any],
]:
    """Validate branches and retain native and compositional target views.

    Correlations are computed from the native target columns.  Aggregation uses
    target columns divided by the sum across *all* native output columns, so an
    auxiliary compartment such as EPIC ``otherCells`` remains in the
    denominator even though it is not itself a consensus target.
    """

    keys = _ordered_keys(estimates, strategy_order)
    if not keys:
        raise ConsensusError("At least one branch estimate is required.")

    first = estimates[keys[0]]
    if not isinstance(first, pd.DataFrame) or first.empty:
        raise ConsensusError("Every branch estimate must be a non-empty DataFrame.")
    if first.index.has_duplicates:
        raise ConsensusError("Branch sample identifiers must be unique.")
    samples = first.index.copy()
    sample_set = set(samples)

    native: dict[BranchKey, np.ndarray] = {}
    aggregation: dict[BranchKey, np.ndarray] = {}
    clipped_by_branch: dict[str, int] = {}
    zero_rows_by_branch: dict[str, list[str]] = {}
    normalization_columns: dict[str, list[str]] = {}

    for key in keys:
        frame = estimates[key]
        if not isinstance(frame, pd.DataFrame) or frame.empty:
            raise ConsensusError(
                f"Branch {key.slug!r} must be a non-empty DataFrame."
            )
        if frame.index.has_duplicates or frame.columns.has_duplicates:
            raise ConsensusError(
                f"Branch {key.slug!r} has duplicate sample or cell-type labels."
            )
        if len(frame.index) != len(samples) or set(frame.index) != sample_set:
            raise ConsensusError(
                f"Branch {key.slug!r} does not contain the same samples "
                "as all branches."
            )

        aligned = frame.reindex(index=samples).copy()
        aligned.columns = pd.Index([str(column) for column in aligned.columns])
        if aligned.columns.has_duplicates:
            raise ConsensusError(
                f"Branch {key.slug!r} has duplicate cell-type labels after "
                "string normalization."
            )
        missing = [
            cell_type for cell_type in cell_types if cell_type not in aligned.columns
        ]
        if missing:
            raise ConsensusError(
                f"Branch {key.slug!r} is missing target cell type(s): "
                + ", ".join(missing)
            )

        try:
            values = aligned.to_numpy(dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ConsensusError(
                f"Branch {key.slug!r} contains non-numeric estimates."
            ) from exc
        if not np.isfinite(values).all():
            raise ConsensusError(
                f"Branch {key.slug!r} contains non-finite estimates."
            )

        material_negative = values < -_NEGATIVE_TOLERANCE
        if material_negative.any():
            minimum = float(values.min())
            raise ConsensusError(
                f"Branch {key.slug!r} contains a materially negative estimate "
                f"({minimum:.6g}); only numerical noise down to "
                f"-{_NEGATIVE_TOLERANCE:g} is clipped."
            )
        negative = values < 0
        clipped_by_branch[key.slug] = int(negative.sum())
        if negative.any():
            values = values.copy()
            values[negative] = 0.0

        target_positions = [
            aligned.columns.get_loc(cell_type) for cell_type in cell_types
        ]
        native[key] = values[:, target_positions].copy()

        row_sums = values.sum(axis=1)
        positive = row_sums > 0
        normalized_all = np.zeros_like(values, dtype=np.float64)
        normalized_all[positive] = values[positive] / row_sums[positive, None]
        aggregation[key] = normalized_all[:, target_positions].copy()
        zero_rows_by_branch[key.slug] = [
            str(sample) for sample in samples[~positive]
        ]
        normalization_columns[key.slug] = list(aligned.columns)

    diagnostics = {
        "clipped_negative_values": clipped_by_branch,
        "zero_sum_rows_by_branch": zero_rows_by_branch,
        "normalization_columns_by_branch": normalization_columns,
    }
    return keys, samples, native, aggregation, diagnostics


def _pearson(left: np.ndarray, right: np.ndarray) -> float | None:
    if left.size < 2:
        return None
    left_centered = left - left.mean()
    right_centered = right - right.mean()
    left_norm = float(np.linalg.norm(left_centered))
    right_norm = float(np.linalg.norm(right_centered))
    if left_norm == 0.0 or right_norm == 0.0:
        return None
    value = float(np.dot(left_centered, right_centered) / (left_norm * right_norm))
    if not np.isfinite(value):
        return None
    return float(np.clip(value, -1.0, 1.0))


def _correlation_matrix(
    keys: Sequence[BranchKey],
    native: Mapping[BranchKey, np.ndarray],
    type_position: int,
) -> np.ndarray:
    size = len(keys)
    matrix = np.full((size, size), np.nan, dtype=np.float64)
    for left_position, left_key in enumerate(keys):
        left = native[left_key][:, type_position]
        for right_position in range(left_position, size):
            right_key = keys[right_position]
            correlation = _pearson(left, native[right_key][:, type_position])
            if correlation is not None:
                matrix[left_position, right_position] = correlation
                matrix[right_position, left_position] = correlation
    return matrix


def _matrix_payload(matrix: np.ndarray) -> list[list[float | None]]:
    return [
        [None if not np.isfinite(value) else float(value) for value in row]
        for row in matrix
    ]


def _pair_labels(
    pairs: Sequence[tuple[int, int]], keys: Sequence[BranchKey]
) -> list[list[str]]:
    return [[keys[left].slug, keys[right].slug] for left, right in pairs]


def _endpoint_weights(
    pairs: Sequence[tuple[int, int]], keys: Sequence[BranchKey]
) -> dict[str, float]:
    if not pairs:
        return {}
    endpoint_weight = 1.0 / (2.0 * len(pairs))
    weights: dict[str, float] = {}
    for left, right in pairs:
        for position in (left, right):
            slug = keys[position].slug
            weights[slug] = weights.get(slug, 0.0) + endpoint_weight
    return weights


def _weighted_values(
    pairs: Sequence[tuple[int, int]],
    keys: Sequence[BranchKey],
    aggregation: Mapping[BranchKey, np.ndarray],
    type_position: int,
) -> np.ndarray:
    weights = _endpoint_weights(pairs, keys)
    sample_count = aggregation[keys[0]].shape[0]
    values = np.zeros(sample_count, dtype=np.float64)
    by_slug = {key.slug: key for key in keys}
    for slug, weight in weights.items():
        values += weight * aggregation[by_slug[slug]][:, type_position]
    return values


def _r_literal_mask(keys: Sequence[BranchKey]) -> np.ndarray:
    size = len(keys)
    references = {key.reference for key in keys}
    mask = np.zeros((size, size), dtype=bool)
    if len(references) > 1:
        for left, left_key in enumerate(keys):
            for right, right_key in enumerate(keys):
                if left_key.reference == right_key.reference:
                    mask[left, right] = True
        return mask

    # Exact special case from DECEPTICON_custom_output.R: for one signature,
    # clear the diagonal and the hard-coded 1-based [4,3]/[3,4] pair only.
    if size < 4:
        raise ConsensusError(
            "r_literal mode with one reference requires at least four branches "
            "to reproduce the R [4,3]/[3,4] mask."
        )
    np.fill_diagonal(mask, True)
    mask[3, 2] = True
    mask[2, 3] = True
    return mask


def _r_literal_pairs(
    correlation: np.ndarray,
    mask: np.ndarray,
) -> tuple[list[tuple[int, int]], list[float]]:
    ranking = np.where(np.isfinite(correlation), correlation, 0.0).copy()
    ranking[mask] = 0.0
    order = np.argsort(-ranking.ravel(order="C"), kind="stable")
    # R: order(...)[1:(2*n)] followed by positions 1, 3 for n == 2.
    selected_flat = order[:4:2]
    if len(selected_flat) != 2:
        raise ConsensusError(
            "r_literal mode requires enough correlation-matrix positions for two pairs."
        )
    size = correlation.shape[0]
    pairs = [divmod(int(position), size) for position in selected_flat]
    values = [float(ranking[left, right]) for left, right in pairs]
    return pairs, values


def _corrected_pairs(
    correlation: np.ndarray,
    keys: Sequence[BranchKey],
    n_pairs: int,
) -> tuple[list[tuple[int, int]], list[float], np.ndarray, int]:
    size = len(keys)
    mask = np.zeros((size, size), dtype=bool)
    candidates: list[tuple[float, int, int]] = []
    for left, left_key in enumerate(keys):
        for right, right_key in enumerate(keys):
            if left == right or left_key.reference == right_key.reference:
                mask[left, right] = True
        for right in range(left + 1, size):
            right_key = keys[right]
            if left_key.reference == right_key.reference:
                continue
            value = correlation[left, right]
            if not np.isfinite(value) or value <= 0.0:
                continue
            candidates.append((float(value), left, right))

    candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
    chosen = candidates[: min(n_pairs, len(candidates))]
    return (
        [(left, right) for _, left, right in chosen],
        [value for value, _, _ in chosen],
        mask,
        len(candidates),
    )


def _reference_balanced_median(
    keys: Sequence[BranchKey],
    aggregation: Mapping[BranchKey, np.ndarray],
    type_position: int,
) -> np.ndarray:
    by_reference: dict[str, list[np.ndarray]] = {}
    for key in keys:
        by_reference.setdefault(key.reference, []).append(
            aggregation[key][:, type_position]
        )
    reference_values = [
        np.median(np.stack(values, axis=0), axis=0)
        for values in by_reference.values()
    ]
    return np.median(np.stack(reference_values, axis=0), axis=0)


def _diagnostic_pairs(
    matrix: np.ndarray,
    keys: Sequence[BranchKey],
    *,
    predicate: Callable[[int, int, float], bool],
) -> list[list[str]]:
    pairs: list[tuple[int, int]] = []
    for left in range(len(keys)):
        for right in range(left, len(keys)):
            if predicate(left, right, matrix[left, right]):
                pairs.append((left, right))
    return _pair_labels(pairs, keys)


def build_consensus(
    estimates: Mapping[BranchKey, pd.DataFrame],
    *,
    cell_types: Sequence[str],
    n_pairs: int = 2,
    mode: ConsensusMode = "r_literal",
    strategy_order: Sequence[BranchKey] | None = None,
) -> ConsensusResult:
    """Build an R-literal or corrected DECEPTICON consensus.

    ``r_literal`` reproduces the R correlation masking, row-major positional
    ranking, repeated endpoint weights, and unclosed primary output.  It uses a
    uniform all-column compositional view for aggregation because the upstream
    custom-output normalization loops are not defined for the 15-branch setup.

    ``corrected`` structurally excludes self/same-reference pairs, selects up
    to ``n_pairs`` unique positive finite pairs, and closes the primary result.
    If no eligible pair exists, it uses a reference-balanced median.
    """

    if isinstance(n_pairs, bool) or not isinstance(n_pairs, (int, np.integer)):
        raise ConsensusError("n_pairs must be a positive integer.")
    if int(n_pairs) <= 0:
        raise ConsensusError("n_pairs must be a positive integer.")
    n_pairs = int(n_pairs)

    mode = str(mode)
    if mode not in _CONSENSUS_MODES:
        raise ConsensusError(
            "mode must be either 'r_literal' or 'corrected'."
        )
    if mode == "r_literal" and n_pairs != 2:
        raise ConsensusError("r_literal mode requires n_pairs=2.")

    ordered_types = [str(value) for value in cell_types]
    if not ordered_types:
        raise ConsensusError("At least one signature cell type is required.")
    if len(set(ordered_types)) != len(ordered_types):
        raise ConsensusError("Signature cell types must be unique.")

    keys, samples, native, aggregation, prep_diagnostics = _prepare_estimates(
        estimates, ordered_types, strategy_order
    )
    unclosed_values = np.empty(
        (len(samples), len(ordered_types)), dtype=np.float64
    )

    correlations: dict[str, list[list[float | None]]] = {}
    eligible_counts: dict[str, int] = {}
    selected_pairs: dict[str, list[list[str]]] = {}
    selected_pair_correlations: dict[str, list[float]] = {}
    endpoint_weights: dict[str, dict[str, float]] = {}
    masked_pairs: dict[str, list[list[str]]] = {}
    undefined_pairs: dict[str, list[list[str]]] = {}
    fallbacks: dict[str, str] = {}

    literal_mask = _r_literal_mask(keys) if mode == "r_literal" else None

    for type_position, cell_type in enumerate(ordered_types):
        correlation = _correlation_matrix(keys, native, type_position)
        correlations[cell_type] = _matrix_payload(correlation)
        undefined_pairs[cell_type] = _diagnostic_pairs(
            correlation,
            keys,
            predicate=lambda _left, _right, value: not np.isfinite(value),
        )

        if mode == "r_literal":
            assert literal_mask is not None
            pairs, pair_values = _r_literal_pairs(correlation, literal_mask)
            mask = literal_mask
            eligible_counts[cell_type] = sum(
                1
                for left in range(len(keys))
                for right in range(left + 1, len(keys))
                if not mask[left, right] and np.isfinite(correlation[left, right])
            )
        else:
            pairs, pair_values, mask, eligible_count = _corrected_pairs(
                correlation, keys, n_pairs
            )
            eligible_counts[cell_type] = eligible_count

        masked_pairs[cell_type] = _diagnostic_pairs(
            mask,
            keys,
            predicate=lambda left, right, _value: bool(mask[left, right]),
        )
        selected_pairs[cell_type] = _pair_labels(pairs, keys)
        selected_pair_correlations[cell_type] = pair_values
        endpoint_weights[cell_type] = _endpoint_weights(pairs, keys)

        if mode == "corrected" and not pairs:
            unclosed_values[:, type_position] = _reference_balanced_median(
                keys, aggregation, type_position
            )
            fallbacks[cell_type] = "no_positive_finite_cross_reference_pairs"
        else:
            unclosed_values[:, type_position] = _weighted_values(
                pairs, keys, aggregation, type_position
            )

    unclosed = pd.DataFrame(
        unclosed_values, index=samples, columns=ordered_types
    )
    row_sums = unclosed_values.sum(axis=1)
    positive = row_sums > 0
    closed_values = np.zeros_like(unclosed_values)
    closed_values[positive] = unclosed_values[positive] / row_sums[positive, None]
    closed = pd.DataFrame(closed_values, index=samples, columns=ordered_types)
    closed_row_sums = closed_values.sum(axis=1)

    primary = unclosed if mode == "r_literal" else closed
    diagnostics: dict[str, Any] = {
        "mode": mode,
        "strategy_order": [key.slug for key in keys],
        "selection_space": "native_target_estimates",
        "aggregation_space": "row_normalized_all_native_columns",
        "primary_space": "unclosed" if mode == "r_literal" else "closed",
        "n_pairs_requested": n_pairs,
        **prep_diagnostics,
        "correlations": correlations,
        "eligible_pair_counts": eligible_counts,
        "selected_pairs": selected_pairs,
        "selected_pair_correlations": selected_pair_correlations,
        "endpoint_weights": endpoint_weights,
        "masked_pairs": masked_pairs,
        "undefined_pairs": undefined_pairs,
        "fallbacks": fallbacks,
        "zero_sum_unclosed_rows": [str(sample) for sample in samples[~positive]],
        "zero_sum_closed_rows": [
            str(sample) for sample in samples[closed_row_sums <= 0]
        ],
    }
    return ConsensusResult(
        unclosed=unclosed,
        closed=closed,
        primary=primary,
        diagnostics=diagnostics,
    )


__all__ = ["ConsensusMode", "ConsensusResult", "build_consensus"]
