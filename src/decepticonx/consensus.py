"""Deterministic correlation-selected consensus deconvolution."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from .exceptions import ConsensusError
from .models import BranchKey


_NEGATIVE_TOLERANCE = 1e-8


@dataclass(frozen=True, slots=True)
class ConsensusResult:
    """Raw cell-type consensus, row-normalised output, and diagnostics."""

    raw: pd.DataFrame
    normalized: pd.DataFrame
    diagnostics: dict[str, Any]


def _ordered_keys(estimates: Mapping[BranchKey, pd.DataFrame]) -> list[BranchKey]:
    keys = list(estimates)
    for key in keys:
        if not isinstance(key, BranchKey):
            raise ConsensusError("Every estimate key must be a BranchKey.")
    return sorted(keys, key=lambda key: (key.method, key.reference))


def _prepare_estimates(
    estimates: Mapping[BranchKey, pd.DataFrame], cell_types: list[str]
) -> tuple[list[BranchKey], pd.Index, dict[BranchKey, np.ndarray], dict[str, Any]]:
    keys = _ordered_keys(estimates)
    if not keys:
        raise ConsensusError("At least one branch estimate is required.")

    first = estimates[keys[0]]
    if not isinstance(first, pd.DataFrame) or first.empty:
        raise ConsensusError("Every branch estimate must be a non-empty DataFrame.")
    if first.index.has_duplicates:
        raise ConsensusError("Branch sample identifiers must be unique.")
    samples = first.index.copy()
    sample_set = set(samples)

    prepared: dict[BranchKey, np.ndarray] = {}
    clipped_by_branch: dict[str, int] = {}
    zero_rows_by_branch: dict[str, list[str]] = {}
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
                f"Branch {key.slug!r} does not contain the same samples as all branches."
            )
        # Missing target types are genuine zero estimates for methods that do
        # not model every signature column.  Non-signature extra columns are
        # intentionally discarded.
        aligned = frame.reindex(index=samples, columns=cell_types, fill_value=0.0)
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

        row_sums = values.sum(axis=1)
        positive = row_sums > 0
        normalised = np.zeros_like(values, dtype=np.float64)
        normalised[positive] = values[positive] / row_sums[positive, None]
        zero_rows_by_branch[key.slug] = [
            str(sample) for sample in samples[~positive]
        ]
        prepared[key] = normalised

    diagnostics = {
        "clipped_negative_values": clipped_by_branch,
        "zero_sum_rows_by_branch": zero_rows_by_branch,
    }
    return keys, samples, prepared, diagnostics


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
    # Floating-point roundoff can put a mathematically exact correlation a
    # handful of ulps outside [-1, 1].
    return float(np.clip(value, -1.0, 1.0))


def build_consensus(
    estimates: Mapping[BranchKey, pd.DataFrame],
    *,
    cell_types: Sequence[str],
    n_pairs: int = 2,
) -> ConsensusResult:
    """Build a DECEPTICON-style consensus from samples-by-types branches.

    Each branch is first clipped for tiny negative numerical noise, restricted
    to the signature's cell types, and row-normalised.  For each cell type we
    rank Pearson correlations between branches, excluding the diagonal and
    branches made from the same reference.  The top ``n_pairs`` *unique
    unordered* pairs are selected.  A pair contributes the mean of its two
    endpoints, matching the intent of the original DECEPTICON R code.

    Constant vectors or too few eligible finite pairs trigger a deterministic
    per-sample median across all branches for that cell type.
    """

    if isinstance(n_pairs, bool) or not isinstance(n_pairs, (int, np.integer)):
        raise ConsensusError("n_pairs must be a positive integer.")
    if int(n_pairs) <= 0:
        raise ConsensusError("n_pairs must be a positive integer.")
    n_pairs = int(n_pairs)

    ordered_types = [str(value) for value in cell_types]
    if not ordered_types:
        raise ConsensusError("At least one signature cell type is required.")
    if len(set(ordered_types)) != len(ordered_types):
        raise ConsensusError("Signature cell types must be unique.")

    keys, samples, branch_values, prep_diagnostics = _prepare_estimates(
        estimates, ordered_types
    )
    raw_values = np.empty((len(samples), len(ordered_types)), dtype=np.float64)
    selected_pairs: dict[str, list[list[str]]] = {}
    pair_correlations: dict[str, list[float]] = {}
    fallbacks: dict[str, str] = {}
    eligible_counts: dict[str, int] = {}

    for type_position, cell_type in enumerate(ordered_types):
        candidates: list[tuple[float, str, str, int, int]] = []
        for left_position, left_key in enumerate(keys):
            left = branch_values[left_key][:, type_position]
            for right_position in range(left_position + 1, len(keys)):
                right_key = keys[right_position]
                # Mask both the diagonal (already excluded by i < j) and
                # branches sharing one reference template.
                if left_key.reference == right_key.reference:
                    continue
                correlation = _pearson(
                    left, branch_values[right_key][:, type_position]
                )
                if correlation is None:
                    continue
                candidates.append(
                    (
                        correlation,
                        left_key.slug,
                        right_key.slug,
                        left_position,
                        right_position,
                    )
                )

        candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
        eligible_counts[cell_type] = len(candidates)
        if len(candidates) < n_pairs:
            stacked = np.stack(
                [branch_values[key][:, type_position] for key in keys], axis=0
            )
            raw_values[:, type_position] = np.median(stacked, axis=0)
            selected_pairs[cell_type] = []
            pair_correlations[cell_type] = []
            reason = "constant_or_undefined_correlation" if not candidates else "insufficient_pairs"
            fallbacks[cell_type] = reason
            continue

        chosen = candidates[:n_pairs]
        pair_means = []
        selected_pairs[cell_type] = []
        pair_correlations[cell_type] = []
        for correlation, left_slug, right_slug, left_position, right_position in chosen:
            left_key = keys[left_position]
            right_key = keys[right_position]
            pair_means.append(
                (
                    branch_values[left_key][:, type_position]
                    + branch_values[right_key][:, type_position]
                )
                / 2.0
            )
            selected_pairs[cell_type].append([left_slug, right_slug])
            pair_correlations[cell_type].append(correlation)
        raw_values[:, type_position] = np.mean(np.stack(pair_means, axis=0), axis=0)

    raw = pd.DataFrame(raw_values, index=samples, columns=ordered_types)
    row_sums = raw_values.sum(axis=1)
    positive = row_sums > 0
    normalized_values = np.zeros_like(raw_values)
    normalized_values[positive] = raw_values[positive] / row_sums[positive, None]
    normalized = pd.DataFrame(
        normalized_values, index=samples, columns=ordered_types
    )

    diagnostics: dict[str, Any] = {
        **prep_diagnostics,
        "branch_order": [key.slug for key in keys],
        "n_pairs_requested": n_pairs,
        "eligible_pair_counts": eligible_counts,
        "selected_pairs": selected_pairs,
        "selected_pair_correlations": pair_correlations,
        "fallbacks": fallbacks,
        "zero_sum_consensus_rows": [str(sample) for sample in samples[~positive]],
    }
    return ConsensusResult(raw=raw, normalized=normalized, diagnostics=diagnostics)


__all__ = ["ConsensusResult", "build_consensus"]
