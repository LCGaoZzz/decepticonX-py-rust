from __future__ import annotations

from collections import OrderedDict

import numpy as np
import pandas as pd

from decepticonx.consensus import build_consensus
from decepticonx.models import BranchKey


SAMPLES = ["s1", "s2", "s3", "s4"]
TYPES = ["A", "B"]


def _frame(a_values) -> pd.DataFrame:
    a = np.asarray(a_values, dtype=float)
    return pd.DataFrame({"A": a, "B": 1.0 - a}, index=SAMPLES)


def test_same_reference_pairs_are_masked_and_pair_endpoints_are_averaged() -> None:
    first = BranchKey("m1", "ref1")
    same_reference = BranchKey("m2", "ref1")
    second = BranchKey("m3", "ref2")
    third = BranchKey("m4", "ref3")
    estimates = {
        first: _frame([0.1, 0.2, 0.7, 0.9]),
        # Perfectly correlated with first, but forbidden because it uses ref1.
        same_reference: _frame([0.2, 0.3, 0.8, 1.0]),
        second: _frame([0.1, 0.2, 0.6, 0.8]),
        third: _frame([0.9, 0.7, 0.3, 0.1]),
    }

    result = build_consensus(estimates, cell_types=TYPES, n_pairs=1)

    selected = result.diagnostics["selected_pairs"]["A"]
    assert selected == [[first.slug, second.slug]]
    assert [first.slug, same_reference.slug] not in selected
    expected_a = (estimates[first]["A"] + estimates[second]["A"]) / 2.0
    np.testing.assert_allclose(result.raw["A"], expected_a)
    np.testing.assert_allclose(result.normalized.sum(axis=1), 1.0)


def test_consensus_is_deterministic_under_mapping_order_and_ties() -> None:
    branches = [
        BranchKey("alpha", "r1"),
        BranchKey("beta", "r2"),
        BranchKey("gamma", "r3"),
    ]
    values = [_frame([0.1, 0.3, 0.7, 0.9])] * 3
    forward = OrderedDict(zip(branches, values, strict=True))
    reverse = OrderedDict(reversed(list(forward.items())))

    first = build_consensus(forward, cell_types=TYPES, n_pairs=2)
    second = build_consensus(reverse, cell_types=TYPES, n_pairs=2)

    pd.testing.assert_frame_equal(first.raw, second.raw)
    pd.testing.assert_frame_equal(first.normalized, second.normalized)
    assert first.diagnostics["selected_pairs"] == second.diagnostics["selected_pairs"]
    assert first.diagnostics["selected_pairs"]["A"] == [
        [branches[0].slug, branches[1].slug],
        [branches[0].slug, branches[2].slug],
    ]


def test_constant_or_insufficient_correlations_fall_back_to_branch_median() -> None:
    branches = [
        BranchKey("m1", "r1"),
        BranchKey("m2", "r2"),
        BranchKey("m3", "r3"),
    ]
    estimates = {
        branches[0]: _frame([0.2, 0.2, 0.2, 0.2]),
        branches[1]: _frame([0.4, 0.4, 0.4, 0.4]),
        branches[2]: _frame([0.8, 0.8, 0.8, 0.8]),
    }

    result = build_consensus(estimates, cell_types=TYPES, n_pairs=2)

    np.testing.assert_allclose(result.raw["A"], 0.4)
    np.testing.assert_allclose(result.raw["B"], 0.6)
    assert result.diagnostics["selected_pairs"] == {"A": [], "B": []}
    assert set(result.diagnostics["fallbacks"]) == set(TYPES)


def test_tiny_negative_noise_is_clipped_before_row_normalisation() -> None:
    branches = {
        BranchKey("m1", "r1"): pd.DataFrame(
            {"A": [-1e-12, 0.25, 0.5, 0.75], "B": [1.0, 0.75, 0.5, 0.25]},
            index=SAMPLES,
        ),
        BranchKey("m2", "r2"): _frame([0.0, 0.2, 0.5, 0.8]),
    }
    result = build_consensus(branches, cell_types=TYPES, n_pairs=1)
    assert result.diagnostics["clipped_negative_values"]["m1__r1"] == 1
    assert (result.raw.to_numpy() >= 0).all()
