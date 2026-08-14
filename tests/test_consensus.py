from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from decepticonx import consensus as consensus_module
from decepticonx.consensus import build_consensus
from decepticonx.exceptions import ConsensusError
from decepticonx.models import BranchKey


SAMPLES = ["s1", "s2", "s3", "s4"]
TYPES = ["A", "B"]
ORACLE = Path(__file__).parent / "oracles" / "r_selector_expected.tsv"


def _frame(a_values, *, index=SAMPLES, total: float = 1.0) -> pd.DataFrame:
    a = np.asarray(a_values, dtype=float)
    return pd.DataFrame({"A": a, "B": total - a}, index=index)


def _selected(result, cell_type: str) -> list[list[str]]:
    return result.diagnostics["selected_pairs"][cell_type]


def _assert_result_contract(result, *, primary: str) -> None:
    """Both names remain available while the primary output is mode-specific."""

    pd.testing.assert_frame_equal(result.raw, result.unclosed)
    pd.testing.assert_frame_equal(result.normalized, result.closed)
    pd.testing.assert_frame_equal(result.primary, getattr(result, primary))
    assert result.diagnostics["selection_space"] == "native_target_estimates"
    assert (
        result.diagnostics["aggregation_space"]
        == "row_normalized_all_native_columns"
    )
    assert result.diagnostics["primary_space"] == primary


def test_corrected_masks_same_reference_and_averages_pair_endpoints() -> None:
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

    result = build_consensus(
        estimates, cell_types=TYPES, n_pairs=1, mode="corrected"
    )

    assert _selected(result, "A") == [[first.slug, second.slug]]
    assert [first.slug, same_reference.slug] not in _selected(result, "A")
    expected_a = (estimates[first]["A"] + estimates[second]["A"]) / 2.0
    np.testing.assert_allclose(result.unclosed["A"], expected_a)
    np.testing.assert_allclose(result.closed.sum(axis=1), 1.0)
    _assert_result_contract(result, primary="closed")


def test_corrected_is_deterministic_under_mapping_order_and_ties() -> None:
    branches = [
        BranchKey("alpha", "r1"),
        BranchKey("beta", "r2"),
        BranchKey("gamma", "r3"),
    ]
    values = [_frame([0.1, 0.3, 0.7, 0.9])] * 3
    forward = OrderedDict(zip(branches, values, strict=True))
    reverse = OrderedDict(reversed(list(forward.items())))

    first = build_consensus(
        forward,
        cell_types=TYPES,
        n_pairs=2,
        mode="corrected",
        strategy_order=branches,
    )
    second = build_consensus(
        reverse,
        cell_types=TYPES,
        n_pairs=2,
        mode="corrected",
        strategy_order=branches,
    )

    pd.testing.assert_frame_equal(first.unclosed, second.unclosed)
    pd.testing.assert_frame_equal(first.closed, second.closed)
    assert first.diagnostics["selected_pairs"] == second.diagnostics["selected_pairs"]
    assert _selected(first, "A") == [
        [branches[0].slug, branches[1].slug],
        [branches[0].slug, branches[2].slug],
    ]


def test_corrected_uses_the_available_positive_pair_instead_of_falling_back() -> None:
    first = BranchKey("m1", "r1")
    second = BranchKey("m2", "r2")
    constant = BranchKey("m3", "r3")
    estimates = {
        first: _frame([0.1, 0.2, 0.7, 0.9]),
        second: _frame([0.1, 0.25, 0.65, 0.85]),
        constant: _frame([0.5, 0.5, 0.5, 0.5]),
    }

    result = build_consensus(
        estimates, cell_types=TYPES, n_pairs=2, mode="corrected"
    )

    assert _selected(result, "A") == [[first.slug, second.slug]]
    assert "A" not in result.diagnostics["fallbacks"]
    expected = (estimates[first]["A"] + estimates[second]["A"]) / 2.0
    np.testing.assert_allclose(result.unclosed["A"], expected)


def test_corrected_no_pair_fallback_is_balanced_by_reference() -> None:
    # With no defined positive correlation, a branch-wise median would be
    # 0.25.  Median-within-reference followed by median-across-references is
    # 0.55, so this fixture detects accidental over-weighting of ref1.
    estimates = {
        BranchKey("m1", "ref1"): _frame([0.1] * 4),
        BranchKey("m2", "ref1"): _frame([0.2] * 4),
        BranchKey("m3", "ref1"): _frame([0.3] * 4),
        BranchKey("m4", "ref2"): _frame([0.9] * 4),
    }

    result = build_consensus(
        estimates, cell_types=TYPES, n_pairs=2, mode="corrected"
    )

    np.testing.assert_allclose(result.unclosed["A"], 0.55)
    np.testing.assert_allclose(result.unclosed["B"], 0.45)
    assert result.diagnostics["selected_pairs"] == {"A": [], "B": []}
    assert set(result.diagnostics["fallbacks"]) == set(TYPES)


def test_r_literal_duplicate_endpoint_accumulates_positional_weight() -> None:
    samples = [f"s{i}" for i in range(6)]
    signal = np.array([0.1, 0.3, 0.6, 0.8, 0.2, 0.7])
    noisy = np.array([0.12, 0.25, 0.67, 0.72, 0.28, 0.64])
    branches = [
        BranchKey("m0", "r0"),
        BranchKey("m1", "r1"),
        BranchKey("m2", "r2"),
    ]
    estimates = {
        branches[0]: _frame(signal, index=samples),
        branches[1]: _frame(signal, index=samples),
        branches[2]: _frame(noisy, index=samples),
    }
    # A backend-only compartment participates in closure even though it is
    # not one of the requested consensus columns.
    estimates[branches[2]]["otherCells"] = 0.5

    result = build_consensus(
        estimates,
        cell_types=TYPES,
        mode="r_literal",
        strategy_order=branches,
    )

    assert _selected(result, "A") == [
        [branches[0].slug, branches[1].slug],
        [branches[0].slug, branches[2].slug],
    ]
    np.testing.assert_allclose(
        result.unclosed["A"], 0.75 * signal + 0.25 * (noisy / 1.5)
    )
    assert result.diagnostics["normalization_columns_by_branch"][branches[2].slug] == [
        "A",
        "B",
        "otherCells",
    ]
    _assert_result_contract(result, primary="unclosed")


def test_r_literal_selector_matches_the_static_base_r_top_pair_oracle() -> None:
    matrix = np.array(
        [
            [0.0, 0.9, 0.1, 0.8],
            [0.9, 0.0, 0.7, 0.2],
            [0.1, 0.7, 0.0, 0.95],
            [0.8, 0.2, 0.95, 0.0],
        ]
    )
    pairs, values = consensus_module._r_literal_pairs(
        matrix, np.zeros_like(matrix, dtype=bool)
    )
    oracle = pd.read_csv(ORACLE, sep="\t")
    expected = oracle.loc[oracle["case"] == "top_pairs"]
    expected_pairs = [
        (int(row.left) - 1, int(row.right) - 1)
        for row in expected.itertuples(index=False)
    ]

    assert pairs == expected_pairs
    np.testing.assert_allclose(values, [0.95, 0.9], rtol=0.0, atol=0.0)


def test_r_literal_matches_masked_zero_oracle_while_corrected_falls_back() -> None:
    oracle = pd.read_csv(ORACLE, sep="\t")
    expected = oracle.loc[oracle["case"] == "masked_negative"]
    branches = [
        BranchKey("m0", "r0"),
        BranchKey("m1", "r1"),
        BranchKey("m2", "r2"),
    ]
    estimates = {
        branches[0]: _frame([2, 1, 0], index=["s1", "s2", "s3"], total=3),
        branches[1]: _frame([0, 2, 1], index=["s1", "s2", "s3"], total=3),
        branches[2]: _frame([1, 0, 2], index=["s1", "s2", "s3"], total=3),
    }

    literal = build_consensus(
        estimates,
        cell_types=TYPES,
        mode="r_literal",
        strategy_order=branches,
    )
    expected_pairs = [
        [branches[int(row.left) - 1].slug, branches[int(row.right) - 1].slug]
        for row in expected.itertuples(index=False)
    ]
    assert _selected(literal, "A") == expected_pairs
    np.testing.assert_allclose(literal.unclosed["A"], [0.5, 1 / 6, 1 / 3])
    np.testing.assert_allclose(literal.unclosed["B"], [0.5, 5 / 6, 2 / 3])

    corrected = build_consensus(
        estimates,
        cell_types=TYPES,
        n_pairs=2,
        mode="corrected",
        strategy_order=branches,
    )
    assert _selected(corrected, "A") == []
    assert _selected(corrected, "B") == []
    np.testing.assert_allclose(corrected.unclosed["A"], 1 / 3)
    np.testing.assert_allclose(corrected.unclosed["B"], 2 / 3)


def test_selection_uses_raw_values_but_aggregation_uses_row_proportions() -> None:
    samples = [f"s{i}" for i in range(6)]
    # Raw top pairs differ from those obtained after per-row closure.  This
    # catches the CIBERSORT-relative/CIBERSORT-ABS collapse at selection time.
    values = np.array(
        [
            [[13, 8], [18, 18], [10, 3], [13, 10], [2, 15], [12, 16]],
            [[3, 14], [16, 13], [20, 9], [8, 19], [5, 8], [8, 11]],
            [[11, 15], [4, 19], [13, 2], [5, 20], [12, 18], [19, 20]],
            [[2, 6], [16, 19], [10, 7], [20, 6], [20, 17], [1, 20]],
        ],
        dtype=float,
    )
    branches = [BranchKey(f"m{i}", f"r{i}") for i in range(4)]
    estimates = {
        key: pd.DataFrame(value, index=samples, columns=TYPES)
        for key, value in zip(branches, values, strict=True)
    }

    result = build_consensus(
        estimates,
        cell_types=TYPES,
        n_pairs=2,
        mode="corrected",
        strategy_order=branches,
    )

    assert _selected(result, "A") == [
        [branches[0].slug, branches[1].slug],
        [branches[1].slug, branches[3].slug],
    ]
    assert _selected(result, "B") == [
        [branches[0].slug, branches[3].slug],
        [branches[0].slug, branches[2].slug],
    ]
    normalized = values / values.sum(axis=2, keepdims=True)
    expected_a = 0.25 * normalized[0, :, 0] + 0.5 * normalized[1, :, 0]
    expected_a += 0.25 * normalized[3, :, 0]
    expected_b = 0.5 * normalized[0, :, 1] + 0.25 * normalized[2, :, 1]
    expected_b += 0.25 * normalized[3, :, 1]
    np.testing.assert_allclose(result.unclosed["A"], expected_a, atol=1e-14)
    np.testing.assert_allclose(result.unclosed["B"], expected_b, atol=1e-14)
    assert not np.allclose(result.unclosed.sum(axis=1), 1.0)
    np.testing.assert_allclose(result.closed.sum(axis=1), 1.0, atol=1e-14)


def test_cibersort_absolute_scale_is_not_erased_before_selection() -> None:
    samples = [f"s{i}" for i in range(6)]
    relative = _frame(
        [0.15, 0.25, 0.4, 0.65, 0.75, 0.55], index=samples
    )
    depth = pd.Series([0.25, 0.5, 1.0, 2.0, 4.0, 8.0], index=samples)
    absolute = relative.mul(depth, axis=0)

    assert not np.allclose(relative, absolute)
    np.testing.assert_allclose(
        absolute.div(absolute.sum(axis=1), axis=0), relative, atol=1e-15
    )


def test_missing_target_cell_type_is_rejected_instead_of_silently_dropped() -> None:
    estimates = {
        BranchKey("m1", "r1"): _frame([0.1, 0.2, 0.3, 0.4]),
        BranchKey("m2", "r2"): pd.DataFrame(
            {"A": [0.2, 0.3, 0.4, 0.5]}, index=SAMPLES
        ),
    }

    with pytest.raises(ConsensusError, match=r"m2__r2.*missing.*B"):
        build_consensus(estimates, cell_types=TYPES, mode="corrected")


def test_ambiguous_branch_slugs_are_rejected() -> None:
    estimates = {
        BranchKey("a__b", "c"): _frame([0.1, 0.2, 0.3, 0.4]),
        BranchKey("a", "b__c"): _frame([0.2, 0.3, 0.4, 0.5]),
    }
    with pytest.raises(ConsensusError, match="duplicate output slugs"):
        build_consensus(estimates, cell_types=TYPES, mode="corrected")


def test_r_literal_requires_exactly_two_pairs() -> None:
    estimates = {
        BranchKey("m1", "r1"): _frame([0.1, 0.2, 0.7, 0.9]),
        BranchKey("m2", "r2"): _frame([0.2, 0.3, 0.6, 0.8]),
    }

    with pytest.raises(ConsensusError, match=r"r_literal.*n_pairs.*2"):
        build_consensus(
            estimates,
            cell_types=TYPES,
            n_pairs=1,
            mode="r_literal",
        )


def test_tiny_negative_noise_is_clipped_before_both_consensus_scales() -> None:
    branches = {
        BranchKey("m1", "r1"): pd.DataFrame(
            {"A": [-1e-12, 0.25, 0.5, 0.75], "B": [1.0, 0.75, 0.5, 0.25]},
            index=SAMPLES,
        ),
        BranchKey("m2", "r2"): _frame([0.0, 0.2, 0.5, 0.8]),
    }
    result = build_consensus(
        branches, cell_types=TYPES, n_pairs=1, mode="corrected"
    )
    assert result.diagnostics["clipped_negative_values"]["m1__r1"] == 1
    assert (result.unclosed.to_numpy() >= 0).all()
    assert (result.closed.to_numpy() >= 0).all()
