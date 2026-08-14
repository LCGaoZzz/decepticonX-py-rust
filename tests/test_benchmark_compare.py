from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

from decepticonx.models import BranchKey


SCRIPT = (
    Path(__file__).parents[1]
    / "benchmarks"
    / "original_comparison"
    / "compare_results.py"
)
SPEC = importlib.util.spec_from_file_location("benchmark_compare_results", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)


def test_align_samples_accepts_only_complete_unique_numeric_suffix_mapping() -> None:
    frame = pd.DataFrame({"A": [10.0, 20.0]}, index=["1", "2"])
    target = pd.Index(["value...1", "value...2"])

    aligned = COMPARE.align_samples(frame, target)

    assert list(aligned.index) == list(target)
    assert aligned["A"].tolist() == [10.0, 20.0]


def test_align_estimate_preserves_only_epic_other_cells_auxiliary() -> None:
    frame = pd.DataFrame(
        {
            "B cell": [0.3, 0.4],
            "T cell": [0.5, 0.4],
            "otherCells": [0.2, 0.2],
            "diagnostic": [99.0, 99.0],
        },
        index=["s2", "s1"],
    )
    samples = pd.Index(["s1", "s2"])

    epic = COMPARE.align_estimate(
        frame, ["B cell", "T cell"], samples, method="epic"
    )
    other = COMPARE.align_estimate(
        frame, ["B cell", "T cell"], samples, method="music"
    )

    assert list(epic.columns) == ["B cell", "T cell", "otherCells"]
    assert list(other.columns) == ["B cell", "T cell"]
    assert epic.loc["s1", "otherCells"] == 0.2


@pytest.mark.parametrize(
    ("source", "target"),
    [
        (["1", "2"], ["value...1", "value...3"]),
        (["1", "2"], ["sample1", "other1"]),
        (["sampleA", "sampleB"], ["value...1", "value...2"]),
    ],
)
def test_align_samples_rejects_incomplete_or_ambiguous_suffixes(
    source: list[str], target: list[str]
) -> None:
    frame = pd.DataFrame({"A": range(len(source))}, index=source)

    with pytest.raises(ValueError):
        COMPARE.align_samples(frame, pd.Index(target))


def test_frozen_provenance_requires_hash_seed_and_expected_solver(
    tmp_path: Path,
) -> None:
    original = tmp_path / "original"
    frozen = tmp_path / "frozen"
    reference_root = original / "custom_signature_matrix"
    reference_root.mkdir(parents=True)
    frozen.mkdir()
    paths = {}
    for name, filename in COMPARE.REFERENCE_FILES.items():
        path = reference_root / filename
        path.write_text(f"{name}\n", encoding="utf-8")
        paths[name] = path

    bulk_hash = "a" * 64
    mapping_hash = "b" * 64
    pd.DataFrame(
        {
            "key": [
                "seed",
                "bulk_sha256",
                "epic_mrna_cell_sha256",
                "epic_mrna_cell_entries",
            ],
            "value": ["7", bulk_hash, mapping_hash, "12"],
        }
    ).to_csv(original / "manifest.tsv", sep="\t", index=False)
    run = {
        "provenance": {
            "seed": 7,
            "bulk_sha256": bulk_hash,
            "epic_solver_requested": "nm",
            "epic_mrna_cell": {
                "mapping_sha256": mapping_hash,
                "entry_count": 12,
                "values_redacted": True,
            },
            "reference_sha256": {
                name: COMPARE.sha256(path) for name, path in paths.items()
            },
        }
    }
    run_path = frozen / "run.json"
    run_path.write_text(json.dumps(run), encoding="utf-8")

    COMPARE.validate_frozen_provenance(original, frozen, paths)

    run["provenance"]["epic_solver_requested"] = "auto"
    run_path.write_text(json.dumps(run), encoding="utf-8")
    COMPARE.validate_frozen_provenance(
        original, frozen, paths, expected_epic_solver="auto"
    )

    run["provenance"]["epic_solver_requested"] = "qp"
    run_path.write_text(json.dumps(run), encoding="utf-8")
    with pytest.raises(ValueError, match="epic_solver='nm'"):
        COMPARE.validate_frozen_provenance(original, frozen, paths)

    run["provenance"]["epic_solver_requested"] = "nm"
    run["provenance"]["epic_mrna_cell"]["mapping_sha256"] = "c" * 64
    run_path.write_text(json.dumps(run), encoding="utf-8")
    with pytest.raises(ValueError, match="EPIC mapping SHA-256 mismatch"):
        COMPARE.validate_frozen_provenance(original, frozen, paths)


def test_mapping_sha256_uses_cross_runtime_fixed_float_format() -> None:
    module_path = (
        Path(__file__).parents[1]
        / "benchmarks"
        / "original_comparison"
        / "run_accelerated_from_refs.py"
    )
    spec = importlib.util.spec_from_file_location(
        "benchmark_run_accelerated_from_refs", module_path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module.mapping_sha256({"z": 1.0, "a": 1.2345678901234567}) == (
        "63f2337a44f06dc7283207363f139ea06739eeb20b54712c1eab3383a3449fdb"
    )


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ([[0.8, -0.1]], "non-negative"),
        ([[0.8, 0.1]], "sum to one"),
    ],
)
def test_truth_frame_rejects_invalid_proportions(
    tmp_path: Path, values: list[list[float]], message: str
) -> None:
    path = tmp_path / "truth.tsv"
    pd.DataFrame(values, index=["sample"], columns=["A", "B"]).to_csv(
        path, sep="\t"
    )

    with pytest.raises(ValueError, match=message):
        COMPARE.truth_frame(path)


def test_stored_consensus_validation_uses_explicit_three_view_contract(
    tmp_path: Path,
) -> None:
    samples = pd.Index(["s1", "s2", "s3", "s4"])
    target_types = ["A", "B"]
    estimates: dict[BranchKey, pd.DataFrame] = {}
    base = pd.Series([0.1, 0.25, 0.7, 0.9], index=samples)
    for method_index, method in enumerate(COMPARE.R_FILENAME_METHOD_ORDER):
        for reference_index, reference in enumerate(COMPARE.REFERENCE_ORDER):
            offset = 0.002 * (method_index + reference_index)
            a = (base + offset).clip(upper=0.99)
            estimates[BranchKey(method, reference)] = pd.DataFrame(
                {"A": a, "B": 1.0 - a}, index=samples
            )

    expected = COMPARE.benchmark_consensus(estimates, target_types)
    expected.primary.to_csv(tmp_path / "consensus.tsv", sep="\t")
    expected.unclosed.to_csv(tmp_path / "consensus_unclosed.tsv", sep="\t")
    expected.closed.to_csv(tmp_path / "consensus_closed.tsv", sep="\t")

    recomputed, rows = COMPARE.checked_stored_consensus(
        tmp_path, estimates, target_types, samples
    )

    pd.testing.assert_frame_equal(recomputed.primary, expected.primary)
    assert [row["scale"] for row in rows] == ["primary", "unclosed", "closed"]
    assert all(row["max_abs"] <= 1e-15 for row in rows)
