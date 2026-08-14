#!/usr/bin/env python3
"""Score a MuSiC2 validation run and reconstruct both consensus semantics."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
from pathlib import Path
import subprocess
import sys
from typing import Sequence

import numpy as np
import pandas as pd

from decepticonx.consensus import build_consensus
from decepticonx.io import load_bulk_expression, prepare_single_cell
from decepticonx.models import BranchKey
from decepticonx.references import build_references


SOURCE_NAMES = (
    "music2_bulk_eset.rds",
    "music2_EMTABsce_healthy.rds",
    "music2_true_proportion.RData",
)
EXPECTED_METHODS = ("cibersort", "cibersort_abs", "deconrnaseq", "music")
EXPECTED_REFERENCES = ("bayesprism", "monocle3", "music2")
TARGET_TYPES = ("acinar", "alpha", "beta", "delta", "ductal", "gamma")
CONSENSUS_ATOL = 1e-12


def read_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError(f"{path} contains duplicate row or column identifiers")
    frame = frame.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(frame.to_numpy(dtype=np.float64)).all():
        raise ValueError(f"{path} contains non-finite values")
    return frame


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def detect_package_commit() -> str:
    repository = Path(__file__).resolve().parents[2]
    expected_package = (repository / "src" / "decepticonx").resolve()
    imported_source = Path(inspect.getfile(build_references)).resolve()
    if not imported_source.is_relative_to(expected_package):
        raise ValueError(
            "the imported decepticonx package is not from this checkout; "
            "install it with `python -m pip install -e .` or pass "
            "--package-commit explicitly"
        )
    try:
        revision = subprocess.run(
            ["git", "-C", str(repository), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        status = subprocess.run(
            [
                "git", "-C", str(repository), "status", "--porcelain",
                "--untracked-files=no",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        raise ValueError(
            "could not determine the package commit; pass --package-commit"
        ) from None
    if status.stdout.strip():
        raise ValueError(
            "the package checkout has tracked modifications; commit them or "
            "pass --package-commit explicitly"
        )
    value = revision.stdout.strip()
    if not value:
        raise ValueError("git returned an empty package commit")
    return value


def close(frame: pd.DataFrame, cell_types: list[str]) -> pd.DataFrame:
    values = frame.loc[:, cell_types].to_numpy(dtype=np.float64)
    if np.any(values < -CONSENSUS_ATOL):
        raise ValueError("cannot close an estimate containing negative values")
    sums = values.sum(axis=1)
    if np.any(sums <= 0):
        raise ValueError("cannot close an estimate with a non-positive row sum")
    return pd.DataFrame(
        values / sums[:, None],
        index=frame.index,
        columns=cell_types,
    )


def aligned_values(
    frame: pd.DataFrame,
    truth: pd.DataFrame,
    sample_ids: pd.Index,
    cell_types: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    expected = truth.loc[sample_ids, cell_types].to_numpy(dtype=np.float64)
    observed = frame.loc[sample_ids, cell_types].to_numpy(dtype=np.float64)
    return observed, expected


def metric_row(
    name: str,
    kind: str,
    frame: pd.DataFrame,
    group: str,
    truth: pd.DataFrame,
    metadata: pd.DataFrame,
    cell_types: list[str],
) -> dict[str, object]:
    if group == "all":
        sample_ids = truth.index
    else:
        sample_ids = metadata.index[metadata["group"].str.lower() == group]
    observed, expected = aligned_values(frame, truth, sample_ids, cell_types)
    error = observed - expected
    sample_rmse = np.sqrt(np.mean(np.square(error), axis=1))
    pearson = float(np.corrcoef(observed.ravel(), expected.ravel())[0, 1])
    row_sums = observed.sum(axis=1)
    return {
        "artifact": name,
        "kind": kind,
        "group": group,
        "samples": len(sample_ids),
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "mae": float(np.mean(np.abs(error))),
        "pearson": pearson,
        "mean_sample_rmse": float(np.mean(sample_rmse)),
        "median_sample_rmse": float(np.median(sample_rmse)),
        "max_sample_rmse": float(np.max(sample_rmse)),
        "row_sum_min": float(np.min(row_sums)),
        "row_sum_max": float(np.max(row_sums)),
        "finite": bool(np.isfinite(observed).all()),
        "minimum": float(np.min(observed)),
    }


def analyze(
    *,
    source_data_dir: Path,
    prepared_dir: Path,
    run_dir: Path,
    output_dir: Path,
    corrected_output_dir: Path,
    package_commit: str | None,
) -> dict[str, object]:
    """Analyze one completed run and write only compact derived artifacts."""

    run_resolved = run_dir.resolve()
    for label, destination in (
        ("output-dir", output_dir),
        ("corrected-output-dir", corrected_output_dir),
    ):
        resolved = destination.resolve()
        if resolved == run_resolved or run_resolved in resolved.parents:
            raise ValueError(f"{label} must be outside the source run directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    corrected_output_dir.mkdir(parents=True, exist_ok=True)

    truth = read_table(prepared_dir / "truth.tsv")
    if len(truth) != 100 or set(truth.columns) != set(TARGET_TYPES):
        raise ValueError("truth must contain exactly 100 samples and six target types")
    truth_values = truth.to_numpy(dtype=np.float64)
    if np.any(truth_values < 0):
        raise ValueError("truth contains negative proportions")
    if np.max(np.abs(truth_values.sum(axis=1) - 1.0)) > CONSENSUS_ATOL:
        raise ValueError("truth rows do not sum to one")
    metadata = pd.read_csv(
        prepared_dir / "sample_metadata.tsv", sep="\t", dtype=str
    ).set_index("sample")
    metadata.index = metadata.index.astype(str)
    if metadata.index.has_duplicates:
        raise ValueError("sample metadata contains duplicate sample IDs")
    if set(metadata.index) != set(truth.index):
        raise ValueError("sample metadata and truth IDs differ")
    if "group" not in metadata.columns:
        raise ValueError("sample metadata is missing the group column")
    group_counts = metadata["group"].str.lower().value_counts().to_dict()
    if group_counts != {"healthy": 50, "t2d": 50}:
        raise ValueError(
            "sample metadata must contain exactly 50 healthy and 50 T2D samples"
        )

    run_payload = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    consensus_diagnostics = run_payload["diagnostics"]["consensus"]
    if consensus_diagnostics.get("mode") != "r_literal":
        raise ValueError("the validation run must use r_literal consensus mode")
    if consensus_diagnostics.get("n_pairs_requested") != 2:
        raise ValueError("the validation run must request exactly two consensus pairs")
    strategy_slugs = list(consensus_diagnostics["strategy_order"])

    estimates: dict[BranchKey, pd.DataFrame] = {}
    keys_by_slug: dict[str, BranchKey] = {}
    for path in sorted((run_dir / "estimates").glob("*.tsv")):
        method, reference = path.stem.split("__", 1)
        key = BranchKey(method=method, reference=reference)
        estimates[key] = read_table(path)
        keys_by_slug[key.slug] = key
    if set(strategy_slugs) != set(keys_by_slug):
        raise ValueError("stored strategy order and estimate files differ")
    expected_slugs = {
        f"{method}__{reference}"
        for method in EXPECTED_METHODS
        for reference in EXPECTED_REFERENCES
    }
    if set(keys_by_slug) != expected_slugs:
        missing = sorted(expected_slugs.difference(keys_by_slug))
        extra = sorted(set(keys_by_slug).difference(expected_slugs))
        raise ValueError(
            f"expected the 12 frozen branches (missing={missing}, extra={extra})"
        )
    backend_diagnostics = run_payload["diagnostics"]["backends"]
    if backend_diagnostics["skipped"]:
        raise ValueError(
            f"the validation run skipped branches: {backend_diagnostics['skipped']}"
        )
    if set(backend_diagnostics["completed"]) != expected_slugs:
        raise ValueError("completed-branch diagnostics do not match the frozen set")
    expected_engines = {
        slug: (
            "rust"
            if slug.startswith(("cibersort__", "cibersort_abs__"))
            else "numpy"
        )
        for slug in expected_slugs
    }
    if backend_diagnostics["engines"] != expected_engines:
        raise ValueError(
            "resolved backend engines differ from the frozen run: "
            f"{backend_diagnostics['engines']}"
        )
    for slug, frame in ((key.slug, value) for key, value in estimates.items()):
        if set(frame.index) != set(truth.index) or len(frame) != len(truth):
            raise ValueError(f"{slug} sample IDs differ from truth")
        if set(frame.columns) != set(TARGET_TYPES) or len(frame.columns) != 6:
            raise ValueError(f"{slug} does not contain exactly the six target types")
    strategy_order = [keys_by_slug[slug] for slug in strategy_slugs]

    saved_primary = read_table(run_dir / "consensus.tsv")
    saved_unclosed = read_table(run_dir / "consensus_unclosed.tsv")
    saved_closed = read_table(run_dir / "consensus_closed.tsv")
    cell_types = list(saved_primary.columns)
    if set(cell_types) != set(truth.columns):
        raise ValueError("consensus and truth cell types differ")
    for label, frame in (
        ("primary consensus", saved_primary),
        ("unclosed consensus", saved_unclosed),
        ("closed consensus", saved_closed),
    ):
        if set(frame.index) != set(truth.index) or len(frame) != len(truth):
            raise ValueError(f"{label} sample IDs differ from truth")
        if set(frame.columns) != set(TARGET_TYPES) or len(frame.columns) != 6:
            raise ValueError(f"{label} does not contain exactly the six target types")

    recomputed_literal = build_consensus(
        estimates,
        cell_types=cell_types,
        n_pairs=2,
        mode="r_literal",
        strategy_order=strategy_order,
    )
    corrected = build_consensus(
        estimates,
        cell_types=cell_types,
        n_pairs=2,
        mode="corrected",
        strategy_order=strategy_order,
    )
    structural_checks = {
        "saved_primary_equals_unclosed_max_abs": float(
            np.max(np.abs(saved_primary.to_numpy() - saved_unclosed.to_numpy()))
        ),
        "r_literal_recompute_unclosed_max_abs": float(
            np.max(
                np.abs(
                    recomputed_literal.unclosed.to_numpy()
                    - saved_unclosed.to_numpy()
                )
            )
        ),
        "r_literal_recompute_closed_max_abs": float(
            np.max(
                np.abs(recomputed_literal.closed.to_numpy() - saved_closed.to_numpy())
            )
        ),
        "r_literal_closed_max_row_sum_error": float(
            np.max(np.abs(saved_closed.sum(axis=1).to_numpy() - 1.0))
        ),
        "corrected_primary_equals_closed_max_abs": float(
            np.max(np.abs(corrected.primary.to_numpy() - corrected.closed.to_numpy()))
        ),
        "corrected_closed_max_row_sum_error": float(
            np.max(np.abs(corrected.closed.sum(axis=1).to_numpy() - 1.0))
        ),
    }
    failed_checks = {
        name: value
        for name, value in structural_checks.items()
        if not np.isfinite(value) or value > CONSENSUS_ATOL
    }
    if failed_checks:
        raise ValueError(
            "consensus reconstruction failed structural checks: "
            f"{failed_checks}"
        )

    corrected.primary.to_csv(corrected_output_dir / "consensus.tsv", sep="\t")
    corrected.unclosed.to_csv(
        corrected_output_dir / "consensus_unclosed.tsv", sep="\t"
    )
    corrected.closed.to_csv(
        corrected_output_dir / "consensus_closed.tsv", sep="\t"
    )
    (corrected_output_dir / "diagnostics.json").write_text(
        json.dumps(corrected.diagnostics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    artifacts: list[tuple[str, str, pd.DataFrame]] = [
        ("r_literal_unclosed", "consensus_unclosed", saved_unclosed),
        ("r_literal_closed", "consensus_closed", saved_closed),
        ("corrected_unclosed", "consensus_unclosed", corrected.unclosed),
        ("corrected_closed", "consensus_closed", corrected.closed),
    ]
    for key in strategy_order:
        artifacts.append((key.slug, "branch_closed", close(estimates[key], cell_types)))

    metric_rows: list[dict[str, object]] = []
    cell_type_rows: list[dict[str, object]] = []
    for name, kind, frame in artifacts:
        for group in ("all", "healthy", "t2d"):
            metric_rows.append(
                metric_row(
                    name, kind, frame, group, truth, metadata, cell_types
                )
            )
            if group == "all":
                sample_ids = truth.index
            else:
                sample_ids = metadata.index[metadata["group"].str.lower() == group]
            observed, expected = aligned_values(
                frame, truth, sample_ids, cell_types
            )
            for position, cell_type in enumerate(cell_types):
                delta = observed[:, position] - expected[:, position]
                cell_type_rows.append(
                    {
                        "artifact": name,
                        "kind": kind,
                        "group": group,
                        "cell_type": cell_type,
                        "rmse": float(np.sqrt(np.mean(np.square(delta)))),
                        "mae": float(np.mean(np.abs(delta))),
                        "pearson": float(
                            np.corrcoef(
                                observed[:, position], expected[:, position]
                            )[0, 1]
                        ),
                    }
                )

    metrics = pd.DataFrame(metric_rows).sort_values(["group", "rmse", "artifact"])
    cell_metrics = pd.DataFrame(cell_type_rows).sort_values(
        ["group", "cell_type", "rmse", "artifact"]
    )
    metrics.to_csv(output_dir / "metrics.tsv", sep="\t", index=False)
    cell_metrics.to_csv(output_dir / "cell_type_metrics.tsv", sep="\t", index=False)

    # Quantify the sole input intervention: the explicit rounding of sparse
    # half-integers in the official SCE for the strict raw-count contract.
    bulk = load_bulk_expression(prepared_dir / "bulk.tsv")
    prepared_rounded = prepare_single_cell(
        prepared_dir / "reference.h5ad",
        bulk,
        cell_type_key="cell_type",
        sample_key="donor",
        counts_layer="auto",
    )
    rounded_refs, _ = build_references(prepared_rounded, bulk)
    prepared_exact = prepare_single_cell(
        prepared_dir / "reference.h5ad",
        bulk,
        cell_type_key="cell_type",
        sample_key="donor",
        counts_layer="X",
        allow_normalized_x=True,
    )
    exact_refs, _ = build_references(prepared_exact, bulk)
    reference_sensitivity: dict[str, dict[str, float]] = {}
    for name, rounded in rounded_refs.items():
        exact = exact_refs[name].loc[rounded.index, rounded.columns]
        difference = rounded.to_numpy() - exact.to_numpy()
        exact_values = exact.to_numpy()
        denominator = float(np.linalg.norm(exact_values))
        reference_sensitivity[name] = {
            "max_abs": float(np.max(np.abs(difference))),
            "mean_abs": float(np.mean(np.abs(difference))),
            "relative_frobenius": float(np.linalg.norm(difference) / denominator),
            "pearson": float(
                np.corrcoef(rounded.to_numpy().ravel(), exact_values.ravel())[0, 1]
            ),
        }

    source_paths = [source_data_dir / name for name in SOURCE_NAMES]
    converted_paths = [
        prepared_dir / "reference.h5ad",
        prepared_dir / "bulk.tsv",
        prepared_dir / "truth.tsv",
    ]
    audit: dict[str, object] = {
        "main_commit": package_commit,
        "branch_count": len(estimates),
        "completed_branches": run_payload["diagnostics"]["backends"]["completed"],
        "skipped_branches": run_payload["diagnostics"]["backends"]["skipped"],
        "engines": run_payload["diagnostics"]["backends"]["engines"],
        "timings_seconds": run_payload["timings_seconds"],
        **structural_checks,
        "corrected_fallbacks": corrected.diagnostics["fallbacks"],
        "reference_rounding_sensitivity": reference_sensitivity,
        "source_hashes_sha256": {
            path.name: sha256_file(path) for path in (*source_paths, *converted_paths)
        },
    }
    (output_dir / "audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print("STRUCTURAL AUDIT")
    print(
        json.dumps(
            {key: value for key, value in audit.items() if key != "source_hashes_sha256"},
            indent=2,
        )
    )
    print("\nTOP OVERALL METRICS")
    print(metrics.loc[metrics["group"] == "all"].head(20).to_string(index=False))
    return audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct both consensus modes and score a completed MuSiC2 "
            "validation run against its known proportions."
        )
    )
    parser.add_argument("--source-data-dir", type=Path, required=True)
    parser.add_argument("--prepared-dir", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--corrected-output-dir",
        type=Path,
        default=None,
        help="corrected consensus output (default: OUTPUT_DIR/corrected_consensus)",
    )
    parser.add_argument(
        "--package-commit",
        default=None,
        help="commit recorded in audit.json (default: auto-detect this checkout)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    corrected_output_dir = args.corrected_output_dir or (
        args.output_dir / "corrected_consensus"
    )
    try:
        analyze(
            source_data_dir=args.source_data_dir,
            prepared_dir=args.prepared_dir,
            run_dir=args.run_dir,
            output_dir=args.output_dir,
            corrected_output_dir=corrected_output_dir,
            package_commit=args.package_commit or detect_package_commit(),
        )
    except (KeyError, OSError, TypeError, ValueError) as exc:
        print(f"analyze_music2_results.py: error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
