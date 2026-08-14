#!/usr/bin/env python3
"""Compare original-R and accelerated outputs on explicit aligned axes."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr, spearmanr

from decepticonx.consensus import ConsensusResult, build_consensus
from decepticonx.models import BranchKey


REFERENCE_ORDER = ("monocle3", "bayesprism", "music2")
REFERENCE_FILES = {
    "monocle3": "Monocle3_base.txt",
    "bayesprism": "BayesPrism_base.txt",
    "music2": "MuSiC2_base.txt",
}
ORIGINAL_PREFIXES = {
    "cibersort": "ciber_sig_",
    "cibersort_abs": "ciber_abs_sig_",
    "epic": "EPIC_sig_",
    "deconrnaseq": "Decon_sig_",
    "music": "music_sig_",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--accelerated-r-refs", type=Path, required=True)
    parser.add_argument("--accelerated-e2e", type=Path)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--label-map", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument(
        "--expected-epic-solver",
        choices=("auto", "nm", "nmf", "qp"),
        default="nm",
        help="required frozen-run EPIC solver (nm is the optimizer-fidelity default)",
    )
    return parser.parse_args()


def load_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    try:
        values = frame.to_numpy(dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"non-numeric table: {path}") from exc
    if frame.empty or not np.isfinite(values).all():
        raise ValueError(f"empty or non-finite table: {path}")
    return frame


def token(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def align_columns(frame: pd.DataFrame, target_types: list[str]) -> pd.DataFrame:
    lookup: dict[str, str] = {}
    for column in frame.columns:
        key = token(column)
        if key in lookup:
            raise ValueError(f"ambiguous normalized column token {key!r}")
        lookup[key] = column
    missing = [cell_type for cell_type in target_types if token(cell_type) not in lookup]
    if missing:
        raise ValueError("missing target cell types: " + ", ".join(missing))
    selected = frame[[lookup[token(cell_type)] for cell_type in target_types]].copy()
    selected.columns = target_types
    return selected


def align_samples(frame: pd.DataFrame, samples: pd.Index) -> pd.DataFrame:
    samples = pd.Index(samples.astype(str))
    if frame.index.has_duplicates:
        raise ValueError("duplicate sample identifiers")
    if samples.has_duplicates:
        raise ValueError("duplicate target sample identifiers")
    if set(frame.index) == set(samples):
        return frame.reindex(samples)

    def unique_suffixes(values: pd.Index, label: str) -> dict[str, str]:
        result: dict[str, str] = {}
        for identifier in values:
            match = re.search(r"([0-9]+)$", identifier)
            if match is None:
                raise ValueError(
                    f"{label} sample {identifier!r} has no terminal numeric suffix"
                )
            suffix = match.group(1)
            if suffix in result:
                raise ValueError(
                    f"{label} samples {result[suffix]!r} and {identifier!r} "
                    f"share numeric suffix {suffix!r}"
                )
            result[suffix] = identifier
        return result

    source_by_suffix = unique_suffixes(frame.index, "source")
    target_by_suffix = unique_suffixes(samples, "target")
    if set(source_by_suffix) != set(target_by_suffix):
        missing = sorted(set(target_by_suffix).difference(source_by_suffix))
        extra = sorted(set(source_by_suffix).difference(target_by_suffix))
        raise ValueError(
            "sample numeric-suffix mismatch; "
            f"missing={missing[:3]}, extra={extra[:3]}"
        )

    # Original DeconRNASeq/MuSiC output uses row names 1..N while the bulk
    # samples are value...1..N.  Only the complete, one-to-one terminal suffix
    # match above is allowed to bridge those identifiers.
    renamed = frame.rename(
        index={
            source: target_by_suffix[suffix]
            for suffix, source in source_by_suffix.items()
        }
    )
    return renamed.reindex(samples)


def load_r_manifest(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if list(frame.columns) != ["key", "value"]:
        raise ValueError(f"R manifest must contain exactly key/value columns: {path}")
    if frame["key"].duplicated().any() or (frame["key"] == "").any():
        raise ValueError(f"R manifest contains an empty or duplicate key: {path}")
    return dict(zip(frame["key"], frame["value"], strict=True))


def checked_sha256(value: object, label: str) -> str:
    normalized = str(value).lower()
    if re.fullmatch(r"[0-9a-f]{64}", normalized) is None:
        raise ValueError(f"{label} is not a SHA-256 digest")
    return normalized


def validate_frozen_provenance(
    original: Path,
    frozen_root: Path,
    reference_paths: Mapping[str, Path],
    *,
    expected_epic_solver: str = "nm",
) -> None:
    manifest = load_r_manifest(original / "manifest.tsv")
    required_manifest = {
        "seed",
        "bulk_sha256",
        "epic_mrna_cell_sha256",
        "epic_mrna_cell_entries",
    }
    missing_manifest = sorted(required_manifest.difference(manifest))
    if missing_manifest:
        raise ValueError(f"R manifest is missing keys: {missing_manifest}")

    run = json.loads((frozen_root / "run.json").read_text(encoding="utf-8"))
    provenance = run.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError("frozen run.json has no provenance object")

    try:
        r_seed = int(manifest["seed"])
        frozen_seed = int(provenance["seed"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("R/frozen provenance contains an invalid seed") from exc
    if r_seed != frozen_seed:
        raise ValueError(f"seed mismatch: R={r_seed}, frozen={frozen_seed}")

    r_bulk = checked_sha256(manifest["bulk_sha256"], "R bulk_sha256")
    frozen_bulk = checked_sha256(
        provenance.get("bulk_sha256"), "frozen bulk_sha256"
    )
    if r_bulk != frozen_bulk:
        raise ValueError(f"bulk SHA-256 mismatch: R={r_bulk}, frozen={frozen_bulk}")

    if provenance.get("epic_solver_requested") != expected_epic_solver:
        raise ValueError(
            "frozen comparison must request "
            f"epic_solver={expected_epic_solver!r}; got "
            f"{provenance.get('epic_solver_requested')!r}"
        )

    frozen_epic_mapping = provenance.get("epic_mrna_cell")
    if not isinstance(frozen_epic_mapping, dict):
        raise ValueError("frozen provenance has no redacted epic_mrna_cell metadata")
    r_mapping_hash = checked_sha256(
        manifest["epic_mrna_cell_sha256"], "R epic_mrna_cell_sha256"
    )
    frozen_mapping_hash = checked_sha256(
        frozen_epic_mapping.get("mapping_sha256"),
        "frozen epic_mrna_cell mapping_sha256",
    )
    try:
        r_mapping_entries = int(manifest["epic_mrna_cell_entries"])
        frozen_mapping_entries = int(frozen_epic_mapping["entry_count"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("R/frozen EPIC mapping entry count is invalid") from exc
    if r_mapping_entries <= 0 or r_mapping_entries != frozen_mapping_entries:
        raise ValueError(
            "EPIC mapping entry-count mismatch: "
            f"R={r_mapping_entries}, frozen={frozen_mapping_entries}"
        )
    if r_mapping_hash != frozen_mapping_hash:
        raise ValueError(
            "EPIC mapping SHA-256 mismatch: "
            f"R={r_mapping_hash}, frozen={frozen_mapping_hash}"
        )
    if frozen_epic_mapping.get("values_redacted") is not True:
        raise ValueError("frozen EPIC mapping provenance must mark values_redacted=true")

    recorded_refs = provenance.get("reference_sha256")
    if not isinstance(recorded_refs, dict) or set(recorded_refs) != set(reference_paths):
        raise ValueError(
            "frozen reference_sha256 must contain exactly: "
            + ", ".join(reference_paths)
        )
    for name, path in reference_paths.items():
        actual = sha256(path)
        recorded = checked_sha256(
            recorded_refs[name], f"frozen reference_sha256[{name!r}]"
        )
        if actual != recorded:
            raise ValueError(
                f"reference SHA-256 mismatch for {name}: "
                f"R file={actual}, frozen provenance={recorded}"
            )


def row_normalize(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame.to_numpy(dtype=np.float64)
    sums = values.sum(axis=1)
    if np.any(sums <= 0):
        raise ValueError("cannot normalize a zero-sum estimate row")
    return pd.DataFrame(values / sums[:, None], index=frame.index, columns=frame.columns)


def safe_pearson(left: np.ndarray, right: np.ndarray) -> float:
    if left.size < 2 or np.std(left) == 0 or np.std(right) == 0:
        return math.nan
    return float(pearsonr(left, right).statistic)


def matrix_metrics(left: pd.DataFrame, right: pd.DataFrame) -> dict[str, float]:
    if left.shape != right.shape or list(left.index) != list(right.index) or list(left.columns) != list(right.columns):
        raise ValueError("matrix axes are not identical")
    a = left.to_numpy(dtype=np.float64)
    b = right.to_numpy(dtype=np.float64)
    difference = a - b
    per_cell = [safe_pearson(a[:, index], b[:, index]) for index in range(a.shape[1])]
    finite_per_cell = [value for value in per_cell if np.isfinite(value)]
    return {
        "max_abs": float(np.max(np.abs(difference))),
        "mae": float(np.mean(np.abs(difference))),
        "rmse": float(np.sqrt(np.mean(np.square(difference)))),
        "pearson_flat": safe_pearson(a.ravel(), b.ravel()),
        "median_per_cell_pearson": (
            float(np.median(finite_per_cell)) if finite_per_cell else math.nan
        ),
        "left_row_sum_min": float(a.sum(axis=1).min()),
        "left_row_sum_median": float(np.median(a.sum(axis=1))),
        "left_row_sum_max": float(a.sum(axis=1).max()),
        "right_row_sum_min": float(b.sum(axis=1).min()),
        "right_row_sum_median": float(np.median(b.sum(axis=1))),
        "right_row_sum_max": float(b.sum(axis=1).max()),
    }


def original_estimates(
    original: Path, target_types: list[str], samples: pd.Index
) -> dict[BranchKey, pd.DataFrame]:
    estimates: dict[BranchKey, pd.DataFrame] = {}
    for method, prefix in ORIGINAL_PREFIXES.items():
        for index, reference in enumerate(REFERENCE_ORDER, start=1):
            path = original / "res" / f"{prefix}{index}.txt"
            frame = align_samples(align_columns(load_table(path), target_types), samples)
            estimates[BranchKey(method=method, reference=reference)] = frame
    if len(estimates) != 15:
        raise AssertionError("expected exactly 15 original branches")
    return estimates


def accelerated_estimates(
    root: Path, target_types: list[str], samples: pd.Index
) -> dict[BranchKey, pd.DataFrame]:
    estimates: dict[BranchKey, pd.DataFrame] = {}
    for method in ORIGINAL_PREFIXES:
        for reference in REFERENCE_ORDER:
            key = BranchKey(method=method, reference=reference)
            path = root / "estimates" / f"{key.slug}.tsv"
            frame = align_samples(align_columns(load_table(path), target_types), samples)
            estimates[key] = frame
    if len(estimates) != 15:
        raise AssertionError("expected exactly 15 accelerated branches")
    return estimates


def write_consensus(
    estimates: Mapping[BranchKey, pd.DataFrame], target_types: list[str], root: Path
) -> ConsensusResult:
    result = build_consensus(estimates, cell_types=target_types, n_pairs=2)
    result.normalized.to_csv(root / "original_common_consensus.tsv", sep="\t")
    result.raw.to_csv(root / "original_common_consensus_raw.tsv", sep="\t")
    (root / "original_common_consensus_diagnostics.json").write_text(
        json.dumps(result.diagnostics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def checked_stored_consensus(
    root: Path,
    estimates: Mapping[BranchKey, pd.DataFrame],
    target_types: list[str],
    samples: pd.Index,
) -> tuple[ConsensusResult, list[dict[str, object]]]:
    recomputed = build_consensus(estimates, cell_types=target_types, n_pairs=2)
    stored = {
        "normalized": align_samples(
            align_columns(load_table(root / "consensus.tsv"), target_types), samples
        ),
        "raw": align_samples(
            align_columns(load_table(root / "consensus_raw.tsv"), target_types), samples
        ),
    }
    expected = {"normalized": recomputed.normalized, "raw": recomputed.raw}
    rows: list[dict[str, object]] = []
    for scale in ("normalized", "raw"):
        metrics = matrix_metrics(expected[scale], stored[scale])
        rows.append({"scale": scale, **metrics})
        if not np.allclose(
            expected[scale].to_numpy(dtype=np.float64),
            stored[scale].to_numpy(dtype=np.float64),
            rtol=1e-12,
            atol=1e-12,
        ):
            raise ValueError(
                f"stored accelerated {scale} consensus does not match "
                f"a recomputation (max_abs={metrics['max_abs']})"
            )
    return recomputed, rows


def reference_metrics(original: Path, accelerated_e2e: Path | None) -> list[dict[str, object]]:
    if accelerated_e2e is None:
        return []
    rows: list[dict[str, object]] = []
    for reference in REFERENCE_ORDER:
        left = load_table(original / "custom_signature_matrix" / REFERENCE_FILES[reference])
        right = load_table(accelerated_e2e / "signatures" / f"{reference}.tsv")
        common_genes = [gene for gene in left.index if gene in right.index]
        common_types = [cell_type for cell_type in left.columns if token(cell_type) in {token(value) for value in right.columns}]
        if not common_genes or not common_types:
            raise ValueError(f"no common axes for reference {reference}")
        left_aligned = align_columns(left.loc[common_genes], common_types)
        right_aligned = align_columns(right.loc[common_genes], common_types)
        right_aligned = right_aligned.reindex(index=left_aligned.index)
        values = matrix_metrics(left_aligned, right_aligned)
        rows.append(
            {
                "reference": reference,
                "original_genes": len(left.index),
                "accelerated_genes": len(right.index),
                "common_genes": len(common_genes),
                "gene_jaccard": len(set(left.index) & set(right.index)) / len(set(left.index) | set(right.index)),
                "common_cell_types": len(common_types),
                **values,
            }
        )
    return rows


def truth_frame(path: Path) -> pd.DataFrame:
    frame = load_table(path)
    frame.index = frame.index.astype(str)
    values = frame.to_numpy(dtype=np.float64)
    if np.any(values < 0):
        raise ValueError("truth proportions must be non-negative")
    row_sums = values.sum(axis=1)
    if not np.allclose(row_sums, 1.0, rtol=0.0, atol=1e-8):
        raise ValueError("truth proportion rows must sum to one")
    return frame


def align_prediction_to_truth(
    prediction: pd.DataFrame, truth: pd.DataFrame, label_map: Mapping[str, str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if set(prediction.columns) != set(label_map):
        raise ValueError("prediction columns do not exactly match label-map keys")
    renamed = prediction.rename(columns=label_map)
    if set(renamed.columns) != set(truth.columns):
        raise ValueError("mapped prediction columns do not exactly match truth columns")
    sample_ids = prediction.index.to_series().str.extract(r"(\d+)$", expand=False)
    if sample_ids.isna().any() or sample_ids.duplicated().any():
        raise ValueError("prediction samples cannot be mapped to truth row IDs")
    renamed.index = sample_ids.to_numpy()
    if set(renamed.index) != set(truth.index):
        raise ValueError("prediction sample IDs do not exactly match truth")
    renamed = row_normalize(renamed.reindex(index=truth.index, columns=truth.columns))
    return renamed, truth.loc[renamed.index, renamed.columns]


def truth_metrics(name: str, prediction: pd.DataFrame, truth: pd.DataFrame) -> dict[str, object]:
    a = prediction.to_numpy(dtype=np.float64)
    b = truth.to_numpy(dtype=np.float64)
    per_cell_pearson = [safe_pearson(a[:, i], b[:, i]) for i in range(a.shape[1])]
    per_cell_spearman = [
        float(spearmanr(a[:, i], b[:, i]).statistic) for i in range(a.shape[1])
    ]
    difference = a - b
    jsd = [
        float(jensenshannon(a[i], b[i], base=2.0) ** 2) for i in range(a.shape[0])
    ]
    l1 = np.abs(difference).sum(axis=1)
    return {
        "run": name,
        "macro_pearson": float(np.nanmean(per_cell_pearson)),
        "macro_spearman": float(np.nanmean(per_cell_spearman)),
        "valid_pearson_cell_types": int(np.isfinite(per_cell_pearson).sum()),
        "valid_spearman_cell_types": int(np.isfinite(per_cell_spearman).sum()),
        "macro_mae": float(np.mean(np.mean(np.abs(difference), axis=0))),
        "macro_rmse": float(np.mean(np.sqrt(np.mean(np.square(difference), axis=0)))),
        "pooled_mae": float(np.mean(np.abs(difference))),
        "pooled_rmse": float(np.sqrt(np.mean(np.square(difference)))),
        "median_sample_l1": float(np.median(l1)),
        "mean_js_divergence": float(np.mean(jsd)),
    }


def truth_metrics_by_cell_type(
    name: str, prediction: pd.DataFrame, truth: pd.DataFrame
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for cell_type in prediction.columns:
        predicted = prediction[cell_type].to_numpy(dtype=np.float64)
        observed = truth[cell_type].to_numpy(dtype=np.float64)
        difference = predicted - observed
        spearman = spearmanr(predicted, observed).statistic
        rows.append(
            {
                "run": name,
                "cell_type": cell_type,
                "pearson": safe_pearson(predicted, observed),
                "spearman": float(spearman),
                "mae": float(np.mean(np.abs(difference))),
                "rmse": float(np.sqrt(np.mean(np.square(difference)))),
            }
        )
    return rows


def bootstrap_delta(
    left: pd.DataFrame,
    right: pd.DataFrame,
    truth: pd.DataFrame,
    *,
    repeats: int,
    seed: int,
) -> dict[str, float]:
    if repeats < 1:
        raise ValueError("--bootstrap must be at least 1")
    rng = np.random.default_rng(seed)
    a = left.to_numpy(dtype=np.float64)
    b = right.to_numpy(dtype=np.float64)
    t = truth.to_numpy(dtype=np.float64)
    deltas = np.empty(repeats, dtype=np.float64)
    for index in range(repeats):
        selected = rng.integers(0, len(t), size=len(t))
        left_rmse = np.sqrt(np.mean(np.square(a[selected] - t[selected])))
        right_rmse = np.sqrt(np.mean(np.square(b[selected] - t[selected])))
        deltas[index] = right_rmse - left_rmse
    return {
        "delta_pooled_rmse_right_minus_left": float(
            np.sqrt(np.mean(np.square(b - t))) - np.sqrt(np.mean(np.square(a - t)))
        ),
        "bootstrap_ci_low": float(np.quantile(deltas, 0.025)),
        "bootstrap_ci_high": float(np.quantile(deltas, 0.975)),
    }


def main() -> int:
    args = arguments()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError(f"output directory must be absent or empty: {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)

    original_reference_paths = {
        name: args.original / "custom_signature_matrix" / filename
        for name, filename in REFERENCE_FILES.items()
    }
    original_refs = {
        name: load_table(path) for name, path in original_reference_paths.items()
    }
    validate_frozen_provenance(
        args.original,
        args.accelerated_r_refs,
        original_reference_paths,
        expected_epic_solver=args.expected_epic_solver,
    )
    target_types = [
        cell_type
        for cell_type in original_refs["monocle3"].columns
        if all(token(cell_type) in {token(value) for value in frame.columns} for frame in original_refs.values())
    ]
    if len(target_types) != 11:
        raise ValueError(f"expected 11 common cell types, found {len(target_types)}")

    first_accelerated = load_table(
        args.accelerated_r_refs / "estimates" / "cibersort__monocle3.tsv"
    )
    samples = first_accelerated.index
    if len(samples) != 100:
        raise ValueError(f"expected 100 samples, found {len(samples)}")

    r_estimates = original_estimates(args.original, target_types, samples)
    accelerated = accelerated_estimates(args.accelerated_r_refs, target_types, samples)
    branch_rows: list[dict[str, object]] = []
    for key in sorted(r_estimates, key=lambda value: (value.method, value.reference)):
        raw = matrix_metrics(r_estimates[key], accelerated[key])
        normalized = matrix_metrics(row_normalize(r_estimates[key]), row_normalize(accelerated[key]))
        branch_rows.append(
            {
                "method": key.method,
                "reference": key.reference,
                **{f"raw_{name}": value for name, value in raw.items()},
                **{f"normalized_{name}": value for name, value in normalized.items()},
            }
        )
    pd.DataFrame(branch_rows).to_csv(args.output / "branch_metrics.tsv", sep="\t", index=False)

    accelerated_e2e: dict[BranchKey, pd.DataFrame] | None = None
    if args.accelerated_e2e is not None:
        accelerated_e2e = accelerated_estimates(
            args.accelerated_e2e, target_types, samples
        )
        e2e_rows: list[dict[str, object]] = []
        for key in sorted(r_estimates, key=lambda value: (value.method, value.reference)):
            raw = matrix_metrics(r_estimates[key], accelerated_e2e[key])
            normalized = matrix_metrics(
                row_normalize(r_estimates[key]), row_normalize(accelerated_e2e[key])
            )
            e2e_rows.append(
                {
                    "method": key.method,
                    "reference": key.reference,
                    **{f"raw_{name}": value for name, value in raw.items()},
                    **{f"normalized_{name}": value for name, value in normalized.items()},
                }
            )
        pd.DataFrame(e2e_rows).to_csv(
            args.output / "branch_metrics_e2e.tsv", sep="\t", index=False
        )

    ref_rows = reference_metrics(args.original, args.accelerated_e2e)
    pd.DataFrame(ref_rows).to_csv(args.output / "reference_metrics.tsv", sep="\t", index=False)

    original_consensus = write_consensus(r_estimates, target_types, args.output)
    solver_consensus, stored_consensus_rows = checked_stored_consensus(
        args.accelerated_r_refs,
        accelerated,
        target_types,
        samples,
    )
    pd.DataFrame(stored_consensus_rows).to_csv(
        args.output / "stored_consensus_validation.tsv", sep="\t", index=False
    )
    consensus_rows = []
    for scale in ("normalized", "raw"):
        metrics = matrix_metrics(
            getattr(original_consensus, scale), getattr(solver_consensus, scale)
        )
        consensus_rows.append({"scale": scale, **metrics})
    pd.DataFrame(consensus_rows).to_csv(
        args.output / "common_consensus_metrics.tsv", sep="\t", index=False
    )
    e2e_consensus: ConsensusResult | None = None
    if args.accelerated_e2e is not None and accelerated_e2e is not None:
        e2e_consensus, e2e_stored_rows = checked_stored_consensus(
            args.accelerated_e2e,
            accelerated_e2e,
            target_types,
            samples,
        )
        pd.DataFrame(e2e_stored_rows).to_csv(
            args.output / "e2e_stored_consensus_validation.tsv",
            sep="\t",
            index=False,
        )
        e2e_consensus_rows = []
        for scale in ("normalized", "raw"):
            metrics = matrix_metrics(
                getattr(original_consensus, scale), getattr(e2e_consensus, scale)
            )
            e2e_consensus_rows.append({"scale": scale, **metrics})
        pd.DataFrame(e2e_consensus_rows).to_csv(
            args.output / "e2e_consensus_metrics.tsv", sep="\t", index=False
        )
    predictions = {
        "r_original_common5": original_consensus.normalized,
        f"accelerated_on_r_refs_{args.expected_epic_solver}": solver_consensus.normalized,
    }
    if e2e_consensus is not None:
        predictions["accelerated_corrected_e2e"] = e2e_consensus.normalized

    truth = truth_frame(args.truth)
    label_map = json.loads(args.label_map.read_text(encoding="utf-8"))
    aligned_predictions: dict[str, pd.DataFrame] = {}
    rows = []
    cell_type_rows: list[dict[str, object]] = []
    truth_aligned: pd.DataFrame | None = None
    for name, prediction in predictions.items():
        aligned, current_truth = align_prediction_to_truth(prediction, truth, label_map)
        aligned_predictions[name] = aligned
        truth_aligned = current_truth
        rows.append(truth_metrics(name, aligned, current_truth))
        cell_type_rows.extend(truth_metrics_by_cell_type(name, aligned, current_truth))
    pd.DataFrame(rows).to_csv(args.output / "truth_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(cell_type_rows).to_csv(
        args.output / "truth_metrics_by_cell_type.tsv", sep="\t", index=False
    )

    if truth_aligned is None:
        raise AssertionError("truth alignment did not run")
    comparisons: list[dict[str, object]] = []
    baseline = aligned_predictions["r_original_common5"]
    for name, prediction in aligned_predictions.items():
        if name == "r_original_common5":
            continue
        comparisons.append(
            {
                "left": "r_original_common5",
                "right": name,
                **bootstrap_delta(
                    baseline,
                    prediction,
                    truth_aligned,
                    repeats=args.bootstrap,
                    seed=args.seed,
                ),
            }
        )
    pd.DataFrame(comparisons).to_csv(args.output / "bootstrap_delta.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
