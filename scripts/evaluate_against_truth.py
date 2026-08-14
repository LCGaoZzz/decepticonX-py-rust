#!/usr/bin/env python3
"""Evaluate a consensus TSV against a samples-by-cell-types truth TSV."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _read(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=None, engine="python", index_col=0)
    frame = frame.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(frame.to_numpy()).all():
        raise ValueError(f"{path} contains NA/Inf")
    return frame


def _sample_suffix(value: object) -> str:
    text = str(value)
    return text.rsplit("...", 1)[-1]


def _canonical_sample_index(frame: pd.DataFrame, *, label: str) -> pd.DataFrame:
    result = frame.copy()
    result.index = pd.Index([_sample_suffix(item) for item in result.index])
    if result.index.has_duplicates:
        duplicates = result.index[result.index.duplicated(keep=False)].unique().tolist()
        raise ValueError(
            f"{label} contains duplicate sample IDs after suffix normalization: "
            f"{duplicates[:5]}"
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("prediction", type=Path)
    parser.add_argument("truth", type=Path)
    parser.add_argument("--label-map", type=Path, help="JSON mapping prediction labels to truth labels")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    prediction = _read(args.prediction)
    truth = _read(args.truth)
    if args.label_map:
        prediction = prediction.rename(columns=json.loads(args.label_map.read_text(encoding="utf-8")))

    prediction = _canonical_sample_index(prediction, label="prediction")
    truth = _canonical_sample_index(truth, label="truth")
    samples = prediction.index.intersection(truth.index, sort=False)
    cell_types = prediction.columns.intersection(truth.columns, sort=False)
    if samples.empty or cell_types.empty:
        raise ValueError("Prediction and truth have no shared samples/cell types")

    prediction = prediction.loc[samples, cell_types]
    truth = truth.loc[samples, cell_types]
    rows: list[dict[str, float | str]] = []
    for cell_type in cell_types:
        pred = prediction[cell_type].to_numpy(float)
        actual = truth[cell_type].to_numpy(float)
        corr = float(np.corrcoef(pred, actual)[0, 1]) if np.std(pred) and np.std(actual) else float("nan")
        rows.append(
            {
                "cell_type": cell_type,
                "pearson": corr,
                "rmse": float(np.sqrt(np.mean((pred - actual) ** 2))),
                "mae": float(np.mean(np.abs(pred - actual))),
            }
        )
    report = pd.DataFrame(rows).set_index("cell_type")
    report.loc["__macro__"] = report.mean(axis=0, skipna=True)
    if args.output:
        report.to_csv(args.output, sep="\t")
    else:
        print(report.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
