#!/usr/bin/env python3
"""Run the five accelerated kernels against frozen original-R references."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import random
import re
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd

from decepticonx.backends import run_backends
from decepticonx.consensus import build_consensus
from decepticonx.io import load_bulk_expression


REFERENCE_FILES = {
    "monocle3": "Monocle3_base.txt",
    "bayesprism": "BayesPrism_base.txt",
    "music2": "MuSiC2_base.txt",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument(
        "--epic-backend", choices=("auto", "rust", "python"), default="rust"
    )
    parser.add_argument(
        "--epic-solver", choices=("auto", "nm", "nmf", "qp"), default="nm"
    )
    parser.add_argument(
        "--epic-mrna-json",
        help=(
            "authorized EPIC mRNA-per-cell mapping as an inline JSON object or "
            "a JSON file path; default is the explicit uniform smoke mapping"
        ),
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def token(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def load_epic_mrna_mapping(argument: str | None) -> dict[str, float]:
    if argument is None:
        return {"default": 1.0}
    stripped = argument.lstrip()
    candidate = None if stripped.startswith("{") else Path(argument).expanduser()
    if candidate is not None and candidate.is_file():
        content = candidate.read_text(encoding="utf-8")
        label = str(candidate)
    else:
        content = argument
        label = "--epic-mrna-json"
    try:
        raw = json.loads(content)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} is not valid JSON") from exc
    if not isinstance(raw, dict) or not raw:
        raise ValueError(f"{label} must contain a non-empty JSON object")
    result: dict[str, float] = {}
    for key, value in raw.items():
        if not isinstance(key, str) or not key or isinstance(value, bool):
            raise ValueError(f"{label} contains an invalid key/value entry")
        try:
            numeric = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{label} value for {key!r} is not numeric") from exc
        if not math.isfinite(numeric) or numeric <= 0:
            raise ValueError(f"{label} value for {key!r} must be finite and positive")
        result[key] = numeric
    return result


def mapping_sha256(mapping: dict[str, float]) -> str:
    # Fixed 17-significant-digit text is stable across Python and R's
    # ``sprintf("%.17g")``.  JSON encoders do not promise identical float
    # spellings across languages, so they are unsuitable for a cross-runtime
    # provenance gate.
    canonical = "\n".join(
        f"{key}\t{float(mapping[key]):.17g}" for key in sorted(mapping)
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def redacted_command() -> list[str]:
    command = [sys.executable]
    redact_next = False
    for argument in sys.argv:
        if redact_next:
            command.append("<redacted>")
            redact_next = False
        elif argument == "--epic-mrna-json":
            command.append(argument)
            redact_next = True
        elif argument.startswith("--epic-mrna-json="):
            command.append("--epic-mrna-json=<redacted>")
        else:
            command.append(argument)
    return command


def canonicalize_references(
    references: dict[str, pd.DataFrame], target_name: str = "monocle3"
) -> tuple[dict[str, pd.DataFrame], list[str]]:
    target = references[target_name]
    if len(target.columns) != 11:
        raise ValueError(f"{target_name} reference must contain exactly 11 columns")
    target_types = list(target.columns)
    target_tokens = [token(value) for value in target_types]
    if any(not value for value in target_tokens) or len(set(target_tokens)) != 11:
        raise ValueError(
            f"{target_name} reference must contain 11 one-to-one normalized labels"
        )

    canonical: dict[str, pd.DataFrame] = {}
    for name, frame in references.items():
        if len(frame.columns) != 11:
            raise ValueError(f"reference {name!r} must contain exactly 11 columns")
        lookup: dict[str, str] = {}
        for column in frame.columns:
            normalized = token(column)
            if not normalized or normalized in lookup:
                raise ValueError(
                    f"reference {name!r} does not have one-to-one normalized labels"
                )
            lookup[normalized] = column
        if set(lookup) != set(target_tokens):
            missing = sorted(set(target_tokens).difference(lookup))
            extra = sorted(set(lookup).difference(target_tokens))
            raise ValueError(
                f"reference {name!r} label-token mismatch; "
                f"missing={missing}, extra={extra}"
            )
        aligned = frame[[lookup[value] for value in target_tokens]].copy()
        aligned.columns = target_types
        canonical[name] = aligned
    return canonical, target_types


def load_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    values = frame.to_numpy(dtype=np.float64)
    if frame.empty or not np.isfinite(values).all():
        raise ValueError(f"{path} is empty or contains non-finite values")
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError(f"{path} contains duplicate row or column identifiers")
    return frame


def reference_path(root: Path, name: str, original_filename: str) -> Path:
    original = root / original_filename
    if original.is_file():
        return original
    stable = root / f"{name}.tsv"
    if stable.is_file():
        return stable
    raise FileNotFoundError(f"missing reference {name!r} in {root}")


def json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def package_versions() -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    for distribution in (
        "decepticonx-py-rust",
        "decepticon-fast",
        "python-cibersort-rs",
        "epic_py",
        "deconrnaseq-py",
        "deconrnaseq-rust",
        "music-py",
        "music-py-rust",
    ):
        try:
            result[distribution] = version(distribution)
        except PackageNotFoundError:
            result[distribution] = None
    return result


def main() -> int:
    args = arguments()
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    if args.output.exists():
        if not args.output.is_dir() or any(args.output.iterdir()):
            raise ValueError(f"output directory must be absent or empty: {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)

    os.environ["RAYON_NUM_THREADS"] = str(args.threads)
    random.seed(args.seed)
    np.random.seed(args.seed)
    epic_mrna_cell = load_epic_mrna_mapping(args.epic_mrna_json)
    started = perf_counter()
    timings: dict[str, float] = {}

    stage = perf_counter()
    bulk = load_bulk_expression(args.bulk)
    timings["bulk_input"] = perf_counter() - stage

    stage = perf_counter()
    reference_paths = {
        name: reference_path(args.references, name, filename)
        for name, filename in REFERENCE_FILES.items()
    }
    loaded_references = {
        name: load_table(path) for name, path in reference_paths.items()
    }
    references, target_types = canonicalize_references(loaded_references)
    timings["reference_input"] = perf_counter() - stage

    estimates = {}
    diagnostics: dict[str, Any] = {}
    method_groups = (
        ("cibersort_pair", ("cibersort", "cibersort_abs")),
        ("epic", ("epic",)),
        ("deconrnaseq", ("deconrnaseq",)),
        ("music", ("music",)),
    )
    for label, methods in method_groups:
        random.seed(args.seed)
        np.random.seed(args.seed)
        stage = perf_counter()
        values, details = run_backends(
            references,
            bulk,
            methods=methods,
            qn=True,
            seed=args.seed,
            threads=args.threads,
            # The default is a uniform smoke-test mapping. Optimizer fidelity needs
            # the same authorized effective mapping used by the R EPIC run.
            epic_mrna_cell=epic_mrna_cell,
            strict=True,
            cibersort_engine="rust",
            deconrnaseq_backend="rust",
            epic_backend=args.epic_backend,
            epic_solver=args.epic_solver,
            music_backend="rust",
        )
        timings[f"deconv_{label}"] = perf_counter() - stage
        estimates.update(values)
        diagnostics[label] = details

    stage = perf_counter()
    consensus = build_consensus(estimates, cell_types=target_types, n_pairs=2)
    timings["consensus"] = perf_counter() - stage
    timings["compute_total"] = perf_counter() - started

    stage = perf_counter()
    signature_dir = args.output / "signatures"
    estimate_dir = args.output / "estimates"
    signature_dir.mkdir()
    estimate_dir.mkdir()
    for name, frame in references.items():
        frame.to_csv(signature_dir / f"{name}.tsv", sep="\t")
    for key, frame in estimates.items():
        frame.to_csv(estimate_dir / f"{key.slug}.tsv", sep="\t")
    consensus.normalized.to_csv(args.output / "consensus.tsv", sep="\t")
    consensus.raw.to_csv(args.output / "consensus_raw.tsv", sep="\t")
    timings["write"] = perf_counter() - stage
    timings["wall_internal"] = perf_counter() - started

    run = {
        "timings_seconds": timings,
        "diagnostics": diagnostics,
        "consensus_diagnostics": consensus.diagnostics,
        "provenance": {
            "kind": "frozen-original-reference solver comparison",
            "python": platform.python_version(),
            "platform": platform.platform(),
            "threads": args.threads,
            "seed": args.seed,
            "bulk_sha256": sha256(args.bulk),
            "reference_sha256": {
                name: sha256(path) for name, path in reference_paths.items()
            },
            "component_versions": package_versions(),
            "epic_backend_requested": args.epic_backend,
            "epic_solver_requested": args.epic_solver,
            "thread_environment": {
                name: os.environ.get(name)
                for name in (
                    "RAYON_NUM_THREADS",
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            },
            "epic_mrna_cell": {
                "mapping_sha256": mapping_sha256(epic_mrna_cell),
                "entry_count": len(epic_mrna_cell),
                "values_redacted": True,
            },
            "command": redacted_command(),
        },
    }
    (args.output / "run.json").write_text(
        json.dumps(run, indent=2, sort_keys=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
