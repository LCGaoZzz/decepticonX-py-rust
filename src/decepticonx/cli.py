"""Command-line interface for :mod:`decepticonx`.

The CLI is deliberately a thin layer over the Python API.  In particular,
``validate`` stops after the two input-contract checks and never imports or
executes an optional deconvolution backend.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import sys
from typing import Any

import pandas as pd

from .backends import probe_backends, probe_runtime
from .exceptions import DecepticonXError, InputValidationError
from .io import load_bulk_expression, prepare_single_cell
from .models import PipelineConfig, config_as_dict
from .pipeline import _configuration, run_decepticonx


def _comma_values(values: Sequence[Any] | None) -> tuple[str, ...] | None:
    """Expand repeatable and comma-separated CLI values."""

    if values is None:
        return None
    flat: list[str] = []
    for value in values:
        nested = value if isinstance(value, (list, tuple)) else (value,)
        for token in nested:
            flat.extend(item.strip() for item in str(token).split(",") if item.strip())
    result = tuple(flat)
    if not result:
        raise InputValidationError("a list option was supplied without any values")
    return result


def _json_object(value: str | Path, *, label: str) -> dict[str, Any]:
    """Read a JSON object from a file path or directly from the command line."""

    token = str(value)
    candidate = Path(token)
    try:
        if candidate.is_file():
            text = candidate.read_text(encoding="utf-8")
        else:
            text = token
        decoded = json.loads(text)
    except (OSError, json.JSONDecodeError) as exc:
        raise InputValidationError(f"{label} must be a JSON object or JSON file: {exc}") from exc
    if not isinstance(decoded, dict):
        raise InputValidationError(f"{label} must decode to a JSON object")
    return decoded


def _add_input_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--cell-type-key",
        default=None,
        help="adata.obs column containing cell-type labels (default: cell_type)",
    )
    parser.add_argument(
        "--sample-key",
        default=None,
        help="optional adata.obs donor/sample column",
    )
    parser.add_argument(
        "--counts-layer",
        default=None,
        help="count source: auto, X, raw.X, or an AnnData layer name",
    )
    parser.add_argument(
        "--gene-symbol-key",
        default=None,
        help="adata.var column mapping source identifiers to human gene symbols",
    )
    parser.add_argument(
        "--exclude-cell-type",
        "--exclude-cell-types",
        action="append",
        nargs="+",
        default=None,
        metavar="TYPE[,TYPE...]",
        help="cell type(s) to omit; repeat or provide a comma-separated list",
    )
    normalized = parser.add_mutually_exclusive_group()
    normalized.add_argument(
        "--allow-normalized-x",
        action="store_true",
        dest="allow_normalized_x",
        default=None,
        help="explicitly permit a non-integer adata.X compatibility run",
    )
    normalized.add_argument(
        "--require-raw-counts",
        action="store_false",
        dest="allow_normalized_x",
        help="require count-like single-cell input (the default)",
    )
    parser.add_argument(
        "--species",
        default=None,
        help="species identifier; human aliases are canonicalized to hs",
    )


def _add_run_options(parser: argparse.ArgumentParser) -> None:
    _add_input_options(parser)
    parser.add_argument(
        "--method",
        "--methods",
        action="append",
        nargs="+",
        default=None,
        metavar="NAME[,NAME...]",
        help=(
            "method(s) to run; default: cibersort,deconrnaseq,music; "
            "cibersort_abs and epic remain explicit opt-ins"
        ),
    )
    parser.add_argument(
        "--reference",
        "--references",
        action="append",
        nargs="+",
        default=None,
        metavar="NAME[,NAME...]",
        help=(
            "reference builder(s); default: bayesprism,music2; "
            "monocle3 remains an explicit opt-in"
        ),
    )
    qn = parser.add_mutually_exclusive_group()
    qn.add_argument(
        "--cibersort-qn",
        "--qn",
        action="store_true",
        dest="cibersort_qn",
        default=None,
        help="enable CIBERSORT quantile normalization (default)",
    )
    qn.add_argument(
        "--no-cibersort-qn",
        "--no-qn",
        action="store_false",
        dest="cibersort_qn",
        help="disable CIBERSORT quantile normalization",
    )
    parser.add_argument("--cibersort-seed", "--seed", type=int, default=None)
    parser.add_argument(
        "--cibersort-engine",
        choices=("rust", "libsvm", "numpy"),
        default=None,
        help="CIBERSORT solver engine (default: rust)",
    )
    parser.add_argument(
        "--epic-backend",
        choices=("auto", "rust", "python"),
        default=None,
        help="EPIC backend; use rust to require its native kernel",
    )
    parser.add_argument(
        "--epic-solver",
        choices=("auto", "nm", "nmf", "qp"),
        default=None,
        help="EPIC optimizer (nm gives strict original-R optimizer parity)",
    )
    parser.add_argument(
        "--deconrnaseq-backend",
        choices=("auto", "rust", "numpy"),
        default=None,
        help="DeconRNASeq backend; use rust to require its native extension",
    )
    parser.add_argument(
        "--music-backend",
        choices=("auto", "rust", "numpy"),
        default=None,
        help="MuSiC backend; use rust to require its native extension",
    )
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--consensus-pairs", "--n-pairs", type=int, default=None)
    parser.add_argument(
        "--consensus-mode",
        choices=("r_literal", "corrected"),
        default=None,
        help=(
            "consensus selector: R positional selection (default) or the "
            "structure-aware corrected selector"
        ),
    )
    strict = parser.add_mutually_exclusive_group()
    strict.add_argument(
        "--strict-backends",
        "--strict",
        action="store_true",
        dest="strict_backends",
        default=None,
        help="fail the run when any requested backend branch fails (default)",
    )
    strict.add_argument(
        "--permissive",
        action="store_false",
        dest="strict_backends",
        help="skip unavailable/failed branches and record them in diagnostics",
    )
    partial = parser.add_mutually_exclusive_group()
    partial.add_argument(
        "--allow-partial-consensus",
        action="store_true",
        dest="allow_partial_consensus",
        default=None,
        help="allow a consensus after permissive mode skips requested branches",
    )
    partial.add_argument(
        "--require-complete-consensus",
        action="store_false",
        dest="allow_partial_consensus",
        help="reject a partial branch set (the default)",
    )
    parser.add_argument(
        "--epic-mrna-cell",
        "--epic-mrna-json",
        "--epic-mrna-cell-json",
        default=None,
        metavar="JSON_OR_FILE",
        help="authorized EPIC mRNA_cell JSON mapping or path to a JSON file",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="decepticonx",
        description="R-free, h5ad-first DECEPTICONx orchestration",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")
    commands = parser.add_subparsers(dest="command", required=True)

    run = commands.add_parser("run", help="run the complete DECEPTICONx workflow")
    run.add_argument("single_cell", help="cells-by-genes .h5ad input")
    run.add_argument("bulk", help="genes-by-samples .csv/.tsv/.txt input")
    run.add_argument(
        "-o",
        "--output-dir",
        "--output",
        required=True,
        help="directory for signatures, estimates, consensus, and run metadata",
    )
    run.add_argument(
        "--config",
        default=None,
        metavar="FILE",
        help="JSON PipelineConfig object; explicit flags override its fields",
    )
    _add_run_options(run)
    run.set_defaults(handler=_run_command)

    validate = commands.add_parser(
        "validate", help="validate h5ad and bulk inputs without running algorithms"
    )
    validate.add_argument("single_cell", help="cells-by-genes .h5ad input")
    validate.add_argument("bulk", help="genes-by-samples .csv/.tsv/.txt input")
    validate.add_argument(
        "--config",
        default=None,
        metavar="FILE",
        help="JSON PipelineConfig object; input flags override its fields",
    )
    _add_input_options(validate)
    validate.set_defaults(handler=_validate_command)

    check = commands.add_parser(
        "check-backends", help="report whether each optional backend is importable"
    )
    check.add_argument(
        "--require-all",
        action="store_true",
        help="return a failing status when one or more backends are unavailable",
    )
    check.add_argument(
        "--require-native",
        action="store_true",
        help="return a failing status unless every component has Rust available",
    )
    check.set_defaults(handler=_check_backends_command)
    return parser


def _config_from_arguments(args: argparse.Namespace, *, run: bool) -> PipelineConfig:
    values: dict[str, Any] = config_as_dict(PipelineConfig())
    if args.config is not None:
        values.update(_json_object(args.config, label="--config"))

    names = [
        "cell_type_key",
        "sample_key",
        "counts_layer",
        "gene_symbol_key",
        "allow_normalized_x",
        "species",
    ]
    for name in names:
        value = getattr(args, name, None)
        if value is not None:
            values[name] = value

    exclusions = _comma_values(args.exclude_cell_type)
    if exclusions is not None:
        values["exclude_cell_types"] = exclusions

    if run:
        methods = _comma_values(args.method)
        reference_names = _comma_values(args.reference)
        if methods is not None:
            values["methods"] = methods
        if reference_names is not None:
            values["references"] = reference_names
        for name in (
            "cibersort_qn",
            "cibersort_seed",
            "threads",
            "consensus_pairs",
            "consensus_mode",
            "strict_backends",
            "allow_partial_consensus",
            "cibersort_engine",
            "epic_backend",
            "epic_solver",
            "deconrnaseq_backend",
            "music_backend",
        ):
            value = getattr(args, name, None)
            if value is not None:
                values[name] = value
        if args.epic_mrna_cell is not None:
            values["epic_mrna_cell"] = _json_object(
                args.epic_mrna_cell, label="--epic-mrna-cell"
            )
    return _configuration(values)


def _run_command(args: argparse.Namespace) -> int:
    config = _config_from_arguments(args, run=True)
    result = run_decepticonx(
        args.single_cell,
        args.bulk,
        config=config,
        output_dir=args.output_dir,
    )
    summary = {
        "ok": True,
        "output_dir": str(Path(args.output_dir).resolve()),
        "references": list(result.signatures),
        "branches": len(result.estimates),
        "samples": int(result.consensus.shape[0]),
        "cell_types": int(result.consensus.shape[1]),
        "total_seconds": result.timings.get("total"),
    }
    print(json.dumps(summary, ensure_ascii=False, default=str))
    return 0


def _validate_command(args: argparse.Namespace) -> int:
    config = _config_from_arguments(args, run=False)
    bulk = load_bulk_expression(args.bulk)
    prepared = prepare_single_cell(
        args.single_cell,
        bulk,
        cell_type_key=config.cell_type_key,
        sample_key=config.sample_key,
        counts_layer=config.counts_layer,
        gene_symbol_key=config.gene_symbol_key,
        exclude_cell_types=config.exclude_cell_types,
        allow_normalized_x=config.allow_normalized_x,
        species=config.species,
    )
    report = {
        "ok": True,
        "bulk": {
            "n_genes": int(bulk.shape[0]),
            "n_samples": int(bulk.shape[1]),
        },
        "single_cell": dict(getattr(prepared, "provenance", {})),
        "warnings": list(getattr(prepared, "warnings", ())),
    }
    print(json.dumps(report, ensure_ascii=False, default=str))
    return 0


def _check_backends_command(args: argparse.Namespace) -> int:
    status = probe_backends()
    runtime = probe_runtime()
    all_available = all(status.values()) and runtime["decepticon_fast"]["available"]
    all_native = all(item["native_available"] for item in runtime.values())
    report = {
        "ok": all_available,
        "backends": status,
        "runtime": runtime,
        "all_native": all_native,
    }
    print(json.dumps(report, ensure_ascii=False, sort_keys=True))
    if args.require_all and not all_available:
        return 1
    if args.require_native and not all_native:
        return 1
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Parse ``argv``, execute one command, and return a process status."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (DecepticonXError, OSError, TypeError, ValueError) as exc:
        print(f"decepticonx: error: {exc}", file=sys.stderr)
        return 2


__all__ = ["build_parser", "main"]
