"""Top-level orchestration for the h5ad-first DECEPTICONx workflow.

This module intentionally contains no implementation of a reference builder or
deconvolution algorithm.  It connects the independently testable input,
reference, backend, and consensus stages and records enough metadata to audit a
run after its tabular outputs have been written.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields
from hashlib import sha256
from importlib import metadata
import json
import math
import platform
from pathlib import Path
import sys
from time import perf_counter
from typing import Any

import pandas as pd

from . import backends, consensus, io, references
from .exceptions import InputValidationError
from .models import BranchKey, DecepticonXResult, PipelineConfig, config_as_dict


_R_FILENAME_METHOD_ORDER = (
    "deconrnaseq",
    "epic",
    "cibersort_abs",
    "cibersort",
    "music",
)


def _configuration(config: PipelineConfig | Mapping[str, Any] | None) -> PipelineConfig:
    """Return a defensive, consistently typed pipeline configuration."""

    if config is None:
        result = PipelineConfig()
    elif isinstance(config, PipelineConfig):
        # PipelineConfig is mutable, so copy it before normalising container
        # fields.  A caller changing their object while a run is active must not
        # be able to alter the provenance or branch selection halfway through.
        result = PipelineConfig(**config_as_dict(config))
    elif isinstance(config, Mapping):
        allowed = {item.name for item in fields(PipelineConfig)}
        unknown = sorted(set(config) - allowed)
        if unknown:
            raise InputValidationError(
                "unknown pipeline configuration field(s): " + ", ".join(unknown)
            )
        result = PipelineConfig(**dict(config))
    else:
        raise TypeError("config must be a PipelineConfig, mapping, or None")

    # JSON configuration files naturally decode tuples as lists.  Canonicalise
    # them here so programmatic and CLI runs have identical provenance.
    def sequence(value: Any, name: str) -> tuple[str, ...]:
        if isinstance(value, str):
            items = value.split(",")
        else:
            try:
                items = list(value)
            except TypeError as exc:
                raise InputValidationError(f"{name} must be a sequence") from exc
        return tuple(str(item).strip() for item in items if str(item).strip())

    try:
        result.methods = sequence(result.methods, "methods")
        result.references = sequence(result.references, "references")
        result.exclude_cell_types = sequence(
            result.exclude_cell_types, "exclude_cell_types"
        )
    except TypeError as exc:
        raise InputValidationError(
            "methods, references, and exclude_cell_types must be sequences"
        ) from exc
    if result.epic_mrna_cell is not None:
        if not isinstance(result.epic_mrna_cell, Mapping):
            raise InputValidationError("epic_mrna_cell must be a JSON object/mapping")
        normalized_mapping: dict[str, float] = {}
        for key, value in result.epic_mrna_cell.items():
            if not isinstance(key, str) or not key or isinstance(value, bool):
                raise InputValidationError(
                    "epic_mrna_cell keys must be non-empty strings and values numeric"
                )
            try:
                numeric = float(value)
            except (TypeError, ValueError) as exc:
                raise InputValidationError(
                    f"epic_mrna_cell value for {key!r} must be numeric"
                ) from exc
            if not math.isfinite(numeric) or numeric <= 0:
                raise InputValidationError(
                    f"epic_mrna_cell value for {key!r} must be finite and positive"
                )
            normalized_mapping[key] = numeric
        result.epic_mrna_cell = normalized_mapping

    if not str(result.cell_type_key).strip():
        raise InputValidationError("cell_type_key cannot be empty")
    if result.sample_key is not None and not str(result.sample_key).strip():
        raise InputValidationError("sample_key cannot be empty")
    if result.gene_symbol_key is not None and not str(result.gene_symbol_key).strip():
        raise InputValidationError("gene_symbol_key cannot be empty")
    if not result.methods:
        raise InputValidationError("at least one deconvolution method is required")
    if not result.references:
        raise InputValidationError("at least one reference builder is required")
    species = str(result.species).strip().lower()
    if species in {"hs", "human", "homo sapiens", "homo_sapiens", "9606"}:
        # Reference builders expose the deliberately narrow first-release
        # spelling ``hs`` even though input validation accepts common aliases.
        result.species = "hs"
    for name in (
        "allow_normalized_x",
        "cibersort_qn",
        "strict_backends",
        "allow_partial_consensus",
    ):
        if type(getattr(result, name)) is not bool:
            raise InputValidationError(f"{name} must be a JSON/Python boolean")
    for name in ("threads", "cibersort_seed", "consensus_pairs"):
        value = getattr(result, name)
        if isinstance(value, bool) or not isinstance(value, int):
            raise InputValidationError(f"{name} must be an integer")
    if result.threads < 1:
        raise InputValidationError("threads must be at least 1")
    if result.consensus_pairs < 1:
        raise InputValidationError("consensus_pairs must be at least 1")

    engine_choices = {
        "cibersort_engine": {"rust", "libsvm", "numpy"},
        "epic_backend": {"auto", "rust", "python"},
        "epic_solver": {"auto", "nm", "nmf", "qp"},
        "deconrnaseq_backend": {"auto", "rust", "numpy"},
        "music_backend": {"auto", "rust", "numpy"},
        "consensus_mode": {"r_literal", "corrected"},
    }
    for name, choices in engine_choices.items():
        value = getattr(result, name)
        if not isinstance(value, str) or value.strip().lower() not in choices:
            raise InputValidationError(
                f"{name} must be one of: {', '.join(sorted(choices))}"
            )
        setattr(result, name, value.strip().lower())
    if result.consensus_mode == "r_literal" and result.consensus_pairs != 2:
        raise InputValidationError(
            "r_literal consensus fixes consensus_pairs at 2 to match optimal_id"
        )
    return result


def _bulk_input(bulk: str | Path | pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load or validate a bulk matrix and describe its origin."""

    if isinstance(bulk, pd.DataFrame):
        frame = io.validate_bulk_expression(bulk)
        origin: dict[str, Any] = {"kind": "dataframe", "path": None}
    elif isinstance(bulk, (str, Path)):
        source = Path(bulk)
        frame = io.load_bulk_expression(source)
        origin = {"kind": "path", "path": str(source.resolve())}
    else:
        raise TypeError("bulk must be a pandas DataFrame or a CSV/TSV path")

    values = frame.to_numpy(copy=False)
    origin.update(
        {
            "n_genes": int(frame.shape[0]),
            "n_samples": int(frame.shape[1]),
            "minimum": float(values.min()),
            "maximum": float(values.max()),
            "sample_names": [str(value) for value in frame.columns],
        }
    )
    return frame, origin


def _reference_result(value: Any) -> tuple[dict[str, pd.DataFrame], Any]:
    """Accept the stable reference return value plus early-development forms."""

    if isinstance(value, tuple) and len(value) == 2:
        signatures, diagnostics = value
    elif hasattr(value, "signatures"):
        signatures = value.signatures
        diagnostics = getattr(value, "diagnostics", {})
    else:
        signatures, diagnostics = value, {}
    if not isinstance(signatures, Mapping) or not signatures:
        raise InputValidationError("reference builders produced no signature matrices")
    return dict(signatures), diagnostics


def _backend_result(value: Any) -> tuple[dict[Any, pd.DataFrame], Any]:
    """Accept the stable backend return value plus a result-object variant."""

    if isinstance(value, tuple) and len(value) == 2:
        estimates, diagnostics = value
    elif hasattr(value, "estimates"):
        estimates = value.estimates
        diagnostics = getattr(value, "diagnostics", {})
    else:
        estimates, diagnostics = value, {}
    if not isinstance(estimates, Mapping):
        raise InputValidationError("deconvolution backends returned an invalid result")
    return dict(estimates), diagnostics


def _consensus_result(
    value: Any,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Any]:
    """Extract primary, unclosed, and closed consensus tables."""

    if all(hasattr(value, name) for name in ("primary", "unclosed", "closed")):
        primary = value.primary
        unclosed = value.unclosed
        closed = value.closed
        diagnostics = getattr(value, "diagnostics", {})
    elif hasattr(value, "normalized") and hasattr(value, "raw"):
        # Compatibility with the initial alpha result object.
        closed = value.normalized
        unclosed = value.raw
        primary = closed
        diagnostics = getattr(value, "diagnostics", {})
    elif isinstance(value, Mapping):
        unclosed = value.get(
            "unclosed", value.get("raw", value.get("consensus_raw"))
        )
        closed = value.get(
            "closed", value.get("normalized", value.get("consensus_normalized"))
        )
        primary = value.get("primary", value.get("consensus", closed))
        diagnostics = value.get("diagnostics", {})
    elif isinstance(value, tuple) and len(value) in (2, 3, 4):
        # Compatibility order: unclosed, closed, optional primary/diagnostics.
        unclosed, closed = value[:2]
        if len(value) == 4:
            primary, diagnostics = value[2:]
        else:
            primary = closed
            diagnostics = value[2] if len(value) == 3 else {}
    else:
        raise InputValidationError("consensus builder returned an invalid result")
    if not all(
        isinstance(frame, pd.DataFrame)
        for frame in (primary, unclosed, closed)
    ):
        raise InputValidationError("consensus builder did not return DataFrame outputs")
    return primary, unclosed, closed, diagnostics


def _cell_types(signatures: Mapping[str, pd.DataFrame]) -> tuple[str, ...]:
    """Require one identifiable cell-type schema across every reference."""

    frames = list(signatures.values())
    for signature in frames:
        if not isinstance(signature, pd.DataFrame):
            raise InputValidationError("every signature matrix must be a pandas DataFrame")
        if signature.columns.has_duplicates:
            raise InputValidationError("signature matrices contain duplicate cell types")
    if not frames:
        raise InputValidationError("reference builders produced no signature matrices")
    ordered = tuple(str(value) for value in frames[0].columns)
    expected = set(ordered)
    if not ordered:
        raise InputValidationError("signature matrices contain no cell types")
    for name, signature in signatures.items():
        actual = {str(value) for value in signature.columns}
        if actual != expected:
            missing = sorted(expected.difference(actual))
            extra = sorted(actual.difference(expected))
            raise InputValidationError(
                f"signature {name!r} has a different cell-type schema "
                f"(missing={missing}, extra={extra}); missing types cannot be "
                "interpreted as zero estimates"
            )
    return ordered


def _strategy_order(
    estimates: Mapping[BranchKey, pd.DataFrame],
    references: tuple[str, ...],
) -> tuple[BranchKey, ...]:
    """Return a stable C-locale equivalent of the R result-file order."""

    available = set(estimates)
    ordered = [
        BranchKey(method=method, reference=reference)
        for method in _R_FILENAME_METHOD_ORDER
        for reference in references
        if BranchKey(method=method, reference=reference) in available
    ]
    remainder = sorted(
        available.difference(ordered), key=lambda key: (key.method, key.reference)
    )
    return tuple((*ordered, *remainder))


def _package_version() -> str:
    try:
        return metadata.version("decepticonx-py-rust")
    except metadata.PackageNotFoundError:
        return "0.1.0"


def _component_versions() -> dict[str, str | None]:
    """Return installed distribution versions without importing backends."""

    distributions = (
        "decepticon-fast",
        "python-cibersort-rs",
        "epic_py",
        "deconrnaseq-py",
        "deconrnaseq-rust",
        "music-py",
        "music-py-rust",
    )
    versions: dict[str, str | None] = {}
    for distribution in distributions:
        try:
            versions[distribution] = metadata.version(distribution)
        except metadata.PackageNotFoundError:
            versions[distribution] = None
    return versions


def _provenance_config(config: PipelineConfig) -> dict[str, Any]:
    """Redact EPIC mRNA values while retaining a reproducibility fingerprint."""

    result = config_as_dict(config)
    mapping = result.get("epic_mrna_cell")
    if mapping is not None:
        encoded = "\n".join(
            f"{key}\t{float(mapping[key]):.17g}" for key in sorted(mapping)
        )
        result["epic_mrna_cell"] = {
            "provided": True,
            "entry_count": len(mapping),
            "sha256": sha256(encoded.encode("utf-8")).hexdigest(),
            "hash_format": "sorted key, tab, %.17g value, newline-separated",
            "values_redacted": True,
        }
    return result


def _redacted_command(arguments: list[str]) -> list[str]:
    """Remove inline EPIC mapping values from recorded command arguments."""

    sensitive = {
        "--epic-mrna-cell",
        "--epic-mrna-json",
        "--epic-mrna-cell-json",
    }
    output: list[str] = []
    redact_next = False
    for argument in arguments:
        if redact_next:
            output.append("<redacted>")
            redact_next = False
            continue
        if argument in sensitive:
            output.append(argument)
            redact_next = True
            continue
        matched = next(
            (flag for flag in sensitive if argument.startswith(flag + "=")),
            None,
        )
        output.append(f"{matched}=<redacted>" if matched else argument)
    return output


def run_decepticonx(
    single_cell: Any,
    bulk: str | Path | pd.DataFrame,
    config: PipelineConfig | Mapping[str, Any] | None = None,
    output_dir: str | Path | None = None,
) -> DecepticonXResult:
    """Run reference construction, five-method deconvolution, and consensus.

    Parameters
    ----------
    single_cell:
        A path to an h5ad file or an in-memory :class:`anndata.AnnData`.
    bulk:
        A strict genes-by-samples CSV/TSV path or an in-memory DataFrame.
    config:
        A :class:`PipelineConfig`, a JSON-like mapping, or ``None`` for defaults.
    output_dir:
        When supplied, write the complete stable result tree before returning.
    """

    started = perf_counter()
    cfg = _configuration(config)
    timings: dict[str, float] = {}

    if any(str(method).strip().lower() == "epic" for method in cfg.methods):
        if cfg.epic_mrna_cell is None:
            raise InputValidationError(
                "EPIC was requested but epic_mrna_cell was not supplied. "
                "Provide an authorized mapping explicitly; no licensed defaults "
                "are bundled or substituted."
            )

    stage = perf_counter()
    bulk_frame, bulk_provenance = _bulk_input(bulk)
    timings["bulk_input"] = perf_counter() - stage

    stage = perf_counter()
    prepared = io.prepare_single_cell(
        single_cell,
        bulk_frame,
        cell_type_key=cfg.cell_type_key,
        sample_key=cfg.sample_key,
        counts_layer=cfg.counts_layer,
        gene_symbol_key=cfg.gene_symbol_key,
        exclude_cell_types=cfg.exclude_cell_types,
        allow_normalized_x=cfg.allow_normalized_x,
        species=cfg.species,
    )
    timings["single_cell_input"] = perf_counter() - stage

    stage = perf_counter()
    signatures, reference_diagnostics = _reference_result(
        references.build_references(
            prepared,
            bulk_frame,
            references=cfg.references,
            species=cfg.species,
        )
    )
    timings["references"] = perf_counter() - stage

    cell_types = _cell_types(signatures)

    stage = perf_counter()
    estimates, backend_diagnostics = _backend_result(
        backends.run_backends(
            signatures,
            bulk_frame,
            methods=cfg.methods,
            qn=cfg.cibersort_qn,
            seed=cfg.cibersort_seed,
            threads=cfg.threads,
            epic_mrna_cell=cfg.epic_mrna_cell,
            strict=cfg.strict_backends,
            cibersort_engine=cfg.cibersort_engine,
            deconrnaseq_backend=cfg.deconrnaseq_backend,
            epic_backend=cfg.epic_backend,
            epic_solver=cfg.epic_solver,
            music_backend=cfg.music_backend,
        )
    )
    timings["backends"] = perf_counter() - stage

    skipped = dict(backend_diagnostics.get("skipped", {})) if isinstance(
        backend_diagnostics, Mapping
    ) else {}
    if skipped and not cfg.allow_partial_consensus:
        raise InputValidationError(
            "one or more requested branches were skipped; set "
            "allow_partial_consensus=True only when a partial consensus is "
            "scientifically intentional (skipped: "
            + ", ".join(sorted(skipped))
            + ")"
        )

    stage = perf_counter()
    primary, unclosed, closed, consensus_diagnostics = _consensus_result(
        consensus.build_consensus(
            estimates,
            cell_types=cell_types,
            n_pairs=cfg.consensus_pairs,
            mode=cfg.consensus_mode,
            # Reference builders return their canonical names in requested
            # order.  Use that realised order so case-normalised configuration
            # spellings cannot perturb the R filename-equivalent tie break.
            strategy_order=_strategy_order(estimates, tuple(signatures)),
        )
    )
    timings["consensus"] = perf_counter() - stage
    timings["total"] = perf_counter() - started

    diagnostics = {
        "input": {
            "warnings": list(getattr(prepared, "warnings", ())),
        },
        "references": reference_diagnostics,
        "backends": backend_diagnostics,
        "consensus": consensus_diagnostics,
    }
    provenance = {
        "package": "decepticonx-py-rust",
        "version": _package_version(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "command": _redacted_command(list(sys.argv)),
        "config": _provenance_config(cfg),
        "component_versions": _component_versions(),
        "bulk": bulk_provenance,
        "single_cell": dict(getattr(prepared, "provenance", {})),
        "cell_types": list(cell_types),
        "reference_names": list(signatures),
        "completed_branches": [
            getattr(key, "slug", str(key)) for key in estimates
        ],
    }

    result = DecepticonXResult(
        signatures=signatures,
        estimates=estimates,
        consensus=primary,
        consensus_unclosed=unclosed,
        consensus_closed=closed,
        diagnostics=diagnostics,
        timings=timings,
        provenance=provenance,
    )
    if output_dir is not None:
        result.write(output_dir)
    return result


__all__ = ["run_decepticonx"]
