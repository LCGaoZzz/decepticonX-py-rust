"""Lazy adapters for the optional deconvolution implementations.

The adapters in this module deliberately contain no fallback algorithms and
no reference data.  Each implementation remains an independently installed
optional dependency; this package only reconciles its input/output contract.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
import importlib
import os
from typing import Any
import warnings

import numpy as np
import pandas as pd

from .exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    InputValidationError,
)
from .models import BranchKey


_BACKEND_MODULES = {
    "cibersort": "python_cibersort",
    "cibersort_abs": "python_cibersort",
    "epic": "epic_py",
    "deconrnaseq": "deconrnaseq",
    "music": "music_py",
}

_REQUIRED_ATTRIBUTES = {
    "python_cibersort": ("cibersort_all",),
    "epic_py": ("EpicReference", "EPIC"),
    "deconrnaseq": ("deconrnaseq",),
    "music_py": ("SCReference", "music_prop"),
}

_INSTALL_HINTS = {
    "python_cibersort": (
        'python -m pip install "git+https://github.com/LCGaoZzz/'
        'python-cibersort-rs.git@cd957af06357d2636ba18af1135bf70b488dfbbb"'
    ),
    "epic_py": (
        "install an authorized epic-py-rust checkout at "
        "2e661cccb5c9d7c8327bd65253e3495a1d3c1991 (see this project's README)"
    ),
    "deconrnaseq": (
        'python -m pip install "git+https://github.com/LCGaoZzz/'
        'deconrnaseq-py.git@a2ec1b9ebe341628e945a6d25c6068cad975f6b5"'
    ),
    "music_py": (
        'python -m pip install "git+https://github.com/LCGaoZzz/'
        'music-py.git@a4464c84f87fbeb61f0fa7032625bad21add52dd"'
    ),
}


def _json_scalar(value: Any) -> Any:
    """Convert a pandas/NumPy scalar into strict JSON-friendly data."""

    if value is None:
        return None
    if isinstance(value, np.generic):
        value = value.item()
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _frame_payload(frame: pd.DataFrame) -> dict[str, Any]:
    """Represent a labelled diagnostic frame without losing sample IDs."""

    columns = [str(column) for column in frame.columns]
    records: list[dict[str, Any]] = []
    for sample, row in frame.iterrows():
        record = {"sample": str(sample)}
        record.update(
            {
                str(column): _json_scalar(row[column])
                for column in frame.columns
            }
        )
        records.append(record)
    return {"columns": columns, "records": records}


def _native_available(method: str, module: Any | None = None) -> bool:
    """Probe the optional native component used by one backend."""

    try:
        if method in {"cibersort", "cibersort_abs"}:
            importlib.import_module("python_cibersort._native")
            return True
        if method == "epic":
            rust_module = importlib.import_module("epic_py._rust")
            return bool(rust_module.rust_available())
        if method == "deconrnaseq":
            backend = module if module is not None else importlib.import_module("deconrnaseq")
            probe = getattr(backend, "rust_available", None)
            return bool(callable(probe) and probe())
        if method == "music":
            backend = module if module is not None else importlib.import_module("music_py")
            probe = getattr(backend, "rust_available", None)
            return bool(callable(probe) and probe())
    except (ImportError, OSError, RuntimeError):
        return False
    return False


def probe_runtime() -> dict[str, dict[str, Any]]:
    """Report package and native-engine availability, including references."""

    status: dict[str, dict[str, Any]] = {}
    try:
        fast = importlib.import_module("decepticon_fast")
        required = ("bayesprism_base", "monocle3_base", "music2_base")
        available = all(callable(getattr(fast, name, None)) for name in required)
        engine = fast.get_backend() if callable(getattr(fast, "get_backend", None)) else None
    except (ImportError, OSError):
        available, engine = False, None
    status["decepticon_fast"] = {
        "available": available,
        "engine": engine,
        "native_available": available and engine == "rust",
    }

    for method in ("cibersort", "cibersort_abs", "epic", "deconrnaseq", "music"):
        module_name = _BACKEND_MODULES[method]
        try:
            module = importlib.import_module(module_name)
            available = all(
                callable(getattr(module, name, None))
                for name in _REQUIRED_ATTRIBUTES[module_name]
            )
        except (ImportError, OSError):
            module, available = None, False
        status[method] = {
            "available": available,
            "native_available": available and _native_available(method, module),
        }
    return status

_METHOD_ALIASES = {
    "cibersort": "cibersort",
    "cibersort_abs": "cibersort_abs",
    "cibersortabs": "cibersort_abs",
    "epic": "epic",
    "deconrnaseq": "deconrnaseq",
    "decon_rna_seq": "deconrnaseq",
    "music": "music",
}


def _load_backend(method: str) -> Any:
    module_name = _BACKEND_MODULES[method]
    try:
        module = importlib.import_module(module_name)
    except (ImportError, OSError) as exc:
        raise BackendUnavailableError(
            f"The optional {method!r} backend ({module_name}) is unavailable. "
            f"Install it with: {_INSTALL_HINTS[module_name]}"
        ) from exc

    missing = [
        name for name in _REQUIRED_ATTRIBUTES[module_name]
        if not callable(getattr(module, name, None))
    ]
    if missing:
        raise BackendUnavailableError(
            f"The installed {module_name!r} package does not expose the required "
            f"API: {', '.join(missing)}. Upgrade/reinstall it with: "
            f"{_INSTALL_HINTS[module_name]}"
        )
    return module


def probe_backends() -> dict[str, bool]:
    """Return installation/API availability for all five method branches.

    Import attempts happen only when this explicit probe is called.  CIBERSORT
    and CIBERSORT-ABS share one package, so their availability is identical.
    """

    by_module: dict[str, bool] = {}
    for module_name, required in _REQUIRED_ATTRIBUTES.items():
        try:
            module = importlib.import_module(module_name)
            by_module[module_name] = all(
                callable(getattr(module, name, None)) for name in required
            )
        except (ImportError, OSError):
            by_module[module_name] = False
    return {
        method: by_module[module_name]
        for method, module_name in _BACKEND_MODULES.items()
    }


def _matrix(value: pd.DataFrame, *, name: str) -> pd.DataFrame:
    if not isinstance(value, pd.DataFrame):
        raise InputValidationError(f"{name} must be a pandas DataFrame")
    if value.empty or value.shape[1] == 0:
        raise InputValidationError(f"{name} must be a non-empty two-dimensional matrix")

    frame = value.copy()
    frame.index = pd.Index([str(item) for item in frame.index])
    frame.columns = pd.Index([str(item) for item in frame.columns])
    if frame.index.has_duplicates:
        raise InputValidationError(f"{name} has duplicate row identifiers")
    if frame.columns.has_duplicates:
        raise InputValidationError(f"{name} has duplicate column identifiers")
    try:
        frame = frame.astype(np.float64)
    except (TypeError, ValueError) as exc:
        raise InputValidationError(f"{name} contains non-numeric values") from exc
    values = frame.to_numpy(copy=False)
    if not np.isfinite(values).all():
        raise InputValidationError(f"{name} contains NaN or infinite values")
    if (values < 0).any():
        raise InputValidationError(f"{name} contains negative expression values")
    return frame


def _prepare_inputs(
    signature: pd.DataFrame,
    bulk: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sig = _matrix(signature, name="signature")
    mix = _matrix(bulk, name="bulk")
    if not set(sig.index).intersection(mix.index):
        raise InputValidationError("signature and bulk have no genes in common")
    return sig, mix


def _finish_output(
    value: Any,
    *,
    samples: Sequence[str],
    cell_types: Sequence[str],
    method: str,
) -> pd.DataFrame:
    if not isinstance(value, pd.DataFrame):
        try:
            value = pd.DataFrame(value)
        except Exception as exc:  # pragma: no cover - defensive package boundary
            raise BackendExecutionError(
                f"{method} returned a non-tabular estimate"
            ) from exc

    frame = value.copy()
    frame.index = pd.Index([str(item) for item in frame.index])
    frame.columns = pd.Index([str(item) for item in frame.columns])
    expected_samples = [str(item) for item in samples]
    expected_types = [str(item) for item in cell_types]

    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise BackendExecutionError(f"{method} returned duplicate labels")
    if set(frame.index) != set(expected_samples):
        missing = sorted(set(expected_samples).difference(frame.index))
        extra = sorted(set(frame.index).difference(expected_samples))
        raise BackendExecutionError(
            f"{method} sample labels do not match the bulk matrix "
            f"(missing={missing}, extra={extra})"
        )
    if set(frame.columns) != set(expected_types):
        missing = sorted(set(expected_types).difference(frame.columns))
        extra = sorted(set(frame.columns).difference(expected_types))
        raise BackendExecutionError(
            f"{method} cell-type labels do not match the signature matrix "
            f"(missing={missing}, extra={extra})"
        )

    frame = frame.reindex(index=expected_samples, columns=expected_types)
    try:
        frame = frame.astype(np.float64)
    except (TypeError, ValueError) as exc:
        raise BackendExecutionError(f"{method} returned non-numeric estimates") from exc
    values = frame.to_numpy(copy=False)
    if not np.isfinite(values).all():
        raise BackendExecutionError(f"{method} returned NaN or infinite estimates")
    if (values < -1e-12).any():
        minimum = float(values.min())
        raise BackendExecutionError(
            f"{method} returned negative estimates (minimum={minimum:.6g})"
        )
    # Solvers can leave harmless round-off immediately below zero.  The public
    # contract is non-negative, so clamp only values within the checked bound.
    if (values < 0).any():
        frame = frame.clip(lower=0.0)
    return frame


def _backend_failure(method: str, exc: BaseException) -> BaseException:
    if isinstance(exc, (BackendUnavailableError, BackendExecutionError)):
        return exc
    module_name = _BACKEND_MODULES[method]
    if isinstance(exc, (ImportError, OSError)):
        return BackendUnavailableError(
            f"The {method} backend could not load a required native component. "
            f"Install/rebuild it with: {_INSTALL_HINTS[module_name]}"
        )
    return BackendExecutionError(f"{method} backend execution failed: {exc}")


def _call_backend(method: str, operation: Callable[[], pd.DataFrame]) -> pd.DataFrame:
    try:
        return operation()
    except Exception as exc:
        error = _backend_failure(method, exc)
        if error is exc:
            raise
        raise error from exc


def _result_frame(result: Any, *, attribute: str, key: str | None = None) -> Any:
    value = getattr(result, attribute, None)
    if value is not None:
        return value
    if key is not None:
        try:
            return result[key]
        except (KeyError, TypeError):
            pass
    raise BackendExecutionError(
        f"backend result has neither attribute {attribute!r} nor key {key!r}"
    )


def _cibersort_table(result: Any, *, method: str) -> pd.DataFrame:
    table = result if isinstance(result, pd.DataFrame) else getattr(result, "table", None)
    if not isinstance(table, pd.DataFrame):
        raise BackendExecutionError(f"{method} did not return a result table")
    table = table.copy()
    diagnostics = {"P-value", "Correlation", "RMSE"}
    diagnostics.update(
        column for column in table.columns
        if str(column).startswith("Absolute score (")
    )
    diagnostic_columns = [
        column for column in table.columns if column in diagnostics
    ]
    estimate = table.drop(columns=diagnostic_columns)
    if diagnostic_columns:
        estimate.attrs["decepticonx_diagnostics"] = {
            "fit": _frame_payload(table.loc[:, diagnostic_columns])
        }
    return estimate


def run_cibersort_pair(
    signature: pd.DataFrame,
    bulk: pd.DataFrame,
    *,
    qn: bool = True,
    seed: int = 0,
    threads: int = 1,
    engine: str = "rust",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run relative and ``sig.score`` absolute CIBERSORT from one fit batch."""

    sig, mix = _prepare_inputs(signature, bulk)
    if threads < 1:
        raise InputValidationError("threads must be at least 1")
    engine = str(engine).strip().lower()
    if engine not in {"rust", "libsvm", "numpy"}:
        raise InputValidationError("CIBERSORT engine must be rust, libsvm, or numpy")

    try:
        backend = _load_backend("cibersort")
        outputs = backend.cibersort_all(
            sig,
            mix,
            perm=0,
            QN=bool(qn),
            seed=int(seed),
            threads=int(threads),
            engine=engine,
        )
        relative_raw = _cibersort_table(outputs["relative"], method="cibersort")
        absolute_raw = _cibersort_table(outputs["sig.score"], method="cibersort_abs")
        relative_details = relative_raw.attrs.get("decepticonx_diagnostics")
        absolute_details = absolute_raw.attrs.get("decepticonx_diagnostics")
        relative = _finish_output(
            relative_raw,
            samples=mix.columns,
            cell_types=sig.columns,
            method="cibersort",
        )
        absolute = _finish_output(
            absolute_raw,
            samples=mix.columns,
            cell_types=sig.columns,
            method="cibersort_abs",
        )
        relative.attrs["decepticonx_engine"] = engine
        absolute.attrs["decepticonx_engine"] = engine
        if relative_details:
            relative.attrs["decepticonx_diagnostics"] = relative_details
        if absolute_details:
            absolute.attrs["decepticonx_diagnostics"] = absolute_details
        return relative, absolute
    except Exception as exc:
        error = _backend_failure("cibersort", exc)
        if error is exc:
            raise
        raise error from exc


def _validated_mrna_mapping(
    mapping: Mapping[str, float] | None,
    cell_types: Sequence[str],
) -> dict[str, float]:
    if not mapping:
        raise BackendUnavailableError(
            "EPIC with a custom signature requires an explicit mRNA_cell mapping "
            "via PipelineConfig.epic_mrna_cell. No EPIC reference or mRNA data are "
            "bundled or loaded automatically because the upstream data are licensed."
        )
    try:
        values = {str(key): float(value) for key, value in mapping.items()}
    except (TypeError, ValueError) as exc:
        raise InputValidationError("EPIC mRNA_cell values must be numeric") from exc
    invalid = [key for key, value in values.items() if not np.isfinite(value) or value <= 0]
    if invalid:
        raise InputValidationError(
            "EPIC mRNA_cell values must be finite and positive; invalid keys: "
            + ", ".join(invalid)
        )
    if "default" not in values:
        required = [str(cell_type) for cell_type in cell_types] + ["otherCells"]
        missing = [cell_type for cell_type in required if cell_type not in values]
        if missing:
            raise BackendUnavailableError(
                "EPIC mRNA_cell has no value for: "
                + ", ".join(missing)
                + ". Supply every custom cell type or an explicit 'default' value; "
                "licensed defaults are never loaded automatically."
            )
    return values


def run_epic(
    signature: pd.DataFrame,
    bulk: pd.DataFrame,
    *,
    mrna_cell: Mapping[str, float] | None,
    threads: int = 1,
    backend_name: str = "auto",
    solver: str = "auto",
) -> pd.DataFrame:
    """Run EPIC with an in-memory custom reference and explicit mRNA values."""

    sig, mix = _prepare_inputs(signature, bulk)
    mrna_values = _validated_mrna_mapping(mrna_cell, sig.columns)
    if threads < 1:
        raise InputValidationError("threads must be at least 1")
    solver_name = str(solver).strip().lower()
    if solver_name not in {"auto", "nm", "nmf", "qp"}:
        raise InputValidationError("EPIC solver must be auto, nm, nmf, or qp")

    def operation() -> pd.DataFrame:
        backend = _load_backend("epic")
        reference = backend.EpicReference(
            ref_profiles=sig,
            sig_genes=list(sig.index),
            ref_profiles_var=None,
        )
        result = backend.EPIC(
            bulk=mix,
            reference=reference,
            mRNA_cell=mrna_values,
            # The original DECEPTICON custom path uses EPIC's default
            # withOtherCells=TRUE. Keep that fit semantics, then expose only
            # the requested signature columns through the common contract.
            withOtherCells=True,
            solver=solver_name,
            backend=backend_name,
            n_threads=int(threads),
        )
        estimate = _result_frame(result, attribute="cellFractions", key="cellFractions")
        if isinstance(estimate, pd.DataFrame):
            estimate = estimate.reindex(columns=list(sig.columns))
        finished = _finish_output(
            estimate,
            samples=mix.columns,
            cell_types=sig.columns,
            method="epic",
        )
        native = _native_available("epic", backend)
        # epic_py only dispatches the vectorized Nelder-Mead paths (``auto``
        # and ``nmf``) to its optional Rust kernel.  The R-faithful ``nm``
        # and active-set ``qp`` solvers always execute in Python, even when
        # backend="rust" was requested.  Mirror epic_py's EPIC_BACKEND
        # override so provenance describes the execution class without
        # claiming that ``auto`` necessarily invoked the native subset.
        environment_backend = os.environ.get("EPIC_BACKEND")
        effective_backend = (
            environment_backend if environment_backend is not None else backend_name
        ).lower()
        native_eligible = solver_name in {"auto", "nmf"}
        if not native_eligible:
            resolved = "python"
            resolution_reason = "solver_python_only"
        elif effective_backend == "python":
            resolved = "python"
            resolution_reason = "python_backend_selected"
        elif effective_backend not in {"auto", "rust"}:
            # epic_py currently falls through to its Python implementation
            # for an unsupported EPIC_BACKEND value.
            resolved = "python"
            resolution_reason = "unsupported_backend_python_fallback"
        elif effective_backend in {"auto", "rust"} and native:
            if solver_name == "auto":
                # epic_py probes every sample with its Python QP solver and
                # dispatches only the unresolved subset to Rust.  That subset
                # may be empty, so calling this a pure Rust execution would be
                # false provenance even though the native kernel is available.
                resolved = "hybrid"
                resolution_reason = "python_qp_probe_rust_for_unresolved"
            else:
                resolved = "rust"
                resolution_reason = "native_kernel_selected"
        else:
            resolved = "python"
            resolution_reason = "native_kernel_unavailable"
        finished.attrs["decepticonx_engine"] = resolved
        details: dict[str, Any] = {
            "backend_requested": backend_name,
            "backend_effective": effective_backend,
            "backend_resolved": resolved,
            "solver_requested": solver_name,
            "native_available": native,
            "native_eligible": native_eligible,
            "backend_resolution_reason": resolution_reason,
        }
        fit_gof = getattr(result, "fit_gof", None)
        if isinstance(fit_gof, pd.DataFrame):
            details["fit_gof"] = _frame_payload(fit_gof)
        finished.attrs["decepticonx_diagnostics"] = details
        return finished

    return _call_backend("epic", operation)


def run_deconrnaseq(
    signature: pd.DataFrame,
    bulk: pd.DataFrame,
    *,
    backend_name: str = "auto",
) -> pd.DataFrame:
    """Run DeconRNASeq's scaled, non-plotting custom-signature path."""

    sig, mix = _prepare_inputs(signature, bulk)

    def operation() -> pd.DataFrame:
        backend = _load_backend("deconrnaseq")
        result = backend.deconrnaseq(
            datasets=mix,
            signatures=sig,
            use_scale=True,
            fig=False,
            backend=backend_name,
        )
        estimate = _result_frame(result, attribute="out_all", key="out.all")
        finished = _finish_output(
            estimate,
            samples=mix.columns,
            cell_types=sig.columns,
            method="deconrnaseq",
        )
        native = _native_available("deconrnaseq", backend)
        finished.attrs["decepticonx_engine"] = (
            backend_name
            if backend_name != "auto"
            else ("rust" if native else "numpy")
        )
        finished.attrs["decepticonx_diagnostics"] = {
            "backend_requested": backend_name,
            "native_available": native,
        }
        return finished

    return _call_backend("deconrnaseq", operation)


def _music_custom_reference(backend: Any, signature: pd.DataFrame) -> Any:
    """Mirror DECEPTICON's custom MuSiC pseudo-reference exactly.

    The R implementation cbinds the signature five times.  Each donor block
    contains one pseudo-cell per cell type, so the Python cells x genes layout
    is five vertically repeated copies of ``signature.T``.
    """

    cell_types = list(signature.columns)
    n_types = len(cell_types)
    counts = np.tile(signature.to_numpy(dtype=np.float64).T, (5, 1))
    obs = pd.DataFrame(
        {
            "sampleID": [donor for donor in range(1, 6) for _ in range(n_types)],
            "subjectname": [
                f"num{donor}" for donor in range(1, 6) for _ in range(n_types)
            ],
            "celltypeID": list(range(1, n_types + 1)) * 5,
            "celltype": cell_types * 5,
        },
        index=[f"cell{index}" for index in range(1, n_types * 5 + 1)],
    )
    return backend.SCReference(
        counts=counts,
        obs=obs,
        gene_names=list(signature.index),
    )


def run_music(
    signature: pd.DataFrame,
    bulk: pd.DataFrame,
    *,
    backend_name: str = "auto",
) -> pd.DataFrame:
    """Run MuSiC using DECEPTICON's five-donor custom-signature construction."""

    sig, mix = _prepare_inputs(signature, bulk)

    def operation() -> pd.DataFrame:
        backend = _load_backend("music")
        reference = _music_custom_reference(backend, sig)
        result = backend.music_prop(
            mix,
            reference,
            clusters="celltype",
            samples="sampleID",
            verbose=False,
            backend=backend_name,
        )
        if callable(getattr(result, "to_pandas", None)):
            estimate = result.to_pandas()["Est.prop.weighted"]
        else:
            estimate = _result_frame(
                result,
                attribute="est_prop_weighted",
                key="Est.prop.weighted",
            )
            if not isinstance(estimate, pd.DataFrame):
                sample_ids = getattr(result, "sample_ids", list(mix.columns))
                cell_types = getattr(result, "cell_types", list(sig.columns))
                estimate = pd.DataFrame(
                    estimate,
                    index=sample_ids,
                    columns=cell_types,
                )
        finished = _finish_output(
            estimate,
            samples=mix.columns,
            cell_types=sig.columns,
            method="music",
        )
        native = _native_available("music", backend)
        finished.attrs["decepticonx_engine"] = (
            backend_name
            if backend_name != "auto"
            else ("rust" if native else "numpy")
        )
        finished.attrs["decepticonx_diagnostics"] = {
            "backend_requested": backend_name,
            "native_available": native,
        }
        return finished

    return _call_backend("music", operation)


def _canonical_methods(methods: Sequence[str]) -> tuple[str, ...]:
    canonical: list[str] = []
    unknown: list[str] = []
    for value in methods:
        token = str(value).strip().lower().replace("-", "_").replace(" ", "")
        method = _METHOD_ALIASES.get(token)
        if method is None:
            unknown.append(str(value))
        elif method not in canonical:
            canonical.append(method)
    if unknown:
        raise InputValidationError(
            "Unknown deconvolution method(s): " + ", ".join(unknown)
        )
    return tuple(canonical)


def run_backends(
    signatures: Mapping[str, pd.DataFrame],
    bulk: pd.DataFrame,
    *,
    methods: Sequence[str],
    qn: bool = True,
    seed: int = 0,
    threads: int = 1,
    epic_mrna_cell: Mapping[str, float] | None = None,
    strict: bool = True,
    cibersort_engine: str = "rust",
    deconrnaseq_backend: str = "auto",
    epic_backend: str = "auto",
    epic_solver: str = "auto",
    music_backend: str = "auto",
) -> tuple[dict[BranchKey, pd.DataFrame], dict[str, Any]]:
    """Run requested method/reference branches behind one stable contract.

    When both CIBERSORT modes are requested, ``cibersort_all`` is invoked only
    once per reference.  With ``strict=False``, unavailable or failed branches
    are omitted and their clear error messages are retained in diagnostics.
    """

    if not isinstance(signatures, Mapping) or not signatures:
        raise InputValidationError("signatures must be a non-empty mapping")
    mix = _matrix(bulk, name="bulk")
    selected = _canonical_methods(methods)
    if threads < 1:
        raise InputValidationError("threads must be at least 1")

    estimates: dict[BranchKey, pd.DataFrame] = {}
    diagnostics: dict[str, Any] = {
        "requested_methods": list(selected),
        "completed": [],
        "skipped": {},
        "warnings": {},
        "engines": {},
        "branch_details": {},
        "backend_calls": {
            "cibersort_all": 0,
            "epic": 0,
            "deconrnaseq": 0,
            "music": 0,
        },
    }

    def record(
        reference_name: str,
        branch_methods: Sequence[str],
        operation: Callable[[], Mapping[str, pd.DataFrame]],
        call_name: str,
    ) -> None:
        keys = [BranchKey(method=method, reference=reference_name) for method in branch_methods]
        caught: list[warnings.WarningMessage] = []
        try:
            # Native/scientific backends often communicate convergence and
            # conditioning problems as warnings while still returning a
            # usable estimate.  Preserve those warnings in run.json instead
            # of letting them disappear in a batch log.
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                values = operation()
            diagnostics["backend_calls"][call_name] += 1
            for key in keys:
                estimates[key] = values[key.method]
                diagnostics["completed"].append(key.slug)
                engine = values[key.method].attrs.get("decepticonx_engine")
                if engine is not None:
                    diagnostics["engines"][key.slug] = engine
                details = values[key.method].attrs.get("decepticonx_diagnostics")
                if details:
                    diagnostics["branch_details"][key.slug] = details
        except (BackendUnavailableError, BackendExecutionError, InputValidationError) as exc:
            if strict:
                raise
            for key in keys:
                diagnostics["skipped"][key.slug] = str(exc)
        finally:
            messages: list[str] = []
            seen: set[str] = set()
            for item in caught:
                message = f"{item.category.__name__}: {item.message}"
                if message not in seen:
                    seen.add(message)
                    messages.append(message)
            if messages:
                for key in keys:
                    diagnostics["warnings"][key.slug] = messages

    for reference, signature in signatures.items():
        reference_name = str(reference)
        if not reference_name:
            raise InputValidationError("signature reference names must be non-empty")

        cib_methods = [
            method for method in ("cibersort", "cibersort_abs")
            if method in selected
        ]
        if cib_methods:
            def run_cib(sig: pd.DataFrame = signature) -> Mapping[str, pd.DataFrame]:
                relative, absolute = run_cibersort_pair(
                    sig,
                    mix,
                    qn=qn,
                    seed=seed,
                    threads=threads,
                    engine=cibersort_engine,
                )
                return {"cibersort": relative, "cibersort_abs": absolute}

            record(reference_name, cib_methods, run_cib, "cibersort_all")

        if "epic" in selected:
            record(
                reference_name,
                ("epic",),
                lambda sig=signature: {
                    "epic": run_epic(
                        sig,
                        mix,
                        mrna_cell=epic_mrna_cell,
                        threads=threads,
                        backend_name=epic_backend,
                        solver=epic_solver,
                    )
                },
                "epic",
            )

        if "deconrnaseq" in selected:
            record(
                reference_name,
                ("deconrnaseq",),
                lambda sig=signature: {
                    "deconrnaseq": run_deconrnaseq(
                        sig,
                        mix,
                        backend_name=deconrnaseq_backend,
                    )
                },
                "deconrnaseq",
            )

        if "music" in selected:
            record(
                reference_name,
                ("music",),
                lambda sig=signature: {
                    "music": run_music(sig, mix, backend_name=music_backend)
                },
                "music",
            )

    return estimates, diagnostics


__all__ = [
    "probe_backends",
    "probe_runtime",
    "run_backends",
    "run_cibersort_pair",
    "run_deconrnaseq",
    "run_epic",
    "run_music",
]
