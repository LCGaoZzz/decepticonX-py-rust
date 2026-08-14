from __future__ import annotations

from types import SimpleNamespace
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

from decepticonx.backends import (
    probe_backends,
    run_backends,
    run_deconrnaseq,
    run_epic,
    run_music,
)
from decepticonx.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    InputValidationError,
)
from decepticonx.models import BranchKey


@pytest.fixture
def matrices() -> tuple[pd.DataFrame, pd.DataFrame]:
    signature = pd.DataFrame(
        [[8.0, 1.0], [2.0, 7.0], [4.0, 3.0]],
        index=["g1", "g2", "g3"],
        columns=["B", "T"],
    )
    bulk = pd.DataFrame(
        [[5.0, 4.0], [6.0, 3.0], [2.0, 7.0]],
        index=["g3", "g1", "g2"],
        columns=["sample_a", "sample_b"],
    )
    return signature, bulk


def test_cibersort_modes_share_one_fit_and_strip_diagnostics(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    bulk_with_extra = pd.concat(
        [
            bulk,
            pd.DataFrame(
                [[1000.0, 1.0]], index=["bulk_only"], columns=bulk.columns
            ),
        ]
    )
    calls: list[tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]] = []

    def cibersort_all(sig, mix, **kwargs):
        calls.append((sig.copy(), mix.copy(), kwargs))
        # Deliberately reverse both axes; the adapter must restore the input
        # sample/cell-type order and remove all diagnostic columns.
        relative = pd.DataFrame(
            {
                "T": [0.8, 0.3],
                "B": [0.2, 0.7],
                "P-value": [9999.0, 9999.0],
                "Correlation": [0.9, 0.8],
                "RMSE": [0.1, 0.2],
            },
            index=["sample_b", "sample_a"],
        )
        absolute = pd.DataFrame(
            {
                "T": [1.6, 0.6],
                "B": [0.4, 1.4],
                "P-value": [9999.0, 9999.0],
                "Correlation": [0.9, 0.8],
                "RMSE": [0.1, 0.2],
                "Absolute score (sig.score)": [2.0, 2.0],
            },
            index=["sample_b", "sample_a"],
        )
        return {
            "relative": SimpleNamespace(table=relative),
            "sig.score": SimpleNamespace(table=absolute),
            "no.sumto1": SimpleNamespace(table=absolute),
        }

    monkeypatch.setitem(
        sys.modules,
        "python_cibersort",
        SimpleNamespace(cibersort_all=cibersort_all),
    )

    estimates, diagnostics = run_backends(
        {"bayesprism": signature},
        bulk_with_extra,
        methods=("CIBERSORT", "CIBERSORT-ABS"),
        qn=True,
        seed=17,
        threads=3,
    )

    assert len(calls) == 1
    passed_signature, passed_bulk, kwargs = calls[0]
    pd.testing.assert_frame_equal(passed_signature, signature)
    # CIBERSORT quantile-normalises the full mixture before intersecting
    # genes, so the adapter must preserve bulk-only rows at this boundary.
    pd.testing.assert_frame_equal(passed_bulk, bulk_with_extra)
    assert kwargs == {
        "perm": 0,
        "QN": True,
        "seed": 17,
        "threads": 3,
        "engine": "rust",
    }
    relative = estimates[BranchKey("cibersort", "bayesprism")]
    absolute = estimates[BranchKey("cibersort_abs", "bayesprism")]
    assert list(relative.index) == list(bulk_with_extra.columns)
    assert list(relative.columns) == list(signature.columns)
    assert list(absolute.index) == list(bulk_with_extra.columns)
    assert list(absolute.columns) == list(signature.columns)
    np.testing.assert_allclose(relative.loc["sample_a"], [0.7, 0.3])
    np.testing.assert_allclose(absolute.loc["sample_a"], [1.4, 0.6])
    assert diagnostics["backend_calls"]["cibersort_all"] == 1


def test_epic_requires_explicit_licensed_mrna_mapping(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    called = False

    def epic(**kwargs):
        nonlocal called
        called = True
        raise AssertionError("EPIC must not run without explicit mRNA_cell values")

    monkeypatch.setitem(
        sys.modules,
        "epic_py",
        SimpleNamespace(EpicReference=lambda **kwargs: kwargs, EPIC=epic),
    )

    with pytest.raises(BackendUnavailableError, match="explicit mRNA_cell mapping") as error:
        run_epic(signature, bulk, mrna_cell=None)
    assert "licensed" in str(error.value)
    assert called is False


def test_epic_uses_only_in_memory_custom_reference_and_explicit_mapping(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    captured: dict[str, object] = {}
    monkeypatch.delenv("EPIC_BACKEND", raising=False)

    class EpicReference:
        def __init__(self, **kwargs):
            captured["reference_kwargs"] = kwargs

    def forbidden_loader(*args, **kwargs):
        raise AssertionError("licensed EPIC data must never be loaded")

    def epic(**kwargs):
        captured["epic_kwargs"] = kwargs
        return SimpleNamespace(
            cellFractions=pd.DataFrame(
                [[0.7, 0.2, 0.1], [0.2, 0.7, 0.1]],
                index=["sample_b", "sample_a"],
                columns=["T", "B", "otherCells"],
            ),
            fit_gof=pd.DataFrame(
                {
                    "convergeCode": [0, 1],
                    "convergeMessage": ["OK", "iteration limit"],
                },
                index=["sample_b", "sample_a"],
            ),
        )

    monkeypatch.setitem(
        sys.modules,
        "epic_py",
        SimpleNamespace(
            EpicReference=EpicReference,
            EPIC=epic,
            load_reference=forbidden_loader,
            load_mrna_cell_default=forbidden_loader,
        ),
    )

    estimate = run_epic(
        signature,
        bulk,
        mrna_cell={"B": 1.5, "T": 2.5, "otherCells": 1.0},
        threads=4,
        backend_name="rust",
        solver="nm",
    )

    ref_kwargs = captured["reference_kwargs"]
    pd.testing.assert_frame_equal(ref_kwargs["ref_profiles"], signature)
    assert ref_kwargs["sig_genes"] == list(signature.index)
    assert ref_kwargs["ref_profiles_var"] is None
    epic_kwargs = captured["epic_kwargs"]
    pd.testing.assert_frame_equal(epic_kwargs["bulk"], bulk)
    assert isinstance(epic_kwargs["reference"], EpicReference)
    assert epic_kwargs["mRNA_cell"] == {"B": 1.5, "T": 2.5, "otherCells": 1.0}
    assert epic_kwargs["withOtherCells"] is True
    assert epic_kwargs["solver"] == "nm"
    assert epic_kwargs["backend"] == "rust"
    assert epic_kwargs["n_threads"] == 4
    assert list(estimate.index) == ["sample_a", "sample_b"]
    assert list(estimate.columns) == ["B", "T", "otherCells"]
    np.testing.assert_allclose(estimate.loc["sample_a"], [0.7, 0.2, 0.1])
    details = estimate.attrs["decepticonx_diagnostics"]
    fit_payload = details["fit_gof"]
    assert estimate.attrs["decepticonx_engine"] == "python"
    assert details["backend_requested"] == "rust"
    assert details["backend_effective"] == "rust"
    assert details["backend_resolved"] == "python"
    assert details["solver_requested"] == "nm"
    assert details["native_eligible"] is False
    assert details["backend_resolution_reason"] == "solver_python_only"
    assert fit_payload["columns"] == ["convergeCode", "convergeMessage"]
    assert fit_payload["records"][1] == {
        "sample": "sample_a",
        "convergeCode": 1,
        "convergeMessage": "iteration limit",
    }


@pytest.mark.parametrize(
    ("solver", "backend_name", "native", "expected", "reason"),
    [
        ("nm", "rust", True, "python", "solver_python_only"),
        ("qp", "rust", True, "python", "solver_python_only"),
        ("nmf", "rust", True, "rust", "native_kernel_selected"),
        (
            "auto",
            "auto",
            True,
            "hybrid",
            "python_qp_probe_rust_for_unresolved",
        ),
        ("nmf", "auto", False, "python", "native_kernel_unavailable"),
        ("auto", "python", True, "python", "python_backend_selected"),
    ],
)
def test_epic_diagnostics_describe_the_possible_execution_path(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
    solver: str,
    backend_name: str,
    native: bool,
    expected: str,
    reason: str,
) -> None:
    signature, bulk = matrices
    monkeypatch.delenv("EPIC_BACKEND", raising=False)
    monkeypatch.setattr(
        "decepticonx.backends._native_available",
        lambda method, module=None: native,
    )
    monkeypatch.setitem(
        sys.modules,
        "epic_py",
        SimpleNamespace(
            EpicReference=lambda **kwargs: kwargs,
            EPIC=lambda **kwargs: SimpleNamespace(
                cellFractions=pd.DataFrame(
                    [[0.7, 0.2, 0.1], [0.2, 0.7, 0.1]],
                    index=["sample_b", "sample_a"],
                    columns=["T", "B", "otherCells"],
                )
            ),
        ),
    )

    estimate = run_epic(
        signature,
        bulk,
        mrna_cell={"default": 1.0},
        backend_name=backend_name,
        solver=solver,
    )

    details = estimate.attrs["decepticonx_diagnostics"]
    assert estimate.attrs["decepticonx_engine"] == expected
    assert details["backend_requested"] == backend_name
    assert details["backend_effective"] == backend_name
    assert details["backend_resolved"] == expected
    assert details["solver_requested"] == solver
    assert details["native_available"] is native
    assert details["native_eligible"] is (solver in {"auto", "nmf"})
    assert details["backend_resolution_reason"] == reason


def test_epic_diagnostics_preserve_request_when_environment_overrides_backend(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    monkeypatch.setenv("EPIC_BACKEND", "python")
    monkeypatch.setattr(
        "decepticonx.backends._native_available",
        lambda method, module=None: True,
    )
    monkeypatch.setitem(
        sys.modules,
        "epic_py",
        SimpleNamespace(
            EpicReference=lambda **kwargs: kwargs,
            EPIC=lambda **kwargs: SimpleNamespace(
                cellFractions=pd.DataFrame(
                    [[0.7, 0.2, 0.1], [0.2, 0.7, 0.1]],
                    index=["sample_b", "sample_a"],
                    columns=["T", "B", "otherCells"],
                )
            ),
        ),
    )

    estimate = run_epic(
        signature,
        bulk,
        mrna_cell={"default": 1.0},
        backend_name="rust",
        solver="nmf",
    )

    details = estimate.attrs["decepticonx_diagnostics"]
    assert details["backend_requested"] == "rust"
    assert details["backend_effective"] == "python"
    assert details["backend_resolved"] == "python"
    assert details["backend_resolution_reason"] == "python_backend_selected"
    assert estimate.attrs["decepticonx_engine"] == "python"


def test_run_backends_forwards_epic_solver(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    captured: dict[str, object] = {}

    def epic_adapter(sig, mix, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame(
            [[0.6, 0.4], [0.3, 0.7]],
            index=mix.columns,
            columns=sig.columns,
        )

    monkeypatch.setattr("decepticonx.backends.run_epic", epic_adapter)
    estimates, _ = run_backends(
        {"music2": signature},
        bulk,
        methods=("epic",),
        epic_mrna_cell={"default": 1.0},
        epic_solver="qp",
    )

    assert BranchKey("epic", "music2") in estimates
    assert captured["solver"] == "qp"


def test_run_backends_records_epic_resolved_engine(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    monkeypatch.delenv("EPIC_BACKEND", raising=False)
    monkeypatch.setattr(
        "decepticonx.backends._native_available",
        lambda method, module=None: True,
    )
    monkeypatch.setitem(
        sys.modules,
        "epic_py",
        SimpleNamespace(
            EpicReference=lambda **kwargs: kwargs,
            EPIC=lambda **kwargs: SimpleNamespace(
                cellFractions=pd.DataFrame(
                    [[0.7, 0.2, 0.1], [0.2, 0.7, 0.1]],
                    index=["sample_b", "sample_a"],
                    columns=["T", "B", "otherCells"],
                )
            ),
        ),
    )

    _, diagnostics = run_backends(
        {"music2": signature},
        bulk,
        methods=("epic",),
        epic_mrna_cell={"default": 1.0},
        epic_backend="rust",
        epic_solver="qp",
    )

    slug = "epic__music2"
    assert diagnostics["engines"][slug] == "python"
    assert diagnostics["branch_details"][slug]["backend_requested"] == "rust"
    assert diagnostics["branch_details"][slug]["backend_resolved"] == "python"


def test_deconrnaseq_preserves_gene_by_sample_direction_and_parameters(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    captured: dict[str, object] = {}

    def deconrnaseq(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(
            out_all=pd.DataFrame(
                [[0.4, 0.6], [0.7, 0.3]],
                index=["sample_b", "sample_a"],
                columns=["T", "B"],
            )
        )

    monkeypatch.setitem(
        sys.modules,
        "deconrnaseq",
        SimpleNamespace(deconrnaseq=deconrnaseq),
    )

    estimate = run_deconrnaseq(signature, bulk, backend_name="numpy")

    pd.testing.assert_frame_equal(captured["datasets"], bulk)
    pd.testing.assert_frame_equal(captured["signatures"], signature)
    assert captured["use_scale"] is True
    assert captured["fig"] is False
    assert captured["backend"] == "numpy"
    assert list(estimate.index) == list(bulk.columns)
    assert list(estimate.columns) == list(signature.columns)
    np.testing.assert_allclose(estimate.loc["sample_a"], [0.3, 0.7])


def test_music_repeats_signature_as_five_donor_screference(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    captured: dict[str, object] = {}

    class SCReference:
        def __init__(self, *, counts, obs, gene_names):
            self.counts = counts
            self.obs = obs
            self.gene_names = gene_names
            captured["reference"] = self

    def music_prop(bulk_mtx, sc_sce, **kwargs):
        captured["bulk"] = bulk_mtx
        captured["music_reference"] = sc_sce
        captured["kwargs"] = kwargs
        return SimpleNamespace(
            est_prop_weighted=np.array([[0.8, 0.2], [0.3, 0.7]]),
            sample_ids=["sample_b", "sample_a"],
            cell_types=["T", "B"],
        )

    monkeypatch.setitem(
        sys.modules,
        "music_py",
        SimpleNamespace(SCReference=SCReference, music_prop=music_prop),
    )

    estimate = run_music(signature, bulk, backend_name="rust")

    reference = captured["reference"]
    expected_block = signature.to_numpy().T
    assert reference.counts.shape == (10, 3)
    for donor in range(5):
        np.testing.assert_allclose(
            reference.counts[donor * 2 : (donor + 1) * 2],
            expected_block,
        )
    assert reference.gene_names == list(signature.index)
    assert reference.obs["sampleID"].tolist() == [
        1, 1, 2, 2, 3, 3, 4, 4, 5, 5
    ]
    assert reference.obs["celltype"].tolist() == ["B", "T"] * 5
    pd.testing.assert_frame_equal(captured["bulk"], bulk)
    assert captured["kwargs"] == {
        "clusters": "celltype",
        "samples": "sampleID",
        "verbose": False,
        "backend": "rust",
    }
    assert list(estimate.index) == list(bulk.columns)
    assert list(estimate.columns) == list(signature.columns)
    np.testing.assert_allclose(estimate.loc["sample_a"], [0.7, 0.3])


@pytest.mark.parametrize(
    "bad_value, message",
    [(np.nan, "NaN or infinite"), (-0.01, "negative estimates")],
)
def test_outputs_must_be_finite_and_nonnegative(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
    bad_value: float,
    message: str,
) -> None:
    signature, bulk = matrices
    values = pd.DataFrame(
        [[bad_value, 0.5], [0.5, 0.5]],
        index=bulk.columns,
        columns=signature.columns,
    )
    monkeypatch.setitem(
        sys.modules,
        "deconrnaseq",
        SimpleNamespace(
            deconrnaseq=lambda **kwargs: SimpleNamespace(out_all=values)
        ),
    )
    with pytest.raises(BackendExecutionError, match=message):
        run_deconrnaseq(signature, bulk)


def test_backend_rejects_a_rank_deficient_common_signature(
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    signature["T"] = signature["B"] * 2.0
    with pytest.raises(InputValidationError, match="rank deficient"):
        run_deconrnaseq(signature, bulk)


def test_backend_rejects_a_target_with_no_mass_on_common_genes(
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    signature["T"] = 0.0
    with pytest.raises(InputValidationError, match="no positive expression"):
        run_deconrnaseq(signature, bulk)


def test_backend_rejects_too_few_common_genes_for_requested_types(
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    one_common_gene = bulk.rename(index={"g2": "x2", "g3": "x3"})
    with pytest.raises(InputValidationError, match="too few genes in common"):
        run_deconrnaseq(signature, one_common_gene)


def test_non_strict_mode_records_epic_authorization_gate(
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    estimates, diagnostics = run_backends(
        {"music2": signature},
        bulk,
        methods=("epic",),
        epic_mrna_cell=None,
        strict=False,
    )
    assert estimates == {}
    assert "explicit mRNA_cell mapping" in diagnostics["skipped"][
        "epic__music2"
    ]


def test_backend_warnings_are_retained_in_diagnostics(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices

    def deconrnaseq(**kwargs):
        warnings.warn("solver reached its iteration limit", RuntimeWarning)
        return SimpleNamespace(
            out_all=pd.DataFrame(
                [[0.4, 0.6], [0.7, 0.3]],
                index=["sample_b", "sample_a"],
                columns=["T", "B"],
            )
        )

    monkeypatch.setitem(
        sys.modules,
        "deconrnaseq",
        SimpleNamespace(deconrnaseq=deconrnaseq),
    )

    estimates, diagnostics = run_backends(
        {"music2": signature},
        bulk,
        methods=("deconrnaseq",),
    )

    assert BranchKey("deconrnaseq", "music2") in estimates
    assert diagnostics["warnings"]["deconrnaseq__music2"] == [
        "RuntimeWarning: solver reached its iteration limit"
    ]
    assert diagnostics["engines"]["deconrnaseq__music2"] == "numpy"
    details = diagnostics["branch_details"]["deconrnaseq__music2"]
    assert details["backend_requested"] == "auto"
    assert details["native_available"] is False
    assert details["input_overlap"] == {
        "signature_gene_count": 3,
        "bulk_gene_count": 3,
        "common_gene_count": 3,
        "signature_overlap_fraction": 1.0,
        "bulk_overlap_fraction": 1.0,
        "common_signature_rank": 2,
        "cell_type_count": 2,
    }


def test_probe_reports_shared_cibersort_package(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(
        sys.modules,
        "python_cibersort",
        SimpleNamespace(cibersort_all=lambda *args, **kwargs: None),
    )
    monkeypatch.setitem(
        sys.modules,
        "epic_py",
        SimpleNamespace(EpicReference=lambda **kwargs: None, EPIC=lambda **kwargs: None),
    )
    monkeypatch.setitem(
        sys.modules,
        "deconrnaseq",
        SimpleNamespace(deconrnaseq=lambda **kwargs: None),
    )
    monkeypatch.setitem(
        sys.modules,
        "music_py",
        SimpleNamespace(SCReference=lambda **kwargs: None, music_prop=lambda *args: None),
    )
    assert probe_backends() == {
        "cibersort": True,
        "cibersort_abs": True,
        "epic": True,
        "deconrnaseq": True,
        "music": True,
    }


def test_missing_package_has_direct_install_hint(
    monkeypatch: pytest.MonkeyPatch,
    matrices: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    signature, bulk = matrices
    real_import = __import__("decepticonx.backends", fromlist=["importlib"]).importlib.import_module

    def missing(name: str):
        if name == "deconrnaseq":
            raise ModuleNotFoundError("No module named 'deconrnaseq'", name=name)
        return real_import(name)

    monkeypatch.setattr("decepticonx.backends.importlib.import_module", missing)
    with pytest.raises(BackendUnavailableError, match="pip install") as error:
        run_deconrnaseq(signature, bulk)
    assert "deconrnaseq-py.git" in str(error.value)
