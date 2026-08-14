"""Opt-in smoke tests against the real optional backend packages.

The regular unit suite uses small contract doubles so contributors can test
the orchestrator without cloning four separately released projects.  Set
``DECEPTICONX_RUN_BACKEND_INTEGRATION=1`` and make the optional packages
importable to exercise this module.  Each unavailable backend is skipped
independently; the public-backend CI job installs/checks out the commits in
``backend-lock.json`` and therefore runs the first three tests.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from decepticonx.backends import (
    run_backends,
    run_deconrnaseq,
    run_epic,
    run_music,
)
from decepticonx.models import BranchKey


pytestmark = pytest.mark.skipif(
    os.environ.get("DECEPTICONX_RUN_BACKEND_INTEGRATION") != "1",
    reason=(
        "real backend smoke tests are opt-in; set "
        "DECEPTICONX_RUN_BACKEND_INTEGRATION=1"
    ),
)


@pytest.fixture(scope="module")
def synthetic_mixture() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Positive, full-rank input with an exactly known three-type mixture."""

    rng = np.random.default_rng(20260814)
    genes = [f"G{index:03d}" for index in range(48)]
    cell_types = ["B", "T", "M"]
    values = rng.uniform(20.0, 80.0, size=(len(genes), len(cell_types)))
    for cell_type in range(len(cell_types)):
        start = cell_type * 12
        values[start : start + 12, cell_type] += 250.0

    signature = pd.DataFrame(values, index=genes, columns=cell_types)
    truth = pd.DataFrame(
        [[0.65, 0.25, 0.10], [0.15, 0.30, 0.55], [0.25, 0.50, 0.25]],
        index=["S1", "S2", "S3"],
        columns=cell_types,
    )
    bulk_values = values @ truth.to_numpy().T
    # Unequal library scales make the CIBERSORT sig.score result observably
    # distinct from its relative result while preserving the true fractions.
    bulk_values *= np.array([1.0, 1.8, 0.7])
    bulk = pd.DataFrame(bulk_values, index=genes, columns=truth.index)
    return signature, bulk, truth


def test_real_cibersort_relative_and_absolute_share_one_fit(
    synthetic_mixture: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend = pytest.importorskip("python_cibersort")
    signature, bulk, truth = synthetic_mixture
    actual = backend.cibersort_all
    calls: list[dict[str, object]] = []

    def observed_call(*args: object, **kwargs: object) -> object:
        calls.append(dict(kwargs))
        return actual(*args, **kwargs)

    monkeypatch.setattr(backend, "cibersort_all", observed_call)
    estimates, diagnostics = run_backends(
        {"synthetic": signature},
        bulk,
        methods=("cibersort", "cibersort_abs"),
        qn=False,
        seed=17,
        threads=1,
        cibersort_engine="numpy",
    )

    relative = estimates[BranchKey("cibersort", "synthetic")]
    absolute = estimates[BranchKey("cibersort_abs", "synthetic")]
    assert len(calls) == 1
    assert calls[0] == {
        "perm": 0,
        "QN": False,
        "seed": 17,
        "threads": 1,
        "engine": "numpy",
    }
    assert diagnostics["backend_calls"]["cibersort_all"] == 1
    assert list(relative.index) == list(truth.index)
    assert list(relative.columns) == list(truth.columns)
    assert np.isfinite(relative.to_numpy()).all()
    assert np.isfinite(absolute.to_numpy()).all()
    assert (relative.to_numpy() >= 0).all()
    assert (absolute.to_numpy() >= 0).all()
    np.testing.assert_allclose(relative.sum(axis=1), 1.0, atol=1e-12)
    np.testing.assert_allclose(relative, truth, atol=2e-3)
    # sig.score changes only the sample scale; normalising it recovers the
    # relative result generated from the same three nu-SVR fits.
    np.testing.assert_allclose(
        absolute.div(absolute.sum(axis=1), axis=0),
        relative,
        rtol=1e-12,
        atol=1e-12,
    )
    assert not np.allclose(absolute.sum(axis=1), 1.0)


def test_real_deconrnaseq_uses_scaled_estimator(
    synthetic_mixture: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend = pytest.importorskip("deconrnaseq")
    signature, bulk, _ = synthetic_mixture
    actual = backend.deconrnaseq
    calls: list[dict[str, object]] = []

    def observed_call(*args: object, **kwargs: object) -> object:
        calls.append(dict(kwargs))
        return actual(*args, **kwargs)

    monkeypatch.setattr(backend, "deconrnaseq", observed_call)
    estimate = run_deconrnaseq(signature, bulk, backend_name="numpy")
    direct = actual(
        datasets=bulk,
        signatures=signature,
        use_scale=True,
        fig=False,
        backend="numpy",
    ).out_all

    assert len(calls) == 1
    assert calls[0]["use_scale"] is True
    assert calls[0]["fig"] is False
    assert calls[0]["backend"] == "numpy"
    np.testing.assert_allclose(estimate, direct, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(estimate.sum(axis=1), 1.0, atol=1e-12)
    assert np.isfinite(estimate.to_numpy()).all()
    assert (estimate.to_numpy() >= 0).all()


def test_real_music_five_donor_path_recovers_exact_mixture(
    synthetic_mixture: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend = pytest.importorskip("music_py")
    signature, bulk, truth = synthetic_mixture
    actual = backend.music_prop
    calls: list[dict[str, object]] = []

    def observed_call(bulk_arg: object, reference: object, **kwargs: object) -> object:
        calls.append(
            {
                **kwargs,
                "reference_shape": tuple(reference.counts.shape),
                "donor_count": int(reference.obs["sampleID"].nunique()),
            }
        )
        return actual(bulk_arg, reference, **kwargs)

    monkeypatch.setattr(backend, "music_prop", observed_call)
    estimate = run_music(signature, bulk, backend_name="numpy")

    assert len(calls) == 1
    assert calls[0]["clusters"] == "celltype"
    assert calls[0]["samples"] == "sampleID"
    assert calls[0]["verbose"] is False
    assert calls[0]["backend"] == "numpy"
    assert calls[0]["reference_shape"] == (15, 48)
    assert calls[0]["donor_count"] == 5
    assert np.isfinite(estimate.to_numpy()).all()
    assert (estimate.to_numpy() >= 0).all()
    np.testing.assert_allclose(estimate.sum(axis=1), 1.0, atol=1e-12)
    np.testing.assert_allclose(estimate, truth, rtol=1e-9, atol=1e-10)


def test_real_epic_when_authorized_package_is_already_installed(
    synthetic_mixture: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    # EPIC is intentionally not fetched by public CI.  This test becomes an
    # executable contract only in an environment with the authorized package.
    pytest.importorskip("epic_py")
    signature, bulk, _ = synthetic_mixture
    estimate = run_epic(
        signature,
        bulk,
        mrna_cell={"default": 1.0},
        backend_name="python",
        solver="qp",
        threads=1,
    )

    assert list(estimate.columns) == [*signature.columns, "otherCells"]
    assert np.isfinite(estimate.to_numpy()).all()
    assert (estimate.to_numpy() >= 0).all()
    np.testing.assert_allclose(estimate.sum(axis=1), 1.0, atol=1e-8)
