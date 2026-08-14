from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from decepticonx.exceptions import InputValidationError
from decepticonx.io import PreparedSingleCell
from decepticonx.references import build_references


@pytest.fixture
def prepared_fixture() -> PreparedSingleCell:
    # The first two cells are A and the final two are B.  Correct raw means:
    #                A     B
    #   g1           2    15
    #   g2           1     0
    # The synthetic remainder makes the library sizes deliberately very
    # different; it must be visible to Monocle3/legacy MuSiC2 but may not alter
    # the corrected MuSiC2 means.
    matrix = np.asarray(
        [
            [1.0, 1.0, 100.0],
            [3.0, 1.0, 200.0],
            [10.0, 0.0, 300.0],
            [20.0, 0.0, 400.0],
        ]
    )
    adata = ad.AnnData(
        matrix,
        obs=pd.DataFrame(
            {"kind": ["A", "A", "B", "B"]},
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["g1", "g2", "__remainder__"]),
    )
    return PreparedSingleCell(
        adata=adata,
        common_genes=("g1", "g2"),
        source_layer="X",
        normalized_input=False,
        cell_type_key="kind",
        sample_key=None,
        warnings=(),
        provenance={"remainder_gene": "__remainder__"},
    )


def _install_fake_fast(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    def bayesprism_base(bulk, sc, labels, species):
        calls["bayesprism"] = (bulk.copy(), sc.copy(), labels.copy(), species)
        means = {
            cell_type: sc.loc[:, labels == cell_type].mean(axis=1)
            for cell_type in pd.unique(labels)
        }
        return pd.DataFrame(means)

    def monocle3_base(sc, metadata):
        calls["monocle3"] = (sc.copy(), metadata.copy())
        labels = metadata.iloc[:, 0]
        means = {
            cell_type: sc.loc[:, labels == cell_type].mean(axis=1)
            for cell_type in sorted(pd.unique(labels))
        }
        return pd.DataFrame(means)

    def music2_base(sc, bulk, labels, expr_low=20, markers=None):
        calls["music2"] = (sc.copy(), bulk.copy(), labels.copy(), markers)
        # Intentionally reproduce a visibly wrong scalar-style result.  The
        # orchestrator must call this only as a compatibility diagnostic and
        # must not expose it as the public reference.
        return pd.DataFrame(
            999.0,
            index=list(reversed(bulk.index)),
            columns=pd.unique(labels),
        )

    fake = SimpleNamespace(
        bayesprism_base=bayesprism_base,
        monocle3_base=monocle3_base,
        music2_base=music2_base,
        get_backend=lambda: "test-double",
    )
    monkeypatch.setitem(__import__("sys").modules, "decepticon_fast", fake)


def test_music2_is_raw_type_mean_in_bulk_nonzero_order(
    monkeypatch: pytest.MonkeyPatch, prepared_fixture: PreparedSingleCell
) -> None:
    calls: dict = {}
    _install_fake_fast(monkeypatch, calls)
    # g3 is not common; g0 is common in neither input; zero-valued common genes
    # would be removed by the same non-zero rule (g2 remains because sample 1
    # is non-zero).  Most importantly, bulk order is g2 then g1.
    bulk = pd.DataFrame(
        [[4.0, 0.0], [8.0, 7.0], [3.0, 2.0]],
        index=["g2", "g1", "g3"],
        columns=["s1", "s2"],
    )

    references, diagnostics = build_references(
        prepared_fixture, bulk, references=("music2",)
    )

    result = references["music2"]
    assert result.index.tolist() == ["g2", "g1"]
    assert result.columns.tolist() == ["A", "B"]
    np.testing.assert_allclose(result.to_numpy(), [[1.0, 0.0], [2.0, 15.0]])

    assert "music2" not in calls
    assert diagnostics["references"]["music2"]["upstream_legacy_call_skipped"]
    assert diagnostics["references"]["music2"]["expression_scale"] == "raw_counts"


def test_remainder_visibility_is_builder_specific(
    monkeypatch: pytest.MonkeyPatch, prepared_fixture: PreparedSingleCell
) -> None:
    calls: dict = {}
    _install_fake_fast(monkeypatch, calls)
    bulk = pd.DataFrame(
        [[1.0], [1.0]], index=["g2", "g1"], columns=["s1"]
    )

    references, _ = build_references(
        prepared_fixture,
        bulk,
        references=("bayesprism", "monocle3", "music2"),
    )

    assert calls["bayesprism"][1].index.tolist() == ["g2", "g1"]
    assert calls["monocle3"][0].index.tolist() == [
        "g2",
        "g1",
        "__remainder__",
    ]
    assert "__remainder__" not in references["monocle3"].index
    assert references["monocle3"].index.tolist() == ["g2", "g1"]


def test_only_hs_species_is_accepted(
    monkeypatch: pytest.MonkeyPatch, prepared_fixture: PreparedSingleCell
) -> None:
    calls: dict = {}
    _install_fake_fast(monkeypatch, calls)
    bulk = pd.DataFrame([[1.0]], index=["g1"], columns=["s1"])
    with pytest.raises(InputValidationError, match="species='hs'"):
        build_references(
            prepared_fixture, bulk, references=("music2",), species="mm"
        )


def test_music2_enforces_the_original_twenty_percent_overlap_gate(
    monkeypatch: pytest.MonkeyPatch, prepared_fixture: PreparedSingleCell
) -> None:
    calls: dict = {}
    _install_fake_fast(monkeypatch, calls)
    prepared = replace(
        prepared_fixture,
        provenance={
            **prepared_fixture.provenance,
            "n_genes_source": 20,
        },
    )
    bulk = pd.DataFrame(
        np.ones((11, 1)),
        index=["g1", "g2", *[f"bulk_only_{index}" for index in range(9)]],
        columns=["s1"],
    )

    with pytest.raises(InputValidationError, match="requires at least 20%"):
        build_references(prepared, bulk, references=("music2",))
