from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from decepticonx.exceptions import InputValidationError
from decepticonx.io import (
    REMAINDER_GENE,
    load_bulk_expression,
    prepare_single_cell,
    validate_bulk_expression,
)


def _bulk(genes: tuple[str, ...] = ("A", "B")) -> pd.DataFrame:
    return pd.DataFrame({"bulk_1": np.arange(1, len(genes) + 1)}, index=genes)


def _adata(x: object, genes: tuple[str, ...] = ("A", "B", "C")) -> ad.AnnData:
    obs = pd.DataFrame(
        {"cell_type": ["T", "B"], "donor": ["d1", "d2"]},
        index=["cell_1", "cell_2"],
    )
    return ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=genes))


def _dense(x: object) -> np.ndarray:
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def test_load_bulk_expression_reads_genes_by_samples(tmp_path) -> None:
    path = tmp_path / "bulk.tsv"
    path.write_text("gene\ts1\ts2\nA\t1\t2.5\nB\t0\t3\n", encoding="utf-8")

    result = load_bulk_expression(path)

    assert result.index.tolist() == ["A", "B"]
    assert result.columns.tolist() == ["s1", "s2"]
    np.testing.assert_allclose(result.to_numpy(), [[1, 2.5], [0, 3]])


def test_load_bulk_expression_accepts_r_table_without_gene_header(tmp_path) -> None:
    path = tmp_path / "bulk.tsv"
    path.write_text("s1\ts2\nA\t1\t2.5\nB\t0\t3\n", encoding="utf-8")

    result = load_bulk_expression(path)

    assert result.index.tolist() == ["A", "B"]
    assert result.columns.tolist() == ["s1", "s2"]
    np.testing.assert_allclose(result.to_numpy(), [[1, 2.5], [0, 3]])


@pytest.mark.parametrize(
    "contents, match",
    [
        ("gene\ts1\ts2\nA\t1\n", "ragged"),
        ("gene\ts1\nA\tNA\n", "non-numeric"),
        ("gene\ts1\nA\tinf\n", "non-finite"),
        ("gene\ts1\nA\t1\nA\t2\n", "duplicate gene"),
        ("gene\ts1\ts1\nA\t1\t2\n", "duplicate sample"),
    ],
)
def test_load_bulk_expression_rejects_malformed_tables(
    tmp_path, contents: str, match: str
) -> None:
    path = tmp_path / "bulk.tsv"
    path.write_text(contents, encoding="utf-8")

    with pytest.raises(InputValidationError, match=match):
        load_bulk_expression(path)


def test_validate_bulk_rejects_dataframe_nan_and_duplicate_names() -> None:
    with pytest.raises(InputValidationError, match="missing, NaN, or infinite"):
        validate_bulk_expression(pd.DataFrame({"s": [1.0, np.nan]}, index=["A", "B"]))
    with pytest.raises(InputValidationError, match="duplicate sample"):
        validate_bulk_expression(pd.DataFrame([[1, 2]], index=["A"], columns=["s", "s"]))
    with pytest.raises(InputValidationError, match="negative"):
        validate_bulk_expression(pd.DataFrame({"s": [-1.0]}, index=["A"]))


def test_auto_prefers_counts_layer_over_normalized_x() -> None:
    adata = _adata(np.array([[0.1, 0.2, 0.3], [1.1, 0.0, 0.4]]))
    counts = sparse.csr_matrix([[1, 2, 7], [4, 0, 6]])
    adata.layers["counts"] = counts

    prepared = prepare_single_cell(adata, _bulk())

    assert prepared.source_layer == "layers['counts']"
    assert prepared.normalized_input is False
    assert prepared.cell_type_key == "cell_type"
    assert prepared.sample_key is None
    assert prepared.common_genes == ("A", "B")
    np.testing.assert_allclose(_dense(prepared.adata.X), [[1, 2, 7], [4, 0, 6]])


def test_noninteger_x_rejected_by_default_and_explicitly_allowed() -> None:
    adata = _adata(sparse.csr_matrix([[0.25, 1.5, 0], [2.1, 0.2, 4.3]]))

    with pytest.raises(InputValidationError, match="not integer-like"):
        prepare_single_cell(adata, _bulk())

    with pytest.warns(RuntimeWarning, match="allow_normalized_x=True"):
        prepared = prepare_single_cell(
            adata, _bulk(), allow_normalized_x=True, sample_key="donor"
        )
    assert prepared.normalized_input is True
    assert prepared.source_layer == "X"
    assert prepared.sample_key == "donor"
    assert prepared.warnings and "normalized space" in prepared.warnings[0]


def test_auto_uses_count_like_raw_before_allowing_normalized_x() -> None:
    adata = _adata(np.array([[0.2, 0.3, 0.5], [1.2, 0.1, 0.0]]))
    raw = _adata(sparse.csr_matrix([[1, 2, 3], [4, 5, 6]]))
    adata.raw = raw

    prepared = prepare_single_cell(adata, _bulk(), allow_normalized_x=True)

    assert prepared.source_layer == "raw.X"
    assert prepared.normalized_input is False
    assert not prepared.warnings


def test_remainder_preserves_full_library_and_duplicate_symbols_are_summed() -> None:
    # Two Ensembl columns map to A. C is not shared with bulk and must land in
    # the remainder. Excluding B cells exercises non-contiguous row selection.
    obs = pd.DataFrame(
        {"cell_type": ["T", "B", "T"]}, index=["c1", "c2", "c3"]
    )
    var = pd.DataFrame(
        {"symbol": ["A", "A", "B", "C"]},
        index=["ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004"],
    )
    source = sparse.csr_matrix([[1, 2, 3, 4], [5, 6, 7, 8], [0, 9, 1, 2]])
    adata = ad.AnnData(X=source, obs=obs, var=var)

    prepared = prepare_single_cell(
        adata,
        _bulk(("A", "B")),
        gene_symbol_key="symbol",
        exclude_cell_types=("B",),
    )

    assert prepared.adata.obs_names.tolist() == ["c1", "c3"]
    assert prepared.adata.var_names.tolist() == ["A", "B", REMAINDER_GENE]
    actual = _dense(prepared.adata.X)
    np.testing.assert_allclose(actual, [[3, 3, 4], [9, 1, 2]])
    np.testing.assert_allclose(actual.sum(axis=1), [10, 12])
    assert prepared.provenance["duplicate_symbol_source_columns"] == 2


def test_remainder_name_is_collision_safe() -> None:
    adata = _adata(np.array([[1, 2, 3], [4, 5, 6]]), ("A", "B", REMAINDER_GENE))

    prepared = prepare_single_cell(adata, _bulk())

    assert prepared.remainder_gene == REMAINDER_GENE + "_1"
    assert prepared.adata.var_names[-1] == REMAINDER_GENE + "_1"


@pytest.mark.parametrize(
    "obs, match",
    [
        (pd.DataFrame({"wrong": ["T", "B"]}, index=["c1", "c2"]), "missing"),
        (pd.DataFrame({"cell_type": ["T", None]}, index=["c1", "c2"]), "missing/empty"),
        (pd.DataFrame({"cell_type": ["T", "  "]}, index=["c1", "c2"]), "missing/empty"),
    ],
)
def test_cell_type_key_must_exist_and_be_complete(obs: pd.DataFrame, match: str) -> None:
    adata = ad.AnnData(
        X=np.array([[1, 2], [3, 4]]), obs=obs, var=pd.DataFrame(index=["A", "B"])
    )

    with pytest.raises(InputValidationError, match=match):
        prepare_single_cell(adata, _bulk())


def test_obs_labels_are_stripped_in_compact_anndata() -> None:
    adata = _adata(np.array([[1, 2, 3], [4, 5, 6]]))
    adata.obs["cell_type"] = [" T ", "B"]
    adata.obs["donor"] = [" d1", "d2 "]

    prepared = prepare_single_cell(adata, _bulk(), sample_key="donor")

    assert prepared.adata.obs["cell_type"].tolist() == ["T", "B"]
    assert prepared.adata.obs["donor"].tolist() == ["d1", "d2"]


def test_duplicate_bulk_gene_names_are_rejected() -> None:
    adata = _adata(np.array([[1, 2, 3], [4, 5, 6]]))
    bulk = pd.DataFrame({"s": [1, 2]}, index=["A", "A"])

    with pytest.raises(InputValidationError, match="duplicate gene"):
        prepare_single_cell(adata, bulk)


def test_ensembl_ids_require_symbol_mapping() -> None:
    adata = _adata(np.array([[1, 2, 3], [4, 5, 6]]), ("ENSG000001", "B", "C"))

    with pytest.raises(InputValidationError, match="gene_symbol_key"):
        prepare_single_cell(adata, _bulk())


def test_sparse_backed_h5ad_is_compacted_without_full_dense_result(tmp_path) -> None:
    adata = _adata(sparse.csr_matrix([[1, 0, 9], [0, 2, 8]]))
    path = tmp_path / "cells.h5ad"
    adata.write_h5ad(path)

    prepared = prepare_single_cell(path, _bulk())

    assert not prepared.adata.isbacked
    assert sparse.issparse(prepared.adata.X)
    np.testing.assert_allclose(_dense(prepared.adata.X).sum(axis=1), [10, 10])
    assert prepared.provenance["source_path"] == str(path.resolve())


def test_nonhuman_species_is_rejected() -> None:
    adata = _adata(np.array([[1, 2, 3], [4, 5, 6]]))

    with pytest.raises(InputValidationError, match="only human"):
        prepare_single_cell(adata, _bulk(), species="mouse")


def test_duplicate_cell_ids_and_zero_library_cells_are_rejected() -> None:
    duplicate = _adata(np.array([[1, 2, 3], [4, 5, 6]]))
    duplicate.obs_names = ["same", "same"]
    with pytest.raises(InputValidationError, match="duplicate cell IDs"):
        prepare_single_cell(duplicate, _bulk())

    zero = _adata(np.array([[1, 2, 3], [0, 0, 0]]))
    with pytest.raises(InputValidationError, match="zero-library"):
        prepare_single_cell(zero, _bulk())
