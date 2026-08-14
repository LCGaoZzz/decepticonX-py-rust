from __future__ import annotations

import json
from types import SimpleNamespace

import pandas as pd
import pytest

from decepticonx.exceptions import InputValidationError
from decepticonx.models import BranchKey, DecepticonXResult, PipelineConfig
from decepticonx import pipeline


def _bulk() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_1": [10.0, 4.0], "sample_2": [3.0, 9.0]},
        index=["G1", "G2"],
    )


def test_pipeline_connects_stages_records_metadata_and_writes(
    monkeypatch, tmp_path
) -> None:
    source_bulk = _bulk()
    validated_bulk = source_bulk.copy()
    prepared = SimpleNamespace(
        warnings=("compatibility warning",),
        provenance={"source_layer": "X", "n_cells_output": 4},
    )
    signature = pd.DataFrame(
        [[8.0, 1.0], [2.0, 7.0]],
        index=["G1", "G2"],
        columns=["B", "T"],
    )
    estimate = pd.DataFrame(
        [[0.7, 0.3], [0.2, 0.8]],
        index=source_bulk.columns,
        columns=signature.columns,
    )
    raw = estimate * 2.0
    normalized = estimate.copy()
    calls: dict[str, object] = {}

    def validate(frame):
        calls["bulk"] = frame
        return validated_bulk

    def prepare(single_cell, bulk, **kwargs):
        calls["prepare"] = (single_cell, bulk, kwargs)
        return prepared

    def build(prepared_value, bulk, **kwargs):
        calls["references"] = (prepared_value, bulk, kwargs)
        return {"music2": signature}, {"built": ["music2"]}

    def run(signatures, bulk, **kwargs):
        calls["backends"] = (signatures, bulk, kwargs)
        return (
            {BranchKey("music", "music2"): estimate},
            {"completed": ["music__music2"]},
        )

    def make_consensus(estimates, **kwargs):
        calls["consensus"] = (estimates, kwargs)
        return SimpleNamespace(
            raw=raw,
            normalized=normalized,
            diagnostics={"selected_pairs": 1},
        )

    monkeypatch.setattr(pipeline.io, "validate_bulk_expression", validate)
    monkeypatch.setattr(pipeline.io, "prepare_single_cell", prepare)
    monkeypatch.setattr(pipeline.references, "build_references", build)
    monkeypatch.setattr(pipeline.backends, "run_backends", run)
    monkeypatch.setattr(pipeline.consensus, "build_consensus", make_consensus)

    output = tmp_path / "result"
    config = PipelineConfig(
        methods=("music",),
        references=("music2",),
        threads=3,
        consensus_pairs=1,
        strict_backends=False,
    )
    result = pipeline.run_decepticonx(
        object(), source_bulk, config=config, output_dir=output
    )

    assert calls["bulk"] is source_bulk
    assert calls["prepare"][1] is validated_bulk
    assert calls["references"][2] == {
        "references": ("music2",),
        "species": "hs",
    }
    assert calls["backends"][2] == {
        "methods": ("music",),
        "qn": True,
        "seed": 0,
        "threads": 3,
        "epic_mrna_cell": None,
        "strict": False,
        "cibersort_engine": "rust",
        "deconrnaseq_backend": "auto",
        "epic_backend": "auto",
        "epic_solver": "auto",
        "music_backend": "auto",
    }
    assert calls["consensus"][1] == {"cell_types": ("B", "T"), "n_pairs": 1}
    pd.testing.assert_frame_equal(result.consensus, normalized)
    pd.testing.assert_frame_equal(result.consensus_raw, raw)
    assert result.diagnostics["references"] == {"built": ["music2"]}
    assert result.diagnostics["backends"]["completed"] == ["music__music2"]
    assert set(result.timings) == {
        "bulk_input",
        "single_cell_input",
        "references",
        "backends",
        "consensus",
        "total",
    }
    assert all(value >= 0 for value in result.timings.values())
    assert result.provenance["single_cell"]["source_layer"] == "X"
    assert result.provenance["bulk"]["n_samples"] == 2
    assert result.provenance["completed_branches"] == ["music__music2"]
    assert (output / "consensus.tsv").is_file()
    assert (output / "signatures" / "music2.tsv").is_file()
    metadata = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert metadata["provenance"]["config"]["threads"] == 3


def test_pipeline_loads_a_bulk_path_and_accepts_mapping_config(
    monkeypatch, tmp_path
) -> None:
    bulk_path = tmp_path / "bulk.tsv"
    frame = _bulk()
    prepared = SimpleNamespace(warnings=(), provenance={})
    signature = pd.DataFrame(
        [[1.0, 2.0], [2.0, 1.0]],
        index=frame.index,
        columns=["B", "T"],
    )
    estimate = pd.DataFrame(
        [[0.5, 0.5], [0.5, 0.5]], index=frame.columns, columns=signature.columns
    )
    loaded: list[object] = []

    monkeypatch.setattr(
        pipeline.io,
        "load_bulk_expression",
        lambda path: loaded.append(path) or frame,
    )
    monkeypatch.setattr(pipeline.io, "prepare_single_cell", lambda *a, **k: prepared)
    monkeypatch.setattr(
        pipeline.references,
        "build_references",
        lambda *a, **k: ({"bayesprism": signature}, {}),
    )
    monkeypatch.setattr(
        pipeline.backends,
        "run_backends",
        lambda *a, **k: ({BranchKey("cibersort", "bayesprism"): estimate}, {}),
    )
    monkeypatch.setattr(
        pipeline.consensus,
        "build_consensus",
        lambda *a, **k: SimpleNamespace(raw=estimate, normalized=estimate, diagnostics={}),
    )

    pipeline.run_decepticonx(
        "cells.h5ad",
        bulk_path,
        config={
            "methods": ["cibersort"],
            "references": ["bayesprism"],
            "exclude_cell_types": ["unknown"],
        },
    )

    assert loaded == [bulk_path]


def test_cell_types_use_first_reference_ordered_intersection() -> None:
    first = pd.DataFrame([[1, 2, 3]], index=["G"], columns=["B", "T", "NK"])
    second = pd.DataFrame([[1, 2, 3]], index=["G"], columns=["T", "B", "M"])
    assert pipeline._cell_types({"first": first, "second": second}) == ("B", "T")


@pytest.mark.parametrize(
    "field,value,match",
    [
        ("allow_normalized_x", "false", "must be a JSON/Python boolean"),
        ("cibersort_qn", 1, "must be a JSON/Python boolean"),
        ("threads", 1.5, "must be an integer"),
    ],
)
def test_configuration_rejects_ambiguous_json_types(field, value, match) -> None:
    with pytest.raises(InputValidationError, match=match):
        pipeline._configuration({field: value})


def test_epic_solver_is_validated_and_canonicalized() -> None:
    assert pipeline._configuration({"epic_solver": " NM "}).epic_solver == "nm"
    with pytest.raises(InputValidationError, match="epic_solver must be one of"):
        pipeline._configuration({"epic_solver": "not-a-solver"})


@pytest.mark.parametrize(
    "mapping",
    [
        {"T": True},
        {"T": "not-a-number"},
        {"T": 0},
        {"T": float("inf")},
        {"": 1.0},
    ],
)
def test_epic_mapping_is_validated_before_provenance_hashing(mapping) -> None:
    with pytest.raises(InputValidationError, match="epic_mrna_cell"):
        pipeline._configuration({"epic_mrna_cell": mapping})


def test_epic_configuration_fails_before_any_input_work(monkeypatch) -> None:
    monkeypatch.setattr(
        pipeline,
        "_bulk_input",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("input work must not start")
        ),
    )
    with pytest.raises(InputValidationError, match="epic_mrna_cell"):
        pipeline.run_decepticonx(
            "cells.h5ad",
            "bulk.tsv",
            config={"methods": ["epic"], "epic_mrna_cell": None},
        )


def test_epic_mapping_and_inline_command_are_redacted() -> None:
    config = PipelineConfig(
        epic_mrna_cell={"T": 2.5, "default": 1.0}, epic_solver="nm"
    )
    recorded = pipeline._provenance_config(config)
    assert recorded["epic_mrna_cell"]["entry_count"] == 2
    assert recorded["epic_mrna_cell"]["values_redacted"] is True
    assert recorded["epic_mrna_cell"]["hash_format"].startswith("sorted key")
    assert recorded["epic_mrna_cell"]["sha256"] == (
        "1360b12f568ec6a9777f7a4252ba9b92f35c589df81a919bf62ab334f91aae71"
    )
    assert recorded["epic_solver"] == "nm"
    assert "2.5" not in json.dumps(recorded)
    command = pipeline._redacted_command(
        ["decepticonx", "run", "--epic-mrna-cell", '{"T":2.5}']
    )
    assert command[-1] == "<redacted>"


def test_result_write_refuses_a_nonempty_output_directory(tmp_path) -> None:
    frame = pd.DataFrame([[1.0]], index=["s1"], columns=["T"])
    result = DecepticonXResult(
        signatures={"music2": frame.T},
        estimates={BranchKey("music", "music2"): frame},
        consensus=frame,
        consensus_raw=frame,
    )
    output = tmp_path / "result"
    result.write(output)
    with pytest.raises(FileExistsError, match="must not already contain"):
        result.write(output)


def test_default_methods_are_runnable_without_licensed_epic_data() -> None:
    assert "epic" not in PipelineConfig().methods
