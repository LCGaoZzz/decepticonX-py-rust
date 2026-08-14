from __future__ import annotations

import json
from types import SimpleNamespace

import pandas as pd

from decepticonx import cli


def test_run_command_builds_config_and_forwards_output(monkeypatch, tmp_path, capsys) -> None:
    captured: dict[str, object] = {}
    consensus = pd.DataFrame([[0.4, 0.6]], index=["s1"], columns=["B", "T"])

    def fake_run(single_cell, bulk, *, config, output_dir):
        captured.update(
            single_cell=single_cell,
            bulk=bulk,
            config=config,
            output_dir=output_dir,
        )
        return SimpleNamespace(
            signatures={"bayesprism": object()},
            estimates={"one": object()},
            consensus=consensus,
            timings={"total": 1.25},
        )

    monkeypatch.setattr(cli, "run_decepticonx", fake_run)
    output = tmp_path / "out"
    status = cli.main(
        [
            "run",
            "cells.h5ad",
            "bulk.tsv",
            "--output-dir",
            str(output),
            "--methods",
            "cibersort,music",
            "--reference",
            "music2",
            "--threads",
            "4",
            "--consensus-pairs",
            "1",
            "--consensus-mode",
            "corrected",
            "--no-cibersort-qn",
            "--permissive",
            "--allow-partial-consensus",
            "--epic-mrna-json",
            '{"B": 1.2, "T": 2.3}',
            "--epic-solver",
            "nm",
        ]
    )

    assert status == 0
    config = captured["config"]
    assert config.methods == ("cibersort", "music")
    assert config.references == ("music2",)
    assert config.threads == 4
    assert config.consensus_pairs == 1
    assert config.consensus_mode == "corrected"
    assert config.cibersort_qn is False
    assert config.strict_backends is False
    assert config.allow_partial_consensus is True
    assert config.epic_mrna_cell == {"B": 1.2, "T": 2.3}
    assert config.epic_solver == "nm"
    assert captured["output_dir"] == str(output)
    report = json.loads(capsys.readouterr().out)
    assert report["ok"] is True
    assert report["branches"] == 1


def test_validate_only_calls_input_stages(monkeypatch, capsys) -> None:
    bulk = pd.DataFrame({"s1": [1.0]}, index=["G1"])
    calls: dict[str, object] = {}

    monkeypatch.setattr(
        cli,
        "load_bulk_expression",
        lambda path: calls.setdefault("load", path) and bulk,
    )

    def prepare(path, frame, **kwargs):
        calls["prepare"] = (path, frame, kwargs)
        return SimpleNamespace(
            warnings=(), provenance={"n_cells_output": 8, "source_layer": "X"}
        )

    monkeypatch.setattr(cli, "prepare_single_cell", prepare)
    monkeypatch.setattr(
        cli,
        "run_decepticonx",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("must not run")),
    )

    status = cli.main(
        [
            "validate",
            "cells.h5ad",
            "bulk.tsv",
            "--cell-type-key",
            "annotation",
            "--allow-normalized-x",
        ]
    )

    assert status == 0
    assert calls["load"] == "bulk.tsv"
    assert calls["prepare"][0] == "cells.h5ad"
    assert calls["prepare"][1] is bulk
    assert calls["prepare"][2]["cell_type_key"] == "annotation"
    assert calls["prepare"][2]["allow_normalized_x"] is True
    report = json.loads(capsys.readouterr().out)
    assert report["ok"] is True
    assert report["single_cell"]["n_cells_output"] == 8


def test_check_backends_reports_status_and_optional_failure(monkeypatch, capsys) -> None:
    status_map = {
        "cibersort": True,
        "cibersort_abs": True,
        "epic": False,
        "deconrnaseq": True,
        "music": True,
    }
    runtime_map = {
        name: {
            "available": name != "epic",
            "native_available": name not in {"epic", "decepticon_fast"},
        }
        for name in ("decepticon_fast", *status_map)
    }
    monkeypatch.setattr(cli, "probe_backends", lambda: status_map)
    monkeypatch.setattr(cli, "probe_runtime", lambda: runtime_map)

    assert cli.main(["check-backends"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report == {
        "ok": False,
        "backends": status_map,
        "runtime": runtime_map,
        "all_native": False,
    }
    assert cli.main(["check-backends", "--require-all"]) == 1
    capsys.readouterr()
    assert cli.main(["check-backends", "--require-native"]) == 1


def test_cli_returns_two_for_invalid_epic_json(capsys, tmp_path) -> None:
    status = cli.main(
        [
            "run",
            "cells.h5ad",
            "bulk.tsv",
            "-o",
            str(tmp_path / "out"),
            "--epic-mrna-cell",
            "not-json",
        ]
    )

    assert status == 2
    assert "must be a JSON object" in capsys.readouterr().err
