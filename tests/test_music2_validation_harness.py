from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "benchmarks" / "music2_validation"


def _load_downloader():
    path = HARNESS / "download_music2_data.py"
    spec = importlib.util.spec_from_file_location("music2_validation_downloader", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_download_manifest_is_commit_pinned_and_matches_frozen_audit() -> None:
    downloader = _load_downloader()
    audit_path = HARNESS / "report" / "2026-08-15" / "audit.json"
    recorded = json.loads(audit_path.read_text(encoding="utf-8"))[
        "source_hashes_sha256"
    ]
    expected = downloader.expected_source_hashes(audit_path)

    assert set(expected) == {source.local_name for source in downloader.SOURCES}
    assert expected == {name: recorded[name] for name in expected}
    for source in downloader.SOURCES:
        assert downloader.UPSTREAM_COMMIT in source.url
        assert source.url.startswith("https://raw.githubusercontent.com/xuranw/MuSiC/")


def test_download_hash_verification_rejects_changed_content(tmp_path: Path) -> None:
    downloader = _load_downloader()
    payload = tmp_path / "object.rds"
    payload.write_bytes(b"commit-pinned test object")
    expected = hashlib.sha256(payload.read_bytes()).hexdigest()

    assert downloader.verify_file(payload, expected) == expected
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        downloader.verify_file(payload, "0" * 64)


def test_source_distribution_prunes_local_validation_data_and_work() -> None:
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8").splitlines()
    assert "prune benchmarks/music2_validation/data" in manifest
    assert "prune benchmarks/music2_validation/work" in manifest


def test_analyzer_refuses_to_write_inside_the_source_run(tmp_path: Path) -> None:
    path = HARNESS / "analyze_music2_results.py"
    spec = importlib.util.spec_from_file_location("music2_validation_analyzer", path)
    assert spec is not None and spec.loader is not None
    analyzer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analyzer)

    run_dir = tmp_path / "run"
    with pytest.raises(ValueError, match="outside the source run directory"):
        analyzer.analyze(
            source_data_dir=tmp_path / "data",
            prepared_dir=tmp_path / "prepared",
            run_dir=run_dir,
            output_dir=run_dir / "analysis",
            corrected_output_dir=tmp_path / "corrected",
            package_commit="0" * 40,
        )


@pytest.mark.parametrize(
    "script",
    [
        "download_music2_data.py",
        "build_music2_h5ad.py",
        "analyze_music2_results.py",
    ],
)
def test_python_harness_cli_help(script: str) -> None:
    completed = subprocess.run(
        [sys.executable, str(HARNESS / script), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "usage:" in completed.stdout.lower()
