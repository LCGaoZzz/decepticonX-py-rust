from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from decepticonx import backends, references
from decepticonx.exceptions import BackendUnavailableError


ROOT = Path(__file__).resolve().parents[1]
FULL_SHA_RE = re.compile(r"(?<![0-9a-fA-F])([0-9a-fA-F]{40})(?![0-9a-fA-F])")
README_FAST_INSTALL_RE = re.compile(
    r"git\+https://github\.com/LCGaoZzz/decepticon-fast\.git@"
    r"([0-9a-fA-F]{40})(?=[\"\s])"
)


def _backend_lock() -> dict[str, object]:
    return json.loads((ROOT / "backend-lock.json").read_text(encoding="utf-8"))


def test_active_documentation_uses_every_locked_backend_revision() -> None:
    components = _backend_lock()["components"]
    assert isinstance(components, dict)
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    notice = (ROOT / "NOTICE.md").read_text(encoding="utf-8")

    for component in components.values():
        assert isinstance(component, dict)
        commit = component["commit"]
        assert isinstance(commit, str)
        assert commit in readme
        assert commit in notice

    reference_builders = components["reference_builders"]
    assert isinstance(reference_builders, dict)
    locked_commit = reference_builders["commit"]
    assert isinstance(locked_commit, str)

    readme_fast_install_shas = README_FAST_INSTALL_RE.findall(readme)
    assert readme_fast_install_shas == [locked_commit]

    notice_reference_builder_rows = [
        line
        for line in notice.splitlines()
        if line.startswith("| Reference builders |")
    ]
    assert len(notice_reference_builder_rows) == 1
    notice_reference_builder_shas = FULL_SHA_RE.findall(
        notice_reference_builder_rows[0]
    )
    assert notice_reference_builder_shas == [locked_commit]


def test_runtime_install_hints_use_locked_backend_revisions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    components = _backend_lock()["components"]
    assert isinstance(components, dict)

    module_by_component = {
        "cibersort": "python_cibersort",
        "epic": "epic_py",
        "deconrnaseq": "deconrnaseq",
        "music": "music_py",
    }
    for component_name, module_name in module_by_component.items():
        component = components[component_name]
        assert isinstance(component, dict)
        commit = component["commit"]
        assert isinstance(commit, str)
        assert FULL_SHA_RE.findall(backends._INSTALL_HINTS[module_name]) == [commit]

    def unavailable(_name: str) -> object:
        raise ImportError("test-only unavailable backend")

    monkeypatch.setattr(references, "import_module", unavailable)
    with pytest.raises(BackendUnavailableError) as error:
        references._load_fast()

    reference_builders = components["reference_builders"]
    assert isinstance(reference_builders, dict)
    locked_commit = reference_builders["commit"]
    assert isinstance(locked_commit, str)
    runtime_reference_builder_shas = FULL_SHA_RE.findall(str(error.value))
    assert runtime_reference_builder_shas == [locked_commit]
