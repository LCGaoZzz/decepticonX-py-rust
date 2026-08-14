from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "evaluate_against_truth.py"
SPEC = importlib.util.spec_from_file_location("evaluate_against_truth", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_sample_suffix_normalization_and_duplicate_gate() -> None:
    frame = pd.DataFrame({"T": [0.2, 0.8]}, index=["value...1", "value...2"])
    normalized = MODULE._canonical_sample_index(frame, label="prediction")
    assert normalized.index.tolist() == ["1", "2"]

    duplicate = pd.DataFrame({"T": [0.2, 0.8]}, index=["a...1", "b...1"])
    with pytest.raises(ValueError, match="duplicate sample IDs"):
        MODULE._canonical_sample_index(duplicate, label="prediction")
