"""Shared result and configuration models."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd


DEFAULT_METHODS = ("cibersort", "deconrnaseq", "music")
DEFAULT_REFERENCES = ("bayesprism", "music2")


@dataclass(frozen=True, order=True, slots=True)
class BranchKey:
    """A deconvolution algorithm applied to one reference template."""

    method: str
    reference: str

    @property
    def slug(self) -> str:
        return f"{self.method}__{self.reference}"


@dataclass(slots=True)
class PipelineConfig:
    """User-visible pipeline configuration."""

    cell_type_key: str = "cell_type"
    sample_key: str | None = None
    counts_layer: str = "auto"
    gene_symbol_key: str | None = None
    exclude_cell_types: tuple[str, ...] = ()
    allow_normalized_x: bool = False
    species: str = "hs"
    methods: tuple[str, ...] = DEFAULT_METHODS
    references: tuple[str, ...] = DEFAULT_REFERENCES
    cibersort_qn: bool = True
    cibersort_seed: int = 0
    threads: int = 1
    consensus_pairs: int = 2
    consensus_mode: str = "r_literal"
    strict_backends: bool = True
    allow_partial_consensus: bool = False
    epic_mrna_cell: Mapping[str, float] | None = None
    cibersort_engine: str = "rust"
    epic_backend: str = "auto"
    epic_solver: str = "auto"
    deconrnaseq_backend: str = "auto"
    music_backend: str = "auto"


@dataclass(slots=True)
class DecepticonXResult:
    """All reference, branch, consensus, and provenance outputs."""

    signatures: dict[str, pd.DataFrame]
    estimates: dict[BranchKey, pd.DataFrame]
    consensus: pd.DataFrame
    consensus_unclosed: pd.DataFrame
    consensus_closed: pd.DataFrame
    diagnostics: dict[str, Any] = field(default_factory=dict)
    timings: dict[str, float] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    def write(self, output_dir: str | Path) -> Path:
        """Write stable, human-readable outputs and return the directory."""

        root = Path(output_dir)
        if root.exists():
            if not root.is_dir() or any(root.iterdir()):
                raise FileExistsError(
                    f"output directory must not already contain files: {root}"
                )
        signatures_dir = root / "signatures"
        estimates_dir = root / "estimates"
        signatures_dir.mkdir(parents=True, exist_ok=True)
        estimates_dir.mkdir(parents=True, exist_ok=True)

        for name, frame in sorted(self.signatures.items()):
            frame.to_csv(signatures_dir / f"{name}.tsv", sep="\t")
        for key, frame in sorted(self.estimates.items()):
            frame.to_csv(estimates_dir / f"{key.slug}.tsv", sep="\t")

        self.consensus.to_csv(root / "consensus.tsv", sep="\t")
        self.consensus_unclosed.to_csv(
            root / "consensus_unclosed.tsv", sep="\t"
        )
        self.consensus_closed.to_csv(root / "consensus_closed.tsv", sep="\t")
        metadata = {
            "diagnostics": self.diagnostics,
            "timings_seconds": self.timings,
            "provenance": self.provenance,
        }
        (root / "run.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False, default=str) + "\n",
            encoding="utf-8",
        )
        return root

    @property
    def consensus_raw(self) -> pd.DataFrame:
        """Backward-compatible alias for the pre-closure consensus."""

        return self.consensus_unclosed

    @property
    def consensus_normalized(self) -> pd.DataFrame:
        """Backward-compatible alias for the row-closed consensus."""

        return self.consensus_closed


def config_as_dict(config: PipelineConfig) -> dict[str, Any]:
    """Convert a configuration to a JSON-friendly dictionary."""

    return asdict(config)
