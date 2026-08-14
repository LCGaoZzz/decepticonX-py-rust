#!/usr/bin/env python3
"""Build the strict-count AnnData input for the MuSiC2 validation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread


UPSTREAM_COMMIT = "0f8e5aa71d6fe76d346c913a823b9f1e6d6bb800"
TARGET_TYPES = ("acinar", "alpha", "beta", "delta", "ductal", "gamma")


def build_h5ad(input_dir: Path, destination: Path) -> dict[str, object]:
    """Convert R-exported Matrix Market and metadata files to AnnData."""

    matrix = mmread(input_dir / "reference_cells_by_genes.mtx").tocsr()
    matrix = matrix.astype(np.float64)
    obs = pd.read_csv(input_dir / "reference_obs.tsv", sep="\t", dtype=str)
    genes = pd.read_csv(input_dir / "reference_genes.tsv", sep="\t", dtype=str)

    required_obs = {"cell_id", "donor", "subject", "cell_type"}
    missing_obs = sorted(required_obs.difference(obs.columns))
    if missing_obs:
        raise ValueError(f"reference_obs.tsv is missing columns: {', '.join(missing_obs)}")
    if "gene" not in genes.columns:
        raise ValueError("reference_genes.tsv is missing the gene column")
    if matrix.shape != (len(obs), len(genes)):
        raise ValueError(
            f"matrix shape {matrix.shape} != metadata shape {(len(obs), len(genes))}"
        )
    if obs["cell_id"].isna().any() or genes["gene"].isna().any():
        raise ValueError("cell IDs and gene symbols must not be missing")
    if obs["cell_id"].duplicated().any() or genes["gene"].duplicated().any():
        raise ValueError("cell IDs and gene symbols must be unique")
    if matrix.data.size and not np.isfinite(matrix.data).all():
        raise ValueError("reference matrix contains non-finite values")
    if matrix.data.size and np.min(matrix.data) < 0:
        raise ValueError("reference matrix contains negative values")
    observed_types = set(obs["cell_type"])
    if observed_types != set(TARGET_TYPES):
        raise ValueError(
            "reference cell types do not match the six-type benchmark contract: "
            f"{sorted(observed_types)}"
        )

    obs = obs.set_index("cell_id")
    var = pd.DataFrame(index=pd.Index(genes["gene"], name="gene"))
    adata = ad.AnnData(X=matrix, obs=obs, var=var)

    # The official object contains a very small number of half-integer values.
    # Preserve those exact values in X and expose the documented intervention as
    # an explicit integer counts layer for the package's strict default contract.
    rounded = matrix.copy()
    rounded.data = np.rint(rounded.data).astype(np.int64)
    rounded.eliminate_zeros()
    adata.layers["counts"] = rounded
    adata.uns["source"] = {
        "repository": "https://github.com/xuranw/MuSiC",
        "commit": UPSTREAM_COMMIT,
        "file": "data/EMTABsce_healthy.rds",
        "selection": ",".join(TARGET_TYPES),
        "official_fractional_values_preserved_in_X": True,
        "rounded_count_layer_is_explicit": True,
    }

    destination.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(destination, compression="gzip")

    fractional = np.abs(matrix.data - np.rint(matrix.data)) > 1e-8
    return {
        "destination": str(destination.resolve()),
        "shape": list(adata.shape),
        "nnz": int(matrix.nnz),
        "fractional_nonzero_entries": int(fractional.sum()),
        "cell_types": obs["cell_type"].value_counts().sort_index().to_dict(),
        "donors": int(obs["donor"].nunique()),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build reference.h5ad from the portable MuSiC2 R export."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="directory produced by export_music2_official.R",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="destination reference.h5ad path",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = build_h5ad(args.input_dir, args.output)
    except (OSError, TypeError, ValueError) as exc:
        print(f"build_music2_h5ad.py: error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
