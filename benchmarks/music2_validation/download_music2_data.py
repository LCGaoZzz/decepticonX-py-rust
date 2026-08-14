#!/usr/bin/env python3
"""Download the three commit-pinned inputs for the MuSiC2 validation.

Only Python's standard library is used so data provenance can be checked
before the benchmark environment (and its optional scientific backends) is
installed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Sequence
from urllib.request import Request, urlopen


UPSTREAM_COMMIT = "0f8e5aa71d6fe76d346c913a823b9f1e6d6bb800"
RAW_ROOT = f"https://raw.githubusercontent.com/xuranw/MuSiC/{UPSTREAM_COMMIT}/data"
HERE = Path(__file__).resolve().parent
DEFAULT_AUDIT = HERE / "report" / "2026-08-15" / "audit.json"


@dataclass(frozen=True)
class Source:
    upstream_name: str
    local_name: str

    @property
    def url(self) -> str:
        return f"{RAW_ROOT}/{self.upstream_name}"


SOURCES = (
    Source("bulk-eset.rds", "music2_bulk_eset.rds"),
    Source("EMTABsce_healthy.rds", "music2_EMTABsce_healthy.rds"),
    Source("true_proportion.RData", "music2_true_proportion.RData"),
)


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of *path*."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_source_hashes(audit_path: Path) -> dict[str, str]:
    """Load and validate the three frozen source hashes from ``audit.json``."""

    try:
        payload = json.loads(audit_path.read_text(encoding="utf-8"))
        recorded = payload["source_hashes_sha256"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
        raise ValueError(f"cannot read source hashes from {audit_path}: {exc}") from exc

    expected_names = {source.local_name for source in SOURCES}
    missing = sorted(expected_names.difference(recorded))
    if missing:
        raise ValueError(
            f"{audit_path} is missing source SHA-256 entries: {', '.join(missing)}"
        )

    hashes: dict[str, str] = {}
    for name in expected_names:
        value = str(recorded[name]).lower()
        if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
            raise ValueError(f"invalid SHA-256 for {name!r} in {audit_path}")
        hashes[name] = value
    return hashes


def verify_file(path: Path, expected_sha256: str) -> str:
    """Verify one file and return its digest, raising on a mismatch."""

    actual = sha256_file(path)
    if actual != expected_sha256:
        raise ValueError(
            f"SHA-256 mismatch for {path}: expected {expected_sha256}, got {actual}"
        )
    return actual


def download_source(source: Source, output_dir: Path, expected_sha256: str) -> Path:
    """Download and atomically install one source after hash verification."""

    destination = output_dir / source.local_name
    if destination.is_file():
        verify_file(destination, expected_sha256)
        print(f"verified existing {destination}")
        return destination
    if destination.exists():
        raise ValueError(f"download destination is not a regular file: {destination}")

    request = Request(source.url, headers={"User-Agent": "decepticonx-music2-validation"})
    temporary: Path | None = None
    try:
        with urlopen(request, timeout=120) as response:  # noqa: S310 - pinned HTTPS URL
            with tempfile.NamedTemporaryFile(
                mode="wb", prefix=f".{source.local_name}.", suffix=".part",
                dir=output_dir, delete=False,
            ) as handle:
                temporary = Path(handle.name)
                shutil.copyfileobj(response, handle, length=1024 * 1024)
        verify_file(temporary, expected_sha256)
        temporary.replace(destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()

    print(f"downloaded and verified {destination}")
    return destination


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download the three official MuSiC2 benchmark objects at upstream "
            f"commit {UPSTREAM_COMMIT} and verify their frozen SHA-256 hashes."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=HERE / "data",
        help="download directory (default: benchmark-local data/)",
    )
    parser.add_argument(
        "--audit",
        type=Path,
        default=DEFAULT_AUDIT,
        help="audit.json containing the expected source_hashes_sha256 values",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        hashes = expected_source_hashes(args.audit)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for source in SOURCES:
            download_source(source, args.output_dir, hashes[source.local_name])
    except (OSError, ValueError) as exc:
        print(f"download_music2_data.py: error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
