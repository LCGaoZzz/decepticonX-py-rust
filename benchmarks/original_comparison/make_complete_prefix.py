#!/usr/bin/env python3
"""Create the mechanically complete-row prefix used by the smoke benchmark."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


KNOWN_SOURCE_SHA256 = "77ca6cf40384dbee4fe1114103b42032e5232bfc3829c22c4449090a576ab2c3"
KNOWN_OUTPUT_SHA256 = "883bf9ae63507dd16094935b82b1803ac6ed03e28da02b96da9fde5e70d24a43"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if digest(args.source) != KNOWN_SOURCE_SHA256:
        raise ValueError("source does not match the supplied truncated bulk file")
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")

    with args.source.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t", strict=True))
    if len(rows) < 2:
        raise ValueError("source has no data rows")
    expected = len(rows[0]) + 1  # R-style header omits the gene-column label.
    incomplete = [index for index, row in enumerate(rows[1:], start=1) if len(row) != expected]
    if incomplete != [len(rows) - 1] or rows[-1][0] != "AFP":
        raise ValueError("the sole incomplete row is not the expected final AFP row")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Match the audited PowerShell-derived artifact byte-for-byte: UTF-8 BOM
    # and CRLF, while preserving every parsed field verbatim.
    with args.output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\r\n")
        writer.writerows(rows[:-1])
    if digest(args.output) != KNOWN_OUTPUT_SHA256:
        raise AssertionError("derived prefix hash differs from the audited artifact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
