#!/usr/bin/env python3
"""Build an explicitly qualified timing table from R and Python/Rust runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-timings", type=Path, required=True)
    parser.add_argument("--accelerated-frozen-run", type=Path, required=True)
    parser.add_argument("--accelerated-e2e-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--r-external-wall", type=float)
    parser.add_argument("--r-rss-kb", type=int)
    parser.add_argument("--accelerated-frozen-external-wall", type=float)
    parser.add_argument("--accelerated-frozen-rss-kb", type=int)
    parser.add_argument("--accelerated-e2e-external-wall", type=float)
    parser.add_argument("--accelerated-e2e-rss-kb", type=int)
    args = parser.parse_args()

    r_frame = pd.read_csv(args.r_timings, sep="\t")
    if set(r_frame["status"]) != {"ok"}:
        failed = r_frame.loc[r_frame["status"] != "ok", ["stage", "error"]]
        raise ValueError(f"R timing file contains failed stages:\n{failed}")
    r = r_frame.set_index("stage")["elapsed_seconds"].to_dict()
    frozen = json.loads(args.accelerated_frozen_run.read_text(encoding="utf-8"))[
        "timings_seconds"
    ]
    e2e = json.loads(args.accelerated_e2e_run.read_text(encoding="utf-8"))[
        "timings_seconds"
    ]

    rows: list[dict[str, object]] = []

    def add(
        stage: str,
        r_seconds: float | None,
        accelerated_seconds: float | None,
        note: str,
    ) -> None:
        speedup = None
        if (
            r_seconds is not None
            and accelerated_seconds is not None
            and accelerated_seconds > 0
        ):
            speedup = r_seconds / accelerated_seconds
        rows.append(
            {
                "stage": stage,
                "r_seconds": r_seconds,
                "accelerated_seconds": accelerated_seconds,
                "r_over_accelerated": speedup,
                "qualification": note,
            }
        )

    add(
        "reference_build_including_single_cell_preprocessing",
        sum(r[name] for name in ("ref_bayesprism", "ref_monocle3", "ref_music2_legacy")),
        e2e["single_cell_input"] + e2e["references"],
        (
            "R legacy builders versus accelerated single-cell loading, bulk-gene "
            "intersection/remainder preparation, and aggregation; compare reference "
            "semantics separately"
        ),
    )
    add(
        "cibersort_pair_three_refs",
        r["deconv_cibersort"] + r["deconv_cibersort_abs"],
        frozen["deconv_cibersort_pair"],
        "R fits relative and ABS separately; accelerated CIBERSORT deliberately shares each fit",
    )
    add(
        "epic_three_refs",
        r["deconv_epic"],
        frozen["deconv_epic"],
        "same frozen refs; solver must be checked in provenance",
    )
    add(
        "deconrnaseq_three_refs",
        r["deconv_deconrnaseq"],
        frozen["deconv_deconrnaseq"],
        "R stock wrapper includes fig=TRUE plotting; accelerated wrapper uses fig=FALSE",
    )
    add(
        "music_three_refs",
        r["music_bulk_eset"] + r["deconv_music"],
        frozen["deconv_music"],
        "same five-pseudo-donor construction; includes R bulk ExpressionSet setup",
    )
    add(
        "backend_common5_total",
        sum(
            r[name]
            for name in (
                "deconv_cibersort",
                "deconv_cibersort_abs",
                "deconv_epic",
                "deconv_deconrnaseq",
                "music_bulk_eset",
                "deconv_music",
            )
        ),
        sum(
            frozen[name]
            for name in (
                "deconv_cibersort_pair",
                "deconv_epic",
                "deconv_deconrnaseq",
                "deconv_music",
            )
        ),
        "internal elapsed; R includes method writes while accelerated writing is separate",
    )
    add(
        "external_process_wall_full_e2e",
        args.r_external_wall,
        args.accelerated_e2e_external_wall,
        "full R workflow versus full accelerated workflow; reference semantics are qualified separately",
    )
    add(
        "external_process_wall_frozen_solver_only",
        None,
        args.accelerated_frozen_external_wall,
        "frozen-reference solver process only; intentionally not compared with the full R workflow",
    )
    if args.r_rss_kb is not None or args.accelerated_e2e_rss_kb is not None:
        rows.append(
            {
                "stage": "maximum_rss_kb_full_e2e",
                "r_seconds": args.r_rss_kb,
                "accelerated_seconds": args.accelerated_e2e_rss_kb,
                "r_over_accelerated": (
                    args.r_rss_kb / args.accelerated_e2e_rss_kb
                    if args.r_rss_kb is not None
                    and args.accelerated_e2e_rss_kb is not None
                    and args.accelerated_e2e_rss_kb > 0
                    else None
                ),
                "qualification": (
                    "full-process maximum RSS from /usr/bin/time -v; peak timing and "
                    "child-process accounting limit cross-runtime interpretation"
                ),
            }
        )
    if args.accelerated_frozen_rss_kb is not None:
        rows.append(
            {
                "stage": "maximum_rss_kb_frozen_solver_only",
                "r_seconds": None,
                "accelerated_seconds": args.accelerated_frozen_rss_kb,
                "r_over_accelerated": None,
                "qualification": (
                    "frozen-reference solver-process maximum RSS only; not an "
                    "end-to-end memory comparison"
                ),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
