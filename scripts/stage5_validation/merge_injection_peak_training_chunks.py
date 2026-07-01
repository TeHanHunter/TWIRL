#!/usr/bin/env python3
"""Merge chunked injection peak-training tables."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
SRC_ROOT = REPO_ROOT / "src"
for path in (SCRIPT_DIR, SRC_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from build_injection_peak_training_table import summarize_peak_table  # noqa: E402


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _chunk_tables(chunk_root: Path, table_name: str) -> list[Path]:
    return sorted(chunk_root.glob(f"chunk_*/{table_name}"))


def merge_peak_chunks(
    *,
    chunk_root: Path,
    out_table: Path,
    table_name: str,
    expect_injections: int,
    id_columns: tuple[str, ...] = ("injection_id",),
) -> dict[str, Any]:
    paths = _chunk_tables(chunk_root, table_name)
    if not paths:
        raise FileNotFoundError(f"no chunk tables named {table_name!r} under {chunk_root}")

    frames: list[pd.DataFrame] = []
    source_rows: list[dict[str, Any]] = []
    for path in paths:
        frame = pd.read_csv(path)
        missing_id_cols = [col for col in id_columns if col not in frame]
        if missing_id_cols:
            raise ValueError(f"missing ID columns in {path}: {', '.join(missing_id_cols)}")
        frames.append(frame)
        unique_ids = int(frame.loc[:, list(id_columns)].drop_duplicates().shape[0])
        source_rows.append(
            {
                "chunk_table": str(path),
                "n_rows": int(len(frame)),
                "n_injections": int(frame["injection_id"].nunique()) if "injection_id" in frame else unique_ids,
                "n_unique_ids": unique_ids,
            }
        )

    merged = pd.concat(frames, ignore_index=True, sort=False)
    chunk_sources = pd.DataFrame(source_rows)
    by_unique_id = int(chunk_sources["n_unique_ids"].sum())
    actual_unique_id = int(merged.loc[:, list(id_columns)].drop_duplicates().shape[0])
    actual_injections = int(merged["injection_id"].nunique()) if "injection_id" in merged else actual_unique_id
    if by_unique_id != actual_unique_id:
        raise ValueError(
            f"chunk ID counts imply {by_unique_id}, merged unique IDs are {actual_unique_id}; "
            f"check for duplicated IDs across chunks using key {','.join(id_columns)}"
        )
    if expect_injections > 0 and actual_injections != expect_injections:
        raise ValueError(f"expected {expect_injections} injections, found {actual_injections}")

    sort_cols = [col for col in (*id_columns, "aperture", "peak_rank") if col in merged.columns]
    if sort_cols:
        merged = merged.sort_values(sort_cols, kind="stable").reset_index(drop=True)

    out_table.parent.mkdir(parents=True, exist_ok=True)
    if out_table.suffix == ".parquet":
        merged.to_parquet(out_table, index=False)
    else:
        merged.to_csv(out_table, index=False)

    summary = summarize_peak_table(merged)
    summary.update(
        {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "chunk_root": str(chunk_root),
            "out_table": str(out_table),
            "n_chunks": int(len(paths)),
            "table_name": table_name,
            "id_columns": list(id_columns),
            "n_unique_ids": actual_unique_id,
        }
    )
    summary_path = out_table.with_name(out_table.stem + "_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    chunk_sources.to_csv(out_table.with_name(out_table.stem + "_chunk_sources.csv"), index=False)
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-root", type=Path, required=True)
    parser.add_argument("--out-table", type=Path, required=True)
    parser.add_argument("--table-name", default="injection_bls_peaks.csv")
    parser.add_argument("--expect-injections", type=int, default=0)
    parser.add_argument(
        "--id-columns",
        default="injection_id",
        help="Comma-separated uniqueness key. Use search_branch,injection_id for merged branch tables.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = merge_peak_chunks(
        chunk_root=args.chunk_root,
        out_table=args.out_table,
        table_name=args.table_name,
        expect_injections=args.expect_injections,
        id_columns=tuple(col.strip() for col in args.id_columns.split(",") if col.strip()),
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
