#!/usr/bin/env python3
"""Summarize progress and partial recall from chunked injected peak tables.

This script is intentionally read-only and stdlib-only. It is safe to run while
the chunked BLS job is active because each chunk writes its CSV and summary only
after completion.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


def _default_repo_root() -> Path:
    script_path = Path(__file__)
    if not script_path.exists():
        return Path.cwd()
    return script_path.resolve().parents[2]


REPO_ROOT = _default_repo_root()
DEFAULT_CHUNK_ROOT = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training/chunks"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "peak_training/chunk_progress"
)
DEFAULT_TOP_KS = (1, 5, 10, 20)
DEFAULT_TMAG_BINS = (
    ("<17", 0.0, 17.0),
    ("17-18", 17.0, 18.0),
    ("18-19", 18.0, 19.0),
    (">19", 19.0, 99.0),
)


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return int(default)


def _as_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(str(value).strip())
    except (TypeError, ValueError):
        return float(default)
    return out


def _parse_manifest_expected_chunks(chunk_root: Path) -> int | None:
    manifest = chunk_root.parent / "chunk_ids" / "manifest.txt"
    if not manifest.exists():
        return None
    for line in manifest.read_text(errors="replace").splitlines():
        key, sep, value = line.partition("=")
        if sep and key.strip() == "n_chunks":
            try:
                return int(value)
            except ValueError:
                return None
    return None


def _chunk_sort_key(path: Path) -> tuple[int, str]:
    name = path.name
    try:
        return int(name.split("_", 1)[1]), name
    except (IndexError, ValueError):
        return 10**9, name


def _tmag_bin(tmag: float) -> str:
    for label, lo, hi in DEFAULT_TMAG_BINS:
        if lo <= tmag < hi:
            return label
    return "unknown"


def _new_tmag_bin_counts() -> dict[str, dict[str, int]]:
    return {label: {"n": 0, **{f"top{k}": 0 for k in DEFAULT_TOP_KS}} for label, *_ in DEFAULT_TMAG_BINS}


def summarize_chunk_progress(chunk_root: Path = DEFAULT_CHUNK_ROOT) -> dict[str, Any]:
    chunk_root = Path(chunk_root)
    chunk_dirs = sorted(
        (path for path in chunk_root.glob("chunk_*") if path.is_dir()),
        key=_chunk_sort_key,
    ) if chunk_root.exists() else []
    completed = [
        path
        for path in chunk_dirs
        if (path / "injection_bls_peaks.csv").exists()
        and (path / "injection_bls_peaks_summary.json").exists()
    ]
    incomplete = [path for path in chunk_dirs if path not in set(completed)]

    injections: set[str] = set()
    first_tmag: dict[str, float] = {}
    hit_by_k: dict[int, set[str]] = {k: set() for k in DEFAULT_TOP_KS}
    match_kind_counts: dict[str, int] = {}
    row_count = 0
    candidate_row_count = 0
    signal_peak_row_count = 0
    chunk_summaries: list[dict[str, Any]] = []

    for chunk_dir in completed:
        csv_path = chunk_dir / "injection_bls_peaks.csv"
        chunk_rows = 0
        chunk_candidates = 0
        with csv_path.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                row_count += 1
                chunk_rows += 1
                injection_id = str(row.get("injection_id", "")).strip()
                if injection_id:
                    injections.add(injection_id)
                    if injection_id not in first_tmag:
                        first_tmag[injection_id] = _as_float(row.get("tmag"))
                if not _as_bool(row.get("is_candidate_peak", "true")):
                    continue
                candidate_row_count += 1
                chunk_candidates += 1
                match_kind = str(row.get("match_kind", "")).strip()
                match_kind_counts[match_kind] = match_kind_counts.get(match_kind, 0) + 1
                if not _as_bool(row.get("is_injected_signal_peak", "")):
                    continue
                signal_peak_row_count += 1
                peak_rank = _as_int(row.get("peak_rank"), default=10**9)
                for k in DEFAULT_TOP_KS:
                    if peak_rank <= k and injection_id:
                        hit_by_k[k].add(injection_id)
        chunk_summaries.append(
            {
                "chunk": chunk_dir.name,
                "rows": chunk_rows,
                "candidate_rows": chunk_candidates,
                "table": str(csv_path),
            }
        )

    n_injections = len(injections)
    recall = {}
    for k, hits in hit_by_k.items():
        n = len(hits)
        recall[f"top{k}"] = {
            "n": n,
            "denom": n_injections,
            "fraction": (n / n_injections if n_injections else None),
        }

    by_tmag = _new_tmag_bin_counts()
    by_tmag.setdefault("unknown", {"n": 0, **{f"top{k}": 0 for k in DEFAULT_TOP_KS}})
    for injection_id in sorted(injections):
        label = _tmag_bin(first_tmag.get(injection_id, float("nan")))
        if label not in by_tmag:
            by_tmag[label] = {"n": 0, **{f"top{k}": 0 for k in DEFAULT_TOP_KS}}
        by_tmag[label]["n"] += 1
        for k, hits in hit_by_k.items():
            if injection_id in hits:
                by_tmag[label][f"top{k}"] += 1

    expected_chunks = _parse_manifest_expected_chunks(chunk_root)
    payload: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "chunk_root": str(chunk_root),
        "expected_chunks": expected_chunks,
        "n_chunk_dirs": len(chunk_dirs),
        "n_completed_chunks": len(completed),
        "n_incomplete_chunks": len(incomplete),
        "incomplete_chunks_preview": [path.name for path in incomplete[:12]],
        "n_rows": row_count,
        "n_candidate_rows": candidate_row_count,
        "n_signal_peak_rows": signal_peak_row_count,
        "n_injections": n_injections,
        "recall": recall,
        "by_tmag": by_tmag,
        "match_kind_counts": dict(sorted(match_kind_counts.items())),
        "chunk_summaries": chunk_summaries,
    }
    return payload


def render_markdown(summary: dict[str, Any]) -> str:
    expected = summary.get("expected_chunks") or "?"
    n_injections = int(summary.get("n_injections", 0))
    lines = [
        "# Injected Peak Chunk Progress",
        "",
        f"Created UTC: `{summary['created_utc']}`",
        "",
        f"- completed chunks: `{summary['n_completed_chunks']}` / `{expected}`",
        f"- chunk dirs: `{summary['n_chunk_dirs']}`",
        f"- injections summarized: `{n_injections:,}`",
        f"- peak rows summarized: `{summary['n_rows']:,}`",
        f"- incomplete preview: `{', '.join(summary['incomplete_chunks_preview']) or 'none'}`",
        "",
        "## Recall",
        "",
    ]
    for label, value in summary.get("recall", {}).items():
        n = int(value.get("n", 0))
        denom = int(value.get("denom", 0))
        frac = value.get("fraction")
        frac_text = "nan" if frac is None else f"{100.0 * float(frac):.1f}%"
        lines.append(f"- {label}: `{n:,}/{denom:,}` ({frac_text})")
    lines.extend(["", "## Tmag Top-20", ""])
    for label, value in summary.get("by_tmag", {}).items():
        n = int(value.get("n", 0))
        top20 = int(value.get("top20", 0))
        frac_text = "nan" if n == 0 else f"{100.0 * top20 / n:.1f}%"
        lines.append(f"- {label}: `{top20:,}/{n:,}` ({frac_text})")
    return "\n".join(lines) + "\n"


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-root", type=Path, default=DEFAULT_CHUNK_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = summarize_chunk_progress(args.chunk_root)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    (args.out_dir / "summary.md").write_text(render_markdown(summary), encoding="utf-8")
    print(render_markdown(summary), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
