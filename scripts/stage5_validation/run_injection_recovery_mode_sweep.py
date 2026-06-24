#!/usr/bin/env python3
"""Run injected-light-curve BLS recovery-mode sweeps.

This is a thin orchestration layer over ``build_s56_pretriage_review_queue``:
it runs BLS on the same injected HDF5 with several aperture/grid settings,
writes one queue per setting, and summarizes strict top-1, exact top-N,
harmonic, and unmatched recovery modes.

Sweep syntax:

```
name:APER1+APER2:n_periods[:duration1+duration2+...]
```

Examples:

```
adp_5k:DET_FLUX_ADP:5000
adp_sml_5k:DET_FLUX_ADP_SML:5000
adp_short_200k:DET_FLUX_ADP:200000:1+2+3+4+5+6+8+10
```
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
import sys
from typing import Any

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
SRC_ROOT = REPO_ROOT / "src"
for path in (SCRIPT_DIR, SRC_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from build_s56_pretriage_review_queue import _finalize_queue, run_injection_bls  # noqa: E402
from summarize_injection_recovery_modes import summarize  # noqa: E402


DEFAULT_DURATIONS_MIN = (3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 30.0)


@dataclass(frozen=True)
class SweepConfig:
    name: str
    apertures: tuple[str, ...]
    n_periods: int
    durations_min: tuple[float, ...] = DEFAULT_DURATIONS_MIN


def _parse_durations(raw: str | None) -> tuple[float, ...]:
    if raw is None or raw == "":
        return DEFAULT_DURATIONS_MIN
    out = tuple(float(part) for part in raw.split("+") if part)
    if not out:
        raise ValueError("duration list cannot be empty")
    if any(value <= 0 for value in out):
        raise ValueError(f"durations must be positive: {raw}")
    return out


def parse_sweep(value: str) -> SweepConfig:
    parts = value.split(":")
    if len(parts) not in (3, 4):
        raise argparse.ArgumentTypeError(
            "sweep must be name:APER1+APER2:n_periods[:duration1+duration2+...]"
        )
    name, apertures_raw, n_periods_raw = parts[:3]
    durations_raw = parts[3] if len(parts) == 4 else None
    name = name.strip()
    apertures = tuple(part.strip() for part in apertures_raw.split("+") if part.strip())
    if not name:
        raise argparse.ArgumentTypeError("sweep name cannot be empty")
    if not apertures:
        raise argparse.ArgumentTypeError("sweep must include at least one aperture")
    try:
        n_periods = int(n_periods_raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid n_periods: {n_periods_raw}") from exc
    if n_periods < 100:
        raise argparse.ArgumentTypeError("n_periods must be >= 100")
    try:
        durations_min = _parse_durations(durations_raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return SweepConfig(name=name, apertures=apertures, n_periods=n_periods, durations_min=durations_min)


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def load_injection_ids(path: Path | None) -> tuple[str, ...] | None:
    if path is None:
        return None
    ids: list[str] = []
    seen: set[str] = set()
    for line in path.read_text().splitlines():
        item = line.strip()
        if not item or item.startswith("#"):
            continue
        # Accept either a one-column text file or a CSV-ish first field.
        item = item.split(",", 1)[0].strip()
        if item and item not in seen:
            ids.append(item)
            seen.add(item)
    return tuple(ids)


def run_sweep(
    *,
    injection_h5: Path,
    out_dir: Path,
    sweeps: tuple[SweepConfig, ...],
    survival_csv: Path | None,
    n_injections: int,
    workers: int,
    reuse: bool,
    injection_ids: tuple[str, ...] | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    out_dir.mkdir(parents=True, exist_ok=True)
    for sweep in sweeps:
        sweep_dir = out_dir / sweep.name
        sweep_dir.mkdir(parents=True, exist_ok=True)
        queue_csv = sweep_dir / "review_queue.csv"
        injected_csv = sweep_dir / "injection_bls_recoveries.csv"
        if reuse and queue_csv.exists() and injected_csv.exists():
            queue = pd.read_csv(queue_csv)
        else:
            injected = run_injection_bls(
                injection_h5=injection_h5,
                apertures=sweep.apertures,
                n_injections=n_injections,
                workers=workers,
                n_periods=sweep.n_periods,
                limit_keys=list(injection_ids) if injection_ids is not None else None,
                cfg_kwargs_extra={"durations_min": sweep.durations_min},
            )
            _write_csv(injected, injected_csv)
            queue = _finalize_queue(pd.DataFrame(), injected)
            _write_csv(queue, sweep_dir / "review_queue_pre_leo.csv")
            _write_csv(queue, queue_csv)
        primary_aperture = sweep.apertures[0]
        summary = summarize(queue_csv, survival_csv, primary_aperture, sweep_dir / "recovery_mode_summary")
        rows.append(
            {
                **asdict(sweep),
                "apertures": "+".join(sweep.apertures),
                "durations_min": "+".join(f"{value:g}" for value in sweep.durations_min),
                "injection_id_filter_n": len(injection_ids) if injection_ids is not None else 0,
                "n_rows": summary["n_rows"],
                "strict_recovered": summary["strict_recovery_status_counts"].get("bls_recovered", 0),
                "top1_recovered": summary["recovery_mode_counts"].get("bls_top1_recovered", 0),
                "exact_topn": summary["recovery_mode_counts"].get("bls_topn_recovered", 0),
                "harmonic_topn": summary["recovery_mode_counts"].get("bls_topn_harmonic_match", 0),
                "unmatched": summary["recovery_mode_counts"].get("bls_peak_mismatch", 0),
                "queue_csv": str(queue_csv),
                "summary_json": str(sweep_dir / "recovery_mode_summary" / "summary.json"),
            }
        )
    overview = pd.DataFrame(rows)
    overview.to_csv(out_dir / "sweep_overview.csv", index=False)
    return overview


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--injection-h5", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sweep", type=parse_sweep, action="append", required=True)
    parser.add_argument("--survival-csv", type=Path, default=None)
    parser.add_argument("--n-injections", type=int, default=0)
    parser.add_argument(
        "--injection-id-file",
        type=Path,
        default=None,
        help="Optional newline-delimited injection_id list. Applied before --n-injections.",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--reuse", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    overview = run_sweep(
        injection_h5=args.injection_h5,
        out_dir=args.out_dir,
        sweeps=tuple(args.sweep),
        survival_csv=args.survival_csv,
        n_injections=args.n_injections,
        workers=args.workers,
        reuse=bool(args.reuse),
        injection_ids=load_injection_ids(args.injection_id_file),
    )
    print(overview.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
