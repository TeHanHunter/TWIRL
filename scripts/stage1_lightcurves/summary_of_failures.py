#!/usr/bin/env python3
"""Scan a tglc-gpu-production run-log directory and report which CCDs took
multiple attempts, which stages failed, and the failure mode (OOM, NFS, etc.).

Looks for:
  - <log_dir>/<run_label>/orbit-<N>_cam<C>_ccd<D>_<stage>.log
  - <log_dir>/<run_label>/orbit-<N>_cam<C>_ccd<D>_summary.json
and:
  - On-disk file counts under tglc-data-dir/orbit-<N>/ffi/cam<C>/ccd<D>/{source,epsf,LC}.

Usage:
  summary_of_failures.py --run-log /pdo/users/tehan/tglc-gpu-production/twirl_logs/s56-gpu-rerun \\
      --tglc-data-dir /pdo/users/tehan/tglc-gpu-production --orbit 119 120
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

STAGES = ("catalogs", "cutouts", "epsfs", "lightcurves")
STAGE_FAIL_PATTERNS = [
    ("OOM", re.compile(r"OutOfMemoryError|out of memory", re.I)),
    ("NFS_IO", re.compile(r"Errno\s*5|Input/output error", re.I)),
    ("traceback", re.compile(r"Traceback \(most recent call last\)")),
    ("RETURN_CODE_NONZERO", re.compile(r"RETURN_CODE\s*[1-9]")),
]


def _detect_failure(log_path: Path) -> str | None:
    if not log_path.exists():
        return None
    try:
        text = log_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return f"READ_ERROR:{e}"
    for label, pat in STAGE_FAIL_PATTERNS:
        if pat.search(text):
            return label
    return None


def _count(d: Path, glob: str) -> int:
    if not d.exists():
        return 0
    return sum(1 for _ in d.glob(glob))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-log", type=Path, required=True,
                    help="Run-label log dir (e.g. .../twirl_logs/s56-gpu-rerun).")
    ap.add_argument("--tglc-data-dir", type=Path, required=True,
                    help="TGLC data root (e.g. /pdo/users/tehan/tglc-gpu-production).")
    ap.add_argument("--orbit", type=int, nargs="+", required=True,
                    help="Orbits to scan (e.g. 119 120).")
    args = ap.parse_args()

    rows: list[dict[str, object]] = []
    for orbit in args.orbit:
        for cam in (1, 2, 3, 4):
            for ccd in (1, 2, 3, 4):
                label = f"cam{cam}_ccd{ccd}"
                ccd_dir = args.tglc_data_dir / f"orbit-{orbit}" / "ffi" / f"cam{cam}" / f"ccd{ccd}"
                source = _count(ccd_dir / "source", "source_*.pkl")
                epsf = _count(ccd_dir / "epsf", "epsf_*.npy")
                lc = _count(ccd_dir / "LC", "*.h5")
                summary_path = args.run_log / f"orbit-{orbit}_{label}_summary.json"
                summary_exists = summary_path.is_file()
                wall_h: float | None = None
                if summary_exists:
                    try:
                        sm = json.loads(summary_path.read_text())
                        wall_h = sm.get("job_wall_seconds", 0) / 3600
                    except Exception:
                        pass
                fail_modes: dict[str, str] = {}
                for stage in STAGES:
                    log_path = args.run_log / f"orbit-{orbit}_{label}_{stage}.log"
                    fm = _detect_failure(log_path)
                    if fm and fm != "RETURN_CODE_NONZERO":
                        fail_modes[stage] = fm
                rows.append({
                    "orbit": orbit,
                    "label": label,
                    "source": source, "epsf": epsf, "lc": lc,
                    "summary": summary_exists,
                    "wall_h": wall_h,
                    "failures": fail_modes,
                })

    print(f"{'orbit':>5} {'ccd':>10} {'source':>6} {'epsf':>5} {'lc':>5} "
          f"{'summary':>7} {'wall_h':>7}  failures")
    n_fail_ccds = 0
    for r in rows:
        wall = f"{r['wall_h']:6.2f}" if r['wall_h'] is not None else "  -   "
        ok = (r['epsf'] == 196 and r['lc'] > 0)
        flag = " " if ok else "*"
        if not ok:
            n_fail_ccds += 1
        fails = ",".join(f"{s}:{m}" for s, m in r['failures'].items()) or "-"
        print(f"{r['orbit']:>5} {r['label']:>10} {r['source']:>6} {r['epsf']:>5} "
              f"{r['lc']:>5} {str(r['summary']):>7} {wall}  {flag}{fails}")
    n_total = len(rows)
    print(f"\n{n_total - n_fail_ccds}/{n_total} CCDs DONE (epsf=196 and LC>0); "
          f"{n_fail_ccds} unfinished or in-progress.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
