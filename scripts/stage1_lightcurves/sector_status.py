#!/usr/bin/env python3
"""Render a real-time markdown dashboard of TWIRL sector production status.

Reads:
  - $REPO/scripts/stage1_lightcurves/sector_queue.txt
  - $ROOT/markers/s<NN>_cutouts_done.flag
  - $ROOT/markers/s<NN>_done.flag
  - $ROOT/orbit-<NN>/ffi/cam*/ccd*/{source,epsf,LC}/ counts
  - $ROOT/twirl_logs/s<NN>-{gpu,prep}/orbit-<NN>_*_summary.json (wall times)
  - $ROOT/hlsp_s00<NN>/ FITS counts (when present)

Writes: $ROOT/STATUS.md (and prints to stdout if --stdout).

Usage:
  sector_status.py                     # write STATUS.md
  sector_status.py --stdout            # also print to stdout
  sector_status.py --watch 60          # refresh every 60 s in a loop
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

DEFAULT_ROOT = Path("/pdo/users/tehan/tglc-gpu-production")
DEFAULT_QUEUE = Path("/pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/sector_queue.txt")
NCCD = 16
NPATCH_PER_CCD = 196
NPATCH_PER_ORBIT = NCCD * NPATCH_PER_CCD  # 3136

# Sectors that started before sector_queue.txt was created. Tracked in the
# dashboard but never picked up by prep_worker (since they're not in the
# queue file). Format: (sector, orbit_1, orbit_2). S57 was added to
# sector_queue.txt for the finalize_worker, so it's no longer here.
PRE_QUEUE_SECTORS = [(56, 119, 120)]


def parse_queue(path: Path) -> list[tuple[int, int, int]]:
    out = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 3:
            out.append((int(parts[0]), int(parts[1]), int(parts[2])))
    return out


def count_files(d: Path, glob: str) -> int:
    if not d.exists():
        return 0
    return sum(1 for _ in d.glob(glob))


def orbit_counts(root: Path, orbit: int) -> dict[str, int]:
    """Sum source/epsf/LC across all 16 CCDs for one orbit. Fast-path skips
    the per-CCD glob when the orbit-N dir doesn't even exist (queued sector)."""
    orbit_root = root / f"orbit-{orbit}" / "ffi"
    if not orbit_root.exists():
        return {"source": 0, "epsf": 0, "lc": 0}
    src = epsf = lc = 0
    for cam in (1, 2, 3, 4):
        for ccd in (1, 2, 3, 4):
            d = orbit_root / f"cam{cam}" / f"ccd{ccd}"
            src += count_files(d / "source", "source_*.pkl")
            epsf += count_files(d / "epsf", "epsf_*.npy")
            lc += count_files(d / "LC", "*.h5")
    return {"source": src, "epsf": epsf, "lc": lc}


def hlsp_count(root: Path, sector: int, only_if_done: bool = True) -> int:
    """rglob over HLSP shards is slow (~19k files in 4-level sharded tree).
    Skip unless the sector marker says DONE. Cache the per-sector count in
    markers/s<NN>_hlsp_count.txt so subsequent dashboard refreshes are O(1).
    """
    if only_if_done and not (root / "markers" / f"s{sector:02d}_done.flag").is_file():
        return 0
    cache = root / "markers" / f"s{sector:02d}_hlsp_count.txt"
    if cache.is_file():
        try:
            return int(cache.read_text().strip())
        except Exception:
            pass
    d = root / f"hlsp_s{sector:04d}"
    if not d.exists():
        return 0
    n = sum(1 for _ in d.rglob("*.fits"))
    try:
        cache.write_text(str(n))
    except Exception:
        pass
    return n


def sector_state(root: Path, sector: int, orbit_1: int, orbit_2: int) -> dict:
    cutouts_flag = root / "markers" / f"s{sector:02d}_cutouts_done.flag"
    done_flag = root / "markers" / f"s{sector:02d}_done.flag"
    o1 = orbit_counts(root, orbit_1)
    o2 = orbit_counts(root, orbit_2)
    total_lc = o1["lc"] + o2["lc"]
    total_src = o1["source"] + o2["source"]
    state = "queued"
    if done_flag.exists():
        state = "DONE"
    elif total_lc > 0:
        # Any LC h5 means we've passed prep into the finalize side, regardless
        # of whether the cutouts_done.flag was written (S56/S57 predate the
        # marker convention).
        state = "finalizing"
    elif cutouts_flag.exists():
        state = "prep done"
    elif total_src > 0:
        state = "prepping"
    return {
        "sector": sector,
        "orbits": (orbit_1, orbit_2),
        "state": state,
        "src": o1["source"] + o2["source"],
        "epsf": o1["epsf"] + o2["epsf"],
        "lc": o1["lc"] + o2["lc"],
        "hlsp": hlsp_count(root, sector),
        "cutouts_flag_mtime": cutouts_flag.stat().st_mtime if cutouts_flag.exists() else 0,
        "done_flag_mtime": done_flag.stat().st_mtime if done_flag.exists() else 0,
    }


STATE_EMOJI = {
    "queued": "⏸",
    "prepping": "🟡",
    "prep done": "🟦",
    "finalizing": "🟠",
    "DONE": "✅",
}


def render(root: Path, queue: list[tuple[int, int, int]]) -> str:
    rows = [sector_state(root, s, o1, o2) for s, o1, o2 in queue]
    lines = []
    now = time.strftime("%Y-%m-%d %H:%M:%S %Z", time.localtime())
    lines.append(f"# TWIRL sector production status — {now}")
    lines.append("")
    n_done = sum(1 for r in rows if r["state"] == "DONE")
    n_prepped = sum(1 for r in rows if r["state"] in ("prep done", "finalizing", "DONE"))
    lines.append(f"**{n_done} / {len(rows)} sectors DONE; {n_prepped} prep complete.**")
    lines.append("")
    lines.append(f"Markers dir: `{root}/markers/`  •  HLSP roots: `{root}/hlsp_s00XX/`")
    lines.append("")
    lines.append("| state | sector | orbits | src/6272 | epsf/6272 | LC | HLSP FITS |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in rows:
        # Icon-only state column so the table aligns naturally; text legend
        # below disambiguates.
        emoji = STATE_EMOJI.get(r["state"], "?")
        sec_tag = f"S{r['sector']}"
        orbits = f"{r['orbits'][0]}+{r['orbits'][1]}"
        # Show counts for every state, including DONE, so the table preserves
        # the production trail (per-sector totals are useful when comparing
        # sectors and spotting under-recovery later).
        src_s = f"{r['src']}/{2*NPATCH_PER_ORBIT}"
        epsf_s = f"{r['epsf']}/{2*NPATCH_PER_ORBIT}"
        lc_s = str(r['lc']) if r['lc'] else "-"
        hlsp_s = str(r['hlsp']) if r['hlsp'] else "-"
        lines.append(f"| {emoji} | {sec_tag} | {orbits} | {src_s} | {epsf_s} | {lc_s} | {hlsp_s} |")
    lines.append("")
    lines.append("Legend: ⏸ queued · 🟡 prepping (cutouts in flight) · 🟦 prep done (waiting for finalize GPU) · 🟠 finalizing (ePSF/LC/detrend/HLSP) · ✅ DONE")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    ap.add_argument("--out", type=Path, help="Output markdown path. Defaults to <root>/STATUS.md.")
    ap.add_argument("--stdout", action="store_true", help="Also write to stdout.")
    ap.add_argument("--watch", type=int, default=0, help="Refresh every N seconds in a loop (0 = one shot).")
    args = ap.parse_args()

    out_path = args.out or (args.root / "STATUS.md")
    queue = PRE_QUEUE_SECTORS + parse_queue(args.queue)

    def render_once():
        text = render(args.root, queue)
        out_path.write_text(text, encoding="utf-8")
        if args.stdout:
            print(text)

    if args.watch <= 0:
        render_once()
        return 0
    while True:
        try:
            render_once()
        except Exception as e:
            print(f"[sector_status] WARN: {e}")
        time.sleep(args.watch)


if __name__ == "__main__":
    raise SystemExit(main())
