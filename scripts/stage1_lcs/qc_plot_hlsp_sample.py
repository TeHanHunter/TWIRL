#!/usr/bin/env python3
"""Plot a random sample (or all) HLSP light curves in a sector for visual QC.

Default mode: 64 random TICs per Tmag bin (faint/bright), one PNG per bin,
8x8 grid of detrended-flux panels. Each panel shows DET_FLUX vs cadence
with QUALITY != 0 cadences masked out, titled with TIC + Tmag + per-LC MAD-RMS.

--all mode: paginate all LCs into multi-page PDF at 8x8 per page.

Usage:
  qc_plot_hlsp_sample.py \\
      --hlsp-root /pdo/users/tehan/tglc-gpu-production/hlsp_s0056 \\
      --output-dir benchmark/qc_s56_gpu \\
      --n-per-bin 64
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.hlsp import iter_hlsp_fits as _iter_hlsp_fits, read_hlsp, quality_mask, tglc_mad_error  # noqa: E402


def _read_lc(path: Path) -> dict | None:
    """Adapt twirl.io.hlsp into this script's panel dict.

    Returns ``cadence``, ``flux`` masked to good cadences, plus a fractional
    MAD-RMS computed via ``tglc_mad_error`` and divided by the median flux.
    Skips LCs with < 10 good cadences (same threshold as before).
    """
    lc = read_hlsp(path)
    if lc is None:
        print(f"[qc] read fail {path.name}")
        return None
    mask = quality_mask(lc, "DET_FLUX")
    if mask.sum() < 10:
        return None
    f = lc.flux["DET_FLUX"][mask]
    cadence = lc.cadenceno[mask]
    med = float(np.median(f))
    sigma = float(np.nanmedian(tglc_mad_error(lc, "DET_FLUX")))
    rms = sigma / med if med > 0 else np.nan
    return dict(path=path, tic=lc.tic, tmag=lc.tmag, cadence=cadence, flux=f, rms=rms)


def _bin_label(tmag: float) -> str:
    if not np.isfinite(tmag):
        return "unknown"
    if tmag < 14:
        return "bright"
    if tmag < 17:
        return "mid"
    return "faint"


def _make_grid_png(lcs: list[dict], out_path: Path, title: str) -> None:
    n = len(lcs)
    if n == 0:
        return
    cols = 8
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.0, rows * 1.5),
                              sharex=False, sharey=False)
    axes = np.atleast_2d(axes).flatten()
    for ax, lc in zip(axes, lcs):
        ax.plot(lc["cadence"], lc["flux"] / np.median(lc["flux"]), ".", ms=0.5)
        ax.set_title(f"TIC {lc['tic']} T={lc['tmag']:.1f} rms={lc['rms']:.3f}",
                     fontsize=6)
        ax.tick_params(labelsize=5)
        ax.set_yticks([])
    for ax in axes[len(lcs):]:
        ax.set_visible(False)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"[qc] wrote {out_path} ({n} LCs)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hlsp-root", type=Path, required=True,
                    help="Root of hlsp_s00XX/ tree (contains shard subdirs).")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--n-per-bin", type=int, default=64,
                    help="Random LCs per Tmag bin in default mode.")
    ap.add_argument("--all", action="store_true",
                    help="Plot every LC in a multi-page PDF (slow).")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    random.seed(args.seed)

    fits_paths = list(_iter_hlsp_fits(args.hlsp_root))
    print(f"[qc] {len(fits_paths)} HLSP FITS found")
    if not fits_paths:
        return 1

    if args.all:
        # Paginate every LC, 64 per page, sorted by TIC for reproducibility.
        page_size = 64
        for page_idx in range(0, len(fits_paths), page_size):
            chunk = fits_paths[page_idx:page_idx + page_size]
            lcs = [r for p in chunk if (r := _read_lc(p)) is not None]
            out = args.output_dir / f"qc_all_p{page_idx // page_size:04d}.png"
            _make_grid_png(lcs, out, title=f"page {page_idx // page_size}")
        return 0

    # Default: read every LC's header for Tmag, bin, sample.
    print("[qc] reading headers to bin by Tmag...")
    by_bin: dict[str, list[Path]] = {}
    for p in fits_paths:
        try:
            tmag = float(fits.getheader(p, ext=0).get("TESSMAG", np.nan))
        except Exception:
            continue
        by_bin.setdefault(_bin_label(tmag), []).append(p)
    for bin_name, paths in by_bin.items():
        sample = random.sample(paths, min(args.n_per_bin, len(paths)))
        lcs = [r for p in sample if (r := _read_lc(p)) is not None]
        _make_grid_png(lcs, args.output_dir / f"qc_{bin_name}.png",
                        title=f"sector QC sample: {bin_name} ({len(lcs)} LCs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
