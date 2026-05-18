#!/usr/bin/env python3
"""Plot random TWIRL-FS S56 detrending examples from a FITS product tree."""
from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402


def _hlsp_path(root: Path, sector: int, tic: int) -> Path:
    s = f"{tic:016d}"
    return (
        root / s[0:4] / s[4:8] / s[8:12] / s[12:16]
        / f"hlsp_twirlfs_tess_ffi_s{sector:04d}-{s}_tess_v01_llc.fits"
    )


def _discover_tics(orbit_roots: list[Path]) -> list[int]:
    tics: set[int] = set()
    for root in orbit_roots:
        for h5 in root.rglob("LC/*.h5"):
            try:
                tics.add(int(h5.stem))
            except ValueError:
                continue
    return sorted(tics)


def _mad_sigma(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    med = np.nanmedian(x)
    return float(1.4826 * np.nanmedian(np.abs(x - med)))


def _finite_ylim(values: np.ndarray, pad_frac: float = 0.08) -> tuple[float, float] | None:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) == 0:
        return None
    lo, hi = np.nanpercentile(finite, [1, 99])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        return None
    pad = pad_frac * (hi - lo)
    return float(lo - pad), float(hi + pad)


def _plot_one(path: Path, out_dir: Path) -> dict[str, object] | None:
    with fits.open(path, memmap=False) as hdul:
        header = hdul[0].header
        data = hdul[1].data
        tic = int(header["TICID"])
        time = np.asarray(data["TIME"], dtype=float)
        sap = np.asarray(data["SAP_FLUX"], dtype=float)
        det = np.asarray(data["DET_FLUX"], dtype=float)
        quality = np.asarray(data["QUALITY"], dtype=int)
        q0 = quality == 0
        neg_q0 = q0 & np.isfinite(sap) & (sap < 0)
        det_mad = _mad_sigma(det[q0])
        det_med = float(np.nanmedian(det[q0]))

        fig, axes = plt.subplots(
            2, 1, figsize=(9.0, 5.4), sharex=True, constrained_layout=True
        )
        fig.patch.set_facecolor("white")
        for ax in axes:
            ax.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.7)
            ax.tick_params(labelsize=9)

        if (~q0).any():
            axes[0].scatter(
                time[~q0], sap[~q0], s=3, c="#b8b8b8", alpha=0.45, linewidths=0,
                label="flagged",
            )
            axes[1].scatter(
                time[~q0], det[~q0], s=3, c="#b8b8b8", alpha=0.45, linewidths=0,
            )
        axes[0].scatter(time[q0], sap[q0], s=3, c="#222222", alpha=0.65, linewidths=0)
        axes[1].scatter(time[q0], det[q0], s=3, c="#2459a6", alpha=0.65, linewidths=0)
        axes[0].axhline(0.0, color="#b22222", linewidth=1.0, alpha=0.8)
        axes[1].axhline(1.0, color="#b22222", linewidth=1.0, alpha=0.8)
        axes[0].set_ylabel("SAP_FLUX")
        axes[1].set_ylabel("DET_FLUX")
        axes[1].set_xlabel("BTJD")

        raw_ylim = _finite_ylim(sap[q0])
        det_ylim = _finite_ylim(det[q0])
        if raw_ylim is not None:
            axes[0].set_ylim(*raw_ylim)
        if det_ylim is not None:
            axes[1].set_ylim(*det_ylim)

        title = (
            f"TIC {tic}  T={float(header.get('TESSMAG', np.nan)):.2f}  "
            f"q0={int(q0.sum())}  neg q0={int(neg_q0.sum())}  "
            f"DET MAD={det_mad:.4f}"
        )
        axes[0].set_title(title, fontsize=11)
        out_png = out_dir / f"tic{tic:010d}_twirl_fs_v1.png"
        fig.savefig(out_png, dpi=160)
        plt.close(fig)

    return {
        "tic": tic,
        "tmag": float(header.get("TESSMAG", np.nan)),
        "path": str(path),
        "plot": str(out_png),
        "n_cadences": int(len(data)),
        "n_q0": int(q0.sum()),
        "neg_sap_q0": int(neg_q0.sum()),
        "finite_det_q0": int(np.isfinite(det[q0]).sum()),
        "det_median_q0": det_med,
        "det_mad_q0": det_mad,
        "pipeline": str(header.get("PIPELINE", "")),
        "method": str(header.get("METHOD", "")),
        "bkspace": float(header.get("BKSPACE", np.nan)),
        "sigclip": float(header.get("SIGCLIP", np.nan)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--orbit-roots", type=Path, nargs="+", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--seed", type=int, default=56)
    parser.add_argument("--min-tmag", type=float, default=18.5)
    parser.add_argument("--min-neg-q0", type=int, default=100)
    args = parser.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tics = _discover_tics(args.orbit_roots)
    rng = random.Random(args.seed)
    rng.shuffle(tics)

    rows: list[dict[str, object]] = []
    for tic in tics:
        path = _hlsp_path(args.root, args.sector, tic)
        if not path.exists():
            continue
        with fits.open(path, memmap=False) as hdul:
            header = hdul[0].header
            data = hdul[1].data
            q0 = np.asarray(data["QUALITY"], dtype=int) == 0
            sap = np.asarray(data["SAP_FLUX"], dtype=float)
            neg_q0 = int((q0 & np.isfinite(sap) & (sap < 0)).sum())
            tmag = float(header.get("TESSMAG", np.nan))
        if np.isfinite(tmag) and tmag < args.min_tmag:
            continue
        if neg_q0 < args.min_neg_q0:
            continue
        row = _plot_one(path, args.out_dir)
        if row is not None:
            rows.append(row)
        if len(rows) >= args.n:
            break

    csv_path = args.out_dir / "random_twirl_fs_v1_sample.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["tic"])
        writer.writeheader()
        writer.writerows(rows)

    md_path = args.out_dir / "summary.md"
    with md_path.open("w") as f:
        f.write("# Random TWIRL-FS v1 Detrending Preview\n\n")
        f.write(f"Seed: `{args.seed}`\n\n")
        f.write(
            f"Selection: random S56 TWIRL-FS targets with `TESSMAG >= {args.min_tmag}` "
            f"and at least `{args.min_neg_q0}` negative quality-zero `SAP_FLUX` cadences.\n\n"
        )
        f.write("| TIC | Tmag | q0 | negative SAP q0 | finite DET q0 | DET q0 MAD | PNG |\n")
        f.write("| --- | ---: | ---: | ---: | ---: | ---: | --- |\n")
        for row in rows:
            png = Path(str(row["plot"])).name
            f.write(
                f"| {row['tic']} | {float(row['tmag']):.3f} | {row['n_q0']} | "
                f"{row['neg_sap_q0']} | {row['finite_det_q0']} | "
                f"{float(row['det_mad_q0']):.4f} | [{png}]({png}) |\n"
            )

    print(f"wrote {len(rows)} plots to {args.out_dir}")
    print(csv_path)
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
