#!/usr/bin/env python3
"""Run LEO-Vetter (TWIRL fork, WD-tuned thresholds) on the top-50
planet-regime candidates from the S56 heuristic-vetter output.

This is the integration smoke test for the wd-host-tuning branch
of LEO-Vetter-twirl. It:

1. Reads ``vetted_planet_candidates.parquet`` (top per-TIC peaks after
   the duration-envelope / period-alias / period-cluster cuts).
2. For each of the top 50 TICs, loads the HLSP FITS, builds a
   ``TCELightCurve``, fills a canonical-WD ``star`` dict, runs LEO's
   ``compute_flux_metrics``, and applies ``check_thresholds_wd``.
3. Writes per-TIC LEO vetting reports under
   ``reports/stage5_validation/leo_vetter_s56_top50/vet_reports/<class>_rank<NN>_tic<TIC>_*.pdf``
   so they can be flipped through in one sortable directory.
4. Writes ``leo_metrics.parquet`` with metrics + FA/FP labels.
5. Prints a per-TIC class summary highlighting WD 1856.

Per-cadence flux error is taken from ``twirl.io.hlsp.tglc_mad_error``
(1.4826 × MAD on the good-quality detrended flux). HLSP ``DET_FLUX_ERR``
is unused: an upstream wrapper issue produces all-NaN errors for many
S56 TICs (including WD 1856).

Canonical host: M_WD=0.6 M_sun, R_WD=0.013 R_sun, Teff=10000 K,
u1=u2=0.05 (low LDC for compact / high-gravity DA atmosphere).
Per-target host parameters via Gaia parallax + Teff is Day-2 work.
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


def _log_row(i: int, n: int, out: dict) -> None:
    tic = int(out.get("tic", -1))
    flag = " *WD1856*" if tic == WD_1856_TIC else ""
    if "error" in out and out["error"]:
        print(f"  [{i:>5}/{n}] TIC {tic:>11}  -> ERROR {out['error']}{flag}")
        return
    label = (
        "PC" if out.get("leo_PC") else
        ("FA" if out.get("leo_FA") else ("FP" if out.get("leo_FP") else "?"))
    )
    print(
        f"  [{i:>5}/{n}] TIC {tic:>11}  P={out.get('P_bls', float('nan')):6.3f}d  "
        f"SDE={out.get('sde_bls', float('nan')):7.1f}  "
        f"MES={out.get('MES', float('nan')):6.2f}  "
        f"SHP={out.get('SHP', float('nan')):5.2f}  "
        f"-> {label}{flag}"
    )


def _log_pace(i: int, n: int, t0: float) -> None:
    elapsed = time.time() - t0
    rate = i / elapsed
    eta = (n - i) / rate
    print(f"  [leo] {i}/{n} done  rate={rate:.2f}/s  eta={eta / 60:.1f}min")

# The data_local/ tree lives in the main repo (not the worktree). Use the
# canonical PycharmProjects/TWIRL path for inputs/outputs so the script works
# from either worktree or main checkout.
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path("/Users/tehan/PycharmProjects/TWIRL")
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

warnings.filterwarnings("ignore", category=RuntimeWarning)

from leo_vetter.main import TCELightCurve  # noqa: E402
from leo_vetter.parameters import derived_parameters  # noqa: E402
from leo_vetter.plots import plot_summary  # noqa: E402
from leo_vetter.wd_thresholds import (  # noqa: E402
    check_thresholds_wd,
    wd_thresholds,
)
from twirl.io.hlsp import read_hlsp, quality_mask, tglc_mad_error  # noqa: E402

WD_1856_TIC = 267574918
SECTOR = 56
import os  # noqa: E402

# All paths env-var override-able so the same script runs locally (Mac) and on PDO.
# Defaults preserve the v1 local layout.
_HLSP_ROOT = Path(os.environ.get(
    "TWIRL_HLSP_ROOT",
    str(DATA_ROOT
        / "data_local/stage2/leo_inputs/sector_0056"
        / "pdo/users/tehan/tglc-gpu-production"
        / f"hlsp_s{SECTOR:04d}"),
))
_BLS_DIR = Path(os.environ.get(
    "TWIRL_BLS_DIR",
    str(DATA_ROOT / "data_local/stage2/bls_first_pass"),
))
_LEO_OUT = Path(os.environ.get(
    "TWIRL_LEO_OUT",
    str(DATA_ROOT / "benchmark" / "leo_vetter_s56_top50"),
))
_TOP_N = int(os.environ.get("TWIRL_TOP_N", "50"))    # 0 = run all survivors
_WORKERS = int(os.environ.get("TWIRL_WORKERS", "1")) # multiprocessing pool size
_VET_REPORT_POLICY = os.environ.get(
    "TWIRL_VET_REPORT_POLICY", "all"
)  # "all" | "pc_only" | "pc_plus_top20"
HLSP_LOCAL_ROOT = _HLSP_ROOT
BLS_PARQUET = _BLS_DIR / f"sector_{SECTOR:04d}" / "vetted_planet_candidates.parquet"
BENCHMARK_DIR = _LEO_OUT
VET_REPORTS_DIR = BENCHMARK_DIR / "vet_reports"


def hlsp_path(tic: int) -> Path:
    s = f"{tic:016d}"
    return (
        HLSP_LOCAL_ROOT
        / s[0:4] / s[4:8] / s[8:12] / s[12:16]
        / f"hlsp_qlp_tess_ffi_s{SECTOR:04d}-{s}_tess_v01_llc.fits"
    )


def load_lc(path: Path):
    """Load HLSP FITS into (time_BJD, raw, flux, flux_err) arrays.

    Uses ``twirl.io.hlsp.tglc_mad_error`` unconditionally — HLSP's
    ``DET_FLUX_ERR`` column is ignored (it ships all-NaN for many S56 TICs
    due to an upstream wrapper issue). MAD-based per-target sigma is the
    TWIRL-canonical estimator.
    """
    lc = read_hlsp(path, columns=("SAP_FLUX", "DET_FLUX"))
    if lc is None:
        raise RuntimeError(f"read_hlsp failed: {path}")
    sigma = tglc_mad_error(lc, aperture="DET_FLUX")
    keep = quality_mask(lc, "DET_FLUX") & np.isfinite(lc.flux["SAP_FLUX"]) & np.isfinite(sigma)
    return (
        lc.time[keep] + 2457000.0,
        lc.flux["SAP_FLUX"][keep].astype(float),
        lc.flux["DET_FLUX"][keep].astype(float),
        sigma[keep].astype(float),
    )


# Canonical WD parameters — see module docstring.
_RHO_WD_G_CM3 = 3.85e5  # 0.6 M_sun / (4/3 pi (0.013 R_sun)^3)
_CANONICAL_WD_STAR = {
    "rad": 0.013, "e_rad": 0.002,
    "mass": 0.6,  "e_mass": 0.05,
    "Teff": 10000.0, "e_Teff": 1000.0,
    "rho": _RHO_WD_G_CM3,  # used only by parameters.get_q, which we override
    "u1": 0.05, "u2": 0.05,
    "Rs": 0.013,  # picked up by unphysical_duration_wd via metrics["Rs"]
}


def wd_star_for(tic: int) -> dict:
    s = dict(_CANONICAL_WD_STAR)
    s["id"] = int(tic)
    s["tic"] = int(tic)
    return s


def class_label(fa: bool, fp: bool) -> str:
    if (not fa) and (not fp):
        return "PC"
    if fa:
        return "FA"
    return "FP"


def render_vet_report(
    tlc: TCELightCurve,
    star: dict,
    rank: int,
    label: str,
    tmag: float,
    P: float,
    out_dir: Path,
) -> Path | None:
    """Render LEO's per-TCE summary plot to a file named for fast triage:
        <CLASS>_rank<NN>_tic<TIC>_T<MAG>_P<P>d.pdf
    Sorting the directory by name groups all PCs first, then FAs, FPs.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fn = (
        f"{label}_rank{rank:02d}_tic{int(tlc.tic):010d}"
        f"_T{tmag:05.2f}_P{P:08.4f}d.pdf"
    )
    out_path = out_dir / fn
    try:
        plot_summary(tlc, star, save_fig=True, save_file=str(out_path))
        return out_path
    except Exception as e:
        return None


def vet_one(row: pd.Series, rank: int) -> dict:
    tic = int(row["tic"])
    P = float(row["period_d"])
    T0 = float(row["t0_bjd"])
    dur_d = float(row["duration_min"]) / (24.0 * 60.0)
    tmag = float(row["tmag"])
    path = hlsp_path(tic)
    if not path.exists():
        return {"tic": tic, "error": "no_hlsp"}
    try:
        time, raw, flux, err = load_lc(path)
        if len(time) < 200:
            return {"tic": tic, "error": "few_cad", "n_cad": len(time)}
        tlc = TCELightCurve(tic, time, raw, flux, err, P, T0, dur_d)
        star = wd_star_for(tic)
        tlc.compute_flux_metrics(star, cap_b=False)
        # Inject Rs so unphysical_duration_wd's chord formula picks it up.
        tlc.metrics["Rs"] = _CANONICAL_WD_STAR["Rs"]
        fa = bool(check_thresholds_wd(tlc.metrics, "FA"))
        fp = bool(check_thresholds_wd(tlc.metrics, "FP"))
        label = class_label(fa, fp)
        # Render PDF only per policy. "all" = every TIC; "pc_only" = PCs only;
        # "pc_plus_top20" = PCs + top-20 of each non-PC class (handled
        # post-hoc by the driver, so vet_one always renders PCs at minimum).
        should_render = (
            _VET_REPORT_POLICY == "all"
            or (label == "PC")
            or (_VET_REPORT_POLICY == "pc_plus_top20" and rank <= 20)
        )
        if should_render:
            render_vet_report(tlc, star, rank, label, tmag, P, VET_REPORTS_DIR)
        m = dict(tlc.metrics)
        m["leo_FA"] = fa
        m["leo_FP"] = fp
        m["leo_PC"] = (not fa) and (not fp)
        m["leo_class"] = label
        return m
    except Exception as e:
        return {"tic": tic, "error": f"{type(e).__name__}: {e}"}


def _vet_one_row(row_dict: dict) -> dict:
    """Multiprocessing-friendly entry. Receives a plain dict (not pd.Series)."""
    row = pd.Series(row_dict)
    out = vet_one(row, rank=int(row["planet_rank"]))
    out["rank_in"] = int(row["planet_rank"])
    out["P_bls"] = float(row["period_d"])
    out["dur_bls_min"] = float(row["duration_min"])
    out["sde_bls"] = float(row["sde_max"])
    out["tmag"] = float(row["tmag"])
    out["n_apertures_agree"] = int(row["n_apertures_agree"])
    return out


def main() -> int:
    print(f"[leo] reading {BLS_PARQUET}")
    df = pd.read_parquet(BLS_PARQUET)
    if _TOP_N > 0:
        df = df.head(_TOP_N).copy()
    else:
        df = df.copy()
    df = df.reset_index(drop=True)
    wd_match = df.index[df["tic"] == WD_1856_TIC]
    wd_msg = f" (WD 1856 at rank {int(wd_match[0]) + 1})" if len(wd_match) else ""
    print(f"[leo] {len(df):,} TICs to vet  workers={_WORKERS}  "
          f"vet_report_policy={_VET_REPORT_POLICY}{wd_msg}")

    BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)
    VET_REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    rows_in = df.to_dict("records")
    rows: list[dict] = []
    t0 = time.time()
    if _WORKERS <= 1:
        for i, r in enumerate(rows_in, 1):
            out = _vet_one_row(r)
            rows.append(out)
            _log_row(i, len(rows_in), out)
            if i % 100 == 0:
                _log_pace(i, len(rows_in), t0)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=_WORKERS) as ex:
            futures = [ex.submit(_vet_one_row, r) for r in rows_in]
            for i, fut in enumerate(as_completed(futures), 1):
                out = fut.result()
                rows.append(out)
                _log_row(i, len(rows_in), out)
                if i % 100 == 0:
                    _log_pace(i, len(rows_in), t0)

    metrics = pd.DataFrame(rows)
    out_path = BENCHMARK_DIR / "leo_metrics.parquet"
    metrics.to_parquet(out_path, compression="zstd")
    print(f"\n[leo] wrote {out_path}  ({len(metrics)} rows, {len(metrics.columns)} cols)")
    print(f"[leo] vet reports in {VET_REPORTS_DIR}")

    if "leo_PC" in metrics.columns:
        n_pc = int(metrics["leo_PC"].sum())
        n_fa = int(metrics["leo_FA"].sum())
        n_fp = int(metrics["leo_FP"].sum())
        n_err = int(metrics["error"].notna().sum()) if "error" in metrics else 0
        print(f"\n[leo] class tally:")
        print(f"  PC (planet candidate):       {n_pc:>3}")
        print(f"  FA (false alarm):            {n_fa:>3}")
        print(f"  FP (astrophysical FP):       {n_fp:>3}")
        print(f"  ERROR (LC/fit failure):      {n_err:>3}")
        wd = metrics[metrics["tic"] == WD_1856_TIC]
        if len(wd):
            r = wd.iloc[0]
            label = "PC" if r.get("leo_PC") else ("FA" if r.get("leo_FA") else "FP")
            print(f"\n[wd1856] LEO label: {label}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
