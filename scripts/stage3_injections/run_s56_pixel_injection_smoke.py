#!/usr/bin/env python3
"""Run one pixel-level TGLC injection smoke test on a PDO source cutout.

This is a calibration prototype, not the production recovery grid. It injects
a finite-exposure BATMAN transit into the target's reconstructed pixel model,
then re-extracts the TGLC light curve for that TIC. The main dense S56
recovery products still use raw-aperture pre-detrend injections because they
are much cheaper.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import importlib.util
import json
import pickle
from pathlib import Path
import sys
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.search.injections import batman_transit_model, choose_observed_epoch  # noqa: E402

WD_RADIUS_REARTH = 0.013 * 109.076
WD_DENSITY_G_CM3 = 3.85e5
G_SI = 6.67430e-11


def _load_script_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_TWIRLFS = None


def _twirlfs_module() -> Any:
    global _TWIRLFS
    if _TWIRLFS is None:
        _TWIRLFS = _load_script_module(
            "build_twirl_hlsp",
            REPO_ROOT / "scripts/stage1_lightcurves/build_twirl_hlsp.py",
        )
    return _TWIRLFS


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def _maybe_int(value: Any) -> int | None:
    if np.ma.is_masked(value):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _a_over_rwd_from_period(period_d: float) -> float:
    rho_kg_m3 = WD_DENSITY_G_CM3 * 1000.0
    period_s = float(period_d) * 86400.0
    return float((G_SI * rho_kg_m3 * period_s**2 / (3.0 * np.pi)) ** (1.0 / 3.0))


def _duration_from_geometry_values(
    *,
    period_d: float,
    radius_rwd: float,
    impact_b: float,
    a_over_rwd: float,
) -> float:
    chord = float(np.sqrt(max((1.0 + radius_rwd) ** 2 - impact_b**2, 1.0e-6)))
    duration_d = float(period_d) / np.pi * chord / max(a_over_rwd, 1.0e-6)
    max_duration_min = min(180.0, 0.20 * float(period_d) * 1440.0)
    return float(np.clip(duration_d * 1440.0, 2.0, max(2.0, max_duration_min)))


def _load_pickle(path: Path) -> Any:
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def _find_target_row(source: Any, tic: int | None, *, min_edge_margin_pix: float) -> dict[str, Any]:
    designations = np.asarray(source.gaia["designation"]).astype(str)
    x_col = f"sector_{int(source.sector)}_x"
    y_col = f"sector_{int(source.sector)}_y"
    candidates: list[dict[str, Any]] = []
    for row in source.tic:
        row_tic = _maybe_int(row["TIC"])
        row_gaia = _maybe_int(row["gaia3"])
        if row_tic is None or row_gaia is None:
            continue
        if tic is not None and row_tic != int(tic):
            continue
        matches = np.flatnonzero(designations == f"Gaia DR3 {row_gaia}")
        if len(matches) != 1:
            continue
        gaia_idx = int(matches[0])
        x = float(source.gaia[x_col][gaia_idx])
        y = float(source.gaia[y_col][gaia_idx])
        tmag = float(source.gaia["tess_mag"][gaia_idx])
        margin = min(x, y, float(source.flux.shape[2] - 1) - x, float(source.flux.shape[1] - 1) - y)
        rec = {
            "tic": row_tic,
            "gaia3": row_gaia,
            "gaia_index": gaia_idx,
            "x": x,
            "y": y,
            "tess_mag": tmag,
            "edge_margin_pix": float(margin),
        }
        if tic is not None:
            if margin < min_edge_margin_pix:
                raise ValueError(
                    f"TIC {tic} is only {margin:.2f} px from the source edge; "
                    f"need >= {min_edge_margin_pix:.2f} px"
                )
            return rec
        if margin >= min_edge_margin_pix:
            candidates.append(rec)
    if tic is not None:
        raise ValueError(f"TIC {tic} was not found in source TIC/Gaia mapping")
    if not candidates:
        raise ValueError(f"no source TICs have edge margin >= {min_edge_margin_pix:.2f} px")
    candidates.sort(key=lambda rec: (abs(float(rec["tess_mag"]) - 17.5), -float(rec["edge_margin_pix"])))
    return candidates[0]


def _subset_source_and_epsf(source: Any, epsf: np.ndarray, max_cadences: int | None) -> tuple[Any, np.ndarray]:
    if max_cadences is None or max_cadences <= 0 or int(max_cadences) >= len(source.time):
        return source, epsf
    keep = np.arange(int(max_cadences), dtype=int)
    out = copy.deepcopy(source)
    out.flux = np.asarray(source.flux[keep], dtype=np.float32)
    out.time = np.asarray(source.time[keep])
    out.quality = np.asarray(source.quality[keep])
    out.cadence = np.asarray(source.cadence[keep])
    return out, np.asarray(epsf[keep], dtype=np.float64)


def _target_pixel_model(
    *,
    source: Any,
    epsf: np.ndarray,
    target: dict[str, Any],
    psf_size: int,
    oversample: int,
    edge_compression: float | None,
) -> np.ndarray:
    from tglc.epsf import make_tglc_design_matrix

    positions = np.array([[target["x"], target["y"]]], dtype=np.float64)
    flux_ratios = np.array(
        [float(source.gaia["tess_flux_ratio"][int(target["gaia_index"])])],
        dtype=np.float64,
    )
    design_matrix, _ = make_tglc_design_matrix(
        source.flux.shape[1:],
        (int(psf_size), int(psf_size)),
        int(oversample),
        positions,
        flux_ratios,
        source.mask.data,
        edge_compression,
    )
    n_psf_pixels = (int(psf_size) * int(oversample) + 1) ** 2
    target_flat = design_matrix[:, :n_psf_pixels] @ np.asarray(epsf[:, :n_psf_pixels], dtype=np.float64).T
    model = target_flat.T.reshape(source.flux.shape)
    return np.asarray(model, dtype=np.float32)


def _write_extracted_h5(
    *,
    source: Any,
    epsf: np.ndarray,
    tic: int,
    output_path: Path,
    psf_size: int,
    oversample: int,
) -> Path:
    from tglc.light_curve import generate_light_curves

    output_path.parent.mkdir(parents=True, exist_ok=True)
    written = False
    for light_curve in generate_light_curves(
        source,
        epsf,
        int(psf_size),
        int(oversample),
        tic_ids=[int(tic)],
    ):
        light_curve.write_hdf5(output_path)
        written = True
    if not written:
        raise RuntimeError(f"TGLC did not generate a light curve for TIC {tic}")
    return output_path


def _primary_raw_flux(path: Path) -> np.ndarray:
    import h5py

    with h5py.File(path, "r") as h5:
        return np.asarray(
            h5["LightCurve"]["AperturePhotometry"]["PrimaryAperture"]["RawFlux"],
            dtype=np.float64,
        )


def _quality_from_h5(path: Path) -> np.ndarray:
    import h5py

    with h5py.File(path, "r") as h5:
        return np.asarray(h5["LightCurve"]["QualityFlag"], dtype=np.int64)


def _write_twirlfs_for_h5(
    h5_path: Path,
    *,
    out_root: Path,
    include_adaptive_columns: bool,
    adaptive_bkspace: float,
    adaptive_gap_split: float,
) -> Path:
    twirlfs = _twirlfs_module()
    from twirl.lightcurves.hlsp_writer import HLSPTarget, hlsp_path, write_twirl_hlsp
    from twirl.lightcurves.tglc_h5_reader import read_tglc_h5

    lc = read_tglc_h5(h5_path)
    ap_det, ap_diag = twirlfs.detrend_apertures(lc, cfg=twirlfs.CANONICAL_DETREND_CONFIG)
    extra_flux_columns = None
    extra_header = None
    if include_adaptive_columns:
        adp_cfg = twirlfs.adaptive_detrend_config(adaptive_bkspace, adaptive_gap_split)
        adp_det, adp_diag = twirlfs.detrend_apertures(lc, cfg=adp_cfg)
        adp_primary_diag = adp_diag.get("Primary", {})
        extra_flux_columns = {
            "DET_FLUX_ADP": adp_det["Primary"][1],
            "DET_FLUX_ADP_ERR": adp_det["Primary"][2],
            "DET_FLUX_ADP_SML": adp_det["Small"][1],
            "DET_FLUX_ADP_LAG": adp_det["Large"][1],
        }
        extra_header = {
            "HASADP": (True, "Includes adaptive detrend compare columns"),
            "ADPMETH": (twirlfs.ADAPTIVE_METHOD_VERSION, "Adaptive compare method"),
            "ADPBKSP": (float(adp_cfg.bkspace_d), "Adaptive spline spacing (d)"),
            "ADPSIG": (float(adp_cfg.sigma_clip), "Adaptive rejection threshold"),
            "ADPGAP": (float(adp_cfg.gap_split_d), "Adaptive gap split threshold"),
            "ADPKNOT": (str(adp_cfg.knot_strategy), "Adaptive knot strategy"),
            "ADPNSEG": (int(adp_primary_diag.get("n_segments", -1)), "Adaptive cotrend segment count"),
            "ADPFIT": (int(adp_primary_diag.get("fit_count", -1)), "Adaptive fit cadence count"),
            "ADPSCAL": (str(adp_primary_diag.get("scale_source", "")), "Adaptive scale source"),
            "ADPCOTS": (str(adp_primary_diag.get("cotrend_status", "")), "Adaptive cotrend status"),
        }

    target = HLSPTarget(
        tic=lc.tic,
        sector=lc.sector,
        cam=lc.cam,
        ccd=lc.ccd,
        tmag=lc.tmag,
        ra=lc.ra,
        dec=lc.dec,
    )
    out_path = hlsp_path(out_root, target)
    orbitid = getattr(lc, "orbitid", np.full(len(lc.time), lc.orbit, dtype=np.int32))
    write_twirl_hlsp(
        out_path,
        target,
        time_btjd=lc.time,
        cadenceno=lc.cadence.astype(np.int32),
        sap_flux=ap_det["Primary"][0],
        det_flux=ap_det["Primary"][1],
        det_flux_err=ap_det["Primary"][2],
        quality=lc.quality.astype(np.int32),
        orbitid=np.asarray(orbitid, dtype=np.int32),
        sap_x=lc.apertures["Primary"].centroid_x.astype(np.float32),
        sap_y=lc.apertures["Primary"].centroid_y.astype(np.float32),
        sap_bkg=lc.background.astype(np.float32),
        sap_bkg_err=lc.background_err.astype(np.float32),
        det_flux_sml=ap_det["Small"][1],
        det_flux_lag=ap_det["Large"][1],
        detrend_config=twirlfs.CANONICAL_DETREND_CONFIG,
        detrend_diagnostics=ap_diag.get("Primary"),
        method_version=twirlfs.METHOD_VERSION,
        extra_flux_columns=extra_flux_columns,
        extra_header=extra_header,
    )
    return out_path


def _phase_error_min(t0_bjd: float, truth_t0_bjd: float, period_d: float) -> float:
    if not all(np.isfinite([t0_bjd, truth_t0_bjd, period_d])) or period_d <= 0:
        return float("nan")
    delta_d = ((t0_bjd - truth_t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d
    return float(delta_d * 1440.0)


def _run_bls_summary(
    hlsp_path: Path,
    *,
    apertures: tuple[str, ...],
    truth_period_d: float,
    truth_t0_bjd: float,
    truth_duration_min: float,
    n_periods: int,
    min_cadences: int,
) -> dict[str, Any]:
    from twirl.io.hlsp import read_hlsp
    from twirl.search.bls import BLSConfig, run_bls_on_lc
    from twirl.search.candidates import result_to_rows

    lc = read_hlsp(hlsp_path, columns=apertures)
    if lc is None:
        raise RuntimeError(f"read_hlsp failed: {hlsp_path}")
    cfg = BLSConfig(apertures=apertures, n_periods=int(n_periods), n_peaks=5, min_cadences=int(min_cadences))
    per_aperture: dict[str, Any] = {}
    best_row: dict[str, Any] | None = None
    recovered_apertures: list[str] = []
    for aperture in apertures:
        result = run_bls_on_lc(lc, cfg=cfg, aperture=aperture)
        rows = result_to_rows(result, run_id="pixel-injection-smoke")
        peak_rows = [r for r in rows if r.get("status") == "ok" and int(r.get("peak_rank", 0)) > 0]
        row = max(peak_rows, key=lambda r: float(r.get("sde", -np.inf))) if peak_rows else rows[0]
        period = float(row.get("period_d", np.nan))
        t0_bjd = float(row.get("t0_bjd", np.nan))
        period_rel_err = (
            abs(period - truth_period_d) / truth_period_d
            if np.isfinite(period) and np.isfinite(truth_period_d) and truth_period_d > 0
            else float("nan")
        )
        t0_phase_err_min = _phase_error_min(t0_bjd, truth_t0_bjd, truth_period_d)
        recovered = (
            row.get("status") == "ok"
            and np.isfinite(period)
            and period_rel_err <= 0.02
            and abs(t0_phase_err_min) <= max(float(truth_duration_min), 20.0)
        )
        if recovered:
            recovered_apertures.append(aperture)
        rec = {
            "status": row.get("status", ""),
            "peak_rank": row.get("peak_rank", ""),
            "period_d": row.get("period_d", ""),
            "t0_bjd": row.get("t0_bjd", ""),
            "duration_min": row.get("duration_min", ""),
            "depth": row.get("depth", ""),
            "depth_snr": row.get("depth_snr", ""),
            "sde": row.get("sde", ""),
            "period_rel_err": period_rel_err,
            "t0_phase_err_min": t0_phase_err_min,
            "recovered": bool(recovered),
        }
        per_aperture[aperture] = rec
        if best_row is None or float(row.get("sde", -np.inf)) > float(best_row.get("sde", -np.inf)):
            best_row = row
    return {
        "hlsp_path": str(hlsp_path),
        "truth_period_d": float(truth_period_d),
        "truth_t0_bjd": float(truth_t0_bjd),
        "truth_duration_min": float(truth_duration_min),
        "apertures": list(apertures),
        "recovered_apertures": recovered_apertures,
        "recovered_any": bool(recovered_apertures),
        "best_aperture": "" if best_row is None else str(best_row.get("aperture", "")),
        "best_status": "" if best_row is None else str(best_row.get("status", "")),
        "best_sde": float("nan") if best_row is None else float(best_row.get("sde", np.nan)),
        "per_aperture": per_aperture,
    }


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-pickle", type=Path, required=True)
    ap.add_argument("--epsf", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--tic", type=int, default=None, help="Target TIC. Omit to choose an interior TIC automatically.")
    ap.add_argument("--min-edge-margin-pix", type=float, default=25.0)
    ap.add_argument("--max-cadences", type=int, default=None, help="Use only the first N cadences for a fast smoke.")
    ap.add_argument("--period-d", type=float, default=0.6)
    ap.add_argument("--radius-rwd", type=float, default=1.0, help="Injected companion radius in WD radii.")
    ap.add_argument("--impact-b", type=float, default=0.2)
    ap.add_argument("--duration-min", type=float, default=None, help="Override geometry-derived duration.")
    ap.add_argument("--t0-tjd", type=float, default=None, help="Transit epoch in TJD/BTJD. Omit to choose observed epoch.")
    ap.add_argument("--min-in-transit", type=int, default=2)
    ap.add_argument("--psf-size", type=int, default=11)
    ap.add_argument("--oversample", type=int, default=2)
    ap.add_argument("--edge-compression", type=float, default=None)
    ap.add_argument("--refit-epsf", action="store_true", help="Refit the ePSF after pixel injection.")
    ap.add_argument("--no-gpu", action="store_true", help="Disable GPU ePSF refit.")
    ap.add_argument("--flux-uncertainty-power", type=float, default=1.4)
    ap.add_argument("--fit-edge-compression", type=float, default=1.0e-4)
    ap.add_argument("--build-twirlfs", action="store_true", help="Write TWIRL-FS FITS for original and injected HDF5 products.")
    ap.add_argument("--no-adaptive-columns", action="store_true", help="Do not add DET_FLUX_ADP columns to TWIRL-FS smoke FITS.")
    ap.add_argument("--adaptive-bkspace", type=float, default=0.3)
    ap.add_argument("--adaptive-gap-split", type=float, default=0.2)
    ap.add_argument("--run-bls", action="store_true", help="Run BLS on the injected TWIRL-FS smoke FITS.")
    ap.add_argument("--bls-apertures", nargs="+", default=["DET_FLUX_ADP", "DET_FLUX"])
    ap.add_argument("--bls-n-periods", type=int, default=5000)
    ap.add_argument("--bls-min-cadences", type=int, default=200)
    ap.add_argument("--overwrite", action="store_true")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"output directory is not empty: {args.out_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    source = _load_pickle(args.source_pickle)
    epsf = np.asarray(np.load(args.epsf), dtype=np.float64)
    source, epsf = _subset_source_and_epsf(source, epsf, args.max_cadences)
    target = _find_target_row(source, args.tic, min_edge_margin_pix=args.min_edge_margin_pix)
    tic = int(target["tic"])

    a_over_rwd = _a_over_rwd_from_period(args.period_d)
    impact_b = float(min(args.impact_b, 0.98 * a_over_rwd, 0.999 * (1.0 + args.radius_rwd)))
    duration_min = (
        float(args.duration_min)
        if args.duration_min is not None
        else _duration_from_geometry_values(
            period_d=args.period_d,
            radius_rwd=args.radius_rwd,
            impact_b=impact_b,
            a_over_rwd=a_over_rwd,
        )
    )
    rng = np.random.default_rng(5606)
    good_for_epoch = np.isfinite(source.time) & np.isfinite(source.flux).any(axis=(1, 2))
    if args.t0_tjd is None:
        t0_tjd, in_transit = choose_observed_epoch(
            source.time,
            period_d=args.period_d,
            duration_min=duration_min,
            rng=rng,
            quality=source.quality,
            finite_mask=good_for_epoch,
            min_in_transit=args.min_in_transit,
        )
    else:
        from twirl.search.injections import box_transit_mask

        t0_tjd = float(args.t0_tjd)
        in_transit = box_transit_mask(
            source.time,
            period_d=args.period_d,
            t0_d=t0_tjd,
            duration_min=duration_min,
        )

    target_model = _target_pixel_model(
        source=source,
        epsf=epsf,
        target=target,
        psf_size=args.psf_size,
        oversample=args.oversample,
        edge_compression=args.edge_compression,
    )
    transit_model = batman_transit_model(
        source.time,
        period_d=args.period_d,
        t0_d=t0_tjd,
        radius_rstar=args.radius_rwd,
        a_over_rstar=a_over_rwd,
        impact_b=impact_b,
        supersample_factor=7,
        exp_time_d=float(getattr(source, "exposure", 200.0)) / 86400.0,
    )

    injected_source = copy.deepcopy(source)
    injected_flux = np.asarray(source.flux, dtype=np.float32) + target_model * (transit_model[:, None, None] - 1.0)
    injected_source.flux = np.asarray(injected_flux, dtype=np.float32)

    injected_epsf = epsf
    if args.refit_epsf:
        from tglc.scripts.epsfs import fit_epsf_for_source

        injected_epsf = fit_epsf_for_source(
            injected_source,
            psf_size=int(args.psf_size),
            oversample_factor=int(args.oversample),
            edge_compression_factor=float(args.fit_edge_compression),
            flux_uncertainty_power=float(args.flux_uncertainty_power),
            use_gpu=not bool(args.no_gpu),
        )
        np.save(args.out_dir / "injected_refit_epsf.npy", injected_epsf)

    original_h5 = _write_extracted_h5(
        source=source,
        epsf=epsf,
        tic=tic,
        output_path=args.out_dir / f"{tic}_original.h5",
        psf_size=args.psf_size,
        oversample=args.oversample,
    )
    injected_h5 = _write_extracted_h5(
        source=injected_source,
        epsf=injected_epsf,
        tic=tic,
        output_path=args.out_dir / f"{tic}_pixel_injected.h5",
        psf_size=args.psf_size,
        oversample=args.oversample,
    )

    twirlfs_original = None
    twirlfs_injected = None
    bls_summary = None
    if args.build_twirlfs or args.run_bls:
        include_adaptive = not bool(args.no_adaptive_columns)
        twirlfs_original = _write_twirlfs_for_h5(
            original_h5,
            out_root=args.out_dir / "twirlfs_original",
            include_adaptive_columns=include_adaptive,
            adaptive_bkspace=args.adaptive_bkspace,
            adaptive_gap_split=args.adaptive_gap_split,
        )
        twirlfs_injected = _write_twirlfs_for_h5(
            injected_h5,
            out_root=args.out_dir / "twirlfs_injected",
            include_adaptive_columns=include_adaptive,
            adaptive_bkspace=args.adaptive_bkspace,
            adaptive_gap_split=args.adaptive_gap_split,
        )
    if args.run_bls:
        bls_summary = _run_bls_summary(
            twirlfs_injected,
            apertures=tuple(str(ap) for ap in args.bls_apertures),
            truth_period_d=float(args.period_d),
            truth_t0_bjd=float(t0_tjd + 2457000.0),
            truth_duration_min=float(duration_min),
            n_periods=int(args.bls_n_periods),
            min_cadences=int(args.bls_min_cadences),
        )
        (args.out_dir / "bls_summary.json").write_text(
            json.dumps(bls_summary, indent=2, sort_keys=True, default=_json_default) + "\n"
        )

    original_primary = _primary_raw_flux(original_h5)
    injected_primary = _primary_raw_flux(injected_h5)
    quality = _quality_from_h5(injected_h5)
    good = (quality == 0) & np.isfinite(original_primary) & np.isfinite(injected_primary)
    in_good = good & in_transit
    oot_good = good & ~in_transit
    delta = injected_primary - original_primary
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_pickle": args.source_pickle,
        "epsf": args.epsf,
        "out_dir": args.out_dir,
        "target": target,
        "injection_model": "pixel_level_batman_quadratic",
        "refit_epsf": bool(args.refit_epsf),
        "used_gpu_for_refit": bool(args.refit_epsf and not args.no_gpu),
        "n_cadences": int(len(source.time)),
        "period_d": float(args.period_d),
        "t0_tjd": float(t0_tjd),
        "duration_min": float(duration_min),
        "radius_rwd": float(args.radius_rwd),
        "radius_rearth": float(args.radius_rwd * WD_RADIUS_REARTH),
        "impact_b": float(impact_b),
        "a_over_rwd": float(a_over_rwd),
        "model_depth": float(1.0 - np.nanmin(transit_model)),
        "n_in_transit": int(np.count_nonzero(in_transit)),
        "n_good_in_transit": int(np.count_nonzero(in_good)),
        "target_model_sum_median": float(np.nanmedian(np.nansum(target_model, axis=(1, 2)))),
        "target_model_sum_in_transit_median": float(np.nanmedian(np.nansum(target_model[in_transit], axis=(1, 2))))
        if np.any(in_transit)
        else float("nan"),
        "primary_delta_in_transit_median": float(np.nanmedian(delta[in_good])) if np.any(in_good) else float("nan"),
        "primary_delta_oot_median": float(np.nanmedian(delta[oot_good])) if np.any(oot_good) else float("nan"),
        "original_h5": original_h5,
        "injected_h5": injected_h5,
        "twirlfs_original": twirlfs_original,
        "twirlfs_injected": twirlfs_injected,
        "bls_summary": bls_summary,
    }
    (args.out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, default=_json_default) + "\n")

    print("[pixel-injection-smoke] complete")
    print(f"  tic: {tic}")
    print(f"  refit_epsf: {bool(args.refit_epsf)}")
    print(f"  n_cadences: {len(source.time):,}")
    print(f"  n_good_in_transit: {manifest['n_good_in_transit']}")
    print(f"  original_h5: {original_h5}")
    print(f"  injected_h5: {injected_h5}")
    if twirlfs_injected is not None:
        print(f"  twirlfs_injected: {twirlfs_injected}")
    if bls_summary is not None:
        print(f"  bls_recovered_any: {bls_summary['recovered_any']}")
        print(f"  bls_summary: {args.out_dir / 'bls_summary.json'}")
    print(f"  manifest: {args.out_dir / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
