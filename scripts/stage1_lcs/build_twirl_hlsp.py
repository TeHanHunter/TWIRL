#!/usr/bin/env python3
"""Generate TWIRL-FS HLSP FITS for a sector from TGLC per-orbit HDF5 files.

Per TIC: read both orbits' HDF5 → merge by time → flux-space detrend each
aperture (Small / Primary / Large) → write a single TWIRL-FS HLSP FITS
mirroring the QLP HLSP schema.

The TGLC HDF5 must have the ``RawFlux`` / ``RawFluxError`` datasets added
by the TWIRL fork patch (branch ``twirl/preserve-negative-flux``). If
those datasets are missing, the reader falls back to a synthesized
flux-from-magnitude path with a warning — usable for sectors not yet
re-extracted, but loses the negative-cadence preservation.

Output layout mirrors QLP HLSP: ``<out_root>/<tic[0:4]>/<tic[4:8]>/...``.
The filename is ``hlsp_twirlfs_tess_ffi_*`` so the TWIRL-FS HLSPs coexist
with QLP's HLSPs in a parallel tree.
"""
from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

import numpy as np  # noqa: E402

from twirl.lightcurves.flux_detrend import (  # noqa: E402
    FluxDetrendConfig,
    flux_space_detrend_result,
)
from twirl.lightcurves.hlsp_writer import (  # noqa: E402
    HLSPTarget,
    hlsp_path,
    write_twirl_hlsp,
)
from twirl.lightcurves.tglc_h5_reader import (  # noqa: E402
    APERTURE_KEYS,
    TGLCLightCurve,
    read_tglc_h5,
)


METHOD_VERSION = "twirl-fs-v2"
ADAPTIVE_METHOD_VERSION = "twirl-fs-v2-adp03q"
CANONICAL_DETREND_CONFIG = FluxDetrendConfig(
    output_mode="subtractive",
    scale_strategy="auto",
)
DETREND_CONFIG = CANONICAL_DETREND_CONFIG


def adaptive_detrend_config(bkspace_d: float, gap_split_d: float) -> FluxDetrendConfig:
    """Return the experimental tighter-spline compare-column config."""
    return FluxDetrendConfig(
        bkspace_d=bkspace_d,
        k=CANONICAL_DETREND_CONFIG.k,
        sigma_clip=CANONICAL_DETREND_CONFIG.sigma_clip,
        max_iter=CANONICAL_DETREND_CONFIG.max_iter,
        edge_pad_d=CANONICAL_DETREND_CONFIG.edge_pad_d,
        output_mode=CANONICAL_DETREND_CONFIG.output_mode,
        scale_strategy=CANONICAL_DETREND_CONFIG.scale_strategy,
        min_scale_abs=CANONICAL_DETREND_CONFIG.min_scale_abs,
        min_scale_snr=CANONICAL_DETREND_CONFIG.min_scale_snr,
        gap_split_d=gap_split_d,
        knot_strategy="quantile",
    )


def discover_orbit_hdf5(orbit_root: Path) -> dict[int, Path]:
    """Walk a per-orbit ``ffi/cam*/ccd*/LC/*.h5`` tree → {tic: path}."""
    out: dict[int, Path] = {}
    for h5 in orbit_root.rglob("LC/*.h5"):
        try:
            tic = int(h5.stem)
        except ValueError:
            continue
        out[tic] = h5
    return out


def merge_orbits(lcs: list[TGLCLightCurve]) -> TGLCLightCurve | None:
    """Concatenate per-orbit ``TGLCLightCurve`` records into one, time-sorted.

    All apertures are merged in lockstep. The PRIMARY aperture's centroid
    is used as the global SAP_X/SAP_Y reference downstream.
    """
    lcs = [lc for lc in lcs if lc is not None]
    if not lcs:
        return None
    lcs.sort(key=lambda lc: lc.orbit)
    head = lcs[0]
    time_arr = np.concatenate([lc.time for lc in lcs])
    cad = np.concatenate([lc.cadence for lc in lcs])
    qual = np.concatenate([lc.quality for lc in lcs])
    cx = np.concatenate([lc.centroid_x for lc in lcs])
    cy = np.concatenate([lc.centroid_y for lc in lcs])
    bg = np.concatenate([lc.background for lc in lcs])
    bg_err = np.concatenate([lc.background_err for lc in lcs])
    orbitid = np.concatenate([
        np.full(len(lc.time), lc.orbit, dtype=np.int32) for lc in lcs
    ])

    order = np.argsort(time_arr, kind="stable")
    apertures = {}
    for ap in APERTURE_KEYS:
        per_ap = [lc.apertures[ap] for lc in lcs]
        rf = np.concatenate([p.raw_flux for p in per_ap])[order]
        rfe = np.concatenate([p.raw_flux_err for p in per_ap])[order]
        rm = np.concatenate([p.raw_magnitude for p in per_ap])[order]
        rme = np.concatenate([p.raw_magnitude_err for p in per_ap])[order]
        px = np.concatenate([p.centroid_x for p in per_ap])[order]
        py = np.concatenate([p.centroid_y for p in per_ap])[order]
        from twirl.lightcurves.tglc_h5_reader import TGLCAperture
        apertures[ap] = TGLCAperture(
            name=ap,
            raw_flux=rf, raw_flux_err=rfe,
            raw_magnitude=rm, raw_magnitude_err=rme,
            centroid_x=px, centroid_y=py,
            flux_was_synthesized=any(p.flux_was_synthesized for p in per_ap),
        )

    merged = TGLCLightCurve(
        tic=head.tic, sector=head.sector, orbit=-1,  # multi-orbit
        cam=head.cam, ccd=head.ccd, tmag=head.tmag,
        ra=head.ra, dec=head.dec, bjd_offset=head.bjd_offset,
        time=time_arr[order], cadence=cad[order], quality=qual[order],
        centroid_x=cx[order], centroid_y=cy[order],
        background=bg[order], background_err=bg_err[order],
        apertures=apertures, path=None,
    )
    merged.orbitid = orbitid[order]  # attach for the writer
    return merged


def detrend_apertures(
    lc: TGLCLightCurve,
    *,
    cfg: FluxDetrendConfig = DETREND_CONFIG,
) -> tuple[dict[str, np.ndarray], dict[str, dict[str, object]]]:
    """Run flux-space detrend on each aperture, return relative-flux arrays.

    Returns ``{aperture_key: (sap_rel, det_rel, det_err_rel)}`` where
    ``sap_rel`` is robust-scale-normalized raw flux (no detrend) and
    ``det_rel`` is the robust-scale subtractive cotrend residual. All
    outputs are recentered to a good-cadence median of 1 to match QLP
    HLSP conventions without using a potentially negative faint-target
    median as the residual denominator.
    """
    out = {}
    diagnostics: dict[str, dict[str, object]] = {}
    for ap_key in APERTURE_KEYS:
        ap = lc.apertures[ap_key]
        result = flux_space_detrend_result(
            lc.time,
            ap.raw_flux,
            quality=lc.quality,
            flux_err=ap.raw_flux_err,
            cfg=cfg,
        )
        diagnostics[ap_key] = {
            "n_segments": result.n_segments,
            "fit_count": result.fit_count,
            "scale_source": result.scale_source,
            "cotrend_status": result.cotrend_status,
        }
        scale = result.scale
        if not np.isfinite(scale) or scale <= 0:
            n = len(ap.raw_flux)
            out[ap_key] = (
                np.full(n, np.nan, dtype=np.float32),
                np.full(n, np.nan, dtype=np.float32),
                np.full(n, np.nan, dtype=np.float32),
            )
            continue
        sap_rel = (ap.raw_flux / scale).astype(np.float32)
        det_rel = result.det_flux
        det_err_rel = result.det_flux_err
        # flux_space_detrend_result returns subtractive residual centered ~1;
        # pin the good-cadence median to exactly 1 with an additive shift.
        m2 = np.nanmedian(det_rel[np.isfinite(det_rel) & (lc.quality == 0)])
        if np.isfinite(m2):
            det_rel = det_rel - m2 + 1.0
        out[ap_key] = (sap_rel, det_rel.astype(np.float32), det_err_rel.astype(np.float32))
    return out, diagnostics


def process_one(args: tuple[int, list[Path], Path, bool, float, float, str]) -> tuple[int, bool, str]:
    (
        tic,
        paths,
        out_root,
        include_adaptive_columns,
        adaptive_bkspace,
        adaptive_gap_split,
        method_version,
    ) = args
    try:
        lcs = [read_tglc_h5(p) for p in paths]
        merged = merge_orbits(lcs)
        if merged is None:
            return tic, False, "no_lc"
        ap_det, ap_diag = detrend_apertures(merged, cfg=CANONICAL_DETREND_CONFIG)
        extra_flux_columns = None
        extra_header = None
        if include_adaptive_columns:
            adp_cfg = adaptive_detrend_config(adaptive_bkspace, adaptive_gap_split)
            adp_det, adp_diag = detrend_apertures(merged, cfg=adp_cfg)
            adp_primary_diag = adp_diag.get("Primary", {})
            extra_flux_columns = {
                "DET_FLUX_ADP": adp_det["Primary"][1],
                "DET_FLUX_ADP_ERR": adp_det["Primary"][2],
                "DET_FLUX_ADP_SML": adp_det["Small"][1],
                "DET_FLUX_ADP_LAG": adp_det["Large"][1],
            }
            extra_header = {
                "HASADP": (True, "Includes adaptive detrend compare columns"),
                "ADPMETH": (ADAPTIVE_METHOD_VERSION, "Adaptive compare method"),
                "ADPBKSP": (float(adp_cfg.bkspace_d), "Adaptive spline spacing (d)"),
                "ADPSIG": (float(adp_cfg.sigma_clip), "Adaptive rejection threshold"),
                "ADPGAP": (float(adp_cfg.gap_split_d), "Adaptive gap split threshold"),
                "ADPKNOT": (str(adp_cfg.knot_strategy), "Adaptive knot strategy"),
                "ADPNSEG": (
                    int(adp_primary_diag.get("n_segments", -1)),
                    "Adaptive cotrend segment count",
                ),
                "ADPFIT": (
                    int(adp_primary_diag.get("fit_count", -1)),
                    "Adaptive fit cadence count",
                ),
                "ADPSCAL": (
                    str(adp_primary_diag.get("scale_source", "")),
                    "Adaptive scale source",
                ),
                "ADPCOTS": (
                    str(adp_primary_diag.get("cotrend_status", "")),
                    "Adaptive cotrend status",
                ),
            }
        target = HLSPTarget(
            tic=merged.tic, sector=merged.sector, cam=merged.cam,
            ccd=merged.ccd, tmag=merged.tmag,
            ra=merged.ra, dec=merged.dec,
        )
        sap_med = ap_det["Primary"][0]
        det_med = ap_det["Primary"][1]
        det_err = ap_det["Primary"][2]
        det_sml = ap_det["Small"][1]
        det_lag = ap_det["Large"][1]
        out_path = hlsp_path(out_root, target)
        write_twirl_hlsp(
            out_path, target,
            time_btjd=merged.time,
            cadenceno=merged.cadence.astype(np.int32),
            sap_flux=sap_med,
            det_flux=det_med,
            det_flux_err=det_err,
            quality=merged.quality.astype(np.int32),
            orbitid=getattr(merged, "orbitid", np.zeros(len(merged.time), dtype=np.int32)),
            sap_x=merged.apertures["Primary"].centroid_x.astype(np.float32),
            sap_y=merged.apertures["Primary"].centroid_y.astype(np.float32),
            sap_bkg=merged.background.astype(np.float32),
            sap_bkg_err=merged.background_err.astype(np.float32),
            det_flux_sml=det_sml,
            det_flux_lag=det_lag,
            detrend_config=CANONICAL_DETREND_CONFIG,
            detrend_diagnostics=ap_diag.get("Primary"),
            method_version=method_version,
            extra_flux_columns=extra_flux_columns,
            extra_header=extra_header,
        )
        return tic, True, "ok"
    except Exception as e:
        return tic, False, f"{type(e).__name__}: {e}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument(
        "--orbit-roots", type=Path, nargs="+", required=True,
        help="One or more orbit-<N>/ffi roots (e.g. orbit-119/ffi orbit-120/ffi)",
    )
    ap.add_argument(
        "--out-root", type=Path, required=True,
        help="Output root; per-TIC FITS are written under "
             "<out_root>/<tic[0:4]>/<tic[4:8]>/<tic[8:12]>/<tic[12:16]>/",
    )
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument(
        "--limit", type=int, default=None,
        help="Process at most N TICs (debug).",
    )
    ap.add_argument(
        "--tic", type=int, default=None,
        help="Process only this TIC (debug).",
    )
    ap.add_argument(
        "--max-failures",
        type=int,
        default=0,
        help=(
            "Allow up to this many per-TIC build failures before returning "
            "nonzero. Default 0 keeps debug runs strict."
        ),
    )
    ap.add_argument(
        "--min-success-frac",
        type=float,
        default=1.0,
        help=(
            "Minimum ok/(ok+fail) fraction required for exit status 0. "
            "Default 1.0 keeps debug runs strict."
        ),
    )
    ap.add_argument(
        "--include-adaptive-columns",
        action="store_true",
        help=(
            "Add experimental tighter-spline compare columns "
            "DET_FLUX_ADP, DET_FLUX_ADP_ERR, DET_FLUX_ADP_SML, "
            "and DET_FLUX_ADP_LAG. DET_FLUX remains the canonical v2 product."
        ),
    )
    ap.add_argument(
        "--adaptive-bkspace",
        type=float,
        default=0.3,
        help=(
            "Spline knot spacing in days for the adaptive compare columns. "
            "Default 0.3 d matches the historical tighter QLP-style spacing."
        ),
    )
    ap.add_argument(
        "--adaptive-gap-split",
        type=float,
        default=0.2,
        help=(
            "Gap threshold in days for adaptive compare columns. Default 0.2 d "
            "splits the S56 orbit-boundary gaps for the tighter 0.3 d spline."
        ),
    )
    ap.add_argument(
        "--method-version",
        default=METHOD_VERSION,
        help="Primary METHOD header value for the canonical DET_FLUX columns.",
    )
    args = ap.parse_args(argv)

    print(
        f"[build-twirl-hlsp] sector={args.sector} workers={args.workers} "
        f"method={args.method_version}"
    )
    if args.include_adaptive_columns:
        print(
            "[build-twirl-hlsp] adaptive compare columns enabled: "
            f"method={ADAPTIVE_METHOD_VERSION} "
            f"bkspace={args.adaptive_bkspace} d "
            f"gap_split={args.adaptive_gap_split} d"
        )
    per_orbit_index: list[dict[int, Path]] = []
    for orbit_root in args.orbit_roots:
        idx = discover_orbit_hdf5(Path(orbit_root))
        per_orbit_index.append(idx)
        print(f"  {orbit_root} -> {len(idx):,} h5 files")

    # Union of TICs across orbits → list of paths per TIC.
    all_tics = set().union(*[i.keys() for i in per_orbit_index])
    if args.tic is not None:
        all_tics = {args.tic} & all_tics
    work: list[tuple[int, list[Path], Path, bool, float, float, str]] = []
    for tic in sorted(all_tics):
        paths = [idx[tic] for idx in per_orbit_index if tic in idx]
        if paths:
            work.append(
                (
                    tic,
                    paths,
                    args.out_root,
                    args.include_adaptive_columns,
                    args.adaptive_bkspace,
                    args.adaptive_gap_split,
                    args.method_version,
                )
            )
    if args.limit:
        work = work[: args.limit]
    print(f"[build-twirl-hlsp] {len(work):,} TICs to process")

    t0 = time.time()
    n_ok, n_fail = 0, 0
    fail_msgs: dict[int, str] = {}
    if args.workers <= 1:
        for w in work:
            tic, ok, msg = process_one(w)
            n_ok += int(ok)
            n_fail += int(not ok)
            if not ok:
                fail_msgs[tic] = msg
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(process_one, w): w[0] for w in work}
            for i, fut in enumerate(as_completed(futures), 1):
                tic, ok, msg = fut.result()
                n_ok += int(ok)
                n_fail += int(not ok)
                if not ok:
                    fail_msgs[tic] = msg
                if i % 500 == 0:
                    elapsed = time.time() - t0
                    rate = i / elapsed
                    eta = (len(work) - i) / rate
                    print(f"  [build-twirl-hlsp] {i:,}/{len(work):,} "
                          f"(ok={n_ok} fail={n_fail}) rate={rate:.1f}/s eta={eta/60:.1f}min")

    dt = time.time() - t0
    print(f"\n[build-twirl-hlsp] done: {n_ok} ok / {n_fail} fail in {dt/60:.1f} min")
    if fail_msgs:
        print(f"  first 5 failures:")
        for tic, msg in list(fail_msgs.items())[:5]:
            print(f"    TIC {tic}: {msg}")
    total = n_ok + n_fail
    success_frac = n_ok / total if total else 0.0
    if n_fail > args.max_failures or success_frac < args.min_success_frac:
        print(
            "[build-twirl-hlsp] failure threshold exceeded: "
            f"n_fail={n_fail} max_failures={args.max_failures} "
            f"success_frac={success_frac:.5f} "
            f"min_success_frac={args.min_success_frac:.5f}"
        )
        return 1
    if n_fail:
        print(
            "[build-twirl-hlsp] continuing despite allowed failures: "
            f"n_fail={n_fail} success_frac={success_frac:.5f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
