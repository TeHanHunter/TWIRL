"""Build native raw/error plus ADP inputs for the S56 harmonic CNN.

The raw TGLC tree remains on PDO.  ``export_tglc_raw_sources`` writes a compact
host-only file there.  ``build_raw_pair_export`` combines that file with ADP
and injection products on ORCD, including the injection-aware uncertainty
channel and shared-grid BLS periodograms.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.lightcurves.tglc_h5_reader import read_tglc_h5
from twirl.search.candidates import compute_sde
from twirl.vetting.harmonic_inputs import (
    CHRONOLOGY_SMALL_CHANNELS,
    CHRONOLOGY_SUPPLEMENTAL_CHANNELS,
    HARMONIC_VIEW_CHANNELS,
    NATIVE_DATASETS,
    PERIODOGRAM_DATASETS,
    PERIODOGRAM_CHANNELS,
    RAW_PAIR_CONTRACT_VERSION,
    injected_raw_uncertainty,
    native_group_path,
)


RAW_SOURCE_CONTRACT_VERSION = "s56_tglc_raw_pair_v1"
DEFAULT_BLS_PERIODS = 4096
DEFAULT_BLS_PERIOD_RANGE_D: tuple[float, float] = (0.12, 15.0)
DEFAULT_BLS_DURATIONS_MIN: tuple[float, ...] = (
    3.0,
    4.0,
    5.0,
    6.0,
    8.0,
    10.0,
    13.0,
    16.0,
    20.0,
    30.0,
)


def _absolute_bjd(time: np.ndarray) -> np.ndarray:
    values = np.asarray(time, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return values + 2457000.0 if finite.size and float(np.nanmedian(finite)) < 1.0e5 else values


def _write_dataset(group: Any, name: str, values: np.ndarray) -> None:
    group.create_dataset(name, data=np.asarray(values), compression="lzf", shuffle=True)


def discover_tglc_paths(
    *,
    raw_root: Path,
    tic: int,
    camera: int | None,
    ccd: int | None,
    orbits: Sequence[int] = (119, 120),
) -> list[Path]:
    detectors = (
        [(int(camera), int(ccd))]
        if camera is not None and ccd is not None
        else [(camera_value, ccd_value) for camera_value in range(1, 5) for ccd_value in range(1, 5)]
    )
    paths = [
        Path(raw_root)
        / f"orbit-{int(orbit)}"
        / "ffi"
        / f"cam{camera_value}"
        / f"ccd{ccd_value}"
        / "LC"
        / f"{int(tic)}.h5"
        for orbit in orbits
        for camera_value, ccd_value in detectors
    ]
    return [path for path in paths if path.exists()]


def compact_adp_detector_lookup(
    compact_adp_h5: Path,
    *,
    tics: Sequence[int],
) -> dict[int, tuple[int, int]]:
    """Read detector coordinates from the exact ADP product used by the model."""

    import h5py

    lookup: dict[int, tuple[int, int]] = {}
    with h5py.File(Path(compact_adp_h5), "r") as h5:
        for tic in tics:
            group_path = f"targets/{int(tic):016d}"
            if group_path not in h5:
                continue
            group = h5[group_path]
            try:
                camera = int(group.attrs["camera"])
                ccd = int(group.attrs["ccd"])
            except (KeyError, TypeError, ValueError):
                continue
            if camera not in range(1, 5) or ccd not in range(1, 5):
                raise ValueError(
                    f"compact ADP detector is invalid for TIC {int(tic)}: cam={camera} ccd={ccd}"
                )
            lookup[int(tic)] = (camera, ccd)
    return lookup


def merge_tglc_raw_paths(paths: Sequence[Path]) -> dict[str, np.ndarray]:
    """Read, concatenate, sort, and deduplicate native raw TGLC cadences."""

    if not paths:
        raise ValueError("no TGLC paths supplied")
    curves = [read_tglc_h5(path) for path in paths]
    for curve in curves:
        for aperture in ("Small", "Primary"):
            if aperture not in curve.apertures:
                raise ValueError(f"{curve.path} is missing {aperture} aperture")
            if curve.apertures[aperture].flux_was_synthesized:
                raise ValueError(f"{curve.path} has synthesized rather than native RawFlux")
    payload = {
        "time": np.concatenate([curve.time_bjd() for curve in curves]).astype(np.float64),
        "cadenceno": np.concatenate([curve.cadence for curve in curves]).astype(np.int64),
        "orbitid": np.concatenate(
            [np.full(len(curve.time), curve.orbit, dtype=np.int32) for curve in curves]
        ),
        "quality": np.concatenate([curve.quality for curve in curves]).astype(np.int32),
        "raw_flux_small": np.concatenate(
            [curve.apertures["Small"].raw_flux for curve in curves]
        ).astype(np.float64),
        "raw_flux_err_small": np.concatenate(
            [curve.apertures["Small"].raw_flux_err for curve in curves]
        ).astype(np.float64),
        "raw_flux_primary": np.concatenate(
            [curve.apertures["Primary"].raw_flux for curve in curves]
        ).astype(np.float64),
        "raw_flux_err_primary": np.concatenate(
            [curve.apertures["Primary"].raw_flux_err for curve in curves]
        ).astype(np.float64),
    }
    order = np.lexsort((payload["time"], payload["cadenceno"]))
    for name in payload:
        payload[name] = payload[name][order]
    _, first = np.unique(payload["cadenceno"], return_index=True)
    first = np.sort(first)
    for name in payload:
        payload[name] = payload[name][first]
    return payload


def align_raw_by_cadence(
    raw: Mapping[str, np.ndarray],
    *,
    cadenceno: np.ndarray,
    time: np.ndarray,
    max_time_error_s: float = 2.0,
) -> dict[str, np.ndarray]:
    """Align native raw arrays to an ADP cadence list, requiring exact coverage."""

    raw_cadence = np.asarray(raw["cadenceno"], dtype=np.int64)
    requested = np.asarray(cadenceno, dtype=np.int64)
    if len(np.unique(raw_cadence)) != len(raw_cadence):
        raise ValueError("raw source contains duplicate cadences")
    lookup = {int(value): index for index, value in enumerate(raw_cadence)}
    missing = [int(value) for value in requested if int(value) not in lookup]
    if missing:
        raise ValueError(f"raw source is missing {len(missing)} ADP cadences; first={missing[:5]}")
    index = np.asarray([lookup[int(value)] for value in requested], dtype=np.int64)
    raw_time = _absolute_bjd(np.asarray(raw["time"], dtype=float)[index])
    adp_time = _absolute_bjd(np.asarray(time, dtype=float))
    delta_s = np.abs(raw_time - adp_time) * 86400.0
    finite = np.isfinite(delta_s)
    if np.any(finite) and float(np.nanmax(delta_s[finite])) > float(max_time_error_s):
        raise ValueError(f"raw/ADP cadence timestamps differ by up to {np.nanmax(delta_s):.3f} s")
    aligned = {name: np.asarray(values)[index] for name, values in raw.items()}
    aligned["_time_delta_s"] = delta_s
    return aligned


def shared_bls_spectrum(
    *,
    time: np.ndarray,
    quality: np.ndarray,
    flux: np.ndarray,
    periods: np.ndarray,
    durations_min: Sequence[float] = DEFAULT_BLS_DURATIONS_MIN,
) -> dict[str, np.ndarray]:
    """Compute production-style BLS power and SDE on a shared grid."""

    from astropy.timeseries import BoxLeastSquares

    time = np.asarray(time, dtype=np.float64)
    quality = np.asarray(quality)
    flux = np.asarray(flux, dtype=np.float64)
    periods = np.asarray(periods, dtype=np.float64)
    invalid = {
        "power": np.full(len(periods), np.nan, dtype=np.float32),
        "sde": np.full(len(periods), np.nan, dtype=np.float32),
    }
    good = (quality == 0) & np.isfinite(time) & np.isfinite(flux)
    if np.count_nonzero(good) < 200:
        return invalid
    t = time[good]
    values = flux[good]
    median = float(np.nanmedian(values))
    if not np.isfinite(median) or abs(median) < 1.0e-12:
        return invalid
    values = values / median
    mad = float(np.nanmedian(np.abs(values - 1.0)))
    if np.isfinite(mad) and mad > 0:
        keep = (values - 1.0) <= 5.0 * 1.4826 * mad
        if np.count_nonzero(keep) >= 200:
            t = t[keep]
            values = values[keep]
    durations_d = np.asarray(durations_min, dtype=float) / 1440.0
    try:
        spectrum = BoxLeastSquares(t, values).power(periods, durations_d, oversample=1)
    except Exception:
        return invalid
    power = np.asarray(spectrum.power, dtype=np.float32)
    return {"power": power, "sde": compute_sde(power.astype(float)).astype(np.float32)}


def shared_bls_periodogram(
    *,
    time: np.ndarray,
    quality: np.ndarray,
    flux: np.ndarray,
    periods: np.ndarray,
    durations_min: Sequence[float] = DEFAULT_BLS_DURATIONS_MIN,
) -> np.ndarray:
    """Backward-compatible shared-grid SDE helper."""

    return shared_bls_spectrum(
        time=time,
        quality=quality,
        flux=flux,
        periods=periods,
        durations_min=durations_min,
    )["sde"]


def _periodogram_payload(
    *,
    time: np.ndarray,
    quality: np.ndarray,
    small: np.ndarray,
    primary: np.ndarray,
    n_periods: int = DEFAULT_BLS_PERIODS,
) -> dict[str, np.ndarray]:
    periods = np.geomspace(
        DEFAULT_BLS_PERIOD_RANGE_D[0], DEFAULT_BLS_PERIOD_RANGE_D[1], int(n_periods)
    ).astype(np.float32)
    small_spectrum = shared_bls_spectrum(
        time=time, quality=quality, flux=small, periods=periods
    )
    primary_spectrum = shared_bls_spectrum(
        time=time, quality=quality, flux=primary, periods=periods
    )
    return {
        "bls_log_period_grid": np.log10(periods).astype(np.float32),
        "bls_power_small": small_spectrum["power"],
        "bls_sde_small": small_spectrum["sde"],
        "bls_power_primary": primary_spectrum["power"],
        "bls_sde_primary": primary_spectrum["sde"],
    }


def export_tglc_raw_sources(
    *,
    training_table: Path,
    raw_root: Path,
    out_h5: Path,
    orbits: Sequence[int] = (119, 120),
    compact_adp_h5: Path | None = None,
) -> dict[str, Any]:
    """Write the compact raw/error host file on the TGLC filesystem."""

    import h5py

    rows = pd.read_csv(training_table, low_memory=False)
    active = (
        _bool_column(rows, "morphology_include_v1")
        | _bool_column(rows, "preserve_include_v1")
        | _bool_column(rows, "harmonic_include_v1")
    )
    rows = rows.loc[active]
    required = rows.loc[:, ["tic", "cam", "ccd"]].copy()
    required["_detector_known"] = required[["cam", "ccd"]].notna().all(axis=1)
    required = required.sort_values("_detector_known", ascending=False).drop_duplicates("tic")
    detector_lookup = (
        compact_adp_detector_lookup(
            compact_adp_h5,
            tics=required["tic"].astype(int).tolist(),
        )
        if compact_adp_h5 is not None
        else {}
    )
    out_h5 = Path(out_h5)
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_h5.with_suffix(out_h5.suffix + ".tmp")
    failures_path = out_h5.with_suffix(".failures.csv")
    failures: list[dict[str, Any]] = []
    written = 0
    with h5py.File(temporary, "w") as output:
        output.attrs["contract_version"] = RAW_SOURCE_CONTRACT_VERSION
        output.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        output.attrs["training_table"] = str(training_table)
        output.attrs["raw_root"] = str(raw_root)
        output.attrs["compact_adp_h5"] = str(compact_adp_h5 or "")
        output.attrs["orbits"] = json.dumps([int(value) for value in orbits])
        output.attrs["time_system"] = "BJD"
        targets = output.create_group("targets")
        for _, row in required.iterrows():
            tic = int(row["tic"])
            compact_detector = detector_lookup.get(tic)
            camera = (
                compact_detector[0]
                if compact_detector is not None
                else int(row["cam"])
                if pd.notna(row["cam"])
                else None
            )
            ccd = (
                compact_detector[1]
                if compact_detector is not None
                else int(row["ccd"])
                if pd.notna(row["ccd"])
                else None
            )
            paths = discover_tglc_paths(
                raw_root=raw_root,
                tic=tic,
                camera=camera,
                ccd=ccd,
                orbits=orbits,
            )
            try:
                detectors = {(path.parents[2].name, path.parents[1].name) for path in paths}
                if len(detectors) > 1:
                    raise ValueError(f"TIC {tic} resolved to multiple detectors: {sorted(detectors)}")
                payload = merge_tglc_raw_paths(paths)
            except Exception as exc:
                failures.append({"tic": tic, "error": str(exc), "paths": "|".join(map(str, paths))})
                continue
            group = targets.create_group(f"{tic:016d}")
            group.attrs["tic"] = tic
            group.attrs["camera"] = int(camera) if camera is not None else -1
            group.attrs["ccd"] = int(ccd) if ccd is not None else -1
            group.attrs["detector_source"] = (
                "compact_adp_attrs" if compact_detector is not None else "training_table_or_discovery"
            )
            group.attrs["source_paths"] = json.dumps([str(path) for path in paths])
            for name, values in payload.items():
                _write_dataset(group, name, values)
            written += 1
            if written % 100 == 0:
                print(f"[raw-source-export] {written}/{len(required)} hosts", flush=True)
    if failures:
        pd.DataFrame(failures).to_csv(failures_path, index=False)
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"raw source export failed for {len(failures)} of {len(required)} TICs")
    temporary.replace(out_h5)
    failures_path.unlink(missing_ok=True)
    return {
        "n_requested": len(required),
        "n_written": written,
        "n_compact_adp_detectors": len(detector_lookup),
        "out_h5": str(out_h5),
    }


def _copy_attrs(source: Any, destination: Any) -> None:
    for key, value in source.attrs.items():
        destination.attrs[str(key)] = value


def _raw_source_payload(group: Any) -> dict[str, np.ndarray]:
    return {name: np.asarray(group[name]) for name in group.keys()}


def _resolve_path(value: Any, *, repo_root: Path) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else Path(repo_root) / path


def _bool_column(rows: pd.DataFrame, name: str) -> pd.Series:
    if name not in rows:
        return pd.Series(False, index=rows.index)
    values = rows[name]
    if values.dtype == bool:
        return values.fillna(False)
    return values.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})


def build_raw_pair_export(
    *,
    training_table: Path,
    raw_source_h5: Path,
    compact_adp_h5: Path,
    injection_pair_h5: Path | None,
    out_h5: Path,
    repo_root: Path,
    n_periods: int = DEFAULT_BLS_PERIODS,
    shard_index: int = 0,
    n_shards: int = 1,
) -> dict[str, Any]:
    """Assemble the final native contract from compact products on ORCD."""

    import h5py

    rows = pd.read_csv(training_table, low_memory=False)
    rows = rows.drop_duplicates("review_id", keep="last") if "review_id" in rows else rows
    active = (
        _bool_column(rows, "morphology_include_v1")
        | _bool_column(rows, "preserve_include_v1")
        | _bool_column(rows, "harmonic_include_v1")
    )
    rows = rows.loc[active].copy()
    rows["native_group_path"] = [native_group_path(row) for row in rows.to_dict("records")]
    if n_shards < 1 or shard_index < 0 or shard_index >= n_shards:
        raise ValueError("shard_index must satisfy 0 <= shard_index < n_shards")
    shard = rows["native_group_path"].map(
        lambda value: int.from_bytes(
            hashlib.sha1(str(value).encode("utf-8")).digest()[:8], "big"
        )
        % int(n_shards)
    )
    rows = rows.loc[shard.eq(int(shard_index))].copy()
    real_rows = rows[~rows["native_group_path"].str.startswith("injections/")].drop_duplicates("tic")
    injection_rows = rows[rows["native_group_path"].str.startswith("injections/")].drop_duplicates("injection_id")

    out_h5 = Path(out_h5)
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_h5.with_suffix(out_h5.suffix + ".tmp")
    failures: list[dict[str, Any]] = []
    with (
        h5py.File(raw_source_h5, "r") as raw_file,
        h5py.File(compact_adp_h5, "r") as adp_file,
        h5py.File(temporary, "w") as output,
    ):
        if str(raw_file.attrs.get("contract_version", "")) != RAW_SOURCE_CONTRACT_VERSION:
            raise ValueError("raw source HDF5 has the wrong contract_version")
        output.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
        output.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        output.attrs["training_table"] = str(training_table)
        output.attrs["raw_source_h5"] = str(raw_source_h5)
        output.attrs["compact_adp_h5"] = str(compact_adp_h5)
        output.attrs["injection_pair_h5"] = str(injection_pair_h5 or "")
        output.attrs["time_system"] = "BJD"
        output.attrs["periodogram_grid"] = "log10_period_d"
        output.attrs["periodogram_n"] = int(n_periods)
        output.attrs["chronology_small_channels"] = json.dumps(CHRONOLOGY_SMALL_CHANNELS)
        output.attrs["chronology_supplemental_channels"] = json.dumps(
            CHRONOLOGY_SUPPLEMENTAL_CHANNELS
        )
        output.attrs["harmonic_view_channels"] = json.dumps(HARMONIC_VIEW_CHANNELS)
        output.attrs["periodogram_channels"] = json.dumps(PERIODOGRAM_CHANNELS)
        output.attrs["shard_index"] = int(shard_index)
        output.attrs["n_shards"] = int(n_shards)
        target_root = output.create_group("targets")
        injection_root = output.create_group("injections")

        for count, (_, row) in enumerate(real_rows.iterrows(), start=1):
            tic = int(row["tic"])
            raw_path = f"targets/{tic:016d}"
            adp_path = f"targets/{tic:016d}"
            try:
                raw = _raw_source_payload(raw_file[raw_path])
                adp = adp_file[adp_path]
                aligned = align_raw_by_cadence(
                    raw,
                    cadenceno=np.asarray(adp["cadenceno"]),
                    time=np.asarray(adp["time"]),
                )
                group = target_root.create_group(f"{tic:016d}")
                _copy_attrs(adp, group)
                group.attrs["raw_source_paths"] = raw_file[raw_path].attrs.get(
                    "source_paths", ""
                )
                group.attrs["raw_adp_time_delta_max_s"] = float(
                    np.nanmax(aligned["_time_delta_s"])
                )
                payload = {
                    "time": _absolute_bjd(np.asarray(adp["time"])),
                    "cadenceno": np.asarray(adp["cadenceno"], dtype=np.int64),
                    "orbitid": np.asarray(adp["orbitid"], dtype=np.int32),
                    "quality": np.asarray(adp["quality"], dtype=np.int32),
                    "raw_flux_small": aligned["raw_flux_small"],
                    "raw_flux_err_small": aligned["raw_flux_err_small"],
                    "raw_flux_primary": aligned["raw_flux_primary"],
                    "raw_flux_err_primary": aligned["raw_flux_err_primary"],
                    "det_flux_adp_sml": np.asarray(adp["DET_FLUX_ADP_SML"]),
                    "det_flux_adp": np.asarray(adp["DET_FLUX_ADP"]),
                }
                payload.update(
                    _periodogram_payload(
                        time=payload["time"],
                        quality=payload["quality"],
                        small=payload["det_flux_adp_sml"],
                        primary=payload["det_flux_adp"],
                        n_periods=n_periods,
                    )
                )
                for name, values in payload.items():
                    _write_dataset(group, name, values)
            except Exception as exc:
                failures.append({"kind": "real", "id": tic, "error": str(exc)})
            if count % 100 == 0:
                print(f"[raw-pair-export] real {count}/{len(real_rows)}", flush=True)

        pair_files: dict[Path, Any] = {}
        canonical_files: dict[Path, Any] = {}
        try:
            for count, (_, row) in enumerate(injection_rows.iterrows(), start=1):
                injection_id = str(row["injection_id"])
                try:
                    pair_path = (
                        _resolve_path(injection_pair_h5, repo_root=repo_root)
                        if injection_pair_h5 is not None
                        else _resolve_path(row["source_h5"], repo_root=repo_root)
                    )
                    pair_file = pair_files.setdefault(pair_path, h5py.File(pair_path, "r"))
                    pair_group_path = (
                        f"injections/{injection_id}"
                        if injection_pair_h5 is not None
                        else str(row["h5_group"])
                    )
                    pair = pair_file[pair_group_path]
                    canonical_path = _resolve_path(
                        pair.attrs.get("source_injection_h5", pair_file.attrs["source_injection_h5"]),
                        repo_root=repo_root,
                    )
                    canonical_file = canonical_files.setdefault(
                        canonical_path, h5py.File(canonical_path, "r")
                    )
                    canonical = canonical_file[f"injections/{injection_id}"]
                    tic = int(pair.attrs["tic"])
                    raw = _raw_source_payload(raw_file[f"targets/{tic:016d}"])
                    aligned = align_raw_by_cadence(
                        raw,
                        cadenceno=np.asarray(pair["cadenceno"]),
                        time=np.asarray(pair["time"]),
                    )
                    model = np.asarray(pair["transit_model"], dtype=float)
                    cadence_s = float(pair.attrs.get("cadence_s", 200.0))
                    group = injection_root.create_group(injection_id)
                    _copy_attrs(pair, group)
                    group.attrs["raw_source_paths"] = raw_file[
                        f"targets/{tic:016d}"
                    ].attrs.get("source_paths", "")
                    group.attrs["raw_adp_time_delta_max_s"] = float(
                        np.nanmax(aligned["_time_delta_s"])
                    )
                    payload = {
                        "time": _absolute_bjd(np.asarray(pair["time"])),
                        "cadenceno": np.asarray(pair["cadenceno"], dtype=np.int64),
                        "orbitid": np.asarray(pair["orbitid"], dtype=np.int32),
                        "quality": np.asarray(pair["quality"], dtype=np.int32),
                        "raw_flux_small": np.asarray(canonical["RAW_FLUX_Small_injected"]),
                        "raw_flux_primary": np.asarray(canonical["RAW_FLUX_Primary_injected"]),
                        "det_flux_adp_sml": np.asarray(pair["DET_FLUX_ADP_SML_injected"]),
                        "det_flux_adp": np.asarray(pair["DET_FLUX_ADP_injected"]),
                        "paired_original_raw_flux_small": np.asarray(canonical["RAW_FLUX_Small_original"]),
                        "paired_original_raw_flux_primary": np.asarray(canonical["RAW_FLUX_Primary_original"]),
                        "paired_original_det_flux_adp_sml": np.asarray(pair["DET_FLUX_ADP_SML_original"]),
                        "paired_original_det_flux_adp": np.asarray(pair["DET_FLUX_ADP_original"]),
                    }
                    for aperture, suffix in (("Small", "small"), ("Primary", "primary")):
                        raw_error = np.asarray(aligned[f"raw_flux_err_{suffix}"], dtype=float)
                        baseline = float(
                            pair.attrs.get(
                                f"injection_baseline_{aperture}",
                                pair.attrs.get(f"raw_baseline_{aperture}", np.nan),
                            )
                        )
                        payload[f"raw_flux_err_{suffix}"] = injected_raw_uncertainty(
                            raw_error,
                            model,
                            source_flux_rate=baseline,
                            cadence_s=cadence_s,
                        )
                        payload[f"paired_original_raw_flux_err_{suffix}"] = raw_error
                    payload.update(
                        _periodogram_payload(
                            time=payload["time"],
                            quality=payload["quality"],
                            small=payload["det_flux_adp_sml"],
                            primary=payload["det_flux_adp"],
                            n_periods=n_periods,
                        )
                    )
                    for name, values in payload.items():
                        _write_dataset(group, name, values)
                except Exception as exc:
                    failures.append({"kind": "injection", "id": injection_id, "error": str(exc)})
                if count % 50 == 0:
                    print(f"[raw-pair-export] injection {count}/{len(injection_rows)}", flush=True)
        finally:
            for handle in (*pair_files.values(), *canonical_files.values()):
                handle.close()

    if failures:
        pd.DataFrame(failures).to_csv(out_h5.with_suffix(".failures.csv"), index=False)
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"final raw-pair export failed for {len(failures)} rows")
    temporary.replace(out_h5)
    return {
        "n_real_targets": len(real_rows),
        "n_injections": len(injection_rows),
        "out_h5": str(out_h5),
        "native_datasets": list(NATIVE_DATASETS),
        "periodogram_datasets": list(PERIODOGRAM_DATASETS),
        "shard_index": int(shard_index),
        "n_shards": int(n_shards),
    }


def merge_raw_pair_shards(*, shard_paths: Sequence[Path], out_h5: Path) -> dict[str, Any]:
    """Merge disjoint native-contract shards and reject duplicate groups."""

    import h5py

    paths = [Path(path) for path in shard_paths]
    if not paths:
        raise ValueError("no raw-pair shards supplied")
    out_h5 = Path(out_h5)
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_h5.with_suffix(out_h5.suffix + ".tmp")
    counts = {"targets": 0, "injections": 0}
    with h5py.File(temporary, "w") as output:
        output.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
        output.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        output.attrs["merged_shards"] = json.dumps([str(path) for path in paths])
        output.attrs["time_system"] = "BJD"
        target_root = output.create_group("targets")
        injection_root = output.create_group("injections")
        roots = {"targets": target_root, "injections": injection_root}
        for path in paths:
            with h5py.File(path, "r") as source:
                if str(source.attrs.get("contract_version", "")) != RAW_PAIR_CONTRACT_VERSION:
                    raise ValueError(f"wrong contract_version in shard {path}")
                for name in (
                    "training_table",
                    "raw_source_h5",
                    "compact_adp_h5",
                    "injection_pair_h5",
                    "periodogram_grid",
                    "periodogram_n",
                    "chronology_small_channels",
                    "chronology_supplemental_channels",
                    "harmonic_view_channels",
                    "periodogram_channels",
                ):
                    if name in source.attrs and name not in output.attrs:
                        output.attrs[name] = source.attrs[name]
                for root_name, destination in roots.items():
                    if root_name not in source:
                        continue
                    for key in source[root_name]:
                        if key in destination:
                            raise ValueError(f"duplicate /{root_name}/{key} while merging {path}")
                        source.copy(source[f"{root_name}/{key}"], destination, name=key)
                        counts[root_name] += 1
    temporary.replace(out_h5)
    return {"out_h5": str(out_h5), "n_shards": len(paths), "counts": counts}


__all__ = [
    "DEFAULT_BLS_DURATIONS_MIN",
    "DEFAULT_BLS_PERIOD_RANGE_D",
    "DEFAULT_BLS_PERIODS",
    "RAW_SOURCE_CONTRACT_VERSION",
    "align_raw_by_cadence",
    "compact_adp_detector_lookup",
    "build_raw_pair_export",
    "discover_tglc_paths",
    "export_tglc_raw_sources",
    "merge_tglc_raw_paths",
    "merge_raw_pair_shards",
    "shared_bls_periodogram",
    "shared_bls_spectrum",
]
