"""Native-cadence inputs for the S56 post-adjudication harmonic CNN."""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


RAW_PAIR_CONTRACT_VERSION = "s56_adp_raw_pair_v1"
A2V1_TEACHER_INPUT_CONTRACT = "s56_A2v1_adp_raw_pair_v1"
RAW_PAIR_APERTURES: tuple[str, str] = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
HARMONIC_FACTORS: tuple[float, ...] = (0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0)
HARMONIC_NAMES: tuple[str, ...] = (
    "p_over_4",
    "p_over_3",
    "p_over_2",
    "p",
    "2p",
    "3p",
    "4p",
)
CHRONOLOGY_SMALL_CHANNELS: tuple[str, ...] = (
    "raw_small_centered_scaled",
    "adp_small_relative",
    "raw_error_small_scaled",
    "relative_time",
    "log_time_gap",
    "relative_cadence_number",
    "log_cadence_gap",
    "relative_orbit_index",
    "orbit_boundary",
    "quality_nonzero",
)
CHRONOLOGY_SUPPLEMENTAL_CHANNELS: tuple[str, ...] = (
    "raw_primary_centered_scaled",
    "adp_primary_relative",
    "raw_error_primary_scaled",
    "raw_primary_minus_small",
    "adp_primary_minus_small",
)
HARMONIC_VIEW_CHANNELS: tuple[str, ...] = (
    "adp_small_relative",
    "adp_primary_relative",
    "adp_primary_minus_small",
    "raw_error_small_scaled",
    "raw_error_primary_scaled",
    "orbital_phase",
    "quality_nonzero",
)
PERIODOGRAM_CHANNELS: tuple[str, ...] = (
    "bls_power_small",
    "bls_sde_small",
    "bls_power_primary",
    "bls_sde_primary",
)
NATIVE_DATASETS: tuple[str, ...] = (
    "time",
    "cadenceno",
    "orbitid",
    "quality",
    "raw_flux_small",
    "raw_flux_err_small",
    "raw_flux_primary",
    "raw_flux_err_primary",
    "det_flux_adp_sml",
    "det_flux_adp",
)
PERIODOGRAM_DATASETS: tuple[str, ...] = (
    "bls_log_period_grid",
    "bls_power_small",
    "bls_sde_small",
    "bls_power_primary",
    "bls_sde_primary",
)
CHANNEL_CONTRACT: Mapping[str, tuple[str, ...]] = {
    "chronology_small_channels": CHRONOLOGY_SMALL_CHANNELS,
    "chronology_supplemental_channels": CHRONOLOGY_SUPPLEMENTAL_CHANNELS,
    "harmonic_view_channels": HARMONIC_VIEW_CHANNELS,
    "periodogram_channels": PERIODOGRAM_CHANNELS,
}


def native_group_path(row: Mapping[str, Any]) -> str:
    """Return the native HDF5 group for a real or injected training row."""

    raw_injected = row.get("is_injected_row", "")
    injected = (
        bool(raw_injected)
        if isinstance(raw_injected, (bool, np.bool_))
        else str(raw_injected).strip().lower() in {"1", "1.0", "true", "t", "yes", "y"}
    )
    if not injected:
        injected = "inject" in str(row.get("source_kind", "")).lower()
    if injected:
        injection_id = str(row.get("injection_id", "")).strip()
        if not injection_id:
            raise ValueError("injected training row has no injection_id")
        return f"injections/{injection_id}"
    try:
        tic = int(float(row["tic"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("real training row has no valid TIC") from exc
    return f"targets/{tic:016d}"


@dataclass(frozen=True)
class NativeLightCurve:
    """One real or injected two-aperture light curve at native cadence."""

    time: np.ndarray
    cadenceno: np.ndarray
    orbitid: np.ndarray
    quality: np.ndarray
    raw_flux_small: np.ndarray
    raw_flux_err_small: np.ndarray
    raw_flux_primary: np.ndarray
    raw_flux_err_primary: np.ndarray
    det_flux_adp_sml: np.ndarray
    det_flux_adp: np.ndarray
    attrs: Mapping[str, Any]
    paired_original: "NativeLightCurve | None" = None
    bls_log_period_grid: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float32))
    bls_power_small: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float32))
    bls_sde_small: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float32))
    bls_power_primary: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float32))
    bls_sde_primary: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float32))

    def validate(self, *, require_errors: bool = True) -> None:
        arrays = {
            "time": self.time,
            "cadenceno": self.cadenceno,
            "orbitid": self.orbitid,
            "quality": self.quality,
            "raw_flux_small": self.raw_flux_small,
            "raw_flux_err_small": self.raw_flux_err_small,
            "raw_flux_primary": self.raw_flux_primary,
            "raw_flux_err_primary": self.raw_flux_err_primary,
            "det_flux_adp_sml": self.det_flux_adp_sml,
            "det_flux_adp": self.det_flux_adp,
        }
        lengths: dict[str, int] = {}
        for name, values in arrays.items():
            arr = np.asarray(values)
            if arr.ndim != 1:
                raise ValueError(f"native dataset {name} must be one-dimensional; got {arr.shape}")
            lengths[name] = len(arr)
        if len(set(lengths.values())) != 1:
            raise ValueError(f"native datasets have inconsistent lengths: {lengths}")
        if lengths["time"] == 0:
            raise ValueError("native light curve is empty")
        finite_time = np.asarray(self.time, dtype=float)
        finite_time = finite_time[np.isfinite(finite_time)]
        if finite_time.size and float(np.nanmedian(finite_time)) < 1.0e5:
            raise ValueError("native time must use absolute BJD")
        if finite_time.size > 1 and np.any(np.diff(finite_time) < 0):
            raise ValueError("native time must be chronological")
        if require_errors:
            for name in ("raw_flux_err_small", "raw_flux_err_primary"):
                values = np.asarray(arrays[name], dtype=float)
                if not np.any(np.isfinite(values) & (values > 0)):
                    raise ValueError(f"native light curve has no finite positive {name}")
        periodogram_lengths = {
            "bls_log_period_grid": len(np.asarray(self.bls_log_period_grid)),
            "bls_power_small": len(np.asarray(self.bls_power_small)),
            "bls_sde_small": len(np.asarray(self.bls_sde_small)),
            "bls_power_primary": len(np.asarray(self.bls_power_primary)),
            "bls_sde_primary": len(np.asarray(self.bls_sde_primary)),
        }
        nonzero = {length for length in periodogram_lengths.values() if length > 0}
        if nonzero and (len(nonzero) != 1 or 0 in periodogram_lengths.values()):
            raise ValueError(f"periodogram datasets have inconsistent lengths: {periodogram_lengths}")


@dataclass(frozen=True)
class NativeChannels:
    small_values: np.ndarray
    small_mask: np.ndarray
    supplemental_values: np.ndarray
    supplemental_mask: np.ndarray


@dataclass(frozen=True)
class HarmonicViews:
    full_values: tuple[np.ndarray, ...]
    full_masks: tuple[np.ndarray, ...]
    primary_values: tuple[np.ndarray, ...]
    primary_masks: tuple[np.ndarray, ...]
    secondary_values: tuple[np.ndarray, ...]
    secondary_masks: tuple[np.ndarray, ...]
    factors: tuple[float, ...] = HARMONIC_FACTORS


def _safe_scale(values: np.ndarray, errors: np.ndarray | None = None) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size:
        median = float(np.nanmedian(finite))
        scale = 1.4826 * float(np.nanmedian(np.abs(finite - median)))
        if np.isfinite(scale) and scale > 0:
            return scale
    if errors is not None:
        err = np.asarray(errors, dtype=float)
        err = err[np.isfinite(err) & (err > 0)]
        if err.size:
            scale = float(np.nanmedian(err))
            if np.isfinite(scale) and scale > 0:
                return scale
    return 1.0


def _normalize_raw_by_orbit(
    flux: np.ndarray,
    error: np.ndarray,
    orbitid: np.ndarray,
    quality: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    flux = np.asarray(flux, dtype=float)
    error = np.asarray(error, dtype=float)
    orbitid = np.asarray(orbitid)
    quality = np.asarray(quality)
    normalized = np.full(len(flux), np.nan, dtype=np.float32)
    normalized_error = np.full(len(flux), np.nan, dtype=np.float32)
    for orbit in np.unique(orbitid):
        in_orbit = orbitid == orbit
        reference = in_orbit & (quality == 0) & np.isfinite(flux)
        if not np.any(reference):
            reference = in_orbit & np.isfinite(flux)
        center = float(np.nanmedian(flux[reference])) if np.any(reference) else 0.0
        scale = _safe_scale(flux[reference], error[in_orbit])
        finite_flux = in_orbit & np.isfinite(flux)
        finite_error = in_orbit & np.isfinite(error) & (error > 0)
        normalized[finite_flux] = ((flux[finite_flux] - center) / scale).astype(np.float32)
        normalized_error[finite_error] = (error[finite_error] / scale).astype(np.float32)
    return normalized, normalized_error


def _relative_det_flux(flux: np.ndarray, orbitid: np.ndarray, quality: np.ndarray) -> np.ndarray:
    flux = np.asarray(flux, dtype=float)
    orbitid = np.asarray(orbitid)
    quality = np.asarray(quality)
    out = np.full(len(flux), np.nan, dtype=np.float32)
    for orbit in np.unique(orbitid):
        in_orbit = orbitid == orbit
        reference = in_orbit & (quality == 0) & np.isfinite(flux)
        if not np.any(reference):
            reference = in_orbit & np.isfinite(flux)
        baseline = float(np.nanmedian(flux[reference])) if np.any(reference) else 1.0
        finite = in_orbit & np.isfinite(flux)
        if np.isfinite(baseline) and abs(baseline) > 1.0e-8:
            out[finite] = (flux[finite] / baseline - 1.0).astype(np.float32)
        else:
            out[finite] = (flux[finite] - baseline).astype(np.float32)
    return out


def _coordinate_channels(
    time: np.ndarray,
    cadenceno: np.ndarray,
    orbitid: np.ndarray,
    quality: np.ndarray,
) -> np.ndarray:
    time = np.asarray(time, dtype=float)
    cadenceno = np.asarray(cadenceno, dtype=float)
    orbitid = np.asarray(orbitid)
    quality = np.asarray(quality)
    relative_time = np.zeros(len(time), dtype=np.float32)
    finite_time = np.isfinite(time)
    if np.any(finite_time):
        low = float(np.nanmin(time[finite_time]))
        span = float(np.nanmax(time[finite_time]) - low)
        if span > 0:
            relative_time[finite_time] = ((time[finite_time] - low) / span).astype(np.float32)
    dt = np.zeros(len(time), dtype=np.float32)
    if len(time) > 1:
        delta = np.diff(time)
        positive = delta[np.isfinite(delta) & (delta > 0)]
        cadence = float(np.nanmedian(positive)) if positive.size else 1.0
        dt[1:] = np.log1p(np.maximum(delta, 0.0) / max(cadence, 1.0e-12)).astype(np.float32)
    cadence_position = np.zeros(len(cadenceno), dtype=np.float32)
    finite_cadence = np.isfinite(cadenceno)
    if np.any(finite_cadence):
        low = float(np.nanmin(cadenceno[finite_cadence]))
        span = float(np.nanmax(cadenceno[finite_cadence]) - low)
        if span > 0:
            cadence_position[finite_cadence] = (
                (cadenceno[finite_cadence] - low) / span
            ).astype(np.float32)
    cadence_gap = np.zeros(len(cadenceno), dtype=np.float32)
    if len(cadenceno) > 1:
        delta = np.diff(cadenceno)
        positive = delta[np.isfinite(delta) & (delta > 0)]
        typical = float(np.nanmedian(positive)) if positive.size else 1.0
        cadence_gap[1:] = np.log1p(
            np.maximum(delta / max(typical, 1.0e-12) - 1.0, 0.0)
        ).astype(np.float32)
    orbit_position = np.zeros(len(orbitid), dtype=np.float32)
    unique_orbits = np.unique(orbitid)
    if len(unique_orbits) > 1:
        orbit_lookup = {
            value: index / float(len(unique_orbits) - 1)
            for index, value in enumerate(unique_orbits)
        }
        orbit_position = np.asarray([orbit_lookup[value] for value in orbitid], dtype=np.float32)
    orbit_boundary = np.zeros(len(time), dtype=np.float32)
    if len(time) > 1:
        orbit_boundary[1:] = (orbitid[1:] != orbitid[:-1]).astype(np.float32)
    quality_nonzero = (quality != 0).astype(np.float32)
    return np.stack(
        [
            relative_time,
            dt,
            cadence_position,
            cadence_gap,
            orbit_position,
            orbit_boundary,
            quality_nonzero,
        ],
        axis=0,
    )


def build_native_channels(lc: NativeLightCurve) -> NativeChannels:
    """Return asymmetric small-aperture and supplemental native channels."""

    lc.validate(require_errors=False)
    raw_small, err_small = _normalize_raw_by_orbit(
        lc.raw_flux_small, lc.raw_flux_err_small, lc.orbitid, lc.quality
    )
    raw_primary, err_primary = _normalize_raw_by_orbit(
        lc.raw_flux_primary, lc.raw_flux_err_primary, lc.orbitid, lc.quality
    )
    det_small = _relative_det_flux(lc.det_flux_adp_sml, lc.orbitid, lc.quality)
    det_primary = _relative_det_flux(lc.det_flux_adp, lc.orbitid, lc.quality)
    shared = _coordinate_channels(lc.time, lc.cadenceno, lc.orbitid, lc.quality)
    small = np.concatenate(
        [np.stack([raw_small, det_small, err_small], axis=0), shared], axis=0
    ).astype(np.float32)
    supplemental = np.stack(
        [
            raw_primary,
            det_primary,
            err_primary,
            raw_primary - raw_small,
            det_primary - det_small,
        ],
        axis=0,
    ).astype(np.float32)
    return NativeChannels(
        small_values=small,
        small_mask=np.isfinite(small),
        supplemental_values=supplemental,
        supplemental_mask=np.isfinite(supplemental),
    )


def orbital_phase(time: np.ndarray, *, period_d: float, t0_bjd: float) -> np.ndarray:
    time = np.asarray(time, dtype=float)
    if not np.isfinite(period_d) or period_d <= 0 or not np.isfinite(t0_bjd):
        return np.full(len(time), np.nan, dtype=float)
    return ((time - float(t0_bjd) + 0.5 * period_d) % period_d) / period_d - 0.5


def _circular_distance(phase: np.ndarray, center: float) -> np.ndarray:
    return np.abs(((phase - center + 0.5) % 1.0) - 0.5)


def _median_cadence_days(time: np.ndarray) -> float:
    delta = np.diff(np.sort(np.asarray(time, dtype=float)))
    delta = delta[np.isfinite(delta) & (delta > 0)]
    return float(np.nanmedian(delta)) if delta.size else 200.0 / 86400.0


def build_harmonic_views(
    lc: NativeLightCurve,
    *,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    factors: Sequence[float] = HARMONIC_FACTORS,
) -> HarmonicViews:
    """Build complete and local unbinned folds for every symmetric factor."""

    lc.validate(require_errors=False)
    channels = build_native_channels(lc)
    det_small = channels.small_values[1]
    det_primary = channels.supplemental_values[1]
    err_small = channels.small_values[2]
    err_primary = channels.supplemental_values[2]
    quality_nonzero = (np.asarray(lc.quality) != 0).astype(np.float32)
    duration_d = float(duration_min) / 1440.0
    cadence_d = _median_cadence_days(lc.time)

    full_values: list[np.ndarray] = []
    full_masks: list[np.ndarray] = []
    primary_values: list[np.ndarray] = []
    primary_masks: list[np.ndarray] = []
    secondary_values: list[np.ndarray] = []
    secondary_masks: list[np.ndarray] = []

    for factor in factors:
        view_period = float(period_d) * float(factor)
        phase = orbital_phase(lc.time, period_d=view_period, t0_bjd=t0_bjd)
        order = np.argsort(np.where(np.isfinite(phase), phase, np.inf), kind="stable")
        values = np.stack(
            [
                det_small,
                det_primary,
                det_primary - det_small,
                err_small,
                err_primary,
                phase.astype(np.float32),
                quality_nonzero,
            ],
            axis=0,
        )[:, order].astype(np.float32)
        masks = np.isfinite(values)
        full_values.append(values)
        full_masks.append(masks)

        half_window_d = min(
            max(4.0 * duration_d, 3.0 * cadence_d),
            0.2 * view_period,
        )
        half_window_phase = half_window_d / max(view_period, 1.0e-12)
        primary_sel = _circular_distance(phase, 0.0) <= half_window_phase
        secondary_sel = _circular_distance(phase, 0.5) <= half_window_phase
        for selected, value_list, mask_list in (
            (primary_sel, primary_values, primary_masks),
            (secondary_sel, secondary_values, secondary_masks),
        ):
            local = values[:, selected[order]]
            # ``order`` is already phase sorted; boolean selection preserves it.
            value_list.append(local.astype(np.float32))
            mask_list.append(np.isfinite(local))

    return HarmonicViews(
        full_values=tuple(full_values),
        full_masks=tuple(full_masks),
        primary_values=tuple(primary_values),
        primary_masks=tuple(primary_masks),
        secondary_values=tuple(secondary_values),
        secondary_masks=tuple(secondary_masks),
        factors=tuple(float(value) for value in factors),
    )


def pad_channel_sequences(
    sequences: Sequence[np.ndarray],
    masks: Sequence[np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pad channel-first sequences without dropping any native samples."""

    if not sequences:
        return (
            np.empty((0, 0, 0), dtype=np.float32),
            np.empty((0, 0, 0), dtype=bool),
            np.empty(0, dtype=np.int32),
        )
    arrays = [np.asarray(value, dtype=np.float32) for value in sequences]
    n_channels = arrays[0].shape[0]
    if any(value.ndim != 2 or value.shape[0] != n_channels for value in arrays):
        raise ValueError("all sequences must have shape (channels, samples) with equal channels")
    lengths = np.asarray([value.shape[1] for value in arrays], dtype=np.int32)
    maximum = int(lengths.max(initial=0))
    values_out = np.zeros((len(arrays), n_channels, maximum), dtype=np.float32)
    mask_out = np.zeros((len(arrays), n_channels, maximum), dtype=bool)
    for index, value in enumerate(arrays):
        length = value.shape[1]
        values_out[index, :, :length] = np.nan_to_num(value, nan=0.0, posinf=0.0, neginf=0.0)
        source_mask = np.isfinite(value) if masks is None else np.asarray(masks[index], dtype=bool)
        if source_mask.shape != value.shape:
            raise ValueError("sequence masks must match their values")
        mask_out[index, :, :length] = source_mask
    return values_out, mask_out, lengths


def injected_raw_uncertainty(
    raw_error: np.ndarray,
    transit_model: np.ndarray,
    *,
    source_flux_rate: float,
    cadence_s: float,
) -> np.ndarray:
    """Propagate source-Poisson plus fixed-floor uncertainty through a transit.

    Raw flux and error are in rate units. The source Poisson variance in rate
    units is ``source_flux_rate / cadence_s``; all remaining variance is kept
    as a cadence-specific background/read floor. The observed noise realization
    itself is deliberately not resampled.
    """

    error = np.asarray(raw_error, dtype=float)
    model = np.asarray(transit_model, dtype=float)
    if error.shape != model.shape:
        raise ValueError("raw_error and transit_model must have identical shape")
    if not np.isfinite(cadence_s) or cadence_s <= 0:
        raise ValueError("cadence_s must be finite and positive")
    if not np.isfinite(source_flux_rate) or source_flux_rate < 0:
        raise ValueError("source_flux_rate must be finite and nonnegative")
    source_variance = float(source_flux_rate) / float(cadence_s)
    floor_variance = np.maximum(np.square(error) - source_variance, 0.0)
    injected_variance = floor_variance + np.clip(model, 0.0, np.inf) * source_variance
    out = np.sqrt(injected_variance)
    out[~np.isfinite(error) | ~np.isfinite(model)] = np.nan
    return out.astype(np.float32)


def _read_group(group: Any, *, paired_prefix: str = "") -> NativeLightCurve:
    def values(name: str, dtype: Any) -> np.ndarray:
        dataset = f"{paired_prefix}{name}"
        if dataset not in group:
            return np.full(len(group["time"]), np.nan, dtype=dtype)
        return np.asarray(group[dataset], dtype=dtype)

    paired = bool(paired_prefix)
    return NativeLightCurve(
        time=np.asarray(group["time"], dtype=np.float64),
        cadenceno=np.asarray(group["cadenceno"], dtype=np.int64),
        orbitid=np.asarray(group["orbitid"], dtype=np.int32),
        quality=np.asarray(group["quality"], dtype=np.int32),
        raw_flux_small=values("raw_flux_small", np.float64),
        raw_flux_err_small=values("raw_flux_err_small", np.float64),
        raw_flux_primary=values("raw_flux_primary", np.float64),
        raw_flux_err_primary=values("raw_flux_err_primary", np.float64),
        det_flux_adp_sml=values("det_flux_adp_sml", np.float64),
        det_flux_adp=values("det_flux_adp", np.float64),
        attrs={str(key): group.attrs[key] for key in group.attrs},
        bls_log_period_grid=(
            np.asarray(group["bls_log_period_grid"], dtype=np.float32)
            if not paired and "bls_log_period_grid" in group
            else np.empty(0, dtype=np.float32)
        ),
        bls_sde_small=(
            np.asarray(group["bls_sde_small"], dtype=np.float32)
            if not paired and "bls_sde_small" in group
            else np.empty(0, dtype=np.float32)
        ),
        bls_power_small=(
            np.asarray(group["bls_power_small"], dtype=np.float32)
            if not paired and "bls_power_small" in group
            else np.empty(0, dtype=np.float32)
        ),
        bls_power_primary=(
            np.asarray(group["bls_power_primary"], dtype=np.float32)
            if not paired and "bls_power_primary" in group
            else np.empty(0, dtype=np.float32)
        ),
        bls_sde_primary=(
            np.asarray(group["bls_sde_primary"], dtype=np.float32)
            if not paired and "bls_sde_primary" in group
            else np.empty(0, dtype=np.float32)
        ),
    )


def read_native_light_curve_from_h5(
    h5: Any,
    *,
    group_path: str,
    require_errors: bool = True,
) -> NativeLightCurve:
    """Read one native light curve from an already-open read-only HDF5 file."""

    contract = str(h5.attrs.get("contract_version", ""))
    if contract != RAW_PAIR_CONTRACT_VERSION:
        raise ValueError(f"unexpected native input contract {contract!r}")
    if str(h5.attrs.get("time_system", "")) != "BJD":
        raise ValueError("native input contract must declare time_system='BJD'")
    for name, expected in CHANNEL_CONTRACT.items():
        observed = tuple(json.loads(str(h5.attrs.get(name, "[]"))))
        if observed != expected:
            raise ValueError(f"native input {name} does not match the v1 channel contract")
    if group_path not in h5:
        raise KeyError(f"missing native input group: {group_path}")
    group = h5[group_path]
    lc = _read_group(group)
    paired_names = [f"paired_original_{name}" for name in NATIVE_DATASETS[4:]]
    paired = None
    if any(name in group for name in paired_names):
        paired = _read_group(group, paired_prefix="paired_original_")
        object.__setattr__(lc, "paired_original", paired)
    lc.validate(require_errors=require_errors)
    if lc.paired_original is not None:
        lc.paired_original.validate(require_errors=require_errors)
    return lc


def read_native_light_curve(
    path: Path,
    *,
    group_path: str,
    require_errors: bool = True,
) -> NativeLightCurve:
    import h5py

    with h5py.File(Path(path), "r") as h5:
        return read_native_light_curve_from_h5(
            h5,
            group_path=group_path,
            require_errors=require_errors,
        )


def verify_raw_pair_contract(
    path: Path,
    *,
    require_errors: bool = True,
    require_periodograms: bool = False,
) -> dict[str, Any]:
    """Verify every group in a native-input HDF5 without loading all at once."""

    import h5py

    path = Path(path)
    failures: list[str] = []
    counts: dict[str, int] = {"targets": 0, "injections": 0}
    if not path.exists():
        return {"passed": False, "failures": [f"missing file: {path}"], "counts": counts}
    with h5py.File(path, "r") as h5:
        contract = str(h5.attrs.get("contract_version", ""))
        if contract != RAW_PAIR_CONTRACT_VERSION:
            failures.append(f"contract_version={contract!r}")
        if str(h5.attrs.get("time_system", "")) != "BJD":
            failures.append("time_system must be BJD")
        for name, expected in CHANNEL_CONTRACT.items():
            try:
                observed = tuple(json.loads(str(h5.attrs.get(name, "[]"))))
            except (TypeError, ValueError, json.JSONDecodeError):
                observed = tuple()
            if observed != expected:
                failures.append(f"{name} does not match channel contract")
        for root_name in ("targets", "injections"):
            if root_name not in h5:
                continue
            for key, group in h5[root_name].items():
                counts[root_name] += 1
                missing = [name for name in NATIVE_DATASETS if name not in group]
                if missing:
                    failures.append(f"/{root_name}/{key}:missing={','.join(missing)}")
                    continue
                lengths = {name: len(group[name]) for name in NATIVE_DATASETS}
                if len(set(lengths.values())) != 1:
                    failures.append(f"/{root_name}/{key}:length_mismatch")
                time = np.asarray(group["time"], dtype=float)
                finite_time = time[np.isfinite(time)]
                if finite_time.size and float(np.nanmedian(finite_time)) < 1.0e5:
                    failures.append(f"/{root_name}/{key}:time_not_absolute_bjd")
                if finite_time.size > 1 and np.any(np.diff(finite_time) < 0):
                    failures.append(f"/{root_name}/{key}:time_not_chronological")
                if require_errors:
                    for name in ("raw_flux_err_small", "raw_flux_err_primary"):
                        error = np.asarray(group[name], dtype=float)
                        if not np.any(np.isfinite(error) & (error > 0)):
                            failures.append(f"/{root_name}/{key}:invalid_{name}")
                periodogram_missing = [name for name in PERIODOGRAM_DATASETS if name not in group]
                if require_periodograms and periodogram_missing:
                    failures.append(
                        f"/{root_name}/{key}:missing={','.join(periodogram_missing)}"
                    )
                if not periodogram_missing:
                    pg_lengths = {name: len(group[name]) for name in PERIODOGRAM_DATASETS}
                    if len(set(pg_lengths.values())) != 1:
                        failures.append(f"/{root_name}/{key}:periodogram_length_mismatch")
                    if require_periodograms:
                        for name in (
                            "bls_power_small",
                            "bls_sde_small",
                            "bls_power_primary",
                            "bls_sde_primary",
                        ):
                            values = np.asarray(group[name], dtype=float)
                            if not np.any(np.isfinite(values)):
                                failures.append(f"/{root_name}/{key}:invalid_{name}")
        if sum(counts.values()) == 0:
            failures.append("no target or injection groups")
    return {"passed": not failures, "failures": failures, "counts": counts}


__all__ = [
    "A2V1_TEACHER_INPUT_CONTRACT",
    "CHRONOLOGY_SMALL_CHANNELS",
    "CHRONOLOGY_SUPPLEMENTAL_CHANNELS",
    "CHANNEL_CONTRACT",
    "HARMONIC_FACTORS",
    "HARMONIC_NAMES",
    "HARMONIC_VIEW_CHANNELS",
    "HarmonicViews",
    "NATIVE_DATASETS",
    "PERIODOGRAM_DATASETS",
    "PERIODOGRAM_CHANNELS",
    "NativeChannels",
    "NativeLightCurve",
    "RAW_PAIR_APERTURES",
    "RAW_PAIR_CONTRACT_VERSION",
    "build_harmonic_views",
    "build_native_channels",
    "injected_raw_uncertainty",
    "native_group_path",
    "orbital_phase",
    "pad_channel_sequences",
    "read_native_light_curve",
    "read_native_light_curve_from_h5",
    "verify_raw_pair_contract",
]
