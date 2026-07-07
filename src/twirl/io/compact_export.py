"""Readers for compact TWIRL light-curve HDF5 exports."""
from __future__ import annotations

from pathlib import Path
from typing import Sequence
from dataclasses import fields

import numpy as np

from twirl.io.hlsp import HLSPLightCurve


def compact_target_key(tic: int) -> str:
    """Return the compact-export target group key for a TIC ID."""

    return f"{int(tic):016d}"


def read_compact_lc_export(
    export_h5: Path,
    *,
    tic: int,
    columns: Sequence[str],
) -> HLSPLightCurve | None:
    """Read one target from a compact S56 light-curve HDF5 export.

    The compact exports store targets under ``/targets/{tic:016d}`` and keep the
    same time, quality, and flux-column naming used by the HLSP FITS reader.
    Returning the shared ``HLSPLightCurve`` dataclass lets downstream vetting
    code run unchanged whether the input came from FITS or compact HDF5.
    """

    import h5py

    path = Path(export_h5)
    key = compact_target_key(tic)
    if not path.exists():
        return None
    with h5py.File(path, "r") as h5:
        if "targets" not in h5 or key not in h5["targets"]:
            return None
        group = h5["targets"][key]
        if "time" not in group or "quality" not in group:
            return None
        time = np.asarray(group["time"], dtype=np.float64)
        n = len(time)
        flux = {
            col: np.asarray(group[col], dtype=np.float64)
            for col in columns
            if col in group
        }
        if not flux:
            return None
        cadenceno = (
            np.asarray(group["cadenceno"], dtype=np.int32)
            if "cadenceno" in group
            else np.arange(n, dtype=np.int32)
        )
        orbitid = (
            np.asarray(group["orbitid"], dtype=np.int16)
            if "orbitid" in group
            else np.zeros(n, dtype=np.int16)
        )
        payload = {
            "tic": int(group.attrs.get("tic", tic)),
            "tmag": float(group.attrs.get("tessmag", np.nan)),
            "sector": int(group.attrs.get("sector", -1)),
            "cam": int(group.attrs.get("camera", -1)),
            "ccd": int(group.attrs.get("ccd", -1)),
            "ra": float(group.attrs.get("ra_obj", np.nan)),
            "dec": float(group.attrs.get("dec_obj", np.nan)),
            "time": time,
            "cadenceno": cadenceno,
            "orbitid": orbitid,
            "quality": np.asarray(group["quality"], dtype=np.int32),
            "flux": flux,
            "path": Path(f"{path}:targets/{key}"),
        }
    accepted = {field.name for field in fields(HLSPLightCurve)}
    return HLSPLightCurve(**{key: value for key, value in payload.items() if key in accepted})


def read_injected_lc_group(
    injection_h5: Path,
    *,
    group_path: str,
    columns: Sequence[str],
) -> HLSPLightCurve | None:
    """Read one pre-detrend injection group as an ``HLSPLightCurve``.

    Injection products store one group per injected signal under ``/injections``.
    Newer multi-aperture products use datasets named ``{aperture}_injected``;
    older single-aperture products use ``flux_injected`` plus an ``aperture``
    attribute. This reader supports both so vetting/report code can consume
    injection products through the same light-curve dataclass used for FITS and
    compact real-LC exports.
    """

    import h5py

    path = Path(injection_h5)
    if not path.exists():
        return None
    with h5py.File(path, "r") as h5:
        if group_path not in h5:
            return None
        group = h5[group_path]
        if "time" not in group or "quality" not in group:
            return None
        time = np.asarray(group["time"], dtype=np.float64)
        n = len(time)
        attrs = group.attrs
        flux: dict[str, np.ndarray] = {}
        for col in columns:
            injected_name = f"{col}_injected"
            if injected_name in group:
                flux[col] = np.asarray(group[injected_name], dtype=np.float64)
            elif "flux_injected" in group and str(attrs.get("aperture", "")) == col:
                flux[col] = np.asarray(group["flux_injected"], dtype=np.float64)
        if not flux:
            return None
        cadenceno = (
            np.asarray(group["cadenceno"], dtype=np.int32)
            if "cadenceno" in group
            else np.arange(n, dtype=np.int32)
        )
        orbitid = (
            np.asarray(group["orbitid"], dtype=np.int16)
            if "orbitid" in group
            else np.zeros(n, dtype=np.int16)
        )
        payload = {
            "tic": int(attrs.get("tic", -1)),
            "tmag": float(attrs.get("tessmag", np.nan)),
            "sector": int(attrs.get("sector", -1)),
            "cam": int(attrs.get("camera", -1)),
            "ccd": int(attrs.get("ccd", -1)),
            "ra": float(attrs.get("ra_obj", attrs.get("ra", np.nan))),
            "dec": float(attrs.get("dec_obj", attrs.get("dec", np.nan))),
            "time": time,
            "cadenceno": cadenceno,
            "orbitid": orbitid,
            "quality": np.asarray(group["quality"], dtype=np.int32),
            "flux": flux,
            "path": Path(f"{path}:{group_path}"),
        }
    accepted = {field.name for field in fields(HLSPLightCurve)}
    return HLSPLightCurve(**{key: value for key, value in payload.items() if key in accepted})
