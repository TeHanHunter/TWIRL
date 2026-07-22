"""Build the S56 WD 1856 independent-extraction evidence for Tier-1 QA.

The builder compares each member of the locked A2v1 ADP pair with an
explicitly mapped flux column from one cadence-level light curve made by an
external extraction.  It does not fetch, invent, or bless a reference product:
callers must supply the external table, its source product, and enough
provenance to audit the independence claim.

All photometric comparisons use the exact common ``CADENCENO`` intersection,
joint effective-quality-zero finite cadences, and the fixed WD 1856 ephemeris.
Effective quality is the native product mask OR the authoritative external
SPOC/QLP overlay.  The output CSV and JSON manifest implement the input contract consumed by
``evaluate_independent_extraction`` in :mod:`twirl.lightcurves.a2v1_tier1_qa`.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd

from twirl.lightcurves.a2v1_qa import (
    WD1856_PERIOD_D,
    WD1856_T0_BJD,
    WD1856_TIC,
    file_sha256,
)
from twirl.lightcurves.a2v1_tier1_qa import _dataframe_content_sha256
from twirl.lightcurves.external_quality import load_external_quality_reference
from twirl.vetting.adp_only import ADP_ONLY_APERTURES


SECTOR = 56
CURRENT_TIME_OFFSET_BJD = 2_457_000.0
CURRENT_EXTRACTOR_FAMILY = "MIT TGLC A2v1"
INDEPENDENT_CONTRACT_VERSION = "s56_a2v1_independent_extraction_v2"
INDEPENDENT_COMPARISON_MODE = "signal_timing_only"
EXPECTED_DETECTOR = (4, 1)


@dataclass(frozen=True)
class IndependentExtractionProvenance:
    """Required human-auditable provenance for the comparison products."""

    current_repository: str
    current_revision: str
    reference_extractor_family: str
    reference_repository: str
    reference_revision: str
    reference_pixel_source: str
    reference_product_aperture: str
    reference_target_id: str
    independence_basis: str

    def validate(self) -> None:
        values = {
            "current_repository": self.current_repository,
            "current_revision": self.current_revision,
            "reference_extractor_family": self.reference_extractor_family,
            "reference_repository": self.reference_repository,
            "reference_revision": self.reference_revision,
            "reference_pixel_source": self.reference_pixel_source,
            "reference_product_aperture": self.reference_product_aperture,
            "reference_target_id": self.reference_target_id,
            "independence_basis": self.independence_basis,
        }
        empty = [name for name, value in values.items() if not str(value).strip()]
        if empty:
            raise ValueError(f"independent-extraction provenance is empty: {empty}")
        family = self.reference_extractor_family.strip().lower()
        if any(token in family for token in ("tglc", "twirl")):
            raise ValueError("reference extractor must be outside the TGLC/TWIRL family")
        if self.reference_repository.strip() == self.current_repository.strip():
            raise ValueError("reference and current repositories must differ")
        target_digits = re.sub(r"\D", "", self.reference_target_id)
        if target_digits != str(WD1856_TIC):
            raise ValueError(
                f"reference_target_id must identify WD 1856 as TIC {WD1856_TIC}"
            )


def _normalize_reference_flux_columns(
    mapping: Mapping[str, str],
) -> dict[str, str]:
    """Validate an exact current-aperture to reference-column mapping."""

    if not isinstance(mapping, Mapping):
        raise TypeError("reference_flux_columns must be a mapping")
    normalized: dict[str, str] = {}
    for current, reference in mapping.items():
        if not isinstance(current, str) or not isinstance(reference, str):
            raise ValueError("reference flux mapping names must be strings")
        current_name = current.strip()
        reference_name = reference.strip()
        if not current_name or not reference_name:
            raise ValueError("reference flux mapping names must be nonempty")
        if current_name in normalized:
            raise ValueError(
                f"duplicate reference flux mapping for current aperture {current_name!r}"
            )
        normalized[current_name] = reference_name

    expected = set(ADP_ONLY_APERTURES)
    observed = set(normalized)
    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise ValueError(
            "reference flux mappings must exactly cover the active current "
            f"apertures; missing={missing}, extra={extra}"
        )
    if len(set(normalized.values())) != len(normalized):
        raise ValueError(
            "reference flux mappings must use a distinct external column for "
            "each active current aperture"
        )
    return {aperture: normalized[aperture] for aperture in ADP_ONLY_APERTURES}


def parse_reference_flux_column_mappings(
    specifications: Sequence[str],
) -> dict[str, str]:
    """Parse repeated ``CURRENT=REFERENCE`` command-line specifications."""

    if isinstance(specifications, (str, bytes)):
        raise ValueError(
            "reference flux mappings must be repeated CURRENT=REFERENCE values"
        )
    mapping: dict[str, str] = {}
    for specification in specifications:
        if not isinstance(specification, str) or specification.count("=") != 1:
            raise ValueError(
                "each reference flux mapping must have the form CURRENT=REFERENCE"
            )
        current, reference = (part.strip() for part in specification.split("=", 1))
        if not current or not reference:
            raise ValueError(
                "each reference flux mapping must have the form CURRENT=REFERENCE"
            )
        if current in mapping:
            raise ValueError(
                f"duplicate reference flux mapping for current aperture {current!r}"
            )
        mapping[current] = reference
    return _normalize_reference_flux_columns(mapping)


def _robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return np.nan
    median = float(np.median(values))
    return float(1.4826 * np.median(np.abs(values - median)))


def _recover_bls_ephemeris(
    *,
    time_bjd: np.ndarray,
    normalized_flux: np.ndarray,
    good: np.ndarray,
    period_relative_half_width: float,
    n_periods: int,
    durations_min: tuple[float, ...],
    oversample: int,
    min_depth_snr: float,
) -> dict[str, float]:
    """Recover a bounded WD 1856 ephemeris from one observed flux series."""

    from astropy.timeseries import BoxLeastSquares

    time = np.asarray(time_bjd, dtype=float)[good]
    flux = np.asarray(normalized_flux, dtype=float)[good]
    if len(time) < 100:
        raise ValueError("bounded BLS requires at least 100 matched good cadences")
    scatter = _robust_mad(flux)
    if not np.isfinite(scatter) or scatter <= 0:
        raise ValueError("bounded BLS cannot estimate a positive robust scatter")
    period_low = WD1856_PERIOD_D * (1.0 - period_relative_half_width)
    period_high = WD1856_PERIOD_D * (1.0 + period_relative_half_width)
    periods = np.linspace(period_low, period_high, n_periods, dtype=float)
    durations_d = np.asarray(durations_min, dtype=float) / 1440.0
    time_origin = float(np.min(time))
    model = BoxLeastSquares(
        time - time_origin,
        flux,
        dy=np.full(len(flux), scatter, dtype=float),
    )
    result = model.power(periods, durations_d, oversample=oversample)
    power = np.asarray(result.power, dtype=float)
    if not np.isfinite(power).any():
        raise ValueError("bounded BLS produced no finite power")
    index = int(np.nanargmax(power))
    recovered = {
        "period_d": float(np.asarray(result.period, dtype=float)[index]),
        "t0_bjd": float(np.asarray(result.transit_time, dtype=float)[index])
        + time_origin,
        "duration_min": float(np.asarray(result.duration, dtype=float)[index]) * 1440.0,
        "depth": float(np.asarray(result.depth, dtype=float)[index]),
        "depth_snr": float(np.asarray(result.depth_snr, dtype=float)[index]),
        "power": float(power[index]),
    }
    if not np.isfinite(list(recovered.values())).all():
        raise ValueError("bounded BLS produced nonfinite ephemeris metrics")
    if recovered["depth"] <= 0 or recovered["depth_snr"] < min_depth_snr:
        raise ValueError(
            "bounded BLS signal is not significant: "
            f"depth={recovered['depth']:.6g}, S/N={recovered['depth_snr']:.6g}"
        )
    return recovered


def _read_reference_table(
    path: Path,
) -> tuple[pd.DataFrame, dict[str, Any], str]:
    """Read a cadence table and retain target/sector FITS header assertions."""

    path = Path(path)
    suffix = path.suffix.lower()
    header: dict[str, Any] = {}
    if suffix == ".parquet":
        table = pd.read_parquet(path)
        source_format = "parquet"
    elif suffix == ".ecsv":
        from astropy.table import Table

        table = Table.read(path, format="ascii.ecsv").to_pandas()
        source_format = "ecsv"
    elif suffix in {".fits", ".fit", ".fts", ".fz"} or path.name.lower().endswith(
        (".fits.gz", ".fit.gz", ".fts.gz")
    ):
        from astropy.io import fits
        from astropy.table import Table

        with fits.open(path, memmap=False) as hdul:
            for hdu in hdul:
                if hdu.header.get("TICID") is not None:
                    header["tic"] = int(hdu.header["TICID"])
                    break
            for hdu in hdul:
                if hdu.header.get("SECTOR") is not None:
                    header["sector"] = int(hdu.header["SECTOR"])
                    break
        table = Table.read(path, format="fits").to_pandas()
        source_format = "fits"
    elif suffix == ".csv" or path.name.lower().endswith(".csv.gz"):
        table = pd.read_csv(
            path, low_memory=False, float_precision="round_trip"
        )
        source_format = "csv"
    else:
        raise ValueError(f"unsupported external cadence-table format: {path}")
    if table.empty:
        raise ValueError("external cadence table is empty")
    return table, header, source_format


def _column_name(table: pd.DataFrame, requested: str) -> str:
    """Resolve one explicitly requested column, allowing only case variation."""

    if requested in table.columns:
        return requested
    matches = [str(name) for name in table.columns if str(name).lower() == requested.lower()]
    if len(matches) != 1:
        raise KeyError(f"external table lacks unambiguous column {requested!r}")
    return matches[0]


def _integer_array(values: pd.Series, *, label: str) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(numeric).all() or not np.equal(numeric, np.floor(numeric)).all():
        raise ValueError(f"{label} must contain finite integers")
    return numeric.astype(np.int64)


def _validate_reference_identity(
    table: pd.DataFrame,
    header: dict[str, Any],
    *,
    source_format: str,
    tic_column: str,
    sector_column: str,
    source_label: str,
) -> dict[str, Any]:
    """Require target/sector identity in the actual external product metadata."""

    if source_format == "fits":
        missing = [name for name in ("tic", "sector") if name not in header]
        if missing:
            raise ValueError(
                f"{source_label} FITS headers lack identity metadata: {missing}"
            )
        observed_tics = [int(header["tic"])]
        observed_sectors = [int(header["sector"])]
        identity_source = {"tic": "FITS:TICID", "sector": "FITS:SECTOR"}
        for requested, expected, label in (
            (tic_column, WD1856_TIC, "TIC"),
            (sector_column, SECTOR, "sector"),
        ):
            matches = [
                str(name)
                for name in table.columns
                if str(name).lower() == requested.lower()
            ]
            if len(matches) > 1:
                raise ValueError(
                    f"{source_label} FITS has ambiguous {label} identity columns"
                )
            if matches:
                observed = sorted(
                    set(
                        _integer_array(
                            table[matches[0]],
                            label=f"{source_label} FITS {label} identity",
                        )
                    )
                )
                if observed != [expected]:
                    raise ValueError(
                        f"{source_label} FITS header and {label} column disagree"
                    )
    else:
        try:
            tic_name = _column_name(table, tic_column)
            sector_name = _column_name(table, sector_column)
        except KeyError as exc:
            raise ValueError(
                f"{source_label} non-FITS table lacks explicit target/sector metadata"
            ) from exc
        observed_tics = sorted(
            set(_integer_array(table[tic_name], label=f"{source_label} TIC identity"))
        )
        observed_sectors = sorted(
            set(
                _integer_array(
                    table[sector_name], label=f"{source_label} sector identity"
                )
            )
        )
        identity_source = {"tic": tic_name, "sector": sector_name}
    if observed_tics != [WD1856_TIC]:
        raise ValueError(
            f"{source_label} does not identify only WD 1856 TIC {WD1856_TIC}: "
            f"{observed_tics}"
        )
    if observed_sectors != [SECTOR]:
        raise ValueError(
            f"{source_label} does not identify only Sector 56: {observed_sectors}"
        )
    return {
        "format": source_format,
        "tic": WD1856_TIC,
        "sector": SECTOR,
        "identity_source": identity_source,
    }


def _reference_arrays(
    table: pd.DataFrame,
    *,
    cadence_column: str,
    time_column: str,
    quality_column: str,
    flux_column: str,
    time_system: str,
) -> dict[str, np.ndarray | str]:
    cadence_name = _column_name(table, cadence_column)
    time_name = _column_name(table, time_column)
    quality_name = _column_name(table, quality_column)
    flux_name = _column_name(table, flux_column)

    cadence = _integer_array(table[cadence_name], label="reference CADENCENO")
    if len(np.unique(cadence)) != len(cadence):
        raise ValueError("reference cadence table contains duplicate CADENCENO values")
    time = pd.to_numeric(table[time_name], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(time).all():
        raise ValueError("reference TIME must be finite")
    normalized_time_system = str(time_system).strip().upper()
    if normalized_time_system == "BTJD":
        time_bjd = time + CURRENT_TIME_OFFSET_BJD
    elif normalized_time_system == "BJD":
        time_bjd = time
    else:
        raise ValueError("reference_time_system must be exactly BJD or BTJD")
    if not np.all((time_bjd > 2_450_000.0) & (time_bjd < 2_500_000.0)):
        raise ValueError("reference timestamps are not plausible BJD values")
    quality = _integer_array(table[quality_name], label="reference QUALITY")
    flux = pd.to_numeric(table[flux_name], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(flux).any():
        raise ValueError("reference flux column contains no finite values")
    return {
        "cadence": cadence,
        "time_bjd": time_bjd,
        "quality": quality,
        "flux": flux,
        "cadence_column": cadence_name,
        "time_column": time_name,
        "quality_column": quality_name,
        "flux_column": flux_name,
        "time_system": normalized_time_system,
    }


def _current_arrays(compact_lc: Path) -> dict[str, Any]:
    compact_lc = Path(compact_lc)
    key = f"{WD1856_TIC:016d}"
    with h5py.File(compact_lc, "r") as h5:
        if int(h5.attrs.get("sector", -1)) != SECTOR:
            raise ValueError("current compact product is not Sector 56")
        time_unit = str(h5.attrs.get("time_unit", "")).replace(" ", "").lower()
        if time_unit != "bjd-2457000":
            raise ValueError(
                "current compact product must declare time_unit='BJD - 2457000'"
            )
        if "targets" not in h5 or key not in h5["targets"]:
            raise KeyError(f"current compact product lacks WD 1856 TIC {WD1856_TIC}")
        group = h5["targets"][key]
        if int(group.attrs.get("tic", -1)) != WD1856_TIC:
            raise ValueError("WD 1856 compact group has the wrong TIC attribute")
        if int(group.attrs.get("sector", -1)) != SECTOR:
            raise ValueError("WD 1856 compact group has the wrong sector attribute")
        camera = int(group.attrs.get("camera", -1))
        ccd = int(group.attrs.get("ccd", -1))
        if (camera, ccd) != EXPECTED_DETECTOR:
            raise ValueError(
                "WD 1856 S56 must be on the locked cam4/ccd1 benchmark detector"
            )
        required = (
            "time",
            "cadenceno",
            "quality",
            "orbitid",
            *ADP_ONLY_APERTURES,
        )
        missing = [name for name in required if name not in group]
        if missing:
            raise KeyError(f"WD 1856 compact group lacks datasets: {missing}")
        time = np.asarray(group["time"], dtype=float)
        cadence_raw = np.asarray(group["cadenceno"])
        quality_raw = np.asarray(group["quality"])
        orbit_raw = np.asarray(group["orbitid"])
        cadence_float = np.asarray(cadence_raw, dtype=float)
        quality_float = np.asarray(quality_raw, dtype=float)
        orbit_float = np.asarray(orbit_raw, dtype=float)
        if (
            not np.isfinite(cadence_float).all()
            or not np.equal(cadence_float, np.floor(cadence_float)).all()
        ):
            raise ValueError("current WD 1856 CADENCENO must contain finite integers")
        if (
            not np.isfinite(quality_float).all()
            or not np.equal(quality_float, np.floor(quality_float)).all()
        ):
            raise ValueError("current WD 1856 QUALITY must contain finite integers")
        if (
            not np.isfinite(orbit_float).all()
            or not np.equal(orbit_float, np.floor(orbit_float)).all()
        ):
            raise ValueError("current WD 1856 ORBITID must contain finite integers")
        cadence = cadence_float.astype(np.int64)
        quality = quality_float.astype(np.int64)
        orbitid = orbit_float.astype(np.int64)
        flux = {
            aperture: np.asarray(group[aperture], dtype=float)
            for aperture in ADP_ONLY_APERTURES
        }
    lengths = {
        len(time),
        len(cadence),
        len(quality),
        len(orbitid),
        *(len(values) for values in flux.values()),
    }
    if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
        raise ValueError("current WD 1856 compact datasets have inconsistent lengths")
    if not np.isfinite(time).all():
        raise ValueError("current WD 1856 TIME contains nonfinite values")
    if len(np.unique(cadence)) != len(cadence):
        raise ValueError("current WD 1856 compact product has duplicate CADENCENO values")
    time_bjd = time + CURRENT_TIME_OFFSET_BJD
    if not np.all((time_bjd > 2_450_000.0) & (time_bjd < 2_500_000.0)):
        raise ValueError("current WD 1856 timestamps are not plausible BJD values")
    if any(not np.isfinite(values).any() for values in flux.values()):
        raise ValueError("one or more locked current apertures contain no finite flux")
    return {
        "cadence": cadence,
        "time_bjd": time_bjd,
        "quality": quality,
        "native_quality": quality.copy(),
        "orbitid": orbitid,
        "flux": flux,
        "camera": camera,
        "ccd": ccd,
    }


def _measure_one_aperture(
    *,
    aperture: str,
    current: dict[str, Any],
    reference: dict[str, Any],
    current_index: np.ndarray,
    reference_index: np.ndarray,
    transit_duration_min: float,
    oot_exclusion_min: float,
    max_time_delta_seconds: float,
    min_in_event_cadences: int,
    min_out_of_event_cadences: int,
    bls_period_relative_half_width: float,
    bls_n_periods: int,
    bls_durations_min: tuple[float, ...],
    bls_oversample: int,
    min_bls_depth_snr: float,
) -> dict[str, Any]:
    current_time = np.asarray(current["time_bjd"])[current_index]
    reference_time = np.asarray(reference["time_bjd"])[reference_index]
    time_delta_s = np.abs(current_time - reference_time) * 86_400.0
    max_delta_s = float(np.max(time_delta_s))
    if max_delta_s > max_time_delta_seconds:
        raise ValueError(
            "exactly matched CADENCENO timestamps differ by "
            f"{max_delta_s:.6g} s, above {max_time_delta_seconds:.6g} s"
        )

    current_flux = np.asarray(current["flux"][aperture], dtype=float)[current_index]
    reference_flux = np.asarray(reference["flux"], dtype=float)[reference_index]
    current_quality = np.asarray(current["quality"])[current_index]
    reference_quality = np.asarray(reference["quality"])[reference_index]
    joint_good = (
        (current_quality == 0)
        & (reference_quality == 0)
        & np.isfinite(current_flux)
        & np.isfinite(reference_flux)
    )
    phase_min = (
        ((current_time - WD1856_T0_BJD + 0.5 * WD1856_PERIOD_D) % WD1856_PERIOD_D)
        - 0.5 * WD1856_PERIOD_D
    ) * 1440.0
    in_event = joint_good & (np.abs(phase_min) <= 0.5 * transit_duration_min)
    out_of_event = joint_good & (np.abs(phase_min) >= oot_exclusion_min)
    n_in_event = int(np.count_nonzero(in_event))
    n_out_of_event = int(np.count_nonzero(out_of_event))
    if n_in_event < min_in_event_cadences:
        raise ValueError(
            f"{aperture} has only {n_in_event} matched in-event cadences; "
            f"need {min_in_event_cadences}"
        )
    if n_out_of_event < min_out_of_event_cadences:
        raise ValueError(
            f"{aperture} has only {n_out_of_event} matched out-of-event cadences; "
            f"need {min_out_of_event_cadences}"
        )

    current_baseline = float(np.median(current_flux[out_of_event]))
    reference_baseline = float(np.median(reference_flux[out_of_event]))
    if not np.isfinite(current_baseline) or current_baseline <= 0:
        raise ValueError(f"{aperture} current out-of-event baseline is not positive")
    if not np.isfinite(reference_baseline) or reference_baseline <= 0:
        raise ValueError("reference out-of-event baseline is not positive")
    current_normalized = current_flux / current_baseline
    reference_normalized = reference_flux / reference_baseline
    current_scatter = _robust_mad(current_normalized[out_of_event]) * 1.0e6
    reference_scatter = _robust_mad(reference_normalized[out_of_event]) * 1.0e6
    current_depth = 1.0 - float(np.median(current_normalized[in_event]))
    reference_depth = 1.0 - float(np.median(reference_normalized[in_event]))
    measured = (current_scatter, reference_scatter, current_depth, reference_depth)
    if not np.isfinite(measured).all() or current_scatter <= 0 or reference_scatter <= 0:
        raise ValueError(f"{aperture} produced invalid scatter or depth measurements")
    current_bls = _recover_bls_ephemeris(
        time_bjd=current_time,
        normalized_flux=current_normalized,
        good=joint_good,
        period_relative_half_width=bls_period_relative_half_width,
        n_periods=bls_n_periods,
        durations_min=bls_durations_min,
        oversample=bls_oversample,
        min_depth_snr=min_bls_depth_snr,
    )
    reference_bls = _recover_bls_ephemeris(
        time_bjd=reference_time,
        normalized_flux=reference_normalized,
        good=joint_good,
        period_relative_half_width=bls_period_relative_half_width,
        n_periods=bls_n_periods,
        durations_min=bls_durations_min,
        oversample=bls_oversample,
        min_depth_snr=min_bls_depth_snr,
    )

    return {
        "tic": WD1856_TIC,
        "sector": SECTOR,
        "aperture": aperture,
        "reference_flux_column": str(reference["flux_column"]),
        "detector": f"cam{current['camera']}_ccd{current['ccd']}",
        "n_current_cadences": int(len(current["cadence"])),
        "n_reference_cadences": int(len(reference["cadence"])),
        "n_common_cadences": int(len(current_index)),
        "current_scatter_ppm": float(current_scatter),
        "reference_scatter_ppm": float(reference_scatter),
        "current_period_d": current_bls["period_d"],
        "reference_period_d": reference_bls["period_d"],
        "current_t0_bjd": current_bls["t0_bjd"],
        "reference_t0_bjd": reference_bls["t0_bjd"],
        "current_depth": float(current_depth),
        "reference_depth": float(reference_depth),
        "current_bls_duration_min": current_bls["duration_min"],
        "reference_bls_duration_min": reference_bls["duration_min"],
        "current_bls_depth": current_bls["depth"],
        "reference_bls_depth": reference_bls["depth"],
        "current_bls_depth_snr": current_bls["depth_snr"],
        "reference_bls_depth_snr": reference_bls["depth_snr"],
        "current_bls_power": current_bls["power"],
        "reference_bls_power": reference_bls["power"],
        "n_common_quality0_finite": int(np.count_nonzero(joint_good)),
        "n_in_event_cadences": n_in_event,
        "n_out_of_event_cadences": n_out_of_event,
        "max_abs_time_delta_s": max_delta_s,
        "transit_duration_min": float(transit_duration_min),
        "oot_exclusion_min": float(oot_exclusion_min),
    }


def _output_preflight(
    *,
    inputs: tuple[Path, ...],
    metrics_csv: Path,
    manifest_json: Path,
    overwrite: bool,
) -> tuple[Path, Path]:
    metrics_csv = Path(metrics_csv)
    manifest_json = Path(manifest_json)
    if metrics_csv.suffix.lower() != ".csv":
        raise ValueError("independent metrics output must be a .csv file")
    if manifest_json.suffix.lower() != ".json":
        raise ValueError("independent manifest output must be a .json file")
    outputs = {metrics_csv.resolve(), manifest_json.resolve()}
    if len(outputs) != 2:
        raise ValueError("metrics CSV and manifest JSON must be different paths")
    input_paths = {Path(path).resolve() for path in inputs}
    collisions = sorted(str(path) for path in input_paths & outputs)
    if collisions:
        raise ValueError(f"independent inputs and outputs collide: {collisions}")
    if not overwrite:
        existing = [str(path) for path in (metrics_csv, manifest_json) if path.exists()]
        if existing:
            raise FileExistsError(f"outputs already exist; pass overwrite=True: {existing}")
    metrics_csv.parent.mkdir(parents=True, exist_ok=True)
    manifest_json.parent.mkdir(parents=True, exist_ok=True)
    return metrics_csv, manifest_json


def build_wd1856_independent_metrics(
    *,
    compact_lc: Path,
    reference_table: Path,
    reference_product: Path,
    cadence_reference_table: Path,
    cadence_reference_manifest: Path,
    metrics_csv: Path,
    manifest_json: Path,
    provenance: IndependentExtractionProvenance,
    reference_time_system: str,
    reference_flux_columns: Mapping[str, str],
    reference_cadence_column: str = "CADENCENO",
    reference_time_column: str = "TIME",
    reference_quality_column: str = "QUALITY",
    reference_tic_column: str = "TICID",
    reference_sector_column: str = "SECTOR",
    transit_duration_min: float = 8.0,
    oot_exclusion_min: float = 30.0,
    max_time_delta_seconds: float = 5.0,
    min_common_cadences: int = 100,
    min_in_event_cadences: int = 4,
    min_out_of_event_cadences: int = 100,
    bls_period_relative_half_width: float = 0.01,
    bls_n_periods: int = 4_001,
    bls_durations_min: tuple[float, ...] = (6.0, 8.0, 10.0),
    bls_oversample: int = 20,
    min_bls_depth_snr: float = 5.0,
    overwrite: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Write hash-bound WD 1856 metrics and the Tier-1 provenance manifest."""

    flux_column_mapping = _normalize_reference_flux_columns(reference_flux_columns)
    compact_lc = Path(compact_lc)
    reference_table = Path(reference_table)
    reference_product = Path(reference_product)
    cadence_reference_table = Path(cadence_reference_table)
    cadence_reference_manifest = Path(cadence_reference_manifest)
    for path in (
        compact_lc,
        reference_table,
        reference_product,
        cadence_reference_table,
        cadence_reference_manifest,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    input_paths = tuple(
        dict.fromkeys(
            path.resolve()
            for path in (
                compact_lc,
                reference_table,
                reference_product,
                cadence_reference_table,
                cadence_reference_manifest,
            )
        )
    )
    initial_input_sha256 = {
        str(path): file_sha256(path) for path in input_paths
    }
    initial_input_size = {str(path): int(path.stat().st_size) for path in input_paths}
    initial_input_mtime_ns = {
        str(path): int(path.stat().st_mtime_ns) for path in input_paths
    }
    provenance.validate()
    floating_parameters = np.asarray(
        [
            transit_duration_min,
            oot_exclusion_min,
            max_time_delta_seconds,
            bls_period_relative_half_width,
            min_bls_depth_snr,
            *bls_durations_min,
        ],
        dtype=float,
    )
    if not np.isfinite(floating_parameters).all() or np.any(floating_parameters <= 0):
        raise ValueError("duration, OOT window, and time tolerance must be finite and positive")
    if not 0 < transit_duration_min < oot_exclusion_min * 2:
        raise ValueError("transit_duration_min must fit inside the OOT exclusion window")
    if not 0.005 <= bls_period_relative_half_width <= 0.05:
        raise ValueError(
            "bounded BLS period half-width must be between 0.005 and 0.05"
        )
    if not bls_durations_min:
        raise ValueError("bounded BLS duration grid must be nonempty")
    if min_bls_depth_snr < 5.0:
        raise ValueError("bounded BLS minimum depth S/N cannot be below 5")
    integer_parameters = (
        min_common_cadences,
        min_in_event_cadences,
        min_out_of_event_cadences,
        bls_n_periods,
        bls_oversample,
    )
    if any(
        isinstance(value, bool)
        or not isinstance(value, (int, np.integer))
        or int(value) <= 0
        for value in integer_parameters
    ):
        raise ValueError("cadence minima and BLS grid sizes must be positive integers")
    if not 1_001 <= bls_n_periods <= 100_001:
        raise ValueError("bounded BLS n_periods must be between 1,001 and 100,001")
    if bls_oversample < 5:
        raise ValueError("bounded BLS oversample must be at least 5")
    metrics_csv, manifest_json = _output_preflight(
        inputs=(
            compact_lc,
            reference_table,
            reference_product,
            cadence_reference_table,
            cadence_reference_manifest,
        ),
        metrics_csv=metrics_csv,
        manifest_json=manifest_json,
        overwrite=overwrite,
    )

    compact_sha256 = initial_input_sha256[str(compact_lc.resolve())]
    reference_product_sha256 = initial_input_sha256[
        str(reference_product.resolve())
    ]
    reference_table_sha256 = initial_input_sha256[str(reference_table.resolve())]
    if compact_sha256 == reference_product_sha256:
        raise ValueError("current compact and external reference products have the same hash")
    current = _current_arrays(compact_lc)
    quality_reference = load_external_quality_reference(
        table_path=cadence_reference_table,
        manifest_path=cadence_reference_manifest,
        sector=SECTOR,
        expected_orbits=(119, 120),
    )
    current_full_overlay = quality_reference.apply(
        sector=SECTOR,
        camera=current["camera"],
        ccd=current["ccd"],
        cadenceno=current["cadence"],
        orbitid=current["orbitid"],
        internal_quality=current["native_quality"],
        context="WD 1856 current compact",
    )
    current["quality"] = current_full_overlay.quality
    current["external_quality"] = current_full_overlay.external_quality
    table, header, table_format = _read_reference_table(reference_table)
    table_identity = _validate_reference_identity(
        table,
        header,
        source_format=table_format,
        tic_column=reference_tic_column,
        sector_column=reference_sector_column,
        source_label="external cadence table",
    )
    if reference_product.resolve() == reference_table.resolve():
        product_identity = dict(table_identity)
    else:
        product_table, product_header, product_format = _read_reference_table(
            reference_product
        )
        product_identity = _validate_reference_identity(
            product_table,
            product_header,
            source_format=product_format,
            tic_column=reference_tic_column,
            sector_column=reference_sector_column,
            source_label="external source product",
        )
    reference_contexts: dict[str, dict[str, Any]] = {}
    for aperture in ADP_ONLY_APERTURES:
        aperture_reference = _reference_arrays(
            table,
            cadence_column=reference_cadence_column,
            time_column=reference_time_column,
            quality_column=reference_quality_column,
            flux_column=flux_column_mapping[aperture],
            time_system=reference_time_system,
        )
        common, current_index, reference_index = np.intersect1d(
            np.asarray(current["cadence"]),
            np.asarray(aperture_reference["cadence"]),
            assume_unique=True,
            return_indices=True,
        )
        if len(common) < min_common_cadences:
            raise ValueError(
                f"{aperture} has only {len(common)} exact common CADENCENO values; "
                f"need {min_common_cadences}"
            )
        current_common_overlay = quality_reference.apply(
            sector=SECTOR,
            camera=current["camera"],
            ccd=current["ccd"],
            cadenceno=current["cadence"][current_index],
            orbitid=current["orbitid"][current_index],
            internal_quality=current["native_quality"][current_index],
            context=f"WD 1856 current common {aperture}",
        )
        reference_common_overlay = quality_reference.apply(
            sector=SECTOR,
            camera=current["camera"],
            ccd=current["ccd"],
            cadenceno=aperture_reference["cadence"][reference_index],
            orbitid=current["orbitid"][current_index],
            internal_quality=aperture_reference["quality"][reference_index],
            context=f"WD 1856 reference common {aperture}",
        )
        effective_reference_quality = np.asarray(
            aperture_reference["quality"], dtype=np.int64
        ).copy()
        effective_reference_quality[reference_index] = (
            reference_common_overlay.quality
        )
        aperture_reference["quality"] = effective_reference_quality
        reference_contexts[aperture] = {
            "reference": aperture_reference,
            "common": common,
            "current_index": current_index,
            "reference_index": reference_index,
            "current_overlay_counts": dict(current_common_overlay.counts),
            "reference_overlay_counts": dict(reference_common_overlay.counts),
        }

    rows = []
    for aperture in ADP_ONLY_APERTURES:
        context = reference_contexts[aperture]
        rows.append(
            _measure_one_aperture(
                aperture=aperture,
                current=current,
                reference=context["reference"],
                current_index=context["current_index"],
                reference_index=context["reference_index"],
                transit_duration_min=transit_duration_min,
                oot_exclusion_min=oot_exclusion_min,
                max_time_delta_seconds=max_time_delta_seconds,
                min_in_event_cadences=min_in_event_cadences,
                min_out_of_event_cadences=min_out_of_event_cadences,
                bls_period_relative_half_width=bls_period_relative_half_width,
                bls_n_periods=bls_n_periods,
                bls_durations_min=bls_durations_min,
                bls_oversample=bls_oversample,
                min_bls_depth_snr=min_bls_depth_snr,
            )
        )
    metrics = pd.DataFrame(rows)
    first_context = reference_contexts[ADP_ONLY_APERTURES[0]]
    reference = first_context["reference"]
    common = first_context["common"]
    resolved_flux_column_mapping = {
        aperture: str(reference_contexts[aperture]["reference"]["flux_column"])
        for aperture in ADP_ONLY_APERTURES
    }
    current_common_counts = {
        tuple(sorted(context["current_overlay_counts"].items()))
        for context in reference_contexts.values()
    }
    reference_common_counts = {
        tuple(sorted(context["reference_overlay_counts"].items()))
        for context in reference_contexts.values()
    }
    if len(current_common_counts) != 1 or len(reference_common_counts) != 1:
        raise RuntimeError("aperture quality-overlay audit counts disagree")
    current_common_audit = dict(next(iter(current_common_counts)))
    reference_common_audit = dict(next(iter(reference_common_counts)))

    metrics_tmp = metrics_csv.with_suffix(metrics_csv.suffix + ".tmp")
    manifest_tmp = manifest_json.with_suffix(manifest_json.suffix + ".tmp")
    if metrics_tmp.resolve() in {
        compact_lc.resolve(),
        reference_table.resolve(),
        reference_product.resolve(),
        cadence_reference_table.resolve(),
        cadence_reference_manifest.resolve(),
    } or manifest_tmp.resolve() in {
        compact_lc.resolve(),
        reference_table.resolve(),
        reference_product.resolve(),
        cadence_reference_table.resolve(),
        cadence_reference_manifest.resolve(),
    }:
        raise ValueError("temporary output path collides with an input")
    try:
        metrics.to_csv(metrics_tmp, index=False)
        roundtrip = pd.read_csv(
            metrics_tmp, low_memory=False, float_precision="round_trip"
        )
        metrics_content_sha256 = _dataframe_content_sha256(roundtrip)
        metrics_file_sha256 = file_sha256(metrics_tmp)
        manifest: dict[str, Any] = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "contract_version": INDEPENDENT_CONTRACT_VERSION,
            "sector": SECTOR,
            "tic": WD1856_TIC,
            "independent": True,
            "comparison_mode": INDEPENDENT_COMPARISON_MODE,
            "independence_basis": provenance.independence_basis.strip(),
            "current_extractor_family": CURRENT_EXTRACTOR_FAMILY,
            "reference_extractor_family": provenance.reference_extractor_family.strip(),
            "current_repository": provenance.current_repository.strip(),
            "reference_repository": provenance.reference_repository.strip(),
            "current_revision": provenance.current_revision.strip(),
            "reference_revision": provenance.reference_revision.strip(),
            "pixel_source": provenance.reference_pixel_source.strip(),
            "reference_target_id": provenance.reference_target_id.strip(),
            "cadence_match_policy": (
                "exact common CADENCENO intersection; joint current/reference "
                "effective QUALITY==0 and finite flux for photometric metrics; "
                "effective quality is native internal OR authoritative external"
            ),
            "scatter_definition": (
                "1.4826*MAD in ppm after each series is normalized by its "
                f"matched out-of-event median; |phase| >= {oot_exclusion_min:g} min"
            ),
            "depth_definition": (
                "1 - median(matched normalized in-event flux) at fixed WD 1856 "
                f"P={WD1856_PERIOD_D:.9f} d, T0={WD1856_T0_BJD:.6f} BJD, "
                f"duration={transit_duration_min:g} min"
            ),
            "ephemeris_recovery_definition": (
                "independent bounded BoxLeastSquares recovery on each matched "
                "quality-zero flux series; published ephemeris only defines the "
                "predeclared search interval, fixed depth window, and later truth "
                "comparison"
            ),
            "current_apertures": list(ADP_ONLY_APERTURES),
            "reference_product_aperture": provenance.reference_product_aperture.strip(),
            "reference_flux_columns": resolved_flux_column_mapping,
            "current_compact_sha256": compact_sha256,
            "current_product_sha256": compact_sha256,
            "reference_product_sha256": reference_product_sha256,
            "reference_table_sha256": reference_table_sha256,
            "metrics_file_sha256": metrics_file_sha256,
            "metrics_content_sha256": metrics_content_sha256,
            "current_compact_path": str(compact_lc.resolve()),
            "reference_table_path": str(reference_table.resolve()),
            "reference_product_path": str(reference_product.resolve()),
            "metrics_path": str(metrics_csv.resolve()),
            "current_time_system": "BJD-2457000 converted to BJD",
            "reference_time_system": str(reference["time_system"]),
            "reference_columns": {
                "cadence": str(reference["cadence_column"]),
                "time": str(reference["time_column"]),
                "quality": str(reference["quality_column"]),
                "flux_by_current_aperture": resolved_flux_column_mapping,
                "tic": str(table_identity["identity_source"]["tic"]),
                "sector": str(table_identity["identity_source"]["sector"]),
            },
            "reference_table_identity": table_identity,
            "reference_product_identity": product_identity,
            "external_quality_overlay": {
                **quality_reference.provenance,
                "applied_before_common_cadence_metrics": True,
                "current_compact_full_audit_counts": dict(
                    current_full_overlay.counts
                ),
                "current_common_audit_counts": current_common_audit,
                "reference_common_audit_counts": reference_common_audit,
            },
            "fixed_ephemeris": {
                "period_d": WD1856_PERIOD_D,
                "t0_bjd": WD1856_T0_BJD,
                "transit_duration_min": float(transit_duration_min),
                "oot_exclusion_min": float(oot_exclusion_min),
            },
            "bounded_bls": {
                "period_relative_half_width": float(
                    bls_period_relative_half_width
                ),
                "n_periods": int(bls_n_periods),
                "durations_min": [float(value) for value in bls_durations_min],
                "oversample": int(bls_oversample),
                "min_depth_snr": float(min_bls_depth_snr),
            },
            "n_current_cadences": int(len(current["cadence"])),
            "n_reference_cadences": int(len(reference["cadence"])),
            "n_common_cadences": int(len(common)),
            "max_time_delta_seconds_allowed": float(max_time_delta_seconds),
            "input_size_bytes": {
                "current_compact": initial_input_size[str(compact_lc.resolve())],
                "reference_table": initial_input_size[str(reference_table.resolve())],
                "reference_product": initial_input_size[str(reference_product.resolve())],
                "cadence_reference_table": initial_input_size[
                    str(cadence_reference_table.resolve())
                ],
                "cadence_reference_manifest": initial_input_size[
                    str(cadence_reference_manifest.resolve())
                ],
            },
            "input_mtime_ns": {
                "current_compact": initial_input_mtime_ns[str(compact_lc.resolve())],
                "reference_table": initial_input_mtime_ns[str(reference_table.resolve())],
                "reference_product": initial_input_mtime_ns[str(reference_product.resolve())],
                "cadence_reference_table": initial_input_mtime_ns[
                    str(cadence_reference_table.resolve())
                ],
                "cadence_reference_manifest": initial_input_mtime_ns[
                    str(cadence_reference_manifest.resolve())
                ],
            },
        }
        manifest_tmp.write_text(
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        quality_reference.assert_unchanged()
        final_input_sha256 = {
            str(path): file_sha256(path) for path in input_paths
        }
        if final_input_sha256 != initial_input_sha256:
            raise ValueError(
                "one or more independent-extraction inputs changed during the build"
            )
        metrics_tmp.replace(metrics_csv)
        manifest_tmp.replace(manifest_json)
    except Exception:
        metrics_tmp.unlink(missing_ok=True)
        manifest_tmp.unlink(missing_ok=True)
        raise
    return roundtrip, manifest


__all__ = [
    "INDEPENDENT_CONTRACT_VERSION",
    "INDEPENDENT_COMPARISON_MODE",
    "IndependentExtractionProvenance",
    "build_wd1856_independent_metrics",
    "parse_reference_flux_column_mappings",
]
