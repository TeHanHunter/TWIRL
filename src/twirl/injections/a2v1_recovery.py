"""Fresh, host-disjoint A2v1 light-curve injection evaluation.

This module deliberately separates the immutable injection schedule from the
expensive per-shard light-curve generation.  A schedule fixes one unique host,
one period-radius cell, and independent random seeds for every injection before
any Slurm task starts.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.io.hlsp import BJDREFI
from twirl.lightcurves.detrend_presets import adp03q_config
from twirl.lightcurves.flux_detrend import flux_space_detrend_result
from twirl.search.injections import inject_batman_transit
from twirl.vetting.harmonic_inputs import injected_raw_uncertainty


def fresh_injection_contract(sector: int) -> str:
    """Return the immutable injection-pair contract for one sector."""

    return f"s{int(sector)}_a2v1_fresh_injection_pair_v1"


def schedule_contract(sector: int) -> str:
    """Return the immutable schedule contract for one sector."""

    return f"s{int(sector)}_a2v1_fresh_injection_schedule_v1"


FRESH_INJECTION_CONTRACT = fresh_injection_contract(56)
SCHEDULE_CONTRACT = schedule_contract(56)
RAW_SOURCE_CONTRACT = "s56_tglc_raw_pair_v1"
ADP_APERTURES: tuple[str, str] = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
T_MAG_EDGES: tuple[float, ...] = (-np.inf, 17.0, 18.0, 19.0, np.inf)
T_MAG_LABELS: tuple[str, ...] = (
    "Tmag < 17",
    "17 <= Tmag < 18",
    "18 <= Tmag < 19",
    "Tmag >= 19",
)
WD_RADIUS_REARTH = 0.013 * 109.076
WD_DENSITY_G_CM3 = 3.85e5
G_SI = 6.67430e-11
MAX_TRANSIT_DUTY_CYCLE = 0.20


@dataclass(frozen=True)
class A2V1RecoveryConfig:
    """Locked configuration for one fresh A2v1 recovery evaluation."""

    name: str = "s56_a2v1_teacher_v1_recovery_v1"
    sector: int = 56
    seed: int = 560201
    n_injections: int = 20_000
    n_shards: int = 40
    period_bins: int = 50
    radius_bins: int = 50
    repeats_per_cell: int = 8
    period_min_d: float = 0.12
    period_max_d: float = 13.0
    radius_min_rearth: float = 0.18
    radius_max_rearth: float = 18.0
    cadence_s: float = 200.0
    min_in_transit: int = 2
    min_good_cadences: int = 200
    max_time_error_s: float = 2.0
    parity_sample_size: int = 100
    parity_median_abs_limit: float = 1.0e-4
    parity_scatter_ratio_tolerance: float = 0.01

    def validate(self) -> None:
        if self.sector < 56:
            raise ValueError("A2v1 recovery sectors must be >= 56")
        if (
            self.n_injections
            != self.period_bins * self.radius_bins * self.repeats_per_cell
        ):
            raise ValueError(
                "n_injections must equal period_bins * radius_bins * repeats_per_cell"
            )
        if self.n_injections % self.n_shards:
            raise ValueError("n_injections must divide evenly across n_shards")
        if not 0 < self.period_min_d < self.period_max_d:
            raise ValueError("invalid period range")
        if not 0 < self.radius_min_rearth < self.radius_max_rearth:
            raise ValueError("invalid radius range")
        if self.min_in_transit < 1 or self.min_good_cadences < 1:
            raise ValueError("cadence requirements must be positive")

    @property
    def rows_per_shard(self) -> int:
        return self.n_injections // self.n_shards


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"cannot serialize {type(value).__name__}")


def load_recovery_config(path: Path) -> A2V1RecoveryConfig:
    """Load and strictly validate a YAML or JSON recovery configuration."""

    path = Path(path)
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text())
    else:
        import yaml

        payload = yaml.safe_load(path.read_text())
    if not isinstance(payload, Mapping):
        raise ValueError(f"configuration must contain a mapping: {path}")
    allowed = set(A2V1RecoveryConfig.__dataclass_fields__)
    unknown = sorted(set(payload) - allowed)
    if unknown:
        raise KeyError(f"unknown recovery configuration fields: {unknown}")
    config = A2V1RecoveryConfig(**dict(payload))
    config.validate()
    return config


def _absolute_bjd(time: np.ndarray) -> np.ndarray:
    values = np.asarray(time, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return (
        values + BJDREFI
        if finite.size and float(np.nanmedian(finite)) < 1.0e5
        else values
    )


def _relative_bjd(time: np.ndarray) -> np.ndarray:
    values = np.asarray(time, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return (
        values - BJDREFI
        if finite.size and float(np.nanmedian(finite)) > 1.0e5
        else values
    )


def _seed(base: int, index: int, stream: int) -> int:
    state = np.random.SeedSequence([int(base), int(index), int(stream)]).generate_state(
        1
    )
    return int(state[0])


def _robust_scatter(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not finite.size:
        return float("nan")
    center = float(np.nanmedian(finite))
    return 1.4826 * float(np.nanmedian(np.abs(finite - center)))


def _sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def compare_adp_compact_products(
    reference_h5: Path,
    active_h5: Path,
    *,
    progress_every: int = 1000,
) -> dict[str, Any]:
    """Require cadence and ADP arrays rebuilt from FITS to match the reference."""

    import h5py

    required = (
        "time",
        "cadenceno",
        "quality",
        "orbitid",
        *ADP_APERTURES,
    )
    mismatch_details: list[dict[str, Any]] = []
    mismatched_targets: set[str] = set()
    n_dataset_comparisons = 0
    with h5py.File(reference_h5, "r") as reference, h5py.File(active_h5, "r") as active:
        reference_keys = set(reference["targets"].keys())
        active_keys = set(active["targets"].keys())
        missing_active = sorted(reference_keys - active_keys)
        extra_active = sorted(active_keys - reference_keys)
        for key in sorted(reference_keys & active_keys):
            reference_group = reference[f"targets/{key}"]
            active_group = active[f"targets/{key}"]
            for name in required:
                n_dataset_comparisons += 1
                if name not in reference_group or name not in active_group:
                    mismatched_targets.add(key)
                    if len(mismatch_details) < 100:
                        mismatch_details.append(
                            {
                                "tic_key": key,
                                "dataset": name,
                                "reason": "missing_dataset",
                                "in_reference": name in reference_group,
                                "in_active": name in active_group,
                            }
                        )
                    continue
                reference_dataset = reference_group[name]
                active_dataset = active_group[name]
                if (
                    reference_dataset.shape != active_dataset.shape
                    or reference_dataset.dtype != active_dataset.dtype
                ):
                    mismatched_targets.add(key)
                    if len(mismatch_details) < 100:
                        mismatch_details.append(
                            {
                                "tic_key": key,
                                "dataset": name,
                                "reason": "schema_mismatch",
                                "reference_shape": reference_dataset.shape,
                                "active_shape": active_dataset.shape,
                                "reference_dtype": str(reference_dataset.dtype),
                                "active_dtype": str(active_dataset.dtype),
                            }
                        )
                    continue
                reference_values = np.asarray(reference_dataset)
                active_values = np.asarray(active_dataset)
                equal = np.array_equal(
                    reference_values,
                    active_values,
                    equal_nan=np.issubdtype(reference_values.dtype, np.inexact),
                )
                if not equal:
                    mismatched_targets.add(key)
                    if len(mismatch_details) < 100:
                        detail: dict[str, Any] = {
                            "tic_key": key,
                            "dataset": name,
                            "reason": "value_mismatch",
                        }
                        if np.issubdtype(reference_values.dtype, np.number):
                            delta = np.abs(
                                reference_values.astype(np.float64)
                                - active_values.astype(np.float64)
                            )
                            finite = delta[np.isfinite(delta)]
                            detail["max_abs_difference"] = (
                                float(np.max(finite)) if finite.size else None
                            )
                        mismatch_details.append(detail)
            if (
                progress_every > 0
                and n_dataset_comparisons % (progress_every * len(required)) == 0
            ):
                print(
                    f"[compact-parity] compared {n_dataset_comparisons // len(required):,} targets",
                    flush=True,
                )
    failures = []
    if missing_active:
        failures.append(
            f"active compact product is missing {len(missing_active)} targets"
        )
    if extra_active:
        failures.append(f"active compact product has {len(extra_active)} extra targets")
    if mismatched_targets:
        failures.append(
            f"active compact product differs for {len(mismatched_targets)} targets"
        )
    return {
        "reference_h5": str(Path(reference_h5).resolve()),
        "active_h5": str(Path(active_h5).resolve()),
        "n_reference_targets": len(reference_keys),
        "n_active_targets": len(active_keys),
        "n_targets_compared": len(reference_keys & active_keys),
        "n_dataset_comparisons": n_dataset_comparisons,
        "n_missing_active_targets": len(missing_active),
        "n_extra_active_targets": len(extra_active),
        "n_mismatched_targets": len(mismatched_targets),
        "missing_active_examples": missing_active[:20],
        "extra_active_examples": extra_active[:20],
        "mismatch_details": mismatch_details,
        "passed": not failures,
        "failures": failures,
    }


def _adp_detrend(
    time: np.ndarray,
    raw_flux: np.ndarray,
    raw_error: np.ndarray,
    quality: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    result = flux_space_detrend_result(
        np.asarray(time, dtype=float),
        np.asarray(raw_flux, dtype=float),
        quality=np.asarray(quality),
        flux_err=np.asarray(raw_error, dtype=float),
        cfg=adp03q_config(),
    )
    detrended = np.asarray(result.det_flux, dtype=float)
    good = (np.asarray(quality) == 0) & np.isfinite(detrended)
    center = float(np.nanmedian(detrended[good])) if np.any(good) else np.nan
    if np.isfinite(center):
        detrended = detrended - center + 1.0
    return detrended.astype(np.float32), {
        "fit_count": int(result.fit_count),
        "n_segments": int(result.n_segments),
        "scale": float(result.scale),
        "scale_source": str(result.scale_source),
        "cotrend_status": str(result.cotrend_status),
    }


def _teacher_tics(path: Path) -> set[int]:
    path = Path(path)
    table = (
        pd.read_parquet(path)
        if path.suffix.lower() == ".parquet"
        else pd.read_csv(path, low_memory=False)
    )
    if "tic" not in table:
        raise KeyError(f"Teacher table has no tic column: {path}")
    values = pd.to_numeric(table["tic"], errors="coerce").dropna().astype(np.int64)
    return set(values.tolist())


def _host_qa_rows(
    *,
    raw_h5: Path,
    adp_h5: Path,
    teacher_tics: set[int],
    prior_evaluation_tics: set[int],
    config: A2V1RecoveryConfig,
) -> pd.DataFrame:
    import h5py

    rows: list[dict[str, Any]] = []
    with h5py.File(raw_h5, "r") as raw_file, h5py.File(adp_h5, "r") as adp_file:
        if str(raw_file.attrs.get("contract_version", "")) != RAW_SOURCE_CONTRACT:
            raise ValueError("raw source HDF5 has the wrong contract_version")
        raw_keys = set(raw_file["targets"].keys())
        adp_keys = set(adp_file["targets"].keys())
        keys = sorted(raw_keys & adp_keys)
        for count, key in enumerate(keys, start=1):
            tic = int(key)
            record: dict[str, Any] = {
                "tic": tic,
                "in_raw_source": True,
                "in_adp_reference": True,
                "teacher_tic_excluded": tic in teacher_tics,
                "prior_evaluation_tic_excluded": tic in prior_evaluation_tics,
                "eligible": False,
                "qa_reason": "",
            }
            adp = adp_file[f"targets/{key}"]
            record.update(
                {
                    "sector": int(adp.attrs.get("sector", config.sector)),
                    "camera": int(adp.attrs.get("camera", -1)),
                    "ccd": int(adp.attrs.get("ccd", -1)),
                    "tessmag": float(adp.attrs.get("tessmag", np.nan)),
                    "source_fits": str(adp.attrs.get("source_fits", "")),
                }
            )
            if int(record["sector"]) != int(config.sector):
                record["qa_reason"] = "wrong_sector"
                rows.append(record)
                continue
            if tic in teacher_tics:
                record["qa_reason"] = "teacher_tic"
                rows.append(record)
                continue
            if tic in prior_evaluation_tics:
                record["qa_reason"] = "prior_evaluation_tic"
                rows.append(record)
                continue
            try:
                raw = raw_file[f"targets/{key}"]
                cadence_raw = np.asarray(raw["cadenceno"], dtype=np.int64)
                cadence_adp = np.asarray(adp["cadenceno"], dtype=np.int64)
                if not np.array_equal(cadence_raw, cadence_adp):
                    raise ValueError("cadenceno_mismatch")
                time_delta_s = (
                    np.abs(
                        _absolute_bjd(np.asarray(raw["time"]))
                        - _absolute_bjd(np.asarray(adp["time"]))
                    )
                    * 86400.0
                )
                if (
                    not np.all(np.isfinite(time_delta_s))
                    or float(np.max(time_delta_s)) > config.max_time_error_s
                ):
                    raise ValueError("time_mismatch")
                quality = np.asarray(adp["quality"], dtype=np.int32)
                if not np.array_equal(
                    np.asarray(raw["quality"], dtype=np.int32), quality
                ):
                    raise ValueError("quality_mismatch")
                arrays = {
                    "raw_small": np.asarray(raw["raw_flux_small"], dtype=float),
                    "raw_err_small": np.asarray(raw["raw_flux_err_small"], dtype=float),
                    "raw_primary": np.asarray(raw["raw_flux_primary"], dtype=float),
                    "raw_err_primary": np.asarray(
                        raw["raw_flux_err_primary"], dtype=float
                    ),
                    "adp_small": np.asarray(adp[ADP_APERTURES[0]], dtype=float),
                    "adp_primary": np.asarray(adp[ADP_APERTURES[1]], dtype=float),
                }
                lengths = {len(values) for values in arrays.values()} | {len(quality)}
                if len(lengths) != 1:
                    raise ValueError("length_mismatch")
                good = quality == 0
                for values in arrays.values():
                    good &= np.isfinite(values)
                good &= arrays["raw_err_small"] > 0
                good &= arrays["raw_err_primary"] > 0
                n_good = int(good.sum())
                if n_good < config.min_good_cadences:
                    raise ValueError("too_few_good_cadences")
                baseline_small = float(np.nanmedian(arrays["raw_small"][good]))
                baseline_primary = float(np.nanmedian(arrays["raw_primary"][good]))
                if not np.isfinite(baseline_small) or baseline_small <= 0:
                    raise ValueError("nonpositive_small_baseline")
                if not np.isfinite(baseline_primary) or baseline_primary <= 0:
                    raise ValueError("nonpositive_primary_baseline")
                if not np.isfinite(record["tessmag"]):
                    raise ValueError("missing_tessmag")
                record.update(
                    {
                        "eligible": True,
                        "qa_reason": "ok",
                        "n_cadences": int(len(quality)),
                        "n_good_cadences": n_good,
                        "max_time_delta_s": float(np.max(time_delta_s)),
                        "baseline_small": baseline_small,
                        "baseline_primary": baseline_primary,
                    }
                )
            except (KeyError, ValueError) as exc:
                record["qa_reason"] = str(exc)
            rows.append(record)
            if count % 1000 == 0:
                print(f"[a2v1-host-qa] {count:,}/{len(keys):,}", flush=True)
        for key in sorted(adp_keys - raw_keys):
            adp = adp_file[f"targets/{key}"]
            rows.append(
                {
                    "tic": int(key),
                    "sector": int(adp.attrs.get("sector", config.sector)),
                    "camera": int(adp.attrs.get("camera", -1)),
                    "ccd": int(adp.attrs.get("ccd", -1)),
                    "tessmag": float(adp.attrs.get("tessmag", np.nan)),
                    "source_fits": str(adp.attrs.get("source_fits", "")),
                    "in_raw_source": False,
                    "in_adp_reference": True,
                    "teacher_tic_excluded": int(key) in teacher_tics,
                    "prior_evaluation_tic_excluded": int(key)
                    in prior_evaluation_tics,
                    "eligible": False,
                    "qa_reason": "missing_raw_source",
                }
            )
    return pd.DataFrame(rows)


def _tmag_bin(values: pd.Series) -> pd.Series:
    return pd.cut(
        pd.to_numeric(values, errors="coerce"),
        T_MAG_EDGES,
        right=False,
        labels=T_MAG_LABELS,
    )


def _assign_cells(
    selected: pd.DataFrame,
    *,
    config: A2V1RecoveryConfig,
) -> pd.DataFrame:
    rng = np.random.default_rng(_seed(config.seed, 0, 23))
    n_cells = config.period_bins * config.radius_bins
    capacity = np.full(n_cells, config.repeats_per_cell, dtype=np.int16)
    assigned: list[dict[str, Any]] = []
    work = selected.copy()
    work["tmag_bin"] = _tmag_bin(work["tessmag"]).astype(str)
    for label in T_MAG_LABELS:
        subset = work.loc[work["tmag_bin"].eq(label)].copy()
        order = rng.permutation(len(subset)) if len(subset) else np.empty(0, dtype=int)
        records = subset.iloc[order].to_dict("records")
        cell_order = np.empty(0, dtype=int)
        cursor = 0
        for record in records:
            if cursor >= len(cell_order):
                available = np.flatnonzero(capacity > 0)
                if not len(available):
                    raise RuntimeError("period-radius grid capacity was exhausted")
                cell_order = rng.permutation(available)
                cursor = 0
            cell = int(cell_order[cursor])
            cursor += 1
            capacity[cell] -= 1
            record["grid_cell_index"] = cell
            assigned.append(record)
    if np.any(capacity != 0):
        raise RuntimeError("host assignment did not fill every period-radius cell")
    out = pd.DataFrame(assigned)
    out["grid_slot"] = out.groupby("grid_cell_index", sort=False).cumcount()
    out["injection_index"] = (
        out["grid_cell_index"] * config.repeats_per_cell + out["grid_slot"]
    ).astype(np.int64)
    out = out.sort_values("injection_index", kind="stable").reset_index(drop=True)
    return out


def _a_over_rwd(period_d: float) -> float:
    rho_kg_m3 = WD_DENSITY_G_CM3 * 1000.0
    period_s = float(period_d) * 86400.0
    return float((G_SI * rho_kg_m3 * period_s**2 / (3.0 * np.pi)) ** (1.0 / 3.0))


def _circle_overlap_depth(radius_ratio: float, impact_b: float) -> float:
    r = float(radius_ratio)
    d = float(impact_b)
    if d >= 1.0 + r:
        return 0.0
    if d <= abs(1.0 - r):
        return float(min(1.0, r * r)) if r <= 1.0 else 1.0
    x = np.clip((d * d + 1.0 - r * r) / (2.0 * d), -1.0, 1.0)
    y = np.clip((d * d + r * r - 1.0) / (2.0 * d * r), -1.0, 1.0)
    area = np.arccos(x) + r * r * np.arccos(y)
    radicand = max((-d + 1.0 + r) * (d + 1.0 - r) * (d - 1.0 + r) * (d + 1.0 + r), 0.0)
    area -= 0.5 * np.sqrt(radicand)
    return float(np.clip(area / np.pi, 0.0, 1.0))


def _physical_parameters(
    row: Mapping[str, Any], config: A2V1RecoveryConfig
) -> dict[str, Any]:
    index = int(row["injection_index"])
    cell = int(row["grid_cell_index"])
    period_bin = cell // config.radius_bins
    radius_bin = cell % config.radius_bins
    period_edges = np.geomspace(
        config.period_min_d, config.period_max_d, config.period_bins + 1
    )
    radius_edges = np.geomspace(
        config.radius_min_rearth, config.radius_max_rearth, config.radius_bins + 1
    )
    rng = np.random.default_rng(_seed(config.seed, index, 31))
    period_d = float(
        np.exp(
            rng.uniform(
                np.log(period_edges[period_bin]), np.log(period_edges[period_bin + 1])
            )
        )
    )
    radius_rearth = float(
        np.exp(
            rng.uniform(
                np.log(radius_edges[radius_bin]), np.log(radius_edges[radius_bin + 1])
            )
        )
    )
    radius_rwd = radius_rearth / WD_RADIUS_REARTH
    a_over_rwd = _a_over_rwd(period_d)
    max_b = max(0.05, 1.0 + radius_rwd)
    impact_b = float(rng.uniform(0.0, min(0.95 * max_b, 0.95 * a_over_rwd)))
    chord = float(np.sqrt(max((1.0 + radius_rwd) ** 2 - impact_b**2, 1.0e-6)))
    duration_d = period_d / np.pi * chord / max(a_over_rwd, 1.0e-6)
    duration_min = float(
        np.clip(
            duration_d * 1440.0,
            2.0,
            max(2.0, min(180.0, MAX_TRANSIT_DUTY_CYCLE * period_d * 1440.0)),
        )
    )
    return {
        "grid_period_bin": period_bin,
        "grid_radius_bin": radius_bin,
        "grid_cell_id": f"p{period_bin:02d}_r{radius_bin:02d}",
        "period_bin_lo_d": float(period_edges[period_bin]),
        "period_bin_hi_d": float(period_edges[period_bin + 1]),
        "radius_bin_lo_rearth": float(radius_edges[radius_bin]),
        "radius_bin_hi_rearth": float(radius_edges[radius_bin + 1]),
        "period_d": period_d,
        "radius_rearth": radius_rearth,
        "radius_rwd": radius_rwd,
        "impact_b": impact_b,
        "a_over_rwd": a_over_rwd,
        "inclination_deg": float(
            np.degrees(np.arccos(np.clip(impact_b / a_over_rwd, 0.0, 1.0)))
        ),
        "duration_min": duration_min,
        "geometric_depth": _circle_overlap_depth(radius_rwd, impact_b),
        "epoch_seed": _seed(config.seed, index, 47),
    }


def build_fresh_injection_schedule(
    *,
    raw_h5: Path,
    adp_h5: Path,
    teacher_table: Path,
    config: A2V1RecoveryConfig,
    out_dir: Path,
    additional_exclusion_tables: Sequence[Path] = (),
    host_overlap_audit_tables: Sequence[Path] = (),
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """QA hosts and write the immutable 20k unique-host injection schedule."""

    config.validate()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    teacher_tics = _teacher_tics(teacher_table)
    prior_evaluation_tics: set[int] = set()
    for path in additional_exclusion_tables:
        prior_evaluation_tics.update(_teacher_tics(path))
    prior_evaluation_tics.difference_update(teacher_tics)
    excluded = teacher_tics | prior_evaluation_tics
    qa = _host_qa_rows(
        raw_h5=Path(raw_h5),
        adp_h5=Path(adp_h5),
        teacher_tics=teacher_tics,
        prior_evaluation_tics=prior_evaluation_tics,
        config=config,
    )
    eligible = qa.loc[qa["eligible"].fillna(False).astype(bool)].copy()
    if len(eligible) < config.n_injections:
        raise ValueError(
            f"only {len(eligible):,} host-disjoint QA-passing TICs are available; "
            f"need {config.n_injections:,} without reuse"
        )
    bright = eligible.loc[pd.to_numeric(eligible["tessmag"], errors="coerce") < 19.0]
    faint = eligible.loc[pd.to_numeric(eligible["tessmag"], errors="coerce") >= 19.0]
    if len(bright) > config.n_injections:
        raise ValueError(
            f"{len(bright):,} eligible hosts have Tmag < 19, which exceeds the "
            f"locked {config.n_injections:,}-host sample; refusing to discard bright hosts"
        )
    need = config.n_injections - len(bright)
    selected = pd.concat(
        [
            bright,
            faint.sample(
                n=need,
                random_state=_seed(config.seed, 0, 59),
                replace=False,
            ),
        ],
        ignore_index=True,
    )
    schedule = _assign_cells(selected, config=config)
    physical = pd.DataFrame(
        [_physical_parameters(row, config) for row in schedule.to_dict("records")]
    )
    schedule = pd.concat([schedule.reset_index(drop=True), physical], axis=1)
    injection_contract = fresh_injection_contract(config.sector)
    locked_schedule_contract = schedule_contract(config.sector)
    schedule["injection_id"] = schedule["injection_index"].map(
        lambda value: f"s{config.sector}a2v1_eval_{int(value):06d}"
    )
    schedule["shard_index"] = (
        schedule["injection_index"] // config.rows_per_shard
    ).astype(np.int16)
    schedule["h5_group"] = "/injections/" + schedule["injection_id"]
    schedule["schedule_contract"] = locked_schedule_contract
    schedule["evaluation_only"] = True
    schedule["teacher_training_excluded"] = True
    schedule["native_input_include"] = True
    schedule["source_kind"] = "injected_validation_holdout"
    schedule["is_injected_row"] = True
    schedule["native_group_path"] = "injections/" + schedule["injection_id"]
    if (
        schedule["tic"].duplicated().any()
        or schedule["injection_id"].duplicated().any()
    ):
        raise RuntimeError("fresh schedule reused a host or injection ID")
    if set(schedule["tic"].astype(int)) & excluded:
        raise RuntimeError("fresh schedule overlaps an excluded host")
    support = (
        schedule.groupby(["grid_period_bin", "grid_radius_bin"], as_index=False)
        .size()
        .rename(columns={"size": "n_injections"})
    )
    if len(support) != config.period_bins * config.radius_bins:
        raise RuntimeError("fresh schedule does not cover every period-radius cell")
    if not support["n_injections"].eq(config.repeats_per_cell).all():
        raise RuntimeError("fresh schedule has the wrong per-cell support")
    qa.to_csv(out_dir / "host_qa.csv.gz", index=False, compression="gzip")
    schedule.to_parquet(
        out_dir / "injection_schedule.parquet", compression="zstd", index=False
    )
    schedule.to_csv(
        out_dir / "injection_schedule.csv.gz", index=False, compression="gzip"
    )
    support.to_csv(out_dir / "period_radius_support.csv", index=False)
    selected_tics = set(schedule["tic"].astype(int))
    host_overlap_audits = []
    for path in host_overlap_audit_tables:
        comparison_tics = _teacher_tics(path)
        host_overlap_audits.append(
            {
                "table": str(Path(path).resolve()),
                "n_comparison_unique_tics": len(comparison_tics),
                "n_selected_host_overlap": len(selected_tics & comparison_tics),
            }
        )
    summary = {
        "created_utc": _utcnow(),
        "contract": locked_schedule_contract,
        "injection_contract": injection_contract,
        "config": asdict(config),
        "raw_h5": str(Path(raw_h5).resolve()),
        "adp_h5": str(Path(adp_h5).resolve()),
        "teacher_table": str(Path(teacher_table).resolve()),
        "additional_exclusion_tables": [
            str(Path(path).resolve()) for path in additional_exclusion_tables
        ],
        "host_overlap_audits": host_overlap_audits,
        "n_teacher_table_rows": int(
            len(pd.read_csv(teacher_table, low_memory=False))
            if Path(teacher_table).suffix.lower() != ".parquet"
            else len(pd.read_parquet(teacher_table))
        ),
        "n_teacher_unique_tics": len(teacher_tics),
        "n_prior_evaluation_unique_tics": len(prior_evaluation_tics),
        "n_total_excluded_unique_tics": len(excluded),
        "n_raw_adp_intersection": int(
            qa["in_raw_source"].fillna(False).astype(bool).sum()
        ),
        "n_adp_without_raw": int(
            (~qa["in_raw_source"].fillna(False).astype(bool)).sum()
        ),
        "n_qa_eligible_after_teacher_exclusion": int(len(eligible)),
        "n_qa_eligible_after_all_exclusions": int(len(eligible)),
        "n_selected": int(len(schedule)),
        "n_selected_unique_tics": int(schedule["tic"].nunique()),
        "n_cells": int(len(support)),
        "cell_support_min": int(support["n_injections"].min()),
        "cell_support_max": int(support["n_injections"].max()),
        "tmag_counts": {
            str(key): int(value)
            for key, value in _tmag_bin(schedule["tessmag"])
            .value_counts(sort=False)
            .items()
        },
        "qa_reason_counts": {
            str(key): int(value)
            for key, value in qa["qa_reason"]
            .fillna("")
            .value_counts()
            .sort_index()
            .items()
        },
    }
    (out_dir / "schedule_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n"
    )
    return schedule, qa, summary


def run_adp_roundtrip_parity(
    *,
    raw_h5: Path,
    adp_h5: Path,
    schedule: pd.DataFrame,
    config: A2V1RecoveryConfig,
    out_dir: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Rerun ADP on unmodified raw flux and compare with A2v1 references."""

    import h5py

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work = schedule.copy()
    work["tmag_bin"] = _tmag_bin(work["tessmag"]).astype(str)
    per_bin = max(1, config.parity_sample_size // len(T_MAG_LABELS))
    samples: list[pd.DataFrame] = []
    for index, label in enumerate(T_MAG_LABELS):
        subset = work.loc[work["tmag_bin"].eq(label)]
        n = min(per_bin, len(subset))
        if n:
            samples.append(
                subset.sample(
                    n=n,
                    random_state=_seed(config.seed, index, 71),
                    replace=False,
                )
            )
    sample = pd.concat(samples, ignore_index=True)
    if len(sample) < min(config.parity_sample_size, len(work)):
        remaining = work.loc[~work["tic"].isin(sample["tic"])].sample(
            n=min(config.parity_sample_size, len(work)) - len(sample),
            random_state=_seed(config.seed, 0, 73),
            replace=False,
        )
        sample = pd.concat([sample, remaining], ignore_index=True)
    rows: list[dict[str, Any]] = []
    with h5py.File(raw_h5, "r") as raw_file, h5py.File(adp_h5, "r") as adp_file:
        for record in sample.to_dict("records"):
            key = f"{int(record['tic']):016d}"
            raw = raw_file[f"targets/{key}"]
            adp = adp_file[f"targets/{key}"]
            quality = np.asarray(adp["quality"], dtype=np.int32)
            time = _relative_bjd(np.asarray(raw["time"], dtype=float))
            cadence_match = np.array_equal(
                np.asarray(raw["cadenceno"], dtype=np.int64),
                np.asarray(adp["cadenceno"], dtype=np.int64),
            )
            for aperture, suffix in (
                (ADP_APERTURES[0], "small"),
                (ADP_APERTURES[1], "primary"),
            ):
                rebuilt, diagnostics = _adp_detrend(
                    time,
                    np.asarray(raw[f"raw_flux_{suffix}"], dtype=float),
                    np.asarray(raw[f"raw_flux_err_{suffix}"], dtype=float),
                    quality,
                )
                reference = np.asarray(adp[aperture], dtype=float)
                good = (quality == 0) & np.isfinite(rebuilt) & np.isfinite(reference)
                median_abs = float(
                    np.nanmedian(np.abs(rebuilt[good] - reference[good]))
                )
                scatter_reference = _robust_scatter(reference[good])
                scatter_rebuilt = _robust_scatter(rebuilt[good])
                ratio = (
                    scatter_rebuilt / scatter_reference
                    if np.isfinite(scatter_reference) and scatter_reference > 1.0e-12
                    else 1.0
                )
                passed = bool(
                    cadence_match
                    and median_abs <= config.parity_median_abs_limit
                    and abs(ratio - 1.0) <= config.parity_scatter_ratio_tolerance
                )
                rows.append(
                    {
                        "tic": int(record["tic"]),
                        "tessmag": float(record["tessmag"]),
                        "tmag_bin": record["tmag_bin"],
                        "aperture": aperture,
                        "cadence_match": cadence_match,
                        "n_compared": int(good.sum()),
                        "median_abs_difference": median_abs,
                        "scatter_reference": scatter_reference,
                        "scatter_rebuilt": scatter_rebuilt,
                        "scatter_ratio": ratio,
                        "passed": passed,
                        **{
                            f"detrend_{key}": value
                            for key, value in diagnostics.items()
                        },
                    }
                )
    metrics = pd.DataFrame(rows)
    metrics.to_csv(out_dir / "adp_roundtrip_parity.csv", index=False)
    summary = {
        "created_utc": _utcnow(),
        "n_hosts": int(metrics["tic"].nunique()),
        "n_aperture_checks": int(len(metrics)),
        "n_passed": int(metrics["passed"].sum()),
        "passed": bool(len(metrics) and metrics["passed"].all()),
        "median_abs_limit": config.parity_median_abs_limit,
        "scatter_ratio_tolerance": config.parity_scatter_ratio_tolerance,
        "max_median_abs_difference": float(metrics["median_abs_difference"].max()),
        "scatter_ratio_min": float(metrics["scatter_ratio"].min()),
        "scatter_ratio_max": float(metrics["scatter_ratio"].max()),
    }
    (out_dir / "adp_roundtrip_parity_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return metrics, summary


def _write_dataset(group: Any, name: str, values: np.ndarray) -> None:
    group.create_dataset(name, data=np.asarray(values), compression="lzf", shuffle=True)


def _inject_at_sampled_epoch(
    *,
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    finite_mask: np.ndarray,
    period_d: float,
    duration_min: float,
    radius_rwd: float,
    a_over_rwd: float,
    impact_b: float,
    baseline: float,
    cadence_s: float,
    min_sampled_cadences: int,
    rng: np.random.Generator,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Choose an epoch using finite-exposure model samples, not cadence centers."""

    time = np.asarray(time, dtype=float)
    good = (
        np.isfinite(time)
        & (np.asarray(quality) == 0)
        & np.asarray(finite_mask, dtype=bool)
    )
    good_indices = np.flatnonzero(good)
    if len(good_indices) < min_sampled_cadences:
        raise ValueError("too few good cadences for sampled-model epoch placement")
    good_time = time[good_indices]
    candidates: list[float] = [
        float(rng.uniform(float(np.min(good_time)), float(np.max(good_time))))
        for _ in range(16)
    ]
    adjacent = good_indices[:-1][np.diff(good_indices) == 1]
    if len(adjacent):
        adjacent = rng.permutation(adjacent)
        candidates.extend(
            float(0.5 * (time[index] + time[index + 1])) for index in adjacent[:256]
        )
    candidates.extend(float(value) for value in rng.permutation(good_time)[:128])
    for t0_d in candidates:
        injected, _, model = inject_batman_transit(
            time,
            flux,
            period_d=period_d,
            t0_d=t0_d,
            duration_min=duration_min,
            radius_rstar=radius_rwd,
            a_over_rstar=a_over_rwd,
            impact_b=impact_b,
            baseline=baseline,
            exp_time_d=cadence_s / 86400.0,
        )
        sampled = good & np.isfinite(model) & (np.asarray(model) < 1.0 - 1.0e-10)
        if int(sampled.sum()) >= int(min_sampled_cadences):
            return (
                t0_d,
                np.asarray(injected, dtype=float),
                np.asarray(model, dtype=float),
            )
    raise ValueError(
        f"could not place finite-exposure transit on {min_sampled_cadences} good cadences"
    )


def _attr_value(value: Any) -> Any:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return ""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, bytes, int, float, bool)):
        return value
    return json.dumps(value, default=_json_default)


def write_fresh_injection_shard(
    *,
    raw_h5: Path,
    adp_h5: Path,
    schedule: pd.DataFrame,
    shard_index: int,
    config: A2V1RecoveryConfig,
    out_h5: Path,
    limit: int | None = None,
    selection_mode: str = "shard",
) -> dict[str, Any]:
    """Generate one self-contained raw/error/ADP injection shard."""

    import h5py

    config.validate()
    contract = fresh_injection_contract(config.sector)
    if shard_index < 0 or shard_index >= config.n_shards:
        raise ValueError("invalid shard_index")
    if selection_mode == "shard":
        rows = schedule.loc[
            pd.to_numeric(schedule["shard_index"], errors="coerce").eq(
                int(shard_index)
            )
        ].copy()
        rows = rows.sort_values("injection_index", kind="stable")
        if limit is not None:
            rows = rows.head(int(limit))
    elif selection_mode == "parameter_spanning":
        if limit is None:
            raise ValueError("parameter_spanning selection requires limit")
        rows = select_parameter_spanning_rows(schedule, n_rows=int(limit))
    else:
        raise ValueError(f"unsupported selection_mode: {selection_mode}")
    expected = (
        config.rows_per_shard
        if limit is None
        else min(config.rows_per_shard, int(limit))
    )
    if len(rows) != expected:
        raise ValueError(
            f"shard {shard_index} has {len(rows)} schedule rows; expected {expected}"
        )
    out_h5 = Path(out_h5)
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_h5.with_suffix(out_h5.suffix + ".tmp")
    manifest_rows: list[dict[str, Any]] = []
    try:
        with (
            h5py.File(raw_h5, "r") as raw_file,
            h5py.File(adp_h5, "r") as adp_file,
            h5py.File(temporary, "w") as output,
        ):
            output.attrs["contract_version"] = contract
            output.attrs["created_utc"] = _utcnow()
            output.attrs["config"] = json.dumps(asdict(config), sort_keys=True)
            output.attrs["source_raw_h5"] = str(Path(raw_h5).resolve())
            output.attrs["source_adp_h5"] = str(Path(adp_h5).resolve())
            output.attrs["source_injection_h5"] = str(out_h5.resolve())
            output.attrs["shard_index"] = int(shard_index)
            output.attrs["n_shards"] = int(config.n_shards)
            output.attrs["n_injections"] = int(len(rows))
            output.attrs["selection_mode"] = str(selection_mode)
            output.attrs["apertures"] = json.dumps(list(ADP_APERTURES))
            root = output.create_group("injections")
            for count, row in enumerate(rows.to_dict("records"), start=1):
                tic = int(row["tic"])
                key = f"{tic:016d}"
                raw = raw_file[f"targets/{key}"]
                adp = adp_file[f"targets/{key}"]
                cadence_raw = np.asarray(raw["cadenceno"], dtype=np.int64)
                cadence_adp = np.asarray(adp["cadenceno"], dtype=np.int64)
                if not np.array_equal(cadence_raw, cadence_adp):
                    raise ValueError(f"TIC {tic}: cadence mismatch during shard build")
                raw_time = _absolute_bjd(np.asarray(raw["time"], dtype=float))
                adp_time = _relative_bjd(np.asarray(adp["time"], dtype=float))
                delta_s = np.abs(raw_time - _absolute_bjd(adp_time)) * 86400.0
                if float(np.nanmax(delta_s)) > config.max_time_error_s:
                    raise ValueError(f"TIC {tic}: time mismatch during shard build")
                quality = np.asarray(adp["quality"], dtype=np.int32)
                orbitid = np.asarray(adp["orbitid"], dtype=np.int32)
                raw_small = np.asarray(raw["raw_flux_small"], dtype=float)
                err_small = np.asarray(raw["raw_flux_err_small"], dtype=float)
                raw_primary = np.asarray(raw["raw_flux_primary"], dtype=float)
                err_primary = np.asarray(raw["raw_flux_err_primary"], dtype=float)
                baseline_small = float(row["baseline_small"])
                baseline_primary = float(row["baseline_primary"])
                finite = (
                    np.isfinite(raw_small)
                    & np.isfinite(raw_primary)
                    & np.isfinite(err_small)
                    & np.isfinite(err_primary)
                )
                epoch_rng = np.random.default_rng(int(row["epoch_seed"]))
                t0_d, injected_small, model = _inject_at_sampled_epoch(
                    time=adp_time,
                    flux=raw_small,
                    quality=quality,
                    finite_mask=finite,
                    period_d=float(row["period_d"]),
                    duration_min=float(row["duration_min"]),
                    radius_rwd=float(row["radius_rwd"]),
                    a_over_rwd=float(row["a_over_rwd"]),
                    impact_b=float(row["impact_b"]),
                    baseline=baseline_small,
                    cadence_s=config.cadence_s,
                    min_sampled_cadences=config.min_in_transit,
                    rng=epoch_rng,
                )
                transit_mask = finite & (quality == 0) & (model < 1.0 - 1.0e-10)
                injected_primary = raw_primary + baseline_primary * (model - 1.0)
                injected_err_small = injected_raw_uncertainty(
                    err_small,
                    model,
                    source_flux_rate=baseline_small,
                    cadence_s=config.cadence_s,
                )
                injected_err_primary = injected_raw_uncertainty(
                    err_primary,
                    model,
                    source_flux_rate=baseline_primary,
                    cadence_s=config.cadence_s,
                )
                det_small, diag_small = _adp_detrend(
                    adp_time, injected_small, injected_err_small, quality
                )
                det_primary, diag_primary = _adp_detrend(
                    adp_time, injected_primary, injected_err_primary, quality
                )
                injection_id = str(row["injection_id"])
                group = root.create_group(injection_id)
                attrs = {
                    **row,
                    "injection_id": injection_id,
                    "sector": config.sector,
                    "tic": tic,
                    "t0_d": t0_d,
                    "t0_bjd": t0_d + BJDREFI,
                    "depth": float(1.0 - np.nanmin(model)),
                    "model_depth": float(1.0 - np.nanmin(model)),
                    "sampled_model_depth": float(1.0 - np.nanmin(model)),
                    "n_good_in_transit": int(
                        np.count_nonzero(transit_mask & (quality == 0) & finite)
                    ),
                    "cadence_s": config.cadence_s,
                    "injection_level": "raw_flux_pre_detrend",
                    "injection_model": "batman_quadratic",
                    "duration_model": "wd_density_batman",
                    "epoch_sampling_policy": "finite_exposure_model_min_sampled_cadences",
                    "sampling_mode": "period_radius_grid",
                    "baseline_source": "raw_median_per_aperture",
                    "injection_baseline_Small": baseline_small,
                    "injection_baseline_Primary": baseline_primary,
                    "raw_baseline_Small": baseline_small,
                    "raw_baseline_Primary": baseline_primary,
                    "aperture": ADP_APERTURES[0],
                    "apertures": json.dumps(list(ADP_APERTURES)),
                    "source_injection_h5": str(out_h5.resolve()),
                    "source_raw_h5": str(Path(raw_h5).resolve()),
                    "source_adp_h5": str(Path(adp_h5).resolve()),
                    "contract_version": contract,
                    "raw_adp_time_delta_max_s": float(np.nanmax(delta_s)),
                }
                for name, value in attrs.items():
                    group.attrs[str(name)] = _attr_value(value)
                for prefix, diagnostics in (
                    (ADP_APERTURES[0], diag_small),
                    (ADP_APERTURES[1], diag_primary),
                ):
                    for name, value in diagnostics.items():
                        group.attrs[f"{prefix}_{name}"] = _attr_value(value)
                datasets = {
                    "time": adp_time,
                    "cadenceno": cadence_adp,
                    "orbitid": orbitid,
                    "quality": quality,
                    "transit_model": model.astype(np.float32),
                    "RAW_FLUX_Small_original": raw_small,
                    "RAW_FLUX_Small_injected": injected_small,
                    "RAW_FLUX_ERR_Small_original": err_small,
                    "RAW_FLUX_ERR_Small_injected": injected_err_small,
                    "RAW_FLUX_Primary_original": raw_primary,
                    "RAW_FLUX_Primary_injected": injected_primary,
                    "RAW_FLUX_ERR_Primary_original": err_primary,
                    "RAW_FLUX_ERR_Primary_injected": injected_err_primary,
                    f"{ADP_APERTURES[0]}_original": np.asarray(adp[ADP_APERTURES[0]]),
                    f"{ADP_APERTURES[0]}_injected": det_small,
                    f"{ADP_APERTURES[1]}_original": np.asarray(adp[ADP_APERTURES[1]]),
                    f"{ADP_APERTURES[1]}_injected": det_primary,
                    "flux_injected": det_small,
                }
                for name, values in datasets.items():
                    _write_dataset(group, name, values)
                manifest_rows.append(
                    {
                        **row,
                        "t0_d": t0_d,
                        "t0_bjd": t0_d + BJDREFI,
                        "depth": float(1.0 - np.nanmin(model)),
                        "model_depth": float(1.0 - np.nanmin(model)),
                        "sampled_model_depth": float(1.0 - np.nanmin(model)),
                        "n_good_in_transit": attrs["n_good_in_transit"],
                        "source_h5": str(out_h5.resolve()),
                        "h5_group": f"/injections/{injection_id}",
                    }
                )
                if count % 50 == 0:
                    print(f"[a2v1-injection-shard] {count:,}/{len(rows):,}", flush=True)
        temporary.replace(out_h5)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    manifest = pd.DataFrame(manifest_rows)
    manifest_path = out_h5.with_name(out_h5.stem + "_manifest.parquet")
    manifest.to_parquet(manifest_path, compression="zstd", index=False)
    summary = {
        "created_utc": _utcnow(),
        "contract": contract,
        "config": asdict(config),
        "shard_index": int(shard_index),
        "n_injections": int(len(manifest)),
        "n_unique_tics": int(manifest["tic"].nunique()),
        "selection_mode": str(selection_mode),
        "out_h5": str(out_h5),
        "manifest": str(manifest_path),
        "sha256": _sha256_file(out_h5),
    }
    out_h5.with_name(out_h5.stem + "_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n"
    )
    return summary


def select_parameter_spanning_rows(
    schedule: pd.DataFrame,
    *,
    n_rows: int,
) -> pd.DataFrame:
    """Select a deterministic smoke sample spanning both grid dimensions."""

    required = {"grid_period_bin", "grid_radius_bin", "injection_index", "tic"}
    missing = sorted(required - set(schedule.columns))
    if missing:
        raise KeyError(f"schedule is missing smoke-selection columns: {missing}")
    if n_rows < 1:
        raise ValueError("n_rows must be positive")

    frame = schedule.copy()
    frame["grid_period_bin"] = pd.to_numeric(
        frame["grid_period_bin"], errors="raise"
    ).astype(int)
    frame["grid_radius_bin"] = pd.to_numeric(
        frame["grid_radius_bin"], errors="raise"
    ).astype(int)
    period_bins = sorted(frame["grid_period_bin"].unique().tolist())
    radius_bins = sorted(frame["grid_radius_bin"].unique().tolist())
    if not period_bins or not radius_bins:
        raise ValueError("schedule has no populated period-radius cells")

    # A coprime stride makes each pass a permutation of the radius bins. The
    # pass offset prevents repeated cells while keeping both marginals even.
    stride = len(radius_bins) - 1 if len(radius_bins) > 1 else 1
    while stride > 1 and np.gcd(stride, len(radius_bins)) != 1:
        stride -= 1
    groups = {
        (int(period_bin), int(radius_bin)): group.sort_values(
            "injection_index", kind="stable"
        )
        for (period_bin, radius_bin), group in frame.groupby(
            ["grid_period_bin", "grid_radius_bin"], sort=False
        )
    }
    selections: list[pd.Series] = []
    use_count: dict[tuple[int, int], int] = {}
    for index in range(int(n_rows)):
        pass_index, period_offset = divmod(index, len(period_bins))
        period_bin = int(period_bins[period_offset])
        radius_offset = (stride * period_offset + pass_index) % len(radius_bins)
        radius_bin = int(radius_bins[radius_offset])
        key = (period_bin, radius_bin)
        cell = groups.get(key)
        if cell is None or cell.empty:
            raise ValueError(f"schedule is missing required smoke cell {key}")
        occurrence = use_count.get(key, 0)
        if occurrence >= len(cell):
            raise ValueError(
                f"smoke selection requests more rows than available in cell {key}"
            )
        selections.append(cell.iloc[occurrence])
        use_count[key] = occurrence + 1
    rows = pd.DataFrame(selections).reset_index(drop=True)
    if rows["injection_index"].duplicated().any() or rows["tic"].duplicated().any():
        raise ValueError("parameter-spanning smoke selection is not unique")
    return rows


__all__ = [
    "ADP_APERTURES",
    "A2V1RecoveryConfig",
    "FRESH_INJECTION_CONTRACT",
    "SCHEDULE_CONTRACT",
    "fresh_injection_contract",
    "schedule_contract",
    "build_fresh_injection_schedule",
    "load_recovery_config",
    "run_adp_roundtrip_parity",
    "select_parameter_spanning_rows",
    "write_fresh_injection_shard",
]
