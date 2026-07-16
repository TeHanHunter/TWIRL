"""Tier-0 integrity and benchmark QA for an A2v1 sector.

The gate combines full-product schema/coverage, complete BLS accounting, a
deterministic ADP light-curve sample, a raw-extraction comparison, and the S56
WD 1856 BLS benchmark. It is a minimum pilot/integrity gate, not the final
science-acceptance contract: population scatter/cadence regression,
injection-preservation thresholds, and a genuinely independent extraction
comparison remain separate Stage 1 exit criteria.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd

from twirl.lightcurves.tglc_h5_reader import read_tglc_h5
from twirl.vetting.adp_only import ADP_ONLY_APERTURES, ADP_ONLY_CONTRACT_VERSION


A2V1_PHOTOMETRIC_QA_VERSION = "a2v1_tier0_integrity_qa_v2"
WD1856_TIC = 267574918
WD1856_PERIOD_D = 1.407939211
WD1856_T0_BJD = 2458779.375083


def _finite_quantile(values: Iterable[float], q: float) -> float:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    return float(np.quantile(array, q)) if array.size else np.nan


def _finite_median(values: Iterable[float]) -> float:
    return _finite_quantile(values, 0.5)


def _robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return np.nan
    median = float(np.median(values))
    return float(1.4826 * np.median(np.abs(values - median)))


def _normalized(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if not finite.size:
        return np.full(values.shape, np.nan, dtype=float)
    median = float(np.median(finite))
    if not np.isfinite(median) or abs(median) < 1.0e-12:
        return np.full(values.shape, np.nan, dtype=float)
    return values / median


def _correlation(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    good = np.isfinite(left) & np.isfinite(right)
    if np.count_nonzero(good) < 20:
        return np.nan
    x = left[good]
    y = right[good]
    if np.std(x) <= 0 or np.std(y) <= 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"unsupported table format: {path}")


def file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def audit_bls_coverage(
    peaks: pd.DataFrame,
    catalog: pd.DataFrame,
    *,
    expected_product_tag: str = "A2v1",
) -> dict[str, Any]:
    """Require one accounted BLS status for every target/aperture pair."""

    config_columns = (
        "bls_n_periods",
        "bls_n_peaks",
        "bls_p_min_d",
        "bls_p_max_cap_d",
        "bls_max_period_fraction",
        "bls_sigma_clip",
        "bls_orbit_edge_trim_d",
    )
    required_columns = {
        "tic",
        "aperture",
        "status",
        "source_product_tag",
        "adp_only_contract_version",
        *config_columns,
    }
    missing_columns = sorted(required_columns - set(peaks.columns))
    expected_tics = {int(tic) for tic in catalog["tic"]}
    expected_apertures = set(ADP_ONLY_APERTURES)
    expected_pairs = {
        (tic, aperture) for tic in expected_tics for aperture in expected_apertures
    }

    observed_pairs: set[tuple[int, str]] = set()
    observed_tics: set[int] = set()
    observed_apertures: set[str] = set()
    failed_pairs: set[tuple[int, str]] = set()
    if not missing_columns:
        work = peaks.copy()
        work["tic"] = pd.to_numeric(work["tic"], errors="coerce")
        work = work.loc[work["tic"].notna()].copy()
        work["tic"] = work["tic"].astype(np.int64)
        work["aperture"] = work["aperture"].fillna("").astype(str)
        observed_tics = {int(tic) for tic in work["tic"]}
        observed_apertures = {value for value in work["aperture"] if value}
        observed_pairs = {
            (int(row.tic), str(row.aperture))
            for row in work[["tic", "aperture"]].drop_duplicates().itertuples(index=False)
            if row.aperture
        }
        for pair, group in work.groupby(["tic", "aperture"], sort=False):
            if not group["status"].fillna("").astype(str).eq("ok").any():
                failed_pairs.add((int(pair[0]), str(pair[1])))

    product_tags = (
        sorted(peaks["source_product_tag"].dropna().astype(str).unique().tolist())
        if "source_product_tag" in peaks
        else []
    )
    contract_versions = (
        sorted(peaks["adp_only_contract_version"].dropna().astype(str).unique().tolist())
        if "adp_only_contract_version" in peaks
        else []
    )
    config_values = {
        column: sorted(peaks[column].dropna().unique().tolist())
        if column in peaks
        else []
        for column in config_columns
    }
    config_consistent = all(len(values) == 1 for values in config_values.values())
    missing_pairs = sorted(expected_pairs - observed_pairs)
    unexpected_pairs = sorted(observed_pairs - expected_pairs)
    unexpected_tics = sorted(observed_tics - expected_tics)
    unexpected_apertures = sorted(observed_apertures - expected_apertures)
    failed_expected_pairs = failed_pairs & expected_pairs
    passed = bool(
        not missing_columns
        and not missing_pairs
        and not unexpected_pairs
        and not unexpected_tics
        and not unexpected_apertures
        and product_tags == [expected_product_tag]
        and contract_versions == [ADP_ONLY_CONTRACT_VERSION]
        and config_consistent
    )
    return {
        "passed": passed,
        "expected_product_tag": expected_product_tag,
        "observed_product_tags": product_tags,
        "expected_contract_version": ADP_ONLY_CONTRACT_VERSION,
        "observed_contract_versions": contract_versions,
        "config_consistent": config_consistent,
        "config_values": config_values,
        "n_expected_tics": len(expected_tics),
        "n_observed_tics": len(observed_tics),
        "n_expected_pairs": len(expected_pairs),
        "n_observed_pairs": len(observed_pairs & expected_pairs),
        "n_missing_pairs": len(missing_pairs),
        "n_unexpected_pairs": len(unexpected_pairs),
        "n_failed_pairs": len(failed_expected_pairs),
        "failed_pair_fraction": (
            len(failed_expected_pairs) / len(expected_pairs) if expected_pairs else np.nan
        ),
        "missing_columns": missing_columns,
        "missing_pair_examples": [list(item) for item in missing_pairs[:20]],
        "unexpected_pair_examples": [list(item) for item in unexpected_pairs[:20]],
        "unexpected_tic_examples": unexpected_tics[:20],
        "unexpected_apertures": unexpected_apertures,
        "failed_pair_examples": [list(item) for item in sorted(failed_expected_pairs)[:20]],
    }


def compact_target_catalog(compact_lc: Path) -> pd.DataFrame:
    """Read target metadata without loading the cadence arrays."""

    rows: list[dict[str, Any]] = []
    with h5py.File(Path(compact_lc), "r") as h5:
        if "targets" not in h5:
            raise KeyError(f"compact export has no /targets group: {compact_lc}")
        sector = int(h5.attrs.get("sector", -1))
        for key, group in h5["targets"].items():
            rows.append(
                {
                    "tic": int(group.attrs.get("tic", int(key))),
                    "sector": int(group.attrs.get("sector", sector)),
                    "camera": int(group.attrs.get("camera", -1)),
                    "ccd": int(group.attrs.get("ccd", -1)),
                    "tmag": float(group.attrs.get("tessmag", np.nan)),
                }
            )
    catalog = pd.DataFrame(rows).sort_values("tic", kind="stable").reset_index(drop=True)
    if catalog.empty:
        raise ValueError(f"compact export has no targets: {compact_lc}")
    if not catalog["tic"].is_unique:
        raise ValueError("compact export has duplicate TIC groups")
    return catalog


def stratified_target_sample(
    catalog: pd.DataFrame,
    *,
    sample_size: int,
    seed: int,
    magnitude_bins: int = 5,
) -> pd.DataFrame:
    """Take a deterministic approximately equal sample across Tmag ranks."""

    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    work = catalog.copy()
    work["tmag"] = pd.to_numeric(work["tmag"], errors="coerce")
    work = work.sort_values(["tmag", "tic"], na_position="last", kind="stable")
    n_take = min(int(sample_size), len(work))
    if n_take == len(work):
        return work.reset_index(drop=True)

    rng = np.random.default_rng(int(seed))
    bins = [part for part in np.array_split(np.arange(len(work)), magnitude_bins) if len(part)]
    allocation = np.full(len(bins), n_take // len(bins), dtype=int)
    allocation[: n_take % len(bins)] += 1
    selected: list[int] = []
    for indices, count in zip(bins, allocation, strict=True):
        count = min(int(count), len(indices))
        selected.extend(rng.choice(indices, size=count, replace=False).tolist())
    if len(selected) < n_take:
        remaining = np.setdiff1d(np.arange(len(work)), np.asarray(selected), assume_unique=False)
        selected.extend(rng.choice(remaining, size=n_take - len(selected), replace=False).tolist())
    return work.iloc[sorted(selected)].reset_index(drop=True)


def audit_compact_adp(compact_lc: Path, sample: pd.DataFrame) -> pd.DataFrame:
    """Measure cadence retention, precision, and aperture agreement."""

    rows: list[dict[str, Any]] = []
    with h5py.File(Path(compact_lc), "r") as h5:
        for record in sample.itertuples(index=False):
            group_path = f"targets/{int(record.tic):016d}"
            if group_path not in h5:
                raise KeyError(f"sampled TIC is absent from compact export: {record.tic}")
            group = h5[group_path]
            required = ("time", "cadenceno", "quality", "orbitid", *ADP_ONLY_APERTURES)
            missing = [name for name in required if name not in group]
            if missing:
                raise KeyError(f"TIC {record.tic} missing compact datasets: {missing}")
            quality = np.asarray(group["quality"], dtype=np.int64)
            orbitid = np.asarray(group["orbitid"], dtype=np.int64)
            q0 = quality == 0
            row: dict[str, Any] = {
                "tic": int(record.tic),
                "sector": int(record.sector),
                "camera": int(record.camera),
                "ccd": int(record.ccd),
                "tmag": float(record.tmag),
                "n_cadences": int(len(quality)),
                "n_quality0": int(np.count_nonzero(q0)),
                "quality0_fraction": float(np.mean(q0)) if len(q0) else np.nan,
                "n_orbits": int(len(np.unique(orbitid[np.isfinite(orbitid)]))),
            }
            normalized: dict[str, np.ndarray] = {}
            for aperture in ADP_ONLY_APERTURES:
                flux = np.asarray(group[aperture], dtype=float)
                good = q0 & np.isfinite(flux)
                values = _normalized(flux)
                normalized[aperture] = values
                suffix = "small" if aperture.endswith("_SML") else "primary"
                row[f"finite_q0_fraction_{suffix}"] = (
                    float(np.count_nonzero(good) / np.count_nonzero(q0))
                    if np.count_nonzero(q0)
                    else np.nan
                )
                row[f"mad_ppm_{suffix}"] = _robust_mad(values[good]) * 1.0e6
                row[f"negative_q0_fraction_{suffix}"] = (
                    float(np.mean(flux[good] < 0)) if np.count_nonzero(good) else np.nan
                )
            both = q0 & np.isfinite(normalized[ADP_ONLY_APERTURES[0]]) & np.isfinite(normalized[ADP_ONLY_APERTURES[1]])
            small = normalized[ADP_ONLY_APERTURES[0]][both]
            primary = normalized[ADP_ONLY_APERTURES[1]][both]
            row["aperture_correlation"] = _correlation(small, primary)
            row["aperture_difference_mad_ppm"] = _robust_mad(small - primary) * 1.0e6
            small_mad = float(row["mad_ppm_small"])
            primary_mad = float(row["mad_ppm_primary"])
            row["primary_to_small_mad_ratio"] = (
                primary_mad / small_mad
                if np.isfinite(small_mad) and small_mad > 0 and np.isfinite(primary_mad)
                else np.nan
            )
            rows.append(row)
    return pd.DataFrame(rows)


def _raw_path(root: Path, orbit: int, camera: int, ccd: int, tic: int) -> Path:
    return (
        Path(root)
        / f"orbit-{int(orbit)}"
        / "ffi"
        / f"cam{int(camera)}"
        / f"ccd{int(ccd)}"
        / "LC"
        / f"{int(tic)}.h5"
    )


def _compare_raw_target(
    payload: tuple[dict[str, Any], str, str, tuple[int, ...]],
) -> list[dict[str, Any]]:
    record, a2v1_root_s, reference_root_s, orbits = payload
    rows: list[dict[str, Any]] = []
    for orbit in orbits:
        a2v1_path = _raw_path(
            Path(a2v1_root_s), orbit, record["camera"], record["ccd"], record["tic"]
        )
        reference_path = _raw_path(
            Path(reference_root_s),
            orbit,
            record["camera"],
            record["ccd"],
            record["tic"],
        )
        base = {
            "tic": int(record["tic"]),
            "tmag": float(record["tmag"]),
            "orbit": int(orbit),
            "camera": int(record["camera"]),
            "ccd": int(record["ccd"]),
            "a2v1_path": str(a2v1_path),
            "reference_path": str(reference_path),
            "separate_file_identity": np.nan,
        }
        if not a2v1_path.exists() or not reference_path.exists():
            rows.extend(
                {**base, "aperture": aperture, "status": "missing_path"}
                for aperture in ("Small", "Primary")
            )
            continue
        separate_file_identity = not os.path.samefile(a2v1_path, reference_path)
        base["separate_file_identity"] = separate_file_identity
        if not separate_file_identity:
            rows.extend(
                {**base, "aperture": aperture, "status": "same_file_identity"}
                for aperture in ("Small", "Primary")
            )
            continue
        try:
            current = read_tglc_h5(a2v1_path)
            reference = read_tglc_h5(reference_path)
        except Exception as exc:  # Keep a row so failures remain auditable.
            rows.extend(
                {
                    **base,
                    "aperture": aperture,
                    "status": f"read_error:{type(exc).__name__}",
                }
                for aperture in ("Small", "Primary")
            )
            continue
        common, current_index, reference_index = np.intersect1d(
            current.cadence,
            reference.cadence,
            assume_unique=False,
            return_indices=True,
        )
        for aperture in ("Small", "Primary"):
            current_aperture = current.apertures[aperture]
            reference_aperture = reference.apertures[aperture]
            if current_aperture.flux_was_synthesized or reference_aperture.flux_was_synthesized:
                rows.append({**base, "aperture": aperture, "status": "synthesized_flux"})
                continue
            current_flux = current_aperture.raw_flux[current_index]
            reference_flux = reference_aperture.raw_flux[reference_index]
            good = (
                (current.quality[current_index] == 0)
                & (reference.quality[reference_index] == 0)
                & np.isfinite(current_flux)
                & np.isfinite(reference_flux)
            )
            current_norm = _normalized(current_flux[good])
            reference_norm = _normalized(reference_flux[good])
            rows.append(
                {
                    **base,
                    "aperture": aperture,
                    "status": "ok" if np.count_nonzero(good) >= 20 else "insufficient_overlap",
                    "n_common_cadences": int(len(common)),
                    "n_quality0_finite": int(np.count_nonzero(good)),
                    "normalized_flux_correlation": _correlation(current_norm, reference_norm),
                    "normalized_difference_mad_ppm": _robust_mad(current_norm - reference_norm) * 1.0e6,
                    "a2v1_negative_fraction": (
                        float(np.mean(current_flux[good] < 0)) if np.count_nonzero(good) else np.nan
                    ),
                }
            )
    return rows


def compare_raw_extractions(
    sample: pd.DataFrame,
    *,
    a2v1_root: Path,
    reference_root: Path,
    orbits: Sequence[int],
    workers: int = 1,
) -> pd.DataFrame:
    """Compare native raw aperture flux against the prior extraction tree."""

    records = sample.loc[:, ["tic", "tmag", "camera", "ccd"]].to_dict("records")
    payloads = [
        (record, str(a2v1_root), str(reference_root), tuple(int(value) for value in orbits))
        for record in records
    ]
    rows: list[dict[str, Any]] = []
    if int(workers) <= 1:
        batches = map(_compare_raw_target, payloads)
    else:
        executor = ThreadPoolExecutor(max_workers=int(workers))
        batches = executor.map(_compare_raw_target, payloads)
    try:
        for batch in batches:
            rows.extend(batch)
    finally:
        if int(workers) > 1:
            executor.shutdown(wait=True)
    return pd.DataFrame(rows)


def audit_wd1856_bls(peaks: pd.DataFrame) -> dict[str, Any]:
    required = {"tic", "aperture", "peak_rank", "period_d", "t0_bjd", "sde"}
    missing = required - set(peaks.columns)
    if missing:
        raise KeyError(f"BLS peak table is missing columns: {sorted(missing)}")
    rows = peaks.loc[
        pd.to_numeric(peaks["tic"], errors="coerce").eq(WD1856_TIC)
        & pd.to_numeric(peaks["peak_rank"], errors="coerce").eq(1)
        & peaks["aperture"].astype(str).isin(ADP_ONLY_APERTURES)
    ].copy()
    diagnostics: list[dict[str, Any]] = []
    for aperture in ADP_ONLY_APERTURES:
        match = rows.loc[rows["aperture"].astype(str).eq(aperture)]
        if match.empty:
            diagnostics.append({"aperture": aperture, "passed": False, "reason": "missing_top1"})
            continue
        row = match.iloc[0]
        period = float(row["period_d"])
        t0 = float(row["t0_bjd"])
        sde = float(row["sde"])
        phase = (t0 - WD1856_T0_BJD) / WD1856_PERIOD_D
        phase_residual_d = abs((phase + 0.5) % 1.0 - 0.5) * WD1856_PERIOD_D
        period_relative_error = abs(period / WD1856_PERIOD_D - 1.0)
        passed = bool(
            np.isfinite(period_relative_error)
            and period_relative_error <= 0.005
            and np.isfinite(sde)
            and sde >= 10.0
            and np.isfinite(phase_residual_d)
            and phase_residual_d * 1440.0 <= 30.0
        )
        diagnostics.append(
            {
                "aperture": aperture,
                "period_d": period,
                "period_relative_error": period_relative_error,
                "t0_bjd": t0,
                "epoch_residual_min": phase_residual_d * 1440.0,
                "sde": sde,
                "passed": passed,
            }
        )
    return {
        "tic": WD1856_TIC,
        "truth_period_d": WD1856_PERIOD_D,
        "truth_t0_bjd": WD1856_T0_BJD,
        "apertures": diagnostics,
        "passed": bool(
            len(diagnostics) == len(ADP_ONLY_APERTURES)
            and all(row["passed"] for row in diagnostics)
        ),
    }


def _schema_gate(schema: Mapping[str, Any], *, sector: int, n_compact_targets: int) -> dict[str, Any]:
    h5 = schema.get("h5", {})
    fits = schema.get("fits", {})
    observed_sector = int(schema.get("sector", -1))
    n_fits = int(fits.get("n_found_unique_tics", fits.get("n_fits", -1)))
    passed = bool(
        observed_sector == int(sector)
        and schema.get("ok", False)
        and schema.get("ok_h5", False)
        and schema.get("ok_fits", False)
        and int(h5.get("n_missing_h5_non_edge", -1)) == 0
        and int(h5.get("n_zero_byte_h5", -1)) == 0
        and int(fits.get("n_missing_fits_non_edge_tics", -1)) == 0
        and int(fits.get("n_bad_checked_fits", -1)) == 0
        and n_fits == int(n_compact_targets)
    )
    return {
        "passed": passed,
        "sector": observed_sector,
        "n_compact_targets": int(n_compact_targets),
        "n_fits_targets": n_fits,
        "n_present_h5": int(h5.get("n_present_h5", -1)),
        "n_missing_h5_non_edge": int(h5.get("n_missing_h5_non_edge", -1)),
        "n_missing_fits_non_edge": int(fits.get("n_missing_fits_non_edge_tics", -1)),
        "n_bad_checked_fits": int(fits.get("n_bad_checked_fits", -1)),
    }


def evaluate_photometric_gates(
    *,
    sector: int,
    schema: Mapping[str, Any],
    n_compact_targets: int,
    bls_coverage: Mapping[str, Any],
    target_metrics: pd.DataFrame,
    raw_comparison: pd.DataFrame,
    wd1856: Mapping[str, Any] | None,
    reference_benchmark: Mapping[str, Any] | None,
    min_raw_comparisons: int = 50,
) -> dict[str, Any]:
    """Evaluate explicit, auditable thresholds and return the gate payload."""

    finite_columns = ["finite_q0_fraction_small", "finite_q0_fraction_primary"]
    finite_values = target_metrics[finite_columns].to_numpy(dtype=float).ravel()
    valid_target_fraction = float(target_metrics[finite_columns].notna().all(axis=1).mean())
    compact_gate = {
        "passed": bool(
            len(target_metrics) >= min(100, n_compact_targets)
            and valid_target_fraction >= 0.95
            and _finite_median(finite_values) >= 0.95
            and _finite_quantile(finite_values, 0.05) >= 0.50
            and _finite_median(target_metrics["quality0_fraction"]) >= 0.50
        ),
        "n_sampled_targets": int(len(target_metrics)),
        "valid_target_fraction": valid_target_fraction,
        "finite_q0_fraction_median": _finite_median(finite_values),
        "finite_q0_fraction_p05": _finite_quantile(finite_values, 0.05),
        "quality0_fraction_median": _finite_median(target_metrics["quality0_fraction"]),
        "mad_ppm_small_median": _finite_median(target_metrics["mad_ppm_small"]),
        "mad_ppm_primary_median": _finite_median(target_metrics["mad_ppm_primary"]),
    }
    aperture_gate = {
        "passed": bool(
            target_metrics["aperture_correlation"].notna().mean() >= 0.90
            and _finite_median(target_metrics["aperture_correlation"]) >= 0.30
            and 0.10 <= _finite_median(target_metrics["primary_to_small_mad_ratio"]) <= 10.0
        ),
        "valid_correlation_fraction": float(target_metrics["aperture_correlation"].notna().mean()),
        "correlation_median": _finite_median(target_metrics["aperture_correlation"]),
        "primary_to_small_mad_ratio_median": _finite_median(
            target_metrics["primary_to_small_mad_ratio"]
        ),
        "difference_mad_ppm_median": _finite_median(
            target_metrics["aperture_difference_mad_ppm"]
        ),
    }
    if raw_comparison.empty:
        valid_raw = raw_comparison
        expected_raw = 0
    else:
        valid_raw = raw_comparison.loc[
            raw_comparison["status"].astype(str).eq("ok")
            & pd.to_numeric(raw_comparison.get("normalized_flux_correlation"), errors="coerce").notna()
        ]
        expected_raw = int(len(raw_comparison))
    separate_identity = valid_raw.get(
        "separate_file_identity", pd.Series(False, index=valid_raw.index, dtype=bool)
    ).fillna(False).astype(bool)
    raw_gate = {
        "passed": bool(
            len(valid_raw) >= int(min_raw_comparisons)
            and (len(valid_raw) / max(expected_raw, 1)) >= 0.50
            and separate_identity.all()
            and _finite_median(valid_raw.get("normalized_flux_correlation", [])) >= 0.50
        ),
        "n_expected_comparisons": expected_raw,
        "n_valid_comparisons": int(len(valid_raw)),
        "valid_comparison_fraction": float(len(valid_raw) / max(expected_raw, 1)),
        "separate_file_identity_fraction": float(separate_identity.mean()),
        "correlation_median": _finite_median(valid_raw.get("normalized_flux_correlation", [])),
        "difference_mad_ppm_median": _finite_median(valid_raw.get("normalized_difference_mad_ppm", [])),
    }
    if int(sector) == 56:
        benchmark_gate = dict(wd1856 or {"passed": False, "reason": "missing_s56_benchmark"})
    else:
        source_wd = (reference_benchmark or {}).get("benchmarks", {}).get("wd1856", {})
        reference_contract = (reference_benchmark or {}).get("contract_version")
        benchmark_gate = {
            "passed": bool(
                (reference_benchmark or {}).get("passed", False)
                and source_wd.get("passed", False)
                and reference_contract == A2V1_PHOTOMETRIC_QA_VERSION
            ),
            "policy": "S57+ requires the current passed S56 A2v1 Tier-0 reference gate",
            "reference_sector": (reference_benchmark or {}).get("sector"),
            "reference_contract_version": (reference_benchmark or {}).get("contract_version"),
            "wd1856": source_wd,
        }
    gates = {
        "schema_and_completeness": _schema_gate(
            schema, sector=sector, n_compact_targets=n_compact_targets
        ),
        "bls_target_aperture_coverage": dict(bls_coverage),
        "sampled_adp_photometry": compact_gate,
        "aperture_consistency": aperture_gate,
        "prior_extraction_tree_comparison": raw_gate,
        "benchmark": benchmark_gate,
    }
    return {"passed": bool(all(gate.get("passed", False) for gate in gates.values())), "gates": gates}


def summarize_by_magnitude(target_metrics: pd.DataFrame, bins: int = 5) -> list[dict[str, Any]]:
    work = target_metrics.sort_values(["tmag", "tic"], na_position="last").copy()
    work["magnitude_bin"] = pd.qcut(work["tmag"], q=min(int(bins), len(work)), duplicates="drop")
    rows: list[dict[str, Any]] = []
    for label, group in work.groupby("magnitude_bin", observed=True):
        rows.append(
            {
                "bin": str(label),
                "n": int(len(group)),
                "tmag_median": _finite_median(group["tmag"]),
                "mad_ppm_small_median": _finite_median(group["mad_ppm_small"]),
                "mad_ppm_primary_median": _finite_median(group["mad_ppm_primary"]),
                "finite_q0_fraction_small_median": _finite_median(group["finite_q0_fraction_small"]),
                "finite_q0_fraction_primary_median": _finite_median(group["finite_q0_fraction_primary"]),
            }
        )
    return rows


def plot_photometric_qa(target_metrics: pd.DataFrame, output: Path) -> None:
    import matplotlib.pyplot as plt

    from twirl.plotting.style import apply_twirl_style, get_ordered_palette

    style = apply_twirl_style("full_page")
    colors = get_ordered_palette(3, "viridis")
    fig, axes = plt.subplots(1, 2, figsize=style["figsize"], constrained_layout=True)
    axes[0].scatter(
        target_metrics["tmag"],
        target_metrics["mad_ppm_small"],
        s=8,
        alpha=0.55,
        color=colors[0],
        linewidths=0,
        label="ADP small",
    )
    axes[0].scatter(
        target_metrics["tmag"],
        target_metrics["mad_ppm_primary"],
        s=8,
        alpha=0.45,
        color=colors[2],
        linewidths=0,
        label="ADP 3x3",
    )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("TESS magnitude")
    axes[0].set_ylabel("Robust scatter (ppm)")
    axes[0].legend(loc="upper left")

    finite = (
        np.isfinite(target_metrics["mad_ppm_small"])
        & np.isfinite(target_metrics["mad_ppm_primary"])
        & (target_metrics["mad_ppm_small"] > 0)
        & (target_metrics["mad_ppm_primary"] > 0)
    )
    points = axes[1].scatter(
        target_metrics.loc[finite, "mad_ppm_small"],
        target_metrics.loc[finite, "mad_ppm_primary"],
        c=target_metrics.loc[finite, "tmag"],
        cmap="viridis",
        s=9,
        alpha=0.65,
        linewidths=0,
    )
    colorbar = fig.colorbar(points, ax=axes[1], pad=0.02, fraction=0.05)
    colorbar.set_label("TESS magnitude")
    if finite.any():
        paired = target_metrics.loc[finite, ["mad_ppm_small", "mad_ppm_primary"]].to_numpy()
        lo = float(np.nanmin(paired))
        hi = float(np.nanmax(paired))
        axes[1].plot([lo, hi], [lo, hi], color="0.25", linewidth=0.9, linestyle="--")
        axes[1].set_xlim(lo, hi)
        axes[1].set_ylim(lo, hi)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_aspect("equal", adjustable="box")
    axes[1].set_xlabel("ADP small robust scatter (ppm)")
    axes[1].set_ylabel("ADP 3x3 robust scatter (ppm)")
    for axis in axes:
        for spine in axis.spines.values():
            spine.set_visible(True)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def run_a2v1_photometric_qa(
    *,
    sector: int,
    orbits: Sequence[int],
    compact_lc: Path,
    schema_summary: Path,
    bls_peaks: Path,
    a2v1_root: Path,
    reference_raw_root: Path,
    out_dir: Path,
    gate_json: Path,
    sample_size: int = 256,
    seed: int = 56,
    workers: int = 8,
    reference_qa_summary: Path | None = None,
) -> dict[str, Any]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    catalog = compact_target_catalog(compact_lc)
    if set(catalog["sector"].astype(int)) != {int(sector)}:
        raise ValueError(f"compact export sector mismatch: {sorted(set(catalog['sector']))}")
    sample = stratified_target_sample(catalog, sample_size=sample_size, seed=seed)
    target_metrics = audit_compact_adp(compact_lc, sample)
    raw_comparison = compare_raw_extractions(
        sample,
        a2v1_root=a2v1_root,
        reference_root=reference_raw_root,
        orbits=orbits,
        workers=workers,
    )
    schema = json.loads(Path(schema_summary).read_text())
    peaks = read_table(bls_peaks)
    bls_coverage = audit_bls_coverage(peaks, catalog)
    reference_benchmark = (
        json.loads(Path(reference_qa_summary).read_text()) if reference_qa_summary else None
    )
    wd1856 = audit_wd1856_bls(peaks) if int(sector) == 56 else None
    evaluated = evaluate_photometric_gates(
        sector=sector,
        schema=schema,
        n_compact_targets=len(catalog),
        bls_coverage=bls_coverage,
        target_metrics=target_metrics,
        raw_comparison=raw_comparison,
        wd1856=wd1856,
        reference_benchmark=reference_benchmark,
    )
    target_path = out_dir / "sampled_adp_photometry.csv"
    raw_path = out_dir / "prior_extraction_tree_comparison.csv"
    figure_path = out_dir / "a2v1_photometric_qa.png"
    target_metrics.to_csv(target_path, index=False)
    raw_comparison.to_csv(raw_path, index=False)
    plot_photometric_qa(target_metrics, figure_path)

    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "contract_version": A2V1_PHOTOMETRIC_QA_VERSION,
        "sector": int(sector),
        "orbits": [int(orbit) for orbit in orbits],
        "passed": bool(evaluated["passed"]),
        "qa_tier": "tier0_integrity_and_benchmark",
        "science_ready": False,
        "science_gate_status": "tier1_not_implemented",
        "gates": evaluated["gates"],
        "benchmarks": {"wd1856": wd1856} if wd1856 else {},
        "magnitude_trends": summarize_by_magnitude(target_metrics),
        "provenance": {
            "compact_lc": str(compact_lc),
            "schema_summary": str(schema_summary),
            "bls_peaks": str(bls_peaks),
            "bls_peaks_sha256": file_sha256(bls_peaks),
            "a2v1_root": str(a2v1_root),
            "reference_raw_root": str(reference_raw_root),
            "reference_qa_summary": str(reference_qa_summary) if reference_qa_summary else None,
            "sample_size": int(sample_size),
            "sample_seed": int(seed),
            "workers": int(workers),
            "detrended_apertures": list(ADP_ONLY_APERTURES),
            "canonical_det_flux_used": False,
        },
        "outputs": {
            "target_metrics": str(target_path),
            "raw_comparison": str(raw_path),
            "figure": str(figure_path),
            "gate_json": str(gate_json),
        },
    }
    gate_json = Path(gate_json)
    gate_json.parent.mkdir(parents=True, exist_ok=True)
    gate_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary
