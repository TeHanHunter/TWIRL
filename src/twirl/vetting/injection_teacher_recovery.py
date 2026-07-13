"""ADP-only injection recovery inputs and Teacher-v1 evaluation helpers."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.io.hlsp import HLSPLightCurve
from twirl.lightcurves.detrend_presets import adp03q_config
from twirl.lightcurves.flux_detrend import flux_space_detrend_result
from twirl.vetting.adp_only import (
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
    classify_period_relation,
)
from twirl.vetting.harmonic_inputs import (
    A2V1_TEACHER_INPUT_CONTRACT,
    injected_raw_uncertainty,
    native_group_path,
)
from twirl.vetting.two_aperture import measure_two_aperture_candidate_metadata


INJECTION_ADP_PAIR_CONTRACT = "s56_predetrend_raw_adp_pair_v1"
TEACHER_RECOVERY_POLICY = "fresh_host_disjoint_a2v1_adp_small_top5"
PEAK_FIELDS: tuple[str, ...] = (
    "peak_rank",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde",
    "log_power",
)


def _bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "1.0", "true", "t", "yes", "y"})
    )


def _absolute_bjd(time: np.ndarray) -> np.ndarray:
    values = np.asarray(time, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return (
        values + 2457000.0
        if finite.size and float(np.nanmedian(finite)) < 1.0e5
        else values
    )


def _stable_shard(value: str, n_shards: int) -> int:
    digest = hashlib.sha1(str(value).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % int(n_shards)


def build_teacher_injection_holdout(
    manifest: pd.DataFrame,
    teacher_training_rows: pd.DataFrame,
    *,
    n_shards: int = 40,
    min_cell_support: int = 4,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Exclude every training injection and every Teacher-v1 host TIC."""

    required = {
        "injection_id",
        "tic",
        "tessmag",
        "grid_period_bin",
        "grid_radius_bin",
    }
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise KeyError(f"injection manifest is missing columns: {missing}")
    if n_shards < 1:
        raise ValueError("n_shards must be positive")
    work = manifest.copy()
    work["injection_id"] = work["injection_id"].fillna("").astype(str)
    if work["injection_id"].eq("").any() or work["injection_id"].duplicated().any():
        raise ValueError("injection manifest IDs must be nonempty and unique")
    work["tic"] = pd.to_numeric(work["tic"], errors="coerce")
    if work["tic"].isna().any():
        raise ValueError("injection manifest contains invalid TIC values")
    work["tic"] = work["tic"].astype(np.int64)

    training = teacher_training_rows.copy()
    injected = (
        _bool_series(training["is_injected_row"])
        if "is_injected_row" in training
        else training.get("source_kind", pd.Series("", index=training.index))
        .fillna("")
        .astype(str)
        .str.contains("inject", case=False)
    )
    training_ids = set(
        training.loc[injected, "injection_id"].dropna().astype(str)
        if "injection_id" in training
        else []
    )
    training_tics = set(
        pd.to_numeric(training.get("tic"), errors="coerce").dropna().astype(np.int64)
    )
    injected_training_tics = set(
        pd.to_numeric(training.loc[injected, "tic"], errors="coerce")
        .dropna()
        .astype(np.int64)
        if "tic" in training
        else []
    )
    work["teacher_training_injection_id"] = work["injection_id"].isin(training_ids)
    work["teacher_training_host_tic"] = work["tic"].isin(training_tics)
    work["teacher_injected_training_host_tic"] = work["tic"].isin(
        injected_training_tics
    )
    work["evaluation_holdout"] = (
        ~work["teacher_training_injection_id"] & ~work["teacher_training_host_tic"]
    )
    retained = work.loc[work["evaluation_holdout"]].copy()
    retained["evaluation_shard"] = retained["injection_id"].map(
        lambda value: _stable_shard(value, n_shards)
    )
    retained["teacher_recovery_policy"] = TEACHER_RECOVERY_POLICY

    cell_columns = ["grid_period_bin", "grid_radius_bin"]
    original_counts = (
        work.groupby(cell_columns, dropna=False).size().rename("n_original")
    )
    retained_counts = (
        retained.groupby(cell_columns, dropna=False).size().rename("n_retained")
    )
    support = (
        pd.concat([original_counts, retained_counts], axis=1).fillna(0).reset_index()
    )
    support[["n_original", "n_retained"]] = support[
        ["n_original", "n_retained"]
    ].astype(int)
    support["retained_fraction"] = support["n_retained"] / support[
        "n_original"
    ].replace(0, np.nan)
    support["meets_min_support"] = support["n_retained"].ge(int(min_cell_support))
    if not support["meets_min_support"].all():
        bad = support.loc[~support["meets_min_support"]]
        raise ValueError(
            f"held-out injection grid has {len(bad)} cells below {min_cell_support} rows"
        )

    tmag = pd.to_numeric(retained["tessmag"], errors="coerce")
    tmag_labels = pd.cut(
        tmag,
        [-np.inf, 17.0, 18.0, 19.0, np.inf],
        right=False,
        labels=["Tmag < 17", "17 <= Tmag < 18", "18 <= Tmag < 19", "Tmag >= 19"],
    )
    tmag_counts = tmag_labels.value_counts(sort=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy": TEACHER_RECOVERY_POLICY,
        "n_manifest": int(len(work)),
        "n_teacher_training_injection_ids": int(len(training_ids)),
        "n_excluded_exact_training_ids": int(
            work["teacher_training_injection_id"].sum()
        ),
        "n_retained": int(len(retained)),
        "n_retained_unique_tics": int(retained["tic"].nunique()),
        "n_period_radius_cells": int(len(support)),
        "n_cells_meeting_min_support": int(support["meets_min_support"].sum()),
        "min_cell_support_required": int(min_cell_support),
        "retained_cell_count_min": int(support["n_retained"].min()),
        "retained_cell_count_median": float(support["n_retained"].median()),
        "retained_cell_count_max": int(support["n_retained"].max()),
        "tmag_bin_counts": {str(key): int(value) for key, value in tmag_counts.items()},
        "n_shards": int(n_shards),
        "n_retained_on_any_teacher_table_host": int(
            retained["teacher_training_host_tic"].sum()
        ),
        "n_retained_on_injected_training_host": int(
            retained["teacher_injected_training_host_tic"].sum()
        ),
    }
    return retained.reset_index(drop=True), support, summary


def subset_raw_source_h5(
    *,
    source_h5: Path,
    tics: Sequence[int],
    out_h5: Path,
) -> dict[str, Any]:
    """Copy only requested host groups from a compact raw-source export."""

    import h5py

    requested = sorted({int(value) for value in tics})
    out_h5 = Path(out_h5)
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_h5.with_suffix(out_h5.suffix + ".tmp")
    missing: list[int] = []
    with h5py.File(source_h5, "r") as source, h5py.File(temporary, "w") as output:
        for key, value in source.attrs.items():
            output.attrs[str(key)] = value
        output.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        output.attrs["subset_source_h5"] = str(source_h5)
        output.attrs["n_requested_tics"] = len(requested)
        targets = output.create_group("targets")
        for index, tic in enumerate(requested, start=1):
            path = f"targets/{tic:016d}"
            if path not in source:
                missing.append(tic)
                continue
            source.copy(source[path], targets, name=f"{tic:016d}")
            if index % 500 == 0:
                print(f"[raw-source-subset] {index:,}/{len(requested):,}", flush=True)
    if missing:
        temporary.unlink(missing_ok=True)
        raise ValueError(
            f"raw source is missing {len(missing)} requested TICs; first={missing[:10]}"
        )
    temporary.replace(out_h5)
    return {
        "source_h5": str(source_h5),
        "out_h5": str(out_h5),
        "n_requested_tics": len(requested),
        "n_written_tics": len(requested),
    }


def align_raw_source_by_time(
    raw: Mapping[str, np.ndarray],
    *,
    injection_time: np.ndarray,
    max_time_error_s: float = 2.0,
) -> tuple[dict[str, np.ndarray], float]:
    """Align compact raw/error arrays to an injection time grid."""

    raw_time = _absolute_bjd(np.asarray(raw["time"], dtype=float))
    requested_time = _absolute_bjd(np.asarray(injection_time, dtype=float))
    if len(raw_time) == len(requested_time):
        direct_delta = np.abs(raw_time - requested_time) * 86400.0
        if (
            np.all(np.isfinite(direct_delta))
            and float(np.max(direct_delta)) <= max_time_error_s
        ):
            return {name: np.asarray(values) for name, values in raw.items()}, float(
                np.max(direct_delta)
            )
    order = np.argsort(raw_time)
    sorted_time = raw_time[order]
    right = np.searchsorted(sorted_time, requested_time, side="left")
    right = np.clip(right, 0, max(len(sorted_time) - 1, 0))
    left = np.clip(right - 1, 0, max(len(sorted_time) - 1, 0))
    right_delta = np.abs(sorted_time[right] - requested_time)
    left_delta = np.abs(sorted_time[left] - requested_time)
    nearest_sorted = np.where(left_delta <= right_delta, left, right)
    index = order[nearest_sorted]
    if len(np.unique(index)) != len(index):
        raise ValueError("raw/injection time alignment is not one-to-one")
    delta_s = np.abs(raw_time[index] - requested_time) * 86400.0
    if not np.all(np.isfinite(delta_s)) or float(np.max(delta_s)) > float(
        max_time_error_s
    ):
        raise ValueError(
            f"raw/injection timestamps differ by up to {np.nanmax(delta_s):.3f} s"
        )
    return {name: np.asarray(values)[index] for name, values in raw.items()}, float(
        np.max(delta_s)
    )


def _adp_detrend(
    *,
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
    diagnostics = {
        "fit_count": int(result.fit_count),
        "n_segments": int(result.n_segments),
        "scale": float(result.scale),
        "scale_source": str(result.scale_source),
        "cotrend_status": str(result.cotrend_status),
    }
    return detrended.astype(np.float32), diagnostics


def rebuild_injected_adp_pair(
    *,
    canonical_injection_h5: Path,
    raw_source_h5: Path,
    injection_ids: Sequence[str],
    out_h5: Path,
) -> dict[str, Any]:
    """Rerun ADP03q on canonical pre-detrend raw-flux injections."""

    import h5py

    ids = [str(value) for value in injection_ids]
    if len(ids) != len(set(ids)):
        raise ValueError("injection_ids contains duplicates")
    out_h5 = Path(out_h5)
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_h5.with_suffix(out_h5.suffix + ".tmp")
    with (
        h5py.File(canonical_injection_h5, "r") as canonical,
        h5py.File(raw_source_h5, "r") as raw_file,
        h5py.File(temporary, "w") as output,
    ):
        output.attrs["contract_version"] = INJECTION_ADP_PAIR_CONTRACT
        output.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        output.attrs["source_injection_h5"] = str(
            Path(canonical_injection_h5).resolve()
        )
        output.attrs["raw_source_h5"] = str(Path(raw_source_h5).resolve())
        output.attrs["aperture"] = ADP_ONLY_APERTURES[0]
        output.attrs["apertures"] = json.dumps(list(ADP_ONLY_APERTURES))
        output.attrs["detrend_preset"] = "twirl-fs-v2-adp03q"
        output.attrs["injection_uncertainty"] = "source_poisson_plus_fixed_floor"
        root = output.create_group("injections")
        for index, injection_id in enumerate(ids, start=1):
            source_path = f"injections/{injection_id}"
            if source_path not in canonical:
                raise KeyError(f"canonical injection product is missing {source_path}")
            source = canonical[source_path]
            tic = int(source.attrs["tic"])
            raw_path = f"targets/{tic:016d}"
            if raw_path not in raw_file:
                raise KeyError(f"compact raw source is missing {raw_path}")
            raw_payload = {
                name: np.asarray(raw_file[raw_path][name])
                for name in raw_file[raw_path]
            }
            time = np.asarray(source["time"], dtype=np.float64)
            aligned, max_delta_s = align_raw_source_by_time(
                raw_payload,
                injection_time=time,
            )
            quality = np.asarray(source["quality"], dtype=np.int32)
            orbitid = np.asarray(source["orbitid"], dtype=np.int32)
            model = np.asarray(source["transit_model"], dtype=np.float32)
            cadence_s = float(source.attrs.get("cadence_s", 200.0))
            destination = root.create_group(injection_id)
            for key, value in source.attrs.items():
                destination.attrs[str(key)] = value
            destination.attrs["source_injection_h5"] = str(
                Path(canonical_injection_h5).resolve()
            )
            destination.attrs["raw_source_h5"] = str(Path(raw_source_h5).resolve())
            destination.attrs["adp_pair_contract_version"] = INJECTION_ADP_PAIR_CONTRACT
            destination.attrs["raw_time_delta_max_s"] = max_delta_s
            destination.attrs["aperture"] = ADP_ONLY_APERTURES[0]
            destination.attrs["apertures"] = json.dumps(list(ADP_ONLY_APERTURES))
            common = {
                "time": time,
                "cadenceno": np.asarray(aligned["cadenceno"], dtype=np.int64),
                "orbitid": orbitid,
                "quality": quality,
                "transit_model": model,
            }
            for name, values in common.items():
                destination.create_dataset(
                    name, data=values, compression="lzf", shuffle=True
                )
            for aperture, suffix, dataset in (
                ("Small", "small", ADP_ONLY_APERTURES[0]),
                ("Primary", "primary", ADP_ONLY_APERTURES[1]),
            ):
                raw_flux = np.asarray(
                    source[f"RAW_FLUX_{aperture}_injected"], dtype=float
                )
                raw_error = np.asarray(aligned[f"raw_flux_err_{suffix}"], dtype=float)
                baseline = float(source.attrs[f"injection_baseline_{aperture}"])
                injected_error = injected_raw_uncertainty(
                    raw_error,
                    model,
                    source_flux_rate=baseline,
                    cadence_s=cadence_s,
                )
                detrended, diagnostics = _adp_detrend(
                    time=time,
                    raw_flux=raw_flux,
                    raw_error=injected_error,
                    quality=quality,
                )
                destination.create_dataset(
                    f"{dataset}_injected",
                    data=detrended,
                    compression="lzf",
                    shuffle=True,
                )
                destination.create_dataset(
                    f"RAW_FLUX_ERR_{aperture}_injected",
                    data=injected_error,
                    compression="lzf",
                    shuffle=True,
                )
                for key, value in diagnostics.items():
                    destination.attrs[f"{dataset}_{key}"] = value
            destination.create_dataset(
                "flux_injected",
                data=np.asarray(destination[f"{ADP_ONLY_APERTURES[0]}_injected"]),
                compression="lzf",
                shuffle=True,
            )
            if index % 50 == 0:
                print(f"[injected-adp-pair] {index:,}/{len(ids):,}", flush=True)
    temporary.replace(out_h5)
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "contract_version": INJECTION_ADP_PAIR_CONTRACT,
        "canonical_injection_h5": str(canonical_injection_h5),
        "raw_source_h5": str(raw_source_h5),
        "out_h5": str(out_h5),
        "n_injections": len(ids),
        "apertures": list(ADP_ONLY_APERTURES),
    }


def normalize_injection_peak_candidates(
    peaks: pd.DataFrame,
    *,
    pair_h5: Path,
    small_peaks_per_injection: int = 5,
) -> pd.DataFrame:
    """Use ADP-small top-N candidates with ADP-primary rank-1 context."""

    frame = peaks.copy()
    required = {"injection_id", "aperture", "status", *PEAK_FIELDS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError(f"injection peak table is missing columns: {missing}")
    apertures = set(frame["aperture"].dropna().astype(str))
    if not apertures.issubset(set(ADP_ONLY_APERTURES)):
        raise ValueError(
            f"canonical/non-ADP peak rows are forbidden: {sorted(apertures)}"
        )
    rank = pd.to_numeric(frame["peak_rank"], errors="coerce")
    small = frame.loc[
        frame["aperture"].astype(str).eq(ADP_ONLY_APERTURES[0])
        & frame["status"].astype(str).eq("ok")
        & rank.between(1, int(small_peaks_per_injection))
    ].copy()
    primary = frame.loc[
        frame["aperture"].astype(str).eq(ADP_ONLY_APERTURES[1])
        & frame["status"].astype(str).eq("ok")
        & rank.eq(1)
    ].copy()
    primary = primary.sort_values("injection_id", kind="stable").drop_duplicates(
        "injection_id", keep="first"
    )
    primary_fields = [column for column in PEAK_FIELDS if column in primary]
    primary = primary[["injection_id", *primary_fields]].rename(
        columns={column: f"adp_{column}" for column in primary_fields}
    )
    out = small.merge(primary, on="injection_id", how="left", validate="many_to_one")
    out["rep_aperture"] = ADP_ONLY_APERTURES[0]
    out["rep_peak_rank"] = pd.to_numeric(out["peak_rank"], errors="coerce").astype(int)
    out["sde_max"] = pd.to_numeric(out["sde"], errors="coerce")
    out["anchor_aperture"] = ADP_ONLY_APERTURES[0]
    out["anchor_period_d"] = pd.to_numeric(out["period_d"], errors="coerce")
    out["anchor_t0_bjd"] = pd.to_numeric(out["t0_bjd"], errors="coerce")
    out["anchor_duration_min"] = pd.to_numeric(out["duration_min"], errors="coerce")
    out["anchor_sde"] = out["sde_max"]
    for column in PEAK_FIELDS:
        if column in out:
            out[f"adp_sml_{column}"] = pd.to_numeric(out[column], errors="coerce")
    relation = classify_period_relation(out["adp_sml_period_d"], out["adp_period_d"])
    out["aperture_period_relation"] = relation
    out["aperture_period_rel_delta"] = (
        np.abs(out["adp_period_d"] - out["adp_sml_period_d"]) / out["adp_sml_period_d"]
    )
    out["aperture_depth_ratio_primary_over_small"] = pd.to_numeric(
        out.get("adp_depth"), errors="coerce"
    ) / pd.to_numeric(out.get("adp_sml_depth"), errors="coerce").replace(0, np.nan)
    out["aperture_disagreement_flag"] = relation.eq("unrelated")
    out["n_apertures_agree"] = np.where(relation.isin({"exact", "harmonic"}), 2, 1)
    out["apertures_agree"] = np.where(
        relation.isin({"exact", "harmonic"}),
        ",".join(ADP_ONLY_APERTURES),
        ADP_ONLY_APERTURES[0],
    )
    out["review_id"] = [
        f"s0056-a2v1-eval-{injection_id}-r{rank_value}"
        for injection_id, rank_value in zip(out["injection_id"], out["rep_peak_rank"])
    ]
    if out["review_id"].duplicated().any():
        raise ValueError("injection candidate review IDs are not unique")
    out["source_kind"] = "injected_validation_holdout"
    out["is_injected_row"] = True
    out["source_h5"] = str(Path(pair_h5))
    out["h5_group"] = "injections/" + out["injection_id"].astype(str)
    out["native_input_include"] = True
    out["native_group_path"] = [
        native_group_path(row) for row in out.to_dict("records")
    ]
    out["source_product_tag"] = "S56_A2v1_fresh_eval_v1"
    out["bls_search_branch"] = "A2v1_ADP_small_top5"
    out["adp_only_contract_version"] = ADP_ONLY_CONTRACT_VERSION
    out["input_contract_version"] = A2V1_TEACHER_INPUT_CONTRACT
    return out.sort_values(
        ["injection_id", "rep_peak_rank"], kind="stable"
    ).reset_index(drop=True)


def _pair_lc(group: Any, pair_h5: Path) -> HLSPLightCurve:
    return HLSPLightCurve(
        tic=int(group.attrs["tic"]),
        tmag=float(group.attrs.get("tessmag", np.nan)),
        sector=int(group.attrs.get("sector", 56)),
        cam=int(group.attrs.get("camera", -1)),
        ccd=int(group.attrs.get("ccd", -1)),
        ra=float("nan"),
        dec=float("nan"),
        time=np.asarray(group["time"], dtype=np.float64),
        cadenceno=np.asarray(group["cadenceno"], dtype=np.int64),
        orbitid=np.asarray(group["orbitid"], dtype=np.int32),
        quality=np.asarray(group["quality"], dtype=np.int32),
        flux={
            aperture: np.asarray(group[f"{aperture}_injected"], dtype=np.float64)
            for aperture in ADP_ONLY_APERTURES
        },
        path=Path(pair_h5),
    )


def enrich_injection_candidate_metadata(
    candidates: pd.DataFrame,
    *,
    pair_h5: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Measure the same fixed-ephemeris scalar metadata used for real rows."""

    import h5py

    rows: list[dict[str, Any]] = []
    with h5py.File(pair_h5, "r") as h5:
        for index, (injection_id, group_rows) in enumerate(
            candidates.groupby("injection_id", sort=True), start=1
        ):
            lc = _pair_lc(h5[f"injections/{injection_id}"], pair_h5)
            for record in group_rows.to_dict("records"):
                anchor = {key: record.get(key, np.nan) for key in PEAK_FIELDS}
                own = {
                    ADP_ONLY_APERTURES[0]: {
                        key: record.get(f"adp_sml_{key}", np.nan) for key in PEAK_FIELDS
                    },
                    ADP_ONLY_APERTURES[1]: {
                        key: record.get(f"adp_{key}", np.nan) for key in PEAK_FIELDS
                    },
                }
                measured = measure_two_aperture_candidate_metadata(
                    lc,
                    anchor_peak=anchor,
                    own_peaks=own,
                    apertures=ADP_ONLY_APERTURES,
                )
                measured["review_id"] = str(record["review_id"])
                measured["metadata_status"] = "ok"
                rows.append(measured)
            if index % 50 == 0:
                print(
                    f"[injection-candidates] metadata {index:,}/"
                    f"{candidates['injection_id'].nunique():,}",
                    flush=True,
                )
    metrics = pd.DataFrame(rows)
    overlap = [
        column for column in metrics if column != "review_id" and column in candidates
    ]
    merged = candidates.drop(columns=overlap).merge(
        metrics, on="review_id", how="left", validate="one_to_one"
    )
    passed = merged["metadata_status"].fillna("missing").eq("ok")
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "pair_h5": str(pair_h5),
        "n_candidates": int(len(merged)),
        "n_injections": int(merged["injection_id"].nunique()),
        "n_metadata_ok": int(passed.sum()),
        "passed": bool(passed.all()),
    }
    return merged, summary


def aggregate_teacher_injection_recovery(
    manifest: pd.DataFrame,
    scored_candidates: pd.DataFrame,
    *,
    preserve_threshold: float = 0.5,
    max_peak_rank: int = 5,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build one end-to-end BLS and Teacher-v1 outcome per injection."""

    scored = scored_candidates.copy()
    required = {
        "injection_id",
        "rep_peak_rank",
        "is_injected_signal_peak",
        "p_preserve",
        "predicted_morphology",
    }
    missing = sorted(required - set(scored.columns))
    if missing:
        raise KeyError(f"teacher score table is missing columns: {missing}")
    scored = scored.loc[
        pd.to_numeric(scored["rep_peak_rank"], errors="coerce").le(int(max_peak_rank))
    ].copy()
    scored["truth_matched_candidate"] = _bool_series(scored["is_injected_signal_peak"])
    scored["preserve_pass"] = pd.to_numeric(scored["p_preserve"], errors="coerce").ge(
        float(preserve_threshold)
    )
    scored["planet_pass"] = scored["predicted_morphology"].astype(str).eq("planet_like")
    rows: list[dict[str, Any]] = []
    for injection_id, group in scored.groupby("injection_id", sort=False):
        matched = group[group["truth_matched_candidate"]]
        rows.append(
            {
                "injection_id": str(injection_id),
                "n_scored_candidates": int(len(group)),
                "n_truth_matched_candidates": int(len(matched)),
                "bls_adp_only_recovered": bool(len(matched)),
                "teacher_v1_preserve_recovered": bool(
                    len(matched) and matched["preserve_pass"].any()
                ),
                "teacher_v1_planet_recovered": bool(
                    len(matched) and matched["planet_pass"].any()
                ),
                "matched_p_preserve_max": (
                    float(pd.to_numeric(matched["p_preserve"], errors="coerce").max())
                    if len(matched)
                    else np.nan
                ),
                "matched_p_planet_max": (
                    float(
                        pd.to_numeric(
                            matched.get("p_planet_like"), errors="coerce"
                        ).max()
                    )
                    if len(matched) and "p_planet_like" in matched
                    else np.nan
                ),
            }
        )
    outcomes = pd.DataFrame(rows)
    result = manifest.merge(
        outcomes, on="injection_id", how="left", validate="one_to_one"
    )
    for column in (
        "bls_adp_only_recovered",
        "teacher_v1_preserve_recovered",
        "teacher_v1_planet_recovered",
    ):
        result[column] = result[column].fillna(False).astype(bool)
    result["n_scored_candidates"] = result["n_scored_candidates"].fillna(0).astype(int)
    result["n_truth_matched_candidates"] = (
        result["n_truth_matched_candidates"].fillna(0).astype(int)
    )
    bls = result["bls_adp_only_recovered"]
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy": TEACHER_RECOVERY_POLICY,
        "preserve_threshold": float(preserve_threshold),
        "max_peak_rank": int(max_peak_rank),
        "n_injections": int(len(result)),
        "n_bls_recovered": int(bls.sum()),
        "bls_recovery_fraction": float(bls.mean()),
        "n_teacher_preserve_recovered": int(
            result["teacher_v1_preserve_recovered"].sum()
        ),
        "teacher_preserve_end_to_end_fraction": float(
            result["teacher_v1_preserve_recovered"].mean()
        ),
        "teacher_preserve_conditional_on_bls": (
            float(result.loc[bls, "teacher_v1_preserve_recovered"].mean())
            if bls.any()
            else np.nan
        ),
        "n_teacher_planet_recovered": int(result["teacher_v1_planet_recovered"].sum()),
        "teacher_planet_end_to_end_fraction": float(
            result["teacher_v1_planet_recovered"].mean()
        ),
    }
    return result, summary


__all__ = [
    "INJECTION_ADP_PAIR_CONTRACT",
    "TEACHER_RECOVERY_POLICY",
    "aggregate_teacher_injection_recovery",
    "align_raw_source_by_time",
    "build_teacher_injection_holdout",
    "enrich_injection_candidate_metadata",
    "normalize_injection_peak_candidates",
    "rebuild_injected_adp_pair",
    "subset_raw_source_h5",
]
