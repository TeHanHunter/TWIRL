"""Real-only S56 label adjudication queue and post-review audit helpers."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from twirl.vetting.label_io import (
    ADJUDICATION_LABEL_COLUMNS,
    latest_label_records,
    normalize_review_queue,
    validate_label_records,
)
from twirl.vetting.two_aperture import TWO_APERTURE_VET_SHEET_VERSION


SIGNAL_LABELS: frozenset[str] = frozenset(
    {
        "planet_like",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "wide_transit_like",
    }
)
CONTROL_LABELS: tuple[str, ...] = ("instrumental_or_systematic", "uncertain")
REPEAT_QUOTAS: Mapping[str, int] = {
    "eclipsing_binary_or_pceb": 4,
    "stellar_variability": 4,
    "planet_like": 3,
    "wide_transit_like": 3,
    "instrumental_or_systematic": 3,
    "uncertain": 3,
}
ORIGINAL_1K = "s56_recovery50_teacher_queue"
NEXT4K_SLICE = "s56_recovery50_teacher_queue_next4k"
ADP_APERTURES = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
ADJUDICATION_VET_SHEET_VERSION = TWO_APERTURE_VET_SHEET_VERSION

PUBLIC_QUEUE_COLUMNS: tuple[str, ...] = (
    "row_id",
    "review_id",
    "tic",
    "sector",
    "cam",
    "ccd",
    "tmag",
    "source_kind",
    "source_bucket",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde_max",
    "aperture_period_rel_delta",
    "aperture_disagreement_flag",
    "tensor_apertures",
    "adp_only_contract_version",
    "twirl_vet_sheet_name",
    "twirl_vet_sheet_pdf_name",
    "candidate_key",
)

HIDDEN_PUBLIC_SUBSTRINGS: tuple[str, ...] = (
    "human_label",
    "prior_label",
    "raw_label",
    "source_cohort",
    "repeat_group",
    "is_repeat",
    "eb_p_",
    "priority_rank",
    "selection_bucket",
    "training_split",
)

ADJUDICATION_LEAKAGE_COLUMNS: frozenset[str] = frozenset(
    {
        "human_label_adjudicated",
        "adjudicated_period_factor",
        "adjudicated_period_status",
        "adjudicated_period_d",
        "adjudicated_notes",
        "adjudicated_labeler",
        "adjudicated_updated_utc",
        "morphology_target_v1",
        "morphology_include_v1",
        "preserve_target_v1",
        "preserve_include_v1",
        "harmonic_target_v1",
        "harmonic_include_v1",
        "note_period_factor",
        "effective_period_factor",
        "effective_period_d",
        "period_factor_source",
        "period_task",
        "variable_period_factor",
        "broad_preserve_only",
        "model_target_policy_version",
        "raw_human_label",
        "raw_human_notes",
        "raw_human_updated_utc",
        "historical_signal_label",
        "source_cohort",
        "labeling_era",
        "repeat_group",
        "repeat_occurrence",
        "is_repeat",
    }
)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value)


def _bool_series(frame: pd.DataFrame, column: str, *, default: bool = False) -> pd.Series:
    if column not in frame:
        return pd.Series(default, index=frame.index, dtype=bool)
    values = frame[column]
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(default).astype(bool)
    return values.fillna("").astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y"})


def _real_rows(frame: pd.DataFrame) -> pd.Series:
    if "is_injected_row" in frame:
        return ~_bool_series(frame, "is_injected_row")
    return frame.get("source_kind", pd.Series("", index=frame.index)).fillna("").astype(str).eq("real_candidate")


def _hash_token(value: str, *, seed: int) -> str:
    return hashlib.sha256(f"{seed}|{value}".encode("utf-8")).hexdigest()


def join_browser_labels(queue_csv: Path, labels_csv: Path) -> pd.DataFrame:
    """Join the latest browser labels through the normalized display identity."""

    queue = normalize_review_queue(read_table(queue_csv))
    for column in ("label", "label_source", "labeler", "notes", "updated_utc"):
        if column in queue:
            queue = queue.rename(columns={column: f"queue_{column}"})
    labels = latest_label_records(Path(labels_csv))
    validate_label_records(queue, labels)
    labels = labels.rename(
        columns={
            "candidate_key": "browser_candidate_key",
            "label": "browser_label",
            "label_source": "browser_label_source",
            "labeler": "browser_labeler",
            "notes": "browser_notes",
            "updated_utc": "browser_updated_utc",
            "period_factor": "browser_period_factor",
            "period_status": "browser_period_status",
        }
    )
    keep = [
        "row_id",
        "browser_candidate_key",
        "browser_label",
        "browser_label_source",
        "browser_labeler",
        "browser_notes",
        "browser_updated_utc",
        "browser_period_factor",
        "browser_period_status",
    ]
    return queue.merge(labels.loc[:, keep], on="row_id", how="left", validate="one_to_one")


def _source_record(
    row: pd.Series,
    *,
    source_uid: str,
    source_cohort: str,
    source_review_id: str,
    origin_queue: str,
    labeling_era: str,
    first_label: str,
    raw_label: str,
    raw_label_source: str,
    raw_labeler: str,
    raw_notes: str,
    raw_updated_utc: str,
    historical_signal_label: str,
    ephemeris_source: str,
) -> dict[str, Any]:
    sector = row.get("sector", 56)
    if not _text(sector):
        sector = 56
    return {
        "source_uid": source_uid,
        "source_cohort": source_cohort,
        "source_review_id": source_review_id,
        "origin_queue": origin_queue,
        "labeling_era": labeling_era,
        "first_human_label": first_label,
        "raw_human_label": raw_label,
        "raw_human_label_source": raw_label_source,
        "raw_human_labeler": raw_labeler,
        "raw_human_notes": raw_notes,
        "raw_human_updated_utc": raw_updated_utc,
        "historical_signal_label": historical_signal_label,
        "ephemeris_source": ephemeris_source,
        "tic": row.get("tic", np.nan),
        "sector": sector,
        "cam": row.get("cam", np.nan),
        "ccd": row.get("ccd", np.nan),
        "tmag": row.get("tmag", np.nan),
        "period_d": row.get("period_d", np.nan),
        "t0_bjd": row.get("t0_bjd", np.nan),
        "duration_min": row.get("duration_min", np.nan),
        "depth": row.get("depth", np.nan),
        "depth_snr": row.get("depth_snr", np.nan),
        "sde_max": row.get("sde_max", row.get("sde", np.nan)),
        "aperture_period_rel_delta": row.get("aperture_period_rel_delta", np.nan),
        "aperture_disagreement_flag": row.get("aperture_disagreement_flag", np.nan),
    }


def _teacher_signal_sources(base: pd.DataFrame, active: pd.DataFrame) -> pd.DataFrame:
    for name, frame in (("base", base), ("active", active)):
        if frame["review_id"].duplicated().any():
            raise ValueError(f"{name} 2k table has duplicate review_id values")
    base_real = base.loc[_real_rows(base)].copy()
    active_real = active.loc[_real_rows(active)].copy()
    historical_ids = set(base_real.loc[base_real["human_label"].isin(SIGNAL_LABELS), "review_id"].astype(str))
    historical_ids.update(active_real.loc[active_real["human_label"].isin(SIGNAL_LABELS), "review_id"].astype(str))
    active_by_id = active_real.set_index(active_real["review_id"].astype(str))
    missing = sorted(historical_ids - set(active_by_id.index))
    if missing:
        raise ValueError(f"historical 2k signal rows missing from active ADP table: {missing[:10]}")
    base_by_id = base_real.set_index(base_real["review_id"].astype(str))
    records: list[dict[str, Any]] = []
    for review_id in sorted(historical_ids):
        row = active_by_id.loc[review_id]
        first = base_by_id.loc[review_id] if review_id in base_by_id.index else row
        first_label = _text(first.get("human_label"))
        current_label = _text(row.get("human_label"))
        historical_label = first_label if first_label in SIGNAL_LABELS else current_label
        origin = _text(row.get("origin_queue"))
        records.append(
            _source_record(
                row,
                source_uid=f"teacher2k:{review_id}",
                source_cohort="teacher2k_historical_signal",
                source_review_id=review_id,
                origin_queue=origin,
                labeling_era="early" if origin == ORIGINAL_1K else "recent",
                first_label=first_label,
                raw_label=current_label,
                raw_label_source=_text(row.get("human_label_source")),
                raw_labeler=_text(row.get("human_labeler")),
                raw_notes=_text(row.get("human_notes")),
                raw_updated_utc=_text(row.get("human_updated_utc")),
                historical_signal_label=historical_label,
                ephemeris_source="current_adp_sml_teacher_anchor",
            )
        )
    return pd.DataFrame(records)


def _browser_signal_sources(
    joined: pd.DataFrame,
    *,
    source_cohort: str,
    source_prefix: str,
    origin_queue: str,
    ephemeris_source: str,
) -> pd.DataFrame:
    selected = joined.loc[joined["browser_label"].fillna("").astype(str).isin(SIGNAL_LABELS)].copy()
    records: list[dict[str, Any]] = []
    for _, row in selected.iterrows():
        source_review_id = _text(row.get("review_id")) or f"tic:{int(float(row['tic']))}"
        label = _text(row.get("browser_label"))
        records.append(
            _source_record(
                row,
                source_uid=f"{source_prefix}:{source_review_id}",
                source_cohort=source_cohort,
                source_review_id=source_review_id,
                origin_queue=origin_queue,
                labeling_era="recent",
                first_label=label,
                raw_label=label,
                raw_label_source=_text(row.get("browser_label_source")),
                raw_labeler=_text(row.get("browser_labeler")),
                raw_notes=_text(row.get("browser_notes")),
                raw_updated_utc=_text(row.get("browser_updated_utc")),
                historical_signal_label=label,
                ephemeris_source=ephemeris_source,
            )
        )
    return pd.DataFrame(records)


def reanchor_deprecated_sources(deprecated: pd.DataFrame, adp_bls: pd.DataFrame) -> pd.DataFrame:
    """Replace every deprecated review ephemeris with current ADP-small BLS top-1."""

    peaks = adp_bls.copy()
    aperture = peaks.get("aperture", pd.Series("", index=peaks.index)).fillna("").astype(str)
    rank = pd.to_numeric(peaks.get("peak_rank", pd.Series(np.nan, index=peaks.index)), errors="coerce")
    status = peaks.get("status", pd.Series("ok", index=peaks.index)).fillna("").astype(str).str.lower()
    peaks = peaks.loc[aperture.eq(ADP_APERTURES[0]) & rank.eq(1) & status.eq("ok")].copy()
    peaks["tic"] = pd.to_numeric(peaks["tic"], errors="raise").astype(np.int64)
    if peaks["tic"].duplicated().any():
        duplicated = peaks.loc[peaks["tic"].duplicated(keep=False), "tic"].astype(int).unique().tolist()
        raise ValueError(f"current ADP-small top-1 table has duplicate TIC rows: {duplicated[:10]}")
    peak_by_tic = peaks.set_index("tic")
    out = deprecated.copy()
    missing: list[int] = []
    for index, row in out.iterrows():
        tic = int(float(row["tic"]))
        if tic not in peak_by_tic.index:
            missing.append(tic)
            continue
        peak = peak_by_tic.loc[tic]
        for column in ("period_d", "t0_bjd", "duration_min", "depth", "depth_snr"):
            out.loc[index, column] = peak.get(column, np.nan)
        out.loc[index, "sde_max"] = peak.get("sde", np.nan)
        out.loc[index, "ephemeris_source"] = "current_adp_sml_bls_top1"
    if missing:
        raise ValueError(f"deprecated signal rows missing current ADP-small BLS top-1: {sorted(missing)[:10]}")
    for column in ("period_d", "t0_bjd", "duration_min"):
        values = pd.to_numeric(out[column], errors="coerce")
        if (~np.isfinite(values)).any() or (values <= 0).any():
            raise ValueError(f"deprecated re-anchoring produced invalid {column}")
    return out


def _quartile(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series("missing", index=values.index, dtype=object)
    finite = np.isfinite(numeric)
    if finite.any():
        ranks = numeric.loc[finite].rank(method="first") - 1.0
        bins = np.floor(4.0 * ranks / max(int(finite.sum()), 1)).clip(0, 3).astype(int)
        out.loc[finite] = bins.map(lambda value: f"q{value + 1}")
    return out


def _stratified_control_sample(frame: pd.DataFrame, *, n: int, seed: int) -> pd.DataFrame:
    if len(frame) < n:
        raise ValueError(f"control stratum has only {len(frame)} rows; need {n}")
    work = frame.copy()
    work["control_period_quartile"] = _quartile(work["period_d"])
    work["control_tmag_quartile"] = _quartile(work["tmag"])
    work["control_sde_quartile"] = _quartile(work["sde_max"])
    disagreement = _bool_series(work, "aperture_disagreement_flag")
    missing_disagreement = work.get(
        "aperture_disagreement_flag", pd.Series(np.nan, index=work.index)
    ).isna()
    rel_delta = pd.to_numeric(
        work.get("aperture_period_rel_delta", pd.Series(np.nan, index=work.index)),
        errors="coerce",
    )
    agreement = pd.Series(np.where(disagreement, "disagree", "agree"), index=work.index, dtype=object)
    agreement.loc[missing_disagreement & ~np.isfinite(rel_delta)] = "unknown"
    agreement.loc[missing_disagreement & np.isfinite(rel_delta)] = np.where(
        rel_delta.loc[missing_disagreement & np.isfinite(rel_delta)] <= 0.02,
        "agree",
        "disagree",
    )
    work["control_aperture_agreement"] = agreement
    work["control_stratum"] = (
        work["control_period_quartile"].astype(str)
        + "|"
        + work["control_tmag_quartile"].astype(str)
        + "|"
        + work["control_sde_quartile"].astype(str)
        + "|"
        + work["control_aperture_agreement"].astype(str)
    )
    work["_row_hash"] = work["review_id"].astype(str).map(lambda value: _hash_token(value, seed=seed))
    work["_stratum_hash"] = work["control_stratum"].map(lambda value: _hash_token(value, seed=seed + 1))
    work = work.sort_values(["control_stratum", "_row_hash"], kind="stable")
    work["_round"] = work.groupby("control_stratum", sort=False).cumcount()
    work = work.sort_values(["_round", "_stratum_hash", "_row_hash"], kind="stable").head(n)
    return work.drop(columns=["_row_hash", "_stratum_hash", "_round"])


def _control_sources(
    active: pd.DataFrame,
    *,
    historical_review_ids: set[str],
    excluded_tics: set[int],
    per_origin_per_label: int,
    seed: int,
) -> pd.DataFrame:
    real = active.loc[_real_rows(active)].copy()
    real["review_id"] = real["review_id"].astype(str)
    real["tic"] = pd.to_numeric(real["tic"], errors="raise").astype(np.int64)
    real = real.loc[~real["review_id"].isin(historical_review_ids) & ~real["tic"].isin(excluded_tics)].copy()
    records: list[dict[str, Any]] = []
    for origin_index, origin in enumerate((ORIGINAL_1K, NEXT4K_SLICE)):
        for label_index, label in enumerate(CONTROL_LABELS):
            pool = real.loc[
                real["origin_queue"].fillna("").astype(str).eq(origin)
                & real["human_label"].fillna("").astype(str).eq(label)
            ].copy()
            chosen = _stratified_control_sample(
                pool,
                n=per_origin_per_label,
                seed=seed + 100 * origin_index + 10 * label_index,
            )
            for _, row in chosen.iterrows():
                review_id = _text(row.get("review_id"))
                record = _source_record(
                    row,
                    source_uid=f"teacher2k:{review_id}",
                    source_cohort="teacher2k_control",
                    source_review_id=review_id,
                    origin_queue=origin,
                    labeling_era="early" if origin == ORIGINAL_1K else "recent",
                    first_label=_text(row.get("human_label")),
                    raw_label=_text(row.get("human_label")),
                    raw_label_source=_text(row.get("human_label_source")),
                    raw_labeler=_text(row.get("human_labeler")),
                    raw_notes=_text(row.get("human_notes")),
                    raw_updated_utc=_text(row.get("human_updated_utc")),
                    historical_signal_label="",
                    ephemeris_source="current_adp_sml_teacher_anchor",
                )
                for column in (
                    "control_period_quartile",
                    "control_tmag_quartile",
                    "control_sde_quartile",
                    "control_aperture_agreement",
                    "control_stratum",
                ):
                    record[column] = row.get(column, "")
                records.append(record)
    return pd.DataFrame(records)


def _select_repeat_sources(
    unique_sources: pd.DataFrame,
    *,
    repeat_quotas: Mapping[str, int],
    seed: int,
) -> pd.DataFrame:
    selected: list[pd.DataFrame] = []
    for label, count in repeat_quotas.items():
        pool = unique_sources.loc[unique_sources["raw_human_label"].fillna("").astype(str).eq(label)].copy()
        if len(pool) < count:
            raise ValueError(f"repeat quota {label}={count} exceeds available rows ({len(pool)})")
        pool["_hash"] = pool["source_uid"].astype(str).map(lambda value: _hash_token(value, seed=seed))
        chosen = pool.sort_values("_hash", kind="stable").head(count).drop(columns="_hash")
        chosen["repeat_reference_label"] = label
        selected.append(chosen)
    out = pd.concat(selected, ignore_index=True, sort=False)
    out["_hash"] = out["source_uid"].astype(str).map(lambda value: _hash_token(value, seed=seed + 7))
    out = out.sort_values("_hash", kind="stable").drop(columns="_hash").reset_index(drop=True)
    out["repeat_group"] = [f"repeat:{index:02d}" for index in range(len(out))]
    return out


def place_blinded_repeats(
    unique_sources: pd.DataFrame,
    repeat_sources: pd.DataFrame,
    *,
    seed: int,
    min_separation: int = 100,
) -> pd.DataFrame:
    """Place repeats in opposite queue ends while shuffling every visible row."""

    n_unique = len(unique_sources)
    n_repeat = len(repeat_sources)
    n_total = n_unique + n_repeat
    if n_repeat != 20 or n_total < 300:
        raise ValueError("the fixed blinded-repeat layout requires 20 repeats and at least 300 total rows")
    repeated_uids = set(repeat_sources["source_uid"].astype(str))
    if len(repeated_uids) != n_repeat:
        raise ValueError("repeat sources must be unique")
    unique_by_uid = unique_sources.set_index("source_uid")
    rng = np.random.default_rng(seed)
    low_limit = 80
    high_start = n_total - 80
    low_positions = rng.choice(np.arange(0, low_limit), size=n_repeat, replace=False)
    high_positions = rng.choice(np.arange(high_start, n_total), size=n_repeat, replace=False)
    slots: dict[int, dict[str, Any]] = {}
    half = n_repeat // 2
    for repeat_index, repeat in repeat_sources.reset_index(drop=True).iterrows():
        uid = str(repeat["source_uid"])
        source = unique_by_uid.loc[uid].to_dict()
        source["source_uid"] = uid
        source["repeat_group"] = str(repeat["repeat_group"])
        source["repeat_reference_label"] = str(repeat["repeat_reference_label"])
        duplicate = dict(source)
        source["repeat_occurrence"] = "original"
        source["is_repeat"] = False
        duplicate["repeat_occurrence"] = "blind_repeat"
        duplicate["is_repeat"] = True
        if repeat_index < half:
            original_position = int(low_positions[repeat_index])
            repeat_position = int(high_positions[repeat_index])
        else:
            original_position = int(high_positions[repeat_index])
            repeat_position = int(low_positions[repeat_index])
        slots[original_position] = source
        slots[repeat_position] = duplicate

    remaining = unique_sources.loc[~unique_sources["source_uid"].astype(str).isin(repeated_uids)].copy()
    remaining["repeat_group"] = ""
    remaining["repeat_reference_label"] = ""
    remaining["repeat_occurrence"] = "single"
    remaining["is_repeat"] = False
    remaining_records = remaining.to_dict("records")
    rng.shuffle(remaining_records)
    open_positions = [position for position in range(n_total) if position not in slots]
    rng.shuffle(open_positions)
    if len(open_positions) != len(remaining_records):
        raise AssertionError("repeat placement did not leave one slot per remaining source")
    for position, record in zip(open_positions, remaining_records):
        slots[position] = record

    placed = pd.DataFrame([slots[position] for position in range(n_total)])
    placed.insert(0, "display_position", np.arange(n_total, dtype=int))
    distances = placed.loc[placed["repeat_group"].ne("")].groupby("repeat_group")["display_position"].agg(
        lambda values: int(values.max() - values.min())
    )
    if len(distances) != n_repeat or (distances < min_separation).any():
        raise ValueError(f"blind-repeat separation failed: {distances.to_dict()}")
    return placed


def _public_and_manifest(placed: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    public = pd.DataFrame(
        {
            "review_id": placed["display_position"].map(lambda value: f"adjudication:{int(value):04d}"),
            "tic": placed["tic"],
            "sector": placed["sector"],
            "cam": placed["cam"],
            "ccd": placed["ccd"],
            "tmag": placed["tmag"],
            "source_kind": "real_candidate",
            "source_bucket": "real_adjudication_blind",
            "period_d": placed["period_d"],
            "t0_bjd": placed["t0_bjd"],
            "duration_min": placed["duration_min"],
            "depth": placed["depth"],
            "depth_snr": placed["depth_snr"],
            "sde_max": placed["sde_max"],
            "aperture_period_rel_delta": placed["aperture_period_rel_delta"],
            "aperture_disagreement_flag": placed["aperture_disagreement_flag"],
            "tensor_apertures": ",".join(ADP_APERTURES),
            "adp_only_contract_version": "s56_adp_pair_v1",
            "twirl_vet_sheet_name": placed["display_position"].map(
                lambda value: f"adjudication_{int(value):04d}_twirl_twoap_current_adp.png"
            ),
            "twirl_vet_sheet_pdf_name": "",
        }
    )
    public = normalize_review_queue(public)
    public = public.loc[:, PUBLIC_QUEUE_COLUMNS]
    hidden = [column for column in public.columns if any(token in column.lower() for token in HIDDEN_PUBLIC_SUBSTRINGS)]
    if hidden:
        raise ValueError(f"public adjudication queue exposes hidden columns: {hidden}")
    manifest = placed.reset_index(drop=True).copy()
    manifest.insert(0, "row_id", public["row_id"].to_numpy())
    manifest.insert(1, "review_id", public["review_id"].to_numpy())
    manifest.insert(2, "candidate_key", public["candidate_key"].to_numpy())
    manifest["review_period_d"] = public["period_d"].to_numpy()
    manifest["review_t0_bjd"] = public["t0_bjd"].to_numpy()
    manifest["review_duration_min"] = public["duration_min"].to_numpy()
    return public, manifest


def verify_queue_contract(
    public: pd.DataFrame,
    manifest: pd.DataFrame,
    *,
    expected_total: int = 343,
    expected_unique: int = 323,
    expected_signal: int = 223,
    expected_controls_per_label: int = 50,
) -> dict[str, Any]:
    failures: list[str] = []
    if len(public) != expected_total:
        failures.append(f"public rows={len(public)} expected={expected_total}")
    if len(manifest) != expected_total:
        failures.append(f"manifest rows={len(manifest)} expected={expected_total}")
    if manifest["source_uid"].nunique() != expected_unique:
        failures.append(f"unique sources={manifest['source_uid'].nunique()} expected={expected_unique}")
    repeat_rows = _bool_series(manifest, "is_repeat")
    if int(repeat_rows.sum()) != expected_total - expected_unique:
        failures.append(f"repeat rows={int(repeat_rows.sum())} expected={expected_total - expected_unique}")
    unique = manifest.loc[~repeat_rows].copy()
    signal_count = int(unique["historical_signal_label"].fillna("").astype(str).isin(SIGNAL_LABELS).sum())
    if signal_count != expected_signal:
        failures.append(f"historical signal sources={signal_count} expected={expected_signal}")
    for label in CONTROL_LABELS:
        count = int(
            (
                unique["source_cohort"].eq("teacher2k_control")
                & unique["raw_human_label"].fillna("").astype(str).eq(label)
            ).sum()
        )
        if count != expected_controls_per_label:
            failures.append(f"control {label}={count} expected={expected_controls_per_label}")
    repeat_counts = (
        manifest.loc[repeat_rows, "repeat_reference_label"].fillna("").astype(str).value_counts().to_dict()
    )
    if repeat_counts != dict(REPEAT_QUOTAS):
        failures.append(f"repeat quotas={repeat_counts} expected={dict(REPEAT_QUOTAS)}")
    distances = manifest.loc[manifest["repeat_group"].fillna("").astype(str).ne("")].groupby("repeat_group")[
        "display_position"
    ].agg(lambda values: int(values.max() - values.min()))
    if len(distances) != sum(REPEAT_QUOTAS.values()) or (distances < 100).any():
        failures.append(f"repeat separation failure: {distances.to_dict()}")
    normalized = normalize_review_queue(public)
    if not normalized["candidate_key"].equals(public["candidate_key"]):
        failures.append("public candidate keys do not reproduce after normalization")
    invalid_ephemeris = (
        ~np.isfinite(pd.to_numeric(public["period_d"], errors="coerce"))
        | ~np.isfinite(pd.to_numeric(public["t0_bjd"], errors="coerce"))
        | ~np.isfinite(pd.to_numeric(public["duration_min"], errors="coerce"))
    )
    if invalid_ephemeris.any():
        failures.append(f"invalid review ephemerides={int(invalid_ephemeris.sum())}")
    return {
        "passed": not failures,
        "failures": failures,
        "n_rows": int(len(public)),
        "n_unique_sources": int(manifest["source_uid"].nunique()),
        "n_repeats": int(repeat_rows.sum()),
        "minimum_repeat_separation": int(distances.min()) if len(distances) else 0,
    }


def _snapshot_file(source: Path, snapshot_dir: Path, *, prefix: str) -> tuple[Path, str]:
    digest = hashlib.sha256(Path(source).read_bytes()).hexdigest()
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    destination = snapshot_dir / f"{prefix}_{digest[:12]}.csv"
    if not destination.exists():
        shutil.copy2(source, destination)
    elif hashlib.sha256(destination.read_bytes()).hexdigest() != digest:
        raise ValueError(f"existing source snapshot has changed: {destination}")
    return destination, digest


def build_real_adjudication_queue(
    *,
    base_training_table: Path,
    active_training_table: Path,
    eb_pilot_queue: Path,
    eb_pilot_labels: Path,
    deprecated_eb_queue: Path,
    deprecated_eb_labels: Path,
    current_adp_bls: Path,
    out_dir: Path,
    seed: int = 56,
    expected_teacher_signals: int = 177,
    expected_pilot_signals: int = 34,
    expected_deprecated_signals: int = 12,
) -> dict[str, Any]:
    """Build the fixed 343-row, real-only S56 adjudication queue."""

    base = read_table(base_training_table)
    active = read_table(active_training_table)
    teacher = _teacher_signal_sources(base, active)
    pilot_joined = join_browser_labels(eb_pilot_queue, eb_pilot_labels)
    pilot = _browser_signal_sources(
        pilot_joined,
        source_cohort="eb_pilot_adp_signal",
        source_prefix="eb_pilot_adp",
        origin_queue="s56_eb_miner_adp_only_pilot100",
        ephemeris_source="current_adp_sml_pilot_queue",
    )
    deprecated_joined = join_browser_labels(deprecated_eb_queue, deprecated_eb_labels)
    deprecated = _browser_signal_sources(
        deprecated_joined,
        source_cohort="deprecated_eb_signal",
        source_prefix="eb_deprecated",
        origin_queue="s56_eb_miner_deprecated_queue",
        ephemeris_source="deprecated_queue_pending_rebuild",
    )
    deprecated = reanchor_deprecated_sources(deprecated, read_table(current_adp_bls))
    actual_signal_counts = (len(teacher), len(pilot), len(deprecated))
    expected_signal_counts = (
        expected_teacher_signals,
        expected_pilot_signals,
        expected_deprecated_signals,
    )
    if actual_signal_counts != expected_signal_counts:
        raise ValueError(f"signal source counts={actual_signal_counts} expected={expected_signal_counts}")
    signals = pd.concat([teacher, pilot, deprecated], ignore_index=True, sort=False)
    if signals["source_uid"].duplicated().any():
        raise ValueError("signal source UIDs overlap across input cohorts")
    signal_tics = pd.to_numeric(signals["tic"], errors="raise").astype(np.int64)
    if signal_tics.duplicated().any():
        raise ValueError("signal TICs overlap across input cohorts")
    controls = _control_sources(
        active,
        historical_review_ids=set(teacher["source_review_id"].astype(str)),
        excluded_tics=set(signal_tics.astype(int)),
        per_origin_per_label=25,
        seed=seed,
    )
    unique_sources = pd.concat([signals, controls], ignore_index=True, sort=False)
    if unique_sources["source_uid"].duplicated().any():
        raise ValueError("adjudication unique-source table contains duplicate source UIDs")
    if len(unique_sources) != 323:
        raise ValueError(f"unique adjudication source count={len(unique_sources)} expected=323")
    repeats = _select_repeat_sources(unique_sources, repeat_quotas=REPEAT_QUOTAS, seed=seed)
    placed = place_blinded_repeats(unique_sources, repeats, seed=seed)
    public, manifest = _public_and_manifest(placed)
    verification = verify_queue_contract(public, manifest)
    if not verification["passed"]:
        raise ValueError(f"adjudication queue verification failed: {verification['failures']}")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue_path = out_dir / "review_queue_real343.csv"
    manifest_path = out_dir / "adjudication_manifest_private.csv"
    source_path = out_dir / "adjudication_unique_sources_private.csv"
    labels_path = out_dir / "human_labels_adjudicated.csv"
    public.to_csv(queue_path, index=False)
    manifest.to_csv(manifest_path, index=False)
    unique_sources.to_csv(source_path, index=False)
    if not labels_path.exists():
        pd.DataFrame(columns=ADJUDICATION_LABEL_COLUMNS).to_csv(labels_path, index=False)
    snapshot_path, snapshot_sha256 = _snapshot_file(
        Path(eb_pilot_labels), out_dir / "source_snapshots", prefix="completed_eb_pilot_labels"
    )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "seed": int(seed),
        "n_display_rows": int(len(public)),
        "n_unique_sources": int(unique_sources["source_uid"].nunique()),
        "n_intentional_repeats": int(_bool_series(manifest, "is_repeat").sum()),
        "signal_source_counts": {
            "teacher2k": int(len(teacher)),
            "eb_pilot_adp": int(len(pilot)),
            "deprecated_eb_reanchored": int(len(deprecated)),
        },
        "control_counts": (
            controls.groupby(["origin_queue", "raw_human_label"]).size().sort_index().astype(int).to_dict()
        ),
        "repeat_quotas": dict(REPEAT_QUOTAS),
        "apertures": list(ADP_APERTURES),
        "vet_sheet_version": ADJUDICATION_VET_SHEET_VERSION,
        "deprecated_ephemeris_policy": "current DET_FLUX_ADP_SML BLS peak_rank=1 only",
        "verification": verification,
        "source_snapshot": str(snapshot_path),
        "source_snapshot_sha256": snapshot_sha256,
        "outputs": {
            "queue": str(queue_path),
            "private_manifest": str(manifest_path),
            "unique_sources": str(source_path),
            "labels": str(labels_path),
        },
    }
    summary["control_counts"] = {
        f"{origin}|{label}": int(value)
        for (origin, label), value in summary["control_counts"].items()
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n"
    )
    return summary


__all__ = [
    "ADJUDICATION_LEAKAGE_COLUMNS",
    "ADJUDICATION_VET_SHEET_VERSION",
    "ADP_APERTURES",
    "CONTROL_LABELS",
    "REPEAT_QUOTAS",
    "SIGNAL_LABELS",
    "build_real_adjudication_queue",
    "join_browser_labels",
    "place_blinded_repeats",
    "read_table",
    "reanchor_deprecated_sources",
    "verify_queue_contract",
]
