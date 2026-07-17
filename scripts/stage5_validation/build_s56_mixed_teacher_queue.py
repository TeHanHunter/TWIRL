#!/usr/bin/env python3
"""Build the S56 mixed teacher pool and first-pass human review queue.

This queue is source-blinded for browser review, but keeps truth/provenance
columns in the CSV for later audit and teacher-table joins.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
for path in (str(SCRIPT_DIR), str(SRC_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)

from build_s56_pretriage_review_queue import (  # noqa: E402
    KNOWN_APERTURES,
    LABEL_COLUMNS,
    _attach_real_review_columns,
    _clean,
    _finalize_queue,
    _read_table,
    annotate_star_parameters,
    load_wd_star_catalog,
    render_leo_reports,
    star_source_counts,
)


DEFAULT_REAL_CANDIDATES = (
    REPO_ROOT / "data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.csv"
)
DEFAULT_INJECTED_CANDIDATES = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "small_pair_200k/review_queue.csv"
)
DEFAULT_INJECTION_H5 = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
    / "pdo_20k_predetrend_batman_periodradius_bright_balanced/injected_lightcurves.h5"
)
DEFAULT_HLSP_ROOT = REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_adp015q_compare"
DEFAULT_STAR_CATALOG = (
    REPO_ROOT / "data_local/catalogs/twirl_master_catalog/twirl_wd_master_catalog_v0_ticmatched.fits"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_mixed_teacher_queue_pdo"

WD_1856_TIC = 267574918
CADENCE_D = 200.0 / 86400.0
LABEL_HEADER = (
    "row_id",
    "candidate_key",
    "tic",
    "sector",
    "label",
    "label_source",
    "labeler",
    "notes",
    "updated_utc",
)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n")


def _value_counts(df: pd.DataFrame, column: str) -> dict[str, int]:
    if column not in df:
        return {}
    return {
        str(key): int(value)
        for key, value in df[column].fillna("").astype(str).value_counts().sort_index().items()
    }


def _finite_numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").replace([np.inf, -np.inf], np.nan)


def _coerce_numeric(frame: pd.DataFrame, columns: tuple[str, ...]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        if column in out:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def _with_cadence_alias_columns(
    frame: pd.DataFrame,
    *,
    period_column: str,
    prefix: str = "",
    tolerance: float,
) -> pd.DataFrame:
    out = frame.copy()
    period = _finite_numeric(out, period_column)
    n_cad = period / CADENCE_D
    nearest = np.rint(n_cad)
    delta = np.abs(n_cad - nearest)
    base = f"{prefix}cadence_" if prefix else "cadence_"
    out[f"{base}period_n"] = n_cad
    out[f"{base}period_nearest_n"] = nearest
    out[f"{base}alias_delta"] = delta
    out[f"{base}alias_flag"] = np.isfinite(delta) & (delta <= float(tolerance))
    return out


def _sample_balanced_cells(
    frame: pd.DataFrame,
    *,
    n: int,
    random_state: int,
    cell_columns: tuple[str, ...],
) -> pd.DataFrame:
    if n <= 0 or frame.empty:
        return frame.head(0).copy()
    if len(frame) <= n:
        return frame.sample(frac=1.0, random_state=random_state).reset_index(drop=True)

    work = frame.copy()
    for column in cell_columns:
        if column not in work:
            work[column] = "missing"
        work[column] = work[column].fillna("missing").astype(str)

    groups = list(work.groupby(list(cell_columns), dropna=False, sort=True))
    per_cell = max(1, int(np.ceil(n / max(len(groups), 1))))
    pieces: list[pd.DataFrame] = []
    for idx, (_, group) in enumerate(groups):
        take = min(len(group), per_cell)
        pieces.append(group.sample(n=take, random_state=random_state + idx))
    selected = pd.concat(pieces, ignore_index=False) if pieces else work.head(0)
    if len(selected) > n:
        selected = selected.sample(n=n, random_state=random_state)
    elif len(selected) < n:
        selected_ids = set(selected.index.tolist())
        remaining = work.loc[[idx for idx in work.index if idx not in selected_ids]]
        fill_n = min(n - len(selected), len(remaining))
        if fill_n > 0:
            fill = remaining.sample(n=fill_n, random_state=random_state + 10_000)
            selected = pd.concat([selected, fill], ignore_index=False)
    return selected.sample(frac=1.0, random_state=random_state + 20_000).reset_index(drop=True)


def _sort_real_candidates(frame: pd.DataFrame) -> pd.DataFrame:
    sort_cols: list[str] = []
    ascending: list[bool] = []
    if "class_rank" in frame:
        sort_cols.append("class_rank")
        ascending.append(True)
    if "sde_max" in frame:
        sort_cols.append("sde_max")
        ascending.append(False)
    if not sort_cols:
        return frame
    return frame.sort_values(sort_cols, ascending=ascending, na_position="last", kind="stable")


def _take_real_bucket(
    frame: pd.DataFrame,
    *,
    name: str,
    mask: pd.Series,
    n: int,
    used_tics: set[int],
    random_state: int,
    random_pick: bool,
    available_counts: dict[str, int],
) -> pd.DataFrame:
    mask = mask.fillna(False).astype(bool)
    available_counts[name] = int(mask.sum())
    available = frame.loc[mask & ~frame["tic"].astype(int).isin(used_tics)].copy()
    if n <= 0 or available.empty:
        return available.head(0)
    if random_pick:
        selected = available.sample(n=min(n, len(available)), random_state=random_state)
    else:
        selected = _sort_real_candidates(available).head(n)
    selected = selected.copy()
    selected["selection_bucket"] = name
    selected["source_bucket"] = name
    used_tics.update(selected["tic"].astype(int).tolist())
    return selected


def select_real_pool(
    path: Path,
    *,
    n_real: int,
    random_state: int,
    alias_tolerance: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if n_real <= 0:
        return pd.DataFrame(), {}
    df = _read_table(path).copy()
    if "tic" not in df:
        raise KeyError(f"real candidate table missing tic: {path}")
    df["tic"] = pd.to_numeric(df["tic"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["tic"]).copy()
    df["tic"] = df["tic"].astype(int)
    df = _coerce_numeric(
        df,
        (
            "period_d",
            "t0_bjd",
            "duration_min",
            "sde_max",
            "depth",
            "depth_snr",
            "tmag",
            "n_apertures_agree",
            "period_cluster_count",
        ),
    )
    df = _with_cadence_alias_columns(
        df,
        period_column="period_d",
        tolerance=alias_tolerance,
    )
    rep_aperture = df.get("rep_aperture", pd.Series("", index=df.index)).fillna("").astype(str)
    finite_ephem = (
        df["period_d"].notna()
        & df["t0_bjd"].notna()
        & df["duration_min"].notna()
        & rep_aperture.isin(KNOWN_APERTURES)
    )
    df = df.loc[finite_ephem].drop_duplicates("tic", keep="first").reset_index(drop=True)
    if len(df) < n_real:
        raise ValueError(f"real candidate table has {len(df):,} usable rows; need {n_real:,}")

    sde = _finite_numeric(df, "sde_max")
    period = _finite_numeric(df, "period_d")
    duration = _finite_numeric(df, "duration_min")
    agree = _finite_numeric(df, "n_apertures_agree")
    cluster = _finite_numeric(df, "period_cluster_count")
    q40 = float(sde.quantile(0.40)) if sde.notna().any() else -np.inf
    q70 = float(sde.quantile(0.70)) if sde.notna().any() else -np.inf
    q90 = float(sde.quantile(0.90)) if sde.notna().any() else -np.inf
    duration_q80 = float(duration.quantile(0.80)) if duration.notna().any() else 30.0
    cluster_q95 = float(cluster.quantile(0.95)) if cluster.notna().any() else np.inf
    vet_class = df.get("vet_class", pd.Series("", index=df.index)).fillna("").astype(str)
    alias_like = (
        df["cadence_alias_flag"].fillna(False).astype(bool)
        | cluster.ge(cluster_q95).fillna(False)
        | (
            df.get("p_alias_pass", pd.Series(True, index=df.index))
            .fillna(True)
            .astype(str)
            .str.lower()
            .isin({"false", "0"})
        )
    )
    non_alias = ~alias_like

    used_tics: set[int] = set()
    available_counts: dict[str, int] = {}
    pieces: list[pd.DataFrame] = []
    quotas: list[tuple[str, pd.Series, int, bool]] = [
        ("real_wd1856_benchmark", df["tic"].eq(WD_1856_TIC), 1, False),
        (
            "real_high_sde_planet_like",
            vet_class.eq("planet_candidate") & non_alias,
            1800,
            False,
        ),
        (
            "real_eb_pceb_like",
            vet_class.isin(["pceb_grid_ceiling", "sub_roche_pceb_suspect"]) & non_alias,
            2200,
            False,
        ),
        (
            "real_variability_broad_duration",
            ((duration.ge(max(15.0, duration_q80))) | period.ge(3.0)) & non_alias,
            1100,
            True,
        ),
        (
            "real_aperture_disagreement",
            agree.le(1.0) & non_alias,
            1000,
            True,
        ),
        ("real_cadence_alias_systematic", alias_like, 700, False),
        ("real_mid_sde_control", sde.ge(q40) & sde.lt(q70) & non_alias, 1400, True),
        ("real_low_sde_control", sde.lt(q40) & non_alias, 800, True),
    ]
    for idx, (name, mask, quota, random_pick) in enumerate(quotas):
        remaining = n_real - sum(len(piece) for piece in pieces)
        if remaining <= 0:
            break
        pieces.append(
            _take_real_bucket(
                df,
                name=name,
                mask=mask,
                n=min(quota, remaining),
                used_tics=used_tics,
                random_state=random_state + idx,
                random_pick=random_pick,
                available_counts=available_counts,
            )
        )

    remaining = n_real - sum(len(piece) for piece in pieces)
    if remaining > 0:
        pieces.append(
            _take_real_bucket(
                df,
                name="real_stratified_fill",
                mask=non_alias,
                n=remaining,
                used_tics=used_tics,
                random_state=random_state + 100,
                random_pick=True,
                available_counts=available_counts,
            )
        )
    remaining = n_real - sum(len(piece) for piece in pieces)
    if remaining > 0:
        pieces.append(
            _take_real_bucket(
                df,
                name="real_alias_overflow_fill",
                mask=pd.Series(True, index=df.index),
                n=remaining,
                used_tics=used_tics,
                random_state=random_state + 200,
                random_pick=True,
                available_counts=available_counts,
            )
        )

    selected = pd.concat([piece for piece in pieces if len(piece)], ignore_index=True)
    if len(selected) != n_real:
        raise ValueError(f"selected {len(selected):,} real rows; expected {n_real:,}")
    selected = selected.sample(frac=1.0, random_state=random_state + 300).reset_index(drop=True)
    selected = _attach_real_review_columns(selected)
    bucket_counts = selected["selection_bucket"].value_counts().to_dict()
    selected["selection_weight"] = selected["selection_bucket"].map(
        lambda bucket: float(available_counts.get(str(bucket), 0) / max(bucket_counts.get(bucket, 1), 1))
    )
    selected["selection_note"] = "stratified_real_bls_vetter_no_ranker"
    summary = {
        "input_rows": int(len(df)),
        "selected_rows": int(len(selected)),
        "sde_quantiles": {"q40": q40, "q70": q70, "q90": q90},
        "duration_q80_min": duration_q80,
        "period_cluster_q95": cluster_q95,
        "selection_bucket_available_counts": available_counts,
        "selection_bucket_counts": _value_counts(selected, "selection_bucket"),
        "cadence_alias_selected": int(selected["cadence_alias_flag"].fillna(False).astype(bool).sum()),
    }
    return selected, summary


def _classify_injection_recovery(frame: pd.DataFrame) -> pd.Series:
    topn = frame.get("topn_recovery_status", pd.Series("", index=frame.index)).fillna("").astype(str)
    recovery = frame.get("recovery_status", pd.Series("", index=frame.index)).fillna("").astype(str)
    mode = pd.Series("inj_peak_mismatch_mid_snr", index=frame.index, dtype=object)
    mode.loc[topn.eq("bls_top1_recovered") | recovery.eq("bls_recovered")] = "inj_bls_top1_recovered"
    mode.loc[topn.isin(["bls_topn_recovered", "bls_topn_harmonic_match"])] = "inj_topn_or_harmonic"

    depth_snr = _finite_numeric(frame, "depth_snr")
    sde = _finite_numeric(frame, "sde_max")
    radius = _finite_numeric(frame, "truth_radius_rearth")
    low_mask = (
        mode.eq("inj_peak_mismatch_mid_snr")
        & (
            depth_snr.le(depth_snr.quantile(0.35)).fillna(False)
            | sde.le(sde.quantile(0.35)).fillna(False)
            | radius.le(radius.quantile(0.25)).fillna(False)
        )
    )
    mode.loc[low_mask] = "inj_low_snr_or_expected_missed"
    return mode


def _take_injected_bucket(
    frame: pd.DataFrame,
    *,
    name: str,
    mask: pd.Series,
    n: int,
    used_ids: set[str],
    random_state: int,
    available_counts: dict[str, int],
) -> pd.DataFrame:
    mask = mask.fillna(False).astype(bool)
    available_counts[name] = int(mask.sum())
    id_series = frame["injection_id"].fillna("").astype(str)
    available = frame.loc[mask & ~id_series.isin(used_ids)].copy()
    if n <= 0 or available.empty:
        return available.head(0)
    selected = _sample_balanced_cells(
        available,
        n=min(n, len(available)),
        random_state=random_state,
        cell_columns=("truth_grid_period_bin", "truth_grid_radius_bin"),
    )
    selected["selection_bucket"] = name
    used_ids.update(selected["injection_id"].fillna("").astype(str).tolist())
    return selected


def select_injected_pool(
    path: Path,
    *,
    n_injected: int,
    random_state: int,
    period_min_d: float,
    period_max_d: float,
    alias_tolerance: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if n_injected <= 0:
        return pd.DataFrame(), {}
    df = _read_table(path).copy()
    if "injection_id" not in df:
        raise KeyError(f"injected candidate table missing injection_id: {path}")
    df = _coerce_numeric(
        df,
        (
            "period_d",
            "t0_bjd",
            "duration_min",
            "sde_max",
            "depth",
            "depth_snr",
            "tmag",
            "truth_period_d",
            "truth_t0_bjd",
            "truth_duration_min",
            "truth_depth",
            "truth_radius_rearth",
            "truth_impact_b",
            "truth_inclination_deg",
            "truth_n_good_in_transit",
        ),
    )
    df = _with_cadence_alias_columns(df, period_column="period_d", tolerance=alias_tolerance)
    df = _with_cadence_alias_columns(
        df,
        period_column="truth_period_d",
        prefix="truth_",
        tolerance=alias_tolerance,
    )
    source_kind = df.get("source_kind", pd.Series("", index=df.index)).fillna("").astype(str)
    if source_kind.ne("").any():
        df = df.loc[source_kind.eq("injection_recovery")].copy()
    good_cadences = df.get("truth_n_good_in_transit", pd.Series(2, index=df.index)).fillna(0).astype(float)
    base = (
        df["period_d"].notna()
        & df["t0_bjd"].notna()
        & df["duration_min"].notna()
        & df.get("rep_aperture", pd.Series("", index=df.index)).fillna("").astype(str).isin(KNOWN_APERTURES)
        & df["truth_period_d"].between(period_min_d, period_max_d, inclusive="both")
        & df["truth_radius_rearth"].notna()
        & df["truth_depth"].notna()
        & df["truth_impact_b"].notna()
        & df["truth_inclination_deg"].notna()
        & good_cadences.ge(2)
    )
    df = df.loc[base].drop_duplicates("injection_id", keep="first").reset_index(drop=True)
    if len(df) < n_injected:
        raise ValueError(f"injected table has {len(df):,} usable rows; need {n_injected:,}")

    df["injection_recovery_selection_mode"] = _classify_injection_recovery(df)
    used_ids: set[str] = set()
    available_counts: dict[str, int] = {}
    pieces: list[pd.DataFrame] = []
    quotas = [
        ("inj_bls_top1_recovered", 300),
        ("inj_topn_or_harmonic", 200),
        ("inj_peak_mismatch_mid_snr", 300),
        ("inj_low_snr_or_expected_missed", 200),
    ]
    for idx, (name, quota) in enumerate(quotas):
        remaining = n_injected - sum(len(piece) for piece in pieces)
        if remaining <= 0:
            break
        pieces.append(
            _take_injected_bucket(
                df,
                name=name,
                mask=df["injection_recovery_selection_mode"].eq(name),
                n=min(quota, remaining),
                used_ids=used_ids,
                random_state=random_state + idx,
                available_counts=available_counts,
            )
        )
    remaining = n_injected - sum(len(piece) for piece in pieces)
    if remaining > 0:
        pieces.append(
            _take_injected_bucket(
                df,
                name="inj_balanced_fill",
                mask=pd.Series(True, index=df.index),
                n=remaining,
                used_ids=used_ids,
                random_state=random_state + 100,
                available_counts=available_counts,
            )
        )

    selected = pd.concat([piece for piece in pieces if len(piece)], ignore_index=True)
    if len(selected) != n_injected:
        raise ValueError(f"selected {len(selected):,} injected rows; expected {n_injected:,}")
    selected = selected.sample(frac=1.0, random_state=random_state + 300).reset_index(drop=True)
    if "truth_source_kind" not in selected:
        selected["truth_source_kind"] = selected.get("source_kind", "injection_recovery")
    if "truth_source_bucket" not in selected:
        selected["truth_source_bucket"] = selected.get("source_bucket", "")
    selected["source_kind"] = "injection_recovery"
    selected["source_bucket"] = selected["selection_bucket"]
    selected["truth_source_kind"] = "injection_recovery"
    selected["recovery_status"] = selected.get("recovery_status", "injection_recovery")
    selected["review_id"] = selected.get("review_id", selected["injection_id"].map(lambda x: f"inj:{x}"))
    bucket_counts = selected["selection_bucket"].value_counts().to_dict()
    selected["selection_weight"] = selected["selection_bucket"].map(
        lambda bucket: float(available_counts.get(str(bucket), 0) / max(bucket_counts.get(bucket, 1), 1))
    )
    selected["selection_note"] = "balanced_injected_batman_predetrend_no_ranker"
    summary = {
        "input_rows": int(len(df)),
        "selected_rows": int(len(selected)),
        "period_request_d": [float(period_min_d), float(period_max_d)],
        "truth_period_range_d": [
            float(selected["truth_period_d"].min()),
            float(selected["truth_period_d"].max()),
        ],
        "truth_radius_range_rearth": [
            float(selected["truth_radius_rearth"].min()),
            float(selected["truth_radius_rearth"].max()),
        ],
        "selection_bucket_available_counts": available_counts,
        "selection_bucket_counts": _value_counts(selected, "selection_bucket"),
        "topn_recovery_counts": _value_counts(selected, "topn_recovery_status"),
        "recovery_status_counts": _value_counts(selected, "recovery_status"),
    }
    if selected["truth_period_d"].min() > period_min_d:
        summary["warning"] = (
            "selected source table does not cover the requested shortest period; "
            f"min truth period is {selected['truth_period_d'].min():.6g} d"
        )
    return selected, summary


def _blind_for_browser(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "truth_source_bucket" not in out:
        out["truth_source_bucket"] = out.get("source_bucket", "")
    if "truth_source_kind" not in out:
        out["truth_source_kind"] = out.get("source_kind", "")
    if "truth_vet_class" not in out:
        out["truth_vet_class"] = out.get("vet_class", "")
    out["source_bucket"] = "review_candidate"
    out["vet_class"] = "review_candidate"
    return out


def _make_review_queue(
    pool: pd.DataFrame,
    *,
    review_real: int,
    review_injected: int,
    random_state: int,
) -> pd.DataFrame:
    source_kind = pool["source_kind"].fillna("").astype(str)
    real = pool.loc[source_kind.eq("real_candidate")].copy()
    injected = pool.loc[source_kind.eq("injection_recovery")].copy()
    if len(real) < review_real:
        raise ValueError(f"not enough real rows for review queue: {len(real)} < {review_real}")
    if len(injected) < review_injected:
        raise ValueError(f"not enough injected rows for review queue: {len(injected)} < {review_injected}")
    pieces = [
        real.sample(n=review_real, random_state=random_state + 10),
        injected.sample(n=review_injected, random_state=random_state + 20),
    ]
    return (
        pd.concat(pieces, ignore_index=True, sort=False)
        .sample(frac=1.0, random_state=random_state + 30)
        .reset_index(drop=True)
    )


def _write_counts(frame: pd.DataFrame, out_path: Path, columns: tuple[str, ...]) -> None:
    if frame.empty:
        pd.DataFrame(columns=[*columns, "n"]).to_csv(out_path, index=False)
        return
    usable = [column for column in columns if column in frame]
    if not usable:
        pd.DataFrame({"n": [len(frame)]}).to_csv(out_path, index=False)
        return
    counts = frame.groupby(usable, dropna=False).size().reset_index(name="n")
    counts.to_csv(out_path, index=False)


def _verify(
    *,
    pool: pd.DataFrame,
    review: pd.DataFrame,
    out_dir: Path,
    n_real_pool: int,
    n_injected_pool: int,
    review_real: int,
    review_injected: int,
    require_leo_reports: bool,
) -> dict[str, Any]:
    failures: list[str] = []
    if len(pool) != n_real_pool + n_injected_pool:
        failures.append(f"pool rows = {len(pool)}; expected {n_real_pool + n_injected_pool}")
    if len(review) != review_real + review_injected:
        failures.append(f"review rows = {len(review)}; expected {review_real + review_injected}")

    for name, frame, expected_real, expected_injected in (
        ("pool", pool, n_real_pool, n_injected_pool),
        ("review", review, review_real, review_injected),
    ):
        source = frame.get("source_kind", pd.Series("", index=frame.index)).fillna("").astype(str)
        n_real = int(source.eq("real_candidate").sum())
        n_inj = int(source.eq("injection_recovery").sum())
        if n_real != expected_real:
            failures.append(f"{name} real_candidate rows = {n_real}; expected {expected_real}")
        if n_inj != expected_injected:
            failures.append(f"{name} injection_recovery rows = {n_inj}; expected {expected_injected}")
        if "review_id" not in frame:
            failures.append(f"{name} missing review_id")
        else:
            dup = int(frame["review_id"].fillna("").astype(str).duplicated().sum())
            if dup:
                failures.append(f"{name} has {dup} duplicate review_id rows")
        for column in ("period_d", "t0_bjd", "duration_min"):
            values = _finite_numeric(frame, column)
            bad = int(values.isna().sum())
            if bad:
                failures.append(f"{name} has {bad} non-finite {column} rows")
        rep = frame.get("rep_aperture", pd.Series("", index=frame.index)).fillna("").astype(str)
        bad_rep = int((~rep.isin(KNOWN_APERTURES)).sum())
        if bad_rep:
            failures.append(f"{name} has {bad_rep} invalid representative apertures")

    for column in ("source_bucket", "vet_class"):
        values = review.get(column, pd.Series("", index=review.index)).fillna("").astype(str)
        bad = int(values.ne("review_candidate").sum())
        if bad:
            failures.append(f"review has {bad} non-blinded visible {column} rows")

    injected = pool.loc[pool.get("source_kind", pd.Series("", index=pool.index)).fillna("").astype(str).eq("injection_recovery")]
    for column in (
        "truth_period_d",
        "truth_radius_rearth",
        "truth_depth",
        "truth_impact_b",
        "truth_inclination_deg",
    ):
        values = _finite_numeric(injected, column)
        bad = int(values.isna().sum())
        if bad:
            failures.append(f"injected pool has {bad} non-finite {column} rows")

    reports_dir = out_dir / "vet_reports"
    if require_leo_reports:
        names = review.get("leo_report_name", pd.Series("", index=review.index)).fillna("").astype(str)
        missing_name = int(names.eq("").sum())
        if missing_name:
            failures.append(f"review has {missing_name} rows without leo_report_name")
        missing_files = [name for name in names[names.ne("")] if not (reports_dir / name).exists()]
        if missing_files:
            failures.append(f"review references {len(missing_files)} missing LEO PDFs; first={missing_files[0]}")

    label_path = out_dir / "human_labels_vetted.csv"
    if not label_path.exists():
        failures.append(f"missing label output placeholder: {label_path}")

    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "passed": not failures,
        "failures": failures,
        "n_pool_rows": int(len(pool)),
        "n_review_rows": int(len(review)),
        "pool_source_kind_counts": _value_counts(pool, "source_kind"),
        "review_source_kind_counts": _value_counts(review, "source_kind"),
        "review_truth_source_kind_counts": _value_counts(review, "truth_source_kind"),
        "pool_selection_bucket_counts": _value_counts(pool, "selection_bucket"),
        "review_selection_bucket_counts": _value_counts(review, "selection_bucket"),
        "leo_class_counts": _value_counts(review, "leo_class"),
        "require_leo_reports": bool(require_leo_reports),
    }


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real-candidates", type=Path, default=DEFAULT_REAL_CANDIDATES)
    parser.add_argument("--injected-candidates", type=Path, default=DEFAULT_INJECTED_CANDIDATES)
    parser.add_argument("--injection-h5", type=Path, default=DEFAULT_INJECTION_H5)
    parser.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    parser.add_argument("--star-catalog", type=Path, default=DEFAULT_STAR_CATALOG)
    parser.add_argument("--star-atmosphere-priority", default="H,He,mixed")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--n-real-pool", type=int, default=9000)
    parser.add_argument("--n-injected-pool", type=int, default=1000)
    parser.add_argument("--n-review", type=int, default=1000)
    parser.add_argument("--review-real", type=int, default=900)
    parser.add_argument("--review-injected", type=int, default=100)
    parser.add_argument("--random-state", type=int, default=5608)
    parser.add_argument("--injected-period-min-d", type=float, default=0.08)
    parser.add_argument("--injected-period-max-d", type=float, default=13.0)
    parser.add_argument("--cadence-alias-tolerance", type=float, default=0.02)
    parser.add_argument("--skip-leo", action="store_true")
    parser.add_argument("--max-leo-reports", type=int, default=1000)
    parser.add_argument("--leo-timeout-s", type=int, default=300)
    parser.add_argument("--leo-workers", type=int, default=8)
    parser.add_argument("--aperture", default="DET_FLUX_ADP015_SML")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.review_real + args.review_injected != args.n_review:
        raise ValueError("--review-real + --review-injected must equal --n-review")

    real, real_summary = select_real_pool(
        args.real_candidates,
        n_real=args.n_real_pool,
        random_state=args.random_state,
        alias_tolerance=args.cadence_alias_tolerance,
    )
    injected, injected_summary = select_injected_pool(
        args.injected_candidates,
        n_injected=args.n_injected_pool,
        random_state=args.random_state + 1000,
        period_min_d=args.injected_period_min_d,
        period_max_d=args.injected_period_max_d,
        alias_tolerance=args.cadence_alias_tolerance,
    )
    pool = _finalize_queue(real, injected)
    pool["teacher_pool_version"] = "s56_mixed_teacher_10k_v1"
    pool["pool_random_state"] = int(args.random_state)
    pool["truth_vet_class"] = pool.get("vet_class", "")

    atmosphere_priority = tuple(
        part.strip() for part in str(args.star_atmosphere_priority).split(",") if part.strip()
    )
    star_records: dict[int, dict[str, Any]] = {}
    if args.star_catalog is not None:
        tics = {int(tic) for tic in pd.to_numeric(pool["tic"], errors="coerce").dropna().astype(int)}
        star_records = load_wd_star_catalog(
            args.star_catalog,
            tics=tics,
            atmosphere_priority=atmosphere_priority,
        )
        pool = annotate_star_parameters(pool, star_records=star_records)

    pool_csv = args.out_dir / "mixed_teacher_pool.csv"
    pool.to_csv(pool_csv, index=False)

    review_unblinded = _make_review_queue(
        pool,
        review_real=args.review_real,
        review_injected=args.review_injected,
        random_state=args.random_state,
    )
    review = _blind_for_browser(review_unblinded)
    pre_leo_csv = args.out_dir / "review_queue_1k_pre_leo.csv"
    review.to_csv(pre_leo_csv, index=False)

    metrics = pd.DataFrame()
    if not args.skip_leo:
        review, metrics = render_leo_reports(
            review,
            out_dir=args.out_dir,
            hlsp_root=args.hlsp_root,
            injection_h5=args.injection_h5,
            aperture=args.aperture,
            max_reports=args.max_leo_reports,
            timeout_s=args.leo_timeout_s,
            overwrite=args.overwrite,
            workers=args.leo_workers,
        )
        metrics.to_csv(args.out_dir / "leo_metrics.csv", index=False)
    else:
        pd.DataFrame().to_csv(args.out_dir / "leo_metrics.csv", index=False)

    review_csv = args.out_dir / "review_queue_1k.csv"
    review.to_csv(review_csv, index=False)

    labels_path = args.out_dir / "human_labels_vetted.csv"
    if not labels_path.exists():
        pd.DataFrame(columns=LABEL_HEADER).to_csv(labels_path, index=False)

    _write_counts(pool, args.out_dir / "pool_bucket_counts.csv", ("source_kind", "selection_bucket"))
    _write_counts(review, args.out_dir / "review_bucket_counts.csv", ("source_kind", "selection_bucket"))
    _write_counts(
        pool.loc[pool["source_kind"].fillna("").astype(str).eq("injection_recovery")],
        args.out_dir / "pool_injected_period_radius_counts.csv",
        ("truth_grid_period_bin", "truth_grid_radius_bin", "selection_bucket"),
    )
    _write_counts(
        review.loc[review["source_kind"].fillna("").astype(str).eq("injection_recovery")],
        args.out_dir / "review_injected_period_radius_counts.csv",
        ("truth_grid_period_bin", "truth_grid_radius_bin", "selection_bucket"),
    )

    require_leo_reports = (not args.skip_leo) and (
        args.max_leo_reports <= 0 or args.max_leo_reports >= len(review)
    )
    verification = _verify(
        pool=pool,
        review=review,
        out_dir=args.out_dir,
        n_real_pool=args.n_real_pool,
        n_injected_pool=args.n_injected_pool,
        review_real=args.review_real,
        review_injected=args.review_injected,
        require_leo_reports=require_leo_reports,
    )
    _write_json(args.out_dir / "verification.json", verification)

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "s56_mixed_teacher_queue",
        "light_curve_product": "S56 TWIRL-FS v2 compare",
        "selection_uses_peak_ranker": False,
        "source_blinding": {
            "browser_visible_source_bucket": "review_candidate",
            "browser_visible_vet_class": "review_candidate",
            "hidden_truth_columns_retained": True,
        },
        "random_state": int(args.random_state),
        "n_real_pool": int(args.n_real_pool),
        "n_injected_pool": int(args.n_injected_pool),
        "n_review": int(len(review)),
        "review_real": int(args.review_real),
        "review_injected": int(args.review_injected),
        "injected_fraction_pool": float(args.n_injected_pool / (args.n_real_pool + args.n_injected_pool)),
        "injected_fraction_review": float(args.review_injected / max(args.n_review, 1)),
        "real_candidates": str(args.real_candidates),
        "injected_candidates": str(args.injected_candidates),
        "injection_h5": str(args.injection_h5),
        "hlsp_root": str(args.hlsp_root),
        "star_catalog": str(args.star_catalog) if args.star_catalog is not None else "",
        "star_atmosphere_priority": list(atmosphere_priority),
        "star_source_counts": star_source_counts(pool),
        "real_summary": real_summary,
        "injected_summary": injected_summary,
        "leo": {
            "skip_leo": bool(args.skip_leo),
            "max_leo_reports": int(args.max_leo_reports),
            "leo_workers": int(args.leo_workers),
            "leo_timeout_s": int(args.leo_timeout_s),
            "n_leo_metrics": int(len(metrics)),
            "n_leo_errors": int(metrics["error"].fillna("").astype(str).ne("").sum()) if "error" in metrics else 0,
            "n_leo_plot_errors": int(metrics["plot_error"].fillna("").astype(str).ne("").sum()) if "plot_error" in metrics else 0,
        },
        "outputs": {
            "mixed_teacher_pool_csv": str(pool_csv),
            "review_queue_1k_pre_leo_csv": str(pre_leo_csv),
            "review_queue_1k_csv": str(review_csv),
            "human_labels_vetted_csv": str(labels_path),
            "verification_json": str(args.out_dir / "verification.json"),
            "summary_json": str(args.out_dir / "summary.json"),
            "vet_reports": str(args.out_dir / "vet_reports"),
        },
        "verification_passed": bool(verification["passed"]),
        "verification_failures": verification["failures"],
    }
    _write_json(args.out_dir / "summary.json", summary)
    print("[mixed-teacher] complete")
    print(f"  pool: {len(pool):,} rows -> {pool_csv}")
    print(f"  review: {len(review):,} rows -> {review_csv}")
    print(f"  labels: {labels_path}")
    print(f"  verification passed: {verification['passed']}")
    if verification["failures"]:
        for failure in verification["failures"]:
            print(f"  failure: {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
