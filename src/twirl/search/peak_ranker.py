"""Injected-truth BLS peak ranker.

This module trains a lightweight model to rank BLS peaks using injection truth.
It is not an astrophysical classifier: the target is only whether a BLS peak
matches the injected ephemeris. Features must therefore be available for real
BLS peaks too; truth columns are used only for labels and audit metrics.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from twirl.vetting.self_training import (
    FeatureConfig,
    LogisticConfig,
    SoftmaxModel,
    append_predictions,
    fit_feature_spec,
    load_model,
    save_model,
    train_softmax_model,
)


DEFAULT_ID_COLUMN = "injection_id"
POSITIVE_LABEL = "signal_peak"
NEGATIVE_LABEL = "background_peak"
DEFAULT_FEATURE_COLUMNS: tuple[str, ...] = (
    "peak_rank",
    "period_d",
    "duration_min",
    "qtran",
    "depth",
    "depth_snr",
    "sde",
    "log_power",
    "n_cad_total",
    "n_cad_quality",
    "n_cad_kept",
    "n_cad_edge_trimmed",
    "n_cad_sigma_clipped",
    "dropout_frac",
    "quality_dropout_frac",
    "baseline_d",
    "n_orbits",
    "tmag",
    "rank_inv",
    "log_period_d",
    "log_duration_min",
    "log_abs_depth",
    "log_depth_snr",
    "log_sde",
)
DEFAULT_CATEGORICAL_COLUMNS: tuple[str, ...] = ("aperture",)


@dataclass(frozen=True)
class PeakRankerConfig:
    """Training/evaluation controls for the peak ranker."""

    id_column: str = DEFAULT_ID_COLUMN
    feature_columns: tuple[str, ...] = DEFAULT_FEATURE_COLUMNS
    categorical_columns: tuple[str, ...] = DEFAULT_CATEGORICAL_COLUMNS
    validation_fraction: float = 0.20
    test_fraction: float = 0.20
    random_state: int = 56
    max_iter: int = 500
    l2: float = 1.0e-3
    score_prefix: str = "ranker"


@dataclass
class PeakRankerResult:
    model: SoftmaxModel
    scored_peaks: pd.DataFrame
    selected_peaks: pd.DataFrame
    summary: dict[str, Any]
    coefficients: pd.DataFrame


def _group_rank_inv(values: pd.Series, group: pd.Series, *, ascending: bool = False) -> pd.Series:
    rank = values.groupby(group).rank(method="min", ascending=ascending)
    return pd.Series(np.where(rank > 0, 1.0 / rank, np.nan), index=values.index)


def _group_max_ratio(values: pd.Series, group: pd.Series) -> pd.Series:
    denom = values.groupby(group).transform("max")
    return pd.Series(np.where(denom > 0, values / denom, np.nan), index=values.index)


def add_peak_features(peaks: pd.DataFrame, *, group_column: str | None = None) -> pd.DataFrame:
    """Add non-truth engineered features to a peak table."""

    df = peaks.copy()
    if "n_cad_quality" not in df.columns and "n_cad_kept" in df.columns:
        df["n_cad_quality"] = df["n_cad_kept"]
    if "n_cad_edge_trimmed" not in df.columns:
        df["n_cad_edge_trimmed"] = 0.0
    if "n_cad_sigma_clipped" not in df.columns:
        df["n_cad_sigma_clipped"] = 0.0
    if "quality_dropout_frac" not in df.columns and "dropout_frac" in df.columns:
        df["quality_dropout_frac"] = df["dropout_frac"]
    for col in ("peak_rank", "period_d", "duration_min", "depth", "depth_snr", "sde", "log_power"):
        if col not in df.columns:
            df[col] = np.nan
    rank = pd.to_numeric(df["peak_rank"], errors="coerce")
    period = pd.to_numeric(df["period_d"], errors="coerce")
    duration = pd.to_numeric(df["duration_min"], errors="coerce")
    depth = pd.to_numeric(df["depth"], errors="coerce")
    depth_snr = pd.to_numeric(df["depth_snr"], errors="coerce")
    sde = pd.to_numeric(df["sde"], errors="coerce")
    log_power = pd.to_numeric(df["log_power"], errors="coerce")

    df["rank_inv"] = np.where(rank > 0, 1.0 / rank, np.nan)
    df["log_period_d"] = np.log10(np.where(period > 0, period, np.nan))
    df["log_duration_min"] = np.log10(np.where(duration > 0, duration, np.nan))
    df["log_abs_depth"] = np.log10(np.where(np.abs(depth) > 0, np.abs(depth), np.nan))
    df["log_depth_snr"] = np.log10(np.where(depth_snr > 0, depth_snr, np.nan))
    df["log_sde"] = np.log10(np.where(sde > 0, sde, np.nan))
    if group_column is not None and group_column in df.columns:
        group = df[group_column].fillna("").astype(str)
        abs_depth = depth.abs()
        df["sde_group_rank_inv"] = _group_rank_inv(sde, group, ascending=False)
        df["depth_snr_group_rank_inv"] = _group_rank_inv(depth_snr, group, ascending=False)
        df["log_power_group_rank_inv"] = _group_rank_inv(log_power, group, ascending=False)
        df["abs_depth_group_rank_inv"] = _group_rank_inv(abs_depth, group, ascending=False)
        df["sde_over_group_max"] = _group_max_ratio(sde, group)
        df["depth_snr_over_group_max"] = _group_max_ratio(depth_snr, group)
        df["abs_depth_over_group_max"] = _group_max_ratio(abs_depth, group)
    return df


def prepare_peak_training_frame(peaks: pd.DataFrame, id_column: str = DEFAULT_ID_COLUMN) -> pd.DataFrame:
    """Filter candidate peaks and attach binary labels."""

    if id_column not in peaks.columns:
        raise KeyError(f"missing ID column: {id_column}")
    if "is_candidate_peak" in peaks.columns:
        mask = peaks["is_candidate_peak"].fillna(False).astype(bool)
        df = peaks.loc[mask].copy()
    else:
        df = peaks.copy()
    if df.empty:
        raise ValueError("peak table has no candidate peak rows")
    if "is_injected_signal_peak" not in df.columns:
        raise KeyError("peak table missing is_injected_signal_peak")
    df = add_peak_features(df, group_column=id_column)
    signal = df["is_injected_signal_peak"].fillna(False).astype(bool)
    df["peak_label"] = np.where(signal, POSITIVE_LABEL, NEGATIVE_LABEL)
    group_has_signal = df.groupby(id_column)["is_injected_signal_peak"].transform("any").astype(bool)
    df["rankable_injection"] = group_has_signal
    return df


def split_injection_ids(
    ids: pd.Series | np.ndarray,
    *,
    validation_fraction: float,
    test_fraction: float,
    random_state: int,
) -> dict[str, set[str]]:
    """Make target-level train/validation/test splits."""

    unique = np.asarray(sorted({str(v) for v in ids}), dtype=object)
    if unique.size < 3:
        return {"train": set(unique.tolist()), "validation": set(), "test": set()}
    rng = np.random.default_rng(random_state)
    rng.shuffle(unique)
    n = unique.size
    n_test = int(round(test_fraction * n))
    n_val = int(round(validation_fraction * n))
    n_test = min(max(n_test, 1 if test_fraction > 0 else 0), n - 2)
    n_val = min(max(n_val, 1 if validation_fraction > 0 else 0), n - n_test - 1)
    test = set(unique[:n_test].tolist())
    validation = set(unique[n_test:n_test + n_val].tolist())
    train = set(unique[n_test + n_val:].tolist())
    return {"train": train, "validation": validation, "test": test}


def _assign_splits(df: pd.DataFrame, cfg: PeakRankerConfig) -> pd.Series:
    splits = split_injection_ids(
        df[cfg.id_column].astype(str),
        validation_fraction=cfg.validation_fraction,
        test_fraction=cfg.test_fraction,
        random_state=cfg.random_state,
    )
    values = df[cfg.id_column].astype(str)
    out = pd.Series("train", index=df.index, dtype=object)
    out.loc[values.isin(splits["validation"])] = "validation"
    out.loc[values.isin(splits["test"])] = "test"
    return out


def _signal_score_column(scored: pd.DataFrame, prefix: str) -> str:
    col = f"{prefix}_p_{POSITIVE_LABEL}"
    if col not in scored.columns:
        raise KeyError(f"missing signal probability column: {col}")
    return col


def _recall_at_k(
    df: pd.DataFrame,
    *,
    id_column: str,
    order_column: str,
    ascending: bool,
    k_values: tuple[int, ...],
) -> dict[int, dict[str, Any]]:
    out: dict[int, dict[str, Any]] = {}
    if df.empty:
        return {
            k: {"n": 0, "denom": 0, "fraction": float("nan"), "rankable_n": 0, "rankable_fraction": float("nan")}
            for k in k_values
        }
    work = df.copy()
    work["_truth"] = work["is_injected_signal_peak"].fillna(False).astype(bool)
    rankable_by_group = work.groupby(id_column)["_truth"].any()
    denom = int(rankable_by_group.size)
    rankable_n = int(rankable_by_group.sum())
    sorted_work = work.sort_values([id_column, order_column], ascending=[True, ascending], kind="stable")
    for k in k_values:
        top = sorted_work.groupby(id_column, group_keys=False).head(k)
        recovered = top.groupby(id_column)["_truth"].any()
        recovered = rankable_by_group.index.to_series().map(recovered).eq(True)
        n_recovered = int(recovered.sum())
        out[int(k)] = {
            "n": n_recovered,
            "denom": denom,
            "fraction": float(n_recovered / denom) if denom else float("nan"),
            "rankable_n": rankable_n,
            "rankable_fraction": float(n_recovered / rankable_n) if rankable_n else float("nan"),
        }
    return out


def _evaluate_split(
    df: pd.DataFrame,
    *,
    id_column: str,
    score_column: str,
    k_values: tuple[int, ...],
) -> dict[str, Any]:
    if df.empty:
        return {
            "n_peak_rows": 0,
            "n_injections": 0,
            "n_rankable_injections": 0,
            "model_recall_at_k": {},
            "bls_rank_recall_at_k": {},
            "sde_recall_at_k": {},
        }
    truth_by_group = df.groupby(id_column)["is_injected_signal_peak"].any()
    return {
        "n_peak_rows": int(len(df)),
        "n_injections": int(truth_by_group.size),
        "n_rankable_injections": int(truth_by_group.sum()),
        "model_recall_at_k": _recall_at_k(
            df, id_column=id_column, order_column=score_column, ascending=False, k_values=k_values
        ),
        "bls_rank_recall_at_k": _recall_at_k(
            df, id_column=id_column, order_column="peak_rank", ascending=True, k_values=k_values
        ),
        "sde_recall_at_k": _recall_at_k(
            df, id_column=id_column, order_column="sde", ascending=False, k_values=k_values
        ),
    }


def _selected_peaks(
    scored: pd.DataFrame,
    *,
    id_column: str,
    score_column: str,
) -> pd.DataFrame:
    idx = (
        scored.sort_values([id_column, score_column], ascending=[True, False], kind="stable")
        .groupby(id_column, sort=False)
        .head(1)
        .index
    )
    keep_cols = list(dict.fromkeys(
        col for col in (
            id_column,
            "split",
            "tic",
            "tmag",
            "aperture",
            "peak_rank",
            "period_d",
            "t0_bjd",
            "duration_min",
            "depth",
            "depth_snr",
            "sde",
            score_column,
            "is_injected_signal_peak",
            "match_kind",
            "truth_period_d",
            "truth_t0_bjd",
            "truth_radius_rearth",
            "truth_impact_b",
        )
        if col in scored.columns
    ))
    return scored.loc[idx, keep_cols].reset_index(drop=True)


def score_peak_table(
    peaks: pd.DataFrame,
    model: SoftmaxModel,
    *,
    score_prefix: str = "ranker",
    candidate_only: bool = True,
    group_column: str | None = None,
) -> pd.DataFrame:
    """Score injected or real BLS peak rows with a trained peak ranker."""

    df = peaks.copy()
    if candidate_only:
        if "status" in df.columns:
            df = df.loc[df["status"].fillna("").astype(str).eq("ok")].copy()
        if "peak_rank" in df.columns:
            df = df.loc[pd.to_numeric(df["peak_rank"], errors="coerce").fillna(0).gt(0)].copy()
    if df.empty:
        return df
    df = add_peak_features(df, group_column=group_column)
    # Models trained on newer tables may expect cadence-cleaning columns that
    # older real-candidate tables do not have. Fill any remaining model
    # features as NaN; the fitted preprocessing will replace them robustly.
    for col in model.feature_spec.numeric_columns:
        if col not in df.columns:
            df[col] = np.nan
    for col in model.feature_spec.categorical_levels:
        if col not in df.columns:
            df[col] = ""
    return append_predictions(df, model, score_prefix)


def select_ranked_ephemerides(
    scored: pd.DataFrame,
    *,
    id_column: str = "tic",
    score_prefix: str = "ranker",
    top_n: int = 1,
) -> pd.DataFrame:
    """Select the top-N ranker-scored ephemerides per target."""

    if id_column not in scored.columns:
        raise KeyError(f"missing ID column: {id_column}")
    score_col = _signal_score_column(scored, score_prefix)
    work = scored.sort_values([id_column, score_col], ascending=[True, False], kind="stable").copy()
    work["ranker_selection_rank"] = work.groupby(id_column).cumcount() + 1
    work = work.loc[work["ranker_selection_rank"].le(int(top_n))]
    keep_cols = list(dict.fromkeys(
        col for col in (
            id_column,
            "ranker_selection_rank",
            "tic",
            "sector",
            "cam",
            "ccd",
            "tmag",
            "aperture",
            "peak_rank",
            "period_d",
            "t0_bjd",
            "duration_min",
            "depth",
            "depth_snr",
            "sde",
            score_col,
            "n_cad_total",
            "n_cad_quality",
            "n_cad_kept",
            "n_cad_edge_trimmed",
            "n_cad_sigma_clipped",
            "dropout_frac",
            "quality_dropout_frac",
            "n_orbits",
            "baseline_d",
            "status",
            "hlsp_path",
            "match_kind",
            "is_injected_signal_peak",
        )
        if col in work.columns
    ))
    return work.loc[:, keep_cols].reset_index(drop=True)


def _coefficient_table(model: SoftmaxModel) -> pd.DataFrame:
    if POSITIVE_LABEL not in model.classes:
        return pd.DataFrame()
    class_index = list(model.classes).index(POSITIVE_LABEL)
    names = ("intercept", *model.feature_spec.feature_names)
    coeff = model.weights[:, class_index]
    return (
        pd.DataFrame({"feature": names, "coefficient": coeff})
        .assign(abs_coefficient=lambda df: df["coefficient"].abs())
        .sort_values("abs_coefficient", ascending=False)
        .reset_index(drop=True)
    )


def train_peak_ranker(
    peaks: pd.DataFrame,
    cfg: PeakRankerConfig | None = None,
    k_values: tuple[int, ...] = (1, 2, 3, 5, 10, 20),
) -> PeakRankerResult:
    """Train and evaluate a BLS peak ranker from injected peak rows."""

    cfg = cfg or PeakRankerConfig()
    df = prepare_peak_training_frame(peaks, id_column=cfg.id_column)
    df["split"] = _assign_splits(df, cfg)
    train_mask = df["split"].eq("train")
    train = df.loc[train_mask].copy()
    if train["peak_label"].nunique() < 2:
        raise ValueError("training split must contain both signal and background peaks")
    train_rankable_by_group = train.groupby(cfg.id_column)["is_injected_signal_peak"].any()

    feature_spec = fit_feature_spec(
        df,
        FeatureConfig(
            feature_columns=cfg.feature_columns,
            categorical_columns=cfg.categorical_columns,
            exclude_columns=(),
            clip_sigma=8.0,
        ),
    )
    model = train_softmax_model(
        train,
        train["peak_label"].astype(str),
        feature_spec,
        LogisticConfig(max_iter=cfg.max_iter, l2=cfg.l2, class_weight_balanced=True),
        metadata={
            "stage": "injection_peak_ranker",
            "id_column": cfg.id_column,
            "n_train_peak_rows": int(len(train)),
            "n_train_injections": int(train[cfg.id_column].nunique()),
            "n_train_rankable_injections": int(train_rankable_by_group.sum()),
            "n_train_background_only_injections": int((~train_rankable_by_group).sum()),
        },
    )
    scored = append_predictions(df, model, cfg.score_prefix)
    score_col = _signal_score_column(scored, cfg.score_prefix)
    selected = _selected_peaks(scored, id_column=cfg.id_column, score_column=score_col)

    summary = {
        "n_peak_rows": int(len(df)),
        "n_candidate_peak_rows": int(len(df)),
        "n_injections": int(df[cfg.id_column].nunique()),
        "n_rankable_injections": int(df.groupby(cfg.id_column)["is_injected_signal_peak"].any().sum()),
        "feature_names": list(model.feature_spec.feature_names),
        "classes": list(model.classes),
        "splits": {},
    }
    for split_name in ("train", "validation", "test", "all"):
        split_df = scored if split_name == "all" else scored.loc[scored["split"].eq(split_name)]
        summary["splits"][split_name] = _evaluate_split(
            split_df,
            id_column=cfg.id_column,
            score_column=score_col,
            k_values=k_values,
        )
    return PeakRankerResult(
        model=model,
        scored_peaks=scored,
        selected_peaks=selected,
        summary=summary,
        coefficients=_coefficient_table(model),
    )


def write_peak_ranker_outputs(result: PeakRankerResult, out_dir: Path, *, write_scored: bool = True) -> None:
    """Write model, scored peaks, selections, coefficients, and summary."""

    import json

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    save_model(result.model, out_dir / "peak_ranker_model.npz")
    if write_scored:
        result.scored_peaks.to_csv(out_dir / "scored_peaks.csv", index=False)
    result.selected_peaks.to_csv(out_dir / "selected_ephemerides.csv", index=False)
    result.coefficients.to_csv(out_dir / "coefficients.csv", index=False)
    (out_dir / "summary.json").write_text(json.dumps(result.summary, indent=2, sort_keys=True) + "\n")


__all__ = [
    "PeakRankerConfig",
    "PeakRankerResult",
    "add_peak_features",
    "load_model",
    "prepare_peak_training_frame",
    "score_peak_table",
    "select_ranked_ephemerides",
    "split_injection_ids",
    "train_peak_ranker",
    "write_peak_ranker_outputs",
]
