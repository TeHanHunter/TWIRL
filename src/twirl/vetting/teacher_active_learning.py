"""Active-learning controls for the S56 harmonic-CNN teacher.

This module intentionally separates three kinds of state:

* legacy human reviews, which remain immutable provenance;
* model scores, which are hidden queue-selection metadata; and
* browser rows, which contain only the evidence needed for a new decision.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.adp_only import ADP_ONLY_APERTURES
from twirl.vetting.label_schema import LABEL_OPTIONS


A2V1_TEACHER_INPUT_CONTRACT = "s56_A2v1_adp_raw_pair_v1"
LEGACY_CANONICAL_REVIEW_CONTRACT = "s56_legacy_canonical_review_v1"
ACTIVE_LEARNING_POLICY_VERSION = "s56_teacher_v1_active_learning_v1"

HARMONIC_FACTORS: tuple[float, ...] = (0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0)
HARMONIC_LABELS: tuple[str, ...] = ("P/4", "P/3", "P/2", "P", "2P", "3P", "4P")
SIGNAL_LABELS: frozenset[str] = frozenset(
    {
        "planet_like",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "wide_transit_like",
    }
)
OTHER_LABELS: frozenset[str] = frozenset(
    {"instrumental_or_systematic", "centroid_contaminant", "uncertain"}
)


@dataclass(frozen=True)
class ActiveLearningQuotas:
    planet_preserve: int = 300
    eclipse_contact: int = 200
    smooth_variable: int = 150
    broad_dip: int = 100
    disagreement_harmonic: int = 150
    stratified_control: int = 100

    @property
    def total(self) -> int:
        return int(sum(self.__dict__.values()))


DEFAULT_ACTIVE_LEARNING_QUOTAS = ActiveLearningQuotas()


def _read_table(path: Path, *, preserve_strings: bool = False) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        if preserve_strings:
            return pd.read_csv(path, dtype=str, keep_default_na=False)
        return pd.read_csv(path, low_memory=False)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"unsupported table format: {path}")


def _write_table(frame: pd.DataFrame, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".csv":
        frame.to_csv(path, index=False)
        return path
    try:
        frame.to_parquet(path, compression="zstd", index=False)
        return path
    except (ImportError, ModuleNotFoundError, ValueError):
        fallback = path.with_suffix(".csv")
        frame.to_csv(fallback, index=False)
        return fallback


def _json(path: Path, payload: Mapping[str, Any]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n")


def _franklin_candidate_key(row: Mapping[str, Any]) -> str:
    """Reproduce the standalone Franklin app's string-preserving key."""

    return "|".join(
        str(row.get(column, "") or "")
        for column in ("review_id", "tic", "sector", "period_d", "t0_bjd")
    )


def _latest_franklin_labels(labels: pd.DataFrame) -> pd.DataFrame:
    if labels.empty:
        return labels.copy()
    work = labels.copy()
    if "updated_utc" in work:
        work = work.sort_values("updated_utc", kind="stable")
    if "row_id" not in work:
        raise KeyError("Franklin labels are missing row_id")
    return work.drop_duplicates("row_id", keep="last")


def effective_human_label(frame: pd.DataFrame) -> pd.Series:
    """Coalesce raw/effective labels, then override only covered adjudications."""

    result = pd.Series("", index=frame.index, dtype=object)
    for column in ("label", "human_label"):
        if column not in frame:
            continue
        values = frame[column].fillna("").astype(str)
        present = values.ne("")
        result.loc[present] = values.loc[present]
    if "human_label_adjudicated" in frame:
        adjudicated = frame["human_label_adjudicated"].fillna("").astype(str)
        covered = adjudicated.ne("")
        result.loc[covered] = adjudicated.loc[covered]
    return result


def _period_match(
    review_period: np.ndarray,
    adp_period: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ratio = review_period / adp_period
    factors = np.asarray(HARMONIC_FACTORS, dtype=float)
    relative = np.abs(ratio[:, None] / factors[None, :] - 1.0)
    relative = np.where(np.isfinite(relative), relative, np.inf)
    nearest = relative.argmin(axis=1)
    return nearest, relative[np.arange(len(relative)), nearest], ratio


def audit_legacy_franklin_queue(
    queue: pd.DataFrame,
    adp_peaks: pd.DataFrame,
    *,
    labels: pd.DataFrame | None = None,
    tolerance: float = 0.02,
    max_peak_rank: int = 10,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Quarantine canonical-sheet reviews and quantify ADP period compatibility."""

    work = queue.copy().reset_index(drop=True)
    work["row_id"] = np.arange(len(work)).astype(str)
    work["candidate_key"] = [
        _franklin_candidate_key(row) for row in work.to_dict("records")
    ]
    tic = pd.to_numeric(work.get("tic"), errors="coerce")
    if tic.isna().any() or tic.duplicated().any():
        raise ValueError("Franklin queue must contain one finite unique TIC per row")
    work["tic"] = tic.astype(np.int64)

    peaks = adp_peaks.copy()
    required = {"tic", "aperture", "peak_rank", "period_d"}
    missing = sorted(required - set(peaks.columns))
    if missing:
        raise KeyError(f"ADP peak table is missing columns: {missing}")
    peaks = peaks.loc[
        peaks["aperture"].astype(str).eq(ADP_ONLY_APERTURES[0])
        & pd.to_numeric(peaks["peak_rank"], errors="coerce").between(1, max_peak_rank)
    ].copy()
    peaks["tic"] = pd.to_numeric(peaks["tic"], errors="coerce").astype("Int64")
    peaks["peak_rank"] = pd.to_numeric(peaks["peak_rank"], errors="coerce").astype("Int64")
    peaks["adp_period_d"] = pd.to_numeric(peaks["period_d"], errors="coerce")
    peaks = peaks.dropna(subset=["tic", "peak_rank", "adp_period_d"])

    expanded = work.loc[:, ["row_id", "tic", "period_d"]].merge(
        peaks.loc[:, ["tic", "peak_rank", "adp_period_d"]],
        on="tic",
        how="left",
        validate="one_to_many",
    )
    review_period = pd.to_numeric(expanded["period_d"], errors="coerce").to_numpy(dtype=float)
    adp_period = expanded["adp_period_d"].to_numpy(dtype=float)
    nearest, relative_error, ratio = _period_match(review_period, adp_period)
    expanded["nearest_factor"] = np.asarray(HARMONIC_LABELS, dtype=object)[nearest]
    expanded["period_rel_error"] = relative_error
    expanded["period_ratio_review_over_adp"] = ratio
    expanded["harmonic_match"] = relative_error <= float(tolerance)
    expanded["direct_match"] = np.isfinite(ratio) & (np.abs(ratio - 1.0) <= float(tolerance))

    rank1 = expanded.loc[expanded["peak_rank"].eq(1)].copy()
    rank1 = rank1.rename(
        columns={
            "adp_period_d": "adp_rank1_period_d",
            "nearest_factor": "adp_rank1_factor",
            "period_rel_error": "adp_rank1_period_rel_error",
            "harmonic_match": "adp_rank1_harmonic_match",
            "direct_match": "adp_rank1_direct_match",
        }
    )
    rank1_columns = [
        "row_id",
        "adp_rank1_period_d",
        "adp_rank1_factor",
        "adp_rank1_period_rel_error",
        "adp_rank1_harmonic_match",
        "adp_rank1_direct_match",
    ]
    best = expanded.sort_values(
        ["row_id", "period_rel_error", "peak_rank"], kind="stable"
    ).drop_duplicates("row_id", keep="first")
    best = best.rename(
        columns={
            "peak_rank": "adp_best_match_peak_rank",
            "adp_period_d": "adp_best_match_period_d",
            "nearest_factor": "adp_best_match_factor",
            "period_rel_error": "adp_best_match_period_rel_error",
            "harmonic_match": "adp_topn_harmonic_match",
        }
    )
    best_columns = [
        "row_id",
        "adp_best_match_peak_rank",
        "adp_best_match_period_d",
        "adp_best_match_factor",
        "adp_best_match_period_rel_error",
        "adp_topn_harmonic_match",
    ]
    work = work.merge(rank1[rank1_columns], on="row_id", how="left", validate="one_to_one")
    work = work.merge(best[best_columns], on="row_id", how="left", validate="one_to_one")

    work["legacy_review_contract"] = LEGACY_CANONICAL_REVIEW_CONTRACT
    work["legacy_canonical_sheet"] = True
    work["active_adp_morphology_eligible"] = False
    work["adp_re_review_required"] = True
    work["legacy_label_status"] = "unlabeled"
    work["human_label"] = ""
    work["human_notes"] = ""
    work["human_labeler"] = ""
    work["human_updated_utc"] = ""

    if labels is not None and not labels.empty:
        latest = _latest_franklin_labels(labels)
        latest["row_id"] = latest["row_id"].astype(str)
        expected = work.set_index("row_id")["candidate_key"]
        joined = latest.loc[:, [
            column
            for column in (
                "row_id", "candidate_key", "label", "notes", "labeler", "updated_utc"
            )
            if column in latest
        ]].copy()
        joined = joined.rename(
            columns={
                "candidate_key": "human_candidate_key",
                "label": "joined_human_label",
                "notes": "joined_human_notes",
                "labeler": "joined_human_labeler",
                "updated_utc": "joined_human_updated_utc",
            }
        )
        if "human_candidate_key" not in joined:
            raise KeyError("Franklin labels are missing candidate_key")
        mismatch = joined.loc[
            joined["human_candidate_key"].astype(str).ne(
                joined["row_id"].map(expected).fillna("").astype(str)
            )
        ]
        if not mismatch.empty:
            raise ValueError(
                f"{len(mismatch)} Franklin labels have candidate-key mismatches"
            )
        work = work.merge(joined, on="row_id", how="left", validate="one_to_one")
        for source, destination in (
            ("joined_human_label", "human_label"),
            ("joined_human_notes", "human_notes"),
            ("joined_human_labeler", "human_labeler"),
            ("joined_human_updated_utc", "human_updated_utc"),
        ):
            if source in work:
                work[destination] = work[source].fillna("").astype(str)
        labeled = work["human_label"].ne("")
        work.loc[labeled, "legacy_label_status"] = "audit_only_pending_adp_review"
        work = work.drop(columns=[column for column in work if column.startswith("joined_")])

    rank1_direct = work["adp_rank1_direct_match"].fillna(False).astype(bool)
    rank1_harmonic = work["adp_rank1_harmonic_match"].fillna(False).astype(bool)
    topn_harmonic = work["adp_topn_harmonic_match"].fillna(False).astype(bool)
    human_label = work["human_label"].fillna("").astype(str)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "legacy_review_contract": LEGACY_CANONICAL_REVIEW_CONTRACT,
        "active_input_apertures": list(ADP_ONLY_APERTURES),
        "n_queue_rows": int(len(work)),
        "n_human_labels_available": int(human_label.ne("").sum()),
        "human_label_counts": {
            str(key): int(value)
            for key, value in human_label[human_label.ne("")].value_counts().sort_index().items()
        },
        "adp_rank1_direct_match_fraction": float(rank1_direct.mean()),
        "adp_rank1_harmonic_match_fraction": float(rank1_harmonic.mean()),
        "adp_topn_harmonic_match_fraction": float(topn_harmonic.mean()),
        "max_peak_rank": int(max_peak_rank),
        "period_tolerance": float(tolerance),
        "n_active_adp_morphology_eligible": 0,
        "policy": (
            "Legacy canonical-sheet labels remain immutable audit provenance and do not "
            "supervise ADP morphology until an ADP re-review or a separately validated "
            "label-transfer gate passes."
        ),
    }
    return work, summary


def run_legacy_franklin_audit(
    *,
    queue_path: Path,
    adp_peaks_path: Path,
    labels_path: Path | None,
    out_dir: Path,
) -> dict[str, Any]:
    queue = _read_table(queue_path, preserve_strings=True)
    peaks = _read_table(adp_peaks_path)
    labels = (
        _read_table(labels_path, preserve_strings=True)
        if labels_path is not None and Path(labels_path).exists()
        else pd.DataFrame()
    )
    rows, summary = audit_legacy_franklin_queue(queue, peaks, labels=labels)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    row_path = _write_table(rows, out_dir / "legacy_queue_adp_compatibility.csv")
    signal = rows["human_label"].isin(SIGNAL_LABELS)
    other = rows["human_label"].isin(OTHER_LABELS)
    signal_path = _write_table(
        rows.loc[signal].copy(), out_dir / "signal_labels_requiring_adp_review.csv"
    )
    control_candidates = rows.loc[other].copy()
    if len(control_candidates) > 100:
        control_candidates = control_candidates.sample(n=100, random_state=56)
    control_path = _write_table(
        control_candidates, out_dir / "other_label_transfer_audit_sample100.csv"
    )
    summary["outputs"] = {
        "compatibility_rows": str(row_path),
        "signal_re_review": str(signal_path),
        "other_transfer_audit": str(control_path),
        "summary": str(out_dir / "summary.json"),
    }
    _json(out_dir / "summary.json", summary)
    notice = (
        "# Franklin S56 Queue Status\n\n"
        "This queue is paused as a legacy canonical-flux review product. Existing labels "
        "must be preserved unchanged, but they are audit-only for the active ADP teacher "
        "until the rows are reviewed on S56-ADP-HV1/A2v1 evidence or pass the documented "
        "transfer gate. Do not continue the remaining canonical-sheet rows.\n"
    )
    (out_dir / "PAUSED_LEGACY_QUEUE.md").write_text(notice)
    return summary


def evaluate_a2v1_transfer_gate(
    compatibility_rows: pd.DataFrame,
    *,
    real_macro_f1: float,
    predicted_class_counts: Mapping[str, int],
    min_macro_f1: float = 0.70,
    min_ephemeris_match_fraction: float = 0.90,
) -> dict[str, Any]:
    """Gate teacher-v1 ranking on A2v1 without using the opened test for tuning."""

    if "a2v1_topn_harmonic_match" in compatibility_rows:
        match = compatibility_rows["a2v1_topn_harmonic_match"]
    elif "adp_topn_harmonic_match" in compatibility_rows:
        match = compatibility_rows["adp_topn_harmonic_match"]
    else:
        raise KeyError("compatibility table has no top-N harmonic-match column")
    if match.dtype != bool:
        match = match.fillna(False).astype(str).str.lower().isin({"1", "true", "yes", "y"})
    match_fraction = float(match.mean()) if len(match) else 0.0
    required_classes = {"planet_like", "eclipse_contact", "smooth_variable", "other"}
    nondegenerate = all(int(predicted_class_counts.get(label, 0)) > 0 for label in required_classes)
    checks = {
        "real_macro_f1": bool(np.isfinite(real_macro_f1) and real_macro_f1 >= min_macro_f1),
        "ephemeris_match": bool(match_fraction >= min_ephemeris_match_fraction),
        "nondegenerate_predictions": bool(nondegenerate),
    }
    return {
        "passed": bool(all(checks.values())),
        "checks": checks,
        "real_macro_f1": float(real_macro_f1),
        "min_real_macro_f1": float(min_macro_f1),
        "ephemeris_match_fraction": match_fraction,
        "min_ephemeris_match_fraction": float(min_ephemeris_match_fraction),
        "predicted_class_counts": {str(k): int(v) for k, v in predicted_class_counts.items()},
        "action": "rank_with_teacher_v1" if all(checks.values()) else "retrain_transfer_checkpoint",
    }


def build_a2v1_label_transfer_table(
    labeled_rows: pd.DataFrame,
    adp_peaks: pd.DataFrame,
    *,
    tolerance: float = 0.02,
    max_peak_rank: int = 10,
    max_scoring_rank: int = 3,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Match existing real human labels to newly searched A2v1 ephemerides."""

    labels = labeled_rows.copy()
    if not any(name in labels for name in ("human_label_adjudicated", "human_label", "label")):
        raise KeyError("labeled rows have no human label column")
    if "tic" not in labels or "period_d" not in labels:
        raise KeyError("labeled rows must contain tic and period_d")
    if "source_kind" in labels:
        labels = labels.loc[
            ~labels["source_kind"].fillna("").astype(str).str.contains("inject", case=False)
        ].copy()
    labels["human_label_transfer"] = effective_human_label(labels)
    labels = labels.loc[labels["human_label_transfer"].isin(LABEL_OPTIONS)].copy()
    labels["tic"] = pd.to_numeric(labels["tic"], errors="coerce").astype("Int64")
    labels["legacy_review_period_d"] = pd.to_numeric(labels["period_d"], errors="coerce")
    labels = labels.dropna(subset=["tic", "legacy_review_period_d"])
    labels = labels.sort_index(kind="stable").drop_duplicates("tic", keep="last")

    peaks = adp_peaks.copy()
    required = {
        "tic", "aperture", "peak_rank", "period_d", "t0_bjd", "duration_min", "status"
    }
    missing = sorted(required - set(peaks.columns))
    if missing:
        raise KeyError(f"A2v1 ADP peaks are missing columns: {missing}")
    if "source_product_tag" not in peaks or set(
        peaks["source_product_tag"].fillna("").astype(str)
    ) != {"A2v1"}:
        raise ValueError("label transfer requires an A2v1-stamped peak table")
    rank = pd.to_numeric(peaks["peak_rank"], errors="coerce")
    peaks = peaks.loc[
        peaks["aperture"].astype(str).eq(ADP_ONLY_APERTURES[0])
        & rank.between(1, int(max_peak_rank))
        & peaks["status"].astype(str).eq("ok")
    ].copy()
    peaks["tic"] = pd.to_numeric(peaks["tic"], errors="coerce").astype("Int64")
    peaks["a2v1_peak_rank"] = pd.to_numeric(peaks["peak_rank"], errors="coerce").astype("Int64")
    peaks["a2v1_period_d"] = pd.to_numeric(peaks["period_d"], errors="coerce")
    peaks["a2v1_t0_bjd"] = pd.to_numeric(peaks["t0_bjd"], errors="coerce")
    peaks["a2v1_duration_min"] = pd.to_numeric(peaks["duration_min"], errors="coerce")
    expanded = labels.merge(
        peaks.loc[
            :,
            [
                "tic",
                "a2v1_peak_rank",
                "a2v1_period_d",
                "a2v1_t0_bjd",
                "a2v1_duration_min",
            ],
        ],
        on="tic",
        how="left",
        validate="one_to_many",
    )
    nearest, relative_error, ratio = _period_match(
        expanded["legacy_review_period_d"].to_numpy(dtype=float),
        expanded["a2v1_period_d"].to_numpy(dtype=float),
    )
    expanded["a2v1_period_factor"] = np.asarray(HARMONIC_LABELS, dtype=object)[nearest]
    expanded["a2v1_period_rel_error"] = relative_error
    expanded["a2v1_period_ratio_legacy_over_new"] = ratio
    expanded["a2v1_topn_harmonic_match"] = relative_error <= float(tolerance)
    best = expanded.sort_values(
        ["tic", "a2v1_period_rel_error", "a2v1_peak_rank"], kind="stable"
    ).drop_duplicates("tic", keep="first")
    best["a2v1_scoring_eligible"] = (
        best["a2v1_topn_harmonic_match"].fillna(False).astype(bool)
        & pd.to_numeric(best["a2v1_peak_rank"], errors="coerce").le(int(max_scoring_rank))
    )
    best["transfer_status"] = np.select(
        [
            best["a2v1_scoring_eligible"],
            best["a2v1_topn_harmonic_match"].fillna(False).astype(bool),
        ],
        ["score_compatible_top3", "compatible_top10_not_scored"],
        default="no_compatible_a2v1_peak",
    )
    match = best["a2v1_topn_harmonic_match"].fillna(False).astype(bool)
    eligible = best["a2v1_scoring_eligible"].fillna(False).astype(bool)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_unique_real_labels": int(len(best)),
        "n_top10_harmonic_matches": int(match.sum()),
        "top10_harmonic_match_fraction": float(match.mean()) if len(best) else 0.0,
        "n_top3_scoring_eligible": int(eligible.sum()),
        "top3_scoring_eligible_fraction": float(eligible.mean()) if len(best) else 0.0,
        "max_peak_rank": int(max_peak_rank),
        "max_scoring_rank": int(max_scoring_rank),
        "period_tolerance": float(tolerance),
        "transfer_status_counts": {
            str(key): int(value)
            for key, value in best["transfer_status"].value_counts().sort_index().items()
        },
    }
    return best.reset_index(drop=True), summary


def evaluate_a2v1_scored_transfer(
    transfer_rows: pd.DataFrame,
    teacher_scores: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Evaluate A2v1 predictions against human labels without truth correction."""

    from twirl.vetting.adjudication_audit import MORPHOLOGY_CLASSES_V1, MORPHOLOGY_TARGET_MAP_V1
    from twirl.vetting.harmonic_training import classification_metrics

    eligible = transfer_rows.loc[
        transfer_rows["a2v1_scoring_eligible"].fillna(False).astype(bool)
    ].copy()
    scores = teacher_scores.copy()
    scores["tic"] = pd.to_numeric(scores["tic"], errors="coerce").astype("Int64")
    score_rank = pd.to_numeric(
        scores.get("rep_peak_rank", scores.get("adp_sml_peak_rank")), errors="coerce"
    ).astype("Int64")
    scores["a2v1_peak_rank"] = score_rank
    probability_columns = [f"p_{label}" for label in MORPHOLOGY_CLASSES_V1]
    missing = sorted(set(probability_columns) - set(scores.columns))
    if missing:
        raise KeyError(f"teacher scores are missing probabilities: {missing}")
    joined = eligible.merge(
        scores.loc[:, ["tic", "a2v1_peak_rank", *probability_columns]],
        on=["tic", "a2v1_peak_rank"],
        how="left",
        validate="one_to_one",
    )
    if joined[probability_columns].isna().any(axis=None):
        missing_rows = int(joined[probability_columns].isna().any(axis=1).sum())
        raise ValueError(f"{missing_rows} A2v1 transfer rows have no teacher score")
    mapped = joined["human_label_transfer"].map(MORPHOLOGY_TARGET_MAP_V1).fillna("")
    active = mapped.ne("")
    truth = np.asarray(
        [MORPHOLOGY_CLASSES_V1.index(value) for value in mapped.loc[active]], dtype=int
    )
    probability = joined.loc[active, probability_columns].to_numpy(dtype=float)
    metrics = classification_metrics(truth, probability, classes=MORPHOLOGY_CLASSES_V1)
    evaluation = joined.loc[active].copy()
    evaluation["morphology_target_transfer"] = mapped.loc[active].to_numpy()
    evaluation["morphology_prediction_transfer"] = np.asarray(
        MORPHOLOGY_CLASSES_V1, dtype=object
    )[probability.argmax(axis=1)]
    gate = evaluate_a2v1_transfer_gate(
        transfer_rows,
        real_macro_f1=float(metrics["macro_f1"]),
        predicted_class_counts=metrics.get("predicted_class_counts", {}),
    )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_evaluated": int(len(evaluation)),
        "human_morphology_metrics": metrics,
        "transfer_gate": gate,
        "truth_policy": "human morphology only; injection truth is not joined",
    }
    return evaluation, summary


def _numeric(frame: pd.DataFrame, name: str, default: float = np.nan) -> pd.Series:
    if name not in frame:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[name], errors="coerce")


def _probability(frame: pd.DataFrame, name: str) -> pd.Series:
    values = _numeric(frame, name, 0.0).fillna(0.0).clip(0.0, 1.0)
    return values


def _selection_scores(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    p_planet = _probability(out, "p_planet_like")
    p_eb = _probability(out, "p_eclipse_contact")
    p_variable = _probability(out, "p_smooth_variable")
    p_other = _probability(out, "p_other")
    p_preserve = _probability(out, "p_preserve")
    period = _numeric(out, "period_d")
    duration = _numeric(out, "duration_min")
    sde = _numeric(out, "sde_max")
    if sde.isna().all():
        sde = _numeric(out, "adp_sml_sde")
    duration_fraction = duration / (period * 1440.0)
    broad_rank = duration_fraction.rank(pct=True).fillna(0.0)
    long_low_sde = ((period >= 1.0) & (sde <= 40.0)).astype(float)
    entropy = _numeric(out, "morphology_entropy", 0.0).fillna(0.0)
    disagreement = _numeric(out, "ensemble_disagreement", 0.0).fillna(0.0)
    rare_harmonic_columns = [
        column
        for column in out.columns
        if column.startswith("p_harmonic_") and column != "p_harmonic_P"
    ]
    rare_harmonic = (
        out.loc[:, rare_harmonic_columns].apply(pd.to_numeric, errors="coerce").max(axis=1).fillna(0.0)
        if rare_harmonic_columns
        else pd.Series(0.0, index=out.index)
    )
    out["_score_planet_preserve"] = 0.55 * p_planet + 0.35 * p_preserve + 0.10 * long_low_sde
    out["_score_eclipse_contact"] = p_eb
    out["_score_smooth_variable"] = p_variable
    out["_score_broad_dip"] = p_preserve * broad_rank * (1.0 - 0.5 * p_eb) * (1.0 - 0.5 * p_variable)
    out["_score_disagreement_harmonic"] = 0.45 * entropy + 0.35 * disagreement + 0.20 * rare_harmonic
    out["_score_stratified_control"] = p_other + (1.0 - p_preserve)
    return out


def _collect_excluded_tics(tables: Iterable[pd.DataFrame]) -> set[int]:
    excluded: set[int] = set()
    for frame in tables:
        if frame is None or frame.empty or "tic" not in frame:
            continue
        values = pd.to_numeric(frame["tic"], errors="coerce").dropna().astype(np.int64)
        excluded.update(int(value) for value in values)
    return excluded


def _pick_rows(
    candidates: pd.DataFrame,
    *,
    score_column: str,
    n_rows: int,
    used_tics: set[int],
    random_control: bool = False,
    seed: int = 56,
) -> pd.DataFrame:
    available = candidates.loc[~candidates["tic"].astype(np.int64).isin(used_tics)].copy()
    available = available.sort_values(score_column, ascending=False, kind="stable")
    available = available.drop_duplicates("tic", keep="first")
    if random_control and len(available) > n_rows:
        sde = _numeric(available, "sde_max")
        if sde.isna().all():
            sde = _numeric(available, "adp_sml_sde")
        try:
            available["_control_sde_bin"] = pd.qcut(
                sde.rank(method="first"), q=4, labels=False, duplicates="drop"
            )
        except ValueError:
            available["_control_sde_bin"] = 0
        parts: list[pd.DataFrame] = []
        remaining = int(n_rows)
        groups = list(available.groupby("_control_sde_bin", dropna=False, sort=True))
        for index, (_, group) in enumerate(groups):
            take = remaining if index == len(groups) - 1 else int(np.ceil(n_rows / max(len(groups), 1)))
            sample_n = min(take, len(group), remaining)
            if sample_n:
                parts.append(group.sample(n=sample_n, random_state=seed + index))
                remaining -= sample_n
        selected = pd.concat(parts, ignore_index=False) if parts else available.head(0)
        if len(selected) < n_rows:
            supplement = available.drop(index=selected.index).head(n_rows - len(selected))
            selected = pd.concat([selected, supplement], ignore_index=False)
        available = selected
    return available.head(int(n_rows)).copy()


def build_active_learning_batch(
    scores: pd.DataFrame,
    *,
    batch_index: int,
    excluded_tables: Sequence[pd.DataFrame] = (),
    quotas: ActiveLearningQuotas = DEFAULT_ACTIVE_LEARNING_QUOTAS,
    seed: int = 56,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Select one deterministic 1k real-TIC batch and a blinded 10% overlap."""

    if quotas.total != 1000:
        raise ValueError(f"active-learning quotas must sum to 1000, got {quotas.total}")
    required = {
        "tic",
        "period_d",
        "t0_bjd",
        "duration_min",
        "p_planet_like",
        "p_eclipse_contact",
        "p_smooth_variable",
        "p_other",
        "p_preserve",
    }
    missing = sorted(required - set(scores.columns))
    if missing:
        raise KeyError(f"teacher score table is missing columns: {missing}")
    candidates = scores.copy()
    candidates["tic"] = pd.to_numeric(candidates["tic"], errors="coerce").astype("Int64")
    candidates = candidates.dropna(subset=["tic", "period_d", "t0_bjd", "duration_min"])
    candidates["tic"] = candidates["tic"].astype(np.int64)
    if "source_kind" in candidates:
        candidates = candidates.loc[
            ~candidates["source_kind"].fillna("").astype(str).str.contains("inject", case=False)
        ].copy()
    excluded_tics = _collect_excluded_tics(excluded_tables)
    candidates = candidates.loc[~candidates["tic"].isin(excluded_tics)].copy()
    candidates = _selection_scores(candidates)

    bucket_specs = (
        ("planet_preserve", quotas.planet_preserve, "_score_planet_preserve", False),
        ("eclipse_contact", quotas.eclipse_contact, "_score_eclipse_contact", False),
        ("smooth_variable", quotas.smooth_variable, "_score_smooth_variable", False),
        ("broad_dip", quotas.broad_dip, "_score_broad_dip", False),
        (
            "disagreement_harmonic",
            quotas.disagreement_harmonic,
            "_score_disagreement_harmonic",
            False,
        ),
        ("stratified_control", quotas.stratified_control, "_score_stratified_control", True),
    )
    used: set[int] = set()
    selected_parts: list[pd.DataFrame] = []
    for offset, (bucket, count, score_column, random_control) in enumerate(bucket_specs):
        available_count = int(
            candidates.loc[~candidates["tic"].astype(np.int64).isin(used), "tic"].nunique()
        )
        chosen = _pick_rows(
            candidates,
            score_column=score_column,
            n_rows=count,
            used_tics=used,
            random_control=random_control,
            seed=seed + 100 * int(batch_index) + offset,
        )
        if len(chosen) != count:
            raise ValueError(
                f"selection bucket {bucket} requested {count} unique TICs but found {len(chosen)}"
            )
        chosen["selection_bucket"] = bucket
        chosen["selection_score"] = chosen[score_column]
        chosen["selection_pool_size"] = available_count
        chosen["selection_fraction"] = float(count) / max(available_count, 1)
        chosen["selection_weight"] = float(available_count) / max(int(count), 1)
        used.update(int(value) for value in chosen["tic"])
        selected_parts.append(chosen)
    hidden = pd.concat(selected_parts, ignore_index=True, sort=False)
    hidden = hidden.sample(
        frac=1.0, random_state=seed + int(batch_index)
    ).reset_index(drop=True)
    hidden["active_learning_policy_version"] = ACTIVE_LEARNING_POLICY_VERSION
    hidden["input_contract_version"] = A2V1_TEACHER_INPUT_CONTRACT
    hidden["batch_index"] = int(batch_index)
    hidden["source_candidate_review_id"] = hidden.get(
        "review_id", pd.Series("", index=hidden.index)
    ).fillna("").astype(str)
    hidden["review_id"] = [
        f"s56_a2v1_teacher_b{int(batch_index):02d}_{index:04d}"
        for index in range(len(hidden))
    ]
    hidden["row_id"] = np.arange(len(hidden), dtype=int)
    hidden["candidate_key"] = [
        _franklin_candidate_key(row) for row in hidden.to_dict("records")
    ]
    hidden["double_review"] = False
    overlap_index = hidden.sample(
        n=100, random_state=seed + 10_000 + int(batch_index)
    ).index
    hidden.loc[overlap_index, "double_review"] = True

    hidden_prefixes = (
        "truth_",
        "source_",
        "selection_",
        "model_",
        "recovery_",
        "native_",
        "member_",
        "_score_",
        "p_",
        "std_p_",
    )
    hidden_columns = {
        column
        for column in hidden.columns
        if column.startswith(hidden_prefixes)
        or column
        in {
            "selection_bucket",
            "selection_score",
            "source_candidate_review_id",
            "active_learning_policy_version",
            "input_contract_version",
            "adp_only_contract_version",
            "bls_search_branch",
            "metadata_status",
            "metadata_error",
            "double_review",
            "morphology_entropy",
            "morphology_margin",
            "ensemble_disagreement",
            "model_version",
            "model_profile",
        }
    }
    queue = hidden.drop(columns=sorted(hidden_columns), errors="ignore").copy()
    for column in ("label", "label_source", "labeler", "notes", "updated_utc"):
        queue[column] = ""
    overlap = queue.loc[overlap_index].copy().reset_index(drop=True)
    overlap["row_id"] = np.arange(len(overlap), dtype=int)
    # Keep candidate_key stable across labelers so disagreements can be joined
    # exactly without exposing the overlap status in either browser queue.
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy_version": ACTIVE_LEARNING_POLICY_VERSION,
        "input_contract_version": A2V1_TEACHER_INPUT_CONTRACT,
        "batch_index": int(batch_index),
        "n_queue_rows": int(len(queue)),
        "n_unique_tics": int(queue["tic"].nunique()),
        "n_double_review": int(len(overlap)),
        "excluded_tics": int(len(excluded_tics)),
        "selection_bucket_counts": {
            str(key): int(value)
            for key, value in hidden["selection_bucket"].value_counts().sort_index().items()
        },
        "selection_weight_policy": (
            "inverse within-bucket selected fraction for audit only; not a prevalence weight"
        ),
    }
    return queue, overlap, hidden, summary


def write_active_learning_batch(
    *,
    scores_path: Path,
    out_dir: Path,
    batch_index: int,
    exclude_paths: Sequence[Path] = (),
) -> dict[str, Any]:
    scores = _read_table(scores_path)
    exclusions = [_read_table(path) for path in exclude_paths if Path(path).exists()]
    queue, overlap, hidden, summary = build_active_learning_batch(
        scores,
        batch_index=batch_index,
        excluded_tables=exclusions,
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue_path = _write_table(queue, out_dir / f"review_queue_batch{batch_index:02d}_1k.csv")
    overlap_path = _write_table(
        overlap, out_dir / f"double_review_queue_batch{batch_index:02d}_100.csv"
    )
    hidden_path = _write_table(
        hidden, out_dir / f"hidden_selection_provenance_batch{batch_index:02d}.parquet"
    )
    summary["outputs"] = {
        "review_queue": str(queue_path),
        "double_review_queue": str(overlap_path),
        "hidden_selection_provenance": str(hidden_path),
        "summary": str(out_dir / "summary.json"),
    }
    _json(out_dir / "summary.json", summary)
    return summary


def teacher_v2_readiness(
    labels: pd.DataFrame,
    *,
    test_metrics: Mapping[str, Any] | None = None,
    unresolved_ephemerides: int = 0,
    join_mismatches: int = 0,
) -> dict[str, Any]:
    """Return morphology-output and student gates from deduplicated real labels."""

    work = labels.copy()
    if "source_kind" in work:
        work = work.loc[
            ~work["source_kind"].fillna("").astype(str).str.contains("inject", case=False)
        ].copy()
    if not any(name in work for name in ("human_label_adjudicated", "human_label", "label")):
        raise KeyError("label table has no human label column")
    key_column = "tic" if "tic" in work else "review_id"
    if key_column not in work:
        raise KeyError("label table has neither tic nor review_id")
    work["_effective_human_label"] = effective_human_label(work)
    work = work.loc[work["_effective_human_label"].isin(LABEL_OPTIONS)].copy()
    work = work.drop_duplicates(key_column, keep="last")
    counts = work["_effective_human_label"].value_counts()
    planet = int(counts.get("planet_like", 0))
    broad = int(counts.get("wide_transit_like", 0))
    n_outputs = 5 if planet >= 50 and broad >= 50 else 4
    metrics = dict(test_metrics or {})
    if "test_metrics" in metrics and isinstance(metrics["test_metrics"], Mapping):
        metrics = dict(metrics["test_metrics"])
    if "morphology_by_source" in metrics and isinstance(
        metrics["morphology_by_source"], Mapping
    ):
        source_metrics = metrics["morphology_by_source"]
        if isinstance(source_metrics.get("real"), Mapping):
            real_metrics = dict(source_metrics["real"])
            calibration = metrics.get("calibration", {})
            if isinstance(calibration, Mapping) and isinstance(calibration.get("real"), Mapping):
                real_metrics["ece"] = calibration["real"].get("ece", np.nan)
            metrics = real_metrics
    per_class = metrics.get("per_class", {}) if isinstance(metrics.get("per_class", {}), Mapping) else {}
    planet_metrics = per_class.get("planet_like", {}) if isinstance(per_class, Mapping) else {}
    checks = {
        "real_planet_support": planet >= 50,
        "real_planet_test_support": int(planet_metrics.get("n", 0)) >= 10,
        "balanced_accuracy": float(metrics.get("balanced_accuracy", np.nan)) >= 0.75,
        "real_planet_recall": float(planet_metrics.get("recall", np.nan)) >= 0.80,
        "eb_recall": float(per_class.get("eclipse_contact", {}).get("recall", np.nan)) >= 0.60,
        "variable_recall": float(per_class.get("smooth_variable", {}).get("recall", np.nan)) >= 0.60,
        "other_recall": float(per_class.get("other", {}).get("recall", np.nan)) >= 0.60,
        "ece": float(metrics.get("ece", np.nan)) <= 0.10,
        "resolved_ephemerides": int(unresolved_ephemerides) == 0,
        "label_joins": int(join_mismatches) == 0,
    }
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "real_label_counts": {str(key): int(value) for key, value in counts.sort_index().items()},
        "teacher_output_count": n_outputs,
        "teacher_target_policy": (
            "s56_harmonic_cnn_v2_five_way"
            if n_outputs == 5
            else "s56_harmonic_cnn_v1_four_way_broad_preserve_only"
        ),
        "teacher_classes": (
            ["planet_like", "eclipse_contact", "smooth_variable", "broad_dip", "other"]
            if n_outputs == 5
            else ["planet_like", "eclipse_contact", "smooth_variable", "other"]
        ),
        "broad_class_promoted": bool(n_outputs == 5),
        "student_ready": bool(all(checks.values())),
        "student_checks": checks,
        "unresolved_ephemerides": int(unresolved_ephemerides),
        "join_mismatches": int(join_mismatches),
    }


def sector_rollout_readiness(
    joined_batch_labels: pd.DataFrame,
    *,
    transfer_gate_passed: bool,
    product_qa_passed: bool,
    expected_rows: int = 1000,
    min_enriched_planets: int = 10,
    min_enrichment_ratio: float = 3.0,
) -> dict[str, Any]:
    """Gate S57+ teacher inference on measured S56 Planet enrichment."""

    required = {"selection_bucket"}
    missing = sorted(required - set(joined_batch_labels.columns))
    if missing:
        raise KeyError(f"joined batch labels are missing columns: {missing}")
    work = joined_batch_labels.copy()
    work["_effective_human_label"] = effective_human_label(work)
    labeled = work["_effective_human_label"].isin(LABEL_OPTIONS)
    enriched = work["selection_bucket"].astype(str).eq("planet_preserve") & labeled
    control = work["selection_bucket"].astype(str).eq("stratified_control") & labeled
    enriched_planets = int(
        work.loc[enriched, "_effective_human_label"].eq("planet_like").sum()
    )
    control_planets = int(
        work.loc[control, "_effective_human_label"].eq("planet_like").sum()
    )
    enriched_n = int(enriched.sum())
    control_n = int(control.sum())
    enriched_rate = enriched_planets / enriched_n if enriched_n else 0.0
    control_rate = control_planets / control_n if control_n else 0.0
    enrichment_ratio = (
        enriched_rate / control_rate
        if control_rate > 0
        else float("inf")
        if enriched_rate > 0
        else 0.0
    )
    checks = {
        "batch_complete": int(labeled.sum()) == int(expected_rows),
        "minimum_enriched_planets": enriched_planets >= int(min_enriched_planets),
        "minimum_enrichment_ratio": enrichment_ratio >= float(min_enrichment_ratio),
        "transfer_gate": bool(transfer_gate_passed),
        "product_qa": bool(product_qa_passed),
    }
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "s57_plus_ready": bool(all(checks.values())),
        "checks": checks,
        "n_labeled": int(labeled.sum()),
        "expected_rows": int(expected_rows),
        "planet_preserve_n": enriched_n,
        "planet_preserve_planets": enriched_planets,
        "planet_preserve_rate": float(enriched_rate),
        "control_n": control_n,
        "control_planets": control_planets,
        "control_rate": float(control_rate),
        "enrichment_ratio": float(enrichment_ratio),
        "min_enriched_planets": int(min_enriched_planets),
        "min_enrichment_ratio": float(min_enrichment_ratio),
    }


__all__ = [
    "A2V1_TEACHER_INPUT_CONTRACT",
    "ACTIVE_LEARNING_POLICY_VERSION",
    "ActiveLearningQuotas",
    "DEFAULT_ACTIVE_LEARNING_QUOTAS",
    "LEGACY_CANONICAL_REVIEW_CONTRACT",
    "audit_legacy_franklin_queue",
    "build_active_learning_batch",
    "build_a2v1_label_transfer_table",
    "evaluate_a2v1_scored_transfer",
    "effective_human_label",
    "evaluate_a2v1_transfer_gate",
    "run_legacy_franklin_audit",
    "sector_rollout_readiness",
    "teacher_v2_readiness",
    "write_active_learning_batch",
]
