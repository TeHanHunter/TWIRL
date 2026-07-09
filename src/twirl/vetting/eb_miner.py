"""Active-learning helpers for finding more S56 EB/PCEB examples."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.adp_only import (
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
    assert_adp_only_search_frame,
    assert_adp_only_training_frame,
    classify_period_relation,
)
from twirl.vetting.recovery50_teacher import (
    candidate_key,
    join_queue_labels,
    json_default,
    read_table,
    write_table,
)


EB_HUMAN_LABEL = "eclipsing_binary_or_pceb"
EB_TARGET_POSITIVE = "eb_pceb"
EB_TARGET_NEGATIVE = "not_eb"
EB_EXCLUDED_LABELS = frozenset({"", "skip"})
DEFAULT_RANDOM_STATE = 56017

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


@dataclass(frozen=True)
class EBTrainingConfig:
    """Controls for the small, intentionally imbalanced EB miner."""

    random_state: int = DEFAULT_RANDOM_STATE
    max_uncertain_negatives: int = 100
    max_negatives_per_label: int = 200
    min_positive_rows: int = 5


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n")


def _safe_int(value: Any, default: int = -1) -> int:
    try:
        if pd.isna(value):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _safe_float(value: Any) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return out if np.isfinite(out) else float("nan")


def _review_id_for_row(row: pd.Series, *, source_pool: str) -> str:
    review_id = str(row.get("review_id", "") or "").strip()
    if review_id:
        return review_id
    tic = _safe_int(row.get("tic"))
    if source_pool.startswith("adp_bls"):
        aperture = str(row.get("aperture", "")).replace("DET_FLUX_", "").lower()
        return f"real:{tic}:adp_bls:{aperture}:peak:{_safe_int(row.get('peak_rank'), 1)}"
    return f"real:{tic}"


def _sheet_name(review_id: Any, branch_name: str = "current_adp") -> str:
    safe = str(review_id).replace("/", "_").replace(":", "_")
    return f"{safe}_twirl_twoap_{branch_name}.png"


def _read_joined_label_source(path: Path, *, source_name: str, priority: int) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    frame = read_table(path).copy()
    if "human_label" not in frame and "label" in frame:
        frame["human_label"] = frame["label"]
    if "human_notes" not in frame and "notes" in frame:
        frame["human_notes"] = frame["notes"]
    if "row_id" not in frame:
        frame.insert(0, "row_id", np.arange(len(frame), dtype=int))
    for col in ("source_kind", "review_id", "candidate_key", "human_label", "human_notes", "human_updated_utc"):
        if col not in frame:
            frame[col] = ""
    frame["eb_label_source_name"] = source_name
    frame["eb_label_source_path"] = str(path)
    frame["eb_label_priority"] = int(priority)
    return frame


def _read_queue_label_source(
    queue_csv: Path,
    labels_csv: Path,
    *,
    source_name: str,
    priority: int,
) -> pd.DataFrame:
    if not queue_csv.exists() or not labels_csv.exists() or labels_csv.stat().st_size == 0:
        return pd.DataFrame()
    frame = join_queue_labels(queue_csv, labels_csv)
    frame["eb_label_source_name"] = source_name
    frame["eb_label_source_path"] = str(labels_csv)
    frame["eb_label_priority"] = int(priority)
    return frame


def _dedupe_labeled_rows(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    out = frame.copy()
    if "review_id" not in out:
        out["review_id"] = ""
    if "candidate_key" not in out:
        out["candidate_key"] = out.apply(candidate_key, axis=1)
    out["review_id"] = out["review_id"].fillna("").astype(str)
    out["candidate_key"] = out["candidate_key"].fillna("").astype(str)
    out["_dedupe_key"] = np.where(out["review_id"].ne(""), out["review_id"], out["candidate_key"])
    out["_human_updated_sort"] = pd.to_datetime(out.get("human_updated_utc", ""), errors="coerce")
    out = out.sort_values(
        ["_dedupe_key", "eb_label_priority", "_human_updated_sort"],
        na_position="first",
        kind="stable",
    )
    out = out.drop_duplicates("_dedupe_key", keep="last")
    return out.drop(columns=["_dedupe_key", "_human_updated_sort"], errors="ignore").reset_index(drop=True)


def add_eb_miner_targets(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach the EB-vs-rest target while preserving the raw human label."""

    out = frame.copy()
    label = out.get("human_label", pd.Series("", index=out.index)).fillna("").astype(str)
    out["eb_miner_target"] = np.where(label.eq(EB_HUMAN_LABEL), EB_TARGET_POSITIVE, EB_TARGET_NEGATIVE)
    out.loc[label.isin(EB_EXCLUDED_LABELS), "eb_miner_target"] = ""
    out["eb_miner_include"] = out["eb_miner_target"].ne("")
    out["main_teacher_target"] = out["eb_miner_target"]
    out["main_teacher_include"] = out["eb_miner_include"]
    out["eb_miner_positive_policy"] = "positive_only_explicit_eclipsing_binary_or_pceb"
    return out


def _cap_eb_training_rows(frame: pd.DataFrame, cfg: EBTrainingConfig) -> pd.DataFrame:
    include = frame["eb_miner_include"].fillna(False).astype(bool)
    work = frame.loc[include].copy()
    positives = work.loc[work["eb_miner_target"].eq(EB_TARGET_POSITIVE)].copy()
    negatives = work.loc[work["eb_miner_target"].eq(EB_TARGET_NEGATIVE)].copy()
    if len(positives) < cfg.min_positive_rows:
        raise ValueError(f"only {len(positives)} EB/PCEB positives; need at least {cfg.min_positive_rows}")

    rng = np.random.default_rng(cfg.random_state)
    pieces = [positives]
    label = negatives.get("human_label", pd.Series("", index=negatives.index)).fillna("").astype(str)
    for raw_label, group in negatives.groupby(label, dropna=False):
        cap = cfg.max_uncertain_negatives if str(raw_label) == "uncertain" else cfg.max_negatives_per_label
        if len(group) > cap:
            take = rng.choice(group.index.to_numpy(), size=int(cap), replace=False)
            group = group.loc[np.sort(take)]
        pieces.append(group)
    out = pd.concat(pieces, ignore_index=False).sample(frac=1.0, random_state=cfg.random_state).reset_index(drop=True)
    if "row_id" in out:
        out["source_row_id"] = out["row_id"]
    out["row_id"] = np.arange(len(out), dtype=int)
    return out


def build_eb_miner_training_table(
    *,
    joined_tables: Sequence[Path],
    queue_label_pairs: Sequence[tuple[Path, Path, str]],
    out_dir: Path,
    cfg: EBTrainingConfig = EBTrainingConfig(),
) -> dict[str, Any]:
    """Build a real-only EB-vs-rest training table from human labels."""

    out_dir.mkdir(parents=True, exist_ok=True)
    pieces: list[pd.DataFrame] = []
    priority = 0
    for path in joined_tables:
        priority += 1
        frame = _read_joined_label_source(path, source_name=path.parent.name, priority=priority)
        if not frame.empty:
            pieces.append(frame)
    for queue_csv, labels_csv, source_name in queue_label_pairs:
        priority += 1
        frame = _read_queue_label_source(queue_csv, labels_csv, source_name=source_name, priority=priority)
        if not frame.empty:
            pieces.append(frame)
    if not pieces:
        raise ValueError("no label sources produced rows")

    joined = pd.concat(pieces, ignore_index=True, sort=False)
    for col in ("source_kind", "human_label", "review_id"):
        if col not in joined:
            joined[col] = ""
        joined[col] = joined[col].fillna("").astype(str)
    real = joined.loc[joined["source_kind"].eq("real_candidate")].copy()
    real = real.loc[real["human_label"].ne("")].copy()
    deduped = _dedupe_labeled_rows(real)
    assert_adp_only_training_frame(deduped)
    labeled = add_eb_miner_targets(deduped)
    audit_path = write_table(labeled, out_dir / "eb_miner_labeled_audit.csv")
    training = _cap_eb_training_rows(labeled, cfg)
    training_path = write_table(training, out_dir / "eb_miner_training_table.csv")

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "out_dir": str(out_dir),
        "n_input_sources": int(len(pieces)),
        "n_joined_labeled_rows": int(len(joined)),
        "n_real_labeled_rows_before_dedupe": int(len(real)),
        "n_real_labeled_rows_after_dedupe": int(len(labeled)),
        "n_training_rows": int(len(training)),
        "target_counts_all_real": labeled["eb_miner_target"].fillna("").astype(str).value_counts().sort_index().to_dict(),
        "target_counts_training": training["eb_miner_target"].fillna("").astype(str).value_counts().sort_index().to_dict(),
        "human_label_counts_all_real": labeled["human_label"].fillna("").astype(str).value_counts().sort_index().to_dict(),
        "human_label_counts_training": training["human_label"].fillna("").astype(str).value_counts().sort_index().to_dict(),
        "source_counts_training": training.get("eb_label_source_name", pd.Series(dtype=str)).fillna("").astype(str).value_counts().sort_index().to_dict(),
        "config": asdict(cfg),
        "outputs": {
            "audit": str(audit_path),
            "training_table": str(training_path),
            "summary": str(out_dir / "summary.json"),
        },
    }
    _write_json(out_dir / "summary.json", summary)
    return summary


def _normalize_candidate_table(
    path: Path,
    *,
    source_pool: str,
    small_peaks_per_tic: int = 3,
    primary_peaks_per_tic: int = 1,
) -> pd.DataFrame:
    frame = read_table(path).copy()
    if frame.empty:
        return frame
    assert_adp_only_search_frame(frame)
    if "tic" not in frame:
        raise KeyError(f"candidate table missing tic: {path}")
    rank = pd.to_numeric(frame["peak_rank"], errors="coerce")
    aperture = frame["aperture"].fillna("").astype(str)
    keep = (
        aperture.eq(ADP_ONLY_APERTURES[0]) & rank.between(1, int(small_peaks_per_tic))
    ) | (
        aperture.eq(ADP_ONLY_APERTURES[1]) & rank.between(1, int(primary_peaks_per_tic))
    )
    out = frame.loc[keep].copy()
    rename = {
        "sde": "sde_max",
        "aperture": "rep_aperture",
        "peak_rank": "rep_peak_rank",
    }
    for src, dst in rename.items():
        if dst not in out and src in out:
            out[dst] = out[src]
    if "sector" not in out:
        out["sector"] = 56
    if "source_bucket" not in out:
        out["source_bucket"] = "real_adp_bls"
    if "selection_bucket" not in out:
        out["selection_bucket"] = "real_adp_bls"
    if "vet_class" not in out:
        out["vet_class"] = "adp_bls"
    if "class_rank" not in out:
        out["class_rank"] = out.get("peak_rank", np.arange(1, len(out) + 1))
    if "blind_rank" not in out:
        out["blind_rank"] = np.arange(1, len(out) + 1)
    if "n_apertures_agree" not in out:
        out["n_apertures_agree"] = 1
    if "apertures_agree" not in out:
        out["apertures_agree"] = out.get("rep_aperture", "")
    for col in (
        "centroid_status",
        "centroid_pass",
        "centroid_delta_pix",
        "centroid_z",
        "n_in_transit",
        "n_oot_band",
    ):
        if col not in out:
            out[col] = ""
    out["tic"] = pd.to_numeric(out["tic"], errors="coerce").astype("Int64")
    out = out.dropna(subset=["tic"]).copy()
    out["tic"] = out["tic"].astype(int)
    for col in ("period_d", "t0_bjd", "duration_min"):
        out[col] = pd.to_numeric(out.get(col), errors="coerce")
    out = out.loc[out["period_d"].gt(0) & out["t0_bjd"].notna() & out["duration_min"].gt(0)].copy()
    out["review_id"] = out.apply(lambda row: _review_id_for_row(row, source_pool=source_pool), axis=1)
    out["source_kind"] = "real_candidate"
    out["recovery_status"] = "real_candidate"
    out["truth_source_kind"] = "real_candidate"
    out["truth_source_bucket"] = out["source_bucket"]
    out["source_pool"] = source_pool
    out["anchor_aperture"] = out["aperture"].astype(str)
    out["anchor_period_d"] = out["period_d"]
    out["anchor_t0_bjd"] = out["t0_bjd"]
    out["anchor_duration_min"] = out["duration_min"]
    out["anchor_sde"] = pd.to_numeric(out.get("sde"), errors="coerce")
    small = out["aperture"].eq(ADP_ONLY_APERTURES[0])
    metric_columns = (
        ("peak_rank", "peak_rank"),
        ("period_d", "period_d"),
        ("t0_bjd", "t0_bjd"),
        ("duration_min", "duration_min"),
        ("depth", "depth"),
        ("depth_snr", "depth_snr"),
        ("sde", "sde"),
        ("log_power", "log_power"),
    )
    source_rank = pd.to_numeric(frame["peak_rank"], errors="coerce")
    for source_aperture, prefix in zip(ADP_ONLY_APERTURES, ("adp_sml", "adp")):
        top1 = frame.loc[
            frame["aperture"].fillna("").astype(str).eq(source_aperture) & source_rank.eq(1)
        ].copy()
        top1 = top1.sort_values("tic", kind="stable").drop_duplicates("tic", keep="first")
        available = [(source, suffix) for source, suffix in metric_columns if source in top1]
        top1 = top1.loc[:, ["tic", *[source for source, _ in available]]].rename(
            columns={source: f"{prefix}_{suffix}" for source, suffix in available}
        )
        out = out.merge(top1, on="tic", how="left", validate="many_to_one")

    # The branch that generated this candidate uses that exact peak; the other
    # branch carries its rank-1 counterpart so combined models see both ADP
    # apertures rather than a systematic missing-metadata pattern.
    for source, suffix in metric_columns:
        values = pd.to_numeric(
            out[source] if source in out else pd.Series(np.nan, index=out.index),
            errors="coerce",
        )
        out.loc[small, f"adp_sml_{suffix}"] = values.loc[small]
        out.loc[~small, f"adp_{suffix}"] = values.loc[~small]
    relation = classify_period_relation(out["adp_sml_period_d"], out["adp_period_d"])
    out["aperture_period_relation"] = relation
    small_period = pd.to_numeric(out["adp_sml_period_d"], errors="coerce")
    primary_period = pd.to_numeric(out["adp_period_d"], errors="coerce")
    out["aperture_period_rel_delta"] = np.abs(primary_period - small_period) / small_period
    small_depth = pd.to_numeric(out["adp_sml_depth"], errors="coerce")
    primary_depth = pd.to_numeric(out["adp_depth"], errors="coerce")
    out["aperture_depth_ratio_primary_over_small"] = primary_depth / small_depth.replace(0, np.nan)
    out["aperture_disagreement_flag"] = relation.eq("unrelated")
    agrees = relation.isin({"exact", "harmonic"})
    out["n_apertures_agree"] = np.where(agrees, 2, 1)
    out["apertures_agree"] = np.where(
        agrees,
        ",".join(ADP_ONLY_APERTURES),
        out["aperture"].astype(str),
    )
    out["candidate_key"] = out.apply(candidate_key, axis=1)
    out["human_label"] = ""
    out["main_teacher_target"] = ""
    out["main_teacher_include"] = False
    return out.reset_index(drop=True)


def _collect_exclusions(paths: Iterable[Path]) -> tuple[set[str], set[int]]:
    review_ids: set[str] = set()
    tics: set[int] = set()
    for path in paths:
        if path is None or not Path(path).exists():
            continue
        try:
            frame = read_table(Path(path))
        except Exception:
            continue
        if "source_kind" in frame:
            frame = frame.loc[frame["source_kind"].fillna("").astype(str).eq("real_candidate")].copy()
        if "review_id" in frame:
            review_ids.update(frame["review_id"].fillna("").astype(str).replace("", np.nan).dropna().tolist())
        if "tic" in frame:
            tic = pd.to_numeric(frame["tic"], errors="coerce").dropna().astype(int)
            tics.update(tic.tolist())
    return review_ids, tics


def build_candidate_scoring_pool(
    *,
    primary_candidates: Path,
    fallback_candidates: Path | None = None,
    exclude_tables: Sequence[Path],
    out_dir: Path,
    small_peaks_per_tic: int = 3,
    primary_peaks_per_tic: int = 1,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    primary = _normalize_candidate_table(
        primary_candidates,
        source_pool="adp_bls",
        small_peaks_per_tic=small_peaks_per_tic,
        primary_peaks_per_tic=primary_peaks_per_tic,
    )
    fallback = pd.DataFrame()
    if fallback_candidates is not None:
        fallback = _normalize_candidate_table(
            fallback_candidates,
            source_pool="adp_bls_fallback",
            small_peaks_per_tic=small_peaks_per_tic,
            primary_peaks_per_tic=primary_peaks_per_tic,
        )
    candidates = pd.concat([primary, fallback], ignore_index=True, sort=False)
    candidates = candidates.drop_duplicates("review_id", keep="first").reset_index(drop=True)
    exclude_review_ids, exclude_tics = _collect_exclusions(exclude_tables)
    keep = ~candidates["review_id"].fillna("").astype(str).isin(exclude_review_ids)
    keep &= ~candidates["tic"].astype(int).isin(exclude_tics)
    candidates = candidates.loc[keep].copy().reset_index(drop=True)
    candidates["row_id"] = np.arange(len(candidates), dtype=int)
    candidates["adp_only_contract_version"] = ADP_ONLY_CONTRACT_VERSION
    assert_adp_only_training_frame(candidates)
    path = write_table(candidates, out_dir / "eb_miner_candidate_pool.csv")
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "primary_candidates": str(primary_candidates),
        "fallback_candidates": str(fallback_candidates) if fallback_candidates is not None else None,
        "n_primary_input": int(len(primary)),
        "n_fallback_input": int(len(fallback)),
        "n_after_exclusions": int(len(candidates)),
        "n_excluded_review_ids": int(len(exclude_review_ids)),
        "n_excluded_tics": int(len(exclude_tics)),
        "small_peaks_per_tic": int(small_peaks_per_tic),
        "primary_peaks_per_tic": int(primary_peaks_per_tic),
        "source_pool_counts": candidates["source_pool"].fillna("").astype(str).value_counts().sort_index().to_dict(),
        "outputs": {"candidate_pool": str(path), "summary": str(out_dir / "candidate_pool_summary.json")},
    }
    _write_json(out_dir / "candidate_pool_summary.json", summary)
    return summary


def aggregate_ensemble_scores(
    score_tables: Sequence[Path],
    *,
    positive_column: str = f"cnn_p_{EB_TARGET_POSITIVE}",
) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    for idx, path in enumerate(score_tables):
        if not Path(path).exists():
            continue
        frame = read_table(Path(path)).copy()
        if positive_column not in frame:
            raise KeyError(f"{path} missing {positive_column}")
        frame["_ensemble_member"] = idx
        pieces.append(frame)
    if not pieces:
        raise ValueError("no score tables supplied")
    scores = pd.concat(pieces, ignore_index=True, sort=False)
    key = "review_id" if "review_id" in scores else "row_id"
    first_cols = [col for col in scores.columns if col not in {positive_column, "_ensemble_member"} and not col.startswith("cnn_")]
    first = scores.sort_values([key, "_ensemble_member"], kind="stable").drop_duplicates(key, keep="first")
    grouped = scores.groupby(key, dropna=False)[positive_column]
    agg = grouped.agg(["mean", "std", "max", "count"]).reset_index()
    agg = agg.rename(
        columns={
            "mean": "eb_p_mean",
            "std": "eb_p_std",
            "max": "eb_p_max",
            "count": "eb_ensemble_n",
        }
    )
    out = first.loc[:, [col for col in first_cols if col in first.columns]].merge(agg, on=key, how="left")
    out["eb_p_std"] = out["eb_p_std"].fillna(0.0)
    return out


def evaluate_eb_miner_release_gate(
    training_summary: dict[str, Any],
    *,
    min_members: int = 3,
    min_heldout_recall: float = 0.40,
    min_weighted_average_precision: float = 0.10,
) -> dict[str, Any]:
    """Require non-degenerate grouped held-out performance before review release."""

    member_rows: list[dict[str, Any]] = []
    total_tp = 0
    total_fn = 0
    ap_weighted_sum = 0.0
    ap_weight = 0
    for member_idx, member in enumerate(training_summary.get("member_summaries", [])):
        evaluation = member.get("eb_binary_eval", {})
        member_tp = 0
        member_fn = 0
        member_ap_sum = 0.0
        member_ap_weight = 0
        for split in ("validation", "test"):
            metrics = evaluation.get(split, {})
            n_positive = _safe_int(metrics.get("n_positive"), 0)
            if n_positive <= 0:
                continue
            member_tp += _safe_int(metrics.get("tp"), 0)
            member_fn += _safe_int(metrics.get("fn"), 0)
            ap = _safe_float(metrics.get("average_precision"))
            if np.isfinite(ap):
                member_ap_sum += ap * n_positive
                member_ap_weight += n_positive
        if member_tp + member_fn <= 0 or member_ap_weight <= 0:
            continue
        recall = member_tp / (member_tp + member_fn)
        average_precision = member_ap_sum / member_ap_weight
        member_rows.append(
            {
                "member_index": int(member_idx),
                "n_heldout_positive": int(member_tp + member_fn),
                "heldout_recall": float(recall),
                "heldout_average_precision": float(average_precision),
            }
        )
        total_tp += member_tp
        total_fn += member_fn
        ap_weighted_sum += member_ap_sum
        ap_weight += member_ap_weight
    heldout_recall = total_tp / (total_tp + total_fn) if total_tp + total_fn else float("nan")
    weighted_ap = ap_weighted_sum / ap_weight if ap_weight else float("nan")
    passed = (
        len(member_rows) >= int(min_members)
        and np.isfinite(heldout_recall)
        and heldout_recall >= float(min_heldout_recall)
        and np.isfinite(weighted_ap)
        and weighted_ap >= float(min_weighted_average_precision)
    )
    return {
        "passed": bool(passed),
        "n_evaluable_members": int(len(member_rows)),
        "heldout_recall": float(heldout_recall),
        "heldout_weighted_average_precision": float(weighted_ap),
        "thresholds": {
            "min_members": int(min_members),
            "min_heldout_recall": float(min_heldout_recall),
            "min_weighted_average_precision": float(min_weighted_average_precision),
        },
        "members": member_rows,
    }


def _is_harmonic(period_a: float, period_b: float, *, tolerance: float) -> bool:
    if not np.isfinite(period_a) or not np.isfinite(period_b) or period_a <= 0 or period_b <= 0:
        return False
    ratio = period_a / period_b
    harmonics = (0.5, 1.0 / 3.0, 2.0, 3.0)
    return any(abs(ratio - h) / h <= tolerance for h in harmonics)


def select_eb_priority_queue(
    scored: pd.DataFrame,
    *,
    n_review: int = 500,
    harmonic_probability_floor: float = 0.80,
    harmonic_tolerance: float = 0.02,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {"tic", "period_d", "t0_bjd", "duration_min", "eb_p_mean", "eb_p_std"}
    missing = sorted(required - set(scored.columns))
    if missing:
        raise KeyError(f"scored candidates missing columns: {missing}")
    work = scored.copy()
    work["eb_p_mean"] = pd.to_numeric(work["eb_p_mean"], errors="coerce").fillna(-np.inf)
    work["eb_p_std"] = pd.to_numeric(work["eb_p_std"], errors="coerce").fillna(0.0)
    if "sde_max" not in work and "sde" in work:
        work["sde_max"] = work["sde"]
    work["sde_max"] = pd.to_numeric(work.get("sde_max", np.nan), errors="coerce")
    work = work.sort_values(["eb_p_mean", "eb_p_std", "sde_max"], ascending=[False, False, False], kind="stable")

    selected_rows: list[pd.Series] = []
    selected_by_tic: dict[int, list[float]] = {}
    duplicate_policy: list[str] = []
    for _, row in work.iterrows():
        tic = _safe_int(row.get("tic"))
        period = _safe_float(row.get("period_d"))
        p_eb = _safe_float(row.get("eb_p_mean"))
        existing = selected_by_tic.get(tic, [])
        policy = "top_per_tic"
        if existing:
            allow = (
                p_eb >= harmonic_probability_floor
                and len(existing) < 2
                and any(_is_harmonic(period, old, tolerance=harmonic_tolerance) for old in existing)
            )
            if not allow:
                continue
            policy = "harmonic_secondary_exception"
        selected_rows.append(row)
        duplicate_policy.append(policy)
        selected_by_tic.setdefault(tic, []).append(period)
        if len(selected_rows) >= n_review:
            break
    if len(selected_rows) < n_review:
        raise ValueError(f"selected only {len(selected_rows)} rows; need {n_review}")
    queue = pd.DataFrame(selected_rows).reset_index(drop=True)
    queue["priority_rank"] = np.arange(1, len(queue) + 1, dtype=int)
    queue["eb_miner_duplicate_policy"] = duplicate_policy
    queue["selection_bucket"] = "eb_miner_top_score"
    queue["source_bucket"] = "real_eb_miner_priority"
    queue["source_kind"] = "real_candidate"
    queue["truth_source_kind"] = "real_candidate"
    queue["twirl_vet_sheet_name"] = queue["review_id"].map(_sheet_name)
    queue["twirl_vet_sheet_pdf_name"] = queue["twirl_vet_sheet_name"].str.replace(".png", ".pdf", regex=False)
    queue["candidate_key"] = queue.apply(candidate_key, axis=1)
    for col in LABEL_HEADER:
        if col not in queue:
            queue[col] = ""
    for col in ("label", "label_source", "labeler", "notes", "updated_utc"):
        queue[col] = ""
    summary = {
        "n_review": int(len(queue)),
        "score_min": float(queue["eb_p_mean"].min()),
        "score_median": float(queue["eb_p_mean"].median()),
        "score_max": float(queue["eb_p_mean"].max()),
        "duplicate_policy_counts": queue["eb_miner_duplicate_policy"].value_counts().sort_index().to_dict(),
        "source_pool_counts": queue.get("source_pool", pd.Series(dtype=str)).fillna("").astype(str).value_counts().sort_index().to_dict(),
    }
    return queue, summary


def write_eb_priority_queue(
    *,
    scored: pd.DataFrame,
    out_dir: Path,
    n_review: int = 500,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    queue, selection_summary = select_eb_priority_queue(scored, n_review=n_review)
    queue_path = out_dir / f"review_queue_eb_priority_{int(n_review)}.csv"
    queue.to_csv(queue_path, index=False)
    labels_path = out_dir / "human_labels_vetted.csv"
    if not labels_path.exists():
        pd.DataFrame(columns=LABEL_HEADER).to_csv(labels_path, index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_path),
        "labels_csv": str(labels_path),
        "selection": selection_summary,
        "verification_passed": len(queue) == n_review and queue["source_kind"].eq("real_candidate").all(),
    }
    _write_json(out_dir / "summary.json", summary)
    return summary
