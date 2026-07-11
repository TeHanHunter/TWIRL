"""Recovery50 teacher/student helpers for S56 human-vetting labels."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.label_io import (
    BASE_LABEL_COLUMNS,
    candidate_key as shared_candidate_key,
    latest_label_records,
    normalize_review_queue,
    validate_label_records,
)

from twirl.io.compact_export import read_compact_lc_export, read_injected_lc_group
from twirl.io.hlsp import BJDREFI, HLSPLightCurve, read_hlsp
from twirl.vetting.adp_only import ADP_ONLY_APERTURES, ADP_ONLY_METADATA_COLUMNS
from twirl.vetting.label_schema import BLS_TRUTH_MATCH_MODES
from twirl.vetting.lightcurve_label_app import find_hlsp_path
from twirl.vetting.self_training import (
    FeatureConfig,
    LogisticConfig,
    append_predictions,
    fit_feature_spec,
    save_model,
    select_pseudo_labels,
    train_softmax_model,
)


TEACHER_POSITIVE_LABEL = "planet_like"
TEACHER_NEGATIVE_LABEL = "instrumental_or_systematic"
TRAINING_LABEL_MAP = {
    # Human convention for the S56 teacher queue: "uncertain" means flat or no
    # obvious transit, so it is a negative example for the main visible-signal
    # teacher while the raw label remains preserved for audit tables.
    "uncertain": TEACHER_NEGATIVE_LABEL,
    "no_visible_signal": TEACHER_NEGATIVE_LABEL,
}
EXCLUDED_MAIN_LABELS = frozenset({"", "skip"})
DEFAULT_APERTURES = ADP_ONLY_APERTURES
DEFAULT_SHAPE_BINS = 65
DEFAULT_SHAPE_WINDOW_DURATIONS = 3.0

SCALAR_METADATA_COLUMNS: tuple[str, ...] = ADP_ONLY_METADATA_COLUMNS

LEAKAGE_PREFIXES = (
    "adjudicated_",
    "compact_",
    "display_",
    "ephemeris_",
    "model_",
    "prior_",
    "pre_revisit_",
    "preserve_",
    "refold_",
    "revisit_",
    "truth_",
    "recovery_status",
    "topn_",
    "strict_top1_recovered",
    "any_exact_or_harmonic_recovered",
    "source_",
    "selection_",
    "teacher_",
    "student_",
    "pseudo_",
    "raw_",
    "repeat_",
    "resolution_",
    "human_",
    "queue_",
    "morphology_",
    "harmonic_",
    "note_period_",
    "effective_period_",
    "variable_period_",
    "broad_",
)
LEAKAGE_EXACT = frozenset(
    {
        "label",
        "label_source",
        "labeler",
        "notes",
        "updated_utc",
        "candidate_key",
        "human_candidate_key",
        "injection_id",
        "h5_group",
        "source_h5",
        "twirl_vet_sheet_name",
        "twirl_vet_sheet_pdf_name",
        "twirl_vet_input",
        "bls_truth_match",
        "truth_signal_present",
        "is_injected_row",
        "is_labeled",
        "main_teacher_include",
        "audit_include",
        "harmonic_suspect",
        "training_split",
        "period_factor",
        "period_status",
        "period_factor_source",
        "period_task",
        "model_target_policy_version",
        "is_repeat",
        "origin_queue",
        "labeling_era",
        "review_id",
        "tic",
        "sector",
        "cam",
        "ccd",
    }
)


@dataclass(frozen=True)
class LabelPolicy:
    """Controls for the first recovery50 teacher target."""

    min_multiclass_count: int = 40
    include_rare_classes: bool = True


def json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".json", ".jsonl"}:
        return pd.read_json(path, lines=suffix == ".jsonl")
    raise ValueError(f"unsupported table format: {path}")


def write_table(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".csv":
        df.to_csv(path, index=False)
        return path
    try:
        df.to_parquet(path, compression="zstd", index=False)
        return path
    except (ImportError, ValueError, ModuleNotFoundError):
        csv_path = path.with_suffix(".csv")
        df.to_csv(csv_path, index=False)
        return csv_path


def latest_labels(labels_csv: Path) -> pd.DataFrame:
    """Return the latest browser label per row_id."""

    return latest_label_records(Path(labels_csv), columns=BASE_LABEL_COLUMNS)


def candidate_key(row: pd.Series) -> str:
    return shared_candidate_key(row)


def join_queue_labels(queue_csv: Path, labels_csv: Path) -> pd.DataFrame:
    """Join a review queue to latest browser labels without mutating either."""

    queue = normalize_review_queue(read_table(queue_csv))
    for col in ("label", "label_source", "labeler", "notes", "updated_utc"):
        if col in queue:
            queue = queue.rename(columns={col: f"queue_{col}"})

    raw_labels = latest_labels(labels_csv)
    validate_label_records(queue, raw_labels)
    labels = raw_labels.rename(
        columns={
            "candidate_key": "human_candidate_key",
            "label": "human_label",
            "label_source": "human_label_source",
            "labeler": "human_labeler",
            "notes": "human_notes",
            "updated_utc": "human_updated_utc",
        }
    )
    labels = labels.drop(columns=[c for c in ("tic", "sector") if c in labels], errors="ignore")
    merged = queue.merge(labels, on="row_id", how="left", validate="one_to_one")
    for col in ("human_candidate_key", "human_label", "human_label_source", "human_labeler", "human_notes", "human_updated_utc"):
        if col not in merged:
            merged[col] = ""
        merged[col] = merged[col].fillna("").astype(str)
    merged["candidate_key_matches_label"] = (
        merged["human_candidate_key"].eq("")
        | merged["candidate_key"].fillna("").astype(str).eq(merged["human_candidate_key"])
    )
    return add_label_roles(merged)


def _as_bool(series: pd.Series | bool, index: pd.Index) -> pd.Series:
    if isinstance(series, bool):
        return pd.Series(series, index=index)
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y"})


def add_label_roles(df: pd.DataFrame, policy: LabelPolicy | None = None) -> pd.DataFrame:
    """Attach explicit human-target and auxiliary truth roles."""

    policy = policy or LabelPolicy()
    out = df.copy()
    label = out.get("human_label", pd.Series("", index=out.index)).fillna("").astype(str)
    training_label = label.replace(TRAINING_LABEL_MAP)
    source_kind = out.get("source_kind", pd.Series("", index=out.index)).fillna("").astype(str)
    topn = out.get("topn_recovery_status", out.get("recovery_status", pd.Series("", index=out.index))).fillna("").astype(str)

    strong_counts = training_label[~label.isin(EXCLUDED_MAIN_LABELS)].value_counts()
    allowed = {TEACHER_POSITIVE_LABEL, TEACHER_NEGATIVE_LABEL}
    if policy.include_rare_classes:
        for rare_label, count in strong_counts.items():
            if rare_label not in allowed and int(count) >= policy.min_multiclass_count:
                allowed.add(str(rare_label))

    out["is_labeled"] = label.ne("")
    out["is_injected_row"] = source_kind.eq("injection_recovery")
    out["truth_signal_present"] = out["is_injected_row"]
    bls_truth_match = topn.isin(BLS_TRUTH_MATCH_MODES)
    for column in out.columns:
        if column.startswith(("topn_exact_recovered_", "topn_harmonic_match_")):
            bls_truth_match |= _as_bool(out[column], out.index)
    if "strict_top1_recovered" in out:
        bls_truth_match |= _as_bool(out["strict_top1_recovered"], out.index)
    out["bls_truth_match"] = bls_truth_match
    out["main_teacher_target"] = np.where(
        (~label.isin(EXCLUDED_MAIN_LABELS)) & training_label.isin(allowed),
        training_label,
        "",
    )
    out["main_teacher_include"] = out["main_teacher_target"].astype(str).ne("")
    out["audit_include"] = label.ne("") & ~label.eq("skip")
    out["label_policy"] = "human_visible_signal_main_truth_aux"
    out["allowed_teacher_classes"] = ",".join(sorted(allowed))
    return out


def add_deterministic_splits(
    df: pd.DataFrame,
    *,
    validation_fraction: float = 0.20,
    test_fraction: float = 0.20,
    random_state: int = 56,
) -> pd.DataFrame:
    """Add stratified train/validation/test splits for main teacher rows."""

    out = df.copy()
    split = pd.Series("unlabeled_or_audit", index=out.index, dtype=object)
    eligible = out["main_teacher_include"].fillna(False).astype(bool)
    if not eligible.any():
        out["training_split"] = split
        return out

    rng = np.random.default_rng(random_state)
    keys = (
        out.loc[eligible, "main_teacher_target"].fillna("").astype(str)
        + "|"
        + out.loc[eligible, "source_kind"].fillna("").astype(str)
    )
    for _, idx_values in keys.groupby(keys).groups.items():
        idx = np.asarray(list(idx_values), dtype=int)
        rng.shuffle(idx)
        n = len(idx)
        n_test = int(round(test_fraction * n))
        n_val = int(round(validation_fraction * n))
        if n >= 3 and test_fraction > 0:
            n_test = max(1, n_test)
        if n - n_test >= 3 and validation_fraction > 0:
            n_val = max(1, n_val)
        n_test = min(n_test, max(0, n - 1))
        n_val = min(n_val, max(0, n - n_test - 1))
        split.loc[idx] = "train"
        if n_test:
            split.loc[idx[:n_test]] = "test"
        if n_val:
            split.loc[idx[n_test:n_test + n_val]] = "validation"
    out["training_split"] = split
    return out


def leakage_columns(columns: Iterable[str]) -> list[str]:
    leaks: list[str] = []
    for col in columns:
        if col in LEAKAGE_EXACT or any(col.startswith(prefix) for prefix in LEAKAGE_PREFIXES):
            leaks.append(col)
    return sorted(set(leaks))


def metadata_feature_columns(df: pd.DataFrame) -> list[str]:
    """Return an allowlisted, leakage-checked scalar feature list."""

    columns = [col for col in SCALAR_METADATA_COLUMNS if col in df.columns]
    columns = [col for col in columns if not col.endswith("_feature_disabled")]
    leaks = leakage_columns(columns)
    if leaks:
        raise ValueError(f"metadata feature allowlist contains leakage columns: {leaks}")
    return columns


def shape_feature_columns(df: pd.DataFrame) -> list[str]:
    columns = [
        col
        for col in df.columns
        if col.startswith("shape_flux_") or col.startswith("shape_mask_") or col.startswith("shape_count_")
    ]
    leaks = leakage_columns(columns)
    if leaks:
        raise ValueError(f"shape feature columns contain leakage columns: {leaks}")
    return sorted(columns)


def _time_to_bjd(time: np.ndarray) -> np.ndarray:
    finite = time[np.isfinite(time)]
    if finite.size and np.nanmedian(finite) < 1.0e5:
        return time + BJDREFI
    return time


def _safe_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


_HALF_PERIOD_RE = re.compile(
    r"(?ix)"
    r"("
    r"\bp\s*/\s*2\b"
    r"|refold\s+(?:at|to|on)\s+p\s*/\s*2"
    r"|half[-\s]*period"
    r"|half\s+of\s+(?:the\s+)?period"
    r")"
)
_DOUBLE_PERIOD_RE = re.compile(
    r"(?ix)"
    r"("
    r"\b2\s*\*?\s*p\b"
    r"|\b2p\b"
    r"|double[-\s]*period"
    r"|twice\s+(?:the\s+)?period"
    r")"
)
_HARMONIC_SUSPECT_RE = re.compile(r"(?i)\b(possible\s+)?harmonic\b|\bwrong\s+period\b|\bwrong\s+fold\b")


def infer_harmonic_refold(notes: Any) -> tuple[bool, float, str]:
    """Infer a simple period refold factor from reviewer notes.

    The return value is ``(harmonic_suspect, refold_factor, reason)``. A finite
    factor means the row can be folded automatically for tensor generation.
    """

    text = "" if notes is None or pd.isna(notes) else str(notes).strip()
    if not text:
        return False, float("nan"), ""
    if _HALF_PERIOD_RE.search(text):
        return True, 0.5, "half_period_note"
    if _DOUBLE_PERIOD_RE.search(text):
        return True, 2.0, "double_period_note"
    if _HARMONIC_SUSPECT_RE.search(text):
        return True, float("nan"), "harmonic_note"
    return False, float("nan"), ""


def select_model_ephemeris(row: pd.Series) -> tuple[float, float, float, str]:
    """Return the ephemeris a model should fold on, with provenance."""

    candidates = (
        ("model_", row.get("model_ephemeris_source", "model")),
        ("display_", row.get("display_ephemeris_source", "display")),
        ("anchor_", "twirl_vet_anchor"),
        ("", "queue"),
    )
    for prefix, source in candidates:
        period_col = f"{prefix}period_d" if prefix else "period_d"
        t0_col = f"{prefix}t0_bjd" if prefix else "t0_bjd"
        duration_col = f"{prefix}duration_min" if prefix else "duration_min"
        period_d = _safe_float(row.get(period_col))
        t0_bjd = _safe_float(row.get(t0_col))
        duration_min = _safe_float(row.get(duration_col))
        if np.isfinite(period_d) and np.isfinite(t0_bjd) and np.isfinite(duration_min) and period_d > 0 and duration_min > 0:
            source_text = str(source or "model").strip() or "model"
            return period_d, t0_bjd, duration_min, source_text
    return float("nan"), float("nan"), float("nan"), "invalid"


def select_bls_ephemeris(row: pd.Series) -> tuple[float, float, float, str]:
    """Return the displayed/BLS ephemeris before human harmonic correction."""

    candidates = (
        ("display_", row.get("display_ephemeris_source", "display")),
        ("anchor_", "twirl_vet_anchor"),
        ("", "queue"),
    )
    for prefix, source in candidates:
        period_col = f"{prefix}period_d" if prefix else "period_d"
        t0_col = f"{prefix}t0_bjd" if prefix else "t0_bjd"
        duration_col = f"{prefix}duration_min" if prefix else "duration_min"
        period_d = _safe_float(row.get(period_col))
        t0_bjd = _safe_float(row.get(t0_col))
        duration_min = _safe_float(row.get(duration_col))
        if np.isfinite(period_d) and np.isfinite(t0_bjd) and np.isfinite(duration_min) and period_d > 0 and duration_min > 0:
            source_text = str(source or "bls").strip() or "bls"
            return period_d, t0_bjd, duration_min, source_text
    return float("nan"), float("nan"), float("nan"), "invalid"


def add_harmonic_ephemeris_annotations(
    df: pd.DataFrame,
    *,
    notes_column: str = "human_notes",
) -> pd.DataFrame:
    """Add reviewer-note harmonic flags and model-fold ephemeris overrides."""

    out = df.copy()
    for col in ("model_period_d", "model_t0_bjd", "model_duration_min"):
        if col not in out:
            out[col] = np.nan
    if "model_ephemeris_source" not in out:
        out["model_ephemeris_source"] = ""

    harmonic_suspect: list[bool] = []
    refold_factor: list[float] = []
    refold_period: list[float] = []
    refold_t0: list[float] = []
    refold_duration: list[float] = []
    ephemeris_status: list[str] = []
    compact_target: list[str] = []
    compact_include: list[bool] = []
    preserve_target: list[str] = []

    for idx, row in out.iterrows():
        notes = row.get(notes_column, row.get("human_notes", ""))
        suspect, factor, reason = infer_harmonic_refold(notes)
        bls_period, bls_t0, bls_duration, bls_source = select_bls_ephemeris(row)
        period = bls_period * factor if suspect and np.isfinite(factor) and np.isfinite(bls_period) else float("nan")
        harmonic_suspect.append(bool(suspect))
        refold_factor.append(float(factor))
        refold_period.append(float(period))
        refold_t0.append(float(bls_t0) if suspect and np.isfinite(period) else float("nan"))
        refold_duration.append(float(bls_duration) if suspect and np.isfinite(period) else float("nan"))

        if suspect and np.isfinite(period):
            out.at[idx, "model_period_d"] = period
            out.at[idx, "model_t0_bjd"] = bls_t0
            out.at[idx, "model_duration_min"] = bls_duration
            out.at[idx, "model_ephemeris_source"] = f"human_{reason}:{bls_source}"
            ephemeris_status.append("human_refold_ready")
        elif suspect:
            ephemeris_status.append("human_harmonic_suspect_unresolved")
        else:
            ephemeris_status.append("display_or_queue")

        label = str(row.get("human_label", "") or "")
        if label in {"planet_like", "wide_transit_like", "eclipsing_binary_or_pceb", "stellar_variability"}:
            preserve_target.append("preserve_signal")
        elif label in {"instrumental_or_systematic", "uncertain", "no_visible_signal"}:
            preserve_target.append("reject_signal")
        else:
            preserve_target.append("")

        if label == "planet_like" and (not suspect or np.isfinite(period)):
            compact_target.append("planet_like")
            compact_include.append(True)
        elif label in {"instrumental_or_systematic", "uncertain", "no_visible_signal", "wide_transit_like", "eclipsing_binary_or_pceb", "stellar_variability"}:
            compact_target.append(label)
            compact_include.append(label != "uncertain")
        else:
            compact_target.append("")
            compact_include.append(False)

    out["harmonic_suspect"] = harmonic_suspect
    out["refold_factor"] = refold_factor
    out["refold_period_d"] = refold_period
    out["refold_t0_bjd"] = refold_t0
    out["refold_duration_min"] = refold_duration
    out["ephemeris_status"] = ephemeris_status
    out["preserve_signal_target"] = preserve_target
    out["compact_morphology_target"] = compact_target
    out["compact_morphology_include"] = compact_include
    return out


def add_display_ephemeris(df: pd.DataFrame, metrics_tables: Sequence[Path] = ()) -> pd.DataFrame:
    """Attach the ephemeris used in vet-sheet displays, falling back to queue timing."""

    out = df.copy()
    index = out.index
    for col in ("period_d", "t0_bjd", "duration_min"):
        out[f"queue_{col}"] = pd.to_numeric(out.get(col, pd.Series(np.nan, index=index)), errors="coerce")

    out["display_period_d"] = out["queue_period_d"]
    out["display_t0_bjd"] = out["queue_t0_bjd"]
    out["display_duration_min"] = out["queue_duration_min"]
    out["display_ephemeris_source"] = "queue"
    out["display_anchor_aperture"] = ""
    out["display_anchor_sde"] = np.nan
    out["display_vet_sheet_name"] = ""
    out["display_vet_sheet_pdf_name"] = ""
    out["display_vet_status"] = ""

    metric_frames: list[pd.DataFrame] = []
    wanted = {
        "review_id",
        "anchor_period_d",
        "anchor_t0_bjd",
        "anchor_duration_min",
        "anchor_aperture",
        "anchor_sde",
        "twirl_vet_sheet_name",
        "twirl_vet_sheet_pdf_name",
        "twirl_vet_status",
    }
    for path in metrics_tables:
        if not path.exists():
            continue
        metrics = read_table(path)
        if "review_id" not in metrics:
            continue
        keep = [col for col in metrics.columns if col in wanted]
        metric = metrics.loc[:, keep].copy()
        for col in wanted:
            if col not in metric:
                metric[col] = np.nan if col.endswith(("_d", "_min", "_sde")) else ""
        metric_frames.append(metric)

    if metric_frames:
        metrics = pd.concat(metric_frames, ignore_index=True)
        metrics = metrics.dropna(subset=["review_id"]).copy()
        metrics["review_id"] = metrics["review_id"].astype(str)
        metrics = metrics.drop_duplicates("review_id", keep="last")
        metrics = metrics.rename(
            columns={
                "anchor_period_d": "_metric_display_period_d",
                "anchor_t0_bjd": "_metric_display_t0_bjd",
                "anchor_duration_min": "_metric_display_duration_min",
                "anchor_aperture": "_metric_display_anchor_aperture",
                "anchor_sde": "_metric_display_anchor_sde",
                "twirl_vet_sheet_name": "_metric_display_vet_sheet_name",
                "twirl_vet_sheet_pdf_name": "_metric_display_vet_sheet_pdf_name",
                "twirl_vet_status": "_metric_display_vet_status",
            }
        )
        out["review_id"] = out["review_id"].astype(str)
        out = out.merge(metrics, on="review_id", how="left", validate="many_to_one")
        metric_period = pd.to_numeric(out["_metric_display_period_d"], errors="coerce")
        metric_t0 = pd.to_numeric(out["_metric_display_t0_bjd"], errors="coerce")
        metric_duration = pd.to_numeric(out["_metric_display_duration_min"], errors="coerce")
        use_anchor = metric_period.gt(0) & metric_t0.notna() & metric_duration.gt(0)
        out.loc[use_anchor, "display_period_d"] = metric_period[use_anchor]
        out.loc[use_anchor, "display_t0_bjd"] = metric_t0[use_anchor]
        out.loc[use_anchor, "display_duration_min"] = metric_duration[use_anchor]
        out.loc[use_anchor, "display_ephemeris_source"] = "twirl_vet_anchor"
        for src, dst in (
            ("_metric_display_anchor_aperture", "display_anchor_aperture"),
            ("_metric_display_anchor_sde", "display_anchor_sde"),
            ("_metric_display_vet_sheet_name", "display_vet_sheet_name"),
            ("_metric_display_vet_sheet_pdf_name", "display_vet_sheet_pdf_name"),
            ("_metric_display_vet_status", "display_vet_status"),
        ):
            out.loc[use_anchor, dst] = out.loc[use_anchor, src]
        out = out.drop(columns=[col for col in out.columns if col.startswith("_metric_display_")])

    denom = out["queue_period_d"].abs()
    out["display_ephemeris_period_rel_delta"] = np.where(
        denom > 0,
        (pd.to_numeric(out["display_period_d"], errors="coerce") - out["queue_period_d"]).abs() / denom,
        np.nan,
    )
    out["display_ephemeris_used_anchor"] = out["display_ephemeris_source"].eq("twirl_vet_anchor")
    return out


def _clean_aperture_name(aperture: str) -> str:
    return aperture.lower().replace("det_flux_", "").replace("det_flux", "flux")


def bin_folded_channel(
    *,
    time: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    n_bins: int = DEFAULT_SHAPE_BINS,
    window_durations: float = DEFAULT_SHAPE_WINDOW_DURATIONS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Fold and median-bin one flux channel in units of transit duration."""

    if period_d <= 0 or duration_min <= 0 or n_bins <= 0:
        values = np.full(n_bins, np.nan, dtype=np.float64)
        return values, np.zeros(n_bins, dtype=bool), np.zeros(n_bins, dtype=np.int16), {"status": "invalid_ephemeris"}

    time_bjd = _time_to_bjd(np.asarray(time, dtype=np.float64))
    flux = np.asarray(flux, dtype=np.float64)
    quality = np.asarray(quality, dtype=np.int32)
    duration_hr = duration_min / 60.0
    phase_hr = (((time_bjd - t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d) * 24.0
    x = phase_hr / duration_hr
    good = (quality == 0) & np.isfinite(x) & np.isfinite(flux)
    if not np.any(good):
        values = np.full(n_bins, np.nan, dtype=np.float64)
        return values, np.zeros(n_bins, dtype=bool), np.zeros(n_bins, dtype=np.int16), {"status": "no_good_cadences"}

    x_good = x[good]
    flux_good = flux[good]
    oot = np.abs(x_good) > 1.5
    baseline = float(np.nanmedian(flux_good[oot])) if np.any(oot) else float(np.nanmedian(flux_good))
    if not np.isfinite(baseline) or abs(baseline) < 1.0e-8:
        baseline = 1.0
    y = flux_good / baseline - 1.0

    centers = np.linspace(-window_durations, window_durations, n_bins)
    step = float(centers[1] - centers[0]) if n_bins > 1 else float(window_durations)
    edges = np.concatenate(
        ([centers[0] - 0.5 * step], 0.5 * (centers[1:] + centers[:-1]), [centers[-1] + 0.5 * step])
    )
    bin_id = np.digitize(x_good, edges) - 1
    values = np.full(n_bins, np.nan, dtype=np.float64)
    counts = np.zeros(n_bins, dtype=np.int16)
    for idx in range(n_bins):
        samples = y[bin_id == idx]
        samples = samples[np.isfinite(samples)]
        if samples.size:
            values[idx] = float(np.nanmedian(samples))
            counts[idx] = min(int(samples.size), np.iinfo(np.int16).max)
    mask = np.isfinite(values)
    in_transit = np.abs(x_good) <= 0.5
    near_oot = (np.abs(x_good) >= 1.5) & (np.abs(x_good) <= window_durations)
    in_median = float(np.nanmedian(y[in_transit])) if np.any(in_transit) else float("nan")
    oot_median = float(np.nanmedian(y[near_oot])) if np.any(near_oot) else float("nan")
    stats = {
        "status": "ok",
        "finite_bin_frac": float(mask.mean()),
        "n_good": int(np.sum(good)),
        "n_in_window": int(np.sum(np.abs(x_good) <= window_durations)),
        "n_in_transit": int(np.sum(in_transit)),
        "folded_in_median": in_median,
        "folded_oot_median": oot_median,
        "folded_depth": float(oot_median - in_median) if np.isfinite(in_median) and np.isfinite(oot_median) else float("nan"),
    }
    return values, mask, counts, stats


def _read_lc_for_row(
    row: pd.Series,
    *,
    apertures: Sequence[str],
    compact_lc_h5: Path | None,
    hlsp_root: Path | None,
    injection_h5_override: Path | None = None,
) -> HLSPLightCurve | None:
    source_kind = str(row.get("source_kind", ""))
    if source_kind == "injection_recovery":
        source_h5 = str(injection_h5_override or row.get("source_h5", "") or "")
        group = str(row.get("h5_group", "") or "")
        if source_h5 and group:
            return read_injected_lc_group(Path(source_h5), group_path=group, columns=apertures)
        return None
    tic = int(row.get("tic"))
    if compact_lc_h5 is not None:
        lc = read_compact_lc_export(compact_lc_h5, tic=tic, columns=apertures)
        if lc is not None:
            return lc
    if hlsp_root is not None:
        path = find_hlsp_path(hlsp_root, tic, int(row.get("sector", 56)))
        if path is not None:
            return read_hlsp(path, columns=apertures)
    return None


def build_folded_shape_features(
    *,
    queue_csv: Path,
    out_dir: Path,
    compact_lc_h5: Path | None,
    hlsp_root: Path | None,
    injection_h5_override: Path | None = None,
    apertures: Sequence[str] = DEFAULT_APERTURES,
    n_bins: int = DEFAULT_SHAPE_BINS,
    window_durations: float = DEFAULT_SHAPE_WINDOW_DURATIONS,
    max_rows: int | None = None,
    progress_every: int = 100,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    queue = read_table(queue_csv)
    if "row_id" not in queue:
        queue.insert(0, "row_id", np.arange(len(queue), dtype=int))
    if max_rows is not None:
        queue = queue.head(max_rows).copy()

    records: list[dict[str, Any]] = []
    status_counts: dict[str, int] = {}
    for n_done, (_, row) in enumerate(queue.iterrows(), start=1):
        rec: dict[str, Any] = {
            "row_id": int(row["row_id"]),
            "review_id": str(row.get("review_id", "")),
            "tic": int(row.get("tic", -1)),
            "source_kind": str(row.get("source_kind", "")),
        }
        period_d, t0_bjd, duration_min, ephemeris_source = select_model_ephemeris(row)
        rec["model_period_d"] = period_d
        rec["model_t0_bjd"] = t0_bjd
        rec["model_duration_min"] = duration_min
        rec["model_ephemeris_source"] = ephemeris_source
        rec["queue_period_d"] = _safe_float(row.get("period_d"))
        rec["queue_t0_bjd"] = _safe_float(row.get("t0_bjd"))
        rec["queue_duration_min"] = _safe_float(row.get("duration_min"))
        lc = _read_lc_for_row(
            row,
            apertures=apertures,
            compact_lc_h5=compact_lc_h5,
            hlsp_root=hlsp_root,
            injection_h5_override=injection_h5_override,
        )
        if lc is None:
            rec["shape_status"] = "missing_light_curve"
            records.append(rec)
            status_counts["missing_light_curve"] = status_counts.get("missing_light_curve", 0) + 1
            continue
        aperture_statuses = []
        for aperture in apertures:
            prefix = _clean_aperture_name(aperture)
            if aperture not in lc.flux:
                aperture_statuses.append(f"{aperture}:missing_aperture")
                continue
            values, mask, counts, stats = bin_folded_channel(
                time=lc.time,
                flux=lc.flux[aperture],
                quality=lc.quality,
                period_d=period_d,
                t0_bjd=t0_bjd,
                duration_min=duration_min,
                n_bins=n_bins,
                window_durations=window_durations,
            )
            aperture_statuses.append(f"{aperture}:{stats['status']}")
            for idx, value in enumerate(values):
                rec[f"shape_flux_{prefix}_bin_{idx:03d}"] = value
                rec[f"shape_mask_{prefix}_bin_{idx:03d}"] = float(mask[idx])
                rec[f"shape_count_{prefix}_bin_{idx:03d}"] = int(counts[idx])
            for key, value in stats.items():
                if key != "status":
                    rec[f"shape_{prefix}_{key}"] = value
        status = "ok" if all(part.endswith(":ok") for part in aperture_statuses) else ";".join(aperture_statuses)
        rec["shape_status"] = status
        records.append(rec)
        status_counts[status] = status_counts.get(status, 0) + 1
        if progress_every > 0 and n_done % progress_every == 0:
            print(f"[shape] processed {n_done:,}/{len(queue):,}", flush=True)

    features = pd.DataFrame(records)
    out_csv = out_dir / "folded_shape_features.csv"
    features.to_csv(out_csv, index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_csv),
        "compact_lc_h5": str(compact_lc_h5) if compact_lc_h5 is not None else "",
        "hlsp_root": str(hlsp_root) if hlsp_root is not None else "",
        "injection_h5_override": str(injection_h5_override) if injection_h5_override is not None else "",
        "out_dir": str(out_dir),
        "n_rows": int(len(features)),
        "apertures": list(apertures),
        "n_bins": int(n_bins),
        "window_durations": float(window_durations),
        "status_counts": {str(k): int(v) for k, v in sorted(status_counts.items())},
        "shape_feature_columns": shape_feature_columns(features),
        "outputs": {"features_csv": str(out_csv)},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    return summary


def _crosstab(df: pd.DataFrame, row: str, col: str, out_path: Path) -> dict[str, dict[str, int]]:
    if row not in df or col not in df:
        return {}
    table = pd.crosstab(df[row].fillna("").astype(str), df[col].fillna("").astype(str))
    table.to_csv(out_path)
    return {
        str(idx): {str(k): int(v) for k, v in vals.items()}
        for idx, vals in table.to_dict(orient="index").items()
    }


def _bin_table(
    df: pd.DataFrame,
    *,
    value_col: str,
    bins: Sequence[float],
    labels: Sequence[str],
) -> pd.DataFrame:
    if value_col not in df or df.empty:
        return pd.DataFrame()
    work = df.copy()
    work["bin"] = pd.cut(pd.to_numeric(work[value_col], errors="coerce"), bins=bins, labels=labels)
    grouped = work.groupby("bin", observed=False)
    rows = []
    for bin_name, part in grouped:
        if len(part) == 0:
            rows.append({"bin": str(bin_name), "n": 0, "n_planet_like": 0, "frac_planet_like": np.nan, "n_certain": 0, "frac_planet_like_certain": np.nan})
            continue
        certain = part[~part["human_label"].isin(["uncertain", "skip"])]
        rows.append(
            {
                "bin": str(bin_name),
                "n": int(len(part)),
                "n_planet_like": int(part["human_label"].eq(TEACHER_POSITIVE_LABEL).sum()),
                "frac_planet_like": float(part["human_label"].eq(TEACHER_POSITIVE_LABEL).mean()),
                "n_certain": int(len(certain)),
                "frac_planet_like_certain": float(certain["human_label"].eq(TEACHER_POSITIVE_LABEL).mean()) if len(certain) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def audit_recovery50_labels(queue_csv: Path, labels_csv: Path, out_dir: Path) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    joined = join_queue_labels(queue_csv, labels_csv)
    joined = add_deterministic_splits(joined)
    joined.to_csv(out_dir / "joined_human_labels.csv", index=False)

    labeled = joined[joined["human_label"].ne("")].copy()
    injected = labeled[labeled["source_kind"].eq("injection_recovery")].copy()
    real = labeled[labeled["source_kind"].eq("real_candidate")].copy()
    certain_injected = injected[~injected["human_label"].isin(["uncertain", "skip"])].copy()

    tables = {
        "label_by_source_kind": _crosstab(labeled, "source_kind", "human_label", out_dir / "label_by_source_kind.csv"),
        "label_by_bls_truth_match": _crosstab(labeled, "bls_truth_match", "human_label", out_dir / "label_by_bls_truth_match.csv"),
        "label_by_topn_recovery_status": _crosstab(labeled, "topn_recovery_status", "human_label", out_dir / "label_by_topn_recovery_status.csv"),
        "label_by_strict_top1_recovered": _crosstab(labeled, "strict_top1_recovered", "human_label", out_dir / "label_by_strict_top1_recovered.csv"),
    }
    if len(labeled):
        median_cols = [
            col
            for col in (
                "period_d",
                "duration_min",
                "depth",
                "depth_snr",
                "sde_max",
                "tmag",
                "truth_period_d",
                "truth_radius_rearth",
                "truth_model_depth",
                "truth_depth",
                "min_twoap_sde",
            )
            if col in labeled
        ]
        if median_cols:
            labeled.groupby(["source_kind", "human_label"], observed=False)[median_cols].median(numeric_only=True).to_csv(
                out_dir / "label_median_parameter_stats.csv"
            )

    trend_specs = {
        "truth_period_d": ((0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, np.inf), ("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-5", "5-10", ">10")),
        "truth_radius_rearth": ((0.0, 1.0, 2.0, 4.0, 8.0, 12.0, np.inf), ("<1", "1-2", "2-4", "4-8", "8-12", ">12")),
        "tmag": ((-np.inf, 16.0, 17.0, 18.0, 19.0, 20.0, np.inf), ("<16", "16-17", "17-18", "18-19", "19-20", ">20")),
        "truth_model_depth": ((0.0, 0.01, 0.03, 0.1, 0.3, 1.0, np.inf), ("<1%", "1-3%", "3-10%", "10-30%", "30-100%", ">100%")),
        "sde_max": ((-np.inf, 10.0, 20.0, 40.0, 80.0, np.inf), ("<10", "10-20", "20-40", "40-80", ">80")),
    }
    trend_outputs: dict[str, str] = {}
    for col, (bins, labels) in trend_specs.items():
        table = _bin_table(injected, value_col=col, bins=bins, labels=labels)
        if not table.empty:
            path = out_dir / f"injected_human_recovery_by_{col}.csv"
            table.to_csv(path, index=False)
            trend_outputs[col] = str(path)

    disagreement = injected[~injected["human_label"].isin([TEACHER_POSITIVE_LABEL, "uncertain", "skip"])].copy()
    disagreement.to_csv(out_dir / "injected_truth_human_disagreement_rows.csv", index=False)

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_csv),
        "labels_csv": str(labels_csv),
        "out_dir": str(out_dir),
        "n_queue_rows": int(len(joined)),
        "n_labeled": int(len(labeled)),
        "n_unlabeled": int(len(joined) - len(labeled)),
        "label_counts": {str(k): int(v) for k, v in labeled["human_label"].value_counts().sort_index().items()},
        "label_by_source_kind": tables["label_by_source_kind"],
        "n_injected_labeled": int(len(injected)),
        "n_real_labeled": int(len(real)),
        "injected_planet_like_fraction_all": float(injected["human_label"].eq(TEACHER_POSITIVE_LABEL).mean()) if len(injected) else np.nan,
        "injected_planet_like_fraction_certain": float(certain_injected["human_label"].eq(TEACHER_POSITIVE_LABEL).mean()) if len(certain_injected) else np.nan,
        "n_injected_certain": int(len(certain_injected)),
        "n_injected_truth_human_disagreements": int(len(disagreement)),
        "tables": tables,
        "trend_outputs": trend_outputs,
        "outputs": {
            "joined": str(out_dir / "joined_human_labels.csv"),
            "disagreement_rows": str(out_dir / "injected_truth_human_disagreement_rows.csv"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    (out_dir / "summary.md").write_text(render_audit_markdown(summary))
    return summary


def render_audit_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# S56 Recovery50 Human Label Audit",
        "",
        f"- labeled rows: `{summary['n_labeled']}` / `{summary['n_queue_rows']}`",
        f"- injected labeled: `{summary['n_injected_labeled']}`",
        f"- real labeled: `{summary['n_real_labeled']}`",
        f"- injected planet-like fraction: `{summary['injected_planet_like_fraction_all']:.3f}`",
        f"- injected planet-like fraction excluding uncertain/skip: `{summary['injected_planet_like_fraction_certain']:.3f}`",
        f"- injected truth/human disagreements: `{summary['n_injected_truth_human_disagreements']}`",
        "",
        "## Label Counts",
        "",
        "| Label | Count |",
        "|---|---:|",
    ]
    for label, count in summary["label_counts"].items():
        lines.append(f"| `{label}` | {count} |")
    lines.append("")
    return "\n".join(lines)


def load_feature_table(
    *,
    training_table: Path,
    shape_features: Path | None = None,
    metrics_tables: Sequence[Path] = (),
) -> pd.DataFrame:
    df = read_table(training_table).copy()
    metric_frames = [
        read_table(Path(path))
        for path in metrics_tables
        if path is not None and Path(path).exists()
    ]
    if metric_frames:
        metrics = pd.concat(metric_frames, ignore_index=True, sort=False)
        keep = [col for col in metrics.columns if col not in df.columns or col in {"row_id", "review_id"}]
        join_cols = [col for col in ("row_id", "review_id") if col in df.columns and col in metrics.columns]
        if join_cols:
            metrics = metrics.drop_duplicates(subset=join_cols, keep="last")
            df = df.merge(metrics.loc[:, keep], on=join_cols, how="left")
    if shape_features is not None and Path(shape_features).exists():
        shape = read_table(Path(shape_features))
        join_cols = [col for col in ("row_id", "review_id") if col in df.columns and col in shape.columns]
        if not join_cols:
            raise KeyError("shape feature table has no join columns in common with training table")
        keep = [col for col in shape.columns if col not in df.columns or col in join_cols]
        df = df.merge(shape.loc[:, keep], on=join_cols, how="left")
    return df


def select_teacher_classes(df: pd.DataFrame, min_class_count: int) -> list[str]:
    include = _as_bool(df["main_teacher_include"], df.index) if "main_teacher_include" in df else pd.Series(False, index=df.index)
    target = df.loc[include, "main_teacher_target"].fillna("").astype(str)
    counts = target.value_counts()
    classes = [TEACHER_POSITIVE_LABEL, TEACHER_NEGATIVE_LABEL]
    classes.extend(sorted(cls for cls, count in counts.items() if cls not in classes and int(count) >= min_class_count))
    return [cls for cls in classes if int(counts.get(cls, 0)) > 0]


def _split_eval_rows(df: pd.DataFrame, classes: Sequence[str]) -> pd.DataFrame:
    work = df[df["main_teacher_target"].isin(classes)].copy()
    work = work[work["training_split"].isin(["train", "validation", "test"])].copy()
    return work


def _confusion_metrics(truth: pd.Series, pred: pd.Series, classes: Sequence[str]) -> dict[str, Any]:
    truth = truth.fillna("").astype(str)
    pred = pred.fillna("").astype(str)
    matrix = pd.crosstab(
        pd.Categorical(truth, categories=list(classes)),
        pd.Categorical(pred, categories=list(classes)),
        dropna=False,
    )
    matrix.index = list(classes)
    matrix.columns = list(classes)
    per_class: dict[str, dict[str, float]] = {}
    recalls = []
    for cls in classes:
        tp = float(matrix.loc[cls, cls])
        actual = float(matrix.loc[cls, :].sum())
        predicted = float(matrix.loc[:, cls].sum())
        recall = tp / actual if actual else np.nan
        precision = tp / predicted if predicted else np.nan
        per_class[cls] = {"precision": precision, "recall": recall, "support": actual, "predicted": predicted}
        if np.isfinite(recall):
            recalls.append(recall)
    accuracy = float(np.mean(truth.to_numpy() == pred.to_numpy())) if len(truth) else np.nan
    return {
        "n": int(len(truth)),
        "accuracy": accuracy,
        "balanced_accuracy": float(np.mean(recalls)) if recalls else np.nan,
        "per_class": per_class,
        "confusion_matrix": {
            str(idx): {str(col): int(value) for col, value in row.items()}
            for idx, row in matrix.to_dict(orient="index").items()
        },
    }


def _evaluate_predictions(scored: pd.DataFrame, classes: Sequence[str], prefix: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for split in ("train", "validation", "test"):
        part = scored[scored["training_split"].eq(split)]
        out[split] = _confusion_metrics(part["main_teacher_target"], part[f"{prefix}_label"], classes)
    for source in ("real_candidate", "injection_recovery"):
        part = scored[scored["source_kind"].eq(source)]
        out[f"source_{source}"] = _confusion_metrics(part["main_teacher_target"], part[f"{prefix}_label"], classes)
    return out


def train_teacher_student_smoke(
    *,
    feature_table: pd.DataFrame,
    out_dir: Path,
    min_class_count: int = 40,
    pseudo_min_confidence: float = 0.98,
    pseudo_min_margin: float = 0.50,
    pseudo_weight: float = 0.25,
    random_state: int = 56,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = feature_table.copy()
    if "training_split" not in df:
        df = add_deterministic_splits(df, random_state=random_state)
    if "main_teacher_include" in df:
        df["main_teacher_include"] = _as_bool(df["main_teacher_include"], df.index)

    classes = select_teacher_classes(df, min_class_count)
    if len(classes) < 2:
        raise ValueError(f"need at least two teacher classes, got {classes}")

    feature_sets = {
        "shape_only": shape_feature_columns(df),
        "metadata_only": metadata_feature_columns(df),
    }
    feature_sets["combined"] = sorted(set(feature_sets["shape_only"]) | set(feature_sets["metadata_only"]))

    profiles: dict[str, Any] = {}
    best_profile = ""
    best_score = -np.inf
    for name, columns in feature_sets.items():
        profile_dir = out_dir / name
        profile_dir.mkdir(parents=True, exist_ok=True)
        if not columns:
            profiles[name] = {"status": "skipped_no_features"}
            continue
        leaks = leakage_columns(columns)
        if leaks:
            raise ValueError(f"{name} contains leakage columns: {leaks}")

        work = _split_eval_rows(df, classes)
        train = work[work["training_split"].eq("train")].copy()
        if train["main_teacher_target"].nunique() < 2:
            profiles[name] = {"status": "skipped_insufficient_train_classes"}
            continue
        feature_cfg = FeatureConfig(
            id_columns=("row_id",),
            feature_columns=tuple(columns),
            categorical_columns=(),
            exclude_columns=(),
            clip_sigma=8.0,
        )
        spec = fit_feature_spec(train, feature_cfg)
        teacher = train_softmax_model(
            train,
            train["main_teacher_target"].astype(str),
            spec,
            LogisticConfig(max_iter=500, l2=1.0e-3, tol=1.0e-7, class_weight_balanced=True),
            metadata={"stage": "recovery50_teacher", "feature_profile": name, "classes": list(classes)},
        )
        scored = append_predictions(work, teacher, "teacher")
        predictions_path = write_table(scored, profile_dir / "teacher_predictions.parquet")
        save_model(teacher, profile_dir / "teacher_model.npz")
        metrics = _evaluate_predictions(scored, teacher.classes, "teacher")
        val_score = metrics.get("validation", {}).get("balanced_accuracy", np.nan)
        if np.isfinite(val_score) and val_score > best_score:
            best_score = float(val_score)
            best_profile = name

        unlabeled = df[~_as_bool(df["main_teacher_include"], df.index)].copy()
        if "human_label" in unlabeled:
            unlabeled = unlabeled[~unlabeled["human_label"].fillna("").astype(str).eq("skip")].copy()
        unlabeled_scored = append_predictions(unlabeled, teacher, "teacher") if len(unlabeled) else unlabeled
        pseudo = select_pseudo_labels(
            unlabeled_scored,
            type(
                "Cfg",
                (),
                {
                    "pseudo_min_confidence": pseudo_min_confidence,
                    "pseudo_min_margin": pseudo_min_margin,
                    "max_pseudo_per_class": 5000,
                },
            )(),
            "teacher",
        ) if len(unlabeled_scored) else pd.DataFrame()
        pseudo_path = write_table(pseudo, profile_dir / "pseudo_labels.parquet")

        student_summary: dict[str, Any] = {"status": "no_pseudo_labels", "n_pseudo": int(len(pseudo))}
        if len(pseudo):
            pseudo_train = pseudo.copy()
            pseudo_train["main_teacher_target"] = pseudo_train["pseudo_label"].astype(str)
            student_train = pd.concat([train, pseudo_train], ignore_index=True)
            student_labels = student_train["main_teacher_target"].astype(str)
            weights = np.concatenate(
                [
                    np.ones(len(train), dtype=np.float64),
                    pseudo_weight * pseudo_train["pseudo_confidence"].to_numpy(dtype=np.float64),
                ]
            )
            student = train_softmax_model(
                student_train,
                student_labels,
                spec,
                LogisticConfig(max_iter=500, l2=1.0e-3, tol=1.0e-7, class_weight_balanced=True),
                sample_weight=weights,
                metadata={
                    "stage": "recovery50_student",
                    "feature_profile": name,
                    "n_human": int(len(train)),
                    "n_pseudo": int(len(pseudo)),
                    "pseudo_weight": float(pseudo_weight),
                },
            )
            student_scored = append_predictions(work, student, "student")
            save_model(student, profile_dir / "student_model.npz")
            student_predictions_path = write_table(student_scored, profile_dir / "student_predictions.parquet")
            student_summary = {
                "status": "trained",
                "n_pseudo": int(len(pseudo)),
                "predictions": str(student_predictions_path),
                "metrics": _evaluate_predictions(student_scored, student.classes, "student"),
            }

        profile_summary = {
            "status": "ok",
            "feature_profile": name,
            "n_features": int(len(columns)),
            "feature_columns": columns,
            "classes": list(teacher.classes),
            "teacher_metrics": metrics,
            "teacher_predictions": str(predictions_path),
            "pseudo_labels": str(pseudo_path),
            "student": student_summary,
        }
        (profile_dir / "summary.json").write_text(json.dumps(profile_summary, indent=2, sort_keys=True, default=json_default) + "\n")
        profiles[name] = profile_summary

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "out_dir": str(out_dir),
        "classes": classes,
        "min_class_count": int(min_class_count),
        "best_profile_by_validation_balanced_accuracy": best_profile,
        "best_validation_balanced_accuracy": best_score if np.isfinite(best_score) else None,
        "profiles": profiles,
        "label_policy": "main human-visible teacher target; injection truth auxiliary/evaluation only",
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    (out_dir / "summary.md").write_text(render_training_markdown(summary))
    return summary


def render_training_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# S56 Recovery50 Teacher/Student Smoke",
        "",
        f"- classes: `{', '.join(summary['classes'])}`",
        f"- best profile: `{summary['best_profile_by_validation_balanced_accuracy']}`",
        f"- best validation balanced accuracy: `{summary['best_validation_balanced_accuracy']}`",
        "",
        "| Profile | Status | Features | Validation Balanced Accuracy | Pseudo Labels |",
        "|---|---|---:|---:|---:|",
    ]
    for name, payload in summary["profiles"].items():
        if payload.get("status") != "ok":
            lines.append(f"| `{name}` | `{payload.get('status')}` | 0 |  |  |")
            continue
        val = payload["teacher_metrics"]["validation"]["balanced_accuracy"]
        pseudo_n = payload["student"].get("n_pseudo", 0)
        lines.append(f"| `{name}` | `ok` | {payload['n_features']} | {val:.3f} | {pseudo_n} |")
    lines.append("")
    return "\n".join(lines)
