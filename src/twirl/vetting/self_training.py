"""Semi-supervised teacher/student triage for TWIRL candidate tables.

This module deliberately starts with candidate-level engineered features, not
raw light-curve tensors. The transparent BLS/dip-search/vetter stack remains
the detection layer; this code learns a triage ranking from human labels,
injections, and high-confidence pseudo-labels with explicit provenance.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp


DEFAULT_ID_COLUMNS: tuple[str, ...] = ("tic",)
DEFAULT_LABEL_COLUMN = "label"
DEFAULT_SOURCE_COLUMN = "label_source"
DEFAULT_CATEGORICAL_COLUMNS: tuple[str, ...] = (
    "vet_class",
    "rep_aperture",
)
DEFAULT_EXCLUDE_COLUMNS: frozenset[str] = frozenset(
    {
        "tic",
        "sector",
        "cam",
        "ccd",
        "label",
        "label_source",
        "labeler",
        "notes",
        "hlsp_path",
        "apertures_agree",
        "bls_run_id",
        "human_label",
        "pseudo_label",
        "teacher_label",
        "student_label",
    }
)


@dataclass
class FeatureConfig:
    """Feature selection and preprocessing controls."""

    id_columns: tuple[str, ...] = DEFAULT_ID_COLUMNS
    feature_columns: tuple[str, ...] | None = None
    categorical_columns: tuple[str, ...] = DEFAULT_CATEGORICAL_COLUMNS
    exclude_columns: tuple[str, ...] = tuple(sorted(DEFAULT_EXCLUDE_COLUMNS))
    max_categorical_levels: int = 16
    clip_sigma: float = 8.0


@dataclass
class LogisticConfig:
    """Small multinomial logistic-regression trainer configuration."""

    max_iter: int = 500
    l2: float = 1.0e-3
    tol: float = 1.0e-7
    class_weight_balanced: bool = True


@dataclass
class SelfTrainingConfig:
    """Teacher/student pseudo-labeling controls."""

    feature: FeatureConfig = field(default_factory=FeatureConfig)
    trainer: LogisticConfig = field(default_factory=LogisticConfig)
    label_column: str = DEFAULT_LABEL_COLUMN
    source_column: str = DEFAULT_SOURCE_COLUMN
    unlabeled_values: tuple[str, ...] = ("", "unknown", "unlabeled", "needs_review")
    pseudo_min_confidence: float = 0.98
    pseudo_min_margin: float = 0.50
    pseudo_weight: float = 0.25
    max_pseudo_per_class: int | None = 5000
    review_queue_size: int = 250


@dataclass
class FeatureSpec:
    """Fitted feature preprocessing state."""

    numeric_columns: tuple[str, ...]
    categorical_levels: dict[str, tuple[str, ...]]
    feature_names: tuple[str, ...]
    medians: np.ndarray
    scales: np.ndarray
    clip_sigma: float = 8.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "numeric_columns": list(self.numeric_columns),
            "categorical_levels": {
                key: list(vals) for key, vals in self.categorical_levels.items()
            },
            "feature_names": list(self.feature_names),
            "medians": self.medians.tolist(),
            "scales": self.scales.tolist(),
            "clip_sigma": float(self.clip_sigma),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "FeatureSpec":
        return cls(
            numeric_columns=tuple(data["numeric_columns"]),
            categorical_levels={
                str(key): tuple(str(v) for v in vals)
                for key, vals in data["categorical_levels"].items()
            },
            feature_names=tuple(data["feature_names"]),
            medians=np.asarray(data["medians"], dtype=np.float64),
            scales=np.asarray(data["scales"], dtype=np.float64),
            clip_sigma=float(data.get("clip_sigma", 8.0)),
        )


@dataclass
class SoftmaxModel:
    """Multiclass linear model with fitted feature state."""

    classes: tuple[str, ...]
    feature_spec: FeatureSpec
    weights: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)

    def predict_proba(self, df: pd.DataFrame) -> np.ndarray:
        x = transform_features(df, self.feature_spec)
        return _softmax(_with_intercept(x) @ self.weights)

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        p = self.predict_proba(df)
        return np.asarray(self.classes, dtype=object)[np.argmax(p, axis=1)]


@dataclass
class SelfTrainingResult:
    teacher: SoftmaxModel
    student: SoftmaxModel
    scored: pd.DataFrame
    pseudo_labels: pd.DataFrame
    review_queue: pd.DataFrame
    summary: dict[str, Any]


def _as_tuple(value: Iterable[str] | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    return tuple(str(v) for v in value)


def config_from_mapping(data: dict[str, Any] | None) -> SelfTrainingConfig:
    """Build a config from a YAML/JSON-compatible mapping."""

    if not data:
        return SelfTrainingConfig()
    feature_data = dict(data.get("feature", {}))
    trainer_data = dict(data.get("trainer", {}))
    feature = FeatureConfig(
        id_columns=_as_tuple(feature_data.get("id_columns")) or DEFAULT_ID_COLUMNS,
        feature_columns=_as_tuple(feature_data.get("feature_columns")),
        categorical_columns=(
            _as_tuple(feature_data.get("categorical_columns"))
            or DEFAULT_CATEGORICAL_COLUMNS
        ),
        exclude_columns=(
            _as_tuple(feature_data.get("exclude_columns"))
            or tuple(sorted(DEFAULT_EXCLUDE_COLUMNS))
        ),
        max_categorical_levels=int(
            feature_data.get("max_categorical_levels", 16)
        ),
        clip_sigma=float(feature_data.get("clip_sigma", 8.0)),
    )
    trainer = LogisticConfig(
        max_iter=int(trainer_data.get("max_iter", 500)),
        l2=float(trainer_data.get("l2", 1.0e-3)),
        tol=float(trainer_data.get("tol", 1.0e-7)),
        class_weight_balanced=bool(
            trainer_data.get("class_weight_balanced", True)
        ),
    )
    return SelfTrainingConfig(
        feature=feature,
        trainer=trainer,
        label_column=str(data.get("label_column", DEFAULT_LABEL_COLUMN)),
        source_column=str(data.get("source_column", DEFAULT_SOURCE_COLUMN)),
        unlabeled_values=tuple(
            str(v).lower() for v in data.get(
                "unlabeled_values",
                ("", "unknown", "unlabeled", "needs_review"),
            )
        ),
        pseudo_min_confidence=float(data.get("pseudo_min_confidence", 0.98)),
        pseudo_min_margin=float(data.get("pseudo_min_margin", 0.50)),
        pseudo_weight=float(data.get("pseudo_weight", 0.25)),
        max_pseudo_per_class=(
            None if data.get("max_pseudo_per_class") is None
            else int(data.get("max_pseudo_per_class", 5000))
        ),
        review_queue_size=int(data.get("review_queue_size", 250)),
    )


def fit_feature_spec(df: pd.DataFrame, cfg: FeatureConfig | None = None) -> FeatureSpec:
    """Choose and fit robust preprocessing for candidate-table features."""

    cfg = cfg or FeatureConfig()
    exclude = set(cfg.exclude_columns)
    numeric_columns: list[str] = []
    categorical_levels: dict[str, tuple[str, ...]] = {}

    if cfg.feature_columns is not None:
        requested = [c for c in cfg.feature_columns if c in df.columns]
        numeric_columns = [
            c for c in requested if pd.api.types.is_numeric_dtype(df[c])
        ]
    else:
        for col in df.columns:
            if col in exclude:
                continue
            if pd.api.types.is_bool_dtype(df[col]) or pd.api.types.is_numeric_dtype(df[col]):
                numeric_columns.append(col)

    for col in cfg.categorical_columns:
        if col not in df.columns:
            continue
        vals = (
            df[col]
            .dropna()
            .astype(str)
            .value_counts()
            .head(cfg.max_categorical_levels)
            .index
            .tolist()
        )
        if vals:
            categorical_levels[col] = tuple(vals)

    raw = _raw_feature_matrix(df, tuple(numeric_columns), categorical_levels)
    med = np.nanmedian(raw, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    centered = raw - med
    mad = np.nanmedian(np.abs(centered), axis=0)
    scale = 1.4826 * mad
    std = np.nanstd(raw, axis=0)
    scale = np.where((scale > 0) & np.isfinite(scale), scale, std)
    scale = np.where((scale > 0) & np.isfinite(scale), scale, 1.0)

    feature_names = list(numeric_columns)
    for col, levels in categorical_levels.items():
        feature_names.extend([f"{col}={level}" for level in levels])

    return FeatureSpec(
        numeric_columns=tuple(numeric_columns),
        categorical_levels=categorical_levels,
        feature_names=tuple(feature_names),
        medians=med.astype(np.float64),
        scales=scale.astype(np.float64),
        clip_sigma=float(cfg.clip_sigma),
    )


def _raw_feature_matrix(
    df: pd.DataFrame,
    numeric_columns: tuple[str, ...],
    categorical_levels: dict[str, tuple[str, ...]],
) -> np.ndarray:
    chunks: list[np.ndarray] = []
    for col in numeric_columns:
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=np.float64)
        chunks.append(vals[:, None])
    for col, levels in categorical_levels.items():
        series = df[col].astype(str).fillna("")
        for level in levels:
            chunks.append((series == level).to_numpy(dtype=np.float64)[:, None])
    if not chunks:
        raise ValueError("no usable numeric or categorical features found")
    return np.hstack(chunks)


def transform_features(df: pd.DataFrame, spec: FeatureSpec) -> np.ndarray:
    raw = _raw_feature_matrix(df, spec.numeric_columns, spec.categorical_levels)
    x = (raw - spec.medians) / spec.scales
    x = np.where(np.isfinite(x), x, 0.0)
    if spec.clip_sigma and spec.clip_sigma > 0:
        x = np.clip(x, -spec.clip_sigma, spec.clip_sigma)
    return x.astype(np.float64)


def _with_intercept(x: np.ndarray) -> np.ndarray:
    return np.hstack([np.ones((x.shape[0], 1), dtype=x.dtype), x])


def _softmax(logits: np.ndarray) -> np.ndarray:
    z = logits - logsumexp(logits, axis=1)[:, None]
    return np.exp(z)


def _encode_labels(labels: pd.Series) -> tuple[np.ndarray, tuple[str, ...]]:
    classes = tuple(sorted(str(v) for v in labels.dropna().unique()))
    if len(classes) < 2:
        raise ValueError("at least two classes are required to train the model")
    index = {cls: i for i, cls in enumerate(classes)}
    y = np.asarray([index[str(v)] for v in labels], dtype=np.int64)
    return y, classes


def train_softmax_model(
    df: pd.DataFrame,
    labels: pd.Series,
    feature_spec: FeatureSpec,
    cfg: LogisticConfig | None = None,
    sample_weight: np.ndarray | None = None,
    metadata: dict[str, Any] | None = None,
) -> SoftmaxModel:
    """Fit a small multiclass logistic-regression model with scipy L-BFGS."""

    cfg = cfg or LogisticConfig()
    y, classes = _encode_labels(labels)
    x = _with_intercept(transform_features(df, feature_spec))
    n, p = x.shape
    k = len(classes)

    sw = (
        np.ones(n, dtype=np.float64)
        if sample_weight is None
        else np.asarray(sample_weight, dtype=np.float64).copy()
    )
    sw = np.where(np.isfinite(sw) & (sw > 0), sw, 0.0)
    if cfg.class_weight_balanced:
        counts = np.bincount(y, minlength=k).astype(np.float64)
        class_weight = n / np.maximum(k * counts, 1.0)
        sw *= class_weight[y]
    sw_sum = float(sw.sum())
    if sw_sum <= 0:
        raise ValueError("sample weights sum to zero")
    sw = sw * (n / sw_sum)

    y_onehot = np.zeros((n, k), dtype=np.float64)
    y_onehot[np.arange(n), y] = 1.0

    l2_mask = np.ones((p, k), dtype=np.float64)
    l2_mask[0, :] = 0.0

    def objective(flat: np.ndarray) -> tuple[float, np.ndarray]:
        w = flat.reshape(p, k)
        logits = x @ w
        log_norm = logsumexp(logits, axis=1)
        loss_vec = log_norm - logits[np.arange(n), y]
        loss = float(np.sum(sw * loss_vec) / n)
        if cfg.l2 > 0:
            loss += 0.5 * cfg.l2 * float(np.sum((w * l2_mask) ** 2))
        prob = np.exp(logits - log_norm[:, None])
        grad = x.T @ ((prob - y_onehot) * sw[:, None]) / n
        if cfg.l2 > 0:
            grad += cfg.l2 * w * l2_mask
        return loss, grad.ravel()

    result = minimize(
        objective,
        np.zeros(p * k, dtype=np.float64),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": int(cfg.max_iter), "ftol": float(cfg.tol)},
    )
    meta = dict(metadata or {})
    meta.update(
        {
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_nit": int(result.nit),
            "optimizer_fun": float(result.fun),
            "trainer": asdict(cfg),
        }
    )
    return SoftmaxModel(
        classes=classes,
        feature_spec=feature_spec,
        weights=result.x.reshape(p, k),
        metadata=meta,
    )


def _merge_labels(
    candidates: pd.DataFrame,
    labels: pd.DataFrame,
    cfg: SelfTrainingConfig,
) -> pd.DataFrame:
    id_cols = list(cfg.feature.id_columns)
    missing = [c for c in id_cols if c not in candidates.columns or c not in labels.columns]
    if missing:
        raise KeyError(f"missing label join columns: {missing}")
    if cfg.label_column not in labels.columns:
        raise KeyError(f"labels table missing {cfg.label_column!r}")
    label_cols = id_cols + [
        c for c in labels.columns
        if c not in id_cols and c in {cfg.label_column, cfg.source_column, "labeler", "notes"}
    ]
    lab = labels[label_cols].copy()
    lab = lab.drop_duplicates(subset=id_cols, keep="last")
    return candidates.merge(lab, on=id_cols, how="left", suffixes=("", "_label"))


def labeled_mask(frame: pd.DataFrame, cfg: SelfTrainingConfig) -> np.ndarray:
    values = frame[cfg.label_column].fillna("").astype(str).str.lower()
    return ~values.isin(set(cfg.unlabeled_values))


def _probability_columns(prefix: str, classes: tuple[str, ...]) -> list[str]:
    return [f"{prefix}_p_{cls}" for cls in classes]


def append_predictions(
    df: pd.DataFrame,
    model: SoftmaxModel,
    prefix: str,
) -> pd.DataFrame:
    out = df.copy()
    prob = model.predict_proba(df)
    cols = _probability_columns(prefix, model.classes)
    for i, col in enumerate(cols):
        out[col] = prob[:, i]
    idx = np.argmax(prob, axis=1)
    part = np.partition(prob, -2, axis=1)
    out[f"{prefix}_label"] = np.asarray(model.classes, dtype=object)[idx]
    out[f"{prefix}_confidence"] = prob[np.arange(prob.shape[0]), idx]
    out[f"{prefix}_margin"] = part[:, -1] - part[:, -2]
    entropy = -np.sum(np.where(prob > 0, prob * np.log(prob), 0.0), axis=1)
    out[f"{prefix}_entropy"] = entropy / np.log(max(prob.shape[1], 2))
    return out


def select_pseudo_labels(
    scored: pd.DataFrame,
    cfg: SelfTrainingConfig,
    prefix: str = "teacher",
) -> pd.DataFrame:
    mask = (
        (scored[f"{prefix}_confidence"] >= cfg.pseudo_min_confidence)
        & (scored[f"{prefix}_margin"] >= cfg.pseudo_min_margin)
    )
    pseudo = scored.loc[mask].copy()
    if cfg.max_pseudo_per_class is not None and not pseudo.empty:
        pseudo = (
            pseudo.sort_values(f"{prefix}_confidence", ascending=False)
            .groupby(f"{prefix}_label", group_keys=False)
            .head(int(cfg.max_pseudo_per_class))
        )
    pseudo["pseudo_label"] = pseudo[f"{prefix}_label"]
    pseudo["pseudo_confidence"] = pseudo[f"{prefix}_confidence"]
    pseudo["pseudo_margin"] = pseudo[f"{prefix}_margin"]
    pseudo["pseudo_source"] = f"{prefix}_self_training"
    return pseudo


def build_review_queue(
    scored: pd.DataFrame,
    is_human_labeled: np.ndarray,
    cfg: SelfTrainingConfig,
    pseudo_index: pd.Index,
) -> pd.DataFrame:
    available = scored.loc[~is_human_labeled].copy()
    if pseudo_index.size:
        available = available.drop(index=pseudo_index, errors="ignore")
    if available.empty:
        return available

    sde = pd.to_numeric(available.get("sde_max", 0.0), errors="coerce").fillna(0.0)
    if float(sde.max()) > float(sde.min()):
        sde_rank = (sde - sde.min()) / (sde.max() - sde.min())
    else:
        sde_rank = pd.Series(0.0, index=available.index)
    entropy = pd.to_numeric(
        available["student_entropy"], errors="coerce"
    ).fillna(0.0)
    confidence = pd.to_numeric(
        available["student_confidence"], errors="coerce"
    ).fillna(0.0)
    # Prioritize uncertain examples, while still surfacing strong candidates.
    available["review_priority"] = 0.65 * entropy + 0.25 * sde_rank + 0.10 * confidence
    keep_cols = [
        c for c in [
            *cfg.feature.id_columns,
            "sector",
            "tmag",
            "vet_class",
            "class_rank",
            "period_d",
            "duration_min",
            "depth",
            "sde_max",
            "n_apertures_agree",
            "student_label",
            "student_confidence",
            "student_margin",
            "student_entropy",
            "review_priority",
        ]
        if c in available.columns
    ]
    return (
        available.sort_values("review_priority", ascending=False)
        .head(cfg.review_queue_size)
        .loc[:, keep_cols]
        .reset_index(drop=True)
    )


def train_teacher_student(
    candidates: pd.DataFrame,
    labels: pd.DataFrame,
    cfg: SelfTrainingConfig | None = None,
) -> SelfTrainingResult:
    """Train a teacher on human labels, pseudo-label, then train a student."""

    cfg = cfg or SelfTrainingConfig()
    merged = _merge_labels(candidates, labels, cfg)
    human = labeled_mask(merged, cfg)
    labeled = merged.loc[human].copy()
    if labeled.empty:
        raise ValueError("no human/injection labels available for teacher training")

    spec = fit_feature_spec(merged, cfg.feature)
    teacher = train_softmax_model(
        labeled,
        labeled[cfg.label_column].astype(str),
        spec,
        cfg.trainer,
        metadata={"stage": "teacher", "n_human": int(len(labeled))},
    )

    teacher_scored = append_predictions(merged, teacher, "teacher")
    unlabeled_scored = teacher_scored.loc[~human].copy()
    pseudo = select_pseudo_labels(unlabeled_scored, cfg, "teacher")

    student_frames = [labeled]
    student_labels = [labeled[cfg.label_column].astype(str)]
    student_weights = [np.ones(len(labeled), dtype=np.float64)]
    if not pseudo.empty:
        pseudo_train = pseudo.copy()
        pseudo_train[cfg.label_column] = pseudo_train["pseudo_label"]
        student_frames.append(pseudo_train)
        student_labels.append(pseudo_train[cfg.label_column].astype(str))
        student_weights.append(
            cfg.pseudo_weight
            * pseudo_train["pseudo_confidence"].to_numpy(dtype=np.float64)
        )

    train_df = pd.concat(student_frames, ignore_index=True)
    train_labels = pd.concat(student_labels, ignore_index=True)
    train_weight = np.concatenate(student_weights)
    student = train_softmax_model(
        train_df,
        train_labels,
        spec,
        cfg.trainer,
        sample_weight=train_weight,
        metadata={
            "stage": "student",
            "n_human": int(len(labeled)),
            "n_pseudo": int(len(pseudo)),
            "pseudo_weight": float(cfg.pseudo_weight),
            "pseudo_min_confidence": float(cfg.pseudo_min_confidence),
            "pseudo_min_margin": float(cfg.pseudo_min_margin),
        },
    )

    scored = append_predictions(teacher_scored, student, "student")
    review = build_review_queue(scored, human, cfg, pseudo.index)
    summary = {
        "n_candidates": int(len(candidates)),
        "n_human_labeled": int(human.sum()),
        "n_unlabeled": int((~human).sum()),
        "n_pseudo": int(len(pseudo)),
        "n_student_train": int(len(train_df)),
        "classes": list(student.classes),
        "feature_names": list(spec.feature_names),
        "pseudo_counts": (
            pseudo["pseudo_label"].value_counts().sort_index().to_dict()
            if not pseudo.empty
            else {}
        ),
    }
    return SelfTrainingResult(
        teacher=teacher,
        student=student,
        scored=scored,
        pseudo_labels=pseudo.reset_index(drop=True),
        review_queue=review,
        summary=summary,
    )


def save_model(model: SoftmaxModel, path: Path) -> None:
    """Save a model as a compressed NPZ with JSON metadata."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "classes": list(model.classes),
        "feature_spec": model.feature_spec.to_dict(),
        "metadata": model.metadata,
    }
    np.savez_compressed(
        path,
        weights=model.weights,
        payload=np.asarray(json.dumps(payload)),
    )


def load_model(path: Path) -> SoftmaxModel:
    with np.load(path, allow_pickle=False) as data:
        weights = np.asarray(data["weights"], dtype=np.float64)
        payload = json.loads(str(data["payload"]))
    return SoftmaxModel(
        classes=tuple(payload["classes"]),
        feature_spec=FeatureSpec.from_dict(payload["feature_spec"]),
        weights=weights,
        metadata=dict(payload.get("metadata", {})),
    )


def write_summary(summary: dict[str, Any], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


__all__ = [
    "FeatureConfig",
    "LogisticConfig",
    "SelfTrainingConfig",
    "SelfTrainingResult",
    "SoftmaxModel",
    "append_predictions",
    "config_from_mapping",
    "fit_feature_spec",
    "load_model",
    "save_model",
    "select_pseudo_labels",
    "train_teacher_student",
    "write_summary",
]
