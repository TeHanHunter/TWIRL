"""Figures for the development-only Teacher v4-SSL full-pool evaluation."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from twirl.plotting.style import apply_twirl_style, get_ordered_palette
from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES


FIGURE_SCHEMA = "twirl_teacher_ssl_fullpool_evaluation_figures_v1"
MODEL_ORDER: tuple[str, ...] = (
    "teacher_v3",
    "fullpool_ssl_frozen_linear_probe",
    "fullpool_ssl_fine_tuned",
)
MODEL_LABELS: Mapping[str, str] = {
    "teacher_v3": "Teacher v3",
    "fullpool_ssl_frozen_linear_probe": "SSL frozen probe",
    "fullpool_ssl_fine_tuned": "SSL fine-tuned",
}
POLICY_LABELS: Mapping[str, str] = {
    "uncertain_as_other": "Uncertain as Other",
    "uncertain_masked": "Uncertain masked",
}
CLASS_LABELS: Mapping[str, str] = {
    "planet_like": "Planet",
    "eclipse_contact": "EB",
    "smooth_variable": "Variable",
    "other": "Other",
}
HUMAN_LABELS: tuple[tuple[str, str], ...] = (
    ("uncertain", "Uncertain"),
    ("instrumental_or_systematic", "Instrument/systematic"),
    ("stellar_variability", "Stellar variability"),
    ("wide_transit_like", "Wide transit-like"),
    ("eclipsing_binary_or_pceb", "EB/PCEB"),
    ("planet_like", "Planet-like"),
)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: Mapping[str, Any]) -> str:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return _file_sha256(path)


def _save_figure(figure: Any, *, stem: Path) -> dict[str, Any]:
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    figure.savefig(png, dpi=300, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    return {
        "png": str(png),
        "png_sha256": _file_sha256(png),
        "pdf": str(pdf),
        "pdf_sha256": _file_sha256(pdf),
    }


def compute_reference_umap(
    *,
    embedding_path: Path,
    output_dir: Path,
    random_state: int = 560064,
    n_neighbors: int = 30,
    min_dist: float = 0.15,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Fit one deterministic UMAP in the coherent fold-0 reference space."""

    try:
        import umap
        from sklearn.decomposition import PCA
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "UMAP plotting requires the optional 'representation' dependencies"
        ) from exc

    embedding_path = Path(embedding_path).expanduser().resolve()
    with np.load(embedding_path) as payload:
        required = {
            "review_id",
            "ssl_observation_id",
            "tic",
            "sector",
            "cv_fold",
            "human_label",
            "morphology_target",
            "period_d",
            "duration_min",
            "embedding",
        }
        missing = sorted(required - set(payload.files))
        if missing:
            raise KeyError(f"reference embedding export lacks {missing}")
        embedding = payload["embedding"].astype(np.float32)
        metadata = {
            name: payload[name].copy()
            for name in required
            if name != "embedding"
        }
    if embedding.ndim != 2 or embedding.shape[1] != 128:
        raise ValueError(f"expected an (n,128) embedding, got {embedding.shape}")
    if not np.isfinite(embedding).all():
        raise ValueError("reference embeddings contain non-finite values")
    norm = np.linalg.norm(embedding, axis=1, keepdims=True)
    normalized = embedding / np.maximum(norm, 1.0e-8)
    pca = PCA(n_components=50, svd_solver="randomized", random_state=int(random_state))
    reduced = pca.fit_transform(normalized)
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=int(n_neighbors),
        min_dist=float(min_dist),
        metric="euclidean",
        random_state=int(random_state),
        transform_seed=int(random_state),
        n_jobs=1,
        low_memory=True,
    )
    coordinates = reducer.fit_transform(reduced)
    if coordinates.shape != (len(embedding), 2) or not np.isfinite(
        coordinates
    ).all():
        raise RuntimeError("UMAP produced an invalid coordinate matrix")
    frame = pd.DataFrame(
        {
            "review_id": metadata["review_id"].astype(str),
            "ssl_observation_id": metadata["ssl_observation_id"].astype(str),
            "tic": metadata["tic"].astype(np.int64),
            "sector": metadata["sector"].astype(np.int16),
            "cv_fold": metadata["cv_fold"].astype(np.int8),
            "human_label": metadata["human_label"].astype(str),
            "morphology_target": metadata["morphology_target"].astype(np.int8),
            "period_d": metadata["period_d"].astype(np.float64),
            "duration_min": metadata["duration_min"].astype(np.float64),
            "umap_1": coordinates[:, 0],
            "umap_2": coordinates[:, 1],
        }
    )
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    coordinates_path = output_dir / "reference_fold_0_umap_coordinates.csv"
    frame.to_csv(
        coordinates_path,
        index=False,
        lineterminator="\n",
        float_format="%.9g",
    )
    summary = {
        "schema_version": FIGURE_SCHEMA,
        "scope": (
            "exploratory_fold_0_reference_encoder; labels_overlay_only; "
            "not_a_performance_estimate"
        ),
        "embedding_path": str(embedding_path),
        "embedding_sha256": _file_sha256(embedding_path),
        "n_rows": int(len(frame)),
        "input_embedding_dim": int(embedding.shape[1]),
        "pca_components": 50,
        "pca_explained_variance_fraction": float(
            pca.explained_variance_ratio_.sum()
        ),
        "umap": {
            "n_neighbors": int(n_neighbors),
            "min_dist": float(min_dist),
            "metric": "euclidean_after_unit_norm_and_pca50",
            "random_state": int(random_state),
        },
        "versions": {
            "umap-learn": importlib.metadata.version("umap-learn"),
            "scikit-learn": importlib.metadata.version("scikit-learn"),
        },
        "coordinates": str(coordinates_path),
        "coordinates_sha256": _file_sha256(coordinates_path),
    }
    _write_json(output_dir / "reference_fold_0_umap.summary.json", summary)
    return frame, summary


def plot_reference_umap(
    *,
    coordinates: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    """Plot label and BLS-period views of one deterministic UMAP."""

    import matplotlib.pyplot as plt

    required = {"human_label", "period_d", "umap_1", "umap_2"}
    missing = sorted(required - set(coordinates.columns))
    if missing:
        raise KeyError(f"UMAP coordinates lack {missing}")
    template = apply_twirl_style("full_page")
    figure, axes = plt.subplots(
        1,
        2,
        figsize=template["figsize"],
        constrained_layout=True,
    )
    left, right = axes
    left.grid(False)
    right.grid(False)
    left.scatter(
        coordinates["umap_1"],
        coordinates["umap_2"],
        s=2.0,
        c="0.82",
        alpha=0.30,
        linewidths=0,
        rasterized=True,
    )
    colors = get_ordered_palette(len(HUMAN_LABELS))
    for (label, display), color in zip(HUMAN_LABELS, colors):
        mask = coordinates["human_label"].astype(str).eq(label)
        if not mask.any():
            continue
        rare = label in {
            "planet_like",
            "eclipsing_binary_or_pceb",
            "wide_transit_like",
        }
        left.scatter(
            coordinates.loc[mask, "umap_1"],
            coordinates.loc[mask, "umap_2"],
            s=12.0 if rare else 4.0,
            c=[color],
            alpha=0.85 if rare else 0.45,
            linewidths=0.25 if rare else 0,
            edgecolors="white" if rare else "none",
            label=f"{display} ({int(mask.sum()):,})",
            rasterized=True,
            zorder=3 if rare else 2,
        )
    left.set_title("Frozen SSL representation by human label")
    left.set_xlabel("UMAP 1")
    left.set_ylabel("UMAP 2")
    left.set_xticks([])
    left.set_yticks([])
    left.legend(
        loc="upper left",
        bbox_to_anchor=(0.01, 0.99),
        borderaxespad=0,
        framealpha=0.92,
        markerscale=1.2,
    )

    period = pd.to_numeric(coordinates["period_d"], errors="coerce").to_numpy(
        dtype=float
    )
    if not np.isfinite(period).all() or (period <= 0).any():
        raise ValueError("UMAP period coloring requires finite positive periods")
    scatter = right.scatter(
        coordinates["umap_1"],
        coordinates["umap_2"],
        s=4.0,
        c=np.log10(period),
        cmap="viridis",
        alpha=0.65,
        linewidths=0,
        rasterized=True,
    )
    right.set_title("Same representation by BLS period")
    right.set_xlabel("UMAP 1")
    right.set_ylabel("UMAP 2")
    right.set_xticks([])
    right.set_yticks([])
    colorbar = figure.colorbar(scatter, ax=right, pad=0.02, fraction=0.05)
    colorbar.set_label(r"$\log_{10}(P / \mathrm{day})$")
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    artifacts = _save_figure(
        figure,
        stem=output_dir / "teacher_ssl_fullpool_reference_umap",
    )
    plt.close(figure)
    return artifacts


def _validate_performance_tables(
    performance: pd.DataFrame,
    per_class: pd.DataFrame,
) -> None:
    expected_pairs = {
        (model, policy)
        for model in MODEL_ORDER
        for policy in POLICY_LABELS
    }
    observed_pairs = set(
        zip(performance["model"].astype(str), performance["label_policy"].astype(str))
    )
    if observed_pairs != expected_pairs:
        raise ValueError(
            "performance table lacks the exact model/policy comparison grid"
        )
    observed_class = set(
        zip(
            per_class["model"].astype(str),
            per_class["label_policy"].astype(str),
            per_class["class"].astype(str),
        )
    )
    expected_class = {
        (model, policy, label)
        for model in MODEL_ORDER
        for policy in POLICY_LABELS
        for label in MORPHOLOGY_CLASSES
    }
    if observed_class != expected_class:
        raise ValueError("per-class table lacks the exact comparison grid")


def plot_development_performance(
    *,
    performance_path: Path,
    per_class_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Plot matched OOF aggregate scores and per-class average precision."""

    import matplotlib.pyplot as plt

    performance_path = Path(performance_path).expanduser().resolve()
    per_class_path = Path(per_class_path).expanduser().resolve()
    performance = pd.read_csv(performance_path)
    per_class = pd.read_csv(per_class_path)
    _validate_performance_tables(performance, per_class)
    template = apply_twirl_style("full_page")
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(template["figsize"][0], 6.0),
        constrained_layout=True,
    )
    colors = get_ordered_palette(len(MODEL_ORDER))
    aggregate_metrics = (
        ("macro_average_precision", "Macro AP"),
        ("macro_f1", "Macro F1"),
        ("balanced_accuracy", "Balanced accuracy"),
    )
    x = np.arange(len(aggregate_metrics), dtype=float)
    width = 0.24
    for column, policy in enumerate(POLICY_LABELS):
        axis = axes[0, column]
        for model_index, (model, color) in enumerate(zip(MODEL_ORDER, colors)):
            row = performance.loc[
                performance["model"].eq(model)
                & performance["label_policy"].eq(policy)
            ].iloc[0]
            values = np.asarray(
                [float(row[name]) for name, _ in aggregate_metrics]
            )
            low = np.asarray(
                [float(row[f"{name}_low"]) for name, _ in aggregate_metrics]
            )
            high = np.asarray(
                [float(row[f"{name}_high"]) for name, _ in aggregate_metrics]
            )
            position = x + (model_index - 1) * width
            axis.bar(
                position,
                values,
                width=width,
                color=color,
                edgecolor="white",
                linewidth=0.5,
                label=MODEL_LABELS[model],
                zorder=2,
            )
            axis.errorbar(
                position,
                values,
                yerr=np.vstack([values - low, high - values]),
                fmt="none",
                ecolor="0.18",
                elinewidth=0.7,
                capsize=1.8,
                capthick=0.7,
                zorder=3,
            )
        axis.set_title(POLICY_LABELS[policy])
        axis.set_xticks(x, [display for _, display in aggregate_metrics])
        axis.tick_params(axis="x", rotation=18)
        axis.set_ylim(0.0, 1.02)
        axis.set_ylabel("Development OOF score" if column == 0 else "")
        if column == 1:
            axis.legend(loc="lower right")

    class_x = np.arange(len(MORPHOLOGY_CLASSES), dtype=float)
    for column, policy in enumerate(POLICY_LABELS):
        axis = axes[1, column]
        for model_index, (model, color) in enumerate(zip(MODEL_ORDER, colors)):
            subset = per_class.loc[
                per_class["model"].eq(model)
                & per_class["label_policy"].eq(policy)
            ].set_index("class")
            values = np.asarray(
                [float(subset.loc[label, "average_precision"]) for label in MORPHOLOGY_CLASSES]
            )
            position = class_x + (model_index - 1) * width
            axis.bar(
                position,
                values,
                width=width,
                color=color,
                edgecolor="white",
                linewidth=0.5,
                zorder=2,
            )
        axis.set_xticks(
            class_x,
            [CLASS_LABELS[label] for label in MORPHOLOGY_CLASSES],
        )
        axis.tick_params(axis="x", rotation=18)
        axis.set_ylim(0.0, 1.02)
        axis.set_ylabel("One-vs-rest average precision" if column == 0 else "")
        axis.set_title(f"Per-class AP — {POLICY_LABELS[policy].lower()}")
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    artifacts = _save_figure(
        figure,
        stem=output_dir / "teacher_ssl_fullpool_development_performance",
    )
    plt.close(figure)
    return {
        **artifacts,
        "performance_input": str(performance_path),
        "performance_input_sha256": _file_sha256(performance_path),
        "per_class_input": str(per_class_path),
        "per_class_input_sha256": _file_sha256(per_class_path),
    }


def make_fullpool_evaluation_figures(
    *,
    embedding_path: Path,
    performance_path: Path,
    per_class_path: Path,
    output_dir: Path,
    random_state: int = 560064,
) -> dict[str, Any]:
    """Build deterministic UMAP and matched-development performance figures."""

    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    coordinates, umap_summary = compute_reference_umap(
        embedding_path=embedding_path,
        output_dir=output_dir,
        random_state=int(random_state),
    )
    umap_artifacts = plot_reference_umap(
        coordinates=coordinates,
        output_dir=output_dir,
    )
    performance_artifacts = plot_development_performance(
        performance_path=performance_path,
        per_class_path=per_class_path,
        output_dir=output_dir,
    )
    summary = {
        "schema_version": FIGURE_SCHEMA,
        "umap": {**umap_summary, "figure": umap_artifacts},
        "development_performance": performance_artifacts,
    }
    _write_json(output_dir / "teacher_ssl_fullpool_figures.provenance.json", summary)
    return summary


__all__ = [
    "compute_reference_umap",
    "make_fullpool_evaluation_figures",
    "plot_development_performance",
    "plot_reference_umap",
]
