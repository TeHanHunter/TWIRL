"""Multi-branch CNN for the S56 post-adjudication morphology teacher.

The model consumes only candidate observables.  Human period corrections are
targets for the auxiliary harmonic head; they never select or alter an input
fold.  All seven folds are generated from the original ADP-small BLS
ephemeris by :mod:`twirl.vetting.harmonic_inputs`.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np
import pandas as pd

from twirl.vetting.adjudication_audit import (
    MORPHOLOGY_CLASSES_V1,
    TRANSIT_HARMONIC_TARGETS_V1,
)
from twirl.vetting.harmonic_inputs import (
    CHRONOLOGY_SMALL_CHANNELS,
    CHRONOLOGY_SUPPLEMENTAL_CHANNELS,
    HARMONIC_VIEW_CHANNELS,
    PERIODOGRAM_CHANNELS,
)


MODEL_VERSION = "s56_harmonic_cnn_v1"
MORPHOLOGY_CLASSES: tuple[str, ...] = MORPHOLOGY_CLASSES_V1
PRESERVE_CLASSES: tuple[str, ...] = ("reject", "preserve")
HARMONIC_CLASSES: tuple[str, ...] = TRANSIT_HARMONIC_TARGETS_V1


@dataclass(frozen=True)
class HarmonicModelConfig:
    small_channels: int = len(CHRONOLOGY_SMALL_CHANNELS)
    supplemental_channels: int = len(CHRONOLOGY_SUPPLEMENTAL_CHANNELS)
    fold_channels: int = len(HARMONIC_VIEW_CHANNELS)
    periodogram_channels: int = len(PERIODOGRAM_CHANNELS)
    metadata_dim: int = 0
    embedding_dim: int = 64
    dropout: float = 0.20
    supplemental_dropout: float = 0.10


@dataclass(frozen=True)
class HarmonicTrainConfig:
    epochs: int = 100
    batch_size: int = 32
    learning_rate: float = 3.0e-4
    weight_decay: float = 1.0e-4
    patience: int = 12
    seed: int = 56
    morphology_loss_weight: float = 1.0
    preserve_loss_weight: float = 0.5
    harmonic_loss_weight: float = 0.25
    class_weight_cap: float = 4.0


MODEL_PROFILES: Mapping[str, frozenset[str]] = {
    # The historical model is evaluated from its existing checkpoint and is
    # deliberately not represented by this architecture.
    "single_period_native_fold": frozenset({"single_fold", "local"}),
    "seven_harmonic_shape": frozenset({"harmonic", "local"}),
    "shape_plus_raw_chronology": frozenset({"harmonic", "local", "chronology"}),
    "shape_plus_periodogram_bls": frozenset({"harmonic", "local", "periodogram", "metadata"}),
    "metadata_only": frozenset({"metadata"}),
    "full_combined": frozenset({"harmonic", "local", "chronology", "periodogram", "metadata"}),
}


def profile_branches(profile: str) -> frozenset[str]:
    if profile not in MODEL_PROFILES:
        raise ValueError(f"unknown harmonic CNN profile {profile!r}; expected {sorted(MODEL_PROFILES)}")
    return MODEL_PROFILES[profile]


def _as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})


def build_sample_loss_weights(
    rows: pd.DataFrame,
    *,
    class_weight_cap: float = 4.0,
) -> np.ndarray:
    """Return per-row morphology weights implementing the v1 loss policy.

    The inverse-square-root class correction is capped.  Real and injected
    Planet-like rows receive equal aggregate weight.  Injected examples in
    every other class are capped at 25% of that class's real aggregate.  Real
    flat and systematic rows contribute equal aggregate weight within Other.
    """

    if "morphology_target_v1" not in rows:
        raise KeyError("morphology_target_v1")
    labels = rows["morphology_target_v1"].fillna("").astype(str)
    include = (
        _as_bool(rows["morphology_include_v1"])
        if "morphology_include_v1" in rows
        else labels.isin(MORPHOLOGY_CLASSES)
    )
    counts = labels[include].value_counts()
    if counts.empty:
        return np.zeros(len(rows), dtype=np.float32)
    largest = float(counts.max())
    class_weight = {
        label: min(float(class_weight_cap), np.sqrt(largest / float(max(count, 1))))
        for label, count in counts.items()
    }
    weights = np.asarray(
        [class_weight.get(label, 0.0) if keep else 0.0 for label, keep in zip(labels, include)],
        dtype=np.float64,
    )
    injected = (
        _as_bool(rows["is_injected_row"]).to_numpy()
        if "is_injected_row" in rows
        else rows.get("source_kind", pd.Series("", index=rows.index)).astype(str).str.contains("inject").to_numpy()
    )

    def equalize(mask_a: np.ndarray, mask_b: np.ndarray) -> None:
        total_a = float(weights[mask_a].sum())
        total_b = float(weights[mask_b].sum())
        if total_a <= 0 or total_b <= 0:
            return
        target = 0.5 * (total_a + total_b)
        weights[mask_a] *= target / total_a
        weights[mask_b] *= target / total_b

    human_label = rows.get("human_label", pd.Series("", index=rows.index)).fillna("").astype(str)
    other = include.to_numpy() & labels.eq("other").to_numpy() & ~injected
    systematic = other & human_label.eq("instrumental_or_systematic").to_numpy()
    flat = other & human_label.isin({"uncertain", "no_visible_signal"}).to_numpy()
    equalize(systematic, flat)

    for label in MORPHOLOGY_CLASSES:
        in_class = include.to_numpy() & labels.eq(label).to_numpy()
        real_mask = in_class & ~injected
        injected_mask = in_class & injected
        real_total = float(weights[real_mask].sum())
        injected_total = float(weights[injected_mask].sum())
        if real_total <= 0 or injected_total <= 0:
            continue
        cap = real_total if label == "planet_like" else 0.25 * real_total
        if injected_total > cap:
            weights[injected_mask] *= cap / injected_total

    positive = weights > 0
    if np.any(positive):
        weights[positive] /= float(np.mean(weights[positive]))
    return weights.astype(np.float32)


def _balanced_group_assignment(
    rows: pd.DataFrame,
    *,
    label_column: str,
    group_column: str,
    n_folds: int,
    seed: int,
) -> pd.Series:
    """Greedily distribute whole groups while balancing label totals."""

    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")
    labels = rows[label_column].fillna("").astype(str)
    groups = rows[group_column].fillna("").astype(str)
    if groups.eq("").any():
        raise ValueError(f"{group_column} contains empty group identifiers")
    classes = sorted(labels.unique())
    table = pd.crosstab(groups, labels).reindex(columns=classes, fill_value=0)
    totals = table.sum(axis=0).to_numpy(dtype=float)
    target = totals / float(n_folds)
    target_size = float(len(rows)) / float(n_folds)
    rarity = (table.to_numpy(dtype=float) / np.maximum(totals, 1.0)).max(axis=1)
    rng = np.random.default_rng(seed)
    tie = rng.random(len(table))
    order = np.lexsort((tie, -table.sum(axis=1).to_numpy(dtype=float), -rarity))
    fold_counts = np.zeros((n_folds, len(classes)), dtype=float)
    fold_sizes = np.zeros(n_folds, dtype=float)
    assignment: dict[str, int] = {}
    for position in order:
        group = str(table.index[position])
        vector = table.iloc[position].to_numpy(dtype=float)
        size = float(vector.sum())
        costs: list[float] = []
        for fold in range(n_folds):
            candidate = fold_counts.copy()
            candidate[fold] += vector
            size_candidate = fold_sizes.copy()
            size_candidate[fold] += size
            class_cost = float(np.sum(np.square(candidate - target) / np.maximum(target, 1.0)))
            size_cost = float(np.sum(np.square(size_candidate - target_size) / max(target_size, 1.0)))
            costs.append(class_cost + 0.10 * size_cost)
        minimum = min(costs)
        choices = [index for index, value in enumerate(costs) if np.isclose(value, minimum)]
        fold = int(rng.choice(choices))
        assignment[group] = fold
        fold_counts[fold] += vector
        fold_sizes[fold] += size
    return groups.map(assignment).astype(np.int16)


def build_grouped_test_and_cv_folds(
    rows: pd.DataFrame,
    *,
    label_column: str = "split_stratum",
    group_column: str = "tic",
    seed: int = 56,
) -> pd.DataFrame:
    """Create a fixed 20% grouped test set and five grouped development folds."""

    work = rows.copy()
    if label_column not in work:
        morphology = work.get("morphology_target_v1", pd.Series("", index=work.index)).fillna("").astype(str)
        preserve_only = work.get("broad_preserve_only", pd.Series(False, index=work.index))
        preserve_only = _as_bool(preserve_only) if isinstance(preserve_only, pd.Series) else preserve_only
        work[label_column] = morphology.where(morphology.ne(""), np.where(preserve_only, "broad_preserve", "auxiliary"))
    outer = _balanced_group_assignment(
        work,
        label_column=label_column,
        group_column=group_column,
        n_folds=5,
        seed=seed,
    )
    test_mask = outer.eq(0)
    result = pd.DataFrame(index=work.index)
    result["fixed_split"] = np.where(test_mask, "test", "development")
    result["cv_fold"] = -1
    development = work.loc[~test_mask]
    if not development.empty:
        dev_folds = _balanced_group_assignment(
            development,
            label_column=label_column,
            group_column=group_column,
            n_folds=5,
            seed=seed + 1,
        )
        result.loc[development.index, "cv_fold"] = dev_folds.astype(int)
    result["cv_fold"] = result["cv_fold"].astype(np.int16)
    return result


def build_harmonic_cnn(
    config: HarmonicModelConfig,
    *,
    profile: str = "full_combined",
) -> Any:
    """Construct the PyTorch multi-branch teacher lazily."""

    branches = profile_branches(profile)
    import torch
    from torch import nn

    embedding_dim = int(config.embedding_dim)

    class ChannelLayerNorm(nn.Module):
        """Normalize channels independently at each cadence, before masking."""

        def __init__(self, channels: int) -> None:
            super().__init__()
            self.norm = nn.LayerNorm(channels)

        def forward(self, values: Any) -> Any:
            return self.norm(values.transpose(1, 2)).transpose(1, 2)

    class ResidualBlock(nn.Module):
        def __init__(self, channels: int, dilation: int) -> None:
            super().__init__()
            padding = dilation * 3
            self.net = nn.Sequential(
                nn.Conv1d(channels, channels, 7, padding=padding, dilation=dilation),
                ChannelLayerNorm(channels),
                nn.GELU(),
                nn.Conv1d(channels, channels, 5, padding=2),
                ChannelLayerNorm(channels),
            )
            self.activation = nn.GELU()

        def forward(self, values: Any, sample_mask: Any) -> Any:
            return self.activation(values + self.net(values) * sample_mask) * sample_mask

    class BranchDropout(nn.Module):
        """Drop an entire supplemental-aperture embedding for selected rows."""

        def __init__(self, probability: float) -> None:
            super().__init__()
            self.probability = float(probability)

        def forward(self, values: Any) -> Any:
            if not self.training or self.probability <= 0:
                return values
            keep_probability = 1.0 - self.probability
            keep = (
                torch.rand((values.shape[0], 1), device=values.device) < keep_probability
            ).to(dtype=values.dtype)
            return values * keep / max(keep_probability, 1.0e-6)

    class MaskedSequenceEncoder(nn.Module):
        def __init__(self, n_channels: int, hidden: int, output: int) -> None:
            super().__init__()
            self.stem = nn.Sequential(
                nn.Conv1d(n_channels, hidden, 9, padding=4),
                ChannelLayerNorm(hidden),
                nn.GELU(),
            )
            self.blocks = nn.ModuleList([ResidualBlock(hidden, 1), ResidualBlock(hidden, 2)])
            self.attention = nn.Conv1d(hidden, 1, 1)
            self.output = nn.Linear(2 * hidden, output)

        def forward(self, values: Any, mask: Any) -> Any:
            sample_mask = mask.any(dim=1, keepdim=True).to(dtype=values.dtype)
            hidden = self.stem(values * mask.to(dtype=values.dtype)) * sample_mask
            for block in self.blocks:
                hidden = block(hidden, sample_mask)
            scores = self.attention(hidden).masked_fill(~sample_mask.bool(), -1.0e4)
            attention = torch.softmax(scores, dim=-1) * sample_mask
            attention = attention / attention.sum(dim=-1, keepdim=True).clamp_min(1.0e-6)
            attended = (hidden * attention).sum(dim=-1)
            maximum = hidden.masked_fill(~sample_mask.bool(), -1.0e4).amax(dim=-1)
            has_data = sample_mask.any(dim=-1)
            maximum = torch.where(has_data, maximum, torch.zeros_like(maximum))
            embedding = self.output(torch.cat([attended, maximum], dim=1))
            return embedding * has_data.to(dtype=embedding.dtype)

    class HarmonicSetEncoder(nn.Module):
        def __init__(self, sequence_encoder: Any, output: int) -> None:
            super().__init__()
            self.sequence_encoder = sequence_encoder
            self.factor_attention = nn.Linear(output, 1)
            self.output = nn.Linear(2 * output, output)

        def forward(self, values: Any, mask: Any, *, single_period: bool = False) -> Any:
            if single_period:
                period_index = 3 if values.shape[1] > 3 else 0
                values = values[:, period_index : period_index + 1]
                mask = mask[:, period_index : period_index + 1]
            batch, views, channels, samples = values.shape
            encoded = self.sequence_encoder(
                values.reshape(batch * views, channels, samples),
                mask.reshape(batch * views, channels, samples),
            ).reshape(batch, views, -1)
            view_mask = mask.any(dim=(2, 3))
            scores = self.factor_attention(encoded).squeeze(-1).masked_fill(~view_mask, -1.0e4)
            attention = torch.softmax(scores, dim=1) * view_mask.to(dtype=encoded.dtype)
            attention = attention / attention.sum(dim=1, keepdim=True).clamp_min(1.0e-6)
            attended = (encoded * attention.unsqueeze(-1)).sum(dim=1)
            maximum = encoded.masked_fill(~view_mask.unsqueeze(-1), -1.0e4).amax(dim=1)
            has_data = view_mask.any(dim=1, keepdim=True)
            maximum = torch.where(has_data, maximum, torch.zeros_like(maximum))
            embedding = self.output(torch.cat([attended, maximum], dim=1))
            return embedding * has_data.to(dtype=embedding.dtype)

    class HarmonicCNN(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.model_version = MODEL_VERSION
            self.profile = profile
            self.config = asdict(config)
            self.small_encoder = MaskedSequenceEncoder(config.small_channels, 64, embedding_dim)
            self.supp_encoder = MaskedSequenceEncoder(config.supplemental_channels, 32, embedding_dim)
            fold_sequence = MaskedSequenceEncoder(config.fold_channels, 48, embedding_dim)
            local_sequence = MaskedSequenceEncoder(config.fold_channels, 48, embedding_dim)
            self.harmonic_encoder = HarmonicSetEncoder(fold_sequence, embedding_dim)
            self.local_encoder = HarmonicSetEncoder(local_sequence, embedding_dim)
            self.periodogram_encoder = MaskedSequenceEncoder(
                config.periodogram_channels, 32, embedding_dim
            )
            self.metadata_encoder = (
                nn.Sequential(
                    nn.Linear(config.metadata_dim, embedding_dim),
                    nn.GELU(),
                    nn.LayerNorm(embedding_dim),
                )
                if config.metadata_dim > 0
                else None
            )
            self.supplement_gate = nn.Sequential(
                nn.Linear(2 * embedding_dim, embedding_dim), nn.Sigmoid()
            )
            self.branch_gate = nn.Sequential(nn.Linear(embedding_dim, embedding_dim), nn.Sigmoid())
            self.supplement_dropout = BranchDropout(config.supplemental_dropout)
            self.fusion = nn.Sequential(
                nn.Linear(5 * embedding_dim, 2 * embedding_dim),
                nn.GELU(),
                nn.LayerNorm(2 * embedding_dim),
                nn.Dropout(config.dropout),
            )
            self.morphology_head = nn.Linear(2 * embedding_dim, len(MORPHOLOGY_CLASSES))
            self.preserve_head = nn.Linear(2 * embedding_dim, len(PRESERVE_CLASSES))
            self.harmonic_head = nn.Linear(2 * embedding_dim, len(HARMONIC_CLASSES))

        def _gate(self, value: Any) -> Any:
            return value * self.branch_gate(value)

        def forward(self, batch: Mapping[str, Any]) -> dict[str, Any]:
            batch_size = int(batch["harmonic_values"].shape[0])
            zero = batch["harmonic_values"].new_zeros((batch_size, embedding_dim))

            if "chronology" in branches:
                small = self.small_encoder(batch["chronology_small"], batch["chronology_small_mask"])
                supplemental = self.supp_encoder(
                    batch["chronology_supplemental"], batch["chronology_supplemental_mask"]
                )
                supplemental = self.supplement_dropout(supplemental)
                supplemental = supplemental * self.supplement_gate(torch.cat([small, supplemental], dim=1))
                chronology = small + supplemental
            else:
                chronology = zero

            harmonic = self.harmonic_encoder(
                batch["harmonic_values"],
                batch["harmonic_mask"],
                single_period="single_fold" in branches,
            ) if ({"harmonic", "single_fold"} & branches) else zero
            if "local" in branches:
                local_values = batch["local_values"]
                local_mask = batch["local_mask"]
                if "single_fold" in branches and local_values.shape[1] >= 11:
                    local_values = torch.cat(
                        [local_values[:, 3:4], local_values[:, 10:11]], dim=1
                    )
                    local_mask = torch.cat(
                        [local_mask[:, 3:4], local_mask[:, 10:11]], dim=1
                    )
                local = self.local_encoder(local_values, local_mask)
            else:
                local = zero
            periodogram = self.periodogram_encoder(
                batch["periodogram_values"], batch["periodogram_mask"]
            ) if "periodogram" in branches else zero
            metadata = (
                self.metadata_encoder(batch["metadata"])
                if "metadata" in branches and self.metadata_encoder is not None
                else zero
            )
            fused = self.fusion(
                torch.cat(
                    [
                        self._gate(chronology),
                        self._gate(harmonic),
                        self._gate(local),
                        self._gate(periodogram),
                        self._gate(metadata),
                    ],
                    dim=1,
                )
            )
            return {
                "embedding": fused,
                "morphology_logits": self.morphology_head(fused),
                "preserve_logits": self.preserve_head(fused),
                "harmonic_logits": self.harmonic_head(fused),
            }

    return HarmonicCNN()


def multitask_loss(
    outputs: Mapping[str, Any],
    batch: Mapping[str, Any],
    *,
    config: HarmonicTrainConfig = HarmonicTrainConfig(),
) -> tuple[Any, dict[str, Any]]:
    """Compute masked, weighted morphology/preserve/harmonic losses."""

    import torch
    from torch.nn import functional as functional

    zero = outputs["morphology_logits"].sum() * 0.0

    def task_loss(logit_name: str, target_name: str, weight_name: str) -> Any:
        target = batch[target_name].long()
        active = target >= 0
        if not torch.any(active):
            return zero
        loss = functional.cross_entropy(outputs[logit_name][active], target[active], reduction="none")
        weight = batch.get(weight_name)
        if weight is None:
            return loss.mean()
        selected = weight[active].to(dtype=loss.dtype)
        return (loss * selected).sum() / selected.sum().clamp_min(1.0e-6)

    morphology = task_loss("morphology_logits", "morphology_target", "morphology_weight")
    preserve = task_loss("preserve_logits", "preserve_target", "preserve_weight")
    harmonic = task_loss("harmonic_logits", "harmonic_target", "harmonic_weight")
    total = (
        config.morphology_loss_weight * morphology
        + config.preserve_loss_weight * preserve
        + config.harmonic_loss_weight * harmonic
    )
    return total, {"morphology": morphology, "preserve": preserve, "harmonic": harmonic}


def probability_outputs(outputs: Mapping[str, Any], *, period_d: Any) -> dict[str, Any]:
    """Convert head logits into calibrated-output-shaped tensors."""

    import torch

    morphology = torch.softmax(outputs["morphology_logits"], dim=1)
    preserve = torch.softmax(outputs["preserve_logits"], dim=1)
    harmonic = torch.softmax(outputs["harmonic_logits"], dim=1)
    factor_values = torch.tensor(
        [0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0],
        dtype=harmonic.dtype,
        device=harmonic.device,
    )
    factor_index = harmonic.argmax(dim=1)
    factor = factor_values[factor_index]
    return {
        "morphology_probability": morphology,
        "p_preserve": preserve[:, 1],
        "harmonic_probability": harmonic,
        "predicted_period_factor": factor,
        "predicted_period_d": period_d * factor,
    }


__all__ = [
    "HARMONIC_CLASSES",
    "HarmonicModelConfig",
    "HarmonicTrainConfig",
    "MODEL_PROFILES",
    "MODEL_VERSION",
    "MORPHOLOGY_CLASSES",
    "PRESERVE_CLASSES",
    "build_grouped_test_and_cv_folds",
    "build_harmonic_cnn",
    "build_sample_loss_weights",
    "multitask_loss",
    "probability_outputs",
    "profile_branches",
]
