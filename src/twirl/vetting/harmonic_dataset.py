"""Dataset and batching helpers for native-cadence S56 harmonic CNN inputs."""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.adjudication_audit import HARMONIC_CNN_TARGET_POLICY
from twirl.vetting.harmonic_cnn import (
    HARMONIC_CLASSES,
    MORPHOLOGY_CLASSES,
    PRESERVE_CLASSES,
    build_grouped_test_and_cv_folds,
    build_sample_loss_weights,
    profile_branches,
)
from twirl.vetting.harmonic_inputs import (
    HARMONIC_FACTORS,
    build_harmonic_views,
    build_native_channels,
    native_group_path,
    pad_channel_sequences,
    read_native_light_curve,
)
from twirl.vetting.recovery50_teacher import leakage_columns


HARMONIC_METADATA_CANDIDATES: tuple[str, ...] = (
    "period_d",
    "duration_min",
    "depth",
    "depth_snr",
    "sde_max",
    "tmag",
    "adp_sml_period_d",
    "adp_sml_duration_min",
    "adp_sml_depth",
    "adp_sml_depth_snr",
    "adp_sml_sde",
    "adp_sml_log_power",
    "adp_sml_own_even_depth",
    "adp_sml_own_odd_depth",
    "adp_sml_own_even_odd_depth_delta",
    "adp_sml_own_even_odd_sigma_delta",
    "adp_sml_trend_ptp",
    "adp_period_d",
    "adp_duration_min",
    "adp_depth",
    "adp_depth_snr",
    "adp_sde",
    "adp_log_power",
    "adp_own_even_depth",
    "adp_own_odd_depth",
    "adp_own_even_odd_depth_delta",
    "adp_own_even_odd_sigma_delta",
    "adp_trend_ptp",
    "aperture_period_rel_delta",
    "aperture_depth_ratio_primary_over_small",
    "aperture_disagreement_flag",
)


@dataclass(frozen=True)
class MetadataNormalization:
    columns: tuple[str, ...]
    center: tuple[float, ...]
    scale: tuple[float, ...]


def candidate_bls_ephemeris(row: Mapping[str, Any]) -> tuple[float, float, float]:
    """Return the untouched queue/ADP-small BLS ephemeris used for all views."""

    def finite_value(*names: str) -> float:
        for name in names:
            try:
                value = float(row.get(name, np.nan))
            except (TypeError, ValueError):
                continue
            if np.isfinite(value):
                return value
        return float("nan")

    # ``period_d`` is the solid-red review period.  Explicitly avoid every
    # model/effective/adjudicated period field here; those are supervision.
    period = finite_value("period_d", "adp_sml_period_d", "display_period_d")
    t0 = finite_value("t0_bjd", "adp_sml_t0_bjd", "display_t0_bjd")
    duration = finite_value("duration_min", "adp_sml_duration_min", "display_duration_min")
    if not np.isfinite(period) or period <= 0:
        raise ValueError("row has no finite positive original BLS period")
    if not np.isfinite(t0):
        raise ValueError("row has no finite original BLS epoch")
    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("row has no finite positive BLS duration")
    return period, t0, duration


def _target_index(value: Any, classes: Sequence[str], include: Any) -> int:
    active = (
        bool(include)
        if isinstance(include, (bool, np.bool_))
        else str(include).strip().lower() in {"1", "true", "t", "yes", "y"}
    )
    if not active:
        return -1
    text = str(value)
    return classes.index(text) if text in classes else -1


def _inverse_sqrt_task_weights(target: np.ndarray, *, cap: float = 4.0) -> np.ndarray:
    out = np.zeros(len(target), dtype=np.float32)
    active = target >= 0
    if not np.any(active):
        return out
    counts = np.bincount(target[active])
    largest = float(counts.max(initial=1))
    for value, count in enumerate(counts):
        if count:
            out[target == value] = min(float(cap), np.sqrt(largest / float(count)))
    out[active] /= float(np.mean(out[active]))
    return out


def _cap_injected_task_weights(
    rows: pd.DataFrame,
    *,
    target_column: str,
    weights: np.ndarray,
) -> np.ndarray:
    """Apply source caps to an auxiliary task without changing its targets."""

    out = np.asarray(weights, dtype=np.float64).copy()
    target = rows[target_column].to_numpy(dtype=int)
    injected = rows.get("is_injected_row", pd.Series(False, index=rows.index))
    if injected.dtype != bool:
        injected = injected.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})
    injected = injected.to_numpy()
    planet = rows["human_label"].fillna("").astype(str).eq("planet_like").to_numpy()
    for value in np.unique(target[target >= 0]):
        in_target = target == value
        real_total = float(out[in_target & ~injected].sum())
        if real_total <= 0:
            continue
        injected_planet = in_target & injected & planet
        injected_nonplanet = in_target & injected & ~planet
        for mask, cap in (
            (injected_planet, real_total),
            (injected_nonplanet, 0.25 * real_total),
        ):
            total = float(out[mask].sum())
            if total > cap and cap >= 0:
                out[mask] *= cap / total
    active = out > 0
    if np.any(active):
        out[active] /= float(np.mean(out[active]))
    return out.astype(np.float32)


def prepare_harmonic_training_rows(rows: pd.DataFrame, *, seed: int = 56) -> pd.DataFrame:
    """Validate targets and attach native groups, task indices, weights, splits."""

    work = rows.copy()
    if "model_target_policy_version" not in work:
        raise KeyError("model_target_policy_version")
    policy = work["model_target_policy_version"].fillna("").astype(str)
    if not policy.eq(HARMONIC_CNN_TARGET_POLICY).all():
        bad = sorted(policy[policy.ne(HARMONIC_CNN_TARGET_POLICY)].unique())
        raise ValueError(f"training rows contain unexpected target policies: {bad}")
    work["morphology_target_index"] = [
        _target_index(value, MORPHOLOGY_CLASSES, include)
        for value, include in zip(work["morphology_target_v1"], work["morphology_include_v1"])
    ]
    work["preserve_target_index"] = [
        _target_index(value, PRESERVE_CLASSES, include)
        for value, include in zip(work["preserve_target_v1"], work["preserve_include_v1"])
    ]
    work["harmonic_target_index"] = [
        _target_index(value, HARMONIC_CLASSES, include)
        for value, include in zip(work["harmonic_target_v1"], work["harmonic_include_v1"])
    ]
    active = (
        work["morphology_target_index"].ge(0)
        | work["preserve_target_index"].ge(0)
        | work["harmonic_target_index"].ge(0)
    )
    work = work.loc[active].copy().reset_index(drop=True)
    work["native_group_path"] = [native_group_path(row) for row in work.to_dict("records")]
    # These placeholders are replaced from each fold's training partition.
    # Held-out label frequencies therefore cannot influence optimization.
    work["morphology_weight"] = work["morphology_target_index"].ge(0).astype(np.float32)
    work["preserve_weight"] = work["preserve_target_index"].ge(0).astype(np.float32)
    work["harmonic_weight"] = work["harmonic_target_index"].ge(0).astype(np.float32)
    morphology = work["morphology_target_v1"].fillna("").astype(str)
    broad = work["broad_preserve_only"]
    if broad.dtype != bool:
        broad = broad.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})
    work["split_stratum"] = morphology.where(
        morphology.ne(""), np.where(broad, "broad_preserve", "auxiliary")
    )
    split = build_grouped_test_and_cv_folds(work, seed=seed)
    work["fixed_split"] = split["fixed_split"].to_numpy()
    work["cv_fold"] = split["cv_fold"].to_numpy()
    return work


def attach_fold_training_weights(
    rows: pd.DataFrame,
    *,
    fit_mask: np.ndarray,
) -> pd.DataFrame:
    """Fit class and source weights from one fold's training rows only."""

    fit_mask = np.asarray(fit_mask, dtype=bool)
    if fit_mask.shape != (len(rows),):
        raise ValueError("fit_mask must have one value per row")
    out = rows.copy()
    for column in ("morphology_weight", "preserve_weight", "harmonic_weight"):
        out[column] = np.float32(0.0)
    fit = out.loc[fit_mask].copy()
    morphology_weight = build_sample_loss_weights(fit)
    preserve_weight = _cap_injected_task_weights(
        fit,
        target_column="preserve_target_index",
        weights=_inverse_sqrt_task_weights(fit["preserve_target_index"].to_numpy()),
    )
    harmonic_weight = _cap_injected_task_weights(
        fit,
        target_column="harmonic_target_index",
        weights=_inverse_sqrt_task_weights(fit["harmonic_target_index"].to_numpy()),
    )
    out.loc[fit_mask, "morphology_weight"] = morphology_weight
    out.loc[fit_mask, "preserve_weight"] = preserve_weight
    out.loc[fit_mask, "harmonic_weight"] = harmonic_weight
    return out


def build_metadata_matrix(
    rows: pd.DataFrame,
    *,
    fit_mask: np.ndarray,
    normalization: MetadataNormalization | None = None,
) -> tuple[np.ndarray, MetadataNormalization]:
    """Build a leakage-audited scalar BLS matrix with train-only scaling."""

    if normalization is None:
        columns = tuple(column for column in HARMONIC_METADATA_CANDIDATES if column in rows)
    else:
        columns = normalization.columns
    leaks = leakage_columns(columns)
    if leaks:
        raise ValueError(f"harmonic metadata selection contains leakage: {leaks}")
    values = np.column_stack(
        [pd.to_numeric(rows.get(column, np.nan), errors="coerce").to_numpy(dtype=float) for column in columns]
    ) if columns else np.empty((len(rows), 0), dtype=float)
    fit_mask = np.asarray(fit_mask, dtype=bool)
    if fit_mask.shape != (len(rows),):
        raise ValueError("fit_mask must have one value per row")
    if normalization is None:
        center = np.nanmedian(values[fit_mask], axis=0) if values.shape[1] else np.empty(0)
        scale = np.nanstd(values[fit_mask], axis=0) if values.shape[1] else np.empty(0)
        center = np.where(np.isfinite(center), center, 0.0)
        scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
        normalization = MetadataNormalization(
            columns=columns,
            center=tuple(float(value) for value in center),
            scale=tuple(float(value) for value in scale),
        )
    center = np.asarray(normalization.center, dtype=float)
    scale = np.asarray(normalization.scale, dtype=float)
    values = np.where(np.isfinite(values), values, center)
    values = (values - center) / scale
    return values.astype(np.float32), normalization


def build_injection_pretraining_rows(
    rows: pd.DataFrame,
    *,
    seed: int = 56,
) -> pd.DataFrame:
    """Build visible-injection, paired-original, and real-negative pretraining rows."""

    work = rows.copy()
    injected = work.get("is_injected_row", pd.Series(False, index=work.index))
    if injected.dtype != bool:
        injected = injected.fillna("").astype(str).str.lower().isin({"1", "true", "t", "yes", "y"})
    visible = (
        injected
        & work["human_label"].fillna("").astype(str).eq("planet_like")
        & work["morphology_target_v1"].fillna("").astype(str).eq("planet_like")
    )
    positives = work.loc[visible].copy()
    if positives.empty:
        raise ValueError("no human-visible injected Planet-like rows for encoder pretraining")
    positives["input_variant"] = "injected"
    positives["pretrain_target"] = 1
    originals = positives.copy()
    originals["input_variant"] = "paired_original"
    originals["pretrain_target"] = 0
    originals["review_id"] = originals["review_id"].astype(str) + ":paired_original"

    real_negative = work.loc[
        ~injected
        & work["human_label"].fillna("").astype(str).isin(
            {"instrumental_or_systematic", "uncertain", "no_visible_signal"}
        )
    ].copy()
    rng = np.random.default_rng(seed)
    selected: list[pd.DataFrame] = []
    desired = len(positives)
    labels = ("instrumental_or_systematic", "uncertain")
    quotas = {labels[0]: desired // 2, labels[1]: desired - desired // 2}
    for label in labels:
        candidates = real_negative[
            real_negative["human_label"].astype(str).eq(label)
            | (label == "uncertain")
            & real_negative["human_label"].astype(str).eq("no_visible_signal")
        ]
        if len(candidates) > quotas[label]:
            chosen = rng.choice(candidates.index.to_numpy(), size=quotas[label], replace=False)
            candidates = candidates.loc[chosen]
        selected.append(candidates)
    negatives = pd.concat(selected, ignore_index=False)
    if len(negatives) < desired:
        remaining = real_negative.drop(index=negatives.index, errors="ignore")
        take = min(desired - len(negatives), len(remaining))
        if take:
            chosen = rng.choice(remaining.index.to_numpy(), size=take, replace=False)
            negatives = pd.concat([negatives, remaining.loc[chosen]], ignore_index=False)
    if len(negatives) != desired:
        raise ValueError(
            f"need {desired} real flat/systematic pretraining negatives; found {len(negatives)}"
        )
    negatives = negatives.copy()
    negatives["input_variant"] = "observed"
    negatives["pretrain_target"] = 0
    result = pd.concat([positives, originals, negatives], ignore_index=True)
    return result.sample(frac=1.0, random_state=seed).reset_index(drop=True)


class HarmonicNativeDataset:
    """Read native HDF5 rows and generate unbinned views on demand."""

    def __init__(
        self,
        rows: pd.DataFrame,
        *,
        native_h5: Path,
        metadata: np.ndarray,
        cache_size: int = 2048,
        profile: str = "full_combined",
    ) -> None:
        self.rows = rows.reset_index(drop=True).copy()
        self.native_h5 = Path(native_h5)
        self.metadata = np.asarray(metadata, dtype=np.float32)
        if self.metadata.shape[0] != len(self.rows):
            raise ValueError("metadata row count does not match training rows")
        self.cache_size = max(0, int(cache_size))
        self._cache: OrderedDict[int, dict[str, Any]] = OrderedDict()
        self.branches = profile_branches(profile)

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, index: int) -> dict[str, Any]:
        index = int(index)
        if index in self._cache:
            self._cache.move_to_end(index)
            return self._cache[index]
        row = self.rows.iloc[index]
        target_payload = {
            "review_id": str(row.get("review_id", "")),
            "tic": int(float(row["tic"])),
            "period_d": np.float32(candidate_bls_ephemeris(row)[0]),
            "metadata": self.metadata[index],
            "morphology_target": int(row["morphology_target_index"]),
            "preserve_target": int(row["preserve_target_index"]),
            "harmonic_target": int(row["harmonic_target_index"]),
            "morphology_weight": np.float32(row["morphology_weight"]),
            "preserve_weight": np.float32(row["preserve_weight"]),
            "harmonic_weight": np.float32(row["harmonic_weight"]),
            "pretrain_target": int(row.get("pretrain_target", -1)),
        }
        if self.branches == {"metadata"}:
            sample = {
                **target_payload,
                "chronology_small": np.zeros((10, 1), dtype=np.float32),
                "chronology_small_mask": np.zeros((10, 1), dtype=bool),
                "chronology_supplemental": np.zeros((5, 1), dtype=np.float32),
                "chronology_supplemental_mask": np.zeros((5, 1), dtype=bool),
                "harmonic_values": tuple(np.zeros((7, 1), dtype=np.float32) for _ in range(7)),
                "harmonic_mask": tuple(np.zeros((7, 1), dtype=bool) for _ in range(7)),
                "local_values": tuple(np.zeros((7, 1), dtype=np.float32) for _ in range(14)),
                "local_mask": tuple(np.zeros((7, 1), dtype=bool) for _ in range(14)),
                "periodogram_values": np.zeros((4, 1), dtype=np.float32),
                "periodogram_mask": np.zeros((4, 1), dtype=bool),
            }
            if self.cache_size:
                self._cache[index] = sample
            return sample
        lc = read_native_light_curve(
            self.native_h5,
            group_path=str(row["native_group_path"]),
            require_errors=True,
        )
        variant = str(row.get("input_variant", "observed"))
        if variant == "paired_original":
            if lc.paired_original is None:
                raise ValueError(f"{row['native_group_path']} has no paired original light curve")
            lc = lc.paired_original
        elif variant not in {"observed", "injected"}:
            raise ValueError(f"unknown native input_variant {variant!r}")
        if float(np.nanmedian(lc.time)) < 1.0e5:
            raise ValueError("native input time must be absolute BJD, not BTJD")
        period, t0, duration = candidate_bls_ephemeris(row)
        chronology = (
            build_native_channels(lc)
            if "chronology" in self.branches
            else None
        )
        harmonic = build_harmonic_views(
            lc,
            period_d=period,
            t0_bjd=t0,
            duration_min=duration,
            factors=(1.0,) if "single_fold" in self.branches else HARMONIC_FACTORS,
        )
        if (
            len(lc.bls_power_small)
            and len(lc.bls_sde_small)
            and len(lc.bls_power_primary)
            and len(lc.bls_sde_primary)
        ):
            periodogram = np.stack(
                [
                    lc.bls_power_small,
                    lc.bls_sde_small,
                    lc.bls_power_primary,
                    lc.bls_sde_primary,
                ],
                axis=0,
            ).astype(np.float32)
            for channel in range(periodogram.shape[0]):
                values = periodogram[channel]
                finite = values[np.isfinite(values)]
                if finite.size:
                    center = float(np.nanmedian(finite))
                    scale = 1.4826 * float(np.nanmedian(np.abs(finite - center)))
                    if not np.isfinite(scale) or scale <= 1.0e-8:
                        scale = float(np.nanstd(finite))
                    if not np.isfinite(scale) or scale <= 1.0e-8:
                        scale = 1.0
                    periodogram[channel] = (values - center) / scale
        else:
            periodogram = np.full((4, 1), np.nan, dtype=np.float32)
        sample = {
            **target_payload,
            "period_d": np.float32(period),
            "chronology_small": (
                chronology.small_values
                if chronology is not None
                else np.zeros((10, 1), dtype=np.float32)
            ),
            "chronology_small_mask": (
                chronology.small_mask
                if chronology is not None
                else np.zeros((10, 1), dtype=bool)
            ),
            "chronology_supplemental": (
                chronology.supplemental_values
                if chronology is not None
                else np.zeros((5, 1), dtype=np.float32)
            ),
            "chronology_supplemental_mask": (
                chronology.supplemental_mask
                if chronology is not None
                else np.zeros((5, 1), dtype=bool)
            ),
            "harmonic_values": harmonic.full_values,
            "harmonic_mask": harmonic.full_masks,
            "local_values": harmonic.primary_values + harmonic.secondary_values,
            "local_mask": harmonic.primary_masks + harmonic.secondary_masks,
            "periodogram_values": periodogram,
            "periodogram_mask": np.isfinite(periodogram),
        }
        if self.cache_size:
            self._cache[index] = sample
            self._cache.move_to_end(index)
            while len(self._cache) > self.cache_size:
                self._cache.popitem(last=False)
        return sample


def _pad_nested_views(samples: Sequence[dict[str, Any]], key: str, mask_key: str) -> tuple[np.ndarray, np.ndarray]:
    n_views = len(samples[0][key])
    flattened = [value for sample in samples for value in sample[key]]
    flattened_masks = [value for sample in samples for value in sample[mask_key]]
    values, masks, _ = pad_channel_sequences(flattened, flattened_masks)
    return (
        values.reshape(len(samples), n_views, values.shape[1], values.shape[2]),
        masks.reshape(len(samples), n_views, masks.shape[1], masks.shape[2]),
    )


def collate_native_samples(samples: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Pad only to the longest sequence in this batch; never truncate."""

    import torch

    if not samples:
        raise ValueError("cannot collate an empty native batch")
    small, small_mask, _ = pad_channel_sequences(
        [sample["chronology_small"] for sample in samples],
        [sample["chronology_small_mask"] for sample in samples],
    )
    supplemental, supplemental_mask, _ = pad_channel_sequences(
        [sample["chronology_supplemental"] for sample in samples],
        [sample["chronology_supplemental_mask"] for sample in samples],
    )
    harmonic, harmonic_mask = _pad_nested_views(samples, "harmonic_values", "harmonic_mask")
    local, local_mask = _pad_nested_views(samples, "local_values", "local_mask")
    periodogram, periodogram_mask, _ = pad_channel_sequences(
        [sample["periodogram_values"] for sample in samples],
        [sample["periodogram_mask"] for sample in samples],
    )

    def tensor(name: str, dtype: Any) -> Any:
        return torch.as_tensor(np.asarray([sample[name] for sample in samples]), dtype=dtype)

    return {
        "review_id": [sample["review_id"] for sample in samples],
        "tic": tensor("tic", torch.long),
        "period_d": tensor("period_d", torch.float32),
        "chronology_small": torch.from_numpy(small),
        "chronology_small_mask": torch.from_numpy(small_mask),
        "chronology_supplemental": torch.from_numpy(supplemental),
        "chronology_supplemental_mask": torch.from_numpy(supplemental_mask),
        "harmonic_values": torch.from_numpy(harmonic),
        "harmonic_mask": torch.from_numpy(harmonic_mask),
        "local_values": torch.from_numpy(local),
        "local_mask": torch.from_numpy(local_mask),
        "periodogram_values": torch.from_numpy(periodogram),
        "periodogram_mask": torch.from_numpy(periodogram_mask),
        "metadata": torch.from_numpy(np.stack([sample["metadata"] for sample in samples])),
        "morphology_target": tensor("morphology_target", torch.long),
        "preserve_target": tensor("preserve_target", torch.long),
        "harmonic_target": tensor("harmonic_target", torch.long),
        "morphology_weight": tensor("morphology_weight", torch.float32),
        "preserve_weight": tensor("preserve_weight", torch.float32),
        "harmonic_weight": tensor("harmonic_weight", torch.float32),
        "pretrain_target": tensor("pretrain_target", torch.long),
    }


__all__ = [
    "HARMONIC_METADATA_CANDIDATES",
    "HarmonicNativeDataset",
    "MetadataNormalization",
    "attach_fold_training_weights",
    "build_metadata_matrix",
    "build_injection_pretraining_rows",
    "candidate_bls_ephemeris",
    "collate_native_samples",
    "native_group_path",
    "prepare_harmonic_training_rows",
]
