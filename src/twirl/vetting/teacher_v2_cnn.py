"""Teacher-v2 CNN with a dedicated injection-trained compact-transit head."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

from twirl.vetting.harmonic_cnn import (
    HarmonicModelConfig,
    HarmonicTrainConfig,
    build_harmonic_cnn,
    multitask_loss,
)
from twirl.vetting.teacher_v2 import COMPACT_CLASSES, TEACHER_V2_MODEL_VERSION


@dataclass(frozen=True)
class TeacherV2TrainConfig(HarmonicTrainConfig):
    compact_loss_weight: float = 1.0


def build_teacher_v2_cnn(
    config: HarmonicModelConfig,
    *,
    profile: str = "full_combined",
) -> Any:
    """Wrap the audited v1 encoder/heads with an independent compact head."""

    from torch import nn

    backbone = build_harmonic_cnn(config, profile=profile)

    class TeacherV2CNN(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.model_version = TEACHER_V2_MODEL_VERSION
            self.profile = profile
            self.config = asdict(config)
            self.backbone = backbone
            self.compact_head = nn.Linear(
                2 * int(config.embedding_dim), len(COMPACT_CLASSES)
            )

        def forward(self, batch: Mapping[str, Any]) -> dict[str, Any]:
            output = self.backbone(batch)
            output["compact_logits"] = self.compact_head(output["embedding"])
            return output

    return TeacherV2CNN()


def teacher_v2_multitask_loss(
    outputs: Mapping[str, Any],
    batch: Mapping[str, Any],
    *,
    config: TeacherV2TrainConfig = TeacherV2TrainConfig(),
) -> tuple[Any, dict[str, Any]]:
    """Combine real-human morphology losses with injection compact loss."""

    from torch.nn import functional as functional

    base_config = HarmonicTrainConfig(
        epochs=config.epochs,
        batch_size=config.batch_size,
        learning_rate=config.learning_rate,
        weight_decay=config.weight_decay,
        patience=config.patience,
        seed=config.seed,
        morphology_loss_weight=config.morphology_loss_weight,
        preserve_loss_weight=config.preserve_loss_weight,
        harmonic_loss_weight=config.harmonic_loss_weight,
        class_weight_cap=config.class_weight_cap,
    )
    base, parts = multitask_loss(outputs, batch, config=base_config)
    target = batch["compact_target"].long()
    active = target >= 0
    if active.any():
        loss = functional.cross_entropy(
            outputs["compact_logits"][active], target[active], reduction="none"
        )
        weight = batch["compact_weight"][active].to(dtype=loss.dtype)
        compact = (loss * weight).sum() / weight.sum().clamp_min(1.0e-6)
    else:
        compact = outputs["compact_logits"].sum() * 0.0
    total = base + float(config.compact_loss_weight) * compact
    return total, {**parts, "compact": compact}


__all__ = [
    "TeacherV2TrainConfig",
    "build_teacher_v2_cnn",
    "teacher_v2_multitask_loss",
]
