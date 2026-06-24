"""TWIRL vetting utilities for Stage 2/5 candidates."""

from twirl.vetting.self_training import (
    SelfTrainingConfig,
    train_teacher_student,
)

__all__ = ["SelfTrainingConfig", "train_teacher_student"]
