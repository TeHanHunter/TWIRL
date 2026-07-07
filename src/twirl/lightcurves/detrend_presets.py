"""Named TWIRL-FS detrending presets.

Keep production method strings and FITS compare-column naming in one place so
S56 validation branches can graduate into the Stage 1 product without changing
their meaning.
"""
from __future__ import annotations

from .flux_detrend import FluxDetrendConfig


TWIRL_FS_V2_METHOD = "twirl-fs-v2"
TWIRL_FS_V2_ADP03Q_METHOD = "twirl-fs-v2-adp03q"
TWIRL_FS_V2_ADP015Q_METHOD = "twirl-fs-v2-adp015q"

# CLI/table-friendly branch names. FITS METHOD strings above keep hyphens.
TWIRL_FS_V2_BRANCH = "twirl_fs_v2"
TWIRL_FS_V2_ADP03Q_BRANCH = "twirl_fs_v2_adp03q"
TWIRL_FS_V2_ADP015Q_BRANCH = "twirl_fs_v2_adp015q"

ADP03Q_COLUMN_TAG = "ADP"
ADP015Q_COLUMN_TAG = "ADP015"


def canonical_config() -> FluxDetrendConfig:
    """Return the canonical TWIRL-FS v2 detrend configuration."""

    return FluxDetrendConfig(
        output_mode="subtractive",
        scale_strategy="auto",
    )


def adaptive_quantile_config(*, bkspace_d: float, gap_split_d: float) -> FluxDetrendConfig:
    """Return a TWIRL-FS adaptive quantile-knot configuration."""

    base = canonical_config()
    return FluxDetrendConfig(
        bkspace_d=float(bkspace_d),
        k=base.k,
        sigma_clip=base.sigma_clip,
        max_iter=base.max_iter,
        edge_pad_d=base.edge_pad_d,
        output_mode=base.output_mode,
        scale_strategy=base.scale_strategy,
        min_scale_abs=base.min_scale_abs,
        min_scale_snr=base.min_scale_snr,
        gap_split_d=float(gap_split_d),
        knot_strategy="quantile",
    )


def adp03q_config() -> FluxDetrendConfig:
    """Return the existing S56 compare-column ADP configuration."""

    return adaptive_quantile_config(bkspace_d=0.3, gap_split_d=0.2)


def adp015q_config() -> FluxDetrendConfig:
    """Return the current S56 BLS-recovery candidate configuration."""

    return adaptive_quantile_config(bkspace_d=0.15, gap_split_d=0.2)


def compare_column_names(tag: str) -> dict[str, str]:
    """Return FITS compare-column names for a primary/small/large branch tag."""

    normalized = str(tag).strip().upper()
    if not normalized:
        raise ValueError("compare column tag must be non-empty")
    return {
        "primary": f"DET_FLUX_{normalized}",
        "primary_err": f"DET_FLUX_{normalized}_ERR",
        "small": f"DET_FLUX_{normalized}_SML",
        "large": f"DET_FLUX_{normalized}_LAG",
    }
