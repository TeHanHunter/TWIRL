"""Scoped Tier-1 QA for A2v1 light curves.

Tier-0 establishes product integrity and the WD 1856 benchmark.  This module
adds the population, cadence, aperture, fixed-injection, and independent-
extraction evidence needed before the active search channels are trusted for
enrichment.  It deliberately distinguishes that bounded scope from promotion
of the complete six-aperture A2v1 product.

The gate is fail-closed: missing or stale evidence is a failure, not a warning.
An ``active_search_pair`` pass may set ``enrichment_ready`` but can never set
``science_ready``.  Only a future ``full_a2v1_product`` configuration with
promotion explicitly enabled can do that.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd

from twirl.io.hlsp import A2V1_APERTURES
from twirl.lightcurves.a2v1_cadence_reference import (
    CADENCE_REFERENCE_BUILDER_VERSION,
    authority_exclusions_sha256,
)
from twirl.lightcurves.a2v1_qa import (
    A2V1_PHOTOMETRIC_QA_VERSION,
    WD1856_PERIOD_D,
    WD1856_T0_BJD,
    WD1856_TIC,
    compact_target_catalog,
    file_sha256,
    stratified_target_sample,
)
from twirl.lightcurves.external_quality import (
    AUTHORITY_EXCLUSION_EXTERNAL_BIT,
    EFFECTIVE_QUALITY_POLICY,
    EXTERNAL_QUALITY_POLICY_CONTRACT,
    load_external_quality_reference,
)
from twirl.vetting.adp_only import ADP_ONLY_APERTURES


TIER1_QA_CONTRACT_VERSION = "a2v1_tier1_science_qa_v1"
TIER1_SCOPES = ("active_search_pair", "full_a2v1_product")
STATUS_ORDER = {"pass": 0, "review": 1, "fail": 2}
TIER0_REQUIRED_GATES = (
    "schema_and_completeness",
    "bls_target_aperture_coverage",
    "sampled_adp_photometry",
    "aperture_consistency",
    "prior_extraction_tree_comparison",
    "benchmark",
)
TIER1_TARGET_REASON_CODES = (
    "aperture_anticorrelation",
    "aperture_correlation_missing",
    "aperture_correlation_low",
    "aperture_mad_ratio",
    "aperture_missing",
    "cadence_loss",
    "duplicate_cadence",
    "finite_flux",
    "orbit_cadence_loss",
    "orbit_mismatch",
    "orbit_missing",
    "scatter_absolute_high",
    "unexpected_cadence",
    "usable_cadence",
)
INJECTION_QUALITY_COUNT_NAMES = (
    "n_cad_total",
    "n_cad_internal_bad",
    "n_cad_external_bad",
    "n_cad_authority_excluded",
    "n_cad_external_only_bad",
    "n_cad_effective_bad",
)


@dataclass(frozen=True)
class Tier1QAConfig:
    """Strict, predeclared thresholds for one Tier-1 evaluation."""

    name: str = "a2v1_tier1_active_search_pair_v1"
    contract_version: str = TIER1_QA_CONTRACT_VERSION
    sector: int = 56
    scope: str = "active_search_pair"
    promotion_enabled: bool = False
    expected_tier0_contract: str = A2V1_PHOTOMETRIC_QA_VERSION
    # Tier-1 may only consume the exact reviewed Tier-0 report and BLS table.
    # The all-zero defaults deliberately make an unreviewed production run
    # impossible.
    expected_tier0_summary_sha256: str = "0" * 64
    expected_tier0_bls_peaks_sha256: str = "0" * 64
    apertures: tuple[str, ...] = ADP_ONLY_APERTURES
    sample_size: int = 0
    sample_seed: int = 560101
    magnitude_edges: tuple[float, ...] = (0.0, 18.0, 19.0, 19.5, 20.0, 20.5, 30.0)
    min_population_targets: int = 512
    min_targets_per_magnitude_bin: int = 30
    min_supported_magnitude_bins: int = 3
    expected_compact_sha256: str = (
        "703193afd1dad78db71411859c69d2724cb40aad1419915ccb5c62495ce10d40"
    )
    # These two all-zero locks are deliberately impossible placeholders until
    # the authoritative S56 QLP-quaternion/SPOC table is built and reviewed.
    expected_cadence_reference_sha256: str = "0" * 64
    expected_cadence_reference_manifest_sha256: str = "0" * 64

    scatter_pass_finite_fraction: float = 0.98
    scatter_review_finite_fraction: float = 0.95
    scatter_pass_slope_min: float = -0.10
    scatter_pass_slope_max: float = 0.60
    scatter_review_slope_min: float = -0.20
    scatter_review_slope_max: float = 0.80
    scatter_high_residual_dex: float = 0.50
    scatter_pass_outlier_fraction: float = 0.05
    scatter_review_outlier_fraction: float = 0.10
    scatter_envelope_reference_tmag: float = 18.0
    scatter_envelope_dex_per_mag: float = 0.25
    scatter_pass_mad_ppm_at_reference: float = 400_000.0
    scatter_review_mad_ppm_at_reference: float = 700_000.0
    scatter_pass_rms5_ppm_at_reference: float = 600_000.0
    scatter_review_rms5_ppm_at_reference: float = 1_000_000.0

    cadence_reference_contract_template: str = "s{sector}_a2v1_cadence_reference_v1"
    cadence_pass_median_loss: float = 0.005
    cadence_review_median_loss: float = 0.01
    cadence_pass_p95_loss: float = 0.02
    cadence_review_p95_loss: float = 0.05
    finite_pass_median: float = 0.99
    finite_review_median: float = 0.95
    finite_pass_p05: float = 0.95
    finite_review_p05: float = 0.90
    usable_pass_median: float = 0.80
    usable_review_median: float = 0.50
    usable_pass_p05: float = 0.50
    usable_review_p05: float = 0.30

    aperture_pass_valid_fraction: float = 0.98
    aperture_review_valid_fraction: float = 0.95
    aperture_moderate_ratio_low: float = 1.0 / 3.0
    aperture_moderate_ratio_high: float = 3.0
    aperture_extreme_ratio_low: float = 0.10
    aperture_extreme_ratio_high: float = 10.0
    aperture_pass_moderate_outlier_fraction: float = 0.05
    aperture_review_moderate_outlier_fraction: float = 0.10
    aperture_pass_extreme_outlier_fraction: float = 0.005
    aperture_review_extreme_outlier_fraction: float = 0.01
    aperture_pass_correlation_valid_fraction: float = 0.98
    aperture_review_correlation_valid_fraction: float = 0.95
    aperture_pass_correlation_median: float = 0.10
    aperture_review_correlation_median: float = 0.0
    aperture_negative_correlation_threshold: float = -0.20
    aperture_pass_negative_correlation_fraction: float = 0.05
    aperture_review_negative_correlation_fraction: float = 0.10

    injection_contract_template: str = "s{sector}_a2v1_fresh_injection_pair_v1"
    injection_selection_mode: str = "full_shard"
    injection_required_shard_indices: tuple[int, ...] = (0, 13, 26, 39)
    # Hashes are ordered exactly like ``injection_required_shard_indices``.
    # Production values remain impossible placeholders until the four source
    # HDF5 shards are rehashed at the accepted evidence location.
    expected_injection_shard_sha256: tuple[str, ...] = ("0" * 64,) * 4
    injection_expected_ids: int = 2_000
    expected_injection_selection_sha256: str = (
        "8ee40ea7e450d5c19027f2e6c0c02279a78102caa246aeac80e032278ede6c5a"
    )
    expected_injection_metadata_sha256: str = (
        "549968899c236387d65700aa1e6d0d84a45ccff63cdf97d0b647f4d19cef1910"
    )
    expected_injection_source_adp_sha256: str = (
        "03f041bc66f09158d7740eaf9116b081d2184b7844720778d74f7ec0cf3b4505"
    )
    expected_injection_parity_report_sha256: str = (
        "f5e98f0e900beda0275d2a1d34ce592e799f5b0690db1490482c34c0f92360f9"
    )
    min_faint_injection_ids: int = 32
    injection_faint_tmag_min: float = 19.0
    min_injection_detectors: int = 8
    min_injection_period_ratio: float = 30.0
    min_injection_duration_ratio: float = 2.0
    min_injection_depth_ratio: float = 5.0
    injection_pass_retention_low: float = 0.85
    injection_pass_retention_high: float = 1.15
    injection_review_retention_low: float = 0.70
    injection_review_retention_high: float = 1.30
    injection_pass_p10: float = 0.70
    injection_review_p10: float = 0.50
    injection_inband_low: float = 0.50
    injection_inband_high: float = 1.50
    injection_pass_inband_fraction: float = 0.95
    injection_review_inband_fraction: float = 0.90

    independent_contract_template: str = "s{sector}_a2v1_independent_extraction_v2"
    # Deliberate impossible locks until the external WD 1856 extraction and
    # its metrics/manifest pair are built and independently reviewed.
    expected_independent_metrics_sha256: str = "0" * 64
    expected_independent_manifest_sha256: str = "0" * 64
    expected_independent_reference_product_sha256: str = "0" * 64
    min_independent_targets: int = 1
    min_independent_detectors: int = 1
    independent_comparison_mode: str = "signal_timing_only"
    independent_max_time_delta_seconds: float = 60.0
    independent_pass_common_fraction: float = 0.90
    independent_review_common_fraction: float = 0.80
    wd1856_pass_period_relative_error: float = 1.0e-4
    wd1856_review_period_relative_error: float = 5.0e-4
    wd1856_pass_epoch_residual_min: float = 5.0
    wd1856_review_epoch_residual_min: float = 15.0

    def validate(self) -> None:
        """Reject ambiguous or internally inconsistent gate configurations."""

        if self.contract_version != TIER1_QA_CONTRACT_VERSION:
            raise ValueError("unsupported Tier-1 contract_version")
        if self.sector != 56:
            raise ValueError("the locked Tier-1 v1 configuration is S56-specific")
        if self.scope not in TIER1_SCOPES:
            raise ValueError(f"unsupported Tier-1 scope: {self.scope}")
        if not self.apertures or len(set(self.apertures)) != len(self.apertures):
            raise ValueError("apertures must be nonempty and unique")
        unknown = sorted(set(self.apertures) - set(A2V1_APERTURES))
        if unknown:
            raise ValueError(f"non-A2v1 apertures requested: {unknown}")
        if self.scope == "active_search_pair" and tuple(self.apertures) != tuple(
            ADP_ONLY_APERTURES
        ):
            raise ValueError("active_search_pair must use the two locked ADP apertures")
        if self.promotion_enabled:
            raise ValueError(
                "promotion is disabled in Tier-1 v1 until the six-channel evaluator is implemented"
            )
        if self.sample_size != 0:
            raise ValueError("Tier-1 v1 requires a full-population scan (sample_size=0)")
        if min(
            self.min_population_targets,
            self.min_targets_per_magnitude_bin,
            self.min_supported_magnitude_bins,
        ) < 1:
            raise ValueError("population coverage counts must be positive")
        for value, label in (
            (self.expected_compact_sha256, "expected_compact_sha256"),
            (
                self.expected_tier0_summary_sha256,
                "expected_tier0_summary_sha256",
            ),
            (
                self.expected_tier0_bls_peaks_sha256,
                "expected_tier0_bls_peaks_sha256",
            ),
            (
                self.expected_cadence_reference_sha256,
                "expected_cadence_reference_sha256",
            ),
            (
                self.expected_cadence_reference_manifest_sha256,
                "expected_cadence_reference_manifest_sha256",
            ),
            (
                self.expected_injection_selection_sha256,
                "expected_injection_selection_sha256",
            ),
            (
                self.expected_injection_metadata_sha256,
                "expected_injection_metadata_sha256",
            ),
            (
                self.expected_injection_source_adp_sha256,
                "expected_injection_source_adp_sha256",
            ),
            (
                self.expected_injection_parity_report_sha256,
                "expected_injection_parity_report_sha256",
            ),
            (
                self.expected_independent_metrics_sha256,
                "expected_independent_metrics_sha256",
            ),
            (
                self.expected_independent_manifest_sha256,
                "expected_independent_manifest_sha256",
            ),
            (
                self.expected_independent_reference_product_sha256,
                "expected_independent_reference_product_sha256",
            ),
        ):
            if (
                value != value.lower()
                or len(value) != 64
                or any(char not in "0123456789abcdef" for char in value)
            ):
                raise ValueError(f"{label} must be a lowercase SHA-256 digest")
        edges = np.asarray(self.magnitude_edges, dtype=float)
        if len(edges) < 3 or not np.all(np.isfinite(edges)) or not np.all(np.diff(edges) > 0):
            raise ValueError("magnitude_edges must be finite and strictly increasing")
        if not np.isfinite(self.scatter_envelope_reference_tmag):
            raise ValueError("scatter envelope reference magnitude must be finite")
        if not 0 <= self.scatter_envelope_dex_per_mag <= 1:
            raise ValueError("invalid scatter envelope slope")
        if not (
            0 < self.scatter_pass_mad_ppm_at_reference
            <= self.scatter_review_mad_ppm_at_reference
        ):
            raise ValueError("invalid absolute MAD envelope")
        if not (
            0 < self.scatter_pass_rms5_ppm_at_reference
            <= self.scatter_review_rms5_ppm_at_reference
        ):
            raise ValueError("invalid absolute RMS envelope")
        if self.injection_selection_mode != "full_shard":
            raise ValueError("Tier-1 v1 requires complete, predeclared injection shards")
        if (
            not self.injection_required_shard_indices
            or len(set(self.injection_required_shard_indices))
            != len(self.injection_required_shard_indices)
            or min(self.injection_required_shard_indices) < 0
        ):
            raise ValueError("invalid required injection shard indices")
        if len(self.expected_injection_shard_sha256) != len(
            self.injection_required_shard_indices
        ):
            raise ValueError("one injection shard hash is required per locked shard index")
        for value in self.expected_injection_shard_sha256:
            if (
                value != value.lower()
                or len(value) != 64
                or any(char not in "0123456789abcdef" for char in value)
            ):
                raise ValueError(
                    "expected_injection_shard_sha256 values must be lowercase SHA-256 digests"
                )
        if min(
            self.injection_expected_ids,
            self.min_faint_injection_ids,
            self.min_injection_detectors,
        ) < 1:
            raise ValueError("injection coverage counts must be positive")
        if self.min_faint_injection_ids > self.injection_expected_ids:
            raise ValueError("faint injection count exceeds the frozen canary size")
        if min(
            self.min_injection_period_ratio,
            self.min_injection_duration_ratio,
            self.min_injection_depth_ratio,
        ) <= 1:
            raise ValueError("injection support ratios must exceed one")
        if min(self.min_independent_targets, self.min_independent_detectors) < 1:
            raise ValueError("independent-extraction coverage counts must be positive")
        if self.independent_comparison_mode != "signal_timing_only":
            raise ValueError(
                "Tier-1 v1 requires independent_comparison_mode="
                "'signal_timing_only'"
            )
        if not 0 < self.independent_max_time_delta_seconds <= 60.0:
            raise ValueError(
                "independent timestamp tolerance must be in (0, 60] seconds"
            )

        for passed, review, label in (
            (self.scatter_pass_finite_fraction, self.scatter_review_finite_fraction, "scatter finite fraction"),
            (self.finite_pass_median, self.finite_review_median, "finite median"),
            (self.finite_pass_p05, self.finite_review_p05, "finite p05"),
            (self.usable_pass_median, self.usable_review_median, "usable median"),
            (self.usable_pass_p05, self.usable_review_p05, "usable p05"),
            (self.aperture_pass_valid_fraction, self.aperture_review_valid_fraction, "aperture valid fraction"),
            (
                self.aperture_pass_correlation_valid_fraction,
                self.aperture_review_correlation_valid_fraction,
                "aperture correlation valid fraction",
            ),
            (self.injection_pass_p10, self.injection_review_p10, "injection p10"),
            (
                self.injection_pass_inband_fraction,
                self.injection_review_inband_fraction,
                "injection in-band fraction",
            ),
            (
                self.independent_pass_common_fraction,
                self.independent_review_common_fraction,
                "independent common fraction",
            ),
        ):
            _validate_pass_review_min(passed, review, label)
        for passed, review, label in (
            (self.cadence_pass_median_loss, self.cadence_review_median_loss, "median cadence loss"),
            (self.cadence_pass_p95_loss, self.cadence_review_p95_loss, "p95 cadence loss"),
            (self.scatter_pass_outlier_fraction, self.scatter_review_outlier_fraction, "scatter outlier fraction"),
            (
                self.aperture_pass_moderate_outlier_fraction,
                self.aperture_review_moderate_outlier_fraction,
                "moderate aperture outlier fraction",
            ),
            (
                self.aperture_pass_extreme_outlier_fraction,
                self.aperture_review_extreme_outlier_fraction,
                "extreme aperture outlier fraction",
            ),
            (
                self.aperture_pass_negative_correlation_fraction,
                self.aperture_review_negative_correlation_fraction,
                "negative aperture correlation fraction",
            ),
            (
                self.wd1856_pass_period_relative_error,
                self.wd1856_review_period_relative_error,
                "WD 1856 period",
            ),
            (
                self.wd1856_pass_epoch_residual_min,
                self.wd1856_review_epoch_residual_min,
                "WD 1856 timing",
            ),
        ):
            _validate_pass_review_max(passed, review, label)
        if not (
            -1 <= self.aperture_negative_correlation_threshold < 0
            and -1
            <= self.aperture_review_correlation_median
            <= self.aperture_pass_correlation_median
            <= 1
        ):
            raise ValueError("invalid aperture correlation thresholds")
        if not (
            self.scatter_review_slope_min
            <= self.scatter_pass_slope_min
            <= self.scatter_pass_slope_max
            <= self.scatter_review_slope_max
        ):
            raise ValueError("scatter slope pass interval must lie inside review interval")
        for review_low, pass_low, pass_high, review_high, label in (
            (
                self.aperture_extreme_ratio_low,
                self.aperture_moderate_ratio_low,
                self.aperture_moderate_ratio_high,
                self.aperture_extreme_ratio_high,
                "aperture ratio",
            ),
            (
                self.injection_review_retention_low,
                self.injection_pass_retention_low,
                self.injection_pass_retention_high,
                self.injection_review_retention_high,
                "injection retention",
            ),
        ):
            if not (0 < review_low <= pass_low <= pass_high <= review_high):
                raise ValueError(f"invalid nested interval for {label}")
        if not (0 < self.injection_inband_low < self.injection_inband_high):
            raise ValueError("invalid injection in-band interval")


def _validate_pass_review_min(passed: float, review: float, label: str) -> None:
    if not (0 <= review <= passed <= 1):
        raise ValueError(f"invalid pass/review minimum thresholds for {label}")


def _validate_pass_review_max(passed: float, review: float, label: str) -> None:
    if not (0 <= passed <= review):
        raise ValueError(f"invalid pass/review maximum thresholds for {label}")


def load_tier1_config(path: Path) -> Tier1QAConfig:
    """Load a strict YAML or JSON Tier-1 configuration."""

    path = Path(path)
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text())
    else:
        import yaml

        payload = yaml.safe_load(path.read_text())
    if not isinstance(payload, Mapping):
        raise ValueError(f"configuration must contain a mapping: {path}")
    allowed = set(Tier1QAConfig.__dataclass_fields__)
    unknown = sorted(set(payload) - allowed)
    if unknown:
        raise KeyError(f"unknown Tier-1 configuration fields: {unknown}")
    values = dict(payload)
    for key in (
        "apertures",
        "magnitude_edges",
        "injection_required_shard_indices",
        "expected_injection_shard_sha256",
    ):
        if key in values:
            values[key] = tuple(values[key])
    config = Tier1QAConfig(**values)
    config.validate()
    return config


def _finite(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(list(values), dtype=float)
    return array[np.isfinite(array)]


def _quantile(values: Iterable[float], q: float) -> float:
    array = _finite(values)
    return float(np.quantile(array, q)) if array.size else np.nan


def _robust_mad(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if not array.size:
        return np.nan
    center = float(np.median(array))
    return float(1.4826 * np.median(np.abs(array - center)))


def _clipped_rms(values: np.ndarray, sigma: float = 5.0, iterations: int = 2) -> float:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    for _ in range(iterations):
        if not array.size:
            return np.nan
        center = float(np.median(array))
        scale = _robust_mad(array)
        if not np.isfinite(scale) or scale <= 0:
            break
        array = array[np.abs(array - center) <= sigma * scale]
    if not array.size:
        return np.nan
    center = float(np.median(array))
    return float(np.sqrt(np.mean((array - center) ** 2)))


def _normalized(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    if not finite.size:
        return np.full(array.shape, np.nan)
    median = float(np.median(finite))
    if not np.isfinite(median) or abs(median) < 1.0e-12:
        return np.full(array.shape, np.nan)
    return array / median


def _correlation(left: np.ndarray, right: np.ndarray) -> float:
    good = np.isfinite(left) & np.isfinite(right)
    if np.count_nonzero(good) < 20:
        return np.nan
    x = np.asarray(left[good], dtype=float)
    y = np.asarray(right[good], dtype=float)
    if np.std(x) <= 0 or np.std(y) <= 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def _status_max(value: float, pass_limit: float, review_limit: float) -> str:
    if not np.isfinite(value):
        return "fail"
    if value <= pass_limit:
        return "pass"
    if value <= review_limit:
        return "review"
    return "fail"


def _status_min(value: float, pass_limit: float, review_limit: float) -> str:
    if not np.isfinite(value):
        return "fail"
    if value >= pass_limit:
        return "pass"
    if value >= review_limit:
        return "review"
    return "fail"


def _status_interval(
    value: float,
    pass_low: float,
    pass_high: float,
    review_low: float,
    review_high: float,
) -> str:
    if not np.isfinite(value):
        return "fail"
    if pass_low <= value <= pass_high:
        return "pass"
    if review_low <= value <= review_high:
        return "review"
    return "fail"


def _status_correlation(value: float, *, passed: float, failed_below: float) -> str:
    """Classify one aperture correlation without rejecting benign low values.

    The population gate still constrains the median and the strongly
    anticorrelated fraction.  At target level, values between the declared
    strong-anticorrelation boundary and the pass threshold are review rather
    than fail so quiet/noisy but non-pathological light curves remain visible.
    """

    if not np.isfinite(value):
        return "fail"
    if value >= passed:
        return "pass"
    if value >= failed_below:
        return "review"
    return "fail"


def _optional_gaia_dr3_source_id(group: h5py.Group) -> int | None:
    """Return an exact Gaia identifier when the compact group records one."""

    candidates = (
        "gaia_dr3_source_id",
        "gaia_source_id",
        "gaiaid",
        "gaia_id",
        "source_id",
    )
    for name in candidates:
        if name not in group.attrs:
            continue
        value = group.attrs[name]
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        # Gaia DR3 identifiers exceed the exact-integer range of float64.
        # Never turn a rounded floating-point attribute into an authoritative
        # identifier.
        if isinstance(value, (float, np.floating)):
            continue
        try:
            source_id = int(str(value).strip())
        except (TypeError, ValueError):
            continue
        if source_id > 0:
            return source_id
    return None


def _worst_status(statuses: Iterable[str]) -> str:
    values = list(statuses)
    return max(values, key=STATUS_ORDER.__getitem__) if values else "fail"


def _safe_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _safe_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_safe_json(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _safe_json(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_strict_json(path: Path, payload: Mapping[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        temporary.write_text(
            json.dumps(_safe_json(payload), indent=2, sort_keys=True, allow_nan=False)
            + "\n"
        )
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _dataframe_content_sha256(frame: pd.DataFrame) -> str:
    payload = {
        "columns": [str(column) for column in frame.columns],
        "records": frame.to_dict(orient="records"),
    }
    encoded = json.dumps(
        _safe_json(payload), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def injection_metadata_sha256(frame: pd.DataFrame) -> str:
    """Hash the immutable per-injection schedule/host metadata."""

    columns = (
        "injection_id",
        "tic",
        "tmag",
        "camera",
        "ccd",
        "period_d",
        "duration_min",
        "model_depth",
        "shard_index",
    )
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise KeyError(f"injection metadata hash is missing columns: {missing}")
    work = frame.loc[:, columns].drop_duplicates().copy()
    if work["injection_id"].astype(str).duplicated().any():
        raise ValueError("injection metadata differs across aperture rows")
    work["injection_id"] = work["injection_id"].astype(str)
    for column in ("tic", "camera", "ccd", "shard_index"):
        work[column] = pd.to_numeric(work[column], errors="raise").astype(np.int64)
    for column in ("tmag", "period_d", "duration_min", "model_depth"):
        work[column] = pd.to_numeric(work[column], errors="raise").astype(float)
    work = work.sort_values("injection_id", kind="stable").reset_index(drop=True)
    encoded = json.dumps(
        work.to_dict(orient="records"),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def audit_compact_population(
    compact_lc: Path,
    *,
    apertures: Sequence[str],
    cadence_reference: pd.DataFrame,
    authority_exclusions: Mapping[
        tuple[int, int, int], Sequence[int]
    ] | None = None,
    sample_size: int = 0,
    seed: int = 560101,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Measure target, aperture, and active-pair metrics from one compact export."""

    catalog = compact_target_catalog(compact_lc)
    reference_required = {
        "sector",
        "orbitid",
        "camera",
        "ccd",
        "cadenceno",
        "spoc_quality",
        "qlp_quality",
        "external_quality",
    }
    missing_reference = sorted(reference_required - set(cadence_reference.columns))
    if missing_reference:
        raise KeyError(f"cadence reference is missing columns: {missing_reference}")
    reference = cadence_reference.copy()
    for column in reference_required:
        reference[column] = pd.to_numeric(reference[column], errors="raise").astype(
            np.int64
        )
    observed_reference_sectors = sorted(set(reference["sector"].astype(int)))
    observed_catalog_sectors = sorted(set(catalog["sector"].astype(int)))
    if observed_reference_sectors != observed_catalog_sectors:
        raise ValueError(
            "cadence reference sector mismatch: "
            f"reference={observed_reference_sectors}, compact={observed_catalog_sectors}"
        )
    normalized_exclusions: dict[tuple[int, int, int], frozenset[int]] = {}
    for raw_key, raw_cadences in (authority_exclusions or {}).items():
        if len(raw_key) != 3:
            raise ValueError(
                "authority-exclusion keys must be (sector, camera, ccd)"
            )
        key = tuple(int(value) for value in raw_key)
        if key[0] not in observed_catalog_sectors:
            raise ValueError(
                f"authority exclusion has unexpected sector/detector key: {key}"
            )
        cadences = frozenset(int(value) for value in raw_cadences)
        if any(value < 0 for value in cadences):
            raise ValueError("authority-exclusion cadences must be nonnegative")
        normalized_exclusions[key] = cadences
    spoc_quality = reference["spoc_quality"].to_numpy(dtype=np.int64)
    qlp_quality = reference["qlp_quality"].to_numpy(dtype=np.int64)
    external_quality = reference["external_quality"].to_numpy(dtype=np.int64)
    if (
        np.any(spoc_quality < 0)
        or not np.isin(qlp_quality, (0, 1)).all()
        or not np.array_equal(
            external_quality, spoc_quality | (qlp_quality << 30)
        )
    ):
        raise ValueError("cadence reference has invalid external quality derivation")
    reference_maps: dict[tuple[int, int], dict[int, tuple[int, int]]] = {}
    reference_orbit_cadences: dict[
        tuple[int, int], dict[int, frozenset[int]]
    ] = {}
    for (camera, ccd), rows in reference.groupby(["camera", "ccd"], sort=False):
        if rows["cadenceno"].duplicated().any():
            raise ValueError(f"cadence reference has duplicate IDs for cam{camera}/ccd{ccd}")
        detector_key = (int(camera), int(ccd))
        reference_maps[detector_key] = {
            int(row.cadenceno): (int(row.external_quality), int(row.orbitid))
            for row in rows.itertuples(index=False)
        }
        reference_orbit_cadences[detector_key] = {
            int(orbitid): frozenset(group["cadenceno"].astype(np.int64))
            for orbitid, group in rows.groupby("orbitid", sort=True)
        }
    sample = (
        stratified_target_sample(catalog, sample_size=sample_size, seed=seed)
        if 0 < int(sample_size) < len(catalog)
        else catalog
    )
    target_rows: list[dict[str, Any]] = []
    aperture_rows: list[dict[str, Any]] = []
    pair_rows: list[dict[str, Any]] = []
    with h5py.File(Path(compact_lc), "r") as h5:
        for record in sample.itertuples(index=False):
            group_path = f"targets/{int(record.tic):016d}"
            group = h5[group_path]
            required = ("quality", "orbitid", "cadenceno", *apertures)
            missing = [name for name in required if name not in group]
            if missing:
                raise KeyError(f"TIC {record.tic} missing Tier-1 datasets: {missing}")
            quality = np.asarray(group["quality"], dtype=np.int64)
            cadence = np.asarray(group["cadenceno"], dtype=np.int64)
            if len(cadence) != len(quality):
                raise ValueError(f"TIC {record.tic}: cadence/quality length mismatch")
            internal_q0 = quality == 0
            orbitid = np.asarray(group["orbitid"], dtype=np.int64)
            if len(orbitid) != len(cadence):
                raise ValueError(f"TIC {record.tic}: cadence/orbit length mismatch")
            detector_key = (int(record.camera), int(record.ccd))
            if detector_key not in reference_maps:
                raise KeyError(
                    f"cadence reference lacks cam{record.camera}/ccd{record.ccd}"
                )
            reference_map = reference_maps[detector_key]
            reference_ids = set(reference_map)
            observed_ids = set(int(value) for value in cadence)
            duplicate_cadences = int(len(cadence) - len(observed_ids))
            missing_cadences = reference_ids - observed_ids
            declared_exclusions = normalized_exclusions.get(
                (
                    int(record.sector),
                    int(record.camera),
                    int(record.ccd),
                ),
                frozenset(),
            )
            if reference_ids & set(declared_exclusions):
                raise ValueError(
                    f"cadence reference includes declared authority exclusions for "
                    f"cam{record.camera}/ccd{record.ccd}"
                )
            unexpected_ids = observed_ids - reference_ids
            authority_excluded_cadences = unexpected_ids & set(
                declared_exclusions
            )
            unexpected_cadences = unexpected_ids - authority_excluded_cadences
            first_index: dict[int, int] = {}
            for index, cadence_id in enumerate(cadence):
                first_index.setdefault(int(cadence_id), int(index))
            external_quality = np.fromiter(
                (
                    (
                        int(reference_map[int(cadence_id)][0])
                        if int(cadence_id) in reference_map
                        else (
                            1 << AUTHORITY_EXCLUSION_EXTERNAL_BIT
                            if int(cadence_id) in authority_excluded_cadences
                            else -1
                        )
                    )
                    for cadence_id in cadence
                ),
                dtype=np.int64,
                count=len(cadence),
            )
            external_known = external_quality >= 0
            external_q0 = external_quality == 0
            q0 = internal_q0 & external_q0
            quality_mask_disagreements = int(
                np.count_nonzero(external_known & (internal_q0 != external_q0))
            )
            orbit_mismatches = sum(
                int(orbitid[first_index[cadence_id]])
                != int(reference_map[cadence_id][1])
                for cadence_id in reference_ids & observed_ids
            )
            expected_by_orbit = reference_orbit_cadences[detector_key]
            expected_orbits = sorted(expected_by_orbit)
            observed_by_orbit: dict[int, set[int]] = {}
            for cadence_id, observed_orbit in zip(cadence, orbitid, strict=True):
                observed_by_orbit.setdefault(int(observed_orbit), set()).add(
                    int(cadence_id)
                )
            orbit_details: dict[str, dict[str, int | float]] = {}
            missing_orbits: list[int] = []
            orbit_loss_fractions: list[float] = []
            n_observed_expected_orbits = 0
            for expected_orbit in expected_orbits:
                expected_orbit_ids = expected_by_orbit[expected_orbit]
                observed_orbit_ids = observed_by_orbit.get(expected_orbit, set())
                observed_expected_ids = expected_orbit_ids & observed_orbit_ids
                missing_orbit_ids = expected_orbit_ids - observed_orbit_ids
                loss_fraction = float(len(missing_orbit_ids) / len(expected_orbit_ids))
                orbit_loss_fractions.append(loss_fraction)
                if observed_expected_ids:
                    n_observed_expected_orbits += 1
                else:
                    missing_orbits.append(expected_orbit)
                orbit_details[str(expected_orbit)] = {
                    "n_expected": int(len(expected_orbit_ids)),
                    "n_observed": int(len(observed_expected_ids)),
                    "n_missing": int(len(missing_orbit_ids)),
                    "missing_fraction": loss_fraction,
                }
            gaia_dr3_source_id = _optional_gaia_dr3_source_id(group)
            target_rows.append(
                {
                    "tic": int(record.tic),
                    "gaia_dr3_source_id": (
                        str(gaia_dr3_source_id)
                        if gaia_dr3_source_id is not None
                        else pd.NA
                    ),
                    "sector": int(record.sector),
                    "camera": int(record.camera),
                    "ccd": int(record.ccd),
                    "detector": f"cam{int(record.camera)}_ccd{int(record.ccd)}",
                    "tmag": float(record.tmag),
                    "n_cadences": int(len(quality)),
                    "n_quality0": int(np.count_nonzero(q0)),
                    "quality0_fraction": float(np.mean(q0)) if len(q0) else np.nan,
                    "internal_quality0_fraction": (
                        float(np.mean(internal_q0)) if len(internal_q0) else np.nan
                    ),
                    "external_quality0_fraction": (
                        float(np.mean(external_q0[external_known]))
                        if np.any(external_known)
                        else np.nan
                    ),
                    "n_internal_quality_nonzero": int(
                        np.count_nonzero(~internal_q0)
                    ),
                    "n_external_quality_nonzero": int(
                        np.count_nonzero(external_known & ~external_q0)
                    ),
                    "n_effective_quality_nonzero": int(np.count_nonzero(~q0)),
                    "n_quality_mask_disagreements": quality_mask_disagreements,
                    "n_orbits": int(len(np.unique(orbitid))),
                    "n_expected_orbits": int(len(expected_orbits)),
                    "n_observed_expected_orbits": int(n_observed_expected_orbits),
                    "missing_orbit_ids": ";".join(str(value) for value in missing_orbits),
                    "max_orbit_missing_cadence_fraction": (
                        float(np.nanmax(orbit_loss_fractions))
                        if orbit_loss_fractions
                        else np.nan
                    ),
                    "orbit_cadence_metrics_json": json.dumps(
                        orbit_details,
                        sort_keys=True,
                        separators=(",", ":"),
                        allow_nan=False,
                    ),
                    "n_expected_cadences": int(len(reference_ids)),
                    "n_missing_cadences": int(len(missing_cadences)),
                    "missing_cadence_fraction": (
                        float(len(missing_cadences) / len(reference_ids))
                        if reference_ids
                        else np.nan
                    ),
                    "n_unexpected_cadences": int(len(unexpected_cadences)),
                    "n_authority_excluded_cadences": int(
                        len(authority_excluded_cadences)
                    ),
                    "n_duplicate_cadences": duplicate_cadences,
                    "n_orbit_mismatches": int(orbit_mismatches),
                }
            )
            normalized: dict[str, np.ndarray] = {}
            per_aperture: dict[str, dict[str, Any]] = {}
            for aperture in apertures:
                flux = np.asarray(group[aperture], dtype=float)
                norm = _normalized(flux)
                normalized[aperture] = norm
                good = q0 & np.isfinite(norm)
                row = {
                    "tic": int(record.tic),
                    "sector": int(record.sector),
                    "camera": int(record.camera),
                    "ccd": int(record.ccd),
                    "detector": f"cam{int(record.camera)}_ccd{int(record.ccd)}",
                    "tmag": float(record.tmag),
                    "aperture": str(aperture),
                    "n_quality0": int(np.count_nonzero(q0)),
                    "n_finite_quality0": int(np.count_nonzero(good)),
                    "finite_quality0_fraction": (
                        float(np.count_nonzero(good) / np.count_nonzero(q0))
                        if np.count_nonzero(q0)
                        else np.nan
                    ),
                    "mad_ppm": _robust_mad(norm[good]) * 1.0e6,
                    "rms5_ppm": _clipped_rms(norm[good]) * 1.0e6,
                }
                aperture_rows.append(row)
                per_aperture[aperture] = row
            if len(apertures) >= 2:
                left_name, right_name = apertures[:2]
                left = normalized[left_name]
                right = normalized[right_name]
                both = q0 & np.isfinite(left) & np.isfinite(right)
                left_mad = float(per_aperture[left_name]["mad_ppm"])
                right_mad = float(per_aperture[right_name]["mad_ppm"])
                pair_rows.append(
                    {
                        "tic": int(record.tic),
                        "sector": int(record.sector),
                        "camera": int(record.camera),
                        "ccd": int(record.ccd),
                        "detector": f"cam{int(record.camera)}_ccd{int(record.ccd)}",
                        "tmag": float(record.tmag),
                        "left_aperture": left_name,
                        "right_aperture": right_name,
                        "mad_ratio": (
                            right_mad / left_mad
                            if np.isfinite(left_mad) and left_mad > 0 and np.isfinite(right_mad)
                            else np.nan
                        ),
                        "correlation": _correlation(left[both], right[both]),
                        "difference_mad_ppm": _robust_mad(left[both] - right[both]) * 1.0e6,
                    }
                )
    targets = pd.DataFrame(target_rows)
    aperture_metrics = pd.DataFrame(aperture_rows)
    pair_metrics = pd.DataFrame(pair_rows)
    return targets, aperture_metrics, pair_metrics


def _robust_linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, np.ndarray]:
    good = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < 3:
        return np.nan, np.nan, np.full(y.shape, np.nan)
    keep = good.copy()
    slope = intercept = np.nan
    for _ in range(3):
        if np.count_nonzero(keep) < 3:
            break
        slope, intercept = np.polyfit(x[keep], y[keep], 1)
        residual = y - (slope * x + intercept)
        center = float(np.nanmedian(residual[keep]))
        scale = _robust_mad(residual[keep])
        if not np.isfinite(scale) or scale <= 0:
            break
        keep = good & (np.abs(residual - center) <= 5.0 * scale)
    residual = y - (slope * x + intercept)
    return float(slope), float(intercept), residual


def evaluate_population_scatter(
    aperture_metrics: pd.DataFrame, config: Tier1QAConfig
) -> tuple[dict[str, Any], pd.DataFrame]:
    required = {"tic", "tmag", "aperture", "mad_ppm", "rms5_ppm", "detector"}
    missing = sorted(required - set(aperture_metrics.columns))
    if missing:
        return {"status": "fail", "reasons": [f"missing columns: {missing}"]}, pd.DataFrame()
    summaries: list[dict[str, Any]] = []
    bin_rows: list[dict[str, Any]] = []
    statuses: list[str] = []
    reasons: list[str] = []
    for aperture in config.apertures:
        rows = aperture_metrics.loc[aperture_metrics["aperture"].astype(str).eq(aperture)].copy()
        n_targets = int(rows["tic"].nunique())
        finite = (
            np.isfinite(pd.to_numeric(rows["tmag"], errors="coerce"))
            & np.isfinite(pd.to_numeric(rows["mad_ppm"], errors="coerce"))
            & (pd.to_numeric(rows["mad_ppm"], errors="coerce") > 0)
            & np.isfinite(pd.to_numeric(rows["rms5_ppm"], errors="coerce"))
            & (pd.to_numeric(rows["rms5_ppm"], errors="coerce") > 0)
        )
        finite_fraction = float(finite.mean()) if len(rows) else 0.0
        work = rows.loc[finite].copy()
        x = work["tmag"].to_numpy(dtype=float)
        y = np.log10(work["mad_ppm"].to_numpy(dtype=float))
        slope, intercept, residual = _robust_linear_fit(x, y)
        outlier_fraction = (
            float(np.mean(residual > config.scatter_high_residual_dex))
            if residual.size
            else np.nan
        )
        work["magnitude_bin"] = pd.cut(
            work["tmag"], bins=list(config.magnitude_edges), include_lowest=True
        )
        supported_bins = 0
        mad_envelope_statuses: list[str] = []
        rms_envelope_statuses: list[str] = []
        for label, group in work.groupby("magnitude_bin", observed=True):
            count = int(len(group))
            supported = count >= config.min_targets_per_magnitude_bin
            tmag_median = _quantile(group["tmag"], 0.5)
            mad_median = _quantile(group["mad_ppm"], 0.5)
            rms_median = _quantile(group["rms5_ppm"], 0.5)
            scale = 10.0 ** (
                config.scatter_envelope_dex_per_mag
                * (tmag_median - config.scatter_envelope_reference_tmag)
            )
            mad_status = _status_max(
                mad_median,
                config.scatter_pass_mad_ppm_at_reference * scale,
                config.scatter_review_mad_ppm_at_reference * scale,
            )
            rms_status = _status_max(
                rms_median,
                config.scatter_pass_rms5_ppm_at_reference * scale,
                config.scatter_review_rms5_ppm_at_reference * scale,
            )
            if supported:
                supported_bins += 1
                mad_envelope_statuses.append(mad_status)
                rms_envelope_statuses.append(rms_status)
            bin_rows.append(
                {
                    "aperture": aperture,
                    "magnitude_bin": str(label),
                    "n": count,
                    "supported": supported,
                    "tmag_median": tmag_median,
                    "mad_ppm_median": mad_median,
                    "mad_ppm_p90": _quantile(group["mad_ppm"], 0.9),
                    "rms5_ppm_median": rms_median,
                    "mad_envelope_status": mad_status,
                    "rms_envelope_status": rms_status,
                }
            )
        finite_status = _status_min(
            finite_fraction,
            config.scatter_pass_finite_fraction,
            config.scatter_review_finite_fraction,
        )
        slope_status = _status_interval(
            slope,
            config.scatter_pass_slope_min,
            config.scatter_pass_slope_max,
            config.scatter_review_slope_min,
            config.scatter_review_slope_max,
        )
        outlier_status = _status_max(
            outlier_fraction,
            config.scatter_pass_outlier_fraction,
            config.scatter_review_outlier_fraction,
        )
        coverage_status = (
            "pass"
            if n_targets >= config.min_population_targets
            and supported_bins >= config.min_supported_magnitude_bins
            else "fail"
        )
        mad_envelope_status = _worst_status(mad_envelope_statuses)
        rms_envelope_status = _worst_status(rms_envelope_statuses)
        aperture_status = _worst_status(
            [
                finite_status,
                slope_status,
                outlier_status,
                coverage_status,
                mad_envelope_status,
                rms_envelope_status,
            ]
        )
        statuses.append(aperture_status)
        if aperture_status == "fail":
            reasons.append(f"{aperture} failed scatter/coverage limits")
        summaries.append(
            {
                "aperture": aperture,
                "status": aperture_status,
                "n_targets": n_targets,
                "finite_metric_fraction": finite_fraction,
                "slope_log10_mad_per_mag": slope,
                "intercept_log10_mad": intercept,
                "high_residual_fraction": outlier_fraction,
                "supported_magnitude_bins": supported_bins,
                "component_status": {
                    "coverage": coverage_status,
                    "finite_metrics": finite_status,
                    "slope": slope_status,
                    "outliers": outlier_status,
                    "absolute_mad_envelope": mad_envelope_status,
                    "absolute_rms_envelope": rms_envelope_status,
                },
            }
        )
    return {
        "status": _worst_status(statuses),
        "reasons": reasons,
        "apertures": summaries,
    }, pd.DataFrame(bin_rows)


def evaluate_cadence_quality(
    target_metrics: pd.DataFrame,
    aperture_metrics: pd.DataFrame,
    config: Tier1QAConfig,
) -> tuple[dict[str, Any], pd.DataFrame]:
    target_required = {
        "tic",
        "sector",
        "detector",
        "n_cadences",
        "n_expected_cadences",
        "missing_cadence_fraction",
        "n_expected_orbits",
        "n_observed_expected_orbits",
        "missing_orbit_ids",
        "max_orbit_missing_cadence_fraction",
        "n_unexpected_cadences",
        "n_authority_excluded_cadences",
        "n_duplicate_cadences",
        "n_orbit_mismatches",
        "quality0_fraction",
        "internal_quality0_fraction",
        "external_quality0_fraction",
        "n_quality_mask_disagreements",
    }
    aperture_required = {"tic", "aperture", "finite_quality0_fraction"}
    missing = sorted(
        (target_required - set(target_metrics.columns))
        | (aperture_required - set(aperture_metrics.columns))
    )
    if missing:
        return {"status": "fail", "reasons": [f"missing columns: {missing}"]}, pd.DataFrame()
    targets = target_metrics.copy()
    if targets["tic"].duplicated().any():
        return {
            "status": "fail",
            "reasons": ["target metrics contain duplicate TIC rows"],
        }, targets
    observed_sectors = sorted(
        set(pd.to_numeric(targets["sector"], errors="coerce").dropna().astype(int))
    )
    if observed_sectors != [config.sector]:
        return {
            "status": "fail",
            "reasons": [f"target sector mismatch: {observed_sectors}"],
        }, targets
    targets["relative_cadence_loss"] = pd.to_numeric(
        targets["missing_cadence_fraction"], errors="coerce"
    )
    median_loss = _quantile(targets["relative_cadence_loss"], 0.5)
    p95_loss = _quantile(targets["relative_cadence_loss"], 0.95)
    orbit_loss = pd.to_numeric(
        targets["max_orbit_missing_cadence_fraction"], errors="coerce"
    )
    p95_orbit_loss = _quantile(orbit_loss, 0.95)
    selected = aperture_metrics.loc[
        aperture_metrics["aperture"].astype(str).isin(config.apertures)
    ].copy()
    duplicate_apertures = selected.duplicated(["tic", "aperture"]).any()
    selected = selected.merge(
        targets.loc[:, ["tic", "quality0_fraction"]],
        on="tic",
        how="inner",
        validate="many_to_one",
    )
    finite_values = selected["finite_quality0_fraction"]
    finite_median = _quantile(finite_values, 0.5)
    finite_p05 = _quantile(finite_values, 0.05)
    selected["usable_cadence_fraction"] = (
        pd.to_numeric(selected["quality0_fraction"], errors="coerce")
        * pd.to_numeric(selected["finite_quality0_fraction"], errors="coerce")
    )
    usable_median = _quantile(selected["usable_cadence_fraction"], 0.5)
    usable_p05 = _quantile(selected["usable_cadence_fraction"], 0.05)
    anomaly_columns = (
        "n_unexpected_cadences",
        "n_duplicate_cadences",
        "n_orbit_mismatches",
    )
    anomaly_totals = {
        column: int(pd.to_numeric(targets[column], errors="coerce").fillna(-1).sum())
        for column in anomaly_columns
    }
    expected_positive = bool(
        (pd.to_numeric(targets["n_expected_cadences"], errors="coerce") > 0).all()
    )
    expected_aperture_rows = len(targets) * len(config.apertures)
    exact_aperture_rows = bool(
        not duplicate_apertures
        and len(selected) == expected_aperture_rows
        and selected["tic"].nunique() == len(targets)
    )
    authoritative_integrity = bool(
        expected_positive
        and exact_aperture_rows
        and all(value == 0 for value in anomaly_totals.values())
    )
    statuses = {
        "median_relative_cadence_loss": _status_max(
            median_loss,
            config.cadence_pass_median_loss,
            config.cadence_review_median_loss,
        ),
        "p95_relative_cadence_loss": _status_max(
            p95_loss, config.cadence_pass_p95_loss, config.cadence_review_p95_loss
        ),
        "p95_max_orbit_cadence_loss": _status_max(
            p95_orbit_loss,
            config.cadence_pass_p95_loss,
            config.cadence_review_p95_loss,
        ),
        "finite_quality0_median": _status_min(
            finite_median, config.finite_pass_median, config.finite_review_median
        ),
        "finite_quality0_p05": _status_min(
            finite_p05, config.finite_pass_p05, config.finite_review_p05
        ),
        "usable_cadence_median": _status_min(
            usable_median, config.usable_pass_median, config.usable_review_median
        ),
        "usable_cadence_p05": _status_min(
            usable_p05, config.usable_pass_p05, config.usable_review_p05
        ),
        "authoritative_reference_integrity": (
            "pass" if authoritative_integrity else "fail"
        ),
    }

    detector_summaries: list[dict[str, Any]] = []
    detector_statuses: list[str] = []
    for detector, detector_targets in targets.groupby("detector", sort=True):
        detector_tics = set(detector_targets["tic"].astype(np.int64))
        detector_apertures = selected.loc[selected["tic"].isin(detector_tics)]
        detector_loss = pd.to_numeric(
            detector_targets["relative_cadence_loss"], errors="coerce"
        )
        detector_orbit_loss = pd.to_numeric(
            detector_targets["max_orbit_missing_cadence_fraction"],
            errors="coerce",
        )
        detector_finite = pd.to_numeric(
            detector_apertures["finite_quality0_fraction"], errors="coerce"
        )
        detector_usable = pd.to_numeric(
            detector_apertures["usable_cadence_fraction"], errors="coerce"
        )
        detector_internal_good = pd.to_numeric(
            detector_targets["internal_quality0_fraction"], errors="coerce"
        )
        detector_external_good = pd.to_numeric(
            detector_targets["external_quality0_fraction"], errors="coerce"
        )
        detector_quality_disagreements = int(
            pd.to_numeric(
                detector_targets["n_quality_mask_disagreements"], errors="coerce"
            )
            .fillna(0)
            .sum()
        )
        detector_anomalies = {
            column: int(
                pd.to_numeric(detector_targets[column], errors="coerce")
                .fillna(-1)
                .sum()
            )
            for column in anomaly_columns
        }
        detector_exact_apertures = bool(
            not detector_apertures.duplicated(["tic", "aperture"]).any()
            and len(detector_apertures)
            == len(detector_targets) * len(config.apertures)
            and detector_apertures["tic"].nunique() == len(detector_targets)
        )
        detector_components = {
            "median_relative_cadence_loss": _status_max(
                _quantile(detector_loss, 0.5),
                config.cadence_pass_median_loss,
                config.cadence_review_median_loss,
            ),
            "p95_relative_cadence_loss": _status_max(
                _quantile(detector_loss, 0.95),
                config.cadence_pass_p95_loss,
                config.cadence_review_p95_loss,
            ),
            "p95_max_orbit_cadence_loss": _status_max(
                _quantile(detector_orbit_loss, 0.95),
                config.cadence_pass_p95_loss,
                config.cadence_review_p95_loss,
            ),
            "finite_quality0_median": _status_min(
                _quantile(detector_finite, 0.5),
                config.finite_pass_median,
                config.finite_review_median,
            ),
            "finite_quality0_p05": _status_min(
                _quantile(detector_finite, 0.05),
                config.finite_pass_p05,
                config.finite_review_p05,
            ),
            "usable_cadence_median": _status_min(
                _quantile(detector_usable, 0.5),
                config.usable_pass_median,
                config.usable_review_median,
            ),
            "usable_cadence_p05": _status_min(
                _quantile(detector_usable, 0.05),
                config.usable_pass_p05,
                config.usable_review_p05,
            ),
            "reference_integrity": (
                "pass"
                if detector_exact_apertures
                and all(value == 0 for value in detector_anomalies.values())
                else "fail"
            ),
        }
        detector_status = _worst_status(detector_components.values())
        detector_statuses.append(detector_status)
        detector_summaries.append(
            {
                "detector": str(detector),
                "status": detector_status,
                "n_targets": int(len(detector_targets)),
                "median_relative_cadence_loss": _quantile(detector_loss, 0.5),
                "p95_relative_cadence_loss": _quantile(detector_loss, 0.95),
                "p95_max_orbit_cadence_loss": _quantile(
                    detector_orbit_loss, 0.95
                ),
                "finite_quality0_fraction_p05": _quantile(detector_finite, 0.05),
                "usable_cadence_fraction_p05": _quantile(detector_usable, 0.05),
                "internal_quality0_fraction_median": _quantile(
                    detector_internal_good, 0.5
                ),
                "external_quality0_fraction_median": _quantile(
                    detector_external_good, 0.5
                ),
                "quality_mask_disagreements_total": detector_quality_disagreements,
                "authority_excluded_cadences_total": int(
                    pd.to_numeric(
                        detector_targets["n_authority_excluded_cadences"],
                        errors="coerce",
                    )
                    .fillna(0)
                    .sum()
                ),
                "reference_anomaly_totals": detector_anomalies,
                "component_status": detector_components,
            }
        )
    statuses["detector_stratified"] = _worst_status(detector_statuses)

    per_target_aperture = (
        selected.groupby("tic", as_index=False)
        .agg(
            n_tier1_apertures=("aperture", "nunique"),
            finite_quality0_fraction_min=("finite_quality0_fraction", "min"),
            usable_cadence_fraction_min=("usable_cadence_fraction", "min"),
        )
    )
    targets = targets.merge(
        per_target_aperture, on="tic", how="left", validate="one_to_one"
    )

    def _target_status(row: pd.Series) -> tuple[str, str]:
        reasons: list[str] = []
        components = [
            _status_max(
                float(row["relative_cadence_loss"]),
                config.cadence_pass_p95_loss,
                config.cadence_review_p95_loss,
            ),
            _status_max(
                float(row["max_orbit_missing_cadence_fraction"]),
                config.cadence_pass_p95_loss,
                config.cadence_review_p95_loss,
            ),
            _status_min(
                float(row["finite_quality0_fraction_min"]),
                config.finite_pass_p05,
                config.finite_review_p05,
            ),
            _status_min(
                float(row["usable_cadence_fraction_min"]),
                config.usable_pass_p05,
                config.usable_review_p05,
            ),
        ]
        if components[0] != "pass":
            reasons.append("cadence_loss")
        if components[1] != "pass":
            reasons.append("orbit_cadence_loss")
        if components[2] != "pass":
            reasons.append("finite_flux")
        if components[3] != "pass":
            reasons.append("usable_cadence")
        missing_orbits = str(row.get("missing_orbit_ids", "")).strip()
        if (
            missing_orbits
            or int(row["n_observed_expected_orbits"]) != int(row["n_expected_orbits"])
        ):
            components.append("fail")
            reasons.append("orbit_missing")
        anomaly_reasons = {
            "n_unexpected_cadences": "unexpected_cadence",
            "n_duplicate_cadences": "duplicate_cadence",
            "n_orbit_mismatches": "orbit_mismatch",
        }
        for column, reason in anomaly_reasons.items():
            if int(row[column]) != 0:
                components.append("fail")
                reasons.append(reason)
        n_tier1_apertures = pd.to_numeric(
            pd.Series([row.get("n_tier1_apertures", np.nan)]), errors="coerce"
        ).iloc[0]
        if (
            not np.isfinite(n_tier1_apertures)
            or int(n_tier1_apertures) != len(config.apertures)
        ):
            components.append("fail")
            reasons.append("aperture_missing")
        return _worst_status(components), ";".join(reasons)

    target_flags = targets.apply(_target_status, axis=1, result_type="expand")
    targets["tier1_cadence_qa_status"] = target_flags[0]
    targets["tier1_cadence_qa_reasons"] = target_flags[1]
    targets["tier1_cadence_qa_pass"] = targets["tier1_cadence_qa_status"].eq("pass")
    status = _worst_status(statuses.values())
    return {
        "status": status,
        "reasons": ["cadence or finite-data limit failed"] if status == "fail" else [],
        "cadence_reference": "authoritative detector/orbit cadence table",
        "cadence_reference_contract": config.cadence_reference_contract_template.format(
            sector=config.sector
        ),
        "median_relative_cadence_loss": median_loss,
        "p95_relative_cadence_loss": p95_loss,
        "p95_max_orbit_cadence_loss": p95_orbit_loss,
        "authoritative_reference_integrity": authoritative_integrity,
        "reference_anomaly_totals": anomaly_totals,
        "quality0_fraction_median": _quantile(targets["quality0_fraction"], 0.5),
        "quality0_fraction_p05": _quantile(targets["quality0_fraction"], 0.05),
        "internal_quality0_fraction_median": _quantile(
            targets["internal_quality0_fraction"], 0.5
        ),
        "external_quality0_fraction_median": _quantile(
            targets["external_quality0_fraction"], 0.5
        ),
        "quality_mask_disagreements_total": int(
            pd.to_numeric(
                targets["n_quality_mask_disagreements"], errors="coerce"
            )
            .fillna(0)
            .sum()
        ),
        "authority_excluded_cadences_total": int(
            pd.to_numeric(
                targets["n_authority_excluded_cadences"], errors="coerce"
            )
            .fillna(0)
            .sum()
        ),
        "finite_quality0_fraction_median": finite_median,
        "finite_quality0_fraction_p05": finite_p05,
        "usable_cadence_fraction_median": usable_median,
        "usable_cadence_fraction_p05": usable_p05,
        "n_target_qa_pass": int(targets["tier1_cadence_qa_pass"].sum()),
        "n_target_qa_review": int(targets["tier1_cadence_qa_status"].eq("review").sum()),
        "n_target_qa_fail": int(targets["tier1_cadence_qa_status"].eq("fail").sum()),
        "detectors": detector_summaries,
        "component_status": statuses,
    }, targets


def evaluate_aperture_outliers(
    pair_metrics: pd.DataFrame, config: Tier1QAConfig
) -> dict[str, Any]:
    required = {"tic", "mad_ratio", "correlation"}
    missing = sorted(required - set(pair_metrics.columns))
    if missing or pair_metrics.empty:
        return {"status": "fail", "reasons": [f"missing aperture-pair evidence: {missing}"]}
    ratio = pd.to_numeric(pair_metrics["mad_ratio"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(ratio) & (ratio > 0)
    valid_fraction = float(np.mean(valid))
    moderate = valid & (
        (ratio < config.aperture_moderate_ratio_low)
        | (ratio > config.aperture_moderate_ratio_high)
    )
    extreme = valid & (
        (ratio < config.aperture_extreme_ratio_low)
        | (ratio > config.aperture_extreme_ratio_high)
    )
    moderate_fraction = float(np.count_nonzero(moderate) / max(np.count_nonzero(valid), 1))
    extreme_fraction = float(np.count_nonzero(extreme) / max(np.count_nonzero(valid), 1))
    correlation = pd.to_numeric(pair_metrics["correlation"], errors="coerce").to_numpy(
        dtype=float
    )
    correlation_valid = np.isfinite(correlation)
    correlation_valid_fraction = float(np.mean(correlation_valid)) if len(correlation) else 0.0
    correlation_median = _quantile(correlation[correlation_valid], 0.5)
    negative_correlation_fraction = (
        float(
            np.mean(
                correlation[correlation_valid]
                < config.aperture_negative_correlation_threshold
            )
        )
        if np.count_nonzero(correlation_valid)
        else np.nan
    )
    statuses = {
        "valid_fraction": _status_min(
            valid_fraction,
            config.aperture_pass_valid_fraction,
            config.aperture_review_valid_fraction,
        ),
        "moderate_outlier_fraction": _status_max(
            moderate_fraction,
            config.aperture_pass_moderate_outlier_fraction,
            config.aperture_review_moderate_outlier_fraction,
        ),
        "extreme_outlier_fraction": _status_max(
            extreme_fraction,
            config.aperture_pass_extreme_outlier_fraction,
            config.aperture_review_extreme_outlier_fraction,
        ),
        "correlation_valid_fraction": _status_min(
            correlation_valid_fraction,
            config.aperture_pass_correlation_valid_fraction,
            config.aperture_review_correlation_valid_fraction,
        ),
        "correlation_median": _status_min(
            correlation_median,
            config.aperture_pass_correlation_median,
            config.aperture_review_correlation_median,
        ),
        "negative_correlation_fraction": _status_max(
            negative_correlation_fraction,
            config.aperture_pass_negative_correlation_fraction,
            config.aperture_review_negative_correlation_fraction,
        ),
    }
    status = _worst_status(statuses.values())
    return {
        "status": status,
        "reasons": ["aperture ratio or correlation limit failed"] if status == "fail" else [],
        "n_targets": int(len(pair_metrics)),
        "valid_ratio_fraction": valid_fraction,
        "mad_ratio_median": _quantile(ratio[valid], 0.5),
        "moderate_outlier_fraction": moderate_fraction,
        "extreme_outlier_fraction": extreme_fraction,
        "correlation_valid_fraction": correlation_valid_fraction,
        "correlation_median": correlation_median,
        "negative_correlation_fraction": negative_correlation_fraction,
        "component_status": statuses,
    }


def attach_target_qa_flags(
    targets: pd.DataFrame,
    aperture_metrics: pd.DataFrame,
    pair_metrics: pd.DataFrame,
    config: Tier1QAConfig,
) -> pd.DataFrame:
    """Attach fail-closed target flags suitable for candidate/teacher joins."""

    required = {"tic", "tmag", "tier1_cadence_qa_status", "tier1_cadence_qa_reasons"}
    missing = sorted(required - set(targets.columns))
    if missing:
        raise KeyError(f"target QA inputs are missing columns: {missing}")
    aperture_required = {"tic", "tmag", "aperture", "mad_ppm", "rms5_ppm"}
    missing = sorted(aperture_required - set(aperture_metrics.columns))
    if missing:
        raise KeyError(f"aperture QA inputs are missing columns: {missing}")
    pair_required = {"tic", "mad_ratio", "correlation"}
    missing = sorted(pair_required - set(pair_metrics.columns))
    if missing:
        raise KeyError(f"aperture-pair QA inputs are missing columns: {missing}")

    aperture = aperture_metrics.loc[
        aperture_metrics["aperture"].astype(str).isin(config.apertures)
    ].copy()
    tmag = pd.to_numeric(aperture["tmag"], errors="coerce")
    scale = 10.0 ** (
        config.scatter_envelope_dex_per_mag
        * (tmag - config.scatter_envelope_reference_tmag)
    )
    mad = pd.to_numeric(aperture["mad_ppm"], errors="coerce")
    rms = pd.to_numeric(aperture["rms5_ppm"], errors="coerce")
    aperture["target_scatter_status"] = [
        _worst_status(
            (
                _status_max(
                    mad_value,
                    config.scatter_pass_mad_ppm_at_reference * scale_value,
                    config.scatter_review_mad_ppm_at_reference * scale_value,
                ),
                _status_max(
                    rms_value,
                    config.scatter_pass_rms5_ppm_at_reference * scale_value,
                    config.scatter_review_rms5_ppm_at_reference * scale_value,
                ),
            )
        )
        for mad_value, rms_value, scale_value in zip(mad, rms, scale, strict=True)
    ]
    scatter_by_target = (
        aperture.groupby("tic", as_index=False)
        .agg(
            n_scatter_apertures=("aperture", "nunique"),
            tier1_scatter_qa_status=(
                "target_scatter_status",
                lambda values: _worst_status(values),
            ),
        )
    )

    pairs = pair_metrics.copy()
    ratio = pd.to_numeric(pairs["mad_ratio"], errors="coerce")
    correlation = pd.to_numeric(pairs["correlation"], errors="coerce")
    pairs["tier1_aperture_ratio_qa_status"] = [
        _status_interval(
            ratio_value,
            config.aperture_moderate_ratio_low,
            config.aperture_moderate_ratio_high,
            config.aperture_extreme_ratio_low,
            config.aperture_extreme_ratio_high,
        )
        for ratio_value in ratio
    ]
    pairs["tier1_aperture_correlation_qa_status"] = [
        _status_correlation(
            correlation_value,
            passed=config.aperture_pass_correlation_median,
            failed_below=config.aperture_negative_correlation_threshold,
        )
        for correlation_value in correlation
    ]
    pairs["tier1_aperture_pair_qa_status"] = [
        _worst_status((ratio_status, correlation_status))
        for ratio_status, correlation_status in zip(
            pairs["tier1_aperture_ratio_qa_status"],
            pairs["tier1_aperture_correlation_qa_status"],
            strict=True,
        )
    ]
    pair_by_target = (
        pairs.groupby("tic", as_index=False)
        .agg(
            n_aperture_pairs=("tic", "size"),
            aperture_correlation=("correlation", "first"),
            tier1_aperture_ratio_qa_status=(
                "tier1_aperture_ratio_qa_status",
                lambda values: _worst_status(values),
            ),
            tier1_aperture_correlation_qa_status=(
                "tier1_aperture_correlation_qa_status",
                lambda values: _worst_status(values),
            ),
            tier1_aperture_pair_qa_status=(
                "tier1_aperture_pair_qa_status",
                lambda values: _worst_status(values),
            ),
        )
    )

    out = targets.merge(scatter_by_target, on="tic", how="left", validate="one_to_one")
    out = out.merge(pair_by_target, on="tic", how="left", validate="one_to_one")

    def _combined(row: pd.Series) -> tuple[str, str]:
        cadence = str(row["tier1_cadence_qa_status"])
        scatter = str(row.get("tier1_scatter_qa_status", "fail"))
        aperture_pair = str(row.get("tier1_aperture_pair_qa_status", "fail"))
        n_scatter = pd.to_numeric(
            pd.Series([row.get("n_scatter_apertures", np.nan)]), errors="coerce"
        ).iloc[0]
        n_pairs = pd.to_numeric(
            pd.Series([row.get("n_aperture_pairs", np.nan)]), errors="coerce"
        ).iloc[0]
        if not np.isfinite(n_scatter) or int(n_scatter) != len(config.apertures):
            scatter = "fail"
        if not np.isfinite(n_pairs) or int(n_pairs) != 1:
            aperture_pair = "fail"
        status = _worst_status((cadence, scatter, aperture_pair))
        reasons = [
            value
            for value in str(row["tier1_cadence_qa_reasons"]).split(";")
            if value and value != "nan"
        ]
        if scatter != "pass":
            reasons.append("scatter_absolute_high")
        ratio_status = str(row.get("tier1_aperture_ratio_qa_status", "fail"))
        correlation_status = str(
            row.get("tier1_aperture_correlation_qa_status", "fail")
        )
        if ratio_status != "pass":
            reasons.append("aperture_mad_ratio")
        correlation_value = pd.to_numeric(
            pd.Series([row.get("aperture_correlation", np.nan)]), errors="coerce"
        ).iloc[0]
        if not np.isfinite(correlation_value):
            reasons.append("aperture_correlation_missing")
        elif correlation_status == "review":
            reasons.append("aperture_correlation_low")
        elif correlation_status == "fail":
            reasons.append("aperture_anticorrelation")
        if not np.isfinite(n_scatter) or int(n_scatter) != len(config.apertures):
            reasons.append("aperture_missing")
        if not np.isfinite(n_pairs) or int(n_pairs) != 1:
            reasons.append("aperture_missing")
        unknown_reasons = sorted(set(reasons) - set(TIER1_TARGET_REASON_CODES))
        if unknown_reasons:
            raise RuntimeError(f"unknown Tier-1 target reason codes: {unknown_reasons}")
        return status, ";".join(dict.fromkeys(reasons))

    combined = out.apply(_combined, axis=1, result_type="expand")
    out["tier1_target_qa_status"] = combined[0]
    out["tier1_target_qa_reasons"] = combined[1]
    out["tier1_target_qa_pass"] = out["tier1_target_qa_status"].eq("pass")
    return out


def summarize_fixed_injection_shards(
    shard_paths: Sequence[Path],
    *,
    sector: int,
    cadence_reference_path: Path,
    cadence_reference_manifest_path: Path,
    apertures: Sequence[str] = ADP_ONLY_APERTURES,
    selected_ids: set[str] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Fit model-weighted depth retention for stored fresh-injection groups."""

    if selected_ids is not None:
        raise ValueError(
            "Tier-1 v1 forbids arbitrary injection-ID selection; provide complete shards"
        )
    quality_reference = load_external_quality_reference(
        table_path=cadence_reference_path,
        manifest_path=cadence_reference_manifest_path,
        sector=int(sector),
    )
    rows: list[dict[str, Any]] = []
    contracts: set[str] = set()
    seen_ids: set[str] = set()
    duplicate_ids: set[str] = set()
    file_hashes: dict[str, str] = {}
    shard_index_hashes: dict[str, str] = {}
    shard_indices: list[int] = []
    source_selection_modes: set[str] = set()
    source_adp_paths: set[str] = set()
    declared_counts: dict[str, int] = {}
    actual_counts: dict[str, int] = {}
    malformed_shards: list[str] = []
    quality_totals = {name: 0 for name in INJECTION_QUALITY_COUNT_NAMES}
    quality_by_shard: dict[str, dict[str, int]] = {}
    n_quality_groups_applied = 0
    for path in sorted(Path(item) for item in shard_paths):
        file_hashes[str(path)] = file_sha256(path)
        with h5py.File(path, "r") as h5:
            contracts.add(str(h5.attrs.get("contract_version", "")))
            if "injections" not in h5:
                malformed_shards.append(f"{path}: missing /injections")
                continue
            try:
                shard_index = int(h5.attrs["shard_index"])
                declared_count = int(h5.attrs["n_injections"])
                source_mode = str(h5.attrs["selection_mode"])
            except (KeyError, TypeError, ValueError) as exc:
                malformed_shards.append(f"{path}: invalid root attributes ({exc})")
                continue
            shard_indices.append(shard_index)
            shard_key = str(shard_index)
            if shard_key in shard_index_hashes:
                malformed_shards.append(
                    f"{path}: duplicate shard_index {shard_index}"
                )
            else:
                shard_index_hashes[shard_key] = file_hashes[str(path)]
            source_selection_modes.add(source_mode)
            source_adp_paths.add(str(h5.attrs.get("source_adp_h5", "")))
            declared_counts[str(path)] = declared_count
            actual_counts[str(path)] = len(h5["injections"])
            shard_quality = {name: 0 for name in INJECTION_QUALITY_COUNT_NAMES}
            shard_quality["n_groups_applied"] = 0
            quality_by_shard[str(path)] = shard_quality
            for injection_id, group in h5["injections"].items():
                if str(injection_id) in seen_ids:
                    duplicate_ids.add(str(injection_id))
                seen_ids.add(str(injection_id))
                camera = int(group.attrs.get("camera", -1))
                ccd = int(group.attrs.get("ccd", -1))
                group_sector = int(group.attrs.get("sector", sector))
                required_datasets = (
                    "quality",
                    "cadenceno",
                    "orbitid",
                    "transit_model",
                )
                missing_datasets = [
                    name for name in required_datasets if name not in group
                ]
                quality_overlay = None
                quality_error = ""
                if missing_datasets:
                    quality_error = f"missing datasets {missing_datasets}"
                else:
                    quality = np.asarray(group["quality"], dtype=np.int64)
                    model = np.asarray(group["transit_model"], dtype=float)
                    try:
                        quality_overlay = quality_reference.apply(
                            sector=group_sector,
                            camera=camera,
                            ccd=ccd,
                            cadenceno=np.asarray(group["cadenceno"]),
                            orbitid=np.asarray(group["orbitid"]),
                            internal_quality=quality,
                            context=f"injection {injection_id}",
                        )
                        if len(model) != len(quality_overlay.quality):
                            raise ValueError(
                                "transit_model length disagrees with cadence arrays"
                            )
                    except (KeyError, TypeError, ValueError) as exc:
                        quality_error = str(exc)
                if quality_overlay is None:
                    malformed_shards.append(
                        f"{path}:/injections/{injection_id}: {quality_error}"
                    )
                    x = np.asarray([], dtype=float)
                else:
                    x = 1.0 - model
                    for name in INJECTION_QUALITY_COUNT_NAMES:
                        value = int(quality_overlay.counts[name])
                        quality_totals[name] += value
                        shard_quality[name] += value
                    n_quality_groups_applied += 1
                    shard_quality["n_groups_applied"] += 1
                for aperture in apertures:
                    original_name = f"{aperture}_original"
                    injected_name = f"{aperture}_injected"
                    base = {
                        "injection_id": str(injection_id),
                        "tic": int(group.attrs.get("tic", -1)),
                        "tmag": float(group.attrs.get("tmag", group.attrs.get("tessmag", np.nan))),
                        "aperture": str(aperture),
                        "camera": camera,
                        "ccd": ccd,
                        "detector": f"cam{camera}_ccd{ccd}",
                        "period_d": float(group.attrs.get("period_d", np.nan)),
                        "duration_min": float(group.attrs.get("duration_min", np.nan)),
                        "model_depth": float(group.attrs.get("model_depth", np.nan)),
                        "radius_rearth": float(group.attrs.get("radius_rearth", np.nan)),
                        "shard_index": shard_index,
                        "source_shard": str(path),
                    }
                    if quality_overlay is None:
                        rows.append(
                            {
                                **base,
                                "status": "quality_overlay_error",
                                "quality_error": quality_error,
                            }
                        )
                        continue
                    if original_name not in group or injected_name not in group:
                        rows.append({**base, "status": "missing_dataset"})
                        continue
                    original = np.asarray(group[original_name], dtype=float)
                    injected = np.asarray(group[injected_name], dtype=float)
                    y = -(injected - original)
                    good = (
                        (quality_overlay.quality == 0)
                        & np.isfinite(x)
                        & np.isfinite(y)
                    )
                    in_transit = good & (x > 1.0e-10)
                    if np.count_nonzero(in_transit) < 2 or np.count_nonzero(good) < 20:
                        rows.append({**base, "status": "insufficient_cadences"})
                        continue
                    design = np.column_stack([np.ones(np.count_nonzero(good)), x[good]])
                    intercept, retention = np.linalg.lstsq(design, y[good], rcond=None)[0]
                    prediction = intercept + retention * x[good]
                    residual = y[good] - prediction
                    rows.append(
                        {
                            **base,
                            "status": "ok",
                            "n_quality0": int(np.count_nonzero(good)),
                            "n_in_transit": int(np.count_nonzero(in_transit)),
                            "depth_retention_fraction": float(retention),
                            "fit_intercept": float(intercept),
                            "fit_residual_mad": _robust_mad(residual),
                            "model_delta_correlation": _correlation(x[good], y[good]),
                        }
                    )
    metrics = pd.DataFrame(rows)
    final_file_hashes = {
        str(path): file_sha256(path)
        for path in sorted(Path(item) for item in shard_paths)
    }
    if final_file_hashes != file_hashes:
        raise ValueError("one or more injection shards changed during summarization")
    quality_reference.assert_unchanged()
    selection_sha256 = hashlib.sha256(
        ("\n".join(sorted(seen_ids)) + "\n").encode()
    ).hexdigest()
    summary = {
        "contracts": sorted(contracts),
        "selection_mode": "full_shard",
        "source_selection_modes": sorted(source_selection_modes),
        "source_adp_h5": sorted(source_adp_paths),
        "shard_indices": sorted(shard_indices),
        "n_shards": len(shard_paths),
        "n_injection_ids": len(seen_ids),
        "n_unique_tics": (
            int(pd.to_numeric(metrics["tic"], errors="coerce").nunique())
            if not metrics.empty
            else 0
        ),
        "n_duplicate_injection_ids": len(duplicate_ids),
        "duplicate_injection_id_examples": sorted(duplicate_ids)[:20],
        "declared_injections_by_shard": declared_counts,
        "actual_injections_by_shard": actual_counts,
        "malformed_shards": malformed_shards,
        "selection_sha256": selection_sha256,
        "metadata_sha256": injection_metadata_sha256(metrics),
        "shard_sha256": file_hashes,
        "shard_index_sha256": shard_index_hashes,
        "metrics_content_sha256": _dataframe_content_sha256(metrics),
        "external_quality_overlay": {
            **quality_reference.provenance,
            "applied_before_retention_metrics": True,
            "n_groups_applied": int(n_quality_groups_applied),
            "aggregate_counts": quality_totals,
            "counts_by_shard": quality_by_shard,
        },
    }
    return metrics, summary


def _positive_support_ratio(values: pd.Series) -> float:
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    good = array[np.isfinite(array) & (array > 0)]
    if not len(good):
        return np.nan
    return float(np.max(good) / np.min(good))


def _fixed_injection_quality_overlay_failures(
    manifest: Mapping[str, Any],
    *,
    sector: int,
    config: Tier1QAConfig,
) -> list[str]:
    """Validate authoritative quality provenance and count arithmetic."""

    overlay = manifest.get("external_quality_overlay")
    if not isinstance(overlay, Mapping):
        return ["fixed-injection external-quality overlay is missing"]
    reasons: list[str] = []
    expected = {
        "policy_contract": EXTERNAL_QUALITY_POLICY_CONTRACT,
        "effective_quality_policy": EFFECTIVE_QUALITY_POLICY,
        "sector": int(sector),
        "cadence_reference_contract_version": (
            config.cadence_reference_contract_template.format(sector=int(sector))
        ),
        "cadence_reference_cadence_authority": "qlp_cam_quat",
        "cadence_reference_quality_authority": "spoc_and_qlp_quality_flags",
        "cadence_reference_table_sha256": (
            config.expected_cadence_reference_sha256
        ),
        "cadence_reference_manifest_sha256": (
            config.expected_cadence_reference_manifest_sha256
        ),
    }
    for name, expected_value in expected.items():
        if overlay.get(name) != expected_value:
            reasons.append(
                f"fixed-injection external-quality overlay has incompatible {name}"
            )
    source_digest = str(
        overlay.get("cadence_reference_source_declaration_sha256", "")
    )
    if re.fullmatch(r"[0-9a-f]{64}", source_digest) is None:
        reasons.append("fixed-injection quality source-declaration hash is invalid")
    if overlay.get("applied_before_retention_metrics") is not True:
        reasons.append(
            "fixed-injection external quality was not applied before retention metrics"
        )

    def parse_counts(value: Any, *, context: str) -> dict[str, int] | None:
        if not isinstance(value, Mapping):
            reasons.append(f"{context} is missing or malformed")
            return None
        parsed: dict[str, int] = {}
        if set(value) != set(INJECTION_QUALITY_COUNT_NAMES):
            reasons.append(f"{context} has wrong count fields")
            return None
        for name in INJECTION_QUALITY_COUNT_NAMES:
            try:
                parsed[name] = int(value[name])
            except (TypeError, ValueError):
                reasons.append(f"{context}.{name} is invalid")
                return None
        if any(value < 0 for value in parsed.values()):
            reasons.append(f"{context} contains negative counts")
        total = parsed["n_cad_total"]
        if any(
            parsed[name] > total
            for name in INJECTION_QUALITY_COUNT_NAMES
            if name != "n_cad_total"
        ):
            reasons.append(f"{context} has counts above n_cad_total")
        if parsed["n_cad_external_only_bad"] > parsed["n_cad_external_bad"]:
            reasons.append(f"{context} has impossible external-only count")
        if parsed["n_cad_authority_excluded"] > parsed["n_cad_external_bad"]:
            reasons.append(f"{context} has impossible authority-excluded count")
        if parsed["n_cad_effective_bad"] != (
            parsed["n_cad_internal_bad"] + parsed["n_cad_external_only_bad"]
        ):
            reasons.append(f"{context} has inconsistent effective-bad arithmetic")
        return parsed

    aggregate = parse_counts(
        overlay.get("aggregate_counts"), context="fixed-injection aggregate quality"
    )
    per_shard = overlay.get("counts_by_shard")
    shard_hashes = manifest.get("shard_sha256")
    if not isinstance(per_shard, Mapping) or not isinstance(shard_hashes, Mapping):
        reasons.append("fixed-injection per-shard quality counts are missing")
        per_shard = {}
        shard_hashes = {}
    elif set(str(value) for value in per_shard) != set(
        str(value) for value in shard_hashes
    ):
        reasons.append("fixed-injection per-shard quality keys disagree with shard hashes")
    summed = {name: 0 for name in INJECTION_QUALITY_COUNT_NAMES}
    summed_groups = 0
    for shard_path, value in per_shard.items():
        if not isinstance(value, Mapping):
            reasons.append(f"fixed-injection shard quality {shard_path} is malformed")
            continue
        counts = parse_counts(
            {name: value.get(name) for name in INJECTION_QUALITY_COUNT_NAMES},
            context=f"fixed-injection shard quality {shard_path}",
        )
        try:
            groups = int(value.get("n_groups_applied", -1))
        except (TypeError, ValueError):
            groups = -1
        if groups < 0:
            reasons.append(
                f"fixed-injection shard quality {shard_path}.n_groups_applied is invalid"
            )
        else:
            summed_groups += groups
        if counts is not None:
            for name in INJECTION_QUALITY_COUNT_NAMES:
                summed[name] += counts[name]
    try:
        declared_groups = int(overlay.get("n_groups_applied", -1))
        n_injections = int(manifest.get("n_injection_ids", -1))
    except (TypeError, ValueError):
        declared_groups = n_injections = -1
    if declared_groups != n_injections or summed_groups != declared_groups:
        reasons.append("fixed-injection quality overlay does not cover every injection")
    if aggregate is not None and summed != aggregate:
        reasons.append("fixed-injection aggregate/per-shard quality counts disagree")
    return reasons


def evaluate_fixed_injections(
    metrics: pd.DataFrame,
    manifest: Mapping[str, Any],
    *,
    sector: int,
    config: Tier1QAConfig,
) -> dict[str, Any]:
    required = {
        "injection_id",
        "tic",
        "tmag",
        "aperture",
        "status",
        "depth_retention_fraction",
        "camera",
        "ccd",
        "detector",
        "period_d",
        "duration_min",
        "model_depth",
        "shard_index",
    }
    missing = sorted(required - set(metrics.columns))
    expected_contract = config.injection_contract_template.format(sector=int(sector))
    observed_contracts = manifest.get("contracts")
    if observed_contracts is None:
        observed = manifest.get("contract_version", manifest.get("contract", ""))
        observed_contracts = [str(observed)] if observed else []
    observed_contracts = sorted(str(value) for value in observed_contracts)
    reasons: list[str] = []
    reasons.extend(
        _fixed_injection_quality_overlay_failures(
            manifest,
            sector=sector,
            config=config,
        )
    )
    if missing:
        reasons.append(f"missing columns: {missing}")
    if observed_contracts != [expected_contract]:
        reasons.append("fresh-injection contract mismatch")
    if str(manifest.get("selection_mode", "")) != config.injection_selection_mode:
        reasons.append("fixed injection selection is not complete-shard evidence")
    if sorted(manifest.get("source_selection_modes", [])) != ["shard"]:
        reasons.append("source injection HDF5 selection_mode is not 'shard'")
    observed_shards = sorted(int(value) for value in manifest.get("shard_indices", []))
    if observed_shards != sorted(config.injection_required_shard_indices):
        reasons.append("frozen injection shard-index set mismatch")
    if int(manifest.get("n_shards", -1)) != len(config.injection_required_shard_indices):
        reasons.append("frozen injection shard count mismatch")
    if int(manifest.get("n_injection_ids", -1)) != config.injection_expected_ids:
        reasons.append("frozen injection ID count mismatch")
    if int(manifest.get("n_unique_tics", -1)) != config.injection_expected_ids:
        reasons.append("frozen injections do not use unique hosts")
    if int(manifest.get("n_duplicate_injection_ids", -1)) != 0:
        reasons.append("duplicate injection IDs were observed across shards")
    if manifest.get("malformed_shards", ["missing manifest field"]):
        reasons.append("one or more injection shards are malformed")
    declared = manifest.get("declared_injections_by_shard", {})
    actual = manifest.get("actual_injections_by_shard", {})
    if not isinstance(declared, Mapping) or not isinstance(actual, Mapping) or declared != actual:
        reasons.append("declared and actual injection shard sizes differ")
    expected_per_shard = config.injection_expected_ids // len(
        config.injection_required_shard_indices
    )
    if (
        config.injection_expected_ids % len(config.injection_required_shard_indices)
        or not isinstance(actual, Mapping)
        or len(actual) != len(config.injection_required_shard_indices)
        or any(int(value) != expected_per_shard for value in actual.values())
    ):
        reasons.append("injection shards do not have the locked equal size")
    if not isinstance(manifest.get("shard_sha256"), Mapping) or len(
        manifest.get("shard_sha256", {})
    ) != len(config.injection_required_shard_indices):
        reasons.append("injection shard hashes are incomplete")
    expected_shard_hashes = {
        str(index): digest
        for index, digest in zip(
            config.injection_required_shard_indices,
            config.expected_injection_shard_sha256,
            strict=True,
        )
    }
    observed_shard_hashes = manifest.get("shard_index_sha256")
    if not isinstance(observed_shard_hashes, Mapping) or {
        str(key): str(value) for key, value in observed_shard_hashes.items()
    } != expected_shard_hashes:
        reasons.append("injection shard files are not pinned by the locked config")
    if str(manifest.get("metrics_content_sha256", "")) != _dataframe_content_sha256(
        metrics
    ):
        reasons.append("injection metrics content hash mismatch")
    if reasons:
        return {
            "status": "fail",
            "reasons": reasons,
            "expected_contract": expected_contract,
            "observed_contracts": observed_contracts,
        }
    selected = metrics.loc[
        metrics["aperture"].astype(str).isin(config.apertures)
    ].copy()
    selected["injection_id"] = selected["injection_id"].astype(str)
    pair_counts = selected.groupby(["injection_id", "aperture"], sort=False).size()
    duplicate_pairs = pair_counts.loc[pair_counts.ne(1)]
    expected_ids = set(selected["injection_id"])
    observed_selection_sha256 = hashlib.sha256(
        ("\n".join(sorted(expected_ids)) + "\n").encode()
    ).hexdigest()
    selection_digest_matches = bool(
        str(manifest.get("selection_sha256", ""))
        == observed_selection_sha256
        == config.expected_injection_selection_sha256
    )
    observed_pairs = set(
        zip(selected["injection_id"], selected["aperture"].astype(str), strict=False)
    )
    expected_pairs = {
        (injection_id, aperture)
        for injection_id in expected_ids
        for aperture in config.apertures
    }
    missing_pairs = expected_pairs - observed_pairs
    work = selected.loc[selected["status"].astype(str).eq("ok")].copy()
    valid_pairs = set(
        zip(work["injection_id"], work["aperture"].astype(str), strict=False)
    )
    finite_retention = np.isfinite(
        pd.to_numeric(work["depth_retention_fraction"], errors="coerce")
    )
    finite_pairs = set(
        zip(
            work.loc[finite_retention, "injection_id"],
            work.loc[finite_retention, "aperture"].astype(str),
            strict=False,
        )
    )
    invalid_pairs = expected_pairs - (valid_pairs & finite_pairs)
    exact_pair_coverage = not duplicate_pairs.size and not missing_pairs and not invalid_pairs
    n_ids = int(len(expected_ids))
    per_injection = selected.drop_duplicates("injection_id").copy()
    try:
        observed_metadata_sha256 = injection_metadata_sha256(selected)
    except (KeyError, TypeError, ValueError):
        observed_metadata_sha256 = ""
    metadata_digest_matches = bool(
        str(manifest.get("metadata_sha256", ""))
        == observed_metadata_sha256
        == config.expected_injection_metadata_sha256
    )
    faint_ids = int(
        per_injection.loc[
            pd.to_numeric(per_injection["tmag"], errors="coerce")
            >= config.injection_faint_tmag_min
        ].shape[0]
    )
    detector = per_injection["detector"].astype(str).str.strip()
    camera = pd.to_numeric(per_injection["camera"], errors="coerce")
    ccd = pd.to_numeric(per_injection["ccd"], errors="coerce")
    valid_detector = pd.Series(
        [
            re.fullmatch(r"cam[1-4]_ccd[1-4]", value) is not None
            and value == f"cam{int(camera_value)}_ccd{int(ccd_value)}"
            if np.isfinite(camera_value) and np.isfinite(ccd_value)
            else False
            for value, camera_value, ccd_value in zip(
                detector, camera, ccd, strict=True
            )
        ],
        index=per_injection.index,
    )
    n_detectors = int(detector.loc[valid_detector].nunique())
    period_ratio = _positive_support_ratio(per_injection["period_d"])
    duration_ratio = _positive_support_ratio(per_injection["duration_min"])
    depth_ratio = _positive_support_ratio(per_injection["model_depth"])
    positive_parameter_support = all(
        (
            np.isfinite(pd.to_numeric(per_injection[column], errors="coerce"))
            & (pd.to_numeric(per_injection[column], errors="coerce") > 0)
        ).all()
        for column in ("period_d", "duration_min", "model_depth")
    )
    model_depth = pd.to_numeric(per_injection["model_depth"], errors="coerce")
    tic = pd.to_numeric(per_injection["tic"], errors="coerce")
    tmag = pd.to_numeric(per_injection["tmag"], errors="coerce")
    shard_index = pd.to_numeric(per_injection["shard_index"], errors="coerce")
    valid_locked_metadata = bool(
        positive_parameter_support
        and (model_depth <= 1).all()
        and np.isfinite(tic).all()
        and (tic > 0).all()
        and (tic == np.floor(tic)).all()
        and np.isfinite(tmag).all()
        and np.isfinite(camera).all()
        and camera.between(1, 4).all()
        and (camera == np.floor(camera)).all()
        and np.isfinite(ccd).all()
        and ccd.between(1, 4).all()
        and (ccd == np.floor(ccd)).all()
        and np.isfinite(shard_index).all()
        and set(shard_index.astype(int))
        == set(config.injection_required_shard_indices)
    )
    metadata_consistent = True
    for column in ("tic", "tmag", "detector", "period_d", "duration_min", "model_depth", "shard_index"):
        if selected.groupby("injection_id")[column].nunique(dropna=False).gt(1).any():
            metadata_consistent = False
            break
    aperture_summaries: list[dict[str, Any]] = []
    statuses: list[str] = []
    for aperture in config.apertures:
        rows = work.loc[work["aperture"].astype(str).eq(aperture)]
        retention = pd.to_numeric(rows["depth_retention_fraction"], errors="coerce")
        median = _quantile(retention, 0.5)
        p10 = _quantile(retention, 0.1)
        finite = retention[np.isfinite(retention)]
        inband = (
            float(
                np.mean(
                    (finite >= config.injection_inband_low)
                    & (finite <= config.injection_inband_high)
                )
            )
            if len(finite)
            else np.nan
        )
        components = {
            "median_retention": _status_interval(
                median,
                config.injection_pass_retention_low,
                config.injection_pass_retention_high,
                config.injection_review_retention_low,
                config.injection_review_retention_high,
            ),
            "p10_retention": _status_min(
                p10, config.injection_pass_p10, config.injection_review_p10
            ),
            "inband_fraction": _status_min(
                inband,
                config.injection_pass_inband_fraction,
                config.injection_review_inband_fraction,
            ),
        }
        status = _worst_status(components.values())
        statuses.append(status)
        aperture_summaries.append(
            {
                "aperture": aperture,
                "status": status,
                "n_rows": int(len(rows)),
                "median_depth_retention": median,
                "p10_depth_retention": p10,
                "inband_fraction": inband,
                "component_status": components,
            }
        )
    coverage_status = (
        "pass"
        if n_ids == config.injection_expected_ids
        and int(per_injection["tic"].nunique()) == config.injection_expected_ids
        and faint_ids >= config.min_faint_injection_ids
        and n_detectors >= config.min_injection_detectors
        and valid_detector.all()
        and positive_parameter_support
        and selection_digest_matches
        and metadata_digest_matches
        and valid_locked_metadata
        and period_ratio >= config.min_injection_period_ratio
        and duration_ratio >= config.min_injection_duration_ratio
        and depth_ratio >= config.min_injection_depth_ratio
        and metadata_consistent
        and exact_pair_coverage
        else "fail"
    )
    statuses.append(coverage_status)
    status = _worst_status(statuses)
    return {
        "status": status,
        "reasons": ["fixed injection coverage or retention failed"] if status == "fail" else [],
        "expected_contract": expected_contract,
        "observed_contracts": observed_contracts,
        "n_injection_ids": n_ids,
        "n_faint_injection_ids": faint_ids,
        "selection_sha256": observed_selection_sha256,
        "selection_digest_matches_locked_config": selection_digest_matches,
        "metadata_sha256": observed_metadata_sha256,
        "metadata_digest_matches_locked_config": metadata_digest_matches,
        "n_unique_tics": int(per_injection["tic"].nunique()),
        "n_detectors": n_detectors,
        "period_support_ratio": period_ratio,
        "duration_support_ratio": duration_ratio,
        "depth_support_ratio": depth_ratio,
        "positive_finite_parameter_support": positive_parameter_support,
        "valid_locked_metadata": valid_locked_metadata,
        "metadata_consistent_across_apertures": metadata_consistent,
        "exact_injection_aperture_coverage": exact_pair_coverage,
        "n_missing_injection_aperture_pairs": len(missing_pairs),
        "n_invalid_injection_aperture_pairs": len(invalid_pairs),
        "n_duplicate_injection_aperture_pairs": int(len(duplicate_pairs)),
        "coverage_status": coverage_status,
        "apertures": aperture_summaries,
    }


def _wd_epoch_residual_minutes(t0_bjd: float) -> float:
    phase = (float(t0_bjd) - WD1856_T0_BJD) / WD1856_PERIOD_D
    return float(abs((phase + 0.5) % 1.0 - 0.5) * WD1856_PERIOD_D * 1440.0)


_INDEPENDENT_QUALITY_COUNT_KEYS = {
    "n_cad_total",
    "n_cad_internal_bad",
    "n_cad_external_bad",
    "n_cad_authority_excluded",
    "n_cad_external_only_bad",
    "n_cad_effective_bad",
}


def _validate_independent_quality_overlay(
    manifest: Mapping[str, Any],
    config: Tier1QAConfig,
    *,
    sector: int,
) -> list[str]:
    """Validate quality-mask provenance and arithmetic for WD 1856 evidence."""

    reasons: list[str] = []
    overlay = manifest.get("external_quality_overlay")
    if not isinstance(overlay, Mapping):
        return ["independent external-quality overlay is missing or malformed"]
    exact_values = {
        "policy_contract": EXTERNAL_QUALITY_POLICY_CONTRACT,
        "effective_quality_policy": EFFECTIVE_QUALITY_POLICY,
        "sector": int(sector),
        "cadence_reference_contract_version": (
            f"s{int(sector)}_a2v1_cadence_reference_v1"
        ),
        "cadence_reference_cadence_authority": "qlp_cam_quat",
        "cadence_reference_quality_authority": "spoc_and_qlp_quality_flags",
        "cadence_reference_table_sha256": (
            config.expected_cadence_reference_sha256
        ),
        "cadence_reference_manifest_sha256": (
            config.expected_cadence_reference_manifest_sha256
        ),
    }
    for field, expected in exact_values.items():
        if overlay.get(field) != expected:
            reasons.append(
                f"independent external-quality overlay has incompatible {field}"
            )
    source_digest = str(
        overlay.get("cadence_reference_source_declaration_sha256", "")
    )
    if not re.fullmatch(r"[0-9a-f]{64}", source_digest):
        reasons.append(
            "independent external-quality source-declaration hash is invalid"
        )
    if overlay.get("applied_before_common_cadence_metrics") is not True:
        reasons.append(
            "independent external quality was not applied before cadence metrics"
        )
    cadence_policy = str(manifest.get("cadence_match_policy", "")).lower()
    if not all(token in cadence_policy for token in ("effective", "internal", "external")):
        reasons.append(
            "independent cadence policy does not declare effective internal/external quality"
        )

    parsed: dict[str, dict[str, int]] = {}
    count_fields = (
        "current_compact_full_audit_counts",
        "current_common_audit_counts",
        "reference_common_audit_counts",
    )
    for field in count_fields:
        raw_counts = overlay.get(field)
        if not isinstance(raw_counts, Mapping) or set(raw_counts) != (
            _INDEPENDENT_QUALITY_COUNT_KEYS
        ):
            reasons.append(f"independent quality audit {field} has wrong keys")
            continue
        counts: dict[str, int] = {}
        valid = True
        for name in _INDEPENDENT_QUALITY_COUNT_KEYS:
            value = raw_counts.get(name)
            if type(value) is not int or value < 0:
                reasons.append(
                    f"independent quality audit {field}.{name} is invalid"
                )
                valid = False
            else:
                counts[name] = int(value)
        if not valid:
            continue
        total = counts["n_cad_total"]
        bounded = all(
            counts[name] <= total
            for name in (
                "n_cad_internal_bad",
                "n_cad_external_bad",
                "n_cad_authority_excluded",
                "n_cad_external_only_bad",
                "n_cad_effective_bad",
            )
        )
        arithmetic = (
            counts["n_cad_external_only_bad"]
            <= counts["n_cad_external_bad"]
            and counts["n_cad_authority_excluded"]
            <= counts["n_cad_external_bad"]
            and counts["n_cad_effective_bad"]
            == counts["n_cad_internal_bad"]
            + counts["n_cad_external_only_bad"]
        )
        if not bounded or not arithmetic:
            reasons.append(f"independent quality audit {field} is inconsistent")
        parsed[field] = counts

    try:
        n_current = int(manifest.get("n_current_cadences", -1))
        n_common = int(manifest.get("n_common_cadences", -1))
    except (TypeError, ValueError):
        n_current, n_common = -1, -1
    full = parsed.get("current_compact_full_audit_counts")
    current_common = parsed.get("current_common_audit_counts")
    reference_common = parsed.get("reference_common_audit_counts")
    if full is not None and full["n_cad_total"] != n_current:
        reasons.append("independent full quality-audit cadence count is inconsistent")
    for field, counts in (
        ("current", current_common),
        ("reference", reference_common),
    ):
        if counts is not None and counts["n_cad_total"] != n_common:
            reasons.append(
                f"independent {field} common quality-audit count is inconsistent"
            )
    if (
        current_common is not None
        and reference_common is not None
        and current_common["n_cad_external_bad"]
        != reference_common["n_cad_external_bad"]
    ):
        reasons.append(
            "independent current/reference common external-quality counts disagree"
        )
    return reasons


def evaluate_independent_extraction(
    metrics: pd.DataFrame,
    manifest: Mapping[str, Any],
    config: Tier1QAConfig,
    *,
    sector: int | None = None,
    catalog_detectors: Mapping[int, str] | None = None,
    current_compact_sha256: str | None = None,
) -> dict[str, Any]:
    sector = config.sector if sector is None else int(sector)
    current_compact_sha256 = (
        config.expected_compact_sha256
        if current_compact_sha256 is None
        else str(current_compact_sha256)
    )
    required_manifest = {
        "contract_version",
        "sector",
        "tic",
        "independent",
        "comparison_mode",
        "independence_basis",
        "current_extractor_family",
        "reference_extractor_family",
        "current_repository",
        "reference_repository",
        "current_revision",
        "reference_revision",
        "pixel_source",
        "cadence_match_policy",
        "scatter_definition",
        "depth_definition",
        "ephemeris_recovery_definition",
        "current_apertures",
        "reference_flux_columns",
        "reference_columns",
        "reference_product_aperture",
        "reference_target_id",
        "current_compact_sha256",
        "current_product_sha256",
        "reference_product_sha256",
        "reference_table_sha256",
        "reference_table_identity",
        "reference_product_identity",
        "fixed_ephemeris",
        "bounded_bls",
        "max_time_delta_seconds_allowed",
        "metrics_file_sha256",
        "metrics_content_sha256",
        "n_current_cadences",
        "n_reference_cadences",
        "n_common_cadences",
        "external_quality_overlay",
    }
    required_metrics = {
        "tic",
        "sector",
        "aperture",
        "reference_flux_column",
        "detector",
        "n_current_cadences",
        "n_reference_cadences",
        "n_common_cadences",
        "current_scatter_ppm",
        "reference_scatter_ppm",
        "current_period_d",
        "reference_period_d",
        "current_t0_bjd",
        "reference_t0_bjd",
        "current_depth",
        "reference_depth",
        "current_bls_duration_min",
        "reference_bls_duration_min",
        "current_bls_depth",
        "reference_bls_depth",
        "current_bls_depth_snr",
        "reference_bls_depth_snr",
        "current_bls_power",
        "reference_bls_power",
        "n_common_quality0_finite",
        "n_in_event_cadences",
        "n_out_of_event_cadences",
        "max_abs_time_delta_s",
    }
    reasons: list[str] = []
    missing_manifest = sorted(required_manifest - set(manifest))
    missing_metrics = sorted(required_metrics - set(metrics.columns))
    if missing_manifest:
        reasons.append(f"independent manifest missing fields: {missing_manifest}")
    if missing_metrics:
        reasons.append(f"independent metrics missing columns: {missing_metrics}")
    if reasons:
        return {"status": "fail", "reasons": reasons}
    expected_contract = config.independent_contract_template.format(sector=sector)
    if str(manifest["contract_version"]) != expected_contract:
        reasons.append("independent-extraction contract mismatch")
    try:
        manifest_sector = int(manifest["sector"])
    except (TypeError, ValueError):
        manifest_sector = -1
    if manifest_sector != sector:
        reasons.append("independent-extraction sector mismatch")
    try:
        manifest_tic = int(manifest["tic"])
    except (TypeError, ValueError):
        manifest_tic = -1
    if manifest_tic != WD1856_TIC:
        reasons.append("independent-extraction target is not WD 1856")
    if str(manifest["comparison_mode"]) != config.independent_comparison_mode:
        reasons.append("independent-extraction comparison mode mismatch")
    reasons.extend(
        _validate_independent_quality_overlay(manifest, config, sector=sector)
    )
    if str(manifest["current_compact_sha256"]) != current_compact_sha256:
        reasons.append("independent evidence is not bound to the current compact product")
    if str(manifest["current_product_sha256"]) != current_compact_sha256:
        reasons.append("independent current-product hash is not the compact product")
    if (
        str(manifest["reference_product_sha256"])
        != config.expected_independent_reference_product_sha256
    ):
        reasons.append("independent reference product is not pinned by the locked config")
    if str(manifest["metrics_content_sha256"]) != _dataframe_content_sha256(metrics):
        reasons.append("independent metrics content hash mismatch")
    if tuple(manifest["current_apertures"]) != tuple(config.apertures):
        reasons.append("independent comparison aperture set mismatch")
    reference_flux_columns = manifest["reference_flux_columns"]
    normalized_reference_flux_columns: dict[str, str] = {}
    if not isinstance(reference_flux_columns, Mapping):
        reasons.append("independent reference_flux_columns is not a mapping")
    else:
        normalized_reference_flux_columns = {
            str(aperture): str(column).strip()
            for aperture, column in reference_flux_columns.items()
        }
        if set(normalized_reference_flux_columns) != set(config.apertures):
            reasons.append(
                "independent reference flux mapping does not cover the locked "
                "aperture set"
            )
        if any(not column for column in normalized_reference_flux_columns.values()):
            reasons.append("independent reference flux mapping has an empty column")
        if len(set(normalized_reference_flux_columns.values())) != len(
            normalized_reference_flux_columns
        ):
            reasons.append(
                "independent reference flux mapping reuses one external column"
            )
    reference_columns = manifest["reference_columns"]
    if not isinstance(reference_columns, Mapping):
        reasons.append("independent reference_columns is not a mapping")
    else:
        nested_mapping = reference_columns.get("flux_by_current_aperture")
        if not isinstance(nested_mapping, Mapping) or {
            str(key): str(value).strip() for key, value in nested_mapping.items()
        } != normalized_reference_flux_columns:
            reasons.append(
                "independent reference column provenance disagrees with its "
                "aperture mapping"
            )
    for key in (
        "current_extractor_family",
        "reference_extractor_family",
        "current_repository",
        "reference_repository",
        "current_revision",
        "reference_revision",
        "pixel_source",
        "cadence_match_policy",
        "scatter_definition",
        "depth_definition",
        "ephemeris_recovery_definition",
        "independence_basis",
        "reference_product_aperture",
        "reference_target_id",
        "current_product_sha256",
        "reference_product_sha256",
        "reference_table_sha256",
    ):
        if not str(manifest[key]).strip():
            reasons.append(f"independent manifest field is empty: {key}")
    current_family = str(manifest["current_extractor_family"]).strip().lower()
    reference_family = str(manifest["reference_extractor_family"]).strip().lower()
    independent = manifest["independent"] is True
    if not independent or current_family == reference_family:
        reasons.append("extractor families are not independent")
    if any(token in reference_family for token in ("tglc", "twirl")):
        reasons.append("reference extraction is from the TGLC/TWIRL family")
    if str(manifest["current_repository"]) == str(manifest["reference_repository"]):
        reasons.append("reference uses the current extraction repository")
    if str(manifest["current_product_sha256"]) == str(manifest["reference_product_sha256"]):
        reasons.append("current and reference product hashes are identical")
    reference_target_digits = re.sub(r"\D", "", str(manifest["reference_target_id"]))
    if reference_target_digits != str(WD1856_TIC):
        reasons.append("independent reference_target_id is not WD 1856")
    for identity_key in ("reference_table_identity", "reference_product_identity"):
        identity = manifest[identity_key]
        if not isinstance(identity, Mapping):
            reasons.append(f"independent {identity_key} is not a mapping")
            continue
        try:
            identity_tic = int(identity.get("tic", -1))
            identity_sector = int(identity.get("sector", -1))
        except (TypeError, ValueError):
            identity_tic, identity_sector = -1, -1
        if identity_tic != WD1856_TIC or identity_sector != sector:
            reasons.append(f"independent {identity_key} has wrong target or sector")
        if not str(identity.get("format", "")).strip() or not isinstance(
            identity.get("identity_source"), Mapping
        ):
            reasons.append(f"independent {identity_key} lacks identity provenance")
    bounded_bls = manifest["bounded_bls"]
    if not isinstance(bounded_bls, Mapping):
        reasons.append("independent bounded_bls definition is not a mapping")
        bounded_bls = {}
    try:
        bls_half_width = float(bounded_bls.get("period_relative_half_width", np.nan))
        bls_n_periods = int(bounded_bls.get("n_periods", -1))
        bls_oversample = int(bounded_bls.get("oversample", -1))
        bls_min_depth_snr = float(bounded_bls.get("min_depth_snr", np.nan))
        bls_durations = np.asarray(bounded_bls.get("durations_min", []), dtype=float)
    except (TypeError, ValueError):
        bls_half_width = np.nan
        bls_n_periods = -1
        bls_oversample = -1
        bls_min_depth_snr = np.nan
        bls_durations = np.asarray([], dtype=float)
    if not (
        np.isfinite(bls_half_width)
        and 0.005 <= bls_half_width <= 0.05
        and bls_n_periods >= 1_001
        and bls_oversample >= 5
        and np.isfinite(bls_min_depth_snr)
        and bls_min_depth_snr >= 5.0
        and bls_durations.size > 0
        and np.isfinite(bls_durations).all()
        and (bls_durations > 0).all()
    ):
        reasons.append("independent bounded BLS definition is invalid")
    ephemeris_definition = str(manifest["ephemeris_recovery_definition"]).lower()
    if "bounded" not in ephemeris_definition or "boxleastsquares" not in ephemeris_definition:
        reasons.append("independent ephemeris recovery is not declared as bounded BLS")
    fixed_ephemeris = manifest["fixed_ephemeris"]
    if not isinstance(fixed_ephemeris, Mapping):
        reasons.append("independent fixed_ephemeris definition is not a mapping")
    else:
        try:
            fixed_period = float(fixed_ephemeris.get("period_d", np.nan))
            fixed_t0 = float(fixed_ephemeris.get("t0_bjd", np.nan))
        except (TypeError, ValueError):
            fixed_period, fixed_t0 = np.nan, np.nan
        if not (
            np.isfinite(fixed_period)
            and np.isfinite(fixed_t0)
            and abs(fixed_period - WD1856_PERIOD_D) <= 1.0e-12
            and abs(fixed_t0 - WD1856_T0_BJD) <= 1.0e-9
        ):
            reasons.append("independent fixed ephemeris is not the locked WD 1856 value")
    if reasons:
        return {"status": "fail", "reasons": reasons}
    work = metrics.copy()
    if work.empty or work.duplicated(["tic", "aperture"]).any():
        return {
            "status": "fail",
            "reasons": ["independent metrics are empty or contain duplicate TIC/aperture rows"],
        }
    metric_sector = pd.to_numeric(work["sector"], errors="coerce")
    if not metric_sector.eq(sector).all():
        reasons.append("independent metric rows have the wrong sector")
    observed_apertures = set(work["aperture"].astype(str))
    if observed_apertures != set(config.apertures):
        reasons.append("independent metric rows do not cover the locked aperture set")
    metric_apertures = work["aperture"].astype(str)
    metric_reference_columns = work["reference_flux_column"].astype(str).str.strip()
    expected_metric_reference_columns = metric_apertures.map(
        normalized_reference_flux_columns
    )
    if (
        expected_metric_reference_columns.isna().any()
        or not metric_reference_columns.eq(expected_metric_reference_columns).all()
    ):
        reasons.append(
            "independent metric aperture mapping disagrees with the manifest"
        )
    tic_numeric = pd.to_numeric(work["tic"], errors="coerce")
    if not np.isfinite(tic_numeric).all() or (tic_numeric <= 0).any():
        reasons.append("independent metric TIC identifiers are invalid")
    metric_tics = set(tic_numeric.astype(int))
    if catalog_detectors is None:
        reasons.append("compact-catalog detector binding was not supplied")
    elif not metric_tics <= set(catalog_detectors):
        reasons.append("independent metrics contain TICs outside the compact catalog")
    detector = work["detector"].astype(str).str.strip()
    valid_detector_format = detector.map(
        lambda value: re.fullmatch(r"cam[1-4]_ccd[1-4]", value) is not None
    )
    if not valid_detector_format.all():
        reasons.append("independent metric detector identifiers are invalid")
    if work.assign(_detector=detector).groupby("tic")["_detector"].nunique().gt(1).any():
        reasons.append("independent detector differs across apertures for one TIC")
    if catalog_detectors is not None:
        detector_matches_catalog = all(
            str(detector_value) == str(catalog_detectors.get(int(tic_value), ""))
            for tic_value, detector_value in zip(tic_numeric, detector, strict=True)
        )
        if not detector_matches_catalog:
            reasons.append("independent detector disagrees with the compact catalog")
    current_n = pd.to_numeric(work["n_current_cadences"], errors="coerce")
    reference_n = pd.to_numeric(work["n_reference_cadences"], errors="coerce")
    common_n = pd.to_numeric(work["n_common_cadences"], errors="coerce")
    valid_counts = (
        np.isfinite(current_n)
        & np.isfinite(reference_n)
        & np.isfinite(common_n)
        & (current_n > 0)
        & (reference_n > 0)
        & (common_n > 0)
        & (common_n <= np.minimum(current_n, reference_n))
    )
    if not valid_counts.all():
        reasons.append("independent cadence counts are invalid")
    work["current_common_fraction"] = common_n / current_n.replace(0, np.nan)
    work["reference_common_fraction"] = common_n / reference_n.replace(0, np.nan)
    work["common_fraction"] = work[
        ["current_common_fraction", "reference_common_fraction"]
    ].min(axis=1)
    current_scatter = pd.to_numeric(work["current_scatter_ppm"], errors="coerce")
    reference_scatter = pd.to_numeric(work["reference_scatter_ppm"], errors="coerce")
    valid_scatter = (
        np.isfinite(current_scatter)
        & np.isfinite(reference_scatter)
        & (current_scatter > 0)
        & (reference_scatter > 0)
    )
    if not valid_scatter.all():
        reasons.append("independent scatter measurements are invalid")
    diagnostic_columns = (
        "current_bls_duration_min",
        "reference_bls_duration_min",
        "current_bls_depth",
        "reference_bls_depth",
        "current_bls_depth_snr",
        "reference_bls_depth_snr",
        "current_bls_power",
        "reference_bls_power",
        "n_common_quality0_finite",
        "n_in_event_cadences",
        "n_out_of_event_cadences",
        "max_abs_time_delta_s",
    )
    diagnostics = {
        column: pd.to_numeric(work[column], errors="coerce")
        for column in diagnostic_columns
    }
    if any(not np.isfinite(values).all() for values in diagnostics.values()):
        reasons.append("independent BLS/cadence diagnostics contain nonfinite values")
    else:
        valid_bls = (
            (diagnostics["current_bls_duration_min"] > 0)
            & (diagnostics["reference_bls_duration_min"] > 0)
            & diagnostics["current_bls_depth"].between(0, 1, inclusive="neither")
            & diagnostics["reference_bls_depth"].between(0, 1, inclusive="neither")
            & (diagnostics["current_bls_depth_snr"] >= bls_min_depth_snr)
            & (diagnostics["reference_bls_depth_snr"] >= bls_min_depth_snr)
            & (diagnostics["current_bls_power"] > 0)
            & (diagnostics["reference_bls_power"] > 0)
        )
        try:
            max_time_delta_allowed = float(
                manifest["max_time_delta_seconds_allowed"]
            )
        except (TypeError, ValueError):
            max_time_delta_allowed = np.nan
        valid_cadence_diagnostics = (
            (diagnostics["n_common_quality0_finite"] > 0)
            & (diagnostics["n_common_quality0_finite"] <= common_n)
            & (diagnostics["n_in_event_cadences"] > 0)
            & (diagnostics["n_out_of_event_cadences"] > 0)
            & (diagnostics["n_in_event_cadences"] <= diagnostics["n_common_quality0_finite"])
            & (diagnostics["n_out_of_event_cadences"] <= diagnostics["n_common_quality0_finite"])
            & (diagnostics["max_abs_time_delta_s"] >= 0)
            & (diagnostics["max_abs_time_delta_s"] <= max_time_delta_allowed)
        )
        if not valid_bls.all():
            reasons.append("independent bounded-BLS diagnostics are invalid")
        if (
            not np.isfinite(max_time_delta_allowed)
            or max_time_delta_allowed <= 0
            or max_time_delta_allowed > config.independent_max_time_delta_seconds
        ):
            reasons.append("independent cadence-match tolerance is invalid")
        elif not valid_cadence_diagnostics.all():
            reasons.append("independent cadence-match diagnostics are invalid")
    if reasons:
        return {"status": "fail", "reasons": reasons}
    work["abs_log_scatter_ratio"] = np.abs(np.log10(current_scatter / reference_scatter))
    n_targets = int(tic_numeric.astype(int).nunique())
    n_detectors = int(detector.nunique())
    common_median = _quantile(work["common_fraction"], 0.5)
    common_min = _quantile(work["common_fraction"], 0.0)
    current_common_min = _quantile(work["current_common_fraction"], 0.0)
    reference_common_min = _quantile(work["reference_common_fraction"], 0.0)
    scatter_median = _quantile(work["abs_log_scatter_ratio"], 0.5)
    scatter_max = _quantile(work["abs_log_scatter_ratio"], 1.0)
    statuses = {
        "coverage": (
            "pass"
            if n_targets >= config.min_independent_targets
            and n_detectors >= config.min_independent_detectors
            else "fail"
        ),
        "common_cadence_fraction": _status_min(
            common_min,
            config.independent_pass_common_fraction,
            config.independent_review_common_fraction,
        ),
    }
    for column in (
        "current_period_d",
        "reference_period_d",
        "current_t0_bjd",
        "reference_t0_bjd",
        "current_depth",
        "reference_depth",
    ):
        work[column] = pd.to_numeric(work[column], errors="coerce")
    wd = work.loc[tic_numeric.eq(WD1856_TIC)].copy()
    wd_summary: dict[str, Any]
    if wd.empty or set(wd["aperture"].astype(str)) != set(config.apertures):
        statuses["wd1856"] = "fail"
        wd_summary = {"status": "fail", "reason": "WD 1856 aperture evidence missing"}
    else:
        aperture_summaries: list[dict[str, Any]] = []
        aperture_statuses: list[str] = []
        for row in wd.itertuples(index=False):
            current_period_error = abs(float(row.current_period_d) / WD1856_PERIOD_D - 1.0)
            reference_period_error = abs(
                float(row.reference_period_d) / WD1856_PERIOD_D - 1.0
            )
            period_error = max(current_period_error, reference_period_error)
            current_epoch_error = _wd_epoch_residual_minutes(float(row.current_t0_bjd))
            reference_epoch_error = _wd_epoch_residual_minutes(float(row.reference_t0_bjd))
            epoch_error = max(current_epoch_error, reference_epoch_error)
            current_depth = float(row.current_depth)
            reference_depth = float(row.reference_depth)
            depth_ratio = reference_depth / current_depth if current_depth > 0 else np.nan
            timing_status = _status_max(
                epoch_error,
                config.wd1856_pass_epoch_residual_min,
                config.wd1856_review_epoch_residual_min,
            )
            period_status = _status_max(
                period_error,
                config.wd1856_pass_period_relative_error,
                config.wd1856_review_period_relative_error,
            )
            physical_depth_status = (
                "pass"
                if 0 < current_depth < 1 and 0 < reference_depth < 1
                else "fail"
            )
            wd_status = _worst_status(
                [timing_status, period_status, physical_depth_status]
            )
            aperture_statuses.append(wd_status)
            aperture_summaries.append(
                {
                    "aperture": str(row.aperture),
                    "status": wd_status,
                    "current_period_relative_error": current_period_error,
                    "reference_period_relative_error": reference_period_error,
                    "current_epoch_residual_min": current_epoch_error,
                    "reference_epoch_residual_min": reference_epoch_error,
                    "current_depth": current_depth,
                    "reference_depth": reference_depth,
                    "reference_to_current_depth_ratio": depth_ratio,
                    "depth_ratio_is_diagnostic_only": True,
                }
            )
        statuses["wd1856"] = _worst_status(aperture_statuses)
        wd_summary = {
            "status": statuses["wd1856"],
            "apertures": aperture_summaries,
        }
    status = _worst_status(statuses.values())
    return {
        "status": status,
        "reasons": ["independent extraction comparison failed"] if status == "fail" else [],
        "n_targets": n_targets,
        "n_detectors": n_detectors,
        "common_cadence_fraction_median": common_median,
        "common_cadence_fraction_min": common_min,
        "current_cadence_coverage_fraction_min": current_common_min,
        "reference_cadence_coverage_fraction_min": reference_common_min,
        "abs_log_scatter_ratio_median": scatter_median,
        "abs_log_scatter_ratio_max": scatter_max,
        "scatter_ratio_is_diagnostic_only": True,
        "comparison_note": (
            "Tier-1 gates on independent signal presence and ephemeris timing. "
            "Depth and scatter ratios are diagnostic because aperture dilution "
            "and decontamination differ between extraction families."
        ),
        "component_status": statuses,
        "wd1856": wd_summary,
        "provenance": dict(manifest),
    }


def evaluate_tier0_prerequisite(
    summary: Mapping[str, Any],
    *,
    sector: int,
    config: Tier1QAConfig,
    tier0_summary_sha256: str,
    current_compact_sha256: str | None = None,
) -> dict[str, Any]:
    current_compact_sha256 = (
        config.expected_compact_sha256
        if current_compact_sha256 is None
        else str(current_compact_sha256)
    )
    reasons: list[str] = []
    if str(tier0_summary_sha256) != config.expected_tier0_summary_sha256:
        reasons.append("Tier-0 summary is not pinned by the locked configuration")
    if int(summary.get("sector", -1)) != int(sector):
        reasons.append("Tier-0 sector mismatch")
    if summary.get("passed") is not True:
        reasons.append("Tier-0 did not pass")
    if str(summary.get("contract_version", "")) != config.expected_tier0_contract:
        reasons.append("Tier-0 contract mismatch")
    if str(summary.get("qa_tier", "")) != "tier0_integrity_and_benchmark":
        reasons.append("Tier-0 report lacks the current tier label")
    if summary.get("science_ready") is not False:
        reasons.append("Tier-0 science_ready must be explicitly false")

    gates = summary.get("gates")
    if not isinstance(gates, Mapping):
        reasons.append("Tier-0 report is missing its gate mapping")
        gates = {}
    observed_gate_names = set(gates)
    expected_gate_names = set(TIER0_REQUIRED_GATES)
    if observed_gate_names != expected_gate_names:
        reasons.append("Tier-0 gate inventory does not match the current contract")
    for gate_name in TIER0_REQUIRED_GATES:
        gate = gates.get(gate_name)
        if not isinstance(gate, Mapping):
            reasons.append(f"Tier-0 gate {gate_name} is missing or malformed")
        elif gate.get("passed") is not True:
            reasons.append(f"Tier-0 gate {gate_name} did not pass")

    benchmarks = summary.get("benchmarks")
    wd1856 = benchmarks.get("wd1856") if isinstance(benchmarks, Mapping) else None
    if not isinstance(wd1856, Mapping) or wd1856.get("passed") is not True:
        reasons.append("Tier-0 WD 1856 benchmark did not pass")
    provenance = summary.get("provenance", {})
    if not isinstance(provenance, Mapping):
        provenance = {}
    tier0_compact_sha256 = str(provenance.get("compact_lc_sha256", ""))
    if (
        current_compact_sha256 != config.expected_compact_sha256
        or tier0_compact_sha256 != current_compact_sha256
    ):
        reasons.append("Tier-0 report is not bound to the locked compact product")
    if tuple(provenance.get("detrended_apertures", ())) != tuple(config.apertures):
        reasons.append("Tier-0 aperture scope mismatch")
    schema_summary_sha256 = str(provenance.get("schema_summary_sha256", ""))
    if not re.fullmatch(r"[0-9a-f]{64}", schema_summary_sha256):
        reasons.append("Tier-0 schema evidence hash is invalid")
    tier0_bls_peaks_sha256 = str(provenance.get("bls_peaks_sha256", ""))
    if tier0_bls_peaks_sha256 != config.expected_tier0_bls_peaks_sha256:
        reasons.append("Tier-0 BLS peak table is not pinned by the locked configuration")
    return {
        "status": "fail" if reasons else "pass",
        "reasons": reasons,
        "contract_version": summary.get("contract_version"),
        "sector": summary.get("sector"),
        "tier0_summary_sha256": str(tier0_summary_sha256),
        "compact_lc_sha256": tier0_compact_sha256,
        "bls_peaks_sha256": tier0_bls_peaks_sha256,
        "required_gates": list(TIER0_REQUIRED_GATES),
    }


def evaluate_tier1_gates(
    *,
    sector: int,
    config: Tier1QAConfig,
    tier0_summary: Mapping[str, Any],
    tier0_summary_sha256: str,
    target_metrics: pd.DataFrame,
    aperture_metrics: pd.DataFrame,
    pair_metrics: pd.DataFrame,
    injection_metrics: pd.DataFrame,
    injection_manifest: Mapping[str, Any],
    independent_metrics: pd.DataFrame,
    independent_manifest: Mapping[str, Any],
    current_compact_sha256: str | None = None,
    catalog_detectors: Mapping[int, str] | None = None,
    cadence_reference_gate: Mapping[str, Any] | None = None,
    injection_source_parity_gate: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    current_compact_sha256 = (
        config.expected_compact_sha256
        if current_compact_sha256 is None
        else str(current_compact_sha256)
    )
    scatter_gate, magnitude_bins = evaluate_population_scatter(aperture_metrics, config)
    cadence_gate, targets_with_loss = evaluate_cadence_quality(
        target_metrics, aperture_metrics, config
    )
    gates = {
        "cadence_reference_prerequisite": dict(cadence_reference_gate)
        if cadence_reference_gate is not None
        else {
            "status": "fail",
            "reasons": ["authoritative cadence-reference manifest was not supplied"],
        },
        "injection_source_parity_prerequisite": dict(
            injection_source_parity_gate
        )
        if injection_source_parity_gate is not None
        else {
            "status": "fail",
            "reasons": ["injection-source parity evidence was not supplied"],
        },
        "tier0_prerequisite": evaluate_tier0_prerequisite(
            tier0_summary,
            sector=sector,
            config=config,
            tier0_summary_sha256=tier0_summary_sha256,
            current_compact_sha256=current_compact_sha256,
        ),
        "population_scatter": scatter_gate,
        "cadence_and_finite_data": cadence_gate,
        "aperture_outliers": evaluate_aperture_outliers(pair_metrics, config),
        "fixed_injection_preservation": evaluate_fixed_injections(
            injection_metrics,
            injection_manifest,
            sector=sector,
            config=config,
        ),
        "independent_extraction": evaluate_independent_extraction(
            independent_metrics,
            independent_manifest,
            config,
            sector=sector,
            catalog_detectors=catalog_detectors,
            current_compact_sha256=current_compact_sha256,
        ),
    }
    targets_with_loss = attach_target_qa_flags(
        targets_with_loss, aperture_metrics, pair_metrics, config
    )
    overall_status = _worst_status(gate["status"] for gate in gates.values())
    passed = overall_status == "pass"
    enrichment_ready = bool(passed and config.scope == "active_search_pair")
    science_ready = False
    return {
        "status": overall_status,
        "passed": passed,
        "enrichment_ready": enrichment_ready,
        "science_ready": science_ready,
        "target_qa": {
            "n_pass": int(targets_with_loss["tier1_target_qa_pass"].sum()),
            "n_review": int(
                targets_with_loss["tier1_target_qa_status"].eq("review").sum()
            ),
            "n_fail": int(targets_with_loss["tier1_target_qa_status"].eq("fail").sum()),
            "candidate_teacher_filter": "tier1_target_qa_pass == True",
        },
        "gates": gates,
    }, magnitude_bins, targets_with_loss


def build_target_eligibility(
    targets: pd.DataFrame, config: Tier1QAConfig
) -> pd.DataFrame:
    """Build the compact sector-observation mask consumed downstream."""

    required = {
        "sector",
        "tic",
        "camera",
        "ccd",
        "detector",
        "tmag",
        "tier1_cadence_qa_status",
        "tier1_scatter_qa_status",
        "tier1_aperture_pair_qa_status",
        "tier1_target_qa_status",
        "tier1_target_qa_reasons",
        "tier1_target_qa_pass",
    }
    missing = sorted(required - set(targets.columns))
    if missing:
        raise KeyError(f"target eligibility inputs are missing columns: {missing}")
    if targets.duplicated(["sector", "tic"]).any():
        raise ValueError("target eligibility requires unique (sector, TIC) rows")
    observed_sectors = sorted(
        set(pd.to_numeric(targets["sector"], errors="coerce").dropna().astype(int))
    )
    if observed_sectors != [config.sector]:
        raise ValueError(f"target eligibility sector mismatch: {observed_sectors}")
    reason_values = {
        reason
        for value in targets["tier1_target_qa_reasons"].fillna("").astype(str)
        for reason in value.split(";")
        if reason
    }
    unknown = sorted(reason_values - set(TIER1_TARGET_REASON_CODES))
    if unknown:
        raise ValueError(f"target eligibility has unknown reason codes: {unknown}")
    status = targets["tier1_target_qa_status"].astype("string")
    if status.isna().any() or not status.isin(STATUS_ORDER).all():
        raise ValueError("target eligibility has invalid target QA status values")
    pass_values = targets["tier1_target_qa_pass"]
    if not pd.api.types.is_bool_dtype(pass_values.dtype) or pass_values.isna().any():
        raise ValueError("target eligibility pass flags must be non-null booleans")
    expected_pass = status.eq("pass")
    if not pass_values.eq(expected_pass).all():
        raise ValueError("target eligibility status and pass flag are inconsistent")

    columns = [
        "sector",
        "tic",
        "gaia_dr3_source_id",
        "camera",
        "ccd",
        "detector",
        "tmag",
        "tier1_cadence_qa_status",
        "tier1_scatter_qa_status",
        "tier1_aperture_pair_qa_status",
        "tier1_target_qa_status",
        "tier1_target_qa_reasons",
        "tier1_target_qa_pass",
    ]
    work = targets.copy()
    if "gaia_dr3_source_id" not in work:
        work["gaia_dr3_source_id"] = pd.NA
    output = work.loc[:, columns].copy()
    output.insert(0, "tier1_scope", config.scope)
    output.insert(0, "tier1_config_name", config.name)
    output.insert(0, "tier1_contract_version", config.contract_version)
    output.insert(
        3,
        "sector_tic_key",
        [
            f"s{int(sector):04d}-tic{int(tic):016d}"
            for sector, tic in zip(output["sector"], output["tic"], strict=True)
        ],
    )
    output["tier1_target_qa_pass"] = expected_pass.astype(bool).to_numpy()
    return output.sort_values(["sector", "tic"], kind="stable").reset_index(drop=True)


def summarize_target_eligibility(eligibility: pd.DataFrame) -> pd.DataFrame:
    """Summarize final target eligibility by detector for QA and plotting."""

    required = {
        "camera",
        "ccd",
        "detector",
        "tier1_target_qa_status",
        "tier1_target_qa_pass",
    }
    missing = sorted(required - set(eligibility.columns))
    if missing:
        raise KeyError(f"detector eligibility inputs are missing columns: {missing}")
    rows: list[dict[str, Any]] = []
    for (camera, ccd, detector), group in eligibility.groupby(
        ["camera", "ccd", "detector"], sort=True
    ):
        status = group["tier1_target_qa_status"].astype(str)
        n_targets = int(len(group))
        n_pass = int(status.eq("pass").sum())
        rows.append(
            {
                "camera": int(camera),
                "ccd": int(ccd),
                "detector": str(detector),
                "n_targets": n_targets,
                "n_pass": n_pass,
                "n_review": int(status.eq("review").sum()),
                "n_fail": int(status.eq("fail").sum()),
                "pass_fraction": float(n_pass / n_targets) if n_targets else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["camera", "ccd"], kind="stable").reset_index(
        drop=True
    )


def _save_figure_atomic(figure: Any, path: Path, *, dpi: int | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        figure.savefig(
            temporary,
            dpi=dpi,
            bbox_inches="tight",
            format=path.suffix.lstrip("."),
        )
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def plot_tier1_diagnostics(
    *,
    targets: pd.DataFrame,
    aperture_metrics: pd.DataFrame,
    pair_metrics: pd.DataFrame,
    injection_metrics: pd.DataFrame,
    config: Tier1QAConfig,
    output_png: Path,
    output_pdf: Path,
) -> None:
    """Render a compact four-panel diagnostic for the bounded Tier-1 gate."""

    import matplotlib.pyplot as plt

    from twirl.plotting.style import apply_twirl_style, get_ordered_palette

    template = apply_twirl_style("full_page")
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(template["figsize"][0], 5.4),
        constrained_layout=True,
    )
    colors = get_ordered_palette(len(config.apertures), "viridis")
    aperture_labels = {
        ADP_ONLY_APERTURES[0]: "ADP small",
        ADP_ONLY_APERTURES[1]: "ADP primary",
    }
    finite_tmag = pd.to_numeric(aperture_metrics["tmag"], errors="coerce")
    finite_tmag = finite_tmag[np.isfinite(finite_tmag)]
    x_min = float(finite_tmag.min()) if len(finite_tmag) else 15.0
    x_max = float(finite_tmag.max()) if len(finite_tmag) else 21.0
    envelope_tmag = np.linspace(x_min, x_max, 256)
    envelope_scale = 10.0 ** (
        config.scatter_envelope_dex_per_mag
        * (envelope_tmag - config.scatter_envelope_reference_tmag)
    )

    for aperture, color in zip(config.apertures, colors, strict=True):
        rows = aperture_metrics.loc[
            aperture_metrics["aperture"].astype(str).eq(aperture)
        ]
        axes[0, 0].scatter(
            rows["tmag"],
            rows["mad_ppm"],
            s=2.0,
            alpha=0.22,
            linewidths=0,
            rasterized=True,
            color=color,
            label=aperture_labels.get(aperture, aperture),
        )
        axes[0, 1].scatter(
            rows["tmag"],
            rows["rms5_ppm"],
            s=2.0,
            alpha=0.22,
            linewidths=0,
            rasterized=True,
            color=color,
            label=aperture_labels.get(aperture, aperture),
        )
    axes[0, 0].plot(
        envelope_tmag,
        config.scatter_pass_mad_ppm_at_reference * envelope_scale,
        color="0.15",
        linewidth=0.9,
        label="pass envelope",
    )
    axes[0, 0].plot(
        envelope_tmag,
        config.scatter_review_mad_ppm_at_reference * envelope_scale,
        color="0.25",
        linewidth=0.9,
        linestyle="--",
        label="review envelope",
    )
    axes[0, 1].plot(
        envelope_tmag,
        config.scatter_pass_rms5_ppm_at_reference * envelope_scale,
        color="0.15",
        linewidth=0.9,
    )
    axes[0, 1].plot(
        envelope_tmag,
        config.scatter_review_rms5_ppm_at_reference * envelope_scale,
        color="0.25",
        linewidth=0.9,
        linestyle="--",
    )
    for axis, ylabel in zip(
        axes[0], ("Robust MAD (ppm)", "5σ-clipped RMS (ppm)"), strict=True
    ):
        axis.set_yscale("log")
        axis.set_xlabel("TESS magnitude")
        axis.set_ylabel(ylabel)
    axes[0, 0].legend(loc="best", frameon=True, ncol=2)

    ratio = pd.to_numeric(pair_metrics["mad_ratio"], errors="coerce")
    correlation = pd.to_numeric(pair_metrics["correlation"], errors="coerce")
    pair_tmag = pd.to_numeric(pair_metrics["tmag"], errors="coerce")
    pair_good = (
        np.isfinite(ratio)
        & (ratio > 0)
        & np.isfinite(correlation)
        & np.isfinite(pair_tmag)
    )
    if pair_good.any():
        pair_scatter = axes[1, 0].scatter(
            ratio[pair_good],
            correlation[pair_good],
            c=pair_tmag[pair_good],
            cmap="viridis",
            s=2.5,
            alpha=0.3,
            linewidths=0,
            rasterized=True,
        )
        pair_colorbar = figure.colorbar(
            pair_scatter,
            ax=axes[1, 0],
            pad=0.02,
            fraction=0.05,
        )
        pair_colorbar.set_label("TESS magnitude")
    axes[1, 0].set_xscale("log")
    for boundary in (
        config.aperture_moderate_ratio_low,
        config.aperture_moderate_ratio_high,
    ):
        axes[1, 0].axvline(boundary, color="0.35", linewidth=0.8, linestyle="--")
    axes[1, 0].axhline(
        config.aperture_pass_correlation_median, color="0.2", linewidth=0.8
    )
    axes[1, 0].axhline(
        config.aperture_negative_correlation_threshold,
        color="0.35",
        linewidth=0.8,
        linestyle="--",
    )
    axes[1, 0].set_xlabel("Primary/small MAD ratio")
    axes[1, 0].set_ylabel("Aperture correlation")

    injection = injection_metrics.copy()
    if "status" in injection:
        injection = injection.loc[injection["status"].astype(str).eq("ok")]
    retention_underflow = 0
    retention_overflow = 0
    for aperture, color in zip(config.apertures, colors, strict=True):
        values = pd.to_numeric(
            injection.loc[
                injection["aperture"].astype(str).eq(aperture),
                "depth_retention_fraction",
            ],
            errors="coerce",
        )
        values = values[np.isfinite(values)]
        if len(values):
            retention_underflow += int((values < 0.0).sum())
            retention_overflow += int((values > 2.0).sum())
            axes[1, 1].hist(
                np.clip(values, 0.0, 2.0),
                bins=np.linspace(0.0, 2.0, 51),
                histtype="step",
                linewidth=1.0,
                color=color,
                label=aperture_labels.get(aperture, aperture),
            )
    axes[1, 1].axvline(1.0, color="0.15", linewidth=0.9)
    for boundary in (config.injection_inband_low, config.injection_inband_high):
        axes[1, 1].axvline(boundary, color="0.35", linewidth=0.8, linestyle="--")
    axes[1, 1].set_xlim(0.0, 2.0)
    axes[1, 1].set_xlabel("Injected-depth retention")
    axes[1, 1].set_ylabel("Injection-aperture rows")
    axes[1, 1].legend(loc="best", frameon=True)

    status_counts = targets["tier1_target_qa_status"].astype(str).value_counts()
    axes[1, 1].text(
        0.98,
        0.96,
        "Targets: "
        + ", ".join(
            f"{name}={int(status_counts.get(name, 0))}"
            for name in ("pass", "review", "fail")
        ),
        transform=axes[1, 1].transAxes,
        ha="right",
        va="top",
        fontsize=template["annotation_size"],
    )
    axes[1, 1].text(
        0.98,
        0.86,
        f"Clipped to edges: <0={retention_underflow}, >2={retention_overflow}",
        transform=axes[1, 1].transAxes,
        ha="right",
        va="top",
        fontsize=template["annotation_size"],
    )
    _save_figure_atomic(figure, output_png, dpi=220)
    _save_figure_atomic(figure, output_pdf)
    plt.close(figure)


def plot_detector_eligibility(
    detector_summary: pd.DataFrame, *, output_png: Path, output_pdf: Path
) -> None:
    """Render the 4x4 target-pass fraction map for one TESS sector."""

    import matplotlib.pyplot as plt

    from twirl.plotting.style import apply_twirl_style

    template = apply_twirl_style("column")
    values = np.full((4, 4), np.nan, dtype=float)
    labels: dict[tuple[int, int], str] = {}
    for row in detector_summary.itertuples(index=False):
        camera = int(row.camera)
        ccd = int(row.ccd)
        values[camera - 1, ccd - 1] = float(row.pass_fraction)
        labels[(camera - 1, ccd - 1)] = f"{int(row.n_pass)}/{int(row.n_targets)}"
    colormap = plt.get_cmap("viridis").copy()
    colormap.set_bad("0.9")
    figure, axis = plt.subplots(figsize=(template["figsize"][0], 3.0))
    image = axis.imshow(values, vmin=0.0, vmax=1.0, cmap=colormap, aspect="equal")
    for (row, column), label in labels.items():
        value = values[row, column]
        text_color = "white" if np.isfinite(value) and value < 0.55 else "black"
        axis.text(column, row, label, ha="center", va="center", color=text_color)
    axis.set_xticks(np.arange(4), labels=["1", "2", "3", "4"])
    axis.set_yticks(np.arange(4), labels=["1", "2", "3", "4"])
    axis.set_xlabel("CCD")
    axis.set_ylabel("Camera")
    colorbar = figure.colorbar(image, ax=axis, pad=0.03, fraction=0.05)
    colorbar.set_label("Tier-1 target pass fraction")
    _save_figure_atomic(figure, output_png, dpi=220)
    _save_figure_atomic(figure, output_pdf)
    plt.close(figure)


def _read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() in {".csv", ".gz"} or path.name.endswith(".csv.gz"):
        # The frozen injection metadata digest is defined on binary float64
        # values.  Pandas' round-trip parser recovers the exact values emitted
        # by ``to_csv``; the default fast parser may move a final bit and make
        # an otherwise identical published evidence table fail its hash lock.
        return pd.read_csv(path, low_memory=False, float_precision="round_trip")
    raise ValueError(f"unsupported table format: {path}")


def evaluate_cadence_reference_evidence(
    table: pd.DataFrame,
    table_path: Path,
    manifest: Mapping[str, Any],
    manifest_path: Path | None = None,
    *,
    sector: int,
    config: Tier1QAConfig,
) -> dict[str, Any]:
    """Validate the provenance envelope around an authoritative cadence table."""

    required = {
        "contract_version",
        "builder_version",
        "sector",
        "cadence_authority",
        "quality_authority",
        "quality_composition",
        "table_sha256",
        "table_columns",
        "n_rows",
        "detectors",
        "orbits",
        "source_file_sha256",
        "sources",
        "n_spoc_authority_files_verified",
        "n_qlp_qflag_files_verified",
        "n_nonzero_spoc_quality",
        "n_nonzero_qlp_quality",
        "n_nonzero_external_quality",
        "n_spoc_rows_excluded_by_quat",
        "authority_exclusions",
        "authority_exclusions_sha256",
    }
    reasons: list[str] = []
    missing = sorted(required - set(manifest))
    if missing:
        reasons.append(f"cadence-reference manifest missing fields: {missing}")
        return {"status": "fail", "reasons": reasons}
    expected_contract = config.cadence_reference_contract_template.format(sector=sector)
    if str(manifest["contract_version"]) != expected_contract:
        reasons.append("cadence-reference contract mismatch")
    if str(manifest["builder_version"]) != CADENCE_REFERENCE_BUILDER_VERSION:
        reasons.append("cadence-reference builder version mismatch")
    try:
        manifest_sector = int(manifest["sector"])
        n_rows = int(manifest["n_rows"])
    except (TypeError, ValueError):
        manifest_sector = -1
        n_rows = -1
    if manifest_sector != sector:
        reasons.append("cadence-reference sector mismatch")
    observed_hash = file_sha256(table_path)
    if str(manifest["table_sha256"]) != observed_hash:
        reasons.append("cadence-reference table hash mismatch")
    if observed_hash != config.expected_cadence_reference_sha256:
        reasons.append("cadence-reference table is not pinned by the locked config")
    manifest_sha256 = file_sha256(manifest_path) if manifest_path is not None else ""
    if manifest_sha256 != config.expected_cadence_reference_manifest_sha256:
        reasons.append("cadence-reference manifest is not pinned by the locked config")
    if n_rows != len(table):
        reasons.append("cadence-reference row-count mismatch")
    if str(manifest["cadence_authority"]) != "qlp_cam_quat":
        reasons.append("cadence reference is not sourced from QLP camera quaternion tables")
    if str(manifest["quality_authority"]) != "spoc_and_qlp_quality_flags":
        reasons.append("quality reference is not sourced from SPOC and QLP flags")
    expected_columns = (
        "sector",
        "orbitid",
        "camera",
        "ccd",
        "cadenceno",
        "spoc_quality",
        "qlp_quality",
        "external_quality",
    )
    expected_quality_composition = {
        "external_quality": "spoc_quality | (qlp_quality << 30)",
        "qlp_quality_raw_values": [0, 1],
        "qlp_quality_external_bit": 30,
    }
    if manifest["quality_composition"] != expected_quality_composition:
        reasons.append("cadence-reference quality composition mismatch")
    if tuple(manifest["table_columns"]) != expected_columns:
        reasons.append("cadence-reference manifest table columns mismatch")
    if tuple(table.columns) != expected_columns:
        reasons.append(
            "cadence-reference table must contain the exact ordered v1 schema"
        )
        observed_detectors: list[str] = []
        observed_orbits: list[int] = []
    else:
        observed_detectors = sorted(
            {
                f"cam{int(camera)}_ccd{int(ccd)}"
                for camera, ccd in zip(table["camera"], table["ccd"], strict=True)
            }
        )
        observed_orbits = sorted(
            set(pd.to_numeric(table["orbitid"], errors="coerce").dropna().astype(int))
        )
        spoc_quality = pd.to_numeric(table["spoc_quality"], errors="coerce")
        qlp_quality = pd.to_numeric(table["qlp_quality"], errors="coerce")
        external_quality = pd.to_numeric(table["external_quality"], errors="coerce")
        spoc_values = spoc_quality.to_numpy(dtype=np.int64, na_value=-1)
        qlp_values = qlp_quality.to_numpy(dtype=np.int64, na_value=-1)
        external_values = external_quality.to_numpy(dtype=np.int64, na_value=-1)
        valid_quality = bool(
            np.isfinite(spoc_quality).all()
            and np.isfinite(qlp_quality).all()
            and np.isfinite(external_quality).all()
            and (spoc_quality >= 0).all()
            and qlp_quality.isin((0, 1)).all()
            and np.array_equal(
                external_values, spoc_values | (qlp_values << 30)
            )
        )
        if not valid_quality:
            reasons.append("cadence-reference external quality derivation is invalid")
        count_fields = {
            "n_nonzero_spoc_quality": "spoc_quality",
            "n_nonzero_qlp_quality": "qlp_quality",
            "n_nonzero_external_quality": "external_quality",
        }
        for manifest_field, column in count_fields.items():
            observed_count = int(
                np.count_nonzero(table[column].to_numpy(dtype=np.int64))
            )
            try:
                declared_count = int(manifest[manifest_field])
            except (TypeError, ValueError):
                declared_count = -1
            if declared_count != observed_count:
                reasons.append(
                    f"cadence-reference {manifest_field} disagrees with the table"
                )
    if sorted(str(value) for value in manifest["detectors"]) != observed_detectors:
        reasons.append("cadence-reference detector list mismatch")
    try:
        manifest_orbits = sorted(int(value) for value in manifest["orbits"])
    except (TypeError, ValueError):
        manifest_orbits = []
    if manifest_orbits != observed_orbits:
        reasons.append("cadence-reference orbit list mismatch")
    if int(sector) == 56:
        expected_detectors = [
            f"cam{camera}_ccd{ccd}"
            for camera in range(1, 5)
            for ccd in range(1, 5)
        ]
        if observed_detectors != expected_detectors:
            reasons.append("S56 cadence reference does not cover all 16 detectors")
        if observed_orbits != [119, 120]:
            reasons.append("S56 cadence reference does not cover exactly orbits 119/120")
    source_hashes = manifest["source_file_sha256"]
    if (
        not isinstance(source_hashes, Mapping)
        or not source_hashes
        or any(
            len(str(value)) != 64
            or any(char not in "0123456789abcdef" for char in str(value).lower())
            for value in source_hashes.values()
        )
    ):
        reasons.append("cadence-reference source-file hashes are missing or invalid")
    sources = manifest["sources"]
    if not isinstance(sources, list) or not sources:
        reasons.append("cadence-reference source inventory is missing or invalid")
    elif isinstance(source_hashes, Mapping):
        source_inventory: dict[str, str] = {}
        for source in sources:
            if not isinstance(source, Mapping):
                reasons.append("cadence-reference source inventory has a non-object")
                continue
            path = str(source.get("path", "")).strip()
            digest = str(source.get("sha256", "")).lower()
            if not path or path in source_inventory:
                reasons.append(
                    "cadence-reference source inventory has an empty/duplicate path"
                )
                continue
            source_inventory[path] = digest
        declared_sources = {
            str(path): str(digest).lower()
            for path, digest in source_hashes.items()
        }
        if source_inventory != declared_sources:
            reasons.append(
                "cadence-reference source inventory disagrees with source hashes"
            )
    try:
        n_spoc_sources = int(manifest["n_spoc_authority_files_verified"])
        n_qlp_sources = int(manifest["n_qlp_qflag_files_verified"])
    except (TypeError, ValueError):
        n_spoc_sources, n_qlp_sources = -1, -1
    if sector == 56 and (n_spoc_sources != 16 or n_qlp_sources != 32):
        reasons.append("S56 cadence reference lacks all 16 SPOC and 32 QLP inputs")
    authority_exclusions = manifest["authority_exclusions"]
    try:
        n_authority_exclusions = int(authority_exclusions["n_rows"])
        n_spoc_excluded = int(manifest["n_spoc_rows_excluded_by_quat"])
    except (KeyError, TypeError, ValueError):
        n_authority_exclusions = n_spoc_excluded = -1
        reasons.append("cadence-reference authority exclusions are malformed")
    if n_authority_exclusions != n_spoc_excluded or n_authority_exclusions < 0:
        reasons.append("cadence-reference authority-exclusion count mismatch")
    try:
        observed_exclusions_sha256 = authority_exclusions_sha256(
            authority_exclusions
        )
    except (TypeError, ValueError):
        observed_exclusions_sha256 = ""
    if str(manifest["authority_exclusions_sha256"]) != (
        observed_exclusions_sha256
    ):
        reasons.append("cadence-reference authority-exclusion hash mismatch")
    return {
        "status": "fail" if reasons else "pass",
        "reasons": reasons,
        "contract_version": manifest.get("contract_version"),
        "table_sha256": observed_hash,
        "manifest_sha256": manifest_sha256,
        "n_rows": int(len(table)),
        "detectors": observed_detectors,
        "orbits": observed_orbits,
        "cadence_authority": manifest.get("cadence_authority"),
        "quality_authority": manifest.get("quality_authority"),
        "authority_exclusions_sha256": observed_exclusions_sha256,
        "n_authority_exclusions": n_authority_exclusions,
    }


def evaluate_injection_source_parity(
    report: Mapping[str, Any],
    report_path: Path,
    injection_manifest: Mapping[str, Any],
    *,
    current_compact_sha256: str,
    config: Tier1QAConfig,
) -> dict[str, Any]:
    """Bind the frozen injection source to the exact compact through audited parity."""

    reasons: list[str] = []
    report_sha256 = file_sha256(report_path)
    if report_sha256 != config.expected_injection_parity_report_sha256:
        reasons.append("injection-source parity report hash mismatch")
    if report.get("passed") is not True or report.get("failures") not in ([], None):
        reasons.append("injection-source input verification did not pass")
    hashes = report.get("hashes", {})
    parity = report.get("compact_rebuild_parity", {})
    if not isinstance(hashes, Mapping) or not isinstance(parity, Mapping):
        reasons.append("injection-source parity report structure is invalid")
        hashes = {}
        parity = {}
    if str(hashes.get("s56_A2v1_adp_pair.h5", "")) != current_compact_sha256:
        reasons.append("parity report does not bind the current compact hash")
    if current_compact_sha256 != config.expected_compact_sha256:
        reasons.append("current compact hash is not the locked Tier-1 product")
    if (
        str(hashes.get("s56_A2v1_adp_pair_rebuilt.h5", ""))
        != config.expected_injection_source_adp_sha256
    ):
        reasons.append("parity report does not bind the frozen injection ADP source")
    try:
        parity_counts_ok = bool(
            parity.get("passed") is True
            and int(parity.get("n_mismatched_targets", -1)) == 0
            and int(parity.get("n_missing_active_targets", -1)) == 0
            and int(parity.get("n_extra_active_targets", -1)) == 0
            and int(parity.get("n_dataset_comparisons", 0)) >= 1
        )
    except (TypeError, ValueError):
        parity_counts_ok = False
    if not parity_counts_ok:
        reasons.append("compact rebuild parity comparison is incomplete or failed")
    source_paths = injection_manifest.get("source_adp_h5", [])
    if not isinstance(source_paths, Sequence) or isinstance(source_paths, (str, bytes)):
        source_paths = []
    source_paths = [str(path) for path in source_paths if str(path)]
    if len(source_paths) != 1:
        reasons.append("injection shards do not name one common ADP source")
    else:
        active_path = str(parity.get("active_h5", ""))
        if Path(source_paths[0]).name != Path(active_path).name:
            reasons.append("injection shard ADP source disagrees with the parity report")
    if Path(str(parity.get("reference_h5", ""))).name != "s56_A2v1_adp_pair.h5":
        reasons.append("parity reference is not the locked S56 compact product")
    return {
        "status": "fail" if reasons else "pass",
        "reasons": reasons,
        "report_sha256": report_sha256,
        "current_compact_sha256": current_compact_sha256,
        "injection_source_adp_sha256": hashes.get("s56_A2v1_adp_pair_rebuilt.h5"),
        "n_dataset_comparisons": parity.get("n_dataset_comparisons"),
        "source_adp_h5": source_paths,
    }


def _sha256_text_mapping(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(_safe_json(payload), sort_keys=True, allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def _write_table_atomic(frame: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        if path.suffix.lower() == ".parquet":
            frame.to_parquet(temporary, compression="zstd", index=False)
        else:
            frame.to_csv(temporary, index=False)
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def run_a2v1_tier1_qa(
    *,
    sector: int,
    config_path: Path,
    tier0_summary_path: Path,
    compact_lc: Path,
    cadence_reference_path: Path,
    cadence_reference_manifest_path: Path,
    injection_source_parity_path: Path,
    independent_metrics_path: Path,
    independent_manifest_path: Path,
    out_dir: Path,
    gate_json: Path,
    injection_metrics_path: Path | None = None,
    injection_manifest_path: Path | None = None,
    injection_shards: Sequence[Path] = (),
) -> dict[str, Any]:
    """Run the scoped Tier-1 gate and write compact, auditable outputs."""

    config_path = Path(config_path)
    tier0_summary_path = Path(tier0_summary_path)
    compact_lc = Path(compact_lc)
    cadence_reference_path = Path(cadence_reference_path)
    cadence_reference_manifest_path = Path(cadence_reference_manifest_path)
    injection_source_parity_path = Path(injection_source_parity_path)
    independent_metrics_path = Path(independent_metrics_path)
    independent_manifest_path = Path(independent_manifest_path)
    out_dir = Path(out_dir)
    gate_json = Path(gate_json)
    injection_metrics_path = (
        Path(injection_metrics_path) if injection_metrics_path is not None else None
    )
    injection_manifest_path = (
        Path(injection_manifest_path) if injection_manifest_path is not None else None
    )
    injection_shards = tuple(Path(path) for path in injection_shards)

    paths = {
        "target_metrics": out_dir / "target_metrics.parquet",
        "target_eligibility": out_dir / "target_eligibility.csv",
        "aperture_metrics": out_dir / "aperture_metrics.parquet",
        "aperture_pair_metrics": out_dir / "aperture_pair_metrics.csv",
        "detector_summary": out_dir / "detector_summary.csv",
        "magnitude_bins": out_dir / "magnitude_bins.csv",
        "injection_metrics": out_dir / "fixed_injection_metrics.csv",
        "injection_manifest": out_dir / "fixed_injection_manifest.json",
        "independent_metrics": out_dir / "independent_extraction_metrics.csv",
        "qa_diagnostics_png": out_dir / "tier1_qa_diagnostics.png",
        "qa_diagnostics_pdf": out_dir / "tier1_qa_diagnostics.pdf",
        "detector_eligibility_png": out_dir / "tier1_detector_eligibility.png",
        "detector_eligibility_pdf": out_dir / "tier1_detector_eligibility.pdf",
        "summary": out_dir / "summary.json",
    }
    input_paths = {
        config_path,
        tier0_summary_path,
        compact_lc,
        cadence_reference_path,
        cadence_reference_manifest_path,
        injection_source_parity_path,
        independent_metrics_path,
        independent_manifest_path,
        *injection_shards,
    }
    if injection_metrics_path is not None:
        input_paths.add(injection_metrics_path)
    if injection_manifest_path is not None:
        input_paths.add(injection_manifest_path)
    input_resolved = {path.resolve() for path in input_paths}
    output_path_list = [*paths.values(), gate_json]
    output_resolved_list = [path.resolve() for path in output_path_list]
    output_resolved = set(output_resolved_list)
    if len(output_resolved) != len(output_resolved_list):
        raise ValueError("Tier-1 output paths collide with one another")
    temporary_resolved_list = [
        path.with_suffix(path.suffix + suffix).resolve()
        for path in output_path_list
        for suffix in (".tmp", ".publish.tmp")
    ]
    temporary_resolved = set(temporary_resolved_list)
    if (
        len(temporary_resolved) != len(temporary_resolved_list)
        or temporary_resolved & output_resolved
        or temporary_resolved & input_resolved
    ):
        raise ValueError("Tier-1 temporary output paths collide with inputs or outputs")
    collisions = sorted(str(path) for path in input_resolved & output_resolved)
    if collisions:
        raise ValueError(f"Tier-1 inputs and outputs collide: {collisions}")

    # Hash every evidence input before any long scan. The same frozen mapping
    # is used for validation and published provenance; live re-hashes are never
    # substituted for the generation that was actually evaluated.
    guarded_input_paths = sorted(
        input_paths, key=lambda path: str(path.resolve())
    )
    initial_input_sha256 = {
        str(path.resolve()): file_sha256(path) for path in guarded_input_paths
    }

    def initial_sha256(path: Path) -> str:
        return initial_input_sha256[str(Path(path).resolve())]

    def assert_inputs_unchanged() -> None:
        final_input_sha256 = {
            str(path.resolve()): file_sha256(path) for path in guarded_input_paths
        }
        if final_input_sha256 != initial_input_sha256:
            raise ValueError(
                "one or more Tier-1 evidence inputs changed during the audit"
            )

    config = load_tier1_config(config_path)
    if int(sector) != config.sector:
        raise ValueError(f"requested sector {sector} does not match locked config S{config.sector}")
    independent_metrics_sha256 = initial_input_sha256[
        str(independent_metrics_path.resolve())
    ]
    independent_manifest_sha256 = initial_input_sha256[
        str(independent_manifest_path.resolve())
    ]
    if independent_metrics_sha256 != config.expected_independent_metrics_sha256:
        raise ValueError(
            "independent metrics file is not pinned by the locked Tier-1 configuration"
        )
    if independent_manifest_sha256 != config.expected_independent_manifest_sha256:
        raise ValueError(
            "independent manifest file is not pinned by the locked Tier-1 configuration"
        )
    compact_sha256 = initial_input_sha256[str(compact_lc.resolve())]
    if compact_sha256 != config.expected_compact_sha256:
        raise ValueError(
            "compact product hash does not match the locked Tier-1 configuration: "
            f"{compact_sha256}"
        )
    tier0_summary_sha256 = initial_input_sha256[str(tier0_summary_path.resolve())]
    tier0_summary = json.loads(tier0_summary_path.read_text())
    tier0_preflight = evaluate_tier0_prerequisite(
        tier0_summary,
        sector=sector,
        config=config,
        tier0_summary_sha256=tier0_summary_sha256,
        current_compact_sha256=compact_sha256,
    )
    if tier0_preflight["status"] != "pass":
        raise ValueError(
            "Tier-0 preflight failed before population scan: "
            + "; ".join(tier0_preflight["reasons"])
        )
    catalog = compact_target_catalog(compact_lc)
    observed_sectors = sorted(set(catalog["sector"].astype(int)))
    if observed_sectors != [int(sector)]:
        raise ValueError(f"compact export sector mismatch: {observed_sectors}")
    catalog_detectors = {
        int(row.tic): f"cam{int(row.camera)}_ccd{int(row.ccd)}"
        for row in catalog.itertuples(index=False)
    }
    cadence_reference = _read_table(cadence_reference_path)
    cadence_reference_manifest = json.loads(
        cadence_reference_manifest_path.read_text()
    )
    cadence_reference_gate = evaluate_cadence_reference_evidence(
        cadence_reference,
        cadence_reference_path,
        cadence_reference_manifest,
        cadence_reference_manifest_path,
        sector=sector,
        config=config,
    )
    if cadence_reference_gate["status"] != "pass":
        raise ValueError(
            "cadence-reference preflight failed before population scan: "
            + "; ".join(cadence_reference_gate["reasons"])
        )
    # Use the central immutable reference loader as the final authority for all
    # downstream masks, including the precomputed-injection branch.
    quality_reference = load_external_quality_reference(
        table_path=cadence_reference_path,
        manifest_path=cadence_reference_manifest_path,
        sector=sector,
        expected_orbits=(119, 120) if int(sector) == 56 else None,
        expected_detectors=(
            tuple((camera, ccd) for camera in range(1, 5) for ccd in range(1, 5))
            if int(sector) == 56
            else None
        ),
    )
    if injection_metrics_path is not None:
        if injection_manifest_path is None:
            raise ValueError("precomputed injection metrics require an injection manifest")
        injection_metrics = _read_table(injection_metrics_path)
        injection_manifest = json.loads(injection_manifest_path.read_text())
        if str(injection_manifest.get("metrics_file_sha256", "")) != file_sha256(
            injection_metrics_path
        ):
            raise ValueError("precomputed injection metrics file hash mismatch")
    else:
        if not injection_shards:
            raise ValueError("provide injection metrics or at least one injection shard")
        injection_metrics, injection_manifest = summarize_fixed_injection_shards(
            injection_shards,
            sector=sector,
            cadence_reference_path=cadence_reference_path,
            cadence_reference_manifest_path=cadence_reference_manifest_path,
            apertures=config.apertures,
        )
    independent_metrics = _read_table(independent_metrics_path)
    independent_manifest = json.loads(independent_manifest_path.read_text())
    if str(independent_manifest.get("metrics_file_sha256", "")) != file_sha256(
        independent_metrics_path
    ):
        raise ValueError("independent metrics file hash mismatch")
    injection_preflight = evaluate_fixed_injections(
        injection_metrics, injection_manifest, sector=sector, config=config
    )
    if injection_preflight["status"] == "fail":
        raise ValueError(
            "fixed-injection preflight failed before population scan: "
            + "; ".join(injection_preflight["reasons"])
        )
    independent_preflight = evaluate_independent_extraction(
        independent_metrics,
        independent_manifest,
        config,
        sector=sector,
        catalog_detectors=catalog_detectors,
        current_compact_sha256=compact_sha256,
    )
    if independent_preflight["status"] == "fail":
        raise ValueError(
            "independent-extraction preflight failed before population scan: "
            + "; ".join(independent_preflight["reasons"])
        )
    injection_source_parity_report = json.loads(
        injection_source_parity_path.read_text()
    )
    injection_source_parity_gate = evaluate_injection_source_parity(
        injection_source_parity_report,
        injection_source_parity_path,
        injection_manifest,
        current_compact_sha256=compact_sha256,
        config=config,
    )
    if injection_source_parity_gate["status"] != "pass":
        raise ValueError(
            "injection-source parity preflight failed before population scan: "
            + "; ".join(injection_source_parity_gate["reasons"])
        )
    target_metrics, aperture_metrics, pair_metrics = audit_compact_population(
        compact_lc,
        apertures=config.apertures,
        cadence_reference=cadence_reference,
        authority_exclusions=quality_reference.authority_exclusions,
        sample_size=config.sample_size,
        seed=config.sample_seed,
    )
    evaluated, magnitude_bins, targets_with_loss = evaluate_tier1_gates(
        sector=sector,
        config=config,
        tier0_summary=tier0_summary,
        tier0_summary_sha256=tier0_summary_sha256,
        target_metrics=target_metrics,
        aperture_metrics=aperture_metrics,
        pair_metrics=pair_metrics,
        injection_metrics=injection_metrics,
        injection_manifest=injection_manifest,
        independent_metrics=independent_metrics,
        independent_manifest=independent_manifest,
        current_compact_sha256=compact_sha256,
        catalog_detectors=catalog_detectors,
        cadence_reference_gate=cadence_reference_gate,
        injection_source_parity_gate=injection_source_parity_gate,
    )
    target_eligibility = build_target_eligibility(targets_with_loss, config)
    detector_summary = summarize_target_eligibility(target_eligibility)
    evaluated["target_qa"].update(
        {
            "observation_key": ["sector", "tic"],
            "reason_code_vocabulary": list(TIER1_TARGET_REASON_CODES),
            "n_with_gaia_dr3_source_id": int(
                target_eligibility["gaia_dr3_source_id"].notna().sum()
            ),
        }
    )

    assert_inputs_unchanged()
    if injection_shards:
        expected_direct_shard_hashes = {
            str(path): str(injection_manifest["shard_sha256"].get(str(path), ""))
            for path in injection_shards
        }
        observed_direct_shard_hashes = {
            str(path): initial_sha256(path) for path in injection_shards
        }
        if observed_direct_shard_hashes != expected_direct_shard_hashes:
            raise ValueError("one or more injection shards changed during the audit")

    out_dir.mkdir(parents=True, exist_ok=True)
    _write_table_atomic(targets_with_loss, paths["target_metrics"])
    _write_table_atomic(target_eligibility, paths["target_eligibility"])
    _write_table_atomic(aperture_metrics, paths["aperture_metrics"])
    _write_table_atomic(pair_metrics, paths["aperture_pair_metrics"])
    _write_table_atomic(detector_summary, paths["detector_summary"])
    _write_table_atomic(magnitude_bins, paths["magnitude_bins"])
    _write_table_atomic(injection_metrics, paths["injection_metrics"])
    _write_table_atomic(independent_metrics, paths["independent_metrics"])
    plot_tier1_diagnostics(
        targets=targets_with_loss,
        aperture_metrics=aperture_metrics,
        pair_metrics=pair_metrics,
        injection_metrics=injection_metrics,
        config=config,
        output_png=paths["qa_diagnostics_png"],
        output_pdf=paths["qa_diagnostics_pdf"],
    )
    plot_detector_eligibility(
        detector_summary,
        output_png=paths["detector_eligibility_png"],
        output_pdf=paths["detector_eligibility_pdf"],
    )
    assert_inputs_unchanged()

    # Publish the complete frozen-canary evidence, not merely a digest of a
    # mutable in-memory mapping.  Rebind the copied manifest to the exact CSV
    # written by this run so it can be independently audited or reused via the
    # precomputed-metrics path.
    published_injection_metrics = _read_table(paths["injection_metrics"])
    published_injection_manifest = dict(injection_manifest)
    published_injection_manifest.update(
        {
            "metrics_file": str(paths["injection_metrics"]),
            "metrics_file_sha256": file_sha256(paths["injection_metrics"]),
            "metrics_content_sha256": _dataframe_content_sha256(
                published_injection_metrics
            ),
            "source_metrics_file": (
                str(injection_metrics_path) if injection_metrics_path else None
            ),
            "source_metrics_file_sha256": (
                initial_sha256(injection_metrics_path)
                if injection_metrics_path is not None
                else None
            ),
            "source_manifest_file": (
                str(injection_manifest_path) if injection_manifest_path else None
            ),
            "source_manifest_file_sha256": (
                initial_sha256(injection_manifest_path)
                if injection_manifest_path is not None
                else None
            ),
        }
    )
    write_strict_json(paths["injection_manifest"], published_injection_manifest)

    provenance = {
        "config": str(config_path),
        "config_sha256": initial_sha256(config_path),
        "tier0_summary": str(tier0_summary_path),
        "tier0_summary_sha256": tier0_summary_sha256,
        "compact_lc": str(compact_lc),
        "compact_lc_sha256": compact_sha256,
        "cadence_reference": str(cadence_reference_path),
        "cadence_reference_sha256": initial_sha256(cadence_reference_path),
        "cadence_reference_manifest": str(cadence_reference_manifest_path),
        "cadence_reference_manifest_sha256": initial_sha256(
            cadence_reference_manifest_path
        ),
        "injection_source_parity": str(injection_source_parity_path),
        "injection_source_parity_sha256": initial_sha256(
            injection_source_parity_path
        ),
        "injection_metrics_input": (
            str(injection_metrics_path) if injection_metrics_path else None
        ),
        "injection_metrics_input_sha256": (
            initial_sha256(injection_metrics_path)
            if injection_metrics_path is not None
            else None
        ),
        "injection_manifest_input": (
            str(injection_manifest_path) if injection_manifest_path else None
        ),
        "injection_manifest_input_sha256": (
            initial_sha256(injection_manifest_path)
            if injection_manifest_path is not None
            else None
        ),
        "injection_shards": [str(path) for path in injection_shards],
        "injection_evidence_sha256": _sha256_text_mapping(injection_manifest),
        "published_injection_manifest_sha256": file_sha256(
            paths["injection_manifest"]
        ),
        "independent_metrics": str(independent_metrics_path),
        "independent_metrics_sha256": independent_metrics_sha256,
        "independent_manifest": str(independent_manifest_path),
        "independent_manifest_sha256": independent_manifest_sha256,
        "output_sha256": {
            name: file_sha256(path)
            for name, path in paths.items()
            if name != "summary"
        },
    }
    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "contract_version": TIER1_QA_CONTRACT_VERSION,
        "config_name": config.name,
        "qa_tier": "tier1_bounded_enrichment_qa",
        "scope": config.scope,
        "sector": int(sector),
        "status": evaluated["status"],
        "passed": evaluated["passed"],
        "enrichment_ready": evaluated["enrichment_ready"],
        "science_ready": evaluated["science_ready"],
        "promotion_enabled": config.promotion_enabled,
        "apertures": list(config.apertures),
        "target_qa": evaluated["target_qa"],
        "gates": evaluated["gates"],
        "thresholds": asdict(config),
        "provenance": provenance,
        "outputs": {name: str(path) for name, path in paths.items()}
        | {"gate_json": str(gate_json)},
    }
    summary_text = (
        json.dumps(_safe_json(summary), indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    )
    summary_temporary = paths["summary"].with_suffix(
        paths["summary"].suffix + ".publish.tmp"
    )
    gate_temporary = gate_json.with_suffix(gate_json.suffix + ".publish.tmp")
    try:
        summary_temporary.write_text(summary_text)
        gate_temporary.write_text(summary_text)
        assert_inputs_unchanged()
        summary_temporary.replace(paths["summary"])
        gate_temporary.replace(gate_json)
        assert_inputs_unchanged()
    except Exception:
        summary_temporary.unlink(missing_ok=True)
        gate_temporary.unlink(missing_ok=True)
        paths["summary"].unlink(missing_ok=True)
        gate_json.unlink(missing_ok=True)
        raise
    return _safe_json(summary)


__all__ = [
    "TIER1_QA_CONTRACT_VERSION",
    "TIER1_TARGET_REASON_CODES",
    "Tier1QAConfig",
    "attach_target_qa_flags",
    "audit_compact_population",
    "build_target_eligibility",
    "evaluate_aperture_outliers",
    "evaluate_cadence_reference_evidence",
    "evaluate_cadence_quality",
    "evaluate_fixed_injections",
    "evaluate_independent_extraction",
    "evaluate_injection_source_parity",
    "evaluate_population_scatter",
    "evaluate_tier0_prerequisite",
    "evaluate_tier1_gates",
    "injection_metadata_sha256",
    "load_tier1_config",
    "plot_detector_eligibility",
    "plot_tier1_diagnostics",
    "run_a2v1_tier1_qa",
    "summarize_target_eligibility",
    "summarize_fixed_injection_shards",
    "write_strict_json",
]
