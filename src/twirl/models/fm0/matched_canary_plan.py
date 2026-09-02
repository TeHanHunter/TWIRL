"""Freeze the identity-only evaluation plan for the matched FM0.3 canary.

The plan reuses the immutable S66--S77 development panel only for probe
training (S66--S71) and validation (S72--S74).  Its old S75--S77 block is not
primary evidence.  A fresh test is selected from the composite release's S77
``temporal_holdout`` role after excluding every component present anywhere in
the development panel.

Only CSV/JSON identity authorities are read.  No shard payload is opened, no
event is injected, no checkpoint is loaded, and no metric is computed here.
One shared schedule is emitted for both matched architectures, preserving one
native 200-second cadence per each of 128 encoder tokens.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from io import StringIO
from pathlib import Path, PurePosixPath
from typing import Any

from .composite_release import (
    COMPOSITE_RELEASE_SCHEMA_VERSION,
    COMPOSITE_RELEASE_STATE,
    EXCLUDED_OVERLAP_ROLE,
    HOLDOUT_ROLE,
    ROLE_INDEX_FIELDS,
    ROLE_INDEX_SCHEMA_VERSION,
    ROLE_ORDER,
    SOURCE_BINDING_FIELDS,
    SOURCE_BINDING_SCHEMA_VERSION,
    CompositeObservation,
)
from .dataset import VIEW_NAMES
from .later_sector_release import (
    LATER_SIX_VIEW_MANIFEST_FIELDS,
    LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION,
)
from .registry import FM0ContractError, sha256_file
from .temporal_panel import DEVELOPMENT_PARTITION
from .temporal_zero_shot import TEMPORAL_COHORTS, load_temporal_panel

MATCHED_CANARY_PLAN_SCHEMA_VERSION = "twirl_fm0_3_matched_canary_plan_v1"
MATCHED_CANARY_PLAN_RECEIPT_SCHEMA_VERSION = (
    "twirl_fm0_3_matched_canary_plan_receipt_v1"
)
MATCHED_CANARY_PLAN_READY_STATE = "FM0_3_MATCHED_CANARY_PLAN_READY"
MATCHED_CANARY_PLAN_SALT = "twirl_fm0_3_matched_canary_plan_v1"

ARCHITECTURES = ("TWIRL-FM0.3.1", "TWIRL-FM0.3.2")
SPLIT_ORDER = ("probe_train", "validation", "fresh_s77_test")
PROBE_TRAIN_SECTORS = tuple(range(66, 72))
VALIDATION_SECTORS = (72, 73, 74)
EXCLUDED_PANEL_SECTORS = (75, 76, 77)
TEST_SECTOR = 77
TARGET_COMPONENTS_PER_COHORT = 240

EVENT_DURATIONS = (1, 3, 9)
EVENT_DEPTHS_TEXT = ("0.01", "0.03", "0.1", "0.3")
EVENT_CELLS = tuple(
    (duration, depth) for duration in EVENT_DURATIONS for depth in EVENT_DEPTHS_TEXT
)
EVENT_CENTER_INDEX_ZERO_BASED = 64
CONTEXT_CADENCES = 128
NOMINAL_CADENCE_SECONDS = 200
REQUIRED_VIEW_INDICES = (2, 3)

SCHEDULE_FIELDS = (
    "schema_version",
    "event_pair_id",
    "split",
    "cohort",
    "selection_rank_one_based",
    "source_authority",
    "source_id",
    "source_partition",
    "source_role",
    "sector",
    "observation_key",
    "leakage_component_id",
    "duration_cadences",
    "fractional_depth",
    "event_center_index_zero_based",
    "context_cadences",
    "nominal_cadence_seconds",
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_COMPONENT = re.compile(r"^leakage_[0-9a-f]{64}$")
_OBSERVATION = re.compile(r"^observation_[0-9a-f]{64}$")


@dataclass(frozen=True, slots=True)
class MatchedCanaryPlanResult:
    """One immutable matched-canary identity and event schedule."""

    root: Path
    schedule_path: Path
    receipt_path: Path
    ready_path: Path
    schedule_sha256: str
    receipt_sha256: str
    receipt: Mapping[str, Any]


@dataclass(frozen=True, slots=True)
class _Identity:
    split: str
    cohort: str
    source_authority: str
    source_id: str
    source_partition: str
    source_role: str
    sector: int
    observation_key: str
    component: str


@dataclass(frozen=True, slots=True)
class _CompositeIdentityAuthority:
    receipt: Mapping[str, Any]
    holdout: tuple[CompositeObservation, ...]
    repeated_components: frozenset[str]


def _digest(value: Any, *, label: str) -> str:
    digest = str(value).strip().lower()
    if _SHA256.fullmatch(digest) is None:
        raise FM0ContractError(f"{label} must be a lowercase SHA-256")
    return digest


def _git_sha(value: Any) -> str:
    sha = str(value).strip().lower()
    if _GIT_SHA.fullmatch(sha) is None:
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    return sha


def _stable_hex(*parts: object) -> str:
    return hashlib.sha256(
        "\x1f".join((MATCHED_CANARY_PLAN_SALT, *(str(part) for part in parts))).encode(
            "utf-8"
        )
    ).hexdigest()


def _identity_order(split: str, cohort: str, identity: str) -> tuple[str, str]:
    return _stable_hex("identity", split, cohort, identity), identity


def _visit_order(split: str, observation: str) -> tuple[str, str]:
    return _stable_hex("visit", split, observation), observation


def _view_presence(value: Any) -> tuple[bool, ...]:
    try:
        parsed = json.loads(str(value))
    except (TypeError, json.JSONDecodeError) as exc:
        raise FM0ContractError("matched-canary view_present_json is invalid") from exc
    if (
        not isinstance(parsed, list)
        or len(parsed) != len(VIEW_NAMES)
        or any(type(item) is not int or item not in (0, 1) for item in parsed)
    ):
        raise FM0ContractError("matched-canary view_present_json is invalid")
    return tuple(bool(item) for item in parsed)


def _required_views_present(value: Any) -> bool:
    present = _view_presence(value)
    return all(present[index] for index in REQUIRED_VIEW_INDICES)


def _check_identity(observation: str, component: str) -> None:
    if _OBSERVATION.fullmatch(observation) is None:
        raise FM0ContractError("matched-canary observation identity is invalid")
    if _COMPONENT.fullmatch(component) is None:
        raise FM0ContractError("matched-canary component identity is invalid")


def _hash_set(values: Sequence[str] | set[str]) -> str:
    return hashlib.sha256(("\n".join(sorted(values)) + "\n").encode()).hexdigest()


def _csv_payload(rows: Sequence[Mapping[str, Any]]) -> bytes:
    stream = StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=list(SCHEDULE_FIELDS), lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        if tuple(row) != SCHEDULE_FIELDS:
            raise FM0ContractError("matched-canary schedule columns drifted")
        writer.writerow(row)
    return stream.getvalue().encode("utf-8")


def _read_schedule(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != SCHEDULE_FIELDS:
            raise FM0ContractError("matched-canary schedule columns drifted")
        rows = [dict(row) for row in reader]
    if not rows:
        raise FM0ContractError("matched-canary schedule is empty")
    return rows


def _bound_file(
    path: str | Path,
    *,
    expected_sha256: str,
    label: str,
    require_read_only: bool,
) -> Path:
    raw = Path(path).expanduser()
    if raw.is_symlink():
        raise FM0ContractError(f"{label} must be materialized")
    source = raw.resolve()
    if not source.is_file() or source.stat().st_size <= 0:
        raise FM0ContractError(f"{label} is missing or empty")
    if require_read_only and source.stat().st_mode & 0o222:
        raise FM0ContractError(f"{label} is writable")
    if sha256_file(source) != _digest(expected_sha256, label=f"{label} hash"):
        raise FM0ContractError(f"{label} hash drifted")
    return source


def _read_exact_csv(
    path: Path, fields: Sequence[str], *, label: str
) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise FM0ContractError(f"{label} columns drifted")
        rows = [dict(row) for row in reader]
    if not rows:
        raise FM0ContractError(f"{label} is empty")
    return rows


def _load_composite_identity_authority(
    root: str | Path,
    *,
    receipt_sha256: str,
    source_bindings_sha256: str,
    role_index_sha256: str,
    require_read_only: bool,
) -> _CompositeIdentityAuthority:
    """Load only accepted identity tables plus S77 manifest metadata.

    This deliberately does not call the full composite validator because that
    validator checks every source shard-directory closure.  Exact accepted
    hashes make the compact identity tables authoritative for this freezer.
    """

    raw_root = Path(root).expanduser()
    if raw_root.is_symlink():
        raise FM0ContractError("composite root must be materialized")
    composite_root = raw_root.resolve()
    if not composite_root.is_dir() or {
        path.name for path in composite_root.iterdir()
    } != {"source_bindings.csv", "role_index.csv", "receipt.json", "READY"}:
        raise FM0ContractError("composite identity root closure drifted")
    if require_read_only and composite_root.stat().st_mode & 0o222:
        raise FM0ContractError("composite identity root is writable")
    receipt_path = _bound_file(
        composite_root / "receipt.json",
        expected_sha256=receipt_sha256,
        label="composite receipt",
        require_read_only=require_read_only,
    )
    source_path = _bound_file(
        composite_root / "source_bindings.csv",
        expected_sha256=source_bindings_sha256,
        label="composite source bindings",
        require_read_only=require_read_only,
    )
    role_path = _bound_file(
        composite_root / "role_index.csv",
        expected_sha256=role_index_sha256,
        label="composite role index",
        require_read_only=require_read_only,
    )
    ready = composite_root / "READY"
    if (
        ready.is_symlink()
        or not ready.is_file()
        or (require_read_only and ready.stat().st_mode & 0o222)
        or ready.read_text(encoding="utf-8").strip() != receipt_sha256
    ):
        raise FM0ContractError("composite READY/receipt binding drifted")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError("composite receipt is invalid") from exc
    if not isinstance(receipt, Mapping):
        raise FM0ContractError("composite receipt must be an object")
    sources = receipt.get("sources")
    selection = receipt.get("selection")
    limits = receipt.get("limits")
    if (
        receipt.get("schema_version") != COMPOSITE_RELEASE_SCHEMA_VERSION
        or receipt.get("release_state") != COMPOSITE_RELEASE_STATE
        or receipt.get("passed") is not True
        or not isinstance(sources, Mapping)
        or not isinstance(selection, Mapping)
        or sources.get("source_bindings_sha256") != source_bindings_sha256
        or selection.get("role_index_sha256") != role_index_sha256
        or selection.get("temporal_holdout_sector") != TEST_SECTOR
        or not isinstance(limits, Mapping)
        or limits.get("identity_only") is not True
        or limits.get("source_shards_opened") is not False
        or limits.get("sealed_rows_selected") != 0
    ):
        raise FM0ContractError("composite receipt identity boundary drifted")

    bindings = _read_exact_csv(
        source_path, SOURCE_BINDING_FIELDS, label="composite source bindings"
    )
    if len(bindings) != 13:
        raise FM0ContractError("composite source-binding count drifted")
    s77 = [row for row in bindings if row["source_id"] == "later_s0077"]
    if (
        len(s77) != 1
        or s77[0]["schema_version"] != SOURCE_BINDING_SCHEMA_VERSION
        or s77[0]["source_kind"] != "later"
        or s77[0]["sector_min"] != "77"
        or s77[0]["sector_max"] != "77"
    ):
        raise FM0ContractError("composite S77 source binding drifted")
    s77_binding = s77[0]
    s77_root_raw = Path(s77_binding["release_root"]).expanduser()
    if s77_root_raw.is_symlink():
        raise FM0ContractError("composite S77 source root must be materialized")
    s77_root = s77_root_raw.resolve()
    if not s77_root.is_dir() or (require_read_only and s77_root.stat().st_mode & 0o222):
        raise FM0ContractError("composite S77 source root is unavailable or writable")
    manifest = _bound_file(
        s77_root / "manifest.csv",
        expected_sha256=s77_binding["manifest_sha256"],
        label="composite S77 manifest",
        require_read_only=require_read_only,
    )

    role_counts = {role: 0 for role in ROLE_ORDER}
    holdout_roles: dict[str, tuple[str, str]] = {}
    repeated_components: set[str] = set()
    with role_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != ROLE_INDEX_FIELDS:
            raise FM0ContractError("composite role-index columns drifted")
        for row in reader:
            observation = str(row["observation_key"])
            component = str(row["leakage_component_id"])
            role = str(row["role"])
            _check_identity(observation, component)
            if (
                row["schema_version"] != ROLE_INDEX_SCHEMA_VERSION
                or role not in role_counts
            ):
                raise FM0ContractError("composite role-index row drifted")
            role_counts[role] += 1
            if role == HOLDOUT_ROLE:
                if row["source_id"] != "later_s0077" or observation in holdout_roles:
                    raise FM0ContractError("composite S77 holdout identity drifted")
                holdout_roles[observation] = (component, row["source_id"])
            elif role == EXCLUDED_OVERLAP_ROLE:
                repeated_components.add(component)
    if role_counts != selection.get("role_counts") or not holdout_roles:
        raise FM0ContractError("composite role-index counts drifted")

    manifest_rows = _read_exact_csv(
        manifest, LATER_SIX_VIEW_MANIFEST_FIELDS, label="composite S77 manifest"
    )
    if len(manifest_rows) != int(s77_binding["n_rows"]):
        raise FM0ContractError("composite S77 manifest count drifted")
    manifest_train_keys = {
        row["observation_key"]
        for row in manifest_rows
        if row["source_partition"] == "poc_train"
    }
    if manifest_train_keys != set(holdout_roles):
        raise FM0ContractError("composite S77 holdout/manifest closure drifted")
    holdout: list[CompositeObservation] = []
    for row in manifest_rows:
        if row["source_partition"] != "poc_train":
            continue
        observation = row["observation_key"]
        component, source_id = holdout_roles[observation]
        relative = PurePosixPath(row["relative_path"])
        _check_identity(observation, component)
        if (
            row["manifest_schema_version"] != LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION
            or row["sector"] != "77"
            or row["leakage_component_id"] != component
            or relative.is_absolute()
            or ".." in relative.parts
            or relative.parts != ("shards", f"{observation}.npz")
        ):
            raise FM0ContractError("composite S77 manifest identity drifted")
        try:
            n_cadences = int(row["n_cadences"])
        except ValueError as exc:
            raise FM0ContractError("composite S77 cadence count is invalid") from exc
        if n_cadences <= 0:
            raise FM0ContractError("composite S77 cadence count is invalid")
        holdout.append(
            CompositeObservation(
                source_id=source_id,
                sector=TEST_SECTOR,
                observation_key=observation,
                component=component,
                relative_path=relative.as_posix(),
                shard_sha256=_digest(
                    row["sha256"], label="composite S77 shard metadata hash"
                ),
                n_cadences=n_cadences,
                view_present_json=json.dumps(
                    [int(value) for value in _view_presence(row["view_present_json"])],
                    separators=(",", ":"),
                ),
            )
        )
    return _CompositeIdentityAuthority(
        receipt=receipt,
        holdout=tuple(holdout),
        repeated_components=frozenset(repeated_components),
    )


def _frozen_analysis_contract() -> dict[str, Any]:
    """Return the exact analysis contract shared by both architectures."""

    return {
        "task": "candidate_centered_clean_vs_single_event_classification",
        "probe": {
            "family": "full_batch_adam_linear_logit",
            "fit_split": "probe_train",
            "fit_sectors": list(PROBE_TRAIN_SECTORS),
            "validation_split": "validation",
            "validation_sectors": list(VALIDATION_SECTORS),
            "validation_role": "diagnostic_only_no_tuning_or_selection",
            "test_split": "fresh_s77_test",
            "test_sector": TEST_SECTOR,
            "standardization_fit_on_probe_train_only": True,
            "test_used_for_fit_tuning_or_threshold": False,
            "one_separate_fitted_weight_vector_per_feature_arm": True,
            "optimizer": "adam",
            "epochs": 400,
            "learning_rate": 0.02,
            "l2_weight": 0.001,
            "initialization_seed": 560306,
            "hyperparameter_tuning": False,
        },
        "center_readout": {
            "candidate_center_index_zero_based": EVENT_CENTER_INDEX_ZERO_BASED,
            "encoder_tokens_required_before_readout": CONTEXT_CADENCES,
            "encoder_token_stride": 1,
            "score_definition": "linear_logit_from_h_cadence_token_64_only",
            "all_128_tokens_exist_before_center_token_selection": True,
            "token_or_window_pooling_before_score": False,
            "off_center_tokens_used_as_probe_inputs": False,
        },
        "arms": {
            "controls": [
                "raw_adp_validity_error_exact_center",
                "quality_only_exact_center",
            ],
            "TWIRL-FM0.3.1": [
                "step0_h_cadence_token64",
                "step2000_h_cadence_token64",
            ],
            "TWIRL-FM0.3.2": [
                "step0_h_cadence_token64",
                "step2000_h_cadence_token64",
            ],
            "same_schedule_fit_and_metric_code_for_every_arm": True,
        },
        "feature_boundaries": {
            "synthetic_support_used_as_input": False,
            "event_duration_used_as_input": False,
            "event_depth_used_as_input": False,
            "event_center_numeric_value_used_as_input": False,
            "bls_or_period_features_used_as_input": False,
            "search_or_candidate_score_used_as_input": False,
            "duration_depth_and_support_role": (
                "injection_generation_and_prespecified_metric_stratification_only"
            ),
        },
        "metrics": {
            "primary": "sample_roc_auc",
            "diagnostic": "paired_clean_injected_ranking_accuracy",
            "report_overall": True,
            "report_by_cohort": True,
            "report_by_duration_depth_cell": True,
            "report_by_depth_aggregated_over_durations": True,
            "confidence_interval": {
                "level": 0.95,
                "method": "deterministic_cluster_bootstrap",
                "replicates": 1000,
                "seed": 560305,
                "resampling_unit": ("whole_leakage_component_with_clean_injected_pair"),
                "same_resample_indices_across_arms_and_paired_deltas": True,
            },
        },
        "primary_test_gates": {
            "gate_split": "fresh_s77_test",
            "raw_control": {
                "minimum_overall_sample_roc_auc_lower_95": 0.90,
            },
            "quality_only_control": {
                "maximum_absolute_overall_sample_roc_auc_from_chance": 0.03,
                "overall_sample_roc_auc_interval_must_contain": 0.50,
            },
            "each_architecture_step2000": {
                "minimum_overall_sample_roc_auc_lower_95": 0.75,
                "minimum_each_cohort_sample_roc_auc_estimate": 0.70,
                "minimum_each_cohort_sample_roc_auc_lower_95_strictly_above": (0.50),
                "paired_step2000_minus_own_step0_auc": {
                    "minimum_estimate": 0.02,
                    "lower_95_strictly_above": 0.0,
                },
                "blocking_depth_aggregates": {
                    "fractional_depths": [0.1, 0.3],
                    "minimum_each_depth_sample_roc_auc_lower_95": 0.80,
                    "aggregate_over_durations": list(EVENT_DURATIONS),
                },
                "nonblocking_depth_aggregates": {
                    "fractional_depths": [0.01, 0.03],
                    "report_with_intervals": True,
                },
            },
        },
        "decision_rule": {
            "one_seed": 560067,
            "raw_and_quality_control_gates_must_pass_before_fm_interpretation": (True),
            "one_seed_pass_authorizes_only": (
                "longer_matched_two_seed_architecture_comparison"
            ),
            "one_seed_pass_authorizes_architecture_promotion": False,
            "one_seed_pass_authorizes_foundation_model_claim": False,
            "architecture_comparison_interpretable_only_if_both_step2000_arms_"
            "pass_all_matched_gates": True,
            "otherwise_architecture_result": "inconclusive",
        },
    }


def _selected_panel_identities(
    rows: Sequence[Mapping[str, str]],
) -> tuple[list[_Identity], Mapping[str, Any], set[str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    component_cohort: dict[str, str] = {}
    component_blocks: dict[str, set[str]] = defaultdict(set)
    seen_observations: set[str] = set()
    all_panel_components: set[str] = set()

    for source in rows:
        row = dict(source)
        observation = str(row.get("observation_key", ""))
        component = str(row.get("leakage_component_id", ""))
        cohort = str(row.get("temporal_cohort", ""))
        try:
            sector = int(str(row.get("sector", "")))
        except ValueError as exc:
            raise FM0ContractError("matched-canary panel sector is invalid") from exc
        _check_identity(observation, component)
        if observation in seen_observations:
            raise FM0ContractError("matched-canary panel observations are not unique")
        seen_observations.add(observation)
        if (
            cohort not in TEMPORAL_COHORTS
            or row.get("source_partition") != DEVELOPMENT_PARTITION
            or sector
            not in (*PROBE_TRAIN_SECTORS, *VALIDATION_SECTORS, *EXCLUDED_PANEL_SECTORS)
        ):
            raise FM0ContractError("matched-canary panel row is outside its authority")
        prior = component_cohort.setdefault(component, cohort)
        if prior != cohort:
            raise FM0ContractError("matched-canary panel component crosses cohorts")
        all_panel_components.add(component)
        if sector in PROBE_TRAIN_SECTORS:
            block = "probe_train"
        elif sector in VALIDATION_SECTORS:
            block = "validation"
        else:
            block = "excluded_panel"
        component_blocks[component].add(block)
        grouped[component].append(row)

    cross_block = {
        component
        for component, blocks in component_blocks.items()
        if "probe_train" in blocks and "validation" in blocks
    }
    selected: list[_Identity] = []
    availability: dict[str, dict[str, int]] = {}
    for split, sectors in (
        ("probe_train", PROBE_TRAIN_SECTORS),
        ("validation", VALIDATION_SECTORS),
    ):
        availability[split] = {}
        for cohort in TEMPORAL_COHORTS:
            usable: dict[str, list[dict[str, str]]] = {}
            for component, component_rows in grouped.items():
                if component in cross_block or component_cohort[component] != cohort:
                    continue
                candidates = [
                    row
                    for row in component_rows
                    if int(row["sector"]) in sectors
                    and _required_views_present(row.get("view_present_json"))
                ]
                if candidates:
                    usable[component] = candidates
            ordered = sorted(
                usable,
                key=lambda component: _identity_order(split, cohort, component),
            )
            availability[split][cohort] = len(ordered)
            if len(ordered) < TARGET_COMPONENTS_PER_COHORT:
                raise FM0ContractError(
                    f"matched-canary {split}/{cohort} has only {len(ordered)} "
                    f"eligible components; {TARGET_COMPONENTS_PER_COHORT} required"
                )
            for component in ordered[:TARGET_COMPONENTS_PER_COHORT]:
                row = min(
                    usable[component],
                    key=lambda value: _visit_order(split, value["observation_key"]),
                )
                selected.append(
                    _Identity(
                        split=split,
                        cohort=cohort,
                        source_authority="temporal_panel",
                        source_id=f"later_s{int(row['sector']):04d}",
                        source_partition=DEVELOPMENT_PARTITION,
                        source_role="development_evaluation",
                        sector=int(row["sector"]),
                        observation_key=row["observation_key"],
                        component=component,
                    )
                )
    audit = {
        "available_components_after_view_and_cross_block_filters": availability,
        "quarantined_train_validation_components": len(cross_block),
        "quarantined_train_validation_components_sha256": _hash_set(cross_block),
        "excluded_panel_sectors": list(EXCLUDED_PANEL_SECTORS),
        "excluded_panel_rows": sum(
            int(row["sector"]) in EXCLUDED_PANEL_SECTORS for row in rows
        ),
    }
    return selected, audit, all_panel_components


def _selected_test_identities(
    holdout_observations: Sequence[CompositeObservation],
    *,
    repeated_components: set[str] | frozenset[str],
    panel_components: set[str],
) -> tuple[list[_Identity], Mapping[str, Any]]:
    grouped: dict[str, list[CompositeObservation]] = defaultdict(list)
    seen_observations: set[str] = set()
    excluded_panel_overlap: set[str] = set()
    excluded_missing_views: set[str] = set()
    for item in holdout_observations:
        _check_identity(item.observation_key, item.component)
        if item.observation_key in seen_observations:
            raise FM0ContractError("composite temporal holdout observations duplicate")
        seen_observations.add(item.observation_key)
        if item.sector != TEST_SECTOR or item.source_id != "later_s0077":
            raise FM0ContractError("composite temporal holdout is not exactly S77")
        if item.component in panel_components:
            excluded_panel_overlap.add(item.component)
            continue
        if not _required_views_present(item.view_present_json):
            excluded_missing_views.add(item.component)
            continue
        grouped[item.component].append(item)

    selected: list[_Identity] = []
    availability: dict[str, int] = {}
    for cohort in TEMPORAL_COHORTS:
        candidates = [
            component
            for component in grouped
            if (component in repeated_components) == (cohort == "repeated")
        ]
        candidates.sort(
            key=lambda component: _identity_order("fresh_s77_test", cohort, component)
        )
        availability[cohort] = len(candidates)
        if len(candidates) < TARGET_COMPONENTS_PER_COHORT:
            raise FM0ContractError(
                f"matched-canary fresh_s77_test/{cohort} has only "
                f"{len(candidates)} eligible components; "
                f"{TARGET_COMPONENTS_PER_COHORT} required"
            )
        for component in candidates[:TARGET_COMPONENTS_PER_COHORT]:
            item = min(
                grouped[component],
                key=lambda value: _visit_order("fresh_s77_test", value.observation_key),
            )
            selected.append(
                _Identity(
                    split="fresh_s77_test",
                    cohort=cohort,
                    source_authority="composite_release",
                    source_id=item.source_id,
                    source_partition="poc_train",
                    source_role=HOLDOUT_ROLE,
                    sector=TEST_SECTOR,
                    observation_key=item.observation_key,
                    component=component,
                )
            )
    return selected, {
        "available_components_after_freshness_and_view_filters": availability,
        "cohort_definition": {
            "repeated": "component_has_earlier_excluded_temporal_overlap_role",
            "new": "component_has_no_earlier_excluded_temporal_overlap_role",
        },
        "excluded_panel_overlap_components": len(excluded_panel_overlap),
        "excluded_panel_overlap_components_sha256": _hash_set(excluded_panel_overlap),
        "excluded_missing_required_view_components": len(excluded_missing_views),
        "excluded_missing_required_view_components_sha256": _hash_set(
            excluded_missing_views
        ),
    }


def _schedule_rows(identities: Sequence[_Identity]) -> list[dict[str, Any]]:
    split_components: dict[str, set[str]] = defaultdict(set)
    rows: list[dict[str, Any]] = []
    for split in SPLIT_ORDER:
        for cohort in TEMPORAL_COHORTS:
            selected = [
                identity
                for identity in identities
                if identity.split == split and identity.cohort == cohort
            ]
            if len(selected) != TARGET_COMPONENTS_PER_COHORT:
                raise FM0ContractError("matched-canary selected identity count drifted")
            selected.sort(
                key=lambda item: _identity_order(split, cohort, item.component)
            )
            for rank, identity in enumerate(selected, start=1):
                duration, depth = EVENT_CELLS[(rank - 1) % len(EVENT_CELLS)]
                pair_hash = _stable_hex(
                    "event_pair",
                    split,
                    cohort,
                    identity.observation_key,
                    duration,
                    depth,
                    EVENT_CENTER_INDEX_ZERO_BASED,
                )
                rows.append(
                    {
                        "schema_version": MATCHED_CANARY_PLAN_SCHEMA_VERSION,
                        "event_pair_id": f"event_pair_{pair_hash}",
                        "split": split,
                        "cohort": cohort,
                        "selection_rank_one_based": str(rank),
                        "source_authority": identity.source_authority,
                        "source_id": identity.source_id,
                        "source_partition": identity.source_partition,
                        "source_role": identity.source_role,
                        "sector": str(identity.sector),
                        "observation_key": identity.observation_key,
                        "leakage_component_id": identity.component,
                        "duration_cadences": str(duration),
                        "fractional_depth": depth,
                        "event_center_index_zero_based": str(
                            EVENT_CENTER_INDEX_ZERO_BASED
                        ),
                        "context_cadences": str(CONTEXT_CADENCES),
                        "nominal_cadence_seconds": str(NOMINAL_CADENCE_SECONDS),
                    }
                )
                split_components[split].add(identity.component)
    for left_index, left in enumerate(SPLIT_ORDER):
        for right in SPLIT_ORDER[left_index + 1 :]:
            if split_components[left] & split_components[right]:
                raise FM0ContractError(
                    f"matched-canary components overlap between {left} and {right}"
                )
    observations = [row["observation_key"] for row in rows]
    components = [row["leakage_component_id"] for row in rows]
    pair_ids = [row["event_pair_id"] for row in rows]
    if (
        len(observations) != len(set(observations))
        or len(components) != len(set(components))
        or len(pair_ids) != len(set(pair_ids))
    ):
        raise FM0ContractError("matched-canary schedule identities are not unique")
    return rows


def _products(
    *,
    temporal_panel_root: str | Path,
    temporal_panel_receipt_sha256: str,
    composite_root: str | Path,
    composite_receipt_sha256: str,
    composite_source_bindings_sha256: str,
    composite_role_index_sha256: str,
    producer_git_sha: str,
    require_read_only: bool,
) -> tuple[bytes, dict[str, Any]]:
    producer = _git_sha(producer_git_sha)
    panel_receipt_hash = _digest(
        temporal_panel_receipt_sha256, label="temporal-panel receipt hash"
    )
    composite_receipt_hash = _digest(
        composite_receipt_sha256, label="composite receipt hash"
    )
    source_bindings_hash = _digest(
        composite_source_bindings_sha256, label="composite source-bindings hash"
    )
    role_index_hash = _digest(
        composite_role_index_sha256, label="composite role-index hash"
    )

    panel_rows, panel_authority = load_temporal_panel(
        temporal_panel_root, receipt_sha256=panel_receipt_hash
    )
    composite = _load_composite_identity_authority(
        composite_root,
        receipt_sha256=composite_receipt_hash,
        source_bindings_sha256=source_bindings_hash,
        role_index_sha256=role_index_hash,
        require_read_only=require_read_only,
    )
    panel_selected, panel_audit, panel_components = _selected_panel_identities(
        panel_rows
    )
    test_selected, test_audit = _selected_test_identities(
        composite.holdout,
        repeated_components=composite.repeated_components,
        panel_components=panel_components,
    )
    schedule_rows = _schedule_rows((*panel_selected, *test_selected))
    schedule_payload = _csv_payload(schedule_rows)
    schedule_hash = hashlib.sha256(schedule_payload).hexdigest()

    count_by_split_cohort: dict[str, dict[str, int]] = {}
    cells_by_split_cohort: dict[str, dict[str, dict[str, int]]] = {}
    for split in SPLIT_ORDER:
        count_by_split_cohort[split] = {}
        cells_by_split_cohort[split] = {}
        for cohort in TEMPORAL_COHORTS:
            subset = [
                row
                for row in schedule_rows
                if row["split"] == split and row["cohort"] == cohort
            ]
            count_by_split_cohort[split][cohort] = len(subset)
            cell_counts = {
                f"duration_{duration}_depth_{depth}": sum(
                    row["duration_cadences"] == str(duration)
                    and row["fractional_depth"] == depth
                    for row in subset
                )
                for duration, depth in EVENT_CELLS
            }
            if set(cell_counts.values()) != {
                TARGET_COMPONENTS_PER_COHORT // len(EVENT_CELLS)
            }:
                raise FM0ContractError("matched-canary event cells are not balanced")
            cells_by_split_cohort[split][cohort] = cell_counts

    panel_receipt_actual = _digest(
        panel_authority.get("receipt_sha256"), label="loaded temporal-panel receipt"
    )
    if panel_receipt_actual != panel_receipt_hash:
        raise FM0ContractError("loaded temporal-panel receipt binding drifted")
    composite_receipt = composite.receipt
    receipt = {
        "schema_version": MATCHED_CANARY_PLAN_RECEIPT_SCHEMA_VERSION,
        "ready_state": MATCHED_CANARY_PLAN_READY_STATE,
        "passed": True,
        "producer_git_sha": producer,
        "input_authorities": {
            "temporal_panel": {
                "root": str(Path(temporal_panel_root).expanduser().resolve()),
                "receipt_sha256": panel_receipt_hash,
                "panel_sha256": panel_authority["panel_sha256"],
                "sector_bindings_sha256": panel_authority["sector_bindings_sha256"],
            },
            "composite_release": {
                "root": str(Path(composite_root).expanduser().resolve()),
                "receipt_sha256": composite_receipt_hash,
                "source_bindings_sha256": source_bindings_hash,
                "role_index_sha256": role_index_hash,
                "source_producer_git_sha": composite_receipt["producer_git_sha"],
            },
        },
        "split_contract": {
            "probe_train_sectors": list(PROBE_TRAIN_SECTORS),
            "validation_sectors": list(VALIDATION_SECTORS),
            "excluded_old_panel_sectors": list(EXCLUDED_PANEL_SECTORS),
            "fresh_test_sector": TEST_SECTOR,
            "fresh_test_source_role": HOLDOUT_ROLE,
            "one_visit_per_component": True,
            "component_spanning_probe_train_and_validation": "quarantine",
            "fresh_test_component_definition": (
                "composite_s77_temporal_holdout_component_absent_from_entire_"
                "frozen_temporal_panel"
            ),
            "pairwise_component_disjoint": True,
            "target_components_per_cohort_per_split": (TARGET_COMPONENTS_PER_COHORT),
            "cohorts": list(TEMPORAL_COHORTS),
            "counts": count_by_split_cohort,
            "panel_audit": panel_audit,
            "fresh_test_audit": test_audit,
        },
        "event_contract": {
            "task": "candidate_centered_single_event_classification",
            "period_defined": False,
            "bls_or_search_features_used": False,
            "duration_cadences": list(EVENT_DURATIONS),
            "fractional_depths": [float(value) for value in EVENT_DEPTHS_TEXT],
            "event_center_index_zero_based": EVENT_CENTER_INDEX_ZERO_BASED,
            "one_clean_and_one_injected_sample_per_schedule_row": True,
            "one_event_cell_per_component": True,
            "balanced_event_cells": True,
            "event_cell_counts": cells_by_split_cohort,
        },
        "cadence_contract": {
            "nominal_cadence_seconds": NOMINAL_CADENCE_SECONDS,
            "context_cadences": CONTEXT_CADENCES,
            "one_encoder_token_per_native_cadence": True,
            "patch_stride": 1,
            "cadence_averaging_or_merging": False,
            "temporal_downsampling": False,
            "representation_pooling": False,
        },
        "architecture_comparison": {
            "variants": list(ARCHITECTURES),
            "shared_schedule_path": "schedule.csv",
            "shared_schedule_sha256": schedule_hash,
            "identical_schedule_across_architectures": True,
            "architecture_specific_identity_or_event_selection": False,
        },
        "analysis_contract": _frozen_analysis_contract(),
        "schedule": {
            "path": "schedule.csv",
            "sha256": schedule_hash,
            "n_rows": len(schedule_rows),
        },
        "limits": {
            "identity_only": True,
            "s77_shard_payloads_opened": False,
            "any_shard_payloads_opened": False,
            "events_injected": False,
            "checkpoints_loaded": False,
            "probe_fitted": False,
            "metrics_computed": False,
            "sealed_test_opened": False,
            "model_training_authorized": False,
        },
        "claim_limit": (
            "Identity-only matched-canary evaluation plan; not an evaluation, "
            "architecture result, model promotion, BLS/search result, sealed-test "
            "result, or training authorization."
        ),
    }
    return schedule_payload, receipt


def _make_tree_read_only(root: Path) -> None:
    for path in root.rglob("*"):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def validate_matched_canary_evaluation_plan(
    root: str | Path,
    *,
    expected_receipt_sha256: str,
    require_read_only: bool = True,
) -> MatchedCanaryPlanResult:
    """Validate the immutable output bundle without opening source shards."""

    output = Path(root).expanduser().resolve()
    if (
        output.is_symlink()
        or not output.is_dir()
        or {path.name for path in output.iterdir()}
        != {"schedule.csv", "receipt.json", "READY"}
    ):
        raise FM0ContractError("matched-canary plan root closure drifted")
    receipt_path = output / "receipt.json"
    schedule_path = output / "schedule.csv"
    ready_path = output / "READY"
    for path in (receipt_path, schedule_path, ready_path):
        if path.is_symlink() or not path.is_file():
            raise FM0ContractError("matched-canary plan artifact is not materialized")
        if require_read_only and path.stat().st_mode & 0o222:
            raise FM0ContractError("matched-canary plan artifact is writable")
    if require_read_only and output.stat().st_mode & 0o222:
        raise FM0ContractError("matched-canary plan root is writable")
    receipt_hash = sha256_file(receipt_path)
    if receipt_hash != _digest(
        expected_receipt_sha256, label="matched-canary receipt hash"
    ):
        raise FM0ContractError("matched-canary receipt hash drifted")
    if ready_path.read_text(encoding="utf-8").strip() != receipt_hash:
        raise FM0ContractError("matched-canary READY/receipt binding drifted")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError("matched-canary receipt is invalid") from exc
    if (
        not isinstance(receipt, Mapping)
        or receipt.get("schema_version") != MATCHED_CANARY_PLAN_RECEIPT_SCHEMA_VERSION
        or receipt.get("ready_state") != MATCHED_CANARY_PLAN_READY_STATE
        or receipt.get("passed") is not True
    ):
        raise FM0ContractError("matched-canary receipt contract drifted")
    schedule_hash = sha256_file(schedule_path)
    if (
        receipt.get("schedule", {}).get("path") != "schedule.csv"
        or receipt.get("schedule", {}).get("sha256") != schedule_hash
        or receipt.get("architecture_comparison", {}).get("shared_schedule_sha256")
        != schedule_hash
    ):
        raise FM0ContractError("matched-canary schedule binding drifted")
    rows = _read_schedule(schedule_path)
    if len(rows) != receipt["schedule"]["n_rows"]:
        raise FM0ContractError("matched-canary schedule count drifted")
    split_components: dict[str, set[str]] = defaultdict(set)
    seen_observations: set[str] = set()
    seen_pairs: set[str] = set()
    for row in rows:
        _check_identity(row["observation_key"], row["leakage_component_id"])
        if (
            row["schema_version"] != MATCHED_CANARY_PLAN_SCHEMA_VERSION
            or row["split"] not in SPLIT_ORDER
            or row["cohort"] not in TEMPORAL_COHORTS
            or int(row["duration_cadences"]) not in EVENT_DURATIONS
            or row["fractional_depth"] not in EVENT_DEPTHS_TEXT
            or int(row["event_center_index_zero_based"])
            != EVENT_CENTER_INDEX_ZERO_BASED
            or int(row["context_cadences"]) != CONTEXT_CADENCES
            or int(row["nominal_cadence_seconds"]) != NOMINAL_CADENCE_SECONDS
        ):
            raise FM0ContractError("matched-canary schedule row contract drifted")
        if (
            row["observation_key"] in seen_observations
            or row["event_pair_id"] in seen_pairs
        ):
            raise FM0ContractError("matched-canary schedule identities duplicate")
        seen_observations.add(row["observation_key"])
        seen_pairs.add(row["event_pair_id"])
        split_components[row["split"]].add(row["leakage_component_id"])
    for left_index, left in enumerate(SPLIT_ORDER):
        for right in SPLIT_ORDER[left_index + 1 :]:
            if split_components[left] & split_components[right]:
                raise FM0ContractError("matched-canary split components overlap")
    limits = receipt.get("limits")
    if limits != {
        "identity_only": True,
        "s77_shard_payloads_opened": False,
        "any_shard_payloads_opened": False,
        "events_injected": False,
        "checkpoints_loaded": False,
        "probe_fitted": False,
        "metrics_computed": False,
        "sealed_test_opened": False,
        "model_training_authorized": False,
    }:
        raise FM0ContractError("matched-canary plan boundary drifted")
    comparison = receipt.get("architecture_comparison", {})
    if (
        tuple(comparison.get("variants", ())) != ARCHITECTURES
        or comparison.get("identical_schedule_across_architectures") is not True
        or comparison.get("architecture_specific_identity_or_event_selection")
        is not False
    ):
        raise FM0ContractError("matched-canary architecture schedule drifted")
    if receipt.get("analysis_contract") != _frozen_analysis_contract():
        raise FM0ContractError("matched-canary analysis contract drifted")
    return MatchedCanaryPlanResult(
        root=output,
        schedule_path=schedule_path,
        receipt_path=receipt_path,
        ready_path=ready_path,
        schedule_sha256=schedule_hash,
        receipt_sha256=receipt_hash,
        receipt=receipt,
    )


def freeze_matched_canary_evaluation_plan(
    output_dir: str | Path,
    *,
    temporal_panel_root: str | Path,
    temporal_panel_receipt_sha256: str,
    composite_root: str | Path,
    composite_receipt_sha256: str,
    composite_source_bindings_sha256: str,
    composite_role_index_sha256: str,
    producer_git_sha: str,
    require_read_only: bool = True,
) -> MatchedCanaryPlanResult:
    """Freeze one shared, identity-only TCN/Conformer evaluation schedule."""

    schedule_payload, receipt = _products(
        temporal_panel_root=temporal_panel_root,
        temporal_panel_receipt_sha256=temporal_panel_receipt_sha256,
        composite_root=composite_root,
        composite_receipt_sha256=composite_receipt_sha256,
        composite_source_bindings_sha256=composite_source_bindings_sha256,
        composite_role_index_sha256=composite_role_index_sha256,
        producer_git_sha=producer_git_sha,
        require_read_only=require_read_only,
    )
    receipt_payload = (
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    receipt_hash = hashlib.sha256(receipt_payload).hexdigest()
    output = Path(output_dir).expanduser().absolute()
    if output.exists():
        return validate_matched_canary_evaluation_plan(
            output,
            expected_receipt_sha256=receipt_hash,
            require_read_only=require_read_only,
        )
    if output.is_symlink() or output.name in {"", ".", ".."}:
        raise FM0ContractError("matched-canary output must be a named directory")
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.partial.", dir=output.parent)
    )
    try:
        (partial / "schedule.csv").write_bytes(schedule_payload)
        (partial / "receipt.json").write_bytes(receipt_payload)
        (partial / "READY").write_text(receipt_hash + "\n", encoding="utf-8")
        _make_tree_read_only(partial)
        os.replace(partial, output)
    except BaseException:
        if partial.exists():
            partial.chmod(0o755)
            for path in partial.rglob("*"):
                path.chmod(0o755 if path.is_dir() else 0o644)
            shutil.rmtree(partial)
        raise
    return validate_matched_canary_evaluation_plan(
        output,
        expected_receipt_sha256=receipt_hash,
        require_read_only=require_read_only,
    )
