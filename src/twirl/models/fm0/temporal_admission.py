"""Fail-closed, label-blind admission of genuinely later FM0 sectors.

The inventory produced here is a control artifact.  It verifies accepted A2v1
receipts and identity authorities, classifies later visits relative to the
complete S56--S64 training-era corpus, and never opens light curves, shards,
search products, labels, embeddings, model scores, or injections.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from statistics import median
from typing import Any

import yaml

from .input_release import FLUX_VIEW_NAMES
from .registry import (
    SPLIT_SALT,
    AliasRegistry,
    FM0ContractError,
    assert_no_search_columns,
    build_alias_registry,
    publish_immutable,
    read_rows,
    sha256_file,
)

CONFIG_SCHEMA_VERSION = "twirl_fm0_2_later_sector_admission_config_v1"
PRODUCTION_CAMPAIGN_ID = "twirl_fm0_2_later_sector_admission_v1"
# Updated only when a separately reviewed production policy version is frozen.
FROZEN_POLICY_SHA256 = (
    "527a7b4d9f9c452f02576eea7f155abe65ad8439057317e8bc10c80b0ed93da3"
)
SECTOR_RECEIPT_SCHEMA_VERSION = "twirl_fm0_later_sector_admission_receipt_v1"
TEMPORAL_SECTOR_RECEIPT_SCHEMA_VERSION = SECTOR_RECEIPT_SCHEMA_VERSION
INVENTORY_SCHEMA_VERSION = "twirl_fm0_later_sector_inventory_v1"
SECTOR_INVENTORY_SCHEMA_VERSION = "twirl_fm0_later_sector_inventory_sector_v1"
QUARANTINE_SCHEMA_VERSION = "twirl_fm0_later_sector_quarantine_v1"
EVIDENCE_SCHEMA_VERSION = "twirl_fm0_later_sector_evidence_v1"
EVIDENCE_METRIC_CONTRACT_VERSION = "twirl_fm0_later_sector_evidence_metrics_v1"
UPSTREAM_CONTROL_SCHEMA_VERSION = "twirl_fm0_later_sector_upstream_control_v1"
OBSERVATION_SCHEMA_VERSION = "twirl_fm0_later_sector_observation_v1"
MINIMUM_REPEATED_HOST_COMPONENTS = 64
MINIMUM_NEW_HOST_COMPONENTS = 256
REQUIRED_EVIDENCE_GATES = (
    "checksum_bound_a2v1_hdf5_provenance",
    "checksum_bound_a2v1_fits_provenance",
    "edge_aware_coverage",
    "hdf5_openability",
    "authoritative_internal_cadence_quality",
    "authoritative_external_cadence_quality",
    "explicit_omissions",
    "fm_channel_mask_finite_numerical_envelope",
    "stable_physical_source_registry_join",
)
EVIDENCE_METRIC_FIELDS = {
    "checksum_bound_a2v1_hdf5_provenance": (
        "n_products",
        "n_checksum_bound_products",
        "n_checksum_mismatches",
    ),
    "checksum_bound_a2v1_fits_provenance": (
        "n_products",
        "n_checksum_bound_products",
        "n_checksum_mismatches",
    ),
    "edge_aware_coverage": (
        "n_expected_observations",
        "n_present_observations",
        "n_edge_omissions",
        "n_non_edge_omissions",
    ),
    "hdf5_openability": (
        "n_hdf5_products",
        "n_hdf5_opened",
        "n_unreadable_hdf5",
    ),
    "authoritative_internal_cadence_quality": (
        "n_observations_checked",
        "n_cadences_checked",
        "n_observations_failed",
    ),
    "authoritative_external_cadence_quality": (
        "n_observations_checked",
        "n_cadences_checked",
        "n_observations_failed",
    ),
    "explicit_omissions": (
        "n_expected_observations",
        "n_present_observations",
        "n_explicit_omissions",
        "n_unexplained_omissions",
    ),
    "fm_channel_mask_finite_numerical_envelope": (
        "n_observations_checked",
        "n_six_view_observations",
        "n_observations_failed",
        "n_nonfinite_active_values",
    ),
    "stable_physical_source_registry_join": (
        "n_identity_rows",
        "n_joined_identity_rows",
        "n_quarantined_identity_rows",
        "n_unmatched_identity_rows",
    ),
}

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_POSITIVE_INTEGER = re.compile(r"^[0-9]+$")

LATER_OBSERVATION_FIELDS = (
    "schema_version",
    "observation_key",
    "product_instance_id",
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "camera",
    "ccd",
    "orbits_json",
    "leakage_component_id",
    "source_partition",
    "host_cohort",
    "cohort_partition",
    "prior_training_era_visit_count",
    "source_sha256",
    "product_state",
    "n_cadences",
    "cadence_retention_fraction",
    "flux_view_names_json",
    "view_present_json",
    "scientific_training_eligible",
)
SECTOR_INVENTORY_FIELDS = (
    "schema_version",
    "sector",
    "expected_orbits_json",
    "admission_receipt_sha256",
    "stage1_acceptance_receipt_sha256",
    "observation_manifest_sha256",
    "alias_authority_sha256",
    "evidence_sha256_json",
    "evidence_metrics_json",
    "evidence_upstream_sha256_json",
    "n_identity_rows",
    "n_eligible_rows",
    "n_repeated_host_rows",
    "n_new_host_rows",
    "n_sealed_identity_rows_excluded",
    "n_quarantined_rows",
    "n_repeated_host_components",
    "n_new_host_components",
    "n_camera_ccd_cells",
    "camera_ccd_coverage_json",
    "cadence_retention_fraction_min",
    "cadence_retention_fraction_median",
)
QUARANTINE_FIELDS = (
    "schema_version",
    "observation_key",
    "sector",
    "gaia_dr3_source_id",
    "tic_id",
    "leakage_component_id",
    "source_partition",
    "reason",
)

_LATER_REQUIRED_FIELDS = frozenset(
    {
        "observation_key",
        "schema_version",
        "product_instance_id",
        "gaia_dr3_source_id",
        "tic_id",
        "sector",
        "camera",
        "ccd",
        "orbits_json",
        "leakage_component_id",
        "source_partition",
        "source_sha256",
        "product_state",
        "n_cadences",
        "cadence_retention_fraction",
        "flux_view_names_json",
        "view_present_json",
        "scientific_training_eligible",
        "quarantined",
    }
)
_LATER_ALLOWED_FIELDS = _LATER_REQUIRED_FIELDS
_BASELINE_MANIFEST_REQUIRED_FIELDS = frozenset(
    {"observation_key", "leakage_component_id", "source_partition"}
)
_BASELINE_SELECTION_REQUIRED_FIELDS = frozenset(
    {
        "gaia_dr3_source_id",
        "tic_id",
        "sector",
        "camera",
        "ccd",
        "leakage_component_id",
        "source_partition",
    }
)


@dataclass(frozen=True)
class VerifiedSectorAdmission:
    sector: int
    receipt_sha256: str
    stage1_acceptance_receipt_sha256: str
    observation_manifest_sha256: str
    alias_authority_sha256: str
    evidence_sha256: Mapping[str, str]
    evidence_metrics: Mapping[str, Mapping[str, int]]
    evidence_upstream_sha256: Mapping[str, Mapping[str, str]]
    alias_rows: tuple[dict[str, str], ...]
    observation_rows: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class BaselineAuthority:
    alias_registry: AliasRegistry
    alias_rows: tuple[dict[str, str], ...]
    component_visit_counts: Mapping[str, int]
    component_partitions: Mapping[str, str]
    component_detector_visits: Mapping[str, tuple[str, ...]]
    observation_keys: frozenset[str]


@dataclass(frozen=True)
class TemporalAdmissionInventory:
    eligible_observations: tuple[dict[str, Any], ...]
    sector_inventory: tuple[dict[str, Any], ...]
    quarantine: tuple[dict[str, Any], ...]
    summary: Mapping[str, Any]


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise FM0ContractError(f"{label} must be a mapping")
    return value


def _sequence(value: Any, *, label: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise FM0ContractError(f"{label} must be a sequence")
    return value


def _boolean(value: Any, *, label: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no", ""}:
            return False
    raise FM0ContractError(f"{label} must be boolean")


def _identifier(value: Any, *, label: str, optional: bool = False) -> str:
    text = "" if value is None else str(value).strip()
    if optional and not text:
        return ""
    if not _POSITIVE_INTEGER.fullmatch(text) or int(text) <= 0:
        raise FM0ContractError(f"{label} must be a positive decimal integer")
    return str(int(text))


def _digest(value: Any, *, label: str) -> str:
    result = str(value).strip().lower()
    if not _SHA256.fullmatch(result):
        raise FM0ContractError(f"{label} must be a lowercase SHA-256 digest")
    return result


def _verify_file(path: str | Path, expected: str, *, label: str) -> Path:
    source = Path(path).expanduser().resolve()
    expected = _digest(expected, label=f"{label} SHA-256")
    if not source.is_file():
        raise FM0ContractError(f"missing {label}: {source}")
    observed = sha256_file(source)
    if observed != expected:
        raise FM0ContractError(
            f"{label} hash mismatch: expected {expected}, observed {observed}"
        )
    return source


def _read_table(
    path: Path, *, label: str
) -> tuple[list[dict[str, Any]], tuple[str, ...]]:
    try:
        rows = read_rows(path)
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not rows:
        raise FM0ContractError(f"{label} is empty: {path}")
    fields = tuple(rows[0])
    if any(tuple(row) != fields for row in rows):
        raise FM0ContractError(f"{label} columns drift between rows")
    return rows, fields


def _read_json(path: Path, *, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    return _mapping(value, label=label)


def _verify_stage1_acceptance_receipt(
    path: Path, *, sector: int, expected_orbits: tuple[int, ...]
) -> None:
    """Require the unchanged accepted A2v1 full-schema validation contract."""

    receipt = _read_json(path, label=f"S{sector} Stage-1 acceptance receipt")
    try:
        receipt_sector = int(receipt.get("sector"))
        observed_orbits = tuple(int(value) for value in receipt.get("orbits", ()))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(
            f"S{sector} Stage-1 acceptance identity is invalid"
        ) from exc
    if receipt_sector != sector or observed_orbits != expected_orbits:
        raise FM0ContractError(f"S{sector} Stage-1 acceptance identity drifted")
    if any(receipt.get(key) is not True for key in ("ok", "ok_h5", "ok_fits")):
        raise FM0ContractError(f"S{sector} Stage-1 acceptance receipt did not pass")

    contract = _mapping(
        receipt.get("expected_contract"),
        label=f"S{sector} Stage-1 expected contract",
    )
    try:
        requested = tuple(int(value) for value in contract.get("requested_orbits", ()))
        observed = tuple(int(value) for value in contract.get("observed_orbits", ()))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(f"S{sector} Stage-1 orbit contract is invalid") from exc
    if (
        contract.get("ok") is not True
        or contract.get("has_expected_rows") is not True
        or contract.get("has_expected_unique_tics") is not True
        or list(contract.get("missing_requested_orbits", ()))
        or requested != expected_orbits
        or observed != expected_orbits
    ):
        raise FM0ContractError(f"S{sector} Stage-1 expected contract did not pass")

    h5 = _mapping(receipt.get("h5"), label=f"S{sector} Stage-1 HDF5 gate")
    fits = _mapping(receipt.get("fits"), label=f"S{sector} Stage-1 FITS gate")
    schema = _mapping(receipt.get("schema"), label=f"S{sector} Stage-1 schema gate")
    try:
        h5_counts = (
            int(h5.get("n_present_h5", 0)),
            int(h5.get("n_missing_h5_non_edge", -1)),
            int(h5.get("n_unreadable_h5", -1)),
            int(h5.get("n_zero_byte_h5", -1)),
        )
        fits_counts = (
            int(fits.get("n_fits", 0)),
            int(fits.get("n_checked_fits", -1)),
            int(fits.get("n_bad_checked_fits", -1)),
            int(fits.get("n_missing_fits_non_edge_tics", -1)),
        )
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(
            f"S{sector} Stage-1 acceptance counts are invalid"
        ) from exc
    if h5_counts[0] <= 0 or h5_counts[1:] != (0, 0, 0):
        raise FM0ContractError(f"S{sector} Stage-1 HDF5 acceptance counts failed")
    if (
        fits_counts[0] <= 0
        or fits_counts[1] != fits_counts[0]
        or fits_counts[2:] != (0, 0)
    ):
        raise FM0ContractError(f"S{sector} Stage-1 FITS acceptance counts failed")
    if (
        schema.get("schema_only") is not True
        or schema.get("check_h5_open") is not True
        or schema.get("expected_method") != "A2v1"
        or schema.get("expected_prodtag") != "A2v1"
    ):
        raise FM0ContractError(f"S{sector} Stage-1 full-schema gate drifted")


def _evidence_metrics(gate: str, raw: Any, *, sector: int) -> dict[str, int]:
    metrics = _mapping(raw, label=f"S{sector} evidence gate {gate} metrics")
    expected_fields = EVIDENCE_METRIC_FIELDS[gate]
    if set(metrics) != set(expected_fields):
        raise FM0ContractError(f"S{sector} evidence gate {gate} metrics drifted")
    values: dict[str, int] = {}
    for field in expected_fields:
        value = metrics[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} metric {field} is invalid"
            )
        values[field] = value

    if gate in {
        "checksum_bound_a2v1_hdf5_provenance",
        "checksum_bound_a2v1_fits_provenance",
    }:
        if (
            values["n_products"] <= 0
            or values["n_checksum_bound_products"] != values["n_products"]
            or values["n_checksum_mismatches"] != 0
        ):
            raise FM0ContractError(f"S{sector} evidence gate {gate} counts failed")
    elif gate == "edge_aware_coverage":
        if (
            values["n_expected_observations"] <= 0
            or values["n_present_observations"] <= 0
            or values["n_non_edge_omissions"] != 0
            or values["n_expected_observations"]
            != values["n_present_observations"]
            + values["n_edge_omissions"]
            + values["n_non_edge_omissions"]
        ):
            raise FM0ContractError(f"S{sector} evidence gate {gate} counts failed")
    elif gate == "hdf5_openability":
        if (
            values["n_hdf5_products"] <= 0
            or values["n_hdf5_opened"] != values["n_hdf5_products"]
            or values["n_unreadable_hdf5"] != 0
        ):
            raise FM0ContractError(f"S{sector} evidence gate {gate} counts failed")
    elif gate in {
        "authoritative_internal_cadence_quality",
        "authoritative_external_cadence_quality",
    }:
        if (
            values["n_observations_checked"] <= 0
            or values["n_cadences_checked"] <= 0
            or values["n_observations_failed"] != 0
        ):
            raise FM0ContractError(f"S{sector} evidence gate {gate} counts failed")
    elif gate == "explicit_omissions":
        if (
            values["n_expected_observations"] <= 0
            or values["n_present_observations"] <= 0
            or values["n_unexplained_omissions"] != 0
            or values["n_expected_observations"]
            != values["n_present_observations"]
            + values["n_explicit_omissions"]
            + values["n_unexplained_omissions"]
        ):
            raise FM0ContractError(f"S{sector} evidence gate {gate} counts failed")
    elif gate == "fm_channel_mask_finite_numerical_envelope":
        if (
            values["n_observations_checked"] <= 0
            or values["n_six_view_observations"] <= 0
            or values["n_six_view_observations"] > values["n_observations_checked"]
            or values["n_observations_failed"] != 0
            or values["n_nonfinite_active_values"] != 0
        ):
            raise FM0ContractError(f"S{sector} evidence gate {gate} counts failed")
    elif gate == "stable_physical_source_registry_join":
        if (
            values["n_identity_rows"] <= 0
            or values["n_joined_identity_rows"] <= 0
            or values["n_unmatched_identity_rows"] != 0
            or values["n_identity_rows"]
            != values["n_joined_identity_rows"]
            + values["n_quarantined_identity_rows"]
            + values["n_unmatched_identity_rows"]
        ):
            raise FM0ContractError(f"S{sector} evidence gate {gate} counts failed")
    else:  # pragma: no cover - guarded by the frozen gate tuple
        raise FM0ContractError(f"unknown evidence gate {gate}")
    return values


def _verify_evidence_upstream_artifacts(
    *,
    evidence_path: Path,
    gate: str,
    sector: int,
    raw: Any,
    expected_metrics: Mapping[str, int],
    minimum_count: int,
    global_paths: set[Path],
    global_digests: set[str],
) -> dict[str, str]:
    bindings = _sequence(raw, label=f"S{sector} evidence gate {gate} upstream")
    if len(bindings) < minimum_count:
        raise FM0ContractError(
            f"S{sector} evidence gate {gate} lacks upstream control artifacts"
        )
    result: dict[str, str] = {}
    for index, raw_binding in enumerate(bindings):
        binding = _mapping(
            raw_binding,
            label=f"S{sector} evidence gate {gate} upstream artifact {index}",
        )
        if set(binding) != {"role", "path", "sha256", "content_class"}:
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream binding drifted"
            )
        role = str(binding.get("role", "")).strip()
        if not role or role in result:
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream role is invalid"
            )
        if binding.get("content_class") != "control_metadata":
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream content is not control metadata"
            )
        source, digest = _resolve_binding(
            evidence_path,
            binding,
            label=f"S{sector} evidence gate {gate} upstream {role}",
        )
        if (
            source == evidence_path
            or source in global_paths
            or digest in global_digests
        ):
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} reuses an upstream control artifact"
            )
        if source.suffix.lower() != ".json" or source.stat().st_size <= 0:
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream control must be nonempty JSON"
            )
        control = _read_json(
            source,
            label=f"S{sector} evidence gate {gate} upstream control {role}",
        )
        if set(control) != {
            "schema_version",
            "gate_name",
            "sector",
            "product_state",
            "passed",
            "metric_contract_version",
            "metrics",
            "producer",
        }:
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream control schema drifted"
            )
        try:
            control_sector = int(control.get("sector"))
        except (TypeError, ValueError) as exc:
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream sector is invalid"
            ) from exc
        if (
            control.get("schema_version") != UPSTREAM_CONTROL_SCHEMA_VERSION
            or control.get("gate_name") != gate
            or control_sector != sector
            or control.get("product_state") != "A2V1_ACCEPTED"
            or control.get("passed") is not True
            or control.get("metric_contract_version")
            != EVIDENCE_METRIC_CONTRACT_VERSION
        ):
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream control identity drifted"
            )
        control_metrics = _evidence_metrics(
            gate,
            control.get("metrics"),
            sector=sector,
        )
        if control_metrics != dict(expected_metrics):
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream metrics differ"
            )

        producer = _mapping(
            control.get("producer"),
            label=f"S{sector} evidence gate {gate} upstream producer",
        )
        if set(producer) != {"git_commit", "code", "config"} or (
            _GIT_SHA.fullmatch(str(producer.get("git_commit", ""))) is None
        ):
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} upstream producer drifted"
            )
        producer_paths: set[Path] = set()
        for producer_role in ("code", "config"):
            producer_binding = _mapping(
                producer.get(producer_role),
                label=(
                    f"S{sector} evidence gate {gate} upstream producer {producer_role}"
                ),
            )
            if set(producer_binding) != {"path", "sha256"}:
                raise FM0ContractError(
                    f"S{sector} evidence gate {gate} producer binding drifted"
                )
            producer_path, _ = _resolve_binding(
                source,
                producer_binding,
                label=(
                    f"S{sector} evidence gate {gate} upstream producer {producer_role}"
                ),
            )
            if (
                producer_path.stat().st_size <= 0
                or producer_path in producer_paths
                or producer_path.suffix.lower() in {".h5", ".hdf5", ".npz", ".npy"}
            ):
                raise FM0ContractError(
                    f"S{sector} evidence gate {gate} producer artifact is invalid"
                )
            producer_paths.add(producer_path)
        global_paths.add(source)
        global_digests.add(digest)
        result[role] = digest
    return dict(sorted(result.items()))


def _load_config(
    path: str | Path, *, expected_sha256: str
) -> tuple[Mapping[str, Any], str]:
    source = _verify_file(
        path,
        expected_sha256,
        label="later-sector frozen config",
    )
    try:
        config = yaml.safe_load(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise FM0ContractError(f"invalid later-sector config: {source}") from exc
    config = _mapping(config, label="later-sector config")
    if config.get("schema_version") != CONFIG_SCHEMA_VERSION:
        raise FM0ContractError("later-sector config schema mismatch")
    observed_hash = _digest(
        expected_sha256,
        label="later-sector config SHA-256",
    )
    if (
        config.get("campaign_id") == PRODUCTION_CAMPAIGN_ID
        and observed_hash != FROZEN_POLICY_SHA256
    ):
        raise FM0ContractError("production later-sector policy hash drifted")
    _validate_config(config)
    return config, observed_hash


def _validate_config(config: Mapping[str, Any]) -> None:
    baseline = _mapping(config.get("baseline"), label="baseline policy")
    identity = _mapping(config.get("identity"), label="identity policy")
    selection = _mapping(config.get("selection"), label="selection policy")
    admission = _mapping(config.get("admission"), label="admission policy")
    inventory = _mapping(config.get("inventory"), label="inventory policy")
    scope = _mapping(config.get("scope"), label="scope policy")
    later_manifest = _mapping(
        config.get("later_observation_manifest"),
        label="later observation manifest policy",
    )
    if scope.get("policy_only") is not True or any(
        scope.get(key) is not False
        for key in (
            "results_recorded",
            "panel_frozen",
            "model_training_authorized",
            "sealed_test_access_authorized",
            "event_retention_authorized",
            "formal_fm0_2_gate_authorized",
            "production_model_claim",
            "foundation_model_claim",
            "prospective_test_claim",
        )
    ):
        raise FM0ContractError(
            "inventory policy cannot authorize a panel or model result"
        )
    if (
        identity.get("split_salt") != SPLIT_SALT
        or int(identity.get("split_modulus", -1)) != 10000
    ):
        raise FM0ContractError("frozen source split authority drifted")
    if list(identity.get("allowed_source_partitions", ())) != [
        "poc_train",
        "poc_development",
    ] or list(identity.get("sealed_source_partitions", ())) != ["poc_sealed_test"]:
        raise FM0ContractError("source partition policy drifted")
    for key in ("manifest_sha256", "alias_authority_sha256", "selection_sha256"):
        _digest(baseline.get(key), label=f"baseline.{key}")
    if int(baseline.get("selection_n_observations", 0)) <= 0:
        raise FM0ContractError("baseline selection row count must be pinned")
    if (
        admission.get("required_receipt_schema_version")
        != SECTOR_RECEIPT_SCHEMA_VERSION
    ):
        raise FM0ContractError("later-sector receipt schema policy drifted")
    if admission.get("required_product_state") != "A2V1_ACCEPTED":
        raise FM0ContractError("later-sector product state policy drifted")
    if admission.get("stage1_acceptance_binding_key") != ("stage1_acceptance_receipt"):
        raise FM0ContractError("Stage-1 acceptance binding key drifted")
    if admission.get("observation_manifest_binding_key") != "observation_manifest":
        raise FM0ContractError("observation manifest binding key drifted")
    if admission.get("alias_authority_binding_key") != "alias_authority":
        raise FM0ContractError("later alias authority binding key drifted")
    required_true = (
        "fail_closed",
        "required_passed",
        "receipt_sha256_supplied_at_runtime",
        "stage1_acceptance_sha256_required",
        "observation_manifest_sha256_required",
        "alias_authority_sha256_required",
        "every_required_evidence_binding_must_pass",
        "distinct_evidence_artifacts_required",
    )
    if any(admission.get(key) is not True for key in required_true):
        raise FM0ContractError("later-sector admission is not fail-closed")
    if tuple(admission.get("required_evidence_gates", ())) != REQUIRED_EVIDENCE_GATES:
        raise FM0ContractError("later-sector required evidence gates drifted")
    if admission.get("required_evidence_schema_version") != EVIDENCE_SCHEMA_VERSION:
        raise FM0ContractError("later-sector evidence schema policy drifted")
    if admission.get("required_evidence_gate_name_field") != "gate_name":
        raise FM0ContractError("later-sector evidence gate-name policy drifted")
    if (
        admission.get("required_evidence_metric_contract_version")
        != EVIDENCE_METRIC_CONTRACT_VERSION
    ):
        raise FM0ContractError("later-sector evidence metric policy drifted")
    if (
        admission.get("required_evidence_upstream_schema_version")
        != UPSTREAM_CONTROL_SCHEMA_VERSION
    ):
        raise FM0ContractError("later-sector upstream control schema drifted")
    if tuple(admission.get("required_evidence_upstream_artifact_fields", ())) != (
        "role",
        "path",
        "sha256",
        "content_class",
    ) or admission.get("required_evidence_upstream_content_class") != (
        "control_metadata"
    ):
        raise FM0ContractError("later-sector evidence upstream policy drifted")
    if int(admission.get("minimum_upstream_artifacts_per_gate", 0)) < 1:
        raise FM0ContractError("later-sector evidence must bind upstream controls")
    if tuple(admission.get("evidence_binding_fields", ())) != (
        "path",
        "sha256",
        "passed",
    ):
        raise FM0ContractError("later-sector evidence binding fields drifted")
    required_false = (
        "deferred_orcd_checkpoint_is_accepted",
        "hdf5_only_sector_is_accepted",
        "diagnostic_product_state_is_accepted",
    )
    if any(admission.get(key) is not False for key in required_false):
        raise FM0ContractError("deferred/diagnostic products cannot be accepted")
    if later_manifest.get("required_schema_version") != OBSERVATION_SCHEMA_VERSION:
        raise FM0ContractError("later observation schema policy drifted")
    if (
        tuple(later_manifest.get("required_flux_view_names", ()))
        != tuple(FLUX_VIEW_NAMES)
        or later_manifest.get("exact_flux_view_order_required") is not True
    ):
        raise FM0ContractError("later observation six-view authority drifted")
    if selection.get("infer_orbits_from_sector_formula") is not False:
        raise FM0ContractError("orbit inference must remain forbidden")
    if int(selection.get("minimum_sectors", 0)) < 5:
        raise FM0ContractError(
            "temporal freeze policy must require at least five sectors"
        )
    if int(inventory.get("minimum_repeated_host_components", -1)) != (
        MINIMUM_REPEATED_HOST_COMPONENTS
    ) or int(inventory.get("minimum_new_host_components", -1)) != (
        MINIMUM_NEW_HOST_COMPONENTS
    ):
        raise FM0ContractError("temporal freeze-readiness cohort floors drifted")
    if (
        inventory.get("count_floors_are_not_panel_freeze_authority") is not True
        or inventory.get("adequacy_thresholds_frozen") is not False
        or inventory.get("panel_freeze_ready_must_remain_false") is not True
    ):
        raise FM0ContractError("unfrozen temporal adequacy policy drifted")
    candidates = tuple(
        int(value) for value in selection.get("ordered_candidate_sectors", ())
    )
    if not candidates or len(set(candidates)) != len(candidates):
        raise FM0ContractError("candidate sectors must be explicit and unique")
    if tuple(sorted(candidates)) != candidates or any(
        right != left + 1 for left, right in pairwise(candidates)
    ):
        raise FM0ContractError(
            "candidate-sector prefix must be chronological and contiguous"
        )
    orbit_map = _mapping(
        selection.get("expected_orbits_by_sector"), label="orbit mapping"
    )
    try:
        orbit_sectors = {int(value) for value in orbit_map}
    except (TypeError, ValueError) as exc:
        raise FM0ContractError("orbit-map sector keys must be integers") from exc
    if orbit_sectors != set(candidates):
        raise FM0ContractError(
            "explicit expected-orbit map must exactly cover candidates"
        )


def _forbidden_columns(
    fields: Iterable[str], config: Mapping[str, Any], *, label: str
) -> None:
    assert_no_search_columns(fields, context=label)
    tokens = [
        str(value).strip().lower() for value in config.get("forbidden_input_tokens", ())
    ]
    bad = [
        str(field)
        for field in fields
        if any(token and token in str(field).strip().lower() for token in tokens)
    ]
    if bad:
        raise FM0ContractError(f"{label} contains forbidden fields: {sorted(set(bad))}")


def _resolve_binding(receipt: Path, binding: Any, *, label: str) -> tuple[Path, str]:
    binding = _mapping(binding, label=f"{label} binding")
    expected = _digest(binding.get("sha256"), label=f"{label} SHA-256")
    text = str(binding.get("path", "")).strip()
    if not text:
        raise FM0ContractError(f"{label} binding lacks a path")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = receipt.parent / path
    return _verify_file(path, expected, label=f"{label} artifact"), expected


def _alias_edges(
    path: Path, config: Mapping[str, Any], *, label: str
) -> tuple[dict[str, str], ...]:
    rows, fields = _read_table(path, label=label)
    _forbidden_columns(fields, config, label=label)
    if set(fields) != {"gaia_dr3_source_id", "tic_id"}:
        raise FM0ContractError(f"{label} must be the edge-only Gaia/TIC authority")
    normalized: set[tuple[str, str]] = set()
    for row in rows:
        normalized.add(
            (
                _identifier(row.get("gaia_dr3_source_id"), label=f"{label} Gaia ID"),
                _identifier(row.get("tic_id"), label=f"{label} TIC ID", optional=True),
            )
        )
    if len(normalized) != len(rows):
        raise FM0ContractError(f"{label} contains duplicate Gaia/TIC edges")
    return tuple(
        {"gaia_dr3_source_id": gaia, "tic_id": tic}
        for gaia, tic in sorted(
            normalized, key=lambda item: (int(item[0]), int(item[1] or 0))
        )
    )


def _component_partitions(registry: AliasRegistry) -> dict[str, str]:
    return {
        str(row["leakage_component_id"]): str(row["source_partition"])
        for row in registry.components
    }


def _load_baseline(
    *,
    config: Mapping[str, Any],
    manifest_path: Path,
    alias_path: Path,
    selection_path: Path,
) -> BaselineAuthority:
    edges = _alias_edges(alias_path, config, label="baseline alias authority")
    registry = build_alias_registry(edges)
    aliases = registry.alias_index()
    partitions = _component_partitions(registry)

    selection_rows, selection_fields = _read_table(
        selection_path, label="baseline corpus selection"
    )
    _forbidden_columns(selection_fields, config, label="baseline corpus selection")
    missing = sorted(_BASELINE_SELECTION_REQUIRED_FIELDS - set(selection_fields))
    if missing:
        raise FM0ContractError(f"baseline corpus selection lacks fields: {missing}")
    baseline_policy = _mapping(config.get("baseline"), label="baseline policy")
    expected_n = int(baseline_policy["selection_n_observations"])
    if len(selection_rows) != expected_n:
        raise FM0ContractError(
            f"baseline selection row count drifted: expected {expected_n}, observed {len(selection_rows)}"
        )
    training_sectors = {
        int(value)
        for value in _sequence(
            baseline_policy.get("training_sectors"), label="training sectors"
        )
    }
    visit_counts: dict[str, int] = {}
    detector_visits: dict[str, list[str]] = {}
    selection_partitions: dict[str, str] = {}
    for row in selection_rows:
        gaia = _identifier(
            row.get("gaia_dr3_source_id"), label="baseline selection Gaia ID"
        )
        tic = _identifier(
            row.get("tic_id"), label="baseline selection TIC ID", optional=True
        )
        alias = aliases.get((gaia, tic))
        if alias is None:
            raise FM0ContractError(
                "baseline selection edge is absent from alias authority"
            )
        component = str(row.get("leakage_component_id", "")).strip()
        partition = str(row.get("source_partition", "")).strip()
        if (
            component != alias["leakage_component_id"]
            or partition != alias["source_partition"]
        ):
            raise FM0ContractError(
                "baseline selection alias/component assignment drifted"
            )
        if bool(alias["quarantined"]):
            raise FM0ContractError(
                "baseline selection contains a quarantined component"
            )
        try:
            sector, camera, ccd = (
                int(row["sector"]),
                int(row["camera"]),
                int(row["ccd"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise FM0ContractError(
                "baseline selection detector placement is invalid"
            ) from exc
        if (
            sector not in training_sectors
            or camera not in range(1, 5)
            or ccd not in range(1, 5)
        ):
            raise FM0ContractError(
                "baseline selection placement is outside frozen scope"
            )
        visit_counts[component] = visit_counts.get(component, 0) + 1
        detector_visits.setdefault(component, []).append(f"cam{camera}/ccd{ccd}")
        prior = selection_partitions.setdefault(component, partition)
        if prior != partition:
            raise FM0ContractError("baseline component crosses source partitions")

    manifest_rows, manifest_fields = _read_table(
        manifest_path, label="baseline input-release manifest"
    )
    _forbidden_columns(manifest_fields, config, label="baseline input-release manifest")
    missing = sorted(_BASELINE_MANIFEST_REQUIRED_FIELDS - set(manifest_fields))
    if missing:
        raise FM0ContractError(f"baseline input manifest lacks fields: {missing}")
    manifest_counts: dict[str, int] = {}
    observation_keys: set[str] = set()
    for row in manifest_rows:
        key = str(row.get("observation_key", "")).strip()
        if not key or key in observation_keys:
            raise FM0ContractError(f"duplicate/empty baseline observation key: {key!r}")
        observation_keys.add(key)
        component = str(row.get("leakage_component_id", "")).strip()
        partition = str(row.get("source_partition", "")).strip()
        if (
            component not in visit_counts
            or partition != selection_partitions[component]
        ):
            raise FM0ContractError(
                "baseline manifest component/partition differs from selection"
            )
        manifest_counts[component] = manifest_counts.get(component, 0) + 1
    if manifest_counts != visit_counts:
        raise FM0ContractError(
            "baseline manifest and corpus selection visit counts differ"
        )
    return BaselineAuthority(
        alias_registry=registry,
        alias_rows=edges,
        component_visit_counts=dict(sorted(visit_counts.items())),
        component_partitions=dict(sorted(partitions.items())),
        component_detector_visits={
            key: tuple(values) for key, values in sorted(detector_visits.items())
        },
        observation_keys=frozenset(observation_keys),
    )


def _expected_orbits(config: Mapping[str, Any], sector: int) -> tuple[int, ...]:
    selection = _mapping(config.get("selection"), label="selection policy")
    mapping = _mapping(
        selection.get("expected_orbits_by_sector"), label="orbit mapping"
    )
    raw = mapping.get(sector, mapping.get(str(sector)))
    if raw is None:
        raise FM0ContractError(f"S{sector} lacks an explicit expected-orbit mapping")
    values = tuple(int(value) for value in _sequence(raw, label=f"S{sector} orbits"))
    if (
        not values
        or len(set(values)) != len(values)
        or any(value <= 0 for value in values)
    ):
        raise FM0ContractError(f"S{sector} orbit mapping is invalid")
    return values


def _normalize_later_row(
    row: Mapping[str, Any], *, sector: int, expected_orbits: tuple[int, ...]
) -> dict[str, Any]:
    if str(row.get("schema_version", "")).strip() != OBSERVATION_SCHEMA_VERSION:
        raise FM0ContractError(f"S{sector} observation manifest schema mismatch")
    try:
        row_sector, camera, ccd = (
            int(row["sector"]),
            int(row["camera"]),
            int(row["ccd"]),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise FM0ContractError(f"S{sector} identity row has invalid placement") from exc
    if row_sector != sector or camera not in range(1, 5) or ccd not in range(1, 5):
        raise FM0ContractError(f"S{sector} identity row placement is outside scope")
    try:
        orbits = tuple(int(value) for value in json.loads(str(row["orbits_json"])))
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise FM0ContractError(
            f"S{sector} identity row has invalid orbits_json"
        ) from exc
    if orbits != expected_orbits:
        raise FM0ContractError(
            f"S{sector} identity row differs from explicit orbit mapping"
        )
    observation_key = str(row.get("observation_key", "")).strip()
    product_instance_id = str(row.get("product_instance_id", "")).strip()
    if not observation_key or not product_instance_id:
        raise FM0ContractError(f"S{sector} observation/product identity is required")
    source_sha256 = _digest(row.get("source_sha256"), label="source_sha256")
    if str(row.get("product_state", "")).strip() != "A2V1_ACCEPTED":
        raise FM0ContractError(f"S{sector} row is not A2V1_ACCEPTED")
    try:
        n_cadences = int(row.get("n_cadences"))
        retention = float(row.get("cadence_retention_fraction"))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(f"S{sector} cadence QA fields are invalid") from exc
    if n_cadences <= 0 or not math.isfinite(retention) or not 0.0 <= retention <= 1.0:
        raise FM0ContractError(f"S{sector} cadence QA fields are outside bounds")
    try:
        flux_view_names = tuple(
            str(value) for value in json.loads(str(row.get("flux_view_names_json", "")))
        )
    except (TypeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"S{sector} flux_view_names_json is invalid") from exc
    if flux_view_names != tuple(FLUX_VIEW_NAMES):
        raise FM0ContractError(f"S{sector} canonical six-view layout drifted")
    try:
        raw_views = json.loads(str(row.get("view_present_json", "")))
    except (TypeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"S{sector} view_present_json is invalid") from exc
    if not isinstance(raw_views, list) or len(raw_views) != 6:
        raise FM0ContractError(f"S{sector} view_present_json must describe six views")
    views: list[bool] = []
    for value in raw_views:
        if isinstance(value, bool):
            views.append(value)
        elif value in (0, 1):
            views.append(bool(value))
        else:
            raise FM0ContractError(f"S{sector} view availability must be boolean")
    return {
        "observation_key": observation_key,
        "product_instance_id": product_instance_id,
        "gaia_dr3_source_id": _identifier(
            row.get("gaia_dr3_source_id"), label="Gaia ID"
        ),
        "tic_id": _identifier(row.get("tic_id"), label="TIC ID"),
        "sector": sector,
        "camera": camera,
        "ccd": ccd,
        "orbits_json": json.dumps(orbits, separators=(",", ":")),
        "leakage_component_id": str(row.get("leakage_component_id", "")).strip(),
        "source_partition": str(row.get("source_partition", "")).strip(),
        "source_sha256": source_sha256,
        "product_state": "A2V1_ACCEPTED",
        "n_cadences": n_cadences,
        "cadence_retention_fraction": repr(retention),
        "flux_view_names_json": json.dumps(flux_view_names, separators=(",", ":")),
        "view_present_json": json.dumps(views, separators=(",", ":")),
        "all_six_views_present": all(views),
        "scientific_training_eligible": _boolean(
            row.get("scientific_training_eligible"),
            label="scientific_training_eligible",
        ),
        "quarantined": _boolean(row.get("quarantined", False), label="quarantined"),
    }


def verify_sector_admission_receipt(
    *,
    sector: int,
    receipt_path: str | Path,
    expected_receipt_sha256: str,
    config: Mapping[str, Any],
) -> VerifiedSectorAdmission:
    """Verify one accepted-sector receipt, evidence, aliases, and identity rows."""

    sector = int(sector)
    receipt_path = _verify_file(
        receipt_path, expected_receipt_sha256, label=f"S{sector} admission receipt"
    )
    receipt = _read_json(receipt_path, label=f"S{sector} admission receipt")
    admission = _mapping(config.get("admission"), label="admission policy")
    schema = receipt.get("receipt_schema_version", receipt.get("schema_version"))
    if schema != admission["required_receipt_schema_version"]:
        raise FM0ContractError(f"S{sector} admission receipt schema mismatch")
    try:
        receipt_sector = int(receipt.get("sector"))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError(f"S{sector} receipt sector is invalid") from exc
    if receipt_sector != sector:
        raise FM0ContractError(
            f"receipt sector S{receipt_sector} does not match the requested accepted sector S{sector}"
        )
    if receipt.get("product_state") != "A2V1_ACCEPTED":
        raise FM0ContractError(
            f"S{sector} receipt product_state is not the requested accepted sector state"
        )
    if receipt.get("passed") is not True:
        raise FM0ContractError(f"S{sector} admission receipt is not passed")

    expected_orbits = _expected_orbits(config, sector)
    stage1_binding = _mapping(
        receipt.get(admission["stage1_acceptance_binding_key"]),
        label=f"S{sector} Stage-1 acceptance binding",
    )
    if set(stage1_binding) != {"path", "sha256"}:
        raise FM0ContractError(f"S{sector} Stage-1 acceptance binding schema drifted")
    stage1_path, stage1_hash = _resolve_binding(
        receipt_path,
        stage1_binding,
        label=f"S{sector} Stage-1 acceptance receipt",
    )
    if stage1_path.suffix.lower() != ".json":
        raise FM0ContractError(f"S{sector} Stage-1 acceptance receipt must be JSON")
    _verify_stage1_acceptance_receipt(
        stage1_path,
        sector=sector,
        expected_orbits=expected_orbits,
    )

    required_gates = tuple(str(value) for value in admission["required_evidence_gates"])
    gates = _mapping(receipt.get("evidence_gates"), label=f"S{sector} evidence gates")
    if set(gates) != set(required_gates):
        raise FM0ContractError(f"S{sector} evidence gates differ from policy")
    evidence_hashes: dict[str, str] = {}
    evidence_metric_values: dict[str, Mapping[str, int]] = {}
    evidence_upstream_hashes: dict[str, Mapping[str, str]] = {}
    evidence_paths: set[Path] = set()
    evidence_digests: set[str] = set()
    upstream_paths: set[Path] = set()
    upstream_digests: set[str] = set()
    for gate in required_gates:
        binding = _mapping(gates[gate], label=f"S{sector} {gate} binding")
        if set(binding) != {"path", "sha256", "passed"}:
            raise FM0ContractError(f"S{sector} evidence gate {gate} binding drifted")
        if binding.get("passed") is not True:
            raise FM0ContractError(f"S{sector} evidence gate {gate} is not passed")
        evidence_path, evidence_hash = _resolve_binding(
            receipt_path, binding, label=f"S{sector} evidence gate {gate}"
        )
        if evidence_path.suffix.lower() != ".json":
            raise FM0ContractError(f"S{sector} evidence gate {gate} must be JSON")
        if evidence_path in evidence_paths or evidence_hash in evidence_digests:
            raise FM0ContractError(
                f"S{sector} evidence gate {gate} reuses another gate artifact"
            )
        evidence_paths.add(evidence_path)
        evidence_digests.add(evidence_hash)
        evidence = _read_json(evidence_path, label=f"S{sector} {gate} evidence")
        if set(evidence) != {
            "schema_version",
            "gate_name",
            "sector",
            "product_state",
            "passed",
            "metrics",
            "upstream_artifacts",
        }:
            raise FM0ContractError(f"S{sector} evidence gate {gate} schema drifted")
        if evidence.get("schema_version") != EVIDENCE_SCHEMA_VERSION:
            raise FM0ContractError(f"S{sector} evidence gate {gate} schema mismatch")
        if evidence.get("gate_name") != gate:
            raise FM0ContractError(f"S{sector} evidence gate {gate} name mismatch")
        if evidence.get("passed") is not True:
            raise FM0ContractError(f"S{sector} JSON evidence {gate} is not passed")
        if evidence.get("product_state") != "A2V1_ACCEPTED":
            raise FM0ContractError(
                f"S{sector} JSON evidence {gate} product_state is not A2V1_ACCEPTED"
            )
        try:
            evidence_sector = int(evidence.get("sector"))
        except (TypeError, ValueError) as exc:
            raise FM0ContractError(
                f"S{sector} JSON evidence {gate} has invalid sector"
            ) from exc
        if evidence_sector != sector:
            raise FM0ContractError(f"S{sector} JSON evidence {gate} has wrong sector")
        metrics = _evidence_metrics(
            gate,
            evidence.get("metrics"),
            sector=sector,
        )
        evidence_metric_values[gate] = metrics
        evidence_upstream_hashes[gate] = _verify_evidence_upstream_artifacts(
            evidence_path=evidence_path,
            gate=gate,
            sector=sector,
            raw=evidence.get("upstream_artifacts"),
            expected_metrics=metrics,
            minimum_count=int(admission["minimum_upstream_artifacts_per_gate"]),
            global_paths=upstream_paths,
            global_digests=upstream_digests,
        )
        evidence_hashes[gate] = evidence_hash

    alias_binding = _mapping(
        receipt.get(admission["alias_authority_binding_key"]),
        label=f"S{sector} alias authority binding",
    )
    if set(alias_binding) != {"path", "sha256"}:
        raise FM0ContractError(f"S{sector} alias authority binding schema drifted")
    alias_path, alias_hash = _resolve_binding(
        receipt_path,
        alias_binding,
        label=f"S{sector} alias authority",
    )
    alias_rows = _alias_edges(alias_path, config, label=f"S{sector} alias authority")
    local_registry = build_alias_registry(alias_rows)
    local_aliases = local_registry.alias_index()

    manifest_binding = _mapping(
        receipt.get(admission["observation_manifest_binding_key"]),
        label=f"S{sector} observation manifest binding",
    )
    if set(manifest_binding) != {"path", "sha256"}:
        raise FM0ContractError(f"S{sector} observation manifest binding schema drifted")
    manifest_path, manifest_hash = _resolve_binding(
        receipt_path,
        manifest_binding,
        label=f"S{sector} observation manifest",
    )
    rows, fields = _read_table(manifest_path, label=f"S{sector} observation manifest")
    _forbidden_columns(fields, config, label=f"S{sector} observation manifest")
    missing = sorted(_LATER_REQUIRED_FIELDS - set(fields))
    extra = sorted(set(fields) - _LATER_ALLOWED_FIELDS)
    if missing or extra:
        raise FM0ContractError(
            f"S{sector} observation schema differs; missing={missing}, extra={extra}"
        )
    normalized: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in rows:
        row = _normalize_later_row(raw, sector=sector, expected_orbits=expected_orbits)
        if row["observation_key"] in seen:
            raise FM0ContractError(f"S{sector} observation keys are duplicated")
        seen.add(row["observation_key"])
        alias = local_aliases.get((row["gaia_dr3_source_id"], row["tic_id"]))
        if alias is None:
            raise FM0ContractError(f"S{sector} observation edge lacks alias authority")
        if (
            row["leakage_component_id"] != alias["leakage_component_id"]
            or row["source_partition"] != alias["source_partition"]
        ):
            raise FM0ContractError(
                f"S{sector} observation alias/component assignment drifted"
            )
        normalized.append(row)
    return VerifiedSectorAdmission(
        sector=sector,
        receipt_sha256=sha256_file(receipt_path),
        stage1_acceptance_receipt_sha256=stage1_hash,
        observation_manifest_sha256=manifest_hash,
        alias_authority_sha256=alias_hash,
        evidence_sha256=dict(sorted(evidence_hashes.items())),
        evidence_metrics=dict(sorted(evidence_metric_values.items())),
        evidence_upstream_sha256=dict(sorted(evidence_upstream_hashes.items())),
        alias_rows=alias_rows,
        observation_rows=tuple(
            sorted(normalized, key=lambda row: row["observation_key"])
        ),
    )


def _validate_sector_lists(
    config: Mapping[str, Any],
    receipt_sectors: tuple[int, ...],
    selected: tuple[int, ...],
) -> None:
    selection = _mapping(config.get("selection"), label="selection policy")
    candidates = tuple(int(value) for value in selection["ordered_candidate_sectors"])
    if len(set(receipt_sectors)) != len(receipt_sectors):
        raise FM0ContractError("receipt sectors contain a duplicate")
    if tuple(sorted(selected)) != selected:
        raise FM0ContractError(
            "selected sectors are not a chronological admitted-sector prefix"
        )
    if not receipt_sectors or receipt_sectors != candidates[: len(receipt_sectors)]:
        raise FM0ContractError(
            "receipt sectors must be a nonempty candidate-sector prefix"
        )
    if not selected or selected != receipt_sectors[: len(selected)]:
        raise FM0ContractError(
            "selected sectors must be a nonempty admitted-sector prefix"
        )
    if any(right != left + 1 for left, right in pairwise(receipt_sectors)):
        raise FM0ContractError(
            "admitted sectors are not contiguous chronological successors"
        )
    configured = tuple(int(value) for value in selection.get("selected_sectors", ()))
    if configured and configured != selected:
        raise FM0ContractError("runtime selected sectors differ from frozen selection")
    first = int(selection["first_allowed_sector"])
    training = tuple(int(value) for value in config["baseline"]["training_sectors"])
    if receipt_sectors[0] != first or first != max(training) + 1:
        raise FM0ContractError(
            "later-sector inventory does not follow the training era"
        )
    for sector in receipt_sectors:
        _expected_orbits(config, sector)


def _combined_alias_registry(
    baseline: BaselineAuthority, admissions: Sequence[VerifiedSectorAdmission]
) -> AliasRegistry:
    edge_map: dict[tuple[str, str], dict[str, str]] = {}
    for row in baseline.alias_rows:
        edge_map[(row["gaia_dr3_source_id"], row["tic_id"])] = row
    baseline_gaia = {edge[0] for edge in edge_map}
    baseline_tic = {edge[1] for edge in edge_map if edge[1]}
    for admission in admissions:
        for row in admission.alias_rows:
            edge = (row["gaia_dr3_source_id"], row["tic_id"])
            if edge not in edge_map and (
                edge[0] in baseline_gaia or edge[1] in baseline_tic
            ):
                raise FM0ContractError(
                    f"new later alias edge touches frozen baseline graph: {edge}"
                )
            edge_map.setdefault(edge, row)
    combined = build_alias_registry(edge_map.values())
    old = baseline.alias_registry.alias_index()
    new = combined.alias_index()
    for edge, authority in old.items():
        observed = new.get(edge)
        if observed is None or any(
            observed[field] != authority[field]
            for field in ("leakage_component_id", "source_partition", "quarantined")
        ):
            raise FM0ContractError(f"later alias closure changed baseline edge {edge}")
    # Every sector authority must already represent its complete connected
    # closure; otherwise component IDs could drift as later sectors arrive.
    for admission in admissions:
        local = build_alias_registry(admission.alias_rows).alias_index()
        for edge, authority in local.items():
            observed = new.get(edge)
            if observed is None or any(
                observed[field] != authority[field]
                for field in ("leakage_component_id", "source_partition", "quarantined")
            ):
                raise FM0ContractError(
                    f"S{admission.sector} alias authority is not a stable full closure"
                )
    return combined


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fields})
    return stream.getvalue().encode("utf-8")


def construct_later_sector_inventory(
    *,
    config_path: str | Path,
    expected_config_sha256: str,
    ordered_sector_receipts: Sequence[tuple[int, str | Path, str]],
    selected_sectors: Sequence[int],
    baseline_manifest_path: str | Path,
    baseline_manifest_sha256: str,
    baseline_alias_authority_path: str | Path,
    baseline_alias_authority_sha256: str,
    baseline_selection_path: str | Path,
    baseline_selection_sha256: str,
) -> TemporalAdmissionInventory:
    """Verify all authorities and construct an identity-only later inventory."""

    config, config_hash = _load_config(
        config_path,
        expected_sha256=expected_config_sha256,
    )
    if any(len(item) != 3 for item in ordered_sector_receipts):
        raise FM0ContractError("sector receipt entries must be (sector, path, SHA-256)")
    receipt_sectors = tuple(int(item[0]) for item in ordered_sector_receipts)
    selected = tuple(int(value) for value in selected_sectors)
    _validate_sector_lists(config, receipt_sectors, selected)

    baseline_policy = _mapping(config.get("baseline"), label="baseline policy")
    bindings = (
        ("manifest_sha256", baseline_manifest_sha256),
        ("alias_authority_sha256", baseline_alias_authority_sha256),
        ("selection_sha256", baseline_selection_sha256),
    )
    for key, runtime_hash in bindings:
        if _digest(runtime_hash, label=key) != _digest(
            baseline_policy[key], label=f"baseline.{key}"
        ):
            raise FM0ContractError(f"runtime {key} differs from frozen config")
    manifest = _verify_file(
        baseline_manifest_path,
        baseline_manifest_sha256,
        label="baseline input manifest",
    )
    aliases = _verify_file(
        baseline_alias_authority_path,
        baseline_alias_authority_sha256,
        label="baseline alias authority",
    )
    selection = _verify_file(
        baseline_selection_path,
        baseline_selection_sha256,
        label="baseline corpus selection",
    )
    baseline = _load_baseline(
        config=config,
        manifest_path=manifest,
        alias_path=aliases,
        selection_path=selection,
    )

    admissions = tuple(
        verify_sector_admission_receipt(
            sector=int(sector),
            receipt_path=receipt_path,
            expected_receipt_sha256=receipt_hash,
            config=config,
        )
        for sector, receipt_path, receipt_hash in ordered_sector_receipts
    )
    combined = _combined_alias_registry(baseline, admissions)
    combined_aliases = combined.alias_index()
    selected_set = set(selected)
    all_rows = [row for admission in admissions for row in admission.observation_rows]
    later_keys = [row["observation_key"] for row in all_rows]
    if len(set(later_keys)) != len(later_keys):
        raise FM0ContractError("later observation keys are duplicated across sectors")
    later_products = [row["product_instance_id"] for row in all_rows]
    if len(set(later_products)) != len(later_products):
        raise FM0ContractError(
            "later product instance IDs are duplicated across sectors"
        )
    if set(later_keys) & baseline.observation_keys:
        raise FM0ContractError("later observation keys overlap the baseline release")
    for row in all_rows:
        alias = combined_aliases.get((row["gaia_dr3_source_id"], row["tic_id"]))
        if (
            alias is None
            or row["leakage_component_id"] != alias["leakage_component_id"]
            or row["source_partition"] != alias["source_partition"]
        ):
            raise FM0ContractError(
                "later observation differs from combined alias authority"
            )

    allowed = set(config["identity"]["allowed_source_partitions"])
    sealed = set(config["identity"]["sealed_source_partitions"])
    eligible: list[dict[str, Any]] = []
    quarantine: list[dict[str, Any]] = []
    sealed_counts = {sector: 0 for sector in receipt_sectors}
    quarantine_counts = {sector: 0 for sector in receipt_sectors}
    identity_counts = {sector: 0 for sector in receipt_sectors}
    for row in all_rows:
        sector = int(row["sector"])
        identity_counts[sector] += 1
        component = str(row["leakage_component_id"])
        partition = str(row["source_partition"])
        alias = combined_aliases[(row["gaia_dr3_source_id"], row["tic_id"])]
        if partition in sealed:
            sealed_counts[sector] += 1
            continue
        reason = ""
        if bool(alias["quarantined"]) or row["quarantined"]:
            reason = "source_registry_quarantined"
        elif partition not in allowed:
            raise FM0ContractError(f"later row has forbidden partition {partition!r}")
        elif not row["scientific_training_eligible"]:
            reason = "scientific_training_ineligible"
        elif not row["all_six_views_present"]:
            reason = "missing_required_six_view_channel"
        if reason:
            quarantine_counts[sector] += 1
            quarantine.append(
                {
                    "schema_version": QUARANTINE_SCHEMA_VERSION,
                    "observation_key": row["observation_key"],
                    "sector": sector,
                    "gaia_dr3_source_id": row["gaia_dr3_source_id"],
                    "tic_id": row["tic_id"],
                    "leakage_component_id": component,
                    "source_partition": partition,
                    "reason": reason,
                }
            )
            continue
        if sector not in selected_set:
            continue
        prior_visits = int(baseline.component_visit_counts.get(component, 0))
        cohort = "repeated_host" if prior_visits else "new_host"
        eligible.append(
            {
                "schema_version": INVENTORY_SCHEMA_VERSION,
                **{
                    field: row[field]
                    for field in (
                        "observation_key",
                        "product_instance_id",
                        "gaia_dr3_source_id",
                        "tic_id",
                        "sector",
                        "camera",
                        "ccd",
                        "orbits_json",
                        "leakage_component_id",
                        "source_partition",
                        "source_sha256",
                        "product_state",
                        "n_cadences",
                        "cadence_retention_fraction",
                        "flux_view_names_json",
                        "view_present_json",
                        "scientific_training_eligible",
                    )
                },
                "host_cohort": cohort,
                "cohort_partition": f"{cohort}_{partition}",
                "prior_training_era_visit_count": prior_visits,
            }
        )
    if not eligible:
        raise FM0ContractError("selected sectors contain no eligible nonsealed rows")
    eligible.sort(
        key=lambda row: (
            row["sector"],
            row["host_cohort"],
            row["leakage_component_id"],
            row["observation_key"],
        )
    )
    quarantine.sort(
        key=lambda row: (
            row["sector"],
            row["leakage_component_id"],
            row["observation_key"],
        )
    )

    repeated_components = {
        row["leakage_component_id"]
        for row in eligible
        if row["host_cohort"] == "repeated_host"
    }
    new_components = {
        row["leakage_component_id"]
        for row in eligible
        if row["host_cohort"] == "new_host"
    }
    detector_transitions: dict[str, int] = {}
    for row in eligible:
        if row["host_cohort"] != "repeated_host":
            continue
        destination = f"cam{row['camera']}/ccd{row['ccd']}"
        for origin in baseline.component_detector_visits[row["leakage_component_id"]]:
            key = f"{origin}->{destination}"
            detector_transitions[key] = detector_transitions.get(key, 0) + 1

    sector_rows: list[dict[str, Any]] = []
    by_sector_admission = {admission.sector: admission for admission in admissions}
    for sector in receipt_sectors:
        rows = [row for row in eligible if row["sector"] == sector]
        repeated = [row for row in rows if row["host_cohort"] == "repeated_host"]
        new = [row for row in rows if row["host_cohort"] == "new_host"]
        detector_cells = sorted({f"cam{row['camera']}/ccd{row['ccd']}" for row in rows})
        retentions = [float(row["cadence_retention_fraction"]) for row in rows]
        admission = by_sector_admission[sector]
        sector_rows.append(
            {
                "schema_version": SECTOR_INVENTORY_SCHEMA_VERSION,
                "sector": sector,
                "expected_orbits_json": json.dumps(
                    _expected_orbits(config, sector), separators=(",", ":")
                ),
                "admission_receipt_sha256": admission.receipt_sha256,
                "stage1_acceptance_receipt_sha256": (
                    admission.stage1_acceptance_receipt_sha256
                ),
                "observation_manifest_sha256": admission.observation_manifest_sha256,
                "alias_authority_sha256": admission.alias_authority_sha256,
                "evidence_sha256_json": json.dumps(
                    admission.evidence_sha256, sort_keys=True, separators=(",", ":")
                ),
                "evidence_metrics_json": json.dumps(
                    admission.evidence_metrics,
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "evidence_upstream_sha256_json": json.dumps(
                    admission.evidence_upstream_sha256,
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "n_identity_rows": identity_counts[sector],
                "n_eligible_rows": len(rows),
                "n_repeated_host_rows": len(repeated),
                "n_new_host_rows": len(new),
                "n_sealed_identity_rows_excluded": sealed_counts[sector],
                "n_quarantined_rows": quarantine_counts[sector],
                "n_repeated_host_components": len(
                    {row["leakage_component_id"] for row in repeated}
                ),
                "n_new_host_components": len(
                    {row["leakage_component_id"] for row in new}
                ),
                "n_camera_ccd_cells": len(detector_cells),
                "camera_ccd_coverage_json": json.dumps(
                    detector_cells,
                    separators=(",", ":"),
                ),
                "cadence_retention_fraction_min": (
                    min(retentions) if retentions else ""
                ),
                "cadence_retention_fraction_median": (
                    median(retentions) if retentions else ""
                ),
            }
        )

    minimum_sectors = int(config["selection"]["minimum_sectors"])
    count_floor_checks = {
        "minimum_selected_sectors": len(selected) >= minimum_sectors,
        "every_selected_sector_has_eligible_rows": all(
            any(row["sector"] == sector for row in eligible) for sector in selected
        ),
        "minimum_repeated_host_components": len(repeated_components)
        >= MINIMUM_REPEATED_HOST_COMPONENTS,
        "minimum_new_host_components": len(new_components)
        >= MINIMUM_NEW_HOST_COMPONENTS,
    }
    camera_ccd_cells = sorted(
        {f"cam{row['camera']}/ccd{row['ccd']}" for row in eligible}
    )
    new_sector_sets: dict[str, set[int]] = {}
    for row in eligible:
        if row["host_cohort"] == "new_host":
            new_sector_sets.setdefault(row["leakage_component_id"], set()).add(
                row["sector"]
            )
    summary = {
        "schema_version": INVENTORY_SCHEMA_VERSION,
        "campaign_id": str(config.get("campaign_id")),
        "config_sha256": config_hash,
        "baseline_manifest_sha256": sha256_file(manifest),
        "baseline_alias_authority_sha256": sha256_file(aliases),
        "baseline_selection_sha256": sha256_file(selection),
        "training_sectors": list(config["baseline"]["training_sectors"]),
        "admitted_sectors": list(receipt_sectors),
        "selected_development_candidate_sectors": list(selected),
        "n_identity_rows": len(all_rows),
        "n_eligible_observations": len(eligible),
        "n_repeated_host_observations": sum(
            row["host_cohort"] == "repeated_host" for row in eligible
        ),
        "n_new_host_observations": sum(
            row["host_cohort"] == "new_host" for row in eligible
        ),
        "n_repeated_host_components_with_prior_and_later_visits": len(
            repeated_components
        ),
        "n_new_host_components": len(new_components),
        "n_new_host_components_with_two_or_more_later_sectors": sum(
            len(value) >= 2 for value in new_sector_sets.values()
        ),
        "n_sealed_identity_rows_excluded": sum(sealed_counts.values()),
        "n_quarantined_rows": len(quarantine),
        "n_camera_ccd_cells": len(camera_ccd_cells),
        "camera_ccd_coverage_json": json.dumps(
            camera_ccd_cells,
            separators=(",", ":"),
        ),
        "eligible_n_cadences_min": min(int(row["n_cadences"]) for row in eligible),
        "eligible_n_cadences_max": max(int(row["n_cadences"]) for row in eligible),
        "eligible_cadence_retention_fraction_min": min(
            float(row["cadence_retention_fraction"]) for row in eligible
        ),
        "eligible_cadence_retention_fraction_max": max(
            float(row["cadence_retention_fraction"]) for row in eligible
        ),
        "repeated_host_detector_transition_visit_pairs": dict(
            sorted(detector_transitions.items())
        ),
        "count_floor_requirements": {
            "minimum_selected_sectors": minimum_sectors,
            "minimum_repeated_host_components": MINIMUM_REPEATED_HOST_COMPONENTS,
            "minimum_new_host_components": MINIMUM_NEW_HOST_COMPONENTS,
        },
        "count_floor_checks": count_floor_checks,
        "count_floor_ready": all(count_floor_checks.values()),
        "adequacy_thresholds_frozen": False,
        "panel_freeze_ready": False,
        "panel_frozen": False,
        "sealed_identity_emitted": False,
        "shards_or_light_curves_opened": False,
        "labels_search_scores_embeddings_or_injections_opened": False,
        "training_authorized": False,
        "sealed_test_access_authorized": False,
        "prospective_claim_authorized": False,
        "claim_limit": str(config.get("claim_limit", "")).strip(),
    }
    return TemporalAdmissionInventory(
        eligible_observations=tuple(eligible),
        sector_inventory=tuple(sector_rows),
        quarantine=tuple(quarantine),
        summary=summary,
    )


def write_later_sector_inventory(
    output_dir: str | Path, inventory: TemporalAdmissionInventory
) -> dict[str, Any]:
    output = Path(output_dir).expanduser().resolve()
    payloads = {
        "later_observations.csv": _csv_bytes(
            inventory.eligible_observations, LATER_OBSERVATION_FIELDS
        ),
        "sector_inventory.csv": _csv_bytes(
            inventory.sector_inventory, SECTOR_INVENTORY_FIELDS
        ),
        "quarantine.csv": _csv_bytes(inventory.quarantine, QUARANTINE_FIELDS),
    }
    summary = dict(inventory.summary)
    summary["outputs_sha256"] = {
        name: hashlib.sha256(payload).hexdigest() for name, payload in payloads.items()
    }
    summary_payload = (
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    for name, payload in payloads.items():
        publish_immutable(output / name, payload)
    publish_immutable(output / "summary.json", summary_payload)
    publish_immutable(
        output / "READY",
        (hashlib.sha256(summary_payload).hexdigest() + "  summary.json\n").encode(
            "utf-8"
        ),
    )
    return summary


def build_later_sector_inventory(
    *,
    config_path: str | Path,
    expected_config_sha256: str,
    ordered_sector_receipts: Sequence[tuple[int, str | Path, str]],
    selected_sectors: Sequence[int],
    baseline_manifest_path: str | Path,
    baseline_manifest_sha256: str,
    baseline_alias_authority_path: str | Path,
    baseline_alias_authority_sha256: str,
    baseline_selection_path: str | Path,
    baseline_selection_sha256: str,
    output_dir: str | Path,
) -> dict[str, Any]:
    inventory = construct_later_sector_inventory(
        config_path=config_path,
        expected_config_sha256=expected_config_sha256,
        ordered_sector_receipts=ordered_sector_receipts,
        selected_sectors=selected_sectors,
        baseline_manifest_path=baseline_manifest_path,
        baseline_manifest_sha256=baseline_manifest_sha256,
        baseline_alias_authority_path=baseline_alias_authority_path,
        baseline_alias_authority_sha256=baseline_alias_authority_sha256,
        baseline_selection_path=baseline_selection_path,
        baseline_selection_sha256=baseline_selection_sha256,
    )
    return write_later_sector_inventory(output_dir, inventory)


__all__ = [
    "INVENTORY_SCHEMA_VERSION",
    "LATER_OBSERVATION_FIELDS",
    "MINIMUM_NEW_HOST_COMPONENTS",
    "MINIMUM_REPEATED_HOST_COMPONENTS",
    "QUARANTINE_FIELDS",
    "SECTOR_INVENTORY_FIELDS",
    "TEMPORAL_SECTOR_RECEIPT_SCHEMA_VERSION",
    "BaselineAuthority",
    "TemporalAdmissionInventory",
    "VerifiedSectorAdmission",
    "build_later_sector_inventory",
    "construct_later_sector_inventory",
    "verify_sector_admission_receipt",
    "write_later_sector_inventory",
]
