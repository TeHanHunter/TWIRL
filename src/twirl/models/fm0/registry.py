"""Identity, split, and immutable-registry contract for TWIRL-FM0.1.

The scientific source is a Gaia DR3 source.  The independence unit used for
all source splits is the connected component of the Gaia--TIC alias graph.
Components containing more than one Gaia source are quarantined rather than
silently assigning aliases from the same physical ambiguity to different
partitions.

This module deliberately contains no BLS, candidate, label, or detector-level
features.  Registry products are audit/control artifacts and are never model
inputs.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import json
import math
from numbers import Integral, Real
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import yaml


REGISTRY_SCHEMA_VERSION = "twirl_fm0_1_registry_v1"
COMPONENT_SCHEMA_VERSION = "twirl_fm0_1_leakage_components_v1"
BUILD_SUMMARY_SCHEMA_VERSION = "twirl_fm0_1_registry_build_summary_v1"

DEFAULT_DESIGN_SHA256 = "94c8672a087884bf8a2c70f5d15315e05de602134af0c6c2073ca1c36232d6f7"
DEFAULT_CONFIG_SHA256 = "7de4e8c9e98c0ce27648a21241acc51766e436a537ba932b54f529fbdbf26d8a"
DEFAULT_FREEZE_RECEIPT_SHA256 = "75092235978c0b582f569be55770ad2d63368e079a6777da0b1c44547f074c25"

SPLIT_SALT = "twirl_fm0_1_source_split_v1"
DIAGNOSTIC_GATES = (
    "immutable_cell_receipts",
    "hdf5_openability",
    "cadence_authority_join",
    "external_quality_join",
    "fm_channel_numeric_gate",
)
DIAGNOSTIC_ADMISSION_RECEIPT_SCHEMA_VERSION = (
    "twirl_fm0_1_diagnostic_admission_receipt_v1"
)
ACCEPTED_SECTORS = frozenset(range(56, 66))
DIAGNOSTIC_SECTORS = frozenset((66, 67))

_DIGITS = re.compile(r"^[0-9]+$")
_FORBIDDEN_TOKENS = (
    "bls",
    "periodogram",
    "period",
    "epoch",
    "duration",
    "depth",
    "power",
    "sde",
    "peak",
    "score",
    "prediction",
    "pseudo",
    "candidate",
    "label",
    "teacher",
    "inject",
    "morphology",
    "fold",
    "event_window",
)


class FM0ContractError(ValueError):
    """Raised when a registry or release violates the frozen FM0.1 contract."""


@dataclass(frozen=True)
class FrozenContract:
    """Byte-bound frozen design and parsed scientific configuration."""

    design_path: Path
    config_path: Path
    freeze_receipt_path: Path
    design_sha256: str
    config_sha256: str
    freeze_receipt_sha256: str
    config: Mapping[str, Any]
    freeze_receipt: Mapping[str, Any]


@dataclass(frozen=True)
class AliasRegistry:
    """Canonical component and Gaia--TIC alias rows."""

    components: tuple[dict[str, Any], ...]
    aliases: tuple[dict[str, Any], ...]

    def alias_index(self) -> dict[tuple[str, str], dict[str, Any]]:
        return {
            (str(row["gaia_dr3_source_id"]), str(row["tic_id"])): row
            for row in self.aliases
        }


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_hash(path: Path, expected: str, label: str) -> str:
    if not path.is_file():
        raise FM0ContractError(f"missing {label}: {path}")
    actual = sha256_file(path)
    if actual != expected:
        raise FM0ContractError(
            f"{label} hash mismatch: expected {expected}, observed {actual}"
        )
    return actual


def load_frozen_contract(
    *,
    design_path: str | Path | None = None,
    config_path: str | Path | None = None,
    freeze_receipt_path: str | Path | None = None,
) -> FrozenContract:
    """Load and verify the byte-exact FM0.1 design authorities.

    The receipt is itself bound to a known hash, then its embedded design and
    config hashes are checked before the few invariants needed by this module
    are validated.  A changed authority requires a new named contract rather
    than a permissive warning.
    """

    root = _repo_root()
    design = Path(design_path or root / "doc/foundation_model_design.md").resolve()
    config = Path(
        config_path or root / "configs/models/twirl_fm0_1_s56_s67_poc.yaml"
    ).resolve()
    receipt = Path(
        freeze_receipt_path
        or root / "reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json"
    ).resolve()

    receipt_hash = _verify_hash(receipt, DEFAULT_FREEZE_RECEIPT_SHA256, "freeze receipt")
    receipt_data = json.loads(receipt.read_text(encoding="utf-8"))
    expected_design = str(receipt_data["design_document"]["sha256"])
    expected_config = str(receipt_data["scientific_config"]["sha256"])
    if expected_design != DEFAULT_DESIGN_SHA256 or expected_config != DEFAULT_CONFIG_SHA256:
        raise FM0ContractError("freeze receipt does not bind the expected FM0.1 authorities")
    design_hash = _verify_hash(design, expected_design, "design document")
    config_hash = _verify_hash(config, expected_config, "scientific config")
    config_data = yaml.safe_load(config.read_text(encoding="utf-8"))

    checks = (
        (config_data.get("campaign_id") == "twirl_fm0_1_s56_s67_poc", "campaign_id"),
        (config_data["identity"]["split_field"] == "leakage_component_id", "split field"),
        (config_data["source_split"]["salt"] == SPLIT_SALT, "split salt"),
        (int(config_data["source_split"]["modulus"]) == 10000, "split modulus"),
        (config_data["windowing"]["length_cadences"] == 2048, "window length"),
        (config_data["windowing"]["evaluation_stride_cadences"] == 1024, "eval stride"),
    )
    failed = [label for passed, label in checks if not passed]
    if failed:
        raise FM0ContractError("frozen config invariant mismatch: " + ", ".join(failed))
    return FrozenContract(
        design_path=design,
        config_path=config,
        freeze_receipt_path=receipt,
        design_sha256=design_hash,
        config_sha256=config_hash,
        freeze_receipt_sha256=receipt_hash,
        config=config_data,
        freeze_receipt=receipt_data,
    )


def _identifier(value: Any, *, name: str, optional: bool = False) -> str:
    """Normalize a positive integer identifier without lossy Gaia floats."""

    if value is None or (isinstance(value, str) and not value.strip()):
        if optional:
            return ""
        raise FM0ContractError(f"{name} is required")
    if isinstance(value, (bool, np.bool_)):
        raise FM0ContractError(f"{name} cannot be boolean")
    if isinstance(value, Integral):
        text = str(int(value))
    elif isinstance(value, Real):
        value = float(value)
        if not math.isfinite(value) or value != math.floor(value):
            raise FM0ContractError(f"{name} must be a positive integer")
        if name == "gaia_dr3_source_id" and abs(value) > 2**53:
            raise FM0ContractError(
                "Gaia identifiers above 2**53 must be supplied as integer/text, not float"
            )
        text = str(int(value))
    else:
        text = str(value).strip()
    if not _DIGITS.fullmatch(text) or int(text) <= 0:
        raise FM0ContractError(f"{name} must be a positive decimal integer")
    return str(int(text))


def _stable_id(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()}"


def deterministic_source_partition(leakage_component_id: str) -> tuple[str, int]:
    """Return the exact frozen 80/10/10 source partition and hash bucket."""

    component = str(leakage_component_id).strip()
    if not component:
        raise FM0ContractError("leakage_component_id is required")
    value = int.from_bytes(
        hashlib.sha256(f"{SPLIT_SALT}:{component}".encode("utf-8")).digest(),
        byteorder="big",
        signed=False,
    ) % 10000
    if value < 8000:
        return "poc_train", value
    if value < 9000:
        return "poc_development", value
    return "poc_sealed_test", value


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, node: str) -> None:
        self.parent.setdefault(node, node)

    def find(self, node: str) -> str:
        parent = self.parent[node]
        if parent != node:
            self.parent[node] = self.find(parent)
        return self.parent[node]

    def union(self, left: str, right: str) -> None:
        lroot, rroot = self.find(left), self.find(right)
        if lroot == rroot:
            return
        keep, move = sorted((lroot, rroot))
        self.parent[move] = keep


def assert_no_search_columns(columns: Iterable[str], *, context: str = "table") -> None:
    """Reject search/candidate/label fields before any registry tensor work."""

    bad = []
    for column in columns:
        normalized = str(column).strip().lower()
        if any(token in normalized for token in _FORBIDDEN_TOKENS):
            bad.append(str(column))
    if bad:
        raise FM0ContractError(f"{context} contains forbidden search/label fields: {sorted(bad)}")


def build_alias_registry(rows: Iterable[Mapping[str, Any]]) -> AliasRegistry:
    """Build canonical connected Gaia--TIC leakage components."""

    normalized: set[tuple[str, str]] = set()
    uf = _UnionFind()
    for raw in rows:
        assert_no_search_columns(raw.keys(), context="alias input")
        gaia = _identifier(raw.get("gaia_dr3_source_id"), name="gaia_dr3_source_id")
        tic = _identifier(raw.get("tic_id"), name="tic_id", optional=True)
        gnode = f"gaia:{gaia}"
        uf.add(gnode)
        if tic:
            tnode = f"tic:{tic}"
            uf.add(tnode)
            uf.union(gnode, tnode)
        normalized.add((gaia, tic))
    if not normalized:
        raise FM0ContractError("alias input is empty")

    nodes_by_root: dict[str, set[str]] = {}
    for node in uf.parent:
        nodes_by_root.setdefault(uf.find(node), set()).add(node)

    component_for_node: dict[str, dict[str, Any]] = {}
    components: list[dict[str, Any]] = []
    for nodes in sorted(nodes_by_root.values(), key=lambda values: sorted(values)):
        gaia_ids = sorted(node.split(":", 1)[1] for node in nodes if node.startswith("gaia:"))
        tic_ids = sorted(node.split(":", 1)[1] for node in nodes if node.startswith("tic:"))
        component_id = _stable_id(
            "leakage", {"gaia_dr3_source_ids": gaia_ids, "tic_ids": tic_ids}
        )
        quarantined = len(gaia_ids) != 1
        partition, bucket = (
            ("quarantine", -1)
            if quarantined
            else deterministic_source_partition(component_id)
        )
        component = {
            "component_schema_version": COMPONENT_SCHEMA_VERSION,
            "leakage_component_id": component_id,
            "gaia_dr3_source_ids_json": json.dumps(gaia_ids, separators=(",", ":")),
            "tic_ids_json": json.dumps(tic_ids, separators=(",", ":")),
            "n_gaia_sources": len(gaia_ids),
            "n_tic_aliases": len(tic_ids),
            "quarantined": quarantined,
            "quarantine_reason": "multiple_gaia_sources" if quarantined else "",
            "source_partition": partition,
            "split_bucket": bucket,
        }
        components.append(component)
        for node in nodes:
            component_for_node[node] = component

    aliases: list[dict[str, Any]] = []
    for gaia, tic in sorted(normalized, key=lambda pair: (int(pair[0]), int(pair[1] or 0))):
        component = component_for_node[f"gaia:{gaia}"]
        aliases.append(
            {
                "registry_schema_version": REGISTRY_SCHEMA_VERSION,
                "physical_source_id": f"gaia_dr3:{gaia}",
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "leakage_component_id": component["leakage_component_id"],
                "quarantined": component["quarantined"],
                "quarantine_reason": component["quarantine_reason"],
                "source_partition": component["source_partition"],
                "split_bucket": component["split_bucket"],
            }
        )
    components.sort(key=lambda row: row["leakage_component_id"])
    return AliasRegistry(tuple(components), tuple(aliases))


def _diagnostic_admission_receipt(
    raw: Mapping[str, Any], *, sector: int
) -> tuple[str, str]:
    """Verify the one immutable admission receipt required by S66--S67.

    Inline gate booleans are deliberately not an authority.  The registry row
    instead binds one JSON receipt by both path and SHA-256, and this function
    verifies the small portion of that receipt schema needed for FM admission.
    """

    legacy_fields = [
        f"{gate}_passed"
        for gate in DIAGNOSTIC_GATES
        if str(raw.get(f"{gate}_passed", "")).strip()
    ]
    if legacy_fields:
        raise FM0ContractError(
            "inline diagnostic gate values are forbidden; use one hash-bound "
            f"diagnostic admission receipt instead: {legacy_fields}"
        )

    path_text = str(raw.get("diagnostic_admission_receipt_path", "")).strip()
    expected_hash = str(
        raw.get("diagnostic_admission_receipt_sha256", "")
    ).strip().lower()
    if not path_text:
        raise FM0ContractError(
            f"S{sector} requires diagnostic_admission_receipt_path"
        )
    if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
        raise FM0ContractError(
            "diagnostic_admission_receipt_sha256 must be a lowercase SHA-256 digest"
        )
    receipt_path = Path(path_text).expanduser().resolve()
    if not receipt_path.is_file():
        raise FM0ContractError(f"missing diagnostic admission receipt: {receipt_path}")
    actual_hash = sha256_file(receipt_path)
    if actual_hash != expected_hash:
        raise FM0ContractError(
            "diagnostic admission receipt hash mismatch: "
            f"expected {expected_hash}, observed {actual_hash}"
        )
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(
            f"invalid diagnostic admission JSON receipt: {receipt_path}"
        ) from exc
    if not isinstance(receipt, dict):
        raise FM0ContractError("diagnostic admission receipt must be a JSON object")
    if (
        receipt.get("receipt_schema_version")
        != DIAGNOSTIC_ADMISSION_RECEIPT_SCHEMA_VERSION
    ):
        raise FM0ContractError("diagnostic admission receipt schema mismatch")
    try:
        receipt_sector = int(receipt.get("sector"))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError("diagnostic admission receipt has invalid sector") from exc
    if receipt_sector != sector:
        raise FM0ContractError(
            f"diagnostic admission receipt sector S{receipt_sector} does not match S{sector}"
        )
    if receipt.get("passed") is not True:
        raise FM0ContractError("diagnostic admission receipt is not passed")
    evidence = receipt.get("evidence_gates")
    if not isinstance(evidence, dict) or set(evidence) != set(DIAGNOSTIC_GATES):
        raise FM0ContractError(
            "diagnostic admission receipt must contain exactly the five evidence gates"
        )
    for gate in DIAGNOSTIC_GATES:
        binding = evidence[gate]
        if not isinstance(binding, dict):
            raise FM0ContractError(
                f"diagnostic evidence gate {gate!r} must be an artifact binding"
            )
        if binding.get("passed") is not True:
            raise FM0ContractError(
                f"diagnostic admission receipt failed evidence gates: {[gate]}"
            )
        evidence_path_text = str(binding.get("path", "")).strip()
        evidence_hash = str(binding.get("sha256", "")).strip().lower()
        if not evidence_path_text or not re.fullmatch(r"[0-9a-f]{64}", evidence_hash):
            raise FM0ContractError(
                f"diagnostic evidence gate {gate!r} lacks a valid path/SHA-256 binding"
            )
        evidence_path = Path(evidence_path_text).expanduser()
        if not evidence_path.is_absolute():
            evidence_path = receipt_path.parent / evidence_path
        evidence_path = evidence_path.resolve()
        if not evidence_path.is_file():
            raise FM0ContractError(
                f"missing diagnostic evidence artifact for {gate!r}: {evidence_path}"
            )
        actual_evidence_hash = sha256_file(evidence_path)
        if actual_evidence_hash != evidence_hash:
            raise FM0ContractError(
                f"diagnostic evidence artifact hash mismatch for {gate!r}: "
                f"expected {evidence_hash}, observed {actual_evidence_hash}"
            )
        if evidence_path.suffix.lower() == ".json":
            try:
                evidence_receipt = json.loads(evidence_path.read_text(encoding="utf-8"))
            except (OSError, UnicodeError, json.JSONDecodeError) as exc:
                raise FM0ContractError(
                    f"invalid JSON diagnostic evidence artifact for {gate!r}"
                ) from exc
            if not isinstance(evidence_receipt, dict):
                raise FM0ContractError(
                    f"JSON diagnostic evidence artifact for {gate!r} must be an object"
                )
            if evidence_receipt.get("passed") is not True:
                raise FM0ContractError(
                    f"JSON diagnostic evidence artifact for {gate!r} is not passed"
                )
            if "sector" in evidence_receipt:
                try:
                    evidence_sector = int(evidence_receipt["sector"])
                except (TypeError, ValueError) as exc:
                    raise FM0ContractError(
                        f"JSON diagnostic evidence artifact for {gate!r} has invalid sector"
                    ) from exc
                if evidence_sector != sector:
                    raise FM0ContractError(
                        f"JSON diagnostic evidence artifact for {gate!r} has sector "
                        f"S{evidence_sector}, expected S{sector}"
                    )
    return str(receipt_path), actual_hash


def build_observation_registry(
    rows: Iterable[Mapping[str, Any]], alias_registry: AliasRegistry
) -> tuple[dict[str, Any], ...]:
    """Bind selected A2v1 product instances to source-sector observations."""

    index = alias_registry.alias_index()
    observations: list[dict[str, Any]] = []
    seen_observations: set[str] = set()
    for raw in rows:
        assert_no_search_columns(raw.keys(), context="observation input")
        gaia = _identifier(raw.get("gaia_dr3_source_id"), name="gaia_dr3_source_id")
        tic = _identifier(raw.get("tic_id"), name="tic_id")
        alias = index.get((gaia, tic))
        if alias is None:
            raise FM0ContractError(f"unregistered Gaia--TIC edge: {gaia}--{tic}")
        sector = int(raw.get("sector"))
        state = str(raw.get("product_state", "")).strip()
        if sector in ACCEPTED_SECTORS:
            if state != "A2V1_ACCEPTED":
                raise FM0ContractError(
                    f"S{sector} requires product_state=A2V1_ACCEPTED, observed {state!r}"
                )
        elif sector in DIAGNOSTIC_SECTORS:
            if state != "ORCD_COMPLETE_DEFERRED":
                raise FM0ContractError(
                    f"S{sector} requires product_state=ORCD_COMPLETE_DEFERRED"
                )
            diagnostic_receipt_path, diagnostic_receipt_sha256 = (
                _diagnostic_admission_receipt(raw, sector=sector)
            )
        else:
            raise FM0ContractError(f"sector {sector} is outside frozen S56--S67 scope")
        if sector not in DIAGNOSTIC_SECTORS:
            diagnostic_receipt_path = ""
            diagnostic_receipt_sha256 = ""
        version = str(raw.get("a2v1_product_version", "")).strip()
        source_hash = str(raw.get("source_sha256", "")).strip().lower()
        if not version:
            raise FM0ContractError("a2v1_product_version is required")
        if not re.fullmatch(r"[0-9a-f]{64}", source_hash):
            raise FM0ContractError("source_sha256 must be a lowercase SHA-256 digest")
        physical = str(alias["physical_source_id"])
        observation_key = _stable_id(
            "observation", {"physical_source_id": physical, "sector": sector}
        )
        if observation_key in seen_observations:
            raise FM0ContractError(
                f"more than one selected product for physical source {physical} in S{sector}"
            )
        seen_observations.add(observation_key)
        product_instance_id = _stable_id(
            "product",
            {
                "observation_key": observation_key,
                "tic_id": tic,
                "a2v1_product_version": version,
                "source_sha256": source_hash,
            },
        )
        observations.append(
            {
                "registry_schema_version": REGISTRY_SCHEMA_VERSION,
                "observation_key": observation_key,
                "product_instance_id": product_instance_id,
                "physical_source_id": physical,
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "sector": sector,
                "a2v1_product_version": version,
                "source_sha256": source_hash,
                "product_state": state,
                "diagnostic_admission_receipt_path": diagnostic_receipt_path,
                "diagnostic_admission_receipt_sha256": diagnostic_receipt_sha256,
                "leakage_component_id": alias["leakage_component_id"],
                "source_partition": alias["source_partition"],
                "quarantined": alias["quarantined"],
            }
        )
    observations.sort(key=lambda row: row["observation_key"])
    return tuple(observations)


def read_rows(path: str | Path) -> list[dict[str, Any]]:
    """Read a small registry input table from CSV or a JSON row list."""

    source = Path(path)
    if source.suffix.lower() == ".json":
        value = json.loads(source.read_text(encoding="utf-8"))
        if not isinstance(value, list) or not all(isinstance(row, dict) for row in value):
            raise FM0ContractError("JSON registry input must be a list of objects")
        return value
    with source.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fields})
    return stream.getvalue().encode("utf-8")


def publish_immutable(path: str | Path, payload: bytes) -> None:
    """Create an immutable artifact; identical retries are idempotent."""

    destination = Path(path)
    if destination.exists():
        if destination.read_bytes() != payload:
            raise FM0ContractError(f"refusing to replace non-identical artifact: {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(payload)


def write_registry_release(
    out_dir: str | Path,
    registry: AliasRegistry,
    observations: Sequence[Mapping[str, Any]] = (),
    *,
    contract: FrozenContract | None = None,
) -> dict[str, Any]:
    """Publish a deterministic local registry release and build summary."""

    contract = contract or load_frozen_contract()
    output = Path(out_dir)
    component_fields = tuple(registry.components[0])
    alias_fields = tuple(registry.aliases[0])
    payloads: dict[str, bytes] = {
        "components.csv": _csv_bytes(registry.components, component_fields),
        "aliases.csv": _csv_bytes(registry.aliases, alias_fields),
    }
    if observations:
        payloads["observations.csv"] = _csv_bytes(observations, tuple(observations[0]))
    hashes = {name: hashlib.sha256(payload).hexdigest() for name, payload in payloads.items()}
    summary = {
        "summary_schema_version": BUILD_SUMMARY_SCHEMA_VERSION,
        "campaign_id": contract.config["campaign_id"],
        "design_sha256": contract.design_sha256,
        "config_sha256": contract.config_sha256,
        "freeze_receipt_sha256": contract.freeze_receipt_sha256,
        "n_components": len(registry.components),
        "n_alias_rows": len(registry.aliases),
        "n_observations": len(observations),
        "n_quarantined_components": sum(bool(row["quarantined"]) for row in registry.components),
        "outputs_sha256": hashes,
        "certifies_full_campaign": False,
    }
    payloads["summary.json"] = (
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    for name, payload in payloads.items():
        publish_immutable(output / name, payload)
    return summary
