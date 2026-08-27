from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
import yaml

from twirl.models.fm0.corpus import (
    CORPUS_SELECTION_FIELDS,
    CORPUS_SELECTION_SCHEMA_VERSION,
)
from twirl.models.fm0.input_release import (
    FLUX_VIEW_NAMES,
    INPUT_RELEASE_SCHEMA_VERSION,
    MANIFEST_COLUMNS,
)
from twirl.models.fm0.registry import (
    FM0ContractError,
    build_alias_registry,
    deterministic_source_partition,
)
from twirl.models.fm0.temporal_admission import (
    TEMPORAL_SECTOR_RECEIPT_SCHEMA_VERSION,
    build_later_sector_inventory,
    construct_later_sector_inventory,
)

_POLICY_PATH = (
    Path(__file__).resolve().parents[1]
    / "configs/models/twirl_fm0_2_later_sector_admission_v1.yaml"
)
_EVIDENCE_SCHEMA_VERSION = "twirl_fm0_later_sector_evidence_v1"
_EVIDENCE_METRIC_CONTRACT_VERSION = "twirl_fm0_later_sector_evidence_metrics_v1"
_LATER_OBSERVATION_SCHEMA_VERSION = "twirl_fm0_later_sector_observation_v1"
_TEST_CAMPAIGN_ID = "twirl_fm0_later_sector_admission_test_fixture_v1"
_UPSTREAM_CONTROL_SCHEMA_VERSION = "twirl_fm0_later_sector_upstream_control_v1"


@dataclass(frozen=True)
class _Host:
    gaia: str
    tics: tuple[str, ...]
    component: str
    partition: str


@dataclass
class _Bundle:
    config_path: Path
    config: dict[str, Any]
    baseline_manifest_path: Path
    baseline_manifest_sha256: str
    baseline_alias_path: Path
    baseline_alias_sha256: str
    baseline_selection_path: Path
    baseline_selection_sha256: str
    receipt_tuples: list[tuple[int, Path, str]]
    receipt_paths: dict[int, Path]
    stage1_acceptance_paths: dict[int, Path]
    manifest_paths: dict[int, Path]
    alias_paths: dict[int, Path]
    selected_sectors: list[int]
    repeated_hosts: list[_Host]
    new_hosts: list[_Host]
    sealed_host: _Host


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str] | None = None,
) -> None:
    if fields is None:
        if not rows:
            raise AssertionError("test CSV rows must not be empty")
        fields = tuple(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {field: row.get(field, "") for field in fields} for row in rows
        )


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _valid_evidence_metrics(gate: str) -> dict[str, int]:
    if gate in {
        "checksum_bound_a2v1_hdf5_provenance",
        "checksum_bound_a2v1_fits_provenance",
    }:
        return {
            "n_products": 100,
            "n_checksum_bound_products": 100,
            "n_checksum_mismatches": 0,
        }
    if gate == "edge_aware_coverage":
        return {
            "n_expected_observations": 100,
            "n_present_observations": 98,
            "n_edge_omissions": 2,
            "n_non_edge_omissions": 0,
        }
    if gate == "hdf5_openability":
        return {
            "n_hdf5_products": 100,
            "n_hdf5_opened": 100,
            "n_unreadable_hdf5": 0,
        }
    if gate in {
        "authoritative_internal_cadence_quality",
        "authoritative_external_cadence_quality",
    }:
        return {
            "n_observations_checked": 100,
            "n_cadences_checked": 409_600,
            "n_observations_failed": 0,
        }
    if gate == "explicit_omissions":
        return {
            "n_expected_observations": 100,
            "n_present_observations": 98,
            "n_explicit_omissions": 2,
            "n_unexplained_omissions": 0,
        }
    if gate == "fm_channel_mask_finite_numerical_envelope":
        return {
            "n_observations_checked": 100,
            "n_six_view_observations": 100,
            "n_observations_failed": 0,
            "n_nonfinite_active_values": 0,
        }
    if gate == "stable_physical_source_registry_join":
        return {
            "n_identity_rows": 100,
            "n_joined_identity_rows": 98,
            "n_quarantined_identity_rows": 2,
            "n_unmatched_identity_rows": 0,
        }
    raise AssertionError(f"unknown evidence gate: {gate}")


def _find_host(
    *,
    namespace: int,
    partitions: set[str],
    n_aliases: int = 1,
) -> _Host:
    for offset in range(10_000):
        gaia = str(namespace * 10_000_000 + offset + 1)
        tics = tuple(
            str(namespace * 100_000_000 + (offset + 1) * 10 + index + 1)
            for index in range(n_aliases)
        )
        registry = build_alias_registry(
            [{"gaia_dr3_source_id": gaia, "tic_id": tic} for tic in tics]
        )
        component = registry.components[0]
        partition = str(component["source_partition"])
        if not component["quarantined"] and partition in partitions:
            return _Host(
                gaia=gaia,
                tics=tics,
                component=str(component["leakage_component_id"]),
                partition=partition,
            )
    raise AssertionError(f"could not find deterministic host in {partitions}")


def _make_hosts(
    count: int,
    *,
    namespace: int,
    partitions: set[str],
    alternate_first: bool = False,
) -> list[_Host]:
    return [
        _find_host(
            namespace=namespace + index,
            partitions=partitions,
            n_aliases=2 if alternate_first and index == 0 else 1,
        )
        for index in range(count)
    ]


def _selection_row(host: _Host, *, sector: int, visit: int) -> dict[str, Any]:
    first_orbit = 2 * sector + 7
    tic = host.tics[0]
    return {
        "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
        "gaia_dr3_source_id": host.gaia,
        "tic_id": tic,
        "sector": sector,
        "camera": visit % 4 + 1,
        "ccd": (visit // 4) % 4 + 1,
        "leakage_component_id": host.component,
        "source_partition": host.partition,
        "orbits_json": json.dumps(
            [first_orbit, first_orbit + 1], separators=(",", ":")
        ),
        "hdf5_paths_json": json.dumps(
            [
                f"orbit-{first_orbit}/ffi/cam1/ccd1/LC/{tic}.h5",
                f"orbit-{first_orbit + 1}/ffi/cam1/ccd1/LC/{tic}.h5",
            ],
            separators=(",", ":"),
        ),
    }


def _baseline_manifest_row(
    selection: Mapping[str, Any], *, index: int
) -> dict[str, Any]:
    key = f"baseline_observation_{index:06d}"
    return {
        "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
        "observation_key": key,
        "product_instance_id": f"baseline_product_{index:06d}",
        "source_sha256": "a" * 64,
        "leakage_component_id": selection["leakage_component_id"],
        "source_partition": selection["source_partition"],
        "product_state": "A2V1_ACCEPTED",
        "relative_path": f"shards/{key}.npz",
        "sha256": "b" * 64,
        "input_source_sha256": "c" * 64,
        "n_cadences": 4096,
        "n_segments": 2,
        "view_present_json": "[1,1,1,1,1,1]",
        "host_visit_offset_cadences": index * 5000,
        "host_visit_gap_cadences": 1000,
        "host_visit_overlaps_previous": False,
        "input_adapter": "a2v1_hdf5_v1",
        "scientific_training_eligible": True,
    }


def _later_manifest_row(
    host: _Host,
    *,
    sector: int,
    role: str,
    index: int,
    tic: str | None = None,
) -> dict[str, Any]:
    first_orbit = 2 * sector + 7
    return {
        "schema_version": _LATER_OBSERVATION_SCHEMA_VERSION,
        "observation_key": f"later_{sector}_{role}_{index:06d}",
        "product_instance_id": f"later_product_{sector}_{role}_{index:06d}",
        "gaia_dr3_source_id": host.gaia,
        "tic_id": tic or host.tics[0],
        "sector": sector,
        "camera": index % 4 + 1,
        "ccd": (index // 4) % 4 + 1,
        "orbits_json": json.dumps(
            [first_orbit, first_orbit + 1], separators=(",", ":")
        ),
        "source_sha256": "d" * 64,
        "product_state": "A2V1_ACCEPTED",
        "leakage_component_id": host.component,
        "source_partition": host.partition,
        "quarantined": False,
        "n_cadences": 4096,
        "cadence_retention_fraction": "0.99",
        "flux_view_names_json": json.dumps(FLUX_VIEW_NAMES, separators=(",", ":")),
        "view_present_json": "[1,1,1,1,1,1]",
        "scientific_training_eligible": True,
    }


def _refresh_config_hashes(bundle: _Bundle) -> None:
    baseline = bundle.config["baseline"]
    baseline["manifest_sha256"] = bundle.baseline_manifest_sha256
    baseline["alias_authority_sha256"] = bundle.baseline_alias_sha256
    baseline["selection_sha256"] = bundle.baseline_selection_sha256
    baseline["selection_n_observations"] = len(
        _read_csv(bundle.baseline_selection_path)
    )
    bundle.config_path.write_text(
        yaml.safe_dump(bundle.config, sort_keys=False), encoding="utf-8"
    )


def _make_bundle(
    tmp_path: Path,
    *,
    sectors: Sequence[int] = (65, 66, 67, 68, 69),
    selected_sectors: Sequence[int] | None = None,
    n_repeated: int = 1,
    n_new: int = 1,
) -> _Bundle:
    repeated_hosts = _make_hosts(
        n_repeated,
        namespace=10,
        partitions={"poc_train", "poc_development"},
        alternate_first=True,
    )
    new_hosts = _make_hosts(
        n_new,
        namespace=1_000,
        partitions={"poc_train", "poc_development"},
    )
    sealed_host = _find_host(
        namespace=9_000,
        partitions={"poc_sealed_test"},
    )

    baseline_alias_rows = [
        {"gaia_dr3_source_id": host.gaia, "tic_id": tic}
        for host in repeated_hosts
        for tic in host.tics
    ]
    baseline_alias_path = tmp_path / "baseline" / "aliases.csv"
    _write_csv(
        baseline_alias_path,
        baseline_alias_rows,
        ("gaia_dr3_source_id", "tic_id"),
    )

    baseline_sectors = list(range(56, 65))
    baseline_selection: list[dict[str, Any]] = []
    # The first component spans the complete era; each additional component has
    # one accepted visit.  This gives the manifest/selection transition audit a
    # deterministic per-component visit-count authority.
    for visit, sector in enumerate(baseline_sectors):
        baseline_selection.append(
            _selection_row(repeated_hosts[0], sector=sector, visit=visit)
        )
    for index, host in enumerate(repeated_hosts[1:], start=1):
        baseline_selection.append(
            _selection_row(
                host,
                sector=baseline_sectors[index % len(baseline_sectors)],
                visit=index + len(baseline_sectors),
            )
        )
    baseline_selection_path = tmp_path / "baseline" / "selection.csv"
    _write_csv(
        baseline_selection_path,
        baseline_selection,
        CORPUS_SELECTION_FIELDS,
    )

    baseline_manifest = [
        _baseline_manifest_row(row, index=index)
        for index, row in enumerate(baseline_selection)
    ]
    baseline_manifest_path = tmp_path / "baseline" / "manifest.csv"
    _write_csv(baseline_manifest_path, baseline_manifest, MANIFEST_COLUMNS)

    config = yaml.safe_load(_POLICY_PATH.read_text(encoding="utf-8"))
    config["campaign_id"] = _TEST_CAMPAIGN_ID
    config["selection"]["ordered_candidate_sectors"] = list(sectors)
    config["selection"]["selected_sectors"] = []
    config["selection"]["expected_orbits_by_sector"] = {
        int(sector): [2 * int(sector) + 7, 2 * int(sector) + 8] for sector in sectors
    }
    config_path = tmp_path / "policy.yaml"
    bundle = _Bundle(
        config_path=config_path,
        config=config,
        baseline_manifest_path=baseline_manifest_path,
        baseline_manifest_sha256=_sha256(baseline_manifest_path),
        baseline_alias_path=baseline_alias_path,
        baseline_alias_sha256=_sha256(baseline_alias_path),
        baseline_selection_path=baseline_selection_path,
        baseline_selection_sha256=_sha256(baseline_selection_path),
        receipt_tuples=[],
        receipt_paths={},
        stage1_acceptance_paths={},
        manifest_paths={},
        alias_paths={},
        selected_sectors=list(selected_sectors or sectors),
        repeated_hosts=repeated_hosts,
        new_hosts=new_hosts,
        sealed_host=sealed_host,
    )
    _refresh_config_hashes(bundle)

    gates = tuple(
        str(value) for value in config["admission"]["required_evidence_gates"]
    )
    for sector in sectors:
        sector = int(sector)
        receipt_dir = tmp_path / f"s{sector:04d}"
        manifest_rows: list[dict[str, Any]] = []
        for index, host in enumerate(repeated_hosts):
            # Exercise a frozen alternate alias for the first repeated host.
            tic = host.tics[-1] if index == 0 else host.tics[0]
            manifest_rows.append(
                _later_manifest_row(
                    host,
                    sector=sector,
                    role="repeated",
                    index=index,
                    tic=tic,
                )
            )
        for index, host in enumerate(new_hosts):
            manifest_rows.append(
                _later_manifest_row(
                    host,
                    sector=sector,
                    role="new",
                    index=index,
                )
            )
        manifest_rows.append(
            _later_manifest_row(
                sealed_host,
                sector=sector,
                role="sealed",
                index=0,
            )
        )
        manifest_path = receipt_dir / "observations.csv"
        _write_csv(manifest_path, manifest_rows)

        # This is the complete edge closure needed to reconstruct every
        # component asserted by the sector manifest.  It intentionally has no
        # derived component or partition columns.
        alias_rows = [
            {"gaia_dr3_source_id": host.gaia, "tic_id": tic}
            for host in (*repeated_hosts, *new_hosts, sealed_host)
            for tic in host.tics
        ]
        alias_path = receipt_dir / "aliases.csv"
        _write_csv(alias_path, alias_rows, ("gaia_dr3_source_id", "tic_id"))

        evidence: dict[str, dict[str, Any]] = {}
        for gate in gates:
            evidence_path = receipt_dir / "evidence" / f"{gate}.json"
            upstream_path = (
                evidence_path.parent / "upstream" / f"{gate}_control_metadata.json"
            )
            producer_dir = upstream_path.parent / "producer"
            producer_code_path = producer_dir / f"{gate}_producer.py"
            producer_config_path = producer_dir / f"{gate}_producer.yaml"
            producer_dir.mkdir(parents=True, exist_ok=True)
            producer_code_path.write_text(
                f"# synthetic producer for {gate}\n",
                encoding="utf-8",
            )
            producer_config_path.write_text(
                f"gate_name: {gate}\nsector: {sector}\n",
                encoding="utf-8",
            )
            metrics = _valid_evidence_metrics(gate)
            _write_json(
                upstream_path,
                {
                    "schema_version": _UPSTREAM_CONTROL_SCHEMA_VERSION,
                    "gate_name": gate,
                    "sector": sector,
                    "product_state": "A2V1_ACCEPTED",
                    "passed": True,
                    "metric_contract_version": _EVIDENCE_METRIC_CONTRACT_VERSION,
                    "metrics": metrics,
                    "producer": {
                        "git_commit": "a" * 40,
                        "code": {
                            "path": str(
                                producer_code_path.relative_to(upstream_path.parent)
                            ),
                            "sha256": _sha256(producer_code_path),
                        },
                        "config": {
                            "path": str(
                                producer_config_path.relative_to(upstream_path.parent)
                            ),
                            "sha256": _sha256(producer_config_path),
                        },
                    },
                },
            )
            _write_json(
                evidence_path,
                {
                    "schema_version": _EVIDENCE_SCHEMA_VERSION,
                    "gate_name": gate,
                    "sector": sector,
                    "product_state": "A2V1_ACCEPTED",
                    "passed": True,
                    "metrics": metrics,
                    "upstream_artifacts": [
                        {
                            "role": f"{gate}_control_metadata",
                            "path": str(
                                upstream_path.relative_to(evidence_path.parent)
                            ),
                            "sha256": _sha256(upstream_path),
                            "content_class": "control_metadata",
                        }
                    ],
                },
            )
            evidence[gate] = {
                "path": str(evidence_path.relative_to(receipt_dir)),
                "sha256": _sha256(evidence_path),
                "passed": True,
            }
        stage1_acceptance_path = receipt_dir / "stage1_acceptance.json"
        _write_json(
            stage1_acceptance_path,
            {
                "sector": sector,
                "orbits": [2 * sector + 7, 2 * sector + 8],
                "ok": True,
                "ok_h5": True,
                "ok_fits": True,
                "expected_contract": {
                    "ok": True,
                    "has_expected_rows": True,
                    "has_expected_unique_tics": True,
                    "missing_requested_orbits": [],
                    "requested_orbits": [2 * sector + 7, 2 * sector + 8],
                    "observed_orbits": [2 * sector + 7, 2 * sector + 8],
                },
                "h5": {
                    "n_present_h5": 100,
                    "n_missing_h5_non_edge": 0,
                    "n_unreadable_h5": 0,
                    "n_zero_byte_h5": 0,
                },
                "fits": {
                    "n_fits": 100,
                    "n_checked_fits": 100,
                    "n_bad_checked_fits": 0,
                    "n_missing_fits_non_edge_tics": 0,
                },
                "schema": {
                    "schema_only": True,
                    "check_h5_open": True,
                    "expected_method": "A2v1",
                    "expected_prodtag": "A2v1",
                },
            },
        )
        receipt_path = receipt_dir / "admission.json"
        _write_json(
            receipt_path,
            {
                "receipt_schema_version": TEMPORAL_SECTOR_RECEIPT_SCHEMA_VERSION,
                "sector": sector,
                "product_state": "A2V1_ACCEPTED",
                "passed": True,
                "observation_manifest": {
                    "path": manifest_path.name,
                    "sha256": _sha256(manifest_path),
                },
                "alias_authority": {
                    "path": alias_path.name,
                    "sha256": _sha256(alias_path),
                },
                "stage1_acceptance_receipt": {
                    "path": stage1_acceptance_path.name,
                    "sha256": _sha256(stage1_acceptance_path),
                },
                "evidence_gates": evidence,
            },
        )
        bundle.receipt_paths[sector] = receipt_path
        bundle.stage1_acceptance_paths[sector] = stage1_acceptance_path
        bundle.manifest_paths[sector] = manifest_path
        bundle.alias_paths[sector] = alias_path
        bundle.receipt_tuples.append((sector, receipt_path, _sha256(receipt_path)))
    return bundle


def _replace_receipt_digest(bundle: _Bundle, sector: int) -> None:
    path = bundle.receipt_paths[sector]
    bundle.receipt_tuples = [
        (item_sector, item_path, _sha256(path) if item_sector == sector else digest)
        for item_sector, item_path, digest in bundle.receipt_tuples
    ]


def _mutate_receipt(
    bundle: _Bundle,
    sector: int,
    mutation: Callable[[dict[str, Any]], None],
) -> None:
    path = bundle.receipt_paths[sector]
    payload = json.loads(path.read_text(encoding="utf-8"))
    mutation(payload)
    _write_json(path, payload)
    _replace_receipt_digest(bundle, sector)


def _refresh_receipt_binding(bundle: _Bundle, sector: int, key: str) -> None:
    bound_paths = {
        "observation_manifest": bundle.manifest_paths[sector],
        "alias_authority": bundle.alias_paths[sector],
        "stage1_acceptance_receipt": bundle.stage1_acceptance_paths[sector],
    }
    bound_path = bound_paths[key]

    def update(payload: dict[str, Any]) -> None:
        payload[key]["sha256"] = _sha256(bound_path)

    _mutate_receipt(bundle, sector, update)


def _mutate_evidence(
    bundle: _Bundle,
    sector: int,
    gate: str,
    mutation: Callable[[dict[str, Any]], None],
) -> None:
    receipt_path = bundle.receipt_paths[sector]
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    evidence_path = receipt_path.parent / receipt["evidence_gates"][gate]["path"]
    payload = json.loads(evidence_path.read_text(encoding="utf-8"))
    mutation(payload)
    _write_json(evidence_path, payload)

    def update(receipt_payload: dict[str, Any]) -> None:
        receipt_payload["evidence_gates"][gate]["sha256"] = _sha256(evidence_path)

    _mutate_receipt(bundle, sector, update)


def _mutate_upstream_control(
    bundle: _Bundle,
    sector: int,
    gate: str,
    mutation: Callable[[dict[str, Any]], None],
) -> None:
    receipt_path = bundle.receipt_paths[sector]
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    evidence_path = receipt_path.parent / receipt["evidence_gates"][gate]["path"]
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    upstream_binding = evidence["upstream_artifacts"][0]
    upstream_path = evidence_path.parent / upstream_binding["path"]
    payload = json.loads(upstream_path.read_text(encoding="utf-8"))
    mutation(payload)
    _write_json(upstream_path, payload)
    upstream_binding["sha256"] = _sha256(upstream_path)
    _write_json(evidence_path, evidence)

    def update(receipt_payload: dict[str, Any]) -> None:
        receipt_payload["evidence_gates"][gate]["sha256"] = _sha256(evidence_path)

    _mutate_receipt(bundle, sector, update)


def _rewrite_producer_artifact(
    bundle: _Bundle,
    sector: int,
    gate: str,
    artifact: str,
    contents: str,
) -> None:
    receipt_path = bundle.receipt_paths[sector]
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    evidence_path = receipt_path.parent / receipt["evidence_gates"][gate]["path"]
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    upstream_binding = evidence["upstream_artifacts"][0]
    upstream_path = evidence_path.parent / upstream_binding["path"]
    upstream = json.loads(upstream_path.read_text(encoding="utf-8"))
    producer_binding = upstream["producer"][artifact]
    artifact_path = upstream_path.parent / producer_binding["path"]
    artifact_path.write_text(contents, encoding="utf-8")
    producer_binding["sha256"] = _sha256(artifact_path)
    _write_json(upstream_path, upstream)
    upstream_binding["sha256"] = _sha256(upstream_path)
    _write_json(evidence_path, evidence)

    def update(receipt_payload: dict[str, Any]) -> None:
        receipt_payload["evidence_gates"][gate]["sha256"] = _sha256(evidence_path)

    _mutate_receipt(bundle, sector, update)


def _mutate_stage1_acceptance(
    bundle: _Bundle,
    sector: int,
    mutation: Callable[[dict[str, Any]], None],
) -> None:
    path = bundle.stage1_acceptance_paths[sector]
    payload = json.loads(path.read_text(encoding="utf-8"))
    mutation(payload)
    _write_json(path, payload)
    _refresh_receipt_binding(bundle, sector, "stage1_acceptance_receipt")


def _build_args(bundle: _Bundle) -> dict[str, Any]:
    return {
        "config_path": bundle.config_path,
        "expected_config_sha256": _sha256(bundle.config_path),
        "ordered_sector_receipts": bundle.receipt_tuples,
        "selected_sectors": bundle.selected_sectors,
        "baseline_manifest_path": bundle.baseline_manifest_path,
        "baseline_manifest_sha256": bundle.baseline_manifest_sha256,
        "baseline_alias_authority_path": bundle.baseline_alias_path,
        "baseline_alias_authority_sha256": bundle.baseline_alias_sha256,
        "baseline_selection_path": bundle.baseline_selection_path,
        "baseline_selection_sha256": bundle.baseline_selection_sha256,
    }


def test_valid_five_sector_inventory_meets_count_floors_but_is_not_freeze_ready(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path, n_repeated=64, n_new=256)
    output = tmp_path / "inventory"
    summary = build_later_sector_inventory(
        **_build_args(bundle),
        output_dir=output,
    )

    assert summary["count_floor_ready"] is True
    assert summary["adequacy_thresholds_frozen"] is False
    assert summary["panel_freeze_ready"] is False
    assert summary["panel_frozen"] is False
    assert len(summary["selected_development_candidate_sectors"]) == 5
    assert summary["n_repeated_host_components_with_prior_and_later_visits"] == 64
    assert summary["n_new_host_components"] == 256
    assert summary["n_sealed_identity_rows_excluded"] == 5
    assert summary["sealed_identity_emitted"] is False
    assert summary["shards_or_light_curves_opened"] is False

    rows = _read_csv(output / "later_observations.csv")
    repeated = [row for row in rows if row["host_cohort"] == "repeated_host"]
    new = [row for row in rows if row["host_cohort"] == "new_host"]
    assert repeated and new
    assert {row["source_partition"] for row in rows} <= {
        "poc_train",
        "poc_development",
    }
    assert not any(
        row["leakage_component_id"] == bundle.sealed_host.component for row in rows
    )
    first = bundle.repeated_hosts[0]
    assert {
        row["tic_id"] for row in repeated if row["gaia_dr3_source_id"] == first.gaia
    } == {first.tics[-1]}


def test_five_sectors_below_cohort_floors_are_not_freeze_ready(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    inventory = construct_later_sector_inventory(**_build_args(bundle))

    assert len(inventory.summary["selected_development_candidate_sectors"]) == 5
    assert inventory.summary["count_floor_checks"]["minimum_selected_sectors"]
    assert not inventory.summary["count_floor_checks"][
        "minimum_repeated_host_components"
    ]
    assert not inventory.summary["count_floor_checks"]["minimum_new_host_components"]
    assert inventory.summary["count_floor_ready"] is False
    assert inventory.summary["adequacy_thresholds_frozen"] is False
    assert inventory.summary["panel_freeze_ready"] is False
    assert inventory.summary["panel_frozen"] is False


def test_fewer_than_five_sectors_are_inventory_only(tmp_path: Path) -> None:
    sectors = (65, 66, 67, 68)
    bundle = _make_bundle(tmp_path, sectors=sectors)
    inventory = construct_later_sector_inventory(**_build_args(bundle))

    assert len(inventory.summary["selected_development_candidate_sectors"]) == 4
    assert not inventory.summary["count_floor_checks"]["minimum_selected_sectors"]
    assert inventory.summary["count_floor_ready"] is False
    assert inventory.summary["adequacy_thresholds_frozen"] is False
    assert inventory.summary["panel_freeze_ready"] is False
    assert inventory.summary["panel_frozen"] is False


def test_expected_config_sha256_is_required_and_rejects_drift(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path)
    missing = _build_args(bundle)
    missing.pop("expected_config_sha256")
    with pytest.raises(TypeError, match="expected_config_sha256"):
        construct_later_sector_inventory(**missing)

    wrong = _build_args(bundle)
    wrong["expected_config_sha256"] = "0" * 64
    with pytest.raises(FM0ContractError, match="config.*hash|hash.*config"):
        construct_later_sector_inventory(**wrong)

    pinned = _build_args(bundle)
    bundle.config_path.write_text(
        bundle.config_path.read_text(encoding="utf-8") + "# byte drift\n",
        encoding="utf-8",
    )
    with pytest.raises(FM0ContractError, match="config.*hash|hash.*config"):
        construct_later_sector_inventory(**pinned)


@pytest.mark.parametrize(
    ("field", "value", "match"),
    (
        ("receipt_schema_version", "wrong", "schema mismatch"),
        ("sector", 66, "does not match|requested accepted sector"),
        ("product_state", "ORCD_COMPLETE_DEFERRED", "product_state|accepted sector"),
        ("passed", False, "not passed"),
    ),
)
def test_receipt_schema_sector_state_and_passed_are_exact(
    tmp_path: Path, field: str, value: Any, match: str
) -> None:
    bundle = _make_bundle(tmp_path)
    _mutate_receipt(bundle, 65, lambda payload: payload.__setitem__(field, value))

    with pytest.raises(FM0ContractError, match=match):
        construct_later_sector_inventory(**_build_args(bundle))


def test_receipt_and_evidence_paths_and_hashes_are_exact(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    wrong_receipts = list(bundle.receipt_tuples)
    sector, path, _ = wrong_receipts[0]
    wrong_receipts[0] = (sector, path, "0" * 64)
    args = _build_args(bundle)
    args["ordered_sector_receipts"] = wrong_receipts
    with pytest.raises(FM0ContractError, match="receipt hash mismatch"):
        construct_later_sector_inventory(**args)

    bundle = _make_bundle(tmp_path / "evidence-hash")
    gate = next(iter(bundle.config["admission"]["required_evidence_gates"]))
    _mutate_receipt(
        bundle,
        65,
        lambda payload: payload["evidence_gates"][gate].__setitem__("sha256", "0" * 64),
    )
    with pytest.raises(FM0ContractError, match="hash mismatch"):
        construct_later_sector_inventory(**_build_args(bundle))

    bundle = _make_bundle(tmp_path / "evidence-path")
    gate = next(iter(bundle.config["admission"]["required_evidence_gates"]))
    _mutate_receipt(
        bundle,
        65,
        lambda payload: payload["evidence_gates"][gate].__setitem__(
            "path", "missing.json"
        ),
    )
    with pytest.raises(FM0ContractError, match="missing .* evidence gate"):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize("mode", ("missing", "extra"))
def test_evidence_gate_set_must_be_exact(tmp_path: Path, mode: str) -> None:
    bundle = _make_bundle(tmp_path)
    gate = next(iter(bundle.config["admission"]["required_evidence_gates"]))

    def mutate(payload: dict[str, Any]) -> None:
        if mode == "missing":
            payload["evidence_gates"].pop(gate)
        else:
            payload["evidence_gates"]["unexpected_gate"] = dict(
                payload["evidence_gates"][gate]
            )

    _mutate_receipt(bundle, 65, mutate)
    with pytest.raises(FM0ContractError, match="evidence gates differ"):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    "mode",
    (
        "missing_schema",
        "wrong_schema",
        "missing_gate_name",
        "wrong_gate_name",
        "wrong_sector",
        "wrong_product_state",
        "not_passed",
        "extra_field",
    ),
)
def test_each_evidence_artifact_has_an_exact_gate_specific_json_contract(
    tmp_path: Path, mode: str
) -> None:
    bundle = _make_bundle(tmp_path)
    gates = list(bundle.config["admission"]["required_evidence_gates"])
    gate = gates[0]

    def mutate(payload: dict[str, Any]) -> None:
        if mode == "missing_schema":
            payload.pop("schema_version")
        elif mode == "wrong_schema":
            payload["schema_version"] = "wrong"
        elif mode == "missing_gate_name":
            payload.pop("gate_name")
        elif mode == "wrong_gate_name":
            payload["gate_name"] = gates[1]
        elif mode == "wrong_sector":
            payload["sector"] = 66
        elif mode == "wrong_product_state":
            payload["product_state"] = "ORCD_COMPLETE_DEFERRED"
        elif mode == "not_passed":
            payload["passed"] = False
        elif mode == "extra_field":
            payload["generic_claim"] = True
        else:  # pragma: no cover
            raise AssertionError(mode)

    _mutate_evidence(bundle, 65, gate, mutate)
    with pytest.raises(
        FM0ContractError,
        match="evidence.*(schema|gate|sector|state|passed|fields)|JSON evidence",
    ):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    "mode",
    (
        "missing_upstream_field",
        "empty_upstream",
        "missing_binding_field",
        "missing_artifact",
        "hash_mismatch",
        "wrong_content_class",
        "duplicate_role",
    ),
)
def test_evidence_upstream_control_artifacts_are_nonempty_and_hash_bound(
    tmp_path: Path, mode: str
) -> None:
    bundle = _make_bundle(tmp_path)
    gate = bundle.config["admission"]["required_evidence_gates"][0]

    def mutate(payload: dict[str, Any]) -> None:
        if mode == "missing_upstream_field":
            payload.pop("upstream_artifacts")
        elif mode == "empty_upstream":
            payload["upstream_artifacts"] = []
        elif mode == "missing_binding_field":
            payload["upstream_artifacts"][0].pop("sha256")
        elif mode == "missing_artifact":
            payload["upstream_artifacts"][0]["path"] = "missing_control.json"
        elif mode == "hash_mismatch":
            payload["upstream_artifacts"][0]["sha256"] = "0" * 64
        elif mode == "wrong_content_class":
            payload["upstream_artifacts"][0]["content_class"] = "science_data"
        elif mode == "duplicate_role":
            payload["upstream_artifacts"].append(dict(payload["upstream_artifacts"][0]))
        else:  # pragma: no cover
            raise AssertionError(mode)

    _mutate_evidence(bundle, 65, gate, mutate)
    with pytest.raises(
        FM0ContractError,
        match="upstream|artifact|binding|hash|content|role|schema",
    ):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    "mode",
    (
        "empty_json",
        "arbitrary_json",
        "extra_field",
        "wrong_gate",
        "wrong_sector",
        "wrong_state",
        "not_passed",
        "wrong_schema_version",
        "wrong_metric_contract_version",
        "wrapper_metric_mismatch",
    ),
)
def test_upstream_control_json_exactly_certifies_its_gate_and_metrics(
    tmp_path: Path, mode: str
) -> None:
    bundle = _make_bundle(tmp_path)
    gates = bundle.config["admission"]["required_evidence_gates"]
    gate = gates[0]

    def mutate(payload: dict[str, Any]) -> None:
        if mode == "empty_json":
            payload.clear()
        elif mode == "arbitrary_json":
            payload.clear()
            payload["arbitrary"] = True
        elif mode == "extra_field":
            payload["generic_claim"] = True
        elif mode == "wrong_gate":
            payload["gate_name"] = gates[1]
        elif mode == "wrong_sector":
            payload["sector"] = 66
        elif mode == "wrong_state":
            payload["product_state"] = "ORCD_COMPLETE_DEFERRED"
        elif mode == "not_passed":
            payload["passed"] = False
        elif mode == "wrong_schema_version":
            payload["schema_version"] = "twirl_fm0_later_sector_upstream_control_v0"
        elif mode == "wrong_metric_contract_version":
            payload["metric_contract_version"] = "generic_metrics_v1"
        elif mode == "wrapper_metric_mismatch":
            payload["metrics"]["n_products"] = 99
            payload["metrics"]["n_checksum_bound_products"] = 99
        else:  # pragma: no cover
            raise AssertionError(mode)

    _mutate_upstream_control(bundle, 65, gate, mutate)
    with pytest.raises(
        FM0ContractError,
        match="upstream|control|schema|gate|sector|state|passed|metric",
    ):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    "mode",
    (
        "missing_producer",
        "producer_extra_field",
        "bad_git_commit",
        "missing_code_binding",
        "missing_config_binding",
        "code_binding_extra_field",
        "bad_code_hash",
        "bad_config_hash",
        "missing_code_file",
        "missing_config_file",
    ),
)
def test_upstream_control_producer_is_exact_and_hash_bound(
    tmp_path: Path, mode: str
) -> None:
    bundle = _make_bundle(tmp_path)
    gate = bundle.config["admission"]["required_evidence_gates"][0]

    def mutate(payload: dict[str, Any]) -> None:
        if mode == "missing_producer":
            payload.pop("producer")
        elif mode == "producer_extra_field":
            payload["producer"]["environment"] = "unknown"
        elif mode == "bad_git_commit":
            payload["producer"]["git_commit"] = "A" * 40
        elif mode == "missing_code_binding":
            payload["producer"].pop("code")
        elif mode == "missing_config_binding":
            payload["producer"].pop("config")
        elif mode == "code_binding_extra_field":
            payload["producer"]["code"]["kind"] = "script"
        elif mode == "bad_code_hash":
            payload["producer"]["code"]["sha256"] = "0" * 64
        elif mode == "bad_config_hash":
            payload["producer"]["config"]["sha256"] = "0" * 64
        elif mode == "missing_code_file":
            payload["producer"]["code"]["path"] = "missing_producer.py"
        elif mode == "missing_config_file":
            payload["producer"]["config"]["path"] = "missing_producer.yaml"
        else:  # pragma: no cover
            raise AssertionError(mode)

    _mutate_upstream_control(bundle, 65, gate, mutate)
    with pytest.raises(
        FM0ContractError,
        match="producer|code|config|commit|hash|missing|upstream|schema",
    ):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize("artifact", ("code", "config"))
def test_upstream_producer_artifacts_must_be_nonempty(
    tmp_path: Path, artifact: str
) -> None:
    bundle = _make_bundle(tmp_path)
    gate = bundle.config["admission"]["required_evidence_gates"][0]
    _rewrite_producer_artifact(bundle, 65, gate, artifact, "")

    with pytest.raises(FM0ContractError, match="producer|code|config|empty"):
        construct_later_sector_inventory(**_build_args(bundle))


def test_cross_gate_upstream_control_reuse_fails(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    gates = bundle.config["admission"]["required_evidence_gates"]
    first, second = gates[:2]
    receipt_path = bundle.receipt_paths[65]
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    first_evidence_path = receipt_path.parent / receipt["evidence_gates"][first]["path"]
    second_evidence_path = (
        receipt_path.parent / receipt["evidence_gates"][second]["path"]
    )
    first_evidence = json.loads(first_evidence_path.read_text(encoding="utf-8"))
    second_evidence = json.loads(second_evidence_path.read_text(encoding="utf-8"))
    second_evidence["upstream_artifacts"] = [
        dict(first_evidence["upstream_artifacts"][0])
    ]
    _write_json(second_evidence_path, second_evidence)

    def update(payload: dict[str, Any]) -> None:
        payload["evidence_gates"][second]["sha256"] = _sha256(second_evidence_path)

    _mutate_receipt(bundle, 65, update)
    with pytest.raises(FM0ContractError, match="upstream|reuse|gate|unique"):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    ("gate", "metric", "value"),
    (
        ("checksum_bound_a2v1_hdf5_provenance", "n_checksum_bound_products", 99),
        ("checksum_bound_a2v1_fits_provenance", "n_checksum_mismatches", 1),
        ("edge_aware_coverage", "n_non_edge_omissions", 1),
        ("hdf5_openability", "n_hdf5_opened", 99),
        ("authoritative_internal_cadence_quality", "n_observations_failed", 1),
        ("authoritative_external_cadence_quality", "n_cadences_checked", 0),
        ("explicit_omissions", "n_unexplained_omissions", 1),
        (
            "fm_channel_mask_finite_numerical_envelope",
            "n_nonfinite_active_values",
            1,
        ),
        ("stable_physical_source_registry_join", "n_unmatched_identity_rows", 1),
    ),
)
def test_gate_specific_evidence_metrics_must_reconcile_and_pass(
    tmp_path: Path, gate: str, metric: str, value: int
) -> None:
    bundle = _make_bundle(tmp_path)
    _mutate_evidence(
        bundle,
        65,
        gate,
        lambda payload: payload["metrics"].__setitem__(metric, value),
    )

    with pytest.raises(FM0ContractError, match="metric|evidence gate|reconcil"):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize("mode", ("missing_metric", "extra_metric"))
def test_gate_specific_metric_fields_are_exact(tmp_path: Path, mode: str) -> None:
    bundle = _make_bundle(tmp_path)
    gate = "hdf5_openability"

    def mutate(payload: dict[str, Any]) -> None:
        if mode == "missing_metric":
            payload["metrics"].pop("n_hdf5_opened")
        else:
            payload["metrics"]["generic_pass_count"] = 100

    _mutate_evidence(bundle, 65, gate, mutate)
    with pytest.raises(FM0ContractError, match="metric|schema|evidence gate"):
        construct_later_sector_inventory(**_build_args(bundle))


def test_evidence_artifacts_must_be_json_and_distinct_per_gate(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path)
    gates = list(bundle.config["admission"]["required_evidence_gates"])
    first, second = gates[:2]

    def reuse(payload: dict[str, Any]) -> None:
        payload["evidence_gates"][second] = dict(payload["evidence_gates"][first])

    _mutate_receipt(bundle, 65, reuse)
    with pytest.raises(FM0ContractError, match="distinct|unique|reus"):
        construct_later_sector_inventory(**_build_args(bundle))

    bundle = _make_bundle(tmp_path / "not-json")
    gate = gates[0]
    receipt_path = bundle.receipt_paths[65]
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    source = receipt_path.parent / receipt["evidence_gates"][gate]["path"]
    non_json = source.with_suffix(".txt")
    non_json.write_bytes(source.read_bytes())

    def replace(payload: dict[str, Any]) -> None:
        payload["evidence_gates"][gate]["path"] = str(
            non_json.relative_to(receipt_path.parent)
        )
        payload["evidence_gates"][gate]["sha256"] = _sha256(non_json)

    _mutate_receipt(bundle, 65, replace)
    with pytest.raises(FM0ContractError, match=r"JSON|\.json"):
        construct_later_sector_inventory(**_build_args(bundle))


def test_stage1_acceptance_receipt_binding_is_required_and_hash_bound(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path)
    _mutate_receipt(
        bundle,
        65,
        lambda payload: payload.pop("stage1_acceptance_receipt"),
    )
    with pytest.raises(FM0ContractError, match="Stage.?1|stage1|acceptance"):
        construct_later_sector_inventory(**_build_args(bundle))

    bundle = _make_bundle(tmp_path / "hash")
    _mutate_receipt(
        bundle,
        65,
        lambda payload: payload["stage1_acceptance_receipt"].__setitem__(
            "sha256", "0" * 64
        ),
    )
    with pytest.raises(FM0ContractError, match="hash mismatch"):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("sector", 66),
        ("orbits", [999, 1000]),
        ("ok", False),
        ("ok_h5", False),
        ("ok_fits", False),
        ("expected_contract", {"ok": False}),
    ),
)
def test_stage1_acceptance_receipt_must_certify_the_exact_accepted_sector(
    tmp_path: Path, field: str, value: Any
) -> None:
    bundle = _make_bundle(tmp_path)
    _mutate_stage1_acceptance(
        bundle,
        65,
        lambda payload: payload.__setitem__(field, value),
    )

    with pytest.raises(
        FM0ContractError,
        match="Stage.?1|stage1|acceptance|orbits|sector|ok",
    ):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize("mode", ("h5_not_opened", "fits_partially_checked"))
def test_stage1_acceptance_requires_open_h5_and_every_fits_checked(
    tmp_path: Path, mode: str
) -> None:
    bundle = _make_bundle(tmp_path)

    def mutate(payload: dict[str, Any]) -> None:
        if mode == "h5_not_opened":
            payload["schema"]["check_h5_open"] = False
        elif mode == "fits_partially_checked":
            payload["fits"]["n_checked_fits"] = payload["fits"]["n_fits"] - 1
        else:  # pragma: no cover
            raise AssertionError(mode)

    _mutate_stage1_acceptance(bundle, 65, mutate)
    with pytest.raises(FM0ContractError, match="HDF5|FITS|Stage.?1|full-schema"):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    ("sectors", "selected", "match"),
    (
        ((65, 66, 68, 69, 70), (65, 66, 68, 69, 70), "contiguous"),
        ((65, 66, 67, 68, 69), (65, 67, 66, 68, 69), "prefix|chronological"),
    ),
)
def test_noncontiguous_and_nonchronological_selections_fail(
    tmp_path: Path,
    sectors: Sequence[int],
    selected: Sequence[int],
    match: str,
) -> None:
    bundle = _make_bundle(tmp_path, sectors=sectors, selected_sectors=selected)
    with pytest.raises(FM0ContractError, match=match):
        construct_later_sector_inventory(**_build_args(bundle))


def test_duplicate_or_implicit_sector_authority_fails(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    duplicate = list(bundle.receipt_tuples)
    duplicate[-1] = duplicate[-2]
    args = _build_args(bundle)
    args["ordered_sector_receipts"] = duplicate
    with pytest.raises(FM0ContractError, match="duplicate|prefix"):
        construct_later_sector_inventory(**args)

    bundle = _make_bundle(tmp_path / "implicit")
    bundle.config["selection"]["expected_orbits_by_sector"].pop(69)
    bundle.config_path.write_text(
        yaml.safe_dump(bundle.config, sort_keys=False), encoding="utf-8"
    )
    with pytest.raises(FM0ContractError, match="explicit .*orbit|expected-orbit"):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    ("case", "match"),
    (
        ("drop_gate", "evidence.*gate|required evidence"),
        ("rename_gate", "evidence.*gate|required evidence"),
        ("disable_fail_closed", "fail.closed"),
        ("candidate_drift", "candidate.*chronological|candidate-sector prefix"),
        ("orbit_map_drift", "explicit orbit mapping|Stage-1 acceptance identity"),
    ),
)
def test_fail_closed_policy_and_sector_authorities_cannot_be_tampered(
    tmp_path: Path, case: str, match: str
) -> None:
    bundle = _make_bundle(tmp_path)
    if case == "drop_gate":
        bundle.config["admission"]["required_evidence_gates"].pop()
    elif case == "rename_gate":
        bundle.config["admission"]["required_evidence_gates"][0] = "renamed_gate"
    elif case == "disable_fail_closed":
        bundle.config["admission"]["fail_closed"] = False
    elif case == "candidate_drift":
        bundle.config["selection"]["ordered_candidate_sectors"][-1] = 70
    elif case == "orbit_map_drift":
        bundle.config["selection"]["expected_orbits_by_sector"][65] = [999, 1000]
    else:  # pragma: no cover
        raise AssertionError(case)
    bundle.config_path.write_text(
        yaml.safe_dump(bundle.config, sort_keys=False), encoding="utf-8"
    )

    with pytest.raises(FM0ContractError, match=match):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize("mode", ("component", "partition"))
def test_sector_alias_authority_rejects_arbitrary_component_or_partition(
    tmp_path: Path, mode: str
) -> None:
    bundle = _make_bundle(tmp_path)
    rows = _read_csv(bundle.manifest_paths[65])
    row = next(item for item in rows if "_new_" in item["observation_key"])
    if mode == "component":
        row["leakage_component_id"] = "leakage_" + "f" * 64
        row["source_partition"] = deterministic_source_partition(
            row["leakage_component_id"]
        )[0]
    else:
        row["source_partition"] = (
            "poc_development" if row["source_partition"] == "poc_train" else "poc_train"
        )
    _write_csv(bundle.manifest_paths[65], rows)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    with pytest.raises(FM0ContractError, match="component|partition|alias authority"):
        construct_later_sector_inventory(**_build_args(bundle))


def test_new_alias_edge_touching_or_merging_baseline_graph_fails(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path)
    rows = _read_csv(bundle.alias_paths[65])
    rows.append(
        {
            "gaia_dr3_source_id": bundle.repeated_hosts[0].gaia,
            "tic_id": bundle.new_hosts[0].tics[0],
        }
    )
    _write_csv(bundle.alias_paths[65], rows, ("gaia_dr3_source_id", "tic_id"))
    _refresh_receipt_binding(bundle, 65, "alias_authority")

    with pytest.raises(FM0ContractError, match="baseline|frozen|component|closure"):
        construct_later_sector_inventory(**_build_args(bundle))


def test_baseline_component_and_visit_count_drift_fail_semantically(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path)
    selection = _read_csv(bundle.baseline_selection_path)
    selection[0]["leakage_component_id"] = "leakage_" + "e" * 64
    _write_csv(bundle.baseline_selection_path, selection, CORPUS_SELECTION_FIELDS)
    manifest = _read_csv(bundle.baseline_manifest_path)
    manifest[0]["leakage_component_id"] = "leakage_" + "e" * 64
    _write_csv(bundle.baseline_manifest_path, manifest, MANIFEST_COLUMNS)
    bundle.baseline_selection_sha256 = _sha256(bundle.baseline_selection_path)
    bundle.baseline_manifest_sha256 = _sha256(bundle.baseline_manifest_path)
    _refresh_config_hashes(bundle)
    with pytest.raises(FM0ContractError, match="baseline.*component|alias.*drift"):
        construct_later_sector_inventory(**_build_args(bundle))

    bundle = _make_bundle(tmp_path / "visit-count")
    manifest = _read_csv(bundle.baseline_manifest_path)
    _write_csv(bundle.baseline_manifest_path, manifest[:-1], MANIFEST_COLUMNS)
    bundle.baseline_manifest_sha256 = _sha256(bundle.baseline_manifest_path)
    _refresh_config_hashes(bundle)
    with pytest.raises(FM0ContractError, match="visit count|per-component|selection"):
        construct_later_sector_inventory(**_build_args(bundle))


def test_forbidden_label_or_search_column_fails_before_inventory(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path)
    rows = _read_csv(bundle.manifest_paths[65])
    for row in rows:
        row["human_label"] = "unknown"
    fields = (*tuple(rows[0].keys())[:-1], "human_label")
    _write_csv(bundle.manifest_paths[65], rows, fields)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    with pytest.raises(FM0ContractError, match="forbidden|policy-forbidden"):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize("value", (False, ""))
def test_later_rows_must_be_scientifically_eligible(tmp_path: Path, value: Any) -> None:
    bundle = _make_bundle(tmp_path)
    rows = _read_csv(bundle.manifest_paths[65])
    rows[0]["scientific_training_eligible"] = value
    _write_csv(bundle.manifest_paths[65], rows)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    inventory = construct_later_sector_inventory(**_build_args(bundle))
    assert not any(
        row["observation_key"] == rows[0]["observation_key"]
        for row in inventory.eligible_observations
    )
    assert any(
        row["observation_key"] == rows[0]["observation_key"]
        and row["reason"] == "scientific_training_ineligible"
        for row in inventory.quarantine
    )


@pytest.mark.parametrize(
    "field", ("product_instance_id", "source_sha256", "product_state")
)
def test_product_identity_and_accepted_state_fields_are_required(
    tmp_path: Path, field: str
) -> None:
    bundle = _make_bundle(tmp_path)
    rows = _read_csv(bundle.manifest_paths[65])
    rows[0][field] = ""
    _write_csv(bundle.manifest_paths[65], rows)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    match = (
        "observation/product identity"
        if field == "product_instance_id"
        else "source_sha256"
        if field == "source_sha256"
        else "A2V1_ACCEPTED"
    )
    with pytest.raises(FM0ContractError, match=match):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    "mode",
    ("missing_schema", "wrong_schema", "legacy_extra_field"),
)
def test_later_observation_manifest_schema_is_exact(tmp_path: Path, mode: str) -> None:
    bundle = _make_bundle(tmp_path)
    rows = _read_csv(bundle.manifest_paths[65])
    if mode == "missing_schema":
        for row in rows:
            row.pop("schema_version")
    elif mode == "wrong_schema":
        rows[0]["schema_version"] = "twirl_fm0_later_sector_observation_v0"
    elif mode == "legacy_extra_field":
        for row in rows:
            row["registry_schema_version"] = "legacy"
    else:  # pragma: no cover
        raise AssertionError(mode)
    _write_csv(bundle.manifest_paths[65], rows)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    with pytest.raises(
        FM0ContractError, match="observation.*schema|schema.*observation"
    ):
        construct_later_sector_inventory(**_build_args(bundle))


@pytest.mark.parametrize(
    "layout",
    (
        list(reversed(FLUX_VIEW_NAMES)),
        list(FLUX_VIEW_NAMES[:-1]),
        [*FLUX_VIEW_NAMES[:-1], "adp015_5x5"],
        "not-json",
    ),
)
def test_later_manifest_declares_the_exact_canonical_six_view_layout(
    tmp_path: Path, layout: Any
) -> None:
    bundle = _make_bundle(tmp_path)
    rows = _read_csv(bundle.manifest_paths[65])
    rows[0]["flux_view_names_json"] = (
        layout if isinstance(layout, str) else json.dumps(layout, separators=(",", ":"))
    )
    _write_csv(bundle.manifest_paths[65], rows)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    with pytest.raises(FM0ContractError, match="flux.*view|view.*layout|canonical"):
        construct_later_sector_inventory(**_build_args(bundle))


def test_one_detector_inventory_exposes_coverage_but_cannot_be_panel_ready(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path, n_repeated=64, n_new=256)
    rows = _read_csv(bundle.manifest_paths[65])
    for row in rows:
        row["camera"] = "1"
        row["ccd"] = "1"
    _write_csv(bundle.manifest_paths[65], rows)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    inventory = construct_later_sector_inventory(**_build_args(bundle))
    sector = next(row for row in inventory.sector_inventory if row["sector"] == 65)
    assert json.loads(sector["camera_ccd_coverage_json"]) == ["cam1/ccd1"]
    assert inventory.summary["count_floor_ready"] is True
    assert inventory.summary["adequacy_thresholds_frozen"] is False
    assert inventory.summary["panel_freeze_ready"] is False


def test_low_retention_inventory_exposes_floor_but_cannot_be_panel_ready(
    tmp_path: Path,
) -> None:
    bundle = _make_bundle(tmp_path, n_repeated=64, n_new=256)
    rows = _read_csv(bundle.manifest_paths[65])
    for row in rows:
        row["cadence_retention_fraction"] = "0.0001"
    _write_csv(bundle.manifest_paths[65], rows)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    inventory = construct_later_sector_inventory(**_build_args(bundle))
    assert inventory.summary[
        "eligible_cadence_retention_fraction_min"
    ] == pytest.approx(0.0001)
    assert inventory.summary["count_floor_ready"] is True
    assert inventory.summary["adequacy_thresholds_frozen"] is False
    assert inventory.summary["panel_freeze_ready"] is False


@pytest.mark.parametrize(
    ("cohort_token", "sector_count_field"),
    (
        ("_repeated_", "n_repeated_host_components"),
        ("_new_", "n_new_host_components"),
    ),
)
def test_weak_per_sector_cohort_coverage_is_visible_but_not_panel_ready(
    tmp_path: Path, cohort_token: str, sector_count_field: str
) -> None:
    bundle = _make_bundle(tmp_path, n_repeated=64, n_new=256)
    rows = _read_csv(bundle.manifest_paths[65])
    retained_one = False
    kept: list[dict[str, str]] = []
    for row in rows:
        if cohort_token not in row["observation_key"]:
            kept.append(row)
        elif not retained_one:
            kept.append(row)
            retained_one = True
    _write_csv(bundle.manifest_paths[65], kept)
    _refresh_receipt_binding(bundle, 65, "observation_manifest")

    inventory = construct_later_sector_inventory(**_build_args(bundle))
    sector = next(row for row in inventory.sector_inventory if row["sector"] == 65)
    assert sector[sector_count_field] == 1
    assert inventory.summary["count_floor_ready"] is True
    assert inventory.summary["adequacy_thresholds_frozen"] is False
    assert inventory.summary["panel_freeze_ready"] is False


def test_immutable_inventory_refuses_nonidentical_overwrite(tmp_path: Path) -> None:
    bundle = _make_bundle(tmp_path)
    output = tmp_path / "inventory"
    first = build_later_sector_inventory(**_build_args(bundle), output_dir=output)
    assert (
        build_later_sector_inventory(**_build_args(bundle), output_dir=output) == first
    )
    (output / "summary.json").write_text("{}\n", encoding="utf-8")

    with pytest.raises(FM0ContractError, match="refusing to replace"):
        build_later_sector_inventory(**_build_args(bundle), output_dir=output)
