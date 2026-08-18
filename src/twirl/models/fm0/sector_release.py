"""Sector-wise construction and merge of the FM0.1 six-view input release.

The source HDF5 archives are deliberately expanded one sector at a time.  This
module applies the same discipline to the derived six-view release: each
sector gets a checksum-bound, immutable partial release, and a later merge
creates one global manifest without rereading or concatenating raw light
curves.  The final merge recomputes host-visit offsets from recorded absolute
visit bounds, so a star observed in several sectors retains its real temporal
spacing without becoming one giant model input.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from .a2v1_adapter import A2V1_HDF5_MANIFEST_FIELDS
from .input_release import (
    BUILD_SUMMARY_SCHEMA_VERSION,
    CADENCE_SECONDS,
    ERROR_VIEW_NAMES,
    FLUX_VIEW_NAMES,
    INPUT_RELEASE_SCHEMA_VERSION,
    MANIFEST_COLUMNS,
    VISIT_TIMING_COLUMNS,
    write_a2v1_hdf5_input_release,
)
from .registry import (
    FM0ContractError,
    FrozenContract,
    load_frozen_contract,
    publish_immutable,
    read_rows,
    sha256_file,
)


SECTOR_INPUT_SCHEMA_VERSION = "twirl_fm0_1_sector_input_release_v1"
SECTOR_INPUT_MERGE_SCHEMA_VERSION = "twirl_fm0_1_sector_input_merge_v1"
_SHA40 = re.compile(r"^[0-9a-f]{40}$")
_REGISTRY_OBSERVATION_FIELDS = (
    "registry_schema_version",
    "observation_key",
    "product_instance_id",
    "physical_source_id",
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "a2v1_product_version",
    "source_sha256",
    "product_state",
    "diagnostic_admission_receipt_path",
    "diagnostic_admission_receipt_sha256",
    "leakage_component_id",
    "source_partition",
    "quarantined",
)


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row[field] for field in fields})
    return stream.getvalue().encode("utf-8")


def _read_exact_csv(path: Path, fields: Sequence[str], *, label: str) -> list[dict[str, str]]:
    rows = read_rows(path)
    if not rows:
        raise FM0ContractError(f"{label} is empty: {path}")
    if tuple(rows[0]) != tuple(fields):
        raise FM0ContractError(f"{label} columns drifted: {path}")
    if any(tuple(row) != tuple(fields) for row in rows):
        raise FM0ContractError(f"{label} row columns drifted: {path}")
    return rows


def _is_quarantined(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1"}


def _sector_name(sector: int) -> str:
    return f"s{int(sector):04d}"


def _relative_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise FM0ContractError(f"unsafe derived shard path: {value!r}")
    return path


def build_sector_input_release(
    *,
    sector: int,
    registry_dir: str | Path,
    hdf5_manifest_path: str | Path,
    output_root: str | Path,
    producer_git_sha: str,
    contract: FrozenContract | None = None,
) -> dict[str, Any]:
    """Publish a BLS-free scientific input sub-release for one sector.

    ``registry_dir`` is the full fixed split.  Keeping it whole makes the
    selected source partition identical in every sector, while the temporary
    per-sector source manifest bounds memory and makes failed sectors
    independently restartable.
    """

    sector = int(sector)
    if sector < 56:
        raise FM0ContractError("FM0.1 sector input releases require S56+")
    if not _SHA40.fullmatch(producer_git_sha):
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    contract = contract or load_frozen_contract()
    registry = Path(registry_dir).resolve()
    observations_path = registry / "observations.csv"
    observations = _read_exact_csv(
        observations_path,
        _REGISTRY_OBSERVATION_FIELDS,
        label="registry observations",
    )
    selected_observations = {
        row["observation_key"]: row
        for row in observations
        if int(row["sector"]) == sector and not _is_quarantined(row["quarantined"])
    }
    if not selected_observations:
        raise FM0ContractError(f"registry has no non-quarantined S{sector} observations")

    hdf5_manifest = Path(hdf5_manifest_path).resolve()
    source_rows = _read_exact_csv(
        hdf5_manifest,
        A2V1_HDF5_MANIFEST_FIELDS,
        label="bound A2v1 HDF5 manifest",
    )
    selected_sources = [
        row for row in source_rows if row["observation_key"] in selected_observations
    ]
    source_keys = {row["observation_key"] for row in selected_sources}
    if source_keys != set(selected_observations):
        missing = sorted(set(selected_observations) - source_keys)
        unexpected = sorted(source_keys - set(selected_observations))
        raise FM0ContractError(
            f"S{sector} source manifest differs from registry; "
            f"missing={len(missing)}, unexpected={len(unexpected)}"
        )
    if len(source_keys) != len(selected_sources):
        raise FM0ContractError(f"S{sector} source manifest has duplicate observation keys")

    root = Path(output_root).resolve()
    final = root / _sector_name(sector)
    partial = final.with_name(final.name + ".partial")
    if final.exists() or partial.exists():
        raise FM0ContractError(f"refusing to overwrite sector input release: {final}")
    partial.mkdir(parents=True)
    selected_sources.sort(key=lambda row: row["observation_key"])
    source_payload = _csv_bytes(selected_sources, A2V1_HDF5_MANIFEST_FIELDS)
    selected_manifest = partial / "a2v1_hdf5_manifest.csv"
    publish_immutable(selected_manifest, source_payload)
    summary = write_a2v1_hdf5_input_release(
        registry_dir=registry,
        hdf5_manifest_path=selected_manifest,
        out_dir=partial,
        contract=contract,
        visit_timing_path=partial / "visit_timing.csv",
    )
    for path in (
        partial / "manifest.csv",
        partial / "summary.json",
        partial / "visit_timing.csv",
    ):
        if not path.is_file() or not path.stat().st_size:
            raise FM0ContractError(f"S{sector} writer omitted {path.name}")
    sector_summary = {
        "schema_version": SECTOR_INPUT_SCHEMA_VERSION,
        "producer_git_sha": producer_git_sha,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": sector,
        "registry_observations_sha256": sha256_file(observations_path),
        "a2v1_hdf5_manifest_sha256": hashlib.sha256(source_payload).hexdigest(),
        "manifest_sha256": sha256_file(partial / "manifest.csv"),
        "summary_sha256": sha256_file(partial / "summary.json"),
        "visit_timing_sha256": sha256_file(partial / "visit_timing.csv"),
        "n_observations": int(summary["n_observations"]),
        "n_cadences": int(summary["n_cadences"]),
        "input_adapter": summary["input_adapter"],
        "scientific_training_eligible": bool(summary["scientific_training_eligible"]),
        "claim_limit": "one sector of a derived input release; not a model result",
    }
    if sector_summary["n_observations"] != len(selected_sources):
        raise FM0ContractError(f"S{sector} writer observation count changed")
    if (
        sector_summary["input_adapter"] != "a2v1_hdf5_quality_aware_v1"
        or not sector_summary["scientific_training_eligible"]
    ):
        raise FM0ContractError(f"S{sector} did not create a scientific A2v1 release")
    publish_immutable(
        partial / "sector_summary.json",
        (json.dumps(sector_summary, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    publish_immutable(partial / "READY", (producer_git_sha + "\n").encode("utf-8"))
    os.replace(partial, final)
    return sector_summary


def _recompute_host_visit_timing(
    manifest_rows: list[dict[str, Any]], timing_rows: Sequence[Mapping[str, str]]
) -> None:
    manifest_by_key = {row["observation_key"]: row for row in manifest_rows}
    if len(manifest_by_key) != len(manifest_rows):
        raise FM0ContractError("duplicate observation keys in input-release merge")
    visits_by_source: dict[str, list[tuple[dict[str, Any], float, float]]] = {}
    seen: set[str] = set()
    for timing in timing_rows:
        key = str(timing["observation_key"])
        if key in seen or key not in manifest_by_key:
            raise FM0ContractError(f"invalid timing observation key: {key!r}")
        seen.add(key)
        source = str(timing["physical_source_id"]).strip()
        try:
            start = float(timing["absolute_visit_start"])
            end = float(timing["absolute_visit_end"])
        except (TypeError, ValueError) as exc:
            raise FM0ContractError(f"invalid timing row for {key}") from exc
        if not source or not math.isfinite(start) or not math.isfinite(end) or end < start:
            raise FM0ContractError(f"invalid visit bounds for {key}")
        visits_by_source.setdefault(source, []).append((manifest_by_key[key], start, end))
    if seen != set(manifest_by_key):
        raise FM0ContractError("input-release timing rows do not cover every manifest row")
    cadence_units_per_day = 86400.0 / CADENCE_SECONDS
    for source, visits in visits_by_source.items():
        visits.sort(key=lambda item: (item[1], str(item[0]["observation_key"])))
        first_start = visits[0][1]
        prior_start: float | None = None
        prior_end: float | None = None
        for row, start, end in visits:
            if prior_start is not None and start <= prior_start:
                raise FM0ContractError(
                    f"non-monotonic sector visit starts for physical source {source}"
                )
            offset = (start - first_start) * cadence_units_per_day
            gap = 0.0 if prior_end is None else (start - prior_end) * cadence_units_per_day
            if not math.isfinite(offset) or not math.isfinite(gap) or offset < 0:
                raise FM0ContractError(f"invalid derived host timing for {row['observation_key']}")
            row["host_visit_offset_cadences"] = repr(offset)
            row["host_visit_gap_cadences"] = repr(gap)
            row["host_visit_overlaps_previous"] = bool(prior_end is not None and gap < 0)
            prior_start = start
            prior_end = end


def merge_sector_input_releases(
    *,
    sector_root: str | Path,
    sectors: Sequence[int],
    registry_dir: str | Path,
    out_dir: str | Path,
    producer_git_sha: str,
    allowed_sector_producer_git_shas: Sequence[str] | None = None,
    contract: FrozenContract | None = None,
) -> dict[str, Any]:
    """Hard-link ready sector releases into one full, integrity-bound release."""

    if not _SHA40.fullmatch(producer_git_sha):
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    allowed = set(allowed_sector_producer_git_shas or (producer_git_sha,))
    if not allowed or any(not _SHA40.fullmatch(value) for value in allowed):
        raise FM0ContractError("approved sector input revisions are invalid")
    contract = contract or load_frozen_contract()
    registry = Path(registry_dir).resolve()
    observations_path = registry / "observations.csv"
    observations = _read_exact_csv(
        observations_path,
        _REGISTRY_OBSERVATION_FIELDS,
        label="registry observations",
    )
    expected_keys = {
        row["observation_key"] for row in observations if not _is_quarantined(row["quarantined"])
    }
    requested_sectors = tuple(sorted({int(value) for value in sectors}))
    if not requested_sectors or any(sector < 56 for sector in requested_sectors):
        raise FM0ContractError("input-release merge sectors must be a non-empty S56+ set")
    expected_sectors = {int(row["sector"]) for row in observations if not _is_quarantined(row["quarantined"])}
    if set(requested_sectors) != expected_sectors:
        raise FM0ContractError(
            f"merge sectors differ from the fixed registry: requested={requested_sectors}, "
            f"registry={tuple(sorted(expected_sectors))}"
        )

    root = Path(sector_root).resolve()
    output = Path(out_dir).resolve()
    partial = output.with_name(output.name + ".partial")
    if output.exists() or partial.exists():
        raise FM0ContractError(f"refusing to overwrite merged input release: {output}")
    partial.mkdir(parents=True)
    (partial / "shards").mkdir()

    manifest_rows: list[dict[str, Any]] = []
    timing_rows: list[dict[str, str]] = []
    source_manifest_rows: list[dict[str, str]] = []
    seen_keys: set[str] = set()
    sector_bindings: dict[str, Any] = {}
    try:
        for sector in requested_sectors:
            sector_dir = root / _sector_name(sector)
            ready = sector_dir / "READY"
            sector_summary_path = sector_dir / "sector_summary.json"
            manifest_path = sector_dir / "manifest.csv"
            timing_path = sector_dir / "visit_timing.csv"
            source_path = sector_dir / "a2v1_hdf5_manifest.csv"
            summary_path = sector_dir / "summary.json"
            required = (ready, sector_summary_path, manifest_path, timing_path, source_path, summary_path)
            if any(not path.is_file() for path in required):
                raise FM0ContractError(f"S{sector} sector input release is incomplete")
            stage_sha = ready.read_text(encoding="utf-8").strip()
            if stage_sha not in allowed:
                raise FM0ContractError(f"S{sector} sector input uses an unapproved revision")
            stage_summary = json.loads(sector_summary_path.read_text(encoding="utf-8"))
            if (
                stage_summary.get("schema_version") != SECTOR_INPUT_SCHEMA_VERSION
                or stage_summary.get("producer_git_sha") != stage_sha
                or int(stage_summary.get("sector", -1)) != sector
                or stage_summary.get("manifest_sha256") != sha256_file(manifest_path)
                or stage_summary.get("summary_sha256") != sha256_file(summary_path)
                or stage_summary.get("visit_timing_sha256") != sha256_file(timing_path)
                or stage_summary.get("a2v1_hdf5_manifest_sha256") != sha256_file(source_path)
                or stage_summary.get("registry_observations_sha256") != sha256_file(observations_path)
                or stage_summary.get("input_adapter") != "a2v1_hdf5_quality_aware_v1"
                or stage_summary.get("scientific_training_eligible") is not True
            ):
                raise FM0ContractError(f"S{sector} sector input binding failed")
            stage_manifest = _read_exact_csv(manifest_path, MANIFEST_COLUMNS, label=f"S{sector} input manifest")
            stage_timing = _read_exact_csv(timing_path, VISIT_TIMING_COLUMNS, label=f"S{sector} visit timing")
            stage_source = _read_exact_csv(source_path, A2V1_HDF5_MANIFEST_FIELDS, label=f"S{sector} source manifest")
            stage_keys = {row["observation_key"] for row in stage_manifest}
            if len(stage_keys) != len(stage_manifest) or stage_keys & seen_keys:
                raise FM0ContractError(f"S{sector} duplicate observation keys in merge")
            if any(int(row["n_cadences"]) <= 0 for row in stage_manifest):
                raise FM0ContractError(f"S{sector} has an empty derived shard")
            for row in stage_manifest:
                if row["input_release_schema_version"] != INPUT_RELEASE_SCHEMA_VERSION:
                    raise FM0ContractError(f"S{sector} input schema drift")
                if row["input_adapter"] != "a2v1_hdf5_quality_aware_v1" or row["scientific_training_eligible"] != "True":
                    raise FM0ContractError(f"S{sector} input eligibility drift")
                relative = _relative_path(row["relative_path"])
                source = sector_dir / relative
                destination = partial / relative
                if not source.is_file() or destination.exists():
                    raise FM0ContractError(f"S{sector} missing or duplicate derived shard: {relative}")
                try:
                    os.link(source, destination)
                except OSError as exc:
                    raise FM0ContractError(
                        f"S{sector} cannot hard-link immutable derived shard: {relative}"
                    ) from exc
            seen_keys.update(stage_keys)
            manifest_rows.extend(stage_manifest)
            timing_rows.extend(stage_timing)
            source_manifest_rows.extend(stage_source)
            sector_bindings[str(sector)] = {
                "producer_git_sha": stage_sha,
                "sector_summary_sha256": sha256_file(sector_summary_path),
                "manifest_sha256": sha256_file(manifest_path),
                "visit_timing_sha256": sha256_file(timing_path),
                "n_observations": len(stage_manifest),
            }

        if seen_keys != expected_keys:
            raise FM0ContractError(
                f"merged release does not match fixed registry; "
                f"missing={len(expected_keys - seen_keys)}, extra={len(seen_keys - expected_keys)}"
            )
        _recompute_host_visit_timing(manifest_rows, timing_rows)
        manifest_rows.sort(key=lambda row: row["observation_key"])
        timing_rows.sort(key=lambda row: row["observation_key"])
        source_manifest_rows.sort(key=lambda row: row["observation_key"])
        if len({row["observation_key"] for row in source_manifest_rows}) != len(source_manifest_rows):
            raise FM0ContractError("duplicate source-manifest key in input-release merge")
        manifest_payload = _csv_bytes(manifest_rows, MANIFEST_COLUMNS)
        timing_payload = _csv_bytes(timing_rows, VISIT_TIMING_COLUMNS)
        source_payload = _csv_bytes(source_manifest_rows, A2V1_HDF5_MANIFEST_FIELDS)
        summary = {
            "summary_schema_version": BUILD_SUMMARY_SCHEMA_VERSION,
            "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
            "campaign_id": contract.config["campaign_id"],
            "design_sha256": contract.design_sha256,
            "config_sha256": contract.config_sha256,
            "freeze_receipt_sha256": contract.freeze_receipt_sha256,
            "registry_observations_sha256": sha256_file(observations_path),
            "hdf5_manifest_sha256": hashlib.sha256(source_payload).hexdigest(),
            "manifest_sha256": hashlib.sha256(manifest_payload).hexdigest(),
            "n_observations": len(manifest_rows),
            "n_cadences": sum(int(row["n_cadences"]) for row in manifest_rows),
            "flux_view_names": list(FLUX_VIEW_NAMES),
            "error_view_names": list(ERROR_VIEW_NAMES),
            "input_adapter": "a2v1_hdf5_quality_aware_v1",
            "scientific_training_eligible": True,
            "host_visit_timing_derivation": "absolute_raw_time_grouped_by_physical_source_id",
            "host_visit_gap_definition": "signed_current_start_minus_previous_end",
            "partial_release": False,
            "certifies_full_campaign": False,
            "sectors": list(requested_sectors),
            "sector_bindings": sector_bindings,
            "claim_limit": "S56-S64 FM0.1 input release only; not a model result",
        }
        publish_immutable(partial / "manifest.csv", manifest_payload)
        publish_immutable(partial / "visit_timing.csv", timing_payload)
        publish_immutable(partial / "a2v1_hdf5_manifest.csv", source_payload)
        publish_immutable(
            partial / "summary.json",
            (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        )
        merge_summary = {
            "schema_version": SECTOR_INPUT_MERGE_SCHEMA_VERSION,
            "producer_git_sha": producer_git_sha,
            "allowed_sector_producer_git_shas": sorted(allowed),
            "sectors": list(requested_sectors),
            "summary_sha256": sha256_file(partial / "summary.json"),
            "manifest_sha256": sha256_file(partial / "manifest.csv"),
            "visit_timing_sha256": sha256_file(partial / "visit_timing.csv"),
            "n_observations": len(manifest_rows),
            "n_cadences": summary["n_cadences"],
            "sector_bindings": sector_bindings,
            "claim_limit": "merge only; an independent full release validation is still required",
        }
        publish_immutable(
            partial / "merge_summary.json",
            (json.dumps(merge_summary, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        )
        publish_immutable(partial / "READY", (producer_git_sha + "\n").encode("utf-8"))
        os.replace(partial, output)
        return merge_summary
    except Exception:
        # The partial directory deliberately remains for forensic inspection;
        # without READY it cannot be mistaken for a consumable release.
        raise


__all__ = [
    "SECTOR_INPUT_MERGE_SCHEMA_VERSION",
    "SECTOR_INPUT_SCHEMA_VERSION",
    "build_sector_input_release",
    "merge_sector_input_releases",
]
