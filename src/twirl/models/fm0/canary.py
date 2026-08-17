"""Deterministic, label-blind S56 canary selection for TWIRL-FM0.1.

The canary validates the real A2v1 input path before a larger release is
materialized.  Selection uses only Gaia/TIC identity, detector geometry, orbit
coverage, and the frozen hash salt.  Search products, labels, and light-curve
values are deliberately absent.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .registry import AliasRegistry, FM0ContractError, publish_immutable


CANARY_SCHEMA_VERSION = "twirl_fm0_1_s56_canary_selection_v1"
CANARY_SALT = "twirl_fm0_1_s56_canary_v1"
S56_ORBITS = frozenset((119, 120))
WD1856_GAIA_DR3_SOURCE_ID = "2146576589564898688"
WD1856_TIC_ID = "267574918"

CANARY_FIELDS = (
    "schema_version",
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "camera",
    "ccd",
    "leakage_component_id",
    "source_partition",
    "rank_sha256",
    "is_benchmark",
    "orbits_json",
    "hdf5_paths_json",
)


def _text_identifier(value: Any, *, name: str) -> str:
    text = str(value).strip()
    if not text or not text.isdigit() or int(text) <= 0:
        raise FM0ContractError(f"{name} must be a positive integer identifier")
    return str(int(text))


def _optional_identifier(value: Any) -> str:
    text = str(value).strip()
    if not text or not text.isdigit() or int(text) <= 0:
        return ""
    return str(int(text))


def _boolean(value: Any) -> bool:
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no", ""}:
            return False
        raise FM0ContractError(f"unparseable boolean value: {value!r}")
    return bool(value)


def _rank(component_id: str, gaia: str, tic: str, camera: int, ccd: int) -> str:
    payload = f"{CANARY_SALT}:{component_id}:{gaia}:{tic}:{camera}:{ccd}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def select_s56_canary(
    observation_rows: Iterable[Mapping[str, Any]],
    alias_registry: AliasRegistry,
    *,
    pdo_root: str | Path,
    per_detector: int = 4,
    benchmark_gaia: str = WD1856_GAIA_DR3_SOURCE_ID,
    benchmark_tic: str = WD1856_TIC_ID,
) -> tuple[tuple[dict[str, Any], ...], tuple[dict[str, str], ...]]:
    """Select balanced two-orbit sources and their complete alias closure."""

    if per_detector <= 0:
        raise FM0ContractError("per_detector must be positive")
    root = Path(pdo_root)
    alias_index = alias_registry.alias_index()
    grouped: dict[tuple[str, str, int, int], set[int]] = {}
    for row in observation_rows:
        try:
            sector = int(row["sector"])
        except (KeyError, TypeError, ValueError) as exc:
            raise FM0ContractError("observation rows require an integer sector") from exc
        if sector != 56 or _boolean(row.get("edge_warn", False)):
            continue
        gaia = _text_identifier(row.get("source_id"), name="source_id")
        tic = _optional_identifier(row.get("tic_id"))
        if not tic:
            continue
        try:
            orbit = int(row["orbit"])
            camera = int(row["camera"])
            ccd = int(row["ccd"])
        except (KeyError, TypeError, ValueError) as exc:
            raise FM0ContractError("observation rows require orbit/camera/ccd") from exc
        if orbit not in S56_ORBITS or camera not in range(1, 5) or ccd not in range(1, 5):
            continue
        alias = alias_index.get((gaia, tic))
        if alias is None or bool(alias["quarantined"]):
            continue
        grouped.setdefault((gaia, tic, camera, ccd), set()).add(orbit)

    candidates: dict[tuple[int, int], list[dict[str, Any]]] = {
        (camera, ccd): [] for camera in range(1, 5) for ccd in range(1, 5)
    }
    for (gaia, tic, camera, ccd), orbits in grouped.items():
        if orbits != S56_ORBITS:
            continue
        alias = alias_index[(gaia, tic)]
        rank = _rank(alias["leakage_component_id"], gaia, tic, camera, ccd)
        paths = [
            str(root / f"orbit-{orbit}" / "ffi" / f"cam{camera}" / f"ccd{ccd}" / "LC" / f"{tic}.h5")
            for orbit in sorted(orbits)
        ]
        candidates[(camera, ccd)].append(
            {
                "schema_version": CANARY_SCHEMA_VERSION,
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "sector": 56,
                "camera": camera,
                "ccd": ccd,
                "leakage_component_id": alias["leakage_component_id"],
                "source_partition": alias["source_partition"],
                "rank_sha256": rank,
                "is_benchmark": gaia == benchmark_gaia and tic == benchmark_tic,
                "orbits_json": json.dumps(sorted(orbits), separators=(",", ":")),
                "hdf5_paths_json": json.dumps(paths, separators=(",", ":")),
            }
        )

    selected: list[dict[str, Any]] = []
    for detector, rows in sorted(candidates.items()):
        rows.sort(key=lambda row: (row["rank_sha256"], int(row["gaia_dr3_source_id"]), int(row["tic_id"])))
        if len(rows) < per_detector:
            raise FM0ContractError(
                f"S56 cam{detector[0]}/ccd{detector[1]} has only {len(rows)} eligible two-orbit sources"
            )
        selected.extend(rows[:per_detector])

    benchmark_key = (str(int(benchmark_gaia)), str(int(benchmark_tic)))
    if not any(
        (row["gaia_dr3_source_id"], row["tic_id"]) == benchmark_key
        for row in selected
    ):
        benchmark_rows = [
            row
            for rows in candidates.values()
            for row in rows
            if (row["gaia_dr3_source_id"], row["tic_id"]) == benchmark_key
        ]
        if len(benchmark_rows) != 1:
            raise FM0ContractError("WD 1856 is missing or ambiguous in the eligible S56 observations")
        selected.append(benchmark_rows[0])

    selected.sort(key=lambda row: (row["camera"], row["ccd"], row["rank_sha256"]))
    selected_components = {row["leakage_component_id"] for row in selected}
    aliases = tuple(
        {
            "gaia_dr3_source_id": str(row["gaia_dr3_source_id"]),
            "tic_id": str(row["tic_id"]),
        }
        for row in alias_registry.aliases
        if row["leakage_component_id"] in selected_components
    )
    return tuple(selected), aliases


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row[field] for field in fields})
    return stream.getvalue().encode("utf-8")


def write_s56_canary_selection(
    out_dir: str | Path,
    selection: Sequence[Mapping[str, Any]],
    aliases: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Publish immutable selection and closure files plus a compact summary."""

    if not selection or not aliases:
        raise FM0ContractError("cannot publish an empty canary selection")
    output = Path(out_dir)
    selection_bytes = _csv_bytes(selection, CANARY_FIELDS)
    alias_fields = ("gaia_dr3_source_id", "tic_id")
    alias_bytes = _csv_bytes(aliases, alias_fields)
    selection_sha = hashlib.sha256(selection_bytes).hexdigest()
    alias_sha = hashlib.sha256(alias_bytes).hexdigest()
    summary = {
        "schema_version": CANARY_SCHEMA_VERSION,
        "sector": 56,
        "selection_salt": CANARY_SALT,
        "n_sources": len(selection),
        "n_benchmark_sources": sum(_boolean(row["is_benchmark"]) for row in selection),
        "n_alias_edges": len(aliases),
        "detector_counts": {
            f"cam{camera}_ccd{ccd}": sum(
                int(row["camera"]) == camera and int(row["ccd"]) == ccd
                for row in selection
            )
            for camera in range(1, 5)
            for ccd in range(1, 5)
        },
        "outputs_sha256": {
            "selection.csv": selection_sha,
            "aliases.csv": alias_sha,
        },
        "claim_limit": "input-path canary only; not a model-performance sample",
    }
    publish_immutable(output / "selection.csv", selection_bytes)
    publish_immutable(output / "aliases.csv", alias_bytes)
    publish_immutable(
        output / "summary.json",
        (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    return summary
