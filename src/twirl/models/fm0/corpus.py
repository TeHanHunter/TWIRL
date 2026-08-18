"""Deterministic, label-blind multi-sector source selection for FM0.1.

The selection is an audit/control artifact.  It uses only Gaia/TIC identity,
sector/orbit/detector placement, edge status, and the frozen alias-component
registry.  It never reads light-curve values, search products, labels, or
candidate information.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

from .registry import AliasRegistry, FM0ContractError, publish_immutable


CORPUS_SELECTION_SCHEMA_VERSION = "twirl_fm0_1_corpus_selection_v1"
CORPUS_SELECTION_FIELDS = (
    "schema_version",
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "camera",
    "ccd",
    "leakage_component_id",
    "source_partition",
    "orbits_json",
    "hdf5_paths_json",
)
CORPUS_AUDIT_FIELDS = (
    "gaia_dr3_source_id",
    "sector",
    "reason",
    "n_complete_products",
    "complete_products_json",
)


def expected_orbits(sector: int) -> tuple[int, int]:
    """Return the two 200-s orbit identifiers for one TESS sector."""

    sector = int(sector)
    if sector < 56:
        raise FM0ContractError("FM0.1 sectors must be at least 56")
    first = 2 * sector + 7
    return first, first + 1


def _identifier(value: Any, *, name: str) -> str:
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


def select_corpus_observations(
    observation_rows: Iterable[Mapping[str, Any]],
    alias_registry: AliasRegistry,
    *,
    sectors: Sequence[int],
) -> tuple[
    tuple[dict[str, Any], ...],
    tuple[dict[str, str], ...],
    tuple[dict[str, Any], ...],
]:
    """Select one unambiguous, complete product for each Gaia source visit.

    A candidate product is a fixed ``(TIC, camera, CCD)`` tuple represented in
    both expected sector orbits.  If a Gaia source/sector has more than one
    complete tuple, the visit is quarantined rather than resolving the
    duplicate through an undocumented preference.  Components containing
    multiple Gaia sources are already quarantined by ``AliasRegistry``.
    """

    sector_set = {int(value) for value in sectors}
    if not sector_set or min(sector_set) < 56:
        raise FM0ContractError("corpus sectors must be a non-empty subset of S56+")
    alias_index = alias_registry.alias_index()
    products: dict[
        tuple[str, int], dict[tuple[str, int, int], set[int]]
    ] = {}
    for row in observation_rows:
        try:
            sector = int(row["sector"])
        except (KeyError, TypeError, ValueError) as exc:
            raise FM0ContractError("observation rows require an integer sector") from exc
        if sector not in sector_set or _boolean(row.get("edge_warn", False)):
            continue
        gaia = _identifier(row.get("source_id"), name="source_id")
        tic = _optional_identifier(row.get("tic_id"))
        if not tic:
            continue
        alias = alias_index.get((gaia, tic))
        if alias is None or bool(alias["quarantined"]):
            continue
        try:
            orbit = int(row["orbit"])
            camera = int(row["camera"])
            ccd = int(row["ccd"])
        except (KeyError, TypeError, ValueError) as exc:
            raise FM0ContractError(
                "observation rows require integer orbit/camera/ccd"
            ) from exc
        if orbit not in expected_orbits(sector):
            continue
        if camera not in range(1, 5) or ccd not in range(1, 5):
            continue
        products.setdefault((gaia, sector), {}).setdefault(
            (tic, camera, ccd), set()
        ).add(orbit)

    selected: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    selected_components: set[str] = set()
    for (gaia, sector), candidates in sorted(
        products.items(), key=lambda item: (item[0][1], int(item[0][0]))
    ):
        required = set(expected_orbits(sector))
        complete = sorted(
            product for product, orbits in candidates.items() if orbits == required
        )
        if len(complete) != 1:
            audit.append(
                {
                    "gaia_dr3_source_id": gaia,
                    "sector": sector,
                    "reason": (
                        "no_complete_two_orbit_product"
                        if not complete
                        else "multiple_complete_two_orbit_products"
                    ),
                    "n_complete_products": len(complete),
                    "complete_products_json": json.dumps(
                        [
                            {"tic_id": item[0], "camera": item[1], "ccd": item[2]}
                            for item in complete
                        ],
                        separators=(",", ":"),
                    ),
                }
            )
            continue
        tic, camera, ccd = complete[0]
        alias = alias_index[(gaia, tic)]
        component = str(alias["leakage_component_id"])
        selected_components.add(component)
        orbits = sorted(required)
        relative_paths = [
            (
                f"orbit-{orbit}/ffi/cam{camera}/ccd{ccd}/LC/{tic}.h5"
            )
            for orbit in orbits
        ]
        selected.append(
            {
                "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "sector": sector,
                "camera": camera,
                "ccd": ccd,
                "leakage_component_id": component,
                "source_partition": alias["source_partition"],
                "orbits_json": json.dumps(orbits, separators=(",", ":")),
                "hdf5_paths_json": json.dumps(relative_paths, separators=(",", ":")),
            }
        )

    if not selected:
        raise FM0ContractError("corpus selection contains no complete observations")
    selected.sort(key=lambda row: (int(row["sector"]), int(row["gaia_dr3_source_id"])))
    aliases = tuple(
        {
            "gaia_dr3_source_id": str(row["gaia_dr3_source_id"]),
            "tic_id": str(row["tic_id"]),
        }
        for row in alias_registry.aliases
        if str(row["leakage_component_id"]) in selected_components
    )
    return tuple(selected), aliases, tuple(audit)


def iter_observation_fits(
    path: str | Path,
    *,
    sectors: Sequence[int],
    chunk_size: int = 250_000,
) -> Iterator[dict[str, Any]]:
    """Yield the small allowlisted geometry view from a large observation FITS."""

    try:
        from astropy.io import fits
        import numpy as np
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError("astropy and numpy are required to read observation FITS") from exc
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    sector_values = np.asarray(sorted({int(value) for value in sectors}), dtype=np.int64)
    with fits.open(Path(path), memmap=True) as hdul:
        data = hdul[1].data
        required = {"source_id", "tic_id", "sector", "orbit", "camera", "ccd", "edge_warn"}
        missing = sorted(required - set(data.names or ()))
        if missing:
            raise FM0ContractError(f"observation FITS lacks required columns: {missing}")
        for start in range(0, len(data), chunk_size):
            chunk = data[start : start + chunk_size]
            mask = np.isin(np.asarray(chunk["sector"], dtype=np.int64), sector_values)
            indices = np.flatnonzero(mask)
            for index in indices:
                yield {
                    "source_id": str(int(chunk["source_id"][index])),
                    "tic_id": str(int(chunk["tic_id"][index])),
                    "sector": int(chunk["sector"][index]),
                    "orbit": int(chunk["orbit"][index]),
                    "camera": int(chunk["camera"][index]),
                    "ccd": int(chunk["ccd"][index]),
                    "edge_warn": bool(chunk["edge_warn"][index]),
                }


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row[field] for field in fields})
    return stream.getvalue().encode("utf-8")


def write_corpus_selection(
    out_dir: str | Path,
    selection: Sequence[Mapping[str, Any]],
    aliases: Sequence[Mapping[str, Any]],
    audit: Sequence[Mapping[str, Any]],
    *,
    sectors: Sequence[int],
    input_authorities: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Publish the immutable selection, alias closure, and ambiguity ledger."""

    if not selection or not aliases:
        raise FM0ContractError("cannot publish an empty corpus selection")
    output = Path(out_dir)
    selection_payload = _csv_bytes(selection, CORPUS_SELECTION_FIELDS)
    alias_payload = _csv_bytes(aliases, ("gaia_dr3_source_id", "tic_id"))
    audit_payload = _csv_bytes(audit, CORPUS_AUDIT_FIELDS)
    outputs_sha256 = {
        "selection.csv": hashlib.sha256(selection_payload).hexdigest(),
        "aliases.csv": hashlib.sha256(alias_payload).hexdigest(),
        "quarantine.csv": hashlib.sha256(audit_payload).hexdigest(),
    }
    counts_by_sector = {
        str(sector): sum(int(row["sector"]) == sector for row in selection)
        for sector in sorted({int(value) for value in sectors})
    }
    visits_by_partition: dict[str, int] = {}
    sources_by_partition: dict[str, set[str]] = {}
    visits_by_source: dict[str, set[int]] = {}
    for row in selection:
        partition = str(row["source_partition"])
        source = str(row["gaia_dr3_source_id"])
        visits_by_partition[partition] = visits_by_partition.get(partition, 0) + 1
        sources_by_partition.setdefault(partition, set()).add(source)
        visits_by_source.setdefault(source, set()).add(int(row["sector"]))
    repeated = [len(values) for values in visits_by_source.values() if len(values) > 1]
    repeated_histogram = {
        str(count): repeated.count(count) for count in sorted(set(repeated))
    }
    summary = {
        "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
        "sectors": sorted({int(value) for value in sectors}),
        "n_observations": len(selection),
        "n_physical_sources": len({str(row["gaia_dr3_source_id"]) for row in selection}),
        "n_leakage_components": len({str(row["leakage_component_id"]) for row in selection}),
        "n_alias_edges": len(aliases),
        "n_quarantined_visits": len(audit),
        "observations_by_sector": counts_by_sector,
        "observations_by_partition": dict(sorted(visits_by_partition.items())),
        "sources_by_partition": {
            key: len(values) for key, values in sorted(sources_by_partition.items())
        },
        "n_repeated_sources": len(repeated),
        "repeated_source_sector_count_histogram": repeated_histogram,
        "selection_policy": "one_complete_two_orbit_tic_detector_product_per_gaia_sector",
        "model_fields_exposed": [],
        "input_authorities": dict(input_authorities or {}),
        "outputs_sha256": outputs_sha256,
        "claim_limit": "BLS-free input selection; not a model result",
    }
    publish_immutable(output / "selection.csv", selection_payload)
    publish_immutable(output / "aliases.csv", alias_payload)
    publish_immutable(output / "quarantine.csv", audit_payload)
    publish_immutable(
        output / "summary.json",
        (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    return summary


__all__ = [
    "CORPUS_AUDIT_FIELDS",
    "CORPUS_SELECTION_FIELDS",
    "CORPUS_SELECTION_SCHEMA_VERSION",
    "expected_orbits",
    "iter_observation_fits",
    "select_corpus_observations",
    "write_corpus_selection",
]
