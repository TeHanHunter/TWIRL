"""Label-blind retained-source inventory for later TWIRL-FM sectors.

The ORCD A2v1 lane retains one immutable output bundle per
``(sector, orbit, camera, CCD)`` cell.  This module verifies that receipt chain,
recovers Gaia--TIC identities without reading search products or labels, and
publishes a checksum-bound inventory grouped across a sector's two orbits.

``source_tic`` files are not, by themselves, a corpus-selection authority:
the same ECSV schema can be produced from an all-TIC catalog, and even the
observation-derived overlays include ``edge_warn`` targets.  A separately
checksum-bound corpus selection is therefore mandatory for source-row
publication.  Observation inventories are geometry metadata only because
they do not bind leakage components or the frozen source partition.  The
overlays are used only to cross-check Gaia--TIC identity.  Unselected HDF5
products are counted and excluded rather than silently entering the model
inventory.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path, PurePosixPath
from typing import Any

from .corpus import (
    CORPUS_SELECTION_FIELDS,
    CORPUS_SELECTION_SCHEMA_VERSION,
    iter_observation_fits,
)
from .orcd_source_admission import (
    SOURCE_READY_STATE,
    SOURCE_RECEIPT_SCHEMA_VERSION,
    verify_retained_sector_source,
)
from .registry import (
    FM0ContractError,
    assert_no_search_columns,
    build_alias_registry,
    deterministic_source_partition,
    publish_immutable,
    sha256_file,
)

INVENTORY_SCHEMA_VERSION = "twirl_fm0_later_retained_source_inventory_v1"
ALIAS_SCHEMA_VERSION = "twirl_fm0_later_retained_alias_v1"
SOURCE_ROW_SCHEMA_VERSION = "twirl_fm0_later_retained_source_row_v1"
SUMMARY_SCHEMA_VERSION = "twirl_fm0_later_retained_source_summary_v1"
IDENTITY_SOURCE_CORPUS_SELECTION = "checksum_bound_corpus_selection_v1"
ALIAS_IDENTITY_SOURCE = "manifest_bound_source_tic_crosscheck_v1"
EXPECTED_SOURCE_GRID_FILES = 196

ALIAS_FIELDS = (
    "alias_schema_version",
    "sector",
    "gaia_dr3_source_id",
    "tic_id",
    "identity_source",
    "target_authority_sha256",
    "identity_ambiguous",
    "has_retained_hdf5",
)
SOURCE_FIELDS = (
    "source_row_schema_version",
    "sector",
    "gaia_dr3_source_id",
    "tic_id",
    "leakage_component_id",
    "source_partition",
    "identity_source",
    "target_authority_sha256",
    "corpus_selection_sha256",
    "corpus_summary_sha256",
    "gaia_tic_alias_authority_sha256",
    "identity_ambiguous",
    "product_state",
    "n_hdf5_products",
    "orbits_json",
    "hdf5_paths_json",
    "hdf5_sha256_json",
    "hdf5_evidence_json",
    "retained_source_sha256",
    "source_row_sha256",
    "source_receipt_path",
    "source_receipt_sha256",
    "model_outcome_blind",
    "panel_admission_authorized",
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_CELL = re.compile(
    r"^s(?P<sector>[0-9]{4})_o(?P<orbit>[0-9]+)_cam(?P<camera>[1-4])_ccd(?P<ccd>[1-4])$"
)
_OVERLAY = re.compile(r"^source_[0-9]+_[0-9]+\.ecsv$")
_LEAKAGE_COMPONENT = re.compile(r"^leakage_[0-9a-f]{64}$")
_SOURCE_PARTITIONS = frozenset({"poc_train", "poc_development", "poc_sealed_test"})


def _make_output_tree_read_only(root: Path) -> None:
    """Remove every write bit before publishing an inventory tree."""

    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(
            f"later-source output is not a materialized directory: {root}"
        )
    for path in sorted(root.rglob("*"), key=lambda value: len(value.parts), reverse=True):
        if path.is_symlink():
            raise FM0ContractError(f"later-source output contains a symlink: {path}")
        path.chmod(0o550 if path.is_dir() else 0o440)
    root.chmod(0o550)


def _assert_output_tree_read_only(root: Path) -> None:
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(
            f"later-source output is not a materialized directory: {root}"
        )
    for path in (root, *root.rglob("*")):
        if path.is_symlink() or path.stat().st_mode & 0o222:
            raise FM0ContractError(f"later-source output is not read-only: {path}")


def _remove_owned_output_tree(root: Path) -> None:
    """Remove exactly one builder-owned tree without following symlinks."""

    if not root.exists() and not root.is_symlink():
        return
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(f"refusing to clean unsafe output path: {root}")
    for directory, names, _files in os.walk(root, topdown=True, followlinks=False):
        directory_path = Path(directory)
        directory_path.chmod(0o700)
        for name in names:
            child = directory_path / name
            if not child.is_symlink():
                child.chmod(0o700)
    shutil.rmtree(root)


def _output_leaf_path(value: str | Path, *, label: str) -> Path:
    """Resolve the parent while preserving the output leaf for symlink checks."""

    requested = Path(value).expanduser()
    if requested.name in {"", ".", ".."}:
        raise FM0ContractError(f"{label} must name a concrete output directory")
    return requested.parent.resolve() / requested.name


def _read_json(path: Path, *, label: str) -> Mapping[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise FM0ContractError(f"{label} must be a materialized file: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not isinstance(value, Mapping):
        raise FM0ContractError(f"{label} must be a JSON object: {path}")
    return value


def _sha256_text(value: Any, *, label: str) -> str:
    text = str(value).strip().lower()
    if _SHA256.fullmatch(text) is None:
        raise FM0ContractError(f"{label} must be a lowercase SHA-256 digest")
    return text


def _positive_identifier(value: Any, *, label: str) -> str:
    if value is None:
        raise FM0ContractError(f"{label} is missing")
    text = str(value).strip()
    if not text or not text.isdigit() or int(text) <= 0:
        raise FM0ContractError(f"{label} must be a positive decimal integer")
    return str(int(text))


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _csv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        if tuple(row) != tuple(fields):
            raise FM0ContractError("later-source inventory row columns drifted")
        writer.writerow({field: row[field] for field in fields})
    return stream.getvalue().encode("utf-8")


def _load_tics(path: Path, *, expected_sha256: str) -> tuple[int, ...]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError(f"missing materialized TIC authority: {path}")
    if sha256_file(path) != expected_sha256:
        raise FM0ContractError(f"TIC authority hash drifted: {path}")
    values: list[int] = []
    try:
        for line_number, raw in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            value = raw.strip()
            if not value or value.startswith("#"):
                continue
            if value != raw or not value.isdigit() or int(value) <= 0:
                raise FM0ContractError(
                    f"invalid TIC authority row at {path}:{line_number}"
                )
            values.append(int(value))
    except UnicodeError as exc:
        raise FM0ContractError(f"invalid TIC authority encoding: {path}") from exc
    if not values or values != sorted(set(values)):
        raise FM0ContractError(
            f"TIC authority must be nonempty, sorted, and unique: {path}"
        )
    return tuple(values)


def _manifest_records(manifest: Mapping[str, Any]) -> dict[str, tuple[int, str]]:
    records = manifest.get("files")
    if not isinstance(records, list) or not records:
        raise FM0ContractError("cell manifest has no file records")
    result: dict[str, tuple[int, str]] = {}
    for raw in records:
        if not isinstance(raw, Mapping):
            raise FM0ContractError("cell manifest contains a non-object file record")
        relative_text = str(raw.get("path", ""))
        relative = PurePosixPath(relative_text)
        if (
            not relative_text
            or relative.is_absolute()
            or ".." in relative.parts
            or relative.parts[0] not in {"source", "source_tic", "metadata"}
        ):
            raise FM0ContractError(f"unsafe cell-manifest path: {relative_text!r}")
        try:
            size = int(raw.get("bytes", -1))
        except (TypeError, ValueError) as exc:
            raise FM0ContractError(
                f"invalid byte count in cell manifest: {relative_text!r}"
            ) from exc
        digest = _sha256_text(
            raw.get("sha256"), label=f"cell-manifest digest for {relative_text}"
        )
        if size < 0 or relative_text in result:
            raise FM0ContractError(
                f"invalid/duplicate cell-manifest path: {relative_text}"
            )
        result[relative_text] = (size, digest)
    return result


def _verify_manifest_file(
    *, root: Path, relative: str, expected_size: int, expected_sha256: str
) -> Path:
    path = root / relative
    if path.is_symlink() or not path.is_file():
        raise FM0ContractError(f"manifest-bound file is not materialized: {path}")
    if path.stat().st_size != expected_size or sha256_file(path) != expected_sha256:
        raise FM0ContractError(f"manifest-bound file drifted: {path}")
    return path


def _output_manifest(path: Path) -> dict[str, str]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError(f"missing output manifest: {path}")
    result: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except UnicodeError as exc:
        raise FM0ContractError(f"invalid output manifest encoding: {path}") from exc
    for line_number, raw in enumerate(lines, start=1):
        fields = raw.split(None, 1)
        if len(fields) != 2:
            raise FM0ContractError(
                f"malformed output manifest row {path}:{line_number}"
            )
        digest = _sha256_text(fields[0], label="output-manifest file digest")
        relative = PurePosixPath(fields[1].strip())
        relative_text = relative.as_posix()
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or len(relative.parts) != 2
            or relative.parts[0] not in {"epsf", "LC"}
            or relative_text in result
        ):
            raise FM0ContractError(
                f"unsafe/duplicate output-manifest path: {relative_text!r}"
            )
        result[relative_text] = digest
    if not result:
        raise FM0ContractError(f"empty output manifest: {path}")
    return result


def _actual_output_paths(retained_root: Path) -> set[str]:
    actual: set[str] = set()
    for directory, pattern in (("epsf", "epsf_*.npy"), ("LC", "*.h5")):
        root = retained_root / directory
        if root.is_symlink() or not root.is_dir():
            raise FM0ContractError(f"retained output directory is missing: {root}")
        for path in root.glob(pattern):
            if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
                raise FM0ContractError(f"invalid retained output: {path}")
            actual.add(f"{directory}/{path.name}")
    return actual


def _overlay_pairs(path: Path) -> set[tuple[str, str]]:
    # Keep Astropy local so importing the finished inventory contract remains
    # lightweight on consumers that only inspect CSV/JSON artifacts.
    from astropy.table import Table

    try:
        table = Table.read(path, format="ascii.ecsv")
    except Exception as exc:  # Astropy exposes several parser exception types.
        raise FM0ContractError(f"cannot read source_tic overlay: {path}") from exc
    if table.colnames != ["TIC", "gaia3"]:
        raise FM0ContractError(
            f"source_tic overlay columns drifted in {path}: {table.colnames!r}"
        )
    pairs: set[tuple[str, str]] = set()
    for row in table:
        tic = _positive_identifier(row["TIC"], label=f"TIC in {path}")
        gaia = _positive_identifier(row["gaia3"], label=f"gaia3 in {path}")
        pairs.add((gaia, tic))
    return pairs


def _inventory_rows(path: Path) -> tuple[list[dict[str, Any]], tuple[str, ...]]:
    suffix = path.suffix.lower()
    if suffix == ".json":
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise FM0ContractError(f"invalid target inventory: {path}") from exc
        if not isinstance(value, list) or not all(
            isinstance(row, Mapping) for row in value
        ):
            raise FM0ContractError("JSON target inventory must be a list of objects")
        rows = [dict(row) for row in value]
        columns = tuple(rows[0]) if rows else ()
    elif suffix == ".csv":
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                columns = tuple(reader.fieldnames or ())
                rows = [dict(row) for row in reader]
        except (OSError, UnicodeError, csv.Error) as exc:
            raise FM0ContractError(f"invalid target inventory: {path}") from exc
    else:
        from astropy.table import Table

        try:
            table = Table.read(path)
        except Exception as exc:
            raise FM0ContractError(
                f"cannot read authoritative target inventory: {path}"
            ) from exc
        columns = tuple(str(name) for name in table.colnames)
        rows = [{name: row[name] for name in table.colnames} for row in table]
    if not rows or not columns:
        raise FM0ContractError(f"authoritative target inventory is empty: {path}")
    if any(tuple(row) != columns for row in rows):
        raise FM0ContractError("target-inventory row columns drifted")
    assert_no_search_columns(columns, context="authoritative target inventory")
    return rows, columns


def _boolean(value: Any, *, label: str) -> bool:
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no", ""}:
            return False
        raise FM0ContractError(f"{label} is not a valid boolean")
    return bool(value)


def _selection_from_observations(
    rows: Sequence[Mapping[str, Any]], *, sector: int, expected_orbits: tuple[int, int]
) -> tuple[dict[tuple[str, str], tuple[int, int]], dict[str, int]]:
    """Reproduce the edge-safe, complete-two-orbit geometry selection."""

    products: dict[str, dict[tuple[str, int, int], set[int]]] = {}
    n_sector_rows = 0
    n_edge_rows = 0
    for row in rows:
        try:
            row_sector = int(row["sector"])
        except (KeyError, TypeError, ValueError) as exc:
            raise FM0ContractError(
                "observation inventory requires integer sector"
            ) from exc
        if row_sector != sector:
            continue
        n_sector_rows += 1
        if _boolean(row.get("edge_warn", False), label="edge_warn"):
            n_edge_rows += 1
            continue
        gaia = _positive_identifier(row.get("source_id"), label="source_id")
        try:
            tic = _positive_identifier(row.get("tic_id"), label="tic_id")
        except FM0ContractError:
            # Missing TIC aliases are visible in the parent catalog but cannot
            # identify an HDF5 product.
            continue
        try:
            orbit = int(row["orbit"])
            camera = int(row["camera"])
            ccd = int(row["ccd"])
        except (KeyError, TypeError, ValueError) as exc:
            raise FM0ContractError(
                "observation inventory requires integer orbit/camera/ccd"
            ) from exc
        if (
            orbit not in expected_orbits
            or camera not in range(1, 5)
            or ccd not in range(1, 5)
        ):
            continue
        products.setdefault(gaia, {}).setdefault((tic, camera, ccd), set()).add(orbit)
    selection: dict[tuple[str, str], tuple[int, int]] = {}
    n_no_complete = 0
    n_multiple_complete = 0
    required = set(expected_orbits)
    for gaia, candidates in products.items():
        complete = [
            product for product, orbits in candidates.items() if orbits == required
        ]
        if not complete:
            n_no_complete += 1
            continue
        if len(complete) != 1:
            n_multiple_complete += 1
            continue
        tic, camera, ccd = complete[0]
        key = (gaia, tic)
        if key in selection and selection[key] != (camera, ccd):
            raise FM0ContractError(
                f"observation inventory selects duplicate product for {key}"
            )
        selection[key] = (camera, ccd)
    if not selection:
        raise FM0ContractError(
            f"observation inventory selects no complete S{sector} products"
        )
    return selection, {
        "n_target_inventory_sector_rows": n_sector_rows,
        "n_target_inventory_edge_rows_excluded": n_edge_rows,
        "n_target_inventory_no_complete_products": n_no_complete,
        "n_target_inventory_multiple_complete_products": n_multiple_complete,
    }


def _selection_from_corpus_rows(
    rows: Sequence[Mapping[str, Any]], *, sector: int, expected_orbits: tuple[int, int]
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, int]]:
    selection: dict[tuple[str, str], dict[str, Any]] = {}
    selected_tic_to_gaia: dict[str, str] = {}
    selected_gaia_to_tic: dict[str, str] = {}
    n_sector_rows = 0
    for row in rows:
        try:
            row_sector = int(row["sector"])
        except (KeyError, TypeError, ValueError) as exc:
            raise FM0ContractError(
                "corpus selection contains an invalid sector"
            ) from exc
        if row_sector != sector:
            continue
        n_sector_rows += 1
        if row.get("schema_version") != CORPUS_SELECTION_SCHEMA_VERSION:
            raise FM0ContractError("corpus-selection schema version drifted")
        gaia = _positive_identifier(
            row.get("gaia_dr3_source_id"), label="gaia_dr3_source_id"
        )
        tic = _positive_identifier(row.get("tic_id"), label="tic_id")
        component = str(row.get("leakage_component_id", "")).strip()
        partition = str(row.get("source_partition", "")).strip()
        if _LEAKAGE_COMPONENT.fullmatch(component) is None:
            raise FM0ContractError(
                "corpus selection has an invalid leakage_component_id"
            )
        if partition not in _SOURCE_PARTITIONS:
            raise FM0ContractError("corpus selection has an invalid source_partition")
        expected_partition, _bucket = deterministic_source_partition(component)
        if partition != expected_partition:
            raise FM0ContractError(
                "corpus-selection source_partition differs from the frozen component split"
            )
        try:
            camera = int(row["camera"])
            ccd = int(row["ccd"])
            orbits = tuple(int(value) for value in json.loads(str(row["orbits_json"])))
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise FM0ContractError(
                "corpus selection has invalid detector/orbit fields"
            ) from exc
        if (
            camera not in range(1, 5)
            or ccd not in range(1, 5)
            or orbits != expected_orbits
        ):
            raise FM0ContractError(
                "corpus selection is not one complete two-orbit product"
            )
        try:
            hdf5_paths = json.loads(str(row["hdf5_paths_json"]))
        except (KeyError, TypeError, json.JSONDecodeError) as exc:
            raise FM0ContractError(
                "corpus selection has invalid hdf5_paths_json"
            ) from exc
        expected_paths = [
            f"orbit-{orbit}/ffi/cam{camera}/ccd{ccd}/LC/{tic}.h5"
            for orbit in expected_orbits
        ]
        if hdf5_paths != expected_paths:
            raise FM0ContractError(
                "corpus-selection hdf5_paths_json differs from its exact "
                "sector/orbit/detector/TIC identity"
            )
        key = (gaia, tic)
        if key in selection:
            raise FM0ContractError(f"duplicate corpus-selection identity: {key}")
        prior_tic = selected_gaia_to_tic.setdefault(gaia, tic)
        prior_gaia = selected_tic_to_gaia.setdefault(tic, gaia)
        if prior_tic != tic or prior_gaia != gaia:
            raise FM0ContractError(
                "corpus selection has an ambiguous selected Gaia--TIC mapping"
            )
        selection[key] = {
            "camera": camera,
            "ccd": ccd,
            "leakage_component_id": component,
            "source_partition": partition,
        }
    if not selection:
        raise FM0ContractError(f"corpus selection contains no S{sector} products")
    return selection, {
        "n_target_inventory_sector_rows": n_sector_rows,
        "n_target_inventory_edge_rows_excluded": 0,
        "n_target_inventory_no_complete_products": 0,
        "n_target_inventory_multiple_complete_products": 0,
    }


def _load_target_selection(
    *,
    path: Path,
    sector: int,
    expected_orbits: tuple[int, int],
    expected_sha256: str,
) -> tuple[dict[tuple[str, str], dict[str, Any]], str, dict[str, int]]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError(f"target inventory must be a materialized file: {path}")
    if sha256_file(path) != expected_sha256:
        raise FM0ContractError("authoritative target-inventory SHA-256 drifted")
    lower_name = path.name.lower()
    if lower_name.endswith((".fits", ".fit", ".fits.gz")):
        # Validate that this really is the geometry authority, then stop.  It
        # cannot bind leakage components or the frozen split, so publishing
        # source rows from it would risk exposing sealed products downstream.
        list(iter_observation_fits(path, sectors=(sector,)))
        raise FM0ContractError(
            "observation inventory is metadata-only; source publication requires "
            "the exact checksum-bound corpus-selection CSV"
        )

    rows, columns = _inventory_rows(path)
    if tuple(columns) == tuple(CORPUS_SELECTION_FIELDS):
        selection, audit = _selection_from_corpus_rows(
            rows, sector=sector, expected_orbits=expected_orbits
        )
        return selection, IDENTITY_SOURCE_CORPUS_SELECTION, audit
    required_observation = {
        "source_id",
        "tic_id",
        "sector",
        "orbit",
        "camera",
        "ccd",
        "edge_warn",
    }
    if required_observation.issubset(columns):
        _selection_from_observations(
            rows, sector=sector, expected_orbits=expected_orbits
        )
        raise FM0ContractError(
            "observation inventory is metadata-only; source publication requires "
            "the exact checksum-bound corpus-selection CSV"
        )
    raise FM0ContractError(
        "target inventory must be the exact frozen corpus-selection schema or "
        "an observation inventory with source_id/tic_id/sector/orbit/camera/ccd/edge_warn"
    )


def _summary_authority(
    summary: Mapping[str, Any], *, key: str, label: str
) -> tuple[str, str]:
    authorities = summary.get("input_authorities")
    if not isinstance(authorities, Mapping):
        raise FM0ContractError("corpus summary has no input_authorities object")
    binding = authorities.get(key)
    if not isinstance(binding, Mapping):
        raise FM0ContractError(f"corpus summary has no {label} authority binding")
    path = str(binding.get("path", "")).strip()
    digest = _sha256_text(binding.get("sha256"), label=f"{label} authority SHA-256")
    if not path:
        raise FM0ContractError(f"corpus summary {label} authority path is empty")
    return path, digest


def _verify_corpus_summary(
    *,
    path: Path,
    expected_sha256: str,
    selection_path: Path,
    selection_sha256: str,
    expected_target_authority_sha256: str,
    sector: int,
    gaia_tic_alias_authority_path: Path | None,
    gaia_tic_alias_authority_sha256: str | None,
) -> tuple[Mapping[str, Any], str, str, bool]:
    expected_summary_sha = _sha256_text(expected_sha256, label="corpus-summary SHA-256")
    summary = _read_json(path, label="corpus summary")
    if sha256_file(path) != expected_summary_sha:
        raise FM0ContractError("corpus-summary SHA-256 drifted")
    try:
        sectors = tuple(int(value) for value in summary.get("sectors", ()))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError("corpus summary has invalid sectors") from exc
    if (
        summary.get("schema_version") != CORPUS_SELECTION_SCHEMA_VERSION
        or sectors != tuple(range(66, 78))
        or sector not in sectors
    ):
        raise FM0ContractError(
            "corpus summary must be the exact S66--S77 corpus-selection authority"
        )
    outputs = summary.get("outputs_sha256")
    if not isinstance(outputs, Mapping):
        raise FM0ContractError("corpus summary has no outputs_sha256 object")
    if (
        selection_path.name != "selection.csv"
        or _sha256_text(
            outputs.get("selection.csv"), label="corpus-summary selection.csv SHA-256"
        )
        != selection_sha256
    ):
        raise FM0ContractError(
            "corpus summary does not bind the supplied selection.csv"
        )
    _observation_path, observation_sha = _summary_authority(
        summary, key="observation_fits", label="observation FITS"
    )
    if observation_sha != expected_target_authority_sha256:
        raise FM0ContractError(
            "corpus summary observation authority differs from the expected target authority"
        )
    summary_alias_path, summary_alias_sha = _summary_authority(
        summary, key="gaia_tic_alias_table", label="Gaia--TIC alias table"
    )
    if (gaia_tic_alias_authority_path is None) != (
        gaia_tic_alias_authority_sha256 is None
    ):
        raise FM0ContractError(
            "Gaia--TIC alias authority path and SHA-256 must be supplied together"
        )
    alias_verified = False
    alias_path = summary_alias_path
    if gaia_tic_alias_authority_path is not None:
        supplied_path = gaia_tic_alias_authority_path.expanduser().resolve()
        supplied_sha = _sha256_text(
            gaia_tic_alias_authority_sha256,
            label="supplied Gaia--TIC alias authority SHA-256",
        )
        if supplied_path.is_symlink() or not supplied_path.is_file():
            raise FM0ContractError(
                f"Gaia--TIC alias authority must be materialized: {supplied_path}"
            )
        if (
            supplied_sha != summary_alias_sha
            or sha256_file(supplied_path) != supplied_sha
        ):
            raise FM0ContractError(
                "supplied Gaia--TIC alias authority differs from the corpus summary"
            )
        alias_verified = True
        alias_path = str(supplied_path)
    return summary, alias_path, summary_alias_sha, alias_verified


def _verify_alias_selection(
    *,
    path: Path,
    selection: Mapping[tuple[str, str], Mapping[str, Any]],
) -> None:
    rows, columns = _inventory_rows(path)
    required = {"gaia_dr3_source_id", "tic_id"}
    if not required.issubset(columns):
        raise FM0ContractError("Gaia--TIC alias authority lacks required columns")
    registry = build_alias_registry(
        {
            "gaia_dr3_source_id": row["gaia_dr3_source_id"],
            "tic_id": row["tic_id"],
        }
        for row in rows
    )
    index = registry.alias_index()
    for (gaia, tic), selected in selection.items():
        alias = index.get((gaia, tic))
        if (
            alias is None
            or bool(alias["quarantined"])
            or alias["leakage_component_id"] != selected["leakage_component_id"]
            or alias["source_partition"] != selected["source_partition"]
        ):
            raise FM0ContractError(
                "corpus selection differs from the supplied Gaia--TIC alias authority"
            )


def _binding_index(
    receipt: Mapping[str, Any], reconstructed: Mapping[str, Any]
) -> dict[str, Mapping[str, Any]]:
    for field in (
        "schema_version",
        "sector",
        "expected_orbits",
        "source_form",
        "source_state",
        "product_state",
        "pdo_sector_accepted",
        "n_cells",
        "n_hdf5_products_declared",
    ):
        if receipt.get(field) != reconstructed.get(field):
            raise FM0ContractError(f"source receipt field drifted: {field}")
    raw_bindings = receipt.get("cell_bindings")
    expected_bindings = reconstructed.get("cell_bindings")
    if not isinstance(raw_bindings, list) or not isinstance(expected_bindings, list):
        raise FM0ContractError("source receipt lacks cell bindings")
    raw_index = {
        str(binding.get("cell")): binding
        for binding in raw_bindings
        if isinstance(binding, Mapping)
    }
    expected_index = {
        str(binding.get("cell")): binding
        for binding in expected_bindings
        if isinstance(binding, Mapping)
    }
    if len(raw_index) != len(raw_bindings) or set(raw_index) != set(expected_index):
        raise FM0ContractError("source receipt cell-binding identities drifted")
    required = (
        "completion_receipt",
        "completion_receipt_sha256",
        "retention_receipt",
        "retention_receipt_sha256",
        "retained_root",
        "n_hdf5_products_declared",
    )
    for cell, expected in expected_index.items():
        raw = raw_index[cell]
        for field in required:
            if raw.get(field) != expected.get(field):
                raise FM0ContractError(
                    f"source receipt {cell} binding drifted: {field}"
                )
    return raw_index


def _cell_inventory(
    *,
    source_root: Path,
    binding: Mapping[str, Any],
    sector: int,
    expected_target_authority_sha256: str,
) -> dict[str, Any]:
    cell = str(binding.get("cell", ""))
    match = _CELL.fullmatch(cell)
    if match is None or int(match.group("sector")) != sector:
        raise FM0ContractError(f"invalid retained cell identity: {cell!r}")
    orbit = int(match.group("orbit"))
    camera = int(match.group("camera"))
    ccd = int(match.group("ccd"))
    retained_root = Path(str(binding["retained_root"])).resolve()
    completion_path = Path(str(binding["completion_receipt"])).resolve()
    completion = _read_json(completion_path, label=f"{cell} completion receipt")

    retained_manifest_path = retained_root / "cell_manifest.json"
    retained_manifest_sha = sha256_file(retained_manifest_path)
    if retained_manifest_sha != _sha256_text(
        completion.get("input_manifest_sha256"), label=f"{cell} input manifest SHA-256"
    ):
        raise FM0ContractError(f"{cell} retained cell manifest is not completion-bound")
    input_root = source_root / "inputs" / f"s{sector:04d}-o{orbit}-cam{camera}-ccd{ccd}"
    if input_root.is_symlink() or not input_root.is_dir():
        raise FM0ContractError(f"{cell} input root is missing: {input_root}")
    input_manifest_path = input_root / "cell_manifest.json"
    if sha256_file(input_manifest_path) != retained_manifest_sha:
        raise FM0ContractError(f"{cell} staged and retained cell manifests differ")
    for ready in (input_root / "READY", retained_root / "READY"):
        if ready.is_symlink() or not ready.is_file():
            raise FM0ContractError(f"{cell} READY authority is missing: {ready}")
        if ready.read_text(encoding="utf-8").strip() != retained_manifest_sha:
            raise FM0ContractError(f"{cell} READY authority drifted: {ready}")

    manifest = _read_json(input_manifest_path, label=f"{cell} input manifest")
    expected_cell = {
        "sector": sector,
        "orbit": orbit,
        "orbit_tag": "o1" if orbit == min(2 * sector + 7, 2 * sector + 8) else "o2",
        "camera": camera,
        "ccd": ccd,
    }
    manifest_cell = manifest.get("cell")
    if (
        manifest.get("schema_version") != 1
        or not isinstance(manifest_cell, Mapping)
        or any(manifest_cell.get(key) != value for key, value in expected_cell.items())
    ):
        raise FM0ContractError(f"{cell} input-manifest identity drifted")
    if (
        _sha256_text(
            manifest.get("target_authority_sha256"), label=f"{cell} target authority"
        )
        != expected_target_authority_sha256
    ):
        raise FM0ContractError(
            f"{cell} target authority differs from requested authority"
        )
    records = _manifest_records(manifest)
    overlay_records = {
        relative: spec
        for relative, spec in records.items()
        if relative.startswith("source_tic/")
    }
    source_records = [
        relative for relative in records if relative.startswith("source/")
    ]
    metadata_records = {
        relative: spec
        for relative, spec in records.items()
        if relative.startswith("metadata/")
    }
    required_metadata = {
        "metadata/requested_tics.txt",
        "metadata/expected_h5_tics.txt",
    }
    if (
        len(source_records) != EXPECTED_SOURCE_GRID_FILES
        or len(overlay_records) != EXPECTED_SOURCE_GRID_FILES
        or set(metadata_records) != required_metadata
        or set(records)
        != set(source_records) | set(overlay_records) | required_metadata
    ):
        raise FM0ContractError(f"{cell} input-manifest payload contract drifted")
    if any(
        _OVERLAY.fullmatch(Path(relative).name) is None for relative in overlay_records
    ):
        raise FM0ContractError(f"{cell} contains an invalid source_tic filename")
    counts = manifest.get("counts")
    if counts != {
        "source": EXPECTED_SOURCE_GRID_FILES,
        "source_tic": EXPECTED_SOURCE_GRID_FILES,
        "metadata": 2,
        "total": (2 * EXPECTED_SOURCE_GRID_FILES) + 2,
    }:
        raise FM0ContractError(f"{cell} input-manifest counts drifted")

    requested_spec = metadata_records["metadata/requested_tics.txt"]
    expected_spec = metadata_records["metadata/expected_h5_tics.txt"]
    requested_path = _verify_manifest_file(
        root=input_root,
        relative="metadata/requested_tics.txt",
        expected_size=requested_spec[0],
        expected_sha256=requested_spec[1],
    )
    expected_path = _verify_manifest_file(
        root=input_root,
        relative="metadata/expected_h5_tics.txt",
        expected_size=expected_spec[0],
        expected_sha256=expected_spec[1],
    )
    requested = set(_load_tics(requested_path, expected_sha256=requested_spec[1]))
    expected = set(_load_tics(expected_path, expected_sha256=expected_spec[1]))
    if not expected.issubset(requested):
        raise FM0ContractError(f"{cell} expected HDF5 TICs are not requested")
    # The retained metadata copies are part of the immutable returned bundle.
    for relative, spec in metadata_records.items():
        _verify_manifest_file(
            root=retained_root,
            relative=relative,
            expected_size=spec[0],
            expected_sha256=spec[1],
        )

    pairs: set[tuple[str, str]] = set()
    for relative, (size, digest) in sorted(overlay_records.items()):
        path = _verify_manifest_file(
            root=input_root,
            relative=relative,
            expected_size=size,
            expected_sha256=digest,
        )
        current = _overlay_pairs(path)
        unrequested = sorted(
            {int(tic) for _, tic in current if int(tic) not in requested}
        )
        if unrequested:
            raise FM0ContractError(
                f"{cell} source_tic overlay contains TICs outside the "
                f"manifest-bound requested authority: {unrequested[:10]}"
            )
        pairs.update(current)

    output_manifest_path = retained_root / "output_manifest.sha256"
    output_manifest_sha = sha256_file(output_manifest_path)
    if output_manifest_sha != _sha256_text(
        completion.get("output_manifest_sha256"),
        label=f"{cell} output manifest SHA-256",
    ):
        raise FM0ContractError(f"{cell} output manifest is not completion-bound")
    output_records = _output_manifest(output_manifest_path)
    actual_paths = _actual_output_paths(retained_root)
    if set(output_records) != actual_paths:
        missing = sorted(set(output_records) - actual_paths)
        extra = sorted(actual_paths - set(output_records))
        raise FM0ContractError(
            f"{cell} output manifest paths drifted; missing={missing[:5]}, extra={extra[:5]}"
        )
    hdf5: list[dict[str, Any]] = []
    actual_tics: set[int] = set()
    for relative, digest in sorted(output_records.items()):
        if not relative.startswith("LC/"):
            continue
        path = retained_root / relative
        tic_text = path.stem
        tic = int(_positive_identifier(tic_text, label=f"HDF5 TIC in {path}"))
        if tic in actual_tics:
            raise FM0ContractError(f"{cell} has duplicate HDF5 TIC {tic}")
        actual_tics.add(tic)
        hdf5.append(
            {
                "tic_id": str(tic),
                "orbit": orbit,
                "camera": camera,
                "ccd": ccd,
                "cell": cell,
                "hdf5_path": str(path),
                "hdf5_sha256": digest,
                "cell_manifest_sha256": retained_manifest_sha,
                "output_manifest_path": str(output_manifest_path),
                "output_manifest_sha256": output_manifest_sha,
            }
        )
    if len(hdf5) != int(binding.get("n_hdf5_products_declared", -1)):
        raise FM0ContractError(f"{cell} HDF5 count differs from the source receipt")
    if not expected.issubset(actual_tics) or not actual_tics.issubset(requested):
        raise FM0ContractError(f"{cell} HDF5 outputs violate their TIC authorities")
    return {
        "cell": cell,
        "orbit": orbit,
        "camera": camera,
        "ccd": ccd,
        "requested_tics": requested,
        "overlay_pairs": pairs,
        "hdf5": hdf5,
    }


def _verify_source_receipt(
    *, source_receipt_path: Path, sector: int, expected_orbits: tuple[int, int]
) -> tuple[Mapping[str, Any], Path, dict[str, Mapping[str, Any]]]:
    source_path = source_receipt_path.expanduser().resolve()
    receipt = _read_json(source_path, label="ORCD source receipt")
    if (
        receipt.get("schema_version") != SOURCE_RECEIPT_SCHEMA_VERSION
        or receipt.get("source_state") != SOURCE_READY_STATE
        or receipt.get("source_form") != "retained_orcd_cells"
        or receipt.get("product_state") != "ORCD_COMPLETE_DEFERRED"
        or receipt.get("pdo_sector_accepted") is not False
        or int(receipt.get("sector", -1)) != sector
        or tuple(int(value) for value in receipt.get("expected_orbits", ()))
        != expected_orbits
        or int(receipt.get("n_cells", -1)) != 32
    ):
        raise FM0ContractError(
            f"S{sector} source receipt is not a retained deferred sector"
        )
    source_root = Path(str(receipt.get("source_root", ""))).expanduser().resolve()
    reconstructed = verify_retained_sector_source(
        sector_root=source_root,
        sector=sector,
        expected_orbits=expected_orbits,
    )
    return receipt, source_root, _binding_index(receipt, reconstructed)


def build_later_source_inventory(
    *,
    source_receipt_path: str | Path,
    source_receipt_sha256: str,
    sector: int,
    expected_orbits: Sequence[int],
    expected_target_authority_sha256: str,
    output_dir: str | Path,
    producer_git_sha: str,
    target_inventory_path: str | Path,
    target_inventory_sha256: str,
    corpus_summary_path: str | Path,
    corpus_summary_sha256: str,
    gaia_tic_alias_authority_path: str | Path | None = None,
    gaia_tic_alias_authority_sha256: str | None = None,
) -> dict[str, Any]:
    """Verify and publish one later-sector retained HDF5 identity inventory.

    The target inventory must use the exact frozen corpus-selection CSV
    schema.  Manifest-bound ``source_tic`` overlays cross-check the selected
    Gaia--TIC identities but never define membership.  Observation inventories
    are rejected here because they cannot bind leakage components or the
    frozen source partition required to keep sealed products closed.
    """

    sector = int(sector)
    orbits = tuple(int(value) for value in expected_orbits)
    if sector not in range(66, 78) or len(orbits) != 2 or len(set(orbits)) != 2:
        raise FM0ContractError(
            "later-source inventory is restricted to S66--S77 and two distinct orbits"
        )
    if orbits != (2 * sector + 7, 2 * sector + 8):
        raise FM0ContractError(f"S{sector} expected orbit identity drifted")
    authority_sha = _sha256_text(
        expected_target_authority_sha256, label="expected target-authority SHA-256"
    )
    if _GIT_SHA.fullmatch(str(producer_git_sha)) is None:
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    inventory_path = Path(target_inventory_path).expanduser().resolve()
    inventory_sha = _sha256_text(
        target_inventory_sha256, label="target-inventory SHA-256"
    )
    expected_source_receipt_sha = _sha256_text(
        source_receipt_sha256, label="source-receipt SHA-256"
    )
    source_path = Path(source_receipt_path).expanduser().resolve()
    if source_path.is_symlink() or not source_path.is_file():
        raise FM0ContractError(f"source receipt must be materialized: {source_path}")
    if sha256_file(source_path) != expected_source_receipt_sha:
        raise FM0ContractError("source-receipt SHA-256 drifted")
    summary_path = Path(corpus_summary_path).expanduser().resolve()
    supplied_alias_path = (
        None
        if gaia_tic_alias_authority_path is None
        else Path(gaia_tic_alias_authority_path)
    )
    corpus_summary, alias_authority_path, alias_authority_sha, alias_verified = (
        _verify_corpus_summary(
            path=summary_path,
            expected_sha256=corpus_summary_sha256,
            selection_path=inventory_path,
            selection_sha256=inventory_sha,
            expected_target_authority_sha256=authority_sha,
            sector=sector,
            gaia_tic_alias_authority_path=supplied_alias_path,
            gaia_tic_alias_authority_sha256=gaia_tic_alias_authority_sha256,
        )
    )
    selection, identity_source, selection_audit = _load_target_selection(
        path=inventory_path,
        sector=sector,
        expected_orbits=orbits,
        expected_sha256=inventory_sha,
    )
    observations_by_sector = corpus_summary.get("observations_by_sector")
    if not isinstance(observations_by_sector, Mapping) or int(
        observations_by_sector.get(str(sector), -1)
    ) != len(selection):
        raise FM0ContractError(
            f"corpus summary S{sector} observation count differs from selection.csv"
        )
    if alias_verified:
        _verify_alias_selection(path=Path(alias_authority_path), selection=selection)

    receipt, source_root, bindings = _verify_source_receipt(
        source_receipt_path=source_path,
        sector=sector,
        expected_orbits=orbits,
    )
    source_receipt_sha = sha256_file(source_path)
    if source_receipt_sha != expected_source_receipt_sha:
        raise FM0ContractError("source receipt changed during inventory construction")
    cells = [
        _cell_inventory(
            source_root=source_root,
            binding=bindings[cell],
            sector=sector,
            expected_target_authority_sha256=authority_sha,
        )
        for cell in sorted(bindings)
    ]
    if {cell["orbit"] for cell in cells} != set(orbits) or len(cells) != 32:
        raise FM0ContractError(
            f"S{sector} cell inventory does not cover both full orbits"
        )
    identity_binding_path = str(inventory_path)
    identity_binding_sha = inventory_sha
    verified_corpus_summary_sha = _sha256_text(
        corpus_summary_sha256, label="corpus-summary SHA-256"
    )

    cell_by_detector = {
        (int(cell["orbit"]), int(cell["camera"]), int(cell["ccd"])): cell
        for cell in cells
    }
    if len(cell_by_detector) != 32:
        raise FM0ContractError(f"S{sector} has duplicate retained detector cells")
    selected_tics = {tic for _, tic in selection}
    overlay_pairs = set().union(*(cell["overlay_pairs"] for cell in cells))
    pairs = {pair for pair in overlay_pairs if pair[1] in selected_tics}
    for (gaia, tic), selected in selection.items():
        camera = int(selected["camera"])
        ccd = int(selected["ccd"])
        for orbit in orbits:
            cell = cell_by_detector[(orbit, camera, ccd)]
            if (gaia, tic) not in cell["overlay_pairs"]:
                raise FM0ContractError(
                    f"selected Gaia--TIC {gaia}--{tic} is absent from the "
                    f"manifest-bound overlay in {cell['cell']}"
                )
        pairs.add((gaia, tic))
    if not pairs:
        raise FM0ContractError(f"S{sector} selected identities have no overlay closure")

    hdf5_index: dict[tuple[str, int, int, int], dict[str, Any]] = {}
    all_hdf5_paths: set[str] = set()
    for cell in cells:
        for hdf5 in cell["hdf5"]:
            path = str(hdf5["hdf5_path"])
            if path in all_hdf5_paths:
                raise FM0ContractError(f"HDF5 path appears in multiple cells: {path}")
            all_hdf5_paths.add(path)
            key = (
                str(hdf5["tic_id"]),
                int(hdf5["orbit"]),
                int(hdf5["camera"]),
                int(hdf5["ccd"]),
            )
            if key in hdf5_index:
                raise FM0ContractError(
                    f"duplicate retained HDF5 product identity: {key}"
                )
            hdf5_index[key] = hdf5
    selected_products: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for pair, selected in selection.items():
        gaia, tic = pair
        camera = int(selected["camera"])
        ccd = int(selected["ccd"])
        products: list[dict[str, Any]] = []
        for orbit in orbits:
            product = hdf5_index.get((tic, orbit, camera, ccd))
            if product is None:
                raise FM0ContractError(
                    f"selected complete product {gaia}--{tic} lacks retained "
                    f"S{sector}/orbit-{orbit}/cam{camera}/ccd{ccd} HDF5"
                )
            products.append(product)
        selected_products[pair] = products
    pairs_by_tic: dict[str, set[str]] = {}
    for gaia, tic in pairs:
        pairs_by_tic.setdefault(tic, set()).add(gaia)
    for gaia, tic in selection:
        if pairs_by_tic.get(tic) != {gaia}:
            raise FM0ContractError(
                "selected TIC has an ambiguous Gaia identity in the "
                "manifest-bound source_tic overlays"
            )

    alias_rows: list[dict[str, Any]] = []
    for gaia, tic in sorted(pairs, key=lambda item: (int(item[0]), int(item[1]))):
        ambiguous = len(pairs_by_tic[tic]) != 1
        alias_rows.append(
            {
                "alias_schema_version": ALIAS_SCHEMA_VERSION,
                "sector": sector,
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "identity_source": ALIAS_IDENTITY_SOURCE,
                "target_authority_sha256": authority_sha,
                "identity_ambiguous": str(ambiguous).lower(),
                "has_retained_hdf5": str(tic in selected_tics).lower(),
            }
        )

    source_rows: list[dict[str, Any]] = []
    selected_paths: set[str] = set()
    for (gaia, tic), selected in sorted(
        selection.items(), key=lambda item: (int(item[0][0]), int(item[0][1]))
    ):
        products = selected_products[(gaia, tic)]
        selected_paths.update(str(product["hdf5_path"]) for product in products)
        portable_evidence = [
            {
                "orbit": int(product["orbit"]),
                "camera": int(product["camera"]),
                "ccd": int(product["ccd"]),
                "cell": str(product["cell"]),
                "hdf5_sha256": str(product["hdf5_sha256"]),
                "cell_manifest_sha256": str(product["cell_manifest_sha256"]),
                "output_manifest_sha256": str(product["output_manifest_sha256"]),
            }
            for product in products
        ]
        retained_source_sha = _canonical_sha256(
            {
                "schema_version": INVENTORY_SCHEMA_VERSION,
                "sector": sector,
                "tic_id": tic,
                "target_authority_sha256": authority_sha,
                "evidence": portable_evidence,
            }
        )
        full_evidence = [
            {
                **portable,
                "hdf5_path": str(product["hdf5_path"]),
                "output_manifest_path": str(product["output_manifest_path"]),
            }
            for portable, product in zip(portable_evidence, products, strict=True)
        ]
        source_row_sha = _canonical_sha256(
            {
                "retained_source_sha256": retained_source_sha,
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "leakage_component_id": selected["leakage_component_id"],
                "source_partition": selected["source_partition"],
                "identity_source": identity_source,
                "identity_binding_sha256": identity_binding_sha,
                "corpus_summary_sha256": verified_corpus_summary_sha,
                "gaia_tic_alias_authority_sha256": alias_authority_sha,
            }
        )
        source_rows.append(
            {
                "source_row_schema_version": SOURCE_ROW_SCHEMA_VERSION,
                "sector": sector,
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "leakage_component_id": selected["leakage_component_id"],
                "source_partition": selected["source_partition"],
                "identity_source": identity_source,
                "target_authority_sha256": authority_sha,
                "corpus_selection_sha256": inventory_sha,
                "corpus_summary_sha256": verified_corpus_summary_sha,
                "gaia_tic_alias_authority_sha256": alias_authority_sha,
                "identity_ambiguous": str(len(pairs_by_tic[tic]) != 1).lower(),
                "product_state": "ORCD_COMPLETE_DEFERRED",
                "n_hdf5_products": len(products),
                "orbits_json": json.dumps(
                    [int(product["orbit"]) for product in products],
                    separators=(",", ":"),
                ),
                "hdf5_paths_json": json.dumps(
                    [str(product["hdf5_path"]) for product in products],
                    separators=(",", ":"),
                ),
                "hdf5_sha256_json": json.dumps(
                    [str(product["hdf5_sha256"]) for product in products],
                    separators=(",", ":"),
                ),
                "hdf5_evidence_json": json.dumps(
                    full_evidence, sort_keys=True, separators=(",", ":")
                ),
                "retained_source_sha256": retained_source_sha,
                "source_row_sha256": source_row_sha,
                "source_receipt_path": str(source_path),
                "source_receipt_sha256": source_receipt_sha,
                "model_outcome_blind": "true",
                "panel_admission_authorized": "false",
            }
        )
    if not source_rows:
        raise FM0ContractError(f"S{sector} target authority selects no retained HDF5")
    source_rows.sort(
        key=lambda row: (int(row["gaia_dr3_source_id"]), int(row["tic_id"]))
    )
    if len(selected_paths) != 2 * len(source_rows):
        raise FM0ContractError(
            "selected HDF5 identity is not exactly two unique products per source row"
        )
    if sha256_file(source_path) != expected_source_receipt_sha:
        raise FM0ContractError("source receipt changed during inventory construction")
    if sha256_file(inventory_path) != inventory_sha:
        raise FM0ContractError("corpus selection changed during inventory construction")
    if sha256_file(summary_path) != verified_corpus_summary_sha:
        raise FM0ContractError("corpus summary changed during inventory construction")
    if (
        alias_verified
        and sha256_file(Path(alias_authority_path)) != alias_authority_sha
    ):
        raise FM0ContractError("Gaia--TIC alias authority changed during construction")
    assert_no_search_columns(ALIAS_FIELDS, context="later alias output")
    assert_no_search_columns(SOURCE_FIELDS, context="later source output")

    alias_payload = _csv_bytes(alias_rows, ALIAS_FIELDS)
    source_payload = _csv_bytes(source_rows, SOURCE_FIELDS)
    output = _output_leaf_path(output_dir, label="later-source output")
    partial = output.with_name(output.name + ".partial")
    if (
        output.exists()
        or output.is_symlink()
        or partial.exists()
        or partial.is_symlink()
    ):
        raise FM0ContractError(
            f"refusing to overwrite later-source inventory: {output}"
        )
    source_rows_by_partition = {
        partition: sum(row["source_partition"] == partition for row in source_rows)
        for partition in sorted(_SOURCE_PARTITIONS)
    }
    summary = {
        "summary_schema_version": SUMMARY_SCHEMA_VERSION,
        "inventory_schema_version": INVENTORY_SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "producer_git_sha": producer_git_sha,
        "producer_code_path": str(Path(__file__).resolve()),
        "producer_code_sha256": sha256_file(Path(__file__).resolve()),
        "sector": sector,
        "expected_orbits": list(orbits),
        "source_receipt_path": str(source_path),
        "source_receipt_sha256": source_receipt_sha,
        "alias_identity_source": ALIAS_IDENTITY_SOURCE,
        "selection_identity_source": identity_source,
        "target_authority_sha256": authority_sha,
        "identity_binding_path": identity_binding_path,
        "identity_binding_sha256": identity_binding_sha,
        "corpus_summary_path": str(summary_path),
        "corpus_summary_sha256": verified_corpus_summary_sha,
        "gaia_tic_alias_authority_path": alias_authority_path,
        "gaia_tic_alias_authority_sha256": alias_authority_sha,
        "gaia_tic_alias_authority_file_verified": alias_verified,
        "n_cells": len(cells),
        "n_hdf5_products_declared": int(receipt["n_hdf5_products_declared"]),
        "n_hdf5_products_selected_unique": len(selected_paths),
        "n_hdf5_products_excluded_not_in_target_inventory": len(
            all_hdf5_paths - selected_paths
        ),
        **selection_audit,
        "n_alias_rows": len(alias_rows),
        "n_source_rows": len(source_rows),
        "source_rows_by_partition": source_rows_by_partition,
        "n_unique_gaia_sources": len({row["gaia_dr3_source_id"] for row in alias_rows}),
        "n_unique_tic_aliases": len({row["tic_id"] for row in alias_rows}),
        "n_ambiguous_tic_aliases": sum(
            len(gaia_ids) != 1 for gaia_ids in pairs_by_tic.values()
        ),
        "aliases_sha256": hashlib.sha256(alias_payload).hexdigest(),
        "sources_sha256": hashlib.sha256(source_payload).hexdigest(),
        "hdf5_hash_authority": "completion_bound_output_manifest_sha256_v1",
        "hdf5_content_opened": False,
        "sealed_hdf5_content_opened": False,
        "model_outcome_blind": True,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "claim_limit": (
            "retained HDF5 provenance and Gaia--TIC identity inventory only; "
            "not accepted A2v1, a six-view release, a temporal panel, or a model result"
        ),
    }
    summary_payload = (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode(
        "utf-8"
    )
    summary_hash = hashlib.sha256(summary_payload).hexdigest()
    owns_partial = False
    owns_final = False
    try:
        partial.parent.mkdir(parents=True, exist_ok=True)
        partial.mkdir()
        owns_partial = True
        publish_immutable(partial / "aliases.csv", alias_payload)
        publish_immutable(partial / "sources.csv", source_payload)
        publish_immutable(partial / "summary.json", summary_payload)
        publish_immutable(
            partial / "READY", summary_hash.encode("ascii") + b"\n"
        )
        load_later_source_rows(
            partial,
            expected_summary_sha256=summary_hash,
            allowed_source_partitions=("poc_train", "poc_development"),
        )
        _make_output_tree_read_only(partial)
        _assert_output_tree_read_only(partial)
        os.replace(partial, output)
        owns_partial = False
        owns_final = True
        _assert_output_tree_read_only(output)
        return summary
    except BaseException:
        if owns_partial and not partial.exists() and output.exists():
            owns_partial = False
            owns_final = True
        cleanup_error: Exception | None = None
        for path, owned in ((partial, owns_partial), (output, owns_final)):
            if not owned:
                continue
            try:
                _remove_owned_output_tree(path)
            except (FM0ContractError, OSError) as exc:
                cleanup_error = exc
        if cleanup_error is not None:
            raise FM0ContractError(
                f"later-source build failed and cleanup also failed: {cleanup_error}"
            ) from cleanup_error
        raise


def _read_exact_csv(
    path: Path, *, fields: Sequence[str], label: str
) -> list[dict[str, str]]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError(f"{label} must be a nonempty materialized file: {path}")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != tuple(fields):
                raise FM0ContractError(f"{label} columns drifted")
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not rows or any(tuple(row) != tuple(fields) for row in rows):
        raise FM0ContractError(f"{label} is empty or has drifting row columns")
    return rows


def _json_list(value: Any, *, label: str) -> list[Any]:
    try:
        result = json.loads(str(value))
    except (TypeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"{label} is not valid JSON") from exc
    if not isinstance(result, list):
        raise FM0ContractError(f"{label} must be a JSON list")
    return result


def load_later_source_rows(
    inventory_dir: str | Path,
    *,
    expected_summary_sha256: str,
    allowed_source_partitions: Sequence[str] = ("poc_train", "poc_development"),
) -> tuple[dict[str, str], ...]:
    """Strictly load non-sealed source rows without opening any HDF5 content.

    The sealed partition is deliberately unavailable through this loader.
    Its rows remain in the immutable audit inventory, but a downstream
    six-view builder must filter them before it can obtain HDF5 paths.
    """

    allowed = tuple(str(value) for value in allowed_source_partitions)
    if (
        not allowed
        or len(set(allowed)) != len(allowed)
        or not set(allowed).issubset(_SOURCE_PARTITIONS)
        or "poc_sealed_test" in allowed
    ):
        raise FM0ContractError(
            "allowed_source_partitions must be a unique nonempty subset of "
            "poc_train/poc_development; sealed rows cannot be loaded"
        )
    root = Path(inventory_dir).expanduser().resolve()
    if root.is_symlink() or not root.is_dir():
        raise FM0ContractError(f"later-source inventory root is invalid: {root}")
    summary_path = root / "summary.json"
    expected_summary_sha = _sha256_text(
        expected_summary_sha256, label="later-source summary SHA-256"
    )
    summary = _read_json(summary_path, label="later-source summary")
    if sha256_file(summary_path) != expected_summary_sha:
        raise FM0ContractError("later-source summary SHA-256 drifted")
    ready = root / "READY"
    if ready.is_symlink() or not ready.is_file():
        raise FM0ContractError("later-source READY authority is missing")
    if ready.read_text(encoding="utf-8").strip() != expected_summary_sha:
        raise FM0ContractError("later-source READY authority drifted")
    try:
        sector = int(summary.get("sector", -1))
    except (TypeError, ValueError) as exc:
        raise FM0ContractError("later-source summary has invalid sector") from exc
    if (
        summary.get("summary_schema_version") != SUMMARY_SCHEMA_VERSION
        or summary.get("inventory_schema_version") != INVENTORY_SCHEMA_VERSION
        or sector not in range(66, 78)
        or summary.get("hdf5_content_opened") is not False
        or summary.get("sealed_hdf5_content_opened") is not False
        or summary.get("panel_admission_authorized") is not False
    ):
        raise FM0ContractError("later-source summary contract drifted")
    aliases_path = root / "aliases.csv"
    sources_path = root / "sources.csv"
    if sha256_file(aliases_path) != _sha256_text(
        summary.get("aliases_sha256"), label="aliases.csv SHA-256"
    ):
        raise FM0ContractError("later-source aliases.csv SHA-256 drifted")
    if sha256_file(sources_path) != _sha256_text(
        summary.get("sources_sha256"), label="sources.csv SHA-256"
    ):
        raise FM0ContractError("later-source sources.csv SHA-256 drifted")
    rows = _read_exact_csv(
        sources_path, fields=SOURCE_FIELDS, label="later-source sources.csv"
    )
    if len(rows) != int(summary.get("n_source_rows", -1)):
        raise FM0ContractError("later-source source-row count drifted")

    summary_selection_sha = _sha256_text(
        summary.get("identity_binding_sha256"), label="corpus selection SHA-256"
    )
    summary_corpus_sha = _sha256_text(
        summary.get("corpus_summary_sha256"), label="corpus summary SHA-256"
    )
    summary_alias_sha = _sha256_text(
        summary.get("gaia_tic_alias_authority_sha256"),
        label="Gaia--TIC alias authority SHA-256",
    )
    target_authority_sha = _sha256_text(
        summary.get("target_authority_sha256"), label="target authority SHA-256"
    )
    source_receipt_sha = _sha256_text(
        summary.get("source_receipt_sha256"), label="source receipt SHA-256"
    )
    expected_orbits = (2 * sector + 7, 2 * sector + 8)
    seen_gaia: set[str] = set()
    seen_tic: set[str] = set()
    seen_hdf5_paths: set[str] = set()
    counts_by_partition = {partition: 0 for partition in _SOURCE_PARTITIONS}
    selected: list[dict[str, str]] = []
    for row in rows:
        gaia = _positive_identifier(
            row["gaia_dr3_source_id"], label="source-row Gaia DR3 source ID"
        )
        tic = _positive_identifier(row["tic_id"], label="source-row TIC ID")
        component = str(row["leakage_component_id"])
        partition = str(row["source_partition"])
        if (
            row["source_row_schema_version"] != SOURCE_ROW_SCHEMA_VERSION
            or int(row["sector"]) != sector
            or _LEAKAGE_COMPONENT.fullmatch(component) is None
            or partition not in _SOURCE_PARTITIONS
            or deterministic_source_partition(component)[0] != partition
            or row["identity_source"] != IDENTITY_SOURCE_CORPUS_SELECTION
            or row["target_authority_sha256"] != target_authority_sha
            or row["corpus_selection_sha256"] != summary_selection_sha
            or row["corpus_summary_sha256"] != summary_corpus_sha
            or row["gaia_tic_alias_authority_sha256"] != summary_alias_sha
            or row["identity_ambiguous"] != "false"
            or row["product_state"] != "ORCD_COMPLETE_DEFERRED"
            or int(row["n_hdf5_products"]) != 2
            or row["source_receipt_sha256"] != source_receipt_sha
            or row["model_outcome_blind"] != "true"
            or row["panel_admission_authorized"] != "false"
        ):
            raise FM0ContractError("later-source source-row contract drifted")
        if gaia in seen_gaia or tic in seen_tic:
            raise FM0ContractError(
                "later-source rows have ambiguous Gaia--TIC identities"
            )
        seen_gaia.add(gaia)
        seen_tic.add(tic)
        orbits = tuple(
            int(value) for value in _json_list(row["orbits_json"], label="orbits_json")
        )
        paths = [
            str(value)
            for value in _json_list(row["hdf5_paths_json"], label="hdf5_paths_json")
        ]
        hashes = [
            _sha256_text(value, label="HDF5 SHA-256")
            for value in _json_list(row["hdf5_sha256_json"], label="hdf5_sha256_json")
        ]
        evidence = _json_list(row["hdf5_evidence_json"], label="hdf5_evidence_json")
        if (
            orbits != expected_orbits
            or len(paths) != 2
            or len(hashes) != 2
            or len(evidence) != 2
        ):
            raise FM0ContractError(
                "later-source row is not one complete two-orbit product"
            )
        portable_evidence: list[dict[str, Any]] = []
        for index, (orbit, hdf5_path, hdf5_sha, raw_evidence) in enumerate(
            zip(orbits, paths, hashes, evidence, strict=True)
        ):
            if not isinstance(raw_evidence, Mapping):
                raise FM0ContractError("later-source HDF5 evidence row is invalid")
            path = Path(hdf5_path)
            match = _CELL.fullmatch(str(raw_evidence.get("cell", "")))
            if (
                not path.is_absolute()
                or tuple(path.parts[-2:]) != ("LC", f"{tic}.h5")
                or hdf5_path in seen_hdf5_paths
                or match is None
                or int(match.group("sector")) != sector
                or int(match.group("orbit")) != orbit
                or int(raw_evidence.get("orbit", -1)) != orbit
                or str(raw_evidence.get("hdf5_path")) != hdf5_path
                or str(raw_evidence.get("hdf5_sha256")) != hdf5_sha
                or index != orbits.index(orbit)
            ):
                raise FM0ContractError("later-source HDF5 evidence identity drifted")
            seen_hdf5_paths.add(hdf5_path)
            portable = {
                key: raw_evidence[key]
                for key in (
                    "orbit",
                    "camera",
                    "ccd",
                    "cell",
                    "hdf5_sha256",
                    "cell_manifest_sha256",
                    "output_manifest_sha256",
                )
            }
            for digest_field in (
                "hdf5_sha256",
                "cell_manifest_sha256",
                "output_manifest_sha256",
            ):
                _sha256_text(portable[digest_field], label=digest_field)
            portable_evidence.append(portable)
        retained_source_sha = _canonical_sha256(
            {
                "schema_version": INVENTORY_SCHEMA_VERSION,
                "sector": sector,
                "tic_id": tic,
                "target_authority_sha256": target_authority_sha,
                "evidence": portable_evidence,
            }
        )
        if row["retained_source_sha256"] != retained_source_sha:
            raise FM0ContractError("later-source retained_source_sha256 drifted")
        source_row_sha = _canonical_sha256(
            {
                "retained_source_sha256": retained_source_sha,
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "leakage_component_id": component,
                "source_partition": partition,
                "identity_source": IDENTITY_SOURCE_CORPUS_SELECTION,
                "identity_binding_sha256": summary_selection_sha,
                "corpus_summary_sha256": summary_corpus_sha,
                "gaia_tic_alias_authority_sha256": summary_alias_sha,
            }
        )
        if row["source_row_sha256"] != source_row_sha:
            raise FM0ContractError("later-source source_row_sha256 drifted")
        counts_by_partition[partition] += 1
        if partition in allowed:
            selected.append(row)
    if len(seen_hdf5_paths) != 2 * len(rows):
        raise FM0ContractError("later-source HDF5 paths are not unique two-per-row")
    summary_counts = summary.get("source_rows_by_partition")
    if not isinstance(summary_counts, Mapping) or any(
        int(summary_counts.get(partition, -1)) != count
        for partition, count in counts_by_partition.items()
    ):
        raise FM0ContractError("later-source partition counts drifted")
    return tuple(selected)


__all__ = [
    "ALIAS_FIELDS",
    "ALIAS_IDENTITY_SOURCE",
    "ALIAS_SCHEMA_VERSION",
    "IDENTITY_SOURCE_CORPUS_SELECTION",
    "INVENTORY_SCHEMA_VERSION",
    "SOURCE_FIELDS",
    "SOURCE_ROW_SCHEMA_VERSION",
    "SUMMARY_SCHEMA_VERSION",
    "build_later_source_inventory",
    "load_later_source_rows",
]
