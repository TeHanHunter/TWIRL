from __future__ import annotations

import importlib.metadata
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict

import numpy as np
from astropy.table import Column, Table


FIRST_200S_SECTOR = 56
FIRST_TESS_SECTOR = 1
FIRST_TESS_ORBIT = 9
DEFAULT_ORBITS_PER_SECTOR = 2
# Official TESS observing pages note that Sectors 97 and 98 each span four
# consecutive orbits instead of the standard two-sector cadence.
NONSTANDARD_ORBITS_PER_SECTOR = {
    97: 4,
    98: 4,
}
TICA_X_MIN = 44.0
TICA_X_MAX = 2093.0
TICA_Y_MIN = 0.0
TICA_Y_MAX = 2049.0
EDGE_WARN_PIXELS = 6.0


ProgressCallback = Callable[[str], None]


@dataclass(frozen=True)
class SectorCoverageConfig:
    min_sector: int = FIRST_200S_SECTOR
    max_sector: int | None = None
    chunk_size: int = 10_000


def _import_tess_point():
    try:
        from tess_stars2px import TESS_Spacecraft_Pointing_Data
    except ImportError as exc:
        raise ImportError(
            "TWIRL sector coverage mapping requires the 'tess-point' package. "
            "Install it with 'python -m pip install tess-point'."
        ) from exc
    return TESS_Spacecraft_Pointing_Data


def get_tess_point_version() -> str:
    try:
        return importlib.metadata.version("tess-point")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def get_available_sectors(min_sector: int = FIRST_200S_SECTOR, max_sector: int | None = None) -> np.ndarray:
    TESS_Spacecraft_Pointing_Data = _import_tess_point()
    scinfo = TESS_Spacecraft_Pointing_Data()
    sectors = np.asarray(scinfo.sectors, dtype=int)
    keep = sectors >= int(min_sector)
    if max_sector is not None:
        keep &= sectors <= int(max_sector)
    return sectors[keep]


def _compute_edge_warn(colpix: float, rowpix: float) -> bool:
    return (
        colpix <= (TICA_X_MIN + EDGE_WARN_PIXELS)
        or colpix >= (TICA_X_MAX - EDGE_WARN_PIXELS)
        or rowpix <= EDGE_WARN_PIXELS
        or rowpix >= (TICA_Y_MAX - EDGE_WARN_PIXELS)
    )


def build_sector_orbit_lookup(max_sector: int) -> Dict[int, np.ndarray]:
    if max_sector < FIRST_TESS_SECTOR:
        raise ValueError(f"max_sector must be >= {FIRST_TESS_SECTOR}")

    current_orbit = FIRST_TESS_ORBIT
    sector_orbits: Dict[int, np.ndarray] = {}
    for sector in range(FIRST_TESS_SECTOR, max_sector + 1):
        n_orbits = NONSTANDARD_ORBITS_PER_SECTOR.get(sector, DEFAULT_ORBITS_PER_SECTOR)
        sector_orbits[sector] = np.arange(current_orbit, current_orbit + n_orbits, dtype=np.int16)
        current_orbit += n_orbits
    return sector_orbits


def expand_sector_observation_table_to_orbits(observation_table: Table) -> Table:
    if "orbit" in observation_table.colnames:
        return observation_table.copy(copy_data=True)

    if "sector" not in observation_table.colnames:
        raise ValueError("observation_table must contain a 'sector' column")

    ordered_names = list(observation_table.colnames)
    insert_at = ordered_names.index("sector") + 1
    ordered_names.insert(insert_at, "orbit")

    sectors = np.asarray(observation_table["sector"], dtype=np.int16)
    if len(sectors) == 0:
        expanded = observation_table.copy(copy_data=True)
        expanded["orbit"] = np.array([], dtype=np.int16)
        return expanded[ordered_names]

    sector_orbit_lookup = build_sector_orbit_lookup(int(sectors.max()))
    orbit_lists = [sector_orbit_lookup[int(sector)] for sector in sectors]
    repeats = np.asarray([len(orbits) for orbits in orbit_lists], dtype=np.int32)

    expanded = Table()
    for name in observation_table.colnames:
        expanded[name] = np.repeat(np.asarray(observation_table[name]), repeats)

    expanded_orbits = np.concatenate(orbit_lists).astype(np.int16)
    expanded["orbit"] = expanded_orbits
    return expanded[ordered_names]


def _predict_chunk_matches_for_sector(
    fpg,
    row_indices: np.ndarray,
    ras: np.ndarray,
    decs: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Return all detector hits for one sector and one target chunk.

    This mirrors the internal tess-point `radec2pix` logic, but keeps track of the
    originating row index for each successful detector hit so the result can be
    merged back into the TWIRL master catalog.
    """

    row_indices = np.asarray(row_indices, dtype=np.int64)
    ras = np.asarray(ras, dtype=float)
    decs = np.asarray(decs, dtype=float)

    vec0s, vec1s, vec2s = fpg.sphereToCart(ras, decs)

    hit_rows: list[int] = []
    hit_cameras: list[int] = []
    hit_ccds: list[int] = []
    hit_cols: list[float] = []
    hit_rows_pix: list[float] = []
    hit_edge_warn: list[bool] = []

    for i, row_index in enumerate(row_indices):
        cur_vec = np.array([vec0s[i], vec1s[i], vec2s[i]], dtype=np.double)
        for camera_zero_based in range(fpg.NCAM):
            cam_vec = np.matmul(fpg.rmat4[camera_zero_based], cur_vec)
            lng, lat = fpg.cartToSphere(cam_vec)
            lng = np.degrees(lng)
            lat = np.degrees(lat)

            if not fpg.star_in_fov(lng, lat):
                continue

            xyfp = fpg.optics_fp(camera_zero_based, lng, lat)
            ccd_zero_based, ccd_pixels, _ = fpg.mm_to_pix(camera_zero_based, xyfp)
            colpix = float(ccd_pixels[0] + 45.0)
            rowpix = float(ccd_pixels[1] + 1.0)

            if not (TICA_X_MIN < colpix < TICA_X_MAX and TICA_Y_MIN < rowpix < TICA_Y_MAX):
                continue

            hit_rows.append(int(row_index))
            hit_cameras.append(int(camera_zero_based + 1))
            hit_ccds.append(int(ccd_zero_based + 1))
            hit_cols.append(colpix)
            hit_rows_pix.append(rowpix)
            hit_edge_warn.append(_compute_edge_warn(colpix, rowpix))

    return {
        "row_index": np.asarray(hit_rows, dtype=np.int32),
        "camera": np.asarray(hit_cameras, dtype=np.int16),
        "ccd": np.asarray(hit_ccds, dtype=np.int16),
        "colpix": np.asarray(hit_cols, dtype=np.float32),
        "rowpix": np.asarray(hit_rows_pix, dtype=np.float32),
        "edge_warn": np.asarray(hit_edge_warn, dtype=bool),
    }


def compute_tess_observation_table(
    table: Table,
    config: SectorCoverageConfig,
    progress: ProgressCallback | None = None,
) -> tuple[Table, Dict[str, Any]]:
    TESS_Spacecraft_Pointing_Data = _import_tess_point()
    scinfo = TESS_Spacecraft_Pointing_Data()

    all_sectors = np.asarray(scinfo.sectors, dtype=int)
    keep = all_sectors >= int(config.min_sector)
    if config.max_sector is not None:
        keep &= all_sectors <= int(config.max_sector)
    sectors = all_sectors[keep]
    sector_orbit_lookup = build_sector_orbit_lookup(int(sectors.max()))

    if len(sectors) == 0:
        raise ValueError(
            f"No tess-point sectors available in the requested range [{config.min_sector}, {config.max_sector}]"
        )

    source_ids = np.asarray(table["source_id"], dtype=np.int64)
    ras = np.asarray(table["ra"], dtype=float)
    decs = np.asarray(table["dec"], dtype=float)
    row_indices = np.arange(len(table), dtype=np.int64)

    chunk_size = max(1, int(config.chunk_size))
    n_chunks = math.ceil(len(table) / chunk_size)

    row_chunks: list[np.ndarray] = []
    sector_chunks: list[np.ndarray] = []
    orbit_chunks: list[np.ndarray] = []
    camera_chunks: list[np.ndarray] = []
    ccd_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    rowpix_chunks: list[np.ndarray] = []
    edge_chunks: list[np.ndarray] = []

    total_matches = 0
    for sector_idx, sector in enumerate(sectors, start=1):
        fpg_index = int(np.where(all_sectors == sector)[0][0])
        fpg = scinfo.fpgObjs[fpg_index]
        sector_orbits = sector_orbit_lookup[int(sector)]

        if progress is not None:
            progress(f"[coverage] sector {sector_idx}/{len(sectors)}: sector={sector} starting")

        for chunk_idx in range(n_chunks):
            start = chunk_idx * chunk_size
            stop = min(len(table), start + chunk_size)
            chunk_hits = _predict_chunk_matches_for_sector(
                fpg=fpg,
                row_indices=row_indices[start:stop],
                ras=ras[start:stop],
                decs=decs[start:stop],
            )
            if len(chunk_hits["row_index"]) == 0:
                continue

            hit_count = len(chunk_hits["row_index"])
            n_orbits = len(sector_orbits)
            match_count = hit_count * n_orbits
            total_matches += match_count

            row_chunks.append(np.repeat(chunk_hits["row_index"], n_orbits))
            sector_chunks.append(np.full(match_count, sector, dtype=np.int16))
            orbit_chunks.append(np.tile(sector_orbits, hit_count))
            camera_chunks.append(np.repeat(chunk_hits["camera"], n_orbits))
            ccd_chunks.append(np.repeat(chunk_hits["ccd"], n_orbits))
            col_chunks.append(np.repeat(chunk_hits["colpix"], n_orbits))
            rowpix_chunks.append(np.repeat(chunk_hits["rowpix"], n_orbits))
            edge_chunks.append(np.repeat(chunk_hits["edge_warn"], n_orbits))

            if progress is not None and ((chunk_idx + 1) % 10 == 0 or chunk_idx + 1 == n_chunks):
                progress(
                    "[coverage] "
                    f"sector {sector} chunk {chunk_idx + 1}/{n_chunks}: "
                    f"matches={total_matches}"
                )

        if progress is not None:
            progress(
                f"[coverage] sector {sector_idx}/{len(sectors)}: sector={sector} finished, total_matches={total_matches}"
            )

    if row_chunks:
        row_index = np.concatenate(row_chunks)
        sectors_out = np.concatenate(sector_chunks)
        orbits_out = np.concatenate(orbit_chunks)
        cameras = np.concatenate(camera_chunks)
        ccds = np.concatenate(ccd_chunks)
        colpix = np.concatenate(col_chunks)
        rowpix = np.concatenate(rowpix_chunks)
        edge_warn = np.concatenate(edge_chunks)
    else:
        row_index = np.array([], dtype=np.int32)
        sectors_out = np.array([], dtype=np.int16)
        orbits_out = np.array([], dtype=np.int16)
        cameras = np.array([], dtype=np.int16)
        ccds = np.array([], dtype=np.int16)
        colpix = np.array([], dtype=np.float32)
        rowpix = np.array([], dtype=np.float32)
        edge_warn = np.array([], dtype=bool)

    observation_table = Table()
    observation_table["catalog_row"] = row_index
    observation_table["source_id"] = source_ids[row_index] if len(row_index) else np.array([], dtype=np.int64)
    observation_table["sector"] = sectors_out
    observation_table["orbit"] = orbits_out
    observation_table["camera"] = cameras
    observation_table["ccd"] = ccds
    observation_table["colpix"] = colpix
    observation_table["rowpix"] = rowpix
    observation_table["edge_warn"] = edge_warn

    metadata = {
        "tess_point_version": get_tess_point_version(),
        "min_sector": int(config.min_sector),
        "max_sector": int(sectors.max()),
        "n_sectors": int(len(sectors)),
        "sectors": [int(sector) for sector in sectors],
        "max_orbit": int(orbits_out.max()) if len(orbits_out) else -1,
        "n_unique_orbits": int(len(np.unique(orbits_out))),
        "nonstandard_orbits_per_sector": {
            str(sector): int(n_orbits) for sector, n_orbits in sorted(NONSTANDARD_ORBITS_PER_SECTOR.items())
        },
        "n_catalog_rows": int(len(table)),
        "n_observation_rows": int(len(observation_table)),
    }
    return observation_table, metadata


def _build_row_json_payloads(
    n_rows: int,
    observation_table: Table,
) -> Dict[str, np.ndarray]:
    has_coverage = np.zeros(n_rows, dtype=bool)
    n_observations = np.zeros(n_rows, dtype=np.int32)
    n_sectors = np.zeros(n_rows, dtype=np.int16)
    n_orbits = np.zeros(n_rows, dtype=np.int16)
    sector_min = np.full(n_rows, -1, dtype=np.int16)
    sector_max = np.full(n_rows, -1, dtype=np.int16)
    orbit_min = np.full(n_rows, -1, dtype=np.int16)
    orbit_max = np.full(n_rows, -1, dtype=np.int16)
    json_payloads: list[str] = ["[]"] * n_rows

    if len(observation_table) == 0:
        return {
            "has_tess_200s_coverage": has_coverage,
            "n_tess_200s_observations": n_observations,
            "n_tess_200s_sectors": n_sectors,
            "n_tess_200s_orbits": n_orbits,
            "tess_200s_sector_min": sector_min,
            "tess_200s_sector_max": sector_max,
            "tess_200s_orbit_min": orbit_min,
            "tess_200s_orbit_max": orbit_max,
            "tess_observations_json": np.asarray(json_payloads, dtype=object),
        }

    row_index = np.asarray(observation_table["catalog_row"], dtype=np.int32)
    sector = np.asarray(observation_table["sector"], dtype=np.int16)
    orbit = np.asarray(observation_table["orbit"], dtype=np.int16)
    camera = np.asarray(observation_table["camera"], dtype=np.int16)
    ccd = np.asarray(observation_table["ccd"], dtype=np.int16)
    colpix = np.asarray(observation_table["colpix"], dtype=float)
    rowpix = np.asarray(observation_table["rowpix"], dtype=float)
    edge_warn = np.asarray(observation_table["edge_warn"], dtype=bool)

    order = np.lexsort((ccd, camera, orbit, sector, row_index))
    row_index = row_index[order]
    sector = sector[order]
    orbit = orbit[order]
    camera = camera[order]
    ccd = ccd[order]
    colpix = colpix[order]
    rowpix = rowpix[order]
    edge_warn = edge_warn[order]

    start = 0
    while start < len(row_index):
        row = int(row_index[start])
        stop = start + 1
        while stop < len(row_index) and row_index[stop] == row:
            stop += 1

        sector_slice = sector[start:stop]
        unique_sectors = np.unique(sector_slice)
        orbit_slice = orbit[start:stop]
        unique_orbits = np.unique(orbit_slice)
        payload = [
            {
                "sector": int(sector[k]),
                "orbit": int(orbit[k]),
                "camera": int(camera[k]),
                "ccd": int(ccd[k]),
                "colpix": round(float(colpix[k]), 3),
                "rowpix": round(float(rowpix[k]), 3),
                "edge_warn": bool(edge_warn[k]),
            }
            for k in range(start, stop)
        ]

        has_coverage[row] = True
        n_observations[row] = stop - start
        n_sectors[row] = len(unique_sectors)
        n_orbits[row] = len(unique_orbits)
        sector_min[row] = int(unique_sectors.min())
        sector_max[row] = int(unique_sectors.max())
        orbit_min[row] = int(unique_orbits.min())
        orbit_max[row] = int(unique_orbits.max())
        json_payloads[row] = json.dumps(payload, separators=(",", ":"))
        start = stop

    return {
        "has_tess_200s_coverage": has_coverage,
        "n_tess_200s_observations": n_observations,
        "n_tess_200s_sectors": n_sectors,
        "n_tess_200s_orbits": n_orbits,
        "tess_200s_sector_min": sector_min,
        "tess_200s_sector_max": sector_max,
        "tess_200s_orbit_min": orbit_min,
        "tess_200s_orbit_max": orbit_max,
        "tess_observations_json": np.asarray(json_payloads, dtype=object),
    }


def attach_tess_observation_columns(table: Table, observation_table: Table) -> Table:
    updated = table.copy(copy_data=True)
    payloads = _build_row_json_payloads(len(updated), observation_table)

    for name, data in payloads.items():
        if name == "tess_observations_json":
            max_len = max(len(value) for value in data) if len(data) else 2
            column = Column(
                data=np.asarray(data, dtype=f"U{max_len}"),
                name=name,
                description=(
                    "JSON array of orbit/sector/camera/CCD/pixel coverage records for the target."
                ),
            )
        elif name == "has_tess_200s_coverage":
            column = Column(
                data=data.astype(bool),
                name=name,
                description="True when at least one 200 s TESS detector hit has been attached.",
            )
        elif name == "n_tess_200s_observations":
            column = Column(
                data=data.astype(np.int32),
                name=name,
                description="Count of attached sector/camera/CCD observation records.",
            )
        elif name == "n_tess_200s_sectors":
            column = Column(
                data=data.astype(np.int16),
                name=name,
                description="Number of unique sectors represented in tess_observations_json.",
            )
        elif name == "n_tess_200s_orbits":
            column = Column(
                data=data.astype(np.int16),
                name=name,
                description="Number of unique orbits represented in tess_observations_json.",
            )
        elif name == "tess_200s_sector_min":
            column = Column(
                data=data.astype(np.int16),
                name=name,
                description="Minimum attached sector. -1 means none attached.",
            )
        elif name == "tess_200s_sector_max":
            column = Column(
                data=data.astype(np.int16),
                name=name,
                description="Maximum attached sector. -1 means none attached.",
            )
        elif name == "tess_200s_orbit_min":
            column = Column(
                data=data.astype(np.int16),
                name=name,
                description="Minimum attached orbit. -1 means none attached.",
            )
        else:
            column = Column(
                data=data.astype(np.int16),
                name=name,
                description="Maximum attached orbit. -1 means none attached.",
            )

        if name in updated.colnames:
            updated.replace_column(name, column)
        else:
            updated.add_column(column)

    return updated


def build_observation_export_table(base_table: Table, observation_table: Table) -> Table:
    export = Table()
    preferred_columns = [
        "source_id",
        "designation",
        "Pwd",
        "is_highconf_wd",
        "is_gaia_quality_ok",
        "tic_id",
        "tic_match_status",
        "ra",
        "dec",
        "pmra",
        "pmdec",
        "phot_g_mean_mag",
        "phot_bp_mean_mag",
        "phot_rp_mean_mag",
        "tmag",
    ]
    row_index = np.asarray(observation_table["catalog_row"], dtype=np.int32)

    export["catalog_row"] = row_index
    for name in preferred_columns:
        if name in base_table.colnames:
            export[name] = base_table[name][row_index]

    export["sector"] = observation_table["sector"]
    export["orbit"] = observation_table["orbit"]
    export["camera"] = observation_table["camera"]
    export["ccd"] = observation_table["ccd"]
    export["colpix"] = observation_table["colpix"]
    export["rowpix"] = observation_table["rowpix"]
    export["edge_warn"] = observation_table["edge_warn"]
    return export


def build_detector_summary_table(observation_export_table: Table) -> Table:
    grouped: dict[tuple[int, int, int, int], dict[str, int]] = {}

    for row in observation_export_table:
        key = (int(row["orbit"]), int(row["sector"]), int(row["camera"]), int(row["ccd"]))
        stats = grouped.setdefault(
            key,
            {
                "n_targets_all": 0,
                "n_targets_unique_tic": 0,
                "n_targets_highconf": 0,
                "n_targets_highconf_unique_tic": 0,
            },
        )
        stats["n_targets_all"] += 1

        if "is_highconf_wd" in observation_export_table.colnames and bool(row["is_highconf_wd"]):
            stats["n_targets_highconf"] += 1
        if "tic_id" in observation_export_table.colnames and int(row["tic_id"]) > 0:
            stats["n_targets_unique_tic"] += 1
            if "is_highconf_wd" in observation_export_table.colnames and bool(row["is_highconf_wd"]):
                stats["n_targets_highconf_unique_tic"] += 1

    summary = Table(names=[
        "orbit",
        "sector",
        "camera",
        "ccd",
        "n_targets_all",
        "n_targets_unique_tic",
        "n_targets_highconf",
        "n_targets_highconf_unique_tic",
    ], dtype=[np.int16, np.int16, np.int16, np.int16, np.int32, np.int32, np.int32, np.int32])

    for (orbit, sector, camera, ccd), stats in sorted(grouped.items()):
        summary.add_row(
            (
                orbit,
                sector,
                camera,
                ccd,
                stats["n_targets_all"],
                stats["n_targets_unique_tic"],
                stats["n_targets_highconf"],
                stats["n_targets_highconf_unique_tic"],
            )
        )

    return summary


def build_sector_summary_table(observation_export_table: Table) -> Table:
    grouped: dict[int, dict[str, Any]] = defaultdict(
        lambda: {
            "source_ids_all": set(),
            "source_ids_unique_tic": set(),
            "source_ids_highconf": set(),
            "source_ids_highconf_unique_tic": set(),
            "detector_keys": set(),
            "n_detector_target_rows": 0,
        }
    )

    has_tic = "tic_id" in observation_export_table.colnames
    has_highconf = "is_highconf_wd" in observation_export_table.colnames

    for row in observation_export_table:
        sector = int(row["sector"])
        source_id = int(row["source_id"])
        detector_key = (int(row["orbit"]), int(row["camera"]), int(row["ccd"]))
        grouped[sector]["source_ids_all"].add(source_id)
        grouped[sector]["detector_keys"].add(detector_key)
        grouped[sector]["n_detector_target_rows"] += 1

        is_unique_tic = has_tic and int(row["tic_id"]) > 0
        is_highconf = has_highconf and bool(row["is_highconf_wd"])

        if is_unique_tic:
            grouped[sector]["source_ids_unique_tic"].add(source_id)
        if is_highconf:
            grouped[sector]["source_ids_highconf"].add(source_id)
        if is_unique_tic and is_highconf:
            grouped[sector]["source_ids_highconf_unique_tic"].add(source_id)

    summary = Table(
        names=[
            "sector",
            "n_unique_targets_all",
            "n_unique_targets_unique_tic",
            "n_unique_targets_highconf",
            "n_unique_targets_highconf_unique_tic",
            "n_detector_target_rows",
            "n_orbit_detector_footprints",
        ],
        dtype=[np.int16, np.int32, np.int32, np.int32, np.int32, np.int32, np.int16],
    )

    for sector, stats in sorted(grouped.items()):
        summary.add_row(
            (
                sector,
                len(stats["source_ids_all"]),
                len(stats["source_ids_unique_tic"]),
                len(stats["source_ids_highconf"]),
                len(stats["source_ids_highconf_unique_tic"]),
                stats["n_detector_target_rows"],
                len(stats["detector_keys"]),
            )
        )

    return summary


def write_detector_target_tables(
    observation_export_table: Table,
    output_directory: Path,
    overwrite: bool = False,
) -> list[Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    written_paths: list[Path] = []

    grouped_indices: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for idx, row in enumerate(observation_export_table):
        key = (int(row["orbit"]), int(row["camera"]), int(row["ccd"]))
        grouped_indices[key].append(idx)

    for (orbit, camera, ccd), indices in sorted(grouped_indices.items()):
        subset = observation_export_table[np.asarray(indices, dtype=np.int32)]
        path = output_directory / f"orbit{orbit:04d}_cam{camera}_ccd{ccd}_twirl_targets.ecsv"
        subset.write(path, format="ascii.ecsv", overwrite=overwrite)
        written_paths.append(path)

    return written_paths
