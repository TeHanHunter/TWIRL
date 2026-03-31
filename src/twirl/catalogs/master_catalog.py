from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
from astropy.table import Column, Table


DEFAULT_BUILD_VERSION = "v0"


def estimate_tess_magnitude(
    g_mag: np.ndarray,
    bp_mag: np.ndarray,
    rp_mag: np.ndarray,
) -> np.ndarray:
    """Estimate TESS magnitude from Gaia photometry.

    Uses the Gaia-to-TESS polynomial adopted in the existing Step 1 plotting
    workflow. Rows with missing BP/RP fall back to a simple G-band offset.
    """

    g_mag = np.asarray(g_mag, dtype=float)
    bp_mag = np.asarray(bp_mag, dtype=float)
    rp_mag = np.asarray(rp_mag, dtype=float)

    color = bp_mag - rp_mag
    tmag = (
        g_mag
        - 0.00522555 * color**3
        + 0.0891337 * color**2
        - 0.633923 * color
        + 0.0324473
    )

    bad_tmag = ~np.isfinite(tmag)
    tmag[bad_tmag] = g_mag[bad_tmag] - 0.430
    return tmag


def preliminary_gaia_quality_mask(table: Table) -> np.ndarray:
    """Return a conservative Gaia-quality convenience flag.

    This is intentionally not the TWIRL parent-sample definition. It is only a
    lightweight builder-time QA flag so downstream code can distinguish WD
    confidence from basic astrometric/photometric cleanliness.
    """

    ruwe = np.asarray(table["ruwe"], dtype=float)
    duplicated_source = np.asarray(table["duplicated_source"], dtype=bool)
    parallax = np.asarray(table["parallax"], dtype=float)
    g_mag = np.asarray(table["phot_g_mean_mag"], dtype=float)
    bp_mag = np.asarray(table["phot_bp_mean_mag"], dtype=float)
    rp_mag = np.asarray(table["phot_rp_mean_mag"], dtype=float)

    return (
        np.isfinite(ruwe)
        & (ruwe <= 1.4)
        & ~duplicated_source
        & np.isfinite(parallax)
        & (parallax > 0.0)
        & np.isfinite(g_mag)
        & np.isfinite(bp_mag)
        & np.isfinite(rp_mag)
    )


def compute_file_sha256(path: Path, chunk_size: int = 1 << 20) -> str:
    """Compute the SHA256 hash for a local file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def build_manifest(
    input_path: Path,
    output_path: Path,
    row_count: int,
    build_version: str,
    highconf_pwd_threshold: float,
    include_input_sha256: bool,
) -> Dict[str, Any]:
    input_stat = input_path.stat()
    built_at = datetime.now(timezone.utc).isoformat()

    manifest = {
        "build_version": build_version,
        "built_at_utc": built_at,
        "seed_catalog_path": str(input_path.resolve()),
        "seed_catalog_size_bytes": input_stat.st_size,
        "seed_catalog_mtime_utc": datetime.fromtimestamp(
            input_stat.st_mtime,
            tz=timezone.utc,
        ).isoformat(),
        "seed_catalog_sha256": None,
        "row_count": row_count,
        "highconf_pwd_threshold": highconf_pwd_threshold,
        "output_path": str(output_path.resolve()),
        "placeholder_columns": [
            "tic_id",
            "tic_match_status",
            "tic_match_sep_arcsec",
            "tmag",
            "has_tess_200s_coverage",
            "n_tess_200s_observations",
            "n_tess_200s_sectors",
            "tess_200s_sector_min",
            "tess_200s_sector_max",
            "tess_observations_json",
        ],
        "notes": [
            "Gaia DR3 remains the authoritative target identifier.",
            "is_highconf_wd tracks Pwd > threshold and is separate from is_gaia_quality_ok.",
            "tess_observations_json is a JSON array placeholder for per-target TESS coverage metadata.",
        ],
    }

    if include_input_sha256:
        manifest["seed_catalog_sha256"] = compute_file_sha256(input_path)

    return manifest


def _replace_or_add_column(table: Table, column: Column) -> None:
    if column.name in table.colnames:
        table.replace_column(column.name, column)
    else:
        table.add_column(column)


def build_master_catalog(
    input_path: Path,
    build_version: str = DEFAULT_BUILD_VERSION,
    highconf_pwd_threshold: float = 0.75,
    include_input_sha256: bool = False,
    output_path: Optional[Path] = None,
) -> Tuple[Table, Dict[str, Any]]:
    """Build the first-pass TWIRL master catalog from the local WD seed file."""

    table = Table.read(input_path, hdu="Joined", unit_parse_strict="silent")
    table.convert_bytestring_to_unicode()

    row_count = len(table)
    pwd = np.asarray(table["Pwd"], dtype=float)
    highconf_mask = np.isfinite(pwd) & (pwd > highconf_pwd_threshold)
    gaia_quality_mask = preliminary_gaia_quality_mask(table)
    tmag = estimate_tess_magnitude(
        g_mag=np.asarray(table["phot_g_mean_mag"], dtype=float),
        bp_mag=np.asarray(table["phot_bp_mean_mag"], dtype=float),
        rp_mag=np.asarray(table["phot_rp_mean_mag"], dtype=float),
    )

    _replace_or_add_column(
        table,
        Column(
            name="is_highconf_wd",
            data=highconf_mask.astype(bool),
            description=f"True when Pwd > {highconf_pwd_threshold:.2f}.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="is_gaia_quality_ok",
            data=gaia_quality_mask.astype(bool),
            description=(
                "Preliminary Gaia quality convenience flag based on RUWE, "
                "duplicate-source status, parallax, and finite Gaia photometry."
            ),
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="tic_id",
            data=np.full(row_count, -1, dtype=np.int64),
            description="Placeholder TIC identifier. -1 means unresolved/unmatched.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="tic_match_status",
            data=np.full(row_count, "unresolved", dtype="U10"),
            description="Placeholder Gaia-to-TIC match status.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="tic_match_sep_arcsec",
            data=np.full(row_count, np.nan, dtype=float),
            description="Angular separation for the adopted Gaia-to-TIC match, in arcsec.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="tmag",
            data=tmag,
            description="Estimated TESS magnitude derived from Gaia photometry.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="has_tess_200s_coverage",
            data=np.zeros(row_count, dtype=bool),
            description="True when at least one 200 s TESS observation has been attached.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="n_tess_200s_observations",
            data=np.zeros(row_count, dtype=np.int32),
            description="Count of attached 200 s TESS sector/orbit/camera/CCD observation records.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="n_tess_200s_sectors",
            data=np.zeros(row_count, dtype=np.int16),
            description="Number of unique sectors represented in tess_observations_json.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="tess_200s_sector_min",
            data=np.full(row_count, -1, dtype=np.int16),
            description="Minimum attached 200 s TESS sector. -1 means none attached.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="tess_200s_sector_max",
            data=np.full(row_count, -1, dtype=np.int16),
            description="Maximum attached 200 s TESS sector. -1 means none attached.",
        ),
    )
    _replace_or_add_column(
        table,
        Column(
            name="tess_observations_json",
            data=np.full(row_count, "[]", dtype="U2"),
            description=(
                "JSON array of per-target TESS observation records, e.g. "
                "[{\"sector\": 56, \"orbit\": 185, \"camera\": 1, \"ccd\": 3}]."
            ),
        ),
    )

    manifest_output_path = output_path if output_path is not None else input_path
    manifest = build_manifest(
        input_path=input_path,
        output_path=manifest_output_path,
        row_count=row_count,
        build_version=build_version,
        highconf_pwd_threshold=highconf_pwd_threshold,
        include_input_sha256=include_input_sha256,
    )

    table.meta["BUILDVER"] = build_version
    table.meta["BLDTIME"] = manifest["built_at_utc"]
    table.meta["SEEDPATH"] = manifest["seed_catalog_path"]
    table.meta["SEEDSIZE"] = manifest["seed_catalog_size_bytes"]
    table.meta["SEEDMTIM"] = manifest["seed_catalog_mtime_utc"]
    table.meta["PWDCUT"] = highconf_pwd_threshold

    if manifest["seed_catalog_sha256"] is not None:
        table.meta["SEEDSHA"] = manifest["seed_catalog_sha256"]

    return table, manifest


def write_master_catalog(
    table: Table,
    manifest: Dict[str, Any],
    output_path: Path,
    manifest_path: Path,
    overwrite: bool = False,
) -> None:
    """Write the built catalog and sidecar manifest to disk."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    table.write(output_path, overwrite=overwrite)
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
