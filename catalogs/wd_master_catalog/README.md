# WD Master Catalog

This directory documents the TWIRL WD master catalog. The built catalog is a
local data product and belongs under:

```text
data_local/catalogs/twirl_master_catalog/
```

The external seed is the Gentile Fusillo et al. (2021) Gaia EDR3 white-dwarf
catalog. Gaia DR3 `source_id` is the authoritative TWIRL target identifier;
TIC identifiers are operational metadata used by the current MIT extraction
stack.

## Current Catalog Design

- one row per seed-catalog white-dwarf target;
- preserve the seed columns rather than silently dropping or renaming them;
- add auditable TWIRL control-plane columns;
- keep `Pwd > 0.75` as the high-confidence reference flag, not an implicit
  final occurrence-rate denominator;
- store multi-visit TESS coverage as a JSON string per target while retaining
  the normalized one-row-per-hit observation table for pipeline work.

## Added Columns

- `is_highconf_wd`: `Pwd > 0.75` convenience flag
- `is_gaia_quality_ok`: provisional Gaia-quality convenience flag
- `tic_id`: TIC identifier, `-1` when unresolved
- `tic_match_status`: Gaia-to-TIC match status
- `tic_match_sep_arcsec`: Gaia-to-TIC angular separation
- `tmag`: TESS-magnitude estimate derived from Gaia photometry
- `has_tess_200s_coverage`: whether coverage metadata is attached
- `n_tess_200s_observations`: number of observation records
- `n_tess_200s_sectors`: number of unique sectors
- `n_tess_200s_orbits`: number of unique orbits
- `tess_200s_sector_min`, `tess_200s_sector_max`: covered sector range, or
  `-1` when absent
- `tess_200s_orbit_min`, `tess_200s_orbit_max`: covered orbit range, or `-1`
  when absent
- `tess_observations_json`: serialized per-target observation records

Example `tess_observations_json` value:

```json
[
  {"sector": 56, "orbit": 119, "camera": 1, "ccd": 3},
  {"sector": 56, "orbit": 120, "camera": 1, "ccd": 3}
]
```

## Build, Coverage, And Promotion

[build_wd_master_catalog.py](../../scripts/stage1_lightcurves/build_wd_master_catalog.py)
writes the FITS master catalog and a JSON provenance manifest.
[map_tess_sector_coverage.py](../../scripts/stage1_lightcurves/map_tess_sector_coverage.py)
populates the coverage columns and writes companion products, including:

- the coverage-enriched master catalog;
- `twirl_wd_tess_observations_v0.fits`, with one row per
  target/orbit/sector/camera/CCD hit;
- detector and sector summaries;
- optional per-orbit/camera/CCD ECSV target tables.

If existing products already contain sector/camera/CCD hits, use
[backfill_tess_orbits.py](../../scripts/stage1_lightcurves/backfill_tess_orbits.py)
to add orbit rows without repeating the sky-geometry calculation.

Once an intermediate catalog state is accepted, promote it with
[promote_master_catalog_version.py](../../scripts/stage1_lightcurves/promote_master_catalog_version.py).
Canonical releases use stable versioned names such as:

- `twirl_wd_master_catalog_v1.fits`
- `twirl_wd_master_catalog_v1_manifest.json`
- `twirl_wd_tess_observations_v1.fits`

Suffix-heavy names such as `*_ticmatched.fits` and `*_tesscoverage.fits` are
intermediate products, not long-term release names. Do not duplicate a large
accepted FITS product solely to make its filename cleaner.

## PDO Gaia-To-TIC Bridge

On PDO,
[export_gaia_dr3_tic_matches.py](../../scripts/stage1_lightcurves/export_gaia_dr3_tic_matches.py)
reads a TWIRL catalog and exports:

- `gaia_dr3_to_tic_matches.csv`, one row per Gaia DR3-to-TIC match;
- `gaia_dr3_to_tic_summary.csv`, one row per Gaia DR3 source;
- `gaia_dr3_to_tic_export_manifest.json`, with provenance and counts.

[merge_tic_summary_into_master_catalog.py](../../scripts/stage1_lightcurves/merge_tic_summary_into_master_catalog.py)
merges the summary conservatively:

- fill `tic_match_status` for Gaia rows represented in the summary;
- fill `tic_id` only for a unique match;
- retain `tic_id = -1` for ambiguous or unresolved matches.

The bridge must not redefine the scientific sample. Unresolved/no-TIC Gaia
targets remain an explicit support gap to characterize, and production
magnitude limits must not silently remove requested WD targets. The A2v1 reuse
path builds requested-TIC overlays from the normalized observation table; a
future Gaia-first emitter may still be required for scientifically important
targets with no usable TIC.

See [sector/orbit-map documentation](../sector_orbit_maps/README.md) and the
[local-data policy](../../doc/local_data.md) for companion products and paths.
