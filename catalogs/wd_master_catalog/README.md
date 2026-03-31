# WD Master Catalog

This directory holds repo-tracked documentation for the TWIRL WD master catalog.
The full built table itself is local-only and should be written under:

```text
data_local/catalogs/twirl_master_catalog/
```

## Current v0 Design

- one row per Gaia DR3 WD seed-catalog target
- preserve all columns from the local Gentile Fusillo et al. (2021) Gaia EDR3 seed file
- add TWIRL control-plane columns rather than dropping or renaming seed columns
- keep Gaia DR3 as the authoritative target identifier
- store multi-visit TESS coverage in a single JSON-string cell per target

## Added Columns In The First Builder

- `is_highconf_wd`: `Pwd > 0.75` convenience flag
- `is_gaia_quality_ok`: provisional Gaia-quality convenience flag
- `tic_id`: placeholder TIC identifier, `-1` when unresolved
- `tic_match_status`: placeholder Gaia-to-TIC match status
- `tic_match_sep_arcsec`: placeholder Gaia-to-TIC angular separation
- `tmag`: estimated TESS magnitude derived from Gaia photometry
- `has_tess_200s_coverage`: summary boolean for attached TESS coverage metadata
- `n_tess_200s_observations`: number of TESS observation records attached
- `n_tess_200s_sectors`: number of unique sectors represented
- `tess_200s_sector_min`: minimum covered sector, `-1` when none
- `tess_200s_sector_max`: maximum covered sector, `-1` when none
- `tess_observations_json`: JSON array placeholder for per-target TESS observation records

Example JSON payload for `tess_observations_json`:

```json
[
  {"sector": 56, "orbit": 185, "camera": 1, "ccd": 3},
  {"sector": 56, "orbit": 186, "camera": 1, "ccd": 3}
]
```

## Builder Outputs

The builder writes two local files:

- a FITS master catalog
- a JSON sidecar manifest with provenance, build version, and seed-file metadata

## PDO Match Export

When run on PDO, the read-only script
`scripts/step1_lcs/export_gaia_dr3_tic_matches.py`
can export Gaia DR3 to TIC matches for the `source_id` values in a local TWIRL catalog.

It writes:

- `gaia_dr3_to_tic_matches.csv`: one row per Gaia DR3 to TIC match
- `gaia_dr3_to_tic_summary.csv`: one row per Gaia DR3 source ID with match counts and status
- `gaia_dr3_to_tic_export_manifest.json`: export metadata and summary counts
