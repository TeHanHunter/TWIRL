# Sector And Orbit Maps

This directory documents the TWIRL sector, orbit, and detector mapping
products. The products themselves are local-only and belong under:

```text
data_local/catalogs/twirl_master_catalog/
```

## Purpose

The products answer a different question from the WD master catalog:

- the master catalog has one row per Gaia DR3 target;
- the observation table has one row per target/orbit/sector/camera/CCD hit.

They support sector-count summaries, orbit/camera/CCD production planning, and
TWIRL target-table exports. They do not replace the full-field MKI TGLC TIC and
Gaia caches written by `tglc catalogs`.

## Expected Local Outputs

The coverage mapper writes:

- a coverage-enriched master-catalog FITS file containing
  `tess_observations_json`;
- a one-row-per-hit observation table;
- a detector summary with one row per orbit/sector/camera/CCD;
- a sector summary;
- optional per-orbit/camera/CCD ECSV target tables.

Run [map_tess_sector_coverage.py](../../scripts/stage1_lightcurves/map_tess_sector_coverage.py)
to build the coverage products. If sector-level products already exist and
only orbit expansion is missing, use
[backfill_tess_orbits.py](../../scripts/stage1_lightcurves/backfill_tess_orbits.py)
instead of rerunning the full `tess-point` geometry pass.

The sector list should contain every requested `Sector >= 56`, including
future scheduled sectors present in the pointing table. Keep observed and
future sectors in the same product so extraction code can apply an explicit,
up-to-date observed-sector cutoff.

Orbit expansion follows the official TESS observing pattern: sectors normally
span two orbits, with the currently represented four-orbit exceptions at
Sectors 97 and 98. Each sector-level detector hit expands to one row per orbit
so Stage 1 can schedule orbit/camera/CCD work directly.

Use
[plot_tess_observation_sky_coverage.py](../../scripts/stage1_lightcurves/plot_tess_observation_sky_coverage.py)
to visualize coverage from the enriched catalog and observation table,
separating observed from future-planned sectors.

The observation table is also the A2v1 target-emission authority: Stage 1
builds lightweight `source_tic` overlays from its detector coordinates so
requested WD TICs can be emitted without regenerating broad all-TIC source
catalogs. Gaia DR3 `source_id` remains the scientific target identifier; TIC
is operational metadata.

See [the WD catalog README](../wd_master_catalog/README.md) for catalog schema
and promotion conventions and [local-data policy](../../doc/local_data.md) for
the expected filesystem layout.
