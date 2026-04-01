# Sector And Detector Maps

This directory holds repo-tracked notes for the TWIRL sector and detector mapping products.

The heavy outputs themselves are local-only and should be written under:

```text
data_local/catalogs/twirl_master_catalog/
```

## Purpose

These products answer a different question from the WD master catalog:

- the master catalog is one row per Gaia DR3 target
- the sector/detector maps are one row per target-sector-camera-CCD hit

They are meant for:

- counting how many sectors each target falls on
- planning which sector/camera/CCD footprints need processing
- exporting TWIRL-side target subsets for later filtering

They are not replacements for the full-field MKI TGLC cached TIC and Gaia catalogs written by
`tglc catalogs`.

## Expected Local Outputs

The coverage mapper should write:

- a master catalog FITS file with attached `tess_observations_json`
- a one-row-per-hit observation table
- a detector summary table with one row per `sector/camera/ccd`
- a sector summary table
- optional per-detector target tables in ECSV format

The sector list should include all requested `Sector >= 56` values, including future scheduled
sectors when they are present in the pointing table. Observed and future sectors should remain in
the same products so later extraction code can apply an up-to-date observed-sector cutoff.

The plotting helper `scripts/step1_lcs/plot_tess_observation_sky_coverage.py` consumes the
coverage-enriched catalog and the observation table to visualize sky coverage split into currently
observed and future-planned sectors.
