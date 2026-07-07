# S56 Missing HDF5 Audit

Generated: 2026-07-07

## Inputs

- S56 target observations: `data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits`
- TIC bridge magnitudes: `data_local/catalogs/twirl_master_catalog/gaia_dr3_to_tic_summary.csv`
- HDF5 product TIC index: `reports/stage1_lightcurves/s56_tic_ids/tic_ids_s0056_lightcurves_orbits119_120.csv`
- Product-index definition: union of TIC IDs with TGLC HDF5 light curves in S56 orbits `119` and `120`.

## Headline

- S56 observation rows with TIC IDs: `63,238` rows, `31,590` unique TICs.
- TGLC HDF5 products: `19,072` unique TICs.
- Missing relative to all TWIRL S56 observation TICs: `12,518` unique TICs.
- Most missing TICs are explained by the TIC-side catalog magnitude cut: `12,430` missing TICs have `unique_tic_tmag > 20`.
- Within the nominal TIC `Tmag <= 20` catalog cut, `88` TICs are missing; `85` of those have `edge_warn=True`, and `3` do not.

## Status Classes

| class | count |
|---|---:|
| `has_h5` | 19,072 |
| `tic_tmag_gt_20_catalog_cut` | 12,430 |
| `within_ticmag20_edge_warn` | 85 |
| `within_ticmag20_no_edge_tmag_boundary` | 3 |

## TIC 1400899528

- Present in the Gaia WD/TIC bridge and S56 observation table: Gaia source `1439208012920230656`, `Pwd=0.910251`, `any_highconf_wd=True`.
- Detector placement: `cam4/ccd4`, orbits `119;120`, `colpix=720.312`, `rowpix=1764.882`, `edge_warn=False`.
- Magnitudes: TWIRL/Gaia-derived target-table `tmag=19.830296`, TIC catalog `unique_tic_tmag=20.013600`.
- Audit classification: `tic_tmag_gt_20_catalog_cut`.

Interpretation: TIC `1400899528` is not absent because it is missing from the Gaia WD catalog, and it is not flagged as a detector-edge target in the TWIRL detector table. The likely drop happens at TGLC catalog construction because the TIC catalog magnitude is just fainter than the `--max-magnitude 20` threshold.

## Missing Within TIC Tmag <= 20 And No Edge Warning

These are the small remaining non-edge cases inside the nominal TIC magnitude cut and should be treated as the next extraction-debug set.

| TIC | Gaia source(s) | max Pwd | TWIRL tmag | TIC Tmag | detector | colpix | rowpix | high-conf WD |
|---:|---|---:|---:|---:|---|---:|---:|---|
| 2004030067 | 1945624724870621568 | 0.214130 | 19.970882 | 20.000000 | cam2/ccd3 | 882.627 | 1825.319 | False |
| 2004659782 | 1950829228802227072 | 0.094979 | 19.930593 | 20.000000 | cam2/ccd2 | 861.880 | 1623.211 | False |
| 2053663055 | 2829622895256028928 | 0.986306 | 20.051444 | 20.000000 | cam1/ccd2 | 360.584 | 1535.355 | True |

## Artifacts

- Full observation-to-HDF5 audit: `s56_observation_to_h5_audit.csv`
- Missing-HDF5 rows only: `s56_missing_h5_all.csv`
- Missing-HDF5 rows with TIC `Tmag <= 20`: `s56_missing_h5_within_ticmag20.csv`
- Julien's ten supplied TICs: `s56_julien_target_h5_audit.csv`
- Machine-readable summary: `summary.json`

## Notes

- This audit uses TIC catalog magnitude for the TGLC catalog-cut classification, not the Gaia-derived TWIRL target-table `tmag`.
- `edge_warn` is taken directly from the S56 TWIRL observation table.
- The three non-edge `Tmag <= 20` rows all sit at the TIC magnitude boundary and need targeted catalog/log checks before being counted as extraction failures.
