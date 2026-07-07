# Julien S56 10-System TWIRL Comparison

Generated: 2026-07-07

## Inputs and provenance

- Julien candidates: TIC/period/RF/SNR/n parsed from the ten supplied S56 PDF filenames.
- Primary TWIRL comparison now matches the active human-label sheet version: `s56_recovery50_teacher_queue/twirl_vet_sheets_fullphase_binmatch`, branch `current_adp`, apertures `DET_FLUX_ADP_SML` and `DET_FLUX_ADP`, `n_periods=20000`, `n_peaks=10`.
- The primary run reads real local S56 HLSP FITS from `data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare`; `injection_h5` and `lc_export_h5` are empty in the summary JSON.
- No injected-light-curve HDF5, `inj:*` review rows, or `injection_*` recovery products were used for the primary comparison.
- Supplementary ADP015 outputs remain in this directory, but they are not the primary human-label comparison in this report.

## Headline

- TWIRL has active-label real S56 two-aperture results for 9 of the 10 TICs.
- In those 9 available light curves, the active-label two-aperture anchor period is either Julien's period itself or a clean integer subharmonic of it.
- TIC `1400899528` is in the Gaia WD catalog, has a unique TIC bridge, and appears in the S56 detector target tables, but no S56 TGLC HDF5 or HLSP product exists for it.

## Per-target summary

| TIC | Julien P (d) | active-label top P (d) | relation | rel err | SDE | aperture flag | status / note |
|---:|---:|---:|---|---:|---:|---|---|
| 267574918 | 1.408300 | 1.407845 | P_J | 0.032% | 65.4 | False | Active human-label two-aperture top period matches Julien period directly. |
| 298666530 | 1.720300 | 0.573353 | P_J/3 | 0.014% | 19.0 | True | Active human-label two-aperture top period is an integer alias/subharmonic of Julien period (P_J/3). Active sheet flags aperture-period disagreement between small and primary apertures. |
| 406539395 | 0.327420 | 0.327479 | P_J | 0.018% | 210.9 | False | Active human-label two-aperture top period matches Julien period directly. |
| 1201294971 | 0.285950 | 0.286004 | P_J | 0.019% | 158.9 | False | Active human-label two-aperture top period matches Julien period directly. |
| 1400899528 | 0.989500 |  |  |  |  |  | In Gaia WD catalog and S56 target table, but no S56 TGLC HDF5/HLSP product exists for this TIC. |
| 1551395445 | 0.866450 | 0.288799 | P_J/3 | 0.006% | 98.8 | False | Active human-label two-aperture top period is an integer alias/subharmonic of Julien period (P_J/3). |
| 1883504789 | 0.787140 | 0.157452 | P_J/5 | 0.016% | 182.5 | False | Active human-label two-aperture top period is an integer alias/subharmonic of Julien period (P_J/5). |
| 1883654820 | 1.616300 | 0.202020 | P_J/8 | 0.009% | 170.7 | True | Active human-label two-aperture top period is an integer alias/subharmonic of Julien period (P_J/8). Active sheet flags aperture-period disagreement between small and primary apertures. |
| 1961800465 | 0.993630 | 0.141931 | P_J/7 | 0.012% | 166.4 | False | Active human-label two-aperture top period is an integer alias/subharmonic of Julien period (P_J/7). |
| 1976964633 | 1.476100 | 0.738262 | P_J/2 | 0.029% | 53.0 | False | Active human-label two-aperture top period is an integer alias/subharmonic of Julien period (P_J/2). |

## TIC 1400899528 audit

- Catalog: present in the Gaia WD/TIC bridge as Gaia source `1439208012920230656`, `tic_match_status=unique_tic`, TIC `1400899528`, catalog Tmag `20.0136`.
- S56 target emission: present in both S56 detector target tables for `orbit 119` and `orbit 120`, `cam4/ccd4`, `Pwd=0.910251`, `is_highconf_wd=True`, target-table `tmag=19.8303`, `edge_warn=False`.
- Product status: absent from `reports/stage1_lightcurves/s56_tic_ids/tic_ids_s0056_lightcurves_orbits119_120.csv`, absent from the local current-ADP HLSP shard, and absent on PDO from both orbit HDF5 paths and current/ADP015 HLSP shards.
- Root cause: not a detector-edge case. It is absent from the generated TGLC `TIC_cam4_ccd4.ecsv` catalogs because the TIC catalog magnitude is just fainter than the `--max-magnitude 20` catalog-construction threshold, even though the TWIRL/Gaia-derived target-table `tmag` is `19.8303`.
- Broader missing-HDF5 audit: `reports/stage1_lightcurves/s56_missing_h5_audit/summary.md` separates all S56 observation TICs into HDF5-present, TIC-magnitude-cut, detector-edge, and remaining boundary cases.

## Artifacts

- Active-label comparison CSV: `julien_twirl_s56_10_comparison_human_label_fullphase_binmatch.csv`
- Active-label two-aperture metrics: `twirl_vet_metrics_human_label_fullphase_binmatch.csv`
- Active-label vet sheets: `twirl_vet_sheets_human_label_fullphase_binmatch/`
- Active-label summary JSON: `twirl_two_aperture_vet_summary_human_label_fullphase_binmatch.json`
- Missing-product audit: `tic1400899528_missing_product_audit.json`
- Supplementary ADP015 comparison CSV: `julien_twirl_s56_10_comparison.csv`

## Notes

- SDE values are TWIRL BLS SDEs, not Julien's RF/SNR metric, so the amplitudes are not one-to-one comparable.
- Several systems show the expected ambiguity between a short BLS subharmonic and Julien's longer period; this is especially strong where each shorter-period event train is coherent in the small aperture.
- In the active-label version, TIC `298666530` and TIC `1883654820` have aperture-period disagreement flags; inspect those sheets before treating the primary-aperture period as equivalent to the small-aperture anchor.
