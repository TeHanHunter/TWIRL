# TGLC production status S56-S93

Snapshot: `pdogpu1 2026-06-02T20:55:55-04:00`.

- Prepared cutouts complete: `33/38` sectors; partial/active prep: `3`.
- ePSF fitted: `8/38` sectors.
- TGLC HDF5 light curves complete: `8/38` sectors.
- Final FITS products exist: `7/38` sectors; current TWIRL-FS v2 sectors: `1`.
- `qc_pause.flag` is present, so normal GPU finalize is intentionally paused.

Active / problematic tail:
- S91: Active prep lane-B from 2026-05-30; orbit 189 has 2870/3136 source pickles and orbit 190 has not started.
- S92: Active prep lane-A from 2026-06-01; orbit 191 is still at 0/3136 source pickles and orbit 192 has not started.
- S66: Stale lease s66.lease from 2026-05-07; source count is 4045/6272 and no cutout-done marker exists.
- S67: Stale lease s67.lease from 2026-05-12; source count is 0/6272 and no cutout-done marker exists.

Artifacts:
- PNG: `tglc_s56_s93_production_status.png`
- PDF: `tglc_s56_s93_production_status.pdf`
