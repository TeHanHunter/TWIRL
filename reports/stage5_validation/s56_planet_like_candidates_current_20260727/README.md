# Sector 56 current Planet-like candidate bundle

This folder contains the newest uniform vetting sheet for every final,
human-accepted `planet_like` morphology candidate in Sector 56. The source
decision snapshot is the accepted `2026-07-24` S56--S62 morphology corpus;
this subset contains **30 unique S56 TICs**:

- `12` from the
  original human-adjudicated S56 set.
- `18`
  from the later model-enriched compact revisit, each subsequently accepted by
  the human morphology pass.
- `13` with aligned periods in both
  current ADP apertures and
  `17` where the current primary
  aperture prefers a different period.

## Priority groups

- **Benchmark (1):** WD 1856 b, the confirmed
  planet and mandatory engineering benchmark.
- **A (12):** current `DET_FLUX_ADP_SML` and
  `DET_FLUX_ADP` period solutions are aligned.
- **B (15):** Planet-like morphology in the small
  search aperture, but the current primary aperture prefers a different
  period.
- **C (2):** the same aperture disagreement plus a
  non-blocking Tier-1 target-QA review warning.

The ordering is a transparent follow-up triage, not an astrophysical
classification or confirmation. Within A, rows are ordered by the smaller of
the two aperture SDE values. Within B/C, rows are ordered by small-aperture SDE.
No model probability enters the ordering. Odd/even and target-QA flags remain
explicit columns in the summary.

## Sheet contract

- Sheet-set version: `s56_s62_a2v1_current_adp_v1`
- Renderer: `S56-ADP-HV2`
- Branch: `current_adp`
- Apertures: `DET_FLUX_ADP_SML, DET_FLUX_ADP`
- Trial periods: `20,000`
- Retained peaks: `10`
- Row ephemeris preserved: `True`

The `morphology_fold_factor` and `morphology_view_period_d` fields describe the
evidence view accepted in morphology review. They do not by themselves claim
the physical orbital period.

## Files

- `s56_planet_like_candidates.xlsx`: sortable summary workbook.
- `s56_planet_like_candidates.csv`: machine-readable candidate table.
- `s56_planet_like_candidates_summary.png` and `.pdf`: one-page visual index.
- `01_TIC...png` through `30_TIC...png`: hash-verified copies of the newest
  current-ADP vetting sheets.
- `manifest.json`: source and output checksums plus selection counts.
