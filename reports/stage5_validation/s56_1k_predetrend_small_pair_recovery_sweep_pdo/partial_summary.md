# S56 1k Pre-Detrend Small-Aperture Recovery Sweep

Date: 2026-06-23

This is the partial local summary while the dense `small_pair_200k` PDO sweep
continues in tmux session `twirl-s56-1k-smallpair-0623`.

## Inputs

- Injection HDF5:
  `data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5`
- Signal-survival table:
  `reports/stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_signal_survival_snr_diagnostics.csv`
- Output directory:
  `reports/stage5_validation/s56_1k_predetrend_small_pair_recovery_sweep_pdo/`

## Completed Sweeps

| sweep | apertures | period grid | top-1 | exact top-N | harmonic top-N | any exact/harmonic | unmatched |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `adp_sml_5k` | `DET_FLUX_ADP_SML` | `5,000` | `82` | `18` | `46` | `146` | `854` |
| `det_sml_5k` | `DET_FLUX_SML` | `5,000` | `81` | `12` | `38` | `131` | `869` |
| `small_pair_5k` | `DET_FLUX_ADP_SML + DET_FLUX_SML` | `5,000` | `88` | `18` | `52` | `158` | `842` |

For comparison, the old medium-ADP `5,000`-period table had `32/1000`
strict top-1 recoveries and `68/1000` exact-or-harmonic/top-N matches.

## Current State

The dense `small_pair_200k` sweep is still running on PDO. Last checked local
operator state: the job had entered the dense pass and reached at least
`250/1000` injected rows. Do not rebuild the final human queue until this
sweep writes its `recovery_mode_summary/summary.json`.

## Interim Interpretation

The full `1,000` sample confirms the small-aperture direction from the
`94`-row priority subset, but the gain is more modest than in the priority
debug set because most of the full sample remains low-SNR. The small-aperture
pair at `5k` improves strict top-1 recovery from `32` to `88` and any
exact/top-N/harmonic recovery from `68` to `158`.

The next decision depends on the dense `200k` pass:

- If `small_pair_200k` materially raises top-1 recovery for the full `1k`,
  rebuild the human-vetting queue from that sweep and render small-aperture
  LEO reports.
- If it does not, use `small_pair_5k` for near-term visual triage and shift
  search work toward systematics suppression / non-BLS matched filters rather
  than more period-grid density.
