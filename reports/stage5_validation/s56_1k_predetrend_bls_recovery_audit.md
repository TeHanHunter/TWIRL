# S56 1k Pre-Detrend BLS Recovery Audit

Date: 2026-06-23

Inputs:

- Queue: `reports/stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_topn5000_pdo/injection_bls_recoveries.csv`
- Signal-survival table: `reports/stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_signal_survival_snr_diagnostics.csv`
- Failure diagnostics: `reports/stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_topn5000_pdo/bls_failure_diagnostics/`
- Aperture-survival report: `reports/stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/aperture_survival_summary/`

## Bottom Line

The low recovery is not random across magnitude. It is mostly the expected loss of low empirical-SNR signals after detrending, with two secondary failures:

1. the current medium `DET_FLUX_ADP` aperture is not the best product for many WDs;
2. among the few SNR-qualified misses, the BLS top peak is usually an alias/harmonic or a stronger systematic period, not the injected period.

Using the `n_periods=5000` top-N accounting, the ADP-only run has:

- strict top-1 recovery: `32 / 1000`
- exact top-N non-top-1 recovery: `11 / 1000`
- harmonic top-N match: `25 / 1000`
- unmatched: `932 / 1000`

With an empirical post-detrend signal-SNR gate at `7`, the unmatched rows split into:

- `904 / 1000` low-SNR unmatched rows
- `28 / 1000` SNR-qualified BLS-ranking-loss rows

So the dominant problem is not "BLS randomly fails." The dominant problem is that most injected signals are not recoverable in the current detrended product.

## Magnitude Dependence

Any exact-or-harmonic period match by Tmag bin:

| Tmag bin | matched / total |
| --- | ---: |
| `<16` | `10 / 11` |
| `16-17` | `10 / 18` |
| `17-18` | `17 / 68` |
| `18-19` | `16 / 219` |
| `19-20` | `15 / 659` |
| `>20` | `0 / 25` |

The steep magnitude trend is mostly a signal-SNR trend. Rows with empirical SNR `>=7` by Tmag bin:

| Tmag bin | SNR-qualified / total | period matched among SNR-qualified |
| --- | ---: | ---: |
| `<16` | `9 / 11` | `9 / 9` |
| `16-17` | `12 / 18` | `9 / 12` |
| `17-18` | `22 / 68` | `12 / 22` |
| `18-19` | `20 / 219` | `10 / 20` |
| `19-20` | `8 / 659` | `3 / 8` |
| `>20` | `0 / 25` | `0 / 0` |

This means we are mostly losing faint targets because they do not survive as high-SNR events in `DET_FLUX_ADP`, not because BLS has a magnitude-dependent bug. The high-SNR sample at `Tmag >= 19` is tiny (`8` rows), so that bin is not enough to diagnose subtle algorithmic magnitude effects.

## Empirical-SNR Dependence

Any exact-or-harmonic period match by empirical post-detrend multi-cadence SNR:

| empirical SNR bin | matched / total |
| --- | ---: |
| `<1` | `6 / 527` |
| `1-3` | `8 / 277` |
| `3-5` | `8 / 96` |
| `5-7` | `3 / 29` |
| `7-10` | `11 / 24` |
| `10-20` | `15 / 30` |
| `>20` | `17 / 17` |

This is the cleanest diagnostic. Once the injected event is very high SNR (`>20`), BLS finds an exact-or-harmonic period every time in this sample. The incomplete recovery at `7-20` SNR is real and should be the target set for BLS/ranking improvements.

## Period Dependence

Any exact-or-harmonic period match by injected period:

| period bin [d] | matched / total |
| --- | ---: |
| `<0.25` | `22 / 153` |
| `0.25-0.5` | `16 / 154` |
| `0.5-1` | `13 / 151` |
| `1-2` | `11 / 142` |
| `2-5` | `5 / 198` |
| `5-10` | `1 / 151` |
| `>10` | `0 / 51` |

This period trend is expected for a 27-day TESS sector plus minute-scale WD transits: longer periods have fewer events, more gap sensitivity, and weaker BLS period coherence. It also argues that the long-period completeness map cannot rely on periodic BLS alone; it needs a dip/single-event branch.

## Aperture/Dataset Issue

The current run used medium `DET_FLUX_ADP`. Signal survival across products shows:

- current `DET_FLUX_ADP`: `71 / 1000` rows with empirical SNR `>=7`, median SNR `0.89`, median depth retention `0.28`
- best single product `DET_FLUX_ADP_SML`: `124 / 1000` rows with empirical SNR `>=7`, median SNR `1.58`, median depth retention `0.50`
- best of all detrended apertures: `137 / 1000` rows with empirical SNR `>=7`

The small apertures can probably recover some of the current losses, but even the best available detrended aperture leaves most injections below the empirical-SNR threshold.

## Interpretation

For the current 1k pre-detrend sample, the hierarchy of failure is:

1. **detrending/product SNR loss dominates**: most faint targets have post-detrend signal SNR below any realistic BLS threshold;
2. **aperture choice is a real but partial fix**: small ADP/canonical apertures improve survival for dozens of cases, not hundreds;
3. **BLS ranking still needs work**: the `28` SNR-qualified misses are the right compact debug set;
4. **long-period signals require a separate branch**: periodic BLS is weak for sparse minute-scale events in a single sector.

The next action should be the priority aperture/BLS sweep on the `94` debug injections. If small apertures convert many of the `66` aperture-recoverable rows, rebuild the 1k human queue on that setting. If they do not, tune the detrending/search statistic before scaling the 10k injection run.
