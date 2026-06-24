# S56 1k Pre-Detrend Small-Pair 200k Recovery

Date: 2026-06-23

## Question

Does the small-aperture pair recover the pre-detrend BATMAN injections better
than the medium-ADP pilot when BLS uses a dense period grid?

## Inputs

- Injection set:
  `data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5`
- Apertures: `DET_FLUX_ADP_SML`, `DET_FLUX_SML`
- BLS periods: `200,000`
- Duration grid: `3, 4, 5, 6, 8, 10, 13, 16, 20, 30 min`
- Execution: restartable PDO chunks, `100` injected rows per chunk, merged after
  all `10` chunks completed.

## Result

| run | strict top-1 | exact top-N only | harmonic top-N | any exact/top-N/harmonic | unmatched |
| --- | ---: | ---: | ---: | ---: | ---: |
| medium ADP, `5k` | `32` | `11` | `25` | `68` | `932` |
| small pair, `5k` | `88` | `18` | `52` | `158` | `842` |
| small pair, `200k` | `137` | `10` | `30` | `177` | `823` |

The dense grid improves strict top-rank recovery, especially for high-SNR rows,
but it does not turn the faint low-SNR population into recoveries.

SNR-bin behavior for the `small_pair_200k` run:

| empirical SNR bin | top-1 | exact top-N | harmonic top-N | unmatched |
| --- | ---: | ---: | ---: | ---: |
| `<1` | `0` | `0` | `6` | `375` |
| `1-3` | `0` | `0` | `11` | `303` |
| `3-5` | `6` | `3` | `6` | `108` |
| `5-7` | `19` | `6` | `7` | `26` |
| `7-10` | `38` | `1` | `0` | `6` |
| `10-20` | `39` | `0` | `0` | `5` |
| `>20` | `35` | `0` | `0` | `0` |

## Interpretation

The small-aperture pair plus dense BLS is the best current pre-human-label
search path for this 1k injection pilot. The improvement is large enough to
rebuild the human-check LEO queue from this table.

The remaining miss population is not mainly a BLS-period-grid issue. Rows with
empirical SNR above `20` are all recovered at top rank, while rows below SNR `3`
are almost never recovered. The next `10k` sample should therefore keep the
small-aperture evidence path, but the faint low-SNR rows should be scored as
completeness-map failures rather than used as confident teacher positives.

## Outputs

- Merged recovery table:
  `small_pair_200k/injection_bls_recoveries.csv`
- Merged review queue before LEO:
  `small_pair_200k/review_queue.csv`
- Recovery summary:
  `small_pair_200k/recovery_mode_summary/summary.json`
- Merge summary:
  `merge_summary.json`
