# S56 Pre-Detrend Detrending Method Audit

Date: 2026-06-23

## Question

Is detrending still suppressing injected WD-transit signals, and is there a better
method that preserves most injected BATMAN signals while removing most long
timescale trends?

## Inputs

- Injection set:
  `data_local/stage3_injections/s56_twirlfs_v2_injection_training/pdo_1k_predetrend_batman_depthgrid_adp_compare/injected_lightcurves.h5`
- Priority subset: `94` rows selected from SNR-qualified BLS-ranking failures
  plus rows newly SNR-qualified in alternate apertures.
- Full confirmation set: all `1,000` pre-detrend BATMAN injections.

The method sweep reads the stored original and injected raw TGLC `RawFlux`
arrays, reruns detrending method variants, then measures the injected-minus-
original truth-window signal. The main metrics are empirical multi-cadence SNR,
truth-window depth retention, residual binned trend 90-10% range, and trend
reduction relative to raw normalized flux.

## Methods Tested

- `constant`: no trend model beyond a per-light-curve normalization.
- `poly2_gap05`, `poly3_gap05`: low-order polynomial fits split across gaps.
- `median06_gap05`, `median10_gap05`, `median20_gap05`: rolling-median
  filter cotrends with 0.6, 1.0, and 2.0 d windows.
- `savgol06_p2_gap05`, `savgol10_p2_gap05`, `savgol20_p2_gap05`:
  Savitzky-Golay order-2 cotrends with 0.6, 1.0, and 2.0 d windows.
- `median03_gap05`, `savgol03_p2_gap05`: shorter 0.3 d local-filter
  baselines, added as a direct test of whether a more local but still
  production-cheap trend model improves the truth-window SNR.
- `pctl60_06_gap05`, `pctl75_*_gap05`, `pctl90_06_gap05`: rolling
  percentile/envelope filters. These are production-cheap non-oracle tests of
  whether using a high-percentile local baseline protects downward transit
  signals better than a median baseline.
- `current_adp03q`: current adaptive comparison product style, using a
  `0.3 d` quantile-knot spline with `0.2 d` gap split.
- `spline_q05_gap02`, `spline_q08_gap02`, `spline_q12_gap05`: quantile-knot
  splines with coarser knot spacing and gap split choices.
- `current_uniform08_gap05`: canonical-like uniform `0.8 d` spline with
  `0.5 d` gap split.
- `spline_uniform12_gap05`, `spline_uniform20_gap05`: looser uniform splines.
- `oracle_adp03q`, `oracle_q05_gap02`: non-production upper-bound tests that
  mask the known injected in-transit cadences from the cotrend fit.

The initial constant/polynomial/spline set was tested on raw `Small`,
`Primary`, and `Large` apertures. The additional filter/oracle methods were
then tested on the full `Small` aperture set, since the first pass showed that
`Small` dominates injected signal SNR.

## Results

### Full 1k set

| product | rows above SNR 7 | median depth retention | median empirical SNR | residual trend ptp | trend reduction |
| --- | ---: | ---: | ---: | ---: | ---: |
| `current_adp03q`, `Small` | `124/1000` | `0.496` | `1.58` | `0.0739` | `0.933` |
| `spline_q05_gap02`, `Small` | `123/1000` | `0.498` | `1.59` | `0.0926` | `0.913` |
| `current_uniform08_gap05`, `Small` | `117/1000` | `0.500` | `1.56` | `0.134` | `0.874` |
| `constant`, `Small` | `81/1000` | `0.504` | `1.07` | `1.07` | `0.000` |
| `current_adp03q`, `Primary` | `71/1000` | `0.284` | `0.89` | `0.0694` | `0.966` |
| `current_adp03q`, `Large` | `46/1000` | `0.167` | `0.53` | `0.0660` | `0.974` |

Full summary and plot:
[summary](s56_1k_predetrend_detrending_method_sweep_full_pdo/summary.md),
[pareto PNG](s56_1k_predetrend_detrending_method_sweep_full_pdo/detrending_method_pareto.png).

### Extended small-aperture sweep

| method | rows above SNR 7 | median depth retention | median empirical SNR | residual trend ptp | trend reduction |
| --- | ---: | ---: | ---: | ---: | ---: |
| `oracle_adp03q`, `Small` | `126/1000` | `0.504` | `1.61` | `0.0739` | `0.933` |
| `oracle_q05_gap02`, `Small` | `125/1000` | `0.504` | `1.60` | `0.0929` | `0.914` |
| `current_adp03q`, `Small` | `124/1000` | `0.496` | `1.58` | `0.0739` | `0.933` |
| `spline_q05_gap02`, `Small` | `123/1000` | `0.498` | `1.59` | `0.0926` | `0.913` |
| `savgol06_p2_gap05`, `Small` | `122/1000` | `0.496` | `1.59` | `0.0743` | `0.933` |
| `median06_gap05`, `Small` | `122/1000` | `0.498` | `1.58` | `0.100` | `0.910` |
| `savgol10_p2_gap05`, `Small` | `121/1000` | `0.498` | `1.58` | `0.0919` | `0.919` |
| `median10_gap05`, `Small` | `120/1000` | `0.501` | `1.55` | `0.157` | `0.849` |
| `median20_gap05`, `Small` | `116/1000` | `0.502` | `1.51` | `0.230` | `0.785` |

Against `current_adp03q`, `oracle_adp03q` adds only `2` new SNR-qualified
rows and loses none; its median SNR gain is `0.008`. This means that masking
the exact injected transit out of the cotrend fit barely changes the result.
The filter methods do not improve the Pareto point: short-window Savitzky-
Golay is close to current ADP, while rolling medians and wider filters leave
more residual trend and recover fewer SNR-qualified rows.

Extended summary and plot:
[summary](s56_1k_predetrend_detrending_method_sweep_extended_small_pdo/summary.md),
[pareto PNG](s56_1k_predetrend_detrending_method_sweep_extended_small_pdo/detrending_method_pareto.png).

### Envelope/local-filter follow-up

| method | rows above SNR 7 | median depth retention | median empirical SNR | residual trend ptp | trend reduction |
| --- | ---: | ---: | ---: | ---: | ---: |
| `oracle_adp03q`, `Small` | `126/1000` | `0.504` | `1.61` | `0.0739` | `0.933` |
| `median03_gap05`, `Small` | `125/1000` | `0.497` | `1.59` | `0.0486` | `0.955` |
| `current_adp03q`, `Small` | `124/1000` | `0.496` | `1.58` | `0.0739` | `0.933` |
| `savgol03_p2_gap05`, `Small` | `123/1000` | `0.486` | `1.57` | `0.0659` | `0.939` |
| `pctl60_06_gap05`, `Small` | `123/1000` | `0.501` | `1.58` | `0.123` | `0.890` |
| `pctl75_03_gap05`, `Small` | `122/1000` | `0.498` | `1.58` | `0.125` | `0.886` |

The shorter rolling median is the only new method that improves the Pareto
point, and the improvement is small: one extra SNR-qualified row relative to
current ADP, with better residual-trend suppression. Percentile/envelope
filters do not produce a signal-preservation jump; the high-percentile
baselines leave more residual low-frequency structure and recover fewer
SNR-qualified rows than the current ADP spline.

Envelope summary and plot:
[summary](s56_1k_predetrend_detrending_method_sweep_envelope_pdo/summary.md),
[pareto PNG](s56_1k_predetrend_detrending_method_sweep_envelope_pdo/detrending_method_pareto.png).

### 94-row priority subset

| product | rows above SNR 7 | median depth retention | median empirical SNR | residual trend ptp | trend reduction |
| --- | ---: | ---: | ---: | ---: | ---: |
| `current_adp03q`, `Small` | `83/94` | `0.648` | `9.19` | `0.0397` | `0.941` |
| `spline_q05_gap02`, `Small` | `82/94` | `0.650` | `9.11` | `0.0507` | `0.922` |
| `current_uniform08_gap05`, `Small` | `76/94` | `0.653` | `8.95` | `0.0805` | `0.878` |
| `constant`, `Small` | `49/94` | `0.656` | `7.18` | `0.580` | `0.000` |

The priority BLS sweep independently points to the same practical next step.
On these `94` rows, medium ADP at `5k` periods found only `7/94` exact or
harmonic matches, while the small-aperture pair with a dense `200k` period grid
found `77/94` and put `73/94` at top rank.

Priority summaries:
[method summary](s56_1k_predetrend_detrending_method_sweep_priority_pdo/summary.md),
[BLS sweep table](s56_1k_predetrend_priority_recovery_sweep_pdo/priority_sweep_summary/priority_sweep_summary.csv).

## Interpretation

Detrending/product SNR is still a real problem for the full faint-heavy `1k`
injection sample. The best tested product still has only `124/1000` rows above
empirical SNR `7`, and the median injected depth retention is about `0.50`.
That is not good enough for a final completeness or human-teacher sample.

However, the method sweep did not find a detrending family that preserves most
signals while still removing most trends. Looser splines, polynomials, most
rolling medians, Savitzky-Golay filters, percentile/envelope filters, and a
constant baseline preserve only marginally more truth-window depth while
leaving more long-timescale structure and yielding fewer SNR-qualified rows.
The oracle transit-masked spline barely improves over current ADP, so the
dominant loss is not simply that the cotrend fit is swallowing the transit. The
0.3 d rolling median is now a useful near-tie or possible follow-up branch, but
its full-set gain is only `125/1000` versus `124/1000` SNR-qualified rows, not
a completeness-changing solution.

The largest win is aperture choice, not a new detrending method. The small raw
aperture carries much more usable injected-signal SNR than the medium or large
apertures. The medium-aperture ADP product used in the first human queue was
therefore a poor default for faint WD injections.

## Decision

Do not scale the current medium-ADP pre-detrend queue to `10k`.

The follow-up full-`1k` recovery-mode sweep on the small-aperture pair
`DET_FLUX_ADP_SML + DET_FLUX_SML` confirmed the priority-subset result. With a
dense `200k` BLS period grid, strict top-rank recovery improves to `137/1000`
and exact/top-N/harmonic recovery improves to `177/1000`, compared with
`32/1000` and `68/1000` for the old medium-ADP `5k` queue. The verified
human-check queue has therefore been rebuilt around small-aperture LEO reports:
`reports/stage5_validation/s56_1k_predetrend_small_pair_200k_review_queue_pdo/`.

Medium and large apertures should remain diagnostic evidence, but they should
not be the default BLS/LEO evidence path for faint injected WD-transit recovery.
The 0.3 d rolling-median small-aperture product is worth a future BLS/LEO
comparison if we decide to add another evidence branch, but it should not
replace the current small-aperture ADP/canonical pair solely on truth-window
metrics.

For science completeness, continue treating the low-SNR faint rows as recovery
map measurements, not as teacher positives. Pixel-level injections remain a
calibration subset, not the dense-grid path.
