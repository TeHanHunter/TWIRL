# S56 Peak Ranker Injection-Recovery Benchmark

Date: `2026-07-01`

## Question

How well does the small injected-truth BLS peak ranker recover the correct
injected period/ephemeris from the list of BLS candidate peaks?

## Data Product

This benchmark uses the coverage-first all-host S56 injected peak table on PDO:

```text
reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_training/s56_allhost_injection_bls_peaks.csv
```

The trained ranker summary is:

```text
reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_ranker_pdo/summary.json
```

The post-BLS gate summary is:

```text
reports/stage5_validation/s56_allhost_predetrend_batman_periodradius_grid_sharded_pdo/peak_training_gate_pdo/summary.json
```

The positive label is `signal_peak`: a BLS peak must match the injected period
or accepted harmonic and the injected transit window. Injection truth columns
are used only for labels/audit, not as model features.

## Model

The ranker is a lightweight logistic/softmax model over BLS and light-curve
features available for real candidates: BLS rank, period, duration, depth,
depth SNR, SDE, log power, cadence-cleaning diagnostics, baseline, Tmag, and
aperture.

It is an ephemeris/peak selector, not an astrophysical classifier.

## Main Result

BLS itself is the dominant limitation. Out of `19,071` evaluated all-host
injections, the correct injected signal appears in the candidate peak set for
`5,085` injections (`26.7%`). The remaining `13,987` injections are not in the
top-20 candidate set, so no peak ranker can recover them without improving the
search/detrending stage.

Failure-mode counts:

| Mode | Count |
|---|---:|
| Top-1 recovered before ranker | `3,803` |
| Ranker-fixable | `1,282` |
| High-observability not in top-20 | `1,360` |
| Low-observability not in top-20 | `12,275` |
| Low-cadence not in top-20 | `352` |

## Held-Out Test Split

The held-out test split has `3,814` injections, of which `1,026` are rankable
because the correct injected signal appears in the candidate peak set.

| Selection | Correct @1 | Correct @2 | Correct @3 | Correct @5 | Correct @10 | Correct @20 |
|---|---:|---:|---:|---:|---:|---:|
| Raw BLS rank | `757` | `768` | `795` | `826` | `867` | `936` |
| SDE rank | `703` | `715` | `738` | `766` | `819` | `881` |
| Ranker | `750` | `788` | `815` | `851` | `913` | `981` |

As a fraction of the rankable test injections:

| Selection | Correct @1 | Correct @2 | Correct @3 | Correct @5 | Correct @10 | Correct @20 |
|---|---:|---:|---:|---:|---:|---:|
| Raw BLS rank | `73.8%` | `74.9%` | `77.5%` | `80.5%` | `84.5%` | `91.2%` |
| SDE rank | `68.5%` | `69.7%` | `71.9%` | `74.7%` | `79.8%` | `85.9%` |
| Ranker | `73.1%` | `76.8%` | `79.4%` | `82.9%` | `89.0%` | `95.6%` |

## Interpretation

The current ranker is useful for ranking a short list of ephemerides, especially
top-2 to top-5, but it is not yet a better single-period selector than raw BLS
rank on the held-out test split. For top-1, it is slightly worse than raw BLS
rank (`750` vs `757`) but better than pure SDE (`750` vs `703`).

The current human triage design, showing the ranker-selected top `3`
ephemerides per target before LEO/human review, is therefore appropriate. A
single forced best-period decision is not yet justified.

The next scientific improvement should focus on increasing the candidate-set
recall of BLS/search itself, because `73.3%` of all-host injections never put
the correct signal into the top-20 peak set. The ranker can only fix the
`~6.7%` ranking-loss population.
