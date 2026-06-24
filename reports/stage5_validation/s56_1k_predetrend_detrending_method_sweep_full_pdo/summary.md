# Pre-Detrend Detrending Method Sweep

Empirical SNR threshold: `7`.

## Best Rows

| method                  | raw_aperture | n    | median_depth_retention | p16_depth_retention | p84_depth_retention | median_multi_snr | n_snr_ge_threshold | frac_snr_ge_threshold | median_det_trend_ptp | median_trend_reduction_frac | median_tmag |
| ----------------------- | ------------ | ---- | ---------------------- | ------------------- | ------------------- | ---------------- | ------------------ | --------------------- | -------------------- | --------------------------- | ----------- |
| current_adp03q          | Small        | 1000 | 0.496                  | 0.308               | 0.7128              | 1.583            | 124                | 0.124                 | 0.07391              | 0.9327                      | 19.44       |
| spline_q05_gap02        | Small        | 1000 | 0.498                  | 0.3103              | 0.7128              | 1.586            | 123                | 0.123                 | 0.09264              | 0.9134                      | 19.44       |
| current_uniform08_gap05 | Small        | 1000 | 0.5002                 | 0.312               | 0.7146              | 1.563            | 117                | 0.117                 | 0.1339               | 0.8738                      | 19.44       |
| constant                | Small        | 1000 | 0.5036                 | 0.3145              | 0.7153              | 1.067            | 81                 | 0.081                 | 1.07                 | 0                           | 19.44       |
| current_adp03q          | Primary      | 1000 | 0.2845                 | 0.1258              | 0.5359              | 0.8915           | 71                 | 0.071                 | 0.0694               | 0.9661                      | 19.44       |
| spline_q05_gap02        | Primary      | 1000 | 0.2861                 | 0.1262              | 0.5394              | 0.8801           | 67                 | 0.067                 | 0.1097               | 0.9437                      | 19.44       |
| current_uniform08_gap05 | Primary      | 1000 | 0.2868                 | 0.1267              | 0.5426              | 0.8481           | 66                 | 0.066                 | 0.178                | 0.9077                      | 19.44       |
| current_adp03q          | Large        | 1000 | 0.1665                 | 0.04748             | 0.3536              | 0.5261           | 46                 | 0.046                 | 0.06598              | 0.9744                      | 19.44       |
| spline_q05_gap02        | Large        | 1000 | 0.1671                 | 0.04754             | 0.3551              | 0.5211           | 42                 | 0.042                 | 0.1158               | 0.9505                      | 19.44       |
| current_uniform08_gap05 | Large        | 1000 | 0.1676                 | 0.04763             | 0.3554              | 0.5093           | 39                 | 0.039                 | 0.1963               | 0.9198                      | 19.44       |
| constant                | Primary      | 1000 | 0.2889                 | 0.1287              | 0.5445              | 0.447            | 29                 | 0.029                 | 2.075                | 0                           | 19.44       |
| constant                | Large        | 1000 | 0.1711                 | 0.04892             | 0.3553              | 0.2123           | 14                 | 0.014                 | 2.457                | 0                           | 19.44       |

## Interpretation Aid

- `n_snr_ge_threshold` is the number of rows whose injected-minus-original truth-window signal reaches the SNR threshold after detrending.
- `median_depth_retention` near `1` means the method preserves the injected BATMAN signal depth.
- `median_det_trend_ptp` is the robust 90-10% range of binned original detrended flux; lower is better.
- `median_trend_reduction_frac` compares that trend range to the raw normalized flux; higher is better.
