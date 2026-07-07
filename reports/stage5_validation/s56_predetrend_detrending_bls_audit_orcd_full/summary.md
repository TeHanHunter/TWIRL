# Raw-Flux Detrending-Strength BLS Audit

This audit starts from stored raw aperture flux, applies each detrending method, reruns BLS, and labels peaks against injected BATMAN truth.

## Branch Summary

| method                  | raw_aperture | n_injections | top1_exact_frac | topn_exact_or_harmonic_frac | median_depth_retention | median_det_trend_ptp | median_sde_rank1 |
| ----------------------- | ------------ | ------------ | --------------- | --------------------------- | ---------------------- | -------------------- | ---------------- |
| spline_q015_gap02       | Small        | 3000         | 0.5113          | 0.571                       | 0.4337                 | 0.02098              | 31.01            |
| spline_q025_gap02       | Small        | 3000         | 0.5087          | 0.569                       | 0.4376                 | 0.02304              | 32.79            |
| spline_q02_gap02        | Small        | 3000         | 0.507           | 0.5693                      | 0.4363                 | 0.02147              | 31               |
| current_adp03q          | Small        | 3000         | 0.505           | 0.5663                      | 0.4386                 | 0.02323              | 33.86            |
| spline_q05_gap02        | Small        | 3000         | 0.4993          | 0.562                       | 0.44                   | 0.03005              | 36.09            |
| spline_q08_gap02        | Small        | 3000         | 0.4947          | 0.556                       | 0.4415                 | 0.03737              | 39.17            |
| current_uniform08_gap05 | Small        | 3000         | 0.4833          | 0.5447                      | 0.4415                 | 0.04795              | 43.99            |
| spline_q015_gap02       | Primary      | 3000         | 0.411           | 0.477                       | 0.3441                 | 0.02424              | 34.21            |
| spline_q02_gap02        | Primary      | 3000         | 0.404           | 0.4673                      | 0.3462                 | 0.02683              | 36.31            |
| spline_q025_gap02       | Primary      | 3000         | 0.3957          | 0.4583                      | 0.3471                 | 0.02885              | 39.24            |
| current_adp03q          | Primary      | 3000         | 0.391           | 0.4537                      | 0.3486                 | 0.02994              | 40.58            |
| spline_q05_gap02        | Primary      | 3000         | 0.3773          | 0.438                       | 0.3509                 | 0.04835              | 43.7             |
| spline_q08_gap02        | Primary      | 3000         | 0.3613          | 0.4263                      | 0.3512                 | 0.0691               | 49.1             |
| current_uniform08_gap05 | Primary      | 3000         | 0.3577          | 0.4207                      | 0.3513                 | 0.08916              | 54.43            |

## Notes

- `top1_exact_frac` is the strict deployable recovery rate: BLS rank 1 matches the injected period and transit window.
- `topn_exact_or_harmonic_frac` includes exact and accepted harmonic matches among all retained peaks.
- `median_depth_retention` measures injected-minus-original depth after detrending; values near 1 preserve the signal.
- This script does not use transit-window masking for deployable branches.
