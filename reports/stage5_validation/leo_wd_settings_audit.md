# WD-Tuned LEO-Vetter Settings Audit

Date: 2026-06-22

Remote review jobs on PDO import LEO-Vetter from:

`/pdo/users/tehan/LEO-Vetter-twirl/leo_vetter/`

The S56 review queue builder uses the WD-specific path:

- `tlc.compute_flux_metrics(star, cap_b=False)`
- `tlc.metrics["Rs"] = 0.013`
- `check_thresholds_wd(tlc.metrics, "FA")`
- `check_thresholds_wd(tlc.metrics, "FP")`

## WD-Specific Thresholds Confirmed On PDO

- `N_transit = 2`: allows long-period single-sector WD candidates with only two observed transits.
- `V_shape = 1e6`: V-shape rejection is effectively disabled for WD hosts.
- `size = 22.4 R_earth`: explicit large-object FP boundary near `2 R_J`.
- `R_comp_max_R_jup = 2.0`: maximum companion radius used by the WD chord-duration envelope.
- `q_ratio_min = 0.4`: loosened duration-ratio threshold for WD/grazing geometries.

## WD Override Behavior Checked

- `vshaped_wd` returns `False` for a WD-like `b=0.9`, `Rp/Rs=8` case.
- `invalid_transits_wd` still rejects `new_N_transit < 2`, but no longer rejects a two-transit candidate solely because pruned `new_MES` is low.
- `large` does not reject the injection-grid ceiling (`16.8 R_earth`, about `1.5 R_J`) and does reject `23.5 R_earth` (`>2 R_J`).
- `non_unique_wd` keeps the primary-vs-noise and primary-vs-tertiary checks but drops the stock positive-feature `MS3` clause.
- `unphysical_duration_wd` uses a WD chord-sum duration envelope rather than the main-sequence point-companion duration formula.

## Remaining Caveat

Several LEO tests remain intentionally inherited from the stock thresholds
(`weak`, `bad_shape`, `chases`, `dmm`, `single_event`, `bad_fit`, `sinusoidal`,
`asymmetric`, `chi`, `data_gapped`, `odd_even`, `secondary`, and optional
`offset`). They are not all independently validated for WD transits yet, so
the current WD tuning is best treated as a permissive first-pass vetter for
human triage rather than a final automated decision rule.

## 2026-06-24 Threshold Smoke Result

The S56 injected-row smoke test in
[summary](s56_10k_predetrend_dense_bls_map_pdo/leo_wd_tuning_smoke/summary.md)
compares possible WD review pass-through presets against BLS exact/top-N/
harmonic recovery on the `1,000` injection rows with LEO metrics.

- Current LEO `PC/FP`: precision `98.7%`, recall `42.9%`.
- `wd_review_high_purity`: precision `91.8%`, recall `50.3%`, adds `13`
  BLS-recovered rows and `7` BLS-missed rows to review.
- `wd_review_balanced`: precision `86.0%`, recall `55.4%`.
- `wd_review_aggressive`: precision `64.9%`, recall `61.6%`.

Interpretation: the physical guard that matters most is primary-event
dominance over secondary/tertiary/positive events. Loosening shape/MES cuts is
reasonable for under-resolved 200 s WD transits, but dropping the primary-
dominance guard quickly admits wrong-ephemeris BLS misses. Implement the first
tune as an explicit `WD_REVIEW` / recovered-but-not-clean class, not as direct
promotion to `PC`.
