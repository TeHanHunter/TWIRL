# S56 Pixel-Level vs Pre-Detrend Injection Feasibility

Date: 2026-06-22

## Decision

For the immediate scientifically useful human-vetting sample, use
raw-aperture light-curve injections before TWIRL-FS detrending. This directly
tests the current risk that canonical or ADP detrending suppresses short WD
occultations, and it reuses the existing BLS, WD-tuned LEO, and browser-vetting
stack.

Pixel-level injection is possible on PDO and now has a working smoke wrapper,
but it remains a calibration-subset project rather than the right way to make
the next 1k/10k human-vetting queue.

## Pixel-Level Feasibility

Read-only PDO probes confirmed the needed raw pixel products exist under the
user-owned S56 production tree:

- TGLC source pickles contain `source.flux` arrays with shape
  `(n_cadences, 150, 150)`.
- Matching ePSF arrays exist under each CCD `epsf/` directory.
- The installed TGLC fork exposes `tglc.scripts.epsfs.fit_epsf_for_source`
  and `tglc.light_curve.generate_light_curves`, so a modified source pickle
  can be refit and re-extracted without changing shared PDO inputs.

Implemented a narrow pixel-level smoke wrapper:

- [pixel injection smoke script](../../scripts/stage3_injections/run_s56_pixel_injection_smoke.py)
- Source tested: PDO `source_7_7.pkl`, ePSF `epsf_7_7.npy`
- Target auto-selected: TIC `389963197`, Tmag `17.50`, edge margin `52 px`
- Reuse-ePSF smoke: `600` cadences, target pixel injection before TGLC
  extraction, primary-aperture median delta `-185.3` in transit and `0.0`
  out of transit
- ePSF-refit smoke: `200` cadences, injection before ePSF refit and TGLC
  extraction, primary-aperture median delta `-116.4` in transit and
  `-3.9e-4` out of transit
- Full-chain reuse-ePSF smoke: pixel injection -> TGLC HDF5 -> TWIRL-FS FITS
  with ADP columns -> BLS summary. A `600`-cadence, `P=0.6 d` case with
  `2` good in-transit cadences produced the expected extracted-aperture
  transit but was not BLS-recovered; the strongest BLS peak landed at
  `0.3628 d` with SDE `3.94`.
- Full-chain short-period control: a `2000`-cadence, `P=0.2 d` case with
  `11` good in-transit cadences also failed the period-recovery criterion.
  Canonical `DET_FLUX` reached SDE `7.30`, but at `0.6253 d` with a
  `30 min` duration box, not the injected `2 min` signal. ADP peaked at
  `0.3747 d`.

The refit path captures extraction, ePSF, crowding, aperture, and centroid
effects. It is the paper-grade calibration layer, but it is too expensive for
the dense human-vetting grid and should be run as a controlled calibration
subset after the raw-aperture detrend/search behavior is understood. The
full-chain BLS smokes show that the immediate failure mode is now search and
detrending sensitivity to very short injected events, not the ability to inject
at pixel level.

## Implemented Pre-Detrend Path

Implemented:

- [pre-detrend injection builder](../../scripts/stage3_injections/make_s56_predetrend_injection_set.py)
- [PDO 1k launcher](../../scripts/stage5_validation/run_s56_1k_predetrend_review_pdo.sh)
- [focused tests](../../tests/test_predetrend_injections.py)

The builder injects finite-exposure BATMAN transits into raw TGLC `RawFlux`,
then reruns canonical `twirl-fs-v2` and ADP `twirl-fs-v2-adp03q` detrending.
The output HDF5 is compatible with the existing Stage 5 review-queue builder.

## Validation

Small PDO smoke tests passed:

- WD 1856 only: `5/5` accepted raw-flux injections, no raw/detrend skips.
- Non-benchmark TIC `47756522`: `10/10` accepted; BLS recovered `4/10`;
  two LEO reports rendered with `0` metric or plot errors.

Completed PDO 1k run:

- input raw targets discovered: `19,072`
- accepted injections: `1,000`
- skipped: `99` nonpositive raw-primary baselines; `0` read, merge, epoch,
  model, or detrend failures
- BLS apertures: `DET_FLUX_ADP`, `DET_FLUX`
- BLS recovered: `33 / 1000`
- representative apertures: `DET_FLUX=741`, `DET_FLUX_ADP=259`
- LEO reports: `1000 / 1000` attempted and rendered
- LEO errors: `0`; LEO plot errors: `0`
- LEO classes: `FA=974`, `PC=25`, `FP=1`

Completed ADP-only rebuild from the same pre-detrend injected HDF5:

- BLS/LEO aperture: `DET_FLUX_ADP` only
- review rows: `1,000`
- BLS recovered: `32 / 1000`
- LEO reports: `1000 / 1000` attempted and rendered
- LEO errors: `0`; LEO plot errors: `0`
- LEO classes: `FA=974`, `PC=25`, `FP=1`

Local browser vetting is serving the ADP-only queue on
`http://127.0.0.1:5005/`.

## BLS And Detrending Diagnostic

The ADP-only pre-detrend queue is not yet a completeness sample. It is a
diagnostic showing that the current detrending/search stack loses many injected
signals before BLS can rank them.

Strict BLS recovery is not random:

- total recovery: `32 / 1000`
- by magnitude: `22 / 97` for `Tmag < 18`, `4 / 684` for `Tmag >= 19`
- by period: `27 / 458` for `P < 1 d`, `0 / 400` for `P >= 2 d`
- deep, short-period subset: `17 / 177` for `model_depth > 0.3` and `P < 1 d`
- deep, long-period subset: `0 / 240` for `model_depth > 0.3` and `P >= 2 d`

The injected signal survives in raw aperture flux but is strongly attenuated
after TWIRL-FS detrending:

- raw primary aperture median depth retention: `0.82`
- `DET_FLUX_ADP` median depth retention: `0.28`
- `DET_FLUX_ADP` median retention for `Tmag >= 19`: `0.23`
- `DET_FLUX_ADP` median retention for deep faint cases: `0.22`

Measured signal SNR after detrending explains most of the low strict recovery:

- empirical multi-transit SNR `< 7`: `929 / 1000` ADP rows
- BLS misses with empirical multi-transit SNR `< 7`: `926 / 968`
- recovery by empirical multi-transit SNR:
  - `< 1`: `0 / 527`
  - `1-3`: `1 / 277`
  - `3-5`: `1 / 96`
  - `5-7`: `1 / 29`
  - `7-10`: `5 / 24`
  - `10-20`: `10 / 30`
  - `> 20`: `14 / 17`

A BLS duration-grid experiment does not rescue the run:

- current `3-30 min` grid: `32 / 1000`
- WD-short `1-10 min` grid: `38 / 1000`
- WD-short `1-6 min` grid: `35 / 1000`

So the immediate blocker is not just that BLS lacks a 1-2 minute duration grid.
The dominant issue is that most faint/long-period injected signals are below
useful empirical SNR after `DET_FLUX_ADP`, and broad non-injection structure or
harmonics often win the BLS ranking for the remaining high-SNR misses. A
200-row `inferred_aperture` baseline smoke slightly improved median ADP depth
retention (`0.31` vs `0.28`) and empirical SNR (`1.14` vs `0.89`) but did not
improve strict recovery (`7 / 200`), so the failure is not just a raw-baseline
normalization mistake. The review-queue builder now preserves strict top-peak
recovery separately from exact top-N recovery and harmonic top-N matches, so
the next rebuilt queue can quantify how many misses are true non-detections
versus BLS ranking/alias failures. The first PDO rebuilds show:

- `n_periods=5000`: strict top-1 `32 / 1000`, exact top-N non-top-1
  `11 / 1000`, harmonic top-N `25 / 1000`, unmatched `932 / 1000`
- `n_periods=200000`: strict top-1 `50 / 1000`, exact top-N non-top-1
  `10 / 1000`, harmonic top-N `17 / 1000`, unmatched `923 / 1000`
- at `n_periods=200000`, empirical SNR `>20` is `17 / 17` top-1 recovered,
  so dense BLS fixes high-SNR grid/ranking misses but not the dominant
  low-SNR population
- aperture signal survival favors the small aperture: `DET_FLUX_ADP_SML`
  has median empirical multi-transit SNR `1.58` with `124` rows above SNR 7,
  versus `0.89` and `71` rows for medium `DET_FLUX_ADP`

Diagnostic artifacts:

- [signal survival CSV](../stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_signal_survival_diagnostics.csv)
- [signal survival + SNR CSV](../stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_signal_survival_snr_diagnostics.csv)
- [BLS duration-grid CSV](../stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/bls_duration_grid_experiment.csv)
- [diagnostic PNG](../stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_adp_bls_recovery_diagnostics.png)
- [SNR diagnostic PNG](../stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_pdo/predet_adp_bls_recovery_snr_diagnostics.png)
- [5k recovery-mode summary](../stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_topn5000_pdo/recovery_mode_summary/summary.json)
- [200k recovery-mode summary](../stage5_validation/s56_1k_predetrend_batman_adp_only_review_queue_topn_pdo/recovery_mode_summary/summary.json)
- [recovery-mode sweep driver](../../scripts/stage5_validation/run_injection_recovery_mode_sweep.py)
- [PDO recovery-mode sweep launcher](../../scripts/stage5_validation/run_s56_predetrend_recovery_mode_sweep_pdo.sh)
- [SNR-stratified human-review subset](../stage5_validation/s56_1k_predetrend_batman_adp_only_snr_stratified_review_queue/review_queue.csv)
- [SNR-stratified local vetting launcher](../../scripts/stage5_validation/run_s56_snr_stratified_vetting_app_local.sh)
- [inferred-baseline smoke](../stage5_validation/s56_200_predetrend_batman_inferred_baseline_bls_smoke/)

## Caveats

This 1k product is a real detrending/search/vetter recovery test, but it is not
yet the final occurrence-rate completeness product.

- It is LC-level before detrending, not pixel-level before extraction.
- The current grid is balanced in period-depth, not radius; tiny-body cases are
  therefore sparse after finite-cadence constraints.
- BLS recovery is intentionally measured through the current stack and is low;
  that is an empirical result, not a UI failure.
- Pixel-level full-chain smokes are currently not recovered by BLS even when
  the extracted aperture flux contains the injected transit. This should be
  treated as a search/detrending diagnostic before scaling to a large pixel
  calibration subset.
- Pixel-level injections are feasible, but should be run next on a calibration
  subset to quantify where raw-aperture injections misestimate
  extraction/crowding losses.
