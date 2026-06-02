# twirl_plan.md

This document is the executable software and survey plan for TWIRL.

## Current Status Snapshot

- TWIRL is a standalone repo and remains separate from `TWIRL_proposal`.
- The core survey is limited to `200 s` TESS FFIs, which means `Sector >= 56`.
- `WD 1856+534` is the benchmark target, and the first fixed benchmark is `Sector 56`, orbits `119` and `120`, `cam4/ccd1`.
- The seed WD catalog is a local external dependency, not a git-tracked repo asset.
- Light curves will be gathered for all WDs in the full parent catalog. The statistical occurrence-rate sample (the locked denominator) will be defined separately after light-curve collection and QA are complete.
- The search is interpretable-first: build a transparent periodic + dip baseline before considering ML triage.
- Stage 1 Sector 56 benchmark (orbits 119 + 120, all 16 CCDs): TGLC lightcurves, detrend, and HLSP FITS production are complete. The cam3/orbit-120 cadence-mismatch blocker has been resolved by `scripts/stage1_lightcurves/fix_tglc_quat_cadence_mismatch.py` (drops the TGLC-only orbit-boundary cadence from every equal-length dataset in each h5; the quat file is left untouched as the authoritative pointing record). After re-running detrend + HLSP at the sector level, **`19,072 / 19,086 = 99.93%`** HLSP FITS are written to `/pdo/users/tehan/tglc-deep-catalogs/hlsp_s0056/` and validation passes with `0/19086 files had issues`. `14` targets remain unrecovered (root cause TBD; pull from `[hlsp skip]` lines in the r2 run log). Next session: Stage 1.6 WD 1856 photometric QA — read `267574918` HLSP FITS, plot stitched LC across orbits 119/120, and record per-aperture RMS/MAD for the `1×1`/`3×3`/`5×5` apertures.
- `2026-04-27`: tglc-mki GPU ePSF path is unblocked. Built `/pdo/users/tehan/twirl-gpu-venv` (`--system-site-packages` + `cupy-cuda11x`) and benchmarked `tglc epsfs` on `pdogpu6`: S56 cam4/ccd1 finished `196/196` in **`30:42`** vs the `4:01:31` CPU baseline → **`~7.9× speedup`** with outputs bit-equivalent to CPU at machine precision. Future Stage 1 ePSF runs should default to this venv on `pdogpu6` (8 GPUs); see progress log §1.4 (`2026-04-27`).
- `2026-04-30`: **S56 GPU production complete and validated.** New tree at `/pdo/users/tehan/tglc-gpu-production/hlsp_s0056/`: `19,068 / 19,086 = 99.91%` HLSP FITS, validation `4 / 19086` issues. WD 1856 parity vs the frozen `tglc-deep-catalogs/hlsp_s0056/` baseline within `~2%` per aperture (medium aperture MAD-RMS `0.1111` vs `0.1087`). The new tree is the science tree of record going forward; the deep-catalogs benchmark stays frozen for reference. Mass-production scale-out to S57+ uses [run_sector_gpu_production.sh](../scripts/stage1_lightcurves/run_sector_gpu_production.sh) + [run_s56_post_lc_chain.sh](../scripts/stage1_lightcurves/run_s56_post_lc_chain.sh) (parameterize per sector). See progress log §1.4 entries for `2026-04-30`.
- `2026-05-01`: Switched S58+ to a two-node pipelined architecture: pdogpu1 prepares catalogs+cutouts (no GPU), pdogpu6 finalizes ePSFs+lightcurves+post-LC chain one sector at a time (sequential for QC). Driver gained `--stages CSV` flag; new prep/finalize scripts and worker loops at `scripts/stage1_lightcurves/`. Projected 39-sector campaign (S58–S96) wall time **`~12 days`** (vs `~32 d` single-node). S57 left undisturbed on pdogpu6 with the older monolithic launcher; new pipeline kicks in at S58.
- `2026-05-06`: **Mass production in flight.** S56–S60 DONE (`5/41` sectors finalized via the new pipeline; HLSP totals: S56 `19,068`, S57 `16,337`, S58 `14,656`, S59 `13,467`, S60 `16,143`). S61 finalizing (currently ~88% of step A; galactic-plane sectors S58–S64 run `~10–12 hr` per finalize wall instead of the `~5–6 hr` budget on sparser fields). pdogpu1 prep stays ahead of pdogpu6 finalize: S62 fully prepped, S63/S64 in flight. Live status at `/pdo/users/tehan/tglc-gpu-production/STATUS.md` rendered every `10 min` by [sector_status.py](../scripts/stage1_lightcurves/sector_status.py) in tmux `twirl-status`. Realistic completion ETA for S61–S94 backlog: **~16 days, finishing ~2026-05-22**.
- `2026-05-06`: Per-sector HLSP QC PDF system shipped. [qc_sector_pdf.py](../scripts/stage1_lightcurves/qc_sector_pdf.py) emits a 7-page PDF per sector (Tmag distribution, galactic-Aitoff sky-coverage map, photometric-precision panel using DMAD recipe overlaid with σ_base from Sullivan+2015, example LCs in 5 Tmag bins with QUALITY-flag + 5σ-clip annotations). Wired into [finalize_sector_gpu.sh](../scripts/stage1_lightcurves/finalize_sector_gpu.sh) as step H so the PDF + companion `.npz` are auto-generated after every HLSP run. Local sync via [sync_qc_pdfs.sh](../scripts/sync_qc_pdfs.sh) (rsync wrapper, PDF-only by default).
- `2026-05-06`: **Pipeline unit-conversion audit complete.** TICA `e-/cadence` → TGLC `e-/s` conversion at `aperture_photometry.py:120` is correct; TESS magnitude zeropoint at `15,000 e-/s @ T=10` matches the TESS Instrument Handbook. Gaia → TESS magnitude polynomial in `ffi.py:202–215` matches Stassun+2019 ApJS 243 Eq. 1; the ePSF design matrix uses dimensionless `tess_flux_ratio` so absolute units cancel in the fit. **Caveat**: HLSP `SAP_FLUX`/`DET_FLUX` are median-normalized relative flux (not e-/s) because QLP `lctools/bin/hlsp.py:258 flux_from_mag` subtracts the median magnitude; `SAP_BKG` is `e-/cadence/pixel`; FITS columns 3/4/10 carry no `TUNIT`. Recovery factor for true e-/s: `15000 * 10**(-0.4*(TESSMAG-10))`. Decision on whether to fork `flux_from_mag` or post-process HLSPs deferred.
- `2026-05-12`–`2026-05-13`: **Pre-talk sprint Week-1 done end-to-end on the existing QLP HLSPs.** Heuristic vetter ([src/twirl/vetting/heuristic.py](../src/twirl/vetting/heuristic.py)) moves WD 1856 from blind rank `1258 / 19,040` to planet-regime rank **`9 / 5,403`** (140× improvement) using three physics-motivated cuts on existing BLS columns — duration upper-envelope at `R_comp_max = 2 R_jup` (chord-sum, M_WD generous at `0.4 M☉`), period-alias rejection at `P > 0.10 d`, and period-cluster pile-up detection. LEO-Vetter integration via [TeHanHunter/LEO-Vetter@wd-host-tuning](https://github.com/TeHanHunter/LEO-Vetter/tree/wd-host-tuning) (new TWIRL fork) — adds WD-host preset + 4 override functions (`vshaped_wd` disabled, `unphysical_duration_wd` chord-sum, `invalid_transits_wd` drops post-pruning MES requirement at 200 s cadence, `non_unique_wd` drops MS3 sig_pos clause); WD 1856 labels as **PC**. Outputs in [reports/stage5_validation/leo_vetter_s56_top50/](../reports/stage5_validation/leo_vetter_s56_top50/) with one PDF per TIC named `<CLASS>_rank<NN>_tic<TIC>_T<MAG>_P<P>d.pdf` for fast triage. BLS v2 config ([src/twirl/search/bls.py:BLSConfig](../src/twirl/search/bls.py)) with `p_min=0.12 d`, denser duration grid `[3,4,5,6,8,10,13,16,20,30]`, `n_peaks=10` was re-run on S56 in `5.87 hr` wall on `pdogpu1`; WD 1856 blind rank improved `1258 → 560` but the heuristic vetter already does the work BLS v2 also does, so the post-vetting rank is essentially unchanged.
- `2026-05-13`: **Plan A — flux-space cotrending — scaffolded and TGLC-side patched.** Root cause of the faint-end cadence loss pinned to [TGLC fork twirl/preserve-negative-flux](https://github.com/TeHanHunter/TESS_Gaia_Light_Curve/tree/twirl/preserve-negative-flux) commit `8235695`: `aperture_photometry.py:114` `flux[flux <= 0] = np.nan` clipped real Poisson + background-subtraction excursions on faint (T≥19) targets before the magnitude conversion, dropping ~50% of `QUALITY=0` cadences for T=19.5+ HLSPs. Patch preserves linear flux through the magnitude conversion and adds `RawFlux`/`RawFluxError` HDF5 datasets per aperture; validated on WD 1856 on `pdogpu1` (magnitude column bit-identical, `RawFlux` keeps 155 previously-NaN'd cadences with values in `[-2922, +22005] e-/s`). Local pipeline modules built and synthetic-tested end-to-end ([src/twirl/lightcurves/{flux_detrend,tglc_h5_reader,hlsp_writer}.py](../src/twirl/lightcurves/), [scripts/stage1_lightcurves/build_twirl_hlsp.py](../scripts/stage1_lightcurves/build_twirl_hlsp.py)); writes TWIRL HLSP FITS in QLP-compatible schema (`hlsp_twirl_tess_ffi_*`) so the existing [src/twirl/io/hlsp.py](../src/twirl/io/hlsp.py) reader and all downstream code (BLS, heuristic vetter, LEO adapter) consume them unchanged. **S56 v3 TWIRL HLSP production succeeded on `pdogpu1` after a `flux_space_detrend` UnboundLocalError fix** (commit `f0e844c`): `19,068 / 19,072 = 99.98%` ok in `7.5 min` wall; 4 failures are upstream corrupt TGLC HDF5s. WD 1856 sanity vs v2 QLP HLSP: `+155 negative-flux cadences preserved`, `+6.3% usable q=0 cadences`, transit shape preserved. Tree at `/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl/`. BLS v3 launched on `pdogpu1` (tmux `bls_s56_v3`) immediately after; v3-vs-v2 comparison is the headline result for the talk.
- `2026-05-17`: **Faint negative-flux detrending experiment shipped.** [compare_flux_detrend_models.py](../scripts/stage1_lightcurves/compare_flux_detrend_models.py) compares divisive, signed-median subtractive, absolute-median subtractive, and robust-auto subtractive variants on synthetic stress tests and a real S56 faint HDF5 sample. The real sample keeps all negative `QUALITY==0` cadences finite for all subtractive variants, but signed-median normalization can invert injected dips when the good-cadence median flux is negative. The current production candidate is now robust-auto subtractive scaling (`1 + (flux - spline) / positive_scale`), wired into [build_twirl_hlsp.py](../scripts/stage1_lightcurves/build_twirl_hlsp.py). Detailed results live in [summary](../reports/stage1_lightcurves/detrend_experiments/s56_faint_sample/summary.md) and progress log §2.6.
- `2026-05-18`: **TWIRL-FS v1 S56 product built.** The initial TWIRL-FS candidate (`twirl-fs-v1`, TWIRL Flux-Spline v1) used robust-auto subtractive scaling with `bkspace_d=0.8`, `sigma_clip=5`, and FITS files named `hlsp_twirlfs_tess_ffi_*`. Full S56 was built on `pdogpu6` at `/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v1`: `19,072 ok / 0 fail` in `17.4 min`; WD 1856 plus seven faint audit targets retained all negative quality-zero cadences as finite detrended measurements. See [methods](twirl_fs_methods.md) and progress log §2.6.
- `2026-05-19`: **TWIRL-FS v1 failed S57 dim-target QA; v2 is the active candidate.** The bad case TIC `1400694779` exposed that v1 fit one spline across the multi-day orbit gap, hit an `LSQUnivariateSpline` knot-condition failure, then silently fell back to a low-order polynomial. S57-S63 relight was stopped and normal S64+ finalize remains paused by `qc_pause.flag`. The candidate fix is `twirl-fs-v2`: independent cotrend fits across gaps larger than `0.5 d`, still using robust-auto subtractive scaling. See [methods](twirl_fs_methods.md) and progress log §2.6.
- `2026-05-22`: **S56 TWIRL-FS v2 rebuilt for collaboration handoff.** Full S56 was rebuilt on `pdogpu6` at `/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2` using the current `twirl-fs-v2` code path. The initial run wrote `19,070` products and logged two transient I/O write failures; both failed TICs were retried successfully and reopen with `METHOD=twirl-fs-v2`, `GAPSPLIT=0.5`, and cotrend diagnostics. This is the current Franklin/Michelle handoff candidate while full product QA remains pending.
- `2026-05-22`: **S56 compare-column product built for Franklin/Michelle.** The handoff tree `/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare` keeps canonical `DET_FLUX` as `twirl-fs-v2` and adds experimental `DET_FLUX_ADP` columns using `twirl-fs-v2-adp03q` (`0.3 d` quantile-knot spline, `0.2 d` adaptive gap split). Full corrected rebuild wrote `19,072 / 19,072` FITS; the 100-target QA sample has all adaptive fits in spline mode and includes all 16 CCDs. Use canonical `DET_FLUX` by default and treat adaptive columns as opt-in comparison inputs.
- `2026-06-02`: **S72-S93 prep is in its final tail, not complete.** On `pdogpu1`, S72-S90 have cutout-done markers and orbits 179-188 are complete. Latest reachable status (`2026-06-01 15:08 EDT`) had S91 orbit 189 still slowly writing (`2649/3136`) and S92 orbit 191 active but still at `0/3136`; S93 was pending. `qc_pause.flag` remains in place, so normal GPU finalize is still intentionally paused while TWIRL-FS v2 product QA remains the gate.
- `2026-05-13`: **TWIRL pivots to a Schwamb-group collaboration.** Michelle Kunimoto brings a well-tuned BLS and LEO-Vetter expertise; her student + Franklin Chen tune LEO-Vetter for WDs in parallel with our `wd-host-tuning` fork. Te Han is the LC producer + data steward (the v3 TWIRL HLSP tree shipped today is the shared survey input) and is offered lead authorship on the **occurrence-rate paper** (verbal — to be locked in writing this week); **catalog paper leadership undecided**. Injection-recovery becomes shared exploratory work with multiple approaches in parallel. See progress log [§2.5](twirl_progress_log.md) for the meeting record and the [Collaboration & Ownership](#collaboration--ownership-2026-05-13) section below for the explicit division of labor and ownership-protection plan.

### Stage 2 note (2026-05-01)

When we open Stage 2 (search), use William Fong's cuvarbase branch + persistent-context M-per-GPU queue architecture as the starting point — he validated S100 cam4/ccd4 equivalence within float roundoff. His numba-vectorized trapezoid fit is also a relevant pattern for any custom CPU-bound search code we write.

## Three-Week Pre-Talk Vertical Slice (2026-05-12 → 2026-06-02)

> **Revised 2026-05-13** to reflect the Schwamb-group collaboration kickoff (see [§Collaboration & Ownership](#collaboration--ownership-2026-05-13) below and progress log [§2.5](twirl_progress_log.md)). Original sprint plan focused on us standing up Stage 2 search + vetting alone; with Michelle's group now bringing a tuned BLS and LEO-Vetter expertise, the sprint reweights toward our actual contribution — light curve quality, the v2-QLP-vs-v3-TWIRL HLSP comparison, and WD 1856 blind recovery as a joint gating milestone.

Time-boxed end-to-end run through Stages 2→5 on the **S56 footprint only**, gated by **WD 1856 blind recovery**. Purpose: surface unknowns now and produce a credible talk for Andrew Vanderburg and Kevin Burdge on ~`2026-06-02`. Stage 1 mass production continues in parallel on its own pdogpu1+pdogpu6 schedule; the sprint does not pull resources from it.

**Why S56 only.** It is the one sector that is fully done, validated, contains the WD 1856 ground truth, and is cadence-aligned. Anything else dilutes the experiment in the time available.

**Non-goals for the sprint** (explicit so they don't creep in):

- Multi-sector phase folding — S56 only.
- Final occurrence-rate posterior — sample too small; any rate quoted is provisional.
- ML triage — heuristic vetter + LEO-Vetter only.
- External follow-up coordination.
- S57+ search runs — even if S57/S58 production finishes in time, freeze the talk results on S56.
- **New non-goal (2026-05-13):** standing up our own dip-search and variable-depth detectors. Deferred until post-talk; Michelle's BLS covers the periodic channel and the WD-1145-style dip-class slide can be motivated with the existing literature.

### Week 1 (2026-05-12 → 2026-05-19): v3 TWIRL HLSPs as the shared survey input

**Done (2026-05-13)**: TGLC negative-flux preservation patch shipped on the [twirl/preserve-negative-flux](https://github.com/TeHanHunter/TESS_Gaia_Light_Curve/tree/twirl/preserve-negative-flux) branch; `flux_space_detrend` UnboundLocalError fix landed (commit `f0e844c`); v3 TWIRL HLSP tree built end-to-end on S56 (`19,068 / 19,072 = 99.98%`); WD 1856 sanity confirms 155 negative-flux cadences preserved and +6.3% usable q=0 cadences with the transit shape intact.

**Still to ship by 2026-05-19:**

- **Hand the v3 TWIRL HLSP tree to Michelle's group as the production search input** (`/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl/`). Schema is QLP-compatible (`hlsp_twirl_tess_ffi_*`), so their reader works without changes.
- **Run our own BLS v3 + heuristic vetter + LEO-Vetter on the v3 tree** as a comparison stack — quantifies what changes vs the v2 QLP run (rank movement for WD 1856, MES recovery for faint hosts, alias-wall behaviour).
- **Build the headline v2-vs-v3 comparison table for the talk**: WD 1856 + ~5 representative faint and bright hosts, side-by-side rank / MES / cadence-count change.

**Week-1 gate**: v3 HLSPs handed off, our own search stack runs end-to-end on v3, comparison table populated.

### Week 2 (2026-05-19 → 2026-05-26): Injection-recovery (shared exploration) + LEO-Vetter tuning experiment

Two parallel work streams; we own the LC-side comparison angle, Schwamb group owns search-side variations. Both are exploratory — multiple approaches in parallel rather than one canonical pipeline. This is now the area where the most methodological learning happens.

- **Our angle — LC-side completeness comparison**: same injected signal grid (`~1,000–2,000` WD-1856-like + `~300–500` WD-1145-like; ranges as in the original sprint) recovered through (a) v2 QLP HLSPs and (b) v3 TWIRL HLSPs with the same downstream search + vetter. The recovered-fraction delta isolates the contribution of the negative-flux preservation as a function of Tmag, depth, and period. This is the cleanest possible methods result for the methods paper.
- **Pixel-level vs LC-level framing for the talk**: LC-level injections are the near-term workhorse because they are cheap enough for dense grids and directly test detrending + search + vetting on the current HLSP products. Pixel-level injections are the paper-grade calibration layer: inject synthetic events into cutout/FFI pixels before ePSF fitting and LC extraction so extraction, deblending, aperture choice, centroid shifts, and crowding failures are included. The talk should present LC-level recovery as the first completeness preview, not the final survey completeness claim, and explicitly name pixel-level injections as the next rigor step.
- **Schwamb-group angle** (parallel): signal-template variations, recovery-statistic alternatives, BLS-vs-other-detector comparisons. Coordination via shared parquet schemas; no need to converge on a single approach.
- **LEO-Vetter tuning experiment**: our [wd-host-tuning](https://github.com/TeHanHunter/LEO-Vetter/tree/wd-host-tuning) fork vs Michelle-student + Franklin's independent tuning. Compare on the same top-N S56 v3 candidate list. Whichever performs better on WD 1856 + the FP-class breakdown wins the talk slide.

**Week-2 gate**: completeness deltas characterized (qualitatively at minimum) for v2 vs v3, and at least one LEO-Vetter tuning passes WD 1856 as PC.

### Week 3 (2026-05-26 → 2026-06-02): WD 1856 joint blind recovery + talk anchor

- **WD 1856 blind-recovery confirmation across the search/vetter combinations** — the gating milestone for the talk. Concretely: WD 1856 (TIC `267574918`) must end up classified as PC by *at least one* (search × vetter) combination on the v3 HLSPs after blind ranking + vetting. The current v2-QLP heuristic-vetter result (rank 9 / 5,403 after vetting) is already proof-of-concept; the v3 run quantifies improvement. If the joint recovery fails, the failure mode becomes the most interesting slide (why → what the search/vetter needs).
- **FP-class characterization** for the Burdge slide — joint with Franklin + Michelle's student: classify the top vetted candidates into PCEBs (WD + M-dwarf), centroid-shift contaminators, systematics ladder peaks, ZZ Ceti pulsators. Count each class. Multi-tuner agreement on these classes is a stronger statement than any single-tuner result.
- **Talk slide outline** (revised for the collaboration framing):
    1. Motivation: WD 1856 + WD 1145, why 200 s FFI matters
    2. The Schwamb-group + Te Han collaboration model (one slide naming who-does-what — keeps the contributions visually distinct)
    3. Light curves: TGLC negative-flux preservation result (Te Han methods slide — v2 vs v3, 155 cadences recovered on WD 1856, +6.3% usable cadences, transit preserved)
    4. Search: BLS rank improvements on v3 (joint with Schwamb-BLS comparison)
    5. Vetting: LEO-Vetter tuning experiment results (joint, two independent tunings)
    6. WD 1856 blind recovery (joint headline)
    7. Injection-recovery preview (LC-level now; pixel-level calibration next)
    8. FP classes
    9. Roadmap: 40-sector campaign, occurrence-rate paper, catalog paper

**Week-3 gate**: a self-contained narrative deck that names the collaborators clearly, attributes contributions explicitly, and holds up under questions from Vanderburg (completeness rigor) and Burdge (FP discrimination, WD-specific eclipsing systems).

### Sprint deliverables

- v3 TWIRL HLSP tree at `/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl/` (**shipped 2026-05-13**).
- v2-vs-v3 HLSP comparison table on a representative TIC sample, with WD 1856 as the headline.
- Side-by-side LEO-Vetter tuning comparison: TWIRL fork vs Michelle-student + Franklin's tuning.
- Rough completeness deltas (v2 vs v3) for the WD-1856-like and WD-1145-like signal regimes.
- Documented WD 1856 blind recovery across (search × vetter) combinations.
- FP-class breakdown of top vetted candidates.
- Talk slide deck with explicit collaborator attribution.

## Collaboration & Ownership (2026-05-13)

This section captures the division of labor and the ownership-protection steps agreed at the Schwamb-group kickoff meeting on `2026-05-13`. Treat as live — update at every meeting with Michelle. Companion meeting note in progress log [§2.5](twirl_progress_log.md).

### Who owns what

| Layer | Owner | Notes |
|---|---|---|
| **TGLC light-curve production** (v3 flux-space-detrended HLSPs at scale) | Te Han (this repo) | Unambiguously ours. Includes the [twirl/preserve-negative-flux](https://github.com/TeHanHunter/TESS_Gaia_Light_Curve/tree/twirl/preserve-negative-flux) TGLC patch and the [flux_space_detrend](../src/twirl/lightcurves/flux_detrend.py) cotrend. Foundation of every downstream result. |
| **Data steward for the survey LCs** | Te Han | Versioned tree, sidecar manifests, public release decisions. |
| **BLS (production search)** | Michelle's group | Their well-tuned implementation is the primary tool. Our [src/twirl/search/bls.py](../src/twirl/search/bls.py) becomes a comparison reference — useful as a sanity-check and for the v2-vs-v3 talk slide, not the headline. |
| **LEO-Vetter tuning for WDs** | Joint experiment: Michelle's student + Franklin Chen one tuning; our [wd-host-tuning](https://github.com/TeHanHunter/LEO-Vetter/tree/wd-host-tuning) fork the other | Two independent tunings produce more learning than one. Talk presents both side-by-side. |
| **Heuristic vetter** | Te Han (this repo) | [src/twirl/vetting/heuristic.py](../src/twirl/vetting/heuristic.py) — physics-motivated, non-rejecting, classifies into `vet_class`. Complementary to LEO. |
| **Injection-recovery** | **Shared exploratory work** — multiple approaches in parallel | LC-level grids provide the fast v2-QLP-vs-v3-TWIRL completeness comparison. Pixel-level injections are the slower calibration layer for extraction, deblending, aperture, centroid, and crowding effects. Their angle: signal-template + statistic variations. |
| **Centroid / pixel-map on-target test** | Te Han (TODO in `ideas.md`) | Per-pixel residual LC cube from `effective_psf.fit_lc` is in the TGLC mitpdo branch already; the build-out is on our side. |
| **Occurrence-rate paper** | **Te Han, lead author** (verbal commitment, to be confirmed by email this week) | Highest-priority ownership prize from this collaboration. |
| **Catalog paper** | **Target: led by a TWIRL student, Te Han as senior/last author.** Michelle's group accommodated with a parallel student-led paper instead. Not yet proposed. | Mentoring credit is the explicit goal — search committees weight student-led papers heavily for faculty applications. Catalog is the natural mentoring vehicle since it sits on the data we produce. Specific TWIRL student to be named before the next meeting (candidate seniority needs to match a catalog-paper workload). |
| **WD-tuned LEO-Vetter methodology paper** | **Proposed: led by Michelle's student** (the parallel-tuning experiment is its natural first-author result) | Offer at the next meeting as a swap for the catalog paper. Frames Michelle's student's tuning work as its own clean first-author deliverable rather than a sub-result inside the catalog paper. Both groups end up with student-led papers. |
| **Methods / pipeline paper** | **Te Han, lead author** (Te Han to propose) | Separate first-author deliverable covering the TGLC negative-flux patch + flux-space cotrend + per-orbit detrend. Decoupling from the occurrence-rate paper protects authorship on at least two papers regardless of how the catalog paper falls. |
| **TGLC fork PR upstream to MIT-Kavli** | Te Han | Preserves Te Han as the patch author of record before any wider adoption. |
| **Discovery paper(s)** — if TWIRL finds new planets / debris systems | **Currently undefined; default convention favors whoever shepherds validation, NOT the survey lead.** Highest-priority item to pre-commit. | See [§Discovery-paper authorship](#discovery-paper-authorship-2026-05-13) below. The first TWIRL discovery should be Te Han lead or TWIRL-student lead (Te Han senior); subsequent discoveries follow a "candidate champion" rule with Te Han as senior author on all TWIRL-originated discoveries. Requires explicit pre-commit with Michelle, Vanderburg, and Burdge *before* any candidate emerges. |

### Strategic trade-off (named explicitly)

The pivot from solo project to group effort is a positive-net trade with a known downside that we should not suppress:

- **Wins:** more productive, broader expertise, higher credibility for the talk and papers, faster convergence on the right BLS + vetter recipe, multiple LEO tunings to compare, shared injection-recovery exploration.
- **Cost:** less individual ownership on the search/vetting layers. Those layers are now group products.

The mitigations in the ownership table above are the concrete response. The v3 HLSP production work and the negative-flux preservation patch remain unambiguously ours; the methods paper makes that contribution a distinct first-author deliverable.

### Items to lock in writing (this week)

1. **Email Michelle** explicitly confirming Te Han leads the occurrence-rate paper — turn the verbal commitment into a written one.
2. **Propose the methods/pipeline paper** in the same email, with a one-paragraph scope (TGLC negative-flux patch + flux-space cotrend + per-orbit detrend + the recovery-fraction-by-Tmag result that quantifies the contribution).
3. **Acknowledge funding** — confirm TGLC compute (PDO) goes through Te Han's affiliation in any joint paper acks.

### Discovery-paper authorship (2026-05-13)

**The most under-defined and most consequential ownership gap in the collaboration.** Currently nothing has been agreed about who leads if TWIRL finds a new WD planet or a new debris-transit system. Defaulting to "whoever shepherds validation leads" — the field convention — is structurally bad for Te Han because validation expertise lives with Vanderburg and Burdge (telescope time, RV mass measurements, spectroscopic FP rejection), not with the survey lead.

**Target structure:**

- **First TWIRL discovery**: Te Han lead, OR a TWIRL student lead with Te Han as senior author (decide at the time based on student readiness). This anchors Te Han as "the person who discovered WD planets with TWIRL" — career-defining for the faculty market.
- **Subsequent discoveries**: "candidate champion" rule — whoever drives validation leads — but Te Han is **always senior author** on TWIRL-originated discoveries. Lets students lead (mentoring credit) and lets collaborators lead (relationship credit) without giving away the survey-author thread.
- **Joint discoveries** (e.g., simultaneous TWIRL + ground-based detection): co-first-authorship with explicit TWIRL contribution statement.

**Why this has to be pre-committed *before* any candidate emerges:**

Once a real candidate is on the table, the political pressure to assign discovery authorship to whoever has telescope time is enormous, and the negotiation becomes retroactive — which always favors the seniority + capability side over the survey side. Raising the question proactively, with no candidate at stake, is the only way to get a fair outcome. Doing it after a candidate emerges marks Te Han as someone trying to grab credit retroactively, which is much worse than not asking at all.

**Critical infrastructure investment: Te Han's independent validation capability.**

Defaulting to Vanderburg or Burdge for follow-up makes the "candidate champion" rule structurally cede discovery authorship to them. Building independent validation capacity — telescope time for follow-up (Magellan, Keck, Gemini proposals), spectroscopy reduction pipeline, RV mass measurement workflow — is the **most important infrastructure investment for the next 12-24 months**, more important than any pipeline code, because it converts the field convention from "works against the survey lead" to "works for the survey lead." Without it, no amount of authorship agreement on paper survives contact with the realities of validation logistics.

**The validation landscape — what Vanderburg actually did for WD 1856, and what Te Han can do:**

**Reality check (corrected 2026-05-13):** the WD 1856 b discovery paper (Vanderburg et al. 2020, *Nature*) was validated using **zero PI-class instruments**. The toolkit was TESS (public), Spitzer IRAC (now JWST, open TAC), Gemini North (facility, open TAC), Hobby-Eberly Telescope (facility), and Gran Telescopio Canarias (facility). Vanderburg first-authored because he recognized the candidate and wrote the proposals — not because he had access to anything we don't.

The Vanderburg-playbook validation chain for a TWIRL discovery — every step facility-accessible:

| Step | Tool | TAC route |
|------|------|-----------|
| Multi-sector photometric confirmation (when applicable) | TESS itself | Public |
| Higher-cadence ground photometric confirmation | **LCO 1m/2m network** (1-2 min cadence) | NSF-funded, free queue access for US researchers — **read the docs this week** |
| Brighter / faster ground photometry | Gemini imaging | Facility, MIT TAC |
| Thermal emission constraint (the planet-vs-BD test) | **JWST** (replaces Spitzer) | Open NASA TAC |
| Spectroscopic companion exclusion | Gemini, Keck, MIKE on Magellan, HET | All facility |
| High-res imaging (background-EB exclusion) | Gemini speckle, Keck AO | Facility |
| RV mass (rare; only if model demands it) | MIKE (Magellan facility, no PI gate) — or ESPRESSO/HARPS via ESO | Mostly facility |
| **Optional faint-end enhancement (T>18 only)** | proto-Lightspeed on Magellan Clay | PI-class — Burdge collaboration only required if read noise matters, i.e. for the very faint regime |

**proto-Lightspeed** ([2026B Magellan list](https://lweb.cfa.harvard.edu/files/tac/2026B_MagInstr.txt), [Burdge group page](https://binaries.mit.edu/proto-lightspeed/)): a high-speed qCMOS imager on Magellan Clay Nasmyth East with 0.29 e⁻ read noise in quiet mode, up to 6600 Hz windowed cadence. PI-class instrument under Kevin Burdge — the official 2026B note is "may be opportunity to collaborate." It's the best-in-class faint-end tool. It is **not** required for validating bright WD candidates (T < ~17), where facility imaging is sufficient. WD 1856 itself (T=16.34) was validated entirely without it.

**Implications for Te Han's first-authorship path:**

1. **For most of the TWIRL Tmag distribution (T < 17), the Vanderburg playbook works directly.** Te Han can first-author discoveries using TESS + LCO + Gemini + JWST + MIKE — all open TAC, no PI gates. Identical structurally to how Vanderburg first-authored WD 1856 b.
2. **proto-Lightspeed is a faint-end add-on, not a gate.** For T > 18 candidates where read noise matters, Burdge collaboration is required and discovery authorship is co-/shared. But that's the minority regime.
3. **LCO is the single most actionable item this week.** Free queue-scheduled time on a global 1m/2m network for US researchers; it's what carries the ground-based photometric confirmation in nearly every TESS planet discovery paper. Te Han should read LCO docs and identify the right observing mode (typically the 1m network with a SDSS-r filter, ~30-60 s exposures).
4. **JWST proposal cycle planning** should start within 6 months of the first TWIRL candidate emerging. The thermal-emission test (planet vs BD discrimination) is the headline validation result for any WD 1856-class candidate. JWST cycle proposals are competitive (~20% acceptance) but no PI gate.
5. **Build the validation pipeline NOW**, before candidates emerge. LCO photometric reduction, MIKE spectroscopic reduction, basic high-res imaging analysis. Each of these is a 1-2 month project. Doing them under no time pressure now is much easier than doing them under candidate-deadline pressure later.

**Burdge collaboration is still strategically valuable, just not mandatory.** A proto-Lightspeed sub-electron-read-noise light curve adds confidence beyond what facility instruments give, especially for the faint regime. Co-PI with Kevin on a Lightspeed proposal *and* independent Te-Han-led LCO+JWST validation isn't an either/or — the strongest discovery papers will use both, with Te Han as lead on the assembled story.

**Critical timing caveat:** Magellan + LCO + Gemini all run ~6-9 month proposal lead times. First TWIRL candidates likely emerge before Te Han's first allocated time. Plan for this: build the *infrastructure* (data reduction pipelines, proposal templates) before candidates emerge; the first candidate's validation may use Director's Discretionary Time (DDT) on facility instruments (faster turnaround), or piggyback on a collaborator's existing time. The *first-discovery authorship rule* must commit Te Han (or a TWIRL student) as lead/co-first/senior on the resulting paper regardless of whose telescope time was used.

**Conversations to have:**

1. **Schwamb-group next meeting**: lock the "first TWIRL discovery → Te Han or TWIRL student lead" rule with Michelle. Frame as "let's settle this now while there are no candidates" — collaborative, not adversarial.
2. **At the 2026-06-02 talk**: raise discovery-authorship with Vanderburg and Burdge. Not "promise me lead authorship" but "what does discovery authorship look like for TWIRL candidates? — want to understand the collaboration model before candidates emerge." Both have been through this many times; they'll respect the early framing.
3. **Internal**: identify TWIRL student(s) who could plausibly lead a discovery paper. Franklin Chen is the obvious candidate but probably too early-PhD to drive validation. The student who leads the catalog paper is the natural candidate for second-discovery lead if and when that happens.

### Catalog-paper negotiation — bring to next meeting

**Goal:** TWIRL student leads catalog paper, Te Han as senior/last author — mentoring credit is explicitly important for faculty-market positioning. Michelle's group offered a parallel student-led paper (WD-tuned LEO-Vetter methodology) so both groups end up with student-led first-author papers.

**Specific TWIRL student to lead:** to be named before the meeting. Candidate seniority needs to match a catalog-paper workload (sample definition, candidate vetting, all the operational pain). Franklin Chen is the obvious candidate given his LEO-tuning involvement, but seniority and readiness need a separate assessment — naming the wrong student is worse than leaving the slot generic.

**Recommended framing for the conversation** (don't open with this — wait until paper authorship comes up naturally):

> "I'd like a TWIRL student to lead the catalog paper. I need to demonstrate mentoring for the faculty market, and catalog is the natural fit since it sits on the data we're producing. For your student — they're doing interesting work on LEO-Vetter WD tuning, and that's its own first-author methods paper: thresholds, override rationale, validation on v3 HLSPs, fork-vs-fork comparison. That gives them a clean first-author paper distinctly theirs, and we end up with a clean portfolio: occurrence rate (Te), catalog (TWIRL student), LEO-WD-tuning (your student), and pipeline methods (Te). Four clean first-author papers, no shared leads, both groups mentoring."

**Risk flags:**

1. **Michelle counters with co-lead on catalog.** Hold the line — shared catalog leads tend to go badly. Better outcome: their student leads LEO paper clean, our student leads catalog clean.
2. **The named TWIRL student isn't ready.** Catalog papers are operationally hard. Real assessment of Franklin (or whichever student) before the meeting.
3. **Order of publication matters.** Set explicit target timelines for both student-led papers so neither student feels deprioritized.
4. **Methods paper out first or in parallel** is doubly important now — fences the pipeline contribution off so it isn't absorbed into the catalog paper's data section.

### Items to bring to the next meeting

- v3 vs v2 HLSP comparison results (WD 1856 + a faint-end sample).
- LEO-Vetter tuning comparison once both tunings have run on the same v3 top-N.
- The fact that the existing TGLC HDF5 has a pre-built per-pixel residual LC cube (`effective_psf.fit_lc` in the mitpdo branch) that the per-pixel pixel-map vetting tool should read — this is a Te Han-side build but uses an existing TGLC artifact.
- Authorship list draft for both the occurrence-rate paper and the methods paper.

## Project Goal

TWIRL will use 200 s TESS FFIs, TGLC-based light-curve extraction, transparent transit/dip searches, optional machine-learning triage, automated vetting, and injection-recovery tests to:

1. search for any transiting or occulting objects around white dwarfs,
2. measure the completeness of that search in the regimes where 200 s data are strongest,
3. derive occurrence-rate constraints or upper limits for those regimes, and
4. produce validated candidate lists and follow-up targets.

The actual execution sequence:

1. generate WD light curves with the MIT-adapted TGLC on MIT PDO machines (and via QLP for newer sectors — see Stage 1.8),
2. build a transparent baseline search for periodic and non-periodic transit-like events,
3. run injection-recovery tests to establish completeness,
4. search the full WD sample,
5. validate resulting candidates,
6. add ML triage only if it materially improves the baseline search ranking or false-positive rate.

## Design Drivers From The Proposal

- Primary input sample: the local external Gentile Fusillo et al. (2021) Gaia EDR3 white dwarf catalogue, recommended at `data_local/catalogs/GaiaEDR3_WD_main.fits`, with Gaia DR3 as the authoritative target identifier and TIC/TESS metadata added when available.
- Photometry source: 200 s TESS FFIs only, which means sectors/orbits from Sector 56 onward.
- Extraction engine: the MIT-adapted TGLC pipeline already used in the QLP environment.
- Search target: any transiting or occulting object around a WD, with first-year emphasis on large, deep, short-duration events in the WD 1856-like regime. Secondary goal (added 2026-05-01): WD 1145+017-style **disintegrating planet candidates** — variable-depth, asymmetric (dust-tail), possibly aperiodic dips. WDs are physically the right host class for these (high surface gravity ⇒ rapid stripping at close orbits; ~3% of DA WDs show circumstellar dust IR excess), and a `Pwd > 0.75` survey of `~360k` WDs is the right population. The search stack therefore needs, in addition to BLS, a depth-allowed-to-vary template detector that can fire on aperiodic dust-tail signatures (cf. BD+05 4868 around metal-rich subgiants); see Stage 2 design notes.
- Science outputs: validated candidates, a ranked follow-up list, and occurrence-rate limits first for large occulters; Earth-size/HZ occurrence work is a longer-term goal only after completeness is demonstrated.
- Companion-vs-host size regime (added `2026-05-07`): for the WD survey we should expect `R_companion >= R_host` to be the *common* case, not the exception. White dwarfs sit at `R_WD ~ 0.013 R_sun ~ 1 R_jup`, so any planet larger than Mars-sized has `R_p > R_WD`, and post-CE M-dwarf or brown-dwarf companions have `R_companion ~ 1-5 R_jup`. This flips the standard main-sequence intuition: transit duration is set by the *larger* body and scales with `R_companion + R_WD`, not just `R_host`. Practical consequences: (1) transit depths can approach `100%` (companion fully occults the WD), as in WD 1856; (2) duration grids and template models for the search stack must allow up to ~hour-scale events for WD+MD post-common-envelope binaries, not just the minute-scale WD-only regime.

## Recommended Repo Layout

The repo should be organized around the pipeline stages rather than papers or notebooks.

```text
doc/
  twirl_plan.md
  twirl_progress_log.md
  mit_tglc_usage_guide.md
  ideas.md
catalogs/
  wd_master_catalog/
  sector_orbit_maps/
configs/
  tglc/
  detection/
  injections/
scripts/
  stage1_lightcurves/   # TGLC pipeline + HLSP build (live)
  stage2_search/        # per-sector BLS + candidate consolidation (live)
  stage3_injections/    # planned, not built yet
  stage4_search/        # planned, not built yet (full-sample search; will reuse stage2_search code with different drivers)
  stage5_validation/    # heuristic vetter + LEO adapter + centroid test (live)
src/
  twirl/
    catalogs/
    io/
    lightcurves/
    plotting/
    search/
    vetting/
    # detection/, injections/ — planned, not built yet
notebooks/
reports/
```

## External Data Policy

Large external inputs and staged survey products should remain local and outside normal git tracking.

Recommended local layout:

```text
data_local/
  catalogs/
    GaiaEDR3_WD_main.fits
    twirl_master_catalog/
  tglc-data/
```

Rules:

- do not commit large raw FITS catalogs to this repo
- do not treat local external data as part of the git-managed source tree
- record local input path, file size, and provenance in derived metadata products
- version-control the code, configs, schemas, and derived table definitions that act on those local files

## Stage 1: Generate WD Light Curves On MIT PDO Machines

This is the foundation for everything else. The MIT fork is no longer a single-target `quick_lc.py` workflow. It is an orbit/camera/CCD production pipeline with separate catalog, cutout, ePSF, and light-curve stages.

### 1.1 Build the WD target table (✓)

The current seed catalogue for this repo is:

- `data_local/catalogs/GaiaEDR3_WD_main.fits` or another user-configured local path
- based on Gentile Fusillo et al. (2021), "A catalogue of white dwarfs in Gaia EDR3"
- this is the main Gaia EDR3 catalogue, not the reduced-proper-motion extension
- this file is a local external dependency and should not be committed to git

Important catalogue facts to treat as part of the survey definition:

- the main catalogue contains about `1.28 million` sources that passed the authors' quality selection
- the paper identifies `Pwd > 0.75` as the high-confidence WD threshold, yielding about `359,073` candidates
- the paper quotes an upper-limit completeness of about `93%` for white dwarfs with `G <= 20`, `Teff > 7000 K`, and `|b| > 20 deg`

TWIRL should build a single master catalog from this file with, at minimum:

- Gaia DR3 designation/source ID
- TIC ID, if available
- RA, Dec
- Gaia `G`, `BP`, `RP`
- Gaia `Pwd`
- TESS magnitude estimate
- WD selection probability / class flag from the input WD catalog
- sector/orbit coverage in TESS
- camera, CCD, and cutout mapping once known

This table becomes the control plane for the whole survey.

The master catalog should also include provenance fields that make later survey definitions reproducible:

- local seed-catalog path used at build time
- seed-catalog file size and modification time, or a stronger file hash when practical
- TWIRL catalog build version
- crossmatch recipe version
- sample-selection recipe version

The parent-sample protocol should be explicit in code and metadata:

- `seed_catalog`: all rows in the local Gentile Fusillo et al. (2021) main catalogue
- `reference_highconf_sample`: rows with `Pwd > 0.75`, used as a default high-confidence comparison sample for QA and pilot work
- `lc_gather_sample`: all WDs in the full parent catalog — light curves are extracted for all of them; no sample cut is applied at this stage
- `twirl_parent_sample`: the final locked denominator for occurrence-rate work, to be defined after light-curve collection and QA are complete; this cut must be frozen before any injection-recovery or occurrence-rate analysis begins
- until `twirl_parent_sample` is frozen, any completeness or occurrence summaries must be labeled provisional
- any exploratory search outside the final locked denominator must be labeled exploratory and excluded from occurrence-rate denominators by default

#### Status

- First Gaia-first WD master catalog and conservative TIC merge path are in place.
- Current milestone facts: `1,280,266` seed rows, `359,073` rows with `Pwd > 0.75`, `1,236,467` unique TIC matches, `43,035` ambiguous TIC matches, and `764` no-bridge cases.
- Detailed dated notes live in the Stage 1.1 [log](twirl_progress_log.md#11-build-the-wd-target-table).

### 1.2 Decide the production sample (✓)

Use two nested samples:

- `search_sample_primary`: the highest-priority WD sample for the first end-to-end run
- `search_sample_full`: the broader sample used for the final occurrence-rate analysis

A practical starting point is:

- use `Pwd > 0.75` as the default high-confidence reference sample for pilot QA and crossmatch studies
- primary sample near `T <= 18`
- full sample extended toward `T <= 20` where recovery remains useful
- include WD 1856+534 as the mandatory smoke-test target
- choose one `Sector >= 56` containing WD 1856+534 as the first end-to-end benchmark; the exact sector can be selected later

This creates two linked but distinct products:

- a statistical survey sample with a locked denominator and documented cuts
- an exploratory event-search layer that can be broader, but does not feed occurrence-rate claims by default

The exact TWIRL WD cut should be frozen only after:

- the Gaia-first target definition and any required MIT TGLC implementation changes are characterized well enough to support production
- benchmark QA is complete
- the baseline search stack is defined well enough to support end-to-end injections

#### Status

- The Stage 1 benchmark and pilot production-sample policy are set: use the Gaia-first parent catalog, keep `Pwd > 0.75` as the default high-confidence reference sample, and preserve unresolved final denominator choices for later occurrence-rate stages.
- Current milestone facts: `1,231,702 / 1,280,266` WDs (`96.2%`) have at least one `Sector >= 56` footprint, and `1,157,242` (`90.4%`) have observed coverage before Sector `100`.
- Detailed dated notes live in the Stage 1.2 [log](twirl_progress_log.md#12-decide-the-production-sample).

### 1.3 Prepare the MIT PDO environment (✓)

The MIT-adapted TGLC expects:

- local TICA FFI files already staged on disk,
- local `pyticdb` access to `tic_82` and `gaia3`,
- a `tglc-data` directory tree organized by orbit/camera/CCD.

For TWIRL, the immediate operational tasks are:

1. freeze the orbit list to process,
2. identify which orbit/camera/CCD combinations contain WD targets,
3. pre-stage the required 200 s FFIs,
4. confirm the `pyticdb` databases resolve on PDO,
5. define batch-job wrappers for one orbit-camera-CCD unit per job.

#### Status

- PDO helper-script environment and `pyticdb` access to `tic_82` are validated for the Gaia-to-TIC export path.
- Detailed dated notes live in the Stage 1.3 [log](twirl_progress_log.md#13-prepare-the-mit-pdo-environment).

### 1.4 Generate catalogs and cutouts (...)

The MIT fork runs in four logical stages:

1. `catalogs`
2. `cutouts`
3. `epsfs`
4. `lightcurves`

For TWIRL, the key configuration change is to lift the TIC magnitude limit well beyond the default value. The MIT fork defaults to bright-star QLP use cases; WD work needs a much deeper target cut.

Initial production assumptions:

- run `--max-magnitude 20`
- keep Gaia catalogs for all relevant field stars
- run per orbit and per CCD, not per target

#### Status

- The fixed Sector `56`, orbit `119/120`, `cam4/ccd1` benchmark has been pushed through `catalogs`, `cutouts`, `epsfs`, and WD-only `lightcurves`; WD 1856 is present in both orbit trees as `267574918.h5`.
- The pre-Sector-67 TICA header fallback is patched in the user-owned TGLC fork, and the current default per-CCD `epsfs` worker count is `--nprocs 16` based on the first one-CCD benchmark.
- The first full-orbit orbit `119` benchmark on `pdogpu1` has completed through HDF5 `lightcurves`: the 15-CCD sweep outside the original `cam4/ccd1` benchmark finished in `33.66 h` wall time with outer CCD concurrency `3` and inner workers `16/16/16/16`, and the final orbit tree now contains `32` catalog files, `3136` cutouts, `3136` ePSFs, and `19086` WD-only `.h5` light curves.
- The clean stopwatch for current performance planning is therefore `33.66 h` for the 15 non-benchmark CCDs; the orbit is now complete across all 16 CCDs, but the original `cam4/ccd1` benchmark was produced earlier in separate runs with different tuning, so it should be treated as a separate seed benchmark rather than folded into one apples-to-apples full-sector wall-time claim.
- Detailed dated notes live in the Stage 1.4 [log](twirl_progress_log.md#14-generate-catalogs-and-cutouts).

### 1.5 Extract and consolidate WD light curves

The full Stage 1 pipeline per orbit/camera/CCD is:

```
tglc catalogs → tglc cutouts → tglc epsfs → tglc lightcurves
  → qlp lctools detrend → qlp lctools hlsp (→ FITS per TIC)
```

`tglc lightcurves` writes one HDF5 file per TIC target. The `detrend` step (run via `qlp lctools detrend`) annotates those HDF5 files with a `bestdmagkey` attribute that `hlsp.py` reads to select the detrended magnitude column. **Without detrending, HLSP generation fails.** The `hlsp` step then converts HDF5 → FITS, requiring the full QLP environment (`lightcurvedb` + `qlp` package, not just TGLC).

For Sector 56 specifically, HLSP generation uses SPOC quality flags (not TICA), because `get_ticaflags` in `hlsp.py` raises `ValueError` for `sector < 67`. These SPOC flags live at `/pdo/qlp-data/spocflags/`.

**Standard recovery rate check** — run after both HDF5 and FITS production:

For each orbit/camera/CCD:
1. Count TIC IDs requested (from TWIRL detector target table)
2. Count HDF5 files produced (`tglc lightcurves` output)
3. Count FITS files produced (`qlp lctools hlsp` output; missing TICs logged as warnings by `generate_qlp_hlsp_fits_file`)
4. Bin the shortfall by TESS magnitude to detect faint-end dropout
5. Cross-check against Gaia-only targets (no TIC bridge): these will always be 0% recovered until the MIT fork is extended to emit Gaia-selected targets directly

This 3-level audit (requested → h5 → FITS) is the standard health check for every production run. Log counts and dropout fractions in the progress log.

After production, build a TWIRL light-curve index with:

- TIC ID and Gaia DR3 ID
- orbit, sector, camera/CCD
- HDF5 path and FITS path
- production status flag (`h5_ok`, `fits_ok`)
- summary quality metrics (RMS, MAD, cadence count)

The index format should explicitly distinguish raw catalog columns, derived TWIRL metadata, crossmatch outputs, processing-status fields, and sample-membership flags.

### 1.6 Run photometric QA before any ML work

Before training or search, measure:

- per-target RMS / MAD versus magnitude
- sector-level failure rates
- missing cadence statistics
- light-curve completeness by target and by orbit
- recovery of WD 1856+534 b in the 200 s products
- agreement or disagreement across the 1x1, 3x3, and 5x5 apertures for benchmark targets
- basic comparison against at least one independent extraction path on the benchmark target set

This step should produce a small QA report and a list of sectors/orbits to reprocess if needed.

### 1.7 Benchmark acceptance criteria

Before scaling beyond the pilot stage, the WD 1856 benchmark should satisfy all of the following at the documentation level:

- the expected transit window is visible in at least one `Sector >= 56` benchmark product
- the event timing is consistent with the published WD 1856 ephemeris within the expected propagated uncertainty
- the signal is present in at least one aperture product and checked across the `1x1`, `3x3`, and `5x5` apertures
- any strong aperture disagreement is documented and investigated as a contamination or extraction issue
- an independent extraction path has been compared on the benchmark set
- the benchmark result is strong enough to serve as a smoke test for the first baseline search implementation

### 1.8 Gaps to close in the current MIT fork (?)

The MIT fork is close to what TWIRL needs, but not identical to the survey requirements.

Known gaps to track early:

- Gaia DR3 should remain the survey-defining target identifier even if the current MIT implementation is still partly TIC-oriented.
- Gaia-to-TIC matching should be audited for metadata and implementation reasons, but it should not define the scientific parent sample by itself.
- If important WD targets lack TIC IDs, the light-curve stage will need to be extended to support Gaia-selected targets directly.
- For current terminology, distinguish `TGLC end-to-end` from `QLP/MAST end-to-end`: Stage 1 extraction is complete at the HDF5 `lightcurves` output, while public-style FITS deliverables require the later QLP detrend + HLSP export path.
- After the full Sector 56 benchmark run, verify that all metadata needed to launch whole-sector TIC light-curve extraction with a magnitude cut is already present and discoverable from the generated products, not hidden in ad hoc scripts.
- After the full Sector 56 benchmark run, confirm that the `cutouts` and `epsf` products under `/pdo/users/tehan/tglc-deep-catalogs/orbit-*/ffi/cam*/ccd*/` are exactly in the MKI TGLC `Manifest` layout so they can be reused by standard `tglc lightcurves` commands without special-case path handling.
- The new output product is decontaminated aperture photometry in HDF5, not the old per-target FITS product with `cal_psf_flux` and `cal_aper_flux`.
- The old `prior`-based single-target workflow is not exposed in the new CLI, so any need for floating-field-star priors must be added explicitly.
- The catalogue completeness quoted by Gentile Fusillo et al. (2021) applies only in specific regimes; TWIRL must not generalize that completeness to the full search sample without its own end-to-end validation.

Questions to answer before mass production:

- What outer CCD concurrency should be used per stage once the first full-orbit benchmark finishes, especially for the DB-heavy `catalogs` stage versus the CPU-heavy `epsfs` stage?
- Why do the WD-only `lightcurves` runs currently write fewer `.h5` files than requested TIC IDs, and is that shortfall astrophysical, metadata-driven, or an extraction/control-plane issue?
- What is the official QLP/LCDB metadata path needed to run detrend and HLSP export without the current local benchmark shims?
- Once the full Sector `56` benchmark is done, should whole-sector production stay WD-only at `lightcurves`, or should TWIRL also preserve a full-TIC sector archive for direct QLP/MAST-style downstream products?
- **QLP ingestion for Sectors 94+**: TGLC has been incorporated into QLP (Petitpas et al. 2025, in prep.) and is producing light curves from Sector 94 onward. Should TWIRL ingest those QLP/MAST products directly for recent sectors rather than re-running MIT PDO extraction? Key questions: format compatibility with the HDF5 archive, whether the WD magnitude limit (`--max-magnitude 20`) was applied, and whether the Gaia-first target list is recoverable from QLP outputs. Resolving this could substantially reduce the PDO compute load for newer sectors.

### Stage 1 Deliverables

- WD master target catalog (✓)
  Current status: v0 catalog built and full Gaia DR3 to TIC export completed.
- orbit/camera/CCD job table
- automated PDO production scripts
- consolidated WD light-curve archive
- QA report and reprocessing list

## Stage 2: Build The Search Stack And Optional WD-Specific Classifier

The proposal points toward Astronet-Triage-like deep-learning search, but the first survey version should not be ML-first. TWIRL should begin with an interpretable search stack and add ML only where it clearly improves ranking or false-positive control.

### 2.1 Define the detection unit

Decide what the search stack scores:

- sector-level candidate windows,
- period-search peaks,
- folded light curves,
- local/global transit views derived from a separate search stage,
- and non-periodic dip clusters that may represent debris or irregular occultations.

For the first version, the cleanest path is:

1. produce a transparent periodic candidate list with a short-duration box/trapezoid search,
2. run a separate dip-search branch for non-periodic or weakly periodic events,
3. pass candidates into automated vetting with LEO-vetter (Kunimoto et al. 2025) — assess input format compatibility with the TWIRL HDF5 archive before committing to it,
4. add a WD-specific classifier later if it materially improves triage.

### 2.1b Minimum baseline search specification

Before adding ML triage, the baseline search should define and document:

- the detrending approach used for minute-scale events, with explicit injection checks that the detrending does not erase the target signals
- how `1x1`, `3x3`, and `5x5` apertures are searched and compared
- the rule for choosing a preferred aperture per candidate
- the candidate statistics recorded by the periodic branch
- the candidate statistics recorded by the non-periodic dip branch
- the multi-sector merge logic for periodic candidates
- the rule for keeping irregular or apparently aperiodic events separate from the periodic survey sample

### 2.2 Build the labeled dataset

Positive examples:

- injected WD transits in real TGLC light curves
- WD 1856+534 b as a benchmark
- synthetic variants covering duration, depth, cadence loss, and crowding
- large-occulting benchmark cases spanning giant planets, brown dwarfs, and compact stellar companions

Negative examples:

- quiet WDs
- eclipsing binaries and blends
- scattered-light / systematics artifacts
- cadence-gap and momentum-dump failures

### 2.3 Train the model

If a classifier is added, the training pipeline should include:

- train/validation/test splits by target, not random cadence
- class-balance control
- calibration of classifier scores
- explicit evaluation versus magnitude, period, depth, and number of sectors

### Stage 2 Deliverables

- reproducible periodic-search baseline
- reproducible dip-search baseline
- reproducible training set definition
- trained detector checkpoint(s), if justified
- evaluation report
- inference script that scores the full TWIRL archive

## Stage 3: Injection-Recovery Tests

This is the completeness backbone of the survey and should run through the same search stack used in Stage 4.

Stage 3 should explicitly keep two injection levels separate:

- **LC-level injections**: inject signal models into existing raw or detrended light-curve products before the relevant downstream step. This is the correct fast path for the S56 talk, v2-vs-TWIRL-FS comparisons, detrend signal-preservation tests, and dense grids over period, duration, depth, Tmag, and sector count. It does **not** measure failures in extraction, deblending, aperture selection, or pixel-level contamination.
- **Pixel-level injections**: inject synthetic occultation signals into the cutout/FFI pixel data before ePSF fitting and aperture extraction, then run the real TGLC/TWIRL light-curve builder and search stack. This is slower and should be run on a calibration subset, but it is the completeness layer needed for publication-grade occurrence rates because it includes extraction systematics, crowding, aperture disagreements, and centroid/on-target behavior.
- **Adopt a two-stage protocol**: use LC-level injections for iteration and first-order completeness surfaces; use pixel-level injections to calibrate where the LC-level approximation breaks, especially near the faint limit, in crowded fields, and for candidates where centroid or aperture behavior drives classification.

### 3.1 Injection design

Inject over a grid in:

- orbital period
- transit depth
- duration / impact parameter
- host magnitude
- number of observed sectors

The grid should be densest near:

- the WD 1856-like large-occulting regime,
- Roche-limit boundary cases,
- deep short-duration transits that are most astrophysically plausible for WD systems,
- and only secondarily the Earth-size/HZ regime, which should not drive the first-year completeness claims.

### 3.2 Recovery protocol

Each injection should pass through:

1. the search stage,
2. the classifier, if used,
3. the vetter,
4. candidate merging rules.

Do not estimate completeness from only one stage; use the true end-to-end recovery fraction.

### 3.3 Completeness products

Produce completeness surfaces as a function of:

- `Tmag`
- period
- depth
- number of sectors
- crowding / contamination metrics

Completeness products should be tagged with the exact sample definition and search-branch definition they apply to. Do not present a completeness surface as survey-wide if it only applies to a provisional cut or one search mode.

### Stage 3 Deliverables

- injection engine
- recovery summaries
- completeness grids for occurrence-rate work

## Stage 4: Search For WD Transiting Objects

After the light curves, detector, and completeness machinery are stable, run the production search.

### 4.1 Full survey inference

Run the periodic and non-periodic search branches across all WD light curves and record:

- candidate ephemerides
- detection scores
- per-sector evidence
- merged multi-sector candidates
- flags for apparently aperiodic or irregular events

### 4.2 Candidate catalog construction

For each candidate, store:

- target identifiers
- discovery sectors
- transit parameters
- detection statistics
- vetter outputs
- links to diagnostic plots and files

### Stage 4 Deliverables

- machine-generated candidate table
- diagnostic plot package
- ranked follow-up target list

## Stage 5: Validate The Targets

Validation should combine automated filtering and external follow-up.

### 5.1 Automated checks

- contamination and crowding review
- centroid and aperture behavior
- odd/even and depth consistency checks where applicable
- ephemeris consistency across sectors
- image-level inspection around predicted transits

### 5.2 Archival time-domain first pass (ZTF + ATLAS)

Before requesting any new telescope time, cross-match surviving candidates against the ZTF and ATLAS public archives. Both surveys provide multi-year optical light curves for most TWIRL targets with `g/r/i` or `c/o` coverage, which is sufficient to:

- rule out obvious stellar variables misclassified as transit candidates (pulsators, eclipsing binaries, cataclysmic variables),
- check for additional dips or flares at other epochs that argue for or against a transit-like interpretation,
- extend the baseline for period refinement on periodic candidates using a single combined WD+ZTF+ATLAS ephemeris fit.

This is a cheap, fully automated pre-filter; candidates that do not survive ZTF/ATLAS consistency do not advance to Stage 5.3 telescope follow-up.

### 5.3 External follow-up

For candidates that survive the archival first pass:

- high-cadence photometry at predicted transit windows
- archival imaging and catalog checks (beyond ZTF/ATLAS)
- spectroscopy or RV constraints where physically meaningful
- high-resolution imaging if blending is a concern

Current planning assumption for MIT-affiliated follow-up:

- treat Magellan as the primary institutional follow-up path to develop early
- treat MMT as a collaborator-driven or backup path rather than the main plan
- identify which MIT-access high-speed photometric instrument is actually schedulable for WD transit windows
- build a rapid-response workflow around short predicted windows and sub-minute cadence requirements
- **WINTER + SPRING (MIT-GO, 2026B call ~mid-May)**: primary near-IR follow-up channel. Y/J/Hs (`0.9–1.7 µm`) gives the depth-vs-band lever arm for planet vs BD vs EB-FP discrimination (cool companions add NIR flux that a giant planet does not). Reach: WINTER `J~18.5` in `960 s` (1°×1° survey FOV); SPRING `J~19.5` in `960 s` (9'×7' precision FOV with 0.4″ pixels). Operationally easy: automated reduction+differencing pipeline, Python queue API. Two-tier plan — WINTER for transit-window re-observation, SPRING for ingress/egress shape on high-priority candidates. The 2026B proposal is a forward-looking commitment for the 40-sector cohort (TWIRL's first candidates land post-Stage-4); pre-deadline TODOs: (a) compute J-band reachability of `twirl_wd_master_catalog_v0_tesscoverage.fits` to quote concrete reachable-target counts, (b) anchor proposal narrative on WD 1856 + a faint S100 candidate.

Follow-up readiness criteria should be satisfied before requesting telescope time for a periodic candidate:

- the event survives the baseline search and automated vetting
- the signal is independently re-extracted or cross-checked in a second reduction path
- image-level and aperture-level contamination checks are complete
- the ephemeris is precise enough for the planned instrument and observing window
- the expected cadence and depth are compatible with the planned facility
- the candidate status, assumptions, and main false-positive modes are documented

For apparently aperiodic or irregular events:

- do not assume repeat-photometry confirmation is possible
- define a separate follow-up path centered on classification, continued monitoring, or archival/context constraints

### 5.4 Final population analysis

Regardless of the number of validated planets, TWIRL should finish with:

- occurrence-rate posteriors or upper limits, first in the large-occulting regime where the 200 s survey is strongest
- candidate catalog publication
- validated discoveries, if present

A publishable null result still requires:

- a locked parent sample
- a documented search definition
- end-to-end completeness tied to that search definition
- a final vetted candidate table or null-candidate statement
- upper limits stated only for the declared sample and event class

### Stage 5 Deliverables

- validated-candidate list
- follow-up status table
- occurrence-rate paper inputs

## Suggested Milestones

### Year 1

- finalize WD sample
- stand up MIT PDO light-curve production
- select and run a first `Sector >= 56` benchmark containing WD 1856+534
- recover WD 1856+534 b with the new 200 s workflow
- stand up the baseline periodic and dip-search branches
- create first injection datasets for large occulters
- identify the concrete MIT-affiliated follow-up path and required collaborators/instruments

### Year 2

- run the full search stack and automated vetting
- publish candidate catalog
- execute first follow-up campaign

### Year 3

- complete validation
- compute occurrence rates with completeness corrections, starting with the large-occulting regime
- publish discoveries and/or upper limits

## Immediate Implementation Priorities For This Repo

1. Rebuild targeted S57 TWIRL-FS v2 products and pass dim-target QA, including TIC `1400694779`, before any production restart.
2. Complete product QA for TWIRL-FS v2 with full-sample cadence retention, RMS/MAD, quality flags, WD 1856 timing, and injected short-event preservation.
3. Run BLS/heuristic/LEO checks on the validated TWIRL-FS v2 S56 tree and compare against the prior QLP/TWIRL S56 products.
4. Resume S57-S63 relight and then S64+ production only after v2 QA passes; leave `qc_pause.flag` in place until then.
5. Build the v2/v3/TWIRL-FS comparison table only after the light-curve product is clearly validated.
6. Define the consolidated HDF5-to-TWIRL index format and attach QA metrics.
7. Audit Gaia-first target support and decide what MIT fork changes are needed for targets without TIC IDs.
8. Build the transparent periodic and dip-search baselines before committing to an ML-heavy workflow.
9. Lock down follow-up readiness and coordination only after the candidate-validation criteria are stable.
