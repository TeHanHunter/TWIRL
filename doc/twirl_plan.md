# twirl_plan.md

This document is the executable software and survey plan for TWIRL.

## Current Status Snapshot

- TWIRL is a standalone repo and remains separate from `TWIRL_proposal`.
- The core survey is limited to `200 s` TESS FFIs, which means `Sector >= 56`.
- `WD 1856+534` is the benchmark target, and the first fixed benchmark is `Sector 56`, orbits `119` and `120`, `cam4/ccd1`.
- The seed WD catalog is a local external dependency, not a git-tracked repo asset.
- Light curves will be gathered for all WDs in the full parent catalog. The statistical occurrence-rate sample (the locked denominator) will be defined separately after light-curve collection and QA are complete.
- The search is interpretable-first: build a transparent periodic + dip baseline before considering ML triage.
- `2026-06-02`: **TWIRL I is the active manuscript target.** It is a framework/overview paper, not the final occurrence-rate paper or a discovery paper, titled *TWIRL I: A Systematic TESS Search for Transiting Planetary Remnants around White Dwarfs*. The clean Overleaf-ready scaffold lives in the independent sibling repo `/Users/tehan/PycharmProjects/twirl-survey-paper`; this pipeline repo remains focused on code, data products, QA, and provenance.
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
- `2026-06-18`: **S56-S93 cutout prep/recovery is complete and size/count QA passed.** On `pdogpu1`, the refreshed S56-S93 completion map shows `38/38` sectors with cutout prep complete; the full source-pickle metadata sweep checked `1,216` orbit/camera/CCD source directories with `0` size/count issues. `qc_pause.flag` remains in place, so normal GPU finalize is still intentionally paused while TWIRL-FS v2 product QA remains the gate.
- `2026-06-03`: **Current-stage talk and QA visuals wrapped for the post-talk checkpoint.** The clean local Keynote/PPTX deck is preserved in `outputs/` and the ignored `reports/exploratory/talks/2026-06-02-current-stage/` archive; the tracked repo now keeps only the reusable scripts and compact report artifacts. New presentation-facing QA products are the S56-S93 production-status map and the rebuilt WD 1856 S56 pixel-map diagnostic. The precision-plot work remains exploratory and should not be used as a final product-QA claim until the normalization/plotting choice is re-run and signed off.
- `2026-06-24`: **ORCD 8xH200 downstream-compute path is usable for TWIRL once compact S56 exports are staged.** The runnable Slurm partition is `pg_mki_aryeh`; the current control-socket probe reaches login host `login007` and sees the H200 node `node4900` with `gpu:h200:8`. PDO remains the Stage 1 TGLC/ePSF production home; the first ORCD pilot should move compact S56 TWIRL-FS v2 light-curve exports, manifests, candidate tables, and recovery outputs to `/orcd/data/mki_aryeh/001/twirl/exports/s56_twirlfs_v2/`, then run CPU/1xH200 smoke tests before larger Stage 3/4 sweeps. Operational details live in [ORCD guide](orcd_h200_usage.md).
- `2026-06-30`: **ORCD S56 staging and baseline environment are now in place.** The `20k` raw-flux pre-detrend BATMAN injection HDF5 plus manifest/labels are staged under the ORCD checkout's `data_local/` tree, and the reusable `twirl-s56` Python environment exists under `/orcd/data/mki_aryeh/001/twirl/envs/`. Next ORCD gate is a `100`-injection peak-table smoke plus H200 visibility smoke before scaling the full peak-table/ranker sequence.
- `2026-06-30`: **ORCD smoke gate passed and the full S56 peak/ranker chain is running CPU-only.** The `100`-injection CPU smoke passed with top-20 signal recall `45/100`; the H200 visibility smoke used exactly `1` H200. The full `20k` injected peak-table build and dependent ranker/apply jobs are now submitted with no GPU/GRES requests, and the ORCD helper refuses routine submissions requesting more than `2` H200s without an explicit override.
- `2026-06-30`: **Coverage-first all-host S56 injection generation is complete on PDO.** The all-host pre-detrend BATMAN run produced `19,072` injections across `77` TIC-grouped shards and covers `19,071 / 19,072` unique exported S56 TICs; the lone missed host is documented as a nonpositive-baseline skip. The sharded robust-BLS peak-table pass is now running on PDO with bounded CPU concurrency.
- `2026-07-01`: **Human triage is ready for the all-host ranker-selected real-candidate queue.** The all-host injected peak table and injected-truth ranker selected `57,204` real S56 ephemerides across `19,068` targets; a shuffled `1,000`-target LEO-backed queue now verifies with `1,000` PDFs, `0` LEO metric errors, and `0` LEO plot errors. The PDO app is live on `pdogpu1:5007`; labels write beside the queue for later teacher-model training.
- `2026-07-01`: **Teacher-queue selection has moved off the peak ranker.** The active human-labeling path is now a blinded `10,000`-row mixed teacher pool: `9,000` stratified real S56 BLS/vetter rows plus `1,000` pre-detrend BATMAN injection-recovery rows, with a random `1,000`-row first-pass queue (`900` real, `100` injected). The builder and PDO/local wrappers are in [script](../scripts/stage5_validation/build_s56_mixed_teacher_queue.py); a local no-LEO dry run verifies row counts, source blinding, hidden truth columns, bucket balance, and finite ephemerides. PDO LEO rendering/sync is the remaining operational step before labeling this queue.
- `2026-07-02`: **Downstream testing/vetting shifts to ORCD by default.** PDO remains the Stage 1/TGLC/HLSP production and compact-export staging home, while raw-flux detrending-strength audits, two-aperture BLS/vetting-sheet production, injection-recovery branch tests, and later ML training should run from compact S56 exports on ORCD CPU/H200 resources as appropriate. The first ORCD ADP+ audit is now treated only as a post-ADP vetting/display diagnostic: it shows residual trends matter, but it is not an independent production search product because it starts from already detrended ADP curves. The production search branch must be chosen from raw-flux re-detrending variants, i.e. the same family of comparison as canonical `DET_FLUX` versus `DET_FLUX_ADP`.
- `2026-07-02`: **Raw-flux detrending-strength audit favors small-aperture search.** ORCD CPU job `17020767` compared `7` fresh raw-flux spline settings across `DET_FLUX_ADP_SML`-style small aperture and default/primary aperture on `3,000` BATMAN injections. The best branch is small aperture with `bkspace_d=0.15 d`, `gap_split_d=0.2 d`, quantile knots (`1534/3000 = 51.1%` strict top-1; `1713/3000 = 57.1%` top-N exact/harmonic), only slightly ahead of current small-aperture ADP (`1515/3000 = 50.5%`) and well ahead of default/primary aperture (`1233/3000 = 41.1%`). Next vetting sheets should therefore search on small aperture, show primary aperture as a comparison, and keep the `0.15 d` setting as a candidate production update pending real-data vetter QA.
- `2026-07-04`: **The `0.15 d` candidate branch is now a concrete S56 product.** The named branch is `twirl-fs-v2-adp015q`, compare FITS columns are `DET_FLUX_ADP015*`, and the two-aperture vetter runs BLS directly on `DET_FLUX_ADP015_SML + DET_FLUX_ADP015` without an additional ADP+ high-pass. Full S56 ADP015 FITS production completed on PDO (`19,072 / 19,072`, zero failures), the compact two-aperture export verifies `19,072` target groups, and a bounded ADP015 two-aperture vet-sheet smoke passed (`10/10 ok`). Full queue rendering/search should now move to ORCD from the compact export, while `twirl-fs-v2` remains canonical until real-data QA and WD 1856 checks pass.
- `2026-06-17`: **Julien joins the active collaboration/follow-up planning.** Meeting notes are recorded in progress log [§2.5](twirl_progress_log.md#25-collaboration-meetings-and-ownership). Immediate implications: compare S56 TWIRL-FS search/vetting results against Julien's SPOC Stage-1 candidate funnel, define what signal classes the current products are sensitive to before first-paper claims, and verify follow-up/funding routes before treating SPECULOOS, MISCOT, LCO 1m, EPRV, or proto-Lightspeed as executable paths.
- `2026-05-13`: **TWIRL pivots to a Schwamb-group collaboration.** Michelle Kunimoto brings a well-tuned BLS and LEO-Vetter expertise; her student + Franklin Chen tune LEO-Vetter for WDs in parallel with our `wd-host-tuning` fork. Te Han is the LC producer + data steward (the v3 TWIRL HLSP tree shipped today is the shared survey input) and is offered lead authorship on the **occurrence-rate paper** (verbal — to be locked in writing this week); **catalog paper leadership undecided**. Injection-recovery becomes shared exploratory work with multiple approaches in parallel. See progress log [§2.5](twirl_progress_log.md) for the meeting record and the [Collaboration & Ownership](#collaboration--ownership-2026-05-13) section below for the explicit division of labor and ownership-protection plan.

### Stage 2 note (2026-05-01)

When we open Stage 2 (search), use William Fong's cuvarbase branch + persistent-context M-per-GPU queue architecture as the starting point — he validated S100 cam4/ccd4 equivalence within float roundoff. His numba-vectorized trapezoid fit is also a relevant pattern for any custom CPU-bound search code we write.

## Publication Roadmap

The active writing target is **TWIRL I**, the survey framework and overview paper. It should establish the WD sample motivation, 200 s TESS/TGLC data-product path, TWIRL-FS light-curve strategy, search/vetting framework, WD 1856 benchmark logic, and completeness plan.

Non-goals for TWIRL I:

- final occurrence-rate posteriors or upper limits
- individual discovery claims
- Earth-size or habitable-zone occurrence claims
- a locked statistical denominator beyond the documented current parent-sample assumptions

Those outputs remain later papers or later sections only after the parent sample, full search definition, end-to-end completeness, and candidate validation are complete.

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

This section captures the division of labor and the ownership-protection steps agreed at the Schwamb-group kickoff meeting on `2026-05-13`. Treat as live — update at every substantive collaboration meeting with Michelle, Julien, or other formal collaborators. Companion meeting notes live in progress log [§2.5](twirl_progress_log.md).

### Julien collaboration note (2026-06-17)

Julien is now part of the active TWIRL collaboration/follow-up planning. Treat this as a collaboration-scope change, not just a one-off consultation. The detailed meeting record is in progress log [§2.5](twirl_progress_log.md#25-collaboration-meetings-and-ownership).

Near-term shared work:

1. Compare validated S56 TWIRL-FS search/vetting outputs against Julien's SPOC Stage-1 candidate funnel.
2. Use Julien's prior vetting experience to stress-test the candidate taxonomy: contaminating blends, SB2 / stellar-companion systems, WD+M-dwarf eclipsers, and planet-like occultation candidates should stay separate.
3. Define what signal classes the current S56 products are actually sensitive to before making first-paper sensitivity or yield claims.
4. Clarify follow-up funding/proposal routes, including whether a GI proposal is realistic and whether "Spitzer" is shorthand for a historical validation model, archival context, or a current JWST-style thermal-IR path.
5. Verify instrument access before committing to any path: SPECULOOS, MISCOT/Avi four-band photometry, LCO 1m, EPRV, and proto-Lightspeed via Kevin Burdge are candidate options, not accepted project infrastructure yet.

First-paper collaboration with Julien is in scope, but paper roles and authorship should be confirmed explicitly before writing or validation commitments depend on them.

### Who owns what

| Layer | Owner | Notes |
|---|---|---|
| **TGLC light-curve production** (v3 flux-space-detrended HLSPs at scale) | Te Han (this repo) | Unambiguously ours. Includes the [twirl/preserve-negative-flux](https://github.com/TeHanHunter/TESS_Gaia_Light_Curve/tree/twirl/preserve-negative-flux) TGLC patch and the [flux_space_detrend](../src/twirl/lightcurves/flux_detrend.py) cotrend. Foundation of every downstream result. |
| **Data steward for the survey LCs** | Te Han | Versioned tree, sidecar manifests, public release decisions. |
| **Julien comparison / follow-up planning** | Te Han + Julien | Compare TWIRL S56 outputs against Julien's SPOC Stage-1 funnel, carry over useful vetting lessons, and jointly evaluate follow-up options. Instrument and funding paths are not accepted until access/proposal details are verified. |
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

#### S56 semi-supervised vetting layer (...)

After the `2026-06-17` Julien meeting, the near-term TWIRL-owned ML path is a
candidate-level semi-supervised **self-training / pseudo-labeling** layer, not
a raw-light-curve discovery engine. The first teacher model is trained on
human-vetted S56 candidates plus injection labels. It then assigns
high-confidence pseudo-labels to the larger S56 candidate pool, and a student
model is trained on the union of human/injection labels and explicitly
provenanced pseudo-labels. Pseudo-labels are never treated as ground truth:
they carry confidence, margin, source model, and iteration metadata, and the
pipeline also emits a human-review queue of uncertain or high-value examples.

Initial implementation lives in [module](../src/twirl/vetting/self_training.py),
[driver](../scripts/stage5_validation/run_self_training_triage.py), and
[config](../configs/detection/self_training_s56.yaml). It intentionally uses
engineered BLS/vetter features so it can run on the existing S56 candidate
tables without blocking on deep-learning infrastructure. Deep models can
replace the estimator later once the labeled sample, injection products, and
ORCD/PDO runtime are stable. The deep-learning target should be an
AstroNet-like multi-aperture raw-light-curve branch combined with ExoMiner-like
diagnostic inputs, not a drop-in stock AstroNet model trained on main-sequence
host priors.

Status (`2026-06-18`): the initial S56 v2 human-vetting template is ready in
[CSV](../reports/stage5_validation/self_training_s56_v2/human_labels_template.csv)
with a companion [guide](../reports/stage5_validation/self_training_s56_v2/labeling_guide.md),
but it is now treated as a bootstrap sheet, not the final review substrate. The
preferred review unit is an end-to-end candidate/recovery object with explicit
provenance: TWIRL-FS v2 canonical `DET_FLUX` light curve, injection truth when
present, BLS recovery columns, WD-tuned LEO-Vetter diagnostics, and then a human
label. The browser vetter ([module](../src/twirl/vetting/lightcurve_label_app.py),
[runner](../scripts/stage5_validation/run_lightcurve_vetting_app.py)) now
presents pre-rendered LEO-Vetter reports first; the plain TWIRL light curve is
collapsed fallback/debug context. Existing LEO outputs cover `51 / 300`
bootstrap-template rows, so that sheet remains only a bootstrap/debug aid. The
active pre-human-triage product is a `pdogpu6` full-S56, multi-aperture run
that injects the same signal into `DET_FLUX_SML`, `DET_FLUX`, and
`DET_FLUX_LAG`, runs BLS recovery across all three apertures, and renders
WD-tuned LEO reports using each row's representative aperture. The final
`100` real + `900` injected review queue passed the pre-label verifier on PDO:
`1,000` rows, `1,000` referenced LEO report files, `0` LEO metric errors, WD
1856 exactly once as `PC`, and finite per-aperture injected-row SDE/recovery
columns. The queue lives under
`reports/stage5_validation/s56_pretriage_review_queue_pdo/`, and the PDO
browser app is live in tmux session `twirl-pretriage-vetting-app`, writing
labels to `human_labels_vetted.csv` beside the queue.

Status (`2026-06-22`): the `1,000`-row queue remains the pilot for confirming
that injected signals are visually recoverable. The next production review
sample is a blinded `10,000`-row S56 queue: `9,000` stratified real BLS/vetter
candidates plus `1,000` LC-level injections. The injected subset now spans
small bodies through giant-planet/brown-dwarf-size occulters using finite-
exposure `batman` models, not box dips. The active sampling default is a
balanced period-depth grid with truth metadata for target depth, finite-
exposure model depth, cadence-sampled model depth, radius, impact parameter,
`a/Rs`, inclination, and grid cell. LEO rendering now uses a plot-only
sanitized copy of each metrics object so pathological fitted models no longer
force fallback sheets. PDO entrypoints are the pure-injection [1k launcher](../scripts/stage5_validation/run_s56_1k_batman_review_pdo.sh)
and the mixed [10k launcher](../scripts/stage5_validation/run_s56_10k_blind_review_pdo.sh).

Status (`2026-06-22`, update): the post-detrend queues above are now treated as
UI and plumbing tests only. The active scientifically useful review pilot is
the raw-flux pre-detrend `1,000`-row queue built by [pre-detrend launcher](../scripts/stage5_validation/run_s56_1k_predetrend_review_pdo.sh):
BATMAN signals are injected into raw TGLC `RawFlux`, canonical and ADP
TWIRL-FS detrending are rerun, then BLS and WD-tuned LEO are applied. The
completed two-column comparison queue has `1000/1000` LEO PDFs and `0` LEO
metric/plot errors, but the human-facing pilot has been rebuilt as ADP-only:
`DET_FLUX_ADP` BLS/LEO for all `1,000` rows, `0` LEO errors, and
`32/1000` BLS recoveries. It is served locally on `127.0.0.1:5005`.
Pixel-level injections now have a working source-pickle/ePSF smoke prototype
on PDO, including an ePSF-refit test and a full pixel -> TGLC HDF5 ->
TWIRL-FS -> BLS smoke path, but should be a calibration subset rather than
the immediate human-vetting queue. The full-chain smokes reproduce the current
search weakness: broad wrong-period BLS peaks outrank injected 2-3 minute WD
signals. See [feasibility report](../reports/stage3_injections/s56_pixel_vs_predetrend_feasibility.md).

Status (`2026-06-22`, diagnostic update): the ADP-only pre-detrend pilot is
now a **detrending/search failure diagnostic**, not a completeness sample.
PDO signal-survival measurements show raw aperture injections retain most of
the BATMAN depth at the truth ephemeris (median retention `0.82`), but
`DET_FLUX_ADP` retains only `0.28` at the median and only `0.23` for
`Tmag >= 19`. Strict BLS recovery is strongly magnitude/period dependent:
`22/97` for `Tmag < 18`, `4/684` for `Tmag >= 19`, and `0/400` for
`P >= 2 d`. Re-running BLS with WD-short duration grids improves only
`32/1000 -> 38/1000`. The sharper SNR diagnostic is that `929/1000`
ADP rows have empirical multi-transit SNR `< 7`; BLS recovery rises to
`14/17` only once SNR exceeds `20`, while high-SNR misses are
`bls_peak_mismatch` cases where harmonics or broad long-period structure
outrank the injected period. An `inferred_aperture` baseline smoke did not
materially change recovery (`7/200`), so the next gate is preserving injected
depth through detrending, using the new exact/top-N/harmonic recovery columns
to isolate BLS ranking failures, and suppressing broad systematics peaks before
scaling the 10k human-vetting or recovery sample. A PDO rebuild with the new
columns showed that the old `32/1000` strict ADP recovery at `n_periods=5000`
contains `11` exact top-N recoveries and `25` harmonic top-N matches; increasing
to `n_periods=200000` raises strict top-1 recovery to `50/1000` and recovers
all empirical-SNR `>20` cases, but still leaves `923/1000` unmatched rows.
Thus the search needs both a denser/refined BLS branch for high-SNR aliases and
a detrending/SNR-preservation branch for the dominant faint low-SNR failures.
Signal-survival diagnostics also show the small ADP aperture is the most
promising next search input (`DET_FLUX_ADP_SML` median empirical multi-SNR
`1.58`, `124` rows above SNR 7) relative to medium ADP (`0.89`, `71` above
SNR 7). The next PDO recovery-mode sweep should therefore compare small,
medium, and large ADP apertures before freezing the human-vetting evidence
aperture.

Until that sweep completes, the immediate human-check queue is the SNR-
stratified subset built from the existing ADP-only 1k product:
[queue](../reports/stage5_validation/s56_1k_predetrend_batman_adp_only_snr_stratified_review_queue/review_queue.csv).
It keeps all `32` strict recoveries, all `11` exact top-N-only cases, all
`25` harmonic top-N cases, and a controlled set of unmatched low/mid/high-SNR
rows, with all `300` rows linked to already-rendered LEO PDFs. Use this queue
for near-term human visual checks; use the full 1k for completeness accounting.
Start the local app with
[launcher](../scripts/stage5_validation/run_s56_snr_stratified_vetting_app_local.sh),
which defaults to `127.0.0.1:5006` and writes labels beside the queue.

Status (`2026-06-23`, diagnostic update): a reusable BLS-failure diagnostic
[script](../scripts/stage5_validation/diagnose_bls_recovery_failures.py) now
confirms the current losses are structured, not random. In the ADP-only
pre-detrend 1k queue, any exact-or-harmonic period match is `10/11` for
`Tmag < 16`, `17/68` for `17 <= Tmag < 18`, only `15/659` for
`19 <= Tmag < 20`, and `0/25` for `Tmag >= 20`. Empirical post-detrend
signal SNR explains the trend more directly: any period match is `6/527` for
SNR `<1`, `11/24` for SNR `7-10`, `15/30` for SNR `10-20`, and `17/17` for
SNR `>20`. Therefore the immediate priority remains signal preservation and
aperture/search comparison, not simply scaling the same ADP setup to more
rows. With a recoverable-SNR gate at `7`, the current decomposition is
`904/1000` low-SNR unmatched rows, `68/1000` exact-or-harmonic matches, and
only `28/1000` SNR-qualified BLS-ranking losses. Treat those `28` as the
near-term search-statistic debugging set; treat the `904` as evidence that the
detrender/aperture product is not preserving enough signal for the faint
injected population.

Status (`2026-06-23`, aperture update): all-aperture signal-survival summary
shows the small detrended apertures are the clear next search inputs. The best
single product is `DET_FLUX_ADP_SML` with `124/1000` rows above empirical SNR
`7`, compared with `71/1000` for current medium `DET_FLUX_ADP`; canonical
`DET_FLUX_SML` is close at `117/1000`. Best-of-all detrended apertures reaches
only `137/1000`, so aperture selection roughly doubles the SNR-qualified set
but does not solve the dominant faint low-SNR problem. The next PDO sweep
should prioritize `DET_FLUX_ADP_SML` and `DET_FLUX_SML`; any `10k` scale-up
should wait until this is confirmed through BLS/LEO, not only truth-window
survival.

Status (`2026-06-23`, priority queue update): the next PDO debug run now has a
compact target list in
[queue dir](../reports/stage5_validation/s56_1k_predetrend_debug_priority_queues/).
It contains `94` unique injected rows: `28` SNR-qualified BLS-ranking failures
and `66` rows newly above empirical SNR `7` under their best aperture. This is
the correct targeted set for search-statistic and aperture experiments before
rerendering a larger human-vetting queue. Recommended apertures in this set are
dominated by `DET_FLUX_ADP_SML` and `DET_FLUX_SML`.

Status (`2026-06-23`, targeted sweep update): the recovery-mode sweep driver
now supports `--injection-id-file`, and
[priority PDO launcher](../scripts/stage5_validation/run_s56_predetrend_priority_recovery_sweep_pdo.sh)
uses the `94`-row priority list by default. Run this targeted sweep before the
full 1k sweep when PDO SSH returns; it should answer whether small-aperture BLS
actually recovers the rows that the truth-window survival diagnostic says are
now SNR-qualified.

Status (`2026-06-23`, sweep-decision update): the priority PDO launcher now
runs a post-sweep verifier
[script](../scripts/stage5_validation/summarize_priority_recovery_sweep.py).
Use its `priority_sweep_summary/summary.json` and delta table to decide the
next branch: if small apertures add many matches relative to `adp_priority_5k`,
run the full 1k small-aperture sweep and rebuild the human queue on that
evidence; if not, the limiting step is not aperture choice and the next work
should focus on detrending/search-statistic changes.

Status (`2026-06-23`, targeted sweep result): the `94`-row priority sweep now
confirms the small-aperture branch. Medium ADP recovered only `7/94`
exact-or-harmonic matches (`2` top-1). `DET_FLUX_ADP_SML` recovered `65/94`,
canonical `DET_FLUX_SML` recovered `55/94`, and the small-aperture pair reached
`77/94` with the dense `200k` period grid (`73/94` top-1). The full `1k`
recovery-mode sweep should therefore prioritize `DET_FLUX_ADP_SML +
DET_FLUX_SML` before any larger `10k` queue is built.

Status (`2026-06-23`, detrending-method update): a direct method sweep over
stored original/injected raw TGLC `RawFlux` arrays finds that detrending is
still a sensitivity bottleneck in the full faint-heavy pilot, but no tested
alternate detrending family preserves most injected signals while also removing
most long-timescale trends. On the full `1,000` injection set, the current ADP
small-aperture product reaches `124/1000` rows above empirical SNR `7`, median
depth retention `0.496`, and trend reduction `0.933`. A fine-grid follow-up
finds the best cheap non-oracle branch is a short rolling median:
`median015_gap05` gives `128/1000` rows above empirical SNR `7`, retention
`0.489`, and trend reduction `0.982`; `0.20-0.30 d` median windows are nearly
tied and preserve slightly more depth with more residual trend. This is useful
for a future BLS/LEO comparison, but it is not a completeness-changing
improvement. Percentile/envelope filters, coarser splines, polynomials,
Savitzky-Golay filters, and constant baselines do not improve usable SNR; an
oracle transit-masked ADP spline reaches only `126/1000`, adding just `2` rows
over current ADP. The practical decision remains to prioritize small apertures
and the search grid first, not to replace the detrending method as the main
fix.
Audit: [report](../reports/stage5_validation/s56_1k_predetrend_detrending_method_audit.md).

Status (`2026-06-23`, active 1k queue): the interrupted dense full-`1k`
small-pair sweep has been replaced by a restartable chunked PDO run and is now
complete. The current best injected-recovery table uses `DET_FLUX_ADP_SML +
DET_FLUX_SML` with a `200k` BLS period grid: `137/1000` strict top-rank
recoveries, `177/1000` exact/top-N/harmonic recoveries, and the expected
high-SNR behavior (`35/35` top-rank recoveries above empirical SNR `20`). A
new LEO-only local review queue has been built from this table with `1,000`
rows, `1,000` full WD-tuned LEO reports, `0` LEO metric/plot errors, shuffled
blinded metadata, and class counts `FA=923`, `PC=72`, `FP=5`. Use this queue
for the immediate visual check; it replaces the older medium-ADP and
SNR-stratified pilot queues. It is still not the final `10k` mixed
real+injected sample.

Status (`2026-06-23`, dense sensitivity map result): the dense pre-detrend
BATMAN period-depth grid completed on `pdogpu6` and was synced locally. The run
uses `10,000` raw-flux injections over a `50 x 50` period-depth grid, reruns
TWIRL-FS detrending, then applies the same `DET_FLUX_ADP_SML + DET_FLUX_SML`,
`200k`-period BLS recovery path used by the current 1k queue. The merged table
has `10000/10000` rows across `100` chunks: `1273/10000` strict top-rank BLS
recoveries and `1683/10000` exact/top-N/harmonic period matches. Recovery is
strongly Tmag-conditioned: exact/top-N/harmonic match fractions are `83%` for
`Tmag < 16`, `57%` for `17 <= Tmag < 18`, `27%` for `18 <= Tmag < 19`, and
`9%` for `19 <= Tmag < 20`. The global recovery rate is therefore dominated by
the S56 faint-target distribution and should not be quoted without a Tmag slice.
Artifacts live under [parameter-space plots](../reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/parameter_space/):
the dense point cloud, raw binned recovery fractions, Tmag-sliced fractions,
and a kernel-smoothed 50% BLS sensitivity boundary. The smoothed boundary gives
the clean visual trend requested for methods inspection; the raw binned plots
remain the audit trail.

Status (`2026-06-24`, duration/radius visualization): the dense map now has a
duration-aware companion visualization that plots the 50% BLS recovery cutoff
as companion radius, not transit depth. The adopted visual model is physically
constrained to the BLS signal proxy `R_p^2 sqrt(duration / period)` plus Tmag,
which avoids over-fitting the correlated period-duration-radius sampling. The
recommended artifacts are now the period-radius 50% boundary split by fixed
duration and the empirical period-radius recovery-fraction map split by
duration and Tmag in
[duration-aware plots](../reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/).
These views replace the literal 3D plane as the presentation path: planet
radius is directly readable, duration is handled by panels, and the raw binned
map remains the audit trail. For bright WDs the 50% boundary lies inside the
injected planet/brown-dwarf radius range; by the faint end the cutoff is often
beyond the injected radius range. The LEO comparison on the `1,000`
LEO-rendered injected rows shows LEO is **not** yet a near-complete proxy for
BLS recovery: LEO `PC/FP` has `98.7%` precision relative to BLS exact/top-N/
harmonic matches, but only `42.9%` recall. The metric split shows those
BLS-recovered LEO-FA cases are lower-MES, fewer-cadence, higher-SHP events
relative to the LEO PC/FP subset, so the next LEO step is targeted WD retuning
and an explicit review/recovered-but-not-clean class, not blind threshold
loosening. The defensible methods wording remains: BLS recovery says whether
the search tracks the injected signal; LEO recovery is a stricter
high-confidence vetting subset.

Status (`2026-06-24`, publication-map and LEO-tuning pass): the publication-
facing sensitivity figure is now the marginalized empirical period-radius map
in
[duration-aware plots](../reports/stage5_validation/s56_10k_predetrend_dense_bls_map_pdo/duration_aware/),
not the duration-split audit grid. It uses only local injected support within
each Tmag bin, labels the 50% recovery contour directly, and marks unsupported
parameter space with a dashed support boundary so white/grey areas are not
misread as zero recovery. The LEO smoke test now supports a conservative
`WD_REVIEW` route: `wd_review_high_purity` raises LEO-vs-BLS recall from
`42.9%` to `50.3%` while keeping `91.8%` precision on the injected 1k smoke
set. Balanced/aggressive relaxations increase recall but lose purity quickly,
so the next implementation should add a review/recovered-but-not-clean class
rather than promote those rows directly to PC.

Status (`2026-06-24`, bright-balanced physical map complete): the publication-
facing period-radius map is now the `20k` pre-detrend BATMAN run on a
period-radius grid with transit-conditioned random impact
parameters/inclinations, equal target sampling across `<17`, `17-18`, `18-19`,
and `>19`, and the same `DET_FLUX_ADP_SML + DET_FLUX_SML`, `200k`-period BLS
recovery path. The final table has `20,000` rows, with panel counts
`5,020`, `4,994`, `5,026`, and `4,960`; exact/top-N/harmonic BLS recovery is
`62.8%`, `45.5%`, `32.3%`, and `17.3%` across those bins. The plot explicitly
labels the 50% empirical recovery contour, the Roche-limit period, the local
injection-support boundary, and local mean total transit duration contours.
Interpret the `P, R_p` surface as BLS recovery marginalized over physically
allowed duration/depth scatter at fixed period and radius. Final report:
[summary](../reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/duration_aware/summary.md).

Operational note (`2026-06-24`): the current `20k` bright-balanced recovery
run stays on `pdogpu1` because it is already chunking against data staged on
PDO. The next scale-up should use the ORCD/H200 path after the compact S56
export is built and copied to project storage.

Operational update (`2026-06-30`): the staged ORCD `20k` injection product is
now explicitly treated as a balanced period-radius recovery/ranker grid, not an
all-host sample (`5,323 / 19,072` unique S56 hosts represented). The companion
all-host path is the coverage-first PDO launcher
[script](../scripts/stage5_validation/run_s56_predetrend_all_host_grid_pdo.sh),
which uses `--target-selection-mode shuffled_cycle` to attempt every discovered
raw S56 TIC before repeats, then runs the host-coverage audit before staging to
ORCD. The full all-host injection product is now built; the sharded robust-BLS
peak table is the active gate, and the all-host handoff wrapper will run the
post-BLS gate, injected-truth peak ranker, and ranker-selected real review queue
after the merged peak table verifies.

Status (`2026-06-24` EOD): wrap checkpoint keeps the S56 raw-flux
pre-detrend injection-recovery path as the active methods branch. The best
current search input is still `DET_FLUX_ADP_SML + DET_FLUX_SML` with the
dense `200k` BLS grid; the bright-balanced `20k` period/radius map is the
current sensitivity visual; and the `1k` LEO-only injection queue is the immediate human
visual-check substrate. Next work should keep three labels separate:
BLS-tracked injection recovery, LEO class, and human object-type label. LEO
tuning should be calibrated against that split with a review/recovered-but-not-
clean bucket rather than by blindly loosening PC thresholds.

Status (`2026-06-23`, `10k` mixed queue result): the mixed scientific queue is
now built from a matching real-candidate search, not the older medium/
large-aperture S56 table. The PDO runner
[script](../scripts/stage5_validation/run_s56_small_pair_stage2_10k_pdo.sh)
searched the S56 TWIRL-FS v2 compare tree using `DET_FLUX_ADP_SML +
DET_FLUX_SML`, then applied the heuristic vetter, centroid enrichment, and a
blinded `9,000` real + `1,000` pre-detrend BATMAN injection LEO queue build.
The strict verifier passed on PDO: `10,000` review rows, `10,000` LEO reports,
`0` LEO metric errors, `0` LEO plot fallbacks, source counts `9000` real +
`1000` injection-recovery, and LEO classes `FA=9905`, `PC=82`, `FP=13`.
Local metadata is synced in
[queue dir](../reports/stage5_validation/s56_10k_predetrend_small_pair_200k_review_queue_pdo/);
the full PDF directory remains on PDO because it is `4.5 GB`.

Status (`2026-07-01`, mixed teacher queue reset): the ranker-selected real
queue is now diagnostic only and is excluded from the teacher-label sampling
loop. The actionable queue is built by
[builder](../scripts/stage5_validation/build_s56_mixed_teacher_queue.py):
`9,000` real rows from the S56 TWIRL-FS v2 BLS/vetter table, stratified across
planet-candidate-like, EB/PCEB-like, broad-duration/variability-like,
single-aperture, cadence-alias/systematic, and control buckets, plus `1,000`
pre-detrend BATMAN injection-recovery rows from the balanced period-radius
grid. The first-pass review queue is a random `1,000` rows with exactly `900`
real and `100` injected rows; browser-visible `source_bucket` and `vet_class`
are blinded as `review_candidate`, while `truth_source_kind`,
`truth_source_bucket`, `truth_vet_class`, injection truth, recovery status,
`selection_bucket`, and `selection_weight` remain in the CSV for audit and
training joins. A local no-LEO dry run passed verification. The injected source
grid currently covers `0.120-13 d`, not the requested `0.08-13 d`; treat the
shortest-period gap as a targeted second-batch item if the first-pass audit
needs it.

Pre-human-labeling path:

1. Use the accepted S56 pilot light-curve product:
   `/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare`.
   The canonical product remains TWIRL-FS v2 `DET_FLUX`; model/training
   inputs must preserve all three aperture channels (`DET_FLUX_SML`,
   `DET_FLUX`, `DET_FLUX_LAG`) because aperture behavior is evidence, not a
   nuisance column.
2. Build a compact full-S56 HDF5 export on PDO from the HLSP FITS tree. This is
   the transfer/training boundary; do not move raw TGLC cutouts or FFI products
   for this step.
3. Generate LC-level injected positives with finite-exposure `batman` transit
   models, injecting the same model into all three apertures for each selected
   target. Store injection truth, split, source target, aperture baselines, and
   original/injected arrays.
4. Run the transparent search/vetting stack before human labels: BLS recovery
   across all three apertures, per-aperture recovery status, representative
   aperture selection, WD-tuned LEO-Vetter metrics/reports, and real-candidate
   provenance from the S56 BLS/vetter table.
5. Produce the browser-ready review queue only after the above steps are
   complete. The active first-pass teacher queue is the mixed
   `10,000`-row pool with a random `1,000`-row review subset built by
   [builder](../scripts/stage5_validation/build_s56_mixed_teacher_queue.py)
   and PDO [runner](../scripts/stage5_validation/run_s56_mixed_teacher_queue_pdo.sh).
   Human labels attach to this queue, not to raw light curves, the old
   `300`-row bootstrap sheet, the injection-only pilot, or ranker-selected
   diagnostic queues.
6. Gate before labeling the mixed real+injected training queue: the pool must
   verify as exactly `9,000` real plus `1,000` injected rows; the first-pass
   queue must verify as exactly `900` real plus `100` injected rows; browser
   visible source fields must be blinded; hidden truth/provenance columns must
   remain present; all rows need finite period, epoch, duration, and
   representative aperture; and every browser-served row must have a referenced
   WD-tuned LEO report. The current builder writes its own `verification.json`.
7. Launch the app only after that verifier passes. The active local entrypoint
   after syncing PDO outputs is
   [launcher](../scripts/stage5_validation/run_s56_mixed_teacher_vetting_app_local.sh),
   which writes labels to `human_labels_vetted.csv` beside the queue. After the
   first `1,000` labels, run the mixed-label audit wrapper to decide whether to
   continue randomly or draw a targeted second batch from deficit buckets.

Only after those gates pass should the teacher/student self-training run use
human labels. The interim table model may use engineered BLS, LEO, centroid,
and aperture features, but stock AstroNet is not the accepted model path.
Future deep learning should use an AstroNet-like multi-aperture folded-light-
curve branch combined with diagnostic features, trained/evaluated on the
TWIRL-specific WD injection and human-label distribution.

### Stage 2 Deliverables

- reproducible periodic-search baseline
- reproducible dip-search baseline
- reproducible training set definition
- trained detector checkpoint(s), if justified
- pseudo-label table and human-review queue for semi-supervised vetting
- evaluation report
- inference script that scores the full TWIRL archive

## Stage 3: Injection-Recovery Tests

This is the completeness backbone of the survey and should run through the same search stack used in Stage 4.

Stage 3 should explicitly keep two injection levels separate:

- **LC-level injections**: inject signal models into existing raw or detrended light-curve products before the relevant downstream step. This is the correct fast path for the S56 talk, v2-vs-TWIRL-FS comparisons, detrend signal-preservation tests, and dense grids over period, duration, depth, Tmag, and sector count. It does **not** measure failures in extraction, deblending, aperture selection, or pixel-level contamination.
- **Pixel-level injections**: inject synthetic occultation signals into the cutout/FFI pixel data before ePSF fitting and aperture extraction, then run the real TGLC/TWIRL light-curve builder and search stack. This is slower and should be run on a calibration subset, but it is the completeness layer needed for publication-grade occurrence rates because it includes extraction systematics, crowding, aperture disagreements, and centroid/on-target behavior.
- **Adopt a two-stage protocol**: use LC-level injections for iteration and first-order completeness surfaces; use pixel-level injections to calibrate where the LC-level approximation breaks, especially near the faint limit, in crowded fields, and for candidates where centroid or aperture behavior drives classification.

The ORCD/H200 path should start from compact S56 TWIRL-FS v2 exports, not raw
TGLC or TICA trees. The current transfer scaffold is [export script](../scripts/stage3_injections/export_s56_lc_training_set.py)
plus [injection builder](../scripts/stage3_injections/make_s56_lc_injection_training_set.py)
and [config](../configs/injections/s56_lc_training_export.yaml): export the
`~19k` S56 HLSP-derived light curves on PDO or another spacious working disk,
build injected positive examples in the WD1856-like, short-deep, and
Roche-boundary regimes, move the HDF5/manifest/candidate tables to
`/orcd/data/mki_aryeh/001/twirl/exports/`, and use H200 jobs for LC-level
injection-recovery, GPU search equivalence, and later raw-LC model experiments.

For detrending validation, the current PDO path is [pre-detrend builder](../scripts/stage3_injections/make_s56_predetrend_injection_set.py):
inject into raw TGLC `RawFlux`, rerun canonical/ADP TWIRL-FS detrending, and
feed the resulting injection HDF5 into the existing BLS/LEO review builder.
This is the default path for the next S56 human-vetting sample until the
pixel-level calibration subset is scaled beyond the current [smoke wrapper](../scripts/stage3_injections/run_s56_pixel_injection_smoke.py).
For host-distribution coverage, use the all-host launcher
[script](../scripts/stage5_validation/run_s56_predetrend_all_host_grid_pdo.sh)
and require the host-coverage audit before calling the sample representative of
S56.
The smoke wrapper can now write TWIRL-FS FITS and BLS summaries, and its first
full-chain tests show non-recovery caused by broad wrong-period BLS peaks. Fix
that search/detrending failure before expanding either raw-aperture or
pixel-level injection counts.

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
- **Julien meeting additions (2026-06-17; not yet verified):** evaluate SPECULOOS, MISCOT/Avi four-band photometry, LCO 1m, EPRV, and proto-Lightspeed/Kevin Burdge as possible follow-up channels. Clarify follow-up funding and whether a GI proposal around TWIRL targets is realistic. Treat Spitzer as a meeting-note keyword to disambiguate, not a current observing option, unless the intended meaning is archival or a JWST successor program.

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
9. Build the S56 self-training vetting loop on top of the verified mixed teacher pool: human labels, injection truth, teacher pseudo-labels, student scores, and targeted follow-up review queues.
10. Stage compact TWIRL-FS exports for ORCD and run S56 equivalence tests before scaling H200 injection-recovery or search jobs.
11. Lock down follow-up readiness and coordination only after the candidate-validation criteria are stable.
