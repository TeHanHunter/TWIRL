# Historical QLP-Ops transition plan (superseded)

> **Historical plan, superseded.** Drafted `2026-04-17` and never approved for
> execution. The accepted Stage 1 product is now A2v1; use the
> [A2v1 production protocol](a2v1_production_protocol.md). This file remains at
> its original path only for provenance and old links. Its open questions and
> commands are not current work.

## Context

**Today's state.** TWIRL's Stage 1 benchmark (Sector 56, orbits 119 + 120, WD-only) is essentially complete: 19,072 / 19,086 HLSP FITS (99.93%) live under `/pdo/users/tehan/tglc-deep-catalogs/hlsp_s0056/`. That run deviates from the QLP Ops "Photometry Roadmap" procedure — it skips `qlp lcdb ingest ffi / lightcurve / qflag`, skips `qlp patools qflag`, skips `qlp lctools qualityflag` + manual review, and uses SPOC flags at HLSP time.

**What this plan is.** A **transition** from the bespoke WD-only pipeline to the QLP-Ops "Photometry Roadmap" procedure for subsequent production. Motivations, per today's decisions:

1. The QLP tooling is already tuned for TGLC, so following the roadmap verbatim should give us **faster production** than our current bespoke wrapping.
2. **LC scope stays WD-only** — we are not producing LCs for every TIC down to Tmag 20. But the roadmap's `tglc catalogs` + `tglc cutouts` + `tglc epsfs` stages are run over the full QLP-default catalog for the CCD; only the final `tglc lightcurves` stage is filtered to WDs. That way the **cutouts and ePSFs become a reusable artifact**: anyone (us, the community) can run `tglc lightcurves` later for any requested TIC on that CCD without redoing cutouts/ePSFs.
3. **Stop at HLSP.** We do not run the Planet Search Roadmap tools (BLS, AstroNet, vetting) here; those are out of scope.
4. Outputs go to a **new folder** (exact path TBD — this is why execution is paused).

**What this plan is not.**

- Not a parallel branch to be maintained indefinitely. The existing `hlsp_s0056/` stays as a validated benchmark-of-record; we do not keep re-running the bespoke path alongside the roadmap path.
- Not a science-scope change. TWIRL's science sample is still the WD subset; the LCs we produce are WD-only.
- **Not an all-TIC HLSP archive.** Only cutouts and ePSFs are produced for all TICs on the CCD; LC extraction and HLSP FITS remain filtered to WDs.

## Proposed output layout (folder name TBD — see blocker)

```
<NEW_ROOT>/                              # e.g. /pdo/users/tehan/qlp_standard_s0056/ — user to decide
    orbit-119/ffi/cam{1..4}/ccd{1..4}/
        catalogs/                        # tglc catalogs — all-TIC (QLP default)
        cutouts/                         # tglc cutouts — all-TIC (shareable artifact)
        epsfs/                           # tglc epsfs    — all-TIC (shareable artifact)
        LC/*.h5                          # tglc lightcurves — WD-only (TWIRL science)
    orbit-120/...                        # same layout
    qflags/cam{C}ccd{D}_qflag.txt        # qlp lctools qualityflag + manual review
    hlsp/                                # sector-level HLSP FITS, WD-only
    logs/
```

The existing `/pdo/users/tehan/tglc-deep-catalogs/hlsp_s0056/` is not modified; it's frozen as the benchmark-of-record.

## Execution sequence (per Photometry Roadmap, Sector 56 orbits 119 + 120)

Activate the QLP 0.13.2 env via `source scripts/activate_qlp_env.sh`.

1. **Orbit setup.** TWIRL-local equivalent of `./new-qlp-orbit.sh` under `<NEW_ROOT>/`.
2. **`tglc catalogs`** per `(orbit, cam, ccd)` — QLP-default catalog (all TICs on the CCD, down to the roadmap's mag limit, typically ~Tmag 16). Reuse sibling-orbit catalog symlinking (CLAUDE.md, MEMORY).
3. **`tglc cutouts`** — all-TIC. This is the shareable artifact.
4. **`tglc epsfs`** — all-TIC. Shareable artifact.
5. **`tglc lightcurves` — WD-only.** Filter the target list down to WDs (reuse the catalog join already used by the bespoke pipeline; the TIC → WD crossmatch lives in `src/twirl/catalogs/`). Only these LCs become h5.
6. **`qlp patools qflag`** — binned quaternions from `/pdo/qlp-data/quats/*.quat.txt`.
7. **Cadence-alignment fix** — run `scripts/stage1_lightcurves/fix_tglc_quat_cadence_mismatch.py` on the new h5 tree before detrend. A fresh `tglc all` will re-introduce the cam3/orbit-120 TGLC-only cadence `699957` (CLAUDE.md).
8. **`qlp lctools detrend -b`** — on the new tree only.
9. **`qlp lctools qualityflag`** — produce `cam{C}ccd{D}_qflag.txt`. **Manual review** via `qflag-plot` (port 5000) on each of the 16 CCD/orbit combinations, plus peer review. Budget a 1–2 day human session; coordinate with QLP ops on the review SOP.
10. **DB ingestion (`qlp lcdb ingest ffi / lightcurve / qflag`) — conditional.** See Open Question Q2 below. If skipped, filesystem outputs are still complete and HLSP still runs.
11. **`qlp lctools hlsp`** at sector level with `--flag-type spoc --flag-source fits` (Sector 56 < 67). Writes into `<NEW_ROOT>/hlsp/`.
12. **Recovery-rate audit** (CLAUDE.md §"Standard recovery rate check"): requested WD TIC IDs → h5 count → FITS count, binned by Tmag.
13. **WD 1856 photometric cross-check.** TIC 267574918 HLSP from the new tree vs the existing `hlsp_s0056/`; overlay, MAD-RMS per aperture, confirm agreement within a few percent of the current 0.109 value (twirl_plan.md §1).

## Blocker — pick before kickoff

**Folder location.** Not yet decided. Options:

- `/pdo/users/tehan/qlp_standard_s0056/` — simplest, keeps everything under our existing user root.
- A coordinated location under `/pdo/qlp-data/` if QLP ops want this to look like a regular orbit in their tree. Requires explicit agreement and probably symlinking from `/pdo/users/tehan/` per CLAUDE.md's PDO write-policy rule ("only edit files under `/pdo/users/tehan/`").
- Something under `/pdo/users/tehan/tglc-deep-catalogs/qlp_standard_s0056/` so it lives next to the benchmark tree.

## Open questions to resolve before kickoff

**Q1. Magnitude limit for cutouts/ePSFs.**
- LCs are WD-only, but cutouts/ePSFs run over the QLP-default catalog for the CCD. The QLP-Ops mainline is typically `--max-magnitude 16`. A brighter limit makes the cutouts/ePSFs smaller and faster; a deeper limit makes them useful for fainter on-demand extractions later.
- Candidate values: 16 (QLP default), 17, 18, or match our WD floor (~20).

**Q2. DB ingestion — in or out?**
- The roadmap writes LCs + qflags into the QLP Postgres DB. Downstream `qlp` tooling (`count_lcs_db.py`, `qlp lctools bls`, etc.) reads from there.
- Since we are *not* running the Planet Search Roadmap here, DB ingestion buys only **operator parity** (same verification tooling) but costs DB coordination.
- Options: (a) skip entirely — filesystem-only; (b) ingest into the shared QLP ops DB under a TWIRL-tagged method (coordinate first); (c) stand up a TWIRL-local DB instance.

**Q3. qflag review human-loop.**
- 16 CCDs × 2 orbits = 32 qflag reviews via `qflag-plot`, plus peer review. Who does the peer review — a QLP ops operator, or an internal TWIRL second pair of eyes?

**Q4. Relationship with existing `hlsp_s0056/`.**
- Recommend **freezing it as the Sector 56 benchmark-of-record**, marked read-only in `doc/twirl_plan.md`. TWIRL science downstream switches to pulling WD TICs out of the new tree once that tree passes the same validation; until then, `hlsp_s0056/` remains the science input.

**Q5. Sharing the cutouts/ePSFs.**
- If the ePSFs are meant to enable community on-demand LC extraction, what's the delivery mechanism? MAST? A public rsync mirror? Internal-only for now? This affects disk layout and naming more than the production run itself, but should be settled before we commit to a directory schema.

## Critical files / artifacts to add (once blocker resolved)

- **`scripts/stage1_lightcurves/run_qlp_standard_orbit_pipeline.py`** — argparse driver; per `(orbit, cam, ccd)` runs steps 2–9 above. Mirrors the existing `run_tglc_orbit_pipeline.py`; key difference is the WD filter applies only at the `tglc lightcurves` stage, not at `cutouts` / `epsfs`.
- **`scripts/stage1_lightcurves/README_qlp_standard.md`** — runbook (commands, envs, output paths, manual qflag review checklist, recovery-rate check).
- **`doc/twirl_plan.md`** — add Stage 1 subsection "1.9 QLP-standard production (WD LCs + shareable cutouts/ePSFs)" with scope and stop-point. Update `## Current Status Snapshot`.
- **`doc/twirl_progress_log.md`** — dated entry recording the decision, rationale, and folder convention.
- **No changes to `src/twirl/`** — operator tooling only.
- **No changes to `scripts/activate_qlp_env.sh`** — reused as-is.

## Risks / things to flag

1. **Cutout/ePSF disk footprint.** The cutouts+ePSFs are the heavy artifact; their size scales with the mag cut. Budget disk under `/pdo/users/tehan/` before starting, not after.
2. **Cadence-alignment is mandatory.** A fresh `tglc all` will re-introduce the cam3/orbit-120 TGLC-only cadence (`699957`). The fix script must run **before** `qlp lctools detrend`, or detrend crashes with the known shape mismatch (CLAUDE.md).
3. **Catalog reuse cross-tree does not apply.** The existing `tglc-deep-catalogs/orbit-{119,120}/…/catalogs/` were built from the WD-filtered input and can't be copied/symlinked into the new tree. Sibling-orbit reuse works *within* the new tree (link orbit-119 catalogs into orbit-120 once built).
4. **`qlp patools qflag` input.** Verify `/pdo/qlp-data/quats/cam{C}_quat.txt` is present for Sector 56 orbits 119 + 120 before starting; gaps will make `qlp patools qflag` fail.
5. **Manual qflag review is a human bottleneck.** 32 reviews + peer review. If the reviewer schedule is not arranged in advance, the pipeline will block at step 9.
6. **No science-claim broadening.** WD-only LC scope is unchanged; this is a production procedure change, not a scope change. State this in `doc/twirl_plan.md` so the distinction is explicit.

## Verification (end-to-end)

Once a single CCD (suggest cam4/ccd1, the WD 1856 CCD) has run through to HLSP:

1. **Cutout/ePSF sanity.** Spot-check cutout cube dimensions + ePSF shape against the existing bespoke run's outputs for the same CCD; they should be equivalent.
2. **WD-LC recovery.** Requested WD TIC IDs for cam4/ccd1 → h5 count → HLSP FITS count, binned by Tmag. Expect ≥99% like the benchmark.
3. **WD 1856 parity.** Stitch TIC 267574918 HLSP FITS from both trees, overlay, compute per-aperture MAD-RMS. Expect agreement to within a few percent of the current 0.109 value.
4. **Cadence check.** `h5ls` a handful of new-tree LCs; confirm `LightCurve/Cadence` ⊆ `camC_quat.txt`.

Only after cam4/ccd1 passes do we scale to the remaining 15 CCDs × 2 orbits.
