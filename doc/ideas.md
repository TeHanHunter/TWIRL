# ideas.md

This file is for brainstorming, referee-style concerns, unresolved decisions, and future directions.

Nothing in this file is automatically a project decision.
Project decisions belong in `doc/twirl_plan.md`, `AGENTS.md`, or another formal repo doc once they are accepted.

## Critical Before Next Stage

These items genuinely block forward progress and should be resolved before Stage 2 begins.

0. **QUALITY-flag vs negative-flux correlation in v3 HLSPs (open, 2026-05-13)**: Inspecting the v3 QC PDF on faint S56 sample LCs (T=18-20), all negative DET_FLUX cadences for TIC 2053458501 have `QUALITY == 2`. WD 1856's negatives (155 in SAP) have `QUALITY == 0`. So for *some* TICs the v3 preservation is invisible downstream because the QC mask correctly drops QUALITY≠0 cadences. Open question: is QUALITY=2 (a) inherited from TICA pointing flags (legit; v2-vs-v3 visually identical because both drop these cadences, just for different reasons), or (b) set by TGLC because flux was negative (over-aggressive; the v3 patch is blocked downstream by an upstream quality assignment). Resolution: grep TICA QUALITY bits for these cadences and compare to TGLC HDF5 QUALITY. If (b), the TGLC fork needs a second patch to stop flagging legitimate negative cadences. Not urgent for the talk — WD 1856 demonstrates the v3 win cleanly without this — but worth resolving for the methods paper and before the survey-wide reprocessing.

1. **QSP detrend signal preservation**: The detrend step (`bkspace_min=0.3` days) needs validation against 5–15 minute WD transit signals. A referee will immediately ask whether the spline erases the signal we are searching for. Injection-recovery tests on the detrend step itself — injecting 2–15 min box transits into real WD light curves before detrending, then checking recovery — must happen before any negative search result can be trusted. (Week 2 of the `2026-05-12 → 2026-06-02` sprint partially answers this — the injection grid pushes 5–30 min signals through the real detrend; see `twirl_plan.md` §Three-Week Pre-Talk Vertical Slice. **2026-05-13 update**: this is now shared exploratory work with Michelle's group — multiple injection-recovery approaches will run in parallel rather than converging on a single canonical method, so the answer will come from a multi-method comparison rather than a single pipeline. Our specific angle is the v2-QLP-HLSP vs v3-TWIRL-HLSP completeness delta, isolating the contribution of the negative-flux preservation patch. **2026-05-17 update**: the first detrend-only smoke test is positive for the robust-auto subtractive model on 10-minute injections, but this is not a full injection-recovery grid and should not be treated as a closed validation item.)
2. **WD 1856 benchmark QA (Stage 1.6–1.7)**: Read `267574918.h5` from orbits `119` and `120`, verify transit timing against the published ephemeris, compare `1x1`/`3x3`/`5x5` apertures, and log the result. This gates everything downstream.
3. **Light-curve consolidation (Stage 1.5)**: Build the TWIRL HDF5 index with target identifiers, orbit, camera/CCD, and summary QA metrics. Without this, Stage 2 has no stable input.
4. **QLP Sector 94+ ingestion decision**: Determine whether TWIRL can consume QLP/MAST products directly for Sectors ≥ 94 (see Stage 1.8). This significantly affects PDO compute planning.
5. **`twirl_parent_sample` freeze trigger**: The locked occurrence-rate denominator must be defined before any injection-recovery or occurrence-rate work begins. Agree on the freeze criteria now so Stage 3 can start cleanly.

---

## Scientific Questions TWIRL Can Definitively Answer

### Primary: Occurrence rate of large transiting companions around WDs
TWIRL's headline deliverable: what fraction of WDs have a giant-planet-to-brown-dwarf-sized companion (R > few R_Earth) at short orbital periods (hours–days)? WD 1856b is the only confirmed case. Even a null result with ~350,000 WDs and documented completeness is the best constraint in existence.

**Why 200s cadence is a threshold requirement, not a convenience:** Transit duration for a WD companion is 2–15 minutes (WD radius ~1 R_Earth; Keplerian speed at P=1 day ~350 km/s). At 200s cadence, a 10-minute transit has ~3 data points. At 30-min cadence (Kepler, earlier TESS), that same transit has zero to one data points — undetectable. Prior systematic WD transit surveys (Faedi+2011, van Sluijs & van Eylen 2018, Barber+2022) were structurally blind to sub-hour transits. TWIRL is the first survey that closes this window at scale.

### Occurrence of disintegrating/debris transiting systems
WD 1145+017 and ~10 other systems show chaotic, quasi-periodic dips from disrupted planetesimals. These may be **more achievable in year 1** than detecting a new giant planet:
- Deep dips (5–60%), irregular, no clean period needed — the dip-search branch is the right tool
- Only ~10 known; this is almost certainly a severe detection-completeness problem, not a true rarity
- ~25–50% of WDs are metal-polluted from accreted rocky material — the statistical connection between transiting debris and atmospheric pollution is the key science
- A TWIRL dip search over 350,000 high-confidence WDs could increase the known sample by 10–100×
- Recent BD+05 4868 (2024–25, MS host) demonstrates that present-generation pipelines can recover comet-like dust-tail signatures end-to-end; same methodology should transfer to WD hosts. Confirmed as Stage-4 secondary goal in `doc/twirl_plan.md` (`Search target` line, 2026-05-01).

### What happens to planetary systems at stellar death?
The occurrence rate tests physical models of post-main-sequence planetary system evolution:
- Tidal migration efficiency during the AGB phase (does the planet spiral in or get ejected?)
- Mass-loss-driven orbital expansion (star loses ~half its mass; planets move outward ~2×)
- Kozai-Lidov secular excitation from outer companions
- Post-AGB dynamical instability and planet-planet scattering timescales

A measured rate (or strong upper limit) discriminates directly between competing theoretical models.

### WD cooling age as a variable
Teff from Gentile Fusillo → cooling age → does the occurrence rate differ for young WDs (hot, Teff > 20,000 K, age < 100 Myr) vs old WDs (cool, Teff < 6,000 K, age > 5 Gyr)? This traces dynamical evolution across Gyr in a single snapshot survey. No other method has this leverage.

### The post-main-sequence habitability question
A WD at Teff ~5000 K has its habitable zone at P ≈ 4–32 hours — directly within the 200s sensitivity window. A 5000 K WD stays in this range for ~8 Gyr. Even a strong upper limit on Earth-sized companions in the WD HZ constrains whether post-main-sequence habitability is possible at statistical rates. This framing has media and broader scientific appeal.

### The most complete picture of the Solar System's fate
The Sun will become a WD in ~5 Gyr. What happens to Earth, Venus, Jupiter? TWIRL's constraints on close-in companions, rocky debris, and atmospheric pollution directly speak to this. No other survey addresses it with this combination of sample size and time resolution.

---

## Things We Have Not Fully Thought Through

### The dip search should be co-equal with the periodic search in year 1
The current plan lists dip search as a "branch." It should be treated as the primary year-1 discovery path alongside the periodic search. The science case is strong, the yield potential is high, and it requires no period measurement. See the debris section above.

### The search statistic for WD transits is probably not BLS
BLS assumes many transits per sector and transit duration much shorter than period. For P=5 days, 5-minute transit, 200s cadence: ~5 transits per sector, 2–3 data points per transit. BLS loses sensitivity rapidly here.
- A **matched filter** convolved with a boxcar of the expected WD transit duration is more sensitive for sparse-transit cases
- **The WD radius is estimatable from Teff + Gaia parallax to ~5–10%** — much better than for FGK stars. This means the expected transit duration can be predicted from the host, reducing the search parameter space. Build this in from the start.
- Before committing to a search statistic, do injection-recovery tests specifically for 2–15 minute transits in 200s data

### Multi-sector phase-folding should be the primary search product
The plan processes orbit-by-orbit for production (correct). But the detection search should operate on multi-sector aggregated light curves, not sector-isolated data. A WD with 30 sectors has ~390 transit opportunities vs ~13 in one sector; sensitivity scales as √N_transits. The coverage map is already built — build the search input format around multi-sector stacks from the beginning.

### TESS CVZ WDs are the highest-priority search targets
WDs in TESS continuous viewing zones (|ecliptic latitude| > ~78°) accumulate 13+ sectors of continuous 200s data. For P < 1 day, that is >365 transit opportunities stacked. These should be sorted to the front of the search queue. The coverage map is already built — extract the CVZ subset and report the count.

### The polluted WD cross-match is a zero-cost science multiplier
Cross-match the master catalog with SDSS/LAMOST/DESI metal-polluted WD catalogs (Koester+2005, Zuckerman+2003, Xu+2023) and AllWISE W4-excess (dusty disk) catalogs to create a "polluted WD" and "IR excess" flag. These are the highest-value dip-search targets: disk geometry implies near edge-on viewing, active accretion implies fresh disruption. The cross-match is one afternoon of work and could focus the dip search where the prior probability is highest. Then ask: do polluted WDs show higher dip rates? Does transit depth correlate with pollution strength?

### WD pulsations as a free byproduct
ZZ Ceti (DAV) WDs pulsate with periods of 100–1200 seconds. At 200s cadence, TWIRL Nyquist-samples periods > 400s — covering the bulk of the DAV instability strip. A systematic ZZ Ceti catalog (Lomb-Scargle at post-detrend stage for all WDs) requires no extra software and is publishable independently. Immediately valuable to the WD asteroseismology community.

### False positives specific to WDs
Main false positives for WD transit searches differ from FGK surveys:
- **WD+M dwarf (PCEB) systems**: very common, eclipses look exactly like giant planet transits. Cross-match with Rebassa-Mansergas+2016 SDSS PCEB catalog and Parsons+2016 before vetting.
- **WD+WD binaries**: double WDs produce deep, short, flat-bottomed transits that mimic WD 1856-like signals. Flag known double WDs from SB2 RV surveys.
- **Centroid displacement as a first-tier vetting step**: TGLC's PSF photometry helps, but centroid shift as a function of phase is the cleanest background-EB diagnostic. Should be first-tier, not a late-stage check.

### What exact search statistic should the periodic branch use?
### How should multi-sector candidate merging be handled?
### How should false alarms be estimated?
### What detrending approach best preserves minute-scale WD events?
### How should the `1x1`, `3x3`, and `5x5` apertures be combined or prioritized?
### At what stage, if any, does ML add real value beyond the transparent baseline search?

---

## Catalog and Sample Questions

- How should Gaia-to-TIC matching be implemented and versioned?
- What should happen to the `764` Gaia WDs with no usable TIC bridge? Characterize their magnitude and sky distribution before deciding whether to extend the MIT fork to emit Gaia-selected targets directly.
- Should the repo support Gaia-first target emission if the MIT fork misses important targets?
- Should pilot QA use only `Pwd > 0.75`, or compare several `Pwd` cuts from the start?
- Which fields from the Gaia seed catalog are mandatory in the TWIRL master catalog?
- **QLP ingestion for Sectors ≥ 94**: Can TWIRL consume QLP/MAST 200-s products directly instead of re-running MIT PDO extraction for newer sectors? Assess: magnitude limit applied, Gaia-WD target coverage, format compatibility with the HDF5 archive. Resolving this could substantially reduce PDO compute load.

---

## Benchmark Questions

- Which `Sector >= 56` should become the official WD 1856 benchmark?
- What exact pass/fail metric should define successful WD 1856 recovery?
- How similar does an independent extraction need to look before the benchmark is considered trustworthy?

---

## Follow-Up Questions

- Which MIT-access instrument is the most realistic high-cadence follow-up path?
- Is Magellan enough on its own, or will collaborator access still be required?
- What cadence and timing precision are required before a periodic candidate is follow-up ready?
- What is the confirmation path for irregular or non-repeatable events?
- How should apparently irregular or aperiodic events be handled scientifically?

---

## Referee-Style Scientific Challenges

The biggest likely referee objections are:

1. The parent sample may not be frozen tightly enough.
2. The sensitivity of `200 s` TESS data may be overstated for very short or shallow WD events.
3. **Completeness may not be measured end to end through the real search and vetting stack** — especially if the detrend erases short signals.
4. Gaia-to-TIC matching may introduce an uncontrolled sample bias.
5. The survey may combine highly heterogeneous targets and sectors too casually.
6. A candidate may be a TESS or TGLC artifact rather than a real transiting object.
7. A detected object may not be uniquely classifiable as a planet.
8. WD 1856 recovery may be treated too casually as proof of global sensitivity.
9. Follow-up may be underpowered if cadence, ephemeris precision, and contamination checks are not explicit.
10. A null result may be overinterpreted if the event class is not clearly defined.
11. **For debris transits**: the connection to atmospheric pollution must be explicitly stated; a survey finding no new disintegrating systems at the current detection threshold is not a contradiction of the ~50% pollution rate.

---

## Concerns About Overreach

- Do not oversell Earth-size or habitable-zone sensitivity early.
- Do not reuse catalog-identification completeness as transit-search completeness.
- Do not treat a benchmark recovery as a substitute for injection-recovery.
- Do not claim a planet when the data only support a transiting or occulting object.
- Do not claim survey-wide limits for a provisional or drifting denominator.
- Do not claim the dip-search is complete for the debris class without characterizing the sensitivity to the range of observed dip depths and durations.

---

## Strategic Ideas

- Keep the statistical survey and exploratory discovery products separate.
- Consider making the first occurrence paper about an observable event class rather than a physical companion class — e.g., "large periodic occultation events" not "giant planets."
- Consider publishing a candidate catalog plus upper limits even if no objects are confirmed.
- Consider a parallel science thread on irregular dips or debris-like events — this could be the faster first-paper path.
- Consider benchmarking MIT TGLC against at least one independent extraction path before scaling up.
- Consider a photometric WD characterization paper as a Stage 1 byproduct: noise floor vs Tmag, contamination rate vs crowding across sectors. Useful to the full WD+TESS community.

---

## Possible Future Additions

- Explicit referee-risk mitigation checklist
- Publication roadmap by paper (methods/instrument paper, debris-transit discovery paper, occurrence-rate paper)
- Fellowship-facing framing ideas
- Instrument-specific follow-up notes
- Candidate triage taxonomy
- WD pulsation catalog design (post-detrend Lomb-Scargle)
- IR excess / dusty disk cross-match design (AllWISE W4)
- Metal-polluted WD cross-match design (SDSS/LAMOST/DESI)

### TODO: robust per-pixel pixel-map vetting tool

The TGLC `mitpdo` branch already builds a per-pixel residual LC cube as a
side product of `tglc.effective_psf.fit_lc` (lines ~199-202): the
`aperture` array is shape `(n_cadences, n_pixels_y, n_pixels_x)` where
each pixel is the post-ePSF-subtraction residual. That IS the prebuilt
pixel map — what's missing is a vetting wrapper that consumes it for
on-target / off-target classification.

The TWIRL-canonical on-target test today (2026-05-13) is the **flux-weighted
centroid shift** at the aperture level
([src/twirl/vetting/centroid_offset.py](../src/twirl/vetting/centroid_offset.py)),
which is a single (Δx, Δy) per target. Cheap (~1 ms per TIC), runs on
the existing HLSP `SAP_X`/`SAP_Y` columns, no dependency on SPOC TPFs,
covers 100% of survey targets at FFI cadence — good Tier-1 filter.

But it's a scalar test, so it misses information available in the full
pixel cube. A robust Tier-2 / paper-grade pixel-map vetting tool should:

1. Use `effective_psf.fit_lc`'s `aperture` cube (or compute it directly
   from the cutouts at `/pdo/users/tehan/tglc-gpu-production/orbit-*/
   ffi/cam*/ccd*/source/source_*.pkl`).
2. For each pixel in the cube, phase-fold its LC at the candidate
   `(P, T0, dur)` and measure the dip amplitude (e.g. fit a box of
   the BLS duration centered on T0 in phase).
3. Build a 2D dip-amplitude map and a 2D significance map (per-pixel
   MAD-based sigma).
4. Compare the **peak of the dip map** to the target's known pixel
   location:
     * On-target: peak aligned with target PSF (within ~0.5 pix), dip
       falls off radially matching the PSF model
     * Off-target / bg-EB: peak displaced from target by > 1 pix
     * Pulsator: dip pattern matches the full PSF rather than a
       single-pixel transit (the variability is intrinsic, so it scales
       with the target's flux distribution rather than appearing as a
       spatially localized signal on top of the PSF)
5. Fit the dip-amplitude map to a 2D PSF model centered on the target
   versus centered on each plausible neighbor (from the Gaia / TIC
   crossmatch already in `source.gaia`); report `pixel_offset_arcsec`
   and which star best matches the dip pattern.

The aperture-depth ladder is NOT a substitute (retracted 2026-05-13;
see progress log §2.4). LEO's pixel `offset` test is the proper
analogue but requires SPOC TPFs that most faint WDs lack.

Time / scope estimate: ~1-2 days to build cleanly + ~few hours to run
on the v2 heuristic-vetter survivors (per-cutout amortizes nicely;
~600 TICs share one cutout in a CCD).

Trigger: do this before the first publication-ready candidate list,
not before the talk on `2026-06-02`.

---

## Maintenance Rule

When a brainstorming item becomes a decision:

- move it into `doc/twirl_plan.md` if it changes current status or priorities
- move it into `AGENTS.md` if it becomes operating protocol
- move it into `doc/twirl_plan.md` if it changes the formal pipeline or scientific plan
