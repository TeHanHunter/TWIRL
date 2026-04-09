# ideas.md

This file is for brainstorming, referee-style concerns, unresolved decisions, and future directions.

Nothing in this file is automatically a project decision.
Project decisions belong in `doc/twirl_plan.md`, `AGENTS.md`, or another formal repo doc once they are accepted.

## Critical Before Next Stage

These items genuinely block forward progress and should be resolved before Stage 2 begins.

1. **WD 1856 benchmark QA (Stage 1.6–1.7)**: Read `267574918.h5` from orbits `119` and `120`, verify transit timing against the published ephemeris, compare `1x1`/`3x3`/`5x5` apertures, and log the result. This gates everything downstream.
2. **Light-curve consolidation (Stage 1.5)**: Build the TWIRL HDF5 index with target identifiers, orbit, camera/CCD, and summary QA metrics. Without this, Stage 2 has no stable input.
3. **QLP Sector 94+ ingestion decision**: Determine whether TWIRL can consume QLP/MAST products directly for Sectors ≥ 94 (see Stage 1.8). This significantly affects PDO compute planning.
4. **`twirl_parent_sample` freeze trigger**: The locked occurrence-rate denominator must be defined before any injection-recovery or occurrence-rate work begins. Agree on the freeze criteria now so Stage 3 can start cleanly.

## Referee-Style Scientific Challenges

The biggest likely referee objections are:

1. The parent sample may not be frozen tightly enough.
2. The sensitivity of `200 s` TESS data may be overstated for very short or shallow WD events.
3. Completeness may not be measured end to end through the real search and vetting stack.
4. Gaia-to-TIC matching may introduce an uncontrolled sample bias.
5. The survey may combine highly heterogeneous targets and sectors too casually.
6. A candidate may be a TESS or TGLC artifact rather than a real transiting object.
7. A detected object may not be uniquely classifiable as a planet.
8. WD 1856 recovery may be treated too casually as proof of global sensitivity.
9. Follow-up may be underpowered if cadence, ephemeris precision, and contamination checks are not explicit.
10. A null result may be overinterpreted if the event class is not clearly defined.

## Current Open Scientific Questions

- What exact WD cuts should define the final TWIRL parent sample?
- What companion or event class should the first occurrence-rate statement target?
- Should the first statistical paper focus on "large periodic occulters" rather than "giant planets" specifically?
- How should apparently irregular or aperiodic events be handled scientifically?
- What is the right balance between broad exploratory search and clean statistical inference?

## Search-Method Questions

- What detrending approach best preserves minute-scale WD events?
- How should the `1x1`, `3x3`, and `5x5` apertures be combined or prioritized?
- What exact search statistic should the periodic branch use?
- What exact statistic or trigger should the dip-search branch use?
- How should multi-sector candidate merging be handled?
- How should false alarms be estimated?
- At what stage, if any, does ML add real value beyond the transparent baseline search?

## Catalog And Sample Questions

- How should Gaia-to-TIC matching be implemented and versioned?
- What should happen to the `764` Gaia WDs with no usable TIC bridge? Characterize their magnitude and sky distribution before deciding whether to extend the MIT fork to emit Gaia-selected targets directly.
- Should the repo support Gaia-first target emission if the MIT fork misses important targets?
- Should pilot QA use only `Pwd > 0.75`, or compare several `Pwd` cuts from the start?
- Which fields from the Gaia seed catalog are mandatory in the TWIRL master catalog?
- **QLP ingestion for Sectors ≥ 94**: Can TWIRL consume QLP/MAST 200-s products directly instead of re-running MIT PDO extraction for newer sectors? Assess: magnitude limit applied, Gaia-WD target coverage, format compatibility with the HDF5 archive.

## Benchmark Questions

- Which `Sector >= 56` should become the official WD 1856 benchmark?
- What exact pass/fail metric should define successful WD 1856 recovery?
- How similar does an independent extraction need to look before the benchmark is considered trustworthy?

## Follow-Up Questions

- Which MIT-access instrument is the most realistic high-cadence follow-up path?
- Is Magellan enough on its own, or will collaborator access still be required?
- What cadence and timing precision are required before a periodic candidate is follow-up ready?
- What is the confirmation path for irregular or non-repeatable events?

## Strategic Ideas

- Keep the statistical survey and exploratory discovery products separate.
- Consider making the first occurrence paper about an observable event class rather than a physical companion class.
- Consider publishing a candidate catalog plus upper limits even if no objects are confirmed.
- Consider a parallel science thread on irregular dips or debris-like events.
- Consider benchmarking MIT TGLC against at least one independent extraction path before scaling up.

## Concerns About Overreach

- Do not oversell Earth-size or habitable-zone sensitivity early.
- Do not reuse catalog-identification completeness as transit-search completeness.
- Do not treat a benchmark recovery as a substitute for injection-recovery.
- Do not claim a planet when the data only support a transiting or occulting object.
- Do not claim survey-wide limits for a provisional or drifting denominator.

## Possible Future Additions

- explicit referee-risk mitigation checklist
- publication roadmap by paper
- fellowship-facing framing ideas
- instrument-specific follow-up notes
- candidate triage taxonomy

## Maintenance Rule

When a brainstorming item becomes a decision:

- move it into `doc/twirl_plan.md` if it changes current status or priorities
- move it into `AGENTS.md` if it becomes operating protocol
- move it into `doc/twirl_plan.md` if it changes the formal pipeline or scientific plan
