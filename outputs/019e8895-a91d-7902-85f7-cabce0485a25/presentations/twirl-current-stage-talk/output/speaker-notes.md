# TWIRL current-stage expert talk speaker notes

## 01. TWIRL current stage

Open by defining the current stage narrowly. TWIRL is not yet presenting a finished survey result or an occurrence-rate claim. The talk is about whether the end-to-end machinery is credible enough to scale: parent catalog, 200 s TESS coverage, MIT TGLC production, TWIRL-FS detrending, BLS search, heuristic vetting, LEO-Vetter collaboration, and pixel-level checks.

The anchor is WD 1856+534 b in Sector 56. It is a known deep, short, white-dwarf transit in the 200 s FFI era. If the pipeline cannot recover that benchmark cleanly, scale-out and completeness claims are premature.

Spend another minute defining the deliverable: this is an editable Keynote-importable PPTX with provenance footers and notes, not a paper draft. Tell the audience that the strongest feedback you want is on the gates: what would make them trust the product enough to scale, and what would still make them stop.

## 02. Table of contents

Use this slide to set expectations for a room that knows both TESS and white-dwarf planets. The first part is literature context, not background filler: it explains why white dwarfs are unusually powerful transit targets and why their signals are easy to overinterpret without strong validation.

The middle of the talk is the TWIRL work product. The S56 benchmark lets us discuss real production, real cadence behavior, real detrending choices, and real BLS/vetting artifacts without implying that the full survey denominator or completeness is finished.

Use the flowchart as a map for questions. If a TESS person asks about extraction, point to the TGLC section. If a white-dwarf planet person asks about validation, point to the WD 1856 and LEO-Vetter sections. This keeps the talk organized while still inviting interruptions.

## 03. Why white dwarfs?

Transition into the science history. The framing is not simply that planets around white dwarfs are interesting. It is that white dwarfs provide a compact, high-contrast stage where planets, debris, accretion, and stellar remnants all intersect.

For this audience, emphasize that the history directly motivates the TWIRL design. WD 1145 shows that transiting debris exists; WD 1856 shows that a planet-scale object can transit a white dwarf; JWST/MIRI now adds a thermal/infrared confirmation axis.

Frame the first section as a shared vocabulary check. The room likely knows the individual systems, but TWIRL depends on connecting them: white-dwarf pollution says planetary material survives, WD 1145 says transiting debris happens, and WD 1856 says planet-scale transits are observable with TESS.

## 04. Small hosts turn modest objects into large photometric events.

Keep this conceptual rather than overclaiming a specific occurrence rate. The key geometric point is that the white dwarf is tiny, so a planet or debris cloud can produce a very deep event. That is why WD 1856 looks so dramatic compared with main-sequence transits.

The same geometry makes cadence central. A signal can be real and still be poorly sampled in a 30-minute FFI product. TWIRL is built around 200 s FFIs because the survey is optimized for large, deep, short-duration events in the WD 1856-like regime.

Pause on the validation burden. The same small-star geometry that gives dramatic transit depths also makes false positives visually tempting. A nearby eclipsing source, a data-product artifact, or a pulsator can look exciting unless the pipeline keeps aperture and pixel evidence attached to every candidate.

## 05. The field moved from circumstellar debris to a planet-scale benchmark.

This is the narrative bridge, now with text large enough to read from the back of the room. Polluted white dwarfs establish that planetary material is present after stellar evolution. WD 1145 gives a transiting debris example. ZTF and related wide-field searches show that debris-transit work is becoming a survey problem, not just a one-object curiosity.

WD 1856 gives a planet-scale object with a much cleaner periodic transit, and JWST/MIRI provides another observational handle. TWIRL is the natural survey question that follows: if 200 s TESS FFIs are now available, can we systematically search a large white-dwarf parent sample for deep short events, while being honest about product quality, vetting, and completeness?

Make the larger history cards do work rather than reading them. Ask the audience to hold two endpoints in mind: evolving debris systems found by K2/ZTF-style searches and stable planet-scale transits like WD 1856. TWIRL needs both the periodic baseline and a future dip/debris branch.

## 06. Debris transits established that white-dwarf systems can be photometrically active.

Use the real publication figure rather than a schematic. The stacked K2 light curves show the transit morphology evolving over the campaign, which is the point TWIRL needs: white-dwarf transit phenomenology includes variable debris as well as stable planet-scale events.

The point for TWIRL is not that WD 1145 is the same search problem as WD 1856. It motivates a broader white-dwarf transit survey that should keep a periodic BLS channel and eventually add a dip-search branch for debris-like or weakly periodic events.

Use the publication figure as the proof object. The point is not to describe every K2 epoch in detail; it is to show that WD 1145 morphology evolves enough that TWIRL should not force all white-dwarf transit behavior into one periodic BLS template.

## 07. WD 1856 is the clean benchmark: deep, short, and repeatable.

Use the all-transits figure to remind the room why WD 1856 is such a strong benchmark. The event is deep and short enough that cadence and aperture choices matter, but it is repeatable enough that a periodic search should recover it.

For TWIRL, this is not just a motivating exoplanet. It is a mandatory technical test. The current first benchmark is Sector 56, orbits 119 and 120, camera 4, CCD 1. A stable recovery in that product is the gate before larger claims.

When the all-transits figure is on screen, orient the audience to the repeated deep events first, then connect the visual depth to the search metric. A good TWIRL product should make this object easy enough to recover that failures become diagnostic rather than ambiguous.

## 08. Ground and Spitzer transits fixed the event shape outside TESS.

This figure is useful for the white-dwarf planet experts because it shows the validation mindset around WD 1856. The discovery was not just a TESS dip; the authors built a consistent multi-instrument transit picture.

TWIRL cannot replicate all of this at the search stage, but it can require that the TESS discovery artifacts are sane before a candidate leaves the pipeline: stable aperture behavior, plausible duration, no obvious alias or cluster pile-up, and a clean pixel/centroid story.

Use this slide to separate discovery from validation. TESS found the event, but the confidence came from compatible follow-up light curves. TWIRL can lower the cost of finding candidates, but the project should still package evidence in a way that makes follow-up decisions rational.

## 09. WD 1856 also teaches that apertures and neighbors are not optional.

Use this slide to connect the literature validation to TWIRL's vetting requirements. In a white-dwarf survey, a deep event can be an astrophysical signal, a contaminant in the aperture, a nearby eclipsing binary, or a product artifact.

That is why TWIRL keeps aperture-to-aperture consistency as a minimum QA item. For WD 1856, the known transit should be visible across the relevant apertures, and the pixel maps should tell the same story before we trust search rankings.

Name the TWIRL translation plainly: every serious candidate should carry small, medium, and large aperture behavior, nearby-source context, and pixel diagnostics. That is the practical inheritance from the WD 1856 literature validation, not a ceremonial citation.

## 10. WD 1856 leaves the formation problem open, which is part of the motivation.

This slide gives the field-level reason the search is interesting to white-dwarf planet experts. WD 1856 raises formation, survival, scattering, and common-envelope questions. The system is not merely another transiting giant planet.

The TWIRL talk does not need to solve the formation channel. It should make clear that a larger sample, or even strong upper limits after a well-calibrated search, would inform whether WD 1856 is rare, representative, or part of a broader class.

This is a good moment to invite the WD experts into the formation question without letting the talk drift. Say that TWIRL's immediate contribution is empirical: either find more systems or build a defensible non-detection after completeness, both of which constrain how rare WD 1856-like systems are.

## 11. JWST/MIRI adds an infrared confirmation axis for WD 1856 b.

Use the MIRI SED figure to update the story beyond the 2020 discovery paper. The important point for this talk is that WD 1856 has become an even stronger benchmark system because independent infrared evidence supports the planet interpretation.

That matters for TWIRL because the benchmark is not arbitrary. If a TESS/TGLC/TWIRL-FS product and search stack can reproduce this system cleanly, the pipeline has passed a relevant, high-value white-dwarf planet test.

Treat the JWST/MIRI slide as a benchmark-strengthening update. The 2025 result does not change the TWIRL pipeline, but it changes the confidence that WD 1856 is the right gate. The target is not merely convenient; it is now even more physically valuable.

## 12. The MIRI result turns WD 1856 from a transit benchmark into a physical planet benchmark.

These panels help close the science-history section. WD 1856 is not only a transit-shape problem; it is a system where atmospheric/thermal constraints are becoming possible. That reinforces why new candidates from TWIRL would be valuable follow-up targets.

At the same time, the slide should not imply that TWIRL can confirm planets by itself. TWIRL is a discovery and triage engine. Strong candidates will need independent extraction, contamination analysis, ephemeris adequacy, and follow-up planning.

Use the two panels to remind the room that discovery value depends on follow-up value. If TWIRL finds a clean analogue, it is not only a TESS candidate; it is potentially an infrared, high-cadence, and atmospheric-characterization target once the ephemeris and contamination story are strong.

## 13. Why TESS and why TGLC?

Transition from science motivation to observability. The key claim is that the cadence and extraction method are the bottlenecks. TWIRL is not simply taking whatever TESS light curve exists; it is using 200 s FFI data and TGLC-style Gaia-prior photometry.

This section should be especially useful to TESS experts: it lays out the survey input, parent sample logic, the TIC/Gaia bridge, and why TGLC's decontaminated aperture products are attractive for crowded white-dwarf fields.

The transition here should be quick. The science case is now established, so the rest of the section asks whether the available TESS product is actually matched to the signal. This is where TESS expertise can most directly improve the project design.

## 14. Most unexplored white dwarfs are a 200 s FFI problem.

Use this as the core TESS motivation slide. The plot makes the proposal problem visual: the well-explored two-minute TOI white dwarfs are a small bright subset, while the 200 s FFI product covers a much larger and fainter Gaia WD population.

This is why TWIRL's core survey is restricted to Sector >= 56. Earlier sectors may be useful for context or special cases, but the production survey input is the 200 s FFI product. That boundary keeps the denominator and completeness problem cleaner.

Make the sampling argument concrete. A 2 minute event can still be poorly represented even at 200 s cadence, while a 10 or 15 minute event becomes much more searchable. This is why the early science target is large and deep events, not small shallow habitable-zone claims.

## 15. Most cataloged WDs have at least one Sector >= 56 footprint.

This is the first survey-scale slide, but keep the caveat clear. Coverage is not the same as a final denominator, and footprint is not the same as a validated light curve. The figure shows why TESS is worth using: the sky coverage is broad enough to support a serious survey.

The high-confidence reference cut is Pwd > 0.75, and Gaia DR3 source_id is the authoritative target identifier. TIC is operational metadata. That distinction matters because the MIT implementation is TIC-oriented while the science sample is Gaia-first.

Add the denominator caution verbally. The coverage map is exciting, but coverage is only the first condition. Each target still needs production, quality flags, cadence retention, light-curve QA, and a defined parent-sample cut before it can enter a statistical denominator.

## 16. The first pass is tuned to deep, short WD events.

This detectability map should be one of the main motivation slides. It explains why the first TWIRL pass is not an Earth-size or habitable-zone occurrence-rate claim: the robust near-term test is large, deep, short events in the WD 1856-like regime.

TIC IDs are still useful operational metadata, especially for MIT TGLC production and TESS bookkeeping, but they should not define the survey sample. The scientific identity remains Gaia-first; the search strategy is constrained by cadence, depth, and product completeness.

Use the detectability figure to keep the scientific claim narrow. The map makes clear that the first credible pass is tuned for deep, short, WD 1856-like events. Earth-like or habitable-zone occurrence-rate language belongs later, after product QA and injection recovery are complete.

## 17. TGLC is attractive because it uses Gaia priors and decontaminated apertures.

The point here is not to sell TGLC abstractly. TWIRL needs a light-curve product that can handle crowded fields and faint WDs, and TGLC's Gaia-prior design is a natural match. The MIT fork is used as a production pipeline, not as a one-target quick extraction.

For downstream search, TWIRL expects HDF5 light curves with decontaminated aperture fluxes in 1x1, 3x3, and 5x5 apertures. That aperture ladder becomes a vetting feature, not just an extraction detail.

This is the TGLC value proposition in one slide. Gaia priors are not just convenient metadata; they help model contamination and produce aperture fluxes that can be compared during vetting. That is why TGLC is better aligned with TWIRL than a black-box one-aperture product.

## 18. The TWIRL section should follow the repo's Stage 1-5 plan.

Use this slide to restructure the TWIRL section around the actual repo plan. Stage 1 is the light-curve production and QA foundation. Stage 2 is the interpretable search stack. Stage 3 is injection recovery. Stage 4 is the full search. Stage 5 is validation and occurrence-rate inputs.

The current talk evidence sits mostly in Stage 1 and the S56 vertical slice through Stages 2 to 5. That distinction matters: S56 is a benchmark and collaboration handoff product, not a final survey denominator. The high-confidence WD reference sample is 359,073 objects, and 339,292 of those have at least one Sector 56-or-later footprint, but the final occurrence-rate denominator is deliberately not frozen yet.

Use the Stage 1-5 roadmap to make the rest of the TWIRL section feel like the project plan rather than a pile of artifacts. Stage 1 is the current production gate, while Stages 2 through 5 are being exercised only through the S56 vertical slice.

## 19. The S56-S94 light-curve scope is a sector-count problem.

This slide replaces the old partial-product dashboard with the catalog-derived scale of the light-curve problem. The bars are the number of unique-TIC WD targets per sector in the local TWIRL sector-summary catalog, from S56 through S94. They are target-sector LC scope, not a claim that every product has been finalized.

The headline numbers are intentionally simple: 2,313,607 unique-TIC target-sector light curves across S56-S94, and 795,096 target-sector entries in the Pwd > 0.75 high-confidence subset. The operational status remains gated: S56 is the handoff product, S57-S60 legacy products exist, S72-S90 prep is near the end, and finalization waits on TWIRL-FS v2 QA.

Spend enough time on the operational assumptions because they define feasibility. The GPU ePSF path, PDO data layout, and HDF5-to-HLSP export are the machinery that lets the project move beyond hand-picked benchmarks. Then repeat the caveat: operational feasibility still needs scientific QA.

## 20. What S56 proves, and what it does not.

This divider is important because the rest of the talk becomes very concrete. Prevent the room from hearing S56 results as final survey claims. The S56 slice is valuable because it contains WD 1856 and because the end-to-end data product exists.

The claims in this section should be phrased as product and benchmark claims: what is built, what is recovered, what failures were found, and what still blocks scale-out.

This divider should reset the standard of evidence. From here forward, every claim is about S56 as a vertical slice. That makes the evidence stronger, not weaker, because it prevents the talk from mixing product-development results with final survey conclusions.

## 21. The current vertical slice runs from catalog to vetting artifact.

Walk through the pipeline as an audit chain. The white-dwarf catalog defines the candidate universe. Sector/orbit/camera/CCD coverage maps decide what TESS data can be used. MIT TGLC produces HDF5 light curves. TWIRL-FS converts those into search-ready HLSP FITS. Search and vetting operate on those products.

This slide is also a collaboration map: the light-curve product and data stewardship are Te Han's responsibility, Franklin owns the BLS/search work shown, and Michelle owns the LEO-Vetter collaboration path.

Use the pipeline loop to make the current state auditable. If someone asks where a number comes from, it should map to a stage: catalog, detector coverage, TGLC HDF5, TWIRL-FS HLSP, BLS table, vetter output, or pixel-map artifact. That is the structure to preserve in future reports.

## 22. The faint-end risk was negative and near-zero background-subtracted flux.

Explain the failure mode at the level needed for both TESS and WD experts. For faint targets, real background-subtracted flux values can be negative or near zero. A magnitude-style conversion or divisive detrending can accidentally turn those cadences into NaNs or unstable ratios.

The TGLC-side patch preserves linear RawFlux and RawFluxError, and TWIRL-FS uses subtractive flux-space detrending. On WD 1856, the RawFlux path preserved 155 previously NaN cadences and increased usable q=0 cadences by about 6.3%, while preserving the transit shape.

Slow down on negative flux because it is easy to underestimate. At the faint end, negative background-subtracted flux is not automatically bad data. Throwing it away can bias cadence retention and search behavior exactly where TWIRL expects many white dwarfs to live.

## 23. The active handoff product is conservative by design.

This is the methods slide for the light-curve product. The canonical search column is DET_FLUX, using robust subtractive flux-space detrending with bkspace_d = 0.8 d and gap splitting at 0.5 d. It is designed to preserve legitimate negative and near-zero flux cadences rather than normalize them away.

The adaptive columns are intentionally opt-in. They can reduce broad residual structure, but they may attenuate long or shallow astrophysical features. That is why the deck presents them as a comparison input, not as a silent replacement for the conservative default.

Explain the conservative product as an intentional compromise. The goal is not the prettiest light curve; it is a search input that preserves short events while avoiding fragile normalization around zero flux. The adaptive columns exist because comparison is useful, but they are not the default scientific product.

## 24. The S56 precision check now uses the expected-flux denominator.

Use this figure to show that the QA is quantitative, not just visual. The important product behavior is how scatter changes with expected flux or magnitude and whether the faint end behaves in a way we understand.

Do not claim that the precision curve alone validates detection completeness. It is one piece of the gate. The product also needs cadence retention statistics, aperture-to-aperture consistency, WD 1856 recovery, independent extraction checks, and injection recovery through the full search and vetting stack.

When presenting the precision plot, point to both panels. The top panel gives an absolute scatter scale; the bottom panel shows behavior relative to an expected-flux trend. Then state the limitation: precision is necessary for search, but it is not the same as recovery completeness.

## 25. The canonical detrend preserves short injected events across the tested grid.

This slide explains why bkspace_d = 0.8 d is the conservative default. The retention sweep tests injected signals and asks whether the detrending removes astrophysical dips along with instrumental structure.

For an expert audience, make clear that this is a product-level smoke test, not the final injection-recovery study. The final completeness calculation must run injections through the actual search, classifier if used, vetter, and candidate-merging path.

Use this figure to justify a parameter choice without overselling it. The canonical setting passes the tested short-event retention grid, which is exactly the regime of the first benchmark. The next question is whether the same choice behaves across sectors and through the actual search stack.

## 26. Longer injected events expose the tradeoff behind adaptive detrending.

This figure helps prevent overconfidence. A tighter detrend can make light curves look cleaner while suppressing longer astrophysical signals. That is acceptable only if the product and search channel are explicitly scoped.

TWIRL's first-year science is optimized for large, deep, short-duration events in the WD 1856-like regime. Longer or debris-like events motivate a separate dip-search branch and dedicated injection tests, not silent reuse of a period-search detrend setting.

The long-event slide is where you show scientific restraint. A method can be excellent for WD 1856-like events and still be wrong for longer or debris-like signals. That is why the project should split search branches instead of pretending one detrend setting optimizes everything.

## 27. Gap splitting fixed a real v1 failure mode at the faint end.

This is a concrete example of why the project is benchmark-gated. TWIRL-FS v1 looked plausible until a dim S57 target revealed that one spline across a multi-day orbit gap could fail and fall back to a weak polynomial.

TWIRL-FS v2 splits gaps larger than 0.5 days and keeps the subtractive flux-space approach. In the 18 dim examples, all negative quality-zero RawFlux cadences stayed finite, and the median binned-trend amplitude reduction was about 0.796.

This dim-target example is a story about finding a failure before scale-out. Make that positive but technical: v1 failed a real dim case, the failure mode was understandable, and v2 encodes the fix. That is exactly how a benchmark-gated pipeline should mature.

## 28. The 100-target S56 v2 compare gallery is the product-level visual audit.

This slide is useful for a TESS-heavy room because it shows that the product is being inspected as a product, not only through the benchmark target. The QA sample includes all 16 CCDs and both canonical and adaptive comparison columns.

Use the gallery to say what remains open: this visual audit supports the handoff, but the full-product QA still has to quantify cadence retention, failure modes, RMS/MAD behavior, quality-flag behavior, and aperture consistency before scale-out claims.

Use the gallery as evidence of breadth. It is not a proof of correctness, but it shows that the handoff was not judged by WD 1856 alone. Invite the TESS experts to identify any patterns they would want quantified before a sector-scale production decision.

## 29. The benchmark event is visible in the S56 search/vet product.

This figure is the bridge from light curves into search. It is the concrete WD 1856 vet sheet from the S56 BLS workflow. The event is visible in the search artifact, which is exactly the kind of object the collaboration can review together.

Be careful with wording. This is a benchmark recovery artifact, not a final blind survey result. The value is that it makes the discussion auditable: experts can point to the plot and ask about cadence masks, apertures, period aliases, and nearby sources.

On the vet sheet, orient the viewer to the raw/folded/search panels and say why this is the collaboration audit surface. The figure is not just a plot; it is the object that lets Franklin, Michelle, Te Han, and external experts critique the same candidate consistently.

## 30. Pixel maps are the next validation layer after aperture agreement.

This slide anchors the centroid/pixel-map branch of the vetting plan. A BLS peak and a folded light curve are not enough if the event could come from a neighboring source. Pixel-level diagnostics ask whether the transit-like behavior is spatially consistent with the WD.

The current TWIRL state includes these diagnostics for WD 1856. The next step is to turn them into a repeatable candidate-vetting module rather than an individual benchmark plot.

The pixel-map slide should connect directly to false positives. A signal that is deep in the aperture but spatially inconsistent with the WD is not a planet candidate. The practical next step is to standardize this diagnostic for every high-priority BLS candidate.

## 31. Search and vetting are now a collaboration interface.

This divider explicitly sets collaboration credit and ownership. It should prevent ambiguity in the room. Franklin's BLS work is the search artifact being showcased. Michelle owns LEO-Vetter and the collaboration path around it. Te Han is responsible for TGLC/TWIRL-FS production and data stewardship.

The technical goal is not to defend one tool. It is to define clean interfaces: search tables, vet sheets, LEO input/output schemas, pixel-map diagnostics, and shared candidate review criteria.

Use this divider to protect collaboration clarity. The talk should not blur roles. Franklin gets explicit credit for the BLS/search work shown. Michelle's ownership of LEO-Vetter is named. Te Han's ownership of production and data stewardship is named. That wording should stay stable.

## 32. Franklin's BLS work gives TWIRL a concrete periodic-search artifact.

Showcase Franklin's work clearly. The BLS search converts the TWIRL-FS light curves into ranked periodic candidates and vet sheets. This is the first discovery engine in the talk because it is interpretable and directly matched to WD 1856-like events.

The current repo logs record a denser S56 BLS v2 configuration with p_min = 0.12 d and durations [3, 4, 5, 6, 8, 10, 13, 16, 20, 30] minutes. WD 1856's blind rank improved from 1258 to 560 in that run, while the heuristic vetter already carried most of the post-vetting gain.

This is Franklin's showcase slide. Give enough technical detail to make the work credible: the BLS search generates ranked candidates and vet sheets, the denser duration grid was tested, and the WD 1856 rank behavior is tracked. Then hand off naturally to vetting.

## 33. The heuristic vetter is a physics-motivated bridge, not the final arbiter.

Explain the heuristic vetter as a transparent bridge between BLS peaks and LEO-Vetter. It uses physics-motivated cuts on duration, period aliases, and period-cluster pile-ups to reduce obvious false positives without pretending to be a full validation engine.

The current benchmark result is strong enough for a talk: WD 1856 moved from blind rank 1258 out of 19,040 to planet-regime rank 9 out of 5,403, about a 140x improvement. That justifies the vetter as an audit layer while LEO-Vetter is tuned for white-dwarf hosts.

Describe the heuristic vetter as a transparent audit layer. It is valuable because experts can understand why a candidate moves up or down. It should not replace LEO-Vetter or human review; it should reduce clutter and make the LEO/human workload more scientifically targeted.

## 34. LEO-Vetter classifies the WD 1856 benchmark as PC in the tuned path.

Use the wording carefully. Michelle owns LEO-Vetter. TWIRL's role is to provide clean light curves, BLS/candidate artifacts, WD-specific metadata, and the questions that need tuning for compact hosts and 200 s cadence.

This slide should be treated as a major proof object. It shows the full LEO-Vetter style report for WD 1856 in the S56 tuned path: raw and detrended views, phase diagram, primary/odd/even/secondary checks, model properties, and the PC classification context.

For the LEO-Vetter slide, keep Michelle's ownership explicit and respectful. TWIRL should supply stable inputs and WD-specific requirements, while Michelle's LEO-Vetter path handles the validation layer. The full WD 1856 report is the key proof object because it shows what a collaborator can audit.

## 35. The handoff should be artifact-based, not memory-based.

This is the collaboration slide requested by the user. Te Han owns the TGLC/TWIRL-FS production tree, product QA, provenance, and data stewardship. Franklin owns the BLS/search work shown, including vet sheets and candidate tables. Michelle owns LEO-Vetter and the collaboration path for WD tuning.

The practical recommendation is to make the interface artifacts stable: HLSP tree and metadata, BLS candidate tables, vet sheets, LEO input/output tables, pixel-map diagnostics, and a shared candidate-review sheet.

Use the collaboration flowchart as an operating plan. The most important phrase is artifact-based handoff. Every stage should leave a table, figure, or metadata record that another collaborator can inspect without relying on a meeting memory or a private notebook.

## 36. The first real candidate list is mostly a false-positive taxonomy exercise.

This slide makes the next scientific work concrete. The top candidates are not automatically planet candidates. They must be sorted into classes: post-common-envelope binaries, WD+WD or compact binaries, centroid-shift contaminants, ZZ Ceti pulsators, systematic ladder peaks, debris-like dips, and possible planet-regime events.

For the expert audience, this is a useful discussion prompt. Ask which failure modes they expect to dominate a 200 s TESS WD search and which ones should block scale-out until the vetting machinery is more mature.

The false-positive taxonomy slide can become a discussion slide if time allows. Ask the room which classes they expect to dominate a 200 s WD search. Their answer can directly shape BLS thresholds, LEO tuning, centroid checks, and follow-up triage.

## 37. Completeness must be measured through the real stack.

This slide protects the scientific claims. TWIRL should not present occurrence rates or Earth-size/habitable-zone claims until completeness is demonstrated. Injection recovery must go through the actual chain: light curve product, search, classifier if any, vetter, and candidate merging.

There are two useful levels. Fast light-curve-level injections can test detection and vetting cheaply. Pixel-level injections are slower but necessary for extraction, crowding, aperture choice, and detrending failures. The talk should invite expert feedback on how to prioritize those layers.

For injections, emphasize that the project needs both speed and realism. Fast light-curve injections are ideal for development, but a paper-grade completeness claim eventually needs to include extraction and crowding failures. This keeps the near-term plan useful without pretending it is the final calibration.

## 38. Scale-out waits for a short list of falsifiable gates.

Use this as the actionable next-steps slide. The gates are TGLC production QA, S56 TWIRL-FS v2 search, BLS candidate tables, LEO-Vetter tuning, centroid and pixel-map vetting, injection recovery, and a separate dip-search branch.

The point is to make the scale decision falsifiable. If the product fails full QA, we do not scale. If WD 1856 is unstable across apertures or vetters, we debug the benchmark. If injections show strong attenuation or vetting loss, we tune before occurrence-rate work.

The gates slide should feel like a work plan. Each gate is falsifiable: a report either exists or it does not; WD 1856 either remains stable or it does not; injections either quantify recovery or they do not. That makes the path to scale-out concrete.

## 39. Scale only after the benchmark evidence is stable.

Close with the decision logic. The survey is compelling because the science payoff is high and TESS 200 s coverage is broad, but the project should scale only after the benchmark evidence is stable.

The final discussion questions for the room are: what would convince you that WD 1856 is recovered well enough; where should independent validation enter the loop; which TESS/WD failure modes are most dangerous; and what minimum injection-recovery design is credible for first occurrence-rate claims?

Close by asking for decisions, not applause. The deck should leave the audience with a small number of technical questions: what evidence is enough for WD 1856, which false positives are most dangerous, how should Michelle's LEO-Vetter interface be packaged, and what injection design is credible.
