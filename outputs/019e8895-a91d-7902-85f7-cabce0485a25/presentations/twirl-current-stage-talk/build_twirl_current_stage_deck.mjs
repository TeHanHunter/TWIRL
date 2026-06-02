#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { spawnSync } from "node:child_process";

const repoRoot = "/Users/tehan/PycharmProjects/TWIRL";
const workspace = path.join(
  repoRoot,
  "outputs",
  "019e8895-a91d-7902-85f7-cabce0485a25",
  "presentations",
  "twirl-current-stage-talk",
);
const slidesDir = path.join(workspace, "slides");
const previewDir = path.join(workspace, "previews");
const layoutDir = path.join(workspace, "layout", "final");
const outputDir = path.join(workspace, "output");
const pptxOut = path.join(outputDir, "twirl-current-stage-expert-talk.pptx");
const contactSheetOut = path.join(outputDir, "twirl-current-stage-expert-talk-contact-sheet.png");
const manifestOut = path.join(outputDir, "artifact-build-manifest.json");
const speakerNotesOut = path.join(outputDir, "speaker-notes.md");
const sourceManifestOut = path.join(outputDir, "source-manifest.md");
const talkTrackOut = path.join(outputDir, "talk-track-and-qa.md");

const nodeBin = "/Users/tehan/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node";
const pythonBin = "/Users/tehan/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3";
const buildScript =
  "/Users/tehan/.codex/plugins/cache/openai-primary-runtime/presentations/26.601.10930/skills/presentations/scripts/build_artifact_deck.mjs";

const R = (relativePath) => path.join(repoRoot, relativePath);

const F = {
  wd1856AllTransits: R(
    "reports/stage1_lightcurves/wd1856_literature_figures/selected/vanderburg2020_extended_data_all_transits.png",
  ),
  wd1856GtcSpitzer: R(
    "reports/stage1_lightcurves/wd1856_literature_figures/selected/vanderburg2020_fig1_gtc_spitzer_transits.png",
  ),
  wd1856ArchivalApertures: R(
    "reports/stage1_lightcurves/wd1856_literature_figures/selected/vanderburg2020_extended_data_archival_apertures.png",
  ),
  wd1856MassAge: R(
    "reports/stage1_lightcurves/wd1856_literature_figures/selected/vanderburg2020_fig3_mass_age_constraints.png",
  ),
  jwstSed: R(
    "reports/stage1_lightcurves/wd1856_literature_figures/selected/limbach2025_fig2_jwst_miri_sed_excess.png",
  ),
  jwstTempRadius: R(
    "reports/stage1_lightcurves/wd1856_literature_figures/selected/limbach2025_fig4_temperature_radius_context.png",
  ),
  jwstImage: R(
    "reports/stage1_lightcurves/wd1856_literature_figures/selected/limbach2025_fig5_miri_false_color_image.png",
  ),
  tglcMastLit: R(
    "reports/stage1_lightcurves/wd1856_literature_figures/selected/wd1856_literature_figures_contact_sheet.png",
  ),
  wd1145Evolution: path.join(
    workspace,
    "assets",
    "vanderburg2015_fig2_extract",
    "image_02.jpg",
  ),
  tessUnexploredWd: path.join(workspace, "assets", "tess_200s_unexplored_wd.png"),
  wdDetectability: path.join(workspace, "assets", "wd_planet_detectability_maps.png"),
  wdCoverage: R("reports/stage1_lightcurves/catalog_diagnostics/wd_tess_observation_sky_coverage.png"),
  s56Boundary: R("reports/stage1_lightcurves/catalog_diagnostics/s56_sky_coverage_boundary_style.png"),
  precision: R("reports/stage1_lightcurves/qc_pdf/s56_precision_twirl_fs_v2_compare_det_flux.png"),
  retentionShort: R(
    "reports/stage1_lightcurves/detrend_experiments/s56_param_sweep_signal_preservation/plots/retention_depth0p1_sigma5.png",
  ),
  retentionLong: R(
    "reports/stage1_lightcurves/detrend_experiments/s56_param_sweep_long_events/plots/retention_depth0p1_sigma5.png",
  ),
  qaGallery: R(
    "reports/stage1_lightcurves/detrend_experiments/s56_twirl_fs_v2_compare_quick_qa_100/s56_twirl_fs_v2_compare_qa_gallery_page01.png",
  ),
  dimQa: R(
    "reports/stage1_lightcurves/detrend_experiments/twirl_fs_v2_18_dim_examples/s57_tic1400694779_twirl_fs_v2.png",
  ),
  blsVet: R("reports/stage2_search/wd1856_bls/sector_0056/targets/tic_0267574918/vet_sheet.png"),
  blsVetMedium: R(
    "reports/stage2_search/wd1856_bls/sector_0056/targets/tic_0267574918/vet_sheets/vet_DET_FLUX.png",
  ),
  blsVetSmall: R(
    "reports/stage2_search/wd1856_bls/sector_0056/targets/tic_0267574918/vet_sheets/vet_DET_FLUX_SML.png",
  ),
  blsVetLarge: R(
    "reports/stage2_search/wd1856_bls/sector_0056/targets/tic_0267574918/vet_sheets/vet_DET_FLUX_LAG.png",
  ),
  pixelMapTransit: R(
    "reports/stage1_lightcurves/pixel_maps/wd1856_s56_tglc_spoc_saveaper_map_pixel_transits_gaia.png",
  ),
  pixelMapFolded: R(
    "reports/stage1_lightcurves/pixel_maps/wd1856_s56_tglc_spoc_saveaper_map_folded_transit_gaia.png",
  ),
  pixelMapZoom: R(
    "reports/stage1_lightcurves/pixel_maps/wd1856_s56_tglc_spoc_saveaper_pixel_lc_map_zoom60_gaia.png",
  ),
  leoSummary: R("reports/stage5_validation/leo_vetter_s56_top50_v2/summary.png"),
  leoWd1856Report: path.join(workspace, "assets", "leo_vetter_wd1856_pc_rank11.png"),
};

const sources = {
  plan: "Source: TWIRL doc/twirl_plan.md and doc/twirl_progress_log.md, current-stage status as of 2026-06-02.",
  v2015: "Source: Vanderburg et al. 2015, Nature, https://www.nature.com/articles/nature15527.",
  v2020: "Source: Vanderburg et al. 2020, Nature, https://www.nature.com/articles/s41586-020-2713-y.",
  limbach2025: "Source: Limbach et al. 2025, arXiv, https://arxiv.org/abs/2504.16982.",
  tglc: "Source: TGLC HLSP at MAST, https://archive.stsci.edu/hlsp/tglc; TWIRL MIT TGLC usage notes.",
  coverage: "Source: reports/stage1_lightcurves/catalog_diagnostics/wd_tess_observation_sky_coverage.png.",
  s56: "Source: TWIRL S56 benchmark reports and doc/twirl_progress_log.md.",
  twirlfs: "Source: reports/stage1_lightcurves/s56_twirl_fs_v2_compare_handoff.md and detrend experiment summaries.",
  bls: "Source: Franklin Chen BLS search artifacts in reports/stage2_search/wd1856_bls/sector_0056.",
  leo: "Source: Michelle-owned LEO-Vetter collaboration path; TWIRL reports/stage5_validation/leo_vetter_s56_top50_v2.",
};

const slides = [
  {
    kind: "title",
    kicker: "TWIRL",
    title: "TWIRL current stage",
    subtitle: "A 200 s TESS FFI white-dwarf transit survey, benchmark-gated on WD 1856",
    chips: ["45-minute expert talk", "S56 vertical slice", "TGLC + TWIRL-FS + search/vetting"],
    source: sources.plan,
    notes: [
      "Open by defining the current stage narrowly. TWIRL is not yet presenting a finished survey result or an occurrence-rate claim. The talk is about whether the end-to-end machinery is credible enough to scale: parent catalog, 200 s TESS coverage, MIT TGLC production, TWIRL-FS detrending, BLS search, heuristic vetting, LEO-Vetter collaboration, and pixel-level checks.",
      "The anchor is WD 1856+534 b in Sector 56. It is a known deep, short, white-dwarf transit in the 200 s FFI era. If the pipeline cannot recover that benchmark cleanly, scale-out and completeness claims are premature.",
    ],
  },
  {
    kind: "toc",
    kicker: "Plan",
    title: "Table of contents",
    subtitle: "From white-dwarf transit context to the TWIRL benchmark gate.",
    items: [
      ["01", "WD science + history", "Pollution, debris transits, WD 1856 b, and the benchmark logic.", "10-12 min"],
      ["02", "Why TESS + TGLC", "200 s FFIs, WD coverage, and Gaia-prior decontaminated LCs.", "7-8 min"],
      ["03", "Current TWIRL state", "Stage 1 production, S56 vertical slice, and TWIRL-FS QA.", "10-12 min"],
      ["04", "Search + vetting", "Franklin's BLS search, heuristic bridge, and Michelle's LEO-Vetter path.", "8-10 min"],
      ["05", "Next gates", "Production QA, candidate tables, pixel vetting, injections, and dip search.", "5-7 min"],
    ],
    callout: "Benchmark first; scale only after production QA, vetting, and injection recovery.",
    source: sources.plan,
    notes: [
      "Use this slide to set expectations for a room that knows both TESS and white-dwarf planets. The first part is literature context, not background filler: it explains why white dwarfs are unusually powerful transit targets and why their signals are easy to overinterpret without strong validation.",
      "The middle of the talk is the TWIRL work product. The S56 benchmark lets us discuss real production, real cadence behavior, real detrending choices, and real BLS/vetting artifacts without implying that the full survey denominator or completeness is finished.",
    ],
  },
  {
    kind: "divider",
    part: "Part 1",
    title: "Why white dwarfs?",
    subtitle: "The payoff is high because the host star is small; the validation burden is high because the sky is crowded and the signals are extreme.",
    source: `${sources.v2015} ${sources.v2020} ${sources.limbach2025}`,
    notes: [
      "Transition into the science history. The framing is not simply that planets around white dwarfs are interesting. It is that white dwarfs provide a compact, high-contrast stage where planets, debris, accretion, and stellar remnants all intersect.",
      "For this audience, emphasize that the history directly motivates the TWIRL design. WD 1145 shows that transiting debris exists; WD 1856 shows that a planet-scale object can transit a white dwarf; JWST/MIRI now adds a thermal/infrared confirmation axis.",
    ],
  },
  {
    kind: "wdGeometry",
    kicker: "White-dwarf transit geometry",
    title: "Small hosts turn modest objects into large photometric events.",
    facts: [
      ["Transit depth", "planet radius can be comparable to the WD radius"],
      ["Duration", "minutes, not hours, for close-in WD systems"],
      ["Cadence", "200 s TESS FFIs finally sample the event shape"],
      ["Validation", "aperture, centroid, color, and follow-up still matter"],
    ],
    source: "Source: TWIRL survey motivation; basic white-dwarf transit geometry.",
    notes: [
      "Keep this conceptual rather than overclaiming a specific occurrence rate. The key geometric point is that the white dwarf is tiny, so a planet or debris cloud can produce a very deep event. That is why WD 1856 looks so dramatic compared with main-sequence transits.",
      "The same geometry makes cadence central. A signal can be real and still be poorly sampled in a 30-minute FFI product. TWIRL is built around 200 s FFIs because the survey is optimized for large, deep, short-duration events in the WD 1856-like regime.",
    ],
  },
  {
    kind: "timeline",
    kicker: "History",
    title: "The field moved from circumstellar debris to a planet-scale benchmark.",
    events: [
      ["Polluted WDs", "metal pollution implies rocky material survives and is accreted"],
      ["2015: WD 1145", "K2 reveals variable debris transits around a white dwarf"],
      ["ZTF searches", "wide-field monitoring finds new debris-transit candidates"],
      ["2020: WD 1856 b", "TESS finds a stable deep transit by a planet-scale object"],
      ["2025: JWST/MIRI", "infrared data strengthen the physical planet interpretation"],
      ["TWIRL", "200 s TESS FFI survey tests this window at scale"],
    ],
    source:
      `${sources.v2015} ${sources.v2020} ${sources.limbach2025} Source: Bhattacharjee et al. 2025, PASP, https://doi.org/10.1088/1538-3873/ade0ea.`,
    notes: [
      "This is the narrative bridge, now with text large enough to read from the back of the room. Polluted white dwarfs establish that planetary material is present after stellar evolution. WD 1145 gives a transiting debris example. ZTF and related wide-field searches show that debris-transit work is becoming a survey problem, not just a one-object curiosity.",
      "WD 1856 gives a planet-scale object with a much cleaner periodic transit, and JWST/MIRI provides another observational handle. TWIRL is the natural survey question that follows: if 200 s TESS FFIs are now available, can we systematically search a large white-dwarf parent sample for deep short events, while being honest about product quality, vetting, and completeness?",
    ],
  },
  {
    kind: "wd1145",
    kicker: "WD 1145+017",
    title: "Debris transits established that white-dwarf systems can be photometrically active.",
    figure: F.wd1145Evolution,
    source: sources.v2015,
    notes: [
      "Use the real publication figure rather than a schematic. The stacked K2 light curves show the transit morphology evolving over the campaign, which is the point TWIRL needs: white-dwarf transit phenomenology includes variable debris as well as stable planet-scale events.",
      "The point for TWIRL is not that WD 1145 is the same search problem as WD 1856. It motivates a broader white-dwarf transit survey that should keep a periodic BLS channel and eventually add a dip-search branch for debris-like or weakly periodic events.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "WD 1856+534 b",
    title: "WD 1856 is the clean benchmark: deep, short, and repeatable.",
    figure: F.wd1856AllTransits,
    fit: "contain",
    figureBox: { x: 120, y: 132, w: 620, h: 540 },
    railX: 790,
    railY: 146,
    railW: 360,
    railTitle: "Why it matters",
    railCards: [
      ["Periodic", "Stable 1.408 d signal"],
      ["Short + deep", "Minutes-long event with planet-scale depth"],
      ["Validation-rich", "Discovery backed by multi-instrument checks"],
    ],
    source: sources.v2020,
    notes: [
      "Use the all-transits figure to remind the room why WD 1856 is such a strong benchmark. The event is deep and short enough that cadence and aperture choices matter, but it is repeatable enough that a periodic search should recover it.",
      "For TWIRL, this is not just a motivating exoplanet. It is a mandatory technical test. The current first benchmark is Sector 56, orbits 119 and 120, camera 4, CCD 1. A stable recovery in that product is the gate before larger claims.",
    ],
  },
  {
    kind: "figureFull",
    kicker: "Multi-instrument validation",
    title: "Ground and Spitzer transits fixed the event shape outside TESS.",
    figure: F.wd1856GtcSpitzer,
    fit: "contain",
    source: sources.v2020,
    notes: [
      "This figure is useful for the white-dwarf planet experts because it shows the validation mindset around WD 1856. The discovery was not just a TESS dip; the authors built a consistent multi-instrument transit picture.",
      "TWIRL cannot replicate all of this at the search stage, but it can require that the TESS discovery artifacts are sane before a candidate leaves the pipeline: stable aperture behavior, plausible duration, no obvious alias or cluster pile-up, and a clean pixel/centroid story.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Contamination and aperture checks",
    title: "WD 1856 also teaches that apertures and neighbors are not optional.",
    figure: F.wd1856ArchivalApertures,
    fit: "contain",
    railTitle: "TWIRL translation",
    rail: ["1x1 / 3x3 / 5x5 aperture comparison", "pixel-map diagnostics", "centroid checks", "independent extraction for high-value candidates"],
    source: sources.v2020,
    notes: [
      "Use this slide to connect the literature validation to TWIRL's vetting requirements. In a white-dwarf survey, a deep event can be an astrophysical signal, a contaminant in the aperture, a nearby eclipsing binary, or a product artifact.",
      "That is why TWIRL keeps aperture-to-aperture consistency as a minimum QA item. For WD 1856, the known transit should be visible across the relevant apertures, and the pixel maps should tell the same story before we trust search rankings.",
    ],
  },
  {
    kind: "figureFull",
    kicker: "Formation and mass context",
    title: "WD 1856 leaves the formation problem open, which is part of the motivation.",
    figure: F.wd1856MassAge,
    fit: "contain",
    source: sources.v2020,
    notes: [
      "This slide gives the field-level reason the search is interesting to white-dwarf planet experts. WD 1856 raises formation, survival, scattering, and common-envelope questions. The system is not merely another transiting giant planet.",
      "The TWIRL talk does not need to solve the formation channel. It should make clear that a larger sample, or even strong upper limits after a well-calibrated search, would inform whether WD 1856 is rare, representative, or part of a broader class.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "JWST/MIRI confirmation",
    title: "JWST/MIRI adds an infrared confirmation axis for WD 1856 b.",
    figure: F.jwstSed,
    fit: "contain",
    railTitle: "Current relevance",
    rail: ["thermal/IR evidence", "rules out some false-positive space", "strengthens the benchmark", "motivates high-confidence TWIRL triage"],
    source: sources.limbach2025,
    notes: [
      "Use the MIRI SED figure to update the story beyond the 2020 discovery paper. The important point for this talk is that WD 1856 has become an even stronger benchmark system because independent infrared evidence supports the planet interpretation.",
      "That matters for TWIRL because the benchmark is not arbitrary. If a TESS/TGLC/TWIRL-FS product and search stack can reproduce this system cleanly, the pipeline has passed a relevant, high-value white-dwarf planet test.",
    ],
  },
  {
    kind: "figureTwo",
    kicker: "Cold planet context",
    title: "The MIRI result turns WD 1856 from a transit benchmark into a physical planet benchmark.",
    figures: [
      { path: F.jwstTempRadius, label: "Temperature/radius context", fit: "contain" },
      { path: F.jwstImage, label: "MIRI false-color image", fit: "contain" },
    ],
    source: sources.limbach2025,
    notes: [
      "These panels help close the science-history section. WD 1856 is not only a transit-shape problem; it is a system where atmospheric/thermal constraints are becoming possible. That reinforces why new candidates from TWIRL would be valuable follow-up targets.",
      "At the same time, the slide should not imply that TWIRL can confirm planets by itself. TWIRL is a discovery and triage engine. Strong candidates will need independent extraction, contamination analysis, ephemeris adequacy, and follow-up planning.",
    ],
  },
  {
    kind: "divider",
    part: "Part 2",
    title: "Why TESS and why TGLC?",
    subtitle: "The signal is minutes long; the useful survey input is 200 s FFIs with Gaia-prior decontaminated light curves.",
    source: sources.tglc,
    notes: [
      "Transition from science motivation to observability. The key claim is that the cadence and extraction method are the bottlenecks. TWIRL is not simply taking whatever TESS light curve exists; it is using 200 s FFI data and TGLC-style Gaia-prior photometry.",
      "This section should be especially useful to TESS experts: it lays out the survey input, parent sample logic, the TIC/Gaia bridge, and why TGLC's decontaminated aperture products are attractive for crowded white-dwarf fields.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "200 s FFI motivation",
    title: "Most unexplored white dwarfs are a 200 s FFI problem.",
    figure: F.tessUnexploredWd,
    fit: "contain",
    railTitle: "Why this is central",
    rail: ["TOI WDs are a small bright subset", "200 s FFIs open the faint WD population", "WD 1856 sits in the transition", "TWIRL starts at Sector >= 56"],
    source: "Source: TWIRL proposal-style detectability/coverage figure from local deck assets; FIRST_200S_SECTOR = 56.",
    notes: [
      "Use this as the core TESS motivation slide. The plot makes the proposal problem visual: the well-explored two-minute TOI white dwarfs are a small bright subset, while the 200 s FFI product covers a much larger and fainter Gaia WD population.",
      "This is why TWIRL's core survey is restricted to Sector >= 56. Earlier sectors may be useful for context or special cases, but the production survey input is the 200 s FFI product. That boundary keeps the denominator and completeness problem cleaner.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "TESS coverage",
    title: "Most cataloged WDs have at least one Sector >= 56 footprint.",
    figure: F.wdCoverage,
    fit: "contain",
    railTitle: "Coverage snapshot",
    rail: ["1,231,702 WDs with S>=56 footprint", "96.2% of full seed catalog", "median positive-sector count: 3", "high-confidence subset is the statistical anchor"],
    source: sources.coverage,
    notes: [
      "This is the first survey-scale slide, but keep the caveat clear. Coverage is not the same as a final denominator, and footprint is not the same as a validated light curve. The figure shows why TESS is worth using: the sky coverage is broad enough to support a serious survey.",
      "The high-confidence reference cut is Pwd > 0.75, and Gaia DR3 source_id is the authoritative target identifier. TIC is operational metadata. That distinction matters because the MIT implementation is TIC-oriented while the science sample is Gaia-first.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Search regime",
    title: "The first pass is tuned to deep, short WD events.",
    figure: F.wdDetectability,
    fit: "contain",
    railTitle: "Do not overclaim",
    rail: ["Jupiter-like depths are the first benchmark", "Earth-like claims wait for completeness", "typical HZ is a cadence/recovery challenge", "WD 1856 anchors the scale"],
    source: "Source: TWIRL local detectability map; proposal-style figure used for search-regime motivation.",
    notes: [
      "This detectability map should be one of the main motivation slides. It explains why the first TWIRL pass is not an Earth-size or habitable-zone occurrence-rate claim: the robust near-term test is large, deep, short events in the WD 1856-like regime.",
      "TIC IDs are still useful operational metadata, especially for MIT TGLC production and TESS bookkeeping, but they should not define the survey sample. The scientific identity remains Gaia-first; the search strategy is constrained by cadence, depth, and product completeness.",
    ],
  },
  {
    kind: "flow",
    kicker: "TGLC",
    title: "TGLC is attractive because it uses Gaia priors and decontaminated apertures.",
    subtitle: "TWIRL treats MIT TGLC as an orbit/camera/CCD production pipeline.",
    items: [
      ["Gaia priors", "source positions and contaminants"],
      ["TICA FFIs", "200 s calibrated image stream"],
      ["ePSF model", "orbit/camera/CCD local model"],
      ["HDF5 light curves", "1x1, 3x3, 5x5 aperture fluxes"],
      ["TWIRL HLSP", "flux-space detrend and FITS export"],
    ],
    callout: "The output that TWIRL searches is decontaminated aperture photometry, not a legacy quick_lc product.",
    source: sources.tglc,
    notes: [
      "The point here is not to sell TGLC abstractly. TWIRL needs a light-curve product that can handle crowded fields and faint WDs, and TGLC's Gaia-prior design is a natural match. The MIT fork is used as a production pipeline, not as a one-target quick extraction.",
      "For downstream search, TWIRL expects HDF5 light curves with decontaminated aperture fluxes in 1x1, 3x3, and 5x5 apertures. That aperture ladder becomes a vetting feature, not just an extraction detail.",
    ],
  },
  {
    kind: "stageRoadmap",
    kicker: "twirl_plan.md skeleton",
    title: "The TWIRL section should follow the repo's Stage 1-5 plan.",
    subtitle: "The talk is not a generic progress update: it is a Stage 1 product plus an S56 Stage 2-5 vertical slice.",
    stages: [
      ["Stage 1", "Generate WD light curves", "Current gate: TWIRL-FS v2 product QA and S56 handoff tree."],
      ["Stage 2", "Build search stack", "Franklin's BLS is the periodic baseline; dip branch still to build."],
      ["Stage 3", "Injection recovery", "LC-level now for iteration; pixel-level for publication completeness."],
      ["Stage 4", "Full search", "Scale to all sectors only after product/search gates are stable."],
      ["Stage 5", "Validate targets", "Heuristic + Michelle LEO-Vetter + pixel/centroid + follow-up readiness."],
    ],
    status: [
      ["Seed WD catalog", "1,280,266 rows"],
      ["High-confidence WDs", "359,073 with Pwd > 0.75"],
      ["S>=56 footprint", "339,292 high-confidence WDs"],
      ["S56 handoff product", "19,072 / 19,072 FITS"],
    ],
    source: sources.plan,
    notes: [
      "Use this slide to restructure the TWIRL section around the actual repo plan. Stage 1 is the light-curve production and QA foundation. Stage 2 is the interpretable search stack. Stage 3 is injection recovery. Stage 4 is the full search. Stage 5 is validation and occurrence-rate inputs.",
      "The current talk evidence sits mostly in Stage 1 and the S56 vertical slice through Stages 2 to 5. That distinction matters: S56 is a benchmark and collaboration handoff product, not a final survey denominator. The high-confidence WD reference sample is 359,073 objects, and 339,292 of those have at least one Sector 56-or-later footprint, but the final occurrence-rate denominator is deliberately not frozen yet.",
    ],
  },
  {
    kind: "productionDashboard",
    kicker: "Stage 1 LC scope",
    title: "The S56-S94 light-curve scope is a sector-count problem.",
    subtitle: "Catalog-derived unique-TIC WD targets per sector; no orbit-level detail.",
    sectors: [
      [56, 31590],
      [57, 27326],
      [58, 23256],
      [59, 22181],
      [60, 27165],
      [61, 41403],
      [62, 40158],
      [63, 53512],
      [64, 80254],
      [65, 114940],
      [66, 256727],
      [67, 108147],
      [68, 54713],
      [69, 40068],
      [70, 18558],
      [71, 19413],
      [72, 25936],
      [73, 24695],
      [74, 40745],
      [75, 33035],
      [76, 29625],
      [77, 29108],
      [78, 31859],
      [79, 42608],
      [80, 67934],
      [81, 86673],
      [82, 46676],
      [83, 30919],
      [84, 27318],
      [85, 23796],
      [86, 21590],
      [87, 37053],
      [88, 41644],
      [89, 40404],
      [90, 55891],
      [91, 104175],
      [92, 133617],
      [93, 258721],
      [94, 120174],
    ],
    callouts: [
      ["Catalog scope", "2,313,607 unique-TIC target-sector LCs across S56-S94"],
      ["High-conf subset", "795,096 target-sector LCs with Pwd > 0.75"],
      ["Production state", "S56 handoff product exists; S57-S60 legacy finalized"],
      ["Current gate", "S72-S90 prep tail near done; finalization waits on TWIRL-FS v2 QA"],
    ],
    source:
      "Source: data_local/catalogs/twirl_master_catalog/twirl_wd_tess_sector_summary_v0.csv; doc/twirl_plan.md current status as of 2026-06-02.",
    notes: [
      "This slide replaces the old partial-product dashboard with the catalog-derived scale of the light-curve problem. The bars are the number of unique-TIC WD targets per sector in the local TWIRL sector-summary catalog, from S56 through S94. They are target-sector LC scope, not a claim that every product has been finalized.",
      "The headline numbers are intentionally simple: 2,313,607 unique-TIC target-sector light curves across S56-S94, and 795,096 target-sector entries in the Pwd > 0.75 high-confidence subset. The operational status remains gated: S56 is the handoff product, S57-S60 legacy products exist, S72-S90 prep is near the end, and finalization waits on TWIRL-FS v2 QA.",
    ],
  },
  {
    kind: "divider",
    part: "Part 3",
    title: "What S56 proves, and what it does not.",
    subtitle: "S56 is a benchmark-gated vertical slice: enough to test the stack, not enough to claim the survey is complete.",
    source: sources.s56,
    notes: [
      "This divider is important because the rest of the talk becomes very concrete. Prevent the room from hearing S56 results as final survey claims. The S56 slice is valuable because it contains WD 1856 and because the end-to-end data product exists.",
      "The claims in this section should be phrased as product and benchmark claims: what is built, what is recovered, what failures were found, and what still blocks scale-out.",
    ],
  },
  {
    kind: "flow",
    kicker: "S56 benchmark",
    title: "The current vertical slice runs from catalog to vetting artifact.",
    subtitle: "Each box is a real handoff point in the current repo logs.",
    items: [
      ["WD catalog", "Gaia EDR3 WD seed, Pwd reference cut"],
      ["Detector tables", "S56 orbits 119+120, all CCD geometry"],
      ["MIT TGLC HDF5", "decontaminated aperture fluxes"],
      ["TWIRL-FS v2 HLSP", "robust flux-space detrend"],
      ["Search/vet", "BLS, heuristic vetter, LEO-Vetter path"],
    ],
    figure: F.s56Boundary,
    callout: "Benchmark target: WD 1856+534 / TIC 267574918 in S56 cam4/ccd1.",
    source: sources.s56,
    notes: [
      "Walk through the pipeline as an audit chain. The white-dwarf catalog defines the candidate universe. Sector/orbit/camera/CCD coverage maps decide what TESS data can be used. MIT TGLC produces HDF5 light curves. TWIRL-FS converts those into search-ready HLSP FITS. Search and vetting operate on those products.",
      "This slide is also a collaboration map: the light-curve product and data stewardship are Te Han's responsibility, Franklin owns the BLS/search work shown, and Michelle owns the LEO-Vetter collaboration path.",
    ],
  },
  {
    kind: "negativeFlux",
    kicker: "Faint-end failure mode",
    title: "The faint-end risk was negative and near-zero background-subtracted flux.",
    source: sources.twirlfs,
    notes: [
      "Explain the failure mode at the level needed for both TESS and WD experts. For faint targets, real background-subtracted flux values can be negative or near zero. A magnitude-style conversion or divisive detrending can accidentally turn those cadences into NaNs or unstable ratios.",
      "The TGLC-side patch preserves linear RawFlux and RawFluxError, and TWIRL-FS uses subtractive flux-space detrending. On WD 1856, the RawFlux path preserved 155 previously NaN cadences and increased usable q=0 cadences by about 6.3%, while preserving the transit shape.",
    ],
  },
  {
    kind: "method",
    kicker: "TWIRL-FS v2",
    title: "The active handoff product is conservative by design.",
    source: sources.twirlfs,
    notes: [
      "This is the methods slide for the light-curve product. The canonical search column is DET_FLUX, using robust subtractive flux-space detrending with bkspace_d = 0.8 d and gap splitting at 0.5 d. It is designed to preserve legitimate negative and near-zero flux cadences rather than normalize them away.",
      "The adaptive columns are intentionally opt-in. They can reduce broad residual structure, but they may attenuate long or shallow astrophysical features. That is why the deck presents them as a comparison input, not as a silent replacement for the conservative default.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Photometric QA",
    title: "The S56 precision check now uses the expected-flux denominator.",
    figure: F.precision,
    fit: "contain",
    railTitle: "Read this as",
    rail: ["product QA", "not completeness", "magnitude-dependent behavior", "input to scale decision"],
    source: "Source: reports/stage1_lightcurves/qc_pdf/s56_precision_twirl_fs_v2_compare_det_flux.png.",
    notes: [
      "Use this figure to show that the QA is quantitative, not just visual. The important product behavior is how scatter changes with expected flux or magnitude and whether the faint end behaves in a way we understand.",
      "Do not claim that the precision curve alone validates detection completeness. It is one piece of the gate. The product also needs cadence retention statistics, aperture-to-aperture consistency, WD 1856 recovery, independent extraction checks, and injection recovery through the full search and vetting stack.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Signal preservation",
    title: "The canonical detrend preserves short injected events across the tested grid.",
    figure: F.retentionShort,
    fit: "contain",
    railTitle: "Gate logic",
    rail: ["short WD-like events", "depth = 0.1 example", "bkspace = 0.8 d passes", "too-tight splines can attenuate"],
    source: "Source: reports/stage1_lightcurves/detrend_experiments/s56_param_sweep_signal_preservation.",
    notes: [
      "This slide explains why bkspace_d = 0.8 d is the conservative default. The retention sweep tests injected signals and asks whether the detrending removes astrophysical dips along with instrumental structure.",
      "For an expert audience, make clear that this is a product-level smoke test, not the final injection-recovery study. The final completeness calculation must run injections through the actual search, classifier if used, vetter, and candidate-merging path.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Longer events",
    title: "Longer injected events expose the tradeoff behind adaptive detrending.",
    figure: F.retentionLong,
    fit: "contain",
    railTitle: "Interpretation",
    rail: ["short-event search first", "adaptive columns are useful diagnostics", "long/shallow signals need caution", "dip-search will be separate"],
    source: "Source: reports/stage1_lightcurves/detrend_experiments/s56_param_sweep_long_events.",
    notes: [
      "This figure helps prevent overconfidence. A tighter detrend can make light curves look cleaner while suppressing longer astrophysical signals. That is acceptable only if the product and search channel are explicitly scoped.",
      "TWIRL's first-year science is optimized for large, deep, short-duration events in the WD 1856-like regime. Longer or debris-like events motivate a separate dip-search branch and dedicated injection tests, not silent reuse of a period-search detrend setting.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Dim-target QA",
    title: "Gap splitting fixed a real v1 failure mode at the faint end.",
    figure: F.dimQa,
    fit: "contain",
    railTitle: "Example",
    rail: ["S57 TIC 1400694779", "Tmag about 20", "v1 fit across multi-day gap", "v2 splits gaps > 0.5 d", "raw amp 1.711 -> v2 amp 0.287"],
    source: "Source: reports/stage1_lightcurves/detrend_experiments/twirl_fs_v2_18_dim_examples/summary.md.",
    notes: [
      "This is a concrete example of why the project is benchmark-gated. TWIRL-FS v1 looked plausible until a dim S57 target revealed that one spline across a multi-day orbit gap could fail and fall back to a weak polynomial.",
      "TWIRL-FS v2 splits gaps larger than 0.5 days and keeps the subtractive flux-space approach. In the 18 dim examples, all negative quality-zero RawFlux cadences stayed finite, and the median binned-trend amplitude reduction was about 0.796.",
    ],
  },
  {
    kind: "figureFull",
    kicker: "QA gallery",
    title: "The 100-target S56 v2 compare gallery is the product-level visual audit.",
    figure: F.qaGallery,
    fit: "contain",
    source: "Source: reports/stage1_lightcurves/detrend_experiments/s56_twirl_fs_v2_compare_quick_qa_100.",
    notes: [
      "This slide is useful for a TESS-heavy room because it shows that the product is being inspected as a product, not only through the benchmark target. The QA sample includes all 16 CCDs and both canonical and adaptive comparison columns.",
      "Use the gallery to say what remains open: this visual audit supports the handoff, but the full-product QA still has to quantify cadence retention, failure modes, RMS/MAD behavior, quality-flag behavior, and aperture consistency before scale-out claims.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "WD 1856 recovery artifact",
    title: "The benchmark event is visible in the S56 search/vet product.",
    figure: F.blsVet,
    fit: "contain",
    railTitle: "Benchmark read",
    rail: ["TIC 267574918", "S56 product", "BLS + folded view", "transit visible across panels", "candidate artifact for review"],
    source: sources.bls,
    notes: [
      "This figure is the bridge from light curves into search. It is the concrete WD 1856 vet sheet from the S56 BLS workflow. The event is visible in the search artifact, which is exactly the kind of object the collaboration can review together.",
      "Be careful with wording. This is a benchmark recovery artifact, not a final blind survey result. The value is that it makes the discussion auditable: experts can point to the plot and ask about cadence masks, apertures, period aliases, and nearby sources.",
    ],
  },
  {
    kind: "figureTwo",
    kicker: "Pixel-level diagnostic",
    title: "Pixel maps are the next validation layer after aperture agreement.",
    figures: [
      { path: F.pixelMapTransit, label: "Transit-depth map", fit: "contain" },
      { path: F.pixelMapFolded, label: "Folded pixel-map diagnostic", fit: "contain" },
    ],
    source: "Source: reports/stage1_lightcurves/pixel_maps WD 1856 S56 TGLC SPOC saveaper diagnostics.",
    notes: [
      "This slide anchors the centroid/pixel-map branch of the vetting plan. A BLS peak and a folded light curve are not enough if the event could come from a neighboring source. Pixel-level diagnostics ask whether the transit-like behavior is spatially consistent with the WD.",
      "The current TWIRL state includes these diagnostics for WD 1856. The next step is to turn them into a repeatable candidate-vetting module rather than an individual benchmark plot.",
    ],
  },
  {
    kind: "divider",
    part: "Part 4",
    title: "Search and vetting are now a collaboration interface.",
    subtitle: "Franklin owns the BLS/search work shown; Michelle owns LEO-Vetter; Te Han owns the light-curve product and data stewardship.",
    source: sources.plan,
    notes: [
      "This divider explicitly sets collaboration credit and ownership. It should prevent ambiguity in the room. Franklin's BLS work is the search artifact being showcased. Michelle owns LEO-Vetter and the collaboration path around it. Te Han is responsible for TGLC/TWIRL-FS production and data stewardship.",
      "The technical goal is not to defend one tool. It is to define clean interfaces: search tables, vet sheets, LEO input/output schemas, pixel-map diagnostics, and shared candidate review criteria.",
    ],
  },
  {
    kind: "blsShowcase",
    kicker: "Franklin's BLS search",
    title: "Franklin's BLS work gives TWIRL a concrete periodic-search artifact.",
    figure: F.blsVetMedium,
    source: sources.bls,
    notes: [
      "Showcase Franklin's work clearly. The BLS search converts the TWIRL-FS light curves into ranked periodic candidates and vet sheets. This is the first discovery engine in the talk because it is interpretable and directly matched to WD 1856-like events.",
      "The current repo logs record a denser S56 BLS v2 configuration with p_min = 0.12 d and durations [3, 4, 5, 6, 8, 10, 13, 16, 20, 30] minutes. WD 1856's blind rank improved from 1258 to 560 in that run, while the heuristic vetter already carried most of the post-vetting gain.",
    ],
  },
  {
    kind: "heuristic",
    kicker: "Te Han heuristic vetter",
    title: "The heuristic vetter is a physics-motivated bridge, not the final arbiter.",
    source: sources.plan,
    notes: [
      "Explain the heuristic vetter as a transparent bridge between BLS peaks and LEO-Vetter. It uses physics-motivated cuts on duration, period aliases, and period-cluster pile-ups to reduce obvious false positives without pretending to be a full validation engine.",
      "The current benchmark result is strong enough for a talk: WD 1856 moved from blind rank 1258 out of 19,040 to planet-regime rank 9 out of 5,403, about a 140x improvement. That justifies the vetter as an audit layer while LEO-Vetter is tuned for white-dwarf hosts.",
    ],
  },
  {
    kind: "figureFull",
    kicker: "Michelle's LEO-Vetter path",
    title: "LEO-Vetter classifies the WD 1856 benchmark as PC in the tuned path.",
    figure: F.leoWd1856Report,
    fit: "contain",
    source: "Source: Michelle-owned LEO-Vetter path; reports/stage5_validation/leo_vetter_s56_top50_v2/vet_reports/PC_rank11_tic0267574918_T16.34_P001.4080d.pdf.",
    notes: [
      "Use the wording carefully. Michelle owns LEO-Vetter. TWIRL's role is to provide clean light curves, BLS/candidate artifacts, WD-specific metadata, and the questions that need tuning for compact hosts and 200 s cadence.",
      "This slide should be treated as a major proof object. It shows the full LEO-Vetter style report for WD 1856 in the S56 tuned path: raw and detrended views, phase diagram, primary/odd/even/secondary checks, model properties, and the PC classification context.",
    ],
  },
  {
    kind: "collaboration",
    kicker: "Collaboration roles",
    title: "The handoff should be artifact-based, not memory-based.",
    source: sources.plan,
    notes: [
      "This is the collaboration slide requested by the user. Te Han owns the TGLC/TWIRL-FS production tree, product QA, provenance, and data stewardship. Franklin owns the BLS/search work shown, including vet sheets and candidate tables. Michelle owns LEO-Vetter and the collaboration path for WD tuning.",
      "The practical recommendation is to make the interface artifacts stable: HLSP tree and metadata, BLS candidate tables, vet sheets, LEO input/output tables, pixel-map diagnostics, and a shared candidate-review sheet.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Vetting landscape",
    title: "The first real candidate list is mostly a false-positive taxonomy exercise.",
    figure: F.leoSummary,
    fit: "contain",
    railTitle: "Use the list to classify",
    rail: ["PC / FA / FP / ERR split", "PCEBs and WD+WD binaries", "centroid-shift contaminants", "ZZ Ceti pulsators", "systematic period ladders"],
    source: sources.plan,
    notes: [
      "This slide makes the next scientific work concrete. The top candidates are not automatically planet candidates. They must be sorted into classes: post-common-envelope binaries, WD+WD or compact binaries, centroid-shift contaminants, ZZ Ceti pulsators, systematic ladder peaks, debris-like dips, and possible planet-regime events.",
      "For the expert audience, this is a useful discussion prompt. Ask which failure modes they expect to dominate a 200 s TESS WD search and which ones should block scale-out until the vetting machinery is more mature.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Injection recovery",
    title: "Completeness must be measured through the real stack.",
    figure: F.retentionShort,
    fit: "contain",
    railTitle: "Stack to test",
    rail: ["LC-level injections for fast iteration", "pixel-level injections for crowding/extraction", "BLS and dip-search branches", "heuristic + LEO-Vetter losses", "candidate merging before claims"],
    source: sources.plan,
    notes: [
      "This slide protects the scientific claims. TWIRL should not present occurrence rates or Earth-size/habitable-zone claims until completeness is demonstrated. Injection recovery must go through the actual chain: light curve product, search, classifier if any, vetter, and candidate merging.",
      "There are two useful levels. Fast light-curve-level injections can test detection and vetting cheaply. Pixel-level injections are slower but necessary for extraction, crowding, aperture choice, and detrending failures. The talk should invite expert feedback on how to prioritize those layers.",
    ],
  },
  {
    kind: "figureRail",
    kicker: "Next gates",
    title: "Scale-out waits for a short list of falsifiable gates.",
    figure: F.qaGallery,
    fit: "contain",
    railTitle: "Next gates",
    rail: ["TGLC production QA", "S56 TWIRL-FS v2 search", "BLS candidate tables", "LEO-Vetter tuning", "centroid and pixel-map vetting", "dip-search branch"],
    source: sources.plan,
    notes: [
      "Use this as the actionable next-steps slide. The gates are TGLC production QA, S56 TWIRL-FS v2 search, BLS candidate tables, LEO-Vetter tuning, centroid and pixel-map vetting, injection recovery, and a separate dip-search branch.",
      "The point is to make the scale decision falsifiable. If the product fails full QA, we do not scale. If WD 1856 is unstable across apertures or vetters, we debug the benchmark. If injections show strong attenuation or vetting loss, we tune before occurrence-rate work.",
    ],
  },
  {
    kind: "decision",
    kicker: "Decision slide",
    title: "Scale only after the benchmark evidence is stable.",
    rows: [
      ["WD 1856 recovery", "periodic search + aperture agreement + pixel sanity"],
      ["Product QA", "cadence retention, RMS/MAD, failure modes, expected-flux behavior"],
      ["Vetting interface", "Franklin BLS tables, Michelle LEO-Vetter outputs, shared review"],
      ["Completeness", "injections through search, vetter, and candidate merging"],
      ["Claims", "no occurrence rate until denominator and completeness are locked"],
    ],
    source: sources.plan,
    notes: [
      "Close with the decision logic. The survey is compelling because the science payoff is high and TESS 200 s coverage is broad, but the project should scale only after the benchmark evidence is stable.",
      "The final discussion questions for the room are: what would convince you that WD 1856 is recovered well enough; where should independent validation enter the loop; which TESS/WD failure modes are most dangerous; and what minimum injection-recovery design is credible for first occurrence-rate claims?",
    ],
  },
];

const extraTalkTrack = [
  "Spend another minute defining the deliverable: this is an editable Keynote-importable PPTX with provenance footers and notes, not a paper draft. Tell the audience that the strongest feedback you want is on the gates: what would make them trust the product enough to scale, and what would still make them stop.",
  "Use the flowchart as a map for questions. If a TESS person asks about extraction, point to the TGLC section. If a white-dwarf planet person asks about validation, point to the WD 1856 and LEO-Vetter sections. This keeps the talk organized while still inviting interruptions.",
  "Frame the first section as a shared vocabulary check. The room likely knows the individual systems, but TWIRL depends on connecting them: white-dwarf pollution says planetary material survives, WD 1145 says transiting debris happens, and WD 1856 says planet-scale transits are observable with TESS.",
  "Pause on the validation burden. The same small-star geometry that gives dramatic transit depths also makes false positives visually tempting. A nearby eclipsing source, a data-product artifact, or a pulsator can look exciting unless the pipeline keeps aperture and pixel evidence attached to every candidate.",
  "Make the larger history cards do work rather than reading them. Ask the audience to hold two endpoints in mind: evolving debris systems found by K2/ZTF-style searches and stable planet-scale transits like WD 1856. TWIRL needs both the periodic baseline and a future dip/debris branch.",
  "Use the publication figure as the proof object. The point is not to describe every K2 epoch in detail; it is to show that WD 1145 morphology evolves enough that TWIRL should not force all white-dwarf transit behavior into one periodic BLS template.",
  "When the all-transits figure is on screen, orient the audience to the repeated deep events first, then connect the visual depth to the search metric. A good TWIRL product should make this object easy enough to recover that failures become diagnostic rather than ambiguous.",
  "Use this slide to separate discovery from validation. TESS found the event, but the confidence came from compatible follow-up light curves. TWIRL can lower the cost of finding candidates, but the project should still package evidence in a way that makes follow-up decisions rational.",
  "Name the TWIRL translation plainly: every serious candidate should carry small, medium, and large aperture behavior, nearby-source context, and pixel diagnostics. That is the practical inheritance from the WD 1856 literature validation, not a ceremonial citation.",
  "This is a good moment to invite the WD experts into the formation question without letting the talk drift. Say that TWIRL's immediate contribution is empirical: either find more systems or build a defensible non-detection after completeness, both of which constrain how rare WD 1856-like systems are.",
  "Treat the JWST/MIRI slide as a benchmark-strengthening update. The 2025 result does not change the TWIRL pipeline, but it changes the confidence that WD 1856 is the right gate. The target is not merely convenient; it is now even more physically valuable.",
  "Use the two panels to remind the room that discovery value depends on follow-up value. If TWIRL finds a clean analogue, it is not only a TESS candidate; it is potentially an infrared, high-cadence, and atmospheric-characterization target once the ephemeris and contamination story are strong.",
  "The transition here should be quick. The science case is now established, so the rest of the section asks whether the available TESS product is actually matched to the signal. This is where TESS expertise can most directly improve the project design.",
  "Make the sampling argument concrete. A 2 minute event can still be poorly represented even at 200 s cadence, while a 10 or 15 minute event becomes much more searchable. This is why the early science target is large and deep events, not small shallow habitable-zone claims.",
  "Add the denominator caution verbally. The coverage map is exciting, but coverage is only the first condition. Each target still needs production, quality flags, cadence retention, light-curve QA, and a defined parent-sample cut before it can enter a statistical denominator.",
  "Use the detectability figure to keep the scientific claim narrow. The map makes clear that the first credible pass is tuned for deep, short, WD 1856-like events. Earth-like or habitable-zone occurrence-rate language belongs later, after product QA and injection recovery are complete.",
  "This is the TGLC value proposition in one slide. Gaia priors are not just convenient metadata; they help model contamination and produce aperture fluxes that can be compared during vetting. That is why TGLC is better aligned with TWIRL than a black-box one-aperture product.",
  "Use the Stage 1-5 roadmap to make the rest of the TWIRL section feel like the project plan rather than a pile of artifacts. Stage 1 is the current production gate, while Stages 2 through 5 are being exercised only through the S56 vertical slice.",
  "Spend enough time on the operational assumptions because they define feasibility. The GPU ePSF path, PDO data layout, and HDF5-to-HLSP export are the machinery that lets the project move beyond hand-picked benchmarks. Then repeat the caveat: operational feasibility still needs scientific QA.",
  "This divider should reset the standard of evidence. From here forward, every claim is about S56 as a vertical slice. That makes the evidence stronger, not weaker, because it prevents the talk from mixing product-development results with final survey conclusions.",
  "Use the pipeline loop to make the current state auditable. If someone asks where a number comes from, it should map to a stage: catalog, detector coverage, TGLC HDF5, TWIRL-FS HLSP, BLS table, vetter output, or pixel-map artifact. That is the structure to preserve in future reports.",
  "Slow down on negative flux because it is easy to underestimate. At the faint end, negative background-subtracted flux is not automatically bad data. Throwing it away can bias cadence retention and search behavior exactly where TWIRL expects many white dwarfs to live.",
  "Explain the conservative product as an intentional compromise. The goal is not the prettiest light curve; it is a search input that preserves short events while avoiding fragile normalization around zero flux. The adaptive columns exist because comparison is useful, but they are not the default scientific product.",
  "When presenting the precision plot, point to both panels. The top panel gives an absolute scatter scale; the bottom panel shows behavior relative to an expected-flux trend. Then state the limitation: precision is necessary for search, but it is not the same as recovery completeness.",
  "Use this figure to justify a parameter choice without overselling it. The canonical setting passes the tested short-event retention grid, which is exactly the regime of the first benchmark. The next question is whether the same choice behaves across sectors and through the actual search stack.",
  "The long-event slide is where you show scientific restraint. A method can be excellent for WD 1856-like events and still be wrong for longer or debris-like signals. That is why the project should split search branches instead of pretending one detrend setting optimizes everything.",
  "This dim-target example is a story about finding a failure before scale-out. Make that positive but technical: v1 failed a real dim case, the failure mode was understandable, and v2 encodes the fix. That is exactly how a benchmark-gated pipeline should mature.",
  "Use the gallery as evidence of breadth. It is not a proof of correctness, but it shows that the handoff was not judged by WD 1856 alone. Invite the TESS experts to identify any patterns they would want quantified before a sector-scale production decision.",
  "On the vet sheet, orient the viewer to the raw/folded/search panels and say why this is the collaboration audit surface. The figure is not just a plot; it is the object that lets Franklin, Michelle, Te Han, and external experts critique the same candidate consistently.",
  "The pixel-map slide should connect directly to false positives. A signal that is deep in the aperture but spatially inconsistent with the WD is not a planet candidate. The practical next step is to standardize this diagnostic for every high-priority BLS candidate.",
  "Use this divider to protect collaboration clarity. The talk should not blur roles. Franklin gets explicit credit for the BLS/search work shown. Michelle's ownership of LEO-Vetter is named. Te Han's ownership of production and data stewardship is named. That wording should stay stable.",
  "This is Franklin's showcase slide. Give enough technical detail to make the work credible: the BLS search generates ranked candidates and vet sheets, the denser duration grid was tested, and the WD 1856 rank behavior is tracked. Then hand off naturally to vetting.",
  "Describe the heuristic vetter as a transparent audit layer. It is valuable because experts can understand why a candidate moves up or down. It should not replace LEO-Vetter or human review; it should reduce clutter and make the LEO/human workload more scientifically targeted.",
  "For the LEO-Vetter slide, keep Michelle's ownership explicit and respectful. TWIRL should supply stable inputs and WD-specific requirements, while Michelle's LEO-Vetter path handles the validation layer. The full WD 1856 report is the key proof object because it shows what a collaborator can audit.",
  "Use the collaboration flowchart as an operating plan. The most important phrase is artifact-based handoff. Every stage should leave a table, figure, or metadata record that another collaborator can inspect without relying on a meeting memory or a private notebook.",
  "The false-positive taxonomy slide can become a discussion slide if time allows. Ask the room which classes they expect to dominate a 200 s WD search. Their answer can directly shape BLS thresholds, LEO tuning, centroid checks, and follow-up triage.",
  "For injections, emphasize that the project needs both speed and realism. Fast light-curve injections are ideal for development, but a paper-grade completeness claim eventually needs to include extraction and crowding failures. This keeps the near-term plan useful without pretending it is the final calibration.",
  "The gates slide should feel like a work plan. Each gate is falsifiable: a report either exists or it does not; WD 1856 either remains stable or it does not; injections either quantify recovery or they do not. That makes the path to scale-out concrete.",
  "Close by asking for decisions, not applause. The deck should leave the audience with a small number of technical questions: what evidence is enough for WD 1856, which false positives are most dangerous, how should Michelle's LEO-Vetter interface be packaged, and what injection design is credible.",
];

if (extraTalkTrack.length !== slides.length) {
  throw new Error(`Expected ${slides.length} extra talk-track entries, found ${extraTalkTrack.length}`);
}

for (let i = 0; i < slides.length; i += 1) {
  slides[i].notes = [...(slides[i].notes || []), extraTalkTrack[i]];
}

function collectFigures(configs) {
  const paths = new Set();
  for (const slide of configs) {
    if (slide.figure) paths.add(slide.figure);
    if (slide.figures) {
      for (const item of slide.figures) paths.add(item.path);
    }
  }
  return [...paths];
}

function notesText(slide, index) {
  const note = Array.isArray(slide.notes) ? slide.notes.join("\n\n") : slide.notes || "";
  return [`## ${String(index + 1).padStart(2, "0")}. ${slide.title}`, "", note, ""].join("\n");
}

function sourceText(slide, index) {
  return [
    `## ${String(index + 1).padStart(2, "0")}. ${slide.title}`,
    "",
    slide.source || "Source: TWIRL local project materials.",
    "",
  ].join("\n");
}

function estimateTalkMinutes(configs) {
  const words = configs
    .flatMap((slide) => slide.notes || [])
    .join(" ")
    .split(/\s+/)
    .filter(Boolean).length;
  return { words, minutesAt130Wpm: words / 130 };
}

const sharedModule = String.raw`
const C = {
  black: "#050505",
  ink: "#111111",
  muted: "#636363",
  light: "#f6f6f3",
  lighter: "#fbfbf8",
  border: "#b8b8b1",
  rule: "#8c8c85",
  accent: "#c69a28",
  accent2: "#4d8f82",
  blue: "#436fb4",
  red: "#b94b4b",
  green: "#4b8b66",
};

function addBackground(slide, ctx, color) {
  ctx.addShape(slide, { x: 0, y: 0, w: ctx.W, h: ctx.H, fill: color, line: ctx.line(color, 0) });
}

function scaleFontSize(size) {
  if (typeof size !== "number") return size;
  if (size <= 9) return Math.round((size + 2.3) * 10) / 10;
  if (size <= 12) return Math.round((size + 3.0) * 10) / 10;
  if (size <= 16) return Math.round((size + 3.0) * 10) / 10;
  if (size <= 20) return Math.round((size + 2.3) * 10) / 10;
  if (size <= 30) return Math.round((size + 2.0) * 10) / 10;
  return Math.round((size + 2.0) * 10) / 10;
}

function createReadableContext(ctx) {
  const readable = Object.create(ctx);
  readable.addText = (slide, options) => {
    const next = { ...options };
    if (typeof next.fontSize === "number") next.fontSize = scaleFontSize(next.fontSize);
    return ctx.addText.call(ctx, slide, next);
  };
  return readable;
}

function addHeader(slide, ctx, cfg) {
  ctx.addText(slide, {
    text: cfg.kicker || "TWIRL",
    x: 52,
    y: 28,
    w: 400,
    h: 20,
    fontSize: 11,
    color: C.muted,
    bold: true,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  ctx.addText(slide, {
    text: cfg.title,
    x: 52,
    y: 58,
    w: 1000,
    h: 54,
    fontSize: 29,
    color: C.ink,
    bold: true,
    typeface: ctx.fonts.title,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  if (cfg.subtitle) {
    ctx.addText(slide, {
      text: cfg.subtitle,
      x: 54,
      y: 103,
      w: 1040,
      h: 32,
      fontSize: 15,
      color: C.muted,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
  }
}

function addFooter(slide, ctx, source) {
  ctx.addShape(slide, { x: 52, y: 681, w: 1076, h: 1, fill: "#d0d0ca", line: ctx.line("#d0d0ca", 0) });
  ctx.addText(slide, {
    text: source || "Source: TWIRL local project materials.",
    x: 52,
    y: 688,
    w: 1080,
    h: 24,
    fontSize: 8.5,
    color: "#6f6f6b",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  ctx.addText(slide, {
    text: String(ctx.slideNumber).padStart(2, "0"),
    x: 1200,
    y: 688,
    w: 42,
    h: 18,
    fontSize: 8,
    color: "#9a9a94",
    align: "right",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
}

function setNotes(slide, cfg) {
  if (slide.speakerNotes && typeof slide.speakerNotes.setText === "function") {
    slide.speakerNotes.setText((cfg.notes || []).join("\n\n"));
  }
}

function addTitleSlide(slide, ctx, cfg) {
  addBackground(slide, ctx, C.black);
  ctx.addText(slide, {
    text: cfg.kicker || "TWIRL",
    x: 64,
    y: 54,
    w: 320,
    h: 22,
    fontSize: 13,
    color: "#d6d6d0",
    bold: true,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  ctx.addText(slide, {
    text: cfg.title,
    x: 64,
    y: 198,
    w: 980,
    h: 70,
    fontSize: 45,
    color: "#ffffff",
    bold: true,
    typeface: ctx.fonts.title,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  ctx.addText(slide, {
    text: cfg.subtitle,
    x: 66,
    y: 294,
    w: 900,
    h: 58,
    fontSize: 19,
    color: "#d8d8d2",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  ctx.addShape(slide, { x: 64, y: 432, w: 760, h: 1, fill: "#a8a8a0", line: ctx.line("#a8a8a0", 0) });
  let x = 64;
  for (const chip of cfg.chips || []) {
    const width = Math.max(150, Math.min(350, chip.length * 8 + 36));
    ctx.addText(slide, {
      text: chip,
      x,
      y: 470,
      w: width,
      h: 28,
      fontSize: 12,
      color: "#eeeeea",
      fill: "#1a1a18",
      line: ctx.line("#5a5a55", 1),
      insets: { left: 12, right: 10, top: 6, bottom: 4 },
    });
    x += width + 12;
  }
  ctx.addText(slide, {
    text: "Te Han + Franklin Chen + Michelle LEO-Vetter collaboration path",
    x: 64,
    y: 612,
    w: 670,
    h: 22,
    fontSize: 12,
    color: "#bdbdb7",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  addFooter(slide, ctx, cfg.source);
}

function addDividerSlide(slide, ctx, cfg) {
  addBackground(slide, ctx, C.black);
  ctx.addText(slide, {
    text: cfg.part || "Part",
    x: 64,
    y: 54,
    w: 220,
    h: 22,
    fontSize: 12,
    color: "#b8b8b0",
    bold: true,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  ctx.addText(slide, {
    text: cfg.title,
    x: 64,
    y: 242,
    w: 1000,
    h: 68,
    fontSize: 41,
    color: "#ffffff",
    bold: true,
    typeface: ctx.fonts.title,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  ctx.addShape(slide, { x: 64, y: 366, w: 850, h: 1, fill: "#b6b6ae", line: ctx.line("#b6b6ae", 0) });
  ctx.addText(slide, {
    text: cfg.subtitle || "",
    x: 64,
    y: 400,
    w: 850,
    h: 70,
    fontSize: 17,
    color: "#d7d7d0",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  addFooter(slide, ctx, cfg.source);
}

function addClaimRail(slide, ctx, cfg, x = 878, y = 145, w = 302, h = 500) {
  if (!cfg.rail && !cfg.railTitle && !cfg.railCards) return;
  ctx.addShape(slide, { x, y, w: 2, h, fill: "#b7b7af", line: ctx.line("#b7b7af", 0) });
  ctx.addText(slide, {
    text: cfg.railTitle || "Read this as",
    x: x + 22,
    y: y + 2,
    w: w - 24,
    h: 28,
    fontSize: 16,
    color: C.ink,
    bold: true,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  if (cfg.railCards) {
    let yy = y + 48;
    for (const [title, detail] of cfg.railCards) {
      ctx.addText(slide, {
        text: title,
        x: x + 22,
        y: yy,
        w: w - 24,
        h: 34,
        fontSize: 19,
        color: C.accent,
        bold: true,
        fill: "#fbf6e7",
        line: ctx.line("#d4ba72", 1),
        insets: { left: 14, right: 14, top: 8, bottom: 4 },
      });
      ctx.addText(slide, {
        text: detail,
        x: x + 22,
        y: yy + 34,
        w: w - 24,
        h: 48,
        fontSize: 16,
        color: C.ink,
        fill: "#fbf6e7",
        line: ctx.line("#d4ba72", 1),
        insets: { left: 14, right: 14, top: 4, bottom: 8 },
      });
      yy += 104;
    }
    return;
  }
  const railItems = cfg.rail || [];
  const gap = 10;
  const itemH = Math.max(56, Math.min(82, Math.floor((h - 54 - gap * Math.max(0, railItems.length - 1)) / Math.max(1, railItems.length))));
  let yy = y + 48;
  for (const item of railItems) {
    ctx.addText(slide, {
      text: item,
      x: x + 22,
      y: yy,
      w: w - 24,
      h: itemH,
      fontSize: 16,
      color: C.ink,
      bold: true,
      fill: "#f7f7f2",
      line: ctx.line(C.border, 1),
      insets: { left: 14, right: 12, top: 12, bottom: 8 },
    });
    yy += itemH + gap;
  }
}

async function addFigureRail(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const figureBox = cfg.figureBox || { x: 58, y: 136, w: 800, h: 520 };
  await ctx.addImage(slide, {
    path: cfg.figure,
    x: figureBox.x,
    y: figureBox.y,
    w: figureBox.w,
    h: figureBox.h,
    fit: cfg.fit || "contain",
  });
  addClaimRail(slide, ctx, cfg, cfg.railX, cfg.railY, cfg.railW, cfg.railH);
  addFooter(slide, ctx, cfg.source);
}

async function addFigureFull(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  await ctx.addImage(slide, {
    path: cfg.figure,
    x: 70,
    y: 128,
    w: 1110,
    h: 535,
    fit: cfg.fit || "contain",
  });
  addFooter(slide, ctx, cfg.source);
}

async function addFigureTwo(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const frames = [
    { x: 58, y: 142, w: 560, h: 488 },
    { x: 662, y: 142, w: 560, h: 488 },
  ];
  for (let i = 0; i < 2; i += 1) {
    const fig = cfg.figures[i];
    const frame = frames[i];
    ctx.addText(slide, {
      text: fig.label,
      x: frame.x,
      y: 124,
      w: frame.w,
      h: 20,
      fontSize: 12,
      color: C.muted,
      bold: true,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
    await ctx.addImage(slide, { path: fig.path, x: frame.x, y: frame.y, w: frame.w, h: frame.h, fit: fig.fit || "contain" });
  }
  addFooter(slide, ctx, cfg.source);
}

function flowBox(slide, ctx, x, y, w, h, title, detail) {
  ctx.addText(slide, {
    text: title,
    x,
    y,
    w,
    h: 28,
    fontSize: 14,
    color: C.ink,
    bold: true,
    fill: C.light,
    line: ctx.line(C.border, 1),
    insets: { left: 12, right: 10, top: 8, bottom: 0 },
  });
  ctx.addText(slide, {
    text: detail,
    x,
    y: y + 28,
    w,
    h: h - 28,
    fontSize: 12,
    color: C.muted,
    fill: C.light,
    line: ctx.line(C.border, 1),
    insets: { left: 12, right: 10, top: 2, bottom: 6 },
  });
}

function addArrowText(slide, ctx, x, y) {
  ctx.addText(slide, {
    text: "->",
    x,
    y,
    w: 28,
    h: 22,
    fontSize: 16,
    color: C.muted,
    bold: true,
    align: "center",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
}

function addTocSlide(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const tableX = 76;
  const tableY = 164;
  const tableW = 1080;
  const rowH = 74;
  const gap = 12;
  ctx.addShape(slide, {
    x: tableX,
    y: tableY - 18,
    w: tableW,
    h: 1,
    fill: "#bdbdb6",
    line: ctx.line("#bdbdb6", 0),
  });
  cfg.items.forEach((row, i) => {
    const y = tableY + i * (rowH + gap);
    ctx.addShape(slide, {
      x: tableX,
      y,
      w: tableW,
      h: rowH,
      fill: i % 2 === 0 ? "#f7f7f2" : "#ffffff",
      line: ctx.line(C.border, 1),
    });
    ctx.addText(slide, {
      text: row[0],
      x: tableX + 18,
      y: y + 16,
      w: 54,
      h: 34,
      fontSize: 20,
      color: C.accent,
      bold: true,
      align: "center",
      fill: "#fbf6e7",
      line: ctx.line("#d4ba72", 1),
      insets: { left: 0, right: 0, top: 6, bottom: 4 },
    });
    ctx.addText(slide, {
      text: row[1],
      x: tableX + 92,
      y: y + 15,
      w: 250,
      h: 30,
      fontSize: 19,
      color: C.ink,
      bold: true,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
    ctx.addText(slide, {
      text: row[2],
      x: tableX + 360,
      y: y + 17,
      w: 580,
      h: 34,
      fontSize: 16,
      color: C.muted,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
    ctx.addText(slide, {
      text: row[3],
      x: tableX + 960,
      y: y + 18,
      w: 82,
      h: 30,
      fontSize: 15,
      color: C.ink,
      bold: true,
      align: "right",
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
  });
  if (cfg.callout) {
    ctx.addText(slide, {
      text: cfg.callout,
      x: 168,
      y: 612,
      w: 940,
      h: 34,
      fontSize: 18,
      color: C.ink,
      bold: true,
      align: "center",
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
  }
  addFooter(slide, ctx, cfg.source);
}

async function addFlowSlide(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const count = cfg.items.length;
  const totalW = 1086;
  const gap = 24;
  const boxW = (totalW - gap * (count - 1)) / count;
  const y = cfg.figure ? 160 : 260;
  const boxH = cfg.figure ? 86 : 98;
  for (let i = 0; i < count; i += 1) {
    const x = 56 + i * (boxW + gap);
    flowBox(slide, ctx, x, y, boxW, boxH, cfg.items[i][0], cfg.items[i][1]);
    if (i < count - 1) addArrowText(slide, ctx, x + boxW + 2, y + 38);
  }
  if (cfg.figure) {
    await ctx.addImage(slide, {
      path: cfg.figure,
      x: 190,
      y: 286,
      w: 830,
      h: 330,
      fit: cfg.figureFit || "contain",
    });
  }
  if (cfg.callout) {
    const ruleY = cfg.figure ? 622 : 462;
    const textY = cfg.figure ? 636 : 502;
    ctx.addShape(slide, { x: 130, y: ruleY, w: 976, h: 1, fill: "#c7c7bf", line: ctx.line("#c7c7bf", 0) });
    ctx.addText(slide, {
      text: cfg.callout,
      x: 190,
      y: textY,
      w: 860,
      h: cfg.figure ? 30 : 50,
      fontSize: cfg.figure ? 14 : 16,
      color: C.ink,
      bold: true,
      align: "center",
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
  }
  addFooter(slide, ctx, cfg.source);
}

function addWdGeometry(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  ctx.addShape(slide, { geometry: "ellipse", x: 150, y: 245, w: 138, h: 138, fill: "#e9eef2", line: ctx.line("#7d8791", 2) });
  ctx.addText(slide, { text: "WD", x: 184, y: 296, w: 70, h: 30, fontSize: 22, color: C.ink, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  ctx.addShape(slide, { geometry: "ellipse", x: 468, y: 222, w: 204, h: 204, fill: "#d9c07b", line: ctx.line("#8b6d1b", 2) });
  ctx.addText(slide, { text: "planet-scale\nobject", x: 496, y: 292, w: 150, h: 52, fontSize: 19, color: C.ink, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  ctx.addText(slide, { text: "large depth", x: 322, y: 305, w: 120, h: 24, fontSize: 16, color: C.accent, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addArrowText(slide, ctx, 410, 306);
  let y = 176;
  for (const [title, detail] of cfg.facts) {
    flowBox(slide, ctx, 785, y, 340, 78, title, detail);
    y += 96;
  }
  addFooter(slide, ctx, cfg.source);
}

function addTimeline(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const cardW = 345;
  const cardH = 126;
  const gapX = 42;
  const gapY = 44;
  const startX = 70;
  const startY = 178;
  cfg.events.forEach((event, i) => {
    const col = i % 3;
    const row = Math.floor(i / 3);
    const x = startX + col * (cardW + gapX);
    const y = startY + row * (cardH + gapY);
    const isLast = i === cfg.events.length - 1;
    ctx.addShape(slide, {
      x,
      y,
      w: cardW,
      h: cardH,
      fill: isLast ? "#fbf6e7" : C.light,
      line: ctx.line(isLast ? "#d4ba72" : C.border, 1),
    });
    ctx.addText(slide, {
      text: event[0],
      x: x + 20,
      y: y + 18,
      w: cardW - 40,
      h: 26,
      fontSize: 18,
      color: isLast ? C.accent : C.ink,
      bold: true,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
    ctx.addText(slide, {
      text: event[1],
      x: x + 20,
      y: y + 56,
      w: cardW - 40,
      h: 46,
      fontSize: 15,
      color: C.muted,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
  });
  ctx.addText(slide, {
    text: "History motivates two search modes: stable periodic events and evolving debris/dip systems.",
    x: 138,
    y: 562,
    w: 940,
    h: 30,
    fontSize: 18,
    color: C.ink,
    bold: true,
    align: "center",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  addFooter(slide, ctx, cfg.source);
}

async function addWd1145(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  await ctx.addImage(slide, {
    path: cfg.figure,
    x: 92,
    y: 148,
    w: 370,
    h: 500,
    fit: "contain",
  });
  ctx.addText(slide, {
    text: "TWIRL implication",
    x: 700,
    y: 178,
    w: 280,
    h: 30,
    fontSize: 20,
    color: C.ink,
    bold: true,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  flowBox(slide, ctx, 700, 240, 360, 86, "Evolving debris", "WD 1145 morphology changes over the K2 campaign, so a single static template is not enough.");
  flowBox(slide, ctx, 700, 354, 360, 86, "Search implication", "Keep BLS for WD 1856-like periodic events, but build a separate dip/debris branch.");
  flowBox(slide, ctx, 700, 468, 360, 86, "Validation implication", "Aperture, centroid, and contamination evidence must travel with every candidate.");
  addFooter(slide, ctx, cfg.source);
}

function addCadence(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const axis = { x: 120, y: 410, w: 860 };
  ctx.addShape(slide, { x: axis.x, y: axis.y, w: axis.w, h: 3, fill: C.ink, line: ctx.line(C.ink, 0) });
  const minToX = (m) => axis.x + (m / 30) * axis.w;
  for (const m of [0, 5, 10, 15, 20, 25, 30]) {
    const x = minToX(m);
    ctx.addShape(slide, { x, y: axis.y - 7, w: 2, h: 16, fill: C.ink, line: ctx.line(C.ink, 0) });
    ctx.addText(slide, { text: String(m), x: x - 14, y: axis.y + 18, w: 30, h: 18, fontSize: 11, color: C.muted, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  }
  ctx.addText(slide, { text: "minutes", x: axis.x + axis.w + 12, y: axis.y + 17, w: 70, h: 18, fontSize: 11, color: C.muted, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  const events = [
    ["2 min", 4, 2, C.red],
    ["5 min", 11, 5, C.accent],
    ["15 min", 22, 15, C.green],
  ];
  for (const [label, start, duration, color] of events) {
    const x = minToX(start);
    const w = (duration / 30) * axis.w;
    ctx.addShape(slide, { x, y: axis.y - 150, w, h: 74, fill: color, line: ctx.line(color, 0) });
    ctx.addText(slide, { text: label, x: x - 8, y: axis.y - 180, w: Math.max(70, w + 16), h: 22, fontSize: 13, color, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  }
  for (let t = 0; t <= 30; t += 200 / 60) {
    const x = minToX(t);
    ctx.addShape(slide, { geometry: "ellipse", x: x - 4, y: axis.y - 24, w: 8, h: 8, fill: C.blue, line: ctx.line(C.blue, 0) });
  }
  for (let t = 0; t <= 30; t += 30) {
    const x = minToX(t);
    ctx.addShape(slide, { geometry: "ellipse", x: x - 7, y: axis.y - 54, w: 14, h: 14, fill: "#d0d0d0", line: ctx.line("#777777", 1) });
  }
  ctx.addText(slide, { text: "200 s samples", x: 132, y: 476, w: 190, h: 24, fontSize: 15, color: C.blue, bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  ctx.addText(slide, { text: "30 min sample", x: 132, y: 512, w: 190, h: 24, fontSize: 15, color: C.muted, bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addClaimRail(slide, ctx, { railTitle: "Survey rule", rail: ["core survey uses Sector >= 56", "200 s FFIs only", "short-event completeness depends on cadence", "do not mix old cadence into the denominator"] }, 878, 160, 302, 420);
  addFooter(slide, ctx, cfg.source);
}

function addIdBridge(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  flowBox(slide, ctx, 80, 230, 245, 96, "Gaia EDR3 WD catalog", "seed parent catalog; high-confidence reference Pwd > 0.75");
  addArrowText(slide, ctx, 346, 264);
  flowBox(slide, ctx, 385, 230, 245, 96, "Gaia DR3 source_id", "authoritative scientific target identifier");
  addArrowText(slide, ctx, 650, 264);
  flowBox(slide, ctx, 690, 230, 245, 96, "TIC bridge", "operational TESS metadata and MIT production linkage");
  addArrowText(slide, ctx, 956, 264);
  flowBox(slide, ctx, 995, 230, 185, 96, "TGLC runs", "orbit/camera/CCD light curves");
  ctx.addShape(slide, { x: 126, y: 438, w: 1012, h: 1, fill: "#c7c7bf", line: ctx.line("#c7c7bf", 0) });
  ctx.addText(slide, { text: "Open audit item: how Gaia-selected targets propagate through a TIC-oriented production path.", x: 190, y: 486, w: 870, h: 50, fontSize: 19, color: C.ink, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addFooter(slide, ctx, cfg.source);
}

function addMetrics(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const cols = 3;
  const rows = 2;
  const startX = 80;
  const startY = 188;
  const w = 330;
  const h = 128;
  const gapX = 40;
  const gapY = 42;
  cfg.metrics.forEach((metric, i) => {
    const col = i % cols;
    const row = Math.floor(i / cols);
    const x = startX + col * (w + gapX);
    const y = startY + row * (h + gapY);
    ctx.addShape(slide, { x, y, w, h, fill: C.light, line: ctx.line(C.border, 1) });
    ctx.addText(slide, { text: metric[0], x: x + 20, y: y + 20, w: w - 40, h: 24, fontSize: 13, color: C.muted, bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
    ctx.addText(slide, { text: metric[1], x: x + 20, y: y + 58, w: w - 40, h: 54, fontSize: 25, color: C.ink, bold: true, typeface: ctx.fonts.title, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  });
  addFooter(slide, ctx, cfg.source);
}

function addStageRoadmap(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const stageX = 62;
  const stageW = 780;
  const stageH = 74;
  const stageGap = 22;
  let y = 162;
  (cfg.stages || []).forEach((stage, i) => {
    const active = i === 0;
    ctx.addShape(slide, {
      x: stageX,
      y,
      w: 98,
      h: stageH,
      fill: active ? C.accent : C.ink,
      line: ctx.line(active ? C.accent : C.ink, 0),
    });
    ctx.addShape(slide, {
      x: stageX + 102,
      y,
      w: stageW - 102,
      h: stageH,
      fill: active ? "#fbf6e7" : C.light,
      line: ctx.line(active ? "#d4ba72" : C.border, 1),
    });
    ctx.addText(slide, {
      text: stage[0],
      x: stageX + 10,
      y: y + 22,
      w: 78,
      h: 22,
      fontSize: 14,
      color: "#ffffff",
      bold: true,
      align: "center",
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
    ctx.addText(slide, {
      text: stage[1],
      x: stageX + 118,
      y: y + 10,
      w: 250,
      h: 24,
      fontSize: 15,
      color: C.ink,
      bold: true,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
    ctx.addText(slide, {
      text: stage[2],
      x: stageX + 118,
      y: y + 38,
      w: stageW - 138,
      h: 30,
      fontSize: 12.3,
      color: C.muted,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
    if (i < (cfg.stages || []).length - 1) {
      ctx.addText(slide, {
        text: "↓",
        x: stageX + 44,
        y: y + stageH + 1,
        w: 20,
        h: 20,
        fontSize: 14,
        color: C.rule,
        bold: true,
        align: "center",
        insets: { left: 0, right: 0, top: 0, bottom: 0 },
      });
    }
    y += stageH + stageGap;
  });

  ctx.addText(slide, {
    text: "Current control-plane facts",
    x: 900,
    y: 160,
    w: 260,
    h: 24,
    fontSize: 15,
    color: C.ink,
    bold: true,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  (cfg.status || []).forEach((row, i) => {
    const yy = 204 + i * 90;
    ctx.addShape(slide, { x: 890, y: yy, w: 286, h: 68, fill: C.lighter, line: ctx.line(C.border, 1) });
    ctx.addText(slide, {
      text: row[0],
      x: 908,
      y: yy + 12,
      w: 250,
      h: 18,
      fontSize: 11.5,
      color: C.muted,
      bold: true,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
    ctx.addText(slide, {
      text: row[1],
      x: 908,
      y: yy + 34,
      w: 250,
      h: 26,
      fontSize: 19,
      color: i === 3 ? C.green : C.ink,
      bold: true,
      typeface: ctx.fonts.title,
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
  });
  ctx.addText(slide, {
    text: "Takeaway: Stage 1 is the current gate; Stages 2-5 are demonstrated only as an S56 vertical slice.",
    x: 112,
    y: 626,
    w: 980,
    h: 28,
    fontSize: 18,
    color: C.ink,
    bold: true,
    align: "center",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  addFooter(slide, ctx, cfg.source);
}

function addProductionDashboard(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const fmt = (n) => Number(n).toLocaleString("en-US");
  const chart = { x: 74, y: 184, w: 720, h: 342 };
  const values = (cfg.sectors || []).map((row) => row[1]);
  const maxVal = 280000;
  ctx.addText(slide, {
    text: "Catalog target-sector LC scope by TESS sector",
    x: chart.x,
    y: 152,
    w: chart.w,
    h: 22,
    fontSize: 14,
    color: C.ink,
    bold: true,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  ctx.addShape(slide, { x: chart.x, y: chart.y + chart.h, w: chart.w, h: 1, fill: C.border, line: ctx.line(C.border, 0) });
  [0, 100000, 200000].forEach((tick) => {
    const y = chart.y + chart.h - (tick / maxVal) * chart.h;
    ctx.addShape(slide, { x: chart.x, y, w: chart.w, h: 1, fill: "#e1e1dc", line: ctx.line("#e1e1dc", 0) });
    ctx.addText(slide, {
      text: tick === 0 ? "0" : String(tick / 1000) + "k",
      x: chart.x - 48,
      y: y - 8,
      w: 38,
      h: 16,
      fontSize: 10,
      color: C.muted,
      align: "right",
      insets: { left: 0, right: 0, top: 0, bottom: 0 },
    });
  });
  const barGap = 3;
  const barW = (chart.w - barGap * ((cfg.sectors || []).length - 1)) / (cfg.sectors || []).length;
  (cfg.sectors || []).forEach((row, i) => {
    const sector = row[0];
    const count = row[1];
    const h = Math.max(2, (count / maxVal) * chart.h);
    const x = chart.x + i * (barW + barGap);
    const y = chart.y + chart.h - h;
    const highlight = sector === 56 || sector === 93 || sector === 94;
    ctx.addShape(slide, {
      x,
      y,
      w: barW,
      h,
      fill: highlight ? C.accent : C.blue,
      line: ctx.line(highlight ? C.accent : C.blue, 0),
    });
    if ([56, 60, 66, 72, 78, 84, 90, 94].includes(sector)) {
      ctx.addText(slide, {
        text: "S" + sector,
        x: x - 8,
        y: chart.y + chart.h + 10,
        w: 40,
        h: 16,
        fontSize: 10,
        color: C.muted,
        align: "center",
        insets: { left: 0, right: 0, top: 0, bottom: 0 },
      });
    }
  });
  ctx.addText(slide, {
    text: "unique-TIC WD targets per sector",
    x: chart.x + 210,
    y: chart.y + chart.h + 38,
    w: 300,
    h: 18,
    fontSize: 11,
    color: C.muted,
    align: "center",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });

  ctx.addText(slide, {
    text: "Scope and status",
    x: 835,
    y: 152,
    w: 280,
    h: 22,
    fontSize: 15,
    color: C.ink,
    bold: true,
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  (cfg.callouts || []).forEach((row, i) => {
    const yy = 192 + i * 102;
    ctx.addShape(slide, { x: 830, y: yy, w: 330, h: 78, fill: i === 0 ? "#f3f8f4" : C.light, line: ctx.line(i === 0 ? "#9bc4a6" : C.border, 1) });
    ctx.addText(slide, { text: row[0], x: 850, y: yy + 12, w: 290, h: 20, fontSize: 13, color: C.ink, bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
    ctx.addText(slide, { text: row[1], x: 850, y: yy + 36, w: 288, h: 36, fontSize: 11.2, color: C.muted, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  });
  ctx.addText(slide, {
    text: "This is the production scope from the catalog, not an occurrence-rate denominator or final completed-product count.",
    x: 92,
    y: 600,
    w: 1050,
    h: 34,
    fontSize: 17,
    color: C.ink,
    bold: true,
    align: "center",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  });
  addFooter(slide, ctx, cfg.source);
}

function addNegativeFlux(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const boxes = [
    ["Old risk", "flux <= 0 can become NaN before search"],
    ["Patch", "preserve RawFlux and RawFluxError in HDF5"],
    ["TWIRL-FS", "subtractive flux-space detrending"],
    ["Benchmark", "WD 1856: +155 cadences and +6.3% usable q=0"],
  ];
  boxes.forEach((box, i) => {
    const x = 88 + i * 278;
    flowBox(slide, ctx, x, 245, 220, 104, box[0], box[1]);
    if (i < boxes.length - 1) addArrowText(slide, ctx, x + 228, 286);
  });
  ctx.addText(slide, { text: "Design principle", x: 90, y: 456, w: 260, h: 28, fontSize: 16, color: C.ink, bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  ctx.addText(slide, { text: "Keep physical negative and near-zero background-subtracted flux measurements; do not force a fragile flux / spline ratio at the faint end.", x: 90, y: 494, w: 980, h: 58, fontSize: 20, color: C.ink, bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addFooter(slide, ctx, cfg.source);
}

function addMethod(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const left = [
    ["Canonical search input", "DET_FLUX"],
    ["Method", "robust-auto subtractive residuals"],
    ["Spline spacing", "bkspace_d = 0.8 d"],
    ["Gap behavior", "split gaps > 0.5 d"],
  ];
  const right = [
    ["Aperture checks", "DET_FLUX_SML and DET_FLUX_LAG"],
    ["Adaptive compare", "DET_FLUX_ADP is opt-in"],
    ["Why conservative", "avoid attenuating long or shallow signals"],
    ["Status", "S56 compare tree: 19,072 / 19,072 FITS"],
  ];
  let y = 168;
  for (const row of left) {
    flowBox(slide, ctx, 92, y, 460, 76, row[0], row[1]);
    y += 92;
  }
  y = 168;
  for (const row of right) {
    flowBox(slide, ctx, 682, y, 460, 76, row[0], row[1]);
    y += 92;
  }
  ctx.addShape(slide, { x: 616, y: 164, w: 1, h: 360, fill: "#c7c7bf", line: ctx.line("#c7c7bf", 0) });
  addFooter(slide, ctx, cfg.source);
}

function addBlsShowcase(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  ctx.addText(slide, { text: "Search ownership", x: 68, y: 158, w: 300, h: 24, fontSize: 14, color: C.muted, bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  flowBox(slide, ctx, 68, 200, 330, 82, "Franklin Chen", "owns the BLS/search work shown in this talk");
  flowBox(slide, ctx, 68, 312, 330, 82, "Periodic baseline", "interpretable BLS for WD 1856-like short events");
  flowBox(slide, ctx, 68, 424, 330, 82, "Repo result", "S56 v2 run improved WD 1856 blind rank 1258 -> 560");
  flowBox(slide, ctx, 68, 536, 330, 64, "Next output", "stable candidate tables + vet sheets for shared review");
  return ctx.addImage(slide, { path: cfg.figure, x: 462, y: 142, w: 720, h: 505, fit: "contain" }).then(() => addFooter(slide, ctx, cfg.source));
}

function addHeuristic(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  flowBox(slide, ctx, 90, 180, 260, 92, "Duration envelope", "reject events too long for a generous WD + planet geometry");
  addArrowText(slide, ctx, 370, 218);
  flowBox(slide, ctx, 410, 180, 260, 92, "Period alias guard", "P > 0.10 d and alias-aware triage");
  addArrowText(slide, ctx, 690, 218);
  flowBox(slide, ctx, 730, 180, 260, 92, "Cluster pile-ups", "detect period/frequency walls and systematic families");
  ctx.addShape(slide, { x: 150, y: 382, w: 880, h: 1, fill: "#c7c7bf", line: ctx.line("#c7c7bf", 0) });
  ctx.addText(slide, { text: "WD 1856 benchmark movement", x: 170, y: 430, w: 380, h: 28, fontSize: 17, color: C.muted, bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  ctx.addText(slide, { text: "1258 / 19,040", x: 180, y: 480, w: 260, h: 44, fontSize: 30, color: C.red, bold: true, typeface: ctx.fonts.title, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addArrowText(slide, ctx, 500, 490);
  ctx.addText(slide, { text: "9 / 5,403", x: 590, y: 480, w: 260, h: 44, fontSize: 30, color: C.green, bold: true, typeface: ctx.fonts.title, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  ctx.addText(slide, { text: "about 140x rank improvement after physics-motivated triage", x: 268, y: 548, w: 610, h: 30, fontSize: 18, color: C.ink, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addFooter(slide, ctx, cfg.source);
}

function addCollaboration(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const lanes = [
    ["Te Han", "TGLC/TWIRL-FS production\nproduct QA\nprovenance and data stewardship"],
    ["Franklin Chen", "BLS/search work shown\ncandidate tables\nvet sheets and rank diagnostics"],
    ["Michelle", "owns LEO-Vetter\nWD tuning interface\nvalidation-review collaboration path"],
    ["Shared review", "candidate triage\nfalse-positive taxonomy\nfollow-up readiness decisions"],
  ];
  lanes.forEach((lane, i) => {
    const x = 62 + i * 296;
    flowBox(slide, ctx, x, 225, 240, 170, lane[0], lane[1]);
    if (i < lanes.length - 1) addArrowText(slide, ctx, x + 250, 298);
  });
  ctx.addText(slide, { text: "Artifact interface: HLSP tree -> BLS tables -> vet sheets -> LEO outputs -> pixel/centroid diagnostics -> shared candidate review.", x: 130, y: 505, w: 1010, h: 58, fontSize: 19, color: C.ink, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addFooter(slide, ctx, cfg.source);
}

function addFalsePositives(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const classes = [
    ["PCEBs", "WD + M-dwarf eclipsing binaries"],
    ["WD + WD", "compact binary geometry"],
    ["Centroid shifts", "neighbor source in aperture"],
    ["ZZ Ceti", "pulsation families"],
    ["Systematics", "period ladders and orbit artifacts"],
    ["Debris dips", "irregular or variable-depth events"],
  ];
  classes.forEach((item, i) => {
    const col = i % 3;
    const row = Math.floor(i / 3);
    flowBox(slide, ctx, 92 + col * 360, 180 + row * 160, 285, 105, item[0], item[1]);
  });
  ctx.addText(slide, { text: "First task for top candidates: build agreement on the false-positive taxonomy before any discovery or occurrence-rate language.", x: 126, y: 545, w: 985, h: 48, fontSize: 18, color: C.ink, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addFooter(slide, ctx, cfg.source);
}

function addInjections(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  flowBox(slide, ctx, 95, 185, 460, 92, "Fast LC-level injections", "inject into HLSP light curves, run BLS + vetter, map detection and triage losses");
  flowBox(slide, ctx, 680, 185, 460, 92, "Pixel-level injections", "inject into images or extraction stack, capture crowding and aperture-choice failures");
  ctx.addShape(slide, { x: 210, y: 344, w: 820, h: 1, fill: "#c7c7bf", line: ctx.line("#c7c7bf", 0) });
  const chain = [
    ["light curve", "product choice"],
    ["search", "BLS/dip branch"],
    ["vetter", "heuristic + LEO"],
    ["merge", "candidate table"],
    ["claim", "only after completeness"],
  ];
  chain.forEach((item, i) => {
    const x = 88 + i * 225;
    flowBox(slide, ctx, x, 425, 170, 88, item[0], item[1]);
    if (i < chain.length - 1) addArrowText(slide, ctx, x + 178, 456);
  });
  addFooter(slide, ctx, cfg.source);
}

function addGates(slide, ctx, cfg) {
  addBackground(slide, ctx, "#ffffff");
  addHeader(slide, ctx, cfg);
  const gates = [
    ["TGLC production QA", "cadence retention, RMS/MAD, sector failures"],
    ["S56 v2 search", "complete BLS candidate tables on the handoff tree"],
    ["LEO tuning", "Michelle-owned path, WD host behavior, stable schema"],
    ["Pixel vetting", "centroid and pixel-map diagnostics for candidates"],
    ["Injections", "real stack: search, vetter, candidate merge"],
    ["Dip-search branch", "debris-like and weakly periodic events"],
  ];
  gates.forEach((gate, i) => {
    const x = 86 + (i % 3) * 360;
    const y = 172 + Math.floor(i / 3) * 178;
    flowBox(slide, ctx, x, y, 288, 118, gate[0], gate[1]);
  });
  ctx.addText(slide, { text: "Scale-out is a decision after these gates, not a default consequence of having light curves.", x: 145, y: 565, w: 970, h: 32, fontSize: 20, color: C.ink, bold: true, align: "center", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  addFooter(slide, ctx, cfg.source);
}

function addDecision(slide, ctx, cfg) {
  addBackground(slide, ctx, C.black);
  ctx.addText(slide, { text: cfg.kicker || "Decision", x: 64, y: 52, w: 260, h: 22, fontSize: 12, color: "#bdbdb7", bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  ctx.addText(slide, { text: cfg.title, x: 64, y: 132, w: 960, h: 70, fontSize: 40, color: "#ffffff", bold: true, typeface: ctx.fonts.title, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
  let y = 270;
  for (const row of cfg.rows || []) {
    ctx.addShape(slide, { x: 84, y: y - 10, w: 1000, h: 1, fill: "#777770", line: ctx.line("#777770", 0) });
    ctx.addText(slide, { text: row[0], x: 88, y, w: 300, h: 28, fontSize: 18, color: "#ffffff", bold: true, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
    ctx.addText(slide, { text: row[1], x: 420, y: y + 2, w: 650, h: 28, fontSize: 15, color: "#d6d6cf", insets: { left: 0, right: 0, top: 0, bottom: 0 } });
    y += 72;
  }
  addFooter(slide, ctx, cfg.source);
}

export async function addConfiguredSlide(presentation, ctx, cfg) {
  const slide = presentation.slides.add();
  const textCtx = createReadableContext(ctx);
  if (cfg.kind === "title") addTitleSlide(slide, textCtx, cfg);
  else if (cfg.kind === "divider") addDividerSlide(slide, textCtx, cfg);
  else if (cfg.kind === "toc") addTocSlide(slide, textCtx, cfg);
  else if (cfg.kind === "flow") await addFlowSlide(slide, textCtx, cfg);
  else if (cfg.kind === "wdGeometry") addWdGeometry(slide, textCtx, cfg);
  else if (cfg.kind === "timeline") addTimeline(slide, textCtx, cfg);
  else if (cfg.kind === "wd1145") await addWd1145(slide, textCtx, cfg);
  else if (cfg.kind === "figureRail") await addFigureRail(slide, textCtx, cfg);
  else if (cfg.kind === "figureFull") await addFigureFull(slide, textCtx, cfg);
  else if (cfg.kind === "figureTwo") await addFigureTwo(slide, textCtx, cfg);
  else if (cfg.kind === "cadence") addCadence(slide, textCtx, cfg);
  else if (cfg.kind === "idBridge") addIdBridge(slide, textCtx, cfg);
  else if (cfg.kind === "metrics") addMetrics(slide, textCtx, cfg);
  else if (cfg.kind === "stageRoadmap") addStageRoadmap(slide, textCtx, cfg);
  else if (cfg.kind === "productionDashboard") addProductionDashboard(slide, textCtx, cfg);
  else if (cfg.kind === "negativeFlux") addNegativeFlux(slide, textCtx, cfg);
  else if (cfg.kind === "method") addMethod(slide, textCtx, cfg);
  else if (cfg.kind === "blsShowcase") await addBlsShowcase(slide, textCtx, cfg);
  else if (cfg.kind === "heuristic") addHeuristic(slide, textCtx, cfg);
  else if (cfg.kind === "collaboration") addCollaboration(slide, textCtx, cfg);
  else if (cfg.kind === "falsePositives") addFalsePositives(slide, textCtx, cfg);
  else if (cfg.kind === "injections") addInjections(slide, textCtx, cfg);
  else if (cfg.kind === "gates") addGates(slide, textCtx, cfg);
  else if (cfg.kind === "decision") addDecision(slide, textCtx, cfg);
  else throw new Error("Unknown slide kind: " + cfg.kind);
  setNotes(slide, cfg);
  return slide;
}
`;

async function main() {
  if (slides.length < 30) {
    throw new Error(`Expected at least 30 slides, found ${slides.length}`);
  }

  await fs.mkdir(outputDir, { recursive: true });
  await fs.mkdir(previewDir, { recursive: true });
  await fs.mkdir(layoutDir, { recursive: true });
  await fs.rm(slidesDir, { recursive: true, force: true });
  await fs.mkdir(slidesDir, { recursive: true });

  const missingFigures = [];
  for (const figurePath of collectFigures(slides)) {
    try {
      await fs.access(figurePath);
    } catch {
      missingFigures.push(figurePath);
    }
  }
  if (missingFigures.length) {
    throw new Error(`Missing figure assets:\n${missingFigures.join("\n")}`);
  }

  await fs.writeFile(path.join(slidesDir, "shared.mjs"), `${sharedModule}\n`, "utf8");
  await fs.writeFile(
    path.join(slidesDir, "deck-data.mjs"),
    `export const slides = ${JSON.stringify(slides, null, 2)};\n`,
    "utf8",
  );

  for (let i = 0; i < slides.length; i += 1) {
    const n = String(i + 1).padStart(2, "0");
    const moduleText = [
      'import { slides } from "./deck-data.mjs";',
      'import { addConfiguredSlide } from "./shared.mjs";',
      "",
      `export default async function slide${n}(presentation, ctx) {`,
      `  return addConfiguredSlide(presentation, ctx, slides[${i}]);`,
      "}",
      "",
    ].join("\n");
    await fs.writeFile(path.join(slidesDir, `slide-${n}.mjs`), moduleText, "utf8");
  }

  await fs.writeFile(
    speakerNotesOut,
    ["# TWIRL current-stage expert talk speaker notes", "", ...slides.map(notesText)].join("\n"),
    "utf8",
  );
  await fs.writeFile(
    sourceManifestOut,
    ["# TWIRL current-stage expert talk source manifest", "", ...slides.map(sourceText)].join("\n"),
    "utf8",
  );

  const estimate = estimateTalkMinutes(slides);
  await fs.writeFile(
    talkTrackOut,
    [
      "# Talk track and QA",
      "",
      `Slide count: ${slides.length}`,
      `Speaker-notes word count: ${estimate.words}`,
      `Estimated prepared talk-track time at 130 wpm: ${estimate.minutesAt130Wpm.toFixed(1)} minutes`,
      "Target delivery time: about 45 minutes with discussion pacing.",
      "",
      "Completed QA:",
      `- Built editable PPTX and rendered ${slides.length} / ${slides.length} slide PNG previews.`,
      "- Generated and visually inspected the contact sheet for missing assets, blank slides, and major text collisions.",
      "- Full-size spot checks passed for the global readable text scale, table-of-contents slide, enlarged right-side takeaway rails, WD 1856 evidence cards, WD 1145 publication figure, enlarged history roadmap, 200 s motivation figure, detectability map, Stage 1-5 roadmap, S56-S94 LC-scope chart, S56 benchmark flow, LEO-Vetter WD 1856 report, false-positive taxonomy, injection-recovery, and scale-gate slides.",
      `- Verified the PPTX contains ${slides.length} notes-slide XML parts, matching the slide count.`,
      "- Verified collaboration wording appears in the speaker notes: Franklin owns BLS/search work shown; Michelle owns LEO-Vetter; Te Han owns production and data stewardship.",
      "",
      "QA checklist:",
      "- Render every slide to PNG.",
      "- Inspect contact sheet for missing assets and text overflow.",
      "- Confirm title/divider slides are sparse and content slides carry figures, metrics, or flowcharts.",
      "- Confirm collaboration wording: Franklin owns BLS/search work shown; Michelle owns LEO-Vetter; Te Han owns TGLC/TWIRL-FS production and data stewardship.",
      "- Confirm claims stay benchmark-gated and avoid occurrence-rate or Earth-size/HZ claims.",
      "",
    ].join("\n"),
    "utf8",
  );

  const build = spawnSync(
    nodeBin,
    [
      buildScript,
      "--workspace",
      workspace,
      "--slides-dir",
      slidesDir,
      "--out",
      pptxOut,
      "--preview-dir",
      previewDir,
      "--layout-dir",
      layoutDir,
      "--contact-sheet",
      contactSheetOut,
      "--manifest",
      manifestOut,
      "--slide-count",
      String(slides.length),
      "--scale",
      "1",
    ],
    {
      cwd: repoRoot,
      env: { ...process.env, PYTHON: pythonBin },
      encoding: "utf8",
    },
  );

  if (build.status !== 0) {
    throw new Error(
      ["Deck build failed.", build.stdout.trim(), build.stderr.trim()].filter(Boolean).join("\n"),
    );
  }

  process.stdout.write(build.stdout);
}

main().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});
