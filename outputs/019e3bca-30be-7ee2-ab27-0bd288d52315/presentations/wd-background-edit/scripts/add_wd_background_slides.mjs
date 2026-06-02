import fs from "node:fs/promises";
import path from "node:path";
import { spawnSync } from "node:child_process";

import {
  ensureArtifactToolWorkspace,
  importArtifactTool,
  saveBlobToFile,
} from "/Users/tehan/.codex/plugins/cache/openai-primary-runtime/presentations/26.601.10930/skills/presentations/scripts/artifact_tool_utils.mjs";

const workspace = "/Users/tehan/PycharmProjects/TWIRL/outputs/019e3bca-30be-7ee2-ab27-0bd288d52315/presentations/wd-background-edit";
const sourcePptx = "/Users/tehan/PycharmProjects/TWIRL/outputs/019e3bca-30be-7ee2-ab27-0bd288d52315/presentations/twirl-wd1856-franklin-flow-edit/output/twirl-wd1856-keynote-bw-talk-franklin-flow.pptx";
const outputDir = path.join(workspace, "output");
const previewDir = path.join(workspace, "previews");
const outputPptx = path.join(outputDir, "twirl-wd1856-keynote-bw-talk-franklin-flow-wd-background.pptx");
const contactSheet = path.join(outputDir, "twirl-wd1856-keynote-bw-talk-franklin-flow-wd-background-contact-sheet.png");
const manifestPath = path.join(outputDir, "wd-background-edit-manifest.json");

const W = 1280;
const H = 720;
const BLACK = "#050505";
const INK = "#111111";
const MID = "#5f5f5f";
const LIGHT = "#d9d9d9";
const PALE = "#f3f3f3";
const WHITE = "#ffffff";
const CLEAR = "#00000000";

function slidesFromPresentation(presentation) {
  if (Array.isArray(presentation.slides?.items)) return presentation.slides.items;
  if (Number.isInteger(presentation.slides?.count) && typeof presentation.slides.getItem === "function") {
    return Array.from({ length: presentation.slides.count }, (_, index) => presentation.slides.getItem(index));
  }
  throw new Error("Could not enumerate presentation slides.");
}

function shape(slide, {
  x = 0,
  y = 0,
  w,
  h,
  geometry = "rect",
  fill = CLEAR,
  lineFill = CLEAR,
  lineWidth = 0,
  name,
}) {
  return slide.shapes.add({
    geometry,
    name,
    position: { left: x, top: y, width: w, height: h },
    fill,
    line: { style: "solid", fill: lineFill, width: lineWidth },
  });
}

function text(slide, value, {
  x,
  y,
  w,
  h,
  size = 24,
  color = INK,
  bold = false,
  align = "left",
  valign = "top",
  face = "Aptos",
  name,
  insets = { left: 0, right: 0, top: 0, bottom: 0 },
}) {
  const box = shape(slide, { x, y, w, h, fill: CLEAR, lineWidth: 0, name });
  box.text = value;
  box.text.fontSize = size;
  box.text.color = color;
  box.text.bold = bold;
  box.text.typeface = face;
  box.text.alignment = align;
  box.text.verticalAlignment = valign;
  box.text.insets = insets;
  return box;
}

function rule(slide, x, y, w, color = INK, thickness = 1.5) {
  return shape(slide, { x, y, w, h: thickness, fill: color, lineWidth: 0 });
}

function circle(slide, cx, cy, r, options = {}) {
  return shape(slide, {
    x: cx - r,
    y: cy - r,
    w: 2 * r,
    h: 2 * r,
    geometry: "ellipse",
    fill: options.fill ?? CLEAR,
    lineFill: options.lineFill ?? INK,
    lineWidth: options.lineWidth ?? 1.5,
    name: options.name,
  });
}

function addSimpleArrow(slide, x1, y1, x2, y2, color = INK) {
  const width = Math.max(1, x2 - x1);
  rule(slide, x1, y1, width, color, 1.2);
  shape(slide, {
    x: x2 - 8,
    y: y2 - 6,
    w: 12,
    h: 12,
    geometry: "triangle",
    fill: color,
    lineWidth: 0,
  });
}

function addSectionSlide(presentation) {
  const slide = presentation.slides.add();
  shape(slide, { x: 0, y: 0, w: W, h: H, fill: BLACK, lineWidth: 0 });
  text(slide, "BACKGROUND", {
    x: 58,
    y: 38,
    w: 240,
    h: 26,
    size: 12,
    color: WHITE,
    bold: true,
    face: "Aptos",
  });
  text(slide, "White dwarfs set the scale.", {
    x: 58,
    y: 266,
    w: 840,
    h: 62,
    size: 38,
    color: WHITE,
    bold: true,
    face: "Aptos Display",
  });
  rule(slide, 58, 405, 780, "#d8d8d8", 1.4);
  text(slide, "The opportunity is large transit signals; the risk is proving where they came from.", {
    x: 58,
    y: 430,
    w: 920,
    h: 52,
    size: 16,
    color: "#dcdcdc",
  });

  circle(slide, 1022, 312, 72, { fill: CLEAR, lineFill: "#dadada", lineWidth: 2 });
  circle(slide, 1022, 312, 17, { fill: WHITE, lineFill: WHITE, lineWidth: 0 });
  circle(slide, 1081, 282, 10, { fill: CLEAR, lineFill: WHITE, lineWidth: 1.5 });
  text(slide, "deep", { x: 948, y: 415, w: 66, h: 28, size: 15, color: WHITE, bold: true, align: "center" });
  text(slide, "short", { x: 1028, y: 415, w: 70, h: 28, size: 15, color: WHITE, bold: true, align: "center" });
  text(slide, "faint", { x: 1105, y: 415, w: 70, h: 28, size: 15, color: WHITE, bold: true, align: "center" });
  return slide;
}

function addWdBasicsSlide(presentation) {
  const slide = presentation.slides.add();
  shape(slide, { x: 0, y: 0, w: W, h: H, fill: WHITE, lineWidth: 0 });
  text(slide, "White dwarfs make small planets look large.", {
    x: 58,
    y: 38,
    w: 980,
    h: 46,
    size: 26,
    bold: true,
    face: "Aptos Display",
  });

  rule(slide, 72, 184, 660, LIGHT, 1.1);
  const xs = [152, 362, 572];
  circle(slide, xs[0], 184, 46, { fill: WHITE, lineFill: INK, lineWidth: 1.7 });
  circle(slide, xs[1], 184, 90, { fill: CLEAR, lineFill: "#8a8a8a", lineWidth: 1.4 });
  circle(slide, xs[2], 184, 18, { fill: INK, lineFill: INK, lineWidth: 0 });
  addSimpleArrow(slide, 210, 184, 302, 184, MID);
  addSimpleArrow(slide, 456, 184, 532, 184, MID);
  text(slide, "main sequence", { x: 82, y: 268, w: 150, h: 24, size: 14, bold: true, align: "center" });
  text(slide, "red giant", { x: 286, y: 268, w: 150, h: 24, size: 14, bold: true, align: "center" });
  text(slide, "white dwarf", { x: 506, y: 268, w: 150, h: 24, size: 14, bold: true, align: "center" });
  text(slide, "compact remnant", { x: 498, y: 295, w: 170, h: 22, size: 12, color: MID, align: "center" });

  text(slide, "radius sets transit depth", { x: 788, y: 126, w: 360, h: 26, size: 16, bold: true });
  rule(slide, 790, 164, 370, INK, 1.2);
  circle(slide, 870, 236, 70, { fill: WHITE, lineFill: INK, lineWidth: 1.5 });
  circle(slide, 918, 236, 11, { fill: INK, lineFill: INK, lineWidth: 0 });
  text(slide, "large star", { x: 788, y: 324, w: 170, h: 24, size: 13, color: MID, align: "center" });
  text(slide, "shallow signal", { x: 788, y: 351, w: 170, h: 24, size: 14, bold: true, align: "center" });
  circle(slide, 1088, 236, 27, { fill: INK, lineFill: INK, lineWidth: 0 });
  circle(slide, 1104, 236, 31, { fill: CLEAR, lineFill: INK, lineWidth: 2.2 });
  text(slide, "white dwarf", { x: 1008, y: 324, w: 190, h: 24, size: 13, color: MID, align: "center" });
  text(slide, "large signal", { x: 1008, y: 351, w: 190, h: 24, size: 14, bold: true, align: "center" });

  const facts = [
    ["small host", "Earth-size scale makes deep eclipses possible"],
    ["fast clock", "events can be minutes, so cadence matters"],
    ["old system", "survival and migration are the science"],
    ["faint target", "photometry and crowding become the limiting systematics"],
  ];
  facts.forEach(([label, body], index) => {
    const y = 445 + index * 48;
    rule(slide, 58, y - 11, 360, "#cfcfcf", 1);
    text(slide, label, { x: 58, y, w: 142, h: 24, size: 16, bold: true });
    text(slide, body, { x: 220, y: y + 2, w: 700, h: 24, size: 14, color: MID });
  });
  text(slide, "depth approx (planet radius / star radius)^2", {
    x: 788,
    y: 472,
    w: 390,
    h: 34,
    size: 15,
    bold: true,
    align: "center",
  });
  text(slide, "White dwarfs trade easy amplitude for hard validation.", {
    x: 788,
    y: 523,
    w: 390,
    h: 44,
    size: 17,
    align: "center",
  });
  return slide;
}

function lane(slide, x, y, w, label, body) {
  rule(slide, x, y, w, INK, 1.2);
  text(slide, label, { x, y: y + 18, w, h: 26, size: 15, bold: true, align: "center" });
  text(slide, body, { x: x + 10, y: y + 58, w: w - 20, h: 72, size: 12, color: MID, align: "center" });
}

function addCurrentEffortsSlide(presentation) {
  const slide = presentation.slides.add();
  shape(slide, { x: 0, y: 0, w: W, h: H, fill: WHITE, lineWidth: 0 });
  text(slide, "Current searches are converging from several angles.", {
    x: 58,
    y: 38,
    w: 1040,
    h: 46,
    size: 26,
    bold: true,
    face: "Aptos Display",
  });

  const y0 = 148;
  lane(slide, 74, y0, 215, "transits", "TESS, ZTF, high-cadence follow-up");
  lane(slide, 333, y0, 215, "infrared", "JWST/MIRI and excess constraints");
  lane(slide, 592, y0, 215, "pollution", "metals and debris trace remnant systems");
  lane(slide, 851, y0, 215, "vetting", "ground + space checks against blends and binaries");

  [181, 440, 699, 958].forEach((cx) => {
    circle(slide, cx, 356, 30, { fill: PALE, lineFill: INK, lineWidth: 1.2 });
    circle(slide, cx, 356, 8, { fill: INK, lineFill: INK, lineWidth: 0 });
  });

  text(slide, "WD 1856+534 b is the benchmark intact case.", {
    x: 74,
    y: 448,
    w: 760,
    h: 34,
    size: 21,
    bold: true,
  });
  rule(slide, 74, 508, 870, INK, 1.2);
  const events = [
    [112, "TESS", "deep periodic transit"],
    [394, "GTC / Spitzer", "same event, different systematics"],
    [676, "JWST / MIRI", "thermal emission confirmation"],
  ];
  events.forEach(([x, label, body]) => {
    circle(slide, x, 508, 8, { fill: INK, lineFill: INK, lineWidth: 0 });
    text(slide, label, { x: x - 52, y: 530, w: 110, h: 24, size: 13, bold: true, align: "center" });
    text(slide, body, { x: x - 88, y: 557, w: 176, h: 40, size: 11, color: MID, align: "center" });
  });
  text(slide, "The evidence base is richer than the confirmed-planet count.", {
    x: 812,
    y: 455,
    w: 315,
    h: 86,
    size: 18,
    bold: true,
    align: "center",
    valign: "mid",
  });
  rule(slide, 815, 555, 308, "#bdbdbd", 1);
  text(slide, "That is why benchmark recovery, contamination checks, and completeness tests matter.", {
    x: 806,
    y: 575,
    w: 330,
    h: 54,
    size: 13,
    color: MID,
    align: "center",
  });
  return slide;
}

function addChallengesSlide(presentation) {
  const slide = presentation.slides.add();
  shape(slide, { x: 0, y: 0, w: W, h: H, fill: WHITE, lineWidth: 0 });
  text(slide, "The hard part is proving the dip is real and local.", {
    x: 58,
    y: 38,
    w: 1040,
    h: 46,
    size: 26,
    bold: true,
    face: "Aptos Display",
  });

  const stages = [
    ["faint WDs", "low flux, negative residuals"],
    ["short events", "cadence retention matters"],
    ["crowded pixels", "blend and aperture disagreement"],
    ["false positives", "EBs, variables, systematics"],
    ["completeness", "injections through the real stack"],
  ];
  stages.forEach(([label, body], index) => {
    const x = 88 + index * 218;
    const y = 196 + index * 22;
    const w = 170 - index * 10;
    shape(slide, { x, y, w, h: 94, fill: index % 2 === 0 ? WHITE : PALE, lineFill: INK, lineWidth: 1.1 });
    text(slide, label, { x: x + 10, y: y + 14, w: w - 20, h: 24, size: 14, bold: true, align: "center" });
    text(slide, body, { x: x + 12, y: y + 45, w: w - 24, h: 44, size: 10.5, color: MID, align: "center" });
    if (index < stages.length - 1) addSimpleArrow(slide, x + w + 16, y + 46, x + 196, y + 46, MID);
  });

  rule(slide, 86, 436, 470, INK, 1.3);
  text(slide, "What TWIRL should contribute", { x: 86, y: 464, w: 420, h: 30, size: 18, bold: true });
  const bullets = [
    "expected-flux photometric QA",
    "WD 1856 recovery across apertures",
    "BLS + dip search on 200 s FFIs",
    "pixel-level localization checks",
    "injection recovery before occurrence claims",
  ];
  bullets.forEach((item, index) => {
    const y = 514 + index * 25;
    shape(slide, { x: 88, y: y + 5, w: 8, h: 8, geometry: "ellipse", fill: INK, lineWidth: 0 });
    text(slide, item, { x: 110, y, w: 450, h: 22, size: 13, color: INK });
  });

  rule(slide, 700, 436, 410, INK, 1.3);
  text(slide, "Discussion prompt", { x: 700, y: 464, w: 360, h: 30, size: 18, bold: true });
  text(slide, "Which failure mode blocks scaling first?", {
    x: 700,
    y: 520,
    w: 380,
    h: 66,
    size: 23,
    bold: true,
  });
  text(slide, "Use WD 1856 as the room's shared calibration target.", {
    x: 700,
    y: 594,
    w: 390,
    h: 42,
    size: 14,
    color: MID,
  });
  return slide;
}

function runContactSheet(previewPaths) {
  const python = process.env.PYTHON || "python";
  const script = "/Users/tehan/.codex/plugins/cache/openai-primary-runtime/presentations/26.601.10930/skills/presentations/scripts/make_contact_sheet.py";
  const result = spawnSync(python, [script, "--output", contactSheet, ...previewPaths], {
    encoding: "utf8",
    maxBuffer: 30 * 1024 * 1024,
  });
  if (result.status !== 0) {
    throw new Error([result.stdout, result.stderr].filter(Boolean).join("\n"));
  }
}

async function main() {
  await fs.mkdir(outputDir, { recursive: true });
  await fs.rm(previewDir, { recursive: true, force: true });
  await fs.mkdir(previewDir, { recursive: true });
  await ensureArtifactToolWorkspace(workspace);
  const { FileBlob, PresentationFile } = await importArtifactTool(workspace);

  const presentation = await PresentationFile.importPptx(await FileBlob.load(sourcePptx));
  const insertedSlides = [
    addSectionSlide(presentation),
    addWdBasicsSlide(presentation),
    addCurrentEffortsSlide(presentation),
    addChallengesSlide(presentation),
  ];
  insertedSlides.forEach((slide, index) => slide.moveTo(1 + index));

  const slides = slidesFromPresentation(presentation);
  const previewPaths = [];
  for (let index = 0; index < slides.length; index += 1) {
    const padded = String(index + 1).padStart(2, "0");
    const previewPath = path.join(previewDir, `slide-${padded}.png`);
    const preview = await presentation.export({ slide: slides[index], format: "png", scale: 0.8 });
    await saveBlobToFile(preview, previewPath);
    previewPaths.push(previewPath);
  }

  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(outputPptx);
  runContactSheet(previewPaths);

  const stat = await fs.stat(outputPptx);
  const manifest = {
    sourcePptx,
    outputPptx,
    outputBytes: stat.size,
    previewDir,
    contactSheet,
    insertedAfterSlide: 1,
    insertedSlides: [
      "Background section title",
      "White dwarf basics and transit geometry",
      "Current search efforts",
      "White dwarf planet search challenges and TWIRL scope",
    ],
    sourcesUsedForSlideClaims: [
      "Vanderburg et al. 2020, Nature, arXiv:2009.07282",
      "Limbach et al. 2025, arXiv:2504.16982",
      "white dwarf pollution/debris disk review literature summarized qualitatively",
    ],
  };
  await fs.writeFile(manifestPath, `${JSON.stringify(manifest, null, 2)}\n`, "utf8");
  console.log(JSON.stringify(manifest, null, 2));
}

main().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});
