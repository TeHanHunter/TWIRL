#!/usr/bin/env node
import fs from "node:fs/promises";
import path from "node:path";
import { createHash } from "node:crypto";
import { createRequire } from "node:module";
import { pathToFileURL } from "node:url";

const require = createRequire(import.meta.url);
const artifactToolUrl = pathToFileURL(require.resolve("@oai/artifact-tool")).href;
const { SpreadsheetFile, Workbook } = await import(artifactToolUrl);

process.on("uncaughtException", (error) => {
  console.error(JSON.stringify({ name: error.name, message: error.message }));
  process.exit(1);
});
process.on("unhandledRejection", (error) => {
  console.error(
    JSON.stringify({
      name: error?.name ?? "UnhandledRejection",
      message: error?.message ?? String(error),
    }),
  );
  process.exit(1);
});

const [inputJson, manifestJson, outputXlsx, previewDir] = process.argv.slice(2);
if (!inputJson || !manifestJson || !outputXlsx || !previewDir) {
  throw new Error(
    "usage: build_s56_planet_like_candidate_workbook.mjs " +
      "<candidates.json> <manifest.json> <output.xlsx> <preview-dir>",
  );
}

const rows = JSON.parse(await fs.readFile(inputJson, "utf8"));
const manifest = JSON.parse(await fs.readFile(manifestJson, "utf8"));
if (!Array.isArray(rows) || rows.length !== 30) {
  throw new Error(`expected 30 candidate rows, found ${rows.length}`);
}

const workbook = Workbook.create();
const overview = workbook.worksheets.add("Overview");
const candidates = workbook.worksheets.add("Candidates");
const provenance = workbook.worksheets.add("Provenance");

const navy = "#3B4D66";
const navyDark = "#26374D";
const paleBlue = "#DCE6F1";
const paleGreen = "#E7F2EC";
const paleYellow = "#FFF4D6";
const paleRed = "#F8E4E4";
const lightGray = "#E3E6EA";
const darkText = "#20252B";

overview.showGridLines = false;
overview.getRange("A1:H1").merge();
overview.getRange("A1").values = [["Sector 56 current Planet-like transit candidates"]];
overview.getRange("A1:H1").format = {
  fill: navyDark,
  font: { bold: true, color: "#FFFFFF", size: 18 },
  verticalAlignment: "center",
};
overview.getRange("A1:H1").format.rowHeight = 34;

overview.getRange("A3:B3").values = [["Candidate set", "Value"]];
overview.getRange("A4:A9").values = [
  ["Human-accepted S56 candidates"],
  ["Original human-adjudicated"],
  ["Model-enriched, human-confirmed"],
  ["Current two-aperture period matches"],
  ["Tier-1 searchable"],
  ["Odd/even review flags"],
];
overview.getRange("B4:B9").formulas = [
  ["=COUNTA('Candidates'!$C$2:$C$31)"],
  ['=COUNTIF(\'Candidates\'!$F$2:$F$31,"human-adjudicated")'],
  ['=COUNTIF(\'Candidates\'!$F$2:$F$31,"model-enriched, human-confirmed")'],
  ['=COUNTIF(\'Candidates\'!$Q$2:$Q$31,"yes")'],
  ['=COUNTIF(\'Candidates\'!$Y$2:$Y$31,"yes")'],
  ['=COUNTIF(\'Candidates\'!$W$2:$W$31,"review")'],
];
overview.getRange("A3:B3").format = {
  fill: navy,
  font: { bold: true, color: "#FFFFFF" },
};
overview.getRange("A3:B9").format.borders = {
  preset: "outside",
  style: "thin",
  color: "#AAB1B9",
};
overview.getRange("A4:A9").format.fill = "#F4F6F8";
overview.getRange("B4:B9").format.numberFormat = "0";

overview.getRange("D3:E3").values = [["Priority group", "Count"]];
overview.getRange("D4:D7").values = [
  ["Benchmark"],
  ["A — aligned current ADP periods"],
  ["B — small-aperture-led"],
  ["C — small-aperture-led + QA review"],
];
overview.getRange("E4:E7").formulas = [
  ['=COUNTIF(\'Candidates\'!$B$2:$B$31,"Benchmark")'],
  ['=COUNTIF(\'Candidates\'!$B$2:$B$31,"A")'],
  ['=COUNTIF(\'Candidates\'!$B$2:$B$31,"B")'],
  ['=COUNTIF(\'Candidates\'!$B$2:$B$31,"C")'],
];
overview.getRange("D3:E3").format = {
  fill: navy,
  font: { bold: true, color: "#FFFFFF" },
};
overview.getRange("D3:E7").format.borders = {
  preset: "outside",
  style: "thin",
  color: "#AAB1B9",
};
overview.getRange("D4").format.fill = paleBlue;
overview.getRange("D5").format.fill = paleGreen;
overview.getRange("D6").format.fill = paleYellow;
overview.getRange("D7").format.fill = paleRed;
overview.getRange("E4:E7").format.numberFormat = "0";

overview.getRange("A11:H11").merge();
overview.getRange("A11").values = [["How to read this workbook"]];
overview.getRange("A11:H11").format = {
  fill: navy,
  font: { bold: true, color: "#FFFFFF" },
};
overview.getRange("A12:H15").merge();
overview.getRange("A12").values = [[
  "This is a follow-up triage of final human-accepted Planet-like morphology " +
    "labels, not an astrophysical confirmation list. Group A requires aligned " +
    "periods in the current DET_FLUX_ADP_SML and DET_FLUX_ADP renderings. " +
    "Groups B/C are small-aperture-led and need contamination and independent-" +
    "extraction checks. The morphology-view period is formula-derived from the " +
    "BLS anchor and review fold factor; the fold is an evidence-view choice, " +
    "not automatically the physical orbital period.",
]];
overview.getRange("A12:H15").format = {
  fill: "#F7F8FA",
  font: { color: darkText, size: 10 },
  wrapText: true,
  verticalAlignment: "top",
  borders: { preset: "outside", style: "thin", color: lightGray },
};

overview.getRange("A17:B17").values = [["Sheet contract", "Value"]];
overview.getRange("A18:A23").values = [
  ["Sheet-set version"],
  ["Renderer"],
  ["Branch"],
  ["Apertures"],
  ["Trial periods"],
  ["Retained peaks"],
];
const sheetSet = manifest.sheet_set;
overview.getRange("B18:B23").values = [
  [sheetSet.sheet_set_version],
  [sheetSet.renderer_version],
  [sheetSet.branch_name],
  [sheetSet.apertures.join(", ")],
  [sheetSet.n_periods],
  [sheetSet.n_peaks],
];
overview.getRange("A17:B17").format = {
  fill: navy,
  font: { bold: true, color: "#FFFFFF" },
};
overview.getRange("A17:B23").format.borders = {
  preset: "outside",
  style: "thin",
  color: "#AAB1B9",
};
overview.getRange("A18:A23").format.fill = "#F4F6F8";
overview.getRange("B22:B23").format.numberFormat = "#,##0";
overview.getRange("A25:H26").merge();
overview.getRange("A25").values = [[
  "Current status: every row is Tier-1 searchable for bounded enrichment. " +
    "Pixel-level localization, archival checks, independent re-extraction, " +
    "false-alarm calibration, and follow-up readiness remain open.",
]];
overview.getRange("A25:H26").format = {
  fill: paleYellow,
  font: { color: darkText, italic: true },
  wrapText: true,
  verticalAlignment: "center",
  borders: { preset: "outside", style: "thin", color: "#D8B35A" },
};
overview.getRange("A1:H26").format.font.name = "Aptos";
overview.getRange("A:A").format.columnWidth = 35;
overview.getRange("B:B").format.columnWidth = 29;
overview.getRange("C:C").format.columnWidth = 3;
overview.getRange("D:D").format.columnWidth = 37;
overview.getRange("E:E").format.columnWidth = 12;
overview.getRange("F:H").format.columnWidth = 12;

const headers = [
  "Priority rank",
  "Group",
  "TIC",
  "Object",
  "Sector",
  "Discovery route",
  "Camera",
  "CCD",
  "Tmag",
  "BLS anchor P [d]",
  "Morphology fold",
  "Morphology view P [d]",
  "Duration [min]",
  "ADP-small depth",
  "ADP-small SDE",
  "ADP-primary SDE",
  "ADP period match",
  "Aperture ΔP/P",
  "Primary/small depth ratio",
  "ADP-small odd-even Δ [σ]",
  "ADP-primary odd-even Δ [σ]",
  "Max odd-even Δ [σ]",
  "Odd-even status",
  "Target QA",
  "Tier-1 searchable",
  "Priority note",
  "Vetting sheet",
  "Sheet SHA-256",
  "Source batch",
  "Review ID",
];
candidates.getRangeByIndexes(0, 0, 1, headers.length).values = [headers];
const values = rows.map((row) => [
  row.priority_rank,
  row.priority_group,
  String(row.tic),
  row.object_name,
  row.sector,
  row.discovery_route,
  row.camera,
  row.ccd,
  row.tmag,
  row.bls_anchor_period_d,
  row.morphology_fold_factor,
  null,
  row.duration_min,
  row.adp_small_depth_fraction,
  row.adp_small_sde,
  row.adp_primary_sde,
  row.aperture_period_match ? "yes" : "no",
  row.aperture_period_rel_delta,
  row.aperture_depth_ratio_primary_over_small,
  row.adp_small_odd_even_sigma,
  row.adp_primary_odd_even_sigma,
  null,
  row.odd_even_review_flag ? "review" : "clear",
  row.tier1_target_qa_status,
  row.tier1_target_searchable ? "yes" : "no",
  row.priority_note,
  row.vet_sheet_file,
  row.vet_sheet_sha256,
  row.source_batch_id,
  row.review_id,
]);
candidates.getRangeByIndexes(1, 0, values.length, headers.length).values = values;
for (let index = 0; index < rows.length; index += 1) {
  const excelRow = index + 2;
  candidates.getRange(`L${excelRow}`).formulas = [[`=J${excelRow}*K${excelRow}`]];
  candidates.getRange(`V${excelRow}`).formulas = [[
    `=MAX(ABS(T${excelRow}),ABS(U${excelRow}))`,
  ]];
}
const table = candidates.tables.add(`A1:AD${rows.length + 1}`, true, "S56PlanetLikeCandidates");
table.style = "TableStyleMedium2";
table.showFilterButton = true;
table.showBandedRows = true;
candidates.freezePanes.freezeRows(1);
candidates.freezePanes.freezeColumns(3);
candidates.showGridLines = false;
candidates.getRange("A1:AD1").format = {
  fill: navyDark,
  font: { bold: true, color: "#FFFFFF", size: 10 },
  wrapText: true,
  verticalAlignment: "center",
};
candidates.getRange("A1:AB1").format.rowHeight = 42;
candidates.getRange(`A2:A${rows.length + 1}`).format.numberFormat = "0";
candidates.getRange(`C2:C${rows.length + 1}`).format.numberFormat = "@";
candidates.getRange(`E2:H${rows.length + 1}`).format.numberFormat = "0";
candidates.getRange(`I2:I${rows.length + 1}`).format.numberFormat = "0.00";
candidates.getRange(`J2:L${rows.length + 1}`).format.numberFormat = "0.000000";
candidates.getRange(`M2:M${rows.length + 1}`).format.numberFormat = "0";
candidates.getRange(`N2:N${rows.length + 1}`).format.numberFormat = "0.0%";
candidates.getRange(`O2:P${rows.length + 1}`).format.numberFormat = "0.0";
candidates.getRange(`R2:R${rows.length + 1}`).format.numberFormat = "0.0000";
candidates.getRange(`S2:S${rows.length + 1}`).format.numberFormat = "0.00";
candidates.getRange(`T2:V${rows.length + 1}`).format.numberFormat = "0.00";
candidates.getRange(`Z2:Z${rows.length + 1}`).format.wrapText = true;
candidates.getRange(`A1:AD${rows.length + 1}`).format.font.name = "Aptos";
candidates.getRange(`A2:AD${rows.length + 1}`).format.rowHeight = 34;
candidates.getRange("A:A").format.columnWidth = 10;
candidates.getRange("B:B").format.columnWidth = 11;
candidates.getRange("C:C").format.columnWidth = 14;
candidates.getRange("D:D").format.columnWidth = 14;
candidates.getRange("E:E").format.columnWidth = 8;
candidates.getRange("F:F").format.columnWidth = 27;
candidates.getRange("G:H").format.columnWidth = 8;
candidates.getRange("I:I").format.columnWidth = 9;
candidates.getRange("J:L").format.columnWidth = 16;
candidates.getRange("M:M").format.columnWidth = 13;
candidates.getRange("N:N").format.columnWidth = 15;
candidates.getRange("O:P").format.columnWidth = 14;
candidates.getRange("Q:Q").format.columnWidth = 14;
candidates.getRange("R:V").format.columnWidth = 15;
candidates.getRange("W:Y").format.columnWidth = 16;
candidates.getRange("Z:Z").format.columnWidth = 55;
candidates.getRange("AA:AA").format.columnWidth = 43;
candidates.getRange("AB:AB").format.columnWidth = 67;
candidates.getRange("AC:AD").format.columnWidth = 34;
candidates.getRange(`B2:B${rows.length + 1}`).conditionalFormats.add(
  "containsText",
  { text: "Benchmark", format: { fill: paleBlue, font: { bold: true } } },
);
candidates.getRange(`B2:B${rows.length + 1}`).conditionalFormats.add(
  "cellIs",
  { operator: "equal", formula: '"A"', format: { fill: paleGreen, font: { bold: true } } },
);
candidates.getRange(`B2:B${rows.length + 1}`).conditionalFormats.add(
  "cellIs",
  { operator: "equal", formula: '"B"', format: { fill: paleYellow, font: { bold: true } } },
);
candidates.getRange(`B2:B${rows.length + 1}`).conditionalFormats.add(
  "cellIs",
  { operator: "equal", formula: '"C"', format: { fill: paleRed, font: { bold: true } } },
);
candidates.getRange(`Q2:Q${rows.length + 1}`).conditionalFormats.add(
  "containsText",
  { text: "no", format: { fill: paleYellow } },
);
candidates.getRange(`W2:W${rows.length + 1}`).conditionalFormats.add(
  "containsText",
  { text: "review", format: { fill: paleRed } },
);
candidates.getRange(`X2:X${rows.length + 1}`).conditionalFormats.add(
  "containsText",
  { text: "review", format: { fill: paleYellow } },
);
candidates.getRange(`Y2:Y${rows.length + 1}`).conditionalFormats.add(
  "containsText",
  { text: "no", format: { fill: paleRed } },
);

provenance.showGridLines = false;
provenance.getRange("A1:D1").merge();
provenance.getRange("A1").values = [["Candidate bundle provenance"]];
provenance.getRange("A1:D1").format = {
  fill: navyDark,
  font: { bold: true, color: "#FFFFFF", size: 16 },
};
provenance.getRange("A3:D3").values = [["Artifact", "Role", "SHA-256", "Path"]];
const sourceRows = Object.entries(manifest.sources).map(([name, entry]) => [
  name,
  "source",
  entry.sha256,
  entry.path,
]);
const outputRows = Object.entries(manifest.outputs)
  .filter(([name]) => name !== path.basename(outputXlsx))
  .map(([name, entry]) => [name, "output", entry.sha256, entry.path]);
const provenanceRows = [...sourceRows, ...outputRows];
provenance.getRangeByIndexes(3, 0, provenanceRows.length, 4).values = provenanceRows;
provenance.getRange("A3:D3").format = {
  fill: navy,
  font: { bold: true, color: "#FFFFFF" },
};
provenance.getRange(`A3:D${provenanceRows.length + 3}`).format.font.name = "Aptos";
provenance.getRange(`A4:D${provenanceRows.length + 3}`).format.rowHeight = 28;
provenance.getRange("A:A").format.columnWidth = 40;
provenance.getRange("B:B").format.columnWidth = 12;
provenance.getRange("C:C").format.columnWidth = 67;
provenance.getRange("D:D").format.columnWidth = 95;
provenance.freezePanes.freezeRows(3);

const overviewCheck = await workbook.inspect({
  kind: "table",
  range: "Overview!A1:H26",
  include: "values,formulas",
  tableMaxRows: 26,
  tableMaxCols: 8,
  maxChars: 8000,
});
const candidateCheck = await workbook.inspect({
  kind: "table",
  range: "Candidates!A1:Y8",
  include: "values,formulas",
  tableMaxRows: 8,
  tableMaxCols: 25,
  maxChars: 8000,
});
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
  maxChars: 3000,
});
console.log(overviewCheck.ndjson);
console.log(candidateCheck.ndjson);
console.log(errors.ndjson);

await fs.mkdir(path.dirname(outputXlsx), { recursive: true });
await fs.mkdir(previewDir, { recursive: true });
for (const [sheetName, range] of [
  ["Overview", "A1:H26"],
  ["Candidates", "A1:Z31"],
  ["Provenance", `A1:D${provenanceRows.length + 3}`],
]) {
  const preview = await workbook.render({
    sheetName,
    range,
    scale: sheetName === "Candidates" ? 0.7 : 1.2,
    format: "png",
  });
  const previewBytes = new Uint8Array(await preview.arrayBuffer());
  await fs.writeFile(path.join(previewDir, `${sheetName.toLowerCase()}.png`), previewBytes);
}
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputXlsx);
const workbookBytes = await fs.readFile(outputXlsx);
manifest.outputs[path.basename(outputXlsx)] = {
  path: outputXlsx,
  sha256: createHash("sha256").update(workbookBytes).digest("hex"),
  size_bytes: workbookBytes.byteLength,
};
await fs.writeFile(manifestJson, `${JSON.stringify(manifest, null, 2)}\n`);
console.log(JSON.stringify({ outputXlsx, previewDir }, null, 2));
