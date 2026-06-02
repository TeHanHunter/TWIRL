
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
