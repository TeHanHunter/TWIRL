#!/usr/bin/env python3
"""Render a self-contained, reader-oriented report for the audited Teacher v2 run."""

from __future__ import annotations

import json
import sys
from pathlib import Path

DEPS = Path(__file__).resolve().parent / "_deps"
if DEPS.exists():
    sys.path.insert(0, str(DEPS))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))

from twirl.plotting.style import apply_twirl_style, get_ordered_palette  # noqa: E402


REPORT_DIR = Path(__file__).resolve().parent
ARTIFACT_DIR = REPORT_DIR.parent / "orcd_final_artifacts"
FIG_DIR = REPORT_DIR / "figures"
PDF_PATH = REPORT_DIR / "s56_teacher_v2_full_report.pdf"

NAVY = "#17365D"
TEAL = "#2F7F6F"
ORANGE = "#B8582F"
BLUE = "#3B6EA5"
SLATE = "#5E6670"
PURPLE = "#7957A6"
LIGHT_BLUE = "#E8F0F8"
LIGHT_TEAL = "#E8F3F0"
LIGHT_ORANGE = "#F9EEE8"
LIGHT_RED = "#F8E8E8"
CLASS_LABELS = ["Planet", "Eclipse", "Variable", "Other"]


def read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def save_figure(fig: plt.Figure, name: str) -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return path


def clean_profile(name: str) -> str:
    return {
        "metadata_only": "Metadata only",
        "single_period_native_fold": "Single fold",
        "seven_harmonic_shape": "Seven-harmonic shape",
        "shape_plus_raw_chronology": "Shape + raw chronology",
        "shape_plus_periodogram_bls": "Shape + periodogram/BLS",
        "full_combined": "Full combined",
    }[name]


def make_pipeline_figure() -> Path:
    apply_twirl_style("full_page")
    fig, ax = plt.subplots(figsize=(7.1, 2.45))
    ax.set_axis_off()
    nodes = [
        (0.02, 0.33, 0.16, 0.38, "A2v1 ADP\nlight curves", "#EDF2F7"),
        (0.23, 0.33, 0.16, 0.38, "ADP-small BLS\n(top 5 peaks)", "#E8F0F8"),
        (0.44, 0.33, 0.18, 0.38, "Teacher v2 ensemble\ncompact score", "#E8F3F0"),
        (0.67, 0.33, 0.14, 0.38, "Frozen threshold\np >= 0.72064", "#F9EEE8"),
        (0.86, 0.33, 0.12, 0.38, "Human\nreview", "#F8E8E8"),
    ]
    for x, y, width, height, label, face in nodes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0.012,rounding_size=0.015",
                linewidth=0.85,
                edgecolor="0.25",
                facecolor=face,
                transform=ax.transAxes,
            )
        )
        ax.text(
            x + width / 2,
            y + height / 2,
            label,
            ha="center",
            va="center",
            fontsize=8,
            transform=ax.transAxes,
        )
    for idx in range(len(nodes) - 1):
        left = nodes[idx]
        right = nodes[idx + 1]
        ax.annotate(
            "",
            xy=(right[0] - 0.01, 0.52),
            xytext=(left[0] + left[2] + 0.01, 0.52),
            xycoords=ax.transAxes,
            arrowprops={"arrowstyle": "->", "color": "0.3", "lw": 1.0},
        )
    ax.text(
        0.5,
        0.10,
        "The teacher is applied only after BLS proposes an ephemeris. It ranks a manageable review set; it does not discover signals that BLS missed.",
        ha="center",
        va="center",
        fontsize=7,
        color="0.28",
        transform=ax.transAxes,
    )
    return save_figure(fig, "teacher_v2_pipeline")


def make_profile_comparison(selection: dict) -> Path:
    ranking = pd.DataFrame(selection["profile_ranking"])
    ranking = ranking.iloc[::-1].reset_index(drop=True)
    labels = [clean_profile(value) for value in ranking["profile"]]
    selected = ranking["profile"].eq(selection["selected_profile"])
    colors_by_row = [TEAL if value else "#AEB8C2" for value in selected]

    apply_twirl_style("full_page")
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 4.4), sharey=True)
    metrics = [
        ("development_injection_recall", "Injection retention\nat fixed 5% workload", "{:.0%}"),
        ("real_morphology_macro_f1", "Human morphology\nmacro F1", "{:.3f}"),
        ("compact_ece", "Compact-score ECE\n(lower is better)", "{:.3f}"),
    ]
    y = np.arange(len(ranking))
    for axis, (column, title, formatting) in zip(axes, metrics, strict=True):
        values = ranking[column].to_numpy(dtype=float)
        bars = axis.barh(y, values, color=colors_by_row, edgecolor="white", linewidth=0.5)
        axis.set_title(title, fontsize=8.2, pad=6)
        axis.set_yticks(y, labels if axis is axes[0] else [])
        axis.grid(axis="x", alpha=0.75)
        axis.set_axisbelow(True)
        if column == "development_injection_recall":
            axis.set_xlim(0.45, 0.76)
        elif column == "real_morphology_macro_f1":
            axis.set_xlim(0.65, 0.84)
        else:
            axis.set_xlim(0.0, 0.21)
        for bar, value in zip(bars, values, strict=True):
            axis.text(
                bar.get_width() + 0.006 * (axis.get_xlim()[1] - axis.get_xlim()[0]),
                bar.get_y() + bar.get_height() / 2,
                formatting.format(value),
                va="center",
                fontsize=6.4,
            )
    axes[0].set_xlabel("fraction retained")
    axes[1].set_xlabel("macro-average F1")
    axes[2].set_xlabel("calibration error")
    fig.text(
        0.01,
        0.02,
        "Development-only architecture comparison. Green marks the frozen Teacher v2 profile. The selection weighted injected compact-transit retention most heavily, then real morphology, compact AP, and calibration.",
        fontsize=6.6,
        color="0.30",
    )
    fig.subplots_adjust(left=0.25, right=0.98, bottom=0.18, top=0.84, wspace=0.40)
    return save_figure(fig, "development_profile_comparison")


def make_training_curves() -> Path:
    histories: list[pd.DataFrame] = []
    for fold in range(5):
        path = ARTIFACT_DIR / "training" / "shape_plus_raw_chronology" / f"fold_{fold}" / "history.csv"
        frame = pd.read_csv(path)
        frame["fold"] = fold + 1
        histories.append(frame)
    palette = get_ordered_palette(5, "viridis")
    apply_twirl_style("full_page")
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.55))
    columns = [
        ("train_loss", "Training loss", "Lower means the training objective is being fit"),
        ("validation_macro_f1", "Validation macro F1", "Higher means more even four-class accuracy"),
        ("validation_compact_ap", "Validation compact AP", "Higher means better compact-transit ranking"),
    ]
    for axis, (column, label, subtitle) in zip(axes, columns, strict=True):
        for frame, color in zip(histories, palette, strict=True):
            axis.plot(frame["epoch"], frame[column], color=color, linewidth=1.0, alpha=0.9, label=f"Fold {frame['fold'].iloc[0]}")
        axis.set_title(label, fontsize=8.2, pad=5)
        axis.set_xlabel("epoch")
        axis.set_ylabel("value")
        axis.text(0.02, 0.03, subtitle, transform=axis.transAxes, fontsize=5.9, color="0.33", va="bottom")
        axis.grid(alpha=0.65)
    axes[0].legend(loc="upper right", ncol=1, fontsize=5.4, frameon=True)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.22, top=0.84, wspace=0.33)
    return save_figure(fig, "selected_profile_training_curves")


def make_calibration_figure(selection: dict) -> Path:
    profile = selection["profile_details"][selection["selected_profile"]]["development_metrics"]
    bins = pd.DataFrame(profile["compact_calibration"]["bins"])
    bins = bins.loc[bins["n"] > 0].copy()
    midpoint = (bins["low"] + bins["high"]) / 2
    ece = profile["compact_calibration"]["ece"]
    apply_twirl_style("full_page")
    fig, ax = plt.subplots(figsize=(4.0, 3.25))
    ax.plot([0.5, 1.0], [0.5, 1.0], color="0.4", linestyle=(0, (3, 2)), linewidth=1.0, label="perfect calibration")
    sizes = 15 + 0.010 * bins["n"].to_numpy(dtype=float)
    ax.scatter(midpoint, bins["accuracy"], s=sizes, color=TEAL, edgecolor="white", linewidth=0.6, zorder=3, label="development bins")
    for x, y, n in zip(midpoint, bins["accuracy"], bins["n"], strict=True):
        ax.annotate(f"n={int(n):,}", (x, y), xytext=(4, 4), textcoords="offset points", fontsize=5.8, color="0.30")
    ax.set_xlim(0.5, 1.0)
    ax.set_ylim(0.5, 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("mean predicted compact-transit probability")
    ax.set_ylabel("observed compact-positive fraction")
    ax.legend(loc="upper left", fontsize=6.2)
    ax.text(0.53, 0.94, f"ECE = {ece:.3f}", fontsize=7.4, color=NAVY, fontweight="bold")
    fig.tight_layout(pad=0.4)
    return save_figure(fig, "selected_profile_calibration")


def draw_matrix(axis: plt.Axes, matrix: list[list[int]], title: str) -> None:
    image = axis.imshow(matrix, cmap="Blues", vmin=0, vmax=max(max(row) for row in matrix))
    for row_index, row in enumerate(matrix):
        for col_index, value in enumerate(row):
            color = "white" if value > 0.55 * image.norm.vmax else "0.15"
            axis.text(col_index, row_index, str(value), ha="center", va="center", color=color, fontsize=7.6)
    axis.set_xticks(range(4), CLASS_LABELS, rotation=32, ha="right")
    axis.set_yticks(range(4), CLASS_LABELS)
    axis.set_xlabel("predicted label")
    axis.set_ylabel("human label")
    axis.set_title(title, fontsize=8.4, pad=6)


def make_human_holdout_figure(human: dict) -> Path:
    v1 = human["teacher_v1"]
    v2 = human["teacher_v2"]
    apply_twirl_style("full_page")
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 3.1), gridspec_kw={"width_ratios": [1, 1, 1.32]})
    draw_matrix(axes[0], v1["confusion_matrix"], "Teacher v1")
    draw_matrix(axes[1], v2["confusion_matrix"], "Teacher v2")
    order = ["planet_like", "eclipse_contact", "smooth_variable", "other"]
    x = np.arange(len(order))
    width = 0.37
    v1_recall = [v1["per_class"][key]["recall"] for key in order]
    v2_recall = [v2["per_class"][key]["recall"] for key in order]
    v1_n = [v1["per_class"][key]["n"] for key in order]
    axes[2].bar(x - width / 2, v1_recall, width, color=SLATE, label="Teacher v1")
    axes[2].bar(x + width / 2, v2_recall, width, color=TEAL, label="Teacher v2")
    axes[2].set_ylim(0, 1.08)
    axes[2].set_xticks(x, [f"{label}\n(n={count})" for label, count in zip(CLASS_LABELS, v1_n, strict=True)])
    axes[2].set_ylabel("recall")
    axes[2].set_title("Recall on the same real human holdout", fontsize=8.4, pad=6)
    axes[2].legend(loc="upper left", fontsize=6.1)
    for index, (left, right) in enumerate(zip(v1_recall, v2_recall, strict=True)):
        axes[2].text(index - width / 2, left + 0.035, f"{left:.0%}", ha="center", fontsize=5.6, color=SLATE)
        axes[2].text(index + width / 2, right + 0.035, f"{right:.0%}", ha="center", fontsize=5.6, color=TEAL)
    fig.text(
        0.01,
        0.01,
        "Rows are human labels and columns are model labels. The high Other count makes accuracy look strong; macro F1 and individual recalls are more informative. Planet results are descriptive only because the holdout has n=2 real planets.",
        fontsize=6.3,
        color="0.30",
    )
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.24, top=0.85, wspace=0.50)
    return save_figure(fig, "same_row_human_holdout")


def make_external_comparison(s56: dict, s57: dict) -> Path:
    def frame(summary: dict) -> pd.DataFrame:
        table = pd.DataFrame(summary["model_comparison"])
        return table.set_index("model").loc[["teacher_v1", "metadata_only", "teacher_v2"]].reset_index()

    left = frame(s56)
    right = frame(s57)
    labels = ["Teacher v1", "Metadata only", "Teacher v2"]
    colors_by_model = [SLATE, ORANGE, TEAL]
    apply_twirl_style("full_page")
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.2), sharey=True)
    for axis, table, title in zip(axes, [left, right], ["S56 locked injection holdout", "S57 external injection set"], strict=True):
        values = table["conditional_on_bls"].to_numpy(dtype=float)
        bars = axis.bar(labels, values, color=colors_by_model, width=0.68)
        axis.set_ylim(0, 1.0)
        axis.set_ylabel("Teacher retention given BLS recovery" if axis is axes[0] else "")
        axis.set_title(title, fontsize=8.4, pad=6)
        axis.tick_params(axis="x", rotation=22)
        for bar, value in zip(bars, values, strict=True):
            axis.text(bar.get_x() + bar.get_width() / 2, value + 0.025, f"{value:.1%}", ha="center", va="bottom", fontsize=6.7)
        axis.axhline(0.80, color="0.35", linestyle=(0, (3, 2)), linewidth=0.9)
        axis.text(2.42, 0.812, "80% gate", fontsize=5.8, color="0.35", ha="right")
    fig.text(
        0.01,
        0.02,
        "Each bar asks: among injections for which BLS put a truth-matched peak in its top five, what fraction did the teacher retain at the frozen S56 review-workload threshold? This isolates teacher triage from BLS detection.",
        fontsize=6.4,
        color="0.30",
    )
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.25, top=0.84, wspace=0.20)
    return save_figure(fig, "s56_s57_teacher_comparison")


def make_tmag_figure() -> Path:
    s56 = pd.read_csv(ARTIFACT_DIR / "s56_locked_holdout_recovery" / "teacher_retention_by_tmag.csv")
    s57 = pd.read_csv(ARTIFACT_DIR / "s57_external_recovery" / "teacher_retention_by_tmag.csv")
    apply_twirl_style("full_page")
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.2), sharey=True)
    width = 0.36
    x = np.arange(len(s56))
    short = ["<17", "17-18", "18-19", ">=19"]
    for axis, table, title in zip(axes, [s56, s57], ["S56", "S57"], strict=True):
        bls = table["bls_top5_fraction"].to_numpy(dtype=float)
        retained = table["teacher_v2_retention_given_bls"].to_numpy(dtype=float)
        bars_a = axis.bar(x - width / 2, bls, width, color=BLUE, label="BLS top-5 recovery")
        bars_b = axis.bar(x + width / 2, retained, width, color=TEAL, label="Teacher v2 retention | BLS")
        axis.set_title(f"{title} by Tmag bin", fontsize=8.4, pad=6)
        axis.set_xticks(x, [f"{label}\n(n={n:,})" for label, n in zip(short, table["n_injections"], strict=True)])
        axis.set_ylim(0, 1.08)
        axis.set_ylabel("fraction" if axis is axes[0] else "")
        for bar in list(bars_a) + list(bars_b):
            height = bar.get_height()
            axis.text(bar.get_x() + bar.get_width() / 2, height + 0.024, f"{height:.0%}", ha="center", va="bottom", fontsize=5.4, rotation=90)
    axes[0].legend(loc="upper left", fontsize=5.8)
    fig.text(
        0.01,
        0.02,
        "BLS recovery declines sharply toward the faintest hosts. Teacher retention after BLS remains substantially higher, but bright-bin map cells are sparse because hosts were not reused; interpret those bins as support-limited.",
        fontsize=6.4,
        color="0.30",
    )
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.24, top=0.84, wspace=0.18)
    return save_figure(fig, "tmag_stratified_recovery")


def image_flowable(path: Path, max_width: float, max_height: float) -> Image:
    with PILImage.open(path) as image:
        width, height = image.size
    scale = min(max_width / width, max_height / height)
    return Image(str(path), width=width * scale, height=height * scale)


def caption(text: str, styles: dict) -> Paragraph:
    return Paragraph(text, styles["caption"])


def add_figure(story: list, path: Path, title: str, text: str, styles: dict, *, height: float = 8.4 * inch) -> None:
    story.append(Paragraph(title, styles["h2"]))
    story.append(image_flowable(path, 7.1 * inch, height))
    story.append(caption(text, styles))
    story.append(PageBreak())


def make_styles() -> dict:
    sample = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "ReportTitle",
            parent=sample["Title"],
            fontName="Times-Bold",
            fontSize=22,
            leading=26,
            textColor=HexColor(NAVY),
            alignment=TA_LEFT,
            spaceAfter=8,
        ),
        "subtitle": ParagraphStyle(
            "Subtitle",
            parent=sample["Normal"],
            fontName="Times-Roman",
            fontSize=11,
            leading=14,
            textColor=HexColor("#364152"),
            spaceAfter=15,
        ),
        "h1": ParagraphStyle(
            "HeadingOne",
            parent=sample["Heading1"],
            fontName="Times-Bold",
            fontSize=15,
            leading=18,
            textColor=HexColor(NAVY),
            spaceBefore=0,
            spaceAfter=8,
        ),
        "h2": ParagraphStyle(
            "HeadingTwo",
            parent=sample["Heading2"],
            fontName="Times-Bold",
            fontSize=12,
            leading=15,
            textColor=HexColor(NAVY),
            spaceBefore=0,
            spaceAfter=6,
        ),
        "body": ParagraphStyle(
            "Body",
            parent=sample["BodyText"],
            fontName="Times-Roman",
            fontSize=9.2,
            leading=12.4,
            textColor=HexColor("#20252B"),
            spaceAfter=7,
        ),
        "small": ParagraphStyle(
            "Small",
            parent=sample["BodyText"],
            fontName="Times-Roman",
            fontSize=7.8,
            leading=10.0,
            textColor=HexColor("#303842"),
            spaceAfter=4,
        ),
        "caption": ParagraphStyle(
            "Caption",
            parent=sample["BodyText"],
            fontName="Times-Italic",
            fontSize=7.4,
            leading=9.4,
            textColor=HexColor("#3E4650"),
            spaceBefore=5,
            spaceAfter=4,
        ),
        "callout": ParagraphStyle(
            "Callout",
            parent=sample["BodyText"],
            fontName="Times-Roman",
            fontSize=9.0,
            leading=11.5,
            textColor=HexColor("#20252B"),
            alignment=TA_LEFT,
        ),
    }


def footer(canvas, document) -> None:
    canvas.saveState()
    canvas.setStrokeColor(HexColor("#AAB2BD"))
    canvas.setLineWidth(0.35)
    canvas.line(0.7 * inch, 0.58 * inch, 7.8 * inch, 0.58 * inch)
    canvas.setFont("Times-Roman", 7)
    canvas.setFillColor(HexColor("#4B5563"))
    canvas.drawString(0.7 * inch, 0.40 * inch, "TWIRL S56 Teacher v2 report | Internal validation")
    canvas.drawRightString(7.8 * inch, 0.40 * inch, f"Page {document.page}")
    canvas.restoreState()


def first_page_footer(canvas, document) -> None:
    canvas.saveState()
    canvas.setFont("Times-Roman", 7)
    canvas.setFillColor(HexColor("#4B5563"))
    canvas.drawString(0.7 * inch, 0.42 * inch, "TWIRL S56 Teacher v2 report | Generated from audited ORCD artifacts")
    canvas.restoreState()


def build_pdf(
    figures: dict[str, Path],
    selection: dict,
    human: dict,
    s56: dict,
    s57: dict,
    roles: dict,
    training_summary: dict,
) -> None:
    styles = make_styles()
    doc = SimpleDocTemplate(
        str(PDF_PATH),
        pagesize=letter,
        rightMargin=0.7 * inch,
        leftMargin=0.7 * inch,
        topMargin=0.64 * inch,
        bottomMargin=0.72 * inch,
        title="TWIRL S56 Teacher v2 full validation report",
        author="TWIRL",
    )
    story: list = []

    story.append(Spacer(1, 0.2 * inch))
    story.append(Paragraph("S56 Teacher v2", styles["title"]))
    story.append(Paragraph("Full validation report: compact-transit triage, human morphology, and external-sector recovery", styles["subtitle"]))
    story.append(Paragraph("Decision in one sentence", styles["h1"]))
    story.append(
        Paragraph(
            "Teacher v2 is a substantially better <b>compact-transit ranking model</b> than Teacher v1, but it is not yet a reliable four-class morphology classifier and is therefore approved only for active-learning enrichment, not automatic student pseudo-labeling.",
            styles["body"],
        )
    )
    callouts = [
        [
            Paragraph("<b>External compact-triage result</b><br/>On untouched S57 injections, Teacher v2 retained <b>77.21%</b> of BLS-recovered truth-matched candidates at the frozen 5% real-TIC workload.", styles["callout"]),
            Paragraph("<b>Relative gain</b><br/>That is <b>+20.96 percentage points</b> over Teacher v1 and <b>+15.27 points</b> over the metadata-only baseline on the same S57 candidate set.", styles["callout"]),
        ],
        [
            Paragraph("<b>Promotion gate missed</b><br/>The predeclared compact-triage target was 80%; 77.21% misses it by <b>2.79 points</b> (about 78 of 2,777 BLS-recovered injections).", styles["callout"]),
            Paragraph("<b>Real morphology warning</b><br/>On the same human holdout, macro F1 fell from <b>0.725</b> in v1 to <b>0.653</b> in v2; eclipse/contact recall is only <b>36.4%</b> (4/11).", styles["callout"]),
        ],
    ]
    table = Table(callouts, colWidths=[3.48 * inch, 3.48 * inch], rowHeights=[0.86 * inch, 0.86 * inch], hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, 0), HexColor(LIGHT_TEAL)),
                ("BACKGROUND", (1, 0), (1, 0), HexColor(LIGHT_BLUE)),
                ("BACKGROUND", (0, 1), (0, 1), HexColor(LIGHT_ORANGE)),
                ("BACKGROUND", (1, 1), (1, 1), HexColor(LIGHT_RED)),
                ("BOX", (0, 0), (-1, -1), 0.45, HexColor("#B0B8C1")),
                ("INNERGRID", (0, 0), (-1, -1), 0.35, HexColor("#B0B8C1")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 0.16 * inch))
    story.append(Paragraph("What this report covers", styles["h1"]))
    story.append(
        Paragraph(
            "The report separates the two scientific questions that were tested: <b>(1) compact-transit triage</b>, evaluated with fresh injected signals after BLS, and <b>(2) morphology classification</b>, evaluated with real human labels. A result can be good for one question and not the other. The recovery-map appendix contains every final S56 and S57 publication-style map: BLS alone, Teacher v1, Teacher v2, and the v2-minus-v1 difference map.",
            styles["body"],
        )
    )
    story.append(Spacer(1, 0.16 * inch))
    story.append(Paragraph("Run identity", styles["h2"]))
    run_table = Table(
        [
            ["Frozen profile", "Shape + raw chronology"],
            ["Model version", selection["model_version"]],
            ["Frozen compact threshold", f"p_compact = {selection['frozen_compact_threshold']:.6f}"],
            ["Review workload", "5% of 6,858 held-out real S56 TICs (342 TICs)"],
            ["Training rows / unique TICs", f"{roles['n_training_rows']:,} / {training_summary['training_verification']['n_unique_tics']:,}"],
            ["Final code checkpoint", "8f1dad40 on codex/s56-teacher-v2"],
        ],
        colWidths=[2.05 * inch, 4.91 * inch],
    )
    run_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, -1), HexColor("#EFF3F7")),
                ("BOX", (0, 0), (-1, -1), 0.4, HexColor("#AAB2BD")),
                ("INNERGRID", (0, 0), (-1, -1), 0.3, HexColor("#C2C8CE")),
                ("FONTNAME", (0, 0), (0, -1), "Times-Bold"),
                ("FONTNAME", (1, 0), (1, -1), "Times-Roman"),
                ("FONTSIZE", (0, 0), (-1, -1), 8.5),
                ("LEADING", (0, 0), (-1, -1), 10.5),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 7),
                ("RIGHTPADDING", (0, 0), (-1, -1), 7),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.append(run_table)
    story.append(PageBreak())

    story.append(Paragraph("How Teacher v2 is evaluated", styles["h1"]))
    story.append(
        Paragraph(
            "The classifier is deliberately downstream of BLS. BLS supplies up to five candidate periods per target. Teacher v2 then ranks those candidate folds using the frozen compact-transit probability. The primary recovery metric therefore has two stages: BLS must first put a truth-matched signal among its top five, then Teacher v2 must retain that candidate. The teacher cannot recover an injection that BLS never proposed.",
            styles["body"],
        )
    )
    story.append(image_flowable(figures["pipeline"], 7.1 * inch, 2.6 * inch))
    story.append(caption("BLS provides the ephemeris proposal. The frozen threshold was chosen before opening the locked S56 holdout and the independent S57 injection set.", styles))
    story.append(Spacer(1, 0.05 * inch))
    story.append(Paragraph("The S57 external funnel", styles["h2"]))
    funnel_table = Table(
        [
            ["Stage", "Count", "Fraction of previous stage", "What it means"],
            ["Fresh S57 injections", f"{s57['n_injections']:,}", "100%", "One unique host TIC per injection; none used to choose the v2 architecture or threshold."],
            ["BLS top-5 truth matches", f"{s57['teacher_v2']['n_bls_top5_recovered']:,}", f"{s57['teacher_v2']['bls_top5_fraction']:.2%}", "BLS proposed the correct period or an allowed harmonic among its first five candidates."],
            ["Teacher v2 retained", f"{s57['teacher_v2']['n_compact_recovered']:,}", f"{s57['teacher_v2']['compact_conditional_on_bls']:.2%}", "The frozen compact score passed the 5% real-TIC workload threshold."],
            ["End-to-end recovery", f"{s57['teacher_v2']['n_compact_recovered']:,}", f"{s57['teacher_v2']['compact_end_to_end_fraction']:.2%}", "Combined BLS plus Teacher v2 recovery, relative to all injected signals."],
        ],
        colWidths=[1.65 * inch, 0.9 * inch, 1.3 * inch, 3.11 * inch],
        repeatRows=1,
    )
    funnel_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), HexColor(NAVY)),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Times-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Times-Roman"),
                ("FONTSIZE", (0, 0), (-1, -1), 7.8),
                ("LEADING", (0, 0), (-1, -1), 9.5),
                ("BOX", (0, 0), (-1, -1), 0.4, HexColor("#9BA5AF")),
                ("INNERGRID", (0, 0), (-1, -1), 0.3, HexColor("#C2C8CE")),
                ("BACKGROUND", (0, 1), (-1, -1), HexColor("#FAFBFC")),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.append(funnel_table)
    story.append(Spacer(1, 0.10 * inch))
    story.append(Paragraph("What was and was not used as a label", styles["h2"]))
    story.append(
        Paragraph(
            "For injected light curves, injection truth defined whether a BLS candidate was a compact-transit positive. For real candidates, the human label remained the morphology target: Planet-like, Eclipse/contact, Smooth variable, or Other. Injection truth, recovery flags, human notes, adjudicated harmonic factors, cohort identity, and paths were excluded from model inputs. Franklin's labels were snapshot as audit data but were not used in the initial v2 fit because A2v1 transfer review was deferred.",
            styles["body"],
        )
    )
    story.append(PageBreak())

    add_figure(
        story,
        figures["profiles"],
        "Architecture comparison on the S56 development data",
        "Each profile was trained as a five-fold TIC-grouped ensemble and evaluated at the same real-candidate workload. The selected profile did not simply have the best morphology score: the selection intentionally prioritized compact-transit retention, because this experiment's immediate use is candidate enrichment.",
        styles,
        height=4.8 * inch,
    )
    add_figure(
        story,
        figures["training"],
        "Selected-profile training trajectories",
        "A fold is one rotation of grouped cross-validation: different TIC groups are held out for validation in each line. Training loss should decline; validation macro F1 and compact average precision show whether that learning generalizes. Variability between folds is expected because rare real classes are unevenly distributed, which is precisely why the model must be judged on the locked holdouts rather than on these curves alone.",
        styles,
        height=4.6 * inch,
    )
    add_figure(
        story,
        figures["calibration"],
        "Compact-score calibration on the S56 development set",
        "Calibration asks whether a score means what it says. Points close to the diagonal imply that candidates assigned, for example, about 80% compact probability were compact positives about 80% of the time in this development mixture. ECE is the average calibration mismatch; 0.015 is low. This is useful for score ranking, but it does not prove equal calibration on real astrophysical candidates.",
        styles,
        height=4.8 * inch,
    )
    add_figure(
        story,
        figures["human"],
        "Same-row real human holdout: morphology remains the weak point",
        "This is a direct, fair v1-versus-v2 comparison on the same 268 real human-labeled rows. Overall accuracy is near 94% mainly because Other has 233 rows. Macro F1 gives every class equal influence and exposes the v2 regression, especially for Eclipse/contact. The two real Planet-like examples are far too few to estimate real-planet completeness.",
        styles,
        height=4.8 * inch,
    )

    story.append(Paragraph("S56 locked injection holdout", styles["h1"]))
    story.append(
        Paragraph(
            f"S56 provides the initial locked test: {s56['n_injections']:,} fresh injections on host TICs excluded from the fitting and development rows. BLS recovered {s56['teacher_v2']['n_bls_top5_recovered']:,} ({s56['teacher_v2']['bls_top5_fraction']:.2%}) in its top five. Teacher v2 kept {s56['teacher_v2']['n_compact_recovered']:,} of those, or {s56['teacher_v2']['compact_conditional_on_bls']:.2%} conditional retention. This is an internal locked test; S57 is the stricter external transfer check.",
            styles["body"],
        )
    )
    add_figure(
        story,
        ARTIFACT_DIR / "s56_locked_holdout_recovery" / "bls_top5_recovery_map" / "period_radius_empirical_recovery_publication.png",
        "S56: ADP-small BLS top-five recovery",
        "This map answers only whether the search proposed a truth-matched period or allowed harmonic among its first five candidates. It is the ceiling for the downstream teacher under this BLS-first workflow. Each panel is a Tmag range; sparse bright panels show limited support because hosts were not reused.",
        styles,
    )
    add_figure(
        story,
        ARTIFACT_DIR / "s56_locked_holdout_recovery" / "teacher_v2_compact_recovery_map" / "period_radius_empirical_recovery_publication.png",
        "S56: ADP BLS plus Teacher v2 compact-transit retention",
        "This uses the same recovered BLS candidate set as the preceding map, then applies the frozen compact threshold. Read it as the end-to-end recovery of the BLS-plus-teacher workflow, not as a standalone neural-network completeness map.",
        styles,
    )
    add_figure(
        story,
        ARTIFACT_DIR / "s56_locked_holdout_recovery" / "teacher_v1_workload_matched_recovery_map" / "period_radius_empirical_recovery_publication.png",
        "S56: Teacher v1 at the matched workload",
        "Teacher v1 is evaluated on the identical BLS candidate set at a workload-matched threshold. This allows map-level comparison without giving either model more human-review capacity.",
        styles,
    )
    add_figure(
        story,
        ARTIFACT_DIR / "s56_locked_holdout_recovery" / "teacher_v2_minus_v1_difference_map" / "teacher_v2_minus_v1_recovery.png",
        "S56: Teacher v2 minus Teacher v1",
        "Positive colors indicate period-radius cells where v2 retains more recovered injected compact transits than v1 at the same review workload. This figure is a difference diagnostic; it is not an absolute recovery map.",
        styles,
    )

    story.append(Paragraph("Workload and magnitude dependence", styles["h1"]))
    story.append(
        Paragraph(
            "Workload is a practical control on the number of real targets sent to a human. Increasing it improves injected retention, but it also increases the review burden and false-positive load. The frozen 5% setting was deliberately conservative for the external comparison; it should be treated as a fixed benchmark, not as the only usable operational setting.",
            styles["body"],
        )
    )
    story.append(image_flowable(REPORT_DIR.parent / "s56_locked_holdout_recovery" / "workload_performance" / "recovery_vs_review_workload.png", 7.1 * inch, 3.7 * inch))
    story.append(caption("S56 locked-holdout compact-transit retention as a function of review workload. The selected 5% reference point preserves 78.55% of BLS-recovered injections; 10%, 20%, and 30% workloads retain 82.35%, 86.73%, and 88.98%, respectively.", styles))
    story.append(Spacer(1, 0.08 * inch))
    story.append(image_flowable(figures["tmag"], 7.1 * inch, 3.25 * inch))
    story.append(caption("Tmag-stratified results. BLS is the dominant faint-end bottleneck. Teacher-v2 retention conditional on BLS stays high compared with BLS recovery itself, while sparse bright-bin support is explicitly flagged rather than made artificially dense through host reuse.", styles))
    story.append(PageBreak())

    story.append(Paragraph("S57 external injection evaluation", styles["h1"]))
    story.append(
        Paragraph(
            f"S57 is the external transfer test. Its {s57['n_injections']:,} injections were never used to choose the architecture, calibration, or operating threshold. The selected profile retained {s57['teacher_v2']['compact_conditional_on_bls']:.2%} of the {s57['teacher_v2']['n_bls_top5_recovered']:,} BLS-recovered injections. It met both relative-improvement requirements but fell short of the absolute 80% compact-triage requirement.",
            styles["body"],
        )
    )
    add_figure(
        story,
        ARTIFACT_DIR / "s57_external_recovery" / "bls_top5_recovery_map" / "period_radius_empirical_recovery_publication.png",
        "S57: ADP-small BLS top-five recovery",
        "This is the fresh external BLS recovery surface. It uses the same map conventions as S56 so the two sectors can be compared visually. The teacher is evaluated only on the candidates represented by this search outcome.",
        styles,
    )
    add_figure(
        story,
        ARTIFACT_DIR / "s57_external_recovery" / "teacher_v2_compact_recovery_map" / "period_radius_empirical_recovery_publication.png",
        "S57: ADP BLS plus frozen Teacher v2",
        "This is the principal external v2 result. The model, calibration, and p_compact threshold were frozen from S56 development before this sector was opened. Its 77.21% conditional retention is therefore a legitimate transfer measurement, though S57 is now diagnostic-only for future model selection because we have inspected it.",
        styles,
    )
    add_figure(
        story,
        ARTIFACT_DIR / "s57_external_recovery" / "teacher_v1_workload_matched_recovery_map" / "period_radius_empirical_recovery_publication.png",
        "S57: Teacher v1 at the matched workload",
        "This reference map shows the prior teacher on precisely the same external BLS candidates. Teacher v1 retains 56.25% conditional on BLS, establishing the large practical gain from v2 for compact-transit ranking.",
        styles,
    )
    add_figure(
        story,
        ARTIFACT_DIR / "s57_external_recovery" / "teacher_v2_minus_v1_difference_map" / "teacher_v2_minus_v1_recovery.png",
        "S57: Teacher v2 minus Teacher v1",
        "The difference map isolates where v2 improves or degrades external recovery. It should be read together with the support masks: a visually strong bright-cell pattern can still be based on only one host.",
        styles,
    )
    add_figure(
        story,
        figures["external"],
        "BLS-conditional compact-transit retention: v2 clears the relative tests",
        "Teacher v2 improved S57 retention by 20.96 percentage points relative to v1 and by 15.27 points relative to metadata-only. Those are real external improvements. The separate 80% absolute gate is stricter and was missed by 2.79 points, so the result supports enrichment but not pseudo-labeling.",
        styles,
        height=4.8 * inch,
    )

    story.append(Paragraph("Interpretation and next decision", styles["h1"]))
    story.append(
        Paragraph(
            "The central result is not that v2 failed. It is that one shared model was asked to do two different jobs, and the evidence points in opposite directions. The injected compact-transit task improved substantially because v2 learned from many BLS-recovered injected positives and raw chronological information. The rarer real morphology classes, especially Eclipse/contact, did not have enough stable support to preserve v1-level performance in the selected compact-optimized architecture.",
            styles["body"],
        )
    )
    next_table = Table(
        [
            ["Use now", "Do not use yet", "Next model goal"],
            [
                Paragraph("<b>Active-learning enrichment.</b><br/>Score BLS candidates, rank by compact probability and uncertainty, then have humans examine a bounded queue. This is exactly what 77% conditional retention is useful for.", styles["small"]),
                Paragraph("<b>Automatic science labels or student pseudo-labels.</b><br/>The real morphology holdout has too few planets and v2 misses too many eclipses for a model output to become a training truth label.", styles["small"]),
                Paragraph("<b>Separate specialists.</b><br/>Choose a compact-transit ensemble for enrichment and a morphology ensemble for EB/variable classification on S56 development. Test the next selection on a newly held-out sector, not S57.", styles["small"]),
            ],
        ],
        colWidths=[2.32 * inch, 2.32 * inch, 2.32 * inch],
    )
    next_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), HexColor(NAVY)),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Times-Bold"),
                ("BACKGROUND", (0, 1), (0, 1), HexColor(LIGHT_TEAL)),
                ("BACKGROUND", (1, 1), (1, 1), HexColor(LIGHT_RED)),
                ("BACKGROUND", (2, 1), (2, 1), HexColor(LIGHT_BLUE)),
                ("BOX", (0, 0), (-1, -1), 0.5, HexColor("#9BA5AF")),
                ("INNERGRID", (0, 0), (-1, -1), 0.35, HexColor("#AAB2BD")),
                ("FONTSIZE", (0, 0), (-1, -1), 8.5),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 7),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
            ]
        )
    )
    story.append(next_table)
    story.append(Spacer(1, 0.12 * inch))
    story.append(Paragraph("Why S57 cannot be reused for the next formal test", styles["h2"]))
    story.append(
        Paragraph(
            "We have now inspected the S57 result and know where v2 succeeds and fails. Any v2.1 architecture or threshold chosen in response to that information would be indirectly tuned to S57, even if no weights were fit there. S57 remains valuable as a diagnostic benchmark; a later sector must become the new untouched external test.",
            styles["body"],
        )
    )
    story.append(PageBreak())

    story.append(Paragraph("Plain-language glossary", styles["h1"]))
    glossary = [
        ("Teacher", "A model trained from labels or controlled injections that ranks candidates for a later review step. It is not a discovery claim or a final astrophysical classifier."),
        ("BLS", "Box Least Squares, a search that tries many periods and box-shaped transit durations. It proposes candidate ephemerides from a light curve."),
        ("Ephemeris", "The predicted timing of repeated events, usually described by a period and an epoch. Folding at the wrong harmonic can make a real signal look like a different morphology."),
        ("Top five BLS candidates", "The five BLS peaks retained for each target. A truth-matched injection counts as BLS-recovered if one of these peaks has the correct period or allowed harmonic and overlap."),
        ("Injection", "A simulated transit added to a real observed light curve. It supplies known ground truth while retaining realistic observing gaps and noise."),
        ("Fresh or external", "Data never used to select the model architecture, thresholds, or calibration. S57 was external for v2, but it cannot be external again for a model designed after seeing its result."),
        ("Frozen model / threshold", "Weights, score calibration, and the p_compact cutoff were fixed before evaluation. Freezing prevents us from quietly choosing a threshold that happens to look best on the test."),
        ("p_compact", "The ensemble's probability-like score that a BLS candidate resembles a compact transit. It is a ranking score, not a probability that the object is a planet."),
        ("Review workload", "The fraction of real TICs allowed through to human inspection. A 5% workload means the score threshold was set so only about 5% of the held-out real target list passes."),
        ("Conditional retention", "Among signals BLS already found, the fraction the teacher kept. It measures the teacher separately from the upstream search."),
        ("End-to-end recovery", "The fraction of all injected signals recovered by both BLS and the teacher. It is always no larger than BLS recovery in this pipeline."),
        ("Recall", "Of all true examples of one class, the fraction the model identifies. Low EB recall means real eclipsing binaries are often being assigned another label."),
        ("Precision", "Of examples assigned a class by the model, the fraction that really belong to that class according to the reference label."),
        ("Accuracy", "The fraction of all labels predicted correctly. It can be misleading when one class, such as Other, is much more common than the rest."),
        ("Macro F1", "A balanced class-performance score: calculate F1 for every class and average them equally. It stops Other from overwhelming the rare classes in the summary."),
        ("Confusion matrix", "A table of human labels versus model labels. Diagonal cells are correct; off-diagonal cells show the specific mistakes."),
        ("Calibration / ECE", "Whether score values match observed frequencies. ECE, expected calibration error, summarizes the mismatch; lower is better."),
        ("Five-fold TIC-grouped ensemble", "Five models are trained with different groups of stars held out. Grouping by TIC prevents one star from appearing in both train and validation. Their averaged scores are more stable than one model's score."),
        ("Pseudo-label", "A model prediction reused as if it were a human label for student training. This is intentionally blocked until the teacher is robust on real data, because wrong pseudo-labels compound."),
    ]
    glossary_rows = [[Paragraph(f"<b>{term}</b>", styles["small"]), Paragraph(definition, styles["small"])] for term, definition in glossary]
    glossary_table = Table(glossary_rows, colWidths=[1.45 * inch, 5.51 * inch], repeatRows=0)
    glossary_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, -1), HexColor("#EFF3F7")),
                ("BOX", (0, 0), (-1, -1), 0.4, HexColor("#AAB2BD")),
                ("INNERGRID", (0, 0), (-1, -1), 0.25, HexColor("#D0D5DA")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    story.append(glossary_table)
    story.append(PageBreak())

    story.append(Paragraph("Provenance and acceptance audit", styles["h1"]))
    story.append(
        Paragraph(
            "The report was generated from the audited ORCD output set for code checkpoint <b>8f1dad40</b>. Artifact integrity and figure QA passed for both S56 and S57. Tests reported by the completed run: 242 passed, 3 skipped locally; 41 passed, 1 skipped on ORCD. The repository has no requested Make targets, so those checks were not applicable.",
            styles["body"],
        )
    )
    provenance_rows = [
        ["Training rows", f"{roles['n_training_rows']:,}"],
        ["Training unique TICs", f"{roles['n_training_rows']:,} rows across {training_summary['training_verification']['n_unique_tics']:,} TICs"],
        ["S56 locked injections", f"{s56['n_injections']:,} unique hosts"],
        ["S57 external injections", f"{s57['n_injections']:,} unique hosts"],
        ["S56 compact result", f"{s56['teacher_v2']['compact_conditional_on_bls']:.2%} conditional on BLS ({s56['teacher_v2']['n_compact_recovered']:,}/{s56['teacher_v2']['n_bls_top5_recovered']:,})"],
        ["S57 compact result", f"{s57['teacher_v2']['compact_conditional_on_bls']:.2%} conditional on BLS ({s57['teacher_v2']['n_compact_recovered']:,}/{s57['teacher_v2']['n_bls_top5_recovered']:,})"],
        ["S57 acceptance", "Relative improvements passed; absolute 80% compact-retention gate failed; human-morphology gate failed."],
    ]
    provenance_table = Table(provenance_rows, colWidths=[2.05 * inch, 4.91 * inch])
    provenance_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, -1), HexColor("#EFF3F7")),
                ("BOX", (0, 0), (-1, -1), 0.4, HexColor("#AAB2BD")),
                ("INNERGRID", (0, 0), (-1, -1), 0.3, HexColor("#D0D5DA")),
                ("FONTNAME", (0, 0), (0, -1), "Times-Bold"),
                ("FONTNAME", (1, 0), (1, -1), "Times-Roman"),
                ("FONTSIZE", (0, 0), (-1, -1), 8.5),
                ("LEADING", (0, 0), (-1, -1), 10.4),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 7),
                ("RIGHTPADDING", (0, 0), (-1, -1), 7),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.append(provenance_table)
    story.append(Spacer(1, 0.16 * inch))
    story.append(Paragraph("Figure inventory", styles["h2"]))
    inventory = [
        "Custom report figures: pipeline, development profile comparison, selected-profile training curves, compact-score calibration, human holdout confusion/recall, S56/S57 model comparison, and Tmag stratification.",
        "S56 recovery maps: BLS top-five, Teacher v2, workload-matched Teacher v1, and v2-minus-v1 difference.",
        "S57 recovery maps: BLS top-five, Teacher v2, workload-matched Teacher v1, and v2-minus-v1 difference.",
        "S56 workload-performance curve: Teacher-v2 retention versus review workload.",
    ]
    for item in inventory:
        story.append(Paragraph(f"- {item}", styles["body"]))
    story.append(Spacer(1, 0.12 * inch))
    story.append(
        Paragraph(
            "Interpret the recovery maps as classifier-validation and active-learning evidence only. They are not an occurrence-rate completeness measurement and do not approve student training or automatic candidate disposition.",
            styles["body"],
        )
    )
    doc.build(story, onFirstPage=first_page_footer, onLaterPages=footer)


def main() -> None:
    if not ARTIFACT_DIR.exists():
        raise FileNotFoundError(f"Missing audited artifact directory: {ARTIFACT_DIR}")
    selection = read_json(ARTIFACT_DIR / "frozen" / "frozen_selection.json")
    human = read_json(ARTIFACT_DIR / "human_holdout" / "same_row_human_holdout_metrics.json")
    s56 = read_json(ARTIFACT_DIR / "s56_locked_holdout_recovery" / "summary.json")
    s57 = read_json(ARTIFACT_DIR / "s57_external_recovery" / "summary.json")
    roles = read_json(ARTIFACT_DIR / "roles" / "role_summary.json")
    training_summary = read_json(
        ARTIFACT_DIR / "training" / "training_summary__shape_plus_raw_chronology.json"
    )
    figures = {
        "pipeline": make_pipeline_figure(),
        "profiles": make_profile_comparison(selection),
        "training": make_training_curves(),
        "calibration": make_calibration_figure(selection),
        "human": make_human_holdout_figure(human),
        "external": make_external_comparison(s56, s57),
        "tmag": make_tmag_figure(),
    }
    build_pdf(figures, selection, human, s56, s57, roles, training_summary)
    manifest = {
        "report_pdf": str(PDF_PATH),
        "selected_profile": selection["selected_profile"],
        "frozen_threshold": selection["frozen_compact_threshold"],
        "s56_conditional_retention": s56["teacher_v2"]["compact_conditional_on_bls"],
        "s57_conditional_retention": s57["teacher_v2"]["compact_conditional_on_bls"],
        "figures": {name: str(path) for name, path in figures.items()},
    }
    with (REPORT_DIR / "report_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(PDF_PATH)


if __name__ == "__main__":
    main()
