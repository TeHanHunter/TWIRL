#!/usr/bin/env python3
"""Build the 2026-07-16 TWIRL repository audit figures and report.

The figure and PDF dependencies live in different project runtimes, so the
builder exposes independent modes::

    PYTHONPATH=src .venv/bin/python scripts/build_repository_audit_report.py figures
    .venv/bin/python scripts/build_repository_audit_report.py markdown
    /Users/tehan/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
        scripts/build_repository_audit_report.py pdf

All outputs are deterministic and use only repository-audit facts.  No
third-party imagery is included.
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path
from textwrap import dedent


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "reports" / "repository_audit_2026-07-16"
CANONICAL_RECOVERY_MAP_DIR = (
    REPO_ROOT
    / "reports"
    / "stage5_validation"
    / "s56_A2v1_teacher_v1_recovery"
    / "full"
    / "final"
    / "bls_top5_recovery_map"
)
CANONICAL_RECOVERY_MAP_STEM = "period_radius_empirical_recovery_publication"
AUDIT_RECOVERY_MAP_STEM = "a2v1_teacher_v1_bls_top5_recovery_surface"

INJECTED = 20_000
BLS_TOP5 = 3_458
TEACHER_RETAINED = 1_527

STAGE_ROWS = (
    {
        "stage": "Stage 1",
        "scope": "Light curves",
        "level": 2,
        "status": "Pilot",
        "gate": "Release manifest and Tier-1\nphotometric QA",
    },
    {
        "stage": "Stage 2",
        "scope": "Search",
        "level": 2,
        "status": "Pilot",
        "gate": "Dip branch, multi-sector\naggregation, and false alarms",
    },
    {
        "stage": "Stage 3",
        "scope": "Injections",
        "level": 1,
        "status": "Prototype",
        "gate": "Frozen A2v1 recovery and\npixel-level calibration subset",
    },
    {
        "stage": "Stage 4",
        "scope": "Survey search",
        "level": 0,
        "status": "Not built",
        "gate": "Survey-wide inference driver\nand release product contract",
    },
    {
        "stage": "Stage 5",
        "scope": "Validation",
        "level": 1,
        "status": "Prototype",
        "gate": "On-target validation, candidate\nmerging, and reproducibility",
    },
)


def _ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sync_canonical_recovery_map(output_dir: Path) -> None:
    """Copy the canonical Teacher-v1 recovery surface into the audit package."""

    output_dir = _ensure_output_dir(output_dir)
    for suffix in ("png", "pdf"):
        source = CANONICAL_RECOVERY_MAP_DIR / (
            f"{CANONICAL_RECOVERY_MAP_STEM}.{suffix}"
        )
        destination = output_dir / f"{AUDIT_RECOVERY_MAP_STEM}.{suffix}"
        if not source.is_file():
            raise FileNotFoundError(f"Missing canonical recovery map: {source}")
        shutil.copy2(source, destination)


def build_figures(output_dir: Path) -> None:
    """Create the two publication-facing audit figures as PNG and PDF."""

    matplotlib_cache = REPO_ROOT / "tmp" / "matplotlib"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))

    import matplotlib.pyplot as plt
    import numpy as np

    from twirl.plotting.style import apply_twirl_style, get_ordered_palette

    output_dir = _ensure_output_dir(output_dir)
    template = apply_twirl_style("full_page")
    palette = get_ordered_palette(4, "viridis")
    neutral = "#b7bcc2"
    dark = "#34383d"

    # Figure 1: engineering maturity and the next explicit gate for each stage.
    fig, (ax, ax_gate) = plt.subplots(
        1,
        2,
        figsize=template["figsize"],
        gridspec_kw={"width_ratios": [1.05, 1.25], "wspace": 0.12},
    )
    y = np.arange(len(STAGE_ROWS))[::-1]
    level_colors = [neutral, palette[1], palette[2], palette[3]]
    for yi, row in zip(y, STAGE_ROWS, strict=True):
        level = int(row["level"])
        ax.plot([0, 3], [yi, yi], color="#e2e4e7", lw=3.0, zorder=1)
        if level > 0:
            ax.plot([0, level], [yi, yi], color=level_colors[level], lw=3.0, zorder=2)
        ax.scatter(
            level,
            yi,
            s=64,
            color=level_colors[level],
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        ax.text(
            level + 0.08,
            yi + 0.20,
            str(row["status"]),
            color=dark,
            fontsize=template["annotation_size"],
            va="center",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{row['stage']}\n{row['scope']}" for row in STAGE_ROWS],
        ha="right",
    )
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(["Not built", "Prototype", "Pilot", "Survey-ready"])
    ax.set_xlim(-0.15, 3.2)
    ax.set_ylim(-0.65, len(STAGE_ROWS) - 0.35)
    ax.set_xlabel("Audit-assessed engineering maturity")
    ax.grid(axis="x", color="#e6e8ea", lw=0.6)
    ax.grid(axis="y", visible=False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(
        0.0,
        1.035,
        "(a) Current readiness",
        transform=ax.transAxes,
        fontsize=template["title_size"],
        va="bottom",
    )

    ax_gate.set_xlim(0, 1)
    ax_gate.set_ylim(ax.get_ylim())
    ax_gate.axis("off")
    ax_gate.axvline(0.02, color="#d5d8dc", lw=0.8)
    for yi, row in zip(y, STAGE_ROWS, strict=True):
        ax_gate.text(
            0.07,
            yi,
            str(row["gate"]),
            va="center",
            ha="left",
            fontsize=template["label_size"],
            linespacing=1.25,
        )
    ax_gate.text(
        0.07,
        1.035,
        "(b) Next evidence gate",
        transform=ax_gate.transAxes,
        fontsize=template["title_size"],
        va="bottom",
    )
    fig.text(
        0.01,
        0.005,
        "Readiness is an engineering audit judgment, not a scientific-performance score.",
        fontsize=template["annotation_size"],
        color="#55595e",
    )
    fig.subplots_adjust(left=0.16, right=0.99, bottom=0.18, top=0.91)
    for suffix in ("png", "pdf"):
        fig.savefig(
            output_dir / f"stage_readiness.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)

    # Figure 2: the exact candidate-retention attrition through BLS and teacher.
    apply_twirl_style("full_page")
    labels = ["Injected signals", "BLS top-five", "Teacher retained"]
    counts = np.asarray([INJECTED, BLS_TOP5, TEACHER_RETAINED], dtype=float)
    fractions = 100.0 * counts / INJECTED
    colors = [neutral, palette[2], palette[3]]

    fig, ax = plt.subplots(figsize=template["figsize"])
    y = np.arange(len(labels))[::-1]
    bars = ax.barh(y, fractions, color=colors, edgecolor="white", height=0.58)
    for bar, count, fraction in zip(bars, counts.astype(int), fractions, strict=True):
        x = bar.get_width()
        inside = x >= 26
        ax.text(
            x - 1.5 if inside else x + 1.4,
            bar.get_y() + bar.get_height() / 2,
            f"{count:,}  ({fraction:.2f}%)",
            ha="right" if inside else "left",
            va="center",
            color="white" if inside else dark,
            fontweight="bold",
            fontsize=template["label_size"],
        )
    conditional = 100.0 * TEACHER_RETAINED / BLS_TOP5
    ax.text(
        37,
        y[2],
        f"Conditional retention\nteacher / BLS top-five = {conditional:.2f}%",
        ha="left",
        va="center",
        fontsize=template["annotation_size"],
        color=dark,
        linespacing=1.35,
        bbox={
            "boxstyle": "round,pad=0.45",
            "facecolor": "#f5f6f7",
            "edgecolor": "#cfd4d8",
            "linewidth": 0.7,
        },
    )
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Fraction of injected signals retained (%)")
    ax.set_xlim(0, 106)
    ax.set_ylim(-0.65, len(labels) - 0.35)
    ax.grid(axis="x", color="#e6e8ea", lw=0.6)
    ax.grid(axis="y", visible=False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(
        0.0,
        -0.25,
        "Candidate-retention diagnostic only: this is not end-to-end survey completeness.",
        transform=ax.transAxes,
        fontsize=template["annotation_size"],
        color="#55595e",
        va="top",
    )
    fig.subplots_adjust(left=0.20, right=0.98, bottom=0.25, top=0.97)
    for suffix in ("png", "pdf"):
        fig.savefig(
            output_dir / f"candidate_retention_attrition.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)
    sync_canonical_recovery_map(output_dir)


def markdown_text() -> str:
    """Return the concise paper-style audit report in Markdown."""

    bls_fraction = 100.0 * BLS_TOP5 / INJECTED
    retained_fraction = 100.0 * TEACHER_RETAINED / INJECTED
    conditional = 100.0 * TEACHER_RETAINED / BLS_TOP5
    return dedent(
        f"""\
        # TWIRL repository audit and near-term execution plan

        **Internal technical report - 2026-07-16**

        ## Abstract

        We audited the TWIRL repository as an executable survey pipeline, with emphasis on reproducibility, production safety, validation semantics, generated-artifact hygiene, and alignment between the Stage 1-5 plan and the code that now exists. The repository has a credible pilot foundation: Stage 1 produces validated A2v1 light curves, Stage 2 supports an interpretable periodic search, and the fast validation suite is healthy. It is not yet a survey-complete system. The critical gaps are a frozen release manifest, science-level photometric QA, the non-periodic dip branch, multi-sector aggregation, false-alarm calibration, frozen-product injection recovery, and on-target end-to-end validation. We recommend keeping teacher development as bounded triage work while moving these baseline survey gates onto the critical path.

        ## Methods

        The scan combined repository-wide inventory, Git/worktree inspection, parser and link checks, shell and Python syntax checks, the project fast test suite, and review of the authoritative plan, production protocol, data conventions, and plotting rules. The audit checkpoint contained 1,964 tracked files (955 under `reports/`, 611 under `outputs/`, 271 under `scripts/`, 49 under `src/`, and 48 under `tests/`) and about 740 MiB of tracked working-tree content. The fast suite completed with 246 passed and 2 skipped tests. Code paths were compared with the declared Stage 1-5 contract. We distinguished a Tier-0 integrity gate (schema, coverage, openability, and benchmark checks) from the still-open Tier-1 science QA gate (photometric precision, cadence loss, aperture outliers, injection preservation, and an independent extraction comparison). Generated data and historical binaries were treated as reproducible artifacts unless they were compact evidence needed to interpret a result.

        ## Results

        ![Stage readiness and next gates](stage_readiness.png)

        **Figure 1.** Audit-assessed engineering maturity and the next evidence gate for each pipeline stage. The readiness scale describes implementation maturity, not scientific performance. Stage 1 and the periodic part of Stage 2 are at pilot maturity; Stage 3 and Stage 5 remain prototypes; the Stage 4 survey-wide inference layer is not yet built.

        The strongest result is separation of the production contract from later analysis: A2v1 generation, validation, compact export, and an interpretable BLS path now have explicit boundaries. The most important weakness is that an integrity-valid sector can still lack the evidence required for scientific readiness. The audit also found a mismatch between the clean forward plan and experimental teacher-v2 and S57 labeling work already present in the repository. These products should be preserved as exploratory evidence, but they should not silently redefine the production baseline or consume the untouched validation holdout.

        Repository hygiene is acceptable at the source level but weak at the artifact-history boundary. At scan time the Git object store was approximately 14 GiB and dominated by historical generated outputs. The on-disk `reports/` tree was approximately 5.5 GiB because it also contained large ignored or untracked review artifacts, while `src/`, `scripts/`, and `tests/` together occupied only about 15 MiB. Nonportable symlinks, stale caches, large third-party reference bundles, and redundant rendered review material should remain outside versioned survey products. History rewriting is not recommended during the active branch stack; compact manifests and an explicit artifact allowlist are safer immediate controls. The local `.venv` reported NumPy 2.2.6 even though `pyproject.toml` requires `numpy>=1.24,<2.0`; the tests passed, so this is environment drift rather than an observed failure.

        ![Candidate-retention attrition](candidate_retention_attrition.png)

        **Figure 2.** Candidate-retention diagnostic for the 20,000-signal teacher-v1 experiment. The top-five BLS stage retained {BLS_TOP5:,} signals ({bls_fraction:.2f}% of injections), and the teacher retained {TEACHER_RETAINED:,} ({retained_fraction:.3f}% of injections; {conditional:.2f}% conditional on top-five BLS recovery). This diagnostic does not include the full search, vetter, candidate-merging, or pixel-level calibration chain and therefore is not end-to-end survey completeness.

        The attrition is scientifically useful because it localizes losses, but it also shows why classifier iteration cannot substitute for baseline search and calibration. The current language should reserve “end-to-end recovery” for a frozen chain that includes search, optional ranking, vetting, candidate merging, and the adopted calibration products. The present teacher result is more accurately a candidate-retention efficiency.

        ### Period-radius candidate-retention surface

        ![A2v1 Teacher-v1 BLS top-five candidate-retention surface](a2v1_teacher_v1_bls_top5_recovery_surface.png)

        **Figure 3.** A2v1 Teacher-v1 BLS top-five candidate-retention surface for 20,000 injected signals, split into four Tmag strata. Recovery fractions are kernel-smoothed; each panel prints its recovered and injected support counts, and the bright bin is comparatively sparse. The map diagnoses where the BLS candidate list retains injected signals. It is not end-to-end survey completeness and should not be used as an occurrence-rate selection function.

        The surface supplies more structure than the aggregate attrition count: candidate retention depends on injected period, companion radius, and target brightness. Its role here is to expose support and search behavior, not to promote a final completeness claim. A frozen survey analysis must propagate the same injections through the adopted vetter, candidate merger, and calibration chain before interpreting the surface as completeness.

        ## Recommendations

        1. **Freeze the survey release contract.** Record an explicit TWIRL-I sector cutoff, target/sample version, A2v1 product tag, search configuration, and checksums in a machine-readable release manifest. Without a cutoff, the `Sector >= 56` denominator expands continuously.
        2. **Complete Tier-1 photometric QA.** Add scatter-versus-magnitude regression, cadence-loss distributions, aperture disagreement thresholds, fixed injection-preservation tests, and a genuinely independent benchmark extraction. Keep Tier-0 integrity status separate.
        3. **Finish the transparent baseline search before new ML gates.** Implement the non-periodic dip branch, multi-sector candidate consolidation, and an empirical false-alarm strategy. These are required even if teacher ranking remains useful.
        4. **Calibrate recovery on frozen products.** Run injections through the adopted search, ranker if used, vetter, and merger; include a bounded pixel-level subset to measure detrending and extraction losses.
        5. **Quarantine exploratory teacher work.** Preserve teacher-v2 and early S57 labels with provenance, mark them exploratory, pause further holdout consumption, and require external-retention and morphology checks before promotion.
        6. **Control repository growth prospectively.** Version compact summaries, manifests, selected figures, and small label tables. Keep review sheets, dependency bundles, shard payloads, and local literature copies in ignored storage. Defer Git history surgery until every active branch and worktree is clean and pushed.
        7. **Rebuild and automate the release environment.** Rebuild the local environment from `pyproject.toml` before release validation. Once the dependency or lock strategy is settled, add pull-request CI for `make test-fast` and `make check-docs`. This is a medium-priority engineering safeguard below the science critical path, not evidence of a current test failure.

        A practical execution order is: (i) finish the A2v1 queue with Tier-0 gates and compact exports, (ii) freeze the release manifest and implement Tier-1 QA, (iii) complete the dip, multi-sector, and false-alarm baseline, (iv) run frozen-product recovery and pixel calibration, and (v) perform survey-wide inference and on-target validation. Bounded human labeling can proceed in parallel only when it does not delay these gates.

        ## Conclusion

        TWIRL is beyond a collection of experiments: it has a testable pilot pipeline and a clear production contract. The next advance should be evidentiary rather than architectural. Freezing the release, separating integrity from science QA, completing the transparent search branches, and measuring recovery through the whole adopted chain will convert the current pilot into a defensible survey framework. The audit figures and this report are reproducible from [`scripts/build_repository_audit_report.py`](../../scripts/build_repository_audit_report.py); the project plan and detailed execution record remain authoritative in [`doc/twirl_plan.md`](../../doc/twirl_plan.md) and [`doc/twirl_progress_log.md`](../../doc/twirl_progress_log.md).
        """
    )


def build_markdown(output_dir: Path) -> None:
    output_dir = _ensure_output_dir(output_dir)
    (output_dir / "repository_audit_report.md").write_text(
        markdown_text(), encoding="utf-8"
    )


def build_pdf(output_dir: Path) -> None:
    """Build a four-page PDF with ReportLab and the generated PNG figures."""

    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
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

    output_dir = _ensure_output_dir(output_dir)
    readiness_path = output_dir / "stage_readiness.png"
    attrition_path = output_dir / "candidate_retention_attrition.png"
    recovery_surface_path = output_dir / f"{AUDIT_RECOVERY_MAP_STEM}.png"
    for path in (readiness_path, attrition_path, recovery_surface_path):
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path.name}; run the 'figures' mode before 'pdf'."
            )

    navy = colors.HexColor("#263746")
    blue = colors.HexColor("#356A8A")
    light_blue = colors.HexColor("#EAF1F5")
    pale = colors.HexColor("#F5F6F7")
    gray = colors.HexColor("#5B6168")
    line = colors.HexColor("#CFD4D8")

    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="AuditTitle",
            parent=styles["Title"],
            fontName="Helvetica-Bold",
            fontSize=20,
            leading=24,
            textColor=navy,
            alignment=TA_LEFT,
            spaceAfter=8,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AuditKicker",
            parent=styles["Normal"],
            fontName="Helvetica-Bold",
            fontSize=7.5,
            leading=9,
            textColor=blue,
            spaceAfter=5,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AuditSection",
            parent=styles["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=12,
            leading=14,
            textColor=navy,
            spaceBefore=7,
            spaceAfter=5,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AuditBody",
            parent=styles["BodyText"],
            fontName="Helvetica",
            fontSize=8.7,
            leading=12.1,
            textColor=colors.HexColor("#25292D"),
            spaceAfter=6,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AuditAbstract",
            parent=styles["AuditBody"],
            fontName="Helvetica",
            fontSize=9.1,
            leading=13,
            leftIndent=11,
            rightIndent=11,
            borderColor=line,
            borderWidth=0.6,
            borderPadding=9,
            backColor=pale,
            spaceBefore=5,
            spaceAfter=10,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AuditCaption",
            parent=styles["BodyText"],
            fontName="Helvetica",
            fontSize=7.4,
            leading=10,
            textColor=gray,
            spaceBefore=3,
            spaceAfter=7,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AuditSmall",
            parent=styles["BodyText"],
            fontName="Helvetica",
            fontSize=7.2,
            leading=9.2,
            textColor=gray,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AuditRecHead",
            parent=styles["BodyText"],
            fontName="Helvetica-Bold",
            fontSize=8,
            leading=10,
            textColor=navy,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AuditRecBody",
            parent=styles["BodyText"],
            fontName="Helvetica",
            fontSize=7.5,
            leading=9.7,
            textColor=colors.HexColor("#25292D"),
        )
    )

    pdf_path = output_dir / "repository_audit_report.pdf"
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=letter,
        rightMargin=0.62 * inch,
        leftMargin=0.62 * inch,
        topMargin=0.58 * inch,
        bottomMargin=0.58 * inch,
        title="TWIRL repository audit and near-term execution plan",
        author="TWIRL project repository audit",
        subject="Repository audit, pipeline readiness, and prioritized next steps",
        invariant=1,
    )

    def para(text: str, style: str = "AuditBody") -> Paragraph:
        return Paragraph(text, styles[style])

    def page_frame(canvas, document) -> None:
        canvas.saveState()
        width, height = letter
        canvas.setStrokeColor(line)
        canvas.setLineWidth(0.5)
        canvas.line(0.62 * inch, height - 0.41 * inch, width - 0.62 * inch, height - 0.41 * inch)
        canvas.setFillColor(gray)
        canvas.setFont("Helvetica", 6.8)
        canvas.drawString(0.62 * inch, height - 0.30 * inch, "TWIRL REPOSITORY AUDIT")
        canvas.drawRightString(width - 0.62 * inch, 0.32 * inch, f"{document.page}")
        canvas.restoreState()

    story = []

    # Page 1: scope and headline results.
    story.extend(
        [
            Spacer(1, 0.14 * inch),
            para("INTERNAL TECHNICAL REPORT | 2026-07-16", "AuditKicker"),
            para("TWIRL repository audit and<br/>near-term execution plan", "AuditTitle"),
            para(
                "A reproducibility, production-safety, and scientific-readiness review of the standalone survey pipeline.",
                "AuditSmall",
            ),
            Spacer(1, 0.20 * inch),
            para("Abstract", "AuditSection"),
            Spacer(1, 0.07 * inch),
            para(
                "We audited TWIRL as an executable survey pipeline, emphasizing reproducibility, production safety, validation semantics, generated-artifact hygiene, and alignment between the Stage 1-5 plan and the code that now exists. The repository has a credible pilot foundation: Stage 1 produces validated A2v1 light curves, Stage 2 supports an interpretable periodic search, and the fast validation suite is healthy. It is not yet a survey-complete system. The critical gaps are a frozen release manifest, science-level photometric QA, the non-periodic dip branch, multi-sector aggregation, false-alarm calibration, frozen-product injection recovery, and on-target end-to-end validation. Teacher development should remain bounded triage work while these baseline survey gates move onto the critical path.",
                "AuditAbstract",
            ),
            para("Audit scope", "AuditSection"),
        ]
    )
    scope_data = [
        [para("Dimension", "AuditRecHead"), para("Evidence reviewed", "AuditRecHead")],
        [para("Pipeline", "AuditRecBody"), para("Stage 1-5 implementation, interfaces, and acceptance gates", "AuditRecBody")],
        [para("Validation", "AuditRecBody"), para("Fast tests, docs, detection sample, parsers, links, and syntax", "AuditRecBody")],
        [para("Operations", "AuditRecBody"), para("Git branches/worktrees, artifact boundaries, local data, and production safety", "AuditRecBody")],
        [para("Science", "AuditRecBody"), para("Claim boundaries, benchmark evidence, recovery semantics, and holdout use", "AuditRecBody")],
    ]
    scope_table = Table(scope_data, colWidths=[1.25 * inch, 5.55 * inch], repeatRows=1)
    scope_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), light_blue),
                ("BOX", (0, 0), (-1, -1), 0.5, line),
                ("INNERGRID", (0, 0), (-1, -1), 0.35, line),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.extend(
        [
            scope_table,
            Spacer(1, 0.12 * inch),
            para("Headline findings", "AuditSection"),
            para(
                "<b>1.</b> A2v1 generation, validation, compact export, and periodic BLS now form a credible pilot path. <b>2.</b> Integrity-valid is not yet science-ready: Tier-1 photometric QA remains open. <b>3.</b> Teacher-v2 and early S57 labels are exploratory evidence and should not silently redefine the baseline or holdout. <b>4.</b> Git history is artifact-heavy; prospective curation is safer than history surgery during the active branch stack.",
            ),
            para(
                "Bottom line: the next advance should be evidentiary rather than architectural.",
                "AuditAbstract",
            ),
            PageBreak(),
        ]
    )

    # Page 2: method and maturity/gate figure.
    story.extend(
        [
            para("Methods", "AuditSection"),
            para(
                "The scan combined repository-wide inventory, Git/worktree inspection, parser and link checks, shell and Python syntax checks, the project fast test suite, and review of the authoritative plan, production protocol, data conventions, and plotting rules. The checkpoint contained 1,964 tracked files and about 740 MiB of tracked content; the fast suite completed with 246 passed and 2 skipped tests. We separated a Tier-0 integrity gate (schema, coverage, openability, and benchmark checks) from the still-open Tier-1 science QA gate (photometric precision, cadence loss, aperture outliers, injection preservation, and an independent extraction comparison).",
            ),
            para("Pipeline readiness", "AuditSection"),
            Image(str(readiness_path), width=6.92 * inch, height=3.99 * inch),
            para(
                "<b>Figure 1.</b> Audit-assessed engineering maturity and the next evidence gate for each pipeline stage. The scale describes implementation maturity, not scientific performance. Stage 1 and the periodic part of Stage 2 are at pilot maturity; Stage 3 and Stage 5 remain prototypes; Stage 4 is not yet built.",
                "AuditCaption",
            ),
            para(
                "The strongest result is separation of the production contract from later analysis. The main weakness is that an integrity-valid sector can still lack the evidence required for scientific readiness. The clean forward plan also diverged from experimental teacher-v2 and S57 labeling work already present. Those products should be preserved with provenance but quarantined from the production baseline and untouched validation claims.",
            ),
            PageBreak(),
        ]
    )

    # Page 3: artifact boundary and injection attrition.
    story.extend(
        [
            para("Results", "AuditSection"),
            para(
                "Repository hygiene is acceptable at the source level but weak at the artifact-history boundary. At scan time the Git object store was approximately 14 GiB and dominated by historical generated outputs. The on-disk reports tree was approximately 5.5 GiB because it also contained large ignored or untracked review artifacts, while source, scripts, and tests together occupied only about 15 MiB. Nonportable symlinks, stale caches, large third-party reference bundles, and redundant rendered material belong outside versioned survey products. The local environment reported NumPy 2.2.6 although pyproject.toml pins NumPy below 2.0; tests passed, so this is environment drift rather than an observed failure.",
            ),
            Image(str(attrition_path), width=6.92 * inch, height=3.99 * inch),
            para(
                f"<b>Figure 2.</b> Candidate-retention diagnostic for the {INJECTED:,}-signal teacher-v1 experiment. Top-five BLS retained {BLS_TOP5:,} signals ({100 * BLS_TOP5 / INJECTED:.2f}% of injections), and the teacher retained {TEACHER_RETAINED:,} ({100 * TEACHER_RETAINED / INJECTED:.3f}% of injections; {100 * TEACHER_RETAINED / BLS_TOP5:.2f}% conditional on top-five BLS recovery). This is not end-to-end survey completeness.",
                "AuditCaption",
            ),
            para(
                "The attrition localizes losses, but it also shows why classifier iteration cannot substitute for baseline search and calibration. Reserve <i>end-to-end recovery</i> for a frozen chain that includes search, optional ranking, vetting, candidate merging, and the adopted calibration products. The present teacher result is a candidate-retention efficiency.",
            ),
            PageBreak(),
        ]
    )

    # Page 4: canonical period-radius candidate-retention evidence.
    story.extend(
        [
            para("Period-radius candidate-retention surface", "AuditSection"),
            para(
                "The aggregate retention counts localize pipeline attrition, while the canonical Teacher-v1 surface shows where the periodic search retains candidates across injected period, companion radius, and target brightness. It is included as diagnostic evidence, not as a final survey selection function.",
            ),
            Image(
                str(recovery_surface_path),
                width=6.92 * inch,
                height=5.61 * inch,
            ),
            para(
                "<b>Figure 3.</b> A2v1 Teacher-v1 BLS top-five candidate-retention surface for 20,000 injected signals, split into four Tmag strata. Fractions are kernel-smoothed; recovered and injected support counts are printed in every panel, and the bright bin is comparatively sparse. This is not end-to-end survey completeness and is not an occurrence-rate selection function.",
                "AuditCaption",
            ),
            para(
                "A frozen survey analysis must propagate the same injections through the adopted vetter, candidate merger, and calibration chain before interpreting this surface as completeness.",
            ),
            PageBreak(),
        ]
    )

    # Page 5: prioritized recommendations and conclusion.
    story.extend([para("Recommendations", "AuditSection")])
    recommendations = [
        ("1 | Freeze the release contract", "Record a TWIRL-I sector cutoff, target/sample version, product tag, search configuration, and checksums in a machine-readable manifest."),
        ("2 | Complete Tier-1 science QA", "Add precision-versus-magnitude regression, cadence-loss distributions, aperture thresholds, fixed injection preservation, and an independent benchmark extraction."),
        ("3 | Finish the transparent baseline", "Implement the dip branch, multi-sector consolidation, and empirical false-alarm strategy before adding new ML gates."),
        ("4 | Calibrate frozen-product recovery", "Run the adopted search, optional ranker, vetter, and merger; include a bounded pixel-level subset to capture extraction and detrending losses."),
        ("5 | Quarantine exploratory teacher work", "Preserve teacher-v2 and early S57 labels with provenance, pause further holdout use, and require external-retention plus morphology checks before promotion."),
        ("6 | Control repository growth", "Version compact summaries, manifests, selected figures, and small label tables; ignore rendered sheets, dependency bundles, shards, and local literature copies."),
        ("7 | Rebuild and automate the environment", "Rebuild from pyproject.toml before release validation. After settling the dependency or lock strategy, add PR CI for make test-fast and make check-docs. This is medium priority below the science gates."),
    ]
    rec_rows = []
    for head, body in recommendations:
        rec_rows.append([para(head, "AuditRecHead"), para(body, "AuditRecBody")])
    rec_table = Table(rec_rows, colWidths=[2.12 * inch, 4.68 * inch])
    rec_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.white),
                ("ROWBACKGROUNDS", (0, 0), (-1, -1), [colors.white, pale]),
                ("BOX", (0, 0), (-1, -1), 0.5, line),
                ("INNERGRID", (0, 0), (-1, -1), 0.35, line),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.extend(
        [
            rec_table,
            Spacer(1, 0.12 * inch),
            para("Execution order", "AuditSection"),
            para(
                "Finish the A2v1 queue with Tier-0 gates and compact exports; freeze the release manifest and implement Tier-1 QA; complete dip, multi-sector, and false-alarm baselines; run frozen-product recovery and pixel calibration; then perform survey-wide inference and on-target validation. Bounded labeling may proceed in parallel only when it does not delay these gates.",
            ),
            para("Conclusion", "AuditSection"),
            para(
                "TWIRL is beyond a collection of experiments: it has a testable pilot pipeline and a clear production contract. Freezing the release, separating integrity from science QA, completing the transparent search branches, and measuring recovery through the whole adopted chain will convert the pilot into a defensible survey framework.",
            ),
            KeepTogether(
                [
                    para("Reproducibility", "AuditSection"),
                    Spacer(1, 0.07 * inch),
                    para(
                        "Figures and this PDF are generated by scripts/build_repository_audit_report.py. The source report is repository_audit_report.md. The project plan and progress log remain the authoritative status records.",
                        "AuditAbstract",
                    ),
                ]
            ),
        ]
    )

    doc.build(story, onFirstPage=page_frame, onLaterPages=page_frame)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        choices=("figures", "markdown", "pdf"),
        help="Artifact family to build; run modes separately for lazy dependencies.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if args.mode == "figures":
        build_figures(output_dir)
    elif args.mode == "markdown":
        build_markdown(output_dir)
    else:
        build_pdf(output_dir)


if __name__ == "__main__":
    main()
