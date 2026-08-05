#!/usr/bin/env python3
"""Package the current human-accepted S56 Planet-like vetting candidates."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.plotting.style import apply_twirl_style  # noqa: E402


DEFAULT_ACCEPTED = (
    REPO_ROOT
    / "reports/stage5_validation/s56_s62_morphology_corpus_teacher_v3_v1"
    / "accepted_signal_rereview.csv"
)
DEFAULT_METRICS = (
    REPO_ROOT
    / "data_local/label_reviews/s56_s62_signal_rereview"
    / "render_provenance_s56_s62_a2v1_current_adp_v1"
    / "s0056_metrics.csv"
)
DEFAULT_ELIGIBILITY = (
    REPO_ROOT
    / "reports/stage1_lightcurves/s56_A2v1_tier1_active_search_pair_v1_20260723"
    / "target_eligibility.csv"
)
DEFAULT_RENDERED_MANIFEST = (
    REPO_ROOT
    / "data_local/label_reviews/s56_s62_signal_rereview"
    / "rendered_sheet_manifest.csv"
)
DEFAULT_SHEET_SET = (
    REPO_ROOT
    / "data_local/label_reviews/s56_s62_signal_rereview"
    / "sheet_set.json"
)
DEFAULT_SOURCE_SHEETS = (
    REPO_ROOT
    / "data_local/label_reviews/s56_s62_signal_rereview"
    / "vet_sheets_s56_s62_a2v1_current_adp_v1"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_planet_like_candidates_current_20260727"
)

BENCHMARK_TIC = 267574918
ORIGIN_LABELS = {
    "s56_adjudicated_real343": "human-adjudicated",
    "s56_compact_planet_revisit_407": "model-enriched, human-confirmed",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _clean_float(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result if np.isfinite(result) else float("nan")


def _font(size: int, *, bold: bool = False) -> ImageFont.ImageFont:
    candidates = (
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
        if bold
        else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    )
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def _priority_group(row: pd.Series) -> str:
    tic = int(row["tic"])
    if tic == BENCHMARK_TIC:
        return "Benchmark"
    if not _as_bool(row["aperture_disagreement_flag_current"]):
        return "A"
    if str(row["tier1_target_qa_status"]).strip().lower() == "pass":
        return "B"
    return "C"


def _priority_note(row: pd.Series) -> str:
    group = str(row["priority_group"])
    if group == "Benchmark":
        return "Confirmed WD 1856 b benchmark; recovered in both current ADP apertures."
    if group == "A":
        return "Current ADP apertures recover aligned periods; prioritize contamination and independent-extraction checks."
    if group == "B":
        return "Planet-like in the small search aperture, but the current primary aperture prefers a different period."
    return (
        "Small-aperture Planet-like morphology with aperture disagreement and a "
        "non-blocking Tier-1 target-QA review warning."
    )


def _build_candidate_table(
    *,
    accepted_path: Path,
    metrics_path: Path,
    eligibility_path: Path,
    rendered_manifest_path: Path,
    source_sheet_dir: Path,
) -> pd.DataFrame:
    accepted = pd.read_csv(accepted_path, dtype={"tic": "int64", "sector": "int64"})
    candidates = accepted[
        accepted["sector"].eq(56) & accepted["human_label"].eq("planet_like")
    ].copy()
    if len(candidates) != 30 or candidates["tic"].nunique() != 30:
        raise RuntimeError(
            "expected the frozen S56 morphology corpus to contain exactly "
            f"30 unique Planet-like candidates, found {len(candidates)} rows "
            f"and {candidates['tic'].nunique()} TICs"
        )

    metrics = pd.read_csv(metrics_path, dtype={"tic": "int64", "sector": "int64"})
    metrics = metrics.add_suffix("_current").rename(
        columns={"review_id_current": "review_id", "tic_current": "tic"}
    )
    candidates = candidates.merge(
        metrics,
        on=["review_id", "tic"],
        how="left",
        validate="one_to_one",
    )
    if candidates["twirl_vet_status_current"].isna().any():
        raise RuntimeError("one or more accepted S56 Planet-like rows lack current render metrics")
    bad_render = ~candidates["twirl_vet_status_current"].eq("ok")
    if bad_render.any():
        raise RuntimeError(
            "current render metrics are not all ok: "
            + ", ".join(candidates.loc[bad_render, "review_id"].astype(str))
        )

    eligibility = pd.read_csv(
        eligibility_path,
        dtype={"tic": "int64", "sector": "int64"},
    )
    eligibility_columns = [
        "tic",
        "tier1_contract_version",
        "tier1_config_name",
        "tier1_target_qa_status",
        "tier1_target_qa_reasons",
        "tier1_target_qa_pass",
        "tier1_target_searchable",
        "tier1_target_searchability_reasons",
        "n_finite_quality0_min",
        "usable_cadence_fraction_min",
    ]
    candidates = candidates.merge(
        eligibility.loc[:, eligibility_columns],
        on="tic",
        how="left",
        validate="many_to_one",
    )
    if candidates["tier1_target_searchable"].isna().any():
        raise RuntimeError("one or more candidates lack Tier-1 target eligibility")
    if not candidates["tier1_target_searchable"].map(_as_bool).all():
        failed = candidates.loc[
            ~candidates["tier1_target_searchable"].map(_as_bool), "tic"
        ]
        raise RuntimeError(
            "one or more accepted Planet-like candidates are not searchable: "
            + ", ".join(failed.astype(str))
        )

    rendered = pd.read_csv(
        rendered_manifest_path,
        dtype={"tic": "int64", "sector": "int64"},
    )
    rendered_columns = [
        "review_id",
        "tic",
        "twirl_vet_sheet_name",
        "width_px",
        "height_px",
        "size_bytes",
        "sha256",
    ]
    rendered = rendered.loc[:, rendered_columns].rename(
        columns={
            "twirl_vet_sheet_name": "source_vet_sheet_name",
            "width_px": "source_sheet_width_px",
            "height_px": "source_sheet_height_px",
            "size_bytes": "source_sheet_size_bytes",
            "sha256": "source_sheet_sha256",
        }
    )
    candidates = candidates.merge(
        rendered,
        on=["review_id", "tic"],
        how="left",
        validate="one_to_one",
    )

    candidates["source_sheet_path"] = candidates["source_vet_sheet_name"].map(
        lambda name: source_sheet_dir / str(name)
    )
    missing = ~candidates["source_sheet_path"].map(Path.exists)
    if missing.any():
        raise FileNotFoundError(
            "missing current vet sheets for: "
            + ", ".join(candidates.loc[missing, "review_id"].astype(str))
        )
    for _, row in candidates.iterrows():
        actual = _sha256(Path(row["source_sheet_path"]))
        if actual != str(row["source_sheet_sha256"]):
            raise RuntimeError(
                f"render-manifest SHA256 mismatch for {row['source_vet_sheet_name']}"
            )

    candidates["discovery_route"] = candidates["source_batch_id"].map(ORIGIN_LABELS)
    if candidates["discovery_route"].isna().any():
        unknown = sorted(candidates.loc[candidates["discovery_route"].isna(), "source_batch_id"].unique())
        raise RuntimeError(f"unrecognized S56 Planet-like source batches: {unknown}")
    candidates["known_benchmark"] = candidates["tic"].eq(BENCHMARK_TIC)
    candidates["object_name"] = np.where(
        candidates["known_benchmark"], "WD 1856 b", ""
    )
    candidates["morphology_fold_factor"] = pd.to_numeric(
        candidates["rereview_period_factor"], errors="coerce"
    ).fillna(1.0)
    candidates["morphology_view_period_d"] = (
        pd.to_numeric(candidates["anchor_period_d_current"], errors="coerce")
        * candidates["morphology_fold_factor"]
    )
    candidates["aperture_period_match"] = ~candidates[
        "aperture_disagreement_flag_current"
    ].map(_as_bool)
    candidates["min_cross_aperture_sde"] = candidates[
        ["adp_sml_sde_current", "adp_sde_current"]
    ].apply(pd.to_numeric, errors="coerce").min(axis=1)
    candidates["max_odd_even_sigma"] = candidates[
        [
            "adp_sml_own_even_odd_sigma_delta_current",
            "adp_own_even_odd_sigma_delta_current",
        ]
    ].apply(pd.to_numeric, errors="coerce").abs().max(axis=1)
    candidates["odd_even_review_flag"] = candidates["max_odd_even_sigma"].ge(3.0)
    candidates["priority_group"] = candidates.apply(_priority_group, axis=1)
    candidates["priority_note"] = candidates.apply(_priority_note, axis=1)

    group_order = {"Benchmark": 0, "A": 1, "B": 2, "C": 3}
    candidates["_group_order"] = candidates["priority_group"].map(group_order)
    candidates["_sort_sde"] = np.where(
        candidates["aperture_period_match"],
        candidates["min_cross_aperture_sde"],
        pd.to_numeric(candidates["adp_sml_sde_current"], errors="coerce"),
    )
    candidates = candidates.sort_values(
        ["_group_order", "odd_even_review_flag", "_sort_sde", "tic"],
        ascending=[True, True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    candidates["priority_rank"] = np.arange(1, len(candidates) + 1)
    return candidates


def _copy_sheets(candidates: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    out = candidates.copy()
    package_names: list[str] = []
    package_hashes: list[str] = []
    for _, row in out.iterrows():
        rank = int(row["priority_rank"])
        tic = int(row["tic"])
        token = str(row["priority_group"]).lower()
        package_name = f"{rank:02d}_TIC{tic}_{token}_current_adp.png"
        destination = out_dir / package_name
        shutil.copy2(Path(row["source_sheet_path"]), destination)
        package_names.append(package_name)
        package_hashes.append(_sha256(destination))
    out["vet_sheet_file"] = package_names
    out["vet_sheet_sha256"] = package_hashes
    if not (out["vet_sheet_sha256"] == out["source_sheet_sha256"]).all():
        raise RuntimeError("one or more copied vet sheets differ from the source rendering")
    return out


def _summary_columns(candidates: pd.DataFrame) -> pd.DataFrame:
    rename = {
        "cam_current": "camera",
        "ccd_current": "ccd",
        "tmag_current": "tmag",
        "anchor_period_d_current": "bls_anchor_period_d",
        "anchor_t0_bjd_current": "bls_anchor_t0_bjd",
        "anchor_duration_min_current": "duration_min",
        "adp_sml_depth_current": "adp_small_depth_fraction",
        "adp_depth_current": "adp_primary_depth_fraction",
        "adp_sml_sde_current": "adp_small_sde",
        "adp_sde_current": "adp_primary_sde",
        "aperture_period_rel_delta_current": "aperture_period_rel_delta",
        "aperture_depth_ratio_primary_over_small_current": (
            "aperture_depth_ratio_primary_over_small"
        ),
        "adp_sml_own_even_odd_sigma_delta_current": (
            "adp_small_odd_even_sigma"
        ),
        "adp_own_even_odd_sigma_delta_current": (
            "adp_primary_odd_even_sigma"
        ),
    }
    superseded = [
        column
        for column in (
            "cam",
            "ccd",
            "tmag",
            "duration_min",
            "aperture_period_rel_delta",
        )
        if column in candidates.columns
    ]
    out = candidates.drop(columns=superseded).rename(columns=rename)
    columns = [
        "priority_rank",
        "priority_group",
        "tic",
        "object_name",
        "sector",
        "discovery_route",
        "source_batch_id",
        "human_label",
        "morphology_review_status",
        "camera",
        "ccd",
        "tmag",
        "bls_anchor_period_d",
        "bls_anchor_t0_bjd",
        "morphology_fold_factor",
        "morphology_view_period_d",
        "rereview_period_status",
        "duration_min",
        "adp_small_depth_fraction",
        "adp_primary_depth_fraction",
        "adp_small_sde",
        "adp_primary_sde",
        "min_cross_aperture_sde",
        "aperture_period_rel_delta",
        "aperture_period_match",
        "aperture_depth_ratio_primary_over_small",
        "adp_small_odd_even_sigma",
        "adp_primary_odd_even_sigma",
        "max_odd_even_sigma",
        "odd_even_review_flag",
        "tier1_target_qa_status",
        "tier1_target_qa_reasons",
        "tier1_target_qa_pass",
        "tier1_target_searchable",
        "n_finite_quality0_min",
        "usable_cadence_fraction_min",
        "priority_note",
        "vet_sheet_file",
        "vet_sheet_sha256",
        "source_vet_sheet_name",
        "source_sheet_sha256",
        "review_id",
        "candidate_key",
        "final_labeler",
        "final_updated_utc",
        "human_label_source",
        "source_label_provenance_sha256",
    ]
    return out.loc[:, columns]


def _render_summary_table(summary: pd.DataFrame, out_dir: Path) -> None:
    apply_twirl_style("full_page")
    fig, ax = plt.subplots(figsize=(11.5, 15.5))
    ax.axis("off")
    counts = summary["priority_group"].value_counts()
    aligned = int(summary["aperture_period_match"].sum())
    human = int(summary["discovery_route"].eq("human-adjudicated").sum())
    model = int(summary["discovery_route"].str.startswith("model-enriched").sum())
    ax.text(
        0.0,
        1.0,
        "Sector 56 current Planet-like transit candidates",
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
        transform=ax.transAxes,
    )
    ax.text(
        0.0,
        0.976,
        (
            f"30 human-accepted morphology candidates: {human} original human-adjudicated, "
            f"{model} model-enriched; {aligned} current two-aperture period matches.\n"
            f"Groups: benchmark={counts.get('Benchmark', 0)}, A={counts.get('A', 0)}, "
            f"B={counts.get('B', 0)}, C={counts.get('C', 0)}."
        ),
        ha="left",
        va="top",
        fontsize=8.5,
        color="0.25",
        transform=ax.transAxes,
    )
    ax.text(
        0.0,
        0.956,
        (
            "Priority is transparent triage, not astrophysical validation: A has aligned "
            "current ADP periods; B is small-aperture-led; C adds a Tier-1 review warning.\n"
            "Morphology fold factors are evidence-view choices, not claimed orbital periods."
        ),
        ha="left",
        va="top",
        fontsize=7.8,
        color="0.30",
        transform=ax.transAxes,
    )

    table_data = []
    for _, row in summary.iterrows():
        route = "Human" if row["discovery_route"] == "human-adjudicated" else "Model→human"
        table_data.append(
            [
                int(row["priority_rank"]),
                str(row["priority_group"]),
                str(int(row["tic"])),
                route,
                f"{float(row['tmag']):.2f}",
                f"{float(row['bls_anchor_period_d']):.6g}",
                f"{float(row['morphology_view_period_d']):.6g}",
                f"{float(row['duration_min']):.0f}",
                f"{100.0 * float(row['adp_small_depth_fraction']):.1f}",
                f"{float(row['adp_small_sde']):.1f}",
                "yes" if _as_bool(row["aperture_period_match"]) else "no",
                str(row["tier1_target_qa_status"]),
            ]
        )
    labels = [
        "Rank",
        "Group",
        "TIC",
        "Route",
        "Tmag",
        "BLS P\n[d]",
        "View P\n[d]",
        "Dur.\n[min]",
        "Small depth\n[%]",
        "Small\nSDE",
        "ADP P\nmatch",
        "Target QA",
    ]
    table = ax.table(
        cellText=table_data,
        colLabels=labels,
        cellLoc="center",
        colLoc="center",
        bbox=[0.0, 0.02, 1.0, 0.90],
        colWidths=[0.045, 0.055, 0.095, 0.11, 0.06, 0.085, 0.085, 0.06, 0.09, 0.06, 0.07, 0.085],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.1)
    header_color = "#3B4D66"
    group_colors = {
        "Benchmark": "#DCE6F1",
        "A": "#E7F2EC",
        "B": "#FFF4D6",
        "C": "#F8E4E4",
    }
    for (row_idx, col_idx), cell in table.get_celld().items():
        cell.set_edgecolor("#D4D7DC")
        cell.set_linewidth(0.35)
        if row_idx == 0:
            cell.set_facecolor(header_color)
            cell.set_text_props(color="white", weight="bold")
            cell.set_height(0.036)
        else:
            group = str(summary.iloc[row_idx - 1]["priority_group"])
            cell.set_facecolor(group_colors[group] if col_idx in {0, 1} else "white")
            cell.set_height(0.027)
            if col_idx in {2, 3}:
                cell.set_text_props(ha="left")

    ax.text(
        0.0,
        0.0,
        (
            "Source: frozen S56–S62 morphology decisions accepted 2026-07-24; "
            "uniform s56_s62_a2v1_current_adp_v1 sheet set (S56-ADP-HV2)."
        ),
        ha="left",
        va="bottom",
        fontsize=7,
        color="0.35",
        transform=ax.transAxes,
    )
    fig.savefig(
        out_dir / "s56_planet_like_candidates_summary.png",
        dpi=220,
        bbox_inches="tight",
    )
    fig.savefig(
        out_dir / "s56_planet_like_candidates_summary.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)


def _render_contact_pages(
    candidates: pd.DataFrame,
    out_dir: Path,
    *,
    page_size: int = 6,
) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    title_font = _font(21, bold=True)
    body_font = _font(17)
    paths: list[Path] = []
    for page_index, start in enumerate(range(0, len(candidates), page_size), start=1):
        page_rows = candidates.iloc[start : start + page_size]
        panels: list[Image.Image] = []
        for _, row in page_rows.iterrows():
            image = Image.open(Path(row["source_sheet_path"])).convert("RGB")
            width = 620
            scale = width / image.width
            image = image.resize(
                (width, int(round(image.height * scale))),
                Image.Resampling.LANCZOS,
            )
            header_height = 62
            panel = Image.new("RGB", (width, header_height + image.height), "white")
            draw = ImageDraw.Draw(panel)
            draw.rectangle(
                [0, 0, panel.width - 1, panel.height - 1],
                outline=(40, 40, 40),
                width=2,
            )
            draw.text(
                (10, 7),
                f"#{int(row['priority_rank']):02d}  TIC {int(row['tic'])}",
                fill=(25, 25, 25),
                font=title_font,
            )
            draw.text(
                (10, 35),
                (
                    f"Group {row['priority_group']} · "
                    f"{row['discovery_route']} · "
                    f"ADP period match={row['aperture_period_match']}"
                ),
                fill=(65, 65, 65),
                font=body_font,
            )
            panel.paste(image, (0, header_height))
            panels.append(panel)

        columns = 3
        rows = 2
        gutter = 12
        width = columns * 620 + (columns - 1) * gutter
        panel_height = max(panel.height for panel in panels)
        canvas = Image.new(
            "RGB",
            (width, rows * panel_height + (rows - 1) * gutter),
            "white",
        )
        for index, panel in enumerate(panels):
            x = (index % columns) * (620 + gutter)
            y = (index // columns) * (panel_height + gutter)
            canvas.paste(panel, (x, y))
        path = out_dir / f"s56_planet_like_contact_page_{page_index:02d}.png"
        canvas.save(path)
        paths.append(path)
    return paths


def _write_readme(summary: pd.DataFrame, out_dir: Path, *, sheet_set: dict[str, Any]) -> None:
    group_counts = summary["priority_group"].value_counts()
    text = f"""# Sector 56 current Planet-like candidate bundle

This folder contains the newest uniform vetting sheet for every final,
human-accepted `planet_like` morphology candidate in Sector 56. The source
decision snapshot is the accepted `2026-07-24` S56--S62 morphology corpus;
this subset contains **{len(summary)} unique S56 TICs**:

- `{int(summary["discovery_route"].eq("human-adjudicated").sum())}` from the
  original human-adjudicated S56 set.
- `{int(summary["discovery_route"].str.startswith("model-enriched").sum())}`
  from the later model-enriched compact revisit, each subsequently accepted by
  the human morphology pass.
- `{int(summary["aperture_period_match"].sum())}` with aligned periods in both
  current ADP apertures and
  `{int((~summary["aperture_period_match"]).sum())}` where the current primary
  aperture prefers a different period.

## Priority groups

- **Benchmark ({group_counts.get("Benchmark", 0)}):** WD 1856 b, the confirmed
  planet and mandatory engineering benchmark.
- **A ({group_counts.get("A", 0)}):** current `DET_FLUX_ADP_SML` and
  `DET_FLUX_ADP` period solutions are aligned.
- **B ({group_counts.get("B", 0)}):** Planet-like morphology in the small
  search aperture, but the current primary aperture prefers a different
  period.
- **C ({group_counts.get("C", 0)}):** the same aperture disagreement plus a
  non-blocking Tier-1 target-QA review warning.

The ordering is a transparent follow-up triage, not an astrophysical
classification or confirmation. Within A, rows are ordered by the smaller of
the two aperture SDE values. Within B/C, rows are ordered by small-aperture SDE.
No model probability enters the ordering. Odd/even and target-QA flags remain
explicit columns in the summary.

## Sheet contract

- Sheet-set version: `{sheet_set["sheet_set_version"]}`
- Renderer: `{sheet_set["renderer_version"]}`
- Branch: `{sheet_set["branch_name"]}`
- Apertures: `{", ".join(sheet_set["apertures"])}`
- Trial periods: `{sheet_set["n_periods"]:,}`
- Retained peaks: `{sheet_set["n_peaks"]}`
- Row ephemeris preserved: `{sheet_set["use_row_ephemeris"]}`

The `morphology_fold_factor` and `morphology_view_period_d` fields describe the
evidence view accepted in morphology review. They do not by themselves claim
the physical orbital period.

## Files

- `s56_planet_like_candidates.xlsx`: sortable summary workbook.
- `s56_planet_like_candidates.csv`: machine-readable candidate table.
- `s56_planet_like_candidates_summary.png` and `.pdf`: one-page visual index.
- `01_TIC...png` through `30_TIC...png`: hash-verified copies of the newest
  current-ADP vetting sheets.
- `manifest.json`: source and output checksums plus selection counts.
"""
    (out_dir / "README.md").write_text(text)


def build_package(
    *,
    accepted_path: Path,
    metrics_path: Path,
    eligibility_path: Path,
    rendered_manifest_path: Path,
    sheet_set_path: Path,
    source_sheet_dir: Path,
    out_dir: Path,
    qa_contact_dir: Path | None,
) -> dict[str, Any]:
    for path in (
        accepted_path,
        metrics_path,
        eligibility_path,
        rendered_manifest_path,
        sheet_set_path,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not source_sheet_dir.is_dir():
        raise FileNotFoundError(source_sheet_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = _build_candidate_table(
        accepted_path=accepted_path,
        metrics_path=metrics_path,
        eligibility_path=eligibility_path,
        rendered_manifest_path=rendered_manifest_path,
        source_sheet_dir=source_sheet_dir,
    )
    candidates = _copy_sheets(candidates, out_dir)
    summary = _summary_columns(candidates)
    csv_path = out_dir / "s56_planet_like_candidates.csv"
    json_path = out_dir / "s56_planet_like_candidates.json"
    summary.to_csv(csv_path, index=False)
    json_path.write_text(
        json.dumps(
            summary.replace({np.nan: None}).to_dict(orient="records"),
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    _render_summary_table(summary, out_dir)
    sheet_set = json.loads(sheet_set_path.read_text())
    _write_readme(summary, out_dir, sheet_set=sheet_set)

    contact_pages: list[Path] = []
    if qa_contact_dir is not None:
        contact_pages = _render_contact_pages(candidates, qa_contact_dir)

    source_paths = {
        "accepted_morphology": accepted_path,
        "current_render_metrics": metrics_path,
        "tier1_target_eligibility": eligibility_path,
        "rendered_sheet_manifest": rendered_manifest_path,
        "sheet_set": sheet_set_path,
    }
    output_paths = [
        csv_path,
        json_path,
        out_dir / "s56_planet_like_candidates_summary.png",
        out_dir / "s56_planet_like_candidates_summary.pdf",
        out_dir / "README.md",
        *[out_dir / name for name in summary["vet_sheet_file"]],
    ]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection": {
            "sector": 56,
            "final_morphology_label": "planet_like",
            "n_candidates": int(len(summary)),
            "n_unique_tics": int(summary["tic"].nunique()),
            "discovery_route_counts": summary["discovery_route"].value_counts().to_dict(),
            "priority_group_counts": summary["priority_group"].value_counts().to_dict(),
            "n_aperture_period_matches": int(summary["aperture_period_match"].sum()),
            "n_tier1_searchable": int(
                summary["tier1_target_searchable"].map(_as_bool).sum()
            ),
            "ranking_contract": (
                "benchmark first; aligned-current-ADP periods before disagreement; "
                "within aligned rows sort by min cross-aperture SDE; within "
                "disagreement rows sort by small-aperture SDE; odd/even review "
                "flags sort after unflagged rows; no model probability used"
            ),
        },
        "sheet_set": sheet_set,
        "sources": {
            name: {
                "path": str(path.relative_to(REPO_ROOT)),
                "sha256": _sha256(path),
            }
            for name, path in source_paths.items()
        },
        "outputs": {
            path.name: {
                "path": str(path.relative_to(REPO_ROOT)),
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in output_paths
        },
        "qa_contact_pages": [str(path) for path in contact_pages],
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest["selection"], indent=2, sort_keys=True))
    return manifest


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted", type=Path, default=DEFAULT_ACCEPTED)
    parser.add_argument("--metrics", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--eligibility", type=Path, default=DEFAULT_ELIGIBILITY)
    parser.add_argument(
        "--rendered-manifest",
        type=Path,
        default=DEFAULT_RENDERED_MANIFEST,
    )
    parser.add_argument("--sheet-set", type=Path, default=DEFAULT_SHEET_SET)
    parser.add_argument("--source-sheet-dir", type=Path, default=DEFAULT_SOURCE_SHEETS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--qa-contact-dir", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    build_package(
        accepted_path=args.accepted.resolve(),
        metrics_path=args.metrics.resolve(),
        eligibility_path=args.eligibility.resolve(),
        rendered_manifest_path=args.rendered_manifest.resolve(),
        sheet_set_path=args.sheet_set.resolve(),
        source_sheet_dir=args.source_sheet_dir.resolve(),
        out_dir=args.out_dir.resolve(),
        qa_contact_dir=args.qa_contact_dir.resolve() if args.qa_contact_dir else None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
