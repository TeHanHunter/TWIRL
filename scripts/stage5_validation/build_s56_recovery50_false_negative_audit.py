#!/usr/bin/env python3
"""Build a vet-sheet audit for human-planet rows missed by a CNN teacher."""
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_teacher import json_default, read_table  # noqa: E402


DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_2k"
DEFAULT_PREDICTIONS = DEFAULT_ROOT / "cnn_teacher_display_anchor/cnn_shape_plus_bls/cnn_predictions.parquet"
DEFAULT_TRAINING_TABLE = DEFAULT_ROOT / "human_training_table/human_vetting_training_table.csv"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "false_negative_planet_audit_display_anchor"
DEFAULT_SHEET_DIRS = (
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/twirl_vet_sheets_fullphase_binmatch",
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/twirl_vet_sheets",
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_next4k/twirl_vet_sheets",
)


def _font(size: int) -> ImageFont.ImageFont:
    for path in (
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ):
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def _find_sheet(row: pd.Series, sheet_dirs: tuple[Path, ...]) -> Path | None:
    names = []
    for col in ("display_vet_sheet_name", "twirl_vet_sheet_name"):
        value = str(row.get(col, "") or "").strip()
        if value and value.lower() != "nan":
            names.append(value)
    tic = str(row.get("tic", "")).split(".")[0]
    if tic:
        names.append(f"real_{tic}_twirl_twoap_current_adp.png")
    for name in names:
        for sheet_dir in sheet_dirs:
            path = sheet_dir / name
            if path.exists():
                return path
    return None


def _contact_sheet(rows: pd.DataFrame, out_png: Path, *, panel_width: int = 430, columns: int = 2) -> None:
    title_font = _font(18)
    body_font = _font(15)
    panels: list[Image.Image] = []
    for _, row in rows.iterrows():
        sheet_path = Path(row["vet_sheet_path"])
        image = Image.open(sheet_path).convert("RGB")
        scale = panel_width / image.width
        image = image.resize((panel_width, int(image.height * scale)), Image.Resampling.LANCZOS)
        header_h = 88
        panel = Image.new("RGB", (panel_width, header_h + image.height), "white")
        draw = ImageDraw.Draw(panel)
        draw.rectangle([0, 0, panel_width - 1, panel.height - 1], outline=(30, 30, 30), width=2)
        label = str(row.get("cnn_label", ""))
        p_planet = float(row.get("cnn_p_planet_like", np.nan))
        model_p = float(row.get("model_period_d", np.nan))
        queue_p = float(row.get("queue_period_d", np.nan))
        draw.text((12, 8), str(row["review_id"]), fill=(20, 20, 20), font=title_font)
        draw.text((12, 34), f"M:{label}  p_planet={p_planet:.3f}", fill=(70, 70, 70), font=body_font)
        draw.text((12, 58), f"P_model={model_p:.6g} d  P_queue={queue_p:.6g} d", fill=(70, 70, 70), font=body_font)
        panel.paste(image, (0, header_h))
        panels.append(panel)

    if not panels:
        out_png.parent.mkdir(parents=True, exist_ok=True)
        Image.new("RGB", (panel_width, 120), "white").save(out_png)
        return

    gutter = 18
    rows_n = int(np.ceil(len(panels) / columns))
    row_heights = [
        max(panel.height for panel in panels[row * columns:(row + 1) * columns])
        for row in range(rows_n)
    ]
    width = columns * panel_width + (columns - 1) * gutter
    height = sum(row_heights) + (rows_n - 1) * gutter
    canvas = Image.new("RGB", (width, height), "white")
    y = 0
    for row in range(rows_n):
        x = 0
        for panel in panels[row * columns:(row + 1) * columns]:
            canvas.paste(panel, (x, y))
            x += panel_width + gutter
        y += row_heights[row] + gutter
    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_png)


def _html_index(rows: pd.DataFrame, out_html: Path, contact_png: Path) -> None:
    lines = [
        "<!doctype html>",
        "<meta charset='utf-8'>",
        "<title>S56 Recovery50 Planet False Negatives</title>",
        "<style>body{font-family:Arial,sans-serif;margin:24px;max-width:1200px} img{max-width:100%} table{border-collapse:collapse} td,th{border:1px solid #bbb;padding:6px 8px;text-align:left}</style>",
        "<h1>S56 Recovery50 Planet False Negatives</h1>",
        f"<p>Contact sheet: <a href='{html.escape(contact_png.name)}'>{html.escape(contact_png.name)}</a></p>",
        "<table>",
        "<tr><th>review_id</th><th>machine</th><th>p_planet</th><th>P_model</th><th>P_queue</th><th>sheet</th></tr>",
    ]
    for _, row in rows.iterrows():
        sheet = Path(row["vet_sheet_path"])
        rel = sheet.name
        lines.append(
            "<tr>"
            f"<td>{html.escape(str(row['review_id']))}</td>"
            f"<td>{html.escape(str(row['cnn_label']))}</td>"
            f"<td>{float(row['cnn_p_planet_like']):.4f}</td>"
            f"<td>{float(row['model_period_d']):.6g}</td>"
            f"<td>{float(row['queue_period_d']):.6g}</td>"
            f"<td><a href='{html.escape(str(sheet))}'>{html.escape(rel)}</a></td>"
            "</tr>"
        )
    lines.extend(["</table>", f"<h2>Contact Sheet</h2><img src='{html.escape(contact_png.name)}'>"])
    out_html.write_text("\n".join(lines) + "\n")


def build_audit(
    *,
    predictions: Path,
    training_table: Path,
    out_dir: Path,
    sheet_dirs: tuple[Path, ...],
    split: str,
    target_label: str,
) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    pred = read_table(predictions)
    table = read_table(training_table)
    join_cols = [col for col in ("row_id", "review_id") if col in pred.columns and col in table.columns]
    keep = [
        col
        for col in table.columns
        if col in join_cols
        or col in {"display_vet_sheet_name", "twirl_vet_sheet_name", "human_label", "display_ephemeris_period_rel_delta"}
    ]
    rows = pred.merge(table.loc[:, keep], on=join_cols, how="left", suffixes=("", "_training"))
    rows = rows[rows["main_teacher_target"].eq(target_label) & ~rows["cnn_label"].eq(target_label)].copy()
    if split != "all":
        rows = rows[rows["cnn_training_split"].eq(split)].copy()
    rows = rows.sort_values("cnn_p_planet_like", ascending=True).reset_index(drop=True)
    rows["vet_sheet_path"] = [
        str(_find_sheet(row, sheet_dirs) or "")
        for _, row in rows.iterrows()
    ]
    missing = rows["vet_sheet_path"].eq("")
    if missing.any():
        missing_ids = ", ".join(rows.loc[missing, "review_id"].astype(str).tolist())
        raise FileNotFoundError(f"missing vet sheets for: {missing_ids}")

    out_csv = out_dir / "heldout_planet_false_negatives.csv"
    out_all_csv = out_dir / "all_split_planet_false_negatives.csv"
    rows.to_csv(out_csv, index=False)
    all_rows = pred[pred["main_teacher_target"].eq(target_label) & ~pred["cnn_label"].eq(target_label)].copy()
    all_rows.to_csv(out_all_csv, index=False)
    contact_png = out_dir / "heldout_planet_false_negative_vet_sheets.png"
    index_html = out_dir / "index.html"
    _contact_sheet(rows, contact_png)
    _html_index(rows, index_html, contact_png)

    summary = {
        "predictions": str(predictions),
        "training_table": str(training_table),
        "split": split,
        "target_label": target_label,
        "n_false_negatives": int(len(rows)),
        "n_all_split_false_negatives": int(len(all_rows)),
        "outputs": {
            "csv": str(out_csv),
            "all_split_csv": str(out_all_csv),
            "contact_sheet": str(contact_png),
            "index_html": str(index_html),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--sheet-dir", type=Path, action="append", default=list(DEFAULT_SHEET_DIRS))
    parser.add_argument("--split", default="test")
    parser.add_argument("--target-label", default="planet_like")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    build_audit(
        predictions=args.predictions,
        training_table=args.training_table,
        out_dir=args.out_dir,
        sheet_dirs=tuple(args.sheet_dir),
        split=args.split,
        target_label=args.target_label,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
