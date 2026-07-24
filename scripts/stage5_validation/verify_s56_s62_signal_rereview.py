#!/usr/bin/env python3
"""Verify queue identity, completion state, and local vet-sheet coverage."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, UnidentifiedImageError

from twirl.vetting.label_io import canonicalize_candidate_key
from twirl.vetting.multisector_signal_review import (
    FINAL_LABELS,
    SIGNAL_REREVIEW_LABEL_RANK,
    tehan_app_candidate_key,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_render_provenance(
    *,
    queue: pd.DataFrame,
    manifest: pd.DataFrame,
    provenance_root: Path,
    sheet_set: dict[str, object],
) -> tuple[list[str], int]:
    if not bool(sheet_set.get("require_render_provenance", False)):
        return [], 0
    errors: list[str] = []
    verified = 0
    expected_apertures = list(sheet_set.get("apertures", []))
    expected_harmonics = [
        float(value) for value in sheet_set.get("harmonic_factors", [])
    ]
    expected_n_periods = int(sheet_set.get("n_periods", -1))
    expected_n_peaks = int(sheet_set.get("n_peaks", -1))
    expected_exports = {
        str(key): str(value)
        for key, value in dict(
            sheet_set.get("lc_export_h5_by_sector", {})
        ).items()
    }
    for sector in sorted(pd.to_numeric(manifest["sector"], errors="raise").unique()):
        sector_int = int(sector)
        metrics_path = provenance_root / f"s{sector_int:04d}_metrics.csv"
        summary_path = provenance_root / f"s{sector_int:04d}_summary.json"
        if not metrics_path.is_file() or not summary_path.is_file():
            errors.append(f"S{sector_int}: missing renderer metrics or summary")
            continue
        metrics = pd.read_csv(metrics_path, dtype=str, keep_default_na=False)
        sector_manifest = manifest.loc[
            pd.to_numeric(manifest["sector"], errors="raise").eq(sector_int)
        ].copy()
        required_metrics = {
            "review_id",
            "tic",
            "sector",
            "twirl_vet_sheet_name",
            "twirl_vet_status",
            "anchor_period_d",
            "anchor_t0_bjd",
            "anchor_duration_min",
            "anchor_source",
            "vet_branch",
            "vet_sheet_version",
        }
        missing_columns = sorted(required_metrics - set(metrics.columns))
        if missing_columns:
            errors.append(
                f"S{sector_int}: renderer metrics missing columns {missing_columns}"
            )
            continue
        if len(metrics) != len(sector_manifest) or metrics["review_id"].duplicated().any():
            errors.append(
                f"S{sector_int}: renderer metrics do not cover the exact sector queue"
            )
            continue
        joined = sector_manifest.merge(
            metrics,
            on="review_id",
            how="left",
            suffixes=("_expected", "_rendered"),
            validate="one_to_one",
        )
        identity_mismatch = (
            joined["twirl_vet_status"].ne("ok")
            | pd.to_numeric(joined["tic_expected"], errors="coerce").ne(
                pd.to_numeric(joined["tic_rendered"], errors="coerce")
            )
            | pd.to_numeric(joined["sector_expected"], errors="coerce").ne(
                pd.to_numeric(joined["sector_rendered"], errors="coerce")
            )
            | joined["twirl_vet_sheet_name_expected"].ne(
                joined["twirl_vet_sheet_name_rendered"]
            )
            | joined["anchor_source"].ne("row_metadata")
            | joined["vet_branch"].ne(str(sheet_set.get("branch_name", "")))
            | joined["vet_sheet_version"].ne(
                str(sheet_set.get("renderer_version", ""))
            )
        )
        for expected_col, rendered_col in (
            ("period_d", "anchor_period_d"),
            ("t0_bjd", "anchor_t0_bjd"),
            ("duration_min", "anchor_duration_min"),
        ):
            expected = pd.to_numeric(joined[expected_col], errors="coerce")
            rendered = pd.to_numeric(joined[rendered_col], errors="coerce")
            identity_mismatch |= ~np.isclose(
                expected,
                rendered,
                rtol=0.0,
                atol=1.0e-10,
                equal_nan=False,
            )
        if identity_mismatch.any():
            errors.append(
                f"S{sector_int}: {int(identity_mismatch.sum())} rendered rows "
                "mismatch the frozen candidate identity or recipe"
            )
            continue
        render_summary = json.loads(summary_path.read_text())
        expected_status = {"ok": int(len(sector_manifest))}
        summary_ok = (
            int(render_summary.get("n_rows", -1)) == len(sector_manifest)
            and render_summary.get("status_counts") == expected_status
            and render_summary.get("branch_name")
            == sheet_set.get("branch_name")
            and render_summary.get("apertures") == expected_apertures
            and bool(render_summary.get("use_row_ephemeris"))
            == bool(sheet_set.get("use_row_ephemeris"))
            and bool(render_summary.get("write_pdf"))
            == bool(sheet_set.get("write_pdf"))
            and [
                float(value)
                for value in render_summary.get("harmonic_factors", [])
            ]
            == expected_harmonics
            and int(render_summary.get("n_periods", -1))
            == expected_n_periods
            and int(render_summary.get("n_peaks", -1)) == expected_n_peaks
            and str(render_summary.get("lc_export_h5", ""))
            == expected_exports.get(str(sector_int), "")
        )
        if not summary_ok:
            errors.append(f"S{sector_int}: renderer summary recipe mismatch")
            continue
        verified += len(sector_manifest)
    return errors, int(verified)


def verify(bundle_dir: Path) -> dict[str, object]:
    bundle_dir = Path(bundle_dir)
    queue_path = bundle_dir / "review_queue_planet_eb.csv"
    labels_path = bundle_dir / "human_labels_final.csv"
    manifest_path = bundle_dir / "required_sheet_manifest.csv"
    sheet_set_path = bundle_dir / "sheet_set.json"
    sheet_set = (
        json.loads(sheet_set_path.read_text())
        if sheet_set_path.is_file()
        else {
            "sheet_root": "vet_sheets",
            "render_provenance_root": "",
            "require_render_provenance": False,
        }
    )
    sheet_root = bundle_dir / str(sheet_set["sheet_root"])
    provenance_root = bundle_dir / str(
        sheet_set.get("render_provenance_root", "")
    )
    queue = pd.read_csv(queue_path, dtype=str, keep_default_na=False)
    labels = pd.read_csv(labels_path, dtype=str, keep_default_na=False)
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    required_queue = {
        "row_id",
        "candidate_key",
        "observation_candidate_key",
        "selected_source_uid",
        "twirl_vet_sheet_name",
    }
    missing_queue = sorted(required_queue - set(queue.columns))
    if missing_queue:
        raise KeyError(f"review queue is missing columns: {missing_queue}")
    required_manifest = {
        "row_id",
        "review_id",
        "sector",
        "tic",
        "period_d",
        "t0_bjd",
        "duration_min",
        "twirl_vet_sheet_name",
    }
    missing_manifest = sorted(required_manifest - set(manifest.columns))
    if missing_manifest:
        raise KeyError(f"sheet manifest is missing columns: {missing_manifest}")
    if queue["row_id"].duplicated().any() or manifest["row_id"].duplicated().any():
        raise ValueError("queue or sheet manifest has duplicate row_id values")
    required_labels = {
        "row_id",
        "candidate_key",
        "tic",
        "sector",
        "label",
        "labeler",
    }
    missing_labels = sorted(required_labels - set(labels.columns))
    if missing_labels:
        raise KeyError(f"final label file is missing columns: {missing_labels}")
    if (
        queue["candidate_key"].eq("").any()
        or queue["selected_source_uid"].eq("").any()
        or queue["candidate_key"].ne(queue["observation_candidate_key"]).any()
    ):
        raise ValueError(
            "review queue has a blank or inconsistent immutable candidate/source "
            "identity"
        )
    expected_order = (
        queue["prior_label"]
        .map(SIGNAL_REREVIEW_LABEL_RANK)
        .fillna(99)
        .astype(int)
        .to_numpy()
    )
    if np.any(expected_order[1:] < expected_order[:-1]):
        raise ValueError("review queue is not ordered Planet-like before EB")
    if set(queue["row_id"]) != set(manifest["row_id"]):
        raise ValueError("sheet manifest does not cover the exact queue")
    queue_binding = queue.loc[
        :,
        [
            "row_id",
            "review_id",
            "sector",
            "tic",
            "period_d",
            "t0_bjd",
            "duration_min",
            "twirl_vet_sheet_name",
        ],
    ].merge(
        manifest.loc[:, list(required_manifest)],
        on="row_id",
        how="left",
        suffixes=("_queue", "_manifest"),
        validate="one_to_one",
    )
    binding_mismatch = pd.Series(False, index=queue_binding.index)
    for column in ("review_id", "twirl_vet_sheet_name"):
        binding_mismatch |= queue_binding[f"{column}_queue"].ne(
            queue_binding[f"{column}_manifest"]
        )
    for column in ("sector", "tic", "period_d", "t0_bjd", "duration_min"):
        binding_mismatch |= ~np.isclose(
            pd.to_numeric(
                queue_binding[f"{column}_queue"], errors="coerce"
            ),
            pd.to_numeric(
                queue_binding[f"{column}_manifest"], errors="coerce"
            ),
            rtol=0.0,
            atol=1.0e-10,
            equal_nan=False,
        )
    if binding_mismatch.any():
        raise ValueError("sheet manifest does not match the exact review queue")
    if labels["row_id"].duplicated().any():
        raise ValueError("final label file has duplicate row_id values")
    if not set(labels["row_id"]).issubset(set(queue["row_id"])):
        raise ValueError("final label file contains unknown row_id values")
    nonblank_labels = labels["label"].astype(str).str.strip().ne("")
    invalid_labels = sorted(
        set(labels.loc[nonblank_labels, "label"]) - FINAL_LABELS
    )
    if invalid_labels:
        raise ValueError(
            f"final label file contains invalid labels: {invalid_labels}"
        )
    if not labels.empty:
        expected = queue.set_index("row_id").apply(
            tehan_app_candidate_key, axis=1
        )
        observed_key = labels["candidate_key"].map(canonicalize_candidate_key)
        expected_key = labels["row_id"].map(expected).map(
            canonicalize_candidate_key
        )
        mismatch = observed_key.ne(expected_key)
        if mismatch.any():
            raise ValueError("final label file candidate_key mismatch")
        expected_identity = queue.set_index("row_id")[["tic", "sector"]]
        for column in ("tic", "sector"):
            observed = pd.to_numeric(labels[column], errors="coerce")
            expected_value = pd.to_numeric(
                labels["row_id"].map(expected_identity[column]),
                errors="coerce",
            )
            if observed.ne(expected_value).any():
                raise ValueError(
                    f"final label file contains a {column} mismatch"
                )
    reviewed = (
        labels["label"].isin(FINAL_LABELS)
        & labels["labeler"].astype(str).str.strip().ne("")
    )
    missing: list[str] = []
    invalid_png: list[str] = []
    hash_rows: list[dict[str, object]] = []
    for _, asset in manifest.iterrows():
        name = str(asset["twirl_vet_sheet_name"])
        path = sheet_root / name
        if not path.is_file():
            missing.append(name)
            continue
        with path.open("rb") as handle:
            signature = handle.read(8)
        if path.stat().st_size < 8 or signature != b"\x89PNG\r\n\x1a\n":
            invalid_png.append(name)
            continue
        try:
            with Image.open(path) as image:
                width, height = image.size
                image_format = image.format
                image.verify()
            if image_format != "PNG" or width <= 0 or height <= 0:
                raise ValueError("invalid PNG format or dimensions")
        except (OSError, SyntaxError, UnidentifiedImageError, ValueError):
            invalid_png.append(name)
            continue
        hash_rows.append(
            {
                "row_id": asset["row_id"],
                "review_id": asset.get("review_id", ""),
                "sector": asset["sector"],
                "tic": asset["tic"],
                "twirl_vet_sheet_name": name,
                "width_px": int(width),
                "height_px": int(height),
                "size_bytes": int(path.stat().st_size),
                "sha256": _sha256(path),
            }
        )
    render_errors, n_render_verified = _verify_render_provenance(
        queue=queue,
        manifest=manifest,
        provenance_root=provenance_root,
        sheet_set=sheet_set,
    )
    pd.DataFrame(hash_rows).to_csv(
        bundle_dir / "rendered_sheet_manifest.csv", index=False
    )
    missing_path = bundle_dir / "missing_sheet_names.txt"
    missing_path.write_text("".join(f"{name}\n" for name in sorted(missing)))
    summary: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_queue_rows": int(len(queue)),
        "n_saved_rows": int(len(labels)),
        "n_reviewed_rows": int(reviewed.sum()),
        "n_pending_rows": int(len(queue) - reviewed.sum()),
        "n_blank_or_unattributed_saved_rows": int((~reviewed).sum()),
        "n_required_sheets": int(len(manifest)),
        "n_present_sheets": int(len(manifest) - len(missing)),
        "n_missing_sheets": int(len(missing)),
        "n_invalid_png": int(len(invalid_png)),
        "n_render_provenance_verified": n_render_verified,
        "render_provenance_errors": render_errors,
        "sheet_set_version": str(sheet_set.get("sheet_set_version", "")),
        "missing_sheet_names": missing,
        "invalid_png_names": invalid_png,
        "review_complete": bool(reviewed.sum() == len(queue)),
        "assets_complete": bool(
            not missing
            and not invalid_png
            and not render_errors
            and (
                not bool(sheet_set.get("require_render_provenance", False))
                or n_render_verified == len(manifest)
            )
        ),
    }
    (bundle_dir / "verification.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Return success while sheets or row reviews are still pending.",
    )
    args = parser.parse_args()
    summary = verify(args.bundle_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    complete = bool(summary["assets_complete"] and summary["review_complete"])
    if not complete and not args.allow_incomplete:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
