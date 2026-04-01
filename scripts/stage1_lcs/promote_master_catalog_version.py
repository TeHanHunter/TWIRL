#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_OUTPUT_DIR = Path("data_local/catalogs/twirl_master_catalog")
DEFAULT_INPUT_CATALOG = DEFAULT_OUTPUT_DIR / "twirl_wd_master_catalog_v0_tesscoverage.fits"
DEFAULT_INPUT_STAGE_MANIFEST = (
    DEFAULT_OUTPUT_DIR / "twirl_wd_master_catalog_v0_tesscoverage_manifest.json"
)
DEFAULT_INPUT_OBSERVATIONS = DEFAULT_OUTPUT_DIR / "twirl_wd_tess_observations_v0.fits"
DEFAULT_INPUT_DETECTOR_SUMMARY = DEFAULT_OUTPUT_DIR / "twirl_wd_tess_detector_summary_v0.csv"
DEFAULT_INPUT_SECTOR_SUMMARY = DEFAULT_OUTPUT_DIR / "twirl_wd_tess_sector_summary_v0.csv"
DEFAULT_INPUT_DETECTOR_DIR = DEFAULT_OUTPUT_DIR / "tess_detector_target_tables_v0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Promote an accepted intermediate TWIRL master-catalog state into a canonical "
            "versioned release such as twirl_wd_master_catalog_v1.fits."
        )
    )
    parser.add_argument(
        "--input-catalog",
        type=Path,
        default=DEFAULT_INPUT_CATALOG,
        help="Coverage-enriched FITS catalog to promote into the next canonical version.",
    )
    parser.add_argument(
        "--input-stage-manifest",
        type=Path,
        default=DEFAULT_INPUT_STAGE_MANIFEST,
        help="Stage manifest that describes how the input catalog was produced.",
    )
    parser.add_argument(
        "--input-observations",
        type=Path,
        default=DEFAULT_INPUT_OBSERVATIONS,
        help="Version-aligned one-row-per-hit observation table to promote alongside the catalog.",
    )
    parser.add_argument(
        "--input-detector-summary",
        type=Path,
        default=DEFAULT_INPUT_DETECTOR_SUMMARY,
        help="Detector summary CSV to promote alongside the catalog.",
    )
    parser.add_argument(
        "--input-sector-summary",
        type=Path,
        default=DEFAULT_INPUT_SECTOR_SUMMARY,
        help="Sector summary CSV to promote alongside the catalog.",
    )
    parser.add_argument(
        "--input-detector-dir",
        type=Path,
        default=DEFAULT_INPUT_DETECTOR_DIR,
        help="Per-detector target-table directory to promote alongside the catalog.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory that will hold the canonical versioned release products.",
    )
    parser.add_argument(
        "--catalog-version",
        type=str,
        default="v1",
        help="Canonical catalog version label, for example v1 or v2.",
    )
    parser.add_argument(
        "--mode",
        choices=("copy", "move"),
        default="copy",
        help="Copy promoted files by default; use move after you are ready to retire the intermediates.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing versioned outputs.",
    )
    return parser.parse_args()


def validate_version(version: str) -> None:
    if not re.fullmatch(r"v\d+", version):
        raise ValueError(f"Invalid catalog version '{version}'. Expected a label like v1 or v2.")


def canonical_paths(output_dir: Path, version: str) -> dict[str, Path]:
    return {
        "catalog": output_dir / f"twirl_wd_master_catalog_{version}.fits",
        "manifest": output_dir / f"twirl_wd_master_catalog_{version}_manifest.json",
        "observations": output_dir / f"twirl_wd_tess_observations_{version}.fits",
        "detector_summary": output_dir / f"twirl_wd_tess_detector_summary_{version}.csv",
        "sector_summary": output_dir / f"twirl_wd_tess_sector_summary_{version}.csv",
        "detector_dir": output_dir / f"tess_detector_target_tables_{version}",
    }


def ensure_output_slot(path: Path, overwrite: bool) -> None:
    if path.exists():
        if not overwrite:
            raise FileExistsError(f"{path} already exists. Use --overwrite to replace it.")
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
    path.parent.mkdir(parents=True, exist_ok=True)


def transfer_path(source: Path, destination: Path, mode: str, overwrite: bool) -> None:
    if not source.exists():
        raise FileNotFoundError(f"Missing input path: {source}")

    ensure_output_slot(destination, overwrite=overwrite)
    if mode == "copy":
        if source.is_dir():
            shutil.copytree(source, destination)
        else:
            shutil.copy2(source, destination)
        return

    shutil.move(str(source), str(destination))


def load_stage_manifest(path: Path) -> dict | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def main() -> None:
    args = parse_args()
    validate_version(args.catalog_version)

    outputs = canonical_paths(args.output_dir, args.catalog_version)
    stage_manifest = load_stage_manifest(args.input_stage_manifest)

    transfer_path(args.input_catalog, outputs["catalog"], mode=args.mode, overwrite=args.overwrite)
    transfer_path(
        args.input_observations,
        outputs["observations"],
        mode=args.mode,
        overwrite=args.overwrite,
    )
    transfer_path(
        args.input_detector_summary,
        outputs["detector_summary"],
        mode=args.mode,
        overwrite=args.overwrite,
    )
    transfer_path(
        args.input_sector_summary,
        outputs["sector_summary"],
        mode=args.mode,
        overwrite=args.overwrite,
    )
    transfer_path(
        args.input_detector_dir,
        outputs["detector_dir"],
        mode=args.mode,
        overwrite=args.overwrite,
    )

    manifest = {
        "promoted_at_utc": datetime.now(timezone.utc).isoformat(),
        "catalog_version": args.catalog_version,
        "promotion_mode": args.mode,
        "source_paths": {
            "catalog": str(args.input_catalog.resolve()),
            "stage_manifest": str(args.input_stage_manifest.resolve()),
            "observations": str(args.input_observations.resolve()),
            "detector_summary": str(args.input_detector_summary.resolve()),
            "sector_summary": str(args.input_sector_summary.resolve()),
            "detector_dir": str(args.input_detector_dir.resolve()),
        },
        "output_paths": {
            key: str(path.resolve()) for key, path in outputs.items()
        },
        "upstream_stage_manifest": stage_manifest,
        "notes": [
            "Canonical TWIRL master-catalog releases should use versioned filenames rather than chained stage suffixes.",
            "Stage-specific filenames remain acceptable as temporary intermediates during a run.",
            "The promoted catalog FITS is expected to already contain TIC and TESS coverage columns before this step.",
        ],
    }
    with outputs["manifest"].open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"[promote] wrote catalog: {outputs['catalog']}")
    print(f"[promote] wrote manifest: {outputs['manifest']}")
    print(f"[promote] wrote observations: {outputs['observations']}")
    print(f"[promote] wrote detector summary: {outputs['detector_summary']}")
    print(f"[promote] wrote sector summary: {outputs['sector_summary']}")
    print(f"[promote] wrote detector table dir: {outputs['detector_dir']}")


if __name__ == "__main__":
    main()
