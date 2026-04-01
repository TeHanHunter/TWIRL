#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.catalogs import DEFAULT_BUILD_VERSION, build_master_catalog, write_master_catalog


DEFAULT_INPUT_CATALOG = Path("data_local/catalogs/GaiaEDR3_WD_main.fits")
DEFAULT_OUTPUT_DIR = Path("data_local/catalogs/twirl_master_catalog")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build the first-pass TWIRL WD master catalog from the local "
            "Gentile Fusillo Gaia EDR3 seed catalog."
        )
    )
    parser.add_argument(
        "--input-catalog",
        type=Path,
        default=DEFAULT_INPUT_CATALOG,
        help="Path to the local Gaia EDR3 WD FITS catalog.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the built master catalog and sidecar manifest.",
    )
    parser.add_argument(
        "--output-name",
        type=str,
        default="twirl_wd_master_catalog_v0.fits",
        help="Filename for the built FITS master catalog.",
    )
    parser.add_argument(
        "--manifest-name",
        type=str,
        default="twirl_wd_master_catalog_v0_manifest.json",
        help="Filename for the JSON sidecar manifest.",
    )
    parser.add_argument(
        "--build-version",
        type=str,
        default=DEFAULT_BUILD_VERSION,
        help="Version tag recorded in the output metadata.",
    )
    parser.add_argument(
        "--highconf-pwd-threshold",
        type=float,
        default=0.75,
        help="Pwd threshold used for the is_highconf_wd convenience flag.",
    )
    parser.add_argument(
        "--input-sha256",
        action="store_true",
        help="Compute and record the SHA256 of the input seed catalog.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_path = args.output_dir / args.output_name
    manifest_path = args.output_dir / args.manifest_name

    table, manifest = build_master_catalog(
        input_path=args.input_catalog,
        build_version=args.build_version,
        highconf_pwd_threshold=args.highconf_pwd_threshold,
        include_input_sha256=args.input_sha256,
        output_path=output_path,
    )
    write_master_catalog(
        table=table,
        manifest=manifest,
        output_path=output_path,
        manifest_path=manifest_path,
        overwrite=args.overwrite,
    )

    print(f"Wrote master catalog: {output_path}")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Rows: {len(table)}")


if __name__ == "__main__":
    main()
