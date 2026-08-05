#!/usr/bin/env python3
"""Fail closed unless a selected native-v2 teacher ensemble is fully provenanced."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_training import (  # noqa: E402
    verify_deployed_teacher_checkpoint_manifest,
    verify_selected_teacher_checkpoint_manifest,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--training-table", type=Path)
    parser.add_argument("--native-h5", type=Path)
    parser.add_argument("--checkpoint-root", type=Path)
    parser.add_argument("--expected-profile")
    args = parser.parse_args()
    source_mode = args.training_table is not None or args.native_h5 is not None
    deployed_mode = args.checkpoint_root is not None
    if source_mode == deployed_mode:
        parser.error(
            "choose exactly one mode: --training-table plus --native-h5, "
            "or --checkpoint-root"
        )
    if source_mode:
        if args.training_table is None or args.native_h5 is None:
            parser.error("source mode requires both --training-table and --native-h5")
        manifest = verify_selected_teacher_checkpoint_manifest(
            manifest_path=args.manifest,
            training_table=args.training_table,
            native_h5=args.native_h5,
        )
    else:
        manifest = verify_deployed_teacher_checkpoint_manifest(
            manifest_path=args.manifest,
            checkpoint_root=args.checkpoint_root,
        )
    if (
        args.expected_profile is not None
        and manifest["selected_profile"] != args.expected_profile
    ):
        raise RuntimeError(
            "selected checkpoint profile mismatch: "
            f"manifest={manifest['selected_profile']!r}, "
            f"expected={args.expected_profile!r}"
        )
    print(
        json.dumps(
            {
                "checkpoint_namespace": manifest["checkpoint_namespace"],
                "input_contract_version": manifest["input_contract_version"],
                "native_h5_sha256": manifest["native_h5_sha256"],
                "training_table_sha256": manifest["training_table_sha256"],
                "selected_profile": manifest["selected_profile"],
                "n_checkpoints": len(manifest["checkpoints"]),
                "passed": True,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
