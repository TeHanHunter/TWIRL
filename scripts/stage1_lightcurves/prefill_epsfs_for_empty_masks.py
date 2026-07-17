#!/usr/bin/env python3
"""Reuse existing ePSFs for cutouts whose saturated-pixel mask is empty.

A2v1 changes the ePSF fit by passing ``source.mask.mask`` into the pixel mask.
For cutouts where that mask has no true pixels, the fit is mathematically the
same as the existing pre-A2v1 ePSF.  This helper pre-fills a new production tree
with links or copies to those old ePSFs, leaving masked cutouts absent so the
normal ``tglc epsfs`` stage recomputes only the cases affected by the new mask.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import pickle
import re
import shutil
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import NamedTuple

import numpy as np


DEFAULT_SOURCE_TGLC_DATA_DIR = Path("/pdo/users/tehan/tglc-gpu-production")
DEFAULT_OUTPUT_TGLC_DATA_DIR = Path("/pdo/users/tehan/tglc-gpu-production-A2v1")
PDO_USER_ROOT = Path("/pdo/users/tehan")
SOURCE_RE = re.compile(r"^source_(?P<x>\d+)_(?P<y>\d+)\.pkl$")


@dataclass(frozen=True)
class Ccd:
    camera: int
    ccd: int


@dataclass
class CcdSummary:
    orbit: int
    camera: int
    ccd: int
    source_dir: str
    old_epsf_dir: str
    output_epsf_dir: str
    n_source_files: int = 0
    n_existing_output_epsfs: int = 0
    n_missing_old_epsfs: int = 0
    n_masked_sources: int = 0
    n_empty_mask_sources: int = 0
    n_links_written: int = 0
    n_copies_written: int = 0
    n_errors: int = 0


class SourceTask(NamedTuple):
    source_file: str
    old_epsf: str
    output_epsf: str
    apply: bool
    copy: bool
    overwrite: bool


class SourceResult(NamedTuple):
    status: str
    message: str = ""


def parse_ccd(value: str) -> Ccd:
    try:
        camera_text, ccd_text = value.split(",", maxsplit=1)
        camera = int(camera_text)
        ccd = int(ccd_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"CCD must have format camera,ccd; got {value!r}"
        ) from exc
    if camera not in {1, 2, 3, 4} or ccd not in {1, 2, 3, 4}:
        raise argparse.ArgumentTypeError(f"Camera/CCD must be in 1..4; got {value!r}")
    return Ccd(camera=camera, ccd=ccd)


def default_ccds() -> list[Ccd]:
    return [Ccd(camera=camera, ccd=ccd) for camera in range(1, 5) for ccd in range(1, 5)]


def enforce_pdo_write_path(path: Path) -> None:
    resolved = Path(path).expanduser().resolve(strict=False)
    if not str(resolved).startswith("/pdo/"):
        return
    allowed = PDO_USER_ROOT.resolve(strict=False)
    if resolved != allowed and allowed not in resolved.parents:
        raise ValueError(f"Refusing to write outside {allowed}: {resolved}")


def parse_source_grid(path: Path) -> tuple[int, int]:
    match = SOURCE_RE.match(path.name)
    if match is None:
        raise ValueError(f"Not a TGLC source pickle name: {path.name}")
    return int(match.group("x")), int(match.group("y"))


def source_has_masked_pixels(source_file: Path) -> bool:
    with source_file.open("rb") as source_pickle:
        source = pickle.load(source_pickle)
    return bool(np.ma.getmaskarray(source.mask).any())


def link_or_copy(old_epsf: Path, output_epsf: Path, *, copy: bool) -> str:
    output_epsf.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_epsf.with_name(output_epsf.name + f".tmp-{os.getpid()}")
    if tmp.exists() or tmp.is_symlink():
        tmp.unlink()
    if copy:
        shutil.copy2(old_epsf, tmp)
        action = "copy"
    else:
        os.symlink(old_epsf.resolve(strict=True), tmp)
        action = "link"
    tmp.rename(output_epsf)
    return action


def process_source_task(task: SourceTask) -> SourceResult:
    source_file = Path(task.source_file)
    old_epsf = Path(task.old_epsf)
    output_epsf = Path(task.output_epsf)
    try:
        if output_epsf.exists() or output_epsf.is_symlink():
            if not task.overwrite:
                return SourceResult("existing_output")
            output_epsf.unlink()

        if not old_epsf.exists():
            return SourceResult("missing_old")

        if source_has_masked_pixels(source_file):
            return SourceResult("masked")

        if task.apply:
            action = link_or_copy(old_epsf, output_epsf, copy=task.copy)
            return SourceResult(action)
        return SourceResult("empty_mask")
    except Exception as exc:  # noqa: BLE001 - batch audit should continue.
        return SourceResult("error", f"{source_file}: {exc}")


def process_ccd(
    *,
    source_tglc_data_dir: Path,
    output_tglc_data_dir: Path,
    orbit: int,
    ccd: Ccd,
    apply: bool,
    copy: bool,
    overwrite: bool,
    nprocs: int,
) -> CcdSummary:
    ccd_rel = Path(f"orbit-{orbit}") / "ffi" / f"cam{ccd.camera}" / f"ccd{ccd.ccd}"
    source_dir = output_tglc_data_dir / ccd_rel / "source"
    old_epsf_dir = source_tglc_data_dir / ccd_rel / "epsf"
    output_epsf_dir = output_tglc_data_dir / ccd_rel / "epsf"
    summary = CcdSummary(
        orbit=orbit,
        camera=ccd.camera,
        ccd=ccd.ccd,
        source_dir=str(source_dir),
        old_epsf_dir=str(old_epsf_dir),
        output_epsf_dir=str(output_epsf_dir),
    )

    if not source_dir.exists():
        return summary

    source_files = sorted(source_dir.glob("source_*.pkl"))
    summary.n_source_files = len(source_files)

    tasks = []
    for source_file in source_files:
        x, y = parse_source_grid(source_file)
        tasks.append(
            SourceTask(
                source_file=str(source_file),
                old_epsf=str(old_epsf_dir / f"epsf_{x}_{y}.npy"),
                output_epsf=str(output_epsf_dir / f"epsf_{x}_{y}.npy"),
                apply=apply,
                copy=copy,
                overwrite=overwrite,
            )
        )

    if nprocs > 1 and len(tasks) > 1:
        with multiprocessing.Pool(nprocs) as pool:
            results = pool.imap_unordered(process_source_task, tasks)
            for result in results:
                update_summary_from_result(summary, result)
    else:
        for task in tasks:
            update_summary_from_result(summary, process_source_task(task))

    return summary


def update_summary_from_result(summary: CcdSummary, result: SourceResult) -> None:
    if result.status == "existing_output":
        summary.n_existing_output_epsfs += 1
    elif result.status == "missing_old":
        summary.n_missing_old_epsfs += 1
    elif result.status == "masked":
        summary.n_masked_sources += 1
    elif result.status == "empty_mask":
        summary.n_empty_mask_sources += 1
    elif result.status == "link":
        summary.n_empty_mask_sources += 1
        summary.n_links_written += 1
    elif result.status == "copy":
        summary.n_empty_mask_sources += 1
        summary.n_copies_written += 1
    elif result.status == "error":
        summary.n_errors += 1
        print(f"[prefill-epsf] ERROR {result.message}", file=sys.stderr)
    else:
        summary.n_errors += 1
        print(f"[prefill-epsf] ERROR unknown status {result}", file=sys.stderr)


def write_summary(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-tglc-data-dir",
        type=Path,
        default=DEFAULT_SOURCE_TGLC_DATA_DIR,
        help="Existing TGLC root containing reusable ePSFs.",
    )
    parser.add_argument(
        "--output-tglc-data-dir",
        type=Path,
        default=DEFAULT_OUTPUT_TGLC_DATA_DIR,
        help="New TGLC root to pre-fill.",
    )
    parser.add_argument("--orbit", type=int, required=True)
    parser.add_argument(
        "--ccd",
        action="append",
        type=parse_ccd,
        help="Camera,CCD pair. Repeatable. Defaults to all 16 CCDs.",
    )
    parser.add_argument("--apply", action="store_true", help="Write links/copies.")
    parser.add_argument("--copy", action="store_true", help="Copy instead of symlink.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing output ePSFs. Do not use while TGLC is running.",
    )
    parser.add_argument(
        "--nprocs",
        type=int,
        default=1,
        help="Number of source-mask readers per CCD.",
    )
    parser.add_argument("--summary-json", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    enforce_pdo_write_path(args.output_tglc_data_dir)
    ccds = args.ccd if args.ccd else default_ccds()

    summaries = [
        process_ccd(
            source_tglc_data_dir=args.source_tglc_data_dir,
            output_tglc_data_dir=args.output_tglc_data_dir,
            orbit=args.orbit,
            ccd=ccd,
            apply=args.apply,
            copy=args.copy,
            overwrite=args.overwrite,
            nprocs=args.nprocs,
        )
        for ccd in ccds
    ]

    totals = {
        "n_ccds": len(summaries),
        "n_source_files": sum(item.n_source_files for item in summaries),
        "n_existing_output_epsfs": sum(item.n_existing_output_epsfs for item in summaries),
        "n_missing_old_epsfs": sum(item.n_missing_old_epsfs for item in summaries),
        "n_masked_sources": sum(item.n_masked_sources for item in summaries),
        "n_empty_mask_sources": sum(item.n_empty_mask_sources for item in summaries),
        "n_links_written": sum(item.n_links_written for item in summaries),
        "n_copies_written": sum(item.n_copies_written for item in summaries),
        "n_errors": sum(item.n_errors for item in summaries),
    }
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_tglc_data_dir": str(args.source_tglc_data_dir),
        "output_tglc_data_dir": str(args.output_tglc_data_dir),
        "orbit": args.orbit,
        "apply": bool(args.apply),
        "copy": bool(args.copy),
        "overwrite": bool(args.overwrite),
        "nprocs": int(args.nprocs),
        "totals": totals,
        "ccds": [asdict(item) for item in summaries],
    }

    for item in summaries:
        print(
            "[prefill-epsf] "
            f"orbit={item.orbit} cam{item.camera}/ccd{item.ccd} "
            f"sources={item.n_source_files} existing={item.n_existing_output_epsfs} "
            f"empty_mask={item.n_empty_mask_sources} masked={item.n_masked_sources} "
            f"links={item.n_links_written} copies={item.n_copies_written} "
            f"missing_old={item.n_missing_old_epsfs} errors={item.n_errors}"
        )
    print("[prefill-epsf] totals " + json.dumps(totals, sort_keys=True))

    if args.summary_json:
        write_summary(args.summary_json, payload)


if __name__ == "__main__":
    main()
