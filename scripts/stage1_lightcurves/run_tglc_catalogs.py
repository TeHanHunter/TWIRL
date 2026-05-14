#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path

from astropy.table import Table


DEFAULT_DETECTOR_SUMMARY = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_tess_detector_summary_v0.csv"
)
DEFAULT_TGLC_DATA_DIR = Path("data_local/tglc-data")
DEFAULT_COMMAND_SCRIPT = Path(
    "data_local/catalogs/twirl_master_catalog/tglc_catalog_jobs_v0.sh"
)
DEFAULT_JOB_TABLE = Path(
    "data_local/catalogs/twirl_master_catalog/tglc_catalog_jobs_v0.csv"
)


@dataclass(frozen=True)
class CatalogJob:
    orbit: int
    sector: int
    camera: int
    ccd: int
    n_targets_all: int
    n_targets_unique_tic: int
    n_targets_highconf: int
    n_targets_highconf_unique_tic: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plan or run MIT `tglc catalogs` jobs from the TWIRL orbit-aware detector summary."
        )
    )
    parser.add_argument(
        "--detector-summary",
        type=Path,
        default=DEFAULT_DETECTOR_SUMMARY,
        help="Orbit-aware TWIRL detector summary CSV.",
    )
    parser.add_argument(
        "--tglc-data-dir",
        type=Path,
        default=DEFAULT_TGLC_DATA_DIR,
        help="Path to the MIT `tglc-data` tree.",
    )
    parser.add_argument(
        "--output-script",
        type=Path,
        default=DEFAULT_COMMAND_SCRIPT,
        help="Optional shell script path for the generated commands.",
    )
    parser.add_argument(
        "--output-job-table",
        type=Path,
        default=DEFAULT_JOB_TABLE,
        help="Optional CSV path for the filtered job table.",
    )
    parser.add_argument(
        "--max-magnitude",
        type=float,
        default=20.0,
        help="`tglc catalogs --max-magnitude` value.",
    )
    parser.add_argument(
        "--nprocs",
        type=int,
        default=16,
        help="`tglc catalogs --nprocs` value.",
    )
    parser.add_argument(
        "--orbit",
        type=int,
        action="append",
        default=None,
        help="Restrict to one or more specific orbits. Repeat to pass multiple values.",
    )
    parser.add_argument(
        "--sector",
        type=int,
        action="append",
        default=None,
        help="Restrict to one or more specific sectors. Repeat to pass multiple values.",
    )
    parser.add_argument(
        "--camera",
        type=int,
        action="append",
        default=None,
        help="Restrict to one or more specific cameras. Repeat to pass multiple values.",
    )
    parser.add_argument(
        "--ccd",
        type=int,
        action="append",
        default=None,
        help="Restrict to one or more specific CCDs. Repeat to pass multiple values.",
    )
    parser.add_argument(
        "--min-targets-all",
        type=int,
        default=1,
        help="Keep only jobs with at least this many TWIRL targets.",
    )
    parser.add_argument(
        "--min-targets-highconf",
        type=int,
        default=0,
        help="Keep only jobs with at least this many `Pwd > 0.75` targets.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on the number of selected jobs.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip jobs whose `TIC_camX_ccdY.ecsv` and `Gaia_camX_ccdY.ecsv` already exist.",
    )
    parser.add_argument(
        "--write-script",
        action="store_true",
        help="Write the generated commands to `--output-script`.",
    )
    parser.add_argument(
        "--write-job-table",
        action="store_true",
        help="Write the filtered job list to `--output-job-table`.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Execute the generated commands sequentially instead of only printing them.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting `--output-script` and `--output-job-table`.",
    )
    return parser.parse_args()


def progress(message: str) -> None:
    print(message, flush=True)


def load_detector_summary(path: Path) -> list[CatalogJob]:
    table = Table.read(path, format="ascii.csv")
    table.convert_bytestring_to_unicode()

    required = [
        "orbit",
        "sector",
        "camera",
        "ccd",
        "n_targets_all",
        "n_targets_unique_tic",
        "n_targets_highconf",
        "n_targets_highconf_unique_tic",
    ]
    missing = [name for name in required if name not in table.colnames]
    if missing:
        raise ValueError(
            f"{path} is missing required columns for orbit-aware job planning: {missing}"
        )

    jobs: list[CatalogJob] = []
    for row in table:
        jobs.append(
            CatalogJob(
                orbit=int(row["orbit"]),
                sector=int(row["sector"]),
                camera=int(row["camera"]),
                ccd=int(row["ccd"]),
                n_targets_all=int(row["n_targets_all"]),
                n_targets_unique_tic=int(row["n_targets_unique_tic"]),
                n_targets_highconf=int(row["n_targets_highconf"]),
                n_targets_highconf_unique_tic=int(row["n_targets_highconf_unique_tic"]),
            )
        )
    return jobs


def catalog_paths_for_job(tglc_data_dir: Path, job: CatalogJob) -> tuple[Path, Path]:
    catalog_dir = tglc_data_dir / f"orbit-{job.orbit}" / "ffi" / "catalogs"
    tic_path = catalog_dir / f"TIC_cam{job.camera}_ccd{job.ccd}.ecsv"
    gaia_path = catalog_dir / f"Gaia_cam{job.camera}_ccd{job.ccd}.ecsv"
    return tic_path, gaia_path


def keep_job(job: CatalogJob, args: argparse.Namespace) -> bool:
    if args.orbit is not None and job.orbit not in set(args.orbit):
        return False
    if args.sector is not None and job.sector not in set(args.sector):
        return False
    if args.camera is not None and job.camera not in set(args.camera):
        return False
    if args.ccd is not None and job.ccd not in set(args.ccd):
        return False
    if job.n_targets_all < int(args.min_targets_all):
        return False
    if job.n_targets_highconf < int(args.min_targets_highconf):
        return False
    return True


def command_for_job(job: CatalogJob, args: argparse.Namespace) -> list[str]:
    return [
        "tglc",
        "catalogs",
        "--orbit",
        str(job.orbit),
        "--ccd",
        f"{job.camera},{job.ccd}",
        "--max-magnitude",
        str(args.max_magnitude),
        "--nprocs",
        str(args.nprocs),
        "--tglc-data-dir",
        str(args.tglc_data_dir),
    ]


def format_command(argv: list[str]) -> str:
    return " ".join(shlex.quote(arg) for arg in argv)


def write_command_script(path: Path, commands: list[tuple[CatalogJob, list[str]]], overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Use --overwrite to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n\n")
        for job, argv in commands:
            handle.write(
                f"# orbit={job.orbit} sector={job.sector} cam={job.camera} ccd={job.ccd} "
                f"targets_all={job.n_targets_all} highconf={job.n_targets_highconf}\n"
            )
            handle.write(format_command(argv))
            handle.write("\n\n")

    path.chmod(0o755)


def write_job_table(path: Path, jobs: list[CatalogJob], overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Use --overwrite to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "orbit",
                "sector",
                "camera",
                "ccd",
                "n_targets_all",
                "n_targets_unique_tic",
                "n_targets_highconf",
                "n_targets_highconf_unique_tic",
            ]
        )
        for job in jobs:
            writer.writerow(
                [
                    job.orbit,
                    job.sector,
                    job.camera,
                    job.ccd,
                    job.n_targets_all,
                    job.n_targets_unique_tic,
                    job.n_targets_highconf,
                    job.n_targets_highconf_unique_tic,
                ]
            )


def main() -> None:
    args = parse_args()

    jobs = load_detector_summary(args.detector_summary)
    progress(f"[tglc-catalogs] loaded {len(jobs)} orbit/camera/ccd jobs from {args.detector_summary}")

    filtered_jobs = [job for job in jobs if keep_job(job, args)]
    progress(f"[tglc-catalogs] {len(filtered_jobs)} jobs remain after filtering")

    if args.skip_existing:
        kept: list[CatalogJob] = []
        skipped = 0
        for job in filtered_jobs:
            tic_path, gaia_path = catalog_paths_for_job(args.tglc_data_dir, job)
            if tic_path.exists() and gaia_path.exists():
                skipped += 1
                continue
            kept.append(job)
        filtered_jobs = kept
        progress(f"[tglc-catalogs] skipped {skipped} jobs with existing TIC+Gaia catalogs")

    if args.limit is not None:
        filtered_jobs = filtered_jobs[: args.limit]
        progress(f"[tglc-catalogs] limited to {len(filtered_jobs)} jobs")

    commands = [(job, command_for_job(job, args)) for job in filtered_jobs]

    if args.write_script:
        write_command_script(args.output_script, commands, overwrite=args.overwrite)
        progress(f"[tglc-catalogs] wrote command script: {args.output_script}")

    if args.write_job_table:
        write_job_table(args.output_job_table, filtered_jobs, overwrite=args.overwrite)
        progress(f"[tglc-catalogs] wrote job table: {args.output_job_table}")

    if not commands:
        progress("[tglc-catalogs] no jobs selected")
        return

    progress("[tglc-catalogs] selected jobs:")
    for job, argv in commands:
        progress(
            "[tglc-catalogs] "
            f"orbit={job.orbit} sector={job.sector} cam={job.camera} ccd={job.ccd} "
            f"targets_all={job.n_targets_all} highconf={job.n_targets_highconf}"
        )
        progress(format_command(argv))

    if args.run:
        for idx, (job, argv) in enumerate(commands, start=1):
            progress(
                "[tglc-catalogs] "
                f"running {idx}/{len(commands)}: orbit={job.orbit} sector={job.sector} "
                f"cam={job.camera} ccd={job.ccd}"
            )
            subprocess.run(argv, check=True)
        progress(f"[tglc-catalogs] completed {len(commands)} jobs")


if __name__ == "__main__":
    main()
