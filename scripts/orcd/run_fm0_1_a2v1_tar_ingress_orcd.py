#!/usr/bin/env python3
"""Transfer checksum-bound FM0.1 sector tars from Blackhole to ORCD."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import time
from typing import Any


BLACKHOLE_ID = "d8ea14bc-dca1-4cbc-92b6-b76d7289b6d2"
ORCD_ID = "ec54b570-cac5-47f7-b2a1-100c2078686f"
SOURCE_ROOT = "/twirl/fm0/source_a2v1_s56_s65_tar_v1"
DEST_ROOT = Path(
    "/orcd/data/mki_aryeh/001/twirl/exports/fm0/source_archives/"
    "a2v1_s56_s65_tar_v1"
)
STATE_PATH = Path(
    "/orcd/data/mki_aryeh/001/twirl/logs/twirl_fm0_1_s56_s65_tar_ingress/"
    "state.json"
)
GLOBUS = Path.home() / ".local/share/twirl-globus-cli/bin/globus"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def run(
    args: list[str], *, check: bool = True, cwd: Path | None = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args, check=check, text=True, capture_output=True, cwd=cwd
    )


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    if temporary.exists():
        raise RuntimeError(f"stale state temporary: {temporary}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def load_state(path: Path, sectors: list[int]) -> dict[str, Any]:
    if path.exists():
        payload = json.loads(path.read_text())
        if payload.get("sectors") != sectors:
            raise RuntimeError("existing ingress state has a different sector range")
        return payload
    return {
        "schema_version": "twirl_fm0_1_a2v1_tar_ingress_v1",
        "sectors": sectors,
        "created_at": utc_now(),
        "items": {},
    }


def source_ready(tag: str) -> bool:
    result = run(
        [
            str(GLOBUS), "ls", f"{BLACKHOLE_ID}:{SOURCE_ROOT}/{tag}.partial/",
            "--filter", "=READY_BLACKHOLE",
        ],
        check=False,
    )
    return result.returncode == 0 and "READY_BLACKHOLE" in result.stdout


def submit(tag: str, item: dict[str, Any]) -> None:
    submission = json.loads(
        run([str(GLOBUS), "task", "generate-submission-id", "-F", "json"]).stdout
    )["value"]
    result = run(
        [
            str(GLOBUS), "transfer", BLACKHOLE_ID, ORCD_ID,
            f"{SOURCE_ROOT}/{tag}.partial/", f"{DEST_ROOT}/{tag}/",
            "--recursive", "--submission-id", submission,
            "--label", f"twirl-fm0-{tag}-tar-ingress",
            "--sync-level", "checksum", "--verify-checksum", "--encrypt",
            "--notify", "off", "-F", "json",
        ]
    )
    item.update(
        {
            "submission_id": submission,
            "task_id": json.loads(result.stdout)["task_id"],
            "submitted_at": utc_now(),
        }
    )


def wait_task(item: dict[str, Any]) -> None:
    task_id = str(item["task_id"])
    while True:
        payload = json.loads(
            run([str(GLOBUS), "task", "show", task_id, "-F", "json"]).stdout
        )
        status = payload["status"]
        if status == "ACTIVE":
            time.sleep(60)
            continue
        if status != "SUCCEEDED":
            raise RuntimeError(f"Globus task {task_id} ended as {status}: {payload}")
        unsafe = (
            "faults", "files_skipped", "subtasks_failed", "subtasks_retrying",
            "subtasks_skipped_errors", "subtasks_canceled", "subtasks_expired",
        )
        if any(int(payload.get(name, 0)) for name in unsafe):
            raise RuntimeError(f"Globus task {task_id} has unsafe counters: {payload}")
        if not payload.get("encrypt_data") or not payload.get("verify_checksum"):
            raise RuntimeError(f"Globus task {task_id} lost transfer protections")
        item["transfer_status"] = status
        item["completed_at"] = utc_now()
        item["bytes_transferred"] = int(payload["bytes_transferred"])
        item["files_transferred"] = int(payload["files_transferred"])
        return


def verify_destination(tag: str, item: dict[str, Any]) -> None:
    root = DEST_ROOT / tag
    archive = f"{tag}_A2v1_raw_hdf5.tar"
    result = run(
        ["sha256sum", "-c", f"{archive}.sha256"], check=False, cwd=root
    )
    if result.returncode != 0:
        raise RuntimeError(f"ORCD archive checksum failed for {tag}: {result.stdout}{result.stderr}")
    ready = root / "READY_ORCD"
    if not ready.exists():
        ready.write_text(json.dumps({"tag": tag, "verified_at": utc_now(), "task_id": item["task_id"]}, sort_keys=True) + "\n")
    item["verified_at"] = utc_now()
    item["ready_orcd"] = str(ready)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--from-sector", type=int, default=56)
    parser.add_argument("--through-sector", type=int, default=65)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    if not args.apply:
        raise SystemExit("refusing transfer without --apply")
    if args.from_sector < 56 or args.through_sector > 65 or args.through_sector < args.from_sector:
        raise SystemExit("sector range must be within S56-S65")
    if not GLOBUS.is_file():
        raise SystemExit(f"missing Globus CLI: {GLOBUS}")
    sectors = list(range(args.from_sector, args.through_sector + 1))
    state = load_state(STATE_PATH, sectors)
    DEST_ROOT.mkdir(parents=True, exist_ok=True)
    for sector in sectors:
        tag = f"s{sector:04d}"
        item = state["items"].setdefault(tag, {})
        while not source_ready(tag):
            time.sleep(args.poll_seconds)
        if "task_id" not in item:
            submit(tag, item)
            atomic_json(STATE_PATH, state)
        if item.get("transfer_status") != "SUCCEEDED":
            wait_task(item)
            atomic_json(STATE_PATH, state)
        if not item.get("verified_at"):
            verify_destination(tag, item)
            atomic_json(STATE_PATH, state)
        print(f"ORCD_TAR_READY {tag}", flush=True)
    state["completed_at"] = utc_now()
    atomic_json(STATE_PATH, state)
    print("ORCD_FM0_TARS_READY S56-S65", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
