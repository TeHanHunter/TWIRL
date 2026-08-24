#!/usr/bin/env python3
"""Build a fail-closed Stage-1-safe admission receipt for one FM0.2 H200 job."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re


def _field(text: str, name: str) -> str:
    match = re.search(r"(?:^|\s)" + re.escape(name) + r"=([^\s]+)", text)
    if not match:
        raise ValueError(f"missing scheduler field: {name}")
    return match.group(1)


def _tres(text: str, name: str, *, zero_if_absent: bool = False) -> int:
    match = re.search(r"(?:^|,)" + re.escape(name) + r"=([0-9]+)", text)
    if not match:
        if zero_if_absent:
            return 0
        raise ValueError(f"missing TRES: {name}")
    return int(match.group(1))


def _memory_mib(value: str) -> int:
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([KMGTP]?)", value)
    if not match:
        raise ValueError(f"unparseable memory: {value}")
    scale = {
        "": 1 / 1024**2,
        "K": 1 / 1024,
        "M": 1,
        "G": 1024,
        "T": 1024**2,
        "P": 1024**3,
    }[match.group(2)]
    return int(float(match.group(1)) * scale)


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_receipt(args: argparse.Namespace) -> dict[str, object]:
    paths = {
        name: Path(getattr(args, name.replace("-", "_")))
        for name in (
            "partition",
            "reclaimable_partition",
            "preemption",
            "node",
            "queue",
            "node_queue",
        )
    }
    raw = {name: path.read_text().strip() for name, path in paths.items()}
    if any(not value for name, value in raw.items() if name not in {"queue", "node_queue"}):
        raise ValueError("required scheduler snapshot is empty")
    if _field(raw["partition"], "PartitionName") != "pg_mki_aryeh" or _field(raw["partition"], "State") != "UP":
        raise ValueError("primary partition is not the exact live authority")
    if _field(raw["reclaimable_partition"], "PartitionName") != "mit_preemptable" or _field(raw["reclaimable_partition"], "State") != "UP":
        raise ValueError("reclaimable partition is not the exact live authority")
    if int(_field(raw["partition"], "PriorityTier")) <= int(_field(raw["reclaimable_partition"], "PriorityTier")):
        raise ValueError("FM partition no longer outranks reclaimable jobs")
    if _field(raw["reclaimable_partition"], "PreemptMode") != "REQUEUE":
        raise ValueError("reclaimable partition no longer requeues")
    if not re.search(r"^PreemptType\s*=\s*preempt/partition_prio$", raw["preemption"], re.MULTILINE):
        raise ValueError("partition-priority preemption is not active")
    if _field(raw["node"], "NodeName") != "node4900":
        raise ValueError("wrong GPU node")
    bad_states = {"DOWN", "DRAIN", "DRAINED", "FAIL", "FAILING", "MAINT", "UNKNOWN"}
    states = set(_field(raw["node"], "State").split("+"))
    if states & bad_states:
        raise ValueError(f"GPU node is ineligible: {sorted(states)}")

    cfg_tres = _field(raw["node"], "CfgTRES")
    alloc_tres = _field(raw["node"], "AllocTRES")
    capacity = {
        "h200": _tres(cfg_tres, "gres/gpu:h200"),
        "cpus": _tres(cfg_tres, "cpu"),
        "memory_mib": int(_field(raw["node"], "RealMemory")),
    }
    allocated = {
        "h200": _tres(alloc_tres, "gres/gpu:h200", zero_if_absent=True),
        "cpus": _tres(alloc_tres, "cpu"),
        "memory_mib": int(_field(raw["node"], "AllocMem")),
    }
    if capacity["h200"] != 8:
        raise ValueError(f"node4900 H200 count drift: {capacity['h200']}")

    stage1_pattern = re.compile(r"^twirl-s[0-9]+-o[0-9]+-c[1-4]d[1-4]-a2v1(?:-p4)?$")
    active_stage1: list[dict[str, object]] = []
    held_stage1: list[dict[str, object]] = []
    runnable_stage1: list[dict[str, object]] = []
    active = {"h200": 0, "cpus": 0, "memory_mib": 0}
    for line in raw["queue"].splitlines():
        fields = line.split("|", 9)
        if len(fields) != 10:
            raise ValueError(f"unparseable partition queue row: {line}")
        jobid, user, state, reason, cpus, memory, gres, name, nodes, dependency = fields
        if name.startswith("twirl-fm0-2-"):
            raise ValueError(f"another FM0.2 job is already live: {jobid}/{state}")
        if not stage1_pattern.fullmatch(name):
            continue
        row = {
            "job_id": jobid,
            "user": user,
            "state": state,
            "reason": reason,
            "cpus": int(cpus),
            "memory": memory,
            "gres": gres,
            "nodes": nodes,
            "dependency": dependency,
        }
        if state in {"RUNNING", "CONFIGURING", "COMPLETING"}:
            match = re.search(r"h200:([0-9]+)", gres.lower())
            if not match or "node4900" not in nodes:
                raise ValueError(f"active Stage-1 resource ambiguity: {jobid}")
            active_stage1.append(row)
            active["h200"] += int(match.group(1))
            active["cpus"] += int(cpus)
            active["memory_mib"] += _memory_mib(memory)
        elif state == "PENDING" and reason == "Dependency":
            held_stage1.append(row)
        elif state == "PENDING":
            runnable_stage1.append(row)
        else:
            raise ValueError(f"unknown Stage-1 scheduler state: {jobid}/{state}")
    if runnable_stage1:
        raise ValueError("runnable Stage-1 job is pending; FM admission denied")
    if active["h200"] > 4 or active["cpus"] > 78:
        raise ValueError("Stage-1 live ceiling already exceeded")

    reclaimable_jobs: list[dict[str, object]] = []
    reclaimable = {"h200": 0, "cpus": 0, "memory_mib": 0}
    for line in raw["node_queue"].splitlines():
        fields = line.split("|", 10)
        if len(fields) != 11:
            raise ValueError(f"unparseable node queue row: {line}")
        jobid, partition, user, state, reason, cpus, memory, gres, name, nodes, dependency = fields
        if partition != "mit_preemptable":
            continue
        if state not in {"RUNNING", "CONFIGURING", "COMPLETING"} or "node4900" not in nodes:
            raise ValueError(f"reclaimable job ambiguity: {jobid}/{state}/{nodes}")
        match = re.search(r"h200:([0-9]+)", gres.lower())
        h200 = int(match.group(1)) if match else 0
        reclaimable_jobs.append({
            "job_id": jobid,
            "partition": partition,
            "user": user,
            "state": state,
            "cpus": int(cpus),
            "memory": memory,
            "gres": gres,
            "name": name,
            "nodes": nodes,
            "dependency": dependency,
        })
        reclaimable["h200"] += h200
        reclaimable["cpus"] += int(cpus)
        reclaimable["memory_mib"] += _memory_mib(memory)

    other = {key: allocated[key] - active[key] - reclaimable[key] for key in allocated}
    if any(value < 0 for value in other.values()):
        raise ValueError("node/queue allocation accounting is inconsistent")
    stage1_reservation = {"h200": 4, "cpus": 64, "memory_mib": 1_572_864}
    requested = {"h200": 1, "cpus": 4, "memory_mib": 32_768}
    needed = {key: other[key] + stage1_reservation[key] + requested[key] for key in capacity}
    if any(needed[key] > capacity[key] for key in needed):
        raise ValueError(f"insufficient capacity after Stage-1 reservation: needed={needed}, capacity={capacity}")
    if args.run_kind == "synthetic_fp32" and args.stop_after_step != 8:
        raise ValueError("synthetic smoke must stop after step 8")
    if args.run_kind == "real_canary" and args.stop_after_step not in {64, 500, 1000, 2000}:
        raise ValueError("real canary stop is not authorized")

    return {
        "schema_version": "twirl_fm0_2_orcd_admission_receipt_v1",
        "passed": True,
        "created_at_epoch": args.now_epoch,
        "expires_at_epoch": args.now_epoch + 300,
        "ttl_seconds": 300,
        "run_kind": args.run_kind,
        "stop_after_step": args.stop_after_step,
        "expected_git_sha": args.expected_git_sha,
        "upstream_input_receipt": {"path": args.input_receipt, "sha256": args.input_receipt_sha256},
        "input_reuse_receipt": {"path": args.input_reuse_receipt, "sha256": args.input_reuse_receipt_sha256},
        "loader_restart_receipt": {"path": args.loader_restart_receipt, "sha256": args.loader_restart_receipt_sha256},
        "partition": "pg_mki_aryeh",
        "gpu_node": "node4900",
        "requested_fm_resources": requested,
        "stage1_reservation": {"lanes": 4, **stage1_reservation},
        "derived": {
            "node_capacity": capacity,
            "node_allocated": allocated,
            "active_stage1_allocated": active,
            "reclaimable_allocated": reclaimable,
            "non_reclaimable_non_stage1_allocated": other,
            "capacity_needed": needed,
            "active_stage1_jobs": active_stage1,
            "dependency_held_stage1_jobs": held_stage1,
            "runnable_stage1_jobs": runnable_stage1,
            "reclaimable_node_jobs": reclaimable_jobs,
        },
        "snapshot_sha256": {name: _digest(path) for name, path in paths.items()},
        "raw_snapshots": raw,
        "sealed_test_access_count": 0,
        "claim_limit": "One gated FM0.2.1 H200 invocation only",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    for name in ("partition", "reclaimable-partition", "preemption", "node", "queue", "node-queue"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--now-epoch", required=True, type=int)
    parser.add_argument("--expected-git-sha", required=True)
    parser.add_argument("--run-kind", required=True, choices=("synthetic_fp32", "real_canary"))
    parser.add_argument("--stop-after-step", required=True, type=int)
    parser.add_argument("--input-receipt", required=True)
    parser.add_argument("--input-receipt-sha256", required=True)
    parser.add_argument("--input-reuse-receipt", required=True)
    parser.add_argument("--input-reuse-receipt-sha256", required=True)
    parser.add_argument("--loader-restart-receipt", required=True)
    parser.add_argument("--loader-restart-receipt-sha256", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.expected_git_sha):
        raise SystemExit("expected Git SHA must be full lowercase hex")
    for value in (
        args.input_receipt_sha256,
        args.input_reuse_receipt_sha256,
        args.loader_restart_receipt_sha256,
    ):
        if not re.fullmatch(r"[0-9a-f]{64}", value):
            raise SystemExit("receipt SHA-256 is invalid")
    payload = build_receipt(args)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
