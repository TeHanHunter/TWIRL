#!/usr/bin/env bash
# Fail-closed local controller for the staged FM0.1 ORCD proof-of-concept.
# It is dry-run by default, never opens an interactive login, never submits a
# GPU array, and never chains a GPU job before its CPU/input gates have passed.
set -euo pipefail
umask 027

DRY_RUN=1
ACTION=""
ADMISSION_TMP=""
readonly PROJECT_ROOT="/orcd/data/mki_aryeh/001/twirl"
readonly PARTITION="pg_mki_aryeh"
readonly RECLAIMABLE_PARTITION="mit_preemptable"
readonly GPU_NODE="node4900"
readonly LOCKED_HOST="tehan@orcd-login.mit.edu"
readonly CONTROL_PATH="/Users/tehan/.ssh/cm/%r@%h:%p"
readonly DESIGN_SHA="94c8672a087884bf8a2c70f5d15315e05de602134af0c6c2073ca1c36232d6f7"
readonly CONFIG_SHA="7de4e8c9e98c0ce27648a21241acc51766e436a537ba932b54f529fbdbf26d8a"
readonly FREEZE_SHA="75092235978c0b582f569be55770ad2d63368e079a6777da0b1c44547f074c25"
readonly ORCD_CONFIG_SHA="582e9e1b47298f23c135f578cc7ddd8e256cea6035a85ea88028030e7b7d96a3"

usage() {
  cat <<'EOF'
Usage: TWIRL_EXPECTED_GIT_SHA=<40-hex> run_twirl_fm0_1_poc_orcd.sh [--execute] ACTION

Actions (run separately; no stage is auto-promoted):
  probe                    Check the existing user-opened socket and live partition.
  deploy                   Create/reuse a clean detached worktree at the exact SHA.
  build-env                Build/reuse the hash-bound FM0.1 environment.
  submit-sector-stage      Expand selected sources from one verified sector tar.
  submit-registry          Submit the CPU leakage-component registry build.
  submit-input-validation  Submit the CPU release build/independent freeze gate.
  submit-loader-smoke      Exercise two frozen shards through the official loader.
  submit-fp32-smoke        Create a <=5-minute all-queue admission receipt, then
                           submit exactly one H200 synthetic mechanics smoke.
  submit-real-train        Create the same fresh Stage-1-safe receipt, then
                           submit one H200 FM0.1.1 real-release training job.
  submit-post-validation   Independently validate the synthetic mechanics artifacts.
  submit-real-post-validation
                           Independently validate one immutable real-data run.
  status                   Read only the named run jobs and receipts.

Default mode only prints the remote commands. --execute is required for any
remote probe or mutation. This script never initiates password, Duo, or
keyboard-interactive authentication.

Required staged inputs for submit-registry/input-validation:
  TWIRL_FM0_ALIASES_TABLE, TWIRL_FM0_SOURCE_ADAPTER, and either
  (TWIRL_FM0_OBSERVATIONS_TABLE + TWIRL_FM0_SOURCE_MANIFEST) for the fixture,
  or TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY for the real A2v1 canary.

TWIRL_FM0_SOURCE_ADAPTER is exactly one of:
  strict_npz_fixture_v1       mechanics-only; never scientific-training eligible
  a2v1_hdf5_quality_aware_v1  checksum-bound A2v1 HDF5 canary/release

Required for submit-sector-stage:
  TWIRL_FM0_SECTOR, TWIRL_FM0_CORPUS_SELECTION,
  TWIRL_FM0_CORPUS_SELECTION_SHA256, TWIRL_FM0_SECTOR_ARCHIVE_DIR,
  TWIRL_FM0_QUALITY_TABLE, TWIRL_FM0_QUALITY_MANIFEST, and
  TWIRL_FM0_SECTOR_STAGE_ROOT. TWIRL_FM0_AFTEROK_JOB_ID is optional and
  serializes the sector behind one successful prior Slurm job.

Required immutable binding for loader and later stages:
  TWIRL_FM0_INPUT_RECEIPT, TWIRL_FM0_INPUT_RECEIPT_SHA256

Additional bindings:
  submit-fp32-smoke: TWIRL_FM0_LOADER_RECEIPT[_SHA256]
  submit-real-train: TWIRL_FM0_LOADER_RECEIPT[_SHA256],
                     TWIRL_FM0_TARGET_STEP (1--20000)
  submit-post-validation: TWIRL_FM0_SMOKE_OUTPUT,
                          TWIRL_FM0_ADMISSION_RECEIPT[_SHA256]
  submit-real-post-validation: TWIRL_FM0_TRAIN_OUTPUT,
                               TWIRL_FM0_RUN_GIT_SHA,
                               TWIRL_FM0_ADMISSION_RECEIPT[_SHA256]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --execute) DRY_RUN=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *)
      [[ -z "${ACTION}" ]] || { echo "Only one action is allowed." >&2; exit 2; }
      ACTION="$1"; shift ;;
  esac
done
[[ -n "${ACTION}" ]] || { usage >&2; exit 2; }
: "${TWIRL_EXPECTED_GIT_SHA:?set the exact reviewed 40-hex commit}"
readonly EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
[[ "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]] || { echo "Expected SHA must be full lowercase hex." >&2; exit 2; }

readonly LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
readonly REMOTE_SOURCE="${PROJECT_ROOT}/code/TWIRL"
readonly REMOTE_REPO="${PROJECT_ROOT}/code/TWIRL_fm0_${EXPECTED_SHA:0:12}"
readonly RUN_ROOT="${PROJECT_ROOT}/reports/stage5_validation/twirl_fm0_1_s56_s67_poc/${EXPECTED_SHA:0:12}"
readonly LAUNCH_DIR="${RUN_ROOT}/launch"
readonly ENV_PREFIX="${PROJECT_ROOT}/envs/twirl-fm0-poc-py3119-torch2110-cu128-v1"

SSH=(
  ssh -o BatchMode=yes -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no -o PubkeyAuthentication=no
  -o HostbasedAuthentication=no -o GSSAPIAuthentication=no
  -o NumberOfPasswordPrompts=0 -o ConnectTimeout=15
  -o ConnectionAttempts=1 -o ControlMaster=no -o "ControlPath=${CONTROL_PATH}"
)

print_cmd() { printf '+'; printf ' %q' "$@"; printf '\n'; }
remote() {
  if [[ "${DRY_RUN}" == 1 ]]; then print_cmd "${SSH[@]}" "${LOCKED_HOST}" "$1";
  else "${SSH[@]}" "${LOCKED_HOST}" "$1"; fi
}
capture() { "${SSH[@]}" "${LOCKED_HOST}" "$1"; }
require_socket() {
  [[ "${DRY_RUN}" == 1 ]] && return
  if ! "${SSH[@]}" -O check "${LOCKED_HOST}" >/dev/null 2>&1; then
    echo "ORCD control socket is absent/stale; open it yourself per doc/orcd_h200_usage.md." >&2
    exit 4
  fi
}
safe_remote_path() {
  local value=$1
  [[ "${value}" == "${PROJECT_ROOT}"/* && "${value}" != *".."* && "${value}" != *","* && "${value}" != *$'\n'* ]]
}
require_hash() { [[ "$2" =~ ^[0-9a-f]{${1}}$ ]] || { echo "Invalid ${1}-hex hash." >&2; exit 2; }; }
base_export() {
  printf '%s' "ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_RUN_ROOT=${RUN_ROOT}"
}
verify_remote_repo="
  [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
  [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
  [[ \$(sha256sum '${REMOTE_REPO}/doc/foundation_model_design.md' | awk '{print \$1}') == '${DESIGN_SHA}' ]]
  [[ \$(sha256sum '${REMOTE_REPO}/configs/models/twirl_fm0_1_s56_s67_poc.yaml' | awk '{print \$1}') == '${CONFIG_SHA}' ]]
  [[ \$(sha256sum '${REMOTE_REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json' | awk '{print \$1}') == '${FREEZE_SHA}' ]]
  [[ \$(sha256sum '${REMOTE_REPO}/configs/orcd/twirl_fm0_1_s56_s67_poc.json' | awk '{print \$1}') == '${ORCD_CONFIG_SHA}' ]]
  test -x '${ENV_PREFIX}/bin/python'
  test -s '${ENV_PREFIX}/build-manifest.json'
"

input_binding() {
  : "${TWIRL_FM0_INPUT_RECEIPT:?set exact immutable input receipt}"
  : "${TWIRL_FM0_INPUT_RECEIPT_SHA256:?set exact input receipt SHA256}"
  safe_remote_path "${TWIRL_FM0_INPUT_RECEIPT}" && [[ "${TWIRL_FM0_INPUT_RECEIPT}" == "${RUN_ROOT}"/* ]] || { echo "Unsafe input receipt path." >&2; exit 2; }
  require_hash 64 "${TWIRL_FM0_INPUT_RECEIPT_SHA256}"
}
loader_binding() {
  : "${TWIRL_FM0_LOADER_RECEIPT:?set exact immutable loader receipt}"
  : "${TWIRL_FM0_LOADER_RECEIPT_SHA256:?set exact loader receipt SHA256}"
  safe_remote_path "${TWIRL_FM0_LOADER_RECEIPT}" && [[ "${TWIRL_FM0_LOADER_RECEIPT}" == "${RUN_ROOT}"/* ]] || { echo "Unsafe loader receipt path." >&2; exit 2; }
  require_hash 64 "${TWIRL_FM0_LOADER_RECEIPT_SHA256}"
}
source_binding() {
  : "${TWIRL_FM0_SOURCE_ADAPTER:?set staged source adapter}"
  case "${TWIRL_FM0_SOURCE_ADAPTER}" in
    strict_npz_fixture_v1)
      : "${TWIRL_FM0_OBSERVATIONS_TABLE:?set staged admitted-observations table}"
      : "${TWIRL_FM0_SOURCE_MANIFEST:?set staged strict-NPZ source manifest}"
      safe_remote_path "${TWIRL_FM0_OBSERVATIONS_TABLE}" && safe_remote_path "${TWIRL_FM0_SOURCE_MANIFEST}" || { echo "Unsafe fixture source path." >&2; exit 2; }
      ;;
    a2v1_hdf5_quality_aware_v1)
      : "${TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY:?set staged A2v1 HDF5 source inventory}"
      safe_remote_path "${TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY}" || { echo "Unsafe A2v1 HDF5 source inventory path." >&2; exit 2; }
      ;;
    *) echo "Unsupported FM0.1 source adapter: ${TWIRL_FM0_SOURCE_ADAPTER}" >&2; exit 2 ;;
  esac
}
submit_once() {
  local record=$1 job_name=$2 wrapper=$3 exports=$4 dependency=${5:-}
  local dependency_option=""
  if [[ -n "${dependency}" ]]; then
    [[ "${dependency}" =~ ^[0-9]+$ ]] || { echo "Invalid afterok job ID." >&2; exit 2; }
    dependency_option="--dependency=afterok:${dependency}"
  fi
  remote "
    set -euo pipefail
    ${verify_remote_repo}
    mkdir -p '${LAUNCH_DIR}' '${PROJECT_ROOT}/logs'
    test ! -e '${record}'
    claim='${record}.claim'
    mkdir -- \"\${claim}\"
    cd '${REMOTE_REPO}'
    job=\$(sbatch --parsable ${dependency_option} --job-name='${job_name}' --export='${exports}' '${wrapper}')
    job=\${job%%;*}
    [[ \${job} =~ ^[0-9]+$ ]]
    printf '%s\\n' \"\${job}\" > '${record}.tmp'
    mv -- '${record}.tmp' '${record}'
    cat '${record}'
  "
}

make_and_submit_admitted_run() {
  local run_kind=$1
  [[ "${run_kind}" == "synthetic" || "${run_kind}" == "real" ]] || { echo "Unknown admitted run kind." >&2; exit 2; }
  if [[ "${run_kind}" == "real" ]]; then
    : "${TWIRL_FM0_TARGET_STEP:?set target optimizer step (1--20000)}"
    [[ "${TWIRL_FM0_TARGET_STEP}" =~ ^[0-9]+$ ]] || { echo "Real target step must be numeric." >&2; exit 2; }
    (( TWIRL_FM0_TARGET_STEP >= 1 && TWIRL_FM0_TARGET_STEP <= 20000 )) || { echo "Invalid real target step." >&2; exit 2; }
  fi
  input_binding; loader_binding; require_socket
  if [[ "${DRY_RUN}" == 1 ]]; then
    echo "+ inspect all ${PARTITION} RUNNING/PENDING jobs and ${GPU_NODE} CfgTRES/AllocTRES"
    echo "+ reserve 4 x (1 H200, 16 CPU, 384 GiB) for Stage 1 before 1 x (1 H200, 4 CPU, 32 GiB) for FM0.1"
    echo "+ publish <=300-second immutable admission receipt and immediately sbatch one non-array H200 ${run_kind} run"
    return
  fi

  capture "
    set -euo pipefail
    ${verify_remote_repo}
    [[ \$(sha256sum '${TWIRL_FM0_INPUT_RECEIPT}' | awk '{print \$1}') == '${TWIRL_FM0_INPUT_RECEIPT_SHA256}' ]]
    [[ \$(sha256sum '${TWIRL_FM0_LOADER_RECEIPT}' | awk '{print \$1}') == '${TWIRL_FM0_LOADER_RECEIPT_SHA256}' ]]
    '${ENV_PREFIX}/bin/python' -c 'import json,sys; i=json.load(open(sys.argv[1])); l=json.load(open(sys.argv[2])); assert i.get(\"schema_version\")==\"twirl_fm0_1_input_release_receipt_v1\" and i.get(\"passed\") is True; assert l.get(\"schema_version\")==\"twirl_fm0_1_loader_smoke_receipt_v1\" and l.get(\"passed\") is True; assert l[\"input_receipt\"][\"sha256\"]==sys.argv[3]' '${TWIRL_FM0_INPUT_RECEIPT}' '${TWIRL_FM0_LOADER_RECEIPT}' '${TWIRL_FM0_INPUT_RECEIPT_SHA256}'
  " >/dev/null

  local tmp remote_epoch local_epoch skew admission_path admission_sha run_output
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/twirl-fm0-admission.XXXXXX")
  ADMISSION_TMP="${tmp}"
  trap 'if [[ -n "${ADMISSION_TMP:-}" ]]; then rm -rf -- "${ADMISSION_TMP}"; fi' EXIT
  remote_epoch=$(capture "date -u +%s")
  local_epoch=$(date -u +%s)
  skew=$(( local_epoch > remote_epoch ? local_epoch - remote_epoch : remote_epoch - local_epoch ))
  (( skew <= 30 )) || { echo "Local/ORCD clock skew exceeds 30 seconds." >&2; exit 5; }
  capture "scontrol show partition -o '${PARTITION}'" > "${tmp}/partition.txt"
  capture "scontrol show partition -o '${RECLAIMABLE_PARTITION}'" > "${tmp}/reclaimable-partition.txt"
  capture "scontrol show config | grep -E '^Preempt(Type|Mode|Parameters|ExemptTime)'" > "${tmp}/preemption.txt"
  capture "scontrol show node -o '${GPU_NODE}'" > "${tmp}/node.txt"
  capture "squeue -p '${PARTITION}' -h -t RUNNING,PENDING,CONFIGURING,COMPLETING -o '%i|%u|%T|%r|%C|%m|%b|%j|%N|%E'" > "${tmp}/queue.txt"
  capture "squeue -w '${GPU_NODE}' -h -t RUNNING,CONFIGURING,COMPLETING -o '%i|%P|%u|%T|%r|%C|%m|%b|%j|%N|%E'" > "${tmp}/node-queue.txt"
  admission_path="${RUN_ROOT}/admission/admission-${run_kind}-${remote_epoch}-${TWIRL_FM0_INPUT_RECEIPT_SHA256:0:12}.json"

  python3 - "${tmp}/partition.txt" "${tmp}/reclaimable-partition.txt" "${tmp}/preemption.txt" \
    "${tmp}/node.txt" "${tmp}/queue.txt" "${tmp}/node-queue.txt" "${tmp}/receipt.json" \
    "${remote_epoch}" "${EXPECTED_SHA}" "${TWIRL_FM0_INPUT_RECEIPT}" "${TWIRL_FM0_INPUT_RECEIPT_SHA256}" <<'PY'
from __future__ import annotations
import hashlib, json, re, sys
from pathlib import Path

partition_text, reclaim_partition_text, preemption_text, node_text, queue_path, node_queue_path, out_path = [Path(x) for x in sys.argv[1:8]]
now, expected_sha, input_path, input_sha = int(sys.argv[8]), sys.argv[9], sys.argv[10], sys.argv[11]
partition_raw = partition_text.read_text().strip()
reclaim_partition_raw = reclaim_partition_text.read_text().strip()
preemption_raw = preemption_text.read_text().strip()
node_raw = node_text.read_text().strip()
if not partition_raw or not reclaim_partition_raw or not preemption_raw or not node_raw:
    raise SystemExit("empty scheduler snapshot")
def field(text, name):
    match = re.search(r"(?:^|\s)" + re.escape(name) + r"=([^\s]+)", text)
    if not match: raise SystemExit(f"missing scheduler field: {name}")
    return match.group(1)
def tres(text, name, *, zero_if_absent=False):
    match = re.search(r"(?:^|,)" + re.escape(name) + r"=([0-9]+)", text)
    if not match:
        if zero_if_absent:
            return 0
        raise SystemExit(f"missing TRES: {name}")
    return int(match.group(1))
def memory_mib(value):
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([KMGTP]?)", value)
    if not match: raise SystemExit(f"unparseable memory: {value}")
    scale = {"": 1/1024**2, "K": 1/1024, "M": 1, "G": 1024, "T": 1024**2, "P": 1024**3}[match.group(2)]
    return int(float(match.group(1)) * scale)
if field(partition_raw, "PartitionName") != "pg_mki_aryeh" or field(partition_raw, "State") != "UP":
    raise SystemExit("partition is not the exact live authority")
if field(reclaim_partition_raw, "PartitionName") != "mit_preemptable" or field(reclaim_partition_raw, "State") != "UP":
    raise SystemExit("reclaimable partition is not the exact live authority")
if int(field(partition_raw, "PriorityTier")) <= int(field(reclaim_partition_raw, "PriorityTier")):
    raise SystemExit("FM partition no longer outranks the reclaimable partition")
if field(reclaim_partition_raw, "PreemptMode") != "REQUEUE":
    raise SystemExit("reclaimable partition no longer uses REQUEUE")
if not re.search(r"^PreemptType\s*=\s*preempt/partition_prio$", preemption_raw, re.MULTILINE):
    raise SystemExit("partition-priority preemption is not active")
if field(node_raw, "NodeName") != "node4900": raise SystemExit("wrong GPU node")
bad_states = {"DOWN", "DRAIN", "DRAINED", "FAIL", "FAILING", "MAINT", "UNKNOWN"}
node_states = set(field(node_raw, "State").split("+"))
if node_states & bad_states: raise SystemExit(f"GPU node is ineligible: {sorted(node_states)}")
cfg_tres, alloc_tres = field(node_raw, "CfgTRES"), field(node_raw, "AllocTRES")
cfg_h200 = tres(cfg_tres, "gres/gpu:h200")
alloc_h200 = tres(alloc_tres, "gres/gpu:h200", zero_if_absent=True)
cfg_cpu, alloc_cpu = tres(cfg_tres, "cpu"), tres(alloc_tres, "cpu")
real_mem, alloc_mem = int(field(node_raw, "RealMemory")), int(field(node_raw, "AllocMem"))
if cfg_h200 != 8: raise SystemExit(f"node4900 H200 count drift: {cfg_h200}")

pattern = re.compile(r"^twirl-s[0-9]+-o[0-9]+-c[1-4]d[1-4]-a2v1(?:-p4)?$")
active, dependency_held, runnable = [], [], []
active_h200 = active_cpu = active_mem = 0
for line in Path(queue_path).read_text().splitlines():
    fields = line.split("|", 9)
    if len(fields) != 10: raise SystemExit(f"unparseable all-partition queue row: {line}")
    jobid, user, state, reason, cpus, mem, gres, name, nodes, dependency = fields
    if name in {"twirl-fm0-fp32-smoke", "twirl-fm0-1-1-real"}:
        raise SystemExit(f"another FM0.1 GPU job is already live: {jobid}/{state}")
    if not pattern.fullmatch(name): continue
    row = {"job_id": jobid, "user": user, "state": state, "reason": reason, "cpus": int(cpus), "memory": mem, "gres": gres, "nodes": nodes, "dependency": dependency}
    if state in {"RUNNING", "CONFIGURING", "COMPLETING"}:
        match = re.search(r"h200:([0-9]+)", gres.lower())
        if not match or "node4900" not in nodes: raise SystemExit(f"active Stage-1 resource ambiguity: {jobid}")
        active.append(row); active_h200 += int(match.group(1)); active_cpu += int(cpus); active_mem += memory_mib(mem)
    elif state == "PENDING" and reason == "Dependency": dependency_held.append(row)
    elif state == "PENDING": runnable.append(row)
    else: raise SystemExit(f"unknown Stage-1 scheduler state: {jobid}/{state}")
if runnable: raise SystemExit("runnable Stage-1 job is pending; FM admission denied")
if active_h200 > 4 or active_cpu > 78: raise SystemExit("Stage-1 live ceiling already exceeded")

reclaimable, reclaim_h200, reclaim_cpu, reclaim_mem = [], 0, 0, 0
for line in node_queue_path.read_text().splitlines():
    fields = line.split("|", 10)
    if len(fields) != 11: raise SystemExit(f"unparseable GPU-node queue row: {line}")
    jobid, partition, user, state, reason, cpus, mem, gres, name, nodes, dependency = fields
    if partition != "mit_preemptable":
        continue
    if state not in {"RUNNING", "CONFIGURING", "COMPLETING"} or "node4900" not in nodes:
        raise SystemExit(f"reclaimable job state/node ambiguity: {jobid}/{state}/{nodes}")
    match = re.search(r"h200:([0-9]+)", gres.lower())
    h200 = int(match.group(1)) if match else 0
    row = {"job_id": jobid, "partition": partition, "user": user, "state": state, "reason": reason, "cpus": int(cpus), "memory": mem, "gres": gres, "name": name, "nodes": nodes, "dependency": dependency}
    reclaimable.append(row)
    reclaim_h200 += h200
    reclaim_cpu += int(cpus)
    reclaim_mem += memory_mib(mem)

other_h200 = alloc_h200-active_h200-reclaim_h200
other_cpu = alloc_cpu-active_cpu-reclaim_cpu
other_mem = alloc_mem-active_mem-reclaim_mem
if min(other_h200, other_cpu, other_mem) < 0:
    raise SystemExit("node/queue allocation accounting is inconsistent after reclaimable jobs")
stage1 = {"lanes": 4, "h200": 4, "cpus": 64, "memory_mib": 1572864}
fm = {"h200": 1, "cpus": 4, "memory_mib": 32768}
needed = {"h200": other_h200+stage1["h200"]+fm["h200"], "cpus": other_cpu+stage1["cpus"]+fm["cpus"], "memory_mib": other_mem+stage1["memory_mib"]+fm["memory_mib"]}
capacity = {"h200": cfg_h200, "cpus": cfg_cpu, "memory_mib": real_mem}
if any(needed[key] > capacity[key] for key in needed):
    raise SystemExit(f"insufficient capacity after Stage-1 reservation: needed={needed}, capacity={capacity}")
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
payload = {
  "schema_version": "twirl_fm0_1_orcd_admission_receipt_v1", "passed": True,
  "created_at_epoch": now, "expires_at_epoch": now+300, "ttl_seconds": 300,
  "expected_git_sha": expected_sha,
  "input_receipt": {"path": input_path, "sha256": input_sha},
  "partition": "pg_mki_aryeh", "gpu_node": "node4900",
  "requested_fm_resources": fm, "stage1_reservation": stage1,
  "derived": {"node_capacity": capacity, "node_allocated": {"h200": alloc_h200, "cpus": alloc_cpu, "memory_mib": alloc_mem}, "reclaimable_allocated": {"h200": reclaim_h200, "cpus": reclaim_cpu, "memory_mib": reclaim_mem}, "non_reclaimable_non_stage1_allocated": {"h200": other_h200, "cpus": other_cpu, "memory_mib": other_mem}, "capacity_needed": needed, "active_stage1_jobs": active, "dependency_held_stage1_jobs": dependency_held, "runnable_stage1_jobs": runnable, "reclaimable_node_jobs": reclaimable},
  "snapshot_sha256": {"partition": digest(partition_text), "reclaimable_partition": digest(reclaim_partition_text), "preemption": digest(preemption_text), "node": digest(node_text), "all_partition_queue": digest(Path(queue_path)), "gpu_node_queue": digest(node_queue_path)},
  "raw_snapshots": {"partition": partition_raw, "reclaimable_partition": reclaim_partition_raw, "preemption": preemption_raw, "node": node_raw, "all_partition_queue": Path(queue_path).read_text(), "gpu_node_queue": node_queue_path.read_text()},
}
Path(out_path).write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
PY
  admission_sha=$(shasum -a 256 "${tmp}/receipt.json" | awk '{print $1}')
  "${SSH[@]}" "${LOCKED_HOST}" "set -euo pipefail; mkdir -p '${RUN_ROOT}/admission'; test ! -e '${admission_path}'; tmp='${admission_path}.tmp'; cat > \"\${tmp}\"; chmod 0444 \"\${tmp}\"; mv -- \"\${tmp}\" '${admission_path}'" < "${tmp}/receipt.json"
  if [[ "${run_kind}" == "synthetic" ]]; then
    run_output="${RUN_ROOT}/model_runs/fp32_synthetic_smoke/${remote_epoch}-${admission_sha:0:12}"
    submit_once "${LAUNCH_DIR}/fp32-smoke.job" "twirl-fm0-fp32-smoke" \
      "scripts/orcd/slurm_twirl_fm0_1_fp32_smoke_h200.sbatch" \
      "$(base_export),TWIRL_FM0_INPUT_RECEIPT=${TWIRL_FM0_INPUT_RECEIPT},TWIRL_FM0_INPUT_RECEIPT_SHA256=${TWIRL_FM0_INPUT_RECEIPT_SHA256},TWIRL_FM0_LOADER_RECEIPT=${TWIRL_FM0_LOADER_RECEIPT},TWIRL_FM0_LOADER_RECEIPT_SHA256=${TWIRL_FM0_LOADER_RECEIPT_SHA256},TWIRL_FM0_ADMISSION_RECEIPT=${admission_path},TWIRL_FM0_ADMISSION_RECEIPT_SHA256=${admission_sha},TWIRL_FM0_SMOKE_OUTPUT=${run_output}"
  else
    run_output="${RUN_ROOT}/model_runs/real_fm0_1_1/${remote_epoch}-${admission_sha:0:12}-step${TWIRL_FM0_TARGET_STEP}"
    submit_once "${LAUNCH_DIR}/real-train-step${TWIRL_FM0_TARGET_STEP}.job" "twirl-fm0-1-1-real" \
      "scripts/orcd/slurm_twirl_fm0_1_real_train_h200.sbatch" \
      "$(base_export),TWIRL_FM0_INPUT_RECEIPT=${TWIRL_FM0_INPUT_RECEIPT},TWIRL_FM0_INPUT_RECEIPT_SHA256=${TWIRL_FM0_INPUT_RECEIPT_SHA256},TWIRL_FM0_LOADER_RECEIPT=${TWIRL_FM0_LOADER_RECEIPT},TWIRL_FM0_LOADER_RECEIPT_SHA256=${TWIRL_FM0_LOADER_RECEIPT_SHA256},TWIRL_FM0_ADMISSION_RECEIPT=${admission_path},TWIRL_FM0_ADMISSION_RECEIPT_SHA256=${admission_sha},TWIRL_FM0_TRAIN_OUTPUT=${run_output},TWIRL_FM0_TARGET_STEP=${TWIRL_FM0_TARGET_STEP}"
  fi
  printf 'admission=%s\nadmission_sha256=%s\nrun_output=%s\n' "${admission_path}" "${admission_sha}" "${run_output}"
  rm -rf -- "${tmp}"
  ADMISSION_TMP=""
  trap - EXIT
}

case "${ACTION}" in
  probe)
    require_socket
    remote "hostname; whoami; sinfo -h -p '${PARTITION}' -o '%P|%D|%t|%G'; scontrol show node -o '${GPU_NODE}'"
    ;;
  deploy)
    require_socket
    remote "set -euo pipefail; test -d '${REMOTE_SOURCE}/.git'; git -C '${REMOTE_SOURCE}' cat-file -e '${EXPECTED_SHA}^{commit}'; if [[ -e '${REMOTE_REPO}' ]]; then [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]; else git -C '${REMOTE_SOURCE}' worktree add --detach '${REMOTE_REPO}' '${EXPECTED_SHA}'; fi; if git -C '${REMOTE_REPO}' sparse-checkout list >/dev/null 2>&1; then git -C '${REMOTE_REPO}' sparse-checkout add reports/stage5_validation/twirl_fm0_1_design_freeze_v1; fi; test -s '${REMOTE_REPO}/reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json'; [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]; git -C '${REMOTE_REPO}' rev-parse HEAD"
    ;;
  build-env)
    require_socket
    remote "set -euo pipefail; [[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]; [[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]; cd '${REMOTE_REPO}'; TWIRL_ORCD_REPO='${REMOTE_REPO}' TWIRL_EXPECTED_GIT_SHA='${EXPECTED_SHA}' bash scripts/orcd/build_orcd_fm0_env.sh"
    ;;
  submit-sector-stage)
    : "${TWIRL_FM0_SECTOR:?set one sector in S56-S65}"
    : "${TWIRL_FM0_CORPUS_SELECTION:?set the immutable corpus selection}"
    : "${TWIRL_FM0_CORPUS_SELECTION_SHA256:?set its SHA-256}"
    : "${TWIRL_FM0_SECTOR_ARCHIVE_DIR:?set the verified sector archive directory}"
    : "${TWIRL_FM0_QUALITY_TABLE:?set the cadence-quality table}"
    : "${TWIRL_FM0_QUALITY_MANIFEST:?set the cadence-quality manifest}"
    : "${TWIRL_FM0_SECTOR_STAGE_ROOT:?set the immutable sector-stage root}"
    [[ "${TWIRL_FM0_SECTOR}" =~ ^(5[6-9]|6[0-5])$ ]] || { echo "Sector must be S56-S65." >&2; exit 2; }
    require_hash 64 "${TWIRL_FM0_CORPUS_SELECTION_SHA256}"
    for path in "${TWIRL_FM0_CORPUS_SELECTION}" "${TWIRL_FM0_SECTOR_ARCHIVE_DIR}" "${TWIRL_FM0_QUALITY_TABLE}" "${TWIRL_FM0_QUALITY_MANIFEST}" "${TWIRL_FM0_SECTOR_STAGE_ROOT}"; do
      safe_remote_path "${path}" || { echo "Unsafe sector-stage path: ${path}" >&2; exit 2; }
    done
    require_socket
    submit_once "${LAUNCH_DIR}/sector-stage-s${TWIRL_FM0_SECTOR}.job" "twirl-fm0-stage-s${TWIRL_FM0_SECTOR}" "scripts/orcd/slurm_twirl_fm0_1_sector_stage_cpu.sbatch" "$(base_export),TWIRL_FM0_SECTOR=${TWIRL_FM0_SECTOR},TWIRL_FM0_CORPUS_SELECTION=${TWIRL_FM0_CORPUS_SELECTION},TWIRL_FM0_CORPUS_SELECTION_SHA256=${TWIRL_FM0_CORPUS_SELECTION_SHA256},TWIRL_FM0_SECTOR_ARCHIVE_DIR=${TWIRL_FM0_SECTOR_ARCHIVE_DIR},TWIRL_FM0_QUALITY_TABLE=${TWIRL_FM0_QUALITY_TABLE},TWIRL_FM0_QUALITY_MANIFEST=${TWIRL_FM0_QUALITY_MANIFEST},TWIRL_FM0_SECTOR_STAGE_ROOT=${TWIRL_FM0_SECTOR_STAGE_ROOT}" "${TWIRL_FM0_AFTEROK_JOB_ID:-}"
    ;;
  submit-registry)
    : "${TWIRL_FM0_ALIASES_TABLE:?set staged aliases table}"
    source_binding
    safe_remote_path "${TWIRL_FM0_ALIASES_TABLE}" || { echo "Unsafe staged alias path." >&2; exit 2; }
    require_socket
    if [[ "${TWIRL_FM0_SOURCE_ADAPTER}" == "strict_npz_fixture_v1" ]]; then
      source_exports="TWIRL_FM0_OBSERVATIONS_TABLE=${TWIRL_FM0_OBSERVATIONS_TABLE},TWIRL_FM0_SOURCE_MANIFEST=${TWIRL_FM0_SOURCE_MANIFEST}"
    else
      source_exports="TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY=${TWIRL_FM0_A2V1_HDF5_SOURCE_INVENTORY}"
    fi
    submit_once "${LAUNCH_DIR}/registry.job" "twirl-fm0-registry" "scripts/orcd/slurm_twirl_fm0_1_registry_cpu.sbatch" "$(base_export),TWIRL_FM0_ALIASES_TABLE=${TWIRL_FM0_ALIASES_TABLE},TWIRL_FM0_SOURCE_ADAPTER=${TWIRL_FM0_SOURCE_ADAPTER},${source_exports}"
    ;;
  submit-input-validation)
    source_binding
    require_socket
    if [[ "${TWIRL_FM0_SOURCE_ADAPTER}" == "strict_npz_fixture_v1" ]]; then
      source_manifest="${TWIRL_FM0_SOURCE_MANIFEST}"
    else
      source_manifest="${RUN_ROOT}/input_build/registry/a2v1_hdf5_manifest.csv"
    fi
    submit_once "${LAUNCH_DIR}/input-validation.job" "twirl-fm0-input" "scripts/orcd/slurm_twirl_fm0_1_input_validation_cpu.sbatch" "$(base_export),TWIRL_FM0_SOURCE_MANIFEST=${source_manifest},TWIRL_FM0_SOURCE_ADAPTER=${TWIRL_FM0_SOURCE_ADAPTER}"
    ;;
  submit-loader-smoke)
    input_binding; require_socket
    submit_once "${LAUNCH_DIR}/loader-smoke.job" "twirl-fm0-loader" "scripts/orcd/slurm_twirl_fm0_1_loader_smoke_cpu.sbatch" "$(base_export),TWIRL_FM0_INPUT_RECEIPT=${TWIRL_FM0_INPUT_RECEIPT},TWIRL_FM0_INPUT_RECEIPT_SHA256=${TWIRL_FM0_INPUT_RECEIPT_SHA256}"
    ;;
  submit-fp32-smoke) make_and_submit_admitted_run synthetic ;;
  submit-real-train) make_and_submit_admitted_run real ;;
  submit-post-validation)
    : "${TWIRL_FM0_SMOKE_OUTPUT:?set exact smoke output}"; : "${TWIRL_FM0_ADMISSION_RECEIPT:?set admission receipt}"; : "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?set admission hash}"
    input_binding; safe_remote_path "${TWIRL_FM0_SMOKE_OUTPUT}" && safe_remote_path "${TWIRL_FM0_ADMISSION_RECEIPT}" || { echo "Unsafe smoke/admission path." >&2; exit 2; }; require_hash 64 "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"
    require_socket
    submit_once "${LAUNCH_DIR}/post-validation.job" "twirl-fm0-postvalidate" "scripts/orcd/slurm_twirl_fm0_1_post_validation_cpu.sbatch" "$(base_export),TWIRL_FM0_SMOKE_OUTPUT=${TWIRL_FM0_SMOKE_OUTPUT},TWIRL_FM0_INPUT_RECEIPT=${TWIRL_FM0_INPUT_RECEIPT},TWIRL_FM0_INPUT_RECEIPT_SHA256=${TWIRL_FM0_INPUT_RECEIPT_SHA256},TWIRL_FM0_ADMISSION_RECEIPT=${TWIRL_FM0_ADMISSION_RECEIPT},TWIRL_FM0_ADMISSION_RECEIPT_SHA256=${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"
    ;;
  submit-real-post-validation)
    : "${TWIRL_FM0_TRAIN_OUTPUT:?set exact real training output}"
    : "${TWIRL_FM0_RUN_GIT_SHA:?set exact training commit SHA}"
    : "${TWIRL_FM0_ADMISSION_RECEIPT:?set admission receipt}"
    : "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256:?set admission hash}"
    input_binding
    safe_remote_path "${TWIRL_FM0_TRAIN_OUTPUT}" && safe_remote_path "${TWIRL_FM0_ADMISSION_RECEIPT}" || { echo "Unsafe training/admission path." >&2; exit 2; }
    require_hash 64 "${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"
    require_hash 40 "${TWIRL_FM0_RUN_GIT_SHA}"
    require_socket
    real_validation_name=$(basename "${TWIRL_FM0_TRAIN_OUTPUT}")
    [[ "${real_validation_name}" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Unsafe real validation name." >&2; exit 2; }
    submit_once "${LAUNCH_DIR}/real-post-validation-${TWIRL_FM0_RUN_GIT_SHA:0:12}-${real_validation_name}.job" "twirl-fm0-real-validate" "scripts/orcd/slurm_twirl_fm0_1_real_post_validation_cpu.sbatch" "$(base_export),TWIRL_FM0_RUN_GIT_SHA=${TWIRL_FM0_RUN_GIT_SHA},TWIRL_FM0_TRAIN_OUTPUT=${TWIRL_FM0_TRAIN_OUTPUT},TWIRL_FM0_INPUT_RECEIPT=${TWIRL_FM0_INPUT_RECEIPT},TWIRL_FM0_INPUT_RECEIPT_SHA256=${TWIRL_FM0_INPUT_RECEIPT_SHA256},TWIRL_FM0_ADMISSION_RECEIPT=${TWIRL_FM0_ADMISSION_RECEIPT},TWIRL_FM0_ADMISSION_RECEIPT_SHA256=${TWIRL_FM0_ADMISSION_RECEIPT_SHA256}"
    ;;
  status)
    require_socket
    remote "set -euo pipefail; squeue -p '${PARTITION}' -h -o '%i|%u|%T|%r|%C|%b|%j|%N'; find '${RUN_ROOT}' -maxdepth 4 -type f \( -name '*.job' -o -name '*.receipt.json' \) -print 2>/dev/null | sort || true"
    ;;
  *) usage >&2; exit 2 ;;
esac
