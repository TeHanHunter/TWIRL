#!/usr/bin/env bash
# Guarded local controller for CPU-only TWIRL-FM S66-S77 full-visit preparation.
#
# The controller reuses only the user-opened ORCD SSH control socket.  It has
# no password, Duo, keyboard-interactive, public-key, hostbased, or GSSAPI
# fallback.  Dry-run is the default, and no action can touch S65 or S78-S80,
# submit an H200 job, freeze a temporal panel, open sealed aperture photometry,
# write a derived sealed shard, or start model training.
set -euo pipefail

DRY_RUN=1
ACTION=""

usage() {
  cat <<'EOF'
Usage: TWIRL_EXPECTED_GIT_SHA=<40-hex> run_twirl_fm0_s66_s77_preparation_orcd.sh [--run] ACTION

Actions:
  probe                  Verify the already-open noninteractive ORCD socket.
  deploy                 Create/reuse one exact clean detached ORCD worktree.
  submit-phase-a         Submit independent chronological S66-S77 quality and
                         source-inventory afterok chains, with durable records.
  submit-phase-a-freeze  Submit the CPU full-validator after both chain tails.
  submit-six-view        Require the completed immutable Phase-A authority
                         record, then submit four chronological sector lanes.
  submit-admission       Require the completed six-view chain, then submit the
                         final CPU admission-v2/full-bundle validation job.
  retry-admission        After one recorded failed admission with no receipt,
                         submit one new immutable repair attempt.
  validate-admission     Read-only terminal job/receipt validation.
  status                 Read-only job-record, Slurm, and artifact status.

Default mode is dry-run. Mutating actions execute only with --run. This
controller never retries, cancels, overwrites, weakens, trains, uses H200,
or accesses a sealed-test light curve/shard.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run)
      DRY_RUN=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -n "${ACTION}" ]]; then
        echo "Only one action may be supplied." >&2
        exit 2
      fi
      ACTION="$1"
      shift
      ;;
  esac
done
if [[ -z "${ACTION}" ]]; then
  usage >&2
  exit 2
fi
case "${ACTION}" in
  probe|deploy|submit-phase-a|submit-phase-a-freeze|submit-six-view|submit-admission|retry-admission|validate-admission|status) ;;
  *) echo "Unknown action: ${ACTION}" >&2; usage >&2; exit 2 ;;
esac

: "${TWIRL_EXPECTED_GIT_SHA:?set TWIRL_EXPECTED_GIT_SHA to the reviewed 40-hex commit}"
readonly EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA}"
if [[ ! "${EXPECTED_SHA}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "TWIRL_EXPECTED_GIT_SHA must be a full lowercase Git SHA." >&2
  exit 2
fi

readonly LOCAL_REPO="${LOCAL_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
readonly MAP_REL="configs/orcd/twirl_fm0_s66_s77_full_visit_preparation_v1.json"
readonly LOCAL_MAP="${LOCAL_REPO}/${MAP_REL}"
[[ -s "${LOCAL_MAP}" ]] || { echo "Missing authority map: ${LOCAL_MAP}" >&2; exit 2; }
readonly MAP_SHA256="$(shasum -a 256 "${LOCAL_MAP}" | awk '{print $1}')"

# The map is trusted only after its exact Git binding is checked in --run mode.
# This parser emits shell-quoted values and rejects any sector/resource drift.
eval "$(python3 - "${LOCAL_MAP}" <<'PY'
from __future__ import annotations

import json
import shlex
import sys
from pathlib import Path

path = Path(sys.argv[1])
data = json.loads(path.read_text(encoding="utf-8"))
scope = data["scope"]
sectors = list(range(66, 78))
assert data["schema_version"] == "twirl_fm0_s66_s77_full_visit_preparation_authority_v1"
assert data["campaign_id"] == "twirl_fm0_s66_s77_full_visit_preparation_v1"
assert scope["ordered_sectors"] == sectors
assert scope["excluded_sectors"] == [65]
assert scope["blocked_sectors"] == [78]
assert scope["untouched_sectors"] == [79, 80]
assert scope["label_blind_sector_wide_hdf5_quality_qa_allowed"] is True
assert scope["sealed_aperture_photometry_or_derived_shard_access_forbidden"] is True
assert scope["model_training_authorized"] is False
assert scope["h200_use_authorized"] is False
orcd = data["orcd"]
assert orcd["host"] == "tehan@orcd-login.mit.edu"
assert orcd["control_socket"] == "/Users/tehan/.ssh/cm/tehan@orcd-login.mit.edu:22"
assert orcd["partition"] == "pg_mki_aryeh"
assert orcd["excluded_node"] == "node4900"
authorities = data["authorities"]
assert sorted(map(int, authorities["source_receipts"])) == sectors
assert sorted(map(int, authorities["hdf5_quality_receipts"])) == sectors
for sector in sectors:
    for name in ("source_receipts", "hdf5_quality_receipts"):
        row = authorities[name][str(sector)]
        assert len(row["sha256"]) == 64 and set(row["sha256"]) <= set("0123456789abcdef")

values = {
    "ORCD_HOST": orcd["host"],
    "CONTROL_PATH": orcd["control_socket"],
    "PROJECT_ROOT": orcd["project_root"],
    "REMOTE_SOURCE": orcd["source_repository"],
    "DETACHED_PREFIX": orcd["detached_repository_prefix"],
    "REMOTE_PYTHON": orcd["python"],
    "QUALITY_TRANSFER_ROOT": authorities["mission_quality_transfer"]["root"],
    "QUALITY_TRANSFER_MANIFEST": authorities["mission_quality_transfer"]["manifest_path"],
    "QUALITY_TRANSFER_SHA256": authorities["mission_quality_transfer"]["manifest_sha256"],
    "CORPUS_SELECTION": authorities["corpus_selection"]["path"],
    "CORPUS_SELECTION_SHA256": authorities["corpus_selection"]["sha256"],
    "CORPUS_SUMMARY": authorities["corpus_selection"]["summary_path"],
    "CORPUS_SUMMARY_SHA256": authorities["corpus_selection"]["summary_sha256"],
    "TARGET_AUTHORITY": authorities["target_observations"]["path"],
    "TARGET_AUTHORITY_SHA256": authorities["target_observations"]["sha256"],
    "ALIAS_AUTHORITY": authorities["gaia_tic_aliases"]["path"],
    "ALIAS_AUTHORITY_SHA256": authorities["gaia_tic_aliases"]["sha256"],
}
values.update({name.upper(): value for name, value in data["outputs"].items()})
for sector in sectors:
    values[f"SOURCE_PATH_{sector}"] = authorities["source_receipts"][str(sector)]["path"]
    values[f"SOURCE_SHA_{sector}"] = authorities["source_receipts"][str(sector)]["sha256"]
    values[f"HDF5_PATH_{sector}"] = authorities["hdf5_quality_receipts"][str(sector)]["path"]
    values[f"HDF5_SHA_{sector}"] = authorities["hdf5_quality_receipts"][str(sector)]["sha256"]
for name, value in values.items():
    print(f"{name}={shlex.quote(str(value))}")
PY
)"

readonly ORCD_HOST CONTROL_PATH PROJECT_ROOT REMOTE_SOURCE DETACHED_PREFIX REMOTE_PYTHON
readonly QUALITY_TRANSFER_ROOT QUALITY_TRANSFER_MANIFEST QUALITY_TRANSFER_SHA256
readonly CORPUS_SELECTION CORPUS_SELECTION_SHA256 CORPUS_SUMMARY CORPUS_SUMMARY_SHA256
readonly TARGET_AUTHORITY TARGET_AUTHORITY_SHA256 ALIAS_AUTHORITY ALIAS_AUTHORITY_SHA256
readonly RUN_ROOT QUALITY_REFERENCE_ROOT SOURCE_INVENTORY_ROOT PHASE_A_AUTHORITY_RECORD
readonly SIX_VIEW_ROOT ADMISSION_RECEIPT LAUNCH_ROOT
readonly REMOTE_REPO="${DETACHED_PREFIX}${EXPECTED_SHA:0:12}"
readonly REMOTE_MAP="${REMOTE_REPO}/${MAP_REL}"

SSH=(
  ssh
  -S "${CONTROL_PATH}"
  -o ControlMaster=no
  -o BatchMode=yes
  -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no
  -o PubkeyAuthentication=no
  -o HostbasedAuthentication=no
  -o GSSAPIAuthentication=no
  -o NumberOfPasswordPrompts=0
  -o ConnectTimeout=15
  -o ConnectionAttempts=1
)

print_cmd() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
}

remote() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    print_cmd "${SSH[@]}" "${ORCD_HOST}" "$1"
  else
    "${SSH[@]}" "${ORCD_HOST}" "$1"
  fi
}

require_socket() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    return
  fi
  if ! "${SSH[@]}" -O check "${ORCD_HOST}" >/dev/null 2>&1; then
    cat >&2 <<EOF
The required user-opened ORCD control socket is unavailable:
  ${CONTROL_PATH}
Open it from your own terminal. This controller has no interactive fallback.
EOF
    exit 4
  fi
}

require_local_exact_commit() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    return
  fi
  [[ "$(git -C "${LOCAL_REPO}" rev-parse HEAD)" == "${EXPECTED_SHA}" ]] || {
    echo "Local HEAD is not TWIRL_EXPECTED_GIT_SHA." >&2; exit 3;
  }
  [[ -z "$(git -C "${LOCAL_REPO}" status --porcelain=v1 --untracked-files=all)" ]] || {
    echo "Local repository must be clean before ORCD execution." >&2; exit 3;
  }
  local committed_map_sha
  committed_map_sha="$(git -C "${LOCAL_REPO}" show "${EXPECTED_SHA}:${MAP_REL}" | shasum -a 256 | awk '{print $1}')"
  [[ "${committed_map_sha}" == "${MAP_SHA256}" ]] || {
    echo "Authority map is not byte-exact at the expected commit." >&2; exit 3;
  }
}

remote_gate() {
  cat <<EOF
[[ \$(git -C '${REMOTE_REPO}' rev-parse HEAD) == '${EXPECTED_SHA}' ]]
! git -C '${REMOTE_REPO}' symbolic-ref -q HEAD >/dev/null 2>&1
[[ -z \$(git -C '${REMOTE_REPO}' status --porcelain=v1 --untracked-files=all) ]]
test -s '${REMOTE_MAP}'
[[ \$(sha256sum '${REMOTE_MAP}' | awk '{print \$1}') == '${MAP_SHA256}' ]]
EOF
}

source_case() {
  local sector path_name sha_name
  for sector in $(seq 66 77); do
    path_name="SOURCE_PATH_${sector}"
    sha_name="SOURCE_SHA_${sector}"
    printf "%s) source_receipt='%s'; source_sha='%s';;\n" \
      "${sector}" "${!path_name}" "${!sha_name}"
  done
}

hdf5_case() {
  local sector path_name sha_name
  for sector in $(seq 66 77); do
    path_name="HDF5_PATH_${sector}"
    sha_name="HDF5_SHA_${sector}"
    printf "%s) hdf5_receipt='%s'; hdf5_sha='%s';;\n" \
      "${sector}" "${!path_name}" "${!sha_name}"
  done
}

readonly REMOTE_GATE="$(remote_gate)"
readonly SOURCE_CASE="$(source_case)"
readonly HDF5_CASE="$(hdf5_case)"

require_socket
require_local_exact_commit

case "${ACTION}" in
  probe)
    remote "set -euo pipefail
      hostname
      whoami
      sinfo -h -p pg_mki_aryeh -o '%P %D %t %G' | head -10"
    ;;

  deploy)
    remote "set -euo pipefail
      if [[ -e '${REMOTE_REPO}' ]]; then
        test -e '${REMOTE_REPO}/.git'
      else
        test -d '${REMOTE_SOURCE}/.git'
        git -C '${REMOTE_SOURCE}' cat-file -e '${EXPECTED_SHA}^{commit}'
        git -C '${REMOTE_SOURCE}' worktree add --detach '${REMOTE_REPO}' '${EXPECTED_SHA}'
      fi
      ${REMOTE_GATE}
      git -C '${REMOTE_REPO}' rev-parse HEAD"
    ;;

  submit-phase-a)
    remote "set -euo pipefail
      ${REMOTE_GATE}
      test -x '${REMOTE_PYTHON}'
      for binding in \
        '${QUALITY_TRANSFER_MANIFEST}=${QUALITY_TRANSFER_SHA256}' \
        '${CORPUS_SELECTION}=${CORPUS_SELECTION_SHA256}' \
        '${CORPUS_SUMMARY}=${CORPUS_SUMMARY_SHA256}' \
        '${TARGET_AUTHORITY}=${TARGET_AUTHORITY_SHA256}' \
        '${ALIAS_AUTHORITY}=${ALIAS_AUTHORITY_SHA256}'; do
        path=\"\${binding%=*}\"; expected=\"\${binding#*=}\"
        test -s \"\${path}\"
        [[ \$(sha256sum \"\${path}\" | awk '{print \$1}') == \"\${expected}\" ]]
      done
      for sector in \$(seq 66 77); do
        case \"\${sector}\" in ${SOURCE_CASE} esac
        test -s \"\${source_receipt}\"
        [[ \$(sha256sum \"\${source_receipt}\" | awk '{print \$1}') == \"\${source_sha}\" ]]
        test ! -e '${QUALITY_REFERENCE_ROOT}/s'\$(printf '%04d' \"\${sector}\")
        test ! -e '${SOURCE_INVENTORY_ROOT}/s'\$(printf '%04d' \"\${sector}\")
      done
      mkdir -p '${LAUNCH_ROOT}' '${RUN_ROOT}/.claims'
      quality_record='${LAUNCH_ROOT}/phase_a_quality.tsv'
      source_record='${LAUNCH_ROOT}/phase_a_source.tsv'
      test ! -e \"\${quality_record}\" && test ! -e \"\${quality_record}.submitting\"
      test ! -e \"\${source_record}\" && test ! -e \"\${source_record}.submitting\"
      mkdir -- '${RUN_ROOT}/.claims/submit-phase-a'

      printf 'sector\tjob_id\tdependency\texpected_git_sha\tauthority_map_sha256\n' > \"\${quality_record}.submitting\"
      previous=''
      for sector in \$(seq 66 77); do
        export_spec='ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_SECTOR='\"\${sector}\"',TWIRL_FM0_QUALITY_TRANSFER_ROOT=${QUALITY_TRANSFER_ROOT},TWIRL_FM0_QUALITY_TRANSFER_MANIFEST_SHA256=${QUALITY_TRANSFER_SHA256},TWIRL_FM0_QUALITY_REFERENCE_ROOT=${QUALITY_REFERENCE_ROOT}'
        if [[ -n \"\${previous}\" ]]; then
          job=\$(sbatch --parsable --dependency=afterok:\"\${previous}\" --export=\"\${export_spec}\" '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_mission_quality_reference_cpu.sbatch')
          dependency=\"afterok:\${previous}\"
        else
          job=\$(sbatch --parsable --export=\"\${export_spec}\" '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_mission_quality_reference_cpu.sbatch')
          dependency=''
        fi
        [[ \"\${job}\" =~ ^[0-9]+$ ]]
        printf '%s\t%s\t%s\t%s\t%s\n' \"\${sector}\" \"\${job}\" \"\${dependency}\" '${EXPECTED_SHA}' '${MAP_SHA256}' >> \"\${quality_record}.submitting\"
        previous=\"\${job}\"
      done
      chmod a-w \"\${quality_record}.submitting\"
      mv -- \"\${quality_record}.submitting\" \"\${quality_record}\"

      printf 'sector\tjob_id\tdependency\texpected_git_sha\tauthority_map_sha256\n' > \"\${source_record}.submitting\"
      previous=''
      for sector in \$(seq 66 77); do
        case \"\${sector}\" in ${SOURCE_CASE} esac
        export_spec='ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_SECTOR='\"\${sector}\"',TWIRL_FM0_SOURCE_RECEIPT='\"\${source_receipt}\"',TWIRL_FM0_SOURCE_RECEIPT_SHA256='\"\${source_sha}\"',TWIRL_FM0_CORPUS_SELECTION=${CORPUS_SELECTION},TWIRL_FM0_CORPUS_SELECTION_SHA256=${CORPUS_SELECTION_SHA256},TWIRL_FM0_CORPUS_SUMMARY=${CORPUS_SUMMARY},TWIRL_FM0_CORPUS_SUMMARY_SHA256=${CORPUS_SUMMARY_SHA256},TWIRL_FM0_GAIA_TIC_ALIAS_AUTHORITY=${ALIAS_AUTHORITY},TWIRL_FM0_GAIA_TIC_ALIAS_AUTHORITY_SHA256=${ALIAS_AUTHORITY_SHA256},TWIRL_FM0_TARGET_AUTHORITY_SHA256=${TARGET_AUTHORITY_SHA256},TWIRL_FM0_SOURCE_INVENTORY_ROOT=${SOURCE_INVENTORY_ROOT}'
        if [[ -n \"\${previous}\" ]]; then
          job=\$(sbatch --parsable --dependency=afterok:\"\${previous}\" --export=\"\${export_spec}\" '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_later_source_inventory_cpu.sbatch')
          dependency=\"afterok:\${previous}\"
        else
          job=\$(sbatch --parsable --export=\"\${export_spec}\" '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_later_source_inventory_cpu.sbatch')
          dependency=''
        fi
        [[ \"\${job}\" =~ ^[0-9]+$ ]]
        printf '%s\t%s\t%s\t%s\t%s\n' \"\${sector}\" \"\${job}\" \"\${dependency}\" '${EXPECTED_SHA}' '${MAP_SHA256}' >> \"\${source_record}.submitting\"
        previous=\"\${job}\"
      done
      chmod a-w \"\${source_record}.submitting\"
      mv -- \"\${source_record}.submitting\" \"\${source_record}\"
      cat \"\${quality_record}\" \"\${source_record}\""
    ;;

  submit-phase-a-freeze)
    remote "set -euo pipefail
      ${REMOTE_GATE}
      quality_record='${LAUNCH_ROOT}/phase_a_quality.tsv'
      source_record='${LAUNCH_ROOT}/phase_a_source.tsv'
      [[ \$(wc -l < \"\${quality_record}\") -eq 13 ]]
      [[ \$(wc -l < \"\${source_record}\") -eq 13 ]]
      quality_tail=\$(tail -1 \"\${quality_record}\" | cut -f2)
      source_tail=\$(tail -1 \"\${source_record}\" | cut -f2)
      [[ \"\${quality_tail}\" =~ ^[0-9]+$ && \"\${source_tail}\" =~ ^[0-9]+$ ]]
      test ! -e '${PHASE_A_AUTHORITY_RECORD}'
      record='${LAUNCH_ROOT}/phase_a_freeze.tsv'
      test ! -e \"\${record}\" && test ! -e \"\${record}.submitting\"
      mkdir -- '${RUN_ROOT}/.claims/submit-phase-a-freeze'
      dependency=afterok:\${quality_tail}:\${source_tail}
      printf 'job_id\tdependency\texpected_git_sha\tauthority_map_sha256\n' > \"\${record}.submitting\"
      job=\$(sbatch --parsable --dependency=\"\${dependency}\" \
        --export='ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256=${MAP_SHA256},TWIRL_FM0_PHASE_A_AUTHORITY_RECORD=${PHASE_A_AUTHORITY_RECORD}' \
        '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_phase_a_authority_freeze_cpu.sbatch')
      [[ \"\${job}\" =~ ^[0-9]+$ ]]
      printf '%s\t%s\t%s\t%s\n' \
        \"\${job}\" \"\${dependency}\" '${EXPECTED_SHA}' '${MAP_SHA256}' >> \"\${record}.submitting\"
      chmod a-w \"\${record}.submitting\"
      mv -- \"\${record}.submitting\" \"\${record}\"
      cat \"\${record}\""
    ;;

  submit-six-view)
    remote "set -euo pipefail
      ${REMOTE_GATE}
      freeze_record='${LAUNCH_ROOT}/phase_a_freeze.tsv'
      [[ \$(wc -l < \"\${freeze_record}\") -eq 2 ]]
      freeze_job=\$(tail -1 \"\${freeze_record}\" | cut -f1)
      phase_a_producer=\$(tail -1 \"\${freeze_record}\" | cut -f3)
      phase_a_map=\$(tail -1 \"\${freeze_record}\" | cut -f4)
      [[ \"\${phase_a_producer}\" =~ ^[0-9a-f]{40}$ ]]
      [[ \"\${phase_a_map}\" == '${MAP_SHA256}' ]]
      state=\$(sacct -nX -j \"\${freeze_job}\" -o State,ExitCode | awk 'NF {print \$1,\$2; exit}')
      [[ \"\${state}\" == 'COMPLETED 0:0' ]]
      test -s '${PHASE_A_AUTHORITY_RECORD}' && test ! -w '${PHASE_A_AUTHORITY_RECORD}'
      phase_sha=\$(sha256sum '${PHASE_A_AUTHORITY_RECORD}' | awk '{print \$1}')
      cd '${REMOTE_REPO}'
      PYTHONPATH='${REMOTE_REPO}/src' '${REMOTE_PYTHON}' scripts/stage5_validation/manage_twirl_fm0_s66_s77_preparation_authorities.py \
        validate-phase-a --authority-map '${REMOTE_MAP}' --authority-map-sha256 '${MAP_SHA256}' --producer-git-sha \"\${phase_a_producer}\" --phase-a-record-sha256 \"\${phase_sha}\" >/dev/null
      for sector in \$(seq 66 77); do
        test ! -e '${SIX_VIEW_ROOT}/s'\$(printf '%04d' \"\${sector}\")
      done
      record='${LAUNCH_ROOT}/six_view.tsv'
      test ! -e \"\${record}\" && test ! -e \"\${record}.submitting\"
      mkdir -- '${RUN_ROOT}/.claims/submit-six-view'
      printf 'sector\tjob_id\tdependency\texpected_git_sha\tauthority_map_sha256\tphase_a_record_sha256\n' > \"\${record}.submitting\"
      previous_0=''
      previous_1=''
      previous_2=''
      previous_3=''
      for sector in \$(seq 66 77); do
        case \"\${sector}\" in ${HDF5_CASE} esac
        quality_dir='${QUALITY_REFERENCE_ROOT}/s'\$(printf '%04d' \"\${sector}\")
        source_dir='${SOURCE_INVENTORY_ROOT}/s'\$(printf '%04d' \"\${sector}\")
        six_dir='${SIX_VIEW_ROOT}/s'\$(printf '%04d' \"\${sector}\")
        quality_sha=\$(sha256sum \"\${quality_dir}/manifest.json\" | awk '{print \$1}')
        source_sha=\$(sha256sum \"\${source_dir}/summary.json\" | awk '{print \$1}')
        export_spec='ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_SECTOR='\"\${sector}\"',TWIRL_FM0_SOURCE_INVENTORY_DIR='\"\${source_dir}\"',TWIRL_FM0_SOURCE_INVENTORY_SUMMARY_SHA256='\"\${source_sha}\"',TWIRL_FM0_QUALITY_REFERENCE_DIR='\"\${quality_dir}\"',TWIRL_FM0_QUALITY_REFERENCE_MANIFEST_SHA256='\"\${quality_sha}\"',TWIRL_FM0_HDF5_QUALITY_RECEIPT='\"\${hdf5_receipt}\"',TWIRL_FM0_HDF5_QUALITY_RECEIPT_SHA256='\"\${hdf5_sha}\"',TWIRL_FM0_SIX_VIEW_RELEASE_DIR='\"\${six_dir}\"''
        lane=\$(( (sector - 66) % 4 ))
        case \"\${lane}\" in
          0) previous=\"\${previous_0}\" ;;
          1) previous=\"\${previous_1}\" ;;
          2) previous=\"\${previous_2}\" ;;
          3) previous=\"\${previous_3}\" ;;
        esac
        if [[ -n \"\${previous}\" ]]; then
          job=\$(sbatch --parsable --dependency=afterok:\"\${previous}\" --export=\"\${export_spec}\" '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_later_six_view_cpu.sbatch')
          dependency=\"afterok:\${previous}\"
        else
          job=\$(sbatch --parsable --export=\"\${export_spec}\" '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_later_six_view_cpu.sbatch')
          dependency=''
        fi
        [[ \"\${job}\" =~ ^[0-9]+$ ]]
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \"\${sector}\" \"\${job}\" \"\${dependency}\" '${EXPECTED_SHA}' '${MAP_SHA256}' \"\${phase_sha}\" >> \"\${record}.submitting\"
        case \"\${lane}\" in
          0) previous_0=\"\${job}\" ;;
          1) previous_1=\"\${job}\" ;;
          2) previous_2=\"\${job}\" ;;
          3) previous_3=\"\${job}\" ;;
        esac
      done
      chmod a-w \"\${record}.submitting\"
      mv -- \"\${record}.submitting\" \"\${record}\"
      cat \"\${record}\""
    ;;

  submit-admission)
    remote "set -euo pipefail
      ${REMOTE_GATE}
      six_record='${LAUNCH_ROOT}/six_view_fast.tsv'
      if [[ ! -e \"\${six_record}\" ]]; then six_record='${LAUNCH_ROOT}/six_view.tsv'; fi
      [[ \$(wc -l < \"\${six_record}\") -eq 13 ]]
      while IFS= read -r six_job; do
        [[ \"\${six_job}\" =~ ^[0-9]+$ ]]
        state=\$(sacct -nX -j \"\${six_job}\" -o State,ExitCode | awk 'NF {print \$1,\$2; exit}')
        [[ \"\${state}\" == 'COMPLETED 0:0' ]]
      done < <(tail -n +2 \"\${six_record}\" | cut -f2)
      six_tail=\$(tail -1 \"\${six_record}\" | cut -f2)
      for sector in \$(seq 66 77); do
        receipt='${SIX_VIEW_ROOT}/s'\$(printf '%04d' \"\${sector}\")'/receipt.json'
        test -s \"\${receipt}\" && test ! -w \"\${receipt}\"
      done
      phase_record='${LAUNCH_ROOT}/phase_a_freeze.tsv'
      [[ \$(wc -l < \"\${phase_record}\") -eq 2 ]]
      phase_a_producer=\$(tail -1 \"\${phase_record}\" | cut -f3)
      phase_a_map=\$(tail -1 \"\${phase_record}\" | cut -f4)
      [[ \"\${phase_a_producer}\" =~ ^[0-9a-f]{40}$ ]]
      [[ \"\${phase_a_map}\" == '${MAP_SHA256}' ]]
      phase_a_record_sha=\$(tail -n +2 \"\${six_record}\" | cut -f6 | sort -u)
      [[ \$(printf '%s\n' \"\${phase_a_record_sha}\" | wc -l) -eq 1 ]]
      [[ \"\${phase_a_record_sha}\" =~ ^[0-9a-f]{64}$ ]]
      [[ \$(sha256sum '${PHASE_A_AUTHORITY_RECORD}' | awk '{print \$1}') == \"\${phase_a_record_sha}\" ]]
      test ! -e '${ADMISSION_RECEIPT}'
      record='${LAUNCH_ROOT}/admission_v2.tsv'
      test ! -e \"\${record}\" && test ! -e \"\${record}.submitting\"
      mkdir -- '${RUN_ROOT}/.claims/submit-admission'
      printf 'job_id\texpected_git_sha\tphase_a_producer_git_sha\tphase_a_record_sha256\tauthority_map_sha256\tsix_view_tail_job\n' > \"\${record}.submitting\"
      job=\$(sbatch --parsable \
        --export='ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_PHASE_A_PRODUCER_GIT_SHA='\"\${phase_a_producer}\"',TWIRL_FM0_PHASE_A_RECORD_SHA256='\"\${phase_a_record_sha}\"',TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256=${MAP_SHA256},TWIRL_FM0_ADMISSION_V2_RECEIPT=${ADMISSION_RECEIPT}' \
        '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_later_admission_v2_cpu.sbatch')
      [[ \"\${job}\" =~ ^[0-9]+$ ]]
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        \"\${job}\" '${EXPECTED_SHA}' \"\${phase_a_producer}\" \"\${phase_a_record_sha}\" '${MAP_SHA256}' \"\${six_tail}\" >> \"\${record}.submitting\"
      chmod a-w \"\${record}.submitting\"
      mv -- \"\${record}.submitting\" \"\${record}\"
      cat \"\${record}\""
    ;;

  retry-admission)
    remote "set -euo pipefail
      ${REMOTE_GATE}
      prior_record='${LAUNCH_ROOT}/admission_v2.tsv'
      [[ \$(wc -l < \"\${prior_record}\") -eq 2 ]]
      prior_job=\$(tail -1 \"\${prior_record}\" | cut -f1)
      [[ \"\${prior_job}\" =~ ^[0-9]+$ ]]
      prior_state=\$(sacct -nX -j \"\${prior_job}\" -o State,ExitCode | awk 'NF {print \$1,\$2; exit}')
      [[ \"\${prior_state}\" == 'FAILED 1:0' ]]
      six_record='${LAUNCH_ROOT}/six_view_fast.tsv'
      if [[ ! -e \"\${six_record}\" ]]; then six_record='${LAUNCH_ROOT}/six_view.tsv'; fi
      [[ \$(wc -l < \"\${six_record}\") -eq 13 ]]
      while IFS= read -r six_job; do
        [[ \"\${six_job}\" =~ ^[0-9]+$ ]]
        state=\$(sacct -nX -j \"\${six_job}\" -o State,ExitCode | awk 'NF {print \$1,\$2; exit}')
        [[ \"\${state}\" == 'COMPLETED 0:0' ]]
      done < <(tail -n +2 \"\${six_record}\" | cut -f2)
      six_tail=\$(tail -1 \"\${six_record}\" | cut -f2)
      phase_record='${LAUNCH_ROOT}/phase_a_freeze.tsv'
      [[ \$(wc -l < \"\${phase_record}\") -eq 2 ]]
      phase_a_producer=\$(tail -1 \"\${phase_record}\" | cut -f3)
      phase_a_map=\$(tail -1 \"\${phase_record}\" | cut -f4)
      [[ \"\${phase_a_producer}\" =~ ^[0-9a-f]{40}$ ]]
      [[ \"\${phase_a_map}\" == '${MAP_SHA256}' ]]
      phase_a_record_sha=\$(tail -n +2 \"\${six_record}\" | cut -f6 | sort -u)
      [[ \$(printf '%s\n' \"\${phase_a_record_sha}\" | wc -l) -eq 1 ]]
      [[ \"\${phase_a_record_sha}\" =~ ^[0-9a-f]{64}$ ]]
      [[ \$(sha256sum '${PHASE_A_AUTHORITY_RECORD}' | awk '{print \$1}') == \"\${phase_a_record_sha}\" ]]
      for sector in \$(seq 66 77); do
        receipt='${SIX_VIEW_ROOT}/s'\$(printf '%04d' \"\${sector}\")'/receipt.json'
        test -s \"\${receipt}\" && test ! -w \"\${receipt}\"
      done
      test ! -e '${ADMISSION_RECEIPT}'
      record='${LAUNCH_ROOT}/admission_v2_retry1.tsv'
      test ! -e \"\${record}\" && test ! -e \"\${record}.submitting\"
      mkdir -- '${RUN_ROOT}/.claims/submit-admission-retry1'
      printf 'job_id\texpected_git_sha\tphase_a_producer_git_sha\tphase_a_record_sha256\tauthority_map_sha256\tsix_view_tail_job\tprior_failed_job\n' > \"\${record}.submitting\"
      job=\$(sbatch --parsable \
        --export='ALL,TWIRL_ORCD_REPO=${REMOTE_REPO},TWIRL_EXPECTED_GIT_SHA=${EXPECTED_SHA},TWIRL_FM0_PHASE_A_PRODUCER_GIT_SHA='\"\${phase_a_producer}\"',TWIRL_FM0_PHASE_A_RECORD_SHA256='\"\${phase_a_record_sha}\"',TWIRL_FM0_PREPARATION_AUTHORITY_MAP_SHA256=${MAP_SHA256},TWIRL_FM0_ADMISSION_V2_RECEIPT=${ADMISSION_RECEIPT}' \
        '${REMOTE_REPO}/scripts/orcd/slurm_twirl_fm0_later_admission_v2_cpu.sbatch')
      [[ \"\${job}\" =~ ^[0-9]+$ ]]
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        \"\${job}\" '${EXPECTED_SHA}' \"\${phase_a_producer}\" \"\${phase_a_record_sha}\" '${MAP_SHA256}' \"\${six_tail}\" \"\${prior_job}\" >> \"\${record}.submitting\"
      chmod a-w \"\${record}.submitting\"
      mv -- \"\${record}.submitting\" \"\${record}\"
      cat \"\${record}\""
    ;;

  validate-admission)
    remote "set -euo pipefail
      ${REMOTE_GATE}
      record='${LAUNCH_ROOT}/admission_v2_retry1.tsv'
      if [[ ! -e \"\${record}\" ]]; then record='${LAUNCH_ROOT}/admission_v2.tsv'; fi
      [[ \$(wc -l < \"\${record}\") -eq 2 ]]
      job=\$(tail -1 \"\${record}\" | cut -f1)
      state=\$(sacct -nX -j \"\${job}\" -o State,ExitCode | awk 'NF {print \$1,\$2; exit}')
      [[ \"\${state}\" == 'COMPLETED 0:0' ]]
      test -s '${ADMISSION_RECEIPT}' && test ! -w '${ADMISSION_RECEIPT}'
      '${REMOTE_PYTHON}' - '${ADMISSION_RECEIPT}' '${EXPECTED_SHA}' <<'PY'
import json, sys
receipt = json.load(open(sys.argv[1], encoding='utf-8'))
assert receipt['schema_version'] == 'twirl_fm0_later_sector_preparation_pool_receipt_v2'
assert receipt['producer_git_sha'] == sys.argv[2]
assert receipt['preparation_pool_sectors'] == list(range(66, 78))
assert receipt['excluded_sectors'] == [65]
assert receipt['s78_included'] is False and receipt['s79_s80_touched'] is False
assert receipt['preparation_pool_admitted'] is True
for key in ('a2v1_accepted', 'scientific_training_eligible', 'temporal_panel_frozen', 'model_training_authorized', 'sealed_aperture_photometry_opened', 'sealed_shards_written'):
    assert receipt[key] is False
PY
      sha256sum '${ADMISSION_RECEIPT}'"
    ;;

  status)
    remote "set -euo pipefail
      printf 'expected_git_sha=%s\nauthority_map_sha256=%s\nrun_root=%s\n' '${EXPECTED_SHA}' '${MAP_SHA256}' '${RUN_ROOT}'
      if [[ -d '${LAUNCH_ROOT}' ]]; then
        find '${LAUNCH_ROOT}' -maxdepth 1 -type f -print -exec sed -n '1,20p' {} \;
      fi
      squeue -h -u tehan -o '%i|%T|%j|%R' | awk -F'|' '\$3 ~ /^twirl-fm0-(quality-ref|source-inventory|phase-a-freeze|six-view|admission-v2)$/ {print}'
      for path in '${PHASE_A_AUTHORITY_RECORD}' '${ADMISSION_RECEIPT}'; do
        if [[ -s \"\${path}\" ]]; then sha256sum \"\${path}\"; else printf 'MISSING %s\n' \"\${path}\"; fi
      done"
    ;;
esac
