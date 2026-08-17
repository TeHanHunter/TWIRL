#!/usr/bin/env bash
# Pull immutable FM0.1 sector tar archives from PDO into Blackhole staging.

set -euo pipefail

if [[ ${1:-} != --apply ]]; then
  echo "usage: $0 --apply [--from-sector 56 --through-sector 65]" >&2
  exit 2
fi

from_sector=56
through_sector=65
shift
while [[ $# -gt 0 ]]; do
  case $1 in
    --from-sector) from_sector=$2; shift 2 ;;
    --through-sector) through_sector=$2; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
if (( from_sector < 56 || through_sector > 65 || through_sector < from_sector )); then
  echo "FM0.1 accepted-sector tar range must be within S56-S65" >&2
  exit 2
fi

pdo_socket=${TWIRL_BLACKHOLE_PDO_SOCKET:-/home/tehan/.ssh/cm/tehan@pdogpu1.mit.edu:22}
pdo_host=${TWIRL_PDO_HOST:-pdogpu1}
pdo_root=${TWIRL_FM0_PDO_TAR_ROOT:-/pdo/users/tehan/twirl-fm0-stage/source_a2v1_s56_s65_tar_v1}
blackhole_root=${TWIRL_FM0_BLACKHOLE_TAR_ROOT:-/globus/tso/twirl/fm0/source_a2v1_s56_s65_tar_v1}
poll_seconds=${TWIRL_FM0_TAR_POLL_SECONDS:-60}

[[ -S $pdo_socket ]] || { echo "missing Blackhole-to-PDO socket: $pdo_socket" >&2; exit 2; }
ssh -S "$pdo_socket" -O check "$pdo_host" >/dev/null
mkdir -p "$blackhole_root"

ssh_options=(
  -S "$pdo_socket" -o BatchMode=yes -o PasswordAuthentication=no
  -o KbdInteractiveAuthentication=no -o NumberOfPasswordPrompts=0
)

for ((sector=from_sector; sector<=through_sector; sector++)); do
  tag=s$(printf '%04d' "$sector")
  source=$pdo_root/$tag
  partial=$blackhole_root/${tag}.partial
  archive_name=${tag}_A2v1_raw_hdf5.tar

  while ! ssh "${ssh_options[@]}" "$pdo_host" test -s "$source/READY"; do
    sleep "$poll_seconds"
  done
  if [[ -s $partial/READY_BLACKHOLE ]]; then
    (
      cd "$partial"
      sha256sum -c "${archive_name}.sha256"
    )
    echo "SKIP_BLACKHOLE_READY $tag"
    continue
  fi
  if [[ -e $partial || -L $partial ]]; then
    echo "colliding partial Blackhole tar stage: $partial" >&2
    exit 2
  fi
  mkdir "$partial"
  rsync -rlt --partial --append-verify --stats \
    -e "ssh -S $pdo_socket -o BatchMode=yes -o PasswordAuthentication=no -o KbdInteractiveAuthentication=no -o NumberOfPasswordPrompts=0" \
    "$pdo_host:$source/" "$partial/"
  (
    cd "$partial"
    sha256sum -c "${archive_name}.sha256"
  )
  cp "$partial/READY" "$partial/READY_BLACKHOLE"
  echo "BLACKHOLE_TAR_READY $tag"
done

echo "BLACKHOLE_FM0_TARS_READY S${from_sector}-S${through_sector}"
