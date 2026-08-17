#!/usr/bin/env bash
# Create immutable, sector-level A2v1 HDF5 tar archives for FM0.1 staging.

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

accepted_root=${TWIRL_FM0_PDO_ACCEPTED_ROOT:-/pdo/users/tehan/tglc-gpu-production-A2v1}
stage_root=${TWIRL_FM0_PDO_TAR_ROOT:-/pdo/users/tehan/twirl-fm0-stage/source_a2v1_s56_s65_tar_v1}
[[ -d $accepted_root && ! -L $accepted_root ]] || {
  echo "missing materialized accepted root: $accepted_root" >&2
  exit 2
}
mkdir -p "$stage_root"
[[ -d $stage_root && ! -L $stage_root ]] || {
  echo "unsafe staging root: $stage_root" >&2
  exit 2
}

for ((sector=from_sector; sector<=through_sector; sector++)); do
  tag=s$(printf '%04d' "$sector")
  orbit1=$((2 * sector + 7))
  orbit2=$((orbit1 + 1))
  final=$stage_root/$tag
  partial=$stage_root/${tag}.partial
  archive_name=${tag}_A2v1_raw_hdf5.tar

  if [[ -s $final/READY ]]; then
    (
      cd "$final"
      sha256sum -c "${archive_name}.sha256"
    )
    echo "SKIP_READY $tag"
    continue
  fi
  if [[ -e $final || -L $final || -e $partial || -L $partial ]]; then
    echo "colliding or partial FM0.1 tar stage: $final $partial" >&2
    exit 2
  fi
  mkdir "$partial"

  list0=$partial/hdf5_files.nul
  list=$partial/hdf5_files.txt
  for orbit in "$orbit1" "$orbit2"; do
    orbit_root=$accepted_root/orbit-$orbit/ffi
    [[ -d $orbit_root && ! -L $orbit_root ]] || {
      echo "missing accepted orbit root: $orbit_root" >&2
      exit 2
    }
    (
      cd "$accepted_root"
      find "orbit-$orbit/ffi" -path '*/LC/*.h5' -type f -print0
    )
  done | sort -z > "$list0"
  tr '\0' '\n' < "$list0" > "$list"
  n_files=$(tr -cd '\0' < "$list0" | wc -c | tr -d ' ')
  [[ $n_files =~ ^[1-9][0-9]*$ ]] || {
    echo "empty HDF5 inventory for $tag" >&2
    exit 2
  }
  source_bytes=$(python3 -c \
    'import os,sys; root=sys.argv[1]; names=open(sys.argv[2], "rb").read().split(b"\0"); print(sum(os.stat(os.path.join(root, os.fsdecode(name))).st_size for name in names if name))' \
    "$accepted_root" "$list0")
  [[ $source_bytes =~ ^[1-9][0-9]*$ ]] || {
    echo "invalid source-byte inventory for $tag" >&2
    exit 2
  }

  archive_partial=$partial/${archive_name}.partial
  stream_digest=$partial/${archive_name}.stream.sha256
  (
    cd "$accepted_root"
    nice -n 10 tar --create --format=pax --numeric-owner --owner=0 --group=0 \
      --mtime=@0 --no-recursion --null --verbatim-files-from \
      --files-from="$list0" --file=-
  ) | tee "$archive_partial" | sha256sum > "$stream_digest"
  digest=$(awk '{print $1}' "$stream_digest")
  [[ $digest =~ ^[0-9a-f]{64}$ ]] || {
    echo "invalid streamed archive digest for $tag" >&2
    exit 2
  }
  mv "$archive_partial" "$partial/$archive_name"
  printf '%s  %s\n' "$digest" "$archive_name" > "$partial/${archive_name}.sha256"
  archive_bytes=$(stat -c '%s' "$partial/$archive_name")
  cat > "$partial/summary.json" <<EOF
{
  "archive": "$archive_name",
  "archive_bytes": $archive_bytes,
  "archive_sha256": "$digest",
  "accepted_sources_modified": false,
  "n_hdf5_files": $n_files,
  "orbits": [$orbit1, $orbit2],
  "schema_version": "twirl_fm0_1_a2v1_sector_tar_v1",
  "sector": $sector,
  "source_bytes": $source_bytes
}
EOF
  rm "$list0" "$stream_digest"
  printf '%s\n' "$digest" > "$partial/READY"
  mv "$partial" "$final"
  echo "PDO_TAR_READY $tag files=$n_files bytes=$archive_bytes sha256=$digest"
done

echo "PDO_FM0_TARS_READY S${from_sector}-S${through_sector}"
