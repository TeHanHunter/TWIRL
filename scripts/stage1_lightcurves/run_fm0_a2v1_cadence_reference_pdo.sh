#!/usr/bin/env bash
# Build one accepted sector's checksum-bound cadence/quality authority for FM0.1.

set -euo pipefail

sector=${1:?usage: $0 SECTOR [OUTPUT_DIR]}
if [[ ! $sector =~ ^[0-9]+$ ]] || (( sector < 56 || sector > 65 )); then
  echo "FM0.1 cadence preparation is bounded to accepted S56-S65" >&2
  exit 2
fi
orbit1=$((2 * sector + 7))
orbit2=$((orbit1 + 1))
tag=s$(printf '%02d' "$sector")

repo=${TWIRL_ROOT:-/pdo/users/tehan/TWIRL}
expected_git_sha=${TWIRL_EXPECTED_GIT_SHA:-}
output=${2:-/pdo/users/tehan/twirl-fm0-stage/cadence_reference/$tag}
python_bin=${PYTHON_BIN:-/sw/qlp-environment/.venv/bin/python}
qlp_root=${TWIRL_QLP_AUTHORITY_ROOT:-/pdo/qlp-data}
spoc_root=${TWIRL_SPOC_FLAG_ROOT:-$qlp_root/spocflags}

[[ $expected_git_sha =~ ^[0-9a-f]{40}$ ]] || {
  echo "TWIRL_EXPECTED_GIT_SHA must be a full Git SHA" >&2
  exit 2
}
[[ -d $repo && ! -L $repo && -s $repo/CODE_REVISION ]] || {
  echo "missing materialized FM0.1 code export or CODE_REVISION: $repo" >&2
  exit 2
}
[[ $(tr -d '[:space:]' < "$repo/CODE_REVISION") == "$expected_git_sha" ]] || {
  echo "FM0.1 code-export revision mismatch" >&2
  exit 2
}
user_root=$(realpath -e -- /pdo/users/tehan)
output=$(realpath -m -- "$output")
if [[ $output != "$user_root"/* || $output == "$repo" || $output == "$repo"/* ]]; then
  echo "output must be under /pdo/users/tehan and outside the code export" >&2
  exit 2
fi
if [[ -e $output || -L $output ]]; then
  echo "refusing to overwrite cadence authority: $output" >&2
  exit 4
fi

quat_args=()
qflag_args=()
spoc_args=()
for orbit in "$orbit1" "$orbit2"; do
  for camera in 1 2 3 4; do
    quat=$qlp_root/orbit-$orbit/ffi/run/cam${camera}_quat.txt
    [[ -s $quat ]] || { echo "missing quaternion authority: $quat" >&2; exit 3; }
    quat_args+=(--quat "$orbit,$camera,$quat")
    for ccd in 1 2 3 4; do
      qflag=$qlp_root/orbit-$orbit/ffi/run/cam${camera}ccd${ccd}_qflag.txt
      [[ -s $qflag ]] || { echo "missing QLP quality authority: $qflag" >&2; exit 3; }
      qflag_args+=(--qlp-qflag "$orbit,$camera,$ccd,$qflag")
    done
  done
done
for camera in 1 2 3 4; do
  for ccd in 1 2 3 4; do
    spoc=$spoc_root/spocffiflag_s${sector}_cam${camera}_ccd${ccd}.txt
    [[ -s $spoc ]] || { echo "missing SPOC quality authority: $spoc" >&2; exit 3; }
    spoc_args+=(--spoc-flag "$camera,$ccd,$spoc")
  done
done

mkdir -p "$output"
cd "$repo"
export PYTHONPATH=$repo/src
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

"$python_bin" -c 'import astropy, numpy, pandas; print("FM0 cadence Python import preflight OK")'

"$python_bin" scripts/stage1_lightcurves/build_a2v1_spoc_quality_table.py \
  --sector "$sector" --expected-orbit "$orbit1" --expected-orbit "$orbit2" \
  "${quat_args[@]}" "${spoc_args[@]}" \
  --output-table "$output/spoc_quality.csv" \
  --output-provenance "$output/spoc_quality.provenance.json"

"$python_bin" scripts/stage1_lightcurves/build_a2v1_cadence_reference.py \
  --sector "$sector" --expected-orbit "$orbit1" --expected-orbit "$orbit2" \
  "${quat_args[@]}" "${qflag_args[@]}" \
  --spoc-quality-table "$output/spoc_quality.csv" \
  --spoc-quality-provenance "$output/spoc_quality.provenance.json" \
  --output-table "$output/cadence_reference.csv" \
  --output-manifest "$output/cadence_reference.manifest.json"

"$python_bin" - "$output" "$expected_git_sha" "$sector" <<'PY'
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
revision = sys.argv[2]
sector = int(sys.argv[3])
names = (
    "spoc_quality.csv",
    "spoc_quality.provenance.json",
    "cadence_reference.csv",
    "cadence_reference.manifest.json",
)
digests = {}
for name in names:
    digest = hashlib.sha256()
    with (root / name).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    digests[name] = digest.hexdigest()
receipt = {
    "schema_version": "twirl_fm0_1_cadence_authority_build_v1",
    "sector": sector,
    "producer_git_sha": revision,
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "outputs_sha256": digests,
}
(root / "fm0_build_receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY
sha256sum "$output"/* > "$output/files.sha256"
printf '%s\n' "$expected_git_sha" > "$output/READY"
echo "FM0_CADENCE_READY sector=$sector output=$output"
