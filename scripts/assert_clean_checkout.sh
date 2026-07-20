#!/usr/bin/env bash
# Refuse production or deployment work from a dirty Git checkout.
set -euo pipefail

REPO="${1:-$(git rev-parse --show-toplevel 2>/dev/null || true)}"
EXPECTED_SHA="${TWIRL_EXPECTED_GIT_SHA:-}"

if [[ -z "${REPO}" || ! -d "${REPO}/.git" ]]; then
  echo "[git-clean] not a Git checkout: ${REPO:-<unset>}" >&2
  exit 2
fi

cd "${REPO}"

status="$(git status --porcelain=v1 --untracked-files=all)"
if [[ -n "${status}" ]]; then
  echo "[git-clean] refusing dirty checkout: ${REPO}" >&2
  printf '%s\n' "${status}" >&2
  exit 3
fi

head_sha="$(git rev-parse HEAD)"
if [[ -n "${EXPECTED_SHA}" && "${head_sha}" != "${EXPECTED_SHA}" ]]; then
  echo "[git-clean] checkout SHA mismatch: ${head_sha} != ${EXPECTED_SHA}" >&2
  exit 4
fi

echo "[git-clean] repo=${REPO}"
echo "[git-clean] sha=${head_sha}"
