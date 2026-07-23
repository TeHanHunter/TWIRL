#!/usr/bin/env bash
# Retired legacy S57 chain. The active teacher candidate-quality contract is
# deliberately S56-only, so this wrapper must not launch partial production.
set -euo pipefail

echo "Blocked: the legacy S57 teacher chain is retired. The current candidate-quality and Tier-1 contracts support S56 only." >&2
exit 64
