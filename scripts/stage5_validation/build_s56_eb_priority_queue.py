#!/usr/bin/env python3
"""Build the first S56 EB/PCEB-priority human-review queue."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
for path in (SCRIPT_DIR, SRC_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from twirl.vetting.eb_miner import evaluate_eb_miner_release_gate, write_eb_priority_queue  # noqa: E402
from twirl.vetting.recovery50_teacher import json_default, read_table  # noqa: E402


DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_eb_miner_adp_only"
DEFAULT_SCORED = DEFAULT_ROOT / "eb_cnn_miner/scored_candidates_ensemble.csv"
DEFAULT_MINER_SUMMARY = DEFAULT_ROOT / "eb_cnn_miner/summary.json"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "review_queue_eb_priority_pilot100"


def _write_launcher(out_dir: Path) -> None:
    app_src = SCRIPT_DIR / "franklin_vetting_app.py"
    if app_src.exists():
        shutil.copy2(app_src, out_dir / "franklin_vetting_app.py")
    launcher = out_dir / "run_eb_priority_vetting.sh"
    launcher.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "ROOT=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")/../../../..\" && pwd)\"\n"
        "cd \"${ROOT}\"\n"
        "exec scripts/stage5_validation/run_s56_eb_priority_vetting_app_local.sh\n"
    )
    launcher.chmod(0o755)


def _write_readme(out_dir: Path, summary: dict) -> None:
    queue_name = Path(summary["queue_csv"]).name
    text = """# S56 EB/PCEB Priority Queue

This queue is an active-learning product from the EB/PCEB miner. Scores are
human-review priorities, not final labels.

## Start

```bash
python3 franklin_vetting_app.py --queue %s --labels-out human_labels_vetted.csv --sheet-root twirl_vet_sheets --check-only
./run_eb_priority_vetting.sh
```

Then open:

```text
http://127.0.0.1:5004/
```

## Label Policy

- `eclipsing_binary_or_pceb`: EB/PCEB morphology, secondary eclipse, odd/even
  mismatch, or stellar-companion-like folded shape.
- Other labels keep their normal S56 meanings.
- `uncertain` still means flat/no obvious useful event.

## Summary

```json
%s
```
""" % (queue_name, json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    (out_dir / "README_EB_priority_vetting.md").write_text(text)


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scored-candidates", type=Path, default=DEFAULT_SCORED)
    ap.add_argument("--miner-summary", type=Path, default=DEFAULT_MINER_SUMMARY)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--n-review", type=int, default=100)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    miner_summary = json.loads(args.miner_summary.read_text())
    release_gate = evaluate_eb_miner_release_gate(miner_summary)
    if not release_gate["passed"]:
        raise RuntimeError(f"EB miner release gate failed: {release_gate}")
    scored = read_table(args.scored_candidates)
    summary = write_eb_priority_queue(scored=scored, out_dir=args.out_dir, n_review=args.n_review)
    summary["scored_candidates"] = str(args.scored_candidates)
    summary["miner_summary"] = str(args.miner_summary)
    summary["release_gate"] = release_gate
    _write_launcher(args.out_dir)
    _write_readme(args.out_dir, summary)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
