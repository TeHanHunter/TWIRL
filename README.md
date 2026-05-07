# TWIRL

**TESS White dwarf Investigation of Remnant pLanets** — a survey pipeline
to search for transiting/occulting objects around white dwarfs using
200 s TESS FFIs (Sector ≥ 56). The benchmark target is **WD 1856+534**.

The repository is the working tree for the production pipeline + planning
documents. Light curves and large catalogs are not committed; see
[`doc/local_data.md`](doc/local_data.md).

## Quick install

```bash
git clone git@github.com:TeHanHunter/TWIRL.git
cd TWIRL
pip install -e .
```

This installs the `twirl` package so downstream code can do, e.g.:

```python
from twirl.io.hlsp import read_hlsp
from twirl.search.bls import run_bls_on_lc, BLSConfig
from twirl.search.consolidate import consolidate_candidates
from twirl.search.diagnostics import plot_vet_sheet
```

The `twirl.search` subpackage exposes the per-sector BLS first pass
(detection layer + cross-aperture consolidation + TOI-style vet sheets).
Stage 2 vetting / completeness work is intended to import from this
package and build on top.

## Pipeline stages (high-level)

1. **Stage 1** — Generate WD light curves on MIT PDO with TGLC.
2. **Stage 2** — Periodic-transit search (BLS + cross-aperture consolidation).
   Heuristic vetter and multi-sector merger are next.
3. **Stage 3** — Injection-recovery completeness mapping.
4. **Stage 4** — Full-sample search.
5. **Stage 5** — Validation and follow-up.

See [`doc/twirl_plan.md`](doc/twirl_plan.md) for the executable plan and
[`doc/twirl_progress_log.md`](doc/twirl_progress_log.md) for dated
execution notes.

## Student work

Stage 2/3 sub-projects live in sibling repositories:

- `twirl-chen` (Franklin) — vetting / triage on top of TWIRL BLS
- `twirl-arimond` (Corina) — injection-recovery completeness using TWIRL BLS
- `twirl-choi` (Emma) — TGLC vs SPOC light-curve benchmarking
- `twirl-le` (Angelina) — pixel-level source-ID for centroid analysis

Each student installs TWIRL via `pip install -e ../TWIRL` and either uses
the `twirl.search` API directly or contributes back via PRs.
