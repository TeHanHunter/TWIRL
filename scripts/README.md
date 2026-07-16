# Script Entrypoints And Lifecycle

`scripts/` contains thin command-line, PDO, and ORCD drivers. Reusable logic
belongs in `src/twirl/`; a script should parse arguments, select configuration,
call package code, and record provenance.

This inventory documents the current tree. It does not move or rename any
entrypoint: active PDO/ORCD jobs and older runbooks may still depend on the
existing paths.

## Lifecycle Labels

- **Active:** accepted for the named current workflow.
- **Compatibility:** retained because a runbook, remote checkout, or older
  product still calls it; do not use it to define a new product.
- **Experimental:** useful for a bounded pilot or comparison, but not a frozen
  survey contract.
- **Deprecated:** has a documented replacement and should remain only as a
  forwarding wrapper until callers migrate.

Names alone do not determine lifecycle. Confirm the live status in
[twirl_plan.md](../doc/twirl_plan.md), the applicable runbook, and the latest
[progress entry](../doc/twirl_progress_log.md) before launching a remote job.

## Current Directory Inventory

| Path | Lifecycle now | Intended Stage ownership |
| --- | --- | --- |
| `stage1_lightcurves/` | Active A2v1 production and QA mixed with compatibility/diagnostic launchers for older QLP and TWIRL-FS products | Stage 1 only: catalogs, detector maps, extraction, A2v1 FITS, indexing, and photometric QA |
| `stage2_search/` | Active transparent periodic-BLS baseline; dip search and multi-sector aggregation remain incomplete | Stage 2: periodic/dip candidate generation and search diagnostics |
| `stage3_injections/` | Active compact export, LC-level injection, and pixel-injection smoke entrypoints | Stage 3: all injection generation, recovery accounting, and completeness drivers |
| `stage4_search/` | Not present yet | Add for frozen full-survey inference, cross-sector merging, ranking, and candidate-table construction |
| `stage5_validation/` | Active but overloaded: human review, LEO/centroid validation, injection experiments, rankers, and teacher-training/inference coexist here | Keep human adjudication and scientific validation in Stage 5; migrate detector-development drivers to Stage 2, injection/recovery drivers to Stage 3, and frozen production-inference drivers to Stage 4 incrementally |
| `orcd/` | Active platform launchers and Slurm wrappers | Infrastructure only; scientific ownership follows the payload's Stage 2–5 workflow |
| top-level wrappers | Mostly environment/compatibility helpers (`activate_qlp_env.sh`, detrend/HLSP wrappers, report sync) | Retain while legacy QLP and remote runbooks require them |

Teacher-v2 drivers and artifacts are experimental historical evidence: the
completed comparison missed its promotion gates. Teacher v1 remains the
bounded active-learning baseline, and neither teacher defines the frozen Stage
4 survey inference contract.

## Current Canonical Stage 1 Path

For a new A2v1 sector with prepared source pickles, use the sequence documented
in [the A2v1 protocol](../doc/a2v1_production_protocol.md):

1. [run_a2v1_reproduction_pdo.sh](stage1_lightcurves/run_a2v1_reproduction_pdo.sh)
   builds observation-table overlays, handles saturated-mask ePSFs, and writes
   orbit HDF5 products.
2. [run_a2v1_hlsp_pdo.sh](stage1_lightcurves/run_a2v1_hlsp_pdo.sh) writes the
   required sector-level ADP/ADP015-only FITS products.
3. [validate_a2v1_product.py](stage1_lightcurves/validate_a2v1_product.py)
   checks HDF5 coverage and FITS schema/completeness.
4. [audit_a2v1_photometric_qa.py](stage1_lightcurves/audit_a2v1_photometric_qa.py)
   performs the current Tier-0 integrity/benchmark QA. Tier-1 science QA is a
   separate open gate and must not be inferred from this script's pass state.

HDF5 completion alone is not an accepted sector product. Older generic GPU,
QLP HLSP, canonical-flux, and TWIRL-FS comparison launchers remain available
for provenance or controlled comparison, but they are not substitutes for the
A2v1 contract.

## Current Canonical Stage 2 Baseline

The transparent periodic baseline is launched by
[stage2_bls_worker.sh](stage2_search/stage2_bls_worker.sh) through
`python -m twirl.search.sector_run`. Supply the validated A2v1 FITS root
explicitly with `--hlsp-root` and load the versioned settings from
[`bls_default.yaml`](../configs/detection/bls_default.yaml). The command must
not fall back to a legacy product root, and its output must record the selected
apertures and BLS/cadence-cleaning configuration. The dip branch,
multi-sector aggregation, and false-alarm calibration remain the next Stage 2
production interfaces to add.

## Incremental Migration Rule

Do not bulk-move the large `stage5_validation/` tree while S56/S57 inference,
review, or remote jobs are active. Migrate one coherent workflow at a time:

1. put shared computation and schemas in the appropriate `src/twirl/` module;
2. add the generic Stage 3 or Stage 4 CLI, driven by versioned config rather
   than another copied `run_sNN_*` launcher;
3. leave the old path as a thin compatibility wrapper with a clear replacement
   message;
4. update PDO/ORCD runbooks, Slurm wrappers, and progress-log links;
5. retire the wrapper only after remote checkouts and active jobs no longer
   reference it.

The target boundary is behavioral:

- search kernels and candidate generation live in `twirl.search`;
- injection models and recovery accounting move toward `twirl.injections`
  rather than expanding Stage 5 scripts;
- automated vetting remains in `twirl.vetting`;
- review/label/adjudication interfaces move toward `twirl.review`, separate
  from model training;
- detector and teacher implementations move toward `twirl.models`; Stage 2
  owns their development/evaluation and Stage 4 owns frozen inference;
- Stage 4 scripts should consume frozen Stage 1 products plus accepted Stage
  2/3 contracts, not encode a new scientific policy ad hoc.

Sector-specific launchers already used by production may remain until a
generic replacement is proven. New workflows should prefer generic CLIs plus
YAML configuration and periodic progress output for long runs.
