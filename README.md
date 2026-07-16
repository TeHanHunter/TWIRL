# TWIRL

**TESS White dwarf Investigation of Remnant pLanets** is a survey pipeline for
transiting or occulting objects around white dwarfs in the `200 s` TESS FFIs
(`Sector >= 56`). The fixed benchmark is **WD 1856+534** in Sector 56, orbits
119 and 120, camera 4/CCD 1.

This is the standalone pipeline repository. Large catalogs and light curves
are local products rather than git assets; see [local-data policy](doc/local_data.md).
The active TWIRL I manuscript is maintained separately in the sibling
`twirl-survey-paper` repository.

## Current checkpoint

As of `2026-07-16`:

- **S56 A2v1 is product-complete and passes Tier-0 QA.** HDF5 coverage and all
  `31,450` ADP/ADP015-only FITS products passed product validation, and WD 1856
  is recovered at BLS rank one in both ADP apertures. Tier-1 science QA remains
  required before survey release.
- **S57 A2v1 is product-complete.** Its `27,213` FITS products and compact ADP
  export passed validation. At audit close, an experimental S57 review queue
  had been created prematurely and contained six human labels. Preserve all
  six, pause further S57 holdout consumption, and do not describe S57 as a
  pristine external holdout.
- **S58-S63 A2v1 production is in progress.** A targeted S58 HDF5 repair passed
  the strengthened openability gate, and the serial queue resumed with S58
  FITS rebuilding. HDF5 completion remains an intermediate checkpoint.
- **S56 active learning is a pilot, not a survey classifier.** The
  `s56_harmonic_cnn_v1` teacher remains the enrichment baseline, and its first
  `1,000`-TIC review batch is partially labeled. Teacher v2 was completed as an
  exploratory comparison but missed promotion gates; student/self-training
  remains blocked.

The authoritative current status and priorities are in
[the project plan](doc/twirl_plan.md); dated runs and benchmarks belong in
[the progress log](doc/twirl_progress_log.md). Files under [reports](reports/README.md)
are evidence snapshots and do not override either document.

## A2v1 product contract

`A2v1` is the current Stage 1 production family. It means:

- no science target magnitude cap (the policy equivalent of
  `--max-magnitude 99`);
- observation-table TIC overlays on reused source pickles;
- saturated-pixel masking when ePSFs are fit;
- ADP and ADP015 flux branches for the `1x1`, `3x3`, and `5x5` apertures;
- sector-level FITS export after both orbit HDF5 trees pass coverage checks.

The operational source-pickle reuse, HDF5, FITS, and validation sequence is in
[the A2v1 production protocol](doc/a2v1_production_protocol.md). Stage 1 runs
on MIT PDO; compact exports support downstream work on ORCD. Raw TGLC/TICA
trees are not moved to ORCD.

## Pipeline status

1. **Stage 1 — catalogs, extraction, A2v1 products, indexing, and photometric
   QA.** S56 passes product validation and Tier-0 integrity/benchmark QA, S57
   is product-validated, and the S58-S63 queue is active. Tier-1 science QA,
   the consolidated archive/index, a frozen release manifest/sector cutoff,
   and the Gaia-first no-TIC bridge remain open.
2. **Stage 2 — transparent search and candidate generation.** The sector-level
   periodic BLS, cross-aperture consolidation, diagnostics, heuristic checks,
   and WD 1856 smoke paths exist. The co-equal non-periodic dip baseline and
   multi-sector aggregation are not yet complete.
3. **Stage 3 — injection/recovery and completeness.** LC-level pre-detrend
   BATMAN workflows, candidate-retention diagnostics, branch comparisons, and
   a pixel-level smoke path exist. Publication-grade end-to-end completeness
   and the representative pixel-level calibration subset remain unfinished.
4. **Stage 4 — frozen full-survey inference and candidate construction.** This
   stage is planned but has no dedicated script directory yet. Current teacher
   scoring and enrichment are active-learning experiments, not a Stage 4
   survey run.
5. **Stage 5 — human review, scientific vetting, validation, and follow-up.**
   S56 review/adjudication, LEO comparison, centroid checks, and the harmonic
   teacher are implemented as pilot infrastructure. Candidate validation and
   follow-up do not begin until the upstream search and completeness gates are
   stable.

See [the script inventory](scripts/README.md) for current entrypoints and the
incremental Stage 1–5 migration policy. Reusable implementation belongs in
`src/twirl/`; scripts should remain thin drivers.

## Install and test

```bash
git clone git@github.com:TeHanHunter/TWIRL.git
cd TWIRL
python -m pip install -e '.[dev]'
make test-fast
```

The base package provides catalog, light-curve, I/O, plotting, search, and
vetting modules. Some production and research workflows additionally require
environment-specific software such as the MIT TGLC/QLP stack or PyTorch/CUDA.
The `dev` extra includes the optional BATMAN injection dependency; use
`pip install -e '.[injections]'` when only that optional workflow is needed.
Follow the PDO or ORCD runbook for remote production jobs.

Example imports:

```python
from twirl.io.hlsp import read_hlsp
from twirl.search.bls import BLSConfig, run_bls_on_lc
from twirl.search.consolidate import consolidate_candidates
```

## Repository map

- `src/twirl/`: importable production code
- `scripts/`: thin local, PDO, and ORCD entrypoints
- `configs/`: reusable run configuration
- `catalogs/`: catalog-product documentation, not large catalog payloads
- `doc/`: authoritative plans, runbooks, methods, and local-data conventions
- `reports/`: compact, dated QA and analysis evidence
- `data_local/`: ignored local catalogs and staged data

Collaboration-specific data, labels, and code may live in sibling repositories
or external handoffs. Do not vendor them into TWIRL; record accepted interfaces
and provenance here.
