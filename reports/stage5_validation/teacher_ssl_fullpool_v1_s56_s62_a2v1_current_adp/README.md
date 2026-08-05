# Teacher v4-SSL full-pool preregistration

**Status:** broad-corpus preprocessing prepared; exact ORCD pool freeze and BLS
search pending.

## Objective

This run replaces the 6,168-observation Teacher v4-SSL infrastructure pilot
with self-supervised pretraining on the broad real S56--S62 A2v1 population.
The source compact exports contain 212,049 sector observations of 149,011
unique TICs. The final pretraining population is intentionally smaller because
it excludes every earlier-sector observation of:

- a TIC in the frozen Teacher v3 fixed test; or
- a TIC reserved by identity for prospective S63 evaluation.

The exact retained count is an output of the checksum-bound ORCD freeze, not a
hand-entered constant. The broad run may therefore be described informally as
the “full 200k” experiment, but all artifacts and results must report its exact
leakage-safe count.

## Frozen boundaries

- Inputs are real A2v1 ADP/ADP-small observations from S56--S62 only.
- The S63 reservation was built from target identities only. No S63 light
  curve was opened.
- Exclusions operate on whole TICs across all seven input sectors.
- BLS uses the locked 50,000-period, ten-peak, two-aperture search.
- S56--S61 require exact compact/reference orbit agreement.
- S62 uses the authoritative cadence reference only for the audited inherited
  orbit-ID mismatch; all corrections are counted and fingerprinted.
- The frozen test, S63 hosts, injections, and human labels are forbidden from
  the SSL objective.

The immutable S63 identity inventory is in
`preregistered/s63_reserved_tics.txt`; its adjacent JSON records the source
observation-table hash and confirms that no S63 light curves were read.

## Fold semantics

The five encoders are independent one-H200 jobs, with at most four running
concurrently. For encoder fold `k`:

- all leakage-safe previously unlabeled TICs are available;
- labeled development TICs assigned to fold `k` are withheld in every sector;
- labeled development TICs in the other four folds are available with their
  labels ignored.

Unseen TICs are not assigned arbitrary pseudo-folds. This preserves the broad
unlabeled sample while keeping each later supervised development evaluation
inductive with respect to its held labeled hosts.

## Pretraining contract

- Architecture/input profile: the existing `s56_harmonic_cnn_v1`
  shape-plus-periodogram encoder and event-preserving VICReg augmentations.
- Independent encoders: five folds, one H200 per fold, at most four concurrent.
- Optimization: `20` epochs, batch size `64`, learning rate `3e-4`, weight
  decay `1e-4`, and seed `560064`.
- Checkpointing: atomic every epoch with Python, NumPy, Torch CPU, and Torch
  CUDA random-number states restored on exact-contract resume.
- Scale language: “full 200k” is shorthand for this experiment family; every
  result must use the exact post-exclusion observation and TIC counts emitted
  by the frozen-pool manifest.

## Execution sequence

1. Freeze and hash the seven compact exports, split registry, and S63
   reservation; emit exact per-sector allowlists.
2. Run an isolated `100`-target S62 BLS canary, then release and
   checksum-merge the full two-aperture BLS products only if it succeeds.
3. Stage byte-identical pool/summary/allowlist authorities to PDO, export
   sharded raw/error channels from the user-owned TGLC tree, transfer those
   shards to ORCD, and build the native seven-harmonic inputs.
4. Run a bounded 10,000-observation, one-epoch throughput smoke.
5. Pretrain five VICReg encoders, one fold per H200, with exact epoch resume.
6. Fine-tune on the existing frozen labels and compare against Teacher v3 on
   matched development support before any promotion decision.

The pretraining checkpoints alone are representation artifacts, not a
science-ready search model. Teacher v3 remains the operational enrichment
baseline until the preregistered downstream comparison is complete.
