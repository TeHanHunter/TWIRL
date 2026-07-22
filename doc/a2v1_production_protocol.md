# A2v1 production protocol

This is the operational protocol for the `A2v1` Stage 1 product on PDO. It
applies to S56 and later sectors with prepared default TGLC source trees,
including S94+. It supersedes the legacy production path whenever the target
product is the no-cap, saturated-mask, ADP/ADP015-only product.

## Product contract

`A2v1` means:

- include every positive TIC in the TWIRL observation table; do not use a
  TIC-magnitude science cut;
- retain prepared source-pickle pixel data through symlinks rather than
  recreating cutouts;
- use `source_tic/source_X_Y.ecsv` overlays, built from observation-table
  detector coordinates, to replace the stale embedded TIC tables;
- reuse an old ePSF only when the source cutout's saturated-pixel mask is
  empty; recompute ePSFs for every non-empty mask with the masked TGLC hook;
- write HDF5 first; the standard full product then adds one sector-level FITS
  product per TIC with only ADP and ADP015 branches for the 1x1, 3x3, and 5x5
  apertures.

The output root is always `/pdo/users/tehan/tglc-gpu-production-A2v1/`:

```text
orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/source/      # symlinks to prepared source pickles
orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/source_tic/  # A2v1 TIC overlays
orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/epsf/        # links or masked refits
orbit-<orbit>/ffi/cam<camera>/ccd<ccd>/LC/          # A2v1 HDF5
hlsp_s<sector4>_A2v1/                                # A2v1-only FITS
twirl_logs/                                          # run and validation records
```

The original source tree and all shared `/pdo/qlp-data/` inputs are read-only.
Never clean or overwrite an existing A2v1 sector while a run is active.

## Required hooks and runtime

Before a production or smoke run, verify the user-owned fork has both hooks:

```bash
grep -q '_source_tic_overlay_path' \
  /pdo/users/tehan/tess-gaia-light-curve-twirl/tglc/scripts/light_curves.py
grep -q 'source.mask.mask' \
  /pdo/users/tehan/tess-gaia-light-curve-twirl/tglc/scripts/epsfs.py
```

Use the GPU venv for `tglc epsfs` and `tglc lightcurves`, keep all BLAS thread
counts at one, and set `HDF5_USE_FILE_LOCKING=FALSE` on PDO storage. The A2v1
wrapper does this itself; direct smoke commands must do the same.

The no-cap policy is enforced by the observation-table overlays. A2v1 normally
does not rerun `tglc catalogs`, so `--max-magnitude 99` is a policy invariant,
not a required CLI argument for the reuse path.

## Pre-flight

For each selected orbit/CCD, require exactly `196` prepared source pickles.
The legacy ePSF state must be either complete (`196` files, permitting
empty-mask reuse) or absent (`0` files, requiring a saturated-mask refit for
every source); reject partial legacy ePSF trees. Confirm TICA FFIs are present
and inspect GPU memory before selecting explicit free devices. For example:

```bash
for orbit in 121 122; do
  base=/pdo/users/tehan/tglc-gpu-production/orbit-$orbit/ffi/cam4/ccd1
  find "$base/source" -name 'source_*.pkl' | wc -l
  find "$base/epsf" -name 'epsf_*.npy' | wc -l
done
nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu \
  --format=csv,noheader
```

Do not proceed if a selected source tree is incomplete. Repair the original
prepared tree first; never regenerate it inside the A2v1 output root.

## One-CCD smoke

Use one orbit and one CCD to prove the reuse chain before a new full sector.
S57's representative smoke is orbit `121`, `cam4/ccd1` (orbit tag `o1`). It
creates only that CCD's A2v1 source links/sidecars/ePSFs/HDF5s.

```bash
export REPO=/pdo/users/tehan/TWIRL
export SOURCE_ROOT=/pdo/users/tehan/tglc-gpu-production
export A2ROOT=/pdo/users/tehan/tglc-gpu-production-A2v1
export FORK=/pdo/users/tehan/tess-gaia-light-curve-twirl
export TGLC_PY=/pdo/users/tehan/twirl-gpu-venv/bin/python
export QLP_PY=/sw/qlp-environment/.venv/bin/python
export LOG=$A2ROOT/twirl_logs/s57_A2v1_smoke_o121_cam4ccd1

export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=NUMEXPR_NUM_THREADS=1
export LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=$REPO/src:$FORK:${PYTHONPATH:-}

$QLP_PY $REPO/scripts/stage1_lightcurves/build_source_tic_overlays.py \
  --source-tglc-data-dir $SOURCE_ROOT \
  --output-tglc-data-dir $A2ROOT \
  --orbit 121 --ccd 4,1 --overlay-from-observations --link-sources \
  --apply --overwrite --summary-json $LOG/overlays.json

$QLP_PY $REPO/scripts/stage1_lightcurves/prefill_epsfs_for_empty_masks.py \
  --source-tglc-data-dir $SOURCE_ROOT \
  --output-tglc-data-dir $A2ROOT \
  --orbit 121 --ccd 4,1 --apply --nprocs 8 --summary-json $LOG/prefill.json

HDF5_USE_FILE_LOCKING=FALSE PYTHONPATH=/sw/qlp-environment/.venv/lib/python3.11/site-packages \
  $TGLC_PY $REPO/scripts/stage1_lightcurves/run_tglc_orbit_pipeline.py \
  --orbit 121 --sector 57 --orbit-tag o1 --ccd 4,1 \
  --tglc-data-dir $A2ROOT --log-dir $LOG --python-bin $TGLC_PY \
  --fork-path $FORK --ld-library-prefix /sw/python-versions/python-3.11.9/lib \
  --stages epsfs,lightcurves --epsfs-nprocs 2 --lightcurves-nprocs 16 \
  --max-parallel-ccd-jobs 1 --gpu-list 0 --run-label s57-a2v1-smoke
```

The smoke passes when the overlay/prefill summaries have no errors, the driver
reports `source=196` and `epsf=196`, and it writes HDF5s for requested targets.
Use a GPU reported free in pre-flight; the `0` above is an example only.

## Full-sector production

After the one-CCD smoke passes, use the generic launcher for both sector orbits
(S57 is the example below):

```bash
TWIRL_A2V1_GPU_LIST=0,1,2,3,6,7 \
TWIRL_A2V1_GPU_MAX_PARALLEL=4 \
TWIRL_A2V1_EPSFS_NPROCS=2 \
HDF5_USE_FILE_LOCKING=FALSE \
bash /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/run_a2v1_reproduction_pdo.sh \
  57 121:o1 122:o2
```

Set `TWIRL_A2V1_PREFILL_EMPTY_MASK_EPSFS=0` only when the prefill has already
completed successfully. On an ePSF OOM, leave completed files intact and
resume the affected CCD with fewer ePSF workers; do not delete the A2v1 tree.

### Sequential sector queue

For a prepared run of multiple sectors, use the gated queue rather than
starting sectors concurrently. It verifies both prepared input counts per orbit,
skips only product stages that already pass their validation gate, and stops at
the first failure:

```bash
TWIRL_A2V1_QUEUE_LABEL=a2v1-s58-s63 \
bash /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/run_a2v1_sector_queue_pdo.sh \
  /pdo/users/tehan/TWIRL/configs/a2v1_production_s58_s63.txt
```

The queue uses an HDF5-only gate that opens every nonzero HDF5 as well as
validating coverage before FITS production. It then validates the full
HDF5-plus-FITS schema before advancing. Its input gate accepts a complete
legacy ePSF tree for empty-mask reuse or a source-only tree that refits all
ePSFs, but fails on partial legacy ePSFs. It is a production chain, not a
status monitor; inspect its persistent queue log manually when needed. Run it
in detached tmux with verbose stage output redirected only to that log, not the
tmux TTY, because a full JSON validation report can otherwise fill the pane's
output buffer and block the queue.

## Final product and validation

HDF5 completion is an intermediate production checkpoint. Every completed A2v1
sector must then build its sector-level ADP/ADP015-only FITS products and pass
the full HDF5-plus-FITS validation report:

```bash
bash /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/run_a2v1_hlsp_pdo.sh \
  57 121 122

cd /pdo/users/tehan/TWIRL
LD_LIBRARY_PATH=/sw/python-versions/python-3.11.9/lib:/pdo/app/anaconda/anaconda2-4.4.0/lib:${LD_LIBRARY_PATH:-} \
HDF5_USE_FILE_LOCKING=FALSE /sw/qlp-environment/.venv/bin/python \
  scripts/stage1_lightcurves/validate_a2v1_product.py \
  --a2v1-root /pdo/users/tehan/tglc-gpu-production-A2v1 \
  --sector 57 --orbits 121 122 --allow-edge-warn-missing \
  --schema-only --fits-workers 8 \
  --summary-json /pdo/users/tehan/tglc-gpu-production-A2v1/twirl_logs/s57_A2v1_validation_full_schema.json
```

The HDF5 contract must report zero missing non-edge rows, zero zero-byte HDF5s,
and zero unreadable HDF5s when run with `--check-h5-open`. A standard complete
product must additionally report zero missing non-edge FITS TICs and zero bad
checked FITS files.
`edge_warn` exclusions remain visible in the report and are accepted only by
the explicit `--allow-edge-warn-missing` option.

The validator contract is strict: the expected target set must be nonempty,
every requested orbit must be represented, and `--allow-missing-fits` may
waive missing FITS only. It must never waive malformed or unreadable FITS.
Queue automation and downstream QA should consume the canonical report name
`s<sector>_A2v1_validation_full_schema.json` rather than a cwd-relative alias.

### QA promotion tiers

Production validation and science QA are separate gates:

- **Product validation** establishes expected HDF5/FITS coverage, openability,
  schema, and explicit edge exclusions.
- **Tier-0 integrity/benchmark QA** verifies compact-export/BLS coverage and
  configuration consistency, finite cadence behavior, aperture availability,
  and the fixed WD 1856 benchmark where applicable. It permits bounded
  diagnostics; it does not make a sector science-ready.
- **Tier-1 science QA** requires predeclared limits for RMS/MAD versus
  magnitude, cadence loss, quality flags, aperture outliers, fixed-injection
  preservation, and a genuinely independent extraction comparison. Only a
  later `full_a2v1_product` Tier-1 pass can promote a sector into a frozen
  survey release; the bounded active-pair scope can authorize enrichment only.

The Tier-1 implementation keeps two scopes explicit:

- `active_search_pair` scans the full compact population in
  `DET_FLUX_ADP_SML` and `DET_FLUX_ADP`, applies fixed magnitude/scatter,
  hash-bound QLP-quaternion/SPOC cadence and quality checks, aperture-ratio,
  fresh-injection preservation, and independent-extraction limits, and may set
  `enrichment_ready=true`. It always leaves `science_ready=false` because it
  does not audit all six A2v1 flux columns.
- `full_a2v1_product` is the later release-promotion scope. It must include all
  ADP/ADP015 `1x1`, `3x3`, and `5x5` channels and detector-stratified
  independent evidence before promotion is enabled.

The current `active_search_pair` implementation, evidence builders, and
contract names are S56-specific. They are not a generic sector-promotion
interface and do not yet define the observation-keyed multi-sector input
contract needed for seven-sector teacher training.

Run the bounded scope with
[`audit_a2v1_tier1_qa.py`](../scripts/stage1_lightcurves/audit_a2v1_tier1_qa.py)
and the locked
[`active-search configuration`](../configs/qa/a2v1_tier1_active_search_pair_v1.yaml).
The command requires the current Tier-0 v2 JSON, the exact compact ADP pair, a
cadence-reference table bound to the original QLP quaternion and SPOC-quality
authorities, the frozen four-shard (`2,000` unique-host) injection canary, and
an independent-extraction metrics table plus manifest. A prior TGLC production
tree is the same extraction family and does not satisfy the independent gate.
The configuration pins both the complete Tier-0 JSON and its BLS peak table;
all six nested Tier-0 gates and the WD 1856 benchmark must pass.

Build the SPOC intermediate, final cadence reference, and official TESSCut
WD 1856 evidence with
[`build_s56_spoc_quality_table.py`](../scripts/stage1_lightcurves/build_s56_spoc_quality_table.py),
[`build_a2v1_cadence_reference.py`](../scripts/stage1_lightcurves/build_a2v1_cadence_reference.py)
and
[`build_wd1856_tesscut_independent_extraction.py`](../scripts/stage1_lightcurves/build_wd1856_tesscut_independent_extraction.py).
The generic
[`build_a2v1_independent_extraction.py`](../scripts/stage1_lightcurves/build_a2v1_independent_extraction.py)
remains available for another genuinely independent, explicitly mapped
reference product. The official TESSCut comparison gates on signal presence
and ephemeris timing; raw depth and scatter ratios are diagnostics because
TESSCut aperture sums and decontaminated A2v1 flux do not share a dilution or
noise transfer function.
The locked configuration deliberately retains impossible authority, shard,
and independent-product hashes until the real S56 artifacts are built and
reviewed; placeholders must not be replaced with inferred values. The fixed-
injection metric is the fitted slope between `1 - transit_model` and the
negative injected-minus-original detrended flux; BLS recovery or teacher
scores do not substitute for this Stage 1 signal-preservation measurement.

Every required gate reports `pass`, `review`, or `fail`; missing evidence,
contract mismatch, or a stale Tier-0 report fails closed. Overall `pass`
requires every gate to pass, while any review result keeps the overall state at
`review`.

The Tier-1 target table is a downstream data contract. Candidate generation,
enrichment review, and teacher-set assembly must join it by exact
`(sector, TIC)` and retain only `tier1_target_qa_pass == True`; a sector-level
pass does not license failed or review-status targets. The published summary
includes the complete fixed-injection manifest and hashes all input and output
evidence.

Custom TWIRL/A2v1 FITS and compact products preserve the TGLC internal quality
flag; they do **not** already contain the sector-level SPOC and QLP flags that
the legacy QLP HLSP writer combines. For Tier-1 and every downstream
enrichment artifact, define a good cadence as
`internal_quality == 0 AND external_quality == 0`, where
`external_quality = spoc_quality | (qlp_quality << 30)`. Regenerate the real
BLS peaks with this hash-bound overlay, and require candidate metadata and
both real and injected teacher-native tensors to use the same cadence table
and manifest. The BLS summary, scoring table, and native HDF5 must retain the
table and manifest checksums; an older internal-only artifact is incompatible.
The S56 native-input contract is therefore v2, while candidate artifacts bind
that input through their separately versioned provenance envelope. A
checkpoint trained under native version 1 must fail compatibility checks
rather than be silently applied to changed quality channels and periodograms.
Encoder caches, final checkpoints, and checkpoint-transfer manifests must also
bind the exact native HDF5 and training-table hashes and live in a distinct
native-v2 namespace. Recovery scoring must verify that selected-checkpoint
manifest before allocating inference work, hash every candidate/native/model
input before and after scoring, and publish its score table and summary
atomically; the native-v1 recovery checkpoint default is retired.

The TWIRL I release manifest must name each accepted sector's product and QA
report checksums, the frozen sector cutoff, catalog/index version, and search
and injection contracts. A sector produced after that cutoff is not silently
added to the publication sample.

Record the run configuration, coverage counts, any GPU retry, and the final
validation JSON in `doc/twirl_progress_log.md`. Record Tier-0 and Tier-1
results separately; neither HDF5 completion nor a Tier-0 pass substitutes for
science-release promotion.
