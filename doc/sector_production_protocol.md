# TWIRL sector production protocol

Procedure for producing TGLC + detrend + HLSP outputs for one TESS sector
(`Sector >= 56`) using the GPU pipeline on PDO. Encodes the working state of
the codebase as of `2026-04-30` (after S56 GPU production validated against the
frozen `tglc-deep-catalogs/hlsp_s0056/` baseline).

## Time budget per sector

| | clean | real-world |
|---|---|---|
| TGLC (catalogs→cutouts→epsfs→LCs) for both orbits | `~18 h` | `~25–30 h` |
| Post-LC chain (cadence-fix + detrend + HLSP + audit) | `~1 h` | `~2 h` |
| **Per sector total** | **`~19 h`** | **`~25–30 h`** |
| **40 sectors, single pdogpu6** | **`~32 days`** | **`~50 days`** |
| **40 sectors, 2 GPU nodes** | **`~16 days`** | **`~25 days`** |
| **40 sectors, 4 GPU nodes** | **`~8 days`** | **`~12 days`** |

The "real-world" column accounts for occasional host reboots, NFS flakes, and
single-file corruptions that the auto-resume + corrupt-file-rm pattern handles
but still cost wall time.

## Output layout (one tree, all sectors)

```
/pdo/users/tehan/tglc-gpu-production/
    orbit-<NNN>/ffi/
        catalogs/                         # Gaia + TIC ecsv per CCD
        cam{1..4}/ccd{1..4}/
            ffi/                          # symlinks into /pdo/qlp-data/tica-delivery/...
            source/                       # cutout pickles
            epsf/                         # ePSF .npy
            LC/                           # per-TIC h5 (WD-only)
        run/                              # qlp.cfg + cam*_quat.txt + cam*ccd*_qflag.txt
    hlsp_s00<NN>/                         # sharded HLSP FITS, sector level
    twirl_logs/
        s<NN>-gpu-rerun/                  # per-CCD stage logs + summary.json
        post_lc_chain/                    # post-LC chain logs
    qc_s<NN>_gpu/                         # QC plot PNGs per sector
```

The frozen `tglc-deep-catalogs/hlsp_s0056/` tree stays untouched as the
benchmark-of-record.

## Pre-flight (5 min)

1. **Sector data presence:**
   ```bash
   ls -d /pdo/qlp-data/tica-delivery/s00${SECTOR}/cam*-ccd*/ | wc -l   # expect 16
   ls /pdo/qlp-data/quats/cam?_quat.txt | wc -l                       # expect 4
   ls /pdo/qlp-data/spocflags/spocffiflag_s${SECTOR}_cam*.txt | wc -l  # only Sector < 67
   ```
   If any are missing, ping QLP ops before starting.

2. **Cadence-range mapping for detrend/HLSP wrappers:**
   ```bash
   /pdo/users/tehan/twirl-gpu-venv/bin/python \
     /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/build_sector_orbit_json.py \
     --sector ${SECTOR} \
     --base-orbit ${SECTOR}:${FIRST_ORBIT_NUMBER} \
     --output /pdo/users/tehan/tglc-gpu-production/sector_orbit.json
   export TWIRL_SECTOR_ORBIT_JSON=/pdo/users/tehan/tglc-gpu-production/sector_orbit.json
   ```
   Add `--orbit-map ${SECTOR}:${OW}:${ABS}` for sectors with non-contiguous
   orbits (e.g. S97/S98 have 4 orbits each).

3. **Disk space:** budget `~350 GB` free under `/pdo/users/tehan/` per sector.
   ```bash
   df -BG /pdo/users/tehan
   ```

4. **GPU node:** `pdogpu6` is the default (8 GPUs, `cupy-cuda11x` venv ready).
   `pdogpu1` has no GPU drivers; pdogpu2–5/7–8 are usually busy with QLP BLS.

## Stage 1: TGLC light curves (per sector)

Run each orbit in its own tmux. The launcher auto-retries on host reboots /
driver crashes (`until` loop, max 10 attempts) and auto-tunes GPU concurrency
via `pick_gpu_config.sh`.

```bash
# Orbit 1 (with Gaia query)
tmux new-session -d -s twirl-s${SECTOR}-o1 \
  'bash /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/run_sector_gpu_production.sh \
       ${SECTOR} ${ORBIT_1} o1; \
   echo DONE_O1; sleep 86400'

# Orbit 2 (sibling-catalog reuse — saves the dense Gaia query)
tmux new-session -d -s twirl-s${SECTOR}-o2 \
  'bash /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/run_sector_gpu_production.sh \
       ${SECTOR} ${ORBIT_2} o2 --reuse-catalogs-from-orbit ${ORBIT_1}; \
   echo DONE_O2; sleep 86400'
```

The driver auto-skips CCDs whose `summary.json` already exists, so a relaunch
is cheap. Per-stage `--replace` is OFF, so partial outputs are reused.

**Monitor (every 20–30 min initially, hourly once stable):**

```bash
ssh pdogpu6 "tmux capture-pane -t twirl-s${SECTOR}-o1 -p | tail -20"
ssh pdogpu6 "for d in /pdo/users/tehan/tglc-gpu-production/orbit-${ORBIT_1}/ffi/cam*/ccd*/; \
  do printf '%-60s ' \$d; \
     for sub in source epsf LC; do n=\$(ls \$d/\$sub/ 2>/dev/null | wc -l); \
       printf '%s=%-5s ' \$sub \$n; done; echo; \
  done"
```

## Failure-mode playbook

| symptom | cause | fix |
|---|---|---|
| `cupy.cuda.memory.OutOfMemoryError` | GPU memory pressure under 4-parallel × 8-epsf | already mitigated — current default is 4-parallel × 4-epsf via `pick_gpu_config.sh`. If still happens, drop `epsfs-nprocs` to 2 |
| `_pickle.UnpicklingError: invalid load key '\x00'` | NFS write left a corrupt source pickle | locate the bad pkl with `python -c "import pickle; pickle.load(open(...))"` over each `source_*.pkl`; `rm` it; the until-loop retry re-cuts only that one |
| `np.load(...) ValueError` on epsf `.npy` | NFS write left a corrupt ePSF | same recipe with `np.load`; `rm` the bad `.npy`; ePSF resume re-fits only that cutout |
| `OSError: [Errno 5] Input/output error` | NFS-side flake mid-write | until-loop will retry; on retry the partial file may still be corrupt — check + `rm` if so |
| `failed to connect to server` (tmux) | host rebooted | `ssh pdogpu6 "uptime"` to confirm; relaunch the same tmux command, auto-resume picks up |
| `Connection closed by UNKNOWN port 65535` | ssh gateway transient | back off 30+ min; the tmux on pdogpu6 keeps running independently |
| Driver "completed cleanly" but a CCD has incomplete LCs | auto-resume threshold bug — already fixed (`summary.json` is now the marker) | should not recur; if it does, manually relaunch with `--no-auto-resume` for the affected CCD |

## Stage 2: post-LC chain (per sector)

```bash
tmux new-session -d -s twirl-s${SECTOR}-postlc \
  'bash /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/run_s56_post_lc_chain.sh; \
   echo DONE_POSTLC; sleep 86400'
```

For sectors other than 56, copy `run_s56_post_lc_chain.sh` to
`run_s${SECTOR}_post_lc_chain.sh` and parameterize the sector + orbit numbers.
Critical bits to preserve:

- All stage outputs use `>> log 2>&1`, NOT `tee` (`tee` deadlocks the
  multiprocessing pool when the parent's pipe blocks).
- Cadence-fix runs serial with `HDF5_USE_FILE_LOCKING=FALSE`. Parallel pool
  also deadlocks on NFS h5py writes.
- `detrend_wrapper.py` and `hlsp_wrapper.py` honor `TWIRL_SECTOR_ORBIT_JSON`,
  so the pre-flight JSON is sufficient for sectors beyond S56.
- HLSP needs `--basedir /pdo/users/tehan/tglc-gpu-production/` and
  `-o /pdo/users/tehan/tglc-gpu-production/hlsp_s00${SECTOR}/`.

## Stage 3: validation

After post-LC chain finishes:

1. **Recovery rate:**
   ```bash
   find /pdo/users/tehan/tglc-gpu-production/hlsp_s00${SECTOR} -name '*.fits' | wc -l
   # Expect close to (sum of LC counts across both orbits' CCDs)
   ```

2. **Per-CCD failure summary:**
   ```bash
   /pdo/users/tehan/twirl-gpu-venv/bin/python \
     /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/summary_of_failures.py \
     --run-log /pdo/users/tehan/tglc-gpu-production/twirl_logs/s${SECTOR}-gpu-rerun \
     --tglc-data-dir /pdo/users/tehan/tglc-gpu-production \
     --orbit ${ORBIT_1} ${ORBIT_2}
   # Expect "32/32 CCDs DONE"
   ```

3. **QC plots (random sample for visual review):**
   ```bash
   /pdo/users/tehan/twirl-gpu-venv/bin/python \
     /pdo/users/tehan/TWIRL/scripts/stage1_lightcurves/qc_plot_hlsp_sample.py \
     --hlsp-root /pdo/users/tehan/tglc-gpu-production/hlsp_s00${SECTOR} \
     --output-dir /pdo/users/tehan/tglc-gpu-production/qc_s${SECTOR}_gpu \
     --n-per-bin 64
   # Then scp the 3 PNGs (faint/mid/bright) down for visual review.
   ```

4. **Headline-target parity check** (where applicable): for sectors that have
   a known reference target (e.g. WD 1856 in S56), load both old and new
   HLSP FITS and compare per-aperture MAD-RMS. Target: within a few percent.

## Documentation hooks

After each sector:

- Append to `doc/twirl_progress_log.md` §1.4: per-CCD wall-time totals,
  HLSP recovery rate, any failure modes encountered, fix applied.
- Update `doc/twirl_plan.md` `## Current Status Snapshot` with the latest
  sector and aggregate count.
- Spot any new failure mode? Add it to the playbook above.
