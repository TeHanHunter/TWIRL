# CLAUDE.md

## Project

TWIRL (TESS White dwarf Investigation of Remnant pLanets) searches for transiting / occulting objects around white dwarfs using 200 s TESS FFIs (Sector ≥ 56). Benchmark target: **WD 1856+534**, Sector 56 orbits 119/120, cam4/ccd1.

## Read First

Before any task, consult `doc/`:

- `twirl_plan.md` — forward plan; `## Current Status Snapshot` and the active sprint section show where things stand
- `twirl_progress_log.md` — dated history; active subsections end with `**Next:**`
- `ideas.md` — open questions / referee risks; starts with `## Critical Before Next Stage`
- `mit_tglc_usage_guide.md` — PDO operational details (envs, per-orbit pipeline, recovery audits, Sector < 67 flag handling)
- `local_data.md` — local data paths (not in git)

No formal test suite. Validation is manual via script outputs and checkpoint files.

## Key Constants

- `FIRST_200S_SECTOR = 56`
- High-confidence WD cut: `Pwd > 0.75` (~359,073 targets)
- Authoritative target identifier: **Gaia DR3 `source_id`** (TIC is secondary)

## Layout

`src/twirl/` is the importable package; `scripts/stageN_*/` are thin argparse CLIs over it. Scripts add `src/` to `sys.path` rather than installing the package editably. `data_local/` is gitignored. `doc/twirl_progress_log.md` is append-only.

## Data Policy

- Never commit anything under `data_local/`.
- No duplicate products with different filenames — use versioned names (`_v0`, `_v1`).
- Every output carries a JSON sidecar with seed-file path, size, mod-time, build version, SHA256, and key assumptions.

## Working Conventions

- Plan first, then edit. Prefer minimal, local changes; don't refactor files unrelated to the task.
- Production code goes in `src/twirl/` and `scripts/`, not notebooks.
- `astropy.table.Table` for structured catalog data; FITS for catalogs, HDF5 for light curves, JSON for metadata.
- Plotting: use `src/twirl/plotting/style.py` templates (`column` 3.4×2.65 in, `full_page` 7.1×4.1 in); no figure titles by default.
- Update `doc/twirl_progress_log.md` after completing a milestone.

## PDO Rules (any `pdo*` host)

- **Read-only outside `/pdo/users/tehan/`.** Never modify, delete, or move files under `/pdo/qlp-data/`, `/sw/`, or other users' dirs. If another user's file is needed in our tree, symlink from `/pdo/users/tehan/` to the upstream path — do not copy or edit upstream.
- **Cap BLAS/OpenMP threads before any multiprocessing pool** (`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=NUMEXPR_NUM_THREADS=1`). Without this, each worker spawns BLAS threads equal to the full core count (128 on pdogpu1) → ~8× oversubscription, 10–100× slowdown. QLP detrend/HLSP wrappers already bake this in.

## End-of-Day Wrap

When the user says **"wrap for the day"**:

1. Update `doc/twirl_plan.md` `## Current Status Snapshot` to current state.
2. Append dated notes to `doc/twirl_progress_log.md` for this session's work.
3. Commit all tracked changes with a concise message.
4. Push to remote.
5. Leave the next concrete step in the progress log as a `**Next:**` line.

## Stop and Ask Before

- Selecting a new benchmark sector
- Broadening science claims beyond the WD 1856-like regime
- Merging content from `TWIRL_proposal`
- Committing to a follow-up instrument
- Changing repo scope
