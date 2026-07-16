# S56 Peak-Ranker Handoff Status

> **Historical snapshot (2026-06-26).** Job states, ETAs, and next actions in
> this file are preserved evidence, not current instructions. Use
> [`doc/twirl_plan.md`](../../doc/twirl_plan.md) and the
> [progress log](../../doc/twirl_progress_log.md) for live status.

Last updated: `2026-06-26`.

## Live PDO Jobs

- `twirl-s56-20k-peak-training`: active current-schema restartable chunked `20k` injected BLS peak-table build on `pdogpu1`.
- `twirl-s56-peak-monitor-handoff`: active monitor on `pdogpu1`; waits for the peak table, then runs the post-BLS gate, trains/applies the injected-truth peak ranker, and builds a ranker-selected real-candidate review queue with LEO rendering skipped.
- The older one-shot build was stopped at `4,850 / 20,000` because it had started before the current cadence-diagnostic schema was deployed. It did not write a final table.
- The replacement chunked build started at `2026-06-26T18:28:28-04:00`, split the `20,000` injections into `80` chunks of `250`, and is running `2` chunks at a time with `8` workers per chunk.
- Latest bounded artifact check at `2026-06-26 23:32 EDT`: direct artifact poll saw `34 / 80` completed chunks, `36` chunk directories, `chunk_034` and `chunk_035` active, no final merged table, and no verification JSON yet.
- Process-health check showed the active chunk workers are CPU-bound rather than dead. Recent chunk-slot cadence is `~17.4 min`; with two chunks in parallel and `64` chunks remaining, the rough ETA to a merged table is `~9.3 hr` plus merge/verification/handoff time. The active chunk directories may stay empty until each chunk finishes because the chunk writer emits the CSV/summary at completion.
- A separate partial-current-schema smoke table from chunks `000`-`005` passed the downstream handoff: `1,500` injections, `60,000` peak rows, injected table verification passed, ranker trained, peak gate wrote recall summaries, a `6,000`-row real candidate sample was scored, and a `20`-row skip-LEO ranker-selected real queue passed verification.
- A group-relative peak-ranker feature smoke on that same partial table was rejected for production because it reduced top-1 injected recovery (`607/1500`) below both raw BLS rank (`629/1500`) and the conservative injected-truth ranker (`644/1500`). The deployed/default ranker remains the conservative absolute-feature model.
- A BLS failure-mode audit on the partial table found `43 / 1,500` high-observability injections where the signal was absent from the top-20 peaks. A targeted short-period branch smoke on those rows (`p_max=2 d`) recovered `15 / 43` in top-20 and `11 / 43` at top-1, so the next BLS iteration should test a short-period branch or per-period-range peak quota before assuming all not-in-top-N misses are below sensitivity.
- A read-only partial recall snapshot from completed chunks `000`-`013` saw `3,500` injections and `140,000` peak rows. Signal recovery was `1361/3500` at top-1 and `1556/3500` at top-20; `195` rows were ranker-fixable top-20-not-top-1 cases and `1944` rows were not in top-20. Top-20 recovery by Tmag bin was `<17: 603/895`, `17-18: 435/868`, `18-19: 321/873`, and `>19: 197/864`. Treat this as a partial progress check only; wait for the merged `20k` verifier before production decisions.
- A later read-only snapshot from completed chunks `000`-`015` saw `4,000` injections and `160,000` peak rows. Signal recovery was `1530/4000` at top-1, `1605/4000` at top-5, `1659/4000` at top-10, and `1742/4000` at top-20. The Tmag split remains steep (`<17: 671/1005`; `17-18: 484/993`; `18-19: 369/1004`; `>19: 218/998` at top-20), supporting the plan to gate search/ranker failures before human triage.
- The reproducible chunk-progress summarizer now writes the latest live partial report under `peak_training/chunk_progress_live/`. At `34/80` chunks it summarizes `8,500` injections and `340,000` peak rows, with top-1/top-5/top-10/top-20 truth recovery of `3240/8500`, `3407/8500`, `3519/8500`, and `3693/8500`. Top-20 recovery remains magnitude dominated (`Tmag < 17`: `1437/2181`; `Tmag > 19`: `450/2109`).
- The current PDO run is still in the table-building phase; do not start ranker production or branch comparison until the merged table and self-verification JSON exist.
- Local code now supports branch-provenance peak tables (`search_branch` plus BLS config columns) and branch-aware chunk merging with `--id-columns search_branch,injection_id`. Do not deploy the updated peak-table builder over the active `20k` PDO run; use it for the next standard+short branch comparison after the current table has merged and verified.
- Local code now also supports an explicit transit-window-overlap truth label for future injected peak tables. A peak must agree in period or supported harmonic and overlap the injected transit window before it is labeled as signal. The verifier is backward-compatible with the active PDO run, but future tables with overlap columns will be rejected if a signal label lacks transit-window overlap.
- Future real-candidate Stage 2 tables now preserve row-level BLS branch/config provenance (`bls_search_branch`, period grid, peak quota, cadence-cleaning settings) so a branch accepted from injection tests can be applied to real S56 targets and later audited after CSV export or ORCD staging.
- The real peak-table verifier remains backward-compatible with the current S56 table, but future branch-specific real tables can require this provenance with `TWIRL_REAL_REQUIRE_BLS_PROVENANCE=1` before ranker application.
- Local code also includes a read-only handoff auditor at `scripts/stage5_validation/audit_s56_peak_handoff_status.py`; run it after syncing or on PDO to write a compact JSON/Markdown gate summary before branch comparison or human-label scale-up.

Monitor log:

```text
/pdo/users/tehan/TWIRL/logs/s56_peak_monitor_handoff.log
```

The monitor reports chunk completion in this form:

```text
chunk progress: <done>/80 complete; active=<chunk names>
```

Future monitor runs also write a read-only partial recall report under:

```text
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training/chunk_progress_live/
```

Set `TWIRL_MONITOR_WRITE_CHUNK_PROGRESS=0` to disable that extra report.

Current peak-build log:

```text
/pdo/users/tehan/TWIRL/logs/s56_20k_peak_training_chunked_pdo.log
```

Stopped stale one-shot log:

```text
/pdo/users/tehan/TWIRL/logs/s56_20k_peak_training_pdo.log
```

## Expected Outputs

After the peak table finishes, the monitor should produce:

- injected peak-table verification: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training/s56_20k_injection_bls_peaks_chunked_verification.json`
- ranker input verification: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_pdo/peak_table_verification.json`
- peak gate: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_gate_pdo/`
- trained ranker: `reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_pdo/`
- ranker-selected real ephemerides: `reports/stage5_validation/s56_ranker_selected_real_candidates_pdo/`
- real review queue without PDFs: `reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo/`
- queue preflight verification: `reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo/verification.json`
- compact ranker-selection audit: `reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo/ranker_selection_summary/summary.md`

Successful smoke outputs live under:

```text
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_partial6_smoke/
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_partial6_smoke/
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_gate_partial6_smoke/
reports/stage5_validation/s56_ranker_selected_real_candidates_partial6_smoke/
reports/stage5_validation/s56_ranker_selected_real_review_queue_partial6_smoke/
```

The current conservative rerun of the same partial smoke is:

```text
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_ranker_partial6_default2_smoke/
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_gate_partial6_default2_smoke/
reports/stage5_validation/s56_ranker_selected_real_candidates_partial6_default2_smoke/
reports/stage5_validation/s56_ranker_selected_real_review_queue_partial6_default2_smoke/
```

The partial BLS failure-mode audit and short-period branch smoke are:

```text
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_gate_partial6_default2_smoke/failure_modes/
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_branch_smokes/short_pmax2_highobs_misses/
```

The real S56 BLS peak table is staged on PDO as CSV because the QLP runtime has
neither `pyarrow` nor `fastparquet`. The staged table is verified as multi-peak
and multi-aperture before ranker application:

```text
data_local/stage2/bls_first_pass_v2/sector_0056/candidates.csv
1,144,023 rows; 19,068 TICs; 20 positive peak ranks; 3 apertures
```

Verification artifact:

```text
reports/stage5_validation/s56_ranker_selected_real_candidates_pdo/real_peak_table_verification.json
```

The injected peak-table handoff now also requires the current cadence-cleaning
diagnostic columns before training the peak ranker. A tiny current-code PDO
schema smoke passed with:

```text
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training/current_code_schema_smoke_verification.json
```

The older `smoke_5_injection_bls_peaks.csv` failed this stricter check, as
expected, because it was produced before the full cadence-diagnostic columns
were added. If the active long one-shot table fails the same check, use the
restartable chunked runner rather than training on the older schema.

The monitor writes the search-gate outputs before attempting ranker/review
handoff. If strict schema verification fails, the recall@K gate remains useful,
but model training should wait for a chunked/current-schema peak table.

Both PDO peak-table runners now self-verify the output table after a long build.
The chunked runner writes:

```text
reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training/s56_20k_injection_bls_peaks_chunked_verification.json
```

This is the fallback artifact to require before training if the older one-shot
table fails the strict injected peak-table verifier.

## Next Checks

```bash
ssh pdogpu1 'cd /pdo/users/tehan/TWIRL && tail -n 80 logs/s56_20k_peak_training_chunked_pdo.log'
ssh pdogpu1 'cd /pdo/users/tehan/TWIRL && tail -n 80 logs/s56_peak_monitor_handoff.log'
PYTHONPATH=src python scripts/stage5_validation/summarize_injection_peak_chunk_progress.py
PYTHONPATH=src python scripts/stage5_validation/audit_s56_peak_handoff_status.py
```

If `s56_ranker_selected_real_review_queue_pdo/verification.json` passes, inspect `ranker_selection_summary/summary.md` for selection-rank, aperture, period, Tmag, and duration coverage, then render WD-tuned LEO PDFs for the selected real-candidate rows before human triage.

After the current `20k` table is merged and verified, the next BLS-method smoke should be a branch-aware comparison:

```bash
bash scripts/stage5_validation/run_s56_short_branch_compare.sh
```

The runner audits the standard table, writes high-observability miss IDs,
builds a `short_pmax2` branch only on those IDs, runs a peak gate on that
branch, and writes a standard-vs-short comparison. The underlying manual steps
are:

The same runner can test an opt-in period-quota variant without changing the
standard table, for example:

```bash
TWIRL_SHORT_BRANCH_NAME=short_pmax2_quota \
TWIRL_SHORT_BRANCH_OUT_ROOT=reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/peak_training_branch_compare/short_pmax2_quota \
TWIRL_SHORT_BRANCH_PERIOD_BIN_EDGES=0.12,0.25,0.5,1.0,2.0 \
TWIRL_SHORT_BRANCH_MAX_PEAKS_PER_PERIOD_BIN=5 \
  bash scripts/stage5_validation/run_s56_short_branch_compare.sh
```

This checks whether period-range coverage rescues high-observability injected
signals that are absent from the standard top-20 list.

If the branch table already exists and you only want to compare it, set:

```bash
TWIRL_SHORT_BRANCH_SKIP_BUILD=1 \
TWIRL_SHORT_BRANCH_PEAK_TABLE=<short-branch-peaks.csv> \
  bash scripts/stage5_validation/run_s56_short_branch_compare.sh
```

The lower-level commands are:

```bash
python scripts/stage5_validation/build_injection_peak_training_table.py \
  --search-branch short_pmax2 \
  --p-max-cap-d 2.0 \
  --n-periods 100000 \
  --limit-keys-file <high-observability-miss-ids.txt> \
  --out-table <short-branch-peaks.csv> \
  --overwrite

python scripts/stage5_validation/merge_injection_peak_training_chunks.py \
  --chunk-root <branch-chunks> \
  --out-table <merged-standard-short.csv> \
  --id-columns search_branch,injection_id

python scripts/stage5_validation/compare_injection_bls_branches.py \
  --baseline-peak-table <standard-peaks.csv> \
  --branch-peak-table <short-branch-peaks.csv> \
  --baseline-label standard \
  --branch-label short_pmax2 \
  --branch-search-branch short_pmax2 \
  --out-dir <branch-comparison-out>
```

Use the comparison outputs, especially `branch_rescued_injections.csv`,
`branch_degraded_injections.csv`, and `branch_comparison_by_cell.csv`, to decide
whether the short-period branch is a production candidate or just a targeted
diagnostic. The comparison defaults to `--comparison-scope intersection`, which
is essential when the short branch is built only on high-observability misses:
missing branch rows for the rest of the `20k` injections are not counted as
degradations.

On ORCD, after the compact products are staged and the standard peak table
exists, the equivalent Slurm wrapper is:

```bash
sbatch scripts/orcd/slurm_s56_short_branch_compare_cpu.sbatch
```

The dedicated PDO runner for the second stage is:

```bash
bash scripts/stage5_validation/run_s56_ranker_selected_real_leo_pdo.sh
```

It consumes `s56_ranker_selected_real_candidates_pdo/selected_ephemerides.csv`, writes the LEO-backed sibling queue at `reports/stage5_validation/s56_ranker_selected_real_leo_queue_pdo/`, verifies referenced reports, and writes its own `ranker_selection_summary/`.

ORCD remains the preferred place for larger downstream training once the compact S56 exports are staged, but non-interactive ORCD commands require the user-authenticated control socket described in [ORCD guide](../../doc/orcd_h200_usage.md).
