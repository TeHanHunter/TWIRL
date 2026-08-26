# TWIRL-FM0.2.1 new-Mac handoff — 2026-08-26

This is a migration snapshot for the completed development-only
`TWIRL-FM0.2.1` step-2000 canary. It records evidence and the next bounded
task; it does not create new scientific or compute authorization.

## Machine and task boundary

- Move only the FM work to the new Mac. Start a fresh local Codex task from a
  clean checkout of the branch below so local tests run on the new Mac.
- Keep A2v1 white-dwarf light-curve production on the original Mac under the
  existing Codex task named `WD LC production`. It may be observed or
  controlled through Codex remote access, but this FM task must not poll,
  modify, or assume ownership of that production workflow.
- The Git handoff carries the durable FM state. Chat visibility or remote
  access does not change which Mac executes commands, so confirm the new task
  is attached to the new-Mac checkout before running local tests.

Read [`AGENTS.md`](../../../AGENTS.md), the
[documentation guide](../../../doc/README.md), and the
[project plan](../../../doc/twirl_plan.md) before acting. The project plan
remains the sole forward-looking authority.

## Git continuation

- remote: `git@github.com:TeHanHunter/TWIRL.git`
- branch: `codex/fm0-2-new-mac-handoff`
- handoff base commit: `8e958a7a` with subject
  `Checkpoint FM0.2 step-2000 handoff`
- handoff tip: the remote branch tip with subject
  `Finalize FM0.2 new-Mac handoff`
- immutable training revision:
  `ddf442aafb8f62966e549e2287abad3474dd556a`
- immutable evaluator revision:
  `83816d07975eebe3825d76dfe7096d22b70376f5`

The handoff commit packages compact reports, figures, receipts, plotting
scripts, and documentation. It does not alter or replace the immutable
training and evaluator revisions above.

If the new Mac has a clean checkout:

```bash
git fetch origin
git switch codex/fm0-2-new-mac-handoff
git pull --ff-only
git log -1 --oneline
```

The earlier rsync copy at `/Users/tehan/Projects/TWIRL` was intentionally dirty
and contains untracked files that this commit now tracks. Do not run `git
clean`, `git reset --hard`, or pull over that copy blindly. Preserve it as a
fallback and make a clean sibling clone instead:

```bash
cd /Users/tehan/Projects
mv TWIRL TWIRL_rsync_snapshot_20260826
git clone --branch codex/fm0-2-new-mac-handoff git@github.com:TeHanHunter/TWIRL.git TWIRL
cd TWIRL
git status --short --branch
```

The scoped Git handoff deliberately excludes unrelated dirty Stage-1,
Stage-2, collaborator, and generated-output work from the old Mac. Keep the
rsync snapshot until those changes receive their own reviewed checkpoint. It
does include the already-referenced frozen A2v1 lock file needed to trace the
FM input-product contract; that static provenance file does not transfer
ownership of the live A2v1 campaign.

Open a fresh local Codex task from this checkout and begin with:

```text
Read AGENTS.md, doc/README.md, doc/twirl_plan.md, and
reports/stage5_validation/
twirl_fm0_2_s56_s64_step2000_evaluation_v1/HANDOFF.md. Continue only the
bounded FM0.2 exact same-seed step-0 CPU representation-evaluation task.
First review and test the checksum-bound implementation locally. Do not touch
A2v1 production, submit more training, open the sealed test, run development
event retention, or apply the formal step-2000 gate.
```

## Result and interpretation

The authorized seed-`560067`, S56--S64, TCN canary stopped exactly at optimizer
step `2000`.

- training job `21303060`: `COMPLETED/0:0` in `1:30:23`, one H200, four CPUs,
  `32 GiB`, BF16;
- strict post-validation job `21327505`: `COMPLETED/0:0`, CPU-only;
- representation evaluator job `21327655`: `COMPLETED/0:0` in `3m19s`, four
  CPUs, `16 GiB`, no GPU;
- compact loss extraction job `21331169`: `COMPLETED/0:0`, CPU-only;
- sealed-test access count: `0` throughout.

Key immutable hashes:

- step-2000 checkpoint:
  `976b5053c857c38b9fbf7a35c9d0605f0023318b2c9bd37a88995992e5aa7bd2`
- step-2000 representation evaluation:
  `3c802741744999fe553e71bc2dccfbef6c309a29a28c30dfaac693a798bb3718`
- post-validation receipt:
  `80d678c4621a22258c5d0e8f0be6f7e08ffa5bf93841c315abf42e9bb006b110`
- evaluator receipt:
  `616ef9a100b7b7a3a1923f81ac19a272cab8b6f4657a4a58c2815688ee6d1191`
- loss-history extract:
  `7d9eb81700f0185a725f814dc693ef8969eeb40953c4c42266b6d803d5909dfe`

The direct VICReg objective repaired the targeted representation pathology:

- `z_window` effective rank: `6.45 -> 19.43 -> 39.40` at steps
  `500 -> 1000 -> 2000`;
- `z_window` constant dimensions: `0`;
- paired-minus-unrelated 95% lower bound at step 2000: `0.05507`;
- trained-minus-random 95% lower bound at step 2000: `0.03311`.

This is still a mixed result, not a canary pass:

- development masked Huber is `0.004869893` versus the frozen zero-predictor
  value `0.004864579`, or `0.109%` worse;
- query-sector-excluded retrieval was not demonstrated and regressed from a
  clustered mean of `0.02778` at step 1000 to `0.02083` at step 2000;
- `h_window` retains one constant dimension;
- the formal scientific go/no-go was not applied;
- development event retention was not run.

The complete interpretation, tables, plots, and receipt inventory are in the
[step-2000 evaluation report](README.md).

## Exact next FM task

Do not submit another training job. The next task is to review, test, and then
execute a dedicated checksum-bound CPU-only representation evaluation of the
exact seed-`560067` initialization checkpoint:

```text
checkpoint_step_00000000.pt
SHA-256 92463070381486f0c6053c190d62da8d0c5c0d31be8072e2bdcd677329ac792c
```

The preserved ORCD run directory is:

```text
/orcd/data/mki_aryeh/001/twirl/reports/stage5_validation/
twirl_fm0_2_s56_s64_objective_canary/ddf442aafb8f/
model_runs/fm0_2_1_canary/seed560067
```

The frozen config requires representation evaluations at
`[0, 500, 1000, 2000]`, but the current controller and Slurm wrapper authorize
only `500|1000|2000`:

- `scripts/orcd/run_twirl_fm0_2_canary_orcd.sh`
- `scripts/orcd/slurm_twirl_fm0_2_representation_health_cpu.sbatch`
- `tests/test_twirl_fm0_2_contract.py`

Do not simply bypass those checks. Review a dedicated step-0 path that:

1. verifies the exact step-0 hash against the passed step-2000 post-validation
   receipt, which already binds the step-0 checkpoint;
2. reuses the identical development population, observation-sector authority,
   evaluator definitions, and evaluator revision used at steps 500--2000;
3. requests exactly four CPUs and `16 GiB`, no GPU, and excludes `node4900`;
4. refuses overwrite and emits a checksum-bound receipt;
5. records zero sealed-test access;
6. labels the exact initialization separately from the matched-random encoder
   control;
7. compares steps `0/500/1000/2000` only after the new receipt passes.

The implementation should receive focused tests and repository review before
submission. ORCD access must use a new user-opened control socket on the new
Mac; never copy an SSH socket or let Codex initiate Duo, password, or
keyboard-interactive authentication.

## Closed actions

Unless the user explicitly makes a new scoped decision:

- no H200 continuation or training beyond step 2000;
- no `FM0.2.2`, second seed, Conformer, or 20,000-step extension;
- no sealed-test access;
- no development event-retention run;
- do not apply the formal step-2000 gate;
- no production, science, or foundation-model claim;
- Stage-1 A2v1 work retains compute priority.

## New-Mac environment check

From the clean clone:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e '.[dev,representation]'
make test-fast
make check-docs
```

The last old-Mac validation completed with `1,157` tests passed and `33`
optional skips. Documentation and whitespace checks also passed.

Git carries only compact evidence. The step-0 and step-2000 checkpoint files
remain on ORCD and must not be copied into the repository.
