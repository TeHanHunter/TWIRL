# Franklin's Teacher v3 training-data handoff

This is the short path to the frozen TWIRL **Teacher v3** data release for an
independent model experiment. It contains Sectors 56--62 only.

## Start here

You need **one CSV plus seven HDF5 files**:

```text
training_table_with_splits.csv
native/
  s56_teacher_v3_native_v2.h5
  s57_teacher_v3_native_v2.h5
  s58_teacher_v3_native_v2.h5
  s59_teacher_v3_native_v2.h5
  s60_teacher_v3_native_v2.h5
  s61_teacher_v3_native_v2.h5
  s62_teacher_v3_native_v2.h5
```

Do **not** download or use the sector FITS products, raw TGLC HDF5 trees,
source pickles, FFIs, or ePSFs. The seven native HDF5 files already contain
the two-aperture, quality-aware light-curve and BLS inputs needed for model
training.

The CSV has the label, fixed split, TIC, and exact HDF5 group for every
observation. Treat it as the authority for row identity and labels.

## PDO handoff location

Franklin has everything needed in this **PDO-only** bundle:

```text
/pdo/users/tehan/twirl_handoffs/franklin_teacher_v3_s56_s62/
```

It contains the CSV at the top level, the seven native files under `native/`,
and a checksum manifest. No ORCD account, source checkout, FITS products, or
raw TGLC files are required. The CSV has already been made portable: its
`native_h5_path` values point to this PDO bundle.

### Quick presence checks

On PDO, first enter the handoff directory:

```bash
handoff=/pdo/users/tehan/twirl_handoffs/franklin_teacher_v3_s56_s62
cd "$handoff"
sha256sum -c SHA256SUMS
```

## What the table means

The table contains `8,181` rows; `8,163` are active training rows. Use the
package helper below to apply the frozen target policy correctly. In
particular, use:

- `teacher_v3_training_include` to select active rows;
- `morphology_target_v1` and the helper-derived integer target, rather than
  trying to remap human labels by hand;
- `fixed_split` to keep the frozen `development` and `test` populations
  separate;
- `cv_fold` (`0`--`4`) for development cross-validation; fixed-test rows have
  `cv_fold = -1`;
- `native_h5_path` and `native_group_path` to open the correct observation.

The split is grouped by TIC. Do not put observations of the same TIC into
both training and evaluation, and do not tune on the fixed `test` population.

The frozen model targets are `planet_like`, `eclipse_contact`,
`smooth_variable`, and `other`. The granular source decisions remain in the
human-label/provenance columns, but use `morphology_target_v1` (or the
helper-derived integer target) for this release's four-class task. They are
human morphology/enrichment labels, not confirmations of planets or eclipsing
binaries. The injected rows are synthetic controls, not real discoveries.

Sector 63 is deliberately absent and remains reserved for prospective
evaluation.

## Minimal Python example

Run this from a TWIRL checkout with its normal dependencies installed. It
loads the frozen rows, holds out one development fold, and reads one native
light curve without touching FITS or raw production files.

```bash
python -m pip install -e .
```

```python
from pathlib import Path

import pandas as pd

from twirl.vetting.harmonic_dataset import prepare_harmonic_training_rows
from twirl.vetting.harmonic_inputs import read_native_light_curve

handoff = Path("/pdo/users/tehan/twirl_handoffs/franklin_teacher_v3_s56_s62")
table_path = handoff / "training_table_with_splits.csv"
native_dir = handoff / "native"

# Apply the frozen target policy and preserve the assigned TIC-grouped split.
source = pd.read_csv(table_path, low_memory=False)
rows = prepare_harmonic_training_rows(source, seed=560062)

held_fold = 0
development = rows[rows["fixed_split"].eq("development")]
train_rows = development[development["cv_fold"].ne(held_fold)].copy()
validation_rows = development[development["cv_fold"].eq(held_fold)].copy()
fixed_test_rows = rows[rows["fixed_split"].eq("test")].copy()

assert set(train_rows["tic"]).isdisjoint(set(validation_rows["tic"]))
assert set(train_rows["tic"]).isdisjoint(set(fixed_test_rows["tic"]))

example = train_rows.iloc[0]
lc = read_native_light_curve(
    Path(example["native_h5_path"]),
    group_path=example["native_group_path"],
)

print(example[["review_id", "tic", "sector", "morphology_target_v1"]])
print(lc.time.shape, lc.raw_flux_small.shape, lc.bls_power_small.shape)
```

Each native HDF5 observation supplies chronological time, cadence number,
orbit ID, effective quality, raw small/primary-aperture fluxes and errors,
ADP small/primary fluxes, and the four BLS profiles on a common period grid.
Real rows are under `targets/`; injected controls are under `injections/`.

## Rules for a comparable experiment

- Keep the supplied fixed test split sealed until the model and evaluation
  choices are frozen.
- Use complete TIC groups for every new split; never split repeated TICs by
  sector observation.
- Preserve `sector`, `tic`, `review_id`, `teacher_v3_observation_id`, and the
  source HDF5/group columns in model outputs so results remain auditable.
- Record whether the experiment includes injected rows and how uncertain or
  broad-preserve-only rows are handled.
- Do not use Teacher v3 output as automatic candidate promotion or as a
  discovery claim. It is a frozen enrichment release.

## More detail

- The release configuration is
  [`configs/models/teacher_v3_s56_s62_a2v1_current_adp.yaml`](../configs/models/teacher_v3_s56_s62_a2v1_current_adp.yaml).
- The data-freezing summary is
  [`native_preparation.summary.json`](../reports/stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/native_preparation/native_preparation.summary.json).
- The main project constraints are in the [project plan](twirl_plan.md).
