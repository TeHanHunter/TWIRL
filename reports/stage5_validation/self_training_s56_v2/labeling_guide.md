# S56 Human Labeling Guide

Fill the `label`, `labeler`, and `notes` columns. Keep `label_source=human` for rows you inspect by eye.

Suggested label vocabulary:

- `planet_like`: WD 1856-like or otherwise transit-like and worth preserving.
- `eclipsing_binary_or_pceb`: likely stellar/PCEB eclipse or Roche-limit contaminant.
- `stellar_variability`: variable-star or coherent non-transit behavior.
- `instrumental_or_systematic`: scattered light, aperture, cadence, or detrending artifact.
- `centroid_contaminant`: signal likely belongs to a nearby source.
- `uncertain`: keep in the human-review pool; do not use as a strong class.
- `skip`: unusable row or duplicate inspection target.

The `source_bucket` column records why each row was selected for this template; it is not a training label.

Browser vetting app:

```bash
PYTHONPATH=src .venv/bin/python scripts/stage5_validation/run_lightcurve_vetting_app.py
```

On PDO, run the same command inside the TWIRL checkout, then from the local
machine tunnel the port:

```bash
ssh -L 5000:localhost:5000 pdogpu6
```

Open `http://localhost:5000/`. The app reads this template and displays the
pre-rendered LEO-Vetter report first when one exists. The TWIRL HLSP light
curve is collapsed by default and is only a fallback/debug context panel. It writes labels to
`human_labels_vetted.csv` so the template remains unchanged. Existing LEO
reports currently cover the prior top-candidate runs, not every row in this
300-object template; rows without a report show the TWIRL light-curve panel
only until the full template is run through LEO.

Bucket counts:

- `planet_candidate_top`: 66
- `sub_roche_pceb_suspect`: 54
- `pceb_grid_ceiling`: 45
- `random_background`: 36
- `planet_centroid_review`: 30
- `high_sde_all_classes`: 30
- `faint_planet_candidate`: 21
- `planet_centroid_pass`: 17
- `wd1856_benchmark`: 1
