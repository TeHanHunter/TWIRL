# Franklin S56 Real-Candidate Vetting

This folder is the lightweight Git handoff for Franklin's S56 real-only vetting
queue. The full data package is produced on ORCD as a separate ZIP because the
rendered vet sheets are several GB and should not be committed to git.

## Download Package

Use the ORCD-built ZIP when it is ready:

```text
/orcd/data/mki_aryeh/001/twirl/code/TWIRL/reports/stage5_validation/s56_franklin_real5k_handoff_png.zip
```

That ZIP is self-contained. After unzipping it, Franklin should have:

- `franklin_review_queue_5k_real.csv`
- `franklin_vetting_app.py`
- `run_franklin_vetting.sh`
- `README_Franklin_vetting.md`
- `reference_examples/`
- `reference_examples.csv`
- `vet_sheets/*.png`
- `franklin_labels_vetted.csv`

## Start Labeling

From inside the unzipped package:

```bash
python3 franklin_vetting_app.py --check-only
./run_franklin_vetting.sh
```

Then open:

```text
http://127.0.0.1:5003/
```

Labels save immediately to `franklin_labels_vetted.csv`.

## Label Policy

- `planet_like`: preserve a compact transit-like signal.
- `wide_transit_like`: preserve a broad/long-duration transit-like signal.
- `eclipsing_binary_or_pceb`: EB/PCEB-like shape, odd/even mismatch, secondary
  eclipse, or stellar-companion morphology.
- `stellar_variability`: astrophysical variability or repeating non-transit
  structure.
- `instrumental_or_systematic`: window-edge leakage, cadence artifacts,
  extraction artifacts, or other non-astrophysical structure.
- `uncertain`: flat/no obvious useful event.
- `skip`: unusable row or app/sheet problem.

If the signal looks folded at the wrong harmonic, keep the preserve/reject
label and add a note such as `half period`, `refold at P/2`, or
`possible harmonic`.

## Current ORCD Jobs

- Render PNG-only fresh sheets: `17503641`
- Package self-contained ZIP after render success: `17503953`

The queue is real-only, excludes injected rows, and excludes prior real review
rows from the current teacher queues.
