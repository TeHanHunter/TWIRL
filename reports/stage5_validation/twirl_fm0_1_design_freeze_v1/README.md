# TWIRL-FM0.1 design freeze v1

This directory binds the frozen scientific design for the S56--S67
`TWIRL-FM0.1` proof-of-concept campaign.

The freeze covers:

- Gaia DR3 source identity and leakage-component grouping;
- the deterministic 80/10/10 train/development/sealed-test source split;
- the six `{raw-relative, ADP, ADP015} x {1x1, 3x3}` stored flux views;
- local and cross-sector elapsed time in nominal `200 s` cadence units;
- the `2,048`-cadence window contract;
- the `TWIRL-FM0.1.1`--`TWIRL-FM0.1.5` comparison ladder;
- the BLS, label, identifier, absolute-time, and detector-metadata exclusion
  boundary;
- diagnostic-only admission of S66--S67 after FM-specific gates; and
- one-H200, Stage-1-priority ORCD execution.

`freeze.json` records SHA-256 hashes for the design document and its
machine-readable configuration. It does not certify an input release, source
split, training job, checkpoint, or model result; those require separate
receipts.
