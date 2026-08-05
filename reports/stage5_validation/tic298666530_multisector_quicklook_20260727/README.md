# TIC 298666530 multi-sector quick-look

## Result

- The target is Gaia DR3 `2154363086096593536`, a high-confidence GF21 white
  dwarf candidate (`Pwd = 0.903566`, `T = 18.326`).
- The 1x1-ADP wide-window refinement gives
  `P = 0.573280500 d`, `T0 = 2460421.88834557 BJD`, and a `12 min` duration.
  Its robust-MAD-weighted combined depth is `2.81%` with depth S/N `19.67`.
- A narrower `0.5733 +/- 0.0005 d` cross-check agrees between apertures:
  `P = 0.573280500 d` in 1x1 ADP and `P = 0.573280350 d` in 3x3 ADP.
  Excluding anomalous S78 changes the 1x1 result by only `0.0000004 d`.
- Independent sector fits find a nearby, depth-S/N >= 5 signal in `14/16`
  sectors in 1x1 ADP, `13/16` in 3x3 ADP, and `12/16` in both. The wide
  3x3-only scan is instead dominated by a `0.5707684 d` systematic peak; it is
  not used as the consensus ephemeris.
- Sector depths and aperture ratios vary dramatically. S78 is the clearest
  outlier, and S77/S83/S84/S86 show weaker or aperture-dependent behavior.
  These products establish recurrent photometric structure near `0.5733 d`,
  not a common physical transit depth or an astrophysical confirmation.

## Pipeline audit

- Yes, TWIRL had encountered the target. The broad S56 search surfaced it but
  selected long `7.3--7.7 d` aliases.
- A target-seeded S56 current-ADP vet recovered `0.573353453 d` in 1x1 ADP
  (`SDE = 19.0`), while 3x3 ADP selected a different alias.
- More importantly, the unseeded S57 enrichment queue independently ranked it
  at `0.573296787 d`; both apertures agreed, with SDE values `107.9` and
  `47.4`. It was not subsequently present in the frozen human-label corpus.
  Exact source rows are preserved in
  [pipeline_discovery_audit.csv](pipeline_discovery_audit.csv).

## Data coverage and provenance

- The master observation table contains `21` planned sectors. The `16`
  delivered/observed sectors through S86 are S56--S60, S73, S75, S77--S80,
  and S82--S86. S117--S121 are forecast-only and were not present in current
  PDO TICA delivery.
- The five existing S56--S60 A2v1 FITS were combined with `11` new,
  target-only sector FITS. New processing used symlinks to the precut source
  pickles; all overlays, ePSFs, HDF5s, FITS, and logs were written only under
  `/pdo/users/tehan/tglc-tic298666530-quicklook-20260727/`.
- Aggregate PDO validation passed `22` ePSFs, `21` HDF5s, and all `11` new
  A2v1 FITS schemas. S79 is detector-edge flagged: orbit 165 did not retain a
  target HDF5 even after a serial retry, so its FITS deliberately contains
  valid orbit 166 only. See
  [the PDO validation report](pdo_product_validation.json) and
  [coverage inventory](coverage_inventory.csv).
- S56-style one-column TIC-ID indexes now accompany every new sector product.
  Their combined manifest, scope warning, and hash validation live under
  `/pdo/users/tehan/tglc-tic298666530-quicklook-20260727/tic_id_indexes/`.
  These indexes correctly contain only TIC 298666530 because the quick-look
  extraction was target-only; they are not full-sector target inventories.
- Exact local FITS paths and SHA-256 hashes are in
  [input_provenance.csv](input_provenance.csv). Cleaning matches the
  transparent per-sector BLS baseline: finite data, internal `QUALITY == 0`,
  per-sector median normalization, and upper-tail clipping. The combined
  refinement additionally uses per-sector robust-MAD weights.

## Outputs

- [Summary figure PNG](tic298666530_multisector_summary.png) and
  [PDF](tic298666530_multisector_summary.pdf)
- [Per-sector folds PNG](tic298666530_sector_folds.png) and
  [PDF](tic298666530_sector_folds.pdf)
- [Sector metrics](sector_metrics.csv),
  [wide-window periodogram](refined_periodogram.csv), and
  [anchored-period cross-check](anchored_period_crosscheck.json)
- [Summary JSON](summary.json) and
  [leave-one-sector-out periods](leave_one_sector_out_periods.csv)

This is an exploratory candidate diagnostic. It does not replace the frozen
multi-sector search/false-alarm contract and is not an astrophysical
confirmation.
