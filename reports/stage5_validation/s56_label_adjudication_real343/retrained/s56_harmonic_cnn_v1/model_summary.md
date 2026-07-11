# S56 Harmonic CNN v1 Model Summary

## Run

- Model contract: `s56_harmonic_cnn_v1` with native input contract
  `s56_adp_raw_pair_v1`.
- ORCD training job: `17655871`, one H200, `02:11:41`, exit code `0:0`.
- Training rows: `2,152`; fixed TIC-grouped test rows: `430`.
- Native input assembly: `1,696` unique real targets plus `456` active
  injected rows, with raw flux/error, two ADP apertures, and two-aperture BLS
  periodograms present.
- The fixed test set was opened once, after development-fold architecture
  selection.
- Checkpoints remain on ORCD under
  `/orcd/data/mki_aryeh/001/twirl/code/TWIRL/reports/stage5_validation/s56_label_adjudication_real343/retrained/s56_harmonic_cnn_v1/`.

## Selected Model

The five-fold development selector chose `shape_plus_periodogram_bls`. This
profile uses the seven native, unbinned harmonic folds and local primary and
secondary windows, together with the two-aperture BLS periodogram and scalar
BLS branch. The full model with raw chronological light curves ranked second.
Thus raw chronology was tested rather than discarded a priori, but it did not
improve development performance with the current sample.

| Development profile | Macro F1 | Balanced accuracy | Real Planet recall | EB recall | Variable recall | Other recall | ECE |
|---|---:|---:|---:|---:|---:|---:|---:|
| Shape + periodogram/BLS | 0.807 | 0.786 | 0.800 | 0.667 | 0.756 | 0.990 | 0.016 |
| Full combined | 0.795 | 0.783 | 0.800 | 0.643 | 0.709 | 0.991 | 0.018 |
| Metadata only | 0.724 | 0.707 | 0.800 | 0.548 | 0.698 | 0.969 | 0.017 |
| Shape + raw chronology | 0.736 | 0.695 | 0.600 | 0.548 | 0.593 | 0.988 | 0.025 |
| Seven-harmonic shape | 0.717 | 0.680 | 0.500 | 0.524 | 0.500 | 0.988 | 0.026 |
| Single-period native fold | 0.696 | 0.677 | 0.100 | 0.571 | 0.605 | 0.987 | 0.024 |

The selected model exceeds metadata-only development macro F1 by `0.084`.

## Locked Test

The four-class morphology evaluation has `428` active rows. The other two
test rows are Broad isolated dips, which intentionally train only the preserve
and harmonic tasks.

| Evaluation subset | n | Accuracy | Macro F1 | Balanced accuracy | ECE |
|---|---:|---:|---:|---:|---:|
| All morphology rows | 428 | 0.909 | 0.757 | 0.750 | 0.048 |
| Real rows | 337 | 0.961 | 0.757 | 0.838 | 0.036 |
| Injected rows | 91 | 0.714 | 0.715 | 0.754 | 0.138 |

| Human morphology | Test support | Correct | Precision | Recall |
|---|---:|---:|---:|---:|
| Planet-like | 64 | 42 | 0.875 | 0.656 |
| Eclipse/contact | 11 | 7 | 0.500 | 0.636 |
| Smooth variable | 22 | 16 | 0.800 | 0.727 |
| Other | 331 | 324 | 0.936 | 0.979 |

Planet-like support is strongly source-imbalanced: `2` real and `62` injected.
Both real examples were recovered, whereas `40/62` injected human-positive
examples were recovered. The apparent real-Planet recall of `1.0` is therefore
not a stable scientific estimate.

The preserve head has balanced accuracy `0.889`; preserve-positive precision
is `0.889` and recall is `0.808`. Both held-out Broad isolated dips were
preserved, but `n=2` is audit evidence only.

The harmonic head has accuracy `0.909` against a `0.792` majority-factor
baseline, exceeding the configured gate by `0.117`. This aggregate result is
not broad harmonic competence: it recovered `58/61` `P` rows and `12/12` `2P`
rows, but `0/4` `P/2` rows, and the other factors had no held-out support.

## Interpretation

All software-defined smoke and promotion gates pass, including calibration,
non-degenerate four-class predictions, and the combined-model comparison.
Scientifically, this remains an active-learning teacher rather than a final
student-label generator because:

- only `12` unique real Planet-like rows exist in the full training table;
- only `2` real Planet-like rows reached the locked test set;
- injected Planet-like examples dominate the Planet test support;
- rare harmonic factors are not represented well enough to validate the
  harmonic head; and
- the raw chronology branch currently overfits relative to the selected
  shape-plus-BLS model.

Human morphology labels are the training targets. Injection truth and recovery
columns remain separate audit data and are not model inputs. Student
pseudo-labeling remains blocked until at least `50` unique real Planet-like
examples are vetted.

## Artifacts

- [Performance PNG](performance_figure/s56_harmonic_cnn_v1_performance.png)
- [Performance PDF](performance_figure/s56_harmonic_cnn_v1_performance.pdf)
- [Figure provenance](performance_figure/figure_summary.json)
- [Development ranking](development_profile_ranking.csv)
- [Locked-test predictions](fixed_test_predictions.csv)
- [Machine-readable summary](summary.json)
- [Injection truth/human audit](injection_truth_human_audit/summary.json)
