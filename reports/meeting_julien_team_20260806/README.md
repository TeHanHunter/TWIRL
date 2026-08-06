# Julien-team meeting package: S56–S62 candidates and model status

Prepared on `2026-08-06` from the frozen S56–S62 human morphology review and the current Teacher v3 / Teacher v4-SSL development records.

## Start here

- [Presentation](TWIRL_S56_S62_candidate_model_update_20260806.pptx): 27 slides with a 12-slide meeting narrative followed by a 14-slide atlas of all `80` final Planet-like or broad/wide-transit observations.
- [All-candidate index](candidate_index_all_509.csv): all `509` reviewed observations across `493` TICs, including the final human morphology label and the matching vet-sheet path.
- [Transit/dip-like index](transit_like_candidate_index_80.csv): the `70` Planet-like plus `10` broad/wide-transit observations shown in the deck appendix.
- `transit_like_candidate_index_80.csv` gives every appendix row a concise
  export filename; the presentation embeds the full `80`-sheet atlas directly.
- `candidate_index_all_509.csv` records the canonical workspace-relative
  source-sheet path for the complete uniform current-ADP review set. Match on
  sector/TIC/label to recover the source sheet for an atlas row.

The generated local symlink collections are intentionally not part of the Git
package: their targets live under `data_local/` and would add no portable
content. Copy the canonical images named in the all-candidate CSV if a
standalone image bundle must leave this machine.

## Suggested presentation route

For a focused 10–15 minute update, present slides 1–12:

1. We now have an auditable seven-sector human-review product: `509` observations / `493` TICs.
2. The concrete candidate discussion set is `80` transit/dip-like morphologies, including five TICs repeated in S56 and S57.
3. WD 1856 remains the benchmark and is recovered in both active apertures.
4. Julien's ten-system comparison already exposed a real target-emission defect; A2v1 repaired the product gap, while the tenth signal-level comparison still needs to be rerun.
5. Teacher v3 improves enrichment on the frozen matched test, but rare positive support, label policy, and incomplete survey branches prevent calling it an all-around classifier.
6. Teacher v4-SSL pretraining is complete across five encoders. Its matched transfer evaluation was still running when the deck was frozen, so the slides quote no transfer metric.
7. The useful meeting decision is to agree on a `10–20` object validation shortlist and one predeclared prospective S63 test.

Use slides 13–27 as a discussion atlas, or open the full-resolution sheets from either CSV while discussing individual systems.

## Exact candidate inventory

| Final human morphology | Observations | Unique TICs |
|---|---:|---:|
| Planet-like | 70 | 65 |
| Eclipse/contact | 419 | 408 |
| Broad/wide transit | 10 | 10 |
| Stellar variability | 7 | 7 |
| Uncertain | 2 | 2 |
| Instrumental/systematic | 1 | 1 |
| **Total** | **509** | **493** |

The rows are an enriched morphology-review set, not an unbiased survey sample or a confirmed-object catalog. Repeated hosts are validation priorities, not multi-sector confirmations.

## Model status represented in the deck

- **Teacher v1:** useful S56 active-learning baseline; fixed-test balanced accuracy `0.750`, macro F1 `0.757`, ECE `0.048`, but only two real Planet-like test examples.
- **Teacher v2:** exploratory injection-retention experiment. It improved conditional retention to `77.21%` but missed its `80%` gate and degraded real-morphology performance, so it was rejected for production use.
- **Teacher v3:** current frozen enrichment release using the v1 architecture and a larger, leakage-safe S56–S62 corpus. On the common `528`-row real non-uncertain fixed-test support, primary versus metadata-only balanced accuracy is `0.779` versus `0.720`, macro F1 is `0.809` versus `0.744`, and ECE is `0.025` versus `0.054`. Planet-like recall is `8/14` and remains descriptive.
- **Teacher v4-SSL:** five VICReg encoders completed `20` epochs each over the broad leakage-safe real-light-curve pool. The frozen model-eligible pool has `175,347` observations. The deck treats the UMAP only as a representation/confound diagnostic and quotes no transfer-performance result.

No teacher is represented as an automatic survey classifier, pseudo-label generator, discovery engine, or completeness measurement.

## Post-freeze update

The deck is an intentionally time-stamped meeting artifact and is not being
rewritten after the fact. The matched Teacher v4-SSL development evaluation
subsequently completed on `2026-08-06`: it did not establish a reliable gain
over Teacher v3 and is not promoted. Teacher v3 remains the frozen enrichment
model for the preregistered S63 run. See the final
[SSL development report](../stage5_validation/teacher_ssl_fullpool_v4_development_performance/README.md)
and the [S63 prospective plan](../stage5_validation/teacher_v3_s63_prospective_v1/README.md).

## Provenance

- Frozen human decisions: [accepted signal rereview](../stage5_validation/s56_s62_morphology_corpus_teacher_v3_v1/accepted_signal_rereview.csv)
- Julien comparison: [summary](../stage5_validation/julien_s56_10_systems_compare/summary.md)
- Teacher v3 performance: [report](../stage5_validation/teacher_v3_s56_s62_a2v1_current_adp/performance/README.md)
- Teacher v4-SSL development diagnostics: [report folder](../stage5_validation/teacher_ssl_fullpool_v4_development_performance/)
- Detailed dated history: [TWIRL progress log](../../doc/twirl_progress_log.md)

Each presentation slide also contains a `[Sources]` block in its speaker notes.
