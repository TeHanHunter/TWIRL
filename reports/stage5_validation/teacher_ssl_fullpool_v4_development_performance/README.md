# Teacher v4-SSL full-pool development evaluation

This report is the final matched development-OOF comparison between the
corrected S56--S62 full-pool Teacher v4-SSL experiment and the frozen Teacher
v3 enrichment release. S63 and the frozen fixed-test tensors were not opened
for this comparison.

## Result

Fine-tuned SSL does not establish a reliable improvement over Teacher v3.

| Label policy | Model | Balanced accuracy | Macro F1 | Macro AP | ECE |
| --- | --- | ---: | ---: | ---: | ---: |
| Uncertain as other | Teacher v3 | 0.801 | 0.789 | 0.805 | 0.005 |
| Uncertain as other | Fine-tuned SSL | 0.775 | 0.789 | 0.796 | 0.006 |
| Uncertain masked | Teacher v3 | 0.822 | 0.802 | 0.853 | 0.033 |
| Uncertain masked | Fine-tuned SSL | 0.821 | 0.808 | 0.803 | 0.031 |

The aggregate bootstrap intervals overlap. On the rarer Planet-like class,
uncertain-masked average precision is `0.763` for Teacher v3 and `0.614` for
fine-tuned SSL. The frozen SSL linear probe is substantially worse under both
label policies, so the representation alone is not a useful classifier in
this setup.

Teacher v4-SSL is therefore archived as development evidence and is not
promoted, used to rank S63, or used to create pseudo-labels. Teacher v3 remains
the operational enrichment model. The UMAP is an exploratory representation
and confound diagnostic only, not a performance estimate.

## Files

- `development_performance.csv`: exact aggregate metrics and TIC-clustered
  bootstrap intervals.
- `development_per_class_metrics.csv`: exact per-class point estimates.
- `teacher_ssl_fullpool_development_performance.{png,pdf}`: aggregate and
  class-level comparison figure.
- `reference_fold_0_umap_coordinates.csv` and
  `teacher_ssl_fullpool_reference_umap.{png,pdf}`: the frozen exploratory
  reference-encoder projection.
- `teacher_ssl_fullpool_figures.provenance.json`: input and figure hashes.
