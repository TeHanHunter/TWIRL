# S56 Real Planet-Like Revisit Queue

- Source rows previously labeled `planet_like` and real: `40`
- Revisit rows written: `40`
- Current-fold rows: `40`
- Half-period refold rows: `0`

Use `wide_transit_like` for broad/long-duration transit-like signals that should be preserved for audit but should not train the compact planet-like class.
Use `eclipsing_binary_or_pceb` when the half-period fold makes the signal look EB/PCEB-like.
The original human labels are not overwritten; new labels belong in `human_labels_revisit.csv`.
