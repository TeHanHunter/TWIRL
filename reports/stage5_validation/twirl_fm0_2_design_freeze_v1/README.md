# TWIRL-FM0.2 objective-canary freeze

This directory binds the gated `TWIRL-FM0.2.1` TCN objective canary. The
scientific change from FM0.1.1 is limited to a scale-matched same-window VICReg
term applied directly to `z_window`. The S56--S64 input bytes, Gaia-component
split, ADP `1x1+3x3` views, model geometry, context, sampler, masks, optimizer,
and `20,000`-step scheduler horizon remain fixed.

The user authorized a real-data canary only after every prestart gate passes.
The first real invocation may stop only at step `64`; exact resumes may then
stop at `500`, `1,000`, and `2,000`. These stops are execution state within one
run contract, not new scientific contracts. No invocation may pass step
`2,000` under this freeze.

The tracked input-reuse receipt binds the already validated FM0.1 S56--S64
release without rebuilding or reassigning any source. The development event
contract freezes its source partitions, deterministic sample selection,
injection grid, metric, baseline, and tolerances before an H200 is used. Its
result is deliberately absent from this directory and remains a post-step-2000
development gate. The evaluator-v2 contract similarly freezes the `h_window`
and `z_window` diagnostics, paired trained-minus-random comparison, genuinely
query-sector-excluded retrieval, train-fit PCA and robust baselines, output
schema, and numeric step-2,000 gate paths.

The sealed `poc_sealed_test` partition remains closed. Nothing here authorizes
FM0.2.2, a 20,000-step training run, production use, a science result, or a
foundation-model claim.
