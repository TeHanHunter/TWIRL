# TWIRL repository audit and near-term execution plan

**Internal technical report - 2026-07-16**

## Abstract

We audited the TWIRL repository as an executable survey pipeline, with emphasis on reproducibility, production safety, validation semantics, generated-artifact hygiene, and alignment between the Stage 1-5 plan and the code that now exists. The repository has a credible pilot foundation: Stage 1 produces validated A2v1 light curves, Stage 2 supports an interpretable periodic search, and the fast validation suite is healthy. It is not yet a survey-complete system. The nearest defensible milestone is an enrichment-ready S56 periodic path: complete the bounded Tier-1 evidence, freeze a confident candidate/label set, and retrain the teacher-v1 baseline while S59-S63 production continues independently. The dip branch, multi-sector merging, false-alarm calibration, frozen-product recovery, and release manifest remain mandatory before a survey-wide science product, but they need not delay this bounded milestone.

## Methods

The scan combined repository-wide inventory, Git/worktree inspection, parser and link checks, shell and Python syntax checks, the project fast test suite, and review of the authoritative plan, production protocol, data conventions, and plotting rules. The integrated audit checkpoint contained 2,147 tracked files (1,056 under `reports/`, 601 under `outputs/`, 329 under `scripts/`, 63 under `src/`, and 59 under `tests/`) and 759.9 MiB of tracked working-tree content. The fast suite completed with 312 passed and 3 skipped tests. Code paths were compared with the declared Stage 1-5 contract. We distinguished Tier-0 integrity (schema, coverage, openability, and benchmark checks) from the bounded Tier-1 enrichment contract (population precision, authoritative cadences, aperture consistency, injection preservation, and an independent extraction comparison) and from the still-future full-product science-release gate. Generated data and historical binaries were treated as reproducible artifacts unless they were compact evidence needed to interpret a result.

## Results

![Stage readiness and next gates](stage_readiness.png)

**Figure 1.** Audit-assessed engineering maturity and the next evidence gate for each pipeline stage. The readiness scale describes implementation maturity, not scientific performance. Stage 1 and the periodic part of Stage 2 are at pilot maturity; Stage 3 and Stage 5 remain prototypes; the Stage 4 survey-wide inference layer is not yet built.

The strongest result is separation of the production contract from later analysis: A2v1 generation, validation, compact export, and an interpretable BLS path now have explicit boundaries. The new fail-closed Tier-1 implementation can qualify the two active S56 ADP channels for bounded enrichment, but it cannot set `science_ready=true`. Its production run remains intentionally blocked until the authoritative cadence product, exact four-shard hashes, and genuinely independent WD 1856 evidence are built and pinned. The audit also found a mismatch between the clean forward plan and experimental teacher-v2 and S57 labeling work already present in the repository. These products should be preserved as exploratory evidence, but they should not silently redefine the production baseline or consume more of the intended validation holdout.

Repository hygiene is acceptable at the source level but weak at the artifact-history boundary. At scan time the Git object store was approximately 14 GiB and dominated by historical generated outputs. The on-disk `reports/` tree was approximately 5.5 GiB because it also contained large ignored or untracked review artifacts, while `src/`, `scripts/`, and `tests/` together occupied only about 15 MiB. Nonportable symlinks, stale caches, large third-party reference bundles, and redundant rendered review material should remain outside versioned survey products. History rewriting is not recommended during the active branch stack; compact manifests and an explicit artifact allowlist are safer immediate controls. The local `.venv` reported NumPy 2.2.6 even though `pyproject.toml` requires `numpy>=1.24,<2.0`; the tests passed, so this is environment drift rather than an observed failure.

![Candidate-retention attrition](candidate_retention_attrition.png)

**Figure 2.** Candidate-retention diagnostic for the 20,000-signal teacher-v1 experiment. The top-five BLS stage retained 3,458 signals (17.29% of injections), and the teacher retained 1,527 (7.635% of injections; 44.16% conditional on top-five BLS recovery). This diagnostic does not include the full search, vetter, candidate-merging, or pixel-level calibration chain and therefore is not end-to-end survey completeness.

The attrition is scientifically useful because it localizes losses, but it also shows why classifier iteration cannot substitute for baseline search and calibration. The current language should reserve “end-to-end recovery” for a frozen chain that includes search, optional ranking, vetting, candidate merging, and the adopted calibration products. The present teacher result is more accurately a candidate-retention efficiency.

### Period-radius candidate-retention surface

![A2v1 Teacher-v1 BLS top-five candidate-retention surface](a2v1_teacher_v1_bls_top5_recovery_surface.png)

**Figure 3.** A2v1 Teacher-v1 BLS top-five candidate-retention surface for 20,000 injected signals, split into four Tmag strata. Recovery fractions are kernel-smoothed; each panel prints its recovered and injected support counts, and the bright bin is comparatively sparse. The map diagnoses where the BLS candidate list retains injected signals. It is not end-to-end survey completeness and should not be used as an occurrence-rate selection function.

The surface supplies more structure than the aggregate attrition count: candidate retention depends on injected period, companion radius, and target brightness. Its role here is to expose support and search behavior, not to promote a final completeness claim. A frozen survey analysis must propagate the same injections through the adopted vetter, candidate merger, and calibration chain before interpreting the surface as completeness.

## Recommendations

1. **Finish the bounded S56 Tier-1 evidence.** Build and pin the original-SPOC-backed cadence reference, the exact four injection shards, and the external WD 1856 comparison; rerun Tier 0 and publish the TIC-level target pass mask. This qualifies enrichment only, never science release.
2. **Freeze a confident S56 candidate and label set.** Apply the target mask, finish the bounded periodic enrichment review, preserve the two-aperture evidence, and merge only adjudicated labels with explicit provenance.
3. **Harden teacher v1 around the data.** Keep the accepted model family; improve versioned label manifests, TIC-grouped and source-separated evaluation, calibration, and bootstrap uncertainty. The small real positive set is the limiting factor, not model capacity.
4. **Keep production independent.** Let the stop-on-failure S59-S63 queue continue and validate each sector, but do not wait for the queue before completing the S56 enrichment milestone. Keep S57 labels and teacher v2 quarantined as exploratory evidence.
5. **Add the deferred survey branches after the periodic path is robust.** Implement the dip branch, multi-sector consolidation, and branch-aware false-alarm strategy, then run the adopted search/ranker/vetter/merger plus a representative pixel-level calibration subset before full-survey enrichment or completeness claims.
6. **Freeze the release and artifact boundary.** Record the TWIRL-I sector cutoff, parent/sample version, accepted products, search contracts, and checksums. Version compact summaries, manifests, selected figures, and small tables; keep rendered sheets, dependency bundles, shards, and local literature outside Git.
7. **Rebuild and automate the release environment.** Rebuild from `pyproject.toml` before release validation and add pull-request CI for `make test-fast` and `make check-docs` once the lock strategy is settled.

A practical execution order is: (i) continue S59-S63 production in parallel, (ii) finish and pin the S56 bounded Tier-1 evidence, (iii) freeze the confident periodic candidate/label set, (iv) retrain and audit teacher v1, (v) add dip, multi-sector, false-alarm, and frozen-chain recovery gates, and only then (vi) perform survey-wide enrichment and science validation.

## Conclusion

TWIRL is beyond a collection of experiments: it has a testable pilot pipeline and a clear production contract. The best next advance is evidentiary rather than architectural: qualify S56 for bounded enrichment, turn that into a confident candidate/label set, and make teacher v1 auditable around the small real training sample. The deferred search branches and end-to-end calibration remain the boundary between this near-term milestone and a science-ready survey. The audit figures and this report are reproducible from [`scripts/build_repository_audit_report.py`](../../scripts/build_repository_audit_report.py); the project plan and detailed execution record remain authoritative in [`doc/twirl_plan.md`](../../doc/twirl_plan.md) and [`doc/twirl_progress_log.md`](../../doc/twirl_progress_log.md).
