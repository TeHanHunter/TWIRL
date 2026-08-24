# TWIRL documentation guide

This page is the entry point to TWIRL documentation. Most work requires only
the [project plan](twirl_plan.md) plus one task-specific reference below.
Do not read the entire documentation tree before starting a bounded task.

## Start here

1. [Project plan](twirl_plan.md) — current status, locked decisions, stage
   gates, and immediate priorities. This is the main human-readable document.
2. [Progress log](twirl_progress_log.md) — dated execution evidence. Read its
   current snapshot first; consult older entries only when tracing a result.
3. [Open questions](ideas.md) — choices that have not yet been decided.

The root [README](../README.md) is a short project introduction. Files under
[reports](../reports/README.md) are immutable or dated evidence for individual
runs; they do not override the plan.

## Read by task

| Task | Read this | Purpose |
| --- | --- | --- |
| Stage 1 A2v1 production or acceptance | [A2v1 production protocol](a2v1_production_protocol.md) | Frozen product recipe, state machine, gates, and campaign resource profile |
| MIT TGLC implementation details | [MIT TGLC guide](mit_tglc_usage_guide.md) | Fork behavior, PDO layout, and compatibility background |
| ORCD access or jobs | [ORCD guide](orcd_h200_usage.md) | Authentication boundary, partition, storage, and generic Slurm policy |
| Local or remote data placement | [Data policy](local_data.md) | Canonical roots, ignored data, and provenance requirements |
| Detrending or light-curve methods | [TWIRL-FS methods](twirl_fs_methods.md) | Flux-space method definitions and assumptions |
| Foundation-model implementation | [FM0.1 frozen design](foundation_model_design.md), [FM0.1 config](../configs/models/twirl_fm0_1_s56_s67_poc.yaml), [freeze receipt](../reports/stage5_validation/twirl_fm0_1_design_freeze_v1/freeze.json), and [FM0.2 objective-canary draft](../configs/models/twirl_fm0_2_s56_s64_objective_canary.yaml) | Frozen proof-of-concept authority plus the current fail-closed repair proposal |
| Scientific or literature claims | [Science background](science_background.md) | Scope and source anchors |
| Publication-facing figures | [Plotting style](plotting_style.md) | Shared visual conventions |

## Authority by topic

- The **plan** decides what the project is doing next.
- A **protocol or runbook** defines how an authorized workflow is executed.
- The **progress log** records what happened and when.
- A **report** binds evidence to named inputs and code.
- The **ideas file** contains unresolved choices only.

The plan decides whether and when a workflow is authorized, paused, or
retired. Within an authorized run, its frozen configuration and protocol
control scientific mechanics and acceptance; an older frozen artifact cannot
override a later plan-level pause. Stop and reconcile any scientific or
acceptance contradiction before execution.

## Historical material

The following files are retained only for provenance and should not be used as
current instructions:

- [pre-A2v1 sector protocol](sector_production_protocol.md)
- [QLP-standard transition plan](qlp_standard_transition_plan.md)
- [archived ideas and progress](archive/README.md)
- [external reference snapshots](reference_snapshots/README.md)

Keep new current decisions in the plan, detailed run events in the progress
log, and generated evidence in reports. Do not create a new planning document
when one of those homes already fits.
