# S56 Active Goal Status

Date: `2026-07-01`

## Objective

Build the S56 path from light curves to injection-aware search, robust BLS peak
ranking, human triage, and later H200-supported ML training.

## Current State

| Requirement | Status | Evidence |
|---|---:|---|
| S56 LC product selected | Ready | Current human/ML path uses the S56 TWIRL-FS v2 compare product with canonical compact-export apertures and ADP display support. |
| Coverage-first all-host injection sample | Ready on PDO | All-host pre-detrend BATMAN run produced `19,072` injections across S56 hosts, with one documented nonpositive-baseline host gap. |
| Balanced recovery/ranker grid | Selected outputs synced | ORCD-selected ephemerides are synced locally and verify with `57,204` rows across `19,068` targets. |
| Robust BLS peak labeling | Ready for all-host PDO path; selected balanced ORCD output available | All-host PDO path already produced the ranker-selected real-candidate LEO queue; balanced ORCD selected ephemerides are now available for comparison/LEO smoke. |
| Human triage queue | Ready | `s56_allhost_ranker_selected_real_leo_queue_pdo/verification.json` passes with `1,000` real candidates, `1,000` LEO PDFs, and `0` LEO metric/plot errors. |
| Human labels | Pending | No `human_labels_vetted.csv` exists yet for the all-host real-candidate queue. |
| Post-label teacher/readiness audit | Pending | PDO monitor `twirl-s56-allhost-label-audit` is waiting for labels and will rebuild label summaries, teacher rows, next-label priorities, and readiness gates. |
| Label-to-ORCD staging | Ready, gated on labels | `stage-labels` now stages only small label/readiness products to ORCD after the post-label readiness gate passes. |
| ORCD selected-output handoff | Selected bundle ready; local LEO smoke artifact not synced | Local `s56_ranker_selected_real_candidates_orcd/` verifies with `57,204` selected rows. The monitor log says the bounded PDO LEO smoke completed, but the `orcd50` queue directory is not present locally. |
| H200 infrastructure | Ready for bounded smokes | Torch H200 tensor smoke and synthetic-label training smoke pass on one H200 with CUDA visible. |
| Real H200 training | Pending | Correctly gated on human-label readiness; no real-label H200 training smoke has been run. |

## Active Monitors

- `twirl-orcd-apply-leo-smoke`: waits for ORCD selected ephemerides, syncs to local/PDO, renders `50` LEO reports on PDO.
- `twirl-s56-allhost-label-audit`: waits for PDO labels and rebuilds teacher/readiness products.
- `twirl-pdo-label-sync`: pulls small PDO label/audit products back locally and refreshes `s56_ml_handoff_readiness/`.

## Next Trigger

1. If human labels appear first: the PDO audit monitor rebuilds teacher rows and readiness; local sync refreshes the ML handoff report. If readiness passes, run `stage-labels`, then submit the real-label H200 smoke.
2. If ORCD artifacts are needed next: inspect the synced selected ephemerides locally and, if needed, sync or rebuild the bounded `orcd50` PDO LEO queue locally.
3. Do not scale H200 training until the real-label readiness gate passes.
