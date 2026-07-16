# S56 ML Handoff Readiness

Created UTC: `2026-07-11T15:23:07.324272+00:00`

## Status

| Gate | Status | Evidence |
|---|---:|---|
| `pdo_leo_queue_verified` | `pass` | passed=True; n_rows=1000; leo_class_counts={'FA': 994, 'FP': 4, 'PC': 2} |
| `human_labels_present` | `pending` | artifact is missing or empty |
| `human_training_readiness` | `pending` | post-label readiness audit has not been run |
| `orcd_selected_ephemerides_synced` | `pass` | passed=True; rows=57204 |
| `orcd_torch_env_documented` | `pass` | torch build info present |
| `h200_tensor_smoke` | `pass` | n_tensor_rows=128; shape=[128, 3, 257]; cuda_available=True; device=cuda:0 |
| `h200_synthetic_train_smoke` | `pass` | synthetic_label_smoke=True; cuda=True; device=cuda:0; split_counts={'test': 26, 'train': 76, 'validation': 26} |
| `h200_real_label_train_smoke_complete` | `pending` | real-label training smoke has not been run |

## Flags

| Flag | Ready |
|---|---:|
| `ready_for_human_triage` | `True` |
| `ready_for_post_label_audit` | `False` |
| `ready_to_submit_real_h200_training_smoke` | `False` |
| `orcd_selected_outputs_ready_for_pdo_leo` | `True` |
| `real_h200_training_smoke_complete` | `False` |

## Next Step

`human_labels_present`: Use the vetting app; the PDO label-audit monitor or run_s56_allhost_real_label_audit_pdo.sh will rebuild readiness after labels appear.
