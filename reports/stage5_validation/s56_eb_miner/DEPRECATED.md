# Deprecated EB Miner Product

Do not use this directory's model, scores, candidate periods, or queue as an
active training/search product. Its candidate ephemerides came from the older
canonical `DET_FLUX*` ranker while the CNN and vet sheets read ADP light
curves. Human review of the first 59 rows found zero EB/PCEB examples.

The labels and sheets are retained only as audit provenance. The replacement
pipeline is under `reports/stage5_validation/s56_eb_miner_adp_only/` and
requires the `s56_adp_pair_v1` contract.
