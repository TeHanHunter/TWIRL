# TWIRL-FM0.2 later-sector admission preflight

Date: `2026-08-26`

## Outcome

The later-sector path is now defined as a separate, fail-closed development
contract. It does not widen the frozen S56--S64 FM0.1/FM0.2 input release.
The [policy](../../../configs/models/twirl_fm0_2_later_sector_admission_v1.yaml),
[inventory builder](../../../scripts/stage5_validation/build_twirl_fm0_later_sector_inventory.py),
and reusable admission logic
([module](../../../src/twirl/models/fm0/temporal_admission.py)) require
an explicitly supplied frozen-policy hash, the unchanged Stage-1 full-schema
acceptance receipt, and distinct gate-specific evidence before a row can enter
a later-sector inventory.

No temporal panel is frozen and no zero-shot model evaluation was run in this
preflight. Only S65 is both later than the training release and accepted, while
the policy requires at least five chronological later sectors plus adequate
new- and repeated-host counts. New training, event retention, the formal
FM0.2 gate, and sealed-test access remain closed.

## Authoritative sector reconciliation

A strict, socket-only read-only probe reached `pdogpu1.mit.edu` as `tehan` and
confirmed `/pdo/users/tehan`. The canonical PDO Stage-1 tree is authoritative.

| Sector block | Observed state | FM admission consequence |
| --- | --- | --- |
| S65 | `228,917` HDF5 products; zero non-edge, unreadable, or zero-byte failures; `114,419` FITS pass the full schema gate | Accepted, but still needs the separate FM cadence/channel/registry admission receipt |
| S66 | All 32 ORCD cells completed; 28 payloads retained; two cam4/ccd4 cells (`3,564` HDF5 total) returned to PDO; both cam4/ccd3 attempt payloads remain intact, but their cleaned source grids must be re-ingressed before retention | Partial/deferred; not admissible |
| S67--S78 | 32/32 input, completion, and retention receipts per sector; every retention receipt says `pdo_return_deferred=true` and `pdo_sector_accepted=false` | Complete compute checkpoints, not accepted sectors |
| S79 | 32 inputs ready; detector compute active | Not admissible |
| S80 | 32 inputs ready; input workflow active | Not admissible |

The operational footprint therefore reaches S80, or 25 sectors from S56, but
the accepted Stage-1 archive remains S56--S65 and the frozen FM release remains
S56--S64. Staged or retained volume is not substituted for accepted temporal
evidence.

## Frozen baseline authorities

The inventory contract binds the real S56--S64 artifacts rather than
reconstructing their schemas from assumptions:

- later-sector policy SHA-256:
  `527a7b4d9f9c452f02576eea7f155abe65ad8439057317e8bc10c80b0ed93da3`;
- input release manifest SHA-256:
  `f9c62ba02ea40ad7b2f457cd4dab55526ee7a453ead8a1a553c4de9f2969598a`;
- label-blind corpus selection SHA-256:
  `3fb2cc3d9cd652999df366ab3cd05f22daa6f23f4ceda66f17c43281e2aaf3b8`;
- Gaia--TIC alias-edge authority SHA-256:
  `03d09fd36df6e24996a4a6688aa78e17f733a80165eb40061b28a821959f8bc8`.

The release manifest binds immutable product instances. The corpus selection
supplies the real sector/camera/CCD visit geometry, and the edge-only alias
authority is rebuilt deterministically under the unchanged
`twirl_fm0_1_source_split_v1` salt. Any later alias edge that changes or merges
a frozen component fails closed.

The production campaign ID is accepted only when the runtime policy hash equals
the in-code frozen hash above. Each sector receipt must also bind the accepted
A2v1 full-schema validation as its Stage-1 acceptance anchor. Every one of the
nine later-admission gates has a distinct JSON artifact with an exact schema,
gate name, sector, accepted state, pass result, gate-specific reconciling
counters, and at least one hash-verified upstream control-metadata artifact; one
generic file cannot stand in for multiple gates. Later observation manifests
use an exact schema and explicitly bind the canonical six-view ordering.

## Why the additional sectors help

The larger sector footprint is valuable first as an evaluation resource. It
provides later spacecraft states and detector placements, many genuinely new
hosts, and repeated hosts that can be queried against their S56--S64 visits.
Those cohorts test whether the learned rank repair transfers across time and
instrument conditions.

More training data may eventually help, but data volume alone does not resolve
the present result: the step-2000 checkpoint learned a much higher-rank pooled
representation while cross-sector retrieval did not improve. The next
discriminating experiment is therefore zero-shot evaluation on a label-blind
later panel, not an immediate larger training run.

## Freeze and execution boundary

The proposed development panel begins with S65 and advances chronologically.
An inventory may be built incrementally, but `panel_frozen` remains false. The
five-sector and `64`/`256` host counts are preliminary count floors, not a panel
freeze authority. The current contract reports `count_floor_ready` while
hard-coding `adequacy_thresholds_frozen=false` and `panel_freeze_ready=false`.
A separate label-blind adequacy policy and freeze can be considered only after
all of the following are true:

- at least five accepted and FM-admitted sectors are present;
- at least 64 repeated-host components have both training-era and later visits;
- at least 256 new-host components are available;
- accepted-sector HDF5/FITS provenance, edge-aware coverage/openability,
  internal and external cadence quality, explicit omissions, FM channel/mask/
  finite/numerical envelopes, and the source-registry join are all hash-bound;
- zero sealed rows are emitted; the inventory builder reads no light curves,
  labels, search results, embeddings, scores, or injections. Numerical and
  channel QA derived upstream from light curves enters only through the
  distinct hash-bound evidence receipts and upstream control metadata, whose
  producer remains subject to a separate label-blind audit.

Once the separate label-blind freeze receipt exists, evaluate the unchanged
step-2000 checkpoint against its exact same-seed step-0 control. Report new and
repeated hosts separately; do not apply the formal canary gate or open the
sealed test.

## Validation

The focused adversarial contract suite passes `100` tests. The complete local
fast suite passes `1,017` tests with `33` optional-environment skips; Ruff,
compilation, CLI-help, diff, and documentation checks pass. A final independent
review found no remaining P1/P2 issue in the policy hash, Stage-1 acceptance
anchor, evidence/upstream provenance, alias closure, manifest schema, cohort
logic, or readiness boundary.

## Next production gate

After the active S80 ingress reaches a stable terminal state, re-ingress only
the two cleaned S66 cam4/ccd3 source grids, retain their intact science
payloads, and return/promote all 30 unpromoted S66 cells one at a time. Then run
the unchanged PDO edge-aware HDF5/openability, sector-FITS, full-schema, and QA
gates before returning S67--S69 in order. S65--S69 can become the first
temporal-panel candidate only after their separate FM admission receipts pass.
