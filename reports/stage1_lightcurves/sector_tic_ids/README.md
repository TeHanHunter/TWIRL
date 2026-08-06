# A2v1 Sector TIC exports

These files provide the complete, one-row-per-TIC target inventories for the
accepted A2v1 S60--S63 products. They are intended for collaborator MAST
crossmatches.

- `tic_ids_s00xx_A2v1_targets.csv` has the simple `tic_id` column used for
  MAST queries.
- `tic_ids_s00xx_A2v1_target_context.csv` adds Gaia DR3 `source_id`, the
  Gentile Fusillo `Pwd`, the `Pwd > 0.75` membership flag, and the two TESS
  orbits.
- The combined S60--S63 files retain a `sector` column, so a TIC appearing in
  more than one sector remains an unambiguous observation-level record.

The selection is every unique positive TIC in the corresponding TWIRL
observation-table sector. It intentionally has **no** post-production `Pwd`
cut: these are the original broad candidate-target lists, not the
high-confidence WD reference sample. The companion context files let a user
apply `is_highconf_wd == 1` (equivalently, `Pwd > 0.75`) when that is the
desired scientific subset.

| Sector | TICs | `Pwd > 0.75` | Lower-`Pwd` candidates |
| --- | ---: | ---: | ---: |
| S60 | 27,165 | 19,652 | 7,513 |
| S61 | 41,403 | 23,501 | 17,902 |
| S62 | 40,158 | 22,070 | 18,088 |
| S63 | 53,512 | 22,685 | 30,827 |

The JSON manifest records the input observation-table hash and the hashes of
every delivered CSV.
