# local_data.md

This note records local-only data conventions for TWIRL. It complements
`doc/twirl_plan.md` and should stay focused on paths, provenance, and what
must not be committed.

## Policy

- Large external inputs and staged survey products belong under `data_local/`
  or another user-configured local path, not in normal git-managed source
  directories.
- Do not commit raw FITS catalogs, staged TICA/TGLC products, or full survey
  light-curve trees.
- Keep code, schemas, configs, compact manifests, and reproducible derived
  table definitions in git.
- Every generated product that depends on local data should record provenance:
  input path, file size, modification time, build version, and important sample
  cuts. Use a hash when practical.

## Recommended Local Layout

```text
data_local/
  catalogs/
    GaiaEDR3_WD_main.fits
    twirl_master_catalog/
  tglc-data/
  stage2/
```

The current seed WD catalog convention is:

```text
data_local/catalogs/GaiaEDR3_WD_main.fits
```

This is the local Gentile Fusillo et al. (2021) main Gaia EDR3 white dwarf
catalog. It is an external input, not a repo asset.

## Production Trees

PDO production products should live under `/pdo/users/tehan/`. Shared PDO trees
such as `/pdo/qlp-data/` are read-only inputs. If a PDO workflow needs a shared
file, stage a user-owned symlink or copy under `/pdo/users/tehan/` rather than
changing the shared source tree.

Accepted master-catalog states should use concise, stable filenames, but do not
duplicate large integrated FITS products just to create cleaner names when the
existing file already represents the accepted state.
