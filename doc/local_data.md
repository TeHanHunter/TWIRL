# local_data.md

This note describes large local-only inputs and staged products that should not be committed to git.

Recommended layout:

```text
data_local/
  catalogs/
    GaiaEDR3_WD_main.fits
    twirl_master_catalog/
  tglc-data/
doc/
  local_data.md
```

Current policy:

- keep large FITS catalogs here, not at the repo root
- keep built local TWIRL catalog products under `data_local/catalogs/twirl_master_catalog/`
- keep staged MIT TGLC products here or in another local path configured by the user
- record the exact local input path, file size, and provenance in derived metadata
- do not commit raw external data into this repository
