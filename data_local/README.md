# Local Data

This directory is for large local-only inputs and staged products that should not be committed to git.

Recommended layout:

```text
data_local/
  README.md
  catalogs/
    GaiaEDR3_WD_main.fits
  tglc-data/
```

Current policy:

- keep large FITS catalogs here, not at the repo root
- keep staged MIT TGLC products here or in another local path configured by the user
- record the exact local input path, file size, and provenance in derived metadata
- do not commit raw external data into this repository
