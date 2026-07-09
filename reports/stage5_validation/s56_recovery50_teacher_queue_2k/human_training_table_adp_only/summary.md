# S56 Labeled 2k ADP-Only Audit

- rows: `2000`
- human display exactly matched ADP small: `1999`
- human display missing ADP small: `1`
- rows with ADP primary supplement: `1843`
- rows with historical canonical supplement: `157`
- active ADP-only rows: `1842`
- excluded pending ADP-only review: `158`

The original queue period frequently differs from the ADP-small BLS anchor, but the
vet sheets used the ADP-small anchor for 1,999/2,000 rows. The lone missing row was
already labeled skip. Historical canonical columns are removed from the active table.
