# Blackhole/Globus transfer probe — 2026-08-07

## Outcome

The new `TESS TSO` collection is accessible by `tehan@mit.edu`, and a
checksum-verified transfer from Blackhole to ORCD completed without faults.
The measured throughput is nevertheless far below the expected rate and is
not yet accepted for the A2v1 production stream.

- Source collection: `TESS TSO`
  (`d8ea14bc-dca1-4cbc-92b6-b76d7289b6d2`), where collection path `/` maps to
  `blackhole:/globus/tso`.
- Destination collection: `MIT ORCD Engaging Collection`
  (`ec54b570-cac5-47f7-b2a1-100c2078686f`).
- Transfer task: `8154fb53-926b-11f1-b41f-02ce27bde401`.
- Payload: eight 512 MiB files, 4 GiB total.
- Result: `SUCCEEDED`, 8/8 files, zero faults, encrypted data path, checksum
  verification requested.
- Elapsed: 139 s from request to completion.
- Effective rate: `30,808,100 B/s = 30.8 MB/s = 246 Mb/s`.
- Main task progress intervals reported `257--271 Mb/s`.
- Globus selected GridFTP `concurrency=4`, `parallelism=8`, and
  `pipelining=20`; the CLI did not under-request transfer concurrency.
- A Blackhole-local 1 GiB `fdatasync` write to `/globus/tso` reached
  `706 MB/s`, so the observed 30.8 MB/s transfer is not explained by that
  single local write test.

An eight-task concurrent probe did not recover aggregate bandwidth. Only a
small subset of tasks transferred concurrently, at roughly `8--10 MB/s` each,
while the others remained active at zero bytes. Redundant active probes were
canceled once this shared-cap behavior was established. All temporary 4 GiB
source and destination payloads were subsequently deleted; no production data
were changed.

At the observed completed-task rate, transfer-only estimates are:

| Prepared-source unit | Bytes | Estimated time |
| --- | ---: | ---: |
| S65 orbit 138 cam4/ccd1 | 101,012,223,860 | 54.6 min |
| S65 orbit 137 cam4/ccd1 | 115,057,047,109 | 62.2 min |
| S66 orbit 139 cam1/ccd1 | 118,824,120,305 | 64.3 min |
| S66 orbit 140 cam1/ccd1 | 101,988,405,375 | 55.2 min |
| Representative two-orbit cam/CCD | about 220 GB | about 2.0 h |

The requested `1 GB/s` target would be about 32 times faster than this probe.

## Access and software state

- `tehan` can log in to `blackhole.mit.edu` through the working `hostess3`
  route, belongs to groups `tso` and `globus`, and can create group-writable
  content under `/globus/tso`.
- Direct TCP connections to Blackhole timed out from the laptop, PDO, and ORCD;
  the `hostess3` proxy route reached it successfully.
- Blackhole now has the documented PDO SSH host blocks. The forwarded TESS RSA
  key is accepted by `hostess3`, after which the server requires an interactive
  authentication step. Automation therefore still needs a user-opened
  Blackhole-to-PDO control socket; the existing laptop-to-PDO socket is not a
  substitute.
- The official Globus CLI `3.42.0` is installed through `pipx` at
  `~/.local/bin/globus`. The MIT identity and both collection data-access
  consents are active; refreshable CLI credentials remove per-transfer browser
  login from the planned controller.
- The current ORCD TWIRL environments contain Astropy and the downstream
  Torch environment, but neither contains CuPy or the MIT/TWIRL TGLC fork.
  A separate A2v1 H200 environment and parity smoke are required before any
  production cell runs there.

## Efficient A2v1 staging unit

Use one `(orbit, camera, CCD)` prepared-source cell as the atomic unit. Do not
transfer both raw TICA FFIs and prepared source pickles. The accepted A2v1
reuse path begins from the immutable 196 prepared `source_*.pkl` files plus
small TIC overlays, so rebuilding cutouts from FFIs on ORCD would duplicate
work and require PDO-only catalog infrastructure.

The intended bounded stream is:

1. PDO prepared source cell to Blackhole by restartable `rsync`.
2. Blackhole to an ORCD `.partial` staging cell through Globus.
3. Count/hash/openability gate, then atomic ready marker.
4. One H200 per orbit cell, with two independent cells active at most.
5. Validate ePSF/HDF5 outputs before Globus egress.
6. Return ePSF and HDF5 to the user-owned PDO A2v1 tree and revalidate there.
7. Retain HDF5 on PDO and ORCD; remove temporary ORCD source/ePSF and
   Blackhole staging only after the returned PDO hashes pass.
8. Build sector FITS and run the unchanged full-product gate on PDO after both
   sector orbits are complete.

Cap live ORCD use at two H200s and 78 CPUs total. A likely first allocation is
two one-H200 orbit-cell workers with no more than 64 CPUs between them, leaving
headroom for bounded validation while staying below the cap.

## Blockers before the first H200 smoke

1. MKI/ORCD should inspect why the managed endpoint pair sustains only about
   `250 Mb/s` despite the automatically selected `4×8` GridFTP settings. Check
   Blackhole gateway NIC/data-node configuration, DTN routing, TCP windows,
   storage path, encryption cost, and endpoint/task concurrency policy.
2. Open and verify a reusable Blackhole-to-PDO control socket from a user
   terminal without copying credentials to Blackhole.
3. Build a version-pinned ORCD A2v1 environment containing the exact TGLC
   fork hooks and CUDA-compatible CuPy.
4. Transfer only a bounded S65 orbit-138 cam4/ccd1 smoke subset first, then run
   one-H200 numerical/product parity checks before enabling the two-H200
   pipeline.

