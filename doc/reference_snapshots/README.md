# External Roadmap Snapshots

**Snapshot date:** `2026-04-16`, from the PDFs' embedded creation metadata.

The repository retains two rendered exports of the MIT `qlp-ops` wiki at their
original paths so existing links do not break:

- [Photometry Roadmap PDF](<../Photometry Roadmap · mit-kavli-institute:qlp-ops Wiki.pdf>)
- [Planet Search Roadmap PDF](<../Planet Search Roadmap · mit-kavli-institute:qlp-ops Wiki.pdf>)

These files preserve operator context as it appeared on the capture date. They
are **snapshots, not live runbooks**. Relative page labels such as “edited
yesterday” or “edited last month” refer to the captured web page and must not be
interpreted relative to the current date.

The PDFs describe the shared QLP photometry and planet-search workflows,
including PDO paths, database ingestion, quality-flag review, BLS, Astronet,
centroid filtering, and TEV delivery. Those procedures are useful for
understanding the upstream environment, but they do not override TWIRL's
current product contract, write boundaries, or scientific search design.

For current TWIRL decisions, use these sources in order:

1. [`twirl_plan.md`](../twirl_plan.md) for the forward-looking project plan.
2. [`a2v1_production_protocol.md`](../a2v1_production_protocol.md) for current
   Stage 1 production and acceptance.
3. [`mit_tglc_usage_guide.md`](../mit_tglc_usage_guide.md) for MIT TGLC usage.
4. These PDFs only for dated QLP-operator context.

Before reusing any command from a snapshot, verify the current QLP wiki,
installed software version, paths, database expectations, and active TWIRL
runbook. In particular, commands shown against shared `/pdo/qlp-data/` do not
authorize TWIRL automation to modify that tree; the repository's current PDO
safety rules still apply.
