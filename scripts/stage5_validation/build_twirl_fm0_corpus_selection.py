#!/usr/bin/env python3
"""Build the BLS-free Gaia-source selection for the FM0.1 corpus."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from twirl.models.fm0.corpus import (
    iter_observation_fits,
    select_corpus_observations,
    write_corpus_selection,
)
from twirl.models.fm0.registry import build_alias_registry, read_rows, sha256_file


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aliases", required=True, type=Path)
    parser.add_argument("--observations-fits", required=True, type=Path)
    parser.add_argument("--from-sector", type=int, default=56)
    parser.add_argument("--through-sector", type=int, default=64)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.from_sector < 56 or args.through_sector < args.from_sector:
        raise SystemExit("invalid FM0.1 sector range")
    sectors = list(range(args.from_sector, args.through_sector + 1))
    aliases = build_alias_registry(read_rows(args.aliases))
    selection, selected_aliases, audit = select_corpus_observations(
        iter_observation_fits(args.observations_fits, sectors=sectors),
        aliases,
        sectors=sectors,
    )
    summary = write_corpus_selection(
        args.out_dir,
        selection,
        selected_aliases,
        audit,
        sectors=sectors,
        input_authorities={
            "gaia_tic_alias_table": {
                "path": str(args.aliases.resolve()),
                "sha256": sha256_file(args.aliases),
            },
            "observation_fits": {
                "path": str(args.observations_fits.resolve()),
                "sha256": sha256_file(args.observations_fits),
            },
        },
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
