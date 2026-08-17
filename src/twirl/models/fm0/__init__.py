"""Experimental TWIRL-FM0.1 registry, input, and model mechanics.

The package implements the frozen proof-of-concept contract and an explicitly
synthetic training smoke.  It does not certify A2v1 products, materialize the
real survey input release, admit an ORCD run, or authorize scientific claims.
"""

from .registry import (
    AliasRegistry,
    FrozenContract,
    build_alias_registry,
    build_observation_registry,
    deterministic_source_partition,
    load_frozen_contract,
)

__all__ = [
    "AliasRegistry",
    "FrozenContract",
    "build_alias_registry",
    "build_observation_registry",
    "deterministic_source_partition",
    "load_frozen_contract",
]
