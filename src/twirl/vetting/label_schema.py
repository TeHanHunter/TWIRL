"""Shared human-label and training-target schema for TWIRL vetting."""
from __future__ import annotations

LABEL_OPTIONS: tuple[str, ...] = (
    "planet_like",
    "eclipsing_binary_or_pceb",
    "stellar_variability",
    "instrumental_or_systematic",
    "centroid_contaminant",
    "uncertain",
    "skip",
)

LABEL_BUTTONS: tuple[tuple[str, str, str], ...] = (
    ("1", "planet_like", "Planet"),
    ("2", "eclipsing_binary_or_pceb", "EB/PCEB"),
    ("3", "stellar_variability", "Variable"),
    ("4", "instrumental_or_systematic", "Systematic"),
    ("5", "uncertain", "Flat/no signal"),
    ("0", "skip", "Skip"),
)

LABEL_KEY_ALIASES: dict[str, str] = {
    "p": "planet_like",
    "e": "eclipsing_binary_or_pceb",
    "v": "stellar_variability",
    "i": "instrumental_or_systematic",
    "f": "uncertain",
    "u": "uncertain",
    "s": "skip",
    "x": "skip",
}

STRONG_LABELS: frozenset[str] = frozenset(
    {
        "planet_like",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "instrumental_or_systematic",
        "centroid_contaminant",
    }
)

AUDIT_LABELS: frozenset[str] = frozenset({"uncertain"})
EXCLUDE_LABELS: frozenset[str] = frozenset({"skip"})

BLS_TRUTH_MATCH_MODES: frozenset[str] = frozenset(
    {
        "bls_top1_recovered",
        "bls_topn_recovered",
        "bls_topn_harmonic_match",
    }
)

OBJECT_TEACHER_LABELS: tuple[str, ...] = (
    "planet_like",
    "eclipsing_binary_or_pceb",
    "stellar_variability",
    "instrumental_or_systematic",
    "centroid_contaminant",
)

VISIBILITY_POSITIVE_LABELS: tuple[str, ...] = ("planet_like",)
VISIBILITY_NEGATIVE_LABELS: tuple[str, ...] = ("instrumental_or_systematic",)


def teacher_target(label: str) -> str:
    """Return the teacher target for a human label, or ``""`` if excluded."""

    label = str(label)
    if label in STRONG_LABELS:
        return label
    if label in AUDIT_LABELS:
        return "uncertain"
    return ""


def is_bls_truth_match(mode: str) -> bool:
    """Return whether a BLS recovery mode corresponds to injected truth."""

    return str(mode) in BLS_TRUTH_MATCH_MODES
