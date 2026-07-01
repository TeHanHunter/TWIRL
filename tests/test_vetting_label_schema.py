from __future__ import annotations

from twirl.vetting import label_schema
from twirl.vetting.lightcurve_label_app import LABEL_BUTTONS, LABEL_KEY_ALIASES, LABEL_OPTIONS


def test_browser_labels_are_schema_labels() -> None:
    assert LABEL_OPTIONS == label_schema.LABEL_OPTIONS
    assert LABEL_BUTTONS == label_schema.LABEL_BUTTONS
    assert LABEL_KEY_ALIASES == label_schema.LABEL_KEY_ALIASES
    assert set(label_schema.STRONG_LABELS).issubset(set(label_schema.LABEL_OPTIONS))
    assert set(label_schema.AUDIT_LABELS).issubset(set(label_schema.LABEL_OPTIONS))
    assert set(label_schema.EXCLUDE_LABELS).issubset(set(label_schema.LABEL_OPTIONS))


def test_teacher_target_and_bls_truth_helpers() -> None:
    assert label_schema.teacher_target("planet_like") == "planet_like"
    assert label_schema.teacher_target("uncertain") == "uncertain"
    assert label_schema.teacher_target("skip") == ""
    assert label_schema.teacher_target("not_a_label") == ""

    assert label_schema.is_bls_truth_match("bls_top1_recovered")
    assert label_schema.is_bls_truth_match("bls_topn_harmonic_match")
    assert not label_schema.is_bls_truth_match("bls_peak_mismatch")
