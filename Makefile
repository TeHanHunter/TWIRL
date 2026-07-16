PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python)
PYTEST := $(PYTHON) -m pytest

.PHONY: test-fast run-detection-sample check-docs

test-fast:
	$(PYTEST) -q

run-detection-sample:
	$(PYTEST) -q tests/test_detection_sample.py

check-docs:
	$(PYTHON) scripts/check_docs.py
