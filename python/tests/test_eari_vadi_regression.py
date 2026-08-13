"""Tests for EARI/VADI validation regression (issue #56)."""

from __future__ import annotations

from pathlib import Path

import pytest
from gut_ibm_tools.hdf5_reader import GutIBMData
from gut_ibm_tools.validation_regression import (
    FISH_TARGETS,
    VALIDATION_TARGETS,
    check_fish_targets,
    check_invariants,
    check_thresholds,
    evaluate_fish_metrics,
    evaluate_metrics,
)


def test_validation_targets_reference_eari_vadi() -> None:
    refs = {
        ref
        for spec in VALIDATION_TARGETS.values()
        for ref in spec.get("references", [])
    }
    assert any("EARI" in r for r in refs)
    assert any("VADI" in r for r in refs)


def test_evaluate_metrics_sample_hdf5(sample_hdf5: Path) -> None:
    metrics = evaluate_metrics(sample_hdf5)
    assert metrics["monochromatic_score"] > 0.7
    assert metrics["comet_tail_ratio"] > 1.0
    assert metrics["resident_retention"] == pytest.approx(2 / 3, rel=1e-3)


def test_check_thresholds_passes_synthetic_spatial(sample_hdf5: Path) -> None:
    metrics = evaluate_metrics(sample_hdf5)
    # Synthetic fixture meets spatial targets; retention is below full-run window.
    spatial_only = {
        k: metrics[k]
        for k in ("monochromatic_score", "comet_tail_ratio", "comet_tail_asymmetry")
    }
    failures = check_thresholds(spatial_only)
    assert failures == []


def test_check_thresholds_fails_low_retention() -> None:
    failures = check_thresholds({"resident_retention": 0.1})
    assert len(failures) == 1
    assert failures[0].metric == "resident_retention"


def test_fish_targets_reference_vadi() -> None:
    refs = {
        ref
        for spec in FISH_TARGETS.values()
        for ref in spec.get("references", [])
    }
    assert any("VADI" in r for r in refs)
    assert any("#25" in r for r in refs)


def test_evaluate_fish_metrics_sample_hdf5(sample_hdf5: Path) -> None:
    metrics = evaluate_fish_metrics(sample_hdf5)
    assert metrics["immunity_hipr_detectable"] == pytest.approx(0.0)
    assert metrics["immunity_hcr_detectable"] == pytest.approx(1.0)
    assert metrics["plasmid_dna_fish_detection_fraction"] > 0.5
    failures = check_fish_targets(metrics)
    assert failures == []


def test_invariants_reject_sample_geometry(sample_hdf5: Path) -> None:
    metrics = evaluate_metrics(sample_hdf5)
    fish_metrics = evaluate_fish_metrics(sample_hdf5)
    with GutIBMData(sample_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert {failure.metric for failure in failures} == {
        "nnd_mean", "mean_exclusion_radius",
    }


def test_evaluate_metrics_requires_multiple_steps(single_step_hdf5: Path) -> None:
    with pytest.raises(ValueError, match="at least two"):
        evaluate_metrics(single_step_hdf5)
