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


def test_invariants_accept_self_consistent_geometry(invariant_hdf5: Path) -> None:
    """A self-consistent synthetic output satisfies all structural invariants."""
    metrics = evaluate_metrics(invariant_hdf5)
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert failures == []


def test_invariants_document_sample_geometry_mismatch(sample_hdf5: Path) -> None:
    """The fixture's agents exceed its declared grid, triggering geometry checks."""
    metrics = evaluate_metrics(sample_hdf5)
    fish_metrics = evaluate_fish_metrics(sample_hdf5)
    with GutIBMData(sample_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert {failure.metric for failure in failures} == {
        "nnd_mean", "mean_exclusion_radius",
    }


def test_invariants_reject_nan_metric(invariant_hdf5: Path) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    metrics["hopkins_statistic"] = float("nan")
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert any(failure.metric == "hopkins_statistic" for failure in failures)


def test_invariants_reject_missing_metric(invariant_hdf5: Path) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    del metrics["hopkins_statistic"]
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert any(failure.metric == "hopkins_statistic" for failure in failures)


def test_invariants_reject_out_of_range_fraction(invariant_hdf5: Path) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    fish_metrics["detection_fraction"] = 1.1
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert any(failure.metric == "detection_fraction" for failure in failures)


def test_invariants_reject_zero_resident_retention(invariant_hdf5: Path) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    metrics["resident_retention"] = 0.0
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert any(failure.metric == "resident_retention" for failure in failures)


def test_invariants_reject_full_resident_retention(invariant_hdf5: Path) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    metrics["resident_retention"] = 1.0
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert any(failure.metric == "resident_retention" for failure in failures)


def test_invariants_reject_excess_detected_agents(invariant_hdf5: Path) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    fish_metrics["n_detected"] = 13.0
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert any(failure.metric == "n_detected" for failure in failures)


def test_invariants_reject_lower_hcr_detection(invariant_hdf5: Path) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    fish_metrics["immunity_hcr_detectable"] = 0.0
    fish_metrics["immunity_hipr_detectable"] = 1.0
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert any(failure.metric == "immunity_detection" for failure in failures)


def test_invariants_reject_lower_rrna_detection(invariant_hdf5: Path) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    fish_metrics["rrna_phylogroup_detection_fraction"] = 0.1
    fish_metrics["immunity_mrna_detection_fraction"] = 0.2
    with GutIBMData(invariant_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert any(failure.metric == "rRNA_detection" for failure in failures)


def test_invariants_reject_missing_summary_intervals(
    invariant_hdf5: Path,
    no_summary_hdf5: Path,
) -> None:
    metrics = evaluate_metrics(invariant_hdf5)
    fish_metrics = evaluate_fish_metrics(invariant_hdf5)
    with GutIBMData(no_summary_hdf5) as data:
        failures = check_invariants(data, metrics, fish_metrics)
    assert len(failures) == 1
    assert failures[0].metric == "summary"
    assert "at least one summary interval" in failures[0].message


def test_evaluate_metrics_requires_multiple_steps(single_step_hdf5: Path) -> None:
    with pytest.raises(ValueError, match="at least two"):
        evaluate_metrics(single_step_hdf5)
