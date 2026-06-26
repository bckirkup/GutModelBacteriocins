"""Tests for EARI/VADI validation regression (issue #56)."""

from __future__ import annotations

from pathlib import Path

import pytest

from gut_ibm_tools.validation_regression import (
    FISH_TARGETS,
    VALIDATION_TARGETS,
    check_fish_targets,
    check_thresholds,
    compare_fish_golden,
    compare_golden,
    evaluate_fish_metrics,
    evaluate_metrics,
    load_golden,
    run_validation,
)


FIXTURES = Path(__file__).parent / "fixtures"


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


def test_compare_golden_exact_match() -> None:
    golden = load_golden(FIXTURES / "sample_hdf5_golden.json")
    metrics = golden["metrics"]
    assert compare_golden(metrics, golden) == []


def test_compare_golden_detects_drift() -> None:
    golden = load_golden(FIXTURES / "sample_hdf5_golden.json")
    metrics = dict(golden["metrics"])
    metrics["monochromatic_score"] += 0.01
    failures = compare_golden(metrics, golden)
    assert any(f.metric == "monochromatic_score" for f in failures)


def test_run_validation_golden_file(sample_hdf5: Path) -> None:
    golden_path = FIXTURES / "sample_hdf5_golden.json"
    metrics, fish_metrics, failures = run_validation(sample_hdf5, golden_path=golden_path)
    assert failures == []
    assert fish_metrics is None
    assert "hopkins_statistic" in metrics


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
    assert metrics["immunity_hipr_detectable"] == 0.0
    assert metrics["immunity_hcr_detectable"] == 1.0
    assert metrics["plasmid_dna_fish_detection_fraction"] > 0.5
    failures = check_fish_targets(metrics)
    assert failures == []


def test_compare_fish_golden_exact_match() -> None:
    golden = load_golden(FIXTURES / "eari_vadi_ci_fish_golden.json")
    metrics = golden["metrics"]
    assert compare_fish_golden(metrics, golden) == []


def test_run_validation_with_fish_golden(sample_hdf5: Path) -> None:
    fish_golden_path = FIXTURES / "sample_hdf5_fish_golden.json"
    _, fish_metrics, failures = run_validation(
        sample_hdf5,
        fish_golden_path=fish_golden_path,
        enforce_fish_targets=True,
    )
    assert fish_metrics is not None
    assert failures == []


def test_ci_golden_fixture_structure() -> None:
    golden = load_golden(FIXTURES / "eari_vadi_ci_golden.json")
    assert golden["scenario"] == "eari_vadi_ci"
    assert "metrics" in golden
    assert "resident_retention" in golden["metrics"]
    assert "comet_tail_ratio" in golden["metrics"]


def test_evaluate_metrics_requires_multiple_steps(single_step_hdf5: Path) -> None:
    with pytest.raises(ValueError, match="at least two"):
        evaluate_metrics(single_step_hdf5)
