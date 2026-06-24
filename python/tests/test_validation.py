"""Tests for gut_ibm_tools.validation metrics."""

from __future__ import annotations

from pathlib import Path

import pytest

from gut_ibm_tools import GutIBMData, validation


def test_validate_spatial_signatures_keys(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        result = validation.validate_spatial_signatures(data, "step_000000")

    expected_keys = {
        "monochromatic_score",
        "comet_tail_ratio",
        "mean_exclusion_radius",
        "hopkins_statistic",
        "nnd_mean",
        "comet_tail_asymmetry",
    }
    assert expected_keys <= set(result)
    assert result["monochromatic_score"] > 0.7
    assert result["comet_tail_ratio"] > 1.0
    assert result["mean_exclusion_radius"] > 0.0


def test_validate_genomic_signatures_retention(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        result = validation.validate_genomic_signatures(data)

    assert "resident_retention" in result
    # Lineages 1 and 2 persist; lineage 3 is replaced at the final step.
    assert result["resident_retention"] == pytest.approx(2 / 3, rel=1e-3)
    assert result["resident_mean_bi_loci"] > result["transient_mean_bi_loci"]


def test_validate_genomic_signatures_single_step(single_step_hdf5: Path) -> None:
    with GutIBMData(single_step_hdf5) as data:
        result = validation.validate_genomic_signatures(data)

    assert result == {"error": -1.0}
