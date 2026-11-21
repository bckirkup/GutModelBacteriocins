"""Tests for gut_ibm_tools.validation metrics."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
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


@pytest.mark.parametrize("dataset_name", ["bacteriocin_BtuB", "bacteriocin_lumped"])
def test_validate_spatial_signatures_accepts_toxin_dataset_names(
    sample_hdf5: Path,
    dataset_name: str,
) -> None:
    if dataset_name == "bacteriocin_lumped":
        with h5py.File(sample_hdf5, "a") as handle:
            for step in handle["grid"].values():
                values = step["bacteriocin_BtuB"][:]
                del step["bacteriocin_BtuB"]
                step.create_dataset(dataset_name, data=values)

    with GutIBMData(sample_hdf5) as data:
        result = validation.validate_spatial_signatures(data, "step_000000")

    assert result["comet_tail_ratio"] > 1.0


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


def test_compute_lineage_retention_skewed_vs_unweighted() -> None:
    # Initial founder distribution: lineage 1 has 90% of cells, lineages 2..10 have ~1.1% each.
    # Final state: lineage 1 survives at 95% share; lineages 2..9 go extinct; lineage 10 has 1/200 cells (0.5%).
    first_ids = np.array([1] * 90 + list(range(2, 11)), dtype=np.int32)
    last_ids = np.array([1] * 190 + [10] * 1 + [99] * 9, dtype=np.int32)

    no_lod = validation.compute_lineage_retention_metrics(first_ids, last_ids, lod_relative_abundance=0.0)
    # Unweighted: 2 of 10 founder lineages survive (0.20), whereas initial-abundance-weighted
    # survival is (90 + 1) / 99 = 91.9%!
    assert no_lod["resident_retention"] == pytest.approx(0.20, abs=1e-9)
    assert no_lod["founder_lineage_persistence"] == pytest.approx(0.20, abs=1e-9)
    assert no_lod["abundance_weighted_retention"] == pytest.approx(91.0 / 99.0, rel=1e-6)
    assert no_lod["final_resident_abundance_fraction"] == pytest.approx(191.0 / 200.0, rel=1e-6)

    # Under a 1% metagenomic/pellet LOD, lineage 10 (0.5%) falls below detection -> 1/10 unweighted, 90/99 weighted
    with_lod = validation.compute_lineage_retention_metrics(first_ids, last_ids, lod_relative_abundance=0.01)
    assert with_lod["resident_retention"] == pytest.approx(0.10, abs=1e-9)
    assert with_lod["abundance_weighted_retention"] == pytest.approx(90.0 / 99.0, rel=1e-6)
    assert with_lod["final_resident_abundance_fraction"] == pytest.approx(190.0 / 200.0, rel=1e-6)


def test_validate_spatial_signatures_optical_depth_window(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        agents = data.get_agents("step_000000")
        z_med = float(sum(agents["z"]) / len(agents["z"]))
        z_max = float(max(agents["z"]))
        sliced = validation.validate_spatial_signatures(
            data, "step_000000", depth_window_m=(z_med, z_max + 1e-6)
        )
    assert "monochromatic_score" in sliced
    assert sliced["comet_tail_ratio"] > 0.0
