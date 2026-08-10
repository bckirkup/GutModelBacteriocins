"""Tests for colony-scale observables."""

from __future__ import annotations

import numpy as np
import pytest
from gut_ibm_tools.colony import ColonyConfig, build_colony_catalog, dbscan_colonies
from gut_ibm_tools.hdf5_reader import GutIBMData


def _agents(points: np.ndarray, ids: np.ndarray | None = None) -> dict[str, np.ndarray]:
    count = len(points)
    return {
        "id": np.arange(count, dtype=np.int64) if ids is None else ids,
        "type": np.zeros(count, dtype=np.int64),
        "n_bi_loci": np.ones(count, dtype=np.int64),
        "lineage_id": np.zeros(count, dtype=np.int64),
        "x": points[:, 0], "y": points[:, 1], "z": points[:, 2],
    }


def test_blobs_recover_size_centroid_and_neighbour_distance() -> None:
    rng = np.random.default_rng(8)
    truth = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    points = np.vstack([center + rng.normal(0.0, 0.15, (20, 3)) for center in truth])
    catalog = build_colony_catalog(_agents(points), ColonyConfig(eps=0.6, min_samples=3))
    assert len(catalog.colonies["colony_id"]) == 2
    np.testing.assert_allclose(np.sort(catalog.colonies["n_members"]), [20, 20])
    assert np.min(catalog.colonies["nn_colony_distance"]) == pytest.approx(10.0, abs=0.4)


def test_ids_and_noise_round_trip() -> None:
    points = np.vstack([np.zeros((4, 3)), [[100.0, 100.0, 100.0]]])
    ids = np.array([91, 4, 77, 12, 300], dtype=np.int64)
    catalog = build_colony_catalog(_agents(points, ids), ColonyConfig(eps=1.0, min_samples=2))
    np.testing.assert_array_equal(catalog.agents["agent_id"], ids)
    assert catalog.agents["colony_id"][-1] == -1
    assert len(catalog.colonies["colony_id"]) == 1


def test_purity_and_mixed_genotypes() -> None:
    points = np.zeros((4, 3))
    agents = _agents(points)
    agents["lineage_id"] = np.array([1, 1, 2, 2])
    catalog = build_colony_catalog(agents, ColonyConfig(eps=1.0, min_samples=2))
    assert catalog.colonies["n_genotypes"][0] == 2
    assert catalog.colonies["dominant_genotype_fraction"][0] == pytest.approx(0.5)
    assert not catalog.colonies["is_clonal"][0]
    agents["lineage_id"][:] = 1
    clonal = build_colony_catalog(agents, ColonyConfig(eps=1.0, min_samples=2))
    assert clonal.colonies["dominant_genotype_fraction"][0] == pytest.approx(1.0)
    assert clonal.colonies["is_clonal"][0]


def test_dbscan_border_points_are_attached_to_core() -> None:
    points = np.array([[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.4, 0.0], [0.4, 0.4, 0.0], [2.0, 0.0, 0.0]])
    labels = dbscan_colonies(points, eps=0.5, min_samples=4)
    assert np.all(labels[:4] == labels[0])
    assert labels[-1] == -1


def test_catalog_reads_sample_hdf5(sample_hdf5) -> None:
    with GutIBMData(sample_hdf5) as data:
        catalog = build_colony_catalog(
            data.get_agents("step_000000"),
            ColonyConfig(eps=30e-6, min_samples=2),
        )
    assert len(catalog.agents["agent_id"]) == 12
    assert catalog.agents["agent_id"].dtype == np.int64


def test_eps_sensitivity_reduces_noise_and_merges_at_largest_factor() -> None:
    rng = np.random.default_rng(13)
    points = np.vstack([
        rng.normal([0.0, 0.0, 0.0], 0.08, (20, 3)),
        rng.normal([1.5, 0.0, 0.0], 0.08, (20, 3)),
    ])
    catalog = build_colony_catalog(
        _agents(points),
        ColonyConfig(eps=0.6, min_samples=3),
    )
    np.testing.assert_array_equal(
        np.diff(catalog.eps_sensitivity["noise_fraction"]) <= 0,
        np.ones(3, dtype=bool),
    )
    assert catalog.eps_sensitivity["n_colonies"][-1] == 1


def test_producer_threshold_flags() -> None:
    points = np.zeros((120, 3))
    agents = _agents(points)
    catalog = build_colony_catalog(agents, ColonyConfig(eps=1.0, min_samples=2))
    assert catalog.colonies["reaches_113"][0]
    assert not catalog.colonies["reaches_527"][0]
    assert not catalog.colonies["reaches_1361"][0]


def test_all_noise_has_well_formed_empty_colony_table() -> None:
    points = np.arange(15, dtype=float).reshape(5, 3) * 100.0
    catalog = build_colony_catalog(_agents(points), ColonyConfig(eps=1.0, min_samples=2))
    assert catalog.colonies["colony_id"].size == 0
    assert catalog.colonies["colony_id"].dtype == np.int64
    assert catalog.colonies["nn_colony_distance"].dtype == float
