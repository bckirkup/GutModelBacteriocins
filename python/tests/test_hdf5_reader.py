"""Tests for gut_ibm_tools.hdf5_reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from gut_ibm_tools import GutIBMData


def test_context_manager_opens_and_closes(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        assert data._file is not None
        assert len(data.steps) == 2
    assert data._file is None


def test_steps_sorted_numerically(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        assert data.steps == ["step_000000", "step_000001"]


def test_get_agents_schema(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        agents = data.get_agents("step_000000")
        assert set(agents) >= {"x", "y", "z", "type", "lineage"}
        assert len(agents["x"]) == 12
        assert agents["type"].dtype == np.int32


def test_get_grid_and_metadata(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        grid = data.get_grid("step_000000")
        assert "bacteriocin" in grid
        assert len(grid["bacteriocin"]) == 20

        meta = data.get_metadata("step_000000")
        assert meta["time"] == 0.0
        assert meta["num_agents"] == 12


def test_get_lineage(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        lineage = data.get_lineage("step_000000")
        assert "num_bi_loci" in lineage
        assert len(lineage["btuB_expression"]) == 12


def test_get_genome(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        genome = data.get_genome("step_000000")
        assert "has_conjugative_plasmid" in genome
        assert len(genome["has_conjugative_plasmid"]) == 12


def test_time_series_num_agents(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        times, counts = data.time_series("num_agents")
        np.testing.assert_allclose(times, [0.0, 3600.0])
        np.testing.assert_array_equal(counts, [12, 12])
