"""Tests for gut_ibm_tools.hdf5_reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from gut_ibm_tools import GutIBMData


def test_context_manager_opens_and_closes(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        assert data._file is not None
        assert len(data.steps) == 2
    assert data._file is None


def test_steps_sorted_numerically(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        assert data.steps == ["step_000000", "step_000001"]


def test_steps_for_returns_layer_specific_steps(mismatched_schedule_hdf5: Path) -> None:
    with GutIBMData(mismatched_schedule_hdf5) as data:
        assert data.steps == ["step_000000", "step_000001", "step_000002"]
        assert data.steps_for("agents") == ["step_000000", "step_000001"]
        assert data.steps_for("summary") == [
            "step_000000",
            "step_000001",
            "step_000002",
        ]


def test_get_agents_schema(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        agents = data.get_agents("step_000000")
        assert set(agents) >= {"x", "y", "z", "type", "lineage"}
        assert len(agents["x"]) == 12
        assert agents["type"].dtype == np.int32


def test_get_grid_and_metadata(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        grid = data.get_grid("step_000000")
        assert "bacteriocin_BtuB" in grid
        assert len(grid["bacteriocin_BtuB"]) == 20

        meta = data.get_metadata("step_000000")
        assert meta["time"] == pytest.approx(0.0)
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
        assert data.get_genome_loci("step_000000") == {}


def test_get_genome_loci_reconstructs_ragged_agents(ragged_genome_hdf5: Path) -> None:
    with GutIBMData(ragged_genome_hdf5) as data:
        loci = data.get_genome_loci("step_000000")
    assert loci == {
        101: (),
        205: ((7, 17), (11, 19)),
        999: ((13, 23),),
    }


def test_time_series_num_agents(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        times, counts = data.time_series("num_agents")
        np.testing.assert_allclose(times, [0.0, 3600.0])
        np.testing.assert_array_equal(counts, [12, 12])
