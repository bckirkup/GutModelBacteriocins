"""Tests for gut_ibm_tools.hdf5_reader."""

from __future__ import annotations

import shutil
from pathlib import Path

import h5py
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
        assert meta["halt_reason_code"] == 0
        summary = data.get_summary("step_000001")
        assert summary["stocks"]["bacteriostatic_live_agents"] == 1
        assert summary["stocks"]["washout_trapped_live_agents"] == 1
        assert summary["halt_reason_code"] == 0
        assert summary["halt_density_cells_per_mL"] == pytest.approx(0.0)
        provenance = data.get_run_provenance()
        assert provenance["completed_total_time"] == 1


def test_get_termination_audit_horizon(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        audit = data.get_termination_audit()
    assert audit["reached_horizon"] is True
    assert audit["cause_code"] == 0
    assert audit["cause"] == "horizon_reached"
    assert audit["detail"] == "test horizon"
    assert audit["step"] == 1
    assert audit["time"] == pytest.approx(3600.0)


def test_get_termination_audit_non_horizon_and_incomplete(
    sample_hdf5: Path, tmp_path: Path
) -> None:
    output = tmp_path / "audit.h5"
    shutil.copy2(sample_hdf5, output)
    with h5py.File(output, "r+") as handle:
        provenance = handle["run_provenance"]
        del provenance["termination_cause_code"]
        provenance.create_dataset(
            "termination_cause_code", data=np.array(3, dtype=np.int32)
        )
        del provenance["termination_cause"]
        provenance.create_dataset("termination_cause", data="stop_requested")
        del provenance["termination_detail"]
        provenance.create_dataset(
            "termination_detail", data="stop requested before next timestep"
        )
    with GutIBMData(output) as data:
        audit = data.get_termination_audit()
    assert audit["reached_horizon"] is False
    assert audit["cause_code"] == 3
    assert audit["cause"] == "stop_requested"
    assert audit["detail"] == "stop requested before next timestep"

    with h5py.File(output, "r+") as handle:
        provenance = handle["run_provenance"]
        del provenance["termination_cause_code"]
        provenance.create_dataset(
            "termination_cause_code", data=np.array(5, dtype=np.int32)
        )
        del provenance["termination_cause"]
        provenance.create_dataset("termination_cause", data="incomplete_unknown")
    with GutIBMData(output) as data:
        audit = data.get_termination_audit()
    assert audit["reached_horizon"] is False
    assert audit["cause_code"] == 5
    assert audit["cause"] == "incomplete_unknown"


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
