"""Shared pytest fixtures for gut_ibm_tools (Spec 4 layered HDF5 schema)."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "integration: requires a built C++ test binary or long-running simulation",
    )


def write_sample_hdf5(path: Path, *, n_agents: int = 12, n_steps: int = 2) -> None:
    """Write a minimal Spec-4-compatible HDF5 file for tests."""
    rng = np.random.default_rng(42)
    nx, ny, nz = 4, 5, 1
    ncells = nx * ny * nz

    with h5py.File(path, "w") as f:
        f.attrs["gutibm_version"] = 4
        f.attrs["nx"] = nx
        f.attrs["ny"] = ny
        f.attrs["nz"] = nz
        f.attrs["grid_dx"] = 5e-6

        for step_idx in range(n_steps):
            step_name = f"step_{step_idx:06d}"

            n_per_type = n_agents // 2
            types = np.array([1] * n_per_type + [2] * n_per_type, dtype=np.int32)
            x = np.concatenate([
                rng.uniform(0, 20e-6, n_per_type),
                rng.uniform(80e-6, 100e-6, n_per_type),
            ])
            y = rng.uniform(0, 20e-6, n_agents)
            z = rng.uniform(0, 10e-6, n_agents)
            lineage_ids = np.array([1, 1, 2, 2, 3, 3] * (n_agents // 6), dtype=np.int64)
            if len(lineage_ids) < n_agents:
                lineage_ids = np.resize(lineage_ids, n_agents)
            if step_idx == n_steps - 1 and n_steps > 1:
                lineage_ids[lineage_ids == 3] = 99

            agents = f.require_group("agents").require_group(step_name)
            agents.create_dataset("id", data=np.arange(n_agents, dtype=np.int64))
            agents.create_dataset("type", data=types)
            agents.create_dataset("state", data=np.zeros(n_agents, dtype=np.int32))
            agents.create_dataset("x", data=x)
            agents.create_dataset("y", data=y)
            agents.create_dataset("z", data=z)
            agents.create_dataset("radius", data=np.full(n_agents, 0.5e-6))
            agents.create_dataset("biomass", data=np.full(n_agents, 1e-15))
            agents.create_dataset("mu_realized", data=np.full(n_agents, 5e-4))
            agents.create_dataset("lineage_id", data=lineage_ids[:n_agents])

            grid = f.require_group("grid").require_group(step_name)
            btub = np.linspace(0.1, 2.0, ncells) ** (1 + 0.2 * step_idx)
            grid.create_dataset("bacteriocin_BtuB", data=btub.reshape(nz, ny, nx))

            summary = f.require_group("summary").require_group(step_name)
            summary.create_dataset("time", data=np.array(step_idx * 3600.0))
            summary.create_dataset("step", data=np.array(step_idx, dtype=np.int32))
            summary.create_dataset("n_total", data=np.array(n_agents, dtype=np.int32))
            summary.create_dataset("num_agents", data=np.array(n_agents, dtype=np.int32))
            summary.create_dataset("num_lineages", data=np.array(3, dtype=np.int32))

            lin = f.require_group("lineage").require_group(step_name)
            lin.create_dataset("btuB_expression", data=rng.uniform(0.2, 1.0, n_agents))
            lin.create_dataset("fepA_expression", data=rng.uniform(0.2, 1.0, n_agents))
            resident_mask = np.isin(lineage_ids[:n_agents], [1, 2])
            n_bi = np.where(resident_mask, 2, 0).astype(np.int32)
            lin.create_dataset("num_bi_loci", data=n_bi)
            lin.create_dataset("generation", data=np.zeros(n_agents, dtype=np.int32))

            genome_rng = np.random.default_rng(9000 + step_idx)
            genome = f.require_group("genome").require_group(step_name)
            has_conj = np.where(resident_mask, 1, 0).astype(np.int32)
            genome.create_dataset("has_conjugative_plasmid", data=has_conj)
            genome.create_dataset(
                "plasmid_cost_amelioration",
                data=genome_rng.uniform(0.0, 0.2, n_agents),
            )
            genome.create_dataset("mutations", data=np.zeros(n_agents, dtype=np.int32))
            genome.create_dataset("parent_id", data=np.arange(n_agents, dtype=np.int64))


@pytest.fixture
def single_step_hdf5(tmp_path: Path) -> Path:
    path = tmp_path / "single_step.h5"
    write_sample_hdf5(path, n_steps=1)
    return path


@pytest.fixture
def sample_hdf5(tmp_path: Path) -> Path:
    path = tmp_path / "sample_output.h5"
    write_sample_hdf5(path)
    return path
