"""Shared pytest fixtures for gut_ibm_tools."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest


def write_sample_hdf5(path: Path, *, n_agents: int = 12, n_steps: int = 2) -> None:
    """Write a minimal GutIBM-compatible HDF5 file for tests."""
    rng = np.random.default_rng(42)

    with h5py.File(path, "w") as f:
        for step_idx in range(n_steps):
            step_name = f"step_{step_idx:06d}"
            step_grp = f.create_group(step_name)

            # Agents: two spatially separated monochromatic patches.
            n_per_type = n_agents // 2
            types = np.array([1] * n_per_type + [2] * n_per_type, dtype=np.int32)
            x = np.concatenate([
                rng.uniform(0, 20e-6, n_per_type),
                rng.uniform(80e-6, 100e-6, n_per_type),
            ])
            y = rng.uniform(0, 20e-6, n_agents)
            z = rng.uniform(0, 10e-6, n_agents)

            atoms = step_grp.create_group("atoms")
            atoms.create_dataset("id", data=np.arange(n_agents, dtype=np.int64))
            atoms.create_dataset("type", data=types)
            atoms.create_dataset("state", data=np.zeros(n_agents, dtype=np.int32))
            atoms.create_dataset("x", data=x)
            atoms.create_dataset("y", data=y)
            atoms.create_dataset("z", data=z)
            atoms.create_dataset("radius", data=np.full(n_agents, 0.5e-6))
            atoms.create_dataset("biomass", data=np.full(n_agents, 1e-15))
            atoms.create_dataset("mu", data=np.full(n_agents, 5e-4))
            # Lineages 1–3 at step 0; lineage 3 lost by final step.
            lineage_ids = np.array([1, 1, 2, 2, 3, 3] * (n_agents // 6), dtype=np.int64)
            if len(lineage_ids) < n_agents:
                lineage_ids = np.resize(lineage_ids, n_agents)
            if step_idx == n_steps - 1 and n_steps > 1:
                lineage_ids[lineage_ids == 3] = 99  # transient immigrant lineage
            atoms.create_dataset("lineage", data=lineage_ids[:n_agents])

            # Grid: downstream-heavy bacteriocin profile (comet-tail analog).
            ncells = 20
            grid = step_grp.create_group("grid")
            bacteriocin = np.linspace(0.1, 2.0, ncells) ** (1 + 0.2 * step_idx)
            grid.create_dataset("bacteriocin", data=bacteriocin)

            meta = step_grp.create_group("metadata")
            meta.create_dataset("time", data=np.array(step_idx * 3600.0))
            meta.create_dataset("step", data=np.array(step_idx, dtype=np.int32))
            meta.create_dataset("num_agents", data=np.array(n_agents, dtype=np.int32))
            meta.create_dataset("num_lineages", data=np.array(3, dtype=np.int32))

            lin = step_grp.create_group("lineage")
            lin.create_dataset("btuB_expression", data=rng.uniform(0.2, 1.0, n_agents))
            lin.create_dataset("fepA_expression", data=rng.uniform(0.2, 1.0, n_agents))
            resident_mask = np.isin(lineage_ids[:n_agents], [1, 2])
            n_bi = np.where(resident_mask, 2, 0).astype(np.int32)
            lin.create_dataset("num_bi_loci", data=n_bi)
            lin.create_dataset("generation", data=np.zeros(n_agents, dtype=np.int32))


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
