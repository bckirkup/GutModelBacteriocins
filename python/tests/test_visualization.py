"""Smoke tests for gut_ibm_tools.visualization (matplotlib optional)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from gut_ibm_tools import GutIBMData, visualization


def test_plot_agent_positions_smoke(sample_hdf5: Path, tmp_path: Path) -> None:
    out = tmp_path / "agents.png"
    with GutIBMData(sample_hdf5) as data:
        visualization.plot_agent_positions(
            data, "step_000000", projection="xy", output=str(out)
        )
    assert out.exists()
    assert out.stat().st_size > 0
