"""Validate GutIBM HDF5 writer output with the Python reader."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from gut_ibm_tools import GutIBMData


def _roundtrip_binary() -> Path | None:
    build = Path(os.environ.get("GUTIBM_BUILD_DIR", "build"))
    exe = build / "tests" / "test_hdf5_roundtrip"
    return exe if exe.exists() else None


@pytest.fixture(scope="module")
def writer_snapshot(tmp_path_factory: pytest.TempPathFactory) -> Path:
    exe = _roundtrip_binary()
    if exe is None:
        pytest.skip("test_hdf5_roundtrip binary not built")

    out_dir = tmp_path_factory.mktemp("hdf5_roundtrip")
    h5_path = out_dir / "serial.h5"
    env = os.environ.copy()
    env["GUTIBM_ROUNDTRIP_H5"] = str(h5_path)
    subprocess.run([str(exe)], check=True, env=env)
    return h5_path


def test_python_reader_roundtrip_schema(writer_snapshot: Path) -> None:
    with GutIBMData(writer_snapshot) as data:
        assert data.steps == ["step_000000", "step_000001", "step_000002"]

        meta0 = data.get_metadata("step_000000")
        assert meta0["step"] == 0
        assert meta0["time"] == pytest.approx(0.0)
        assert meta0["num_agents"] == 12

        agents2 = data.get_agents("step_000002")
        meta2 = data.get_metadata("step_000002")
        assert len(agents2["x"]) == meta2["num_agents"]
        assert meta2["step"] == 2
        assert meta2["time"] == pytest.approx(120.0)

        grid2 = data.get_grid("step_000002")
        assert "bacteriocin" in grid2
        assert "carbon" in grid2
        assert len(grid2["bacteriocin"]) > 0

        lineage2 = data.get_lineage("step_000002")
        assert len(lineage2["num_bi_loci"]) == len(agents2["x"])

        times, counts = data.time_series("num_agents")
        np.testing.assert_allclose(times, [0.0, 60.0, 120.0])
        assert counts[0] == 12
        assert counts[-1] == meta2["num_agents"]
