"""
HDF5 reader for GutIBM output – compatible with nufeb_tools format.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np


class GutIBMData:
    """Lazy-loading wrapper around a GutIBM HDF5 output file."""

    def __init__(self, filepath: str | Path) -> None:
        self.path = Path(filepath)
        self._file: h5py.File | None = None

    def open(self) -> None:
        self._file = h5py.File(self.path, "r")

    def close(self) -> None:
        if self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self) -> "GutIBMData":
        self.open()
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    @property
    def steps(self) -> list[str]:
        """Return sorted list of step group names."""
        assert self._file is not None
        return sorted(
            [k for k in self._file.keys() if k.startswith("step_")],
            key=lambda s: int(s.split("_")[1]),
        )

    def get_agents(self, step: str) -> dict[str, np.ndarray]:
        """Return agent arrays for a given step."""
        assert self._file is not None
        grp = self._file[f"{step}/atoms"]
        return {name: np.array(ds) for name, ds in grp.items()}

    def get_grid(self, step: str) -> dict[str, np.ndarray]:
        """Return chemical grid arrays for a given step."""
        assert self._file is not None
        grp = self._file[f"{step}/grid"]
        return {name: np.array(ds) for name, ds in grp.items()}

    def get_metadata(self, step: str) -> dict[str, Any]:
        """Return metadata scalars for a given step."""
        assert self._file is not None
        grp = self._file[f"{step}/metadata"]
        return {name: np.array(ds).item() for name, ds in grp.items()}

    def get_lineage(self, step: str) -> dict[str, np.ndarray]:
        """Return lineage-tracking arrays for a given step."""
        assert self._file is not None
        grp = self._file[f"{step}/lineage"]
        return {name: np.array(ds) for name, ds in grp.items()}

    def time_series(self, field: str = "num_agents") -> tuple[np.ndarray, np.ndarray]:
        """Extract a scalar time series from metadata."""
        times, vals = [], []
        for step in self.steps:
            meta = self.get_metadata(step)
            times.append(meta["time"])
            vals.append(meta.get(field, 0))
        return np.array(times), np.array(vals)
