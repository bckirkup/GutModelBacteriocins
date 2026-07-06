"""
HDF5 reader for GutIBM output (Spec 4 layered schema).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from .path_utils import validate_input_path


class GutIBMData:
    """Lazy-loading wrapper around a GutIBM HDF5 output file."""

    def __init__(self, filepath: str | Path) -> None:
        self.path = validate_input_path(filepath)
        self._file: h5py.File | None = None
        self._nx: int = 0
        self._ny: int = 0
        self._nz: int = 0

    def open(self) -> None:
        self._file = h5py.File(self.path, "r")
        self._nx = int(self._file.attrs.get("nx", 0))
        self._ny = int(self._file.attrs.get("ny", 0))
        self._nz = int(self._file.attrs.get("nz", 0))

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
        """Return sorted list of step group names from summary or agents layers."""
        assert self._file is not None
        layer = "summary" if "summary" in self._file else "agents"
        if layer not in self._file:
            return []
        return sorted(
            [k for k in self._file[layer].keys() if k.startswith("step_")],
            key=lambda s: int(s.split("_")[1]),
        )

    def get_summary(self, step: str) -> dict[str, Any]:
        """Return summary scalars and nested groups for a step."""
        assert self._file is not None
        grp = self._file[f"summary/{step}"]

        def read_scalar(name: str) -> Any:
            return np.array(grp[name]).item()

        out: dict[str, Any] = {
            "time": read_scalar("time"),
            "step": read_scalar("step"),
            "n_total": read_scalar("n_total") if "n_total" in grp else read_scalar("num_agents"),
            "num_agents": read_scalar("num_agents") if "num_agents" in grp else read_scalar("n_total"),
        }
        if "num_lineages" in grp:
            out["num_lineages"] = read_scalar("num_lineages")
        if "events" in grp:
            out["events"] = {name: np.array(ds).item() for name, ds in grp["events"].items()}
        if "chem" in grp:
            out["chem"] = {name: np.array(ds).item() for name, ds in grp["chem"].items()}
        if "spatial" in grp:
            out["spatial"] = {name: np.array(ds).item() for name, ds in grp["spatial"].items()}
        return out

    def get_agents(self, step: str) -> dict[str, np.ndarray]:
        """Return agent arrays for a given step."""
        assert self._file is not None
        grp = self._file[f"agents/{step}"]
        out = {name: np.array(ds) for name, ds in grp.items()}
        if "mu_realized" in out and "mu" not in out:
            out["mu"] = out["mu_realized"]
        if "lineage_id" in out and "lineage" not in out:
            out["lineage"] = out["lineage_id"]
        return out

    def get_grid(self, step: str) -> dict[str, np.ndarray]:
        """Return chemical grid arrays for a step (3D datasets flattened to 1D)."""
        assert self._file is not None
        path = f"grid/{step}"
        if path not in self._file:
            return {}
        grp = self._file[path]
        out: dict[str, np.ndarray] = {}
        for name, ds in grp.items():
            arr = np.array(ds)
            out[name] = arr.ravel()
        return out

    def get_metadata(self, step: str) -> dict[str, Any]:
        """Compatibility alias for summary scalars."""
        summary = self.get_summary(step)
        return {
            "time": summary["time"],
            "step": summary["step"],
            "num_agents": summary["num_agents"],
            "num_lineages": summary.get("num_lineages", 0),
        }

    def get_lineage(self, step: str) -> dict[str, np.ndarray]:
        """Return lineage-tracking arrays for a given step."""
        assert self._file is not None
        grp = self._file[f"lineage/{step}"]
        return {name: np.array(ds) for name, ds in grp.items()}

    def get_genome(self, step: str) -> dict[str, np.ndarray]:
        """Return per-agent genome arrays for a given step (if present)."""
        assert self._file is not None
        path = f"genome/{step}"
        if path not in self._file:
            return {}
        grp = self._file[path]
        return {name: np.array(ds) for name, ds in grp.items()}

    def time_series(self, field: str = "num_agents") -> tuple[np.ndarray, np.ndarray]:
        """Extract a scalar time series from summary layers."""
        times, vals = [], []
        for step in self.steps:
            meta = self.get_metadata(step)
            times.append(meta["time"])
            vals.append(meta.get(field, 0))
        return np.array(times), np.array(vals)
