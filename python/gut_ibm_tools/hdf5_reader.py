"""
HDF5 reader for GutIBM output (Spec 4 layered schema).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Self

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

    def __enter__(self) -> Self:
        self.open()
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    @property
    def steps(self) -> list[str]:
        """Return sorted list of step group names from summary or agents layers."""
        assert self._file is not None
        layer = "summary" if "summary" in self._file else "agents"
        return self.steps_for(layer)

    def steps_for(self, layer: str) -> list[str]:
        """Return sorted step group names present in one output layer."""
        assert self._file is not None
        if layer not in self._file:
            return []
        return sorted(
            [k for k in self._file[layer] if k.startswith("step_")],
            key=lambda s: int(s.split("_")[1]),
        )

    @property
    def domain_diagonal(self) -> float:
        """Return the physical domain diagonal from HDF5 metadata."""
        assert self._file is not None
        legacy_dx = float(self._file.attrs.get("grid_dx", 0.0))
        grid_dx = tuple(
            float(self._file.attrs.get(name, legacy_dx))
            for name in ("grid_dx_x", "grid_dx_y", "grid_dx_z")
        )
        lengths = []
        for axis, cells in enumerate(
            (("x", self._nx), ("y", self._ny), ("z", self._nz))
        ):
            axis_name, cell_count = cells
            lo = self._file.attrs.get(f"domain_lo_{axis_name}")
            hi = self._file.attrs.get(f"domain_hi_{axis_name}")
            if lo is not None and hi is not None:
                lengths.append(float(hi) - float(lo))
            else:
                lengths.append(float(cell_count) * grid_dx[axis])
        return float(np.linalg.norm(lengths))

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
        if "stocks" in grp:
            out["stocks"] = {name: np.array(ds).item() for name, ds in grp["stocks"].items()}
        if "mechanics" in grp:
            out["mechanics"] = {
                name: np.array(ds).item() for name, ds in grp["mechanics"].items()
            }
        for name in ("halt_reason_code", "halt_density_cells_per_mL"):
            if name in grp:
                out[name] = read_scalar(name)
        return out

    def get_run_provenance(self) -> dict[str, Any]:
        """Return run-level provenance and termination metadata."""
        assert self._file is not None
        if "run_provenance" not in self._file:
            return {}
        group = self._file["run_provenance"]
        result: dict[str, Any] = {}
        for name, dataset in group.items():
            value = np.array(dataset)
            value = value.item() if value.ndim == 0 else value
            if isinstance(value, bytes):
                value = value.decode("utf-8")
            result[name] = value
        return result

    def get_termination_audit(self) -> dict[str, Any]:
        """Return the authoritative horizon and termination assessment."""
        provenance = self.get_run_provenance()
        code = int(provenance.get("termination_cause_code", 5))
        cause = provenance.get("termination_cause", "incomplete_unknown")
        detail = provenance.get("termination_detail", "")
        if isinstance(cause, bytes):
            cause = cause.decode("utf-8")
        if isinstance(detail, bytes):
            detail = detail.decode("utf-8")
        return {
            "reached_horizon": code == 0,
            "cause_code": code,
            "cause": cause,
            "detail": detail,
            "step": provenance.get("termination_step"),
            "time": provenance.get("termination_time"),
        }

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
        metadata = {
            "time": summary["time"],
            "step": summary["step"],
            "num_agents": summary["num_agents"],
            "num_lineages": summary.get("num_lineages", 0),
        }
        for name in ("halt_reason_code", "halt_density_cells_per_mL"):
            if name in summary:
                metadata[name] = summary[name]
        return metadata

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

    def get_genome_loci(
        self,
        step: str,
    ) -> dict[int, tuple[tuple[int, int], ...]]:
        """Return BI-locus identities grouped by agent ID when offsets exist.

        Older genome layers lack the per-agent ID and offset datasets.  Those
        layers return an empty mapping so callers can use count-based fallback
        genotypes.
        """
        genome = self.get_genome(step)
        required = {"id", "bi_offset", "bi_count", "bi_toxin_id", "bi_immunity_id"}
        if not required.issubset(genome):
            return {}
        ids = genome["id"]
        offsets = genome["bi_offset"]
        counts = genome["bi_count"]
        toxin_ids = genome["bi_toxin_id"]
        immunity_ids = genome["bi_immunity_id"]
        if len(ids) != len(offsets) or len(ids) != len(counts):
            return {}
        if len(set(ids.tolist())) != len(ids):
            return {}
        if len(toxin_ids) != len(immunity_ids):
            return {}
        loci: dict[int, tuple[tuple[int, int], ...]] = {}
        previous_end = 0
        for agent_id, offset, count in zip(ids, offsets, counts):
            start = int(offset)
            stop = start + int(count)
            if (
                start < previous_end
                or count < 0
                or stop > len(toxin_ids)
                or stop > len(immunity_ids)
            ):
                return {}
            loci[int(agent_id)] = tuple(
                sorted(
                    (int(toxin), int(immunity))
                    for toxin, immunity in zip(
                        toxin_ids[start:stop], immunity_ids[start:stop]
                    )
                )
            )
            previous_end = stop
        if previous_end != len(toxin_ids):
            return {}
        return loci

    def time_series(self, field: str = "num_agents") -> tuple[np.ndarray, np.ndarray]:
        """Extract a scalar time series from summary layers."""
        times, vals = [], []
        for step in self.steps:
            meta = self.get_metadata(step)
            times.append(meta["time"])
            vals.append(meta.get(field, 0))
        return np.array(times), np.array(vals)
