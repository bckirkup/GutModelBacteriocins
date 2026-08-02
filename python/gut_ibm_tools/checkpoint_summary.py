"""Extract per-checkpoint summary statistics from Spec-4 GutIBM HDF5 files."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

from . import analysis
from .hdf5_reader import GutIBMData
from .path_utils import validate_input_path

DEAD_STATE = 3
MIXING_RADIUS_M = 10e-6
MIN_AGENTS_FOR_SPATIAL = 2

# Prefer grid volume means; fall back to summary/chem scalars when grid absent.
GRID_CHEM_NAMES = (
    "carbon",
    "iron",
    "oxygen",
    "acetate",
    "ethanolamine",
    "b12",
    "mucin",
    "siderophore",
    "ai2",
    "bacteriocin_BtuB",
    "bacteriocin_FepA",
    "bacteriocin_CirA",
    "bacteriocin_FhuA",
)

EVENT_KEYS = (
    "sos_inductions",
    "phage_inductions",
    "colicin_kills",
    "cdi_kills",
    "washout_deaths",
    "boundary_deaths",
    "starvation_deaths",
    "divisions",
    "conjugation_transfers",
    "mutations",
)


def _alive_mask(agents: dict[str, np.ndarray]) -> np.ndarray:
    if "state" not in agents:
        return np.ones(len(agents.get("x", [])), dtype=bool)
    return agents["state"] != DEAD_STATE


def _positions(agents: dict[str, np.ndarray], mask: np.ndarray) -> np.ndarray:
    return np.column_stack([agents["x"][mask], agents["y"][mask], agents["z"][mask]])


def _radius_of_gyration(positions: np.ndarray) -> float:
    if len(positions) == 0:
        return float("nan")
    com = positions.mean(axis=0)
    squared = np.sum((positions - com) ** 2, axis=1)
    return float(np.sqrt(np.mean(squared)))


def _per_type_com_stats(
    positions: np.ndarray,
    types: np.ndarray,
) -> dict[str, float]:
    out: dict[str, float] = {}
    for type_id in sorted(np.unique(types)):
        tid = int(type_id)
        pts = positions[types == type_id]
        prefix = f"type{tid}"
        out[f"n_{prefix}"] = float(len(pts))
        if len(pts) == 0:
            continue
        com = pts.mean(axis=0)
        std = pts.std(axis=0)
        out[f"{prefix}_com_x"] = float(com[0])
        out[f"{prefix}_com_y"] = float(com[1])
        out[f"{prefix}_com_z"] = float(com[2])
        out[f"{prefix}_std_x"] = float(std[0])
        out[f"{prefix}_std_y"] = float(std[1])
        out[f"{prefix}_std_z"] = float(std[2])
        out[f"{prefix}_rg"] = _radius_of_gyration(pts)
        out[f"{prefix}_mean_z"] = float(pts[:, 2].mean())
    return out


def _growth_stats(agents: dict[str, np.ndarray], mask: np.ndarray) -> dict[str, float]:
    out: dict[str, float] = {}
    if "type" not in agents:
        return out
    types = agents["type"][mask]
    mu = agents.get("mu_realized", agents.get("mu"))
    mu_max = agents.get("mu_max")
    for type_id in sorted(np.unique(types)):
        tid = int(type_id)
        tmask = types == type_id
        if mu is not None:
            out[f"type{tid}_mean_mu"] = float(np.mean(mu[mask][tmask]))
        if mu_max is not None:
            out[f"type{tid}_mean_mu_max"] = float(np.mean(mu_max[mask][tmask]))
    if mu is not None and np.any(mask):
        out["mean_mu"] = float(np.mean(mu[mask]))
    return out


def _crypt_stats(agents: dict[str, np.ndarray], mask: np.ndarray) -> dict[str, float]:
    if "in_crypt" not in agents or "type" not in agents:
        return {}
    crypt = agents["in_crypt"][mask]
    types = agents["type"][mask]
    out: dict[str, float] = {"n_in_crypt": float(np.sum(crypt != 0))}
    for type_id in sorted(np.unique(types)):
        tid = int(type_id)
        tmask = types == type_id
        out[f"type{tid}_n_in_crypt"] = float(np.sum(crypt[tmask] != 0))
    return out


def _spatial_metrics(positions: np.ndarray, types: np.ndarray) -> dict[str, float]:
    out: dict[str, float] = {
        "rg_global": _radius_of_gyration(positions),
        "com_x": float(positions[:, 0].mean()) if len(positions) else float("nan"),
        "com_y": float(positions[:, 1].mean()) if len(positions) else float("nan"),
        "com_z": float(positions[:, 2].mean()) if len(positions) else float("nan"),
        "mean_z": float(positions[:, 2].mean()) if len(positions) else float("nan"),
    }
    if len(positions) < MIN_AGENTS_FOR_SPATIAL:
        return out
    out["hopkins"] = analysis.hopkins_statistic(
        positions, rng=np.random.default_rng(0)
    )
    if len(np.unique(types)) >= 2:
        out["monochromatic_score"] = analysis.monochromatic_patch_score(
            positions, types, radius=MIXING_RADIUS_M
        )
    return out


def _grid_chem_means(volumes: dict[str, np.ndarray]) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, arr in volumes.items():
        flat = np.asarray(arr, dtype=float).ravel()
        if flat.size == 0:
            continue
        out[f"grid_mean_{name}"] = float(np.mean(flat))
        out[f"grid_max_{name}"] = float(np.max(flat))
        out[f"grid_min_{name}"] = float(np.min(flat))
    return out


def _oxygen_depth_profile(volumes: dict[str, np.ndarray]) -> dict[str, float]:
    """Mean/min/max oxygen along z (axis 0 in writer layout may vary — try last axis)."""
    if "oxygen" not in volumes:
        return {}
    arr = np.asarray(volumes["oxygen"], dtype=float)
    if arr.ndim != 3:
        return {
            "oxygen_mean": float(np.mean(arr)),
            "oxygen_min": float(np.min(arr)),
            "oxygen_max": float(np.max(arr)),
        }
    # Writer stores (nx, ny, nz) — profile along depth axis = last.
    depth_mean = arr.mean(axis=(0, 1))
    return {
        "oxygen_mean": float(np.mean(arr)),
        "oxygen_min": float(np.min(arr)),
        "oxygen_max": float(np.max(arr)),
        "oxygen_z_mean_lo": float(depth_mean[0]),
        "oxygen_z_mean_hi": float(depth_mean[-1]),
    }


def _flatten_summary_scalars(summary: dict[str, Any]) -> dict[str, float | int]:
    out: dict[str, float | int] = {
        "step": int(summary["step"]),
        "time": float(summary["time"]),
        "n_total": int(summary.get("n_total", summary.get("num_agents", 0))),
    }
    if "num_lineages" in summary:
        out["num_lineages"] = int(summary["num_lineages"])
    if "dt" in summary:
        out["dt"] = float(summary["dt"])

    events = summary.get("events") or {}
    for key in EVENT_KEYS:
        if key in events:
            out[f"event_{key}"] = int(events[key])

    chem = summary.get("chem") or {}
    for key, value in chem.items():
        out[f"chem_{key}"] = float(value)

    spatial = summary.get("spatial") or {}
    for key, value in spatial.items():
        out[f"summary_spatial_{key}"] = float(value)

    for array_key in ("n_by_type", "n_in_crypt", "n_by_state"):
        if array_key not in summary:
            continue
        arr = np.asarray(summary[array_key]).ravel()
        for i, value in enumerate(arr):
            out[f"{array_key}_{i}"] = int(value) if array_key.startswith("n_") else float(value)

    for array_key in ("mean_z_by_type", "mean_mu_by_type"):
        if array_key not in summary:
            continue
        arr = np.asarray(summary[array_key]).ravel()
        for i, value in enumerate(arr):
            out[f"{array_key}_{i}"] = float(value)

    return out


def _agent_layer_stats(agents: dict[str, np.ndarray]) -> dict[str, float]:
    if "x" not in agents or "y" not in agents or "z" not in agents:
        return {"n_alive": 0.0}
    mask = _alive_mask(agents)
    out: dict[str, float] = {
        "n_agents": float(len(agents["x"])),
        "n_alive": float(np.sum(mask)),
        "n_dead": float(np.sum(~mask)),
    }
    if not np.any(mask):
        return out

    positions = _positions(agents, mask)
    types = agents["type"][mask] if "type" in agents else np.zeros(np.sum(mask), dtype=np.int32)
    out.update(_per_type_com_stats(positions, types))
    out.update(_spatial_metrics(positions, types))
    out.update(_growth_stats(agents, mask))
    out.update(_crypt_stats(agents, mask))

    if "state" in agents:
        states = agents["state"][mask]
        for state_id, count in zip(*np.unique(states, return_counts=True), strict=True):
            out[f"n_state_{int(state_id)}"] = float(count)
    return out


def extract_checkpoint_row(filepath: str | Path) -> dict[str, Any]:
    """Open one checkpoint, extract a flat summary row, then close the file."""
    path = validate_input_path(filepath)
    row: dict[str, Any] = {
        "file": path.name,
        "path": str(path),
    }
    with GutIBMData(path) as data:
        steps = data.steps
        if not steps:
            raise ValueError(f"no step groups found in {path}")
        step = steps[-1]
        row["step_group"] = step
        nx, ny, nz = data.grid_shape
        row["nx"] = nx
        row["ny"] = ny
        row["nz"] = nz
        row["grid_dx"] = data.grid_dx

        if data.has_layer("summary"):
            summary = data.get_summary(step)
            row.update(_flatten_summary_scalars(summary))

        try:
            agents = data.get_agents(step)
        except (KeyError, AssertionError, OSError):
            agents = {}
        if agents:
            row.update(_agent_layer_stats(agents))

        volumes = data.get_grid_volumes(step)
        if volumes:
            row["grid_species"] = ",".join(sorted(volumes))
            row.update(_grid_chem_means(volumes))
            row.update(_oxygen_depth_profile(volumes))
            # Record which expected chem names were present for autodiscovery metadata.
            present = [name for name in GRID_CHEM_NAMES if name in volumes]
            row["known_grid_species"] = ",".join(present)

    return row


def extract_many(
    paths: list[Path],
    *,
    progress: bool = True,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Extract rows for many checkpoints; return (rows, skip_messages)."""
    rows: list[dict[str, Any]] = []
    skips: list[str] = []
    total = len(paths)
    for index, path in enumerate(paths, start=1):
        if progress:
            print(f"checkpoint {index}/{total}: {path.name}", file=sys.stderr, flush=True)
        try:
            rows.append(extract_checkpoint_row(path))
        except (OSError, ValueError, KeyError, TypeError) as exc:
            message = f"{path.name}: {exc}"
            skips.append(message)
            print(f"warning: skipping {message}", file=sys.stderr, flush=True)
    return rows, skips
