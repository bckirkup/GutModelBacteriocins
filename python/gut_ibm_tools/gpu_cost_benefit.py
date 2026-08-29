"""GPU cost/benefit benchmark harness.

This module only assembles configurations, executes one requested arm, and
merges already-produced result files.  The merge path never starts GutIBM.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import subprocess
import time
from pathlib import Path
from statistics import median
from typing import Any

import h5py

from .batch_config import apply_overrides
from .path_utils import (
    PathValidationError,
    validate_input_path,
    validate_input_path_within,
    validate_path_syntax,
    write_json_file,
)

ARM_MATRIX: dict[str, dict[str, Any]] = {
    "A1": {
        "axis": "A",
        "metabolism.uptake_limit": "none",
        "gpu_enabled": False,
    },
    "A2": {
        "axis": "A",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
    },
    "A3": {
        "axis": "A",
        "metabolism.uptake_limit": "delivery",
        "gpu_enabled": False,
    },
    "A4": {
        "axis": "A",
        "metabolism.uptake_limit": "none",
        "gpu_enabled": True,
    },
    "A5": {
        "axis": "A",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": True,
    },
    "A6": {
        "axis": "A",
        "metabolism.uptake_limit": "delivery",
        "gpu_enabled": True,
        "status": "blocked",
        "blocked_reason": (
            'parser refusal: gpu_enabled=true cannot be combined with '
            'metabolism.uptake_limit="delivery": CUDA parity is not '
            "implemented yet"
        ),
    },
    # Keep FMM off for the non-FMM B baseline so B4 measures its declared
    # corrected-series-plus-FMM variant rather than duplicating B2.
    "B1": {
        "axis": "B",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "pre_fix_duplicated_reflection",
        "image_series_relative_tolerance": 0.0,
        "image_series_max_shells": 3,
        "use_fmm": False,
    },
    "B2": {
        "axis": "B",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        "use_fmm": False,
    },
    "B3": {
        "axis": "B",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-6,
        "use_fmm": False,
    },
    "B4": {
        "axis": "B",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        "use_fmm": True,
    },
    "C1": {
        "axis": "C",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        "toxin.lumen_transfer_length": "inf",
        "toxin.lumen_transfer_basis": "effective",
    },
    "C2": {
        "axis": "C",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        "toxin.lumen_transfer_length": 100.0e-6,
        "toxin.lumen_transfer_basis": "effective",
    },
    "C3": {
        "axis": "C",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        "toxin.lumen_transfer_length": 100.0e-6,
        "toxin.lumen_transfer_basis": "free",
    },
    "C4": {
        "axis": "C",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        "toxin.lumen_transfer_length": 300.0e-6,
        "toxin.lumen_transfer_basis": "effective",
    },
    "C5": {
        "axis": "C",
        "metabolism.uptake_limit": "sherwood",
        "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        "toxin.lumen_transfer_length": 30.0e-6,
        "toxin.lumen_transfer_basis": "effective",
    },
}

SCALES = ("1e5", "1e6")
SEEDS = (55, 56, 57)


def _without_comments(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _without_comments(item)
            for key, item in value.items()
            if not key.startswith("_")
        }
    if isinstance(value, list):
        return [_without_comments(item) for item in value]
    return value


def _read_json(path: Path) -> dict[str, Any]:
    with validate_input_path(path).open(encoding="utf-8") as handle:
        return _without_comments(json.load(handle))


def _read_declared_json(root: Path, declared: str) -> dict[str, Any]:
    with validate_input_path_within(root, declared).open(
        encoding="utf-8"
    ) as handle:
        return _without_comments(json.load(handle))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    write_json_file(path, payload)


def _scale_config(base: dict[str, Any], scale: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in scale.items():
        if key in {"hdf5_file", "hdf5"}:
            continue
        result[key] = copy.deepcopy(value)
    result["profile_steps"] = True
    result["hdf5_file"] = "benchmark-output.h5"
    result["hdf5"] = {
        "enabled": True,
        "schedule": {"summary": 1, "provenance": 1},
        "grid_species": [],
    }
    return result


def generate_configs(
    base_path: Path,
    output_dir: Path,
    scale_paths: dict[str, Path],
) -> dict[str, Any]:
    """Generate all arm/scale configs and a manifest from one base config."""
    base = _read_json(base_path)
    arms: dict[str, Any] = {}
    for scale, scale_path in scale_paths.items():
        scale_base = _scale_config(base, _read_json(scale_path))
        for arm, definition in ARM_MATRIX.items():
            config = copy.deepcopy(scale_base)
            overrides = {
                key: value
                for key, value in definition.items()
                if key not in {"axis", "status", "blocked_reason"}
            }
            apply_overrides(config, overrides)
            arm_dir = output_dir / arm
            config_path = arm_dir / f"{scale}.json"
            _write_json(config_path, config)
            arms.setdefault(arm, {
                "axis": definition["axis"],
                "overrides": overrides,
                "status": definition.get("status", "runnable"),
                "blocked_reason": definition.get("blocked_reason"),
            })
            arms[arm].setdefault("configs", {})[scale] = str(config_path.resolve())
    manifest = {
        "schema_version": 1,
        "scales": list(scale_paths),
        "seeds": list(SEEDS),
        "arms": arms,
        "baselines": {"A": "A3", "B": "B2", "C": "C1"},
        "reporting": {
            "missing_records": "explicit",
            "raw_result_files_retained": True,
            "simulation_execution": "runner_only",
        },
    }
    _write_json(output_dir / "manifest.json", manifest)
    return manifest


def _hdf5_value(value: Any) -> Any:
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def _read_hdf5_group(group: h5py.Group) -> dict[str, Any]:
    values: dict[str, Any] = {}
    for name, item in group.items():
        if isinstance(item, h5py.Group):
            values[name] = _read_hdf5_group(item)
        else:
            values[name] = _hdf5_value(item[()])
    return values


def _read_provenance(path: Path) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        if "run_provenance" not in handle:
            return {}
        return _read_hdf5_group(handle["run_provenance"])


def _validated_binary_path(binary: Path) -> Path:
    candidate = validate_path_syntax(binary).resolve()
    if candidate.name != "gut_ibm":
        raise PathValidationError("benchmark binary must be named gut_ibm")
    if not candidate.is_relative_to(Path.cwd().resolve()):
        raise PathValidationError(
            f"benchmark binary must be under the working directory: {binary}"
        )
    return candidate


def _validated_mpi_launcher(mpirun: str | None) -> str | None:
    if mpirun is None:
        return None
    launcher = validate_path_syntax(mpirun).name
    if launcher == "mpirun":
        return "mpirun"
    if launcher == "mpiexec":
        return "mpiexec"
    if launcher == "mpiexec.hydra":
        return "mpiexec.hydra"
    raise PathValidationError(f"unsupported MPI launcher: {mpirun}")


def _validated_mpi_ranks(mpi_ranks: int) -> str:
    if not 1 <= mpi_ranks <= 4096:
        raise ValueError("mpi_ranks must be between 1 and 4096")
    return str(mpi_ranks)


def _run_pass(
    pass_name: str,
    profile_steps: bool,
    config: dict[str, Any],
    config_path: Path,
    hdf5_path: Path,
    binary: Path,
    *,
    mpirun: str | None,
    mpi_ranks: int,
) -> dict[str, Any]:
    config["profile_steps"] = profile_steps
    config["hdf5_file"] = str(hdf5_path.resolve())
    _write_json(config_path, config)
    config_argument = str(config_path)
    if mpirun is None:
        command = ["./gut_ibm", config_argument]
    else:
        command = [
            mpirun,
            "-np",
            _validated_mpi_ranks(mpi_ranks),
            "./gut_ibm",
            config_argument,
        ]
    started = time.monotonic()
    completed = subprocess.run(command, cwd=binary.parent, check=False)
    wall_seconds = time.monotonic() - started
    status = "completed" if completed.returncode == 0 else "failed"
    return {
        "name": pass_name,
        "profile_steps": profile_steps,
        "config_file": str(config_path.resolve()),
        "hdf5_file": str(hdf5_path.resolve()),
        "status": status,
        "completion_status": status,
        "exit_code": completed.returncode,
        "wall_seconds": wall_seconds,
        "provenance": (
            _read_provenance(hdf5_path) if hdf5_path.is_file() else {}
        ),
    }


def _pass_pair_overhead(passes: dict[str, Any]) -> tuple[float | None, float | None]:
    cost = passes.get("cost", {})
    profile = passes.get("profile", {})
    cost_wall = cost.get("wall_seconds")
    profile_wall = profile.get("wall_seconds")
    if not all(isinstance(value, (int, float))
               for value in (cost_wall, profile_wall)):
        return None, None
    difference = float(profile_wall) - float(cost_wall)
    if float(cost_wall) <= 0.0:
        return difference, None
    return difference, difference / float(cost_wall)


def run_one_arm(
    manifest_path: Path,
    arm: str,
    scale: str,
    seed: int,
    binary: Path,
    output_dir: Path,
    *,
    mpirun: str | None = None,
    mpi_ranks: int = 1,
) -> Path:
    """Run one arm/scale/seed with separate cost and profile passes."""
    manifest = _read_json(manifest_path)
    arm_info = manifest["arms"][arm]
    result_path = output_dir / f"{arm}_{scale}_seed{seed}.json"
    config = _read_declared_json(
        manifest_path.parent, arm_info["configs"][scale]
    )
    config["seed"] = seed
    record: dict[str, Any] = {
        "schema_version": 2,
        "arm": arm,
        "axis": arm_info["axis"],
        "scale": scale,
        "seed": seed,
        "host": {"hostname": os.uname().nodename},
        "backend": {"gpu_requested": bool(config.get("gpu_enabled", False))},
    }
    if arm_info["status"] == "blocked":
        blocked_passes = {
            name: {
                "name": name,
                "profile_steps": profile_steps,
                "status": "blocked",
                "completion_status": "not_run",
                "exit_code": None,
                "wall_seconds": None,
                "provenance": {},
            }
            for name, profile_steps in (("cost", False), ("profile", True))
        }
        record.update({
            "status": "blocked",
            "completion_status": "not_run",
            "blocked_reason": arm_info["blocked_reason"],
            "passes": blocked_passes,
            "profiling_wall_difference_seconds": None,
            "profiling_wall_difference_fraction": None,
        })
        _write_json(result_path, record)
        return result_path

    binary = _validated_binary_path(binary)
    mpirun = _validated_mpi_launcher(mpirun)
    passes: dict[str, Any] = {}
    for pass_name, profile_steps in (("cost", False), ("profile", True)):
        pass_config = copy.deepcopy(config)
        pass_config_path = output_dir / (
            f"{arm}_{scale}_seed{seed}.{pass_name}.config.json"
        )
        pass_hdf5_path = output_dir / (
            f"{arm}_{scale}_seed{seed}.{pass_name}.h5"
        )
        passes[pass_name] = _run_pass(
            pass_name,
            profile_steps,
            pass_config,
            pass_config_path,
            pass_hdf5_path,
            binary,
            mpirun=mpirun,
            mpi_ranks=mpi_ranks,
        )
    pass_statuses = [item["status"] for item in passes.values()]
    status = "completed" if all(
        item == "completed" for item in pass_statuses
    ) else "failed"
    difference, fraction = _pass_pair_overhead(passes)
    record.update({
        "status": status,
        "completion_status": status,
        "passes": passes,
        "profiling_wall_difference_seconds": difference,
        "profiling_wall_difference_fraction": fraction,
    })
    _write_json(result_path, record)
    return result_path


def _completed_pass_pair(record: dict[str, Any]) -> bool:
    passes = record.get("passes", {})
    return (
        record.get("status") == "completed"
        and all(
            passes.get(name, {}).get("status") == "completed"
            for name in ("cost", "profile")
        )
    )


def _numeric_pass_values(
    records: list[dict[str, Any]], key: str
) -> list[float]:
    values = []
    for record in records:
        value = record.get(key)
        if isinstance(value, (int, float)):
            values.append(float(value))
    return values


def _median_numeric_value(
    records: list[dict[str, Any]], key: str
) -> float | None:
    values = _numeric_pass_values(records, key)
    return median(values) if values else None


def merge_results(manifest_path: Path, result_paths: list[Path]) -> dict[str, Any]:
    """Merge raw records, preserving explicit missing arm/scale/seed entries."""
    manifest = _read_json(manifest_path)
    records: dict[tuple[str, str, int], dict[str, Any]] = {}
    for path in result_paths:
        record = _read_json(path)
        key = (record["arm"], record["scale"], int(record["seed"]))
        if key in records:
            raise ValueError(f"duplicate benchmark result: {key}")
        record["raw_result_file"] = str(Path(path).resolve())
        records[key] = record
    merged: dict[str, Any] = {
        "schema_version": 2,
        "scales": manifest["scales"],
        "seeds": manifest["seeds"],
        "baselines": manifest["baselines"],
        "reporting": manifest.get("reporting", {}),
        "arms": {},
    }
    for arm, arm_info in manifest["arms"].items():
        arm_result = {
            "axis": arm_info["axis"],
            "status": arm_info["status"],
            "blocked_reason": arm_info.get("blocked_reason"),
            "scales": {},
        }
        for scale in manifest["scales"]:
            entries = []
            for seed in manifest["seeds"]:
                record = records.get((arm, scale, int(seed)))
                if record is not None:
                    entries.append(record)
                elif arm_info["status"] == "blocked":
                    entries.append({
                        "status": "blocked",
                        "arm": arm,
                        "scale": scale,
                        "seed": seed,
                        "blocked_reason": arm_info["blocked_reason"],
                    })
                else:
                    entries.append({
                        "status": "missing",
                        "arm": arm,
                        "scale": scale,
                        "seed": seed,
                    })
            completed = [item for item in entries if _completed_pass_pair(item)]
            cost_passes = [item["passes"]["cost"] for item in completed]
            profile_passes = [item["passes"]["profile"] for item in completed]
            phase_totals: dict[str, list[float]] = {}
            for item in profile_passes:
                profile = item.get("provenance", {}).get("step_profile", {})
                for phase, value in profile.items():
                    if phase == "step_count" or not isinstance(value, (int, float)):
                        continue
                    phase_totals.setdefault(phase, []).append(float(value))
            cost_walls = [
                item["wall_seconds"] for item in cost_passes
                if isinstance(item.get("wall_seconds"), (int, float))
            ]
            steps_per_second = []
            for item in cost_passes:
                steps = item.get("provenance", {}).get("termination_step")
                wall = item.get("wall_seconds")
                if (isinstance(steps, (int, float))
                        and isinstance(wall, (int, float))
                        and wall > 0.0):
                    steps_per_second.append(float(steps) / float(wall))
            kernel_counts = []
            direct_counts = []
            table_builds = []
            table_evictions = []
            for item in profile_passes:
                provenance = item.get("provenance", {})
                for target, key in (
                    (kernel_counts, "green_function_kernel_evaluations"),
                    (direct_counts, "robin_direct_mode_evaluations"),
                    (table_builds, "robin_tables_built"),
                    (table_evictions, "robin_table_evictions"),
                ):
                    value = provenance.get(key)
                    if isinstance(value, (int, float)):
                        target.append(float(value))
            arm_result["scales"][scale] = {
                "records": entries,
                "available_count": len(completed),
                "missing_count": sum(item.get("status") == "missing"
                                     for item in entries),
                "cost_summary": {
                    "wall_seconds_median": (
                        median(cost_walls) if cost_walls else None
                    ),
                    "steps_per_second_median": (
                        median(steps_per_second) if steps_per_second else None
                    ),
                },
                "attribution_summary": {
                    "step_profile_median": {
                        phase: median(values)
                        for phase, values in phase_totals.items()
                    },
                    "kernel_evaluations_median": (
                        median(kernel_counts) if kernel_counts else None
                    ),
                    "robin_direct_mode_evaluations_median": (
                        median(direct_counts) if direct_counts else None
                    ),
                    "robin_tables_built_median": (
                        median(table_builds) if table_builds else None
                    ),
                    "robin_table_evictions_median": (
                        median(table_evictions) if table_evictions else None
                    ),
                },
                "profiling_overhead": {
                    "wall_seconds_median": _median_numeric_value(
                        completed, "profiling_wall_difference_seconds"
                    ),
                    "fraction_median": _median_numeric_value(
                        completed, "profiling_wall_difference_fraction"
                    ),
                },
            }
        merged["arms"][arm] = arm_result
    return merged


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--base", type=Path, required=True)
    generate.add_argument("--output-dir", type=Path, required=True)
    generate.add_argument("--scale", action="append", nargs=2, metavar=("NAME", "PATH"),
                          required=True)
    run = subparsers.add_parser("run")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--arm", required=True)
    run.add_argument("--scale", required=True)
    run.add_argument("--seed", type=int, required=True)
    run.add_argument("--binary", type=Path, required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--mpirun")
    run.add_argument("--mpi-ranks", type=int, default=1)
    merge = subparsers.add_parser("merge")
    merge.add_argument("--manifest", type=Path, required=True)
    merge.add_argument("--output", type=Path, required=True)
    merge.add_argument("results", nargs="+", type=Path)
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    if args.command == "generate":
        generate_configs(args.base, args.output_dir, dict(args.scale))
    elif args.command == "run":
        run_one_arm(args.manifest, args.arm, args.scale, args.seed,
                    args.binary, args.output_dir, mpirun=args.mpirun,
                    mpi_ranks=args.mpi_ranks)
    else:
        merged = merge_results(args.manifest, args.results)
        _write_json(args.output, merged)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
