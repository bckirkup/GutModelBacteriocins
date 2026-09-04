"""GPU cost/benefit benchmark harness.

The harness assembles configurations, executes one requested arm, and merges
already-produced result files.  The merge and report paths never start GutIBM
or read simulation HDF5 files.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import platform
import socket
import subprocess
import time
from pathlib import Path
from statistics import median
from typing import Any

import h5py

from .path_utils import (
    PathValidationError,
    prepare_output_file,
    validate_input_path,
    validate_input_path_within,
    validate_path_syntax,
    write_json_file,
)

UPTAKE_LIMIT_KEY = "metabolism.uptake_limit"
LUMEN_TRANSFER_LENGTH_KEY = "toxin.lumen_transfer_length"
LUMEN_TRANSFER_BASIS_KEY = "toxin.lumen_transfer_basis"
DRIFT_CORRECTION_KEY = "qssa.drift_correction"
MPI_LAUNCHER = "/usr/bin/mpirun"

_A_SCALES = ("s1", "s2")
_CPU_SCALES = ("s0",)
_E_SCALES = ("p1",)
_E_SEEDS = (55, 56, 57, 58, 59)

ARM_MATRIX: dict[str, dict[str, Any]] = {
    "A1": {"axis": "A", UPTAKE_LIMIT_KEY: "none", "gpu_enabled": False,
           "accepted_placements": {"host"}, "scales": _A_SCALES},
    "A2": {"axis": "A", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
           "accepted_placements": {"host"}, "scales": _A_SCALES},
    # host_forced_delivery requires a declined device request, unavailable
    # when GPU delivery is disabled.
    "A3": {"axis": "A", UPTAKE_LIMIT_KEY: "delivery", "gpu_enabled": False,
           "accepted_placements": {"host"}, "scales": _A_SCALES},
    "A4": {"axis": "A", UPTAKE_LIMIT_KEY: "none", "gpu_enabled": True,
           "accepted_placements": {"device"}, "scales": _A_SCALES},
    "A5": {"axis": "A", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": True,
           "accepted_placements": {"device"}, "scales": _A_SCALES},
    "A6": {"axis": "A", UPTAKE_LIMIT_KEY: "delivery", "gpu_enabled": True,
           "accepted_placements": {"device_delivery"}, "scales": _A_SCALES},
    "E1": {
        "axis": "E", UPTAKE_LIMIT_KEY: "none", "gpu_enabled": False,
        "accepted_placements": {"host"}, "scales": _E_SCALES,
        "outcome_only": True,
    },
    "E2": {
        "axis": "E", UPTAKE_LIMIT_KEY: "none", "gpu_enabled": True,
        "accepted_placements": {"device"}, "scales": _E_SCALES,
        "outcome_only": True,
    },
    "E3": {
        "axis": "E", UPTAKE_LIMIT_KEY: "delivery", "gpu_enabled": False,
        "accepted_placements": {"host"}, "scales": _E_SCALES,
        "outcome_only": True,
    },
    "E4": {
        "axis": "E", UPTAKE_LIMIT_KEY: "delivery", "gpu_enabled": True,
        "accepted_placements": {"device_delivery"}, "scales": _E_SCALES,
        "outcome_only": True,
    },
    "E2r": {
        "axis": "E", UPTAKE_LIMIT_KEY: "none", "gpu_enabled": True,
        "accepted_placements": {"device"}, "scales": _E_SCALES,
        "outcome_only": True,
    },
    "B1": {
        "axis": "B", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "pre_fix_duplicated_reflection",
        "image_series_relative_tolerance": 0.0,
        "image_series_max_shells": 3, "use_fmm": False,
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "B2": {
        "axis": "B", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10, "use_fmm": False,
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "B3": {
        "axis": "B", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-6, "use_fmm": False,
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "B4": {
        "axis": "B", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10, "use_fmm": True,
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "C1": {
        "axis": "C", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        LUMEN_TRANSFER_LENGTH_KEY: "inf",
        LUMEN_TRANSFER_BASIS_KEY: "effective",
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "C2": {
        "axis": "C", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        LUMEN_TRANSFER_LENGTH_KEY: 100.0e-6,
        LUMEN_TRANSFER_BASIS_KEY: "effective",
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "C3": {
        "axis": "C", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        LUMEN_TRANSFER_LENGTH_KEY: 100.0e-6,
        LUMEN_TRANSFER_BASIS_KEY: "free",
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "C4": {
        "axis": "C", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        LUMEN_TRANSFER_LENGTH_KEY: 300.0e-6,
        LUMEN_TRANSFER_BASIS_KEY: "effective",
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "C5": {
        "axis": "C", UPTAKE_LIMIT_KEY: "sherwood", "gpu_enabled": False,
        "image_series_mode": "corrected",
        "image_series_relative_tolerance": 1.0e-10,
        LUMEN_TRANSFER_LENGTH_KEY: 30.0e-6,
        LUMEN_TRANSFER_BASIS_KEY: "effective",
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
    "D1": {"axis": "D", DRIFT_CORRECTION_KEY: False, "gpu_enabled": False,
           "accepted_placements": {"host"}, "scales": _CPU_SCALES},
    "D2": {"axis": "D", DRIFT_CORRECTION_KEY: True, "gpu_enabled": False,
           "accepted_placements": {"host"}, "scales": _CPU_SCALES},
    "D3": {
        "axis": "D", DRIFT_CORRECTION_KEY: True, "gpu_enabled": True,
        "status": "blocked",
        "blocked_reason": (
            "corrected sealed kernel routes to host fallback by design (#394)"
        ),
        "accepted_placements": {"host"}, "scales": _CPU_SCALES,
    },
}

SCALES = ("s0", "s1", "s2", "p1")
SEEDS = (55, 56, 57)


def _without_comments(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _without_comments(item)
            for key, item in value.items()
            if not str(key).startswith("_")
        }
    if isinstance(value, list):
        return [_without_comments(item) for item in value]
    return value


def _read_json(path: Path) -> dict[str, Any]:
    with validate_input_path(path).open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"JSON object required: {path}")
    return _without_comments(payload)


def _read_declared_json(root: Path, declared: str) -> dict[str, Any]:
    with validate_input_path_within(root, declared).open(
        encoding="utf-8"
    ) as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"JSON object required: {declared}")
    return _without_comments(payload)


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if hasattr(value, "tolist"):
        return _json_safe(value.tolist())
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    write_json_file(path, _json_safe(payload), indent=2, allow_external=True)


def _scale_config(
    base: dict[str, Any], scale: dict[str, Any], scale_name: str | None = None
) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in scale.items():
        if key not in {"hdf5_file", "hdf5"}:
            result[key] = copy.deepcopy(value)
    if scale_name == "p1":
        result["total_time"] = 21600
        result["profile_steps"] = False
        result["hdf5"] = {
            "enabled": True,
            "schedule": {
                "summary": 1, "provenance": 1, "agents": 0, "grid": 0,
                "lineage": 0, "genome": 0,
            },
        }
        result["hdf5_file"] = "gpu-precision-p1.h5"
    else:
        result["profile_steps"] = True
        result["hdf5_file"] = (
            f"benchmark-{scale_name}.h5" if scale_name else "benchmark-output.h5"
        )
        result["hdf5"] = {
            "enabled": True,
            "schedule": {"summary": 1, "provenance": 1, "grid_species": []},
        }
    return result


def _arm_overrides(definition: dict[str, Any]) -> dict[str, Any]:
    excluded = {
        "axis", "status", "blocked_reason", "accepted_placements", "scales",
        "outcome_only",
    }
    return {
        key: copy.deepcopy(value)
        for key, value in definition.items()
        if key not in excluded
    }


def generate_configs(
    base_path: Path, output_dir: Path, scale_paths: dict[str, Path]
) -> dict[str, Any]:
    """Generate matrix configurations and a manifest."""
    base = _read_json(base_path)
    arms: dict[str, Any] = {}
    for scale, scale_path in scale_paths.items():
        if scale not in SCALES:
            raise ValueError(f"unknown benchmark scale: {scale}")
        scale_base = _scale_config(base, _read_json(scale_path), scale)
        for arm, definition in ARM_MATRIX.items():
            if scale not in definition["scales"]:
                continue
            config = copy.deepcopy(scale_base)
            overrides = _arm_overrides(definition)
            # Strict JSON parsing accepts dotted root keys, not nested fix
            # objects.  Keep these keys literal in generated JSON.
            config.update(overrides)
            config["hdf5_file"] = f"{arm}_{scale}_seed{{seed}}.h5"
            arm_dir = output_dir / arm
            config_path = arm_dir / f"{scale}.json"
            _write_json(config_path, config)
            arm_entry = arms.setdefault(
                arm,
                {
                    "axis": definition["axis"],
                    "overrides": overrides,
                    "status": definition.get("status", "runnable"),
                    "blocked_reason": definition.get("blocked_reason"),
                    "outcome_only": definition.get("outcome_only", False),
                    "accepted_placements": sorted(
                        definition["accepted_placements"]
                    ),
                    "scales": list(definition["scales"]),
                    "configs": {},
                },
            )
            arm_entry["configs"][scale] = str(config_path.resolve())
    manifest = {
        "schema_version": 2,
        "scales": list(scale_paths),
        "seeds": list(_E_SEEDS if "p1" in scale_paths else SEEDS),
        "arms": arms,
        "baselines": {"A": "A3", "B": "B2", "C": "C1", "D": "D1"},
        "reporting": {
            "missing_records": "explicit",
            "invalid_placement_records": "explicit",
            "raw_result_files_retained": True,
            "simulation_execution": "runner_only",
            "primary_gpu_metric": "chemistry_s",
        },
    }
    _write_json(output_dir / "manifest.json", manifest)
    return manifest


def _hdf5_value(value: Any) -> Any:
    if hasattr(value, "shape") and getattr(value, "shape", ()) == ():
        value = value.item()
    elif hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, list) and len(value) == 1:
        value = value[0]
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
    if not candidate.is_file():
        raise PathValidationError(f"benchmark binary is not a regular file: {binary}")
    if not os.access(candidate, os.X_OK):
        raise PathValidationError(f"benchmark binary is not executable: {binary}")
    return candidate


def _validated_mpi_launcher(mpirun: str | None) -> str | None:
    if mpirun is None:
        return None
    launcher = validate_path_syntax(mpirun).name
    if launcher in {"mpirun", "mpiexec", "mpiexec.hydra"}:
        return str(validate_path_syntax(mpirun))
    raise PathValidationError(f"unsupported MPI launcher: {mpirun}")


def _validated_mpi_ranks(mpi_ranks: int) -> str:
    if not 1 <= mpi_ranks <= 4096:
        raise ValueError("mpi_ranks must be between 1 and 4096")
    return str(mpi_ranks)


def _read_cpu_model() -> str | None:
    try:
        with Path("/proc/cpuinfo").open(encoding="utf-8") as handle:
            for line in handle:
                if line.lower().startswith("model name"):
                    _, _, value = line.partition(":")
                    return value.strip() or None
    except (OSError, UnicodeError):
        return platform.processor() or None
    return platform.processor() or None


def _read_mem_total_kib() -> int | None:
    try:
        with Path("/proc/meminfo").open(encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    fields = line.split()
                    return int(fields[1]) if len(fields) >= 2 else None
    except (OSError, ValueError, IndexError):
        return None
    return None


def _probe_gpu() -> dict[str, Any]:
    result: dict[str, Any] = {
        "name": None, "driver_version": None, "memory_total": None
    }
    try:
        completed = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=name,driver_version,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True, text=True, check=False, timeout=5
        )
    except (OSError, subprocess.SubprocessError):
        return result
    if completed.returncode != 0 or not completed.stdout.strip():
        return result
    fields = [field.strip() for field in completed.stdout.splitlines()[0].split(",")]
    if fields:
        result["name"] = fields[0] or None
    if len(fields) > 1:
        result["driver_version"] = fields[1] or None
    if len(fields) > 2:
        result["memory_total"] = fields[2] or None
    return result


def _host_identity() -> dict[str, Any]:
    gpu = _probe_gpu()
    cpu_model = _read_cpu_model()
    logical_cores = os.cpu_count()
    mem_total_kib = _read_mem_total_kib()
    aws_batch = {
        "job_id": os.environ.get("AWS_BATCH_JOB_ID"),
        "jq_name": os.environ.get("AWS_BATCH_JQ_NAME"),
        "ce_name": os.environ.get("AWS_BATCH_CE_NAME"),
    }
    instance_type = (
        os.environ.get("EC2_INSTANCE_TYPE")
        or os.environ.get("AWS_INSTANCE_TYPE")
        or os.environ.get("AWS_BATCH_INSTANCE_TYPE")
    )
    identity_text = json.dumps(
        {
            "cpu_model": cpu_model,
            "logical_cores": logical_cores,
            "mem_total_kib": mem_total_kib,
            "gpu_name": gpu.get("name") or "absent",
        },
        sort_keys=True, separators=(",", ":")
    )
    return {
        "hostname": socket.gethostname(),
        "cpu_model": cpu_model,
        "logical_cores": logical_cores,
        "mem_total_kib": mem_total_kib,
        "gpu": gpu,
        "aws_batch": aws_batch,
        "instance_type": instance_type,
        "host_fingerprint": hashlib.sha256(
            identity_text.encode("utf-8")
        ).hexdigest(),
    }


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
    config_path.parent.mkdir(parents=True, exist_ok=True)
    hdf5_path.parent.mkdir(parents=True, exist_ok=True)
    _write_json(config_path, config)
    config_argument = str(config_path.resolve())
    launcher = mpirun is not None
    if launcher:
        _validated_mpi_launcher(mpirun)
    command = [str(binary), config_argument]
    if launcher:
        command = [
            MPI_LAUNCHER, "-np", _validated_mpi_ranks(mpi_ranks),
            str(binary), config_argument
        ]
    stdout_path = hdf5_path.with_suffix(".stdout")
    stderr_path = hdf5_path.with_suffix(".stderr")
    started = time.monotonic()
    completed = subprocess.run(
        command, cwd=binary.parent, check=False,
        capture_output=True, text=True
    )
    wall_seconds = time.monotonic() - started
    stdout_path.write_text(completed.stdout or "", encoding="utf-8")
    stderr_path.write_text(completed.stderr or "", encoding="utf-8")
    status = "completed" if completed.returncode == 0 else "failed"
    provenance = _read_provenance(hdf5_path) if hdf5_path.is_file() else {}
    return {
        "name": pass_name, "profile_steps": profile_steps,
        "config_file": str(config_path.resolve()),
        "hdf5_file": str(hdf5_path.resolve()),
        "stdout_file": str(stdout_path.resolve()),
        "stderr_file": str(stderr_path.resolve()),
        "status": status, "completion_status": status,
        "exit_code": completed.returncode, "wall_seconds": wall_seconds,
        "provenance": provenance,
    }


def _pass_pair_overhead(passes: dict[str, Any]) -> tuple[float | None, float | None]:
    cost_wall = passes.get("cost", {}).get("wall_seconds")
    profile_wall = passes.get("profile", {}).get("wall_seconds")
    if not all(isinstance(value, (int, float)) for value in (cost_wall, profile_wall)):
        return None, None
    difference = float(profile_wall) - float(cost_wall)
    return difference, (
        difference / float(cost_wall) if float(cost_wall) > 0.0 else None
    )


def _pass_plan(arm_info: dict[str, Any]) -> tuple[tuple[str, bool], ...]:
    if arm_info.get("outcome_only"):
        return (("outcome", False),)
    return (("cost", False), ("profile", True))


def _placement_validation(
    pass_record: dict[str, Any],
    arm_info: dict[str, Any],
    gpu_requested: bool,
) -> dict[str, Any]:
    provenance = pass_record.get("provenance", {})
    observed = provenance.get("chemistry_placement")
    gpu_compiled = provenance.get("gpu_compiled")
    accepted = sorted(arm_info.get("accepted_placements", []))
    failures: list[str] = []
    if observed not in accepted:
        failures.append(f"chemistry placement {observed!r} is not in {accepted!r}")
    if gpu_requested and gpu_compiled != 1:
        failures.append(f"gpu_compiled must be 1 for GPU arm, got {gpu_compiled!r}")
    validation = {
        "observed_placement": observed,
        "gpu_compiled": gpu_compiled,
        "accepted_placements": accepted,
        "valid": not failures,
    }
    if failures:
        validation["reason"] = "; ".join(failures)
        pass_record["placement_status"] = "invalid_placement"
    else:
        pass_record["placement_status"] = "accepted"
    pass_record["placement_validation"] = validation
    return validation


def _current_arm_definition(
    arm: str, manifest_definition: dict[str, Any]
) -> dict[str, Any]:
    return ARM_MATRIX.get(arm, manifest_definition)


def _gpu_requested_for_arm(arm_definition: dict[str, Any]) -> bool:
    if "gpu_enabled" in arm_definition:
        return bool(arm_definition["gpu_enabled"])
    return bool(arm_definition.get("overrides", {}).get("gpu_enabled", False))


def _revalidate_record_placement(
    record: dict[str, Any], arm_definition: dict[str, Any]
) -> None:
    if record.get("status") in {"blocked", "missing", "not_applicable"}:
        return
    passes = record.get("passes", {})
    if not isinstance(passes, dict):
        return
    gpu_requested = _gpu_requested_for_arm(arm_definition)
    stored_status = record.get("status")
    stored_pass_statuses = {
        name: pass_record.get("placement_status")
        for name, pass_record in passes.items()
        if isinstance(pass_record, dict)
    }
    invalid_passes: list[str] = []
    for pass_name, pass_record in passes.items():
        if not isinstance(pass_record, dict):
            continue
        if pass_record.get("status") != "completed":
            continue
        validation = _placement_validation(
            pass_record, arm_definition, gpu_requested
        )
        if not validation["valid"]:
            invalid_passes.append(pass_name)
    pass_statuses = [
        pass_record.get("status")
        for pass_record in passes.values()
        if isinstance(pass_record, dict)
    ]
    if invalid_passes:
        status = "invalid_placement"
        record["invalid_placement_passes"] = invalid_passes
        record["invalid_placement_reason"] = (
            "provenance gate failed for pass(es): " + ", ".join(invalid_passes)
        )
    elif all(item == "completed" for item in pass_statuses):
        status = "completed"
        record.pop("invalid_placement_passes", None)
        record.pop("invalid_placement_reason", None)
    else:
        status = "failed"
    rederived_pass_statuses = {
        name: pass_record.get("placement_status")
        for name, pass_record in passes.items()
        if isinstance(pass_record, dict)
    }
    if stored_status != status or stored_pass_statuses != rederived_pass_statuses:
        record["placement_revalidated"] = {
            "stored": stored_status,
            "rederived": status,
        }
    record["status"] = status
    record["completion_status"] = (
        "not_run" if status == "blocked" else status
    )


def run_one_arm(
    manifest_path: Path, arm: str, scale: str, seed: int, binary: Path,
    output_dir: Path, *, mpirun: str | None = None, mpi_ranks: int = 1
) -> Path:
    """Run one arm/scale/seed with timing or outcome-only passes."""
    manifest = _read_json(manifest_path)
    arm_info = manifest["arms"][arm]
    if scale not in arm_info.get("configs", {}):
        raise ValueError(f"scale {scale!r} is not configured for arm {arm}")
    result_path = output_dir / f"{arm}_{scale}_seed{seed}.json"
    prepare_output_file(result_path)
    config = _read_declared_json(
        manifest_path.parent, arm_info["configs"][scale]
    )
    config["seed"] = seed
    gpu_requested = bool(config.get("gpu_enabled", False))
    record: dict[str, Any] = {
        "schema_version": 3, "arm": arm, "axis": arm_info["axis"],
        "scale": scale, "seed": seed, "host": _host_identity(),
        "backend": {"gpu_requested": gpu_requested},
        "accepted_placements": sorted(arm_info.get("accepted_placements", [])),
    }
    if arm_info.get("status") == "blocked":
        blocked_passes = {
            name: {
                "name": name, "profile_steps": profile_steps,
                "status": "blocked", "completion_status": "not_run",
                "exit_code": None, "wall_seconds": None, "provenance": {},
            }
            for name, profile_steps in _pass_plan(arm_info)
        }
        record.update({
            "status": "blocked", "completion_status": "not_run",
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
    for pass_name, profile_steps in _pass_plan(arm_info):
        pass_config = copy.deepcopy(config)
        pass_config_path = output_dir / (
            f"{arm}_{scale}_seed{seed}.{pass_name}.config.json"
        )
        pass_hdf5_path = output_dir / f"{arm}_{scale}_seed{seed}.{pass_name}.h5"
        passes[pass_name] = _run_pass(
            pass_name, profile_steps, pass_config, pass_config_path,
            pass_hdf5_path, binary, mpirun=mpirun, mpi_ranks=mpi_ranks
        )
    invalid_passes = []
    for pass_name, pass_record in passes.items():
        if pass_record["status"] == "completed":
            validation = _placement_validation(
                pass_record, arm_info, gpu_requested
            )
            if not validation["valid"]:
                invalid_passes.append(pass_name)
    pass_statuses = [item["status"] for item in passes.values()]
    if invalid_passes:
        status = "invalid_placement"
        record["invalid_placement_passes"] = invalid_passes
        record["invalid_placement_reason"] = (
            "provenance gate failed for pass(es): " + ", ".join(invalid_passes)
        )
    elif all(item == "completed" for item in pass_statuses):
        status = "completed"
    else:
        status = "failed"
    difference, fraction = _pass_pair_overhead(passes)
    record.update({
        "status": status, "completion_status": status, "passes": passes,
        "profiling_wall_difference_seconds": difference,
        "profiling_wall_difference_fraction": fraction,
    })
    _write_json(result_path, record)
    return result_path


def _completed_pass_pair(record: dict[str, Any]) -> bool:
    passes = record.get("passes", {})
    return bool(passes) and record.get("status") == "completed" and all(
        isinstance(pass_record, dict)
        and pass_record.get("status") == "completed"
        for pass_record in passes.values()
    )


def _cost_pass(record: dict[str, Any]) -> dict[str, Any]:
    passes = record.get("passes", {})
    return passes.get("cost") or passes.get("outcome") or {}


def _profile_pass(record: dict[str, Any]) -> dict[str, Any]:
    passes = record.get("passes", {})
    return passes.get("profile") or passes.get("outcome") or {}


def _numeric_pass_values(records: list[dict[str, Any]], key: str) -> list[float]:
    return [
        float(record[key]) for record in records
        if isinstance(record.get(key), (int, float))
    ]


def _median_numeric_value(
    records: list[dict[str, Any]], key: str
) -> float | None:
    values = _numeric_pass_values(records, key)
    return median(values) if values else None


def _record_fingerprint(record: dict[str, Any]) -> str | None:
    fingerprint = record.get("host", {}).get("host_fingerprint")
    return fingerprint if isinstance(fingerprint, str) else None


_COUNTER_KEYS = (
    "green_function_kernel_evaluations", "robin_direct_mode_evaluations",
    "robin_tables_built", "robin_table_evictions"
)


def _phase_totals(profile_passes: list[dict[str, Any]]) -> dict[str, list[float]]:
    phase_totals: dict[str, list[float]] = {}
    for item in profile_passes:
        profile = item.get("provenance", {}).get("step_profile", {})
        for phase, value in profile.items():
            if phase != "step_count" and isinstance(value, (int, float)):
                phase_totals.setdefault(phase, []).append(float(value))
    return phase_totals


def _steps_per_second(cost_passes: list[dict[str, Any]]) -> list[float]:
    values = []
    for item in cost_passes:
        steps = item.get("provenance", {}).get("termination_step")
        wall = item.get("wall_seconds")
        if (
            isinstance(steps, (int, float))
            and isinstance(wall, (int, float))
            and wall > 0
        ):
            values.append(float(steps) / float(wall))
    return values


def _counter_values(
    profile_passes: list[dict[str, Any]]
) -> dict[str, list[float]]:
    values = {key: [] for key in _COUNTER_KEYS}
    for item in profile_passes:
        provenance = item.get("provenance", {})
        for key in _COUNTER_KEYS:
            value = provenance.get(key)
            if isinstance(value, (int, float)):
                values[key].append(float(value))
    return values


def _summary_for_entries(entries: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [item for item in entries if _completed_pass_pair(item)]
    cost_passes = [_cost_pass(item) for item in completed]
    profile_passes = [_profile_pass(item) for item in completed]
    cost_walls = [
        item["wall_seconds"] for item in cost_passes
        if isinstance(item.get("wall_seconds"), (int, float))
    ]
    phase_totals = _phase_totals(profile_passes)
    steps_per_second = _steps_per_second(cost_passes)
    counter_values = _counter_values(profile_passes)
    return {
        "available_count": len(completed),
        "missing_count": sum(item.get("status") == "missing" for item in entries),
        "invalid_placement_count": sum(
            item.get("status") == "invalid_placement" for item in entries
        ),
        "blocked_count": sum(item.get("status") == "blocked" for item in entries),
        "not_applicable_count": sum(
            item.get("status") == "not_applicable" for item in entries
        ),
        "host_fingerprints": sorted({
            fingerprint for item in entries
            if (fingerprint := _record_fingerprint(item)) is not None
        }),
        "cost_summary": {
            "wall_seconds_median": median(cost_walls) if cost_walls else None,
            "steps_per_second_median": (
                median(steps_per_second) if steps_per_second else None
            ),
        },
        "attribution_summary": {
            "step_profile_median": {
                phase: median(values) for phase, values in phase_totals.items()
            },
            **{
                f"{key}_median": median(values) if values else None
                for key, values in counter_values.items()
            },
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


def _build_comparable_groups(
    records: dict[tuple[str, str, int], dict[str, Any]]
) -> list[dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}
    for (arm, scale, seed), record in records.items():
        if not _completed_pass_pair(record):
            continue
        fingerprint = _record_fingerprint(record)
        if fingerprint is None:
            continue
        group = groups.setdefault(
            fingerprint, {"host_fingerprint": fingerprint, "arms": {}, "records": []}
        )
        group["arms"].setdefault(arm, set()).add(scale)
        group["records"].append({
            "arm": arm, "scale": scale, "seed": seed,
            "raw_result_file": record.get("raw_result_file"),
        })
    return [
        {
            "host_fingerprint": fingerprint,
            "arms": [
                {"arm": arm, "scales": sorted(scales)}
                for arm, scales in sorted(group["arms"].items())
            ],
            "records": sorted(
                group["records"],
                key=lambda item: (item["arm"], item["scale"], item["seed"])
            ),
        }
        for fingerprint, group in sorted(groups.items())
    ]


def merge_results(manifest_path: Path, result_paths: list[Path]) -> dict[str, Any]:
    """Merge raw records, preserving every missing or invalid entry."""
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
        "schema_version": 3, "scales": manifest["scales"],
        "seeds": manifest["seeds"], "baselines": manifest["baselines"],
        "reporting": manifest.get("reporting", {}), "arms": {},
    }
    for arm, arm_info in manifest["arms"].items():
        current_arm_info = _current_arm_definition(arm, arm_info)
        arm_result = {
            "axis": current_arm_info["axis"],
            "status": current_arm_info.get("status", "runnable"),
            "blocked_reason": current_arm_info.get("blocked_reason"),
            "accepted_placements": sorted(
                current_arm_info.get("accepted_placements", [])
            ),
            "scales": {},
        }
        applicable_scales = set(
            arm_info.get("scales", arm_info.get("configs", {}).keys())
        )
        for scale in manifest["scales"]:
            entries: list[dict[str, Any]] = []
            for seed in manifest["seeds"]:
                record = records.get((arm, scale, int(seed)))
                if record is not None:
                    _revalidate_record_placement(record, current_arm_info)
                    entries.append(record)
                elif scale not in applicable_scales:
                    entries.append({
                        "status": "not_applicable", "arm": arm, "scale": scale,
                        "seed": seed, "reason": "arm is declared only at "
                        + ", ".join(sorted(applicable_scales)),
                    })
                elif current_arm_info.get("status") == "blocked":
                    entries.append({
                        "status": "blocked", "arm": arm, "scale": scale,
                        "seed": seed,
                        "blocked_reason": current_arm_info["blocked_reason"],
                    })
                else:
                    entries.append({
                        "status": "missing", "arm": arm, "scale": scale,
                        "seed": seed, "reason": "raw result was not supplied",
                    })
            arm_result["scales"][scale] = {
                "records": entries, **_summary_for_entries(entries)
            }
        merged["arms"][arm] = arm_result
    merged["comparable_arm_groups"] = _build_comparable_groups(records)
    merged["invalid_placement_count"] = sum(
        scale_data["invalid_placement_count"]
        for arm_data in merged["arms"].values()
        for scale_data in arm_data["scales"].values()
    )
    return merged


def _record_chemistry_seconds(record: dict[str, Any]) -> float | None:
    if not _completed_pass_pair(record):
        return None
    value = (
        _profile_pass(record).get("provenance", {})
        .get("step_profile", {}).get("chemistry_s")
    )
    return float(value) if isinstance(value, (int, float)) else None


def _in_comparable_group(
    merged: dict[str, Any], records: list[dict[str, Any]], fingerprint: str
) -> bool:
    required_keys = {
        (record.get("arm"), record.get("scale"), int(record["seed"]))
        for record in records
    }
    for group in merged.get("comparable_arm_groups", []):
        if group.get("host_fingerprint") != fingerprint:
            continue
        group_keys = {
            (record.get("arm"), record.get("scale"), int(record["seed"]))
            for record in group.get("records", [])
        }
        if required_keys <= group_keys:
            return True
    return False


def _same_group_ratio(
    merged: dict[str, Any], arm: str, scale: str, seed: int, baseline: str
) -> str | None:
    current = merged["arms"].get(arm, {}).get("scales", {}).get(scale, {})
    base = merged["arms"].get(baseline, {}).get("scales", {}).get(scale, {})
    current_record = next(
        (item for item in current.get("records", []) if item.get("seed") == seed),
        None,
    )
    if current_record is None:
        return None
    fingerprint = _record_fingerprint(current_record)
    if fingerprint is None:
        return None
    base_record = next(
        (
            item for item in base.get("records", [])
            if item.get("seed") == seed
            and _record_fingerprint(item) == fingerprint
        ),
        None,
    )
    if base_record is None or not _in_comparable_group(
        merged, [current_record, base_record], fingerprint
    ):
        return None
    numerator = _record_chemistry_seconds(current_record)
    denominator = _record_chemistry_seconds(base_record or {})
    if numerator is None or denominator in (None, 0.0):
        return None
    ratio = numerator / denominator
    return (
        f"chemistry_s ratio {arm}/{baseline}={ratio:.6g} "
        f"(scale={scale}, seeds=1, host={fingerprint})"
    )


def _completed_records(scale_data: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        record for record in scale_data.get("records", [])
        if _completed_pass_pair(record)
    ]


def _aggregate_ratio(
    merged: dict[str, Any], arm: str, scale: str, baseline: str
) -> dict[str, Any]:
    current = merged["arms"].get(arm, {}).get("scales", {}).get(scale, {})
    base = merged["arms"].get(baseline, {}).get("scales", {}).get(scale, {})
    current_records = _completed_records(current)
    base_records = _completed_records(base)
    result = {
        "arm": arm,
        "baseline": baseline,
        "scale": scale,
        "arm_seed_count": len(current_records),
        "baseline_seed_count": len(base_records),
    }
    if not current_records or not base_records:
        result["reason"] = (
            "completed records unavailable "
            f"(arm seeds={len(current_records)}, "
            f"baseline seeds={len(base_records)})"
        )
        return result
    current_seeds = {int(record["seed"]) for record in current_records}
    base_seeds = {int(record["seed"]) for record in base_records}
    if current_seeds != base_seeds:
        result["reason"] = (
            "completed seed sets differ "
            f"(arm={sorted(current_seeds)}, baseline={sorted(base_seeds)})"
        )
        return result
    fingerprints = {
        _record_fingerprint(record)
        for record in (*current_records, *base_records)
    }
    if None in fingerprints or len(fingerprints) != 1:
        current_fingerprints = sorted(
            {_record_fingerprint(item) for item in current_records},
            key=str,
        )
        base_fingerprints = sorted(
            {_record_fingerprint(item) for item in base_records},
            key=str,
        )
        result["reason"] = (
            "completed records do not share one host fingerprint "
            f"(arm={current_fingerprints}, baseline={base_fingerprints})"
        )
        return result
    fingerprint = next(iter(fingerprints))
    if not _in_comparable_group(
        merged, [*current_records, *base_records], fingerprint
    ):
        result["reason"] = (
            "completed records are not contained in one comparable host group"
        )
        return result
    current_chemistry = [
        value for record in current_records
        if (value := _record_chemistry_seconds(record)) is not None
    ]
    base_chemistry = [
        value for record in base_records
        if (value := _record_chemistry_seconds(record)) is not None
    ]
    current_total = [
        float(value) for record in current_records
        if isinstance(
            value := _cost_pass(record).get("wall_seconds"),
            (int, float),
        )
    ]
    base_total = [
        float(value) for record in base_records
        if isinstance(
            value := _cost_pass(record).get("wall_seconds"),
            (int, float),
        )
    ]
    if not current_chemistry or not base_chemistry:
        result["reason"] = "chemistry_s is absent from completed profile passes"
        return result
    chemistry_denominator = median(base_chemistry)
    if not chemistry_denominator:
        result["reason"] = "baseline chemistry_s median is zero"
        return result
    result["host_fingerprint"] = fingerprint
    result["chemistry_ratio"] = median(current_chemistry) / chemistry_denominator
    if current_total and base_total and median(base_total):
        result["total_wall_ratio"] = median(current_total) / median(base_total)
    else:
        result["total_wall_ratio"] = None
    return result


def _report_record_row(
    merged: dict[str, Any],
    arm: str,
    scale: str,
    record: dict[str, Any],
    baseline: str | None,
) -> str:
    status = record.get("status", "unknown")
    reason = (
        record.get("blocked_reason")
        or record.get("invalid_placement_reason")
        or record.get("reason")
        or ""
    )
    status_text = f"{status}: {reason}" if reason else status
    revalidation = record.get("placement_revalidated")
    if isinstance(revalidation, dict):
        status_text += (
            " (placement revalidated: "
            f"stored={revalidation.get('stored')}, "
            f"rederived={revalidation.get('rederived')})"
        )
    chemistry = _record_chemistry_seconds(record)
    total = _cost_pass(record).get("wall_seconds")
    precision = record.get("precision_summary") or "not recorded; see raw provenance"
    raw_path = record.get("raw_result_file", "—")
    ratio = (
        _same_group_ratio(merged, arm, scale, int(record["seed"]), baseline)
        if baseline and arm != baseline
        else None
    )
    cells = [
        arm, scale, str(record.get("seed", "—")), status_text,
        f"{chemistry:.6g}" if chemistry is not None else "—",
        f"{float(total):.6g}" if isinstance(total, (int, float)) else "—",
        str(precision).replace("|", "\\|"),
        str(raw_path).replace("|", "\\|"), ratio or "—",
    ]
    return "| " + " | ".join(cells) + " |"


def _aggregate_report_row(
    merged: dict[str, Any],
    arm: str,
    scale: str,
    baseline: str,
) -> str:
    aggregate = _aggregate_ratio(merged, arm, scale, baseline)
    seed_counts = (
        f"{aggregate['arm_seed_count']}/{aggregate['baseline_seed_count']}"
    )
    label = (
        f"{arm}/{baseline} (scale={scale}, "
        f"seeds={seed_counts}, arms={arm}/{baseline})"
    )
    if "chemistry_ratio" in aggregate:
        chemistry_ratio = f"{label}={aggregate['chemistry_ratio']:.6g}"
        total_ratio = aggregate.get("total_wall_ratio")
        total_text = (
            f"{label}={total_ratio:.6g}"
            if total_ratio is not None else "unavailable"
        )
        status_text = "available"
    else:
        chemistry_ratio = "unavailable"
        total_text = "unavailable"
        status_text = f"unavailable: {aggregate['reason']}"
    cells = [
        arm, scale, seed_counts, f"{arm}/{baseline}",
        chemistry_ratio, total_text, status_text,
    ]
    return "| " + " | ".join(cells) + " |"


def render_report(merged: dict[str, Any]) -> str:
    """Render Markdown using only a merged-results object."""
    lines = [
        "# GPU cost/benefit benchmark report", "",
        "Chemistry phase is the primary GPU quantity; total wall time is secondary.",
        "",
        (
            "| Arm | Scale | Seed | Status / reason | Cost: chemistry_s (s) | "
            "Cost: total wall (s) | Precision / observations | Raw result file | "
            "Ratio |"
        ),
        "|---|---|---:|---|---:|---:|---|---|---|",
    ]
    baselines = merged.get("baselines", {})
    for arm, arm_data in merged.get("arms", {}).items():
        baseline = baselines.get(arm_data.get("axis"))
        for scale, scale_data in arm_data.get("scales", {}).items():
            for record in scale_data.get("records", []):
                lines.append(
                    _report_record_row(merged, arm, scale, record, baseline)
                )
    lines.extend([
        "", "## Comparable host groups", "",
        "Ratios are emitted only for records sharing one host fingerprint.", "",
    ])
    for group in merged.get("comparable_arm_groups", []):
        arm_names = ", ".join(item["arm"] for item in group.get("arms", []))
        lines.append(f"- `{group['host_fingerprint']}`: {arm_names}")
    lines.extend([
        "", "## Aggregate ratios", "",
        (
            "Chemistry is primary; total wall time is explicitly secondary. "
            "Ratios require matching completed seed sets on one host fingerprint."
        ),
        "",
        (
            "| Arm | Scale | Seed counts (arm/baseline) | Arms | "
            "Chemistry ratio (primary) | Total wall ratio (secondary) | "
            "Status / reason |"
        ),
        "|---|---|---:|---|---:|---:|---|",
    ])
    for arm, arm_data in merged.get("arms", {}).items():
        baseline = baselines.get(arm_data.get("axis"))
        if not baseline or arm == baseline:
            continue
        for scale in arm_data.get("scales", {}):
            lines.append(_aggregate_report_row(merged, arm, scale, baseline))
    return "\n".join(lines) + "\n"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--base", type=Path, required=True)
    generate.add_argument("--output-dir", type=Path, required=True)
    generate.add_argument(
        "--scale", action="append", nargs=2, metavar=("NAME", "PATH"), required=True
    )
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
    report = subparsers.add_parser("report")
    report.add_argument(
        "--merged", "--input", dest="merged", type=Path, required=True
    )
    report.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    if args.command == "generate":
        generate_configs(args.base, args.output_dir, dict(args.scale))
    elif args.command == "run":
        run_one_arm(
            args.manifest, args.arm, args.scale, args.seed, args.binary,
            args.output_dir, mpirun=args.mpirun, mpi_ranks=args.mpi_ranks
        )
    elif args.command == "merge":
        _write_json(args.output, merge_results(args.manifest, args.results))
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            render_report(_read_json(args.merged)), encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
