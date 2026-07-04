"""Parse batch runner JSON and expand parameter sweeps into job definitions."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from itertools import product
from pathlib import Path
from typing import Any

from .path_utils import PathValidationError, validate_input_path

MANIFEST_VERSION = 1
DEFAULT_BINARY = "build/gut_ibm"
DEFAULT_MPI_RANKS = 1
DEFAULT_MPIRUN = "mpirun"
DEFAULT_ON_FAIL = "continue"
ON_FAIL_VALUES = frozenset({"continue", "abort"})


class BatchConfigError(ValueError):
    """Raised when batch configuration is invalid."""


@dataclass(frozen=True)
class ValidateConfig:
    golden: str | None = None
    fish_golden: str | None = None
    check_targets: bool = False
    check_fish_targets: bool = False


@dataclass(frozen=True)
class BatchSettings:
    output_dir: Path
    base_config: Path
    binary: Path
    mpi_ranks: int
    mpirun: str
    mpirun_args: tuple[str, ...]
    on_fail: str
    env: dict[str, str]
    validate: ValidateConfig | None


@dataclass(frozen=True)
class JobSpec:
    job_id: str
    overrides: dict[str, Any] = field(default_factory=dict)


def _load_json_object(path: Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise BatchConfigError(f"batch config must be a JSON object: {path}")
    return payload


def _resolve_repo_path(path: str | Path, *, base_dir: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return (base_dir / candidate).resolve()


def _parse_validate_block(raw: Any) -> ValidateConfig | None:
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise BatchConfigError("'validate' must be an object")
    return ValidateConfig(
        golden=raw.get("golden"),
        fish_golden=raw.get("fish_golden"),
        check_targets=bool(raw.get("check_targets", False)),
        check_fish_targets=bool(raw.get("check_fish_targets", False)),
    )


def parse_batch_config(
    batch_path: str | Path,
    *,
    repo_root: Path | None = None,
    require_binary: bool = True,
) -> tuple[BatchSettings, list[JobSpec]]:
    """Load batch JSON and return settings plus expanded job specifications."""
    validated = validate_input_path(batch_path)
    base_dir = repo_root.resolve() if repo_root is not None else Path.cwd().resolve()
    payload = _load_json_object(validated)

    if "sweep" in payload and "runs" in payload:
        raise BatchConfigError("batch config cannot contain both 'sweep' and 'runs'")

    for required in ("output_dir", "base_config"):
        if required not in payload:
            raise BatchConfigError(f"missing required key: {required}")

    on_fail = str(payload.get("on_fail", DEFAULT_ON_FAIL))
    if on_fail not in ON_FAIL_VALUES:
        raise BatchConfigError(f"on_fail must be one of {sorted(ON_FAIL_VALUES)}")

    env_raw = payload.get("env", {})
    if not isinstance(env_raw, dict):
        raise BatchConfigError("'env' must be an object")
    env = {str(key): str(value) for key, value in env_raw.items()}

    mpirun_args_raw = payload.get("mpirun_args", [])
    if not isinstance(mpirun_args_raw, list):
        raise BatchConfigError("'mpirun_args' must be an array")
    mpirun_args = tuple(str(arg) for arg in mpirun_args_raw)

    settings = BatchSettings(
        output_dir=_resolve_repo_path(payload["output_dir"], base_dir=base_dir),
        base_config=_resolve_repo_path(payload["base_config"], base_dir=base_dir),
        binary=_resolve_repo_path(payload.get("binary", DEFAULT_BINARY), base_dir=base_dir),
        mpi_ranks=int(payload.get("mpi_ranks", DEFAULT_MPI_RANKS)),
        mpirun=str(payload.get("mpirun", DEFAULT_MPIRUN)),
        mpirun_args=mpirun_args,
        on_fail=on_fail,
        env=env,
        validate=_parse_validate_block(payload.get("validate")),
    )

    if settings.mpi_ranks < 1:
        raise BatchConfigError("mpi_ranks must be >= 1")

    if not settings.base_config.is_file():
        raise BatchConfigError(f"base_config not found: {settings.base_config}")
    if require_binary and not settings.binary.is_file():
        raise BatchConfigError(f"binary not found: {settings.binary}")

    if "runs" in payload:
        jobs = _expand_explicit_runs(payload["runs"])
    elif "sweep" in payload:
        jobs = _expand_sweep(payload["sweep"])
    else:
        raise BatchConfigError("batch config must contain 'sweep' or 'runs'")

    return settings, jobs


def _expand_explicit_runs(runs_raw: Any) -> list[JobSpec]:
    if not isinstance(runs_raw, list) or not runs_raw:
        raise BatchConfigError("'runs' must be a non-empty array")
    jobs: list[JobSpec] = []
    seen: set[str] = set()
    for index, entry in enumerate(runs_raw):
        if not isinstance(entry, dict):
            raise BatchConfigError(f"runs[{index}] must be an object")
        job_id = entry.get("id")
        if not isinstance(job_id, str) or not job_id.strip():
            raise BatchConfigError(f"runs[{index}] requires a non-empty 'id'")
        if job_id in seen:
            raise BatchConfigError(f"duplicate run id: {job_id}")
        seen.add(job_id)
        overrides = entry.get("overrides", {})
        if not isinstance(overrides, dict):
            raise BatchConfigError(f"runs[{index}].overrides must be an object")
        jobs.append(JobSpec(job_id=_sanitize_job_id(job_id), overrides=dict(overrides)))
    return jobs


def _expand_sweep(sweep_raw: Any) -> list[JobSpec]:
    if not isinstance(sweep_raw, dict) or not sweep_raw:
        raise BatchConfigError("'sweep' must be a non-empty object")
    keys = sorted(sweep_raw.keys())
    value_lists: list[list[Any]] = []
    for key in keys:
        values = sweep_raw[key]
        if not isinstance(values, list) or not values:
            raise BatchConfigError(f"sweep.{key} must be a non-empty array")
        value_lists.append(values)

    jobs: list[JobSpec] = []
    seen: set[str] = set()
    for combo in product(*value_lists):
        overrides = dict(zip(keys, combo, strict=True))
        job_id = _job_id_from_params(overrides)
        if job_id in seen:
            raise BatchConfigError(f"duplicate sweep job id: {job_id}")
        seen.add(job_id)
        jobs.append(JobSpec(job_id=job_id, overrides=overrides))
    return jobs


def _format_param_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        text = f"{value:g}"
        return text.replace(".", "p").replace("-", "m")
    text = str(value)
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", text)


def _job_id_from_params(params: dict[str, Any]) -> str:
    parts = [f"{key}={_format_param_value(params[key])}" for key in sorted(params)]
    return _sanitize_job_id("_".join(parts))


def _sanitize_job_id(job_id: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._=-]+", "_", job_id.strip())
    if not sanitized:
        raise BatchConfigError("job id is empty after sanitization")
    return sanitized


def load_simulation_template(base_config: Path) -> dict[str, Any]:
    """Load base simulation JSON as a mutable dict."""
    with open(base_config, encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise BatchConfigError(f"base_config must be a JSON object: {base_config}")
    return payload


def _set_dot_path(config: dict[str, Any], dotted_key: str, value: Any) -> None:
    parts = dotted_key.split(".")
    if len(parts) == 1:
        config[parts[0]] = value
        return

    cursor: Any = config
    for part in parts[:-1]:
        if isinstance(cursor, list):
            if not part.isdigit():
                raise BatchConfigError(
                    f"list index must be numeric in dot path '{dotted_key}'",
                )
            index = int(part)
            cursor = cursor[index]
        elif isinstance(cursor, dict):
            if part not in cursor:
                raise BatchConfigError(f"dot path segment not found: {part}")
            cursor = cursor[part]
        else:
            raise BatchConfigError(f"cannot traverse dot path '{dotted_key}'")

    last = parts[-1]
    if isinstance(cursor, list):
        if not last.isdigit():
            raise BatchConfigError(
                f"list index must be numeric in dot path '{dotted_key}'",
            )
        cursor[int(last)] = value
    elif isinstance(cursor, dict):
        cursor[last] = value
    else:
        raise BatchConfigError(f"cannot set dot path '{dotted_key}'")


def apply_overrides(config: dict[str, Any], overrides: dict[str, Any]) -> None:
    """Merge job overrides into a simulation config dict (in place)."""
    for key, value in overrides.items():
        if "." in key:
            _set_dot_path(config, key, value)
        else:
            config[key] = value


def strip_metadata_keys(config: dict[str, Any]) -> dict[str, Any]:
    """Return config without underscore-prefixed documentation keys."""
    return {key: value for key, value in config.items() if not str(key).startswith("_")}


def build_job_config(
    base_config: Path,
    overrides: dict[str, Any],
    *,
    hdf5_file: Path,
) -> dict[str, Any]:
    """Build per-job simulation JSON from template + overrides."""
    config = strip_metadata_keys(load_simulation_template(base_config))
    apply_overrides(config, overrides)
    config["hdf5_file"] = str(hdf5_file.resolve())
    return config


def resolve_batch_path(batch_path: str | Path) -> Path:
    """Validate batch config path and return resolved location."""
    try:
        return validate_input_path(batch_path)
    except PathValidationError as exc:
        raise BatchConfigError(str(exc)) from exc
