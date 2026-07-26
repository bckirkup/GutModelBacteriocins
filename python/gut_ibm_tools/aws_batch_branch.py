"""Branch a new AWS Batch job from an immutable midstream restart artifact.

Science contract (Tier 2): restored state includes agents + lineage + genome +
chemical grid + clock. RNG is reseeded from the job ``seed`` (not serialized).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .batch_config import (
    BatchConfigError,
    load_simulation_template,
    strip_metadata_keys,
)
from .path_utils import PathValidationError, validate_input_path, validate_path_syntax

AWS_CLI = "aws"
S3_SCHEME = "s3://"
HDF5_PLACEHOLDER = "output.h5"

S3GetText = Callable[[str], str]
S3Put = Callable[[str, str], None]
Submit = Callable[[list[str]], str]


@dataclass(frozen=True)
class BranchResult:
    input_uri: str
    checkpoint_uri: str
    output_uri: str
    checkpoint_prefix: str
    job_id: str | None
    config: dict[str, Any]


def _default_s3_get_text(s3_uri: str) -> str:
    result = subprocess.run(
        [AWS_CLI, "s3", "cp", s3_uri, "-"],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def _default_s3_put(content: str, s3_uri: str) -> None:
    subprocess.run(
        [AWS_CLI, "s3", "cp", "-", s3_uri],
        input=content,
        text=True,
        check=True,
    )


def _default_submit(argv: list[str]) -> str:
    result = subprocess.run(
        argv,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def _require_s3_uri(label: str, value: str) -> str:
    if not value.startswith(S3_SCHEME):
        raise BatchConfigError(f"{label} must be an s3:// URI, got: {value}")
    return value.rstrip("/")


def resolve_restart_uri(
    from_uri: str,
    *,
    s3_get_text: S3GetText = _default_s3_get_text,
) -> tuple[str, dict[str, Any] | None]:
    """Resolve ``--from`` to a concrete restart HDF5 URI and optional latest.json."""
    uri = from_uri.strip()
    if not uri.startswith(S3_SCHEME):
        raise BatchConfigError(f"--from must be an s3:// URI, got: {uri}")
    if uri.endswith(".h5"):
        return uri, None
    if uri.endswith(("latest.json", "/latest.json")):
        payload = json.loads(s3_get_text(uri))
        restart_uri = str(payload.get("uri") or "")
        if not restart_uri.startswith(S3_SCHEME):
            raise BatchConfigError(f"latest.json missing usable uri field: {uri}")
        return restart_uri, payload
    # Treat as checkpoint prefix.
    latest = uri.rstrip("/") + "/latest.json"
    payload = json.loads(s3_get_text(latest))
    restart_uri = str(payload.get("uri") or "")
    if not restart_uri.startswith(S3_SCHEME):
        raise BatchConfigError(f"latest.json missing usable uri field: {latest}")
    return restart_uri, payload


def apply_overlay(
    base: dict[str, Any],
    *,
    overlays: dict[str, Any],
    seed: int | None,
    total_time: float | None,
) -> dict[str, Any]:
    """Merge branch overlays onto a simulation config (strip metadata keys)."""
    cfg = strip_metadata_keys(dict(base))
    for key, value in overlays.items():
        cfg[key] = value
    if seed is not None:
        cfg["seed"] = seed
    if total_time is not None:
        cfg["total_time"] = total_time
    cfg["hdf5_file"] = HDF5_PLACEHOLDER
    cfg["gpu_enabled"] = True
    # entry.sh injects restart.directory under WORK; keep enabled for further forks.
    restart = cfg.get("restart")
    if not isinstance(restart, dict):
        restart = {}
    restart["enabled"] = True
    restart.setdefault("interval_steps", 60)
    restart.setdefault("directory", "restart")
    cfg["restart"] = restart
    return cfg


def _load_overlay_file(path: str | Path) -> dict[str, Any]:
    validated = validate_input_path(path)
    return strip_metadata_keys(load_simulation_template(validated))


def _parse_set_overrides(items: list[str]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for item in items:
        if "=" not in item:
            raise BatchConfigError(f"--set expects key=value, got: {item}")
        key, raw = item.split("=", 1)
        key = key.strip()
        raw = raw.strip()
        if not key:
            raise BatchConfigError(f"--set missing key in: {item}")
        try:
            out[key] = json.loads(raw)
        except json.JSONDecodeError:
            out[key] = raw
    return out


def _sanitize_job_name(stem: str) -> str:
    cleaned = "".join(ch if (ch.isalnum() or ch in "-_") else "-" for ch in stem)
    cleaned = cleaned.strip("-_") or "gutibm-branch"
    return cleaned[:128]


def _extract_job_id(stdout: str) -> str | None:
    stdout = stdout.strip()
    if not stdout:
        return None
    try:
        payload = json.loads(stdout)
    except json.JSONDecodeError:
        return None
    job_id = payload.get("jobId")
    return str(job_id) if job_id is not None else None


def branch_job(
    *,
    from_uri: str,
    out_prefix: str,
    job_queue: str,
    job_definition: str,
    overlay_path: str | Path | None = None,
    set_overrides: list[str] | None = None,
    seed: int | None = None,
    total_time: float | None = None,
    job_name: str | None = None,
    dry_run: bool = False,
    require_gpu: bool = True,
    s3_get_text: S3GetText = _default_s3_get_text,
    s3_put: S3Put = _default_s3_put,
    submit: Submit = _default_submit,
) -> BranchResult:
    """Build + optionally submit a single Batch job forked from a restart URI."""
    restart_uri, latest = resolve_restart_uri(from_uri, s3_get_text=s3_get_text)
    out_root = _require_s3_uri("--out-prefix", out_prefix)
    input_uri = f"{out_root}/input.json"
    output_uri = f"{out_root}/output.h5.gz"
    ckpt_prefix = f"{out_root}/ckpt/"

    base: dict[str, Any] = {}
    if overlay_path is not None:
        base = _load_overlay_file(overlay_path)
    overlays = _parse_set_overrides(set_overrides or [])
    config = apply_overlay(base, overlays=overlays, seed=seed, total_time=total_time)
    provenance = {
        "parent_restart_uri": restart_uri,
        "parent_latest": latest,
        "fidelity": "tier2_agents_grid",
        "rng": "reseeded_from_config_seed",
    }
    config["_branch_provenance"] = provenance

    if not dry_run:
        s3_put(json.dumps(config, indent=2) + "\n", input_uri)

    environment = [
        {"name": "INPUT_S3_URI", "value": input_uri},
        {"name": "OUTPUT_S3_URI", "value": output_uri},
        {"name": "CHECKPOINT_S3_PREFIX", "value": ckpt_prefix},
        {"name": "CHECKPOINT_S3_URI", "value": restart_uri},
        {"name": "STATUS_S3_URI", "value": f"{out_root}/status.json"},
        {"name": "MPI_RANKS", "value": "1"},
    ]
    if require_gpu:
        environment.append({"name": "REQUIRE_GPU", "value": "1"})

    job_id: str | None = None
    if not dry_run:
        name = _sanitize_job_name(job_name or Path(out_root).name or "gutibm-branch")
        argv = [
            AWS_CLI,
            "batch",
            "submit-job",
            "--job-name",
            name,
            "--job-queue",
            job_queue,
            "--job-definition",
            job_definition,
            "--container-overrides",
            json.dumps({"environment": environment}),
        ]
        job_id = _extract_job_id(submit(argv))

    return BranchResult(
        input_uri=input_uri,
        checkpoint_uri=restart_uri,
        output_uri=output_uri,
        checkpoint_prefix=ckpt_prefix,
        job_id=job_id,
        config=config,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fork a Batch job from an immutable restart (latest.json or step_*.h5)."
        ),
    )
    parser.add_argument(
        "--from",
        dest="from_uri",
        required=True,
        help="s3://.../latest.json, s3://.../ckpt/, or s3://.../step_NNNNNN.h5",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="s3://bucket/campaign/fork_name (input/output/ckpt land underneath)",
    )
    parser.add_argument("--job-queue", required=True)
    parser.add_argument("--job-definition", required=True)
    parser.add_argument(
        "--overlay",
        default=None,
        help="Simulation JSON whose parameters become the fork base",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        help="Overlay key=value (JSON value or string); repeatable",
    )
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--total-time", type=float, default=None)
    parser.add_argument("--job-name", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--no-require-gpu",
        action="store_true",
        help="Do not set REQUIRE_GPU=1 on the forked job",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        if args.overlay is not None:
            validate_path_syntax(args.overlay)
        result = branch_job(
            from_uri=args.from_uri,
            out_prefix=args.out_prefix,
            job_queue=args.job_queue,
            job_definition=args.job_definition,
            overlay_path=args.overlay,
            set_overrides=args.set,
            seed=args.seed,
            total_time=args.total_time,
            job_name=args.job_name,
            dry_run=args.dry_run,
            require_gpu=not args.no_require_gpu,
        )
    except (BatchConfigError, PathValidationError, json.JSONDecodeError, OSError) as exc:
        print(f"aws batch branch error: {exc}", file=sys.stderr)
        return 2

    print(f"parent_restart: {result.checkpoint_uri}")
    print(f"input:          {result.input_uri}")
    print(f"output:         {result.output_uri}")
    print(f"ckpt_prefix:    {result.checkpoint_prefix}")
    if args.dry_run:
        print("dry-run: no upload, no submit")
        print(json.dumps(result.config, indent=2))
    elif result.job_id is not None:
        print(f"submitted:      {result.job_id}")
    else:
        print("submitted:      (no jobId in response)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
