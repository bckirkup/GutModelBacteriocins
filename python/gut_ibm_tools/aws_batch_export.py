"""Export a batch_*.json manifest into an S3 array-job tree and submit it.

Phase 3 of ``docs/AWS_BATCH.md``: expand a ``batch_runner`` manifest into one
per-index ``input.json`` under an S3 prefix (identical expansion to
``batch_runner --dry-run``, since both call :func:`parse_batch_config` /
:func:`build_job_config`), then submit a single AWS Batch array job. Because the
campaign job definition requests ``GPU=1`` on single-GPU instances, Batch places
one run per instance automatically.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from .batch_config import (
    BatchConfigError,
    BatchSettings,
    JobSpec,
    build_job_config,
    parse_batch_config,
    resolve_batch_path,
)
from .path_utils import PathValidationError

AWS = "aws"
S3_SCHEME = "s3://"
INPUT_FILE_NAME = "input.json"
GPU_ENABLED_KEY = "gpu_enabled"
HDF5_FILE_KEY = "hdf5_file"
# entry.sh overwrites hdf5_file with its own work-dir path, so the exported value
# is a stable placeholder rather than a leaked local absolute path.
HDF5_PLACEHOLDER = "output.h5"

ENV_INPUT_PREFIX = "INPUT_S3_PREFIX"
ENV_OUTPUT_PREFIX = "OUTPUT_S3_PREFIX"
ENV_CHECKPOINT_PREFIX = "CHECKPOINT_S3_PREFIX"
ENV_MPI_RANKS = "MPI_RANKS"

# Callable that uploads ``content`` to ``s3_uri`` (injected in tests).
S3Put = Callable[[str, str], None]
# Callable that runs an ``aws batch submit-job`` argv and returns stdout.
Submit = Callable[[list[str]], str]


@dataclass(frozen=True)
class ExportResult:
    array_size: int
    input_uris: list[str]
    configs: list[dict[str, Any]]
    job_id: str | None


def _default_s3_put(content: str, s3_uri: str) -> None:
    # Fixed argv, no shell; content is uploaded via stdin.
    subprocess.run(
        [AWS, "s3", "cp", "-", s3_uri],
        input=content,
        text=True,
        check=True,
    )


def _default_submit(argv: list[str]) -> str:
    # Fixed argv (assembled in _submit_argv), no shell.
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


def build_array_job_config(settings: BatchSettings, job: JobSpec) -> dict[str, Any]:
    """Per-index simulation config, matching ``batch_runner`` expansion.

    Reuses :func:`build_job_config` so sweep/override application is identical to
    the local runner, then pins ``gpu_enabled`` and a stable ``hdf5_file`` name.
    """
    config = build_job_config(
        settings.base_config,
        job.overrides,
        hdf5_file=Path(HDF5_PLACEHOLDER),
    )
    config[HDF5_FILE_KEY] = HDF5_PLACEHOLDER
    config[GPU_ENABLED_KEY] = True
    return config


def _input_uri(input_prefix: str, index: int) -> str:
    return f"{input_prefix}/{index}/{INPUT_FILE_NAME}"


def _container_overrides(
    *,
    input_prefix: str,
    output_prefix: str,
    checkpoint_prefix: str | None,
) -> dict[str, Any]:
    environment = [
        {"name": ENV_INPUT_PREFIX, "value": input_prefix},
        {"name": ENV_OUTPUT_PREFIX, "value": output_prefix},
        {"name": ENV_MPI_RANKS, "value": "1"},
    ]
    if checkpoint_prefix is not None:
        environment.append({"name": ENV_CHECKPOINT_PREFIX, "value": checkpoint_prefix})
    return {"environment": environment}


def _sanitize_job_name(stem: str) -> str:
    cleaned = "".join(ch if (ch.isalnum() or ch in "-_") else "-" for ch in stem)
    cleaned = cleaned.strip("-_") or "gutibm-array"
    return cleaned[:128]


def _submit_argv(
    *,
    job_name: str,
    job_queue: str,
    job_definition: str,
    array_size: int,
    overrides: dict[str, Any],
) -> list[str]:
    return [
        AWS,
        "batch",
        "submit-job",
        "--job-name",
        job_name,
        "--job-queue",
        job_queue,
        "--job-definition",
        job_definition,
        "--array-properties",
        f"size={array_size}",
        "--container-overrides",
        json.dumps(overrides),
    ]


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


def export_array(
    batch_path: str | Path,
    *,
    input_prefix: str,
    output_prefix: str,
    checkpoint_prefix: str | None = None,
    job_queue: str,
    job_definition: str,
    dry_run: bool = False,
    repo_root: Path | None = None,
    s3_put: S3Put = _default_s3_put,
    submit: Submit = _default_submit,
) -> ExportResult:
    """Expand ``batch_path`` to an S3 array tree and (unless dry-run) submit it."""
    resolved = resolve_batch_path(batch_path)
    settings, jobs = parse_batch_config(
        resolved,
        repo_root=repo_root,
        require_binary=False,
    )
    in_prefix = _require_s3_uri("--input-prefix", input_prefix)
    out_prefix = _require_s3_uri("--output-prefix", output_prefix)
    ckpt_prefix = (
        _require_s3_uri("--checkpoint-prefix", checkpoint_prefix)
        if checkpoint_prefix is not None
        else None
    )

    configs: list[dict[str, Any]] = []
    input_uris: list[str] = []
    for index, job in enumerate(jobs):
        config = build_array_job_config(settings, job)
        uri = _input_uri(in_prefix, index)
        configs.append(config)
        input_uris.append(uri)
        if not dry_run:
            s3_put(json.dumps(config, indent=2) + "\n", uri)

    array_size = len(jobs)
    job_id: str | None = None
    if not dry_run:
        argv = _submit_argv(
            job_name=_sanitize_job_name(resolved.stem),
            job_queue=job_queue,
            job_definition=job_definition,
            array_size=array_size,
            overrides=_container_overrides(
                input_prefix=in_prefix,
                output_prefix=out_prefix,
                checkpoint_prefix=ckpt_prefix,
            ),
        )
        job_id = _extract_job_id(submit(argv))

    return ExportResult(
        array_size=array_size,
        input_uris=input_uris,
        configs=configs,
        job_id=job_id,
    )


def _print_summary(result: ExportResult, *, dry_run: bool) -> None:
    print(f"array size: {result.array_size}")
    for uri in result.input_uris:
        print(f"  {uri}")
    if dry_run:
        print("dry-run: no uploads, no job submitted")
    elif result.job_id is not None:
        print(f"submitted array job: {result.job_id}")
    else:
        print("submitted array job (no jobId in response)")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export a batch_*.json manifest into an S3 array job and submit it.",
    )
    parser.add_argument("batch_config", help="Path to batch runner JSON manifest")
    parser.add_argument(
        "--input-prefix",
        required=True,
        help="S3 prefix for per-index inputs, e.g. s3://bucket/campaign/jobs",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="S3 prefix for per-index outputs (passed to entry.sh)",
    )
    parser.add_argument(
        "--checkpoint-prefix",
        default=None,
        help="Optional S3 prefix for Spot-resilient checkpoints",
    )
    parser.add_argument("--job-queue", required=True, help="Batch job queue name")
    parser.add_argument(
        "--job-definition",
        required=True,
        help="Batch job definition name",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Expand and print the array without uploading or submitting",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        result = export_array(
            args.batch_config,
            input_prefix=args.input_prefix,
            output_prefix=args.output_prefix,
            checkpoint_prefix=args.checkpoint_prefix,
            job_queue=args.job_queue,
            job_definition=args.job_definition,
            dry_run=args.dry_run,
        )
    except (BatchConfigError, PathValidationError, json.JSONDecodeError) as exc:
        print(f"aws batch export error: {exc}", file=sys.stderr)
        return 2
    _print_summary(result, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
