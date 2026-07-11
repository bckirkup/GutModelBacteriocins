"""CLI batch runner for GutIBM parameter scans and experiment grids."""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from .batch_config import (
    BatchConfigError,
    BatchSettings,
    build_job_config,
    parse_batch_config,
    resolve_batch_path,
)
from .batch_manifest import (
    JOB_STATUS_DONE,
    JOB_STATUS_FAILED,
    JOB_STATUS_INTERRUPTED,
    JOB_STATUS_PENDING,
    JOB_STATUS_RUNNING,
    BatchManifest,
    JobRecord,
    create_manifest,
    format_status_line,
    job_paths,
    load_manifest,
    manifest_exists,
    merge_resume_manifest,
    pending_jobs,
    save_manifest,
)
from .path_utils import PathValidationError, prepare_output_file
from .validation_regression import run_validation

EXIT_INTERRUPTED = 130


class _InterruptState:
    def __init__(self) -> None:
        self.requested = False
        self.current_job: JobRecord | None = None
        self.manifest: BatchManifest | None = None
        self.output_dir: Path | None = None


_INTERRUPT = _InterruptState()


def _handle_signal(signum: int, _frame: Any) -> None:
    _INTERRUPT.requested = True
    signal_name = "SIGINT" if signum == signal.SIGINT else "SIGTERM"
    print(f"\n{signal_name} received; finishing current job then stopping...", file=sys.stderr)


def _register_signal_handlers() -> None:
    signal.signal(signal.SIGINT, _handle_signal)
    signal.signal(signal.SIGTERM, _handle_signal)


def _format_job_label(job: JobRecord) -> str:
    if not job.params:
        return job.job_id
    parts = [f"{key}={job.params[key]}" for key in sorted(job.params)]
    return " ".join(parts)


def _print_dry_run(manifest: BatchManifest, output_dir: Path) -> None:
    print(f"Batch: {len(manifest.jobs)} jobs -> {output_dir}")
    for index, job in enumerate(manifest.jobs, start=1):
        paths = job_paths(output_dir, job.job_id)
        print(
            f"{index:4d}  {job.job_id}  "
            f"config={paths['config_path']}  hdf5={paths['hdf5_path']}",
        )


def _write_json_output(path: Path, payload: dict[str, Any]) -> None:
    out = prepare_output_file(path)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def _write_job_config(settings: BatchSettings, job: JobRecord, output_dir: Path) -> Path:
    paths = job_paths(output_dir, job.job_id)
    paths["job_dir"].mkdir(parents=True, exist_ok=True)
    config = build_job_config(
        settings.base_config,
        job.params,
        hdf5_file=paths["hdf5_path"],
    )
    config_path = paths["config_path"]
    _write_json_output(config_path, config)
    job.config_path = str(config_path.resolve())
    job.hdf5_path = str(paths["hdf5_path"].resolve())
    job.log_path = str(paths["log_path"].resolve())
    return config_path


def _build_command(settings: BatchSettings, config_path: Path) -> list[str]:
    command = [
        settings.mpirun,
        *settings.mpirun_args,
        "-np",
        str(settings.mpi_ranks),
        str(settings.binary),
        str(config_path),
    ]
    return command


def _run_simulation(
    settings: BatchSettings,
    config_path: Path,
    log_path: Path,
) -> subprocess.CompletedProcess[str]:
    command = _build_command(settings, config_path)
    env = {**os.environ, **settings.env}
    with open(log_path, "w", encoding="utf-8") as log_handle:
        return subprocess.run(
            command,
            cwd=settings.binary.parent,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            env=env,
            check=False,
        )


def _maybe_validate(
    settings: BatchSettings,
    job: JobRecord,
) -> list[str]:
    if settings.validate is None or job.hdf5_path is None:
        return []
    h5_path = Path(job.hdf5_path)
    if not h5_path.is_file():
        return [f"HDF5 output missing: {h5_path}"]

    golden = settings.validate.golden
    fish_golden = settings.validate.fish_golden
    _, _, failures = run_validation(
        h5_path,
        golden_path=Path(golden) if golden else None,
        fish_golden_path=Path(fish_golden) if fish_golden else None,
        check_targets=settings.validate.check_targets,
        enforce_fish_targets=settings.validate.check_fish_targets,
    )
    return [str(failure) for failure in failures]


def _status_suffix(
    job: JobRecord,
    *,
    exit_code: int,
    duration_s: float,
) -> str:
    if job.status == JOB_STATUS_DONE and not job.validation_failures:
        return f"OK ({duration_s:.1f}s)"
    if job.status == JOB_STATUS_DONE and job.validation_failures:
        return f"VALIDATION_FAIL ({duration_s:.1f}s)"
    if job.status == JOB_STATUS_INTERRUPTED:
        return "INTERRUPTED"
    return f"FAIL (exit {exit_code})"


def _run_single_job(
    settings: BatchSettings,
    manifest: BatchManifest,
    job: JobRecord,
    output_dir: Path,
    *,
    index: int,
    total: int,
) -> int:
    from datetime import datetime, timezone

    def utc_now() -> str:
        return datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    job.status = JOB_STATUS_RUNNING
    job.started_at = utc_now()
    save_manifest(manifest, output_dir)

    config_path = _write_job_config(settings, job, output_dir)
    paths = job_paths(output_dir, job.job_id)
    if paths["hdf5_path"].exists():
        paths["hdf5_path"].unlink()

    started = time.monotonic()
    result = _run_simulation(settings, config_path, paths["log_path"])
    duration_s = time.monotonic() - started
    job.duration_s = duration_s
    job.finished_at = utc_now()
    job.exit_code = result.returncode

    if _INTERRUPT.requested:
        job.status = JOB_STATUS_INTERRUPTED
    elif result.returncode == 0:
        job.validation_failures = _maybe_validate(settings, job)
        job.status = (
            JOB_STATUS_FAILED if job.validation_failures else JOB_STATUS_DONE
        )
    else:
        job.status = JOB_STATUS_FAILED

    save_manifest(manifest, output_dir)

    label = _format_job_label(job)
    suffix = _status_suffix(job, exit_code=result.returncode, duration_s=duration_s)
    print(f"[{index}/{total}] {label} ... {suffix}", file=sys.stderr)
    if job.validation_failures:
        for failure in job.validation_failures:
            print(f"  validation: {failure}", file=sys.stderr)

    if job.status == JOB_STATUS_FAILED and settings.on_fail == "abort":
        return 1
    if job.status == JOB_STATUS_INTERRUPTED:
        return EXIT_INTERRUPTED
    return 0 if job.status == JOB_STATUS_DONE else 1


def _prepare_manifest(
    batch_path: Path,
    settings: BatchSettings,
    jobs: list,
    *,
    resume: bool,
) -> BatchManifest:
    output_dir = settings.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if resume:
        if not manifest_exists(output_dir):
            raise BatchConfigError(
                f"--resume requires existing manifest in {output_dir}",
            )
        manifest = load_manifest(output_dir)
        return merge_resume_manifest(manifest, jobs)

    if manifest_exists(output_dir):
        raise BatchConfigError(
            f"manifest already exists in {output_dir}; use --resume or a new output_dir",
        )
    return create_manifest(
        batch_config=batch_path,
        output_dir=output_dir,
        jobs=jobs,
    )


def run_batch(
    batch_path: str | Path,
    *,
    resume: bool = False,
    dry_run: bool = False,
    status_only: bool = False,
    repo_root: Path | None = None,
) -> int:
    resolved = resolve_batch_path(batch_path)
    settings, jobs = parse_batch_config(
        resolved,
        repo_root=repo_root,
        require_binary=not dry_run,
    )
    output_dir = settings.output_dir

    if status_only:
        if not manifest_exists(output_dir):
            raise BatchConfigError(f"manifest not found in {output_dir}")
        manifest = load_manifest(output_dir)
        print(format_status_line(manifest))
        summary = manifest.summary
        if summary[JOB_STATUS_FAILED] or summary[JOB_STATUS_INTERRUPTED]:
            return 1
        if summary[JOB_STATUS_PENDING] or summary[JOB_STATUS_RUNNING]:
            return 1
        return 0

    manifest = _prepare_manifest(resolved, settings, jobs, resume=resume)

    if dry_run:
        _print_dry_run(manifest, output_dir)
        return 0

    _register_signal_handlers()
    _INTERRUPT.manifest = manifest
    _INTERRUPT.output_dir = output_dir

    to_run = pending_jobs(manifest)
    total = len(manifest.jobs)
    failures = 0

    for job in to_run:
        if _INTERRUPT.requested:
            break
        _INTERRUPT.current_job = job
        index = manifest.jobs.index(job) + 1
        rc = _run_single_job(
            settings,
            manifest,
            job,
            output_dir,
            index=index,
            total=total,
        )
        if rc == EXIT_INTERRUPTED:
            print(
                f"interrupted after {index}/{total} (resume with --resume)",
                file=sys.stderr,
            )
            return EXIT_INTERRUPTED
        if rc != 0:
            failures += 1
            if settings.on_fail == "abort":
                return 1

    if _INTERRUPT.requested:
        print(
            f"interrupted ({format_status_line(manifest)}; resume with --resume)",
            file=sys.stderr,
        )
        return EXIT_INTERRUPTED

    print(format_status_line(manifest))
    return 1 if failures else 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run GutIBM batch experiments from a JSON manifest.",
    )
    parser.add_argument(
        "batch_config",
        help="Path to batch runner JSON configuration",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from existing batch_manifest.json in output_dir",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Expand jobs and print paths without running simulations",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Print manifest summary and exit",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        return run_batch(
            args.batch_config,
            resume=args.resume,
            dry_run=args.dry_run,
            status_only=args.status,
        )
    except (BatchConfigError, PathValidationError, json.JSONDecodeError) as exc:
        print(f"batch runner error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
