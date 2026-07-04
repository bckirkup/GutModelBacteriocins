"""Incremental batch manifest for resumable parameter-scan runs."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .batch_config import MANIFEST_VERSION, BatchConfigError, JobSpec
from .path_utils import prepare_output_file

MANIFEST_NAME = "batch_manifest.json"
JOB_STATUS_DONE = "done"
JOB_STATUS_FAILED = "failed"
JOB_STATUS_INTERRUPTED = "interrupted"
JOB_STATUS_PENDING = "pending"
JOB_STATUS_RUNNING = "running"
RESUMABLE_STATUSES = frozenset({
    JOB_STATUS_PENDING,
    JOB_STATUS_FAILED,
    JOB_STATUS_INTERRUPTED,
})


class ManifestError(ValueError):
    """Raised when manifest state is invalid."""


@dataclass
class JobRecord:
    job_id: str
    params: dict[str, Any] = field(default_factory=dict)
    status: str = JOB_STATUS_PENDING
    config_path: str | None = None
    hdf5_path: str | None = None
    log_path: str | None = None
    exit_code: int | None = None
    duration_s: float | None = None
    started_at: str | None = None
    finished_at: str | None = None
    validation_failures: list[str] = field(default_factory=list)


@dataclass
class BatchManifest:
    version: int
    batch_config: str
    output_dir: str
    created_at: str
    updated_at: str
    jobs: list[JobRecord] = field(default_factory=list)

    @property
    def summary(self) -> dict[str, int]:
        counts = {
            JOB_STATUS_DONE: 0,
            JOB_STATUS_FAILED: 0,
            JOB_STATUS_INTERRUPTED: 0,
            JOB_STATUS_PENDING: 0,
            JOB_STATUS_RUNNING: 0,
        }
        for job in self.jobs:
            counts[job.status] = counts.get(job.status, 0) + 1
        counts["total"] = len(self.jobs)
        return counts

    def job_by_id(self, job_id: str) -> JobRecord | None:
        for job in self.jobs:
            if job.job_id == job_id:
                return job
        return None


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def manifest_path(output_dir: Path) -> Path:
    return output_dir / MANIFEST_NAME


def manifest_exists(output_dir: Path) -> bool:
    return manifest_path(output_dir).is_file()


def _record_from_dict(raw: dict[str, Any]) -> JobRecord:
    validation_failures = raw.get("validation_failures", [])
    if not isinstance(validation_failures, list):
        raise ManifestError("validation_failures must be an array")
    return JobRecord(
        job_id=str(raw["job_id"]),
        params=dict(raw.get("params", {})),
        status=str(raw.get("status", JOB_STATUS_PENDING)),
        config_path=raw.get("config_path"),
        hdf5_path=raw.get("hdf5_path"),
        log_path=raw.get("log_path"),
        exit_code=raw.get("exit_code"),
        duration_s=raw.get("duration_s"),
        started_at=raw.get("started_at"),
        finished_at=raw.get("finished_at"),
        validation_failures=[str(item) for item in validation_failures],
    )


def load_manifest(output_dir: Path) -> BatchManifest:
    path = manifest_path(output_dir)
    if not path.is_file():
        raise ManifestError(f"manifest not found: {path}")
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ManifestError("manifest must be a JSON object")
    jobs_raw = payload.get("jobs", [])
    if not isinstance(jobs_raw, list):
        raise ManifestError("manifest.jobs must be an array")
    return BatchManifest(
        version=int(payload.get("version", MANIFEST_VERSION)),
        batch_config=str(payload["batch_config"]),
        output_dir=str(payload["output_dir"]),
        created_at=str(payload.get("created_at", _utc_now())),
        updated_at=str(payload.get("updated_at", _utc_now())),
        jobs=[_record_from_dict(item) for item in jobs_raw],
    )


def create_manifest(
    *,
    batch_config: Path,
    output_dir: Path,
    jobs: list[JobSpec],
) -> BatchManifest:
    now = _utc_now()
    return BatchManifest(
        version=MANIFEST_VERSION,
        batch_config=str(batch_config.resolve()),
        output_dir=str(output_dir.resolve()),
        created_at=now,
        updated_at=now,
        jobs=[
            JobRecord(job_id=job.job_id, params=dict(job.overrides))
            for job in jobs
        ],
    )


def save_manifest(manifest: BatchManifest, output_dir: Path) -> Path:
    manifest.updated_at = _utc_now()
    path = manifest_path(output_dir)
    out = prepare_output_file(path)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(asdict(manifest), handle, indent=2)
        handle.write("\n")
    return path


def pending_jobs(manifest: BatchManifest) -> list[JobRecord]:
    return [job for job in manifest.jobs if job.status in RESUMABLE_STATUSES]


def format_status_line(manifest: BatchManifest) -> str:
    summary = manifest.summary
    return (
        f"{summary[JOB_STATUS_DONE]}/{summary['total']} done, "
        f"{summary[JOB_STATUS_FAILED]} failed, "
        f"{summary[JOB_STATUS_INTERRUPTED]} interrupted, "
        f"{summary[JOB_STATUS_PENDING]} pending"
    )


def job_paths(output_dir: Path, job_id: str) -> dict[str, Path]:
    job_dir = output_dir / "jobs" / job_id
    return {
        "job_dir": job_dir,
        "config_path": job_dir / "input.json",
        "hdf5_path": job_dir / "output.h5",
        "log_path": job_dir / "run.log",
    }


def ensure_manifest_matches_jobs(
    manifest: BatchManifest,
    jobs: list[JobSpec],
) -> None:
    expected_ids = [job.job_id for job in jobs]
    actual_ids = [job.job_id for job in manifest.jobs]
    if expected_ids != actual_ids:
        raise BatchConfigError(
            "manifest job list does not match batch config; "
            "use a fresh output_dir or delete batch_manifest.json",
        )


def merge_resume_manifest(
    existing: BatchManifest,
    jobs: list[JobSpec],
) -> BatchManifest:
    ensure_manifest_matches_jobs(existing, jobs)
    return existing
