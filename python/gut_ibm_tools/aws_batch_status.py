"""Live AWS Batch job status for GutIBM (Phase Observability).

Combines ``describe-jobs``, optional S3 ``status.json`` heartbeats from
``deploy/aws/entry.sh``, and optional CloudWatch log tails. No live AWS calls
in unit tests — inject callables.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

AWS = "aws"
STATUS_FILE_NAME = "status.json"
HOST_EC2_PREFIX = "Host EC2"
DEFAULT_MU_WARN = 1.0e-5
DEFAULT_STALE_HEARTBEAT_S = 900.0
DEFAULT_MEM_WARN_MB = 4096.0
DEFAULT_GPU_WARN_MB = 1024.0

# Callable: argv -> stdout text (injected in tests).
AwsRunner = Callable[[list[str]], str]
# Callable: s3_uri -> JSON text or raises FileNotFoundError.
S3GetText = Callable[[str], str]

PROGRESS_RE = re.compile(
    r"Step\s+(?P<step>\d+)"
    r".*?\bt=(?P<t>[0-9.eE+-]+)s"
    r".*?\bglobal_agents=(?P<agents>\d+)"
    r".*?\bmu_avg=(?P<mu>[0-9.eE+-]+)"
    r"(?:.*?\bpct=(?P<pct>[0-9.eE+-]+))?"
    r"(?:.*?\brate=(?P<rate>[0-9.eE+-]+))?"
    r"(?:.*?\beta_s=(?P<eta>[0-9.eE+-]+))?",
)


@dataclass(frozen=True)
class UsefulnessWarning:
    code: str
    message: str


@dataclass
class JobStatusReport:
    job_id: str
    batch_status: str | None = None
    status_reason: str | None = None
    attempts: int = 0
    log_stream: str | None = None
    started_at_ms: int | None = None
    stopped_at_ms: int | None = None
    status_json: dict[str, Any] | None = None
    status_uri: str | None = None
    log_tail: list[str] = field(default_factory=list)
    warnings: list[UsefulnessWarning] = field(default_factory=list)
    spot_reclaim_hint: bool = False


def _default_aws_run(argv: list[str]) -> str:
    result = subprocess.run(argv, capture_output=True, text=True, check=True)
    return result.stdout


def _default_s3_get(s3_uri: str) -> str:
    result = subprocess.run(
        [AWS, "s3", "cp", s3_uri, "-"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise FileNotFoundError(s3_uri)
    return result.stdout


def parse_progress_line(line: str) -> dict[str, Any] | None:
    """Parse one gut_ibm stdout progress line into a dict (or None)."""
    match = PROGRESS_RE.search(line.strip())
    if match is None:
        return None
    groups = match.groupdict()
    out: dict[str, Any] = {
        "step": int(groups["step"]),
        "sim_time_s": float(groups["t"]),
        "global_agents": int(groups["agents"]),
        "mu_avg": float(groups["mu"]),
    }
    if groups.get("pct") is not None:
        out["pct"] = float(groups["pct"])
    if groups.get("rate") is not None:
        out["rate"] = float(groups["rate"])
    if groups.get("eta") is not None:
        out["eta_s"] = float(groups["eta"])
    return out


def resolve_status_uri(
    *,
    status_uri: str | None = None,
    checkpoint_prefix: str | None = None,
    output_prefix: str | None = None,
    array_index: str | int | None = None,
) -> str | None:
    """Resolve where entry.sh wrote status.json."""
    if status_uri:
        return status_uri.rstrip("/")
    if checkpoint_prefix:
        prefix = checkpoint_prefix.rstrip("/")
        if array_index is not None and f"/{array_index}" not in prefix + "/":
            prefix = f"{prefix}/{array_index}"
        return f"{prefix}/{STATUS_FILE_NAME}"
    if output_prefix is not None and array_index is not None:
        return f"{output_prefix.rstrip('/')}/{array_index}/{STATUS_FILE_NAME}"
    if output_prefix:
        return f"{output_prefix.rstrip('/')}/{STATUS_FILE_NAME}"
    return None


def evaluate_usefulness(
    status: dict[str, Any] | None,
    *,
    mu_warn: float = DEFAULT_MU_WARN,
    stale_heartbeat_s: float = DEFAULT_STALE_HEARTBEAT_S,
    mem_warn_mb: float = DEFAULT_MEM_WARN_MB,
    gpu_warn_mb: float = DEFAULT_GPU_WARN_MB,
    now: datetime | None = None,
) -> list[UsefulnessWarning]:
    """Triage hints only — not auto-cancel criteria."""
    warnings: list[UsefulnessWarning] = []
    if not status:
        return warnings

    agents = status.get("global_agents")
    if isinstance(agents, (int, float)) and agents <= 1:
        warnings.append(
            UsefulnessWarning(
                "population_collapse",
                f"global_agents={agents} — run may have hit the population stop",
            )
        )

    mu_avg = status.get("mu_avg")
    if isinstance(mu_avg, (int, float)) and mu_avg < mu_warn:
        warnings.append(
            UsefulnessWarning(
                "low_mu",
                f"mu_avg={mu_avg} below warn threshold {mu_warn} (washout risk)",
            )
        )

    pct = status.get("pct")
    wall = status.get("wall_elapsed_s")
    if (
        isinstance(pct, (int, float))
        and isinstance(wall, (int, float))
        and wall > 3600
        and pct < 0.5
    ):
        warnings.append(
            UsefulnessWarning(
                "slow_progress",
                f"only {pct}% complete after {wall:.0f}s wall time",
            )
        )

    updated = status.get("updated_at")
    if isinstance(updated, str) and updated:
        try:
            stamp = datetime.strptime(updated, "%Y-%m-%dT%H:%M:%SZ").replace(
                tzinfo=timezone.utc
            )
            ref = now or datetime.now(timezone.utc)
            age = (ref - stamp).total_seconds()
            if age > stale_heartbeat_s:
                warnings.append(
                    UsefulnessWarning(
                        "stale_heartbeat",
                        f"status.json age {age:.0f}s > {stale_heartbeat_s:.0f}s "
                        "(process may be hung, OOM-killed, or Spot-killed without flush)",
                    )
                )
        except ValueError:
            pass

    if status.get("spot_interruption"):
        warnings.append(
            UsefulnessWarning(
                "spot_interruption",
                "status.json reports spot_interruption=true",
            )
        )

    if status.get("memory_pressure"):
        warnings.append(
            UsefulnessWarning(
                "memory_pressure",
                "status.json reports memory_pressure=true "
                "(graceful stop; resize instance before resume)",
            )
        )
    else:
        mem_free = status.get("mem_effective_free_mb")
        if isinstance(mem_free, (int, float)) and mem_free < mem_warn_mb:
            warnings.append(
                UsefulnessWarning(
                    "low_memory",
                    f"mem_effective_free_mb={mem_free} below warn {mem_warn_mb}",
                )
            )
        gpu_free = status.get("gpu_free_mb")
        if isinstance(gpu_free, (int, float)) and gpu_free < gpu_warn_mb:
            warnings.append(
                UsefulnessWarning(
                    "low_gpu_memory",
                    f"gpu_free_mb={gpu_free} below warn {gpu_warn_mb}",
                )
            )

    return warnings


def _parse_describe_job(payload: dict[str, Any], job_id: str) -> JobStatusReport:
    jobs = payload.get("jobs") or []
    if not jobs:
        return JobStatusReport(job_id=job_id, batch_status=None)
    job = jobs[0]
    attempts = job.get("attempts") or []
    reason = job.get("statusReason")
    spot_hint = isinstance(reason, str) and reason.startswith(HOST_EC2_PREFIX)
    if not spot_hint:
        for attempt in attempts:
            att_reason = attempt.get("statusReason")
            if isinstance(att_reason, str) and att_reason.startswith(HOST_EC2_PREFIX):
                spot_hint = True
                break
    container = job.get("container") or {}
    return JobStatusReport(
        job_id=job.get("jobId", job_id),
        batch_status=job.get("status"),
        status_reason=reason,
        attempts=len(attempts) if attempts else int(job.get("attempt") or 0),
        log_stream=container.get("logStreamName"),
        started_at_ms=job.get("startedAt"),
        stopped_at_ms=job.get("stoppedAt"),
        spot_reclaim_hint=spot_hint,
    )


def fetch_job_status(
    job_id: str,
    *,
    region: str = "us-east-1",
    status_uri: str | None = None,
    checkpoint_prefix: str | None = None,
    output_prefix: str | None = None,
    array_index: str | int | None = None,
    tail_logs: int = 0,
    log_group: str = "/aws/batch/job",
    aws_run: AwsRunner | None = None,
    s3_get: S3GetText | None = None,
    mu_warn: float = DEFAULT_MU_WARN,
) -> JobStatusReport:
    """Build a JobStatusReport for one Batch job id (child or standalone)."""
    runner = aws_run or _default_aws_run
    getter = s3_get or _default_s3_get

    raw = runner(
        [
            AWS,
            "batch",
            "describe-jobs",
            "--jobs",
            job_id,
            "--region",
            region,
            "--output",
            "json",
        ]
    )
    payload = json.loads(raw) if raw.strip() else {"jobs": []}
    report = _parse_describe_job(payload, job_id)

    uri = resolve_status_uri(
        status_uri=status_uri,
        checkpoint_prefix=checkpoint_prefix,
        output_prefix=output_prefix,
        array_index=array_index,
    )
    report.status_uri = uri
    if uri:
        try:
            text = getter(uri)
            report.status_json = json.loads(text)
        except (json.JSONDecodeError, OSError):
            report.status_json = None

    if tail_logs > 0 and report.log_stream:
        try:
            log_raw = runner(
                [
                    AWS,
                    "logs",
                    "get-log-events",
                    "--region",
                    region,
                    "--log-group-name",
                    log_group,
                    "--log-stream-name",
                    report.log_stream,
                    "--limit",
                    str(tail_logs),
                    "--output",
                    "json",
                ]
            )
            events = json.loads(log_raw).get("events") or []
            report.log_tail = [str(e.get("message", "")).rstrip() for e in events]
        except (subprocess.CalledProcessError, json.JSONDecodeError, OSError):
            report.log_tail = []

    report.warnings = evaluate_usefulness(report.status_json, mu_warn=mu_warn)
    if report.spot_reclaim_hint and not any(
        w.code == "spot_interruption" for w in report.warnings
    ):
        report.warnings.append(
            UsefulnessWarning(
                "batch_host_ec2",
                f"Batch statusReason suggests Spot reclaim: {report.status_reason}",
            )
        )
    return report


def format_report(report: JobStatusReport) -> str:
    """One-screen human summary."""
    lines = [
        f"job_id:          {report.job_id}",
        f"batch_status:    {report.batch_status or '(unknown)'}",
        f"status_reason:   {report.status_reason or '-'}",
        f"attempts:        {report.attempts}",
        f"log_stream:      {report.log_stream or '-'}",
        f"spot_reclaim:    {report.spot_reclaim_hint}",
    ]
    if report.status_uri:
        lines.append(f"status_uri:      {report.status_uri}")
    sj = report.status_json
    if sj:
        lines.extend(
            [
                f"sim_state:       {sj.get('state')}",
                f"pct:             {sj.get('pct')}",
                f"sim_time_s:      {sj.get('sim_time_s')}",
                f"eta_s:           {sj.get('eta_s')}",
                f"global_agents:   {sj.get('global_agents')}",
                f"mu_avg:          {sj.get('mu_avg')}",
                f"wall_elapsed_s:  {sj.get('wall_elapsed_s')}",
                f"updated_at:      {sj.get('updated_at')}",
                f"resume:          {sj.get('resume_from_checkpoint')}",
                f"spot_flag:       {sj.get('spot_interruption')}",
                f"mem_free_mb:     {sj.get('mem_effective_free_mb')}",
                f"gpu_free_mb:     {sj.get('gpu_free_mb')}",
                f"mem_pressure:    {sj.get('memory_pressure')}",
            ]
        )
    else:
        lines.append("status.json:     (not found)")
    if report.warnings:
        lines.append("warnings:")
        for warn in report.warnings:
            lines.append(f"  - [{warn.code}] {warn.message}")
    if report.log_tail:
        lines.append("log_tail:")
        for msg in report.log_tail[-20:]:
            lines.append(f"  {msg}")
    return "\n".join(lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize GutIBM AWS Batch job progress and Spot hints.",
    )
    parser.add_argument("job_id", help="Batch job id (child id for array jobs)")
    parser.add_argument("--region", default="us-east-1")
    parser.add_argument("--status-uri", default=None, help="Explicit s3://…/status.json")
    parser.add_argument(
        "--checkpoint-prefix",
        default=None,
        help="S3 checkpoint prefix (status.json lives beside checkpoint.h5)",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="S3 output prefix (fallback location for status.json)",
    )
    parser.add_argument(
        "--array-index",
        default=None,
        help="Array index when resolving prefixes",
    )
    parser.add_argument(
        "--tail-logs",
        type=int,
        default=0,
        help="Optional number of CloudWatch log events to append",
    )
    parser.add_argument(
        "--mu-warn",
        type=float,
        default=DEFAULT_MU_WARN,
        help="Warn when status mu_avg is below this threshold",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        report = fetch_job_status(
            args.job_id,
            region=args.region,
            status_uri=args.status_uri,
            checkpoint_prefix=args.checkpoint_prefix,
            output_prefix=args.output_prefix,
            array_index=args.array_index,
            tail_logs=args.tail_logs,
            mu_warn=args.mu_warn,
        )
    except (subprocess.CalledProcessError, json.JSONDecodeError, OSError) as exc:
        print(f"aws batch status error: {exc}", file=sys.stderr)
        return 2
    print(format_report(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
