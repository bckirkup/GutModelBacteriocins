"""Unit tests for AWS Batch status / progress parsing (no live AWS)."""

from __future__ import annotations

import json
from datetime import UTC, datetime

import pytest

from gut_ibm_tools.aws_batch_status import (
    evaluate_usefulness,
    fetch_job_status,
    format_report,
    parse_progress_line,
    resolve_status_uri,
)

PROGRESS_LINE = (
    "Step 10  t=600s  dt=60s  global_agents=50  local_agents=50  "
    "mu_avg=0.0005  pct=10  rate=1.2  eta_s=4500"
)


def test_parse_progress_line_fields() -> None:
    parsed = parse_progress_line(PROGRESS_LINE)
    assert parsed is not None
    assert parsed["step"] == 10
    assert parsed["sim_time_s"] == pytest.approx(600.0)
    assert parsed["global_agents"] == 50
    assert parsed["mu_avg"] == pytest.approx(0.0005)
    assert parsed["pct"] == pytest.approx(10.0)
    assert parsed["rate"] == pytest.approx(1.2)
    assert parsed["eta_s"] == pytest.approx(4500.0)


def test_parse_progress_line_rejects_noise() -> None:
    assert parse_progress_line("Downloading input: s3://x") is None


def test_resolve_status_uri_prefers_explicit() -> None:
    uri = resolve_status_uri(
        status_uri="s3://b/status.json",
        checkpoint_prefix="s3://b/ckpt/",
        array_index=0,
    )
    assert uri == "s3://b/status.json"


def test_resolve_status_uri_checkpoint_index() -> None:
    uri = resolve_status_uri(checkpoint_prefix="s3://b/ckpt", array_index=3)
    assert uri == "s3://b/ckpt/3/status.json"


def test_evaluate_usefulness_population_and_spot() -> None:
    status = {
        "global_agents": 1,
        "mu_avg": 1e-8,
        "pct": 0.1,
        "wall_elapsed_s": 7200,
        "spot_interruption": True,
        "updated_at": "2026-01-01T00:00:00Z",
    }
    now = datetime(2026, 1, 1, 1, 0, 0, tzinfo=UTC)
    codes = {w.code for w in evaluate_usefulness(status, now=now)}
    assert "population_collapse" in codes
    assert "low_mu" in codes
    assert "slow_progress" in codes
    assert "spot_interruption" in codes
    assert "stale_heartbeat" in codes


def test_evaluate_usefulness_memory_pressure() -> None:
    status = {
        "memory_pressure": True,
        "mem_effective_free_mb": 100,
        "gpu_free_mb": 50,
        "updated_at": "2026-07-25T12:00:00Z",
    }
    codes = {w.code for w in evaluate_usefulness(status)}
    assert "memory_pressure" in codes
    # Pressure flag supersedes low_memory soft warnings.
    assert "low_memory" not in codes


def test_evaluate_usefulness_low_memory_soft_warn() -> None:
    status = {
        "memory_pressure": False,
        "mem_effective_free_mb": 1500,
        "gpu_free_mb": 200,
        "updated_at": "2026-07-25T12:00:00Z",
    }
    codes = {w.code for w in evaluate_usefulness(status)}
    assert "low_memory" in codes
    assert "low_gpu_memory" in codes
    assert "memory_pressure" not in codes


def test_mem_warn_threshold_sensitivity() -> None:
    status = {
        "memory_pressure": False,
        "mem_effective_free_mb": 3000,
        "gpu_free_mb": 8000,
        "updated_at": "2026-07-25T12:00:00Z",
    }
    loose = {w.code for w in evaluate_usefulness(status, mem_warn_mb=2000)}
    tight = {w.code for w in evaluate_usefulness(status, mem_warn_mb=4000)}
    assert "low_memory" not in loose
    assert "low_memory" in tight


def test_resolve_status_uri_output_prefix_sensitivity() -> None:
    bare = resolve_status_uri(output_prefix="s3://b/out")
    indexed = resolve_status_uri(output_prefix="s3://b/out", array_index=2)
    assert bare == "s3://b/out/status.json"
    assert indexed == "s3://b/out/2/status.json"
    assert bare != indexed

def test_fetch_job_status_combines_batch_and_s3() -> None:
    describe = {
        "jobs": [
            {
                "jobId": "job-1",
                "status": "RUNNING",
                "statusReason": None,
                "attempt": 1,
                "attempts": [{"statusReason": "Host EC2 terminated"}],
                "container": {"logStreamName": "gutibm/default/abc"},
                "startedAt": 1_700_000_000_000,
            }
        ]
    }
    status = {
        "state": "running",
        "pct": 25.0,
        "sim_time_s": 1000.0,
        "eta_s": 3000.0,
        "global_agents": 120,
        "mu_avg": 0.001,
        "wall_elapsed_s": 600,
        "updated_at": "2026-07-25T12:00:00Z",
        "spot_interruption": False,
        "resume_from_checkpoint": True,
    }

    def aws_run(argv: list[str]) -> str:
        assert "describe-jobs" in argv
        return json.dumps(describe)

    def s3_get(uri: str) -> str:
        assert uri.endswith("status.json")
        return json.dumps(status)

    report = fetch_job_status(
        "job-1",
        checkpoint_prefix="s3://bucket/ckpt",
        array_index=0,
        aws_run=aws_run,
        s3_get=s3_get,
    )
    assert report.batch_status == "RUNNING"
    assert report.spot_reclaim_hint is True
    assert report.status_json is not None
    assert report.status_json["pct"] == pytest.approx(25.0)
    text = format_report(report)
    assert "pct:" in text
    assert "batch_host_ec2" in text or "Host EC2" in text
