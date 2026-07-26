"""Tests for midstream branch / re-init helper."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from gut_ibm_tools.aws_batch_branch import (
    apply_overlay,
    branch_job,
    resolve_restart_uri,
)
from gut_ibm_tools.batch_config import BatchConfigError

REPO_ROOT = Path(__file__).resolve().parents[2]
BURNIN = (
    REPO_ROOT
    / "experiments/diversity_campaign/stage3_campaign/3a_burnin.json"
)
KD = (
    REPO_ROOT
    / "experiments/diversity_campaign/stage3_campaign/3_kd_sweep_1e-7.json"
)


class _Recorder:
    def __init__(self) -> None:
        self.gets: list[str] = []
        self.puts: list[tuple[str, str]] = []
        self.submits: list[list[str]] = []
        self.latest = {
            "step": "step_000120",
            "uri": "s3://bucket/burnin/ckpt/step_000120.h5",
            "size_bytes": 12345,
            "sha256": "abc",
            "fidelity": "tier2_agents_grid",
        }

    def s3_get_text(self, uri: str) -> str:
        self.gets.append(uri)
        if uri.endswith("latest.json"):
            return json.dumps(self.latest)
        raise AssertionError(f"unexpected get: {uri}")

    def s3_put(self, content: str, uri: str) -> None:
        self.puts.append((content, uri))

    def submit(self, argv: list[str]) -> str:
        self.submits.append(argv)
        return json.dumps({"jobId": "branch-job-1"})


def test_resolve_restart_uri_from_latest() -> None:
    rec = _Recorder()
    uri, payload = resolve_restart_uri(
        "s3://bucket/burnin/ckpt/latest.json",
        s3_get_text=rec.s3_get_text,
    )
    assert uri == "s3://bucket/burnin/ckpt/step_000120.h5"
    assert payload is not None
    assert payload["step"] == "step_000120"


def test_resolve_restart_uri_from_prefix() -> None:
    rec = _Recorder()
    uri, _ = resolve_restart_uri(
        "s3://bucket/burnin/ckpt",
        s3_get_text=rec.s3_get_text,
    )
    assert uri.endswith("step_000120.h5")
    assert rec.gets[0].endswith("/latest.json")


def test_resolve_restart_uri_direct_h5() -> None:
    uri, payload = resolve_restart_uri("s3://bucket/burnin/ckpt/step_000060.h5")
    assert uri.endswith("step_000060.h5")
    assert payload is None


def test_apply_overlay_sets_seed_and_restart() -> None:
    base = {"total_time": 100.0, "seed": 1, "kd_corrinoid_btuB": 1e-9}
    cfg = apply_overlay(base, overlays={"kd_corrinoid_btuB": 1e-7}, seed=2001, total_time=604800)
    assert cfg["seed"] == 2001
    assert cfg["total_time"] == pytest.approx(604800)
    assert cfg["kd_corrinoid_btuB"] == pytest.approx(1e-7)
    assert cfg["restart"]["enabled"] is True
    assert cfg["gpu_enabled"] is True


def test_branch_job_dry_run_and_sensitivity(tmp_path: Path) -> None:
    rec = _Recorder()
    result = branch_job(
        from_uri="s3://bucket/burnin/ckpt/latest.json",
        out_prefix="s3://bucket/forks/kd1e7_seed2001",
        job_queue="gutibm-gpu-campaign",
        job_definition="gutibm-cuda-campaign",
        overlay_path=BURNIN,
        set_overrides=["kd_corrinoid_btuB=1e-7"],
        seed=2001,
        total_time=604800.0,
        dry_run=True,
        s3_get_text=rec.s3_get_text,
        s3_put=rec.s3_put,
        submit=rec.submit,
    )
    assert result.job_id is None
    assert not rec.puts
    assert not rec.submits
    assert result.checkpoint_uri.endswith("step_000120.h5")
    assert result.config["seed"] == 2001
    assert result.config["kd_corrinoid_btuB"] == pytest.approx(1e-7)

    # Sensitivity: different seed must change config fingerprint.
    other = branch_job(
        from_uri="s3://bucket/burnin/ckpt/latest.json",
        out_prefix="s3://bucket/forks/kd1e7_seed2002",
        job_queue="gutibm-gpu-campaign",
        job_definition="gutibm-cuda-campaign",
        overlay_path=BURNIN,
        set_overrides=["kd_corrinoid_btuB=1e-7"],
        seed=2002,
        dry_run=True,
        s3_get_text=rec.s3_get_text,
        s3_put=rec.s3_put,
        submit=rec.submit,
    )
    assert other.config["seed"] != result.config["seed"]


def test_branch_job_submits_with_checkpoint_uri() -> None:
    rec = _Recorder()
    result = branch_job(
        from_uri="s3://bucket/burnin/ckpt/latest.json",
        out_prefix="s3://bucket/forks/run1",
        job_queue="gutibm-gpu-campaign",
        job_definition="gutibm-cuda-campaign",
        overlay_path=KD if KD.exists() else BURNIN,
        seed=3001,
        dry_run=False,
        s3_get_text=rec.s3_get_text,
        s3_put=rec.s3_put,
        submit=rec.submit,
    )
    assert result.job_id == "branch-job-1"
    assert len(rec.puts) == 1
    assert rec.puts[0][1].endswith("/input.json")
    env = json.loads(rec.submits[0][rec.submits[0].index("--container-overrides") + 1])
    names = {e["name"]: e["value"] for e in env["environment"]}
    assert names["CHECKPOINT_S3_URI"].endswith("step_000120.h5")
    assert names["REQUIRE_GPU"] == "1"
    assert names["CHECKPOINT_S3_PREFIX"].endswith("/ckpt/")


def test_resolve_rejects_non_s3() -> None:
    with pytest.raises(BatchConfigError):
        resolve_restart_uri("/tmp/latest.json")
