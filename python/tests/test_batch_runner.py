"""Tests for batch runner configuration, manifest, and CLI."""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path

import pytest

from gut_ibm_tools.batch_config import (
    BatchConfigError,
    apply_overrides,
    build_job_config,
    parse_batch_config,
)
from gut_ibm_tools.batch_manifest import (
    JOB_STATUS_DONE,
    create_manifest,
    format_status_line,
    load_manifest,
    pending_jobs,
    save_manifest,
)
from gut_ibm_tools.batch_runner import main as batch_main

FIXTURES = Path(__file__).parent / "fixtures"


def _write_base_config(path: Path) -> None:
    payload = {
        "_comment": "minimal template",
        "total_time": 60,
        "bio_dt": 60,
        "seed": 1,
        "domain_x": 0.0005,
        "domain_y": 0.0005,
        "domain_z": 0.0001,
        "grid_dx": 5e-6,
        "initial_strains": [
            {
                "type": 1,
                "count": 10,
                "mu_max": 5.5e-4,
                "plasmids": ["ColE1"],
                "conjugative": True,
            }
        ],
        "hdf5": {
            "enabled": False
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _write_fake_binary(path: Path) -> None:
    path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _write_batch_config(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


@pytest.fixture
def batch_workspace(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "workspace"
    root.mkdir()
    _write_base_config(root / "base_config.json")
    _write_fake_binary(root / "gut_ibm")
    monkeypatch.chdir(root)
    return root


def test_sweep_expansion_produces_six_jobs(batch_workspace: Path) -> None:
    batch_path = batch_workspace / "batch.json"
    _write_batch_config(
        batch_path,
        {
            "output_dir": "out",
            "base_config": "base_config.json",
            "binary": "gut_ibm",
            "sweep": {
                "seed": [1, 2, 3],
                "kd_colicinE_btuB": [1e-10, 2e-10],
            },
        },
    )
    _, jobs = parse_batch_config(batch_path, repo_root=batch_workspace)
    assert len(jobs) == 6
    ids = [job.job_id for job in jobs]
    assert len(set(ids)) == 6


def test_explicit_runs_preserve_ids(batch_workspace: Path) -> None:
    batch_path = batch_workspace / "batch.json"
    _write_batch_config(
        batch_path,
        {
            "output_dir": "out",
            "base_config": "base_config.json",
            "binary": "gut_ibm",
            "runs": [
                {"id": "baseline", "overrides": {}},
                {"id": "short", "overrides": {"total_time": 120}},
            ],
        },
    )
    _, jobs = parse_batch_config(batch_path, repo_root=batch_workspace)
    assert [job.job_id for job in jobs] == ["baseline", "short"]
    assert jobs[1].overrides["total_time"] == 120


def test_dot_path_override_merges_nested_field(batch_workspace: Path) -> None:
    base = batch_workspace / "base_config.json"
    config = json.loads(base.read_text(encoding="utf-8"))
    apply_overrides(config, {"initial_strains.0.count": 99})
    assert config["initial_strains"][0]["count"] == 99


def test_build_job_config_sets_absolute_hdf5(batch_workspace: Path) -> None:
    hdf5_path = batch_workspace / "jobs" / "run_a" / "output.h5"
    config = build_job_config(
        batch_workspace / "base_config.json",
        {"seed": 42},
        hdf5_file=hdf5_path,
    )
    assert config["seed"] == 42
    assert config["hdf5_file"] == str(hdf5_path.resolve())
    assert "_comment" not in config


def test_sweep_and_runs_mutually_exclusive(batch_workspace: Path) -> None:
    batch_path = batch_workspace / "batch.json"
    _write_batch_config(
        batch_path,
        {
            "output_dir": "out",
            "base_config": "base_config.json",
            "binary": "gut_ibm",
            "sweep": {"seed": [1]},
            "runs": [{"id": "a", "overrides": {}}],
        },
    )
    with pytest.raises(BatchConfigError, match="both 'sweep' and 'runs'"):
        parse_batch_config(batch_path, repo_root=batch_workspace)


def test_manifest_resume_skips_done_jobs(batch_workspace: Path) -> None:
    from gut_ibm_tools.batch_config import JobSpec

    output_dir = batch_workspace / "out"
    jobs = [
        JobSpec(job_id="a", overrides={"seed": 1}),
        JobSpec(job_id="b", overrides={"seed": 2}),
        JobSpec(job_id="c", overrides={"seed": 3}),
        JobSpec(job_id="d", overrides={"seed": 4}),
    ]
    manifest = create_manifest(
        batch_config=batch_workspace / "batch.json",
        output_dir=output_dir,
        jobs=jobs,
    )
    manifest.jobs[0].status = JOB_STATUS_DONE
    manifest.jobs[1].status = JOB_STATUS_DONE
    save_manifest(manifest, output_dir)

    loaded = load_manifest(output_dir)
    remaining = pending_jobs(loaded)
    assert [job.job_id for job in remaining] == ["c", "d"]
    assert "2/4 done" in format_status_line(loaded)


def test_dry_run_cli(batch_workspace: Path, capsys: pytest.CaptureFixture[str]) -> None:
    batch_path = batch_workspace / "batch.json"
    _write_batch_config(
        batch_path,
        {
            "output_dir": "out",
            "base_config": "base_config.json",
            "binary": "gut_ibm",
            "sweep": {"seed": [1, 2]},
        },
    )
    rc = batch_main([str(batch_path), "--dry-run"])
    captured = capsys.readouterr()
    assert rc == 0
    assert "2 jobs" in captured.out
    assert "seed=1" in captured.out or "seed=2" in captured.out


def test_status_cli_requires_manifest(batch_workspace: Path) -> None:
    batch_path = batch_workspace / "batch.json"
    _write_batch_config(
        batch_path,
        {
            "output_dir": "out",
            "base_config": "base_config.json",
            "binary": "gut_ibm",
            "sweep": {"seed": [1]},
        },
    )
    assert batch_main([str(batch_path), "--status"]) == 2


def test_status_cli_reports_pending(batch_workspace: Path, capsys: pytest.CaptureFixture[str]) -> None:
    from gut_ibm_tools.batch_config import JobSpec

    batch_path = batch_workspace / "batch.json"
    output_dir = batch_workspace / "out"
    _write_batch_config(
        batch_path,
        {
            "output_dir": "out",
            "base_config": "base_config.json",
            "binary": "gut_ibm",
            "sweep": {"seed": [1, 2]},
        },
    )
    manifest = create_manifest(
        batch_config=batch_path,
        output_dir=output_dir,
        jobs=[JobSpec(job_id="seed=1", overrides={"seed": 1})],
    )
    save_manifest(manifest, output_dir)

    rc = batch_main([str(batch_path), "--status"])
    captured = capsys.readouterr()
    assert rc == 1
    assert "pending" in captured.out


def test_fixture_batch_sweep_expands_to_four_jobs(batch_workspace: Path) -> None:
    fixture = json.loads((FIXTURES / "batch_sweep.json").read_text(encoding="utf-8"))
    batch_path = batch_workspace / "batch.json"
    _write_batch_config(batch_path, fixture)
    _, jobs = parse_batch_config(batch_path, repo_root=batch_workspace)
    assert len(jobs) == 4


@pytest.mark.integration
def test_integration_single_job_batch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    build_dir = Path(os.environ.get("GUTIBM_BUILD_DIR", "build"))
    binary = build_dir / "gut_ibm"
    if not binary.is_file():
        pytest.skip("gut_ibm binary not built")

    root = tmp_path / "integration"
    root.mkdir()
    monkeypatch.chdir(root)
    _write_base_config(root / "base_config.json")
    batch_path = root / "batch.json"
    _write_batch_config(
        batch_path,
        {
            "output_dir": "out",
            "base_config": "base_config.json",
            "binary": str(binary.resolve()),
            "mpi_ranks": 1,
            "mpirun_args": ["--allow-run-as-root"],
            "runs": [{"id": "smoke", "overrides": {"total_time": 60}}],
        },
    )
    rc = batch_main([str(batch_path)])
    assert rc == 0
    manifest = load_manifest(root / "out")
    assert manifest.jobs[0].status == JOB_STATUS_DONE
    assert Path(manifest.jobs[0].hdf5_path).is_file()
