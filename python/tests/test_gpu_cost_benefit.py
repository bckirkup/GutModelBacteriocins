"""Tests for the GPU cost/benefit benchmark harness."""

from __future__ import annotations

import json
from pathlib import Path

import gut_ibm_tools.gpu_cost_benefit as benchmark
import pytest
from gut_ibm_tools.gpu_cost_benefit import (
    ARM_MATRIX,
    generate_configs,
    merge_results,
    run_one_arm,
)
from gut_ibm_tools.path_utils import PathValidationError


def _write_config(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


def _flatten_config(value: dict) -> dict[str, object]:
    flattened: dict[str, object] = {}

    def visit(current: dict, prefix: str = "") -> None:
        for key, item in current.items():
            path = f"{prefix}.{key}" if prefix else key
            if isinstance(item, dict):
                visit(item, path)
            else:
                flattened[path] = item

    visit(value)
    return flattened


def test_generated_arms_change_only_their_declared_axis(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.chdir(tmp_path)
    base_path = tmp_path / "base.json"
    scale_path = tmp_path / "scale.json"
    _write_config(base_path, {"total_time": 60, "initial_strains": []})
    _write_config(scale_path, {"domain_x": 0.001, "initial_strains": []})
    manifest = generate_configs(
        base_path, tmp_path / "generated", {"1e5": scale_path}
    )

    baseline = {
        "A": {
            "metabolism.uptake_limit": "delivery",
            "gpu_enabled": False,
            "image_series_mode": "corrected",
            "image_series_relative_tolerance": 1.0e-10,
            "image_series_max_shells": None,
            "use_fmm": False,
            "toxin.lumen_transfer_length": "inf",
            "toxin.lumen_transfer_basis": "effective",
        },
        "B": {
            "metabolism.uptake_limit": "sherwood",
            "gpu_enabled": False,
            "image_series_mode": "corrected",
            "image_series_relative_tolerance": 1.0e-10,
            "image_series_max_shells": None,
            "use_fmm": False,
        },
        "C": {
            "metabolism.uptake_limit": "sherwood",
            "gpu_enabled": False,
            "image_series_mode": "corrected",
            "image_series_relative_tolerance": 1.0e-10,
            "image_series_max_shells": None,
            "use_fmm": False,
            "toxin.lumen_transfer_length": "inf",
            "toxin.lumen_transfer_basis": "effective",
        },
    }
    axis_keys = {
        "A": {"metabolism.uptake_limit", "gpu_enabled"},
        "B": {
            "image_series_mode",
            "image_series_relative_tolerance",
            "image_series_max_shells",
            "use_fmm",
        },
        "C": {"toxin.lumen_transfer_length", "toxin.lumen_transfer_basis"},
    }
    baseline_arms = {"A": "A3", "B": "B2", "C": "C1"}
    generated = tmp_path / "generated"
    for arm, info in manifest["arms"].items():
        values = dict(baseline[info["axis"]])
        values.update(info["overrides"])
        changed_keys = {
            key for key in values
            if values.get(key) != baseline[info["axis"]].get(key)
        }
        assert changed_keys <= axis_keys[info["axis"]]
        assert ARM_MATRIX[arm]["axis"] == info["axis"]

        baseline_path = generated / baseline_arms[info["axis"]] / "1e5.json"
        arm_path = generated / arm / "1e5.json"
        baseline_config = _flatten_config(json.loads(baseline_path.read_text()))
        arm_config = _flatten_config(json.loads(arm_path.read_text()))
        all_keys = set(baseline_config) | set(arm_config)
        actual_changes = {
            key for key in all_keys
            if baseline_config.get(key) != arm_config.get(key)
        }
        expected_changes = {
            key for key, value in info["overrides"].items()
            if baseline_config.get(key) != value
        }
        assert actual_changes == expected_changes
        assert actual_changes <= axis_keys[info["axis"]]


def test_partial_merge_marks_missing_records_without_filling_defaults(
    tmp_path: Path,
) -> None:
    manifest = {
        "scales": ["1e5"],
        "seeds": [55, 56],
        "baselines": {"A": "A3", "B": "B2", "C": "C1"},
        "arms": {
            "A1": {
                "axis": "A",
                "status": "runnable",
                "configs": {"1e5": str(tmp_path / "A1.json")},
            },
            "A6": {
                "axis": "A",
                "status": "blocked",
                "blocked_reason": "parser refusal",
                "configs": {"1e5": str(tmp_path / "A6.json")},
            },
        },
    }
    manifest_path = tmp_path / "manifest.json"
    result_path = tmp_path / "A1_1e5_seed55.json"
    _write_config(manifest_path, manifest)
    _write_config(
        result_path,
        {
            "arm": "A1",
            "scale": "1e5",
            "seed": 55,
            "status": "completed",
            "wall_seconds": 999.0,
            "passes": {
                "cost": {
                    "status": "completed",
                    "wall_seconds": 10.0,
                    "provenance": {"termination_step": 20},
                },
                "profile": {
                    "status": "completed",
                    "wall_seconds": 11.0,
                    "provenance": {"step_profile": {"biology_s": 3.0}},
                },
            },
            "profiling_wall_difference_seconds": 1.0,
            "profiling_wall_difference_fraction": 0.1,
        },
    )

    merged = merge_results(manifest_path, [result_path])
    records = merged["arms"]["A1"]["scales"]["1e5"]["records"]
    assert records[0]["status"] == "completed"
    assert records[0]["raw_result_file"] == str(result_path.resolve())
    assert merged["arms"]["A1"]["scales"]["1e5"]["cost_summary"] == {
        "wall_seconds_median": 10.0,
        "steps_per_second_median": 2.0,
    }
    assert merged["arms"]["A1"]["scales"]["1e5"][
        "attribution_summary"
    ]["step_profile_median"] == {"biology_s": 3.0}
    assert records[1]["status"] == "missing"
    assert "wall_seconds" not in records[1]
    assert merged["arms"]["A6"]["status"] == "blocked"
    assert merged["arms"]["A6"]["blocked_reason"] == "parser refusal"
    assert all(
        record["status"] == "blocked"
        for record in merged["arms"]["A6"]["scales"]["1e5"]["records"]
    )


def test_runner_separates_cost_and_profile_passes(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    base_path = tmp_path / "base.json"
    scale_path = tmp_path / "scale.json"
    _write_config(base_path, {"total_time": 60, "initial_strains": []})
    _write_config(scale_path, {"domain_x": 0.001, "initial_strains": []})
    generated = tmp_path / "generated"
    generate_configs(base_path, generated, {"1e5": scale_path})
    manifest_path = generated / "manifest.json"
    output_dir = tmp_path / "results"
    calls = []

    def fake_run(command, cwd, check):
        calls.append((command, cwd, check))
        return benchmark.subprocess.CompletedProcess(command, 0)

    clock = iter((10.0, 12.0, 20.0, 23.0))
    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)
    monkeypatch.setattr(benchmark.time, "monotonic", lambda: next(clock))
    monkeypatch.setattr(
        benchmark,
        "_read_provenance",
        lambda path: {
            "termination_step": 4,
            "step_profile": {"biology_s": 1.5},
        },
    )

    result_path = run_one_arm(
        manifest_path,
        "A1",
        "1e5",
        55,
        tmp_path / "gut_ibm",
        output_dir,
    )
    result = json.loads(result_path.read_text(encoding="utf-8"))
    assert len(calls) == 2
    cost_config = json.loads(
        (output_dir / "A1_1e5_seed55.cost.config.json").read_text()
    )
    profile_config = json.loads(
        (output_dir / "A1_1e5_seed55.profile.config.json").read_text()
    )
    assert cost_config["profile_steps"] is False
    assert profile_config["profile_steps"] is True
    assert result["passes"]["cost"]["profile_steps"] is False
    assert result["passes"]["profile"]["profile_steps"] is True
    assert result["profiling_wall_difference_seconds"] == pytest.approx(1.0)
    assert result["profiling_wall_difference_fraction"] == pytest.approx(0.5)


def test_runner_rejects_manifest_config_outside_manifest_directory(
    tmp_path: Path,
) -> None:
    manifest_dir = tmp_path / "manifest"
    manifest_dir.mkdir()
    outside_config = tmp_path / "outside.json"
    _write_config(outside_config, {"total_time": 60})
    manifest_path = manifest_dir / "manifest.json"
    _write_config(
        manifest_path,
        {
            "arms": {
                "A1": {
                    "axis": "A",
                    "status": "runnable",
                    "configs": {"1e5": str(outside_config)},
                }
            }
        },
    )

    with pytest.raises(PathValidationError, match="outside trusted root"):
        run_one_arm(
            manifest_path,
            "A1",
            "1e5",
            55,
            tmp_path / "gut_ibm",
            tmp_path / "results",
        )
