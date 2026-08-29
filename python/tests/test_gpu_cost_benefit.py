"""Tests for the GPU cost/benefit benchmark harness."""

from __future__ import annotations

import json
from pathlib import Path

from gut_ibm_tools.gpu_cost_benefit import (
    ARM_MATRIX,
    generate_configs,
    merge_results,
)


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
        {"arm": "A1", "scale": "1e5", "seed": 55, "status": "completed"},
    )

    merged = merge_results(manifest_path, [result_path])
    records = merged["arms"]["A1"]["scales"]["1e5"]["records"]
    assert records[0]["status"] == "completed"
    assert records[0]["raw_result_file"] == str(result_path.resolve())
    assert records[1]["status"] == "missing"
    assert "wall_seconds" not in records[1]
    assert merged["arms"]["A6"]["status"] == "blocked"
    assert merged["arms"]["A6"]["blocked_reason"] == "parser refusal"
    assert all(
        record["status"] == "blocked"
        for record in merged["arms"]["A6"]["scales"]["1e5"]["records"]
    )
