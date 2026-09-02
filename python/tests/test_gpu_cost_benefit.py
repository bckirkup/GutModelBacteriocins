"""Tests for the GPU cost/benefit benchmark harness."""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path

import h5py
import pytest

import gut_ibm_tools.gpu_cost_benefit as benchmark
from gut_ibm_tools.gpu_cost_benefit import (
    ARM_MATRIX,
    generate_configs,
    merge_results,
    render_report,
    run_one_arm,
)
from gut_ibm_tools.path_utils import PathValidationError


def _write_config(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _make_binary(path: Path) -> Path:
    path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def _base_configs(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    base = tmp_path / "base.json"
    _write_config(
        base,
        {
            "total_time": 60,
            "bio_dt": 60,
            "initial_strains": [],
            "domain_x": 0.0006,
            "domain_y": 0.0006,
            "domain_z": 0.0001,
            "grid_dx": 2e-6,
        },
    )
    scales = {}
    for scale in benchmark.SCALES:
        path = tmp_path / f"{scale}.json"
        _write_config(path, {"total_time": 600, "initial_strains": []})
        scales[scale] = path
    return base, scales


def test_matrix_matches_reconciled_scopes() -> None:
    assert benchmark.SCALES == ("s0", "s1", "s2")
    for arm in ("A1", "A2"):
        assert ARM_MATRIX[arm]["accepted_placements"] == {"host"}
    assert ARM_MATRIX["A3"]["accepted_placements"] == {"host_forced_delivery"}
    for arm in ("A4", "A5"):
        assert ARM_MATRIX[arm]["accepted_placements"] == {"device"}
    assert ARM_MATRIX["A6"]["accepted_placements"] == {"device_delivery"}
    assert "blocked_reason" not in ARM_MATRIX["A6"]
    assert ARM_MATRIX["D3"]["status"] == "blocked"
    assert ARM_MATRIX["D3"]["accepted_placements"] == {"host"}
    assert ARM_MATRIX["D3"]["blocked_reason"] == (
        "corrected sealed kernel routes to host fallback by design (#394)"
    )
    assert all("mixed" not in item["accepted_placements"] for item in ARM_MATRIX.values())
    assert all(ARM_MATRIX[arm]["scales"] == ("s1", "s2") for arm in ("A1", "A6"))
    assert all(ARM_MATRIX[arm]["scales"] == ("s0",) for arm in ("B1", "C1", "D1"))


def test_generated_configs_have_parser_compatible_shape(tmp_path: Path) -> None:
    base, scales = _base_configs(tmp_path)
    generated = tmp_path / "generated"
    manifest = generate_configs(base, generated, scales)
    assert set(manifest["scales"]) == set(benchmark.SCALES)
    expected_hdf5 = {
        "enabled": True,
        "schedule": {"summary": 1, "provenance": 1, "grid_species": []},
    }
    for arm, info in manifest["arms"].items():
        for scale, declared in info["configs"].items():
            config = json.loads(Path(declared).read_text(encoding="utf-8"))
            assert config["hdf5"] == expected_hdf5
            assert "grid_species" not in config["hdf5"]
            assert "metabolism" not in config
            for key, value in info["overrides"].items():
                assert config[key] == value
            assert scale in info["scales"]
            assert set(info["accepted_placements"]) == ARM_MATRIX[arm][
                "accepted_placements"
            ]


def test_runner_accepts_expected_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base, scales = _base_configs(tmp_path)
    generated = tmp_path / "generated"
    generate_configs(base, generated, scales)
    monkeypatch.setattr(
        benchmark, "_host_identity", lambda: {"host_fingerprint": "host-a"}
    )
    binary = _make_binary(tmp_path / "gut_ibm")

    def fake_run_pass(
        pass_name, profile_steps, config, config_path, hdf5_path, binary, **kwargs
    ):
        with h5py.File(hdf5_path, "w") as handle:
            provenance = handle.create_group("run_provenance")
            provenance.create_dataset("chemistry_placement", data="host")
            provenance.create_dataset("gpu_compiled", data=0)
        return {
            "name": pass_name,
            "status": "completed",
            "completion_status": "completed",
            "profile_steps": profile_steps,
            "wall_seconds": 1.0,
            "provenance": benchmark._read_provenance(hdf5_path),
        }

    monkeypatch.setattr(benchmark, "_run_pass", fake_run_pass)
    result = run_one_arm(
        generated / "manifest.json", "D1", "s0", 55, binary, tmp_path / "results"
    )
    payload = json.loads(result.read_text(encoding="utf-8"))
    assert payload["status"] == "completed"
    assert all(
        item["placement_status"] == "accepted"
        for item in payload["passes"].values()
    )


@pytest.mark.parametrize(
    ("placement", "gpu_compiled", "expected_reason"),
    [
        ("device", 0, "gpu_compiled must be 1"),
        ("mixed", 0, "chemistry placement"),
    ],
)
def test_runner_rejects_invalid_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    placement: str,
    gpu_compiled: int,
    expected_reason: str,
) -> None:
    base, scales = _base_configs(tmp_path)
    generated = tmp_path / "generated"
    generate_configs(base, generated, scales)
    monkeypatch.setattr(
        benchmark, "_host_identity", lambda: {"host_fingerprint": "host-a"}
    )
    binary = _make_binary(tmp_path / "gut_ibm")

    def fake_run_pass(
        pass_name, profile_steps, config, config_path, hdf5_path, binary, **kwargs
    ):
        with h5py.File(hdf5_path, "w") as handle:
            provenance = handle.create_group("run_provenance")
            provenance.create_dataset("chemistry_placement", data=placement)
            provenance.create_dataset("gpu_compiled", data=gpu_compiled)
        return {
            "name": pass_name,
            "status": "completed",
            "completion_status": "completed",
            "profile_steps": profile_steps,
            "wall_seconds": 1.0,
            "provenance": benchmark._read_provenance(hdf5_path),
        }

    monkeypatch.setattr(benchmark, "_run_pass", fake_run_pass)
    arm = "A4" if placement == "device" else "A1"
    result = run_one_arm(
        generated / "manifest.json", arm, "s1", 55, binary, tmp_path / "results"
    )
    payload = json.loads(result.read_text(encoding="utf-8"))
    assert payload["status"] == "invalid_placement"
    assert payload["invalid_placement_passes"] == ["cost", "profile"]
    assert expected_reason in payload["invalid_placement_reason"] or any(
        expected_reason in item["placement_validation"].get("reason", "")
        for item in payload["passes"].values()
    )
    assert all(
        item["placement_validation"]["accepted_placements"]
        for item in payload["passes"].values()
    )


def test_merge_excludes_invalid_records_and_groups_hosts(tmp_path: Path) -> None:
    manifest = {
        "scales": ["s1"],
        "seeds": [55, 56],
        "baselines": {"A": "A3"},
        "arms": {
            arm: {
                "axis": "A",
                "status": "runnable",
                "scales": ["s1"],
                "accepted_placements": ["host"],
                "configs": {"s1": str(tmp_path / f"{arm}.json")},
            }
            for arm in ("A1", "A3")
        },
    }
    manifest_path = tmp_path / "manifest.json"
    _write_config(manifest_path, manifest)

    def record(arm: str, seed: int, status: str, wall: float) -> Path:
        path = tmp_path / f"{arm}_{seed}.json"
        _write_config(
            path,
            {
                "arm": arm,
                "scale": "s1",
                "seed": seed,
                "status": status,
                "host": {"host_fingerprint": "same-host"},
                "passes": {
                    "cost": {
                        "status": "completed",
                        "wall_seconds": wall,
                        "provenance": {"termination_step": 2},
                    },
                    "profile": {
                        "status": "completed",
                        "wall_seconds": wall,
                        "provenance": {
                            "step_profile": {"chemistry_s": wall / 2}
                        },
                    },
                },
            },
        )
        return path

    valid = record("A1", 55, "completed", 10.0)
    invalid = record("A1", 56, "invalid_placement", 1.0)
    baseline = record("A3", 55, "completed", 20.0)
    merged = merge_results(manifest_path, [valid, invalid, baseline])
    scale_data = merged["arms"]["A1"]["scales"]["s1"]
    assert scale_data["available_count"] == 1
    assert scale_data["invalid_placement_count"] == 1
    assert scale_data["cost_summary"]["wall_seconds_median"] == 10.0
    assert merged["invalid_placement_count"] == 1
    assert merged["comparable_arm_groups"][0]["host_fingerprint"] == "same-host"
    assert {item["arm"] for item in merged["comparable_arm_groups"][0]["records"]} == {
        "A1",
        "A3",
    }


def test_report_has_status_rows_cost_precision_and_same_host_ratio() -> None:
    record = {
        "arm": "A1",
        "scale": "s1",
        "seed": 55,
        "status": "completed",
        "host": {"host_fingerprint": "same"},
        "raw_result_file": "/results/A1.json",
        "passes": {
            "cost": {"status": "completed", "wall_seconds": 10.0},
            "profile": {
                "status": "completed",
                "provenance": {"step_profile": {"chemistry_s": 2.0}},
            },
        },
    }
    baseline = {
        **record,
        "arm": "A3",
        "raw_result_file": "/results/A3.json",
        "passes": {
            "cost": {"status": "completed", "wall_seconds": 20.0},
            "profile": {
                "status": "completed",
                "provenance": {"step_profile": {"chemistry_s": 4.0}},
            },
        },
    }
    merged = {
        "baselines": {"A": "A3"},
        "arms": {
            arm: {
                "axis": "A",
                "scales": {"s1": {"records": [item]}},
            }
            for arm, item in (("A1", record), ("A3", baseline))
        },
        "comparable_arm_groups": [
            {
                "host_fingerprint": "same",
                "arms": [{"arm": "A1"}, {"arm": "A3"}],
                "records": [],
            }
        ],
    }
    report = render_report(merged)
    assert "Cost: chemistry_s" in report
    assert "Precision / observations" in report
    assert "/results/A1.json" in report
    assert "scale=s1" in report
    assert "seeds=1" in report
    assert "A1/A3" in report


def test_runner_rejects_config_outside_manifest_directory(tmp_path: Path) -> None:
    manifest_dir = tmp_path / "manifest"
    manifest_dir.mkdir()
    outside_config = tmp_path / "outside.json"
    _write_config(outside_config, {"total_time": 60})
    manifest_path = manifest_dir / "manifest.json"
    _write_config(
        manifest_path,
        {
            "arms": {
                "D1": {
                    "axis": "D",
                    "status": "runnable",
                    "scales": ["s0"],
                    "accepted_placements": ["host"],
                    "configs": {"s0": str(outside_config)},
                }
            }
        },
    )
    with pytest.raises(PathValidationError, match="outside trusted root"):
        run_one_arm(
            manifest_path, "D1", "s0", 55, tmp_path / "gut_ibm", tmp_path / "results"
        )


@pytest.mark.integration
def test_generated_config_runs_when_binary_is_available(tmp_path: Path) -> None:
    binary_value = os.environ.get("GUTIBM_TEST_BINARY")
    if not binary_value or not Path(binary_value).is_file():
        pytest.skip("GUTIBM_TEST_BINARY is not available")
    base, scales = _base_configs(tmp_path)
    generated = tmp_path / "generated"
    generate_configs(base, generated, {"s0": scales["s0"]})
    result = run_one_arm(
        generated / "manifest.json",
        "D1",
        "s0",
        55,
        Path(binary_value),
        tmp_path / "results",
    )
    assert json.loads(result.read_text(encoding="utf-8"))["status"] == "completed"
