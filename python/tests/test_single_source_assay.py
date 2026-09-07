"""End-to-end tests for experiments/single_source_transport_assay_v1.

The assay generator and analyzer live outside the ``gut_ibm_tools`` package, so
they are loaded from their on-disk path.  Synthetic outputs are built here with
a known radial profile per amplitude arm: the tests check that the analyzer
passes only when its predeclared conditions hold and blocks (rather than
quietly adapting) on missing outputs, missing exact provenance timing, missing
paired snapshot times, and a moved source.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
ASSAY = REPO / "experiments" / "single_source_transport_assay_v1"
GENERATED = ASSAY / "generated"

LYSIS_TIME_S = 240.0
LYSIS_STEP = 4
GRID_STEPS = (6, 8, 10)
DECAY_UM = {0.0: 10.0, 15.0: 8.0, 60.0: 6.0}
SOURCE_FRACTION = 0.5


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def analyze_assay():
    return load_module("assay_analyzer_under_test", ASSAY / "analyze_assay.py")


@pytest.fixture(scope="module")
def contract():
    return json.loads((ASSAY / "assay_contract.json").read_text())


def resolved_config(cfg: dict) -> dict:
    schedule = cfg["hdf5"]["schedule"]
    resolved = {
        "seed": int(cfg["seed"]),
        "bacteriocin.mucin_charge.amplitude": float(cfg["bacteriocin.mucin_charge.amplitude"]),
        "kd_b12_btuB": float(cfg["kd_corrinoid_btuB"]),
        "kd_colicinE_btuB": float(cfg["kd_colicinE_btuB"]),
        "b12.initial_conc": float(cfg["b12_initial_conc"]),
        "burst_release_tau": float(cfg["burst_release_tau"]),
        "hdf5.schedule.grid_species": list(schedule["grid_species"]),
    }
    for key, value in schedule.items():
        if key != "grid_species":
            resolved[f"hdf5.schedule.{key}"] = int(value)
    return resolved


def radial_field(
    cfg: dict, source_m: tuple[float, float, float], decay_um: float, scale: float
):
    dx = float(cfg["grid_dx"])
    shape = tuple(round(float(cfg[f"domain_{axis}"]) / dx) for axis in ("x", "y", "z"))
    centres = [(np.arange(n) + 0.5) * dx for n in shape]
    dxs = centres[0] - source_m[0]
    dys = centres[1] - source_m[1]
    dzs = centres[2] - source_m[2]
    r = np.sqrt(
        dxs[:, None, None] ** 2 + dys[None, :, None] ** 2 + dzs[None, None, :] ** 2
    )
    r_um = np.maximum(r * 1.0e6, 0.5 * dx * 1.0e6)
    return scale * np.exp(-r_um / decay_um) / r_um


def write_output(
    path: Path,
    cfg: dict,
    sha: str,
    *,
    producer: bool,
    source_m: tuple[float, float, float],
    exact_timing: bool = True,
    grid_steps: tuple[int, ...] = GRID_STEPS,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    bio_dt = float(cfg["bio_dt"])
    amplitude = float(cfg["bacteriocin.mucin_charge.amplitude"])
    field = (
        radial_field(cfg, source_m, DECAY_UM[amplitude], 1.0 + amplitude)
        if producer
        else np.zeros(
            tuple(
                round(float(cfg[f"domain_{axis}"]) / float(cfg["grid_dx"]))
                for axis in ("x", "y", "z")
            )
        )
    )
    with h5py.File(path, "w") as h5:
        prov = h5.create_group("run_provenance")
        prov.create_dataset("git_sha", data=np.bytes_(sha))
        prov.create_dataset("chemistry_placement", data=np.bytes_("device_delivery"))
        prov.create_dataset("mpi_rank_count", data=np.int64(1))
        prov.create_dataset(
            "resolved_config", data=np.bytes_(json.dumps(resolved_config(cfg)))
        )
        for step in grid_steps:
            key = f"step_{step:06d}"
            h5.create_dataset(f"summary/{key}/time", data=np.float64(step * bio_dt))
            h5.create_dataset(f"grid/{key}/bacteriocin_BtuB", data=field)
        if not producer:
            h5.create_group("provenance")
            return
        group = h5.create_group(f"provenance/step_{LYSIS_STEP:06d}")
        group.create_dataset("cause", data=np.array([5], dtype=np.int64))
        group.create_dataset("strain", data=np.array([1], dtype=np.int64))
        group.create_dataset("x", data=np.array([source_m[0]], dtype=float))
        group.create_dataset("y", data=np.array([source_m[1]], dtype=float))
        group.create_dataset("z", data=np.array([source_m[2]], dtype=float))
        if exact_timing:
            group.create_dataset("event_step", data=np.array([LYSIS_STEP], dtype=np.int64))
            group.create_dataset(
                "event_time_s", data=np.array([LYSIS_TIME_S], dtype=float)
            )


def build_results(
    root: Path,
    runs: list[dict],
    sha: str,
    *,
    exact_timing: bool = True,
    skip_index: int | None = None,
    drop_grid_step_for_index: int | None = None,
    move_source_for_index: int | None = None,
) -> None:
    for run in runs:
        index = int(run["array_index"])
        if index == skip_index:
            continue
        cfg = json.loads((ASSAY / run["input_relpath"]).read_text())
        domain = [float(cfg[f"domain_{axis}"]) for axis in ("x", "y", "z")]
        source = tuple(SOURCE_FRACTION * length for length in domain)
        if index == move_source_for_index:
            source = (source[0] + 4.0e-6, source[1], source[2])
        steps = GRID_STEPS
        if index == drop_grid_step_for_index:
            steps = GRID_STEPS[:-1]
        write_output(
            root / str(index) / "output.h5",
            cfg,
            sha,
            producer=run["arm"] == "producer",
            source_m=source,
            exact_timing=exact_timing,
            grid_steps=steps,
        )


def run_analyzer(work: Path, results: Path) -> tuple[int, dict]:
    proc = subprocess.run(
        [
            sys.executable,
            str(ASSAY / "analyze_assay.py"),
            "--results-root",
            str(results),
            "--output-dir",
            "analysis",
        ],
        cwd=work,
        capture_output=True,
        text=True,
        check=False,
    )
    gate_path = work / "analysis" / "assay_gate.json"
    assert gate_path.is_file(), proc.stdout + proc.stderr
    return proc.returncode, json.loads(gate_path.read_text())


@pytest.fixture(scope="module")
def runs(analyze_assay):
    manifest_runs, sha = analyze_assay.load_runs(Path("/nonexistent"))
    return manifest_runs, sha


def test_generated_assay_package_shape(runs, contract):
    manifest_runs, _sha = runs
    assert len(manifest_runs) == 12
    manifest = json.loads((GENERATED / "manifest.json").read_text())
    assert manifest["mpi_ranks"] == 1
    producers = [r for r in manifest_runs if r["arm"] == "producer"]
    nulls = [r for r in manifest_runs if r["arm"] == "toxin_free_null"]
    assert len(producers) == 9
    assert len(nulls) == 3
    amplitudes = contract["design"]["axes"]["bacteriocin.mucin_charge.amplitude"]
    assert sorted({r["amplitude"] for r in producers}) == list(amplitudes)
    assert sorted({r["seed"] for r in producers}) == sorted(contract["design"]["seeds"])
    for run in manifest_runs:
        cfg = json.loads((ASSAY / run["input_relpath"]).read_text())
        assert cfg["seed"] == run["seed"]
        assert cfg["bacteriocin.mucin_charge.amplitude"] == pytest.approx(run["amplitude"])
        assert cfg["hdf5"]["schedule"]["provenance"] == 1
        assert cfg["hdf5"]["schedule"]["grid_species"] == ["bacteriocin_BtuB"]


def test_deployment_generation_refuses_mismatched_sha():
    proc = subprocess.run(
        [
            sys.executable,
            str(ASSAY / "prepare_assay.py"),
            "--deployment",
            "--execution-source-sha",
            "0" * 40,
            "--image-digest",
            "repo@sha256:" + "1" * 64,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    assert "REFUSED" in proc.stdout + proc.stderr


def test_analyzer_passes_on_ordered_synthetic_arms(tmp_path, runs, contract):
    manifest_runs, sha = runs
    results = tmp_path / "results"
    build_results(results, manifest_runs, sha)
    code, gate = run_analyzer(tmp_path, results)
    assert gate["blockers"] == []
    assert gate["status"] == "PASS"
    assert code == 0
    assert gate[contract["gate"]["id"]] is True
    assert gate["C_transport_gate"] is False
    assert gate["checks"]["n_paired_snapshot_times"] == len(GRID_STEPS) * len(
        contract["design"]["seeds"]
    )
    assert gate["checks"]["nulls_toxin_free"] is True


def test_analyzer_fails_when_amplitudes_share_one_profile(tmp_path, runs, monkeypatch):
    manifest_runs, sha = runs
    monkeypatch.setitem(DECAY_UM, 15.0, DECAY_UM[0.0])
    results = tmp_path / "results"
    build_results(results, manifest_runs, sha)
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "FAIL"
    assert code == 1
    assert gate["checks"]["intrinsic_radii_ordered_0_gt_15_gt_60"] is False
    assert gate["blockers"] == []


def test_analyzer_blocks_on_missing_output(tmp_path, runs):
    manifest_runs, sha = runs
    results = tmp_path / "results"
    build_results(results, manifest_runs, sha, skip_index=3)
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED"
    assert code == 2
    assert gate["checks"]["all_outputs_returned"] is False
    assert any("no returned output file" in b for b in gate["blockers"])


def test_analyzer_blocks_on_legacy_provenance_without_exact_timing(tmp_path, runs):
    manifest_runs, sha = runs
    results = tmp_path / "results"
    build_results(results, manifest_runs, sha, exact_timing=False)
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED"
    assert code == 2
    assert any("event_time_s" in b for b in gate["blockers"])


def test_analyzer_blocks_on_incomplete_paired_times(tmp_path, runs):
    manifest_runs, sha = runs
    results = tmp_path / "results"
    build_results(results, manifest_runs, sha, drop_grid_step_for_index=3)
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED"
    assert code == 2
    assert gate["checks"]["paired_times_complete"] is False
    assert any("not shared by all amplitudes" in b for b in gate["blockers"])


def test_analyzer_blocks_when_source_moves_between_arms(tmp_path, runs):
    manifest_runs, sha = runs
    results = tmp_path / "results"
    build_results(results, manifest_runs, sha, move_source_for_index=3)
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED"
    assert code == 2
    assert gate["checks"]["source_identical_across_amplitudes"] is False


def test_analyzer_blocks_on_foreign_execution_sha(tmp_path, runs):
    manifest_runs, _sha = runs
    results = tmp_path / "results"
    build_results(results, manifest_runs, "f" * 40)
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED"
    assert code == 2
    assert gate["checks"]["all_runs_authenticated"] is False
