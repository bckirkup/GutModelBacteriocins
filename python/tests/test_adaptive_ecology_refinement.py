"""End-to-end tests for experiments/adaptive_ecology_refinement_v1.

The refinement generator, preflight, command printer, and analyzer live outside
the ``gut_ibm_tools`` package, so they are exercised by subprocess (or, where
git state must be stubbed, by module import).  Synthetic outputs are built with
h5py; the tests check that the analyzer blocks — rather than adapting — on
missing outputs, wrong provenance, a missing intrinsic gate, and a missing or
mismatched approval, and that the predeclared selection rule is applied exactly
as written.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
PKG = REPO / "experiments" / "adaptive_ecology_refinement_v1"
GENERATED = PKG / "generated"

SEEDS = [20260911, 20260913, 20260917]
AMPLITUDES = [15, 20, 30]
N_TYPE1 = 2000
T_END_S = 21600.0
BIO_DT = 60.0
# Agent snapshots inside and before the 7200 s analysis window.
AGENT_STEPS = (60, 240, 300, 360)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def run_tool(script: str, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(PKG / script), *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def manifest() -> dict:
    return json.loads((GENERATED / "manifest.json").read_text())


def resolved_config(entry: dict, cfg: dict) -> dict:
    schedule = cfg["hdf5"]["schedule"]
    resolved = {
        "seed": int(entry["seed"]),
        "bacteriocin.mucin_charge.amplitude": float(entry["amplitude"]),
        "kd_b12_btuB": 1.0e-4,
        "kd_colicinE_btuB": 5.0e-7,
        "b12.initial_conc": 1.0e-3,
        "burst_release_tau": 300.0,
        "hdf5.schedule.grid_species": list(schedule["grid_species"]),
    }
    for key, value in schedule.items():
        if key != "grid_species":
            resolved[f"hdf5.schedule.{key}"] = int(value)
    return resolved


def type_array(n1: int, n2: int) -> np.ndarray:
    return np.array([1] * max(n1, 0) + [2] * max(n2, 0), dtype=np.int64)


def write_output(
    path: Path,
    entry: dict,
    sha: str,
    *,
    slope_per_h: float = 0.0,
    final_n2: int = 60,
    kills: float = 0.0,
    lysis: float = 1.0,
    placement: str = "device_delivery",
    mpi_ranks: int = 1,
    termination: str = "horizon_reached",
) -> None:
    """Write one synthetic output whose window slope and final n_type2 are exact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((PKG / entry["input_relpath"]).read_text())
    ratio_final = math.log10((N_TYPE1 + 0.5) / (final_n2 + 0.5))
    with h5py.File(path, "w") as h5:
        prov = h5.create_group("run_provenance")
        prov.create_dataset("git_sha", data=np.bytes_(sha))
        prov.create_dataset("chemistry_placement", data=np.bytes_(placement))
        prov.create_dataset("mpi_rank_count", data=np.int64(mpi_ranks))
        prov.create_dataset("termination_cause", data=np.bytes_(termination))
        prov.create_dataset(
            "resolved_config", data=np.bytes_(json.dumps(resolved_config(entry, cfg)))
        )
        for step in AGENT_STEPS:
            key = f"step_{step:06d}"
            t = step * BIO_DT
            h5.create_dataset(f"summary/{key}/time", data=np.float64(t))
            ratio = ratio_final + slope_per_h * (t - T_END_S) / 3600.0
            n2 = max(1, round((N_TYPE1 + 0.5) * 10.0 ** (-ratio) - 0.5))
            h5.create_dataset(f"agents/{key}/type", data=type_array(N_TYPE1, n2))
        events = h5.create_group(f"summary/step_{AGENT_STEPS[-1]:06d}/events")
        events.create_dataset(
            "cumulative_mortality_colicin", data=np.float64(kills * lysis)
        )
        events.create_dataset("cumulative_mortality_lysis", data=np.float64(lysis))
        events.create_dataset(
            "cumulative_divisions_by_type", data=np.array([0, 500], dtype=np.int64)
        )


def build_results(
    root: Path,
    man: dict,
    *,
    plan: dict[float, dict] | None = None,
    skip: set[int] | None = None,
    sha: str | None = None,
    placement: str = "device_delivery",
) -> None:
    """Write the 12 outputs. ``plan`` maps amplitude -> {slope, final_n2, kills}."""
    skip = skip or set()
    exec_sha = sha or man["execution_source_sha"]
    for entry in man["runs"]:
        index = int(entry["array_index"])
        if index in skip:
            continue
        if entry["arm"] == "producer":
            spec = (plan or {}).get(float(entry["amplitude"]), {})
            slope = spec.get("slope", 0.3)
            final_n2 = spec.get("final_n2", 1200)
            kills = spec.get("kills", 10.0)
        else:
            slope, final_n2, kills = 0.0, 60, 0.0
        write_output(
            root / str(index) / "output.h5",
            entry,
            exec_sha,
            slope_per_h=slope,
            final_n2=final_n2,
            kills=kills,
            lysis=10.0,
            placement=placement,
        )


def run_analyzer(work: Path, results: Path, *extra: str) -> tuple[int, dict]:
    proc = subprocess.run(
        [
            sys.executable,
            str(PKG / "analyze_refinement.py"),
            "--results-root",
            str(results),
            "--results-dir",
            "analysis",
            *extra,
        ],
        cwd=work,
        capture_output=True,
        text=True,
        check=False,
    )
    gate_path = work / "analysis" / "refinement_gate.json"
    assert gate_path.is_file(), proc.stdout + proc.stderr
    return proc.returncode, json.loads(gate_path.read_text())


def intrinsic_gate_file(work: Path, man: dict, *, mismatch: str | None = None) -> Path:
    data = {
        "single_source_transport_gate": True,
        "execution_source_sha": man["execution_source_sha"],
        "container_image_digest": man["container_image_digest"],
    }
    if mismatch == "sha":
        data["execution_source_sha"] = "0" * 40
    if mismatch == "digest":
        data["container_image_digest"] = "other@sha256:" + "2" * 64
    path = work / "intrinsic_gate.json"
    path.write_text(json.dumps(data))
    return path


def approval_file(work: Path, man: dict, amplitude: float, **overrides) -> Path:
    data = {
        "authorize_c_transport_gate": True,
        "execution_source_sha": man["execution_source_sha"],
        "container_image_digest": man["container_image_digest"],
        "jobs": 12,
        "selected_mucin_charge_amplitude": amplitude,
    }
    data.update(overrides)
    path = work / "approval.json"
    path.write_text(json.dumps(data))
    return path


def qualifying_plan() -> dict[float, dict]:
    """Amplitude 20 qualifies on every criterion; 15 and 30 do not."""
    return {
        15.0: {"slope": 0.45, "final_n2": 1200, "kills": 9.0},  # slope window miss
        20.0: {"slope": 0.30, "final_n2": 1200, "kills": 10.0},
        30.0: {"slope": 0.35, "final_n2": 1500, "kills": 8.0},
    }


def test_expected_run_defaults_match_contract_identity():
    sys.path.insert(0, str(REPO / "python"))
    from gut_ibm_tools.transport_metrics import ExpectedRun

    expected = ExpectedRun(
        execution_source_sha="x", seed=1, amplitude=15.0,
        chemistry_placement="device_delivery",
    )
    assert expected.kd_corrinoid_btuB == pytest.approx(1.0e-4)
    assert expected.kd_colicinE_btuB == pytest.approx(5.0e-7)
    assert expected.b12_initial_conc == pytest.approx(1.0e-3)
    assert expected.burst_release_tau_s == pytest.approx(300.0)
    assert expected.grid_species == ("bacteriocin_BtuB",)


def test_planning_generation_shape_and_order():
    proc = run_tool("prepare_refinement.py", "--clean")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    man = manifest()
    assert man["schema_version"] == 2
    assert man["refinement_id"] == "adaptive_ecology_refinement_v1"
    assert man["mpi_ranks"] == 1
    assert man["gpu_required"] is True
    assert man["attempt_timeout_s"] == 7200
    assert man["chemistry_placement"] == "device_delivery"
    assert man["gate"] == "C_transport_gate"
    expected_ids = [
        f"C_amp{amp}_producer_s{seed}" for amp in AMPLITUDES for seed in SEEDS
    ] + [f"C_amp15_plasmid_free_null_s{seed}" for seed in SEEDS]
    assert [r["run_id"] for r in man["runs"]] == expected_ids
    for index, run in enumerate(man["runs"]):
        assert run["array_index"] == index
        assert run["input_relpath"] == f"generated/jobs/{index}/input.json"
        assert run["output_relpath"] == f"generated/results/{index}/output.h5.gz"
        path = GENERATED / "jobs" / str(index) / "input.json"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == run["input_sha256"]
        cfg = json.loads(path.read_text())
        assert cfg["metabolism.uptake_limit"] == "delivery"
        assert cfg["fixes"] == ["metabolism", "bacteriocin", "receptor", "mechanics"]
        plasmids = [s["plasmids"] for s in cfg["initial_strains"]]
        if run["arm"] == "producer":
            assert plasmids == [["ColE1"], []]
        else:
            assert plasmids == [[], []]


def test_real_identity_without_deployment_is_refused():
    proc = run_tool(
        "prepare_refinement.py", "--execution-source-sha", "0" * 40
    )
    assert proc.returncode != 0
    assert "REFUSED" in proc.stdout + proc.stderr


def test_deployment_generation_refuses_mismatched_sha():
    proc = run_tool(
        "prepare_refinement.py",
        "--deployment",
        "--execution-source-sha", "0" * 40,
        "--image-digest", "repo@sha256:" + "1" * 64,
    )
    assert proc.returncode != 0
    assert "REFUSED" in proc.stdout + proc.stderr


def test_planning_preflight_passes_with_explicit_warnings():
    proc = run_tool("preflight_refinement.py")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    result = json.loads(proc.stdout)
    assert result["status"] == "PASS"
    assert result["mode"] == "planning"
    assert result["runs_checked"] == 12
    assert result["warnings"]


def mutate_package(tmp_path: Path, mutate) -> Path:
    """Copy the package to a scratch dir and apply a mutation to the copy."""
    dest = tmp_path / "pkg"
    dest.mkdir()
    for name in (
        "refinement_contract.json",
        "refinement_decision_record.json",
        "prepare_refinement.py",
        "preflight_refinement.py",
        "aws_commands_refinement.py",
        "analyze_refinement.py",
        "README.md",
    ):
        shutil.copy2(PKG / name, dest / name)
    shutil.copytree(GENERATED, dest / "generated")
    mutate(dest)
    return dest


def rehash_manifest_inputs(pkg: Path) -> None:
    man_path = pkg / "generated" / "manifest.json"
    man = json.loads(man_path.read_text())
    for run in man["runs"]:
        run["input_sha256"] = hashlib.sha256(
            (pkg / run["input_relpath"]).read_bytes()
        ).hexdigest()
    man_path.write_text(json.dumps(man, indent=2))


def mutate_input(pkg: Path, index: int, edit) -> None:
    path = pkg / "generated" / "jobs" / str(index) / "input.json"
    cfg = json.loads(path.read_text())
    edit(cfg)
    path.write_text(json.dumps(cfg, indent=2))
    rehash_manifest_inputs(pkg)


@pytest.mark.parametrize(
    "mutate, needle",
    [
        (
            lambda pkg: mutate_input(pkg, 0, lambda c: c.__setitem__("bacteriocin.mucin_charge.amplitude", 21)),
            "amplitude",
        ),
        (
            lambda pkg: mutate_input(pkg, 0, lambda c: c.__setitem__("kd_corrinoid_btuB", 2e-4)),
            "kd_corrinoid_btuB",
        ),
        (
            lambda pkg: mutate_input(pkg, 0, lambda c: c.__setitem__("b12_initial_conc", 2e-3)),
            "b12_initial_conc",
        ),
        (
            lambda pkg: mutate_input(pkg, 0, lambda c: c["hdf5"]["schedule"].__setitem__("agents", 5)),
            "HDF5 schedule drift",
        ),
        (
            lambda pkg: mutate_input(pkg, 9, lambda c: c["initial_strains"][0].__setitem__("plasmids", ["ColE1"])),
            "plasmid identity drift",
        ),
        (
            lambda pkg: mutate_manifest(pkg, lambda m: m.__setitem__("chemistry_placement", "host_forced_delivery")),
            "device_delivery",
        ),
        (
            lambda pkg: mutate_manifest(pkg, lambda m: m["runs"][0].__setitem__("run_id", "SS_amp15_producer_s20260911")),
            "only ecological C_amp",
        ),
        (
            lambda pkg: mutate_manifest(pkg, lambda m: m["runs"].pop()),
            "expected 12",
        ),
        (
            lambda pkg: mutate_manifest(pkg, lambda m: m["runs"][0].__setitem__("input_sha256", "0" * 64)),
            "input SHA-256 mismatch",
        ),
    ],
)
def test_preflight_fails_closed_on_mutated_tree(tmp_path, mutate, needle):
    pkg = mutate_package(tmp_path, mutate)
    proc = subprocess.run(
        [sys.executable, str(pkg / "preflight_refinement.py")],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 1
    result = json.loads(proc.stdout)
    assert result["status"] == "FAIL"
    assert any(needle in error for error in result["errors"]), result["errors"]


def mutate_manifest(pkg: Path, edit) -> None:
    man_path = pkg / "generated" / "manifest.json"
    man = json.loads(man_path.read_text())
    edit(man)
    man_path.write_text(json.dumps(man, indent=2))


def test_preflight_fails_on_device_placement(tmp_path):
    pkg = mutate_package(
        tmp_path,
        lambda pkg: mutate_manifest(
            pkg, lambda m: m.__setitem__("chemistry_placement", "device")
        ),
    )
    proc = subprocess.run(
        [sys.executable, str(pkg / "preflight_refinement.py")],
        capture_output=True, text=True, check=False,
    )
    result = json.loads(proc.stdout)
    assert result["status"] == "FAIL"
    assert any("device_delivery" in error for error in result["errors"])


@pytest.fixture()
def aws():
    return load_module(
        "refinement_aws_under_test", PKG / "aws_commands_refinement.py"
    )


def aws_argv(man: dict, **overrides) -> list[str]:
    args = {
        "--execution-source-sha": "a" * 40,
        "--image-uri": "repo@sha256:" + "1" * 64,
        "--input-prefix": "s3://bucket/adaptive-ecology-refinement-v1/inputs",
        "--output-prefix": "s3://bucket/adaptive-ecology-refinement-v1/outputs",
        "--job-queue": "queue",
        "--job-definition": "definition",
    }
    args.update(overrides)
    return [item for pair in args.items() for item in pair]


def stub_deployment(aws, monkeypatch, man: dict | None, head: str = "a" * 40):
    monkeypatch.setattr(aws, "git", lambda *a: head)
    if man is not None:
        monkeypatch.setattr(aws, "load_manifest", lambda: man)
    monkeypatch.setattr(
        aws.subprocess,
        "run",
        lambda *a, **k: SimpleNamespace(returncode=0, stdout="", stderr=""),
    )


@pytest.mark.parametrize(
    "overrides, needle",
    [
        ({"--execution-source-sha": "deadbeef"}, "full lowercase 40-hex"),
        ({"--image-uri": "mutable:latest"}, "image URI must be immutable"),
        (
            {"--input-prefix": "s3://bucket/single-source-assay/inputs"},
            "adaptive-ecology-refinement",
        ),
        (
            {"--input-prefix": "s3://bucket/adaptive-ecology-refinement/stage_c/inputs"},
            "stale",
        ),
        (
            {"--input-prefix": "s3://bucket/adaptive-ecology-refinement-v1/outputs",
             "--output-prefix": "s3://bucket/adaptive-ecology-refinement-v1/outputs"},
            "isolated",
        ),
        ({"--job-name": "gutibm-single-source-assay"}, "job name"),
    ],
)
def test_command_generator_refuses_bad_requests(aws, monkeypatch, capsys, overrides, needle):
    stub_deployment(aws, monkeypatch, manifest())
    monkeypatch.setattr(sys, "argv", ["aws_commands_refinement.py"] + aws_argv(manifest(), **overrides))
    with pytest.raises(SystemExit) as excinfo:
        aws.main()
    assert needle in str(excinfo.value)


def test_command_generator_refuses_sha_not_at_head(aws, monkeypatch):
    stub_deployment(aws, monkeypatch, manifest(), head="b" * 40)
    monkeypatch.setattr(
        sys, "argv", ["aws_commands_refinement.py"] + aws_argv(manifest())
    )
    with pytest.raises(SystemExit) as excinfo:
        aws.main()
    assert "differs from git HEAD" in str(excinfo.value)


def test_command_generator_refuses_mismatched_manifest_identity(aws, monkeypatch):
    man = manifest()
    man["execution_source_sha"] = "a" * 40
    stub_deployment(aws, monkeypatch, man)
    monkeypatch.setattr(
        sys, "argv", ["aws_commands_refinement.py"] + aws_argv(man)
    )
    with pytest.raises(SystemExit) as excinfo:
        aws.main()
    assert "not byte-for-byte equal" in str(excinfo.value)


def test_command_generator_prints_but_never_executes(aws, monkeypatch, capsys):
    man = manifest()
    man["execution_source_sha"] = "a" * 40
    man["container_image_digest"] = "repo@sha256:" + "1" * 64
    stub_deployment(aws, monkeypatch, man)
    monkeypatch.setattr(
        sys, "argv", ["aws_commands_refinement.py"] + aws_argv(man)
    )
    assert aws.main() == 0
    out = capsys.readouterr().out
    assert "# REVIEW ONLY. This program did not execute AWS commands." in out
    assert "aws s3 cp" in out and "--include '*/input.json'" in out
    assert "aws batch submit-job" in out and "size=12" in out
    assert "attemptDurationSeconds=7200" in out


def test_analyzer_ready_for_approval_when_amplitude_qualifies(tmp_path):
    man = manifest()
    results = tmp_path / "results"
    build_results(results, man, plan=qualifying_plan())
    gate_file = intrinsic_gate_file(tmp_path, man)
    code, gate = run_analyzer(tmp_path, results, "--intrinsic-gate-file", str(gate_file))
    assert gate["status"] == "READY_FOR_APPROVAL"
    assert code == 1
    assert gate["selected_amplitude"] == pytest.approx(20.0)
    assert gate["selection_basis"] == "qualified"
    assert gate["C_transport_gate"] is False
    assert gate["stage_d_release"] == "BLOCKED"
    selection = json.loads((tmp_path / "analysis" / "amplitude_selection.json").read_text())
    by_amp = {r["amplitude"]: r for r in selection["per_amplitude"]}
    assert by_amp[15.0]["qualifying"] is False
    assert by_amp[20.0]["qualifying"] is True
    assert by_amp[30.0]["qualifying"] is True
    # Graded response: the kills knob orders the medians exactly as injected.
    assert by_amp[20.0]["median_kills_per_lysis"] > by_amp[30.0]["median_kills_per_lysis"]
    # Integer agent counts quantize the series; the recovered slope is within
    # ~1e-3 of the injected 0.30 log10/h.
    assert by_amp[20.0]["median_delta_slope_log10_ratio_per_h"] == pytest.approx(0.30, abs=1e-3)


def test_selection_tie_rule_prefers_lower_amplitude(tmp_path):
    man = manifest()
    plan = {
        15.0: {"slope": 0.30, "final_n2": 1200, "kills": 9.5},
        20.0: {"slope": 0.30, "final_n2": 1300, "kills": 10.0},
        30.0: {"slope": 0.30, "final_n2": 1400, "kills": 7.0},
    }
    results = tmp_path / "results"
    build_results(results, man, plan=plan)
    run_analyzer(tmp_path, results)
    selection = json.loads((tmp_path / "analysis" / "amplitude_selection.json").read_text())
    assert selection["basis"] == "qualified"
    # 15 is within 10% of the maximum median kills (10.0) -> lowest wins.
    assert selection["selected_amplitude"] == pytest.approx(15.0)
    assert selection["amplitudes_within_tolerance"] == [15.0, 20.0]


def test_analyzer_blocks_on_missing_output(tmp_path):
    man = manifest()
    results = tmp_path / "results"
    build_results(results, man, plan=qualifying_plan(), skip={4})
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED_MISSING_OUTPUTS"
    assert code == 2
    assert gate["C_transport_gate"] is False
    missing = json.loads((tmp_path / "analysis" / "missing_outputs.json").read_text())
    assert [r["run_id"] for r in missing] == [man["runs"][4]["run_id"]]


def test_analyzer_blocks_on_foreign_execution_sha(tmp_path):
    man = manifest()
    results = tmp_path / "results"
    build_results(results, man, plan=qualifying_plan(), sha="f" * 40)
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED_AUTHENTICATION"
    assert code == 2
    assert gate["authentication_problems"]


def test_analyzer_blocks_on_host_forced_delivery(tmp_path):
    man = manifest()
    results = tmp_path / "results"
    build_results(
        results, man, plan=qualifying_plan(), placement="host_forced_delivery"
    )
    code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED_AUTHENTICATION"
    assert code == 2


def test_analyzer_records_fallback_when_no_amplitude_qualifies(tmp_path):
    man = manifest()
    plan = {
        15.0: {"slope": 0.45, "final_n2": 800, "kills": 9.0},   # positive but out of window; sus >= 500
        20.0: {"slope": -0.10, "final_n2": 1200, "kills": 10.0},  # negative slope
        30.0: {"slope": 0.50, "final_n2": 300, "kills": 8.0},
    }
    results = tmp_path / "results"
    build_results(results, man, plan=plan)
    gate_file = intrinsic_gate_file(tmp_path, man)
    code, gate = run_analyzer(tmp_path, results, "--intrinsic-gate-file", str(gate_file))
    assert gate["status"] == "BLOCKED_NO_QUALIFYING_AMPLITUDE"
    assert code == 2
    assert gate["selection_basis"] == "fallback_reproduction_of_prior_effect"
    assert gate["selected_amplitude"] == pytest.approx(15.0)
    assert gate["C_transport_gate"] is False


def test_analyzer_blocks_without_intrinsic_gate(tmp_path):
    man = manifest()
    results = tmp_path / "results"
    build_results(results, man, plan=qualifying_plan())
    _code, gate = run_analyzer(tmp_path, results)
    assert gate["status"] == "BLOCKED_MISSING_INTRINSIC_GATE"
    assert gate["single_source_transport_gate"] is None
    _code, gate = run_analyzer(
        tmp_path, results,
        "--intrinsic-gate-file", str(intrinsic_gate_file(tmp_path, man, mismatch="sha")),
    )
    assert gate["status"] == "BLOCKED_MISSING_INTRINSIC_GATE"
    assert gate["single_source_transport_gate"] is False


def test_gate_passes_only_with_matching_approval(tmp_path):
    man = manifest()
    results = tmp_path / "results"
    build_results(results, man, plan=qualifying_plan())
    gate_file = intrinsic_gate_file(tmp_path, man)
    good = approval_file(tmp_path, man, 20.0)
    code, gate = run_analyzer(
        tmp_path, results,
        "--intrinsic-gate-file", str(gate_file),
        "--approval-file", str(good),
    )
    assert gate["status"] == "PASS"
    assert code == 0
    assert gate["C_transport_gate"] is True
    assert gate["stage_d_release"] == "RELEASED"

    bad_amp = tmp_path / "approval_bad_amp.json"
    bad_amp.write_text(json.dumps({
        "authorize_c_transport_gate": True,
        "execution_source_sha": man["execution_source_sha"],
        "container_image_digest": man["container_image_digest"],
        "jobs": 12,
        "selected_mucin_charge_amplitude": 15.0,
    }))
    code, gate = run_analyzer(
        tmp_path, results,
        "--intrinsic-gate-file", str(gate_file),
        "--approval-file", str(bad_amp),
    )
    assert gate["status"] != "PASS"
    assert gate["C_transport_gate"] is False
    assert gate["stage_d_release"] == "BLOCKED"
    assert any("amplitude" in p for p in gate["approval_problems"])

    bad_sha = tmp_path / "approval_bad_sha.json"
    bad_sha.write_text(json.dumps({
        "authorize_c_transport_gate": True,
        "execution_source_sha": "0" * 40,
        "container_image_digest": man["container_image_digest"],
        "jobs": 12,
        "selected_mucin_charge_amplitude": 20.0,
    }))
    code, gate = run_analyzer(
        tmp_path, results,
        "--intrinsic-gate-file", str(gate_file),
        "--approval-file", str(bad_sha),
    )
    assert gate["status"] != "PASS"
    assert any("execution_source_sha" in p for p in gate["approval_problems"])
