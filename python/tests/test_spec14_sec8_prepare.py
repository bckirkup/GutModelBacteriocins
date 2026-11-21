"""Unit tests for experiments/spec14_sec8/prepare.py (#464)."""
from __future__ import annotations

import importlib.util
import json
import pathlib

import pytest


def _load_prepare_module():
    repo_root = pathlib.Path(__file__).resolve().parents[2]
    script_path = repo_root / "experiments" / "spec14_sec8" / "prepare.py"
    spec = importlib.util.spec_from_file_location("spec14_sec8_prepare", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_prepare_generates_carriage_matched_and_non_btub_controls(
    tmp_path: pathlib.Path,
) -> None:
    prepare = _load_prepare_module()
    expected_treatments = {
        "producer",
        "null",
        "kill_null",
        "immunity_null",
        "colIa_producer",
        "colIa_kill_null",
        "colB_producer",
        "colB_kill_null",
    }
    assert expected_treatments.issubset(set(prepare.TREATMENTS.keys()))

    prepare.main(["--output-dir", str(tmp_path), "--conditions", "vitro"])

    for treatment in expected_treatments:
        arm_file = tmp_path / f"vitro_{treatment}_s{prepare.SEEDS[0]}" / "input.json"
        assert arm_file.is_file()
        cfg = json.loads(arm_file.read_text())
        strains = cfg["initial_strains"]
        assert len(strains) == 2

        if treatment == "producer":
            assert strains[0]["plasmids"] == ["ColE2"]
            assert strains[1]["plasmids"] == []
        elif treatment == "null":
            assert strains[0]["plasmids"] == []
            assert strains[1]["plasmids"] == []
        elif treatment == "kill_null":
            assert strains[0]["plasmids"] == ["ColE2"]
            assert strains[1]["plasmids"] == []
            assert cfg["kill_rate_colicin"] == pytest.approx(0.0)
            assert cfg["kill_rate_microcin"] == pytest.approx(0.0)
            assert cfg["sos_cross_induction_rate"] == pytest.approx(0.0)
        elif treatment == "immunity_null":
            assert strains[0]["plasmids"] == ["ColE2"]
            assert strains[1]["plasmids"] == ["ColE2"]
        elif treatment.startswith("colIa"):
            assert strains[0]["plasmids"] == ["ColIa"]
            assert "ColE1" not in strains[0]["plasmids"]
        elif treatment.startswith("colB"):
            assert strains[0]["plasmids"] == ["ColB"]
            assert "ColE1" not in strains[0]["plasmids"]


def test_prepare_decoupling_cli_flags(tmp_path: pathlib.Path) -> None:
    prepare = _load_prepare_module()
    prepare.main([
        "--output-dir", str(tmp_path),
        "--conditions", "vitro",
        "--treatments", "producer", "kill_null",
        "--no-b12-growth-essential",
        "--btub-competitor-conc", "1e-6",
        "--b12-free-fraction", "0.05",
    ])
    arm_file = tmp_path / f"vitro_producer_s{prepare.SEEDS[0]}" / "input.json"
    cfg = json.loads(arm_file.read_text())
    assert cfg["metabolism.b12_growth_essential"] is False
    assert cfg["receptor.btub_competitor_conc"] == pytest.approx(1.0e-6)
    assert cfg["receptor.b12_free_fraction"] == pytest.approx(0.05)
