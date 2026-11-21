"""Tests for Spec 13 gradient post-purge readout (#467) and oxygen BC presets (#470)."""

from __future__ import annotations

import importlib.util
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def _load_module(name: str, rel_path: str):
    path = REPO / rel_path
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _extract_known_parser_keys() -> set[str]:
    """Extract recognized scalar and nested config keys from C++ parser sources."""
    keys: set[str] = {"_comment", "hdf5", "layer2", "layer3", "initial_strains"}
    for rel in ("src/io/input_parser.cpp", "src/io/config_json.cpp"):
        text = (REPO / rel).read_text()
        for match in re.finditer(r'"([A-Za-z0-9_.]+)"', text):
            keys.add(match.group(1))
    return keys


def test_spec13_gradient_strict_post_purge_and_72h_endpoint():
    readout = _load_module(
        "spec13_run_readout", "analysis/spec13_gradient/run_readout.py"
    )
    ts_path = REPO / "analysis/spec13_gradient/out/gradient_seed1001_timeseries.csv"
    metrics = readout.tail_metrics(ts_path)

    post = metrics["post_purge_window_gt_48h"]
    assert post["sample_count_per_region"] == 143
    assert math.isclose(
        post["time_range_s"][0], 173400.0, rel_tol=0.0, abs_tol=1e-9
    )
    assert math.isclose(
        post["time_range_s"][1], 258600.0, rel_tol=0.0, abs_tol=1e-9
    )
    assert 1.0e6 <= post["stool_cfu_g"] <= 1.0e8
    assert post["flatness"] <= 4.0

    # Post-purge window mean (~0.9945) differs from 71h50m (0.9025) and 72h prov (0.8975).
    assert post["occupancy"]["cecum"] > 0.99
    assert math.isclose(
        metrics["last_csv_sample"]["cecum"]["occ_frac"],
        0.9025,
        rel_tol=0.0,
        abs_tol=1e-9,
    )
    ep_cecum = metrics["endpoint_72h_prov"]["regions"]["cecum"]
    assert math.isclose(
        ep_cecum["occupancy_frac"], 0.8975, rel_tol=0.0, abs_tol=1e-9
    )
    assert ep_cecum["mean_cells_per_occupied_patch"] > 15.0


def test_spec13_occupancy_sensitivity_generator(tmp_path: Path):
    prep = _load_module(
        "spec13_occ_sens_prepare",
        "experiments/spec13_occupancy_sensitivity/prepare.py",
    )
    written = prep.generate_all(output_dir=tmp_path, seeds=(1001, 1002))
    assert len(written) == 12

    seen_p0: set[float] = set()
    for cfg_path in written:
        data = json.loads(cfg_path.read_text())
        assert data["layer3"]["enabled"] is True
        assert data["layer2"]["n_patches"] == data["layer3"]["n_patches_per_region"]
        seen_p0.add(data["layer3"]["reattach_p0"])
    assert len(seen_p0) == 3
    assert any(math.isclose(v, 0.005, rel_tol=0.0, abs_tol=1e-12) for v in seen_p0)
    assert any(math.isclose(v, 0.015, rel_tol=0.0, abs_tol=1e-12) for v in seen_p0)
    assert any(math.isclose(v, 0.05, rel_tol=0.0, abs_tol=1e-12) for v in seen_p0)


def test_oxygen_boundary_presets_match_parser_invariants():
    known_keys = _extract_known_parser_keys()
    preset_dir = REPO / "examples/oxygen_boundary_presets"
    preset_files = sorted(preset_dir.glob("*.json"))
    assert len(preset_files) == 4

    modes_seen: set[str] = set()
    for pfile in preset_files:
        cfg = json.loads(pfile.read_text())
        for key in cfg:
            assert key in known_keys, (
                f"Unrecognized top-level key {key} in {pfile.name}"
            )
        assert cfg["oxygen.enabled"] is True
        mode = cfg["oxygen.epithelial_boundary"]
        modes_seen.add(mode)
        if mode in ("robin", "flux"):
            assert cfg["oxygen.z_gradient"] is False
            if mode == "robin":
                assert cfg["oxygen.epithelial_transfer_coeff"] > 0.0
        if cfg.get("oxygen.respiration_driver") == "funded":
            assert cfg.get("uptake_limit") == "delivery"
            assert cfg.get("oxygen.delivery_uptake_enabled") is True
        maint_factor = cfg["oxygen.anaerobic_maintenance_factor"]
        assert 3.4 <= maint_factor <= 8.0

    assert modes_seen == {"dirichlet", "robin"}
