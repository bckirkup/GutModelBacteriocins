"""Tests for GutIBM checkpoint analysis CLI and helpers."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from gut_ibm_tools.checkpoint_analyze import main as analyze_main
from gut_ibm_tools.checkpoint_scan import discover_checkpoints
from gut_ibm_tools.checkpoint_summary import extract_checkpoint_row
from gut_ibm_tools.hdf5_reader import GutIBMData

from conftest import write_sample_hdf5

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")


def _install_checkpoints(dest: Path, n_files: int = 2) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for step_idx in range(n_files):
        # One closed restart file per step (single step group inside).
        path = dest / f"step_{step_idx:06d}.h5"
        write_sample_hdf5(path, n_agents=12, n_steps=1)
        # Overwrite summary step index/time to match filename.
        with h5py.File(path, "a") as handle:
            grp = handle["summary/step_000000"]
            del grp["step"]
            del grp["time"]
            grp.create_dataset("step", data=np.array(step_idx * 30, dtype=np.int32))
            grp.create_dataset("time", data=np.array(float(step_idx * 1800.0)))
            events = grp["events"]
            del events["colicin_kills"]
            events.create_dataset(
                "colicin_kills", data=np.array(step_idx * 5, dtype=np.int32)
            )
        paths.append(path)
    # Incomplete multipart junk must be ignored by discovery.
    (dest / "step_000099.h5.ABCDEF01").write_bytes(b"not-a-real-hdf5")
    return paths


def test_discover_checkpoints_sorts_and_ignores_junk(tmp_path: Path) -> None:
    ckpt = tmp_path / "ckpt"
    _install_checkpoints(ckpt, n_files=3)
    # Out-of-order create already sorted by name pattern; add high step early name.
    late = ckpt / "step_000100.h5"
    shutil.copy(ckpt / "step_000000.h5", late)

    found = discover_checkpoints(ckpt)
    names = [p.name for p in found]
    assert names == [
        "step_000000.h5",
        "step_000001.h5",
        "step_000002.h5",
        "step_000100.h5",
    ]
    assert not any(n.endswith(".ABCDEF01") for n in names)


def test_extract_checkpoint_row_golden(tmp_path: Path) -> None:
    path = tmp_path / "step_000000.h5"
    write_sample_hdf5(path, n_agents=12, n_steps=1)
    row = extract_checkpoint_row(path)

    assert row["n_alive"] == pytest.approx(12.0)
    assert row["n_type1"] == pytest.approx(6.0)
    assert row["n_type2"] == pytest.approx(6.0)
    assert row["event_colicin_kills"] == 0
    assert row["chem_mean_carbon"] == pytest.approx(1.0)
    assert "grid_mean_bacteriocin_BtuB" in row
    assert "hopkins" in row
    assert "monochromatic_score" in row
    assert row["nx"] == 4
    assert row["grid_dx"] == pytest.approx(5e-6)


def test_extract_mixing_sensitivity(tmp_path: Path) -> None:
    """Segregated clusters → higher monochromatic score than mixed layout."""
    segregated = tmp_path / "seg.h5"
    mixed = tmp_path / "mix.h5"
    write_sample_hdf5(segregated, n_agents=12, n_steps=1)

    # Rebuild mixed positions in-place: interleave types spatially.
    write_sample_hdf5(mixed, n_agents=12, n_steps=1)
    with h5py.File(mixed, "a") as handle:
        agents = handle["agents/step_000000"]
        rng = np.random.default_rng(7)
        for axis in ("x", "y", "z"):
            del agents[axis]
            agents.create_dataset(axis, data=rng.uniform(0, 100e-6, 12))

    seg_row = extract_checkpoint_row(segregated)
    mix_row = extract_checkpoint_row(mixed)
    assert seg_row["monochromatic_score"] > mix_row["monochromatic_score"]


def test_summarize_cli_writes_csv_and_meta(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    ckpt = tmp_path / "ckpt"
    _install_checkpoints(ckpt, n_files=2)

    rc = analyze_main(
        [
            "--summarize",
            "--checkpoint-dir",
            str(ckpt),
            "--output",
            "summary.csv",
            "--quiet",
        ]
    )
    assert rc == 0
    frame = pd.read_csv(tmp_path / "summary.csv")
    assert len(frame) == 2
    assert list(frame["step"]) == [0, 30]
    assert "n_alive" in frame.columns
    assert "event_colicin_kills" in frame.columns
    meta = json.loads((tmp_path / "summary_meta.json").read_text(encoding="utf-8"))
    assert meta["n_written"] == 2
    assert meta["n_skipped"] == 0


def test_summarize_skips_corrupt_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    ckpt = tmp_path / "ckpt"
    _install_checkpoints(ckpt, n_files=1)
    bad = ckpt / "step_000050.h5"
    bad.write_bytes(b"corrupt")

    rc = analyze_main(
        [
            "--summarize",
            "--checkpoint-dir",
            str(ckpt),
            "--output",
            "summary.csv",
            "--quiet",
        ]
    )
    assert rc == 0
    frame = pd.read_csv(tmp_path / "summary.csv")
    assert len(frame) == 1
    meta = json.loads((tmp_path / "summary_meta.json").read_text(encoding="utf-8"))
    assert meta["n_skipped"] == 1


def test_snapshot_and_timeseries_cli(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    ckpt = tmp_path / "ckpt"
    paths = _install_checkpoints(ckpt, n_files=2)

    rc_sum = analyze_main(
        [
            "--summarize",
            "--checkpoint-dir",
            str(ckpt),
            "--output",
            "summary.csv",
            "--quiet",
        ]
    )
    assert rc_sum == 0

    rc_snap = analyze_main(
        [
            "--snapshot",
            "--checkpoint",
            str(paths[0]),
            "--output",
            "snap0",
            "--quiet",
        ]
    )
    assert rc_snap == 0
    snap_dir = tmp_path / "snap0"
    assert (snap_dir / "agents_xz.png").is_file()
    assert (snap_dir / "density_vs_depth.png").is_file()
    assert (snap_dir / "heatmap_bacteriocin_BtuB.png").is_file()
    assert (snap_dir / "heatmap_carbon.png").is_file()

    rc_ts = analyze_main(
        [
            "--timeseries",
            "--summary",
            "summary.csv",
            "--output",
            "ts",
            "--quiet",
        ]
    )
    assert rc_ts == 0
    assert (tmp_path / "ts" / "population.png").is_file()


def test_get_grid_volumes_preserves_shape(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        volumes = data.get_grid_volumes("step_000000")
        assert volumes["carbon"].shape == (4, 5, 1)
        flat = data.get_grid("step_000000")
        assert flat["carbon"].shape == (20,)
        assert data.has_layer("summary")
        assert data.grid_shape == (4, 5, 1)
