"""Tests for the GPU-vs-CPU outcome comparator."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from gut_ibm_tools.gpu_precision import (
    compare,
    load_run,
    mann_whitney_exact_two_sided,
)

STEPS = [0, 1, 2, 3]
BASE_OXYGEN = [1.0, 1.0, 1.0, 1.0]
ZERO_CLIP = [0.0, 0.0, 0.0, 0.0]


def _write_run(
    path: Path,
    populations: list[int],
    *,
    clip: list[float] | None = None,
    oxygen: list[float] | None = None,
) -> None:
    clip_values = clip or ZERO_CLIP
    oxygen_values = oxygen or BASE_OXYGEN
    with h5py.File(path, "w") as handle:
        for index, step in enumerate(STEPS):
            group = handle.create_group(f"summary/step_{step:06d}")
            group.create_dataset("n_total", data=np.int32(populations[index]))
            group.create_dataset("num_lineages", data=np.int32(3))
            group.create_dataset(
                "mean_realized_fermentation_fraction", data=np.float64(0.5)
            )
            group.create_dataset("halt_reason_code", data=np.int32(0))
            events = group.create_group("events")
            for name in (
                "divisions",
                "sos_inductions",
                "mortality_colicin",
                "mortality_lysis",
                "mortality_cdi",
                "conjugation_transfers",
                "outflow_washout",
                "outflow_boundary",
            ):
                events.create_dataset(f"cumulative_{name}", data=np.int32(index))
            chemistry = group.create_group("chem")
            chemistry.create_dataset("mean_carbon", data=np.float64(1.0))
            chemistry.create_dataset("mean_iron", data=np.float64(1.0))
            chemistry.create_dataset(
                "mean_oxygen", data=np.float64(oxygen_values[index])
            )
            for receptor in ("BtuB", "FepA", "CirA", "FhuA"):
                chemistry.create_dataset(
                    f"max_toxin_{receptor}", data=np.float64(1.0e-6)
                )
            stocks = group.create_group("stocks")
            stocks.create_dataset("bacteriostatic_live_agents", data=np.int32(1))
            stocks.create_dataset(
                "washout_trapped_live_agents", data=np.int32(0)
            )
            nutrient_flux = group.create_group("nutrient_flux")
            nutrient_flux.create_dataset(
                "reaction_clip_cumulative",
                data=np.array([clip_values[index], 0.0]),
            )
            nutrient_flux.create_dataset(
                "uptake_shortfall_cumulative", data=np.array([0.0, 0.0])
            )
            nutrient_flux.create_dataset(
                "maintenance_shortfall_cumulative", data=np.array([0.0, 0.0])
            )


def _load_pair(
    tmp_path: Path,
    tag: str,
    host_values: dict[int, list[int]],
    device_values: dict[int, list[int]],
    *,
    device_clip: list[float] | None = None,
    device_oxygen: list[float] | None = None,
    repeat_values: list[int] | None = None,
):
    host_runs = []
    device_runs = []
    for seed, values in host_values.items():
        path = tmp_path / f"{tag}_host_{seed}.h5"
        _write_run(path, values)
        host_runs.append(load_run(path, "E1", seed))
    for seed, values in device_values.items():
        path = tmp_path / f"{tag}_device_{seed}.h5"
        _write_run(
            path,
            values,
            clip=device_clip,
            oxygen=device_oxygen,
        )
        device_runs.append(load_run(path, "E2", seed))
    repeats = []
    if repeat_values is not None:
        seed = min(device_values)
        path = tmp_path / f"{tag}_repeat_{seed}.h5"
        _write_run(
            path,
            repeat_values,
            clip=device_clip,
            oxygen=device_oxygen,
        )
        repeats.append((device_runs[0], load_run(path, "E2r", seed)))
    return compare(host_runs, device_runs, repeats)


def _verdicts(record: dict) -> dict[str, dict]:
    return {item["estimand"]: item for item in record["verdicts"]}


def test_backend_shift_inside_seed_spread_is_interchangeable(
    tmp_path: Path,
) -> None:
    host = {
        55: [10, 40, 70, 100],
        56: [10, 45, 75, 110],
        57: [10, 50, 80, 120],
        58: [10, 55, 85, 130],
        59: [10, 60, 90, 140],
    }
    device = {seed: values[:-1] + [values[-1] + 2]
              for seed, values in host.items()}
    record = _load_pair(
        tmp_path, "inside", host, device, repeat_values=device[55]
    )
    assert _verdicts(record)["n_total"]["verdict"] == "interchangeable"


def test_backend_shift_beyond_seed_spread_is_flagged(tmp_path: Path) -> None:
    host = {
        55: [10, 40, 70, 100],
        56: [10, 45, 75, 110],
        57: [10, 50, 80, 120],
        58: [10, 55, 85, 130],
        59: [10, 60, 90, 140],
    }
    device = {seed: values[:-1] + [values[-1] + 500]
              for seed, values in host.items()}
    record = _load_pair(
        tmp_path, "beyond", host, device, repeat_values=device[55]
    )
    assert _verdicts(record)["n_total"]["verdict"] == "flagged"


def test_device_reaction_clip_surge_is_an_invariant_difference(
    tmp_path: Path,
) -> None:
    values = {seed: [10, 40, 70, 100] for seed in range(55, 60)}
    record = _load_pair(
        tmp_path,
        "clip",
        values,
        values,
        device_clip=[0.0, 1.0, 2.0, 3.0],
        repeat_values=values[55],
    )
    verdicts = _verdicts(record)
    assert verdicts["reaction_clip_cumulative_total"]["verdict"] == (
        "invariant_differs"
    )
    assert "reaction_clip_cumulative_total" in record["flagged_estimands"]


def test_nonidentical_device_repeat_fails_reproducibility(tmp_path: Path) -> None:
    values = {seed: [10, 40, 70, 100] for seed in range(55, 60)}
    record = _load_pair(
        tmp_path,
        "repeat",
        values,
        values,
        repeat_values=[10, 40, 70, 103],
    )
    assert record["reproducibility"]["reproducible"] is False


def test_nonfinite_trajectory_is_flagged_even_with_finite_final_value(
    tmp_path: Path,
) -> None:
    values = {seed: [10, 40, 70, 100] for seed in range(55, 60)}
    record = _load_pair(
        tmp_path,
        "nonfinite",
        values,
        values,
        device_oxygen=[1.0, 1.0, np.nan, 1.0],
        repeat_values=values[55],
    )
    verdicts = _verdicts(record)
    assert verdicts["mean_oxygen"]["verdict"] == "non_finite"
    assert "mean_oxygen" in record["flagged_estimands"]


def test_mann_whitney_exact_two_sided_floor_and_interleaving() -> None:
    separated = mann_whitney_exact_two_sided(
        [1.0, 2.0, 3.0, 4.0, 5.0],
        [6.0, 7.0, 8.0, 9.0, 10.0],
    )
    assert separated["p_two_sided"] == pytest.approx(2.0 / 252.0)
    assert separated["p_floor"] == pytest.approx(2.0 / 252.0)

    interleaved = mann_whitney_exact_two_sided(
        [1.0, 3.0, 5.0, 7.0, 9.0],
        [2.0, 4.0, 6.0, 8.0, 10.0],
    )
    assert interleaved["p_two_sided"] > separated["p_two_sided"]
