"""Regression tests for physical comet-tail geometry and periodic metrics."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest
from gut_ibm_tools import GutIBMData, validation
from gut_ibm_tools.analysis import (
    comet_tail_asymmetry_index,
    comet_tail_index,
)


def _write_grid_file(
    path: Path,
    *,
    nx: int = 40,
    ny: int = 8,
    nz: int = 4,
    with_domain: bool = True,
) -> None:
    with h5py.File(path, "w") as handle:
        handle.attrs["nx"] = nx
        handle.attrs["ny"] = ny
        handle.attrs["nz"] = nz
        handle.attrs["grid_dx"] = 1.0
        if with_domain:
            handle.attrs["domain_lo_x"] = 0.0
            handle.attrs["domain_lo_y"] = 0.0
            handle.attrs["domain_lo_z"] = 0.0
            handle.attrs["domain_hi_x"] = float(nx)
            handle.attrs["domain_hi_y"] = float(ny)
            handle.attrs["domain_hi_z"] = float(nz)
        grid = handle.require_group("grid").require_group("step_000000")
        grid.create_dataset(
            "marker",
            data=np.arange(nx * ny * nz).reshape(nz, ny, nx),
        )


def _physical_field(
    centers: np.ndarray,
    *,
    reference: np.ndarray,
    x_scale: float,
    y_scale: float = 1.5,
    z_scale: float = 0.8,
) -> np.ndarray:
    displacement = centers - reference
    x_sigma = np.where(displacement[:, 0] >= 0, x_scale, x_scale / 6)
    return np.exp(
        -0.5 * (
            (displacement[:, 0] / x_sigma) ** 2
            + (displacement[:, 1] / y_scale) ** 2
            + (displacement[:, 2] / z_scale) ** 2
        )
    )


def test_grid_cell_centers_preserve_flattened_order(tmp_path: Path) -> None:
    path = tmp_path / "geometry.h5"
    _write_grid_file(path, nx=3, ny=2, nz=2)

    with GutIBMData(path) as data:
        assert data.grid_shape() == (3, 2, 2)
        centers = data.grid_cell_centers()
        marker = data.get_grid("step_000000")["marker"]
        periods = data.grid_periods()

    expected = np.array([
        [0.5, 0.5, 0.5],
        [1.5, 0.5, 0.5],
        [2.5, 0.5, 0.5],
        [0.5, 1.5, 0.5],
        [1.5, 1.5, 0.5],
        [2.5, 1.5, 0.5],
        [0.5, 0.5, 1.5],
        [1.5, 0.5, 1.5],
        [2.5, 0.5, 1.5],
        [0.5, 1.5, 1.5],
        [1.5, 1.5, 1.5],
        [2.5, 1.5, 1.5],
    ])
    np.testing.assert_allclose(marker, np.arange(12))
    np.testing.assert_allclose(centers, expected)
    assert np.all(centers > 0)
    assert np.all(centers < np.array([3.0, 2.0, 2.0]))
    np.testing.assert_allclose(centers[0], [0.5, 0.5, 0.5])
    np.testing.assert_allclose(centers[-1], [2.5, 1.5, 1.5])
    np.testing.assert_allclose(periods, [3.0, 2.0, 2.0])


def test_grid_cell_centers_use_zero_origin_spacing_fallback(tmp_path: Path) -> None:
    path = tmp_path / "legacy_geometry.h5"
    _write_grid_file(path, nx=2, ny=2, nz=1, with_domain=False)
    with h5py.File(path, "a") as handle:
        handle.attrs["grid_dx_x"] = 2.0
        handle.attrs["grid_dx_y"] = 3.0
        handle.attrs["grid_dx_z"] = 4.0

    with GutIBMData(path) as data:
        centers = data.grid_cell_centers()
        periods = data.grid_periods()

    np.testing.assert_allclose(centers, [
        [1.0, 1.5, 2.0],
        [3.0, 1.5, 2.0],
        [1.0, 4.5, 2.0],
        [3.0, 4.5, 2.0],
    ])
    np.testing.assert_allclose(periods, [4.0, 6.0, 4.0])


def test_symmetric_three_dimensional_field_has_no_tail(tmp_path: Path) -> None:
    path = tmp_path / "symmetric.h5"
    _write_grid_file(path)
    with GutIBMData(path) as data:
        centers = data.grid_cell_centers()

    reference = np.array([20.0, 4.0, 2.0])
    displacement = centers - reference
    field = np.exp(-0.5 * (
        (displacement[:, 0] / 4.0) ** 2
        + (displacement[:, 1] / 1.5) ** 2
        + (displacement[:, 2] / 0.8) ** 2
    ))
    ratio = comet_tail_index(
        centers,
        field,
        reference=reference,
        period=40.0,
    )
    asymmetry = comet_tail_asymmetry_index(
        centers,
        field,
        reference=reference,
        period=40.0,
    )
    assert ratio == pytest.approx(1.0, rel=1e-12)
    assert asymmetry == pytest.approx(1.0, rel=1e-12)


def test_real_geometry_detects_plume_hidden_by_fabricated_line(tmp_path: Path) -> None:
    path = tmp_path / "plume.h5"
    _write_grid_file(path)
    with GutIBMData(path) as data:
        centers = data.grid_cell_centers()

    reference = np.array([20.0, 4.0, 2.0])
    field = _physical_field(centers, reference=reference, x_scale=12.0)
    real_ratio = comet_tail_index(
        centers,
        field,
        reference=reference,
        period=40.0,
    )
    real_asymmetry = comet_tail_asymmetry_index(
        centers,
        field,
        reference=reference,
        period=40.0,
    )
    fabricated = np.column_stack([
        np.linspace(0.0, 1e-3, len(field)),
        np.zeros(len(field)),
        np.zeros(len(field)),
    ])
    fabricated_ratio = comet_tail_index(fabricated, field)
    fabricated_asymmetry = comet_tail_asymmetry_index(fabricated, field)

    assert real_ratio > 1.5
    assert real_asymmetry > 1.0
    assert fabricated_ratio == pytest.approx(1.0, abs=0.1)
    assert fabricated_asymmetry == pytest.approx(1.0, abs=0.1)


def test_periodic_wrap_matches_translated_configuration() -> None:
    edge_positions = np.array([
        [37.5, 0.0, 0.0],
        [38.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [1.5, 0.0, 0.0],
    ])
    middle_positions = np.array([
        [17.5, 0.0, 0.0],
        [18.5, 0.0, 0.0],
        [20.5, 0.0, 0.0],
        [21.5, 0.0, 0.0],
    ])
    concentrations = np.array([1.0, 2.0, 3.0, 4.0])
    edge_reference = np.array([39.5, 0.0, 0.0])
    middle_reference = np.array([19.5, 0.0, 0.0])

    edge_ratio = comet_tail_index(
        edge_positions,
        concentrations,
        reference=edge_reference,
        period=40.0,
    )
    middle_ratio = comet_tail_index(
        middle_positions,
        concentrations,
        reference=middle_reference,
        period=40.0,
    )
    edge_asymmetry = comet_tail_asymmetry_index(
        edge_positions,
        concentrations,
        reference=edge_reference,
        period=40.0,
    )
    middle_asymmetry = comet_tail_asymmetry_index(
        middle_positions,
        concentrations,
        reference=middle_reference,
        period=40.0,
    )

    assert edge_ratio == pytest.approx(middle_ratio)
    assert edge_asymmetry == pytest.approx(middle_asymmetry)


@pytest.mark.parametrize("elongated_axis", [1, 2])
def test_transverse_elongation_is_not_an_x_comet_tail(
    tmp_path: Path,
    elongated_axis: int,
) -> None:
    path = tmp_path / f"transverse_{elongated_axis}.h5"
    _write_grid_file(path)
    with GutIBMData(path) as data:
        centers = data.grid_cell_centers()

    reference = np.array([20.0, 4.0, 2.0])
    displacement = centers - reference
    scales = np.array([3.0, 10.0, 10.0])
    scales[3 - elongated_axis] = 0.8
    field = np.exp(-0.5 * np.sum((displacement / scales) ** 2, axis=1))
    ratio = comet_tail_index(
        centers,
        field,
        reference=reference,
        period=40.0,
    )
    asymmetry = comet_tail_asymmetry_index(
        centers,
        field,
        reference=reference,
        period=40.0,
    )

    assert ratio == pytest.approx(1.0, rel=1e-12)
    assert asymmetry == pytest.approx(1.0, rel=1e-12)


def test_validation_rejects_mismatched_toxin_length(
    sample_hdf5: Path,
    tmp_path: Path,
) -> None:
    path = tmp_path / "mismatched_toxin.h5"
    path.write_bytes(sample_hdf5.read_bytes())
    with h5py.File(path, "a") as handle:
        dataset = handle["grid/step_000000/bacteriocin_BtuB"]
        values = dataset[:].ravel()[:-1]
        del handle["grid/step_000000/bacteriocin_BtuB"]
        handle["grid/step_000000"].create_dataset(
            "bacteriocin_BtuB",
            data=values,
        )

    with GutIBMData(path) as data, pytest.raises(
        ValueError, match="flattened toxin field length"
    ):
        validation.validate_spatial_signatures(data, "step_000000")


def test_toxin_reference_prefers_biomass_weighted_producers() -> None:
    grid_positions = np.column_stack([
        np.arange(-1.0, 10.0, 2.0),
        np.zeros(6),
        np.zeros(6),
    ])
    agents = {
        "x": np.array([1.0, 3.0, 9.0]),
        "y": np.zeros(3),
        "z": np.zeros(3),
        "biomass": np.array([1.0, 3.0, 4.0]),
        "n_bi_loci": np.array([1, 1, 0]),
    }
    concentrations = np.array([1.0, 1.0, 10.0, 10.0, 10.0, 10.0])

    reference = validation._toxin_source_reference(agents, grid_positions)
    all_agent_reference = np.array([5.75, 0.0, 0.0])
    producer_ratio = comet_tail_index(
        grid_positions,
        concentrations,
        reference=reference,
    )
    all_agent_ratio = comet_tail_index(
        grid_positions,
        concentrations,
        reference=all_agent_reference,
    )

    assert reference[0] == pytest.approx(2.5)
    assert not np.isclose(all_agent_reference[0], reference[0])
    assert producer_ratio == pytest.approx(10.0)
    assert all_agent_ratio == pytest.approx(10.0 / 5.5)
    assert producer_ratio - all_agent_ratio > 5.0


@pytest.mark.parametrize(
    "agents",
    [
        {
            "x": np.array([1.0, 3.0, 9.0]),
            "y": np.zeros(3),
            "z": np.zeros(3),
            "biomass": np.array([1.0, 3.0, 4.0]),
            "n_bi_loci": np.array([0, 0, 0]),
        },
        {
            "x": np.array([1.0, 3.0, 9.0]),
            "y": np.zeros(3),
            "z": np.zeros(3),
            "biomass": np.array([1.0, 3.0, 4.0]),
        },
    ],
    ids=["no-producers", "missing-producer-field"],
)
def test_toxin_reference_falls_back_to_all_agents(
    agents: dict[str, np.ndarray],
) -> None:
    grid_positions = np.column_stack([
        np.arange(-1.0, 10.0, 2.0),
        np.zeros(6),
        np.zeros(6),
    ])
    concentrations = np.array([1.0, 1.0, 10.0, 10.0, 10.0, 10.0])

    reference = validation._toxin_source_reference(agents, grid_positions)
    all_agent_ratio = comet_tail_index(
        grid_positions,
        concentrations,
        reference=reference,
    )

    assert reference[0] == pytest.approx(5.75)
    assert all_agent_ratio == pytest.approx(10.0 / 5.5)


def test_toxin_reference_uses_domain_midpoint_without_agents() -> None:
    grid_positions = np.column_stack([
        np.arange(0.5, 10.0, 1.0),
        np.zeros(10),
        np.zeros(10),
    ])
    concentrations = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 4.0, 4.0, 4.0, 4.0, 4.0])

    reference = validation._toxin_source_reference({}, grid_positions)
    midpoint_ratio = comet_tail_index(
        grid_positions,
        concentrations,
        reference=reference,
    )

    assert reference == pytest.approx([5.0, 0.0, 0.0])
    assert midpoint_ratio == pytest.approx(4.0)
