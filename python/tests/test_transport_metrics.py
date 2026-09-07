"""Tests for authenticated source-intrinsic transport metrics."""

from __future__ import annotations

import json

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from gut_ibm_tools.transport_metrics import (
    Box,
    ExpectedRun,
    TransportInputError,
    authenticate_run,
    paired_times,
    profile_radii,
    read_lysis_sources,
    release_weight,
    shell_profile,
)

SHA = "a" * 40
SEED = 20260911
DX = 2.0e-6
N = 50
BOX = Box(Lx=N * DX, Ly=N * DX, Lz=N * DX)
CENTER = (0.5 * BOX.Lx, 0.5 * BOX.Ly, 0.5 * BOX.Lz)
R_MAX_UM = 40.0
SHELL_UM = 2.0


def resolved_config(**overrides) -> dict:
    cfg = {
        "seed": SEED,
        "bio_dt": 60.0,
        "grid_dx": DX,
        "domain_x": BOX.Lx,
        "domain_y": BOX.Ly,
        "domain_z": BOX.Lz,
        "bacteriocin.mucin_charge.amplitude": 15.0,
        "kd_b12_btuB": 1.0e-4,
        "kd_colicinE_btuB": 5.0e-7,
        "b12.initial_conc": 1.0e-3,
        "burst_release_tau": 300.0,
        "hdf5.schedule.summary": 1,
        "hdf5.schedule.agents": 10,
        "hdf5.schedule.grid": 60,
        "hdf5.schedule.lineage": 0,
        "hdf5.schedule.genome": 0,
        "hdf5.schedule.provenance": 10,
        "hdf5.schedule.grid_species": ["bacteriocin_BtuB"],
        "plasmid_overrides": {},
    }
    cfg.update(overrides)
    return cfg


def write_run(path, *, config=None, sha=SHA, placement="device", ranks=1, events=None):
    """Write a minimal output file with run provenance and a provenance layer."""
    with h5py.File(path, "w") as h:
        group = h.create_group("run_provenance")
        group.create_dataset("resolved_config", data=json.dumps(config or resolved_config()))
        group.create_dataset("git_sha", data=sha)
        group.create_dataset("chemistry_placement", data=placement)
        group.create_dataset("mpi_rank_count", data=np.int32(ranks))
        for step, records in (events or {}).items():
            prov = h.create_group(f"provenance/step_{step:06d}")
            for name, values in records.items():
                prov.create_dataset(name, data=np.asarray(values))
    return path


def lysis_records(*, event_time_s=1200.0, event_step=20, exact=True):
    records = {
        "cause": np.array([5], dtype=np.int32),
        "strain": np.array([1], dtype=np.int32),
        "x": np.array([CENTER[0]]),
        "y": np.array([CENTER[1]]),
        "z": np.array([CENTER[2]]),
    }
    if exact:
        records["event_time_s"] = np.array([event_time_s])
        records["event_step"] = np.array([event_step], dtype=np.int64)
    return records


def exponential_field(decay_um: float) -> np.ndarray:
    """Concentration ~ exp(-r / decay) / r around the domain centre."""
    axis = (np.arange(N, dtype=float) + 0.5) * DX
    dz = axis[:, None, None] - CENTER[2]
    dy = axis[None, :, None] - CENTER[1]
    dx = axis[None, None, :] - CENTER[0]
    r_um = np.sqrt(dz**2 + dy**2 + dx**2) * 1.0e6
    r_um = np.maximum(r_um, 0.5 * DX * 1.0e6)
    return np.exp(-r_um / decay_um) / r_um


def profile_for(decay_um: float):
    return shell_profile(exponential_field(decay_um), DX, CENTER, BOX, R_MAX_UM, SHELL_UM)


class TestSensitivity:
    """Graded response: steeper decay must give smaller profile radii."""

    def test_radii_order_with_decay_length(self):
        decays = [2.0, 5.0, 12.0, 30.0]
        radii = [profile_radii(profile_for(d)) for d in decays]
        r50 = [r[0] for r in radii]
        r90 = [r[1] for r in radii]

        assert r50 == sorted(r50)
        assert r90 == sorted(r90)
        assert r50[-1] - r50[0] > 5.0
        assert r90[-1] - r90[0] > 5.0
        assert all(a < b for a, b in zip(r50, r90))
        assert all(0.0 < r <= R_MAX_UM for r in r50 + r90)

    def test_shell_means_decay_monotonically(self):
        profile = profile_for(8.0)
        means = profile.means
        assert np.all(np.isfinite(means))
        assert np.all(np.diff(means) < 0.0)
        assert profile.complete

    def test_radii_are_independent_of_source_placement(self):
        """The metric is intrinsic: an off-centre source with equal support agrees."""
        offset = (CENTER[0] + 10 * DX, CENTER[1] - 6 * DX, CENTER[2])
        axis = (np.arange(N, dtype=float) + 0.5) * DX
        dz = axis[:, None, None] - offset[2]
        dy = axis[None, :, None] - offset[1]
        dx = axis[None, None, :] - offset[0]
        dy -= BOX.Ly * np.round(dy / BOX.Ly)
        dx -= BOX.Lx * np.round(dx / BOX.Lx)
        r_um = np.maximum(np.sqrt(dz**2 + dy**2 + dx**2) * 1.0e6, 0.5 * DX * 1.0e6)
        field = np.exp(-r_um / 8.0) / r_um

        moved = profile_radii(shell_profile(field, DX, offset, BOX, R_MAX_UM, SHELL_UM))
        centred = profile_radii(profile_for(8.0))
        assert moved[0] == pytest.approx(centred[0], rel=0.02)
        assert moved[1] == pytest.approx(centred[1], rel=0.02)


class TestBoundaryHandling:
    def test_support_beyond_domain_is_refused(self):
        near_wall = (CENTER[0], CENTER[1], 10 * DX)
        with pytest.raises(TransportInputError, match="exceeds the in-domain support"):
            shell_profile(exponential_field(8.0), DX, near_wall, BOX, R_MAX_UM, SHELL_UM)

    def test_support_reports_the_binding_constraint(self):
        profile = profile_for(8.0)
        assert profile.support_um == pytest.approx(50.0)

    def test_non_integral_shell_count_is_refused(self):
        with pytest.raises(TransportInputError, match="whole shells"):
            shell_profile(exponential_field(8.0), DX, CENTER, BOX, 41.0, SHELL_UM)

    def test_empty_shell_yields_undefined_radii(self):
        profile = shell_profile(exponential_field(8.0), DX, CENTER, BOX, 4.0, 0.25)
        assert not profile.complete
        assert all(np.isnan(r) for r in profile_radii(profile))


class TestExactEventTiming:
    def test_exact_event_time_is_used(self, tmp_path):
        path = write_run(tmp_path / "exact.h5", events={30: lysis_records(event_time_s=1200.0)})
        with h5py.File(path, "r") as h:
            sources = read_lysis_sources(h)
        assert len(sources) == 1
        assert sources[0].time_s == pytest.approx(1200.0)
        assert sources[0].step == 20

    def test_missing_event_time_blocks_analysis(self, tmp_path):
        path = write_run(tmp_path / "legacy.h5", events={30: lysis_records(exact=False)})
        with h5py.File(path, "r") as h, pytest.raises(TransportInputError, match="event_time_s"):
            read_lysis_sources(h)

    def test_release_weight_decays_and_prunes(self, tmp_path):
        path = write_run(tmp_path / "w.h5", events={30: lysis_records(event_time_s=1200.0)})
        with h5py.File(path, "r") as h:
            source = read_lysis_sources(h)[0]
        assert release_weight(source, 1200.0, 300.0) == pytest.approx(1.0)
        assert release_weight(source, 1500.0, 300.0) == pytest.approx(np.exp(-1.0))
        assert release_weight(source, 1100.0, 300.0) == pytest.approx(0.0)
        assert release_weight(source, 1200.0 + 1600.0, 300.0) == pytest.approx(0.0)

    def test_export_step_delay_would_have_inflated_the_weight(self, tmp_path):
        """The export step is up to one interval late; that is not the event time."""
        path = write_run(
            tmp_path / "delay.h5", events={30: lysis_records(event_time_s=1200.0, event_step=20)}
        )
        with h5py.File(path, "r") as h:
            source = read_lysis_sources(h)[0]
        exact = release_weight(source, 2600.0, 300.0)
        export_time_substitute = np.exp(-(2600.0 - 1800.0) / 300.0)
        assert export_time_substitute / exact > 5.0


class TestAuthentication:
    def expected(self, **kwargs) -> ExpectedRun:
        base = {"execution_source_sha": SHA, "seed": SEED, "amplitude": 15.0}
        base.update(kwargs)
        return ExpectedRun(**base)

    def test_matching_run_has_no_violations(self, tmp_path):
        path = write_run(tmp_path / "ok.h5")
        with h5py.File(path, "r") as h:
            assert authenticate_run(h, self.expected()) == []

    @pytest.mark.parametrize(
        ("kwargs", "needle"),
        [
            ({"sha": "b" * 40}, "git_sha"),
            ({"placement": "host"}, "chemistry_placement"),
            ({"placement": "host_forced_delivery"}, "chemistry_placement"),
            ({"ranks": 2}, "mpi_rank_count"),
        ],
    )
    def test_runtime_identity_violations(self, tmp_path, kwargs, needle):
        path = write_run(tmp_path / "identity.h5", **kwargs)
        with h5py.File(path, "r") as h:
            violations = authenticate_run(h, self.expected())
        assert any(needle in v for v in violations)

    @pytest.mark.parametrize(
        ("key", "value", "needle"),
        [
            ("seed", 1234, "seed"),
            ("bacteriocin.mucin_charge.amplitude", 60.0, "amplitude"),
            ("kd_b12_btuB", 1.0e-6, "kd_b12_btuB"),
            ("b12.initial_conc", 1.0e-5, "b12.initial_conc"),
            ("burst_release_tau", 60.0, "burst_release_tau"),
            ("hdf5.schedule.provenance", 20, "hdf5.schedule.provenance"),
            ("hdf5.schedule.grid", 0, "hdf5.schedule.grid"),
            ("hdf5.schedule.grid_species", ["carbon"], "grid_species"),
        ],
    )
    def test_scientific_setting_violations(self, tmp_path, key, value, needle):
        path = write_run(tmp_path / "cfg.h5", config=resolved_config(**{key: value}))
        with h5py.File(path, "r") as h:
            violations = authenticate_run(h, self.expected())
        assert any(needle in v for v in violations)

    def test_cole1_overrides_are_authenticated(self, tmp_path):
        config = resolved_config(
            plasmid_overrides={"ColE1": {"burst_size": 1.0e3, "retardation": 2.0}}
        )
        path = write_run(tmp_path / "override.h5", config=config)
        with h5py.File(path, "r") as h:
            violations = authenticate_run(h, self.expected())
        assert any("burst_size" in v for v in violations)
        assert any("retardation" in v for v in violations)

    def test_matching_cole1_override_is_accepted(self, tmp_path):
        config = resolved_config(plasmid_overrides={"ColE1": {"burst_size": 1.0e5}})
        path = write_run(tmp_path / "same.h5", config=config)
        with h5py.File(path, "r") as h:
            assert authenticate_run(h, self.expected()) == []

    def test_genome_layer_locus_is_authenticated(self, tmp_path):
        path = write_run(tmp_path / "genome.h5")
        with h5py.File(path, "a") as h:
            group = h.create_group("genome/step_000010")
            group.create_dataset("bi_pI", data=np.array([6.5]))
            group.create_dataset("bi_diff_coeff", data=np.array([3.5e-11]))
        with h5py.File(path, "r") as h:
            violations = authenticate_run(h, self.expected())
        assert any("bi_pI" in v for v in violations)
        assert any("bi_diff_coeff" in v for v in violations)

    def test_absent_resolved_config_is_an_error(self, tmp_path):
        path = tmp_path / "bare.h5"
        with h5py.File(path, "w") as h:
            h.create_group("run_provenance")
        with h5py.File(path, "r") as h, pytest.raises(TransportInputError):
            authenticate_run(h, self.expected())


class TestPairedTimes:
    def test_identical_times_pair_completely(self):
        pairs = paired_times({0: [600.0, 1200.0], 15: [600.0, 1200.0], 60: [600.0, 1200.0]})
        assert pairs.complete
        np.testing.assert_allclose(pairs.common_s, [600.0, 1200.0])

    def test_missing_time_blocks_instead_of_shrinking_support(self):
        pairs = paired_times({0: [600.0, 1200.0], 15: [600.0], 60: [600.0, 1200.0]})
        assert not pairs.complete
        np.testing.assert_allclose(pairs.common_s, [600.0])
        assert pairs.missing_s[15] == [1200.0]
        assert pairs.missing_s[0] == []

    def test_no_arms_pair_to_nothing(self):
        pairs = paired_times({})
        assert pairs.common_s == []
        assert pairs.complete
