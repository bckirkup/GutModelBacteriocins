"""Authenticated, source-intrinsic bacteriocin transport metrics.

Two facilities are provided:

Input authentication
    :func:`authenticate_run` compares an output file against the exact
    scientific and runtime settings a transport analysis claims to have
    analyzed (execution source commit, chemistry placement, rank count, seed,
    mucin-charge amplitude, receptor dissociation constants, corrinoid
    concentration, burst release timescale, ColE1 identity, and the HDF5
    grid/provenance schedules).  It returns violations instead of raising so a
    caller can report every mismatch of a campaign in one pass.

Source-intrinsic radial profiles
    Radii computed against the *nearest* of many simultaneously active sources
    measure source packing: sparse sources own larger Voronoi cells and
    therefore report larger radii at identical diffusivity.  The profile here
    is instead measured around a single known source, on fixed shells, over a
    radial support that is fully contained in the domain, so it depends on the
    release and transport law alone.

Predeclared definitions used by both the assay and any gate built on it:

shell averaging
    Shells are fixed-width spherical shells centred on the source.  A shell
    value is the unweighted arithmetic mean of the concentrations of all voxel
    centres whose radius falls in ``[r_lo, r_hi)`` (measure ``dr``), so shells
    are not weighted by how many voxels a radius happens to contain.
radial support
    Only shells with ``r_hi <= support`` are used, where ``support`` is the
    largest radius whose full sphere around the source lies inside the domain
    (half the shorter periodic span in x/y, distance to the nearer z wall).
    A requested ``r_max`` above the support is a hard error, not a silent trim.
boundary handling
    x and y use minimum-image distances (periodic); z is non-periodic.  No
    image sources are added: the support rule keeps the profile away from the
    z walls instead of modelling reflections.
profile r50/r90
    Shell means are converted to a radial mass density ``m(r) * r**2`` and
    accumulated over the support.  ``rXX`` is the radius at which the
    normalized cumulative mass first reaches ``XX/100``, linearly interpolated
    between shell centres (and from ``F(0) = 0`` inside the first shell).
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

CAUSE_LYSIS = 5
PRODUCER_STRAIN = 1
COLE1_LIBRARY_PI = 9.0
COLE1_LIBRARY_DIFF_COEFF = 4.0e-11
COLE1_LIBRARY_BURST_SIZE = 1.0e5


class TransportInputError(RuntimeError):
    """An input cannot support the requested transport analysis."""


def _decode(raw: Any) -> str:
    if isinstance(raw, bytes):
        return raw.decode()
    if hasattr(raw, "tobytes"):
        return raw.tobytes().decode()
    return str(raw)


def read_string(h5, path: str) -> str | None:
    """Return a string dataset, or None when it is absent."""
    if path not in h5:
        return None
    return _decode(h5[path][()])


def read_resolved_config(h5) -> dict:
    """Return the runtime resolved configuration document."""
    raw = read_string(h5, "run_provenance/resolved_config")
    if raw is None:
        raise TransportInputError("run_provenance/resolved_config is absent")
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        raise TransportInputError(f"resolved_config is not valid JSON: {exc}") from exc


def _scalar(h5, path: str) -> float | None:
    if path not in h5:
        return None
    return float(np.asarray(h5[path][()]).reshape(-1)[0])


@dataclass(frozen=True)
class ExpectedRun:
    """The settings a transport analysis requires of one run.

    Defaults carry shared transport settings and physical GPU placement
    (``device``); callers exercising ecological device delivery must request
    ``device_delivery`` explicitly. Per-run commit, seed, and amplitude have no defaults.
    """

    execution_source_sha: str
    seed: int
    amplitude: float
    kd_corrinoid_btuB: float = 1.0e-4
    kd_colicinE_btuB: float = 5.0e-7
    b12_initial_conc: float = 1.0e-3
    burst_release_tau_s: float = 300.0
    cole1_pI: float = COLE1_LIBRARY_PI
    cole1_diff_coeff: float = COLE1_LIBRARY_DIFF_COEFF
    cole1_burst_size: float = COLE1_LIBRARY_BURST_SIZE
    mpi_ranks: int = 1
    chemistry_placement: str = "device"
    grid_species: tuple[str, ...] = ("bacteriocin_BtuB",)
    hdf5_schedule: Mapping[str, int] = field(
        default_factory=lambda: {
            "summary": 1,
            "agents": 10,
            "grid": 60,
            "lineage": 0,
            "genome": 0,
            "provenance": 10,
        }
    )
    rtol: float = 1.0e-9


def _close(observed: float | None, expected: float, rtol: float) -> bool:
    if observed is None:
        return False
    return abs(observed - expected) <= rtol * abs(expected)


def authenticate_run(h5, expected: ExpectedRun) -> list[str]:
    """Return every way ``h5`` fails to be the run ``expected`` describes.

    An empty list means the file's runtime provenance and resolved
    configuration match the claimed execution image and scientific settings
    exactly, so a transport measurement taken from it can be attributed to
    those settings.
    """
    violations: list[str] = []

    git_sha = read_string(h5, "run_provenance/git_sha")
    if git_sha != expected.execution_source_sha:
        violations.append(
            f"run_provenance/git_sha {git_sha!r} != execution_source_sha "
            f"{expected.execution_source_sha!r}"
        )

    placement = read_string(h5, "run_provenance/chemistry_placement")
    if placement != expected.chemistry_placement:
        violations.append(
            f"chemistry_placement {placement!r} != {expected.chemistry_placement!r}"
        )

    ranks = _scalar(h5, "run_provenance/mpi_rank_count")
    if ranks is None or int(ranks) != expected.mpi_ranks:
        violations.append(f"mpi_rank_count {ranks} != {expected.mpi_ranks}")

    cfg = read_resolved_config(h5)

    if int(cfg.get("seed", -1)) != int(expected.seed):
        violations.append(f"resolved seed {cfg.get('seed')} != {expected.seed}")

    checks: list[tuple[str, float]] = [
        ("bacteriocin.mucin_charge.amplitude", expected.amplitude),
        ("kd_b12_btuB", expected.kd_corrinoid_btuB),
        ("kd_colicinE_btuB", expected.kd_colicinE_btuB),
        ("b12.initial_conc", expected.b12_initial_conc),
        ("burst_release_tau", expected.burst_release_tau_s),
    ]
    for key, want in checks:
        observed = cfg.get(key)
        value = float(observed) if isinstance(observed, (int, float)) else None
        if not _close(value, want, expected.rtol):
            violations.append(f"resolved {key} {observed} != {want}")

    for key, want in expected.hdf5_schedule.items():
        observed = cfg.get(f"hdf5.schedule.{key}")
        if observed is None or int(observed) != int(want):
            violations.append(f"resolved hdf5.schedule.{key} {observed} != {want}")

    species = tuple(cfg.get("hdf5.schedule.grid_species") or ())
    if species != tuple(expected.grid_species):
        violations.append(f"resolved grid_species {species} != {tuple(expected.grid_species)}")

    violations.extend(_cole1_violations(h5, cfg, expected))
    return violations


def _cole1_violations(h5, cfg: Mapping[str, Any], expected: ExpectedRun) -> list[str]:
    """Authenticate the ColE1 release identity.

    ``pI``, ``diff_coeff`` and ``burst_size`` are plasmid-library constants of
    the execution source, so an authenticated commit plus the absence of a
    conflicting ``plasmid_overrides`` entry pins them.  When the run also
    emitted the genome layer, the producer's realized locus is compared
    directly.
    """
    violations: list[str] = []
    overrides = (cfg.get("plasmid_overrides") or {}).get("ColE1") or {}
    override_expectations = {
        "diff_coeff": expected.cole1_diff_coeff,
        "burst_size": expected.cole1_burst_size,
    }
    for key, value in overrides.items():
        if key == "retardation":
            violations.append(
                "plasmid_overrides.ColE1.retardation overrides the mucin-charge "
                "amplitude under test"
            )
            continue
        want = override_expectations.get(key)
        if want is None or not _close(float(value), want, expected.rtol):
            violations.append(f"plasmid_overrides.ColE1.{key} {value} != {want}")

    if "genome" not in h5:
        return violations
    for step_key in sorted(h5["genome"].keys()):
        group = h5["genome"][step_key]
        if "bi_pI" not in group or group["bi_pI"].shape[0] == 0:
            continue
        pI = np.asarray(group["bi_pI"][()], dtype=float)
        diff = np.asarray(group["bi_diff_coeff"][()], dtype=float)
        if not np.all(np.isclose(pI, expected.cole1_pI, rtol=expected.rtol, atol=0.0)):
            violations.append(f"genome/{step_key} bi_pI {sorted(set(pI))} != {expected.cole1_pI}")
        if not np.all(np.isclose(diff, expected.cole1_diff_coeff, rtol=expected.rtol, atol=0.0)):
            violations.append(
                f"genome/{step_key} bi_diff_coeff {sorted(set(diff))} != "
                f"{expected.cole1_diff_coeff}"
            )
        break
    return violations


@dataclass(frozen=True)
class LysisSource:
    """One producer lysis, at the simulation clock value of the lysis itself."""

    time_s: float
    step: int
    x: float
    y: float
    z: float
    strain: int

    @property
    def position(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


def _step_index(name: str) -> int:
    return int(name.rsplit("_", 1)[1])


def sorted_steps(group) -> list[str]:
    return sorted(group.keys(), key=_step_index)


def read_lysis_sources(h5, producer_strain: int = PRODUCER_STRAIN) -> list[LysisSource]:
    """Return producer lysis sources stamped with exact event times.

    Raises when the provenance layer lacks ``event_time_s``/``event_step``:
    the provenance-export step of a buffered event is an upper bound on the
    event time (up to one export interval late), and the release window decays
    on a 300 s timescale, so an export-time substitute is not a source time.
    """
    if "provenance" not in h5:
        raise TransportInputError("no provenance layer: cannot locate lysis sources")

    sources: list[LysisSource] = []
    for step_key in sorted_steps(h5["provenance"]):
        group = h5["provenance"][step_key]
        if "cause" not in group or group["cause"].shape[0] == 0:
            continue
        if "event_time_s" not in group or "event_step" not in group:
            raise TransportInputError(
                f"provenance/{step_key} has no event_time_s/event_step: this output "
                "predates exact provenance event timing and cannot be used for "
                "transport analysis"
            )
        cause = np.asarray(group["cause"][()])
        strain = np.asarray(group["strain"][()])
        event_time = np.asarray(group["event_time_s"][()], dtype=float)
        event_step = np.asarray(group["event_step"][()], dtype=float)
        x = np.asarray(group["x"][()], dtype=float)
        y = np.asarray(group["y"][()], dtype=float)
        z = np.asarray(group["z"][()], dtype=float)
        selected = (cause == CAUSE_LYSIS) & (strain == producer_strain)
        for i in np.flatnonzero(selected):
            sources.append(
                LysisSource(
                    time_s=float(event_time[i]),
                    step=int(event_step[i]),
                    x=float(x[i]),
                    y=float(y[i]),
                    z=float(z[i]),
                    strain=int(strain[i]),
                )
            )
    return sources


def release_weight(source: LysisSource, t_s: float, tau_s: float, prune_tau: float = 5.0) -> float:
    """Exponential release-window weight of a source at time ``t_s``."""
    age = t_s - source.time_s
    if age < -1e-12 or age > prune_tau * tau_s + 1e-12:
        return 0.0
    return math.exp(-age / tau_s)


@dataclass(frozen=True)
class Box:
    """Domain extents; x and y are periodic, z is not."""

    Lx: float
    Ly: float
    Lz: float


def radial_support_m(source: Sequence[float], box: Box) -> float:
    """Largest radius whose full sphere around ``source`` is inside the domain."""
    return min(0.5 * box.Lx, 0.5 * box.Ly, source[2], box.Lz - source[2])


def voxel_radius_m(shape: tuple[int, int, int], dx: float, source: Sequence[float], box: Box):
    """Distance from ``source`` to every voxel centre of a ``(nz, ny, nx)`` grid."""
    nz, ny, nx = shape
    xs = (np.arange(nx, dtype=float) + 0.5) * dx
    ys = (np.arange(ny, dtype=float) + 0.5) * dx
    zs = (np.arange(nz, dtype=float) + 0.5) * dx
    dxs = xs - source[0]
    dys = ys - source[1]
    dzs = zs - source[2]
    dxs -= box.Lx * np.round(dxs / box.Lx)
    dys -= box.Ly * np.round(dys / box.Ly)
    return np.sqrt(
        dzs[:, None, None] ** 2 + dys[None, :, None] ** 2 + dxs[None, None, :] ** 2
    )


@dataclass(frozen=True)
class ShellProfile:
    """Fixed-width shell means around one source, inside the radial support."""

    centers_um: np.ndarray
    means: np.ndarray
    counts: np.ndarray
    support_um: float
    r_max_um: float
    shell_width_um: float

    @property
    def complete(self) -> bool:
        """True when every shell in the support contained at least one voxel."""
        return bool(np.all(self.counts > 0))


def shell_profile(
    conc: np.ndarray,
    dx_m: float,
    source: Sequence[float],
    box: Box,
    r_max_um: float,
    shell_width_um: float,
) -> ShellProfile:
    """Shell-mean concentration profile around a single source.

    Raises when ``r_max_um`` exceeds the radial support: a profile that leaves
    the domain would mix transport with boundary truncation, and silently
    trimming it would make the analyzed support depend on where the source
    happens to sit.
    """
    support_um = radial_support_m(source, box) * 1.0e6
    if r_max_um > support_um + 1e-9:
        raise TransportInputError(
            f"requested radial support {r_max_um} um exceeds the in-domain support "
            f"{support_um:.3f} um for source {tuple(source)}"
        )
    n_shells = round(r_max_um / shell_width_um)
    if n_shells < 2 or abs(n_shells * shell_width_um - r_max_um) > 1e-9:
        raise TransportInputError(
            f"radial support {r_max_um} um is not at least two whole shells of "
            f"{shell_width_um} um"
        )

    radius_um = voxel_radius_m(conc.shape, dx_m, source, box) * 1.0e6
    edges = np.linspace(0.0, r_max_um, n_shells + 1)
    means = np.full(n_shells, np.nan, dtype=float)
    counts = np.zeros(n_shells, dtype=int)
    flat_c = conc.ravel()
    flat_r = radius_um.ravel()
    for i in range(n_shells):
        in_shell = (flat_r >= edges[i]) & (flat_r < edges[i + 1])
        counts[i] = int(in_shell.sum())
        if counts[i]:
            means[i] = float(flat_c[in_shell].mean())
    return ShellProfile(
        centers_um=0.5 * (edges[:-1] + edges[1:]),
        means=means,
        counts=counts,
        support_um=support_um,
        r_max_um=r_max_um,
        shell_width_um=shell_width_um,
    )


def profile_radii(profile: ShellProfile, quantiles: Sequence[float] = (0.5, 0.9)) -> list[float]:
    """Radii enclosing ``quantiles`` of the profile's radial mass.

    The radial mass density of a spherically averaged profile is
    ``m(r) * r**2`` (the ``4 * pi * dr`` factor is constant across shells and
    cancels in the normalization), so this is a property of the profile shape
    over a fixed support, independent of how many sources the domain holds.
    """
    if not profile.complete:
        return [float("nan")] * len(quantiles)
    weights = np.maximum(profile.means, 0.0) * profile.centers_um**2
    total = float(weights.sum())
    if total <= 0.0:
        return [float("nan")] * len(quantiles)
    cumulative = np.cumsum(weights) / total
    radii: list[float] = []
    for q in quantiles:
        index = int(np.searchsorted(cumulative, q, side="left"))
        index = min(index, len(cumulative) - 1)
        r_hi = float(profile.centers_um[index])
        f_hi = float(cumulative[index])
        r_lo = float(profile.centers_um[index - 1]) if index else 0.0
        f_lo = float(cumulative[index - 1]) if index else 0.0
        if f_hi <= f_lo:
            radii.append(r_hi)
        else:
            radii.append(r_lo + (q - f_lo) * (r_hi - r_lo) / (f_hi - f_lo))
    return radii


@dataclass(frozen=True)
class PairedTimes:
    """Snapshot times shared by every arm, and what each arm lacks."""

    common_s: list[float]
    missing_s: dict[Any, list[float]]

    @property
    def complete(self) -> bool:
        return not any(self.missing_s.values())


def paired_times(times_by_arm: Mapping[Any, Sequence[float]], tol_s: float = 1e-9) -> PairedTimes:
    """Intersect snapshot times across arms.

    Amplitude arms are compared only at times all of them observed.  A time
    missing from one arm is reported rather than dropped from the union, so a
    caller can block instead of comparing arms over different supports.
    """
    if not times_by_arm:
        return PairedTimes(common_s=[], missing_s={})
    union: list[float] = []
    for times in times_by_arm.values():
        for t in times:
            if not any(abs(t - u) <= tol_s for u in union):
                union.append(float(t))
    union.sort()

    def holds(times: Sequence[float], t: float) -> bool:
        return any(abs(t - other) <= tol_s for other in times)

    common = [t for t in union if all(holds(times, t) for times in times_by_arm.values())]
    missing = {
        arm: [t for t in union if not holds(times, t)] for arm, times in times_by_arm.items()
    }
    return PairedTimes(common_s=common, missing_s=missing)
