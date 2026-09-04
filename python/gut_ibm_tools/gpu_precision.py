"""GPU-vs-CPU scientific outcome comparison.

The comparison answers one question: does running GutIBM on the device change
the scientific answer?  Seed-matched host/device runs of the same
configuration are compared on the estimands the model reports, and the
backend difference is judged against the seed-to-seed spread within a backend.

Paired host/device trajectories diverge once a field difference of order 1e-7
flips a single Bernoulli draw, after which the two runs are different
realizations of the same process.  Divergence of a pair is therefore reported
descriptively and is not a defect.  The defect criterion is distributional: a
backend that shifts a reported median by more than the within-backend seed
spread is not interchangeable with the other.

This module reads simulation HDF5 output and writes comparison records.  It
never starts a simulation.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Any

import h5py

from .path_utils import validate_input_path, write_json_file, write_text_file

SCHEMA_VERSION = 1

#: Relative difference at or below which two values are called identical.
IDENTICAL_RELATIVE_TOLERANCE = 1.0e-12

#: Divergence-onset thresholds for the descriptive paired profile.
DIVERGENCE_THRESHOLDS = (1.0e-12, 1.0e-6, 1.0e-1)
RUN_SPEC_METAVAR = "ARM:SEED:PATH"


@dataclass(frozen=True)
class Estimand:
    """One comparable quantity recorded once per summary step."""

    name: str
    dataset: str
    #: ``stochastic`` quantities are judged distributionally; ``invariant``
    #: quantities are judged absolutely on both backends.
    kind: str
    #: Index into a per-species or per-type vector dataset, when applicable.
    index: int | None = None
    #: Arms whose uptake mode makes the quantity meaningful; empty means all.
    arms: tuple[str, ...] = ()
    #: ``dynamical`` quantities carry state or draw-derived counts, so any
    #: repeat mismatch means the runs are different realizations.  An
    #: ``accounting`` quantity is a diagnostic accumulator that no state
    #: variable reads back, so last-bit repeat noise cannot change the
    #: trajectory.
    repeat_class: str = "dynamical"


STOCHASTIC_ESTIMANDS: tuple[Estimand, ...] = (
    Estimand("n_total", "n_total", "stochastic"),
    Estimand("num_lineages", "num_lineages", "stochastic"),
    Estimand("cumulative_divisions", "events/cumulative_divisions",
             "stochastic"),
    Estimand("cumulative_sos_inductions", "events/cumulative_sos_inductions",
             "stochastic"),
    Estimand("cumulative_mortality_colicin",
             "events/cumulative_mortality_colicin", "stochastic"),
    Estimand("cumulative_mortality_lysis", "events/cumulative_mortality_lysis",
             "stochastic"),
    Estimand("cumulative_mortality_cdi", "events/cumulative_mortality_cdi",
             "stochastic"),
    Estimand("cumulative_conjugation_transfers",
             "events/cumulative_conjugation_transfers", "stochastic"),
    Estimand("cumulative_outflow_washout", "events/cumulative_outflow_washout",
             "stochastic"),
    Estimand("cumulative_outflow_boundary",
             "events/cumulative_outflow_boundary", "stochastic"),
    Estimand("mean_realized_fermentation_fraction",
             "mean_realized_fermentation_fraction", "stochastic"),
    Estimand("max_toxin_BtuB", "chem/max_toxin_BtuB", "stochastic"),
    Estimand("max_toxin_FepA", "chem/max_toxin_FepA", "stochastic"),
    Estimand("max_toxin_CirA", "chem/max_toxin_CirA", "stochastic"),
    Estimand("max_toxin_FhuA", "chem/max_toxin_FhuA", "stochastic"),
    Estimand("mean_carbon", "chem/mean_carbon", "stochastic"),
    Estimand("mean_iron", "chem/mean_iron", "stochastic"),
    Estimand("mean_oxygen", "chem/mean_oxygen", "stochastic"),
    Estimand("bacteriostatic_live_agents",
             "stocks/bacteriostatic_live_agents", "stochastic"),
    Estimand("washout_trapped_live_agents",
             "stocks/washout_trapped_live_agents", "stochastic"),
)

#: Accounting-health channels.  A backend difference in any of these is a
#: defect rather than a difference of realization, because they report the
#: health of the numerics rather than a biological draw.
INVARIANT_ESTIMANDS: tuple[Estimand, ...] = (
    Estimand("halt_reason_code", "halt_reason_code", "invariant"),
    Estimand("reaction_clip_cumulative_total",
             "nutrient_flux/reaction_clip_cumulative", "invariant",
             repeat_class="accounting"),
    Estimand("uptake_shortfall_cumulative_total",
             "nutrient_flux/uptake_shortfall_cumulative", "invariant",
             repeat_class="accounting"),
    Estimand("maintenance_shortfall_cumulative_total",
             "nutrient_flux/maintenance_shortfall_cumulative", "invariant",
             repeat_class="accounting"),
    Estimand("delivery_infeasible_cumulative_total",
             "nutrient_flux/delivery_infeasible_cumulative", "invariant",
             repeat_class="accounting"),
    Estimand("delivery_retry_events_cumulative_total",
             "nutrient_flux/delivery_retry_events_cumulative", "invariant",
             repeat_class="accounting"),
    Estimand("delivery_reduction_cumulative_total",
             "nutrient_flux/delivery_reduction_cumulative", "invariant",
             repeat_class="accounting"),
)

ALL_ESTIMANDS = STOCHASTIC_ESTIMANDS + INVARIANT_ESTIMANDS


class PrecisionInputError(RuntimeError):
    """Raised when a run file cannot supply a comparable trajectory."""


@dataclass
class RunSeries:
    """Per-step trajectories extracted from one run."""

    path: str
    arm: str
    seed: int
    steps: list[int] = field(default_factory=list)
    times: dict[int, float] = field(default_factory=dict)
    series: dict[str, list[float]] = field(default_factory=dict)

    def final(self, name: str) -> float:
        values = self.series[name]
        if not values:
            raise PrecisionInputError(f"empty series for {name} in {self.path}")
        return values[-1]


def _step_groups(handle: h5py.File) -> list[str]:
    if "summary" not in handle:
        raise PrecisionInputError("no summary group; summary output was off")
    group = handle["summary"]
    names = [name for name in group if name.startswith("step_")]
    if not names:
        raise PrecisionInputError("summary group contains no steps")
    return sorted(names, key=lambda name: int(name.removeprefix("step_")))


def _scalar(node: Any) -> float:
    """Reduce a recorded dataset to one comparable number.

    Vector datasets (per species, per receptor) are summed: the comparison
    asks whether the backend moved the reported quantity, and a per-species
    breakdown is recoverable from the raw files when a flagged estimand needs
    attributing.
    """
    value = node[()]
    if hasattr(value, "size"):
        if value.size == 0:
            return 0.0
        flattened = value.reshape(-1)
        if value.size > 1:
            return float(flattened.sum())
        return float(flattened[0])
    return float(value)


def load_run(path: Path, arm: str, seed: int) -> RunSeries:
    """Extract every available estimand trajectory from one run file."""
    resolved = validate_input_path(path)
    run = RunSeries(path=str(resolved), arm=arm, seed=seed)
    with h5py.File(resolved, "r") as handle:
        names = _step_groups(handle)
        run.steps = [int(name.removeprefix("step_")) for name in names]
        for name, step in zip(names, run.steps):
            node = handle[f"summary/{name}"].get("time")
            if node is not None:
                run.times[step] = _scalar(node)
        for estimand in ALL_ESTIMANDS:
            values: list[float] = []
            for name in names:
                node = handle[f"summary/{name}"].get(estimand.dataset)
                if node is None:
                    break
                values.append(_scalar(node))
            if len(values) == len(names):
                run.series[estimand.name] = values
    if not run.series:
        raise PrecisionInputError(f"no comparable estimand in {resolved}")
    return run


def export_run(path: Path, arm: str, seed: int, output: Path) -> None:
    """Export one loaded run as newline-delimited JSON."""
    run = load_run(path, arm, seed)
    header = {
        "schema_version": SCHEMA_VERSION,
        "kind": "header",
        "arm": run.arm,
        "seed": run.seed,
        "source": Path(run.path).name,
        "steps": run.steps,
        "times": [run.times.get(step) for step in run.steps],
    }
    lines = [json.dumps(header)]
    lines.extend(
        json.dumps({"kind": "series", "name": name, "values": values})
        for name, values in run.series.items()
    )
    write_text_file(output, "\n".join(lines) + "\n", allow_external=True)


def load_export(path: Path) -> RunSeries:
    """Load and validate a newline-delimited precision export."""
    resolved = validate_input_path(path)
    try:
        lines = resolved.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as error:
        raise PrecisionInputError(
            f"could not read precision export {resolved}: {error}"
        ) from error
    if not lines:
        raise PrecisionInputError(f"missing header line in {resolved}")
    try:
        header = json.loads(lines[0])
    except json.JSONDecodeError as error:
        raise PrecisionInputError(
            f"unparseable JSON in {resolved} line 1: {error.msg}"
        ) from error
    if not isinstance(header, dict) or header.get("kind") != "header":
        raise PrecisionInputError(f"missing header line in {resolved}")
    schema_version = header.get("schema_version")
    if type(schema_version) is not int or schema_version != SCHEMA_VERSION:
        raise PrecisionInputError(
            f"unexpected schema_version in {resolved}: {schema_version!r}"
        )
    required = {"arm", "seed", "source", "steps", "times"}
    missing = sorted(required - set(header))
    if missing:
        raise PrecisionInputError(
            f"missing header fields in {resolved}: {', '.join(missing)}"
        )
    arm = header["arm"]
    seed = header["seed"]
    steps = header["steps"]
    times = header["times"]
    if not isinstance(arm, str) or type(seed) is not int:
        raise PrecisionInputError(f"invalid arm or seed in {resolved}")
    if (
        not isinstance(steps, list)
        or any(type(step) is not int for step in steps)
        or len(set(steps)) != len(steps)
    ):
        raise PrecisionInputError(f"invalid steps in {resolved}")
    if not isinstance(times, list) or len(times) != len(steps):
        raise PrecisionInputError(
            f"times length differs from steps length in {resolved}"
        )
    run_times: dict[int, float] = {}
    for step, value in zip(steps, times):
        if value is not None:
            if not isinstance(value, (int, float)) or isinstance(value, bool):
                raise PrecisionInputError(
                    f"invalid time for step {step} in {resolved}"
                )
            run_times[step] = float(value)

    series: dict[str, list[float]] = {}
    for line_number, line in enumerate(lines[1:], 2):
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as error:
            raise PrecisionInputError(
                f"unparseable JSON in {resolved} line {line_number}: {error.msg}"
            ) from error
        if (
            not isinstance(payload, dict)
            or payload.get("kind") != "series"
            or not isinstance(payload.get("name"), str)
        ):
            raise PrecisionInputError(
                f"unexpected record in {resolved} line {line_number}"
            )
        name = payload["name"]
        values = payload.get("values")
        if name in series:
            raise PrecisionInputError(
                f"duplicate series {name!r} in {resolved}"
            )
        if not isinstance(values, list) or len(values) != len(steps):
            raise PrecisionInputError(
                f"series {name!r} length differs from steps length in {resolved}"
            )
        if any(
            not isinstance(value, (int, float)) or isinstance(value, bool)
            for value in values
        ):
            raise PrecisionInputError(
                f"invalid values for series {name!r} in {resolved}"
            )
        series[name] = [float(value) for value in values]
    if not series:
        raise PrecisionInputError(f"no series in precision export {resolved}")
    return RunSeries(
        path=str(resolved),
        arm=arm,
        seed=seed,
        steps=steps,
        times=run_times,
        series=series,
    )


def relative_difference(lhs: float, rhs: float) -> float:
    scale = max(1.0, abs(lhs), abs(rhs))
    return abs(lhs - rhs) / scale


def bitwise_equal(lhs: float, rhs: float) -> bool:
    """Compare two doubles by their IEEE-754 representation.

    The reproducibility control needs literal bit equality, not numerical
    closeness, so the comparison is made on the packed bytes.
    """
    return struct.pack("<d", lhs) == struct.pack("<d", rhs)


def ulp_distance(lhs: float, rhs: float) -> int:
    """Return the monotone IEEE-754 representation distance between doubles."""
    lhs_bits = struct.unpack("<q", struct.pack("<d", lhs))[0]
    rhs_bits = struct.unpack("<q", struct.pack("<d", rhs))[0]
    if lhs_bits < 0:
        lhs_bits = 0x8000000000000000 - lhs_bits
    if rhs_bits < 0:
        rhs_bits = 0x8000000000000000 - rhs_bits
    return abs(lhs_bits - rhs_bits)


def repeat_report(first: RunSeries, second: RunSeries) -> dict[str, Any]:
    """Compare two runs that should be bit-identical.

    Used for the same-seed, same-backend reproducibility control.  A failure
    here makes every paired divergence number uninterpretable, so it is
    evaluated before any verdict is formed.
    """
    shared = sorted(set(first.series) & set(second.series))
    repeat_classes = {
        estimand.name: estimand.repeat_class for estimand in ALL_ESTIMANDS
    }
    mismatches: list[dict[str, Any]] = []
    for name in shared:
        lhs = first.series[name]
        rhs = second.series[name]
        limit = min(len(lhs), len(rhs))
        mismatch_steps: list[int] = []
        max_relative = 0.0
        max_ulp = 0
        for index in range(limit):
            if not bitwise_equal(lhs[index], rhs[index]):
                mismatch_steps.append(first.steps[index])
                max_relative = max(
                    max_relative,
                    relative_difference(lhs[index], rhs[index]),
                )
                max_ulp = max(max_ulp, ulp_distance(lhs[index], rhs[index]))
        if mismatch_steps:
            mismatches.append({
                "estimand": name,
                "repeat_class": repeat_classes.get(name, "dynamical"),
                "first_mismatch_step": mismatch_steps[0],
                "mismatch_step_count": len(mismatch_steps),
                "max_relative_difference": max_relative,
                "max_ulp_distance": max_ulp,
            })
    dynamical_mismatches = [
        item for item in mismatches if item["repeat_class"] == "dynamical"
    ]
    accounting_mismatches = [
        item for item in mismatches if item["repeat_class"] == "accounting"
    ]
    return {
        "arm": first.arm,
        "seed": first.seed,
        "compared_estimands": len(shared),
        "step_count": min(len(first.steps), len(second.steps)),
        "reproducible": not mismatches,
        "dynamical_reproducible": not dynamical_mismatches,
        "accounting_reproducible": not accounting_mismatches,
        "mismatches": mismatches,
    }


def _run_label(run: RunSeries) -> str:
    return f"{run.arm}:{run.seed} ({run.path})"


def _comparison_runs(
    host_runs: list[RunSeries],
    device_runs: list[RunSeries],
    repeats: list[tuple[RunSeries, RunSeries]],
) -> list[RunSeries]:
    runs = [*host_runs, *device_runs]
    for first, second in repeats:
        runs.extend((first, second))
    if not runs:
        raise PrecisionInputError("no runs supplied for comparison")
    empty = [run for run in runs if not run.steps]
    if empty:
        labels = ", ".join(_run_label(run) for run in empty)
        raise PrecisionInputError(f"run has no steps: {labels}")
    return runs


def _comparison_metadata(runs: list[RunSeries]) -> dict[str, Any]:
    common_steps = set(runs[0].steps)
    for run in runs[1:]:
        common_steps.intersection_update(run.steps)
    if not common_steps:
        labels = ", ".join(_run_label(run) for run in runs)
        raise PrecisionInputError(
            f"no common comparison step across runs: {labels}"
        )
    comparison_step = max(common_steps)
    last_steps = {
        f"{run.arm}:{run.seed}": max(run.steps)
        for run in runs
    }
    comparison_time = next(
        (
            run.times[comparison_step]
            for run in runs
            if comparison_step in run.times
        ),
        None,
    )
    return {
        "comparison_step": comparison_step,
        "comparison_time_seconds": comparison_time,
        "run_last_steps": last_steps,
        "truncated_runs_present": len(set(last_steps.values())) > 1,
    }


def _value_at_step(run: RunSeries, name: str, step: int) -> float:
    try:
        index = run.steps.index(step)
    except ValueError as error:
        raise PrecisionInputError(
            f"run {_run_label(run)} has no step {step}"
        ) from error
    values = run.series.get(name)
    if values is None or index >= len(values):
        raise PrecisionInputError(
            f"run {_run_label(run)} has no value for {name} at step {step}"
        )
    return values[index]


def divergence_profile(
    host: RunSeries,
    device: RunSeries,
    comparison_step: int | None = None,
) -> dict[str, Any]:
    """First step at which a seed-matched pair separates, per threshold."""
    shared = sorted(set(host.series) & set(device.series))
    shared_steps = sorted(set(host.steps) & set(device.steps))
    if comparison_step is not None:
        shared_steps = [
            step for step in shared_steps if step <= comparison_step
        ]
    onsets: dict[str, dict[str, int | None]] = {}
    for name in shared:
        per_threshold: dict[str, int | None] = {}
        for threshold in DIVERGENCE_THRESHOLDS:
            onset: int | None = None
            for step in shared_steps:
                if relative_difference(
                    _value_at_step(host, name, step),
                    _value_at_step(device, name, step),
                ) > threshold:
                    onset = step
                    break
            per_threshold[f"{threshold:g}"] = onset
        onsets[name] = per_threshold
    return {
        "seed": host.seed,
        "host_arm": host.arm,
        "device_arm": device.arm,
        "compared_steps": len(shared_steps),
        "comparison_step": comparison_step,
        "onset_step": onsets,
    }


def _quartile(values: list[float], upper: bool) -> float:
    ordered = sorted(values)
    count = len(ordered)
    half = count // 2
    part = ordered[count - half:] if upper else ordered[:half]
    return median(part) if part else median(ordered)


def interquartile_range(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return _quartile(values, True) - _quartile(values, False)


def mann_whitney_exact_two_sided(lhs: list[float],
                                 rhs: list[float]) -> dict[str, Any]:
    """Exact two-sided Mann-Whitney U by enumeration.

    Reported as corroboration only.  With five samples per backend the
    smallest attainable two-sided p-value is 2/C(10,5) = 0.0079, so this test
    cannot establish agreement; the effect size against seed spread is the
    primary statistic.
    """
    n_lhs = len(lhs)
    n_rhs = len(rhs)
    if n_lhs == 0 or n_rhs == 0:
        return {"u": None, "p_two_sided": None, "p_floor": None}
    combined = lhs + rhs
    total = len(combined)

    def statistic(selection: tuple[int, ...]) -> float:
        chosen = [combined[index] for index in selection]
        others = [combined[index] for index in range(total)
                  if index not in selection]
        count = 0.0
        for left in chosen:
            for right in others:
                if left > right:
                    count += 1.0
                elif left >= right:
                    count += 0.5
        return count

    observed = statistic(tuple(range(n_lhs)))
    extreme = 0
    considered = 0
    for selection in itertools.combinations(range(total), n_lhs):
        value = statistic(selection)
        considered += 1
        if abs(value - n_lhs * n_rhs / 2.0) >= abs(
                observed - n_lhs * n_rhs / 2.0):
            extreme += 1
    return {
        "u": observed,
        "p_two_sided": extreme / considered,
        "p_floor": 2.0 / considered,
    }


def _accounting_noise_entries(
    repeat_reports: list[dict],
) -> list[dict]:
    entries = []
    for report in repeat_reports:
        for mismatch in report["mismatches"]:
            if mismatch["repeat_class"] == "accounting":
                entries.append({
                    "arm": report["arm"],
                    "seed": report["seed"],
                    **mismatch,
                })
    return entries


def _accounting_noise_floors(
    accounting_noise: list[dict],
) -> dict[str, float]:
    floors: dict[str, float] = {}
    for entry in accounting_noise:
        name = entry["estimand"]
        floors[name] = max(
            floors.get(name, 0.0),
            entry["max_relative_difference"],
        )
    return floors


def _accounting_noise_text(
    accounting_noise: list[dict],
) -> str:
    details = "; ".join(
        f"{entry['estimand']} (max ULP "
        f"{entry['max_ulp_distance']}, max relative difference "
        f"{entry['max_relative_difference']:.6g})"
        for entry in accounting_noise
    )
    return (
        "Accounting repeat noise was observed in: "
        f"{details}. Verdicts stand; these accounting counters are judged "
        "against their measured repeat noise, and a host/device difference "
        "below that floor is not attributable to the backend."
    )


def distributional_verdict(name: str, kind: str, host_values: list[float],
                           device_values: list[float],
                           paired_identical: bool,
                           repeat_noise: float | None = None,
                           ) -> dict[str, Any]:
    """Judge one estimand: backend shift against within-backend seed spread."""
    host_median = median(host_values)
    device_median = median(device_values)
    shift = abs(device_median - host_median)
    spread = max(interquartile_range(host_values),
                 interquartile_range(device_values))
    repeat_noise_floor = max(
        IDENTICAL_RELATIVE_TOLERANCE,
        repeat_noise if repeat_noise is not None else 0.0,
    )
    if paired_identical:
        verdict = "identical"
    elif kind == "invariant":
        verdict = (
            "invariant_match"
            if relative_difference(device_median, host_median)
            <= repeat_noise_floor
            else "invariant_differs"
        )
    elif shift <= spread:
        verdict = "interchangeable"
    else:
        verdict = "flagged"
    return {
        "estimand": name,
        "kind": kind,
        "host_median": host_median,
        "device_median": device_median,
        "host_values": host_values,
        "device_values": device_values,
        "shift": shift,
        "seed_spread": spread,
        "shift_over_spread": (shift / spread) if spread > 0.0 else None,
        "repeat_noise_floor": repeat_noise_floor,
        "verdict": verdict,
        "mann_whitney": mann_whitney_exact_two_sided(host_values,
                                                     device_values),
    }


def compare(host_runs: list[RunSeries], device_runs: list[RunSeries],
            repeats: list[tuple[RunSeries, RunSeries]]) -> dict[str, Any]:
    """Full comparison record for one uptake mode."""
    runs = _comparison_runs(host_runs, device_runs, repeats)
    metadata = _comparison_metadata(runs)
    comparison_step = metadata["comparison_step"]
    repeat_reports = [repeat_report(first, second) for first, second in repeats]
    reproducible = all(report["reproducible"] for report in repeat_reports)
    dynamical_reproducible = all(
        report["dynamical_reproducible"] for report in repeat_reports
    )
    accounting_reproducible = all(
        report["accounting_reproducible"] for report in repeat_reports
    )
    accounting_noise = _accounting_noise_entries(repeat_reports)
    accounting_noise_floors = _accounting_noise_floors(accounting_noise)

    host_by_seed = {run.seed: run for run in host_runs}
    device_by_seed = {run.seed: run for run in device_runs}
    paired_seeds = sorted(set(host_by_seed) & set(device_by_seed))
    profiles = [
        divergence_profile(
            host_by_seed[seed], device_by_seed[seed], comparison_step
        )
        for seed in paired_seeds
    ] if dynamical_reproducible else []

    verdicts: list[dict[str, Any]] = []
    available = sorted(
        set.intersection(*[set(run.series) for run in host_runs + device_runs])
    ) if host_runs and device_runs else []
    by_name = {estimand.name: estimand for estimand in ALL_ESTIMANDS}
    if dynamical_reproducible:
        for name in available:
            host_values = [
                _value_at_step(run, name, comparison_step)
                for run in host_runs
            ]
            device_values = [
                _value_at_step(run, name, comparison_step)
                for run in device_runs
            ]
            noise_floor = accounting_noise_floors.get(name)
            paired_tolerance = max(
                IDENTICAL_RELATIVE_TOLERANCE,
                noise_floor if noise_floor is not None else 0.0,
            )
            paired_identical = all(
                relative_difference(
                    _value_at_step(
                        host_by_seed[seed], name, comparison_step
                    ),
                    _value_at_step(
                        device_by_seed[seed], name, comparison_step
                    ),
                )
                <= paired_tolerance
                for seed in paired_seeds
            ) if paired_seeds else False
            verdicts.append(distributional_verdict(
                name, by_name[name].kind, host_values, device_values,
                paired_identical, noise_floor))

    non_finite = [
        {"path": run.path, "estimand": name}
        for run in host_runs + device_runs
        for name, values in run.series.items()
        if any(not math.isfinite(value) for value in values)
    ]
    # A non-finite value anywhere in a trajectory is a defect even when the
    # final step it is judged on happens to be finite again, so it forces the
    # estimand into the flagged set rather than resting on its verdict.
    non_finite_estimands = {item["estimand"] for item in non_finite}
    for item in verdicts:
        if item["estimand"] in non_finite_estimands:
            item["verdict"] = "non_finite"
    flagged = [item["estimand"] for item in verdicts
               if item["verdict"] in {"flagged", "invariant_differs",
                                      "non_finite"}]
    flagged.extend(sorted(non_finite_estimands - set(flagged)))
    if not dynamical_reproducible:
        interpretation = (
            "Reproducibility control FAILED: the device is not bit-identical "
            "across repeats of one seed at this scale, so paired divergence "
            "cannot be attributed to the backend and no distributional "
            "verdict should be quoted."
        )
    elif accounting_noise:
        interpretation = _accounting_noise_text(accounting_noise)
    else:
        interpretation = (
            "Paired divergence is descriptive: a stochastic nonlinear model "
            "separates once a field difference flips one draw. The verdict is "
            "distributional; 'interchangeable' means the backend shift in the "
            "median is within the within-backend seed spread, which is "
            "reported alongside so the resolution of the measurement is "
            "visible."
        )
    record = {
        "schema_version": SCHEMA_VERSION,
        "host_arm": host_runs[0].arm if host_runs else None,
        "device_arm": device_runs[0].arm if device_runs else None,
        "host_seeds": sorted(host_by_seed),
        "device_seeds": sorted(device_by_seed),
        "paired_seeds": paired_seeds,
        "reproducibility": {
            "reproducible": reproducible,
            "dynamical_reproducible": dynamical_reproducible,
            "accounting_reproducible": accounting_reproducible,
            "reports": repeat_reports,
        },
        "dynamical_reproducible": dynamical_reproducible,
        "accounting_reproducible": accounting_reproducible,
        "accounting_repeat_noise": accounting_noise,
        "divergence_profiles": profiles,
        "verdicts": verdicts,
        "flagged_estimands": flagged,
        "non_finite": non_finite,
        "interpretation": interpretation,
    }
    record.update(metadata)
    if metadata["truncated_runs_present"]:
        record["horizon_warning"] = (
            "WARNING: supplied runs have different trajectory depths; "
            f"all verdicts use common step {comparison_step}, not each "
            "run's own final step."
        )
    else:
        record["horizon_warning"] = None
    return record


def _verdict_table(record: dict[str, Any]) -> list[str]:
    lines = [
        ("| estimand | kind | host median | device median | shift |"
         " seed spread | shift/spread | verdict |"),
        "|---|---|---|---|---|---|---|---|",
    ]
    for item in record["verdicts"]:
        ratio = item["shift_over_spread"]
        lines.append(
            f"| `{item['estimand']}` | {item['kind']} "
            f"| {item['host_median']:.6g} | {item['device_median']:.6g} "
            f"| {item['shift']:.6g} | {item['seed_spread']:.6g} "
            f"| {'—' if ratio is None else f'{ratio:.3g}'} "
            f"| {item['verdict']} |"
        )
    return lines


def _repeat_table(record: dict) -> list[str]:
    mismatches = [
        mismatch
        for report in record["reproducibility"]["reports"]
        for mismatch in report["mismatches"]
    ]
    lines = [
        (
            "| estimand | class | mismatch steps | max ULP |"
            " max relative difference |"
        ),
        "|---|---|---:|---:|---:|",
    ]
    lines.extend(
        f"| `{item['estimand']}` | {item['repeat_class']} "
        f"| {item['mismatch_step_count']} | {item['max_ulp_distance']} "
        f"| {item['max_relative_difference']:.6g} |"
        for item in mismatches
    )
    return lines


def render_markdown(record: dict[str, Any]) -> str:
    lines = [
        (f"# GPU-vs-CPU outcome comparison: {record['host_arm']} vs "
         f"{record['device_arm']}"),
        "",
        f"Paired seeds: {record['paired_seeds']}",
        (f"Comparison step: {record['comparison_step']} "
         f"(time={record['comparison_time_seconds']})"),
        "",
    ]
    if record["truncated_runs_present"]:
        lines.extend([
            "",
            "## WARNING: TRUNCATED RUNS PRESENT",
            "",
            f"**{record['horizon_warning']}**",
            "",
            f"Run last steps: `{record['run_last_steps']}`",
            "",
        ])
    lines.extend([
        "## Reproducibility control",
        "",
        ("Device repeats bit-identical: "
         f"{record['reproducibility']['reproducible']}"),
        ("Dynamical reproducibility: "
         f"{record['dynamical_reproducible']}"),
        ("Accounting reproducibility: "
         f"{record['accounting_reproducible']}"),
        "",
        "### Repeat mismatch details",
        "",
    ])
    mismatch_rows = _repeat_table(record)
    if len(mismatch_rows) == 2:
        lines.append("No repeat mismatches recorded.")
    else:
        lines.extend(mismatch_rows)
    if record["accounting_repeat_noise"]:
        lines.extend([
            "",
            "### Accounting repeat-noise caveat",
            "",
            record["interpretation"],
        ])
    lines.extend([
        "",
        "## Distributional verdict",
        "",
    ])
    lines.extend(_verdict_table(record))
    lines.extend([
        "",
        "## Interpretation",
        "",
        record["interpretation"],
        "",
    ])
    if record["flagged_estimands"]:
        lines.extend([
            "## Flagged",
            "",
            ("These estimands moved by more than the within-backend seed "
             "spread, or differ on an invariant channel:"),
            "",
        ])
        lines.extend(f"- `{name}`" for name in record["flagged_estimands"])
        lines.append("")
    if record["non_finite"]:
        lines.extend(["## Non-finite values", ""])
        lines.extend(
            f"- `{item['estimand']}` in `{item['path']}`"
            for item in record["non_finite"]
        )
        lines.append("")
    return "\n".join(lines)


def _parse_run_spec(spec: str) -> tuple[str, int, Path]:
    parts = spec.split(":", 2)
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            f"run spec must be {RUN_SPEC_METAVAR}, got {spec}")
    arm, seed, path = parts
    try:
        seed_value = int(seed)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"seed must be an integer in {spec}") from error
    return arm, seed_value, Path(path)


def _load_spec_run(path: Path, arm: str, seed: int) -> RunSeries:
    if path.suffix.lower() in {".json", ".jsonl"}:
        run = load_export(path)
        if run.arm != arm or run.seed != seed:
            raise PrecisionInputError(
                f"run spec identity {arm}:{seed} does not match "
                f"export {run.arm}:{run.seed} in {path}"
            )
        return run
    return load_run(path, arm, seed)


def _export_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Export one HDF5 run for CloudWatch transport."
    )
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    export_run(args.hdf5, args.arm, args.seed, args.output)
    return 0


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == "export":
        return _export_main(arguments[1:])
    parser = argparse.ArgumentParser(
        description="Compare host and device runs on reported estimands.")
    parser.add_argument("--host", action="append", default=[],
                        type=_parse_run_spec, metavar=RUN_SPEC_METAVAR)
    parser.add_argument("--device", action="append", default=[],
                        type=_parse_run_spec, metavar=RUN_SPEC_METAVAR)
    parser.add_argument("--repeat", action="append", default=[],
                        type=_parse_run_spec, metavar=RUN_SPEC_METAVAR,
                        help="device repeat of an already supplied seed")
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--markdown-out", type=Path)
    args = parser.parse_args(argv)

    if not args.host or not args.device:
        parser.error("at least one --host and one --device run are required")

    host_runs = [
        _load_spec_run(path, arm, seed) for arm, seed, path in args.host
    ]
    device_runs = [
        _load_spec_run(path, arm, seed) for arm, seed, path in args.device
    ]
    device_by_seed = {run.seed: run for run in device_runs}
    repeats: list[tuple[RunSeries, RunSeries]] = []
    for arm, seed, path in args.repeat:
        original = device_by_seed.get(seed)
        if original is None:
            parser.error(f"--repeat seed {seed} has no matching --device run")
        repeats.append((original, _load_spec_run(path, arm, seed)))

    record = compare(host_runs, device_runs, repeats)
    write_json_file(args.json_out, record, indent=2, allow_external=True)
    if args.markdown_out:
        write_text_file(
            args.markdown_out,
            render_markdown(record),
            allow_external=True,
        )
    if record["truncated_runs_present"]:
        print(record["horizon_warning"])
    print(json.dumps({
        "comparison_step": record["comparison_step"],
        "comparison_time_seconds": record["comparison_time_seconds"],
        "run_last_steps": record["run_last_steps"],
        "truncated_runs_present": record["truncated_runs_present"],
        "reproducible": record["reproducibility"]["reproducible"],
        "dynamical_reproducible": record["dynamical_reproducible"],
        "accounting_reproducible": record["accounting_reproducible"],
        "accounting_repeat_noise": record["accounting_repeat_noise"],
        "flagged": record["flagged_estimands"],
        "paired_seeds": record["paired_seeds"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
