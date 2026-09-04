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
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Any

import h5py

from .path_utils import validate_input_path, write_json_file

SCHEMA_VERSION = 1

#: Relative difference at or below which two values are called identical.
IDENTICAL_RELATIVE_TOLERANCE = 1.0e-12

#: Divergence-onset thresholds for the descriptive paired profile.
DIVERGENCE_THRESHOLDS = (1.0e-12, 1.0e-6, 1.0e-1)


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
             "nutrient_flux/reaction_clip_cumulative", "invariant"),
    Estimand("uptake_shortfall_cumulative_total",
             "nutrient_flux/uptake_shortfall_cumulative", "invariant"),
    Estimand("maintenance_shortfall_cumulative_total",
             "nutrient_flux/maintenance_shortfall_cumulative", "invariant"),
    Estimand("delivery_infeasible_cumulative_total",
             "nutrient_flux/delivery_infeasible_cumulative", "invariant"),
    Estimand("delivery_retry_events_cumulative_total",
             "nutrient_flux/delivery_retry_events_cumulative", "invariant"),
    Estimand("delivery_reduction_cumulative_total",
             "nutrient_flux/delivery_reduction_cumulative", "invariant"),
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


def relative_difference(lhs: float, rhs: float) -> float:
    scale = max(1.0, abs(lhs), abs(rhs))
    return abs(lhs - rhs) / scale


def bitwise_equal(lhs: float, rhs: float) -> bool:
    """Compare two doubles by their IEEE-754 representation.

    The reproducibility control needs literal bit equality, not numerical
    closeness, so the comparison is made on the packed bytes.
    """
    return struct.pack("<d", lhs) == struct.pack("<d", rhs)


def repeat_report(first: RunSeries, second: RunSeries) -> dict[str, Any]:
    """Compare two runs that should be bit-identical.

    Used for the same-seed, same-backend reproducibility control.  A failure
    here makes every paired divergence number uninterpretable, so it is
    evaluated before any verdict is formed.
    """
    shared = sorted(set(first.series) & set(second.series))
    mismatches: list[dict[str, Any]] = []
    for name in shared:
        lhs = first.series[name]
        rhs = second.series[name]
        limit = min(len(lhs), len(rhs))
        for index in range(limit):
            if not bitwise_equal(lhs[index], rhs[index]):
                mismatches.append({
                    "estimand": name,
                    "step": first.steps[index],
                    "first": lhs[index],
                    "second": rhs[index],
                    "relative_difference": relative_difference(
                        lhs[index], rhs[index]),
                })
                break
    return {
        "arm": first.arm,
        "seed": first.seed,
        "compared_estimands": len(shared),
        "step_count": min(len(first.steps), len(second.steps)),
        "reproducible": not mismatches,
        "mismatches": mismatches,
    }


def divergence_profile(host: RunSeries, device: RunSeries) -> dict[str, Any]:
    """First step at which a seed-matched pair separates, per threshold."""
    shared = sorted(set(host.series) & set(device.series))
    onsets: dict[str, dict[str, int | None]] = {}
    for name in shared:
        lhs = host.series[name]
        rhs = device.series[name]
        limit = min(len(lhs), len(rhs))
        per_threshold: dict[str, int | None] = {}
        for threshold in DIVERGENCE_THRESHOLDS:
            onset: int | None = None
            for index in range(limit):
                if relative_difference(lhs[index], rhs[index]) > threshold:
                    onset = host.steps[index]
                    break
            per_threshold[f"{threshold:g}"] = onset
        onsets[name] = per_threshold
    return {
        "seed": host.seed,
        "host_arm": host.arm,
        "device_arm": device.arm,
        "compared_steps": min(len(host.steps), len(device.steps)),
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
                elif not (left < right):
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


def distributional_verdict(name: str, kind: str, host_values: list[float],
                           device_values: list[float],
                           paired_identical: bool) -> dict[str, Any]:
    """Judge one estimand: backend shift against within-backend seed spread."""
    host_median = median(host_values)
    device_median = median(device_values)
    shift = abs(device_median - host_median)
    spread = max(interquartile_range(host_values),
                 interquartile_range(device_values))
    if paired_identical:
        verdict = "identical"
    elif kind == "invariant":
        verdict = "invariant_match" if shift <= 0.0 else "invariant_differs"
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
        "verdict": verdict,
        "mann_whitney": mann_whitney_exact_two_sided(host_values,
                                                     device_values),
    }


def compare(host_runs: list[RunSeries], device_runs: list[RunSeries],
            repeats: list[tuple[RunSeries, RunSeries]]) -> dict[str, Any]:
    """Full comparison record for one uptake mode."""
    repeat_reports = [repeat_report(first, second) for first, second in repeats]
    reproducible = all(report["reproducible"] for report in repeat_reports)

    host_by_seed = {run.seed: run for run in host_runs}
    device_by_seed = {run.seed: run for run in device_runs}
    paired_seeds = sorted(set(host_by_seed) & set(device_by_seed))
    profiles = [
        divergence_profile(host_by_seed[seed], device_by_seed[seed])
        for seed in paired_seeds
    ]

    verdicts: list[dict[str, Any]] = []
    available = sorted(
        set.intersection(*[set(run.series) for run in host_runs + device_runs])
    ) if host_runs and device_runs else []
    by_name = {estimand.name: estimand for estimand in ALL_ESTIMANDS}
    for name in available:
        host_values = [run.final(name) for run in host_runs]
        device_values = [run.final(name) for run in device_runs]
        paired_identical = all(
            relative_difference(host_by_seed[seed].final(name),
                                device_by_seed[seed].final(name))
            <= IDENTICAL_RELATIVE_TOLERANCE
            for seed in paired_seeds
        ) if paired_seeds else False
        verdicts.append(distributional_verdict(
            name, by_name[name].kind, host_values, device_values,
            paired_identical))

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
    return {
        "schema_version": SCHEMA_VERSION,
        "host_arm": host_runs[0].arm if host_runs else None,
        "device_arm": device_runs[0].arm if device_runs else None,
        "host_seeds": sorted(host_by_seed),
        "device_seeds": sorted(device_by_seed),
        "paired_seeds": paired_seeds,
        "reproducibility": {
            "reproducible": reproducible,
            "reports": repeat_reports,
        },
        "divergence_profiles": profiles,
        "verdicts": verdicts,
        "flagged_estimands": flagged,
        "non_finite": non_finite,
        "interpretation": (
            "Paired divergence is descriptive: a stochastic nonlinear model "
            "separates once a field difference flips one draw. The verdict is "
            "distributional; 'interchangeable' means the backend shift in the "
            "median is within the within-backend seed spread, which is "
            "reported alongside so the resolution of the measurement is "
            "visible."
            if reproducible else
            "Reproducibility control FAILED: the device is not bit-identical "
            "across repeats of one seed at this scale, so paired divergence "
            "cannot be attributed to the backend and no distributional "
            "verdict should be quoted."
        ),
    }


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


def render_markdown(record: dict[str, Any]) -> str:
    lines = [
        (f"# GPU-vs-CPU outcome comparison: {record['host_arm']} vs "
         f"{record['device_arm']}"),
        "",
        f"Paired seeds: {record['paired_seeds']}",
        "",
        "## Reproducibility control",
        "",
        ("Device repeats bit-identical: "
         f"{record['reproducibility']['reproducible']}"),
        "",
        "## Distributional verdict",
        "",
    ]
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
            "run spec must be ARM:SEED:PATH, got " + spec)
    arm, seed, path = parts
    try:
        seed_value = int(seed)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"seed must be an integer in {spec}") from error
    return arm, seed_value, Path(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare host and device runs on reported estimands.")
    parser.add_argument("--host", action="append", default=[],
                        type=_parse_run_spec, metavar="ARM:SEED:PATH")
    parser.add_argument("--device", action="append", default=[],
                        type=_parse_run_spec, metavar="ARM:SEED:PATH")
    parser.add_argument("--repeat", action="append", default=[],
                        type=_parse_run_spec, metavar="ARM:SEED:PATH",
                        help="device repeat of an already supplied seed")
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--markdown-out", type=Path)
    args = parser.parse_args(argv)

    if not args.host or not args.device:
        parser.error("at least one --host and one --device run are required")

    host_runs = [load_run(path, arm, seed) for arm, seed, path in args.host]
    device_runs = [load_run(path, arm, seed) for arm, seed, path in args.device]
    device_by_seed = {run.seed: run for run in device_runs}
    repeats: list[tuple[RunSeries, RunSeries]] = []
    for arm, seed, path in args.repeat:
        original = device_by_seed.get(seed)
        if original is None:
            parser.error(f"--repeat seed {seed} has no matching --device run")
        repeats.append((original, load_run(path, arm, seed)))

    record = compare(host_runs, device_runs, repeats)
    write_json_file(args.json_out, record, indent=2, allow_external=True)
    if args.markdown_out:
        args.markdown_out.write_text(render_markdown(record), encoding="utf-8")
    print(json.dumps({
        "reproducible": record["reproducibility"]["reproducible"],
        "flagged": record["flagged_estimands"],
        "paired_seeds": record["paired_seeds"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
