#!/usr/bin/env python3
"""Analyze the single-source transport assay against its predeclared contract.

Every definition this program applies — shell averaging, radial support,
boundary handling, r50/r90, the analysis window, the paired-time rule, and the
gate criteria — is read from ``assay_contract.json``, which is committed before
the runs exist.  Nothing here is chosen after seeing the results.

Exit status: 0 only when the gate passes; 2 when the assay is blocked (missing
outputs, failed authentication, missing paired times), 1 when it ran and the
predeclared ordering failed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import shutil
import sys
import tempfile
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from itertools import pairwise
from pathlib import Path

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover
    raise SystemExit(f"numpy required: {exc}") from exc

try:
    import h5py
except ImportError as exc:  # pragma: no cover
    raise SystemExit(f"h5py required: {exc}") from exc

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO / "python"))
from gut_ibm_tools.path_utils import prepare_output_directory
from gut_ibm_tools.transport_metrics import (
    Box,
    ExpectedRun,
    TransportInputError,
    authenticate_run,
    paired_times,
    profile_radii,
    radial_support_m,
    read_lysis_sources,
    shell_profile,
    sorted_steps,
)

CONTRACT = json.loads((ROOT / "assay_contract.json").read_text())
DECL = CONTRACT["predeclared_analysis"]
GATE = CONTRACT["gate"]
R_MAX_UM = float(DECL["radial_support_um"])
SHELL_UM = float(DECL["shell_width_um"])
NEAR_FIELD_UM = tuple(DECL["near_field_um"])
WINDOW_S = 1500.0
BURST_TAU_S = 300.0
ORDER_MARGIN_R50_UM = float(GATE["order_margin_r50_um"])
ORDER_MARGIN_R90_UM = float(GATE["order_margin_r90_um"])
NULL_ABS_TOL = float(GATE["null_abs_tol"])
POSITION_TOL_M = float(GATE["position_tol_m"])
EVENT_TIME_TOL_S = float(GATE["event_time_tol_s"])
MONOTONIC_ABS_TOL = float(GATE["monotonic_abs_tol"])
AMPLITUDES = [float(v) for v in CONTRACT["design"]["axes"]["bacteriocin.mucin_charge.amplitude"]]
ENDPOINT_AMPLITUDES = [float(v) for v in CONTRACT["design"]["endpoint_amplitudes"]]
STATUS_BLOCKED = "BLOCKED"
STATUS_FAIL = "FAIL"
STATUS_PASS = "PASS"


@contextmanager
def open_h5(path: Path) -> Iterator[h5py.File]:
    if path.suffix != ".gz":
        with h5py.File(path, "r") as h5:
            yield h5
        return
    with tempfile.TemporaryDirectory() as tmpdir:
        raw = Path(tmpdir) / "run.h5"
        with gzip.open(path, "rb") as src, raw.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        with h5py.File(raw, "r") as h5:
            yield h5


def find_output(results_root: Path, index: int) -> Path | None:
    for name in ("output.h5.gz", "output.h5"):
        candidate = results_root / str(index) / name
        if candidate.exists():
            return candidate
    return None


def snapshot_time(h5, step_key: str, bio_dt: float) -> float:
    summary = h5.get("summary")
    if summary is not None and step_key in summary and "time" in summary[step_key]:
        return float(np.asarray(summary[step_key]["time"][()]).reshape(-1)[0])
    return float(int(step_key.rsplit("_", 1)[1]) * bio_dt)


def analyze_producer(run: dict, path: Path, execution_sha: str) -> dict:
    """Radial profiles around the run's single release source."""
    report = {
        "run_id": run["run_id"],
        "seed": run["seed"],
        "amplitude": run["amplitude"],
        "array_index": run["array_index"],
        "arm": run["arm"],
        "violations": [],
        "snapshots": [],
    }
    with open_h5(path) as h5:
        expected = ExpectedRun(
            execution_source_sha=execution_sha,
            seed=int(run["seed"]),
            amplitude=float(run["amplitude"]),
            hdf5_schedule=run["hdf5_schedule"],
            chemistry_placement=CONTRACT["runtime"]["chemistry_placement"],
        )
        report["violations"] = authenticate_run(h5, expected)

        sources = read_lysis_sources(h5)
        report["n_lysis_events"] = len(sources)
        if len(sources) != 1:
            report["violations"].append(
                f"expected exactly one producer lysis source, found {len(sources)}"
            )
            return report
        source = sources[0]
        box = Box(Lx=run["Lx"], Ly=run["Ly"], Lz=run["Lz"])
        support_um = radial_support_m(source.position, box) * 1.0e6
        report["source"] = {
            "event_time_s": source.time_s,
            "event_step": source.step,
            "x": source.x,
            "y": source.y,
            "z": source.z,
            "radial_support_um": support_um,
        }
        if support_um < R_MAX_UM - 1e-9:
            report["violations"].append(
                f"source radial support {support_um:.3f} um < predeclared {R_MAX_UM} um"
            )
            return report

        for step_key in sorted_steps(h5["grid"]):
            t_s = snapshot_time(h5, step_key, run["bio_dt"])
            if t_s <= source.time_s + 1e-9 or t_s > source.time_s + WINDOW_S + 1e-9:
                continue
            conc = np.asarray(h5["grid"][step_key]["bacteriocin_BtuB"][()], dtype=float)
            profile = shell_profile(
                conc, run["dx"], source.position, box, R_MAX_UM, SHELL_UM
            )
            r50, r90 = profile_radii(profile)
            near = profile.means[
                (profile.centers_um >= NEAR_FIELD_UM[0])
                & (profile.centers_um < NEAR_FIELD_UM[1])
            ]
            report["snapshots"].append(
                {
                    "run_id": run["run_id"],
                    "seed": run["seed"],
                    "amplitude": run["amplitude"],
                    "t_s": t_s,
                    "source_age_s": t_s - source.time_s,
                    "release_weight": math.exp(-(t_s - source.time_s) / BURST_TAU_S),
                    "r50_um": r50,
                    "r90_um": r90,
                    "near_field_0_10_mean": float(np.mean(near)),
                    "profile_complete": profile.complete,
                    "shell_means": [float(v) for v in profile.means],
                }
            )
        return report


def analyze_null(run: dict, path: Path, execution_sha: str) -> dict:
    report = {
        "run_id": run["run_id"],
        "seed": run["seed"],
        "arm": run["arm"],
        "violations": [],
    }
    with open_h5(path) as h5:
        report["violations"] = authenticate_run(
            h5,
            ExpectedRun(
                execution_source_sha=execution_sha,
                seed=int(run["seed"]),
                amplitude=float(run["amplitude"]),
                hdf5_schedule=run["hdf5_schedule"],
                chemistry_placement=CONTRACT["runtime"]["chemistry_placement"],
            ),
        )
        max_abs = 0.0
        n_snapshots = 0
        for step_key in sorted_steps(h5["grid"]):
            conc = np.asarray(h5["grid"][step_key]["bacteriocin_BtuB"][()], dtype=float)
            max_abs = max(max_abs, float(np.max(np.abs(conc))))
            n_snapshots += 1
        report["n_grid_snapshots"] = n_snapshots
        report["max_abs_bacteriocin_BtuB"] = max_abs
        report["toxin_free"] = max_abs <= NULL_ABS_TOL
        report["n_lysis_events"] = len(read_lysis_sources(h5))
        return report


def load_runs(results_root: Path) -> tuple[list[dict], str]:
    manifest = json.loads((ROOT / "generated" / "manifest.json").read_text())
    runs = []
    for entry in manifest["runs"]:
        cfg = json.loads((ROOT / entry["input_relpath"]).read_text())
        runs.append(
            {
                **entry,
                "bio_dt": float(cfg["bio_dt"]),
                "dx": float(cfg["grid_dx"]),
                "Lx": float(cfg["domain_x"]),
                "Ly": float(cfg["domain_y"]),
                "Lz": float(cfg["domain_z"]),
                "hdf5_schedule": {
                    key: int(value)
                    for key, value in cfg["hdf5"]["schedule"].items()
                    if key != "grid_species"
                },
                "output_path": find_output(results_root, entry["array_index"]),
            }
        )
    return runs, manifest["execution_source_sha"]


def endpoint_decreasing(values_by_amp: dict[float, float], margin: float) -> bool:
    values = [values_by_amp[a] for a in ENDPOINT_AMPLITUDES]
    return all(
        np.isfinite(a) and np.isfinite(b) and a - b >= margin
        for a, b in pairwise(values)
    )


def endpoint_increasing(values_by_amp: dict[float, float]) -> bool:
    values = [values_by_amp[a] for a in ENDPOINT_AMPLITUDES]
    return all(np.isfinite(a) and np.isfinite(b) and b > a for a, b in pairwise(values))


def nonincreasing(values: list[float], tolerance: float) -> bool:
    return all(
        np.isfinite(a) and np.isfinite(b) and b <= a + tolerance
        for a, b in pairwise(values)
    )


def nondecreasing(values: list[float], tolerance: float) -> bool:
    return all(
        np.isfinite(a) and np.isfinite(b) and b + tolerance >= a
        for a, b in pairwise(values)
    )


def order_checks(producer_reports: list[dict]) -> tuple[dict, list[str]]:
    """Endpoint margins plus interpolation monotonicity at paired times, per seed."""
    by_seed_amp: dict[int, dict[float, dict[float, dict]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    for report in producer_reports:
        for snap in report["snapshots"]:
            by_seed_amp[int(snap["seed"])][float(snap["amplitude"])][
                float(snap["t_s"])
            ] = snap

    detail: dict[str, dict] = {}
    blockers: list[str] = []
    expected = set(AMPLITUDES)
    for seed in CONTRACT["design"]["seeds"]:
        by_amp = by_seed_amp.get(int(seed), {})
        observed = set(by_amp)
        pairs = paired_times({amp: sorted(by_amp.get(amp, {})) for amp in AMPLITUDES})
        seed_detail = {
            "paired_times_s": pairs.common_s,
            "missing_times_s": {str(amp): times for amp, times in pairs.missing_s.items()},
            "n_amplitude_arms": len(observed),
            "missing_amplitudes": sorted(expected - observed),
            "times": [],
        }
        if observed != expected:
            blockers.append(
                f"seed {seed}: amplitude arms {sorted(observed)}, expected {AMPLITUDES}"
            )
        if not pairs.complete:
            blockers.append(
                f"seed {seed}: snapshot times are not shared by all amplitudes ({pairs.missing_s})"
            )
        if not pairs.common_s:
            blockers.append(f"seed {seed}: no snapshot time is shared by all amplitudes")

        for t_s in pairs.common_s:
            snaps = [by_amp[amp][t_s] for amp in AMPLITUDES]
            r50 = [snap["r50_um"] for snap in snaps]
            r90 = [snap["r90_um"] for snap in snaps]
            near = [snap["near_field_0_10_mean"] for snap in snaps]
            r50_by_amp = dict(zip(AMPLITUDES, r50))
            r90_by_amp = dict(zip(AMPLITUDES, r90))
            near_by_amp = dict(zip(AMPLITUDES, near))
            finite = all(np.isfinite(value) for value in (*r50, *r90, *near)) and all(
                snap["profile_complete"] for snap in snaps
            )
            seed_detail["times"].append({
                "t_s": t_s,
                "amplitudes": AMPLITUDES,
                "r50_um": r50,
                "r90_um": r90,
                "near_field_0_10_mean": near,
                "finite_and_complete": finite,
                "endpoint_r50_ordered_0_gt_15_gt_60": endpoint_decreasing(r50_by_amp, ORDER_MARGIN_R50_UM),
                "endpoint_r90_ordered_0_gt_15_gt_60": endpoint_decreasing(r90_by_amp, ORDER_MARGIN_R90_UM),
                "endpoint_near_field_ordered_0_lt_15_lt_60": endpoint_increasing(near_by_amp),
                "interpolation_r50_nonincreasing": nonincreasing(r50, MONOTONIC_ABS_TOL),
                "interpolation_r90_nonincreasing": nonincreasing(r90, MONOTONIC_ABS_TOL),
                "interpolation_near_field_nondecreasing": nondecreasing(near, MONOTONIC_ABS_TOL),
            })
        detail[str(seed)] = seed_detail
    return detail, blockers


def position_checks(producer_reports: list[dict]) -> tuple[dict, list[str]]:
    """All five amplitude arms of a seed must share source position and time."""
    by_seed: dict[int, list[dict]] = defaultdict(list)
    for report in producer_reports:
        if "source" in report:
            by_seed[int(report["seed"])].append(report)
    detail: dict[str, dict] = {}
    blockers: list[str] = []
    for seed in CONTRACT["design"]["seeds"]:
        reports = sorted(by_seed.get(int(seed), []), key=lambda r: float(r["amplitude"]))
        amplitudes = [float(r["amplitude"]) for r in reports]
        positions = [(r["source"]["x"], r["source"]["y"], r["source"]["z"]) for r in reports]
        times = [float(r["source"]["event_time_s"]) for r in reports]
        if positions:
            spread = [max(abs(p[axis] - positions[0][axis]) for p in positions) for axis in range(3)]
            max_position_spread = max(spread)
            max_time_spread = max(abs(t - times[0]) for t in times)
        else:
            max_position_spread = math.inf
            max_time_spread = math.inf
        complete = amplitudes == AMPLITUDES
        identical_position = complete and max_position_spread <= POSITION_TOL_M
        identical_time = complete and max_time_spread <= EVENT_TIME_TOL_S
        detail[str(seed)] = {
            "n_arms": len(reports), "amplitudes": amplitudes, "positions_m": positions,
            "max_spread_m": max_position_spread, "identical_across_amplitudes": identical_position,
            "event_times_s": times, "max_event_time_spread_s": max_time_spread,
            "event_time_identical_across_amplitudes": identical_time,
        }
        if not complete:
            blockers.append(f"seed {seed}: source records for amplitudes {amplitudes}, expected {AMPLITUDES}")
        if complete and not identical_position:
            blockers.append(f"seed {seed}: source position differs across amplitudes by {max_position_spread:.3e} m")
        if complete and not identical_time:
            blockers.append(f"seed {seed}: source event time differs across amplitudes by {max_time_spread:.3e} s")
    return detail, blockers


def write_snapshot_csv(path: Path, producer_reports: list[dict]) -> int:
    rows = []
    for report in producer_reports:
        for snap in report["snapshots"]:
            row = {key: value for key, value in snap.items() if key != "shell_means"}
            rows.append(row)
    if not rows:
        path.write_text("")
        return 0
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def jsonable(value):
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root", type=Path, default=ROOT / "generated" / "results"
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "analysis")
    args = parser.parse_args()
    output_dir = prepare_output_directory(args.output_dir)

    runs, execution_sha = load_runs(args.results_root)
    blockers: list[str] = []
    producer_reports: list[dict] = []
    null_reports: list[dict] = []

    for run in runs:
        if run["output_path"] is None:
            blockers.append(f"{run['run_id']}: no returned output file")
            continue
        try:
            if run["arm"] == "producer":
                producer_reports.append(
                    analyze_producer(run, run["output_path"], execution_sha)
                )
            else:
                null_reports.append(
                    analyze_null(run, run["output_path"], execution_sha)
                )
        except TransportInputError as exc:
            blockers.append(f"{run['run_id']}: {exc}")

    for report in producer_reports + null_reports:
        for violation in report["violations"]:
            blockers.append(f"{report['run_id']}: {violation}")

    positions, position_blockers = position_checks(producer_reports)
    orders, order_blockers = order_checks(producer_reports)
    blockers.extend(position_blockers)
    blockers.extend(order_blockers)

    nulls_toxin_free = len(null_reports) == len(CONTRACT["design"]["seeds"]) and all(
        report.get("toxin_free") for report in null_reports
    )
    if not nulls_toxin_free:
        blockers.append("toxin-free nulls did not hold |bacteriocin_BtuB| <= tolerance")

    endpoint_ordered = bool(orders) and all(
        time_detail["finite_and_complete"]
        and time_detail["endpoint_r50_ordered_0_gt_15_gt_60"]
        and time_detail["endpoint_r90_ordered_0_gt_15_gt_60"]
        and time_detail["endpoint_near_field_ordered_0_lt_15_lt_60"]
        for seed_detail in orders.values() for time_detail in seed_detail["times"]
    )
    interpolation_monotonic = bool(orders) and all(
        time_detail["finite_and_complete"]
        and time_detail["interpolation_r50_nonincreasing"]
        and time_detail["interpolation_r90_nonincreasing"]
        and time_detail["interpolation_near_field_nondecreasing"]
        for seed_detail in orders.values() for time_detail in seed_detail["times"]
    )
    ordered = endpoint_ordered and interpolation_monotonic
    n_paired_times = sum(len(detail["times"]) for detail in orders.values())

    if blockers:
        status = STATUS_BLOCKED
    elif ordered and n_paired_times > 0:
        status = STATUS_PASS
    else:
        status = STATUS_FAIL

    n_rows = write_snapshot_csv(output_dir / "assay_snapshots.csv", producer_reports)
    metrics = {
        "assay_id": CONTRACT["assay_id"],
        "execution_source_sha": execution_sha,
        "results_root": str(args.results_root),
        "predeclared_analysis": DECL,
        "producer_runs": producer_reports,
        "null_runs": null_reports,
        "source_positions": positions,
        "order_checks_by_seed": orders,
        "n_snapshot_rows": n_rows,
    }
    (output_dir / "assay_metrics.json").write_text(
        json.dumps(jsonable(metrics), indent=2, allow_nan=False) + "\n"
    )

    gate = {
        "status": status,
        GATE["id"]: status == STATUS_PASS,
        # This assay measures the transport law; promoting a Stage C gate or
        # launching Stage D remains a separate, recorded maintainer decision.
        "C_transport_gate": False,
        "checks": {
            "all_outputs_returned": all(run["output_path"] is not None for run in runs),
            "all_runs_authenticated": not any(
                report["violations"] for report in producer_reports + null_reports
            ),
            "single_source_per_producer_run": all(
                report.get("n_lysis_events") == 1 for report in producer_reports
            ),
            "source_identical_across_amplitudes": all(
                detail["identical_across_amplitudes"] for detail in positions.values()
            ),
            "source_time_identical_across_amplitudes": all(
                detail["event_time_identical_across_amplitudes"] for detail in positions.values()
            ),
            "paired_times_complete": not order_blockers,
            "n_paired_snapshot_times": n_paired_times,
            "endpoint_original_criteria_0_gt_15_gt_60": endpoint_ordered,
            "adaptive_interpolation_monotonic_0_15_20_30_60": interpolation_monotonic,
            "intrinsic_radii_ordered_0_gt_15_gt_60": endpoint_ordered,
            "nulls_toxin_free": nulls_toxin_free,
        },
        "blockers": blockers,
        "criteria": GATE["criteria"],
    }
    (output_dir / "assay_gate.json").write_text(
        json.dumps(jsonable(gate), indent=2, allow_nan=False) + "\n"
    )

    print(json.dumps({"status": status, "blockers": blockers[:10]}, indent=2))
    if status == STATUS_PASS:
        return 0
    return 2 if status == STATUS_BLOCKED else 1


if __name__ == "__main__":
    raise SystemExit(main())
