#!/usr/bin/env python3
"""Analyze the adaptive ecology refinement against its predeclared contract.

Every definition applied here — the shared calendar window, the per-seed
producer-null slope contrast, the selection rule, and the gate precedence — is
read from ``refinement_contract.json``, committed before the runs exist.
Missing or unauthenticated output blocks the gate; absent values are never
inferred as zero and a missing seed is never dropped from a median.

Exit status: 0 on PASS, 1 on READY_FOR_APPROVAL (analysis complete, awaiting
the recorded operator approval), 2 on any BLOCKED status.
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
from collections.abc import Iterator
from contextlib import contextmanager
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
from gut_ibm_tools.transport_metrics import ExpectedRun, authenticate_run

CONTRACT = json.loads((ROOT / "refinement_contract.json").read_text())
DECL = CONTRACT["predeclared_analysis"]
RULE = CONTRACT["selection_rule"]
GATE = CONTRACT["gate"]
SEEDS = [int(s) for s in CONTRACT["design"]["seeds"]]
CANDIDATES = [float(a) for a in RULE["candidates"]]
WINDOW_S = float(DECL["window_s"])
PRODUCER_TYPE = 1
GRID_SPECIES = tuple(CONTRACT["runtime"]["hdf5_schedule"]["grid_species"])
HDF5_SCHEDULE = {
    key: int(value)
    for key, value in CONTRACT["runtime"]["hdf5_schedule"].items()
    if key != "grid_species"
}
STATUS_PASS = "PASS"
STATUS_READY = "READY_FOR_APPROVAL"


@contextmanager
def open_h5(path: Path) -> Iterator[h5py.File]:
    """Open a plain or gzipped HDF5 output; the gz stream is staged to a temp file."""
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


def val(dataset) -> object:
    """Scalar value of an HDF5 dataset, decoding bytes."""
    arr = np.asarray(dataset[()])
    item = arr.reshape(-1)[0] if arr.size else None
    if isinstance(item, (bytes, np.bytes_)):
        return item.decode(errors="replace")
    if isinstance(item, np.generic):
        return item.item()
    return item


def step_keys(group) -> list[str]:
    return sorted(group.keys(), key=lambda name: int(name.rsplit("_", 1)[1]))


def slope(times: list[float], values: list[float]) -> float | None:
    if len(times) < 2 or max(times) == min(times):
        return None
    return float(
        np.polyfit(np.asarray(times, dtype=float), np.asarray(values, dtype=float), 1)[0]
    )


def window_metrics(series: list[tuple], t_end: float | None, width_s: float) -> dict:
    """Slope/tail on series points with t in [t_end - width_s, t_end]."""
    if not series or t_end is None:
        return {
            "slope_log10_ratio_per_h": None,
            "tail_median_log10_ratio": None,
            "samples": 0,
        }
    window = [
        point
        for point in series
        if t_end - width_s - 1e-12 <= point[0] <= t_end + 1e-12
    ]
    per_second = slope([p[0] for p in window], [p[3] for p in window])
    return {
        "slope_log10_ratio_per_h": None if per_second is None else per_second * 3600.0,
        "tail_median_log10_ratio": (
            float(np.median([p[3] for p in window])) if window else None
        ),
        "samples": len(window),
    }


def get_event(last_summary, names: list[str]):
    events = last_summary.get("events")
    if events is None:
        return None
    for name in names:
        if name in events:
            return val(events[name])
    return None


def producer_divisions_from_events(last_summary) -> int | None:
    """Integer cumulative producer divisions from cumulative_divisions_by_type[1].

    Returns None when the dataset is absent (legacy outputs) or too short; a
    composition-weighted proxy is never substituted.
    """
    events = last_summary.get("events")
    if events is None or "cumulative_divisions_by_type" not in events:
        return None
    raw = events["cumulative_divisions_by_type"]
    data = raw[()] if hasattr(raw, "__getitem__") else raw
    try:
        values = list(data)
    except TypeError:
        values = [data]
    if len(values) <= PRODUCER_TYPE:
        return None
    return int(values[PRODUCER_TYPE])


def find_output(results_root: Path, index: int) -> Path | None:
    for name in ("output.h5.gz", "output.h5"):
        candidate = results_root / str(index) / name
        if candidate.exists():
            return candidate
    return None


def composition_series(h5, bio_dt: float) -> list[tuple]:
    """(t, n_type1, n_type2, log10 ratio) per agent snapshot."""
    series = []
    summary = h5["summary"]
    for key in step_keys(h5["agents"]):
        agents = h5["agents"][key]
        types = np.asarray(agents["type"][()])
        step = int(key.rsplit("_", 1)[1])
        t = step * float(bio_dt)
        if key in summary and "time" in summary[key]:
            t = float(val(summary[key]["time"]))
        n1 = int((types == 1).sum())
        n2 = int((types == 2).sum())
        series.append((t, n1, n2, math.log10((n1 + 0.5) / (n2 + 0.5))))
    return series


def expected_run(entry: dict, execution_sha: str) -> ExpectedRun:
    return ExpectedRun(
        execution_source_sha=execution_sha,
        seed=int(entry["seed"]),
        amplitude=float(entry["amplitude"]),
        chemistry_placement=CONTRACT["runtime"]["chemistry_placement"],
        hdf5_schedule=HDF5_SCHEDULE,
        grid_species=GRID_SPECIES,
    )


def analyze_run(
    entry: dict, cfg: dict, manifest: dict, results_root: Path
) -> dict:
    path = find_output(results_root, int(entry["array_index"]))
    row = {
        "run_id": entry["run_id"],
        "arm": entry["arm"],
        "seed": int(entry["seed"]),
        "amplitude": float(entry["amplitude"]),
        "output_status": "missing",
        "output_path": str(path) if path else None,
        "analysis_error": None,
        "execution_source_sha_observed": None,
        "execution_source_sha_match": None,
        "chemistry_placement": None,
        "mpi_rank_count": None,
        "termination_cause": None,
        "t_end_s": None,
        "authentication_violations": [],
        "n_type1": None,
        "n_type2": None,
        "mortality_colicin": None,
        "mortality_lysis": None,
        "producer_divisions": None,
        "kills_per_lysis": None,
        "realized_lysis_per_producer_division": None,
        "slope_log10_ratio_per_h": None,
        "tail_median_log10_ratio": None,
        "window_samples": 0,
    }
    if path is None:
        return row
    try:
        with open_h5(path) as h5:
            for required in ("run_provenance", "summary", "agents"):
                if required not in h5:
                    raise ValueError(f"missing /{required}")
            provenance = h5["run_provenance"]
            observed = (
                str(val(provenance["git_sha"])) if "git_sha" in provenance else None
            )
            placement = (
                str(val(provenance["chemistry_placement"]))
                if "chemistry_placement" in provenance
                else None
            )
            ranks = (
                int(val(provenance["mpi_rank_count"]))
                if "mpi_rank_count" in provenance
                else None
            )
            termination = (
                str(val(provenance["termination_cause"]))
                if "termination_cause" in provenance
                else "missing"
            )
            try:
                row["authentication_violations"] = authenticate_run(
                    h5, expected_run(entry, manifest["execution_source_sha"])
                )
            except Exception as exc:  # noqa: BLE001 - authentication reports, never raises
                row["authentication_violations"] = [f"{type(exc).__name__}: {exc}"]
            series = composition_series(h5, cfg["bio_dt"])
            if not series:
                raise ValueError("no agent snapshots")
            t_end = series[-1][0]
            last = h5["summary"][step_keys(h5["summary"])[-1]]
            kills = get_event(
                last, ["cumulative_mortality_colicin", "mortality_colicin"]
            )
            lysis = get_event(last, ["cumulative_mortality_lysis", "mortality_lysis"])
            divisions = producer_divisions_from_events(last)
            final_types = np.asarray(h5["agents"][step_keys(h5["agents"])[-1]]["type"][()])
            own_end = window_metrics(series, t_end, WINDOW_S)
            row.update(
                {
                    "output_status": (
                        "complete" if termination == "horizon_reached" else "terminated"
                    ),
                    "execution_source_sha_observed": observed,
                    "execution_source_sha_match": observed
                    == manifest["execution_source_sha"],
                    "chemistry_placement": placement,
                    "mpi_rank_count": ranks,
                    "termination_cause": termination,
                    "t_end_s": t_end,
                    "n_type1": int((final_types == 1).sum()),
                    "n_type2": int((final_types == 2).sum()),
                    "mortality_colicin": kills,
                    "mortality_lysis": lysis,
                    "producer_divisions": divisions,
                    "kills_per_lysis": (
                        kills / lysis if kills is not None and lysis else None
                    ),
                    "realized_lysis_per_producer_division": (
                        lysis / divisions
                        if lysis is not None and divisions
                        else None
                    ),
                    "slope_log10_ratio_per_h": own_end[
                        "slope_log10_ratio_per_h"
                    ],
                    "tail_median_log10_ratio": own_end[
                        "tail_median_log10_ratio"
                    ],
                    "window_samples": own_end["samples"],
                    "_series": series,
                }
            )
    except Exception as exc:  # noqa: BLE001 - a bad output is data, not a crash
        row["output_status"] = "invalid"
        row["analysis_error"] = f"{type(exc).__name__}: {exc}"
    return row


def paired_contrast(producer: dict, null: dict, manifest: dict) -> dict:
    """Producer minus the plasmid-free null of its own seed, on t_common."""
    t_common = min(producer["t_end_s"], null["t_end_s"])
    treated = window_metrics(producer.get("_series") or [], t_common, WINDOW_S)
    control = window_metrics(null.get("_series") or [], t_common, WINDOW_S)
    slopes = (
        treated["slope_log10_ratio_per_h"],
        control["slope_log10_ratio_per_h"],
    )
    tails = (
        treated["tail_median_log10_ratio"],
        control["tail_median_log10_ratio"],
    )
    return {
        "run_id": producer["run_id"],
        "control_run_id": null["run_id"],
        "seed": producer["seed"],
        "amplitude": producer["amplitude"],
        "contrast": "producer_minus_same_seed_plasmid_free_null",
        "t_common_s": t_common,
        "window_s": WINDOW_S,
        "slope_log10_ratio_per_h_treatment": treated["slope_log10_ratio_per_h"],
        "slope_log10_ratio_per_h_control": control["slope_log10_ratio_per_h"],
        "delta_slope_log10_ratio_per_h": (
            None if None in slopes else slopes[0] - slopes[1]
        ),
        "tail_median_log10_ratio_treatment": tails[0],
        "tail_median_log10_ratio_control": tails[1],
        "delta_tail_median_log10_ratio": None if None in tails else tails[0] - tails[1],
        "samples_treatment": treated["samples"],
        "samples_control": control["samples"],
        "execution_source_sha_match_both": bool(
            producer.get("execution_source_sha_match")
            and null.get("execution_source_sha_match")
        ),
        "container_image_digest": manifest.get("container_image_digest"),
    }


def median(values: list[float | None]) -> float | None:
    """Median over all seeds; any missing value makes the median missing."""
    if any(value is None for value in values):
        return None
    ordered = sorted(float(v) for v in values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def amplitude_record(
    amplitude: float, producers: list[dict], pairs: list[dict]
) -> dict:
    """Apply the predeclared selection-rule criteria to one candidate amplitude."""
    deltas = {
        int(pair["seed"]): pair["delta_slope_log10_ratio_per_h"]
        for pair in pairs
        if float(pair["amplitude"]) == amplitude
    }
    finals = {
        int(row["seed"]): row["n_type2"]
        for row in producers
        if float(row["amplitude"]) == amplitude
    }
    kills = {
        int(row["seed"]): row["kills_per_lysis"]
        for row in producers
        if float(row["amplitude"]) == amplitude
    }
    delta_list = [deltas.get(seed) for seed in SEEDS]
    final_list = [finals.get(seed) for seed in SEEDS]
    kills_list = [kills.get(seed) for seed in SEEDS]
    record = {
        "amplitude": amplitude,
        "delta_slope_log10_ratio_per_h": deltas,
        "median_delta_slope_log10_ratio_per_h": median(delta_list),
        "final_susceptible_counts": finals,
        "median_final_susceptible": median(final_list),
        "kills_per_lysis": kills,
        "median_kills_per_lysis": median(kills_list),
        "criteria": [],
        "qualifying": False,
        "non_qualifying_reason": None,
    }

    def criterion(description: str, passed: bool, reason: str | None = None) -> bool:
        record["criteria"].append(
            {"criterion": description, "pass": bool(passed), "reason": reason}
        )
        return passed

    missing = [seed for seed, value in zip(SEEDS, delta_list) if value is None]
    if missing:
        criterion(
            "all three seed-level producer-null delta slopes are positive",
            False,
            f"delta slope missing for seeds {missing}; absent values are never "
            "inferred as zero and seeds are never dropped",
        )
    else:
        bad = [seed for seed, value in deltas.items() if not value > 0.0]
        criterion(
            "all three seed-level producer-null delta slopes are positive",
            not bad,
            f"non-positive delta slopes for seeds {sorted(bad)}" if bad else None,
        )
    med = record["median_delta_slope_log10_ratio_per_h"]
    criterion(
        f"median delta slope in [{RULE['median_slope_min']}, "
        f"{RULE['median_slope_max']}] log10 per hour inclusive",
        med is not None
        and RULE["median_slope_min"] <= med <= RULE["median_slope_max"],
        None if med is not None else "median delta slope is missing",
    )
    low = [
        seed
        for seed, count in finals.items()
        if count is None or count < RULE["per_seed_final_susceptible_min"]
    ]
    med_final = record["median_final_susceptible"]
    criterion(
        f"median final susceptible count >= "
        f"{RULE['median_final_susceptible_min']} and no seed below "
        f"{RULE['per_seed_final_susceptible_min']}",
        not low
        and med_final is not None
        and med_final >= RULE["median_final_susceptible_min"],
        (
            f"final susceptible count missing or below "
            f"{RULE['per_seed_final_susceptible_min']} for seeds {sorted(low)}"
            if low
            else None
        ),
    )
    if any(value is None for value in kills_list):
        criterion(
            "median kills per producer lysis available for all three seeds",
            False,
            "kills_per_lysis missing for at least one seed",
        )
    record["qualifying"] = all(item["pass"] for item in record["criteria"])
    if not record["qualifying"]:
        reasons = [c["reason"] for c in record["criteria"] if not c["pass"]]
        record["non_qualifying_reason"] = "; ".join(r for r in reasons if r)
    return record


def select_amplitude(records: list[dict]) -> dict:
    """Resolve the predeclared selection rule over the candidate amplitudes."""
    qualifying = [r for r in records if r["qualifying"]]
    outcome = {"selected_amplitude": None, "basis": None, "per_amplitude": records}
    if qualifying:
        best = max(r["median_kills_per_lysis"] for r in qualifying)
        tolerance = float(RULE["tie_relative_tolerance"])
        tied = [
            r
            for r in qualifying
            if r["median_kills_per_lysis"] >= best * (1.0 - tolerance)
        ]
        chosen = min(tied, key=lambda r: r["amplitude"])
        outcome.update(
            {
                "selected_amplitude": chosen["amplitude"],
                "basis": "qualified",
                "max_median_kills_per_lysis": best,
                "tie_relative_tolerance": tolerance,
                "amplitudes_within_tolerance": sorted(
                    r["amplitude"] for r in tied
                ),
            }
        )
        return outcome
    first = next((r for r in records if r["amplitude"] == CANDIDATES[0]), None)
    deltas = list((first or {}).get("delta_slope_log10_ratio_per_h", {}).values())
    med_final = (first or {}).get("median_final_susceptible")
    reproduces = (
        first is not None
        and len(deltas) == len(SEEDS)
        and all(value is not None and value > 0.0 for value in deltas)
        and med_final is not None
        and med_final >= 500
    )
    if reproduces:
        outcome.update(
            {
                "selected_amplitude": CANDIDATES[0],
                "basis": "fallback_reproduction_of_prior_effect",
            }
        )
    else:
        outcome["basis"] = "no_qualifying_amplitude"
    return outcome


def intrinsic_gate(path: Path | None, manifest: dict) -> tuple[bool | None, list[str]]:
    """True only when the supplied file proves the intrinsic gate at this image."""
    if path is None:
        return None, ["no --intrinsic-gate-file supplied"]
    try:
        data = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return False, [f"intrinsic gate file unreadable: {exc}"]
    problems = []
    if data.get("single_source_transport_gate") is not True:
        problems.append("intrinsic file does not show single_source_transport_gate true")
    if data.get("execution_source_sha") != manifest.get("execution_source_sha"):
        problems.append("intrinsic gate execution_source_sha differs from manifest")
    if data.get("container_image_digest") != manifest.get("container_image_digest"):
        problems.append("intrinsic gate container_image_digest differs from manifest")
    return not problems, problems


def approval_status(
    path: Path | None, manifest: dict, selected: float | None
) -> tuple[bool | None, list[str]]:
    """True only when the operator approval binds this exact gate decision."""
    if path is None:
        return None, ["no --approval-file supplied"]
    try:
        data = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return False, [f"approval file unreadable: {exc}"]
    problems = []
    if data.get("authorize_c_transport_gate") is not True:
        problems.append("approval lacks authorize_c_transport_gate=true")
    if data.get("execution_source_sha") != manifest.get("execution_source_sha"):
        problems.append("approval execution_source_sha differs from manifest")
    if data.get("container_image_digest") != manifest.get("container_image_digest"):
        problems.append("approval container_image_digest differs from manifest")
    if data.get("jobs") != int(CONTRACT["design"]["jobs"]):
        problems.append(f"approval jobs != {CONTRACT['design']['jobs']}")
    if selected is None or not math.isclose(
        float(data.get("selected_mucin_charge_amplitude", float("nan"))),
        float(selected),
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        problems.append(
            "approval selected_mucin_charge_amplitude does not equal the "
            "amplitude the selection rule resolved"
        )
    return not problems, problems


def public_row(row: dict) -> dict:
    return {key: value for key, value in row.items() if not key.startswith("_")}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


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
        "--results-root",
        type=Path,
        default=ROOT / "generated" / "results",
        help="directory holding <array_index>/output.h5[.gz] outputs",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=ROOT / "analysis",
        help="directory the analysis artifacts are written into",
    )
    parser.add_argument("--intrinsic-gate-file", type=Path)
    parser.add_argument("--approval-file", type=Path)
    args = parser.parse_args()
    output_dir = prepare_output_directory(args.results_dir)

    manifest = json.loads((ROOT / "generated" / "manifest.json").read_text())
    rows = []
    for entry in manifest["runs"]:
        cfg = json.loads((ROOT / entry["input_relpath"]).read_text())
        rows.append(analyze_run(entry, cfg, manifest, args.results_root))

    public = [public_row(row) for row in rows]
    csv_rows = []
    for item in public:
        csv_row = dict(item)
        csv_row["authentication_violation_count"] = len(
            item["authentication_violations"]
        )
        csv_row["authentication_violations"] = "; ".join(
            item["authentication_violations"]
        )
        csv_rows.append(csv_row)
    write_csv(output_dir / "run_metrics.csv", csv_rows)
    (output_dir / "run_metrics.json").write_text(
        json.dumps(jsonable(public), indent=2, allow_nan=False) + "\n"
    )

    missing_rows = [
        public_row(row)
        for row in rows
        if row["output_status"] not in ("complete", "terminated")
    ]
    (output_dir / "missing_outputs.json").write_text(
        json.dumps(jsonable(missing_rows), indent=2, allow_nan=False) + "\n"
    )

    complete = [
        row for row in rows if row["output_status"] in ("complete", "terminated")
    ]
    pairs = []
    for row in complete:
        if row["arm"] != "producer":
            continue
        null = next(
            (
                candidate
                for candidate in complete
                if candidate["arm"] == "plasmid_free_null"
                and candidate["seed"] == row["seed"]
            ),
            None,
        )
        if null is not None:
            pairs.append(paired_contrast(row, null, manifest))
    write_csv(output_dir / "paired_metrics.csv", pairs)
    (output_dir / "paired_metrics.json").write_text(
        json.dumps(jsonable(pairs), indent=2, allow_nan=False) + "\n"
    )

    producers = [row for row in complete if row["arm"] == "producer"]
    selection = select_amplitude(
        [amplitude_record(amp, producers, pairs) for amp in CANDIDATES]
    )
    (output_dir / "amplitude_selection.json").write_text(
        json.dumps(jsonable(selection), indent=2, allow_nan=False) + "\n"
    )

    gate = build_gate(manifest, rows, selection, args)
    (output_dir / "refinement_gate.json").write_text(
        json.dumps(jsonable(gate), indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(jsonable({"status": gate["status"], "gate": gate}), indent=2))
    if gate["status"] == STATUS_PASS:
        return 0
    return 1 if gate["status"] == STATUS_READY else 2


def build_gate(manifest: dict, rows: list[dict], selection: dict, args) -> dict:
    """Apply the contract's status precedence, first blocking condition wins."""
    unreadable = [
        row["run_id"]
        for row in rows
        if row["output_status"] not in ("complete", "terminated")
    ]
    auth_problems = []
    for row in rows:
        if row["output_status"] not in ("complete", "terminated"):
            continue
        for violation in row["authentication_violations"]:
            auth_problems.append(f"{row['run_id']}: {violation}")
        if row["execution_source_sha_match"] is not True:
            auth_problems.append(
                f"{row['run_id']}: run_provenance/git_sha does not match the "
                "manifest execution source"
            )
        if row["chemistry_placement"] != "device_delivery":
            auth_problems.append(
                f"{row['run_id']}: chemistry_placement "
                f"{row['chemistry_placement']!r} != 'device_delivery'"
            )
        if row["mpi_rank_count"] != 1:
            auth_problems.append(
                f"{row['run_id']}: mpi_rank_count {row['mpi_rank_count']!r} != 1"
            )
    intrinsic_ok, intrinsic_problems = intrinsic_gate(args.intrinsic_gate_file, manifest)
    selected = selection.get("selected_amplitude")
    approval_ok, approval_problems = (
        (None, [])
        if intrinsic_ok is not True
        else approval_status(args.approval_file, manifest, selected)
    )
    qualified = selection.get("basis") == "qualified" and selected is not None

    if unreadable:
        status = "BLOCKED_MISSING_OUTPUTS"
    elif auth_problems:
        status = "BLOCKED_AUTHENTICATION"
    elif intrinsic_ok is not True:
        status = "BLOCKED_MISSING_INTRINSIC_GATE"
    elif not qualified:
        status = "BLOCKED_NO_QUALIFYING_AMPLITUDE"
    elif approval_ok is True:
        status = STATUS_PASS
    else:
        status = STATUS_READY
    passed = status == STATUS_PASS
    return {
        "gate": GATE["id"],
        "status": status,
        "single_source_transport_gate": intrinsic_ok,
        "C_transport_gate": passed,
        "stage_d_release": "RELEASED" if passed else "BLOCKED",
        "execution_source_sha": manifest.get("execution_source_sha"),
        "container_image_digest": manifest.get("container_image_digest"),
        "selected_amplitude": selected,
        "selection_basis": selection.get("basis"),
        "runs_total": len(rows),
        "runs_readable": len(rows) - len(unreadable),
        "unreadable_runs": unreadable,
        "authentication_problems": auth_problems,
        "intrinsic_gate_problems": intrinsic_problems,
        "approval_problems": approval_problems,
        "fallback_note": (
            "a fallback-retained amplitude never satisfies this gate"
            if selection.get("basis") == "fallback_reproduction_of_prior_effect"
            else None
        ),
    }


if __name__ == "__main__":
    raise SystemExit(main())
