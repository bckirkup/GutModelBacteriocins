"""Spec 13 Phase-1 occupancy sweep.

Sweeps the SPEC13_MULTISCALE Layer-2 physiological table at segment scale
and asks the falsifiable question: does any combination produce stationary
occupancy in the 0.1-1% band with segment mean 1e4-1e5 CFU/mL?

Design:
  - full factorial over the swept physiological table (3 levels per axis)
    at 3 values of the discriminating fragment:single establishment ratio,
    one seed per cell
  - each arm runs to 24 h; the readout is the tail-6h stationary mean of
    type-stratified occupancy and segment mean CFU/mL
  - arms landing in the admissible band are then re-run at 5 further seeds
    (stage B) before any "admissible region" claim
"""

import csv
import itertools
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
BINARY = REPO / "build-serial" / "gut_ibm"
CONFIG_DIR = HERE / "sweep_configs"
RUN_DIR = HERE / "sweep_runs"
RESULTS = HERE / "results.csv"

DT = 60.0
HORIZON_S = 86400.0  # 24 h
TAIL_S = 21600.0  # stationary window = last 6 h
N_PATCHES = 1000
SUMMARY_STEPS = 6  # CSV row every 6 min

OCC_LO, OCC_HI = 1e-3, 1e-2  # 0.1-1% occupancy band
MEAN_LO, MEAN_HI = 1e4, 1e5  # CFU/mL band

FIXED = {
    "establish_prob_single": 0.01,
    "crypt_fraction": 0.10,
    "crypt_single_cell_loss_per_h": 0.02,
    "crypt_supply_mult": 2.0,
    "proximal_fraction": 0.45,
    "proximal_single_cell_loss_per_h": 0.10,
    "proximal_supply_mult": 1.0,
    "distal_fraction": 0.45,
    "distal_single_cell_loss_per_h": 0.20,
    "distal_supply_mult": 0.5,
}

# The physiological table (each axis 3 levels).
GRID = {
    "contraction_rate_per_min": [0.2, 0.6, 1.0],
    "exposed_disruption_prob": [0.3, 0.5, 0.7],
    "agent_loss_fraction": [0.1, 0.3, 0.5],
    "crypt_disruption_prob": [0.05, 0.10, 0.15],
    "transit_half_life_h": [2.0, 4.0, 8.0],
    "reattach_prob_per_transit": [0.01, 0.05, 0.10],
}
RATIOS = [1.0, 10.0, 30.0, 100.0, 300.0, 1000.0]


def base_config(seed: int, point: dict) -> dict:
    l2 = {
        "enabled": True,
        "n_patches": N_PATCHES,
        "patch_lateral_m": 300e-6,
        "patch_depth_m": 150e-6,
        "mucus_thickness_m": 200e-6,
        "occupancy_min_agents": 1,
        "initial_patch_fraction": 0.05,
        "initial_founder_cells": 10,
        "initial_pool_cells": 0,
        "bloom_factor": 10.0,
        "bloom_sustain_s": 3600.0,
        "spatial_packing_fraction": 0.6,
        "cell_radius_m": 0.5e-6,
        "segment_dysbiosis_threshold": 0.0,
        "segment_guard_window_s": 1800.0,
        "invalid_patch_fraction_stop": 0.5,
        "audit_patch_index": -1,
        "summary_interval_steps": SUMMARY_STEPS,
        "checkpoint_interval_steps": 0,
        "checkpoint_final": False,
    }
    l2.update(FIXED)
    l2.update(point)
    # The spec sweeps one "exposed disruption" axis shared by both
    # exposed patch types; pseudo-keys do not reach the parser.
    l2["proximal_disruption_prob"] = point["exposed_disruption_prob"]
    l2["distal_disruption_prob"] = point["exposed_disruption_prob"]
    l2["transit_half_life_s"] = point["transit_half_life_h"] * 3600.0
    l2.pop("exposed_disruption_prob", None)
    l2.pop("transit_half_life_h", None)
    return {
        "seed": seed,
        "total_time": HORIZON_S,
        "bio_dt": DT,
        "output_interval": 3600,
        "initial_strains": [
            {"type": 1, "count": 100, "mu_max": 5e-4, "plasmids": ["ColE1"]}
        ],
        "layer2": l2,
    }


def run_arm(name: str, point: dict, seed: int) -> dict:
    cfg_path = CONFIG_DIR / f"{name}.json"
    out_dir = RUN_DIR / name
    out_dir.mkdir(parents=True, exist_ok=True)
    ts_path = out_dir / "ts.csv"
    prov_path = out_dir / "prov.json"
    cfg = base_config(seed, point)
    cfg["layer2"]["timeseries_file"] = str(ts_path)
    cfg["layer2"]["provenance_file"] = str(prov_path)
    cfg_path.write_text(json.dumps(cfg, indent=1))

    t0 = time.time()
    proc = subprocess.run(
        [str(BINARY), str(cfg_path)],
        capture_output=True,
        text=True,
        timeout=3600,
        check=False,
    )
    wall = time.time() - t0
    row = {"name": name, "seed": seed, "wall_s": round(wall, 2)}
    row.update(point)
    row["exit_code"] = proc.returncode
    if proc.returncode != 0 or not ts_path.exists():
        row["error"] = (proc.stderr or proc.stdout)[-300:]
        return row

    tail = []
    with ts_path.open() as fh:
        for rec in csv.DictReader(fh):
            if float(rec["time_s"]) >= HORIZON_S - TAIL_S:
                tail.append(rec)
    if not tail:
        row["error"] = "no tail rows"
        return row

    def mean(col):
        return sum(float(r[col]) for r in tail) / len(tail)

    row["occ_crypt"] = mean("occ_crypt")
    row["occ_proximal"] = mean("occ_exposed_proximal")
    row["occ_distal"] = mean("occ_exposed_distal")
    row["occ_total"] = (
        cfg["layer2"]["crypt_fraction"] * row["occ_crypt"]
        + cfg["layer2"]["proximal_fraction"] * row["occ_proximal"]
        + cfg["layer2"]["distal_fraction"] * row["occ_distal"]
    )
    row["segment_mean_cfu_ml"] = mean("segment_mean_cfu_ml")
    row["segment_cfu_cm2"] = mean("segment_cfu_cm2")
    row["pool_cells"] = mean("pool_cells")
    row["patch_cells"] = mean("patch_cells")
    row["density_crypt"] = mean("density_occupied_crypt_cfu_ml")
    row["density_proximal"] = mean("density_occupied_proximal_cfu_ml")
    row["density_distal"] = mean("density_occupied_distal_cfu_ml")
    row["n_occupied"] = mean("n_occupied_total")
    row["admissible"] = (
        OCC_LO <= row["occ_total"] <= OCC_HI
        and MEAN_LO <= row["segment_mean_cfu_ml"] <= MEAN_HI
    )
    return row


def all_points():
    keys = list(GRID)
    for ratio in RATIOS:
        for combo in itertools.product(*(GRID[k] for k in keys)):
            point = dict(zip(keys, combo))
            point["establish_ratio"] = ratio
            yield point


def main() -> int:
    if not BINARY.exists():
        print(f"missing binary: {BINARY}", file=sys.stderr)
        return 2

    arms = [(f"arm_{i:04d}", p) for i, p in enumerate(all_points())]
    print(f"{len(arms)} arms")

    seed_base = 20260929
    rows = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {
            pool.submit(run_arm, name, p, seed_base + k): name
            for k, (name, p) in enumerate(arms)
        }
        for done, fut in enumerate(as_completed(futures), start=1):
            row = fut.result()
            rows.append(row)
            if done % 100 == 0 or row.get("admissible"):
                print(
                    f"[{done}/{len(arms)}] {row['name']} "
                    f"occ={row.get('occ_total', float('nan')):.4f} "
                    f"mean={row.get('segment_mean_cfu_ml', float('nan')):.3e} "
                    f"adm={row.get('admissible')}",
                    flush=True,
                )

    keys = [
        "name",
        "seed",
        "wall_s",
        "exit_code",
        "occ_total",
        "occ_crypt",
        "occ_proximal",
        "occ_distal",
        "n_occupied",
        "segment_mean_cfu_ml",
        "segment_cfu_cm2",
        "density_crypt",
        "density_proximal",
        "density_distal",
        "pool_cells",
        "patch_cells",
        "admissible",
        "error",
    ]
    ordered = sorted(rows, key=lambda r: r["name"])
    extra = sorted({k for r in rows for k in r if k not in keys})
    with RESULTS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys + extra)
        w.writeheader()
        w.writerows(ordered)
    print(f"\nwrote {RESULTS}")
    adm = [r["name"] for r in ordered if r.get("admissible")]
    print(f"admissible arms ({len(adm)}):", adm if adm else "NONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
