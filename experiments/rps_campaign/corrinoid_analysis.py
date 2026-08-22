"""Corrinoid/BtuB competition sweep analysis.

Reads the 12 runs of the kd_b12_btuB sweep (4 values x 3 seeds) on the
flux-0.3x epithelial delivery arm and reports, at matched simulated time,
target-side receptor occupancy, expected colicin kills (summed per-step
Bernoulli probability, a lower-variance estimator than realized deaths),
realized colicin kills, producer lysis cost and population trajectories.

Expected layout under the run root:
```
corrinoid/results_<arm>/output.h5
corrinoid/results_<arm>/<arm>_trace.{steps,agents}.csv
```
The analysis summary is written to ``corrinoid/corrinoid_analysis.json``.
"""
import argparse
import csv
import json
import os
import statistics
from pathlib import Path

import h5py

KDS = [("kd1e6", 1e-6), ("kd1e5", 1e-5), ("kd1e4", 1e-4), ("kd1e3", 1e-3)]
SEEDS = [1001, 2027, 3313]
B12 = 1.0e-3
KD_TOX = 5.0e-7


def summary_at(h5, step):
    """Cumulative counters and stocks at the last snapshot <= step."""
    keys = sorted(k for k in h5["summary"])
    chosen = None
    for k in keys:
        if int(h5[f"summary/{k}/step"][()][0]) <= step:
            chosen = k
    g = h5[f"summary/{chosen}"]
    n_by_type = g["n_by_type"][()].tolist()
    return {
        "step": int(g["step"][()][0]),
        "time": float(g["time"][()][0]),
        "colicin": int(g["events/cumulative_mortality_colicin"][()][0]),
        "lysis": int(g["events/cumulative_mortality_lysis"][()][0]),
        "sos": int(g["events/cumulative_sos_inductions"][()][0]),
        "divisions": int(g["events/cumulative_divisions"][()][0]),
        "n1": n_by_type[1],
        "n2": n_by_type[2],
    }


def prov(h5):
    p = h5["run_provenance"]
    return {
        "halt_code": int(p["halt_reason_code"][()][0]),
        # halt_time is 0 for full-horizon runs; termination_time is the
        # simulated time actually reached in either case.
        "halt_time": float(p["termination_time"][()][0]),
        "halt_density": float(p["halt_density_cells_per_mL"][()][0]),
    }


def steps_csv(path, step):
    exp = 0.0
    pos_steps = 0
    realized = 0
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if int(row["step"]) > step:
                break
            exp += float(row["sum_type2_kill_prob"])
            if int(row["toxin_positive_type2"]) > 0:
                pos_steps += 1
            realized = int(row["mortality_colicin"])
    return exp, pos_steps, realized


def agents_csv(path, step):
    occ, tox, haz = [], [], []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if int(row["step"]) > step:
                break
            if int(row["type"]) != 2:
                continue
            occ.append(float(row["btuB_occupancy"]))
            tox.append(float(row["btuB_toxin"]))
            haz.append(float(row["btuB_hazard"]))
    return occ, tox, haz


def pct(values, q):
    if not values:
        return 0.0
    ordered = sorted(values)
    idx = min(len(ordered) - 1, int(q * len(ordered)))
    return ordered[idx]


def main(run_root: Path):
    root = run_root / "corrinoid"
    runs = {}
    halt_times = []
    for tag, _kd in KDS:
        for seed in SEEDS:
            stem = f"corr_{tag}_s{seed}"
            with h5py.File(root / f"results_{stem}" / "output.h5", "r") as h5:
                runs[stem] = {"prov": prov(h5)}
            halt_times.append(runs[stem]["prov"]["halt_time"])
    matched_time = min(halt_times)
    matched_step = int(matched_time / 60.0)
    print(f"matched simulated time: {matched_time} s (step {matched_step})")

    rows = []
    for tag, kd in KDS:
        agg = {"kd": kd, "factor": 1.0 + B12 / kd, "seeds": []}
        agg["apparent_kd"] = KD_TOX * agg["factor"]
        for seed in SEEDS:
            stem = f"corr_{tag}_s{seed}"
            d = root / f"results_{stem}"
            with h5py.File(d / "output.h5", "r") as h5:
                s = summary_at(h5, matched_step)
            exp, pos_steps, realized = steps_csv(
                d / f"{stem}_trace.steps.csv", matched_step)
            occ, tox, haz = agents_csv(
                d / f"{stem}_trace.agents.csv", matched_step)
            s.update({
                "seed": seed,
                "expected_kills": exp,
                "positive_toxin_steps": pos_steps,
                "trace_realized": realized,
                "occ_p50": pct(occ, 0.50), "occ_p90": pct(occ, 0.90),
                "occ_p99": pct(occ, 0.99), "occ_max": max(occ or [0.0]),
                "tox_p99": pct(tox, 0.99), "tox_max": max(tox or [0.0]),
                "haz_p99": pct(haz, 0.99),
                "halt_time": runs[stem]["prov"]["halt_time"],
                "halt_code": runs[stem]["prov"]["halt_code"],
            })
            agg["seeds"].append(s)
        rows.append(agg)

    def mean(arm, key):
        return statistics.mean(s[key] for s in arm["seeds"])

    hdr = ("kd_b12_btuB  factor  apparent_Kd   occ_p99   occ_max   "
           "E[kills]  kills  lysis  SOS  div    n1    n2   halt_h")
    print("\n" + hdr)
    for arm in rows:
        print(f"{arm['kd']:11.0e}  {arm['factor']:6.0f}  "
              f"{arm['apparent_kd']:.3e}  "
              f"{mean(arm, 'occ_p99'):.3e} {mean(arm, 'occ_max'):.3e}  "
              f"{mean(arm, 'expected_kills'):8.3f}  "
              f"{mean(arm, 'colicin'):5.1f}  {mean(arm, 'lysis'):5.1f}  "
              f"{mean(arm, 'sos'):3.0f}  {mean(arm, 'divisions'):5.0f}  "
              f"{mean(arm, 'n1'):4.0f}  {mean(arm, 'n2'):4.0f}  "
              f"{mean(arm, 'halt_time') / 3600.0:5.1f}")

    print("\nper-seed detail")
    for arm in rows:
        for s in arm["seeds"]:
            print(f"  kd={arm['kd']:.0e} seed={s['seed']} "
                  f"E[kills]={s['expected_kills']:.3f} "
                  f"kills={s['colicin']} lysis={s['lysis']} sos={s['sos']} "
                  f"div={s['divisions']} n1={s['n1']} n2={s['n2']} "
                  f"occ_p99={s['occ_p99']:.3e} occ_max={s['occ_max']:.3e} "
                  f"tox_max={s['tox_max']:.3e} "
                  f"pos_steps={s['positive_toxin_steps']} "
                  f"halt={s['halt_time'] / 3600.0:.1f}h/{s['halt_code']}")

    out = root / "corrinoid_analysis.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({"matched_time_s": matched_time, "arms": rows}, fh, indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        default=os.environ.get("GUTIBM_CAMPAIGN_ROOT", "."),
        help="campaign root containing corrinoid/ (default: "
        "GUTIBM_CAMPAIGN_ROOT or current directory)",
    )
    main(Path(parser.parse_args().run_root).expanduser().resolve())
