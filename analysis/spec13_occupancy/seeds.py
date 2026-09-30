"""Stage B: seed-robustness of the admissible cells.

Re-runs every stage-A admissible arm at 5 additional seeds and reports,
per cell, how many of the 6 total seeds land in the 0.1-1% occupancy /
1e4-1e5 CFU/mL band.  A cell counts toward the admissible region only if
the band is seed-typical (>=4/6 seeds admissible); cells admissible at
fewer seeds are reported as marginal.
"""

import csv
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.csv"
SEED_RESULTS = HERE / "results_seeds.csv"
sys.path.insert(0, str(HERE))
from sweep import run_arm

SEEDS = [777001, 777002, 777003, 777004, 777005]
POINT_KEYS = [
    "contraction_rate_per_min",
    "exposed_disruption_prob",
    "agent_loss_fraction",
    "crypt_disruption_prob",
    "transit_half_life_h",
    "reattach_prob_per_transit",
    "establish_ratio",
]


def main() -> int:
    rows = list(csv.DictReader(RESULTS.open()))
    adm = [r for r in rows if r.get("admissible") == "True"]
    print(f"{len(adm)} admissible cells x {len(SEEDS)} extra seeds")

    jobs = []
    for r in adm:
        point = {k: float(r[k]) for k in POINT_KEYS}
        jobs.append((r["name"], point))
    # Deduplicate cells (names are unique already).
    out = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        futs = {}
        for name, point in jobs:
            for s in SEEDS:
                f = pool.submit(run_arm, f"{name}_s{s % 1000}", point, s)
                futs[f] = (name, s)
        for done, f in enumerate(as_completed(futs), start=1):
            row = f.result()
            name, s = futs[f]
            out.append(
                (
                    name,
                    s,
                    row.get("occ_total"),
                    row.get("segment_mean_cfu_ml"),
                    row.get("admissible"),
                )
            )
            if done % 500 == 0:
                print(f"{done}/{len(futs)}", flush=True)

    tally = {}
    for name, s, occ, mean, admissible in out:
        tally.setdefault(name, []).append(bool(admissible))
    robust = sum(1 for flags in tally.values() if sum(flags) >= 4)
    marginal = sum(1 for flags in tally.values() if 0 < sum(flags) < 4)
    dead = sum(1 for flags in tally.values() if sum(flags) == 0)

    with SEED_RESULTS.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "seed", "occ_total", "segment_mean_cfu_ml", "admissible"])
        for name, s, occ, mean, admissible in sorted(out):
            w.writerow([name, s, occ, mean, admissible])

    print(f"robust (>=4/6 incl. stage-A): {robust}")
    print(f"marginal (1-3/6): {marginal}")
    print(f"not admissible at any seed: {dead}")
    print(f"wrote {SEED_RESULTS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
