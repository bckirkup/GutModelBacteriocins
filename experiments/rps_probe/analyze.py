#!/usr/bin/env python3
"""btuB-null resistant-strain probe: analysis at population scale.

Asserted contrasts, not per-agent invariance (#334 passed per agent and was
inert in a real run):

  A1  With a ColE1 producer present, the btuB-null strain (type 3) ends the
      horizon at a higher share of the non-producer population than the
      sensitive strain (type 2).
  A2  In the null arm (identical founders, no ColE1 anywhere) that contrast
      disappears -- otherwise whatever separates type 2 from type 3 is not
      colicin.
  B   The resistance cost, measured as the type3/type2 ratio with no producer
      present, is monotone non-increasing as corrinoid falls.

Every number is labeled with the run's termination cause and delivery
rationing factor; a crowded or early-terminated arm is said to be one.
"""
import itertools
import json
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("/home/ubuntu/gutibm-campaign/rps-probe")


def keys_sorted(group):
    return sorted(group.keys(), key=lambda v: int(v.rsplit("_", 1)[1]))


def scalar(dataset):
    return float(np.asarray(dataset[()]).reshape(-1)[0])


def text(dataset):
    raw = dataset[()]
    if isinstance(raw, bytes):
        return raw.decode()
    return str(np.asarray(raw).reshape(-1)[0])


def species(flux):
    return [bytes(row[row != 0]).decode() for row in flux["species_names"][:]]


def analyze(arm):
    directory = ROOT / arm
    cfg = json.loads((directory / "input.json").read_text())
    with h5py.File(directory / "output.h5", "r") as handle:
        summary_steps = keys_sorted(handle["summary"])
        last = handle["summary"][summary_steps[-1]]
        flux = last["nutrient_flux"]
        carbon = species(flux).index("carbon")
        rationing = min(
            float(handle["summary"][key]["nutrient_flux"][
                "delivery_rationing_factor_interval"][carbon])
            for key in summary_steps[1:])

        agent_steps = keys_sorted(handle["agents"])
        series = []
        for key in agent_steps:
            types = np.asarray(handle["agents"][key]["type"][()])
            series.append({int(t): int((types == t).sum())
                           for t in (1, 2, 3)})
        final = series[-1]
        events = last["events"]
        return {
            "arm": arm,
            "seed": cfg["seed"],
            "b12": cfg.get("b12_initial_conc", 1.0e-3),
            "producer": bool(cfg["initial_strains"][0]["plasmids"]),
            "termination_cause": text(
                handle["run_provenance"]["termination_cause"]),
            "git_sha": text(handle["run_provenance"]["git_sha"])[:12],
            "rationing_factor_c": rationing,
            "n1": final[1], "n2": final[2], "n3": final[3],
            "n_total": final[1] + final[2] + final[3],
            "divisions": int(scalar(events["cumulative_divisions"])),
            "mortality_colicin": int(scalar(
                events["cumulative_mortality_colicin"])),
            "mortality_lysis": int(scalar(events["cumulative_mortality_lysis"])),
            "outflow_boundary": int(scalar(
                events["cumulative_outflow_boundary"])),
            "series": series,
        }


def ratio(row):
    return row["n3"] / max(row["n2"], 1)


def comparable(row):
    """An arm can only report an R/S ratio if both strains are still there.

    In the 3c166c1 run one type was extinct in 8 of 12 arms, and the R/S
    "ratio" was then n3/1 -- an arbitrary number set by whichever founders
    survived the early crash, not by treatment.  Averaging those across seeds
    produced a monotone B sequence that passed the assertion while measuring
    nothing.  A ratio computed against an extinct denominator is not evidence,
    so it is excluded here and the shortfall is reported rather than averaged
    away.
    """
    return row["n2"] > 0 and row["n3"] > 0


def verdict(passed, comparable_n, required):
    if comparable_n < required:
        return f"INSUFFICIENT ({comparable_n}/{required} arms with both strains)"
    return "PASS" if passed else "FAIL"


def main():
    arms = sorted(p.name for p in ROOT.iterdir()
                  if (p / "output.h5").exists())
    rows = [analyze(a) for a in arms]
    (ROOT / "metrics.json").write_text(json.dumps(rows, indent=1) + "\n")

    print(f"{'arm':32s} {'term':18s} {'ration':>7s} "
          f"{'C':>6s} {'S':>6s} {'R':>6s} {'R/S':>7s}  flag")
    for row in rows:
        flag = "" if comparable(row) else "one strain extinct: R/S not usable"
        print(f"{row['arm']:32s} {row['termination_cause']:18s} "
              f"{row['rationing_factor_c']:7.3f} "
              f"{row['n1']:6d} {row['n2']:6d} {row['n3']:6d} "
              f"{ratio(row):7.2f}  {flag}")

    ok = True
    with_producer = [r for r in rows
                     if r["arm"].startswith("A_three") and comparable(r)]
    null = [r for r in rows if r["arm"].startswith("A_null") and comparable(r)]
    rp = float(np.mean([ratio(r) for r in with_producer])) if with_producer else float("nan")
    rn = float(np.mean([ratio(r) for r in null])) if null else float("nan")
    a1 = bool(with_producer) and rp > 1.0
    a2 = bool(with_producer and null) and rn < rp
    ok &= a1 and a2
    print(f"\nA1 R/S > 1 with producer: "
          f"{verdict(a1, len(with_producer), 1)} (R/S={rp:.2f})")
    print(f"A2 contrast vanishes without producer: "
          f"{verdict(a2, min(len(with_producer), len(null)), 1)} "
          f"(null R/S={rn:.2f})")

    cost = {}
    counted = {}
    for row in rows:
        if row["arm"].startswith("B_"):
            counted[row["b12"]] = counted.get(row["b12"], 0) + 1
            if comparable(row):
                cost.setdefault(row["b12"], []).append(ratio(row))
    if counted:
        levels = sorted(counted, reverse=True)
        print("\nB resistance cost, no producer (R/S by corrinoid):")
        for level in levels:
            usable = cost.get(level, [])
            mean = float(np.mean(usable)) if usable else float("nan")
            print(f"   b12={level:.0e}  R/S={mean:.3f}  "
                  f"usable={len(usable)}/{counted[level]}")
        usable_levels = [b for b in levels if len(cost.get(b, [])) >= 2]
        means = [float(np.mean(cost[b])) for b in usable_levels]
        monotone = (len(usable_levels) == len(levels)
                    and all(b <= a + 1e-9
                            for a, b in itertools.pairwise(means)))
        ok &= monotone
        print(f"B cost non-increasing as corrinoid falls: "
              f"{verdict(monotone, len(usable_levels), len(levels))}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
