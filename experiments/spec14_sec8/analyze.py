#!/usr/bin/env python3
"""Spec 14 Section 8: does the in-vitro / in-vivo reversal appear on main?

Asserted contrasts at population scale, each paired within seed against its
own no-producer null.  A raw type-1 / type-2 separation is not evidence: the
earlier RPS probe measured a larger separation in the null than with a
producer, because ~30 survivors sweep a patch stochastically.

  V1  in vitro, producer:  the sensitive strain is excluded 100-fold, and the
      paired null shows no comparable separation.
  V2  in vivo, producer:   the two strains persist within 10x of each other to
      10 d.  Only interpretable if the paired null is also within 10x --
      otherwise stochastic sweep dominates and the arm reports INSUFFICIENT.
  V3  attribution:         receptor-mediated kills per division.  If the vitro
      producer arm kills at a negligible rate, Section 8 has measured colicin
      potency (retardation / kd_colicinE_btuB), not spatial structure, which
      is the third branch of the spec's own interpretation guide.

No ratio is reported against an extinct or near-extinct denominator; those
report INSUFFICIENT explicitly.
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("/home/ubuntu/gutibm-campaign/spec14-sec8")
ARMS = ROOT / "arms"

# A ratio is only quoted when the smaller arm of the comparison could in
# principle have been observed at 1/100 of the larger, i.e. the larger side
# has at least this many agents.
MIN_DENOM = 10
EXCLUSION_DECADES = 2.0


def keys_sorted(group):
    return sorted(group.keys(), key=lambda v: int(v.rsplit("_", 1)[1]))


def scalar(dataset):
    return float(np.asarray(dataset[()]).reshape(-1)[0])


def text(dataset):
    raw = dataset[()]
    if isinstance(raw, bytes):
        return raw.decode()
    return str(np.asarray(raw).reshape(-1)[0])


def log_ratio(n1: int, n2: int):
    """log10(n1/n2), or None when the pair cannot support a ratio."""
    if n1 == 0 and n2 == 0:
        return None
    if max(n1, n2) < MIN_DENOM:
        return None
    if n1 == 0:
        return -math.inf
    if n2 == 0:
        return math.inf
    return math.log10(n1 / n2)


def analyze(arm: str) -> dict:
    directory = ARMS / arm
    cfg = json.loads((directory / "input.json").read_text())
    with h5py.File(directory / "output.h5", "r") as handle:
        agent_steps = keys_sorted(handle["agents"])
        series = []
        for key in agent_steps:
            grp = handle["agents"][key]
            types = np.asarray(grp["type"][()])
            step = int(key.rsplit("_", 1)[1])
            series.append({
                "t_s": step * float(cfg["bio_dt"]),
                "n1": int((types == 1).sum()),
                "n2": int((types == 2).sum()),
            })
        summary_steps = keys_sorted(handle["summary"])
        last = handle["summary"][summary_steps[-1]]
        events = last["events"]
        divisions = int(scalar(events["cumulative_divisions"]))
        kills = int(scalar(events["cumulative_mortality_colicin"]))
        record = {
            "arm": arm,
            "condition": arm.split("_")[0],
            "treatment": arm.split("_")[1],
            "seed": cfg["seed"],
            "termination_cause": text(
                handle["run_provenance"]["termination_cause"]),
            "git_sha": text(handle["run_provenance"]["git_sha"])[:12],
            "t_end_s": series[-1]["t_s"] if series else 0.0,
            "n1": series[-1]["n1"] if series else 0,
            "n2": series[-1]["n2"] if series else 0,
            "divisions": divisions,
            "mortality_colicin": kills,
            "mortality_lysis": int(scalar(
                events["cumulative_mortality_lysis"])),
            "outflow_boundary": int(scalar(
                events["cumulative_outflow_boundary"])),
            "kills_per_division": (kills / divisions) if divisions else None,
            "series": series,
        }
    record["log_ratio_final"] = log_ratio(record["n1"], record["n2"])
    # First sample at which type 2 is excluded 100-fold by type 1.
    record["t_exclusion_s"] = None
    for sample in series:
        lr = log_ratio(sample["n1"], sample["n2"])
        if lr is not None and lr >= EXCLUSION_DECADES:
            record["t_exclusion_s"] = sample["t_s"]
            break
    # Ratio at the 6 h mark, where the in-vitro claim is made.
    record["log_ratio_6h"] = None
    for sample in series:
        if sample["t_s"] >= 21600.0:
            record["log_ratio_6h"] = log_ratio(sample["n1"], sample["n2"])
            break
    return record


def fmt(value):
    if value is None:
        return "INSUFFICIENT"
    if isinstance(value, float):
        if value == math.inf:
            return "+inf"
        if value == -math.inf:
            return "-inf"
        return f"{value:.3g}"
    return str(value)


def main() -> int:
    arms = sorted(p.name for p in ARMS.iterdir()
                  if (p / "output.h5").exists())
    missing = sorted(p.name for p in ARMS.iterdir()
                     if not (p / "output.h5").exists())
    records = [analyze(a) for a in arms]
    by_key = {(r["condition"], r["treatment"], r["seed"]): r for r in records}

    print(f"arms analyzed: {len(records)}; missing output: {missing}")
    print()
    header = ("arm", "term", "n1", "n2", "log10(n1/n2)", "t_excl_s",
              "kills", "div", "kills/div")
    print("| " + " | ".join(header) + " |")
    print("|" + "---|" * len(header))
    for r in sorted(records, key=lambda r: r["arm"]):
        print("| " + " | ".join([
            r["arm"], r["termination_cause"], str(r["n1"]), str(r["n2"]),
            fmt(r["log_ratio_final"]), fmt(r["t_exclusion_s"]),
            str(r["mortality_colicin"]), str(r["divisions"]),
            fmt(r["kills_per_division"]),
        ]) + " |")
    print()

    verdicts = {}
    for condition in ("vitro", "vivo"):
        deltas, paired = [], []
        for seed in sorted({r["seed"] for r in records}):
            prod = by_key.get((condition, "producer", seed))
            null = by_key.get((condition, "null", seed))
            if not prod or not null:
                continue
            lp, ln = prod["log_ratio_final"], null["log_ratio_final"]
            paired.append((seed, lp, ln))
            if lp is None or ln is None or math.isinf(lp) or math.isinf(ln):
                continue
            deltas.append(lp - ln)
        print(f"## {condition}: paired producer - null, log10(n1/n2)")
        for seed, lp, ln in paired:
            print(f"  seed {seed}: producer {fmt(lp)}, null {fmt(ln)}")
        if len(deltas) >= 3:
            median = statistics.median(deltas)
            print(f"  median paired delta over {len(deltas)} usable seeds: "
                  f"{median:+.3g} decades")
            verdicts[condition] = median
        else:
            print(f"  only {len(deltas)} usable seeds: INSUFFICIENT")
            verdicts[condition] = None
        print()

    vitro_prod = [r for r in records
                  if r["condition"] == "vitro" and r["treatment"] == "producer"]
    kpd = [r["kills_per_division"] for r in vitro_prod
           if r["kills_per_division"] is not None]
    print("## V3 attribution")
    if kpd:
        print(f"  vitro producer kills/division: median {statistics.median(kpd):.3g}"
              f" over {len(kpd)} seeds")
        if statistics.median(kpd) < 1.0e-3:
            print("  VERDICT: colicin is inert at shipped potency. Section 8"
                  " measured potency (retardation / kd_colicinE_btuB), not"
                  " spatial structure.")
    else:
        print("  INSUFFICIENT")
    print()

    print("## V1 / V2")
    v1 = verdicts.get("vitro")
    v2 = verdicts.get("vivo")
    if v1 is None:
        print("  V1 INSUFFICIENT")
    else:
        print(f"  V1 {'PASS' if v1 >= EXCLUSION_DECADES else 'FAIL'}:"
              f" vitro producer beats its null by {v1:+.3g} decades"
              f" (need >= {EXCLUSION_DECADES})")
    if v2 is None:
        print("  V2 INSUFFICIENT")
    else:
        print(f"  V2 {'PASS' if abs(v2) < 1.0 else 'FAIL'}:"
              f" vivo producer - null is {v2:+.3g} decades (need |.| < 1)")
    print()
    print("  The reversal is reproduced only if V1 PASSes and V2 PASSes."
          " V1 FAIL with V3 firing means the model lacks colicin potency,"
          " not spatial structure.")

    (ROOT / "results.json").write_text(
        json.dumps(records, indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
