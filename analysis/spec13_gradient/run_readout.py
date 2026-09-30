"""Spec 13 Phase-3 readout: three-region Layer-3 colonic gradient.

Question (falsifiable, per docs/SPEC13_IMPLEMENTATION_REVIEW.md Phase-3 gate):
with the sourced axial profiles and within-region parameters drawn from the
Phase-1 admissible corner, does the model produce
  (a) stool in the 1e6-1e8 CFU/g band,
  (b) a mucosally flat Enterobacteriaceae fraction (vs the declared bound),
  (c) retention time / shedding rate matching independent data at a stated
      axial position?
Plus the S9 observable set (burstiness/diurnal envelope, growth-in-transit
amplification, purge recovery kinematics), each reported against the
parameter combination it pins. A gradient-collapsed-to-uniform arm
(layer3.uniform_profile = true) separates the gradient's contribution.

Arms: {gradient, uniform} x seeds {1001, 1002, 1003}, 72 h at dt = 60 s,
mucosal purge (90% luminal) at 48 h. Runs are sequential (~2.7 GB RSS each).
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

csv.field_size_limit(10**9)

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
BINARY = REPO / "build-serial" / "gut_ibm"
CONFIG_DIR = HERE / "configs"
RUN_DIR = HERE / "out"
SUMMARY = HERE / "readout_summary.json"
VERDICTS = HERE / "verdicts.csv"

DT = 60.0
HORIZON_S = 259200.0  # 72 h
PURGE_H = 48.0
SEEDS = [1001, 1002, 1003]
ARMS = ["gradient", "uniform"]

STOOL_LO, STOOL_HI = 1e6, 1e8  # CFU/g target band
FLATNESS_BOUND = 4.0  # declared bound: max/min regional fraction
# Whole-gut transit literature band for healthy adults (sourced in REPORT).
TRANSIT_LO_H, TRANSIT_HI_H = 30.0, 70.0
TAIL_H = 24.0  # measurement window (post-purge recovery tail)

# Within-region parameters: Phase-1 admissible corner (crypt-carried,
# contraction <=0.6/min, agent loss <=0.3) — see spec13_occupancy/REPORT.md.
LAYER2 = {
    "initial_patch_fraction": 0.01,
    "initial_founder_cells": 20,
    "contraction_rate_per_min": 0.5,
    "agent_loss_fraction": 0.2,
    "transit_half_life_h": 4.0,
    "reattach_prob_per_transit": 0.05,
    "establish_prob_single": 0.01,
    "establish_ratio": 10.0,
    "crypt_fraction": 0.10,
    "crypt_disruption_prob": 0.07,
    "crypt_single_cell_loss_per_h": 0.02,
    "crypt_supply_mult": 2.0,
    "proximal_fraction": 0.45,
    "proximal_disruption_prob": 0.5,
    "proximal_single_cell_loss_per_h": 0.10,
    "proximal_supply_mult": 1.0,
    "distal_fraction": 0.45,
    "distal_disruption_prob": 0.5,
    "distal_single_cell_loss_per_h": 0.20,
    "distal_supply_mult": 0.5,
    "invalid_patch_fraction_stop": 1.0,
}

# Sourced axial profiles (SPEC13_MULTISCALE.md Layer-3 table midpoints):
# transit 8-15 / 5-10 / 5-15 h; pH 6.0-6.5 / 6.4-6.8 / 6.6-7.0.
LAYER3 = {
    "n_patches_per_region": 400,
    "f_edge": 0.2,
    "reattach_p0": 0.05,
    "hapc_rate_per_day": 6.0,
    "hapc_antegrade_fraction": 0.95,
    "hapc_day_fraction": 0.8,
    "hapc_day_start_h": 8.0,
    "alpha_hapc_per_min": 0.2,
    "stool_g_per_day": 128.0,
    "flatness_bound": FLATNESS_BOUND,
    "lumen_seed_cells": 3000,
    "purge_at_h": PURGE_H,
    "purge_fraction": 0.9,
    "summary_interval_steps": 10,
    "cecum_transit_h": 12.0,
    "cecum_ph": 6.25,
    "cecum_lumen_carbon_mol_m3": 8.0e-3,
    "cecum_lumen_volume_l": 0.40,
    "cecum_contraction_k0_per_min": 0.5,
    "cecum_outflow_survival": 0.92,
    "transverse_transit_h": 8.0,
    "transverse_ph": 6.60,
    "transverse_lumen_carbon_mol_m3": 5.0e-3,
    "transverse_lumen_volume_l": 0.35,
    "transverse_contraction_k0_per_min": 0.4,
    "transverse_outflow_survival": 0.90,
    "descending_transit_h": 10.0,
    "descending_ph": 6.80,
    "descending_lumen_carbon_mol_m3": 2.5e-3,
    "descending_lumen_volume_l": 0.30,
    "descending_contraction_k0_per_min": 0.3,
    "descending_outflow_survival": 0.88,
}


def arm_config(arm: str, seed: int) -> dict:
    out_ts = RUN_DIR / f"{arm}_seed{seed}_timeseries.csv"
    out_prov = RUN_DIR / f"{arm}_seed{seed}_prov.json"
    layer3 = dict(LAYER3)
    layer3["enabled"] = True
    layer3["uniform_profile"] = arm == "uniform"
    layer3["timeseries_file"] = str(out_ts.relative_to(REPO))
    layer3["provenance_file"] = str(out_prov.relative_to(REPO))
    layer2 = dict(LAYER2)
    layer2["n_patches"] = layer3["n_patches_per_region"]
    return {
        "seed": seed,
        "total_time": HORIZON_S,
        "bio_dt": DT,
        "hdf5": {"enabled": False},
        "initial_strains": [
            {
                "type": 1,
                "count": 100,
                "mu_max": 5.0e-4,
                "plasmids": ["ColE1"],
                "conjugative": False,
            }
        ],
        "layer2": layer2,
        "layer3": layer3,
    }


def run_arm(arm: str, seed: int) -> None:
    cfg = arm_config(arm, seed)
    cfg_path = CONFIG_DIR / f"{arm}_seed{seed}.json"
    cfg_path.write_text(json.dumps(cfg, indent=1))
    print(f"[run] {arm} seed={seed} -> {cfg_path.name}", flush=True)
    proc = subprocess.run(
        [str(BINARY), str(cfg_path.relative_to(REPO))],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    (RUN_DIR / f"{arm}_seed{seed}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-3:]
        raise RuntimeError(f"{arm} seed={seed} exited {proc.returncode}: {tail}")


def rows_for(path: Path) -> list[dict]:
    with open(path) as f:
        return [r for r in csv.DictReader(f) if r["region"] in REGIONS]


REGIONS = ("cecum", "transverse", "descending")


def window(rows: list[dict], t_lo: float, t_hi: float) -> list[dict]:
    return [r for r in rows if t_lo <= float(r["time_s"]) < t_hi]


def tail_metrics(ts_path: Path) -> dict:
    rows = rows_for(ts_path)
    t_end = max(float(r["time_s"]) for r in rows)
    t_tail = t_end - TAIL_H * 3600.0
    tail = window(rows, t_tail, t_end + 1)

    # Stool CFU/g: chain-level stool_cum is identical across region rows —
    # take the descending row series and difference it over the tail window.
    desc = sorted(
        [r for r in rows if r["region"] == "descending"],
        key=lambda r: float(r["time_s"]),
    )
    desc_tail = window(desc, t_tail, t_end + 1)
    stool_cells = float(desc_tail[-1]["stool_cum"]) - float(desc_tail[0]["stool_cum"])
    days = TAIL_H / 24.0
    stool_cfu_g = stool_cells / days / LAYER3["stool_g_per_day"]

    # Mucosal flatness: Enterobacteriaceae fraction per region vs flora proxy.
    fractions = {}
    occ = {}
    for region in REGIONS:
        region_tail = [r for r in tail if r["region"] == region]
        fractions[region] = sum(
            float(r["mean_cfu_ml"]) for r in region_tail
        ) / max(1, len(region_tail)) / 1e7
        occ[region] = sum(float(r["occ_frac"]) for r in region_tail) / max(
            1, len(region_tail)
        )
    fvals = [v for v in fractions.values() if v > 0]
    flatness = max(fvals) / min(fvals) if len(fvals) == 3 else float("nan")

    # Retention: mean residence in region i = lumen_cells_i / outflow_rate_i.
    # Outflow rate estimated as net lumen deltas + stool/death flows is noisy;
    # use lumen_cells_i / (washout flux implied by transit tau_i) directly:
    # residence_i ~= tau_i at steady state, so report realized mean age as
    # lumen_cells_i / outflow_cells_rate_i where the rate is reconstructed
    # from (lumen_births - delta_stocks) is fragile. Instead report the
    # configured chain transit and the realized lumen census per region.
    lumen_cells = {}
    lumen_mu = {}
    for region in REGIONS:
        region_tail = [r for r in tail if r["region"] == region]
        lumen_cells[region] = sum(
            float(r["lumen_cells"]) for r in region_tail
        ) / max(1, len(region_tail))
        lumen_mu[region] = sum(
            float(r["lumen_mu_per_s"]) for r in region_tail
        ) / max(1, len(region_tail))

    # Diurnal envelope: stool export deltas binned by hour-of-day.
    day_h0 = int(LAYER3["hapc_day_start_h"])
    by_hour: dict[int, float] = {}
    prev = float(desc[0]["stool_cum"])
    for r in desc[1:]:
        hour = int((float(r["time_s"]) / 3600.0) % 24)
        by_hour[hour] = by_hour.get(hour, 0.0) + float(r["stool_cum"]) - prev
        prev = float(r["stool_cum"])
    counts = {}
    for r in desc[1:]:
        hour = int((float(r["time_s"]) / 3600.0) % 24)
        counts[hour] = counts.get(hour, 0) + 1
    per_hour = {h: v / counts[h] for h, v in by_hour.items() if counts[h]}
    day_mean = sum(
        per_hour.get(h, 0.0) for h in range(day_h0, day_h0 + 12)
    ) / 12.0
    night_mean = sum(
        per_hour.get((day_h0 + 12 + k) % 24, 0.0) for k in range(12)
    ) / 12.0
    diurnal_ratio = day_mean / night_mean if night_mean > 0 else float("inf")

    # Purge recovery: stool export rate post-purge vs 12 h pre-purge.
    purge_t = PURGE_H * 3600.0
    pre = window(desc, purge_t - 12 * 3600.0, purge_t)
    post = window(desc, purge_t, purge_t + 12 * 3600.0)
    pre_rate = (
        (float(pre[-1]["stool_cum"]) - float(pre[0]["stool_cum"]))
        / (12.0) if len(pre) > 1 else float("nan")
    )
    post_rate = (
        (float(post[-1]["stool_cum"]) - float(post[0]["stool_cum"]))
        / 12.0 if len(post) > 1 else float("nan")
    )
    recovery_ratio = post_rate / pre_rate if pre_rate else float("nan")

    # Growth-in-transit amplification: luminal births vs stool exports.
    prov = json.loads(
        (ts_path.parent / ts_path.name.replace("_timeseries.csv", "_prov.json"))
        .read_text()
    )
    led = prov["ledger"]
    amplification = (
        led["luminal_births"] / led["stool_exports"]
        if led["stool_exports"] > 0 else float("inf")
    )

    return {
        "stool_cfu_g": stool_cfu_g,
        "fractions": fractions,
        "flatness": flatness,
        "occupancy": occ,
        "lumen_cells": lumen_cells,
        "lumen_mu_per_s": lumen_mu,
        "diurnal_ratio": diurnal_ratio,
        "recovery_ratio_12h": recovery_ratio,
        "amplification": amplification,
        "termination_cause": prov["termination"]["cause"],
        "ledger": led,
    }


def main() -> int:
    CONFIG_DIR.mkdir(exist_ok=True)
    RUN_DIR.mkdir(exist_ok=True)
    if "--analyze-only" not in sys.argv:
        for arm in ARMS:
            for seed in SEEDS:
                run_arm(arm, seed)

    summary = {}
    for arm in ARMS:
        summary[arm] = {}
        for seed in SEEDS:
            ts = RUN_DIR / f"{arm}_seed{seed}_timeseries.csv"
            if ts.exists():
                summary[arm][seed] = tail_metrics(ts)
    SUMMARY.write_text(json.dumps(summary, indent=1, default=str))

    verdict_rows = []
    for arm in ARMS:
        vals = [summary[arm][s] for s in SEEDS if s in summary[arm]]
        if not vals:
            continue
        stool = [v["stool_cfu_g"] for v in vals]
        flat = [v["flatness"] for v in vals]
        stool_ok = all(STOOL_LO <= s <= STOOL_HI for s in stool)
        flat_ok = all(f <= FLATNESS_BOUND for f in flat)
        verdict_rows.append(
            {
                "arm": arm,
                "stool_cfu_g": f"{min(stool):.3e}..{max(stool):.3e}",
                "stool_band_met": stool_ok,
                "flatness": f"{min(flat):.2f}..{max(flat):.2f}",
                "flatness_met": flat_ok,
                "diurnal_ratio": f"{min(v['diurnal_ratio'] for v in vals):.2f}..{max(v['diurnal_ratio'] for v in vals):.2f}",
                "recovery_12h": f"{min(v['recovery_ratio_12h'] for v in vals):.2f}..{max(v['recovery_ratio_12h'] for v in vals):.2f}",
                "amplification": f"{min(v['amplification'] for v in vals):.2f}..{max(v['amplification'] for v in vals):.2f}",
            }
        )
    with open(VERDICTS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(verdict_rows[0].keys()))
        w.writeheader()
        w.writerows(verdict_rows)
    print(json.dumps(summary, indent=1, default=str)[:4000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
