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

Issue #467 (GAP-05) hygiene:
- Preserves legacy 145-sample tail keys (`47h50m-71h50m` inclusive) so
  historical `verdicts.csv` comparisons remain invariant.
- Separates strictly post-purge tail-window statistics (`t > 48 h`, 143 samples
  `173,400-258,600 s`, where `purged_cum > 0`) from the last CSV sample
  (`71h50m`, `258,600 s`) and exact `72 h` (`259,200 s`) `prov.json` endpoints.
- Reports per-patch and per-occupied-patch cell densities alongside count-
  thresholded occupancy (`N >= 1`, `N >= 10`, `N >= 50` cells/patch, reading
  direct threshold columns when exported or homogeneous-Poisson reference
  expectations from `mean_cells_per_patch` on legacy summaries).
"""

import csv
import importlib.util
import io
import json
import math
import subprocess
import sys
import types
from pathlib import Path

csv.field_size_limit(10**9)

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

_REPO_PYTHON = REPO / "python"
sys.path.insert(0, str(_REPO_PYTHON))
try:
    from gut_ibm_tools.path_utils import (
        validate_input_path,
        write_json_file,
        write_text_file,
    )
except ModuleNotFoundError as error:
    if error.name != "h5py":
        raise
    package = types.ModuleType("gut_ibm_tools")
    package.__path__ = [str(_REPO_PYTHON / "gut_ibm_tools")]
    sys.modules["gut_ibm_tools"] = package
    module_path = _REPO_PYTHON / "gut_ibm_tools" / "path_utils.py"
    spec = importlib.util.spec_from_file_location(
        "gut_ibm_tools.path_utils", module_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {module_path}") from error
    path_utils = importlib.util.module_from_spec(spec)
    sys.modules["gut_ibm_tools.path_utils"] = path_utils
    spec.loader.exec_module(path_utils)
    validate_input_path = path_utils.validate_input_path
    write_json_file = path_utils.write_json_file
    write_text_file = path_utils.write_text_file

BINARY = REPO / "build-serial" / "gut_ibm"
CONFIG_DIR = HERE / "configs"
RUN_DIR = HERE / "out"
SUMMARY = HERE / "readout_summary.json"
VERDICTS = HERE / "verdicts.csv"

DT = 60.0
HORIZON_S = 259200.0  # 72 h
PURGE_H = 48.0
FLORA_DENSITY_CFU_ML = 1.0e8
SEEDS = [1001, 1002, 1003]
ARMS = ["gradient", "uniform"]
COUNT_THRESHOLDS = (1, 10, 50)

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

REGIONS = ("cecum", "transverse", "descending")


def arm_config(
    arm: str,
    seed: int,
    *,
    reattach_p0: float | None = None,
    run_dir: Path = RUN_DIR,
) -> dict:
    out_ts = run_dir / f"{arm}_seed{seed}_timeseries.csv"
    out_prov = run_dir / f"{arm}_seed{seed}_prov.json"
    layer3 = dict(LAYER3)
    layer3["enabled"] = True
    layer3["uniform_profile"] = arm == "uniform"
    if reattach_p0 is not None:
        layer3["reattach_p0"] = float(reattach_p0)
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
    write_json_file(cfg_path, cfg, indent=1)
    print(f"[run] {arm} seed={seed} -> {cfg_path.name}", flush=True)
    proc = subprocess.run(
        [str(BINARY), str(cfg_path.relative_to(REPO))],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    write_text_file(
        RUN_DIR / f"{arm}_seed{seed}.log", proc.stdout + proc.stderr
    )
    if proc.returncode != 0:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-3:]
        raise RuntimeError(f"{arm} seed={seed} exited {proc.returncode}: {tail}")


def rows_for(path: Path) -> list[dict]:
    with open(validate_input_path(path)) as f:
        return [r for r in csv.DictReader(f) if r["region"] in REGIONS]


def window(rows: list[dict], t_lo: float, t_hi: float) -> list[dict]:
    return [r for r in rows if t_lo <= float(r["time_s"]) < t_hi]


def poisson_tail_ge(lam: float, k_min: int) -> float:
    """Return P(X >= k_min) for X ~ Poisson(lam)."""
    if k_min <= 0:
        return 1.0
    if lam <= 0.0:
        return 0.0
    if lam >= 100.0 and k_min <= 50:
        # Complement sum in log space avoids exp(-lam) underflow.
        cdf_lt = 0.0
        for k in range(k_min):
            log_pmf = -lam + k * math.log(lam) - math.lgamma(k + 1)
            cdf_lt += math.exp(log_pmf)
        return max(0.0, min(1.0, 1.0 - cdf_lt))
    term = math.exp(-lam)
    cdf_lt = term
    for k in range(1, k_min):
        term *= lam / k
        cdf_lt += term
    return max(0.0, min(1.0, 1.0 - cdf_lt))


def _region_threshold_metrics(
    region_rows: list[dict],
    n_patches: int,
) -> dict:
    """Summarize patch cell density and count-thresholded occupancy over rows."""
    n_rows = max(1, len(region_rows))
    occ_ge_1 = sum(float(r["occ_frac"]) for r in region_rows) / n_rows
    mean_cells_patch = (
        sum(float(r["patch_cells"]) / max(1, n_patches) for r in region_rows)
        / n_rows
    )
    mean_cells_occupied = (
        sum(
            float(r["patch_cells"])
            / max(1.0, float(r["occ_frac"]) * n_patches)
            for r in region_rows
        )
        / n_rows
    )

    sample_row = region_rows[0] if region_rows else {}
    has_direct_cols = all(
        f"occ_ge_{thr}_frac" in sample_row for thr in (10, 50)
    )
    if has_direct_cols:
        occ_ge_10 = (
            sum(float(r["occ_ge_10_frac"]) for r in region_rows) / n_rows
        )
        occ_ge_50 = (
            sum(float(r["occ_ge_50_frac"]) for r in region_rows) / n_rows
        )
        threshold_source = "timeseries_direct"
    else:
        occ_ge_10 = (
            sum(
                poisson_tail_ge(float(r["patch_cells"]) / max(1, n_patches), 10)
                for r in region_rows
            )
            / n_rows
        )
        occ_ge_50 = (
            sum(
                poisson_tail_ge(float(r["patch_cells"]) / max(1, n_patches), 50)
                for r in region_rows
            )
            / n_rows
        )
        threshold_source = "poisson_reference_from_mean_cells_per_patch"

    return {
        "occ_ge_1_frac": occ_ge_1,
        "occ_ge_10_frac": occ_ge_10,
        "occ_ge_50_frac": occ_ge_50,
        "threshold_source": threshold_source,
        "mean_cells_per_patch": mean_cells_patch,
        "mean_cells_per_occupied_patch": mean_cells_occupied,
    }


def _strict_post_purge_metrics(
    rows: list[dict],
    desc: list[dict],
    prov: dict,
) -> dict:
    """Compute strictly post-purge (`time_s > purge_at_s`) window & 72h endpoint metrics."""
    purge_at_s = float(
        prov.get("layer3", {}).get("purge_at_s", PURGE_H * 3600.0)
    )
    n_patches = int(
        prov.get("layer3", {}).get(
            "n_patches_per_region", LAYER3["n_patches_per_region"]
        )
    )
    stool_g_per_day = float(
        prov.get("layer3", {}).get("stool_g_per_day", LAYER3["stool_g_per_day"])
    )

    post_rows = [r for r in rows if float(r["time_s"]) > purge_at_s]
    desc_post = [r for r in desc if float(r["time_s"]) > purge_at_s]

    if len(desc_post) >= 2:
        span_days = (
            float(desc_post[-1]["time_s"]) - float(desc_post[0]["time_s"])
        ) / 86400.0
        stool_delta = float(desc_post[-1]["stool_cum"]) - float(
            desc_post[0]["stool_cum"]
        )
        post_stool_cfu_g = (
            stool_delta / span_days / stool_g_per_day
            if span_days > 0.0
            else float("nan")
        )
        time_range_s = [
            float(desc_post[0]["time_s"]),
            float(desc_post[-1]["time_s"]),
        ]
    else:
        post_stool_cfu_g = float("nan")
        time_range_s = []

    post_fractions: dict[str, float] = {}
    post_occ: dict[str, float] = {}
    post_density: dict[str, dict] = {}
    for region in REGIONS:
        r_post = [r for r in post_rows if r["region"] == region]
        n_r = max(1, len(r_post))
        post_fractions[region] = (
            sum(float(r["mean_cfu_ml"]) for r in r_post)
            / n_r
            / FLORA_DENSITY_CFU_ML
        )
        post_occ[region] = sum(float(r["occ_frac"]) for r in r_post) / n_r
        post_density[region] = _region_threshold_metrics(r_post, n_patches)

    fvals = [v for v in post_fractions.values() if v > 0.0]
    post_flatness = (
        max(fvals) / min(fvals) if len(fvals) == len(REGIONS) else float("nan")
    )

    last_t = max(float(r["time_s"]) for r in rows)
    last_csv: dict[str, dict] = {}
    for region in REGIONS:
        r_reg = sorted(
            [r for r in rows if r["region"] == region],
            key=lambda x: float(x["time_s"]),
        )
        last_row = r_reg[-1]
        occ_val = float(last_row["occ_frac"])
        p_cells = float(last_row["patch_cells"])
        last_csv[region] = {
            "time_s": last_t,
            "occ_frac": occ_val,
            "patch_cells": int(round(p_cells)),
            "mean_cells_per_patch": p_cells / max(1, n_patches),
            "mean_cells_per_occupied_patch": p_cells
            / max(1.0, occ_val * n_patches),
            "mean_cfu_ml": float(last_row["mean_cfu_ml"]),
            "fraction_at_1e8": float(last_row["mean_cfu_ml"])
            / FLORA_DENSITY_CFU_ML,
        }

    endpoint_72h_regions: dict[str, dict] = {}
    for reg_obs in prov.get("observables", {}).get("regions", []):
        rname = reg_obs["name"]
        occ_val = float(reg_obs["occupancy_frac"])
        p_cells = float(reg_obs["patch_cells"])
        endpoint_72h_regions[rname] = {
            "occupancy_frac": occ_val,
            "patch_cells": int(round(p_cells)),
            "mean_cells_per_patch": p_cells / max(1, n_patches),
            "mean_cells_per_occupied_patch": p_cells
            / max(1.0, occ_val * n_patches),
            "mean_cfu_ml": float(reg_obs["mean_cfu_ml"]),
            "fraction_at_1e8": float(reg_obs["mean_cfu_ml"])
            / FLORA_DENSITY_CFU_ML,
        }

    return {
        "post_purge_window_gt_48h": {
            "purge_at_s": purge_at_s,
            "sample_count_per_region": len(desc_post),
            "time_range_s": time_range_s,
            "stool_cfu_g": post_stool_cfu_g,
            "fractions": post_fractions,
            "flatness": post_flatness,
            "occupancy": post_occ,
            "patch_density_and_thresholds": post_density,
        },
        "last_csv_sample": last_csv,
        "endpoint_72h_prov": {
            "wall_time_s": float(prov.get("wall_time_s", HORIZON_S)),
            "cumulative_stool_cfu_per_g": float(
                prov.get("observables", {}).get("stool_cfu_per_g", float("nan"))
            ),
            "regions": endpoint_72h_regions,
        },
    }


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
    stool_cells = float(desc_tail[-1]["stool_cum"]) - float(
        desc_tail[0]["stool_cum"]
    )
    days = TAIL_H / 24.0
    stool_cfu_g = stool_cells / days / LAYER3["stool_g_per_day"]

    # Mucosal flatness: Enterobacteriaceae fraction per region vs flora proxy.
    # Denominator = shipped layer3 flora_density_cfu_ml default (1e8 cells/mL,
    # FISH-scale total community; docs/gutibm_density_reconciliation.md).
    fractions = {}
    occ = {}
    for region in REGIONS:
        region_tail = [r for r in tail if r["region"] == region]
        fractions[region] = (
            sum(float(r["mean_cfu_ml"]) for r in region_tail)
            / max(1, len(region_tail))
            / FLORA_DENSITY_CFU_ML
        )
        occ[region] = sum(float(r["occ_frac"]) for r in region_tail) / max(
            1, len(region_tail)
        )
    fvals = [v for v in fractions.values() if v > 0]
    flatness = max(fvals) / min(fvals) if len(fvals) == 3 else float("nan")

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
    counts: dict[int, int] = {}
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
        (float(pre[-1]["stool_cum"]) - float(pre[0]["stool_cum"])) / 12.0
        if len(pre) > 1
        else float("nan")
    )
    post_rate = (
        (float(post[-1]["stool_cum"]) - float(post[0]["stool_cum"])) / 12.0
        if len(post) > 1
        else float("nan")
    )
    recovery_ratio = post_rate / pre_rate if pre_rate else float("nan")

    # Growth-in-transit amplification: luminal births vs stool exports.
    prov_path = ts_path.parent / ts_path.name.replace(
        "_timeseries.csv", "_prov.json"
    )
    prov = json.loads(validate_input_path(prov_path).read_text())
    led = prov["ledger"]
    amplification = (
        led["luminal_births"] / led["stool_exports"]
        if led["stool_exports"] > 0
        else float("inf")
    )

    extra = _strict_post_purge_metrics(rows, desc, prov)

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
        **extra,
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
    write_json_file(SUMMARY, summary, indent=1)

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
    header = list(verdict_rows[0].keys())
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=header)
    w.writeheader()
    w.writerows(verdict_rows)
    write_text_file(VERDICTS, buf.getvalue())
    print(json.dumps(summary, indent=1, default=str)[:4000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
