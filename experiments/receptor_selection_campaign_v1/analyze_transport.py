#!/usr/bin/env python3
"""Source-centred Stage C bacteriocin transport profiles.

Builds radial profiles around active producer-lysis burst sources using the
model's exponential release window (tau=300 s, prune after 5 tau = 1500 s),
nearest-source voxel assignment, and active-source-weight normalization.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import shutil
import tempfile
from collections import defaultdict
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
CAUSE_LYSIS = 5
PRODUCER_STRAIN = 1
BURST_TAU_S = 300.0
BURST_PRUNE_S = 5.0 * BURST_TAU_S  # 1500 s
BIN_EDGES_UM = np.array([0.0, 10.0, 25.0, 50.0, 75.0, 100.0], dtype=float)
BIN_LABELS = ["0-10", "10-25", "25-50", "50-75", "75-100"]
D_FREE = 4.0e-11
PI_COLE1 = 9.0
PH = 7.0
R_MIN = 1.2
DZ_HALF = 1.35
WIDTH = 1.0
EXPECTED_TRANSPORT = {
    0: {"retardation": 1.200, "D_eff": 3.333e-11},
    15: {"retardation": 13.456, "D_eff": 2.973e-12},
    60: {"retardation": 50.225, "D_eff": 7.964e-13},
}
NULL_ABS_TOL = 1e-15
RETARDATION_RTOL = 1e-3
DEFF_RTOL = 1e-3


def retardation_from_pI(amplitude: float) -> float:
    exponent = (DZ_HALF - (PI_COLE1 - PH)) / WIDTH
    return R_MIN + amplitude / (1.0 + 10.0**exponent)


def open_h5(path: Path):
    if path.suffix != ".gz":
        return h5py.File(path, "r"), None
    tmp = tempfile.NamedTemporaryFile(suffix=".h5", delete=False)
    tmp.close()
    with gzip.open(path, "rb") as src, open(tmp.name, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return h5py.File(tmp.name, "r"), Path(tmp.name)


def step_index(name: str) -> int:
    return int(name.rsplit("_", 1)[1])


def steps(group) -> list[str]:
    return sorted(group.keys(), key=step_index)


def min_image_delta(delta: float, length: float) -> float:
    return delta - length * round(delta / length)


def voxel_centers(nx: int, ny: int, nz: int, dx: float):
    # HDF5 grid is stored (nz, ny, nx); centers at (i+0.5)*dx.
    xs = (np.arange(nx, dtype=float) + 0.5) * dx
    ys = (np.arange(ny, dtype=float) + 0.5) * dx
    zs = (np.arange(nz, dtype=float) + 0.5) * dx
    zz, yy, xx = np.meshgrid(zs, ys, xs, indexing="ij")
    return xx, yy, zz


def load_lysis_sources(h, bio_dt: float):
    """Return list of (t_s, x, y, z) for producer lysis events."""
    out = []
    if "provenance" not in h:
        return out
    for sk in steps(h["provenance"]):
        g = h["provenance"][sk]
        if "cause" not in g or g["cause"].shape[0] == 0:
            continue
        cause = np.asarray(g["cause"][()])
        strain = np.asarray(g["strain"][()])
        x = np.asarray(g["x"][()], dtype=float)
        y = np.asarray(g["y"][()], dtype=float)
        z = np.asarray(g["z"][()], dtype=float)
        t = float(step_index(sk) * bio_dt)
        mask = (cause == CAUSE_LYSIS) & (strain == PRODUCER_STRAIN)
        for i in np.flatnonzero(mask):
            out.append((t, float(x[i]), float(y[i]), float(z[i])))
    return out


def active_sources(sources, t_grid: float):
    active = []
    for t_i, x, y, z in sources:
        age = t_grid - t_i
        if age < -1e-12 or age > BURST_PRUNE_S + 1e-12:
            continue
        w = math.exp(-age / BURST_TAU_S)
        active.append((w, x, y, z, age))
    return active


def nearest_source_distances(xx, yy, zz, active, Lx: float, Ly: float):
    n = xx.size
    best_d = np.full(n, np.inf, dtype=float)
    flat_x = xx.ravel()
    flat_y = yy.ravel()
    flat_z = zz.ravel()
    for _w, sx, sy, sz, _age in active:
        dx = min_image_delta_vec(flat_x - sx, Lx)
        dy = min_image_delta_vec(flat_y - sy, Ly)
        dz = flat_z - sz  # non-periodic
        d = np.sqrt(dx * dx + dy * dy + dz * dz)
        best_d = np.minimum(best_d, d)
    return best_d.reshape(xx.shape)


def min_image_delta_vec(delta: np.ndarray, length: float) -> np.ndarray:
    return delta - length * np.round(delta / length)


def profile_snapshot(conc, xx, yy, zz, active, Lx: float, Ly: float):
    w_sum = float(sum(a[0] for a in active))
    if w_sum <= 0.0:
        return None
    dist_m = nearest_source_distances(xx, yy, zz, active, Lx, Ly)
    dist_um = dist_m * 1.0e6
    c_norm = conc / w_sum

    bin_means = []
    bin_counts = []
    for lo, hi in zip(BIN_EDGES_UM[:-1], BIN_EDGES_UM[1:]):
        mask = (dist_um >= lo) & (dist_um < hi)
        count = int(mask.sum())
        bin_counts.append(count)
        bin_means.append(float(np.mean(c_norm[mask])) if count else float("nan"))

    # Concentration-mass radii: radius enclosing 50%/90% of normalized
    # concentration mass, with distance measured to the nearest active source.
    flat_c = np.maximum(c_norm.ravel(), 0.0)
    flat_d = dist_um.ravel()
    total = float(flat_c.sum())
    if total <= 0.0:
        r50 = float("nan")
        r90 = float("nan")
    else:
        order = np.argsort(flat_d)
        cum = np.cumsum(flat_c[order]) / total
        d_sorted = flat_d[order]
        r50 = float(d_sorted[np.searchsorted(cum, 0.50, side="left")])
        r90 = float(d_sorted[np.searchsorted(cum, 0.90, side="left")])
    # Diagnostic shape radii from shell-mean vs r (measure dr); not used by gate.
    r50_shape, r90_shape = radial_mean_mass_radii(c_norm, dist_um)

    return {
        "n_active_sources": len(active),
        "active_source_weight_sum": w_sum,
        "bin_means": bin_means,
        "bin_counts": bin_counts,
        "near_field": bin_means[0],
        "intermediate": bin_means[1],
        "far_field": bin_means[3] if len(bin_means) > 3 else float("nan"),
        # far-field requested as 50–100 µm = bins 50-75 and 75-100 combined mean
        "far_field_50_100": _combined_mean(c_norm, dist_um, 50.0, 100.0),
        "r50_um": r50,
        "r90_um": r90,
        "r50_shape_um": r50_shape,
        "r90_shape_um": r90_shape,
    }


def _combined_mean(c_norm, dist_um, lo, hi):
    mask = (dist_um >= lo) & (dist_um < hi)
    if not np.any(mask):
        return float("nan")
    return float(np.mean(c_norm[mask]))


def radial_mean_mass_radii(c_norm: np.ndarray, dist_um: np.ndarray, r_max_um: float = 100.0, n_bins: int = 200):
    """r50/r90 of shell-mean concentration versus radius (measure dr)."""
    edges = np.linspace(0.0, r_max_um, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    means = np.zeros(n_bins, dtype=float)
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (dist_um >= lo) & (dist_um < hi)
        if np.any(mask):
            means[i] = float(np.mean(c_norm[mask]))
    weights = np.maximum(means, 0.0)
    total = float(weights.sum())
    if total <= 0.0:
        return float("nan"), float("nan")
    cum = np.cumsum(weights) / total
    r50 = float(centers[min(int(np.searchsorted(cum, 0.50, side="left")), n_bins - 1)])
    r90 = float(centers[min(int(np.searchsorted(cum, 0.90, side="left")), n_bins - 1)])
    return r50, r90


def find_output(results_root: Path, index: int) -> Path | None:
    for name in ("output.h5.gz", "output.h5"):
        p = results_root / str(index) / name
        if p.exists():
            return p
    return None


def stage_c_runs(results_root: Path):
    manifest = json.loads((ROOT / "generated" / "campaign_manifest.json").read_text())
    runs = []
    for entry in manifest["runs"]:
        if entry.get("stage") != "C":
            continue
        cfg = json.loads((ROOT / entry["input_relpath"]).read_text())
        camp = cfg["_campaign"]
        runs.append(
            {
                "array_index": entry["array_index"],
                "run_id": camp["run_id"],
                "arm": camp["arm"],
                "seed": int(cfg["seed"]),
                "amplitude": float(camp["axes"]["amplitude"]),
                "bio_dt": float(cfg["bio_dt"]),
                "dx": float(cfg["grid_dx"]),
                "Lx": float(cfg["domain_x"]),
                "Ly": float(cfg["domain_y"]),
                "Lz": float(cfg["domain_z"]),
                "input_path": ROOT / entry["input_relpath"],
                "output_path": find_output(results_root, entry["array_index"]),
                "burst_release_tau": float(cfg.get("burst_release_tau", BURST_TAU_S)),
                "burst_size_cfg": cfg.get("burst_size"),  # may be absent (plasmid default)
            }
        )
    return runs


def analyze_null(run) -> dict:
    path = run["output_path"]
    if path is None:
        return {
            "run_id": run["run_id"],
            "seed": run["seed"],
            "status": "missing_output",
            "max_abs_bacteriocin_BtuB": None,
            "toxin_free": False,
        }
    h, tmp = open_h5(path)
    try:
        max_abs = 0.0
        n_grids = 0
        for sk in steps(h["grid"]):
            conc = np.asarray(h["grid"][sk]["bacteriocin_BtuB"][()], dtype=float)
            max_abs = max(max_abs, float(np.max(np.abs(conc))))
            n_grids += 1
        return {
            "run_id": run["run_id"],
            "seed": run["seed"],
            "status": "ok",
            "n_grid_snapshots": n_grids,
            "max_abs_bacteriocin_BtuB": max_abs,
            "toxin_free": max_abs <= NULL_ABS_TOL,
        }
    finally:
        h.close()
        if tmp:
            tmp.unlink(missing_ok=True)


def analyze_producer(run) -> dict:
    path = run["output_path"]
    base = {
        "run_id": run["run_id"],
        "seed": run["seed"],
        "amplitude": run["amplitude"],
        "array_index": run["array_index"],
    }
    if path is None:
        return {**base, "status": "missing_output", "snapshots": [], "n_excluded_no_source": 0}

    h, tmp = open_h5(path)
    try:
        sources = load_lysis_sources(h, run["bio_dt"])
        grid_keys = steps(h["grid"])
        sample = np.asarray(h["grid"][grid_keys[0]]["bacteriocin_BtuB"][()])
        nz, ny, nx = sample.shape
        xx, yy, zz = voxel_centers(nx, ny, nz, run["dx"])

        # Burst identity check from resolved config when present.
        burst_tau = run["burst_release_tau"]
        burst_size = None
        if "run_provenance" in h and "resolved_config" in h["run_provenance"]:
            raw = h["run_provenance"]["resolved_config"][()]
            if isinstance(raw, bytes):
                raw = raw.decode()
            elif hasattr(raw, "tobytes"):
                raw = raw.tobytes().decode()
            resolved = json.loads(str(raw))
            burst_tau = float(resolved.get("burst_release_tau", burst_tau))
            # ColE1 burst_size lives on the plasmid, not always a top-level key.
            # Prefer explicit key if present.
            if "burst_size" in resolved:
                burst_size = float(resolved["burst_size"])

        snaps = []
        n_excluded = 0
        for sk in grid_keys:
            if "time" in h["summary"].get(sk, {}):
                t_grid = float(np.asarray(h["summary"][sk]["time"][()]).reshape(-1)[0])
            else:
                t_grid = float(step_index(sk) * run["bio_dt"])
            active = active_sources(sources, t_grid)
            if not active:
                n_excluded += 1
                continue
            conc = np.asarray(h["grid"][sk]["bacteriocin_BtuB"][()], dtype=float)
            prof = profile_snapshot(conc, xx, yy, zz, active, run["Lx"], run["Ly"])
            if prof is None:
                n_excluded += 1
                continue
            row = {
                "run_id": run["run_id"],
                "seed": run["seed"],
                "amplitude": run["amplitude"],
                "grid_step": step_index(sk),
                "t_grid_s": t_grid,
                "n_active_sources": prof["n_active_sources"],
                "active_source_weight_sum": prof["active_source_weight_sum"],
                "near_field_0_10": prof["near_field"],
                "intermediate_10_25": prof["intermediate"],
                "far_field_50_100": prof["far_field_50_100"],
                "r50_um": prof["r50_um"],
                "r90_um": prof["r90_um"],
                "r50_shape_um": prof["r50_shape_um"],
                "r90_shape_um": prof["r90_shape_um"],
            }
            for label, mean, count in zip(BIN_LABELS, prof["bin_means"], prof["bin_counts"]):
                row[f"bin_{label}_mean"] = mean
                row[f"bin_{label}_n_voxels"] = count
            # Drop zero-mass snapshots (active sources but undefined radii).
            if math.isnan(prof["r50_um"]) or math.isnan(prof["r90_um"]):
                n_excluded += 1
                continue
            snaps.append(row)

        return {
            **base,
            "status": "ok",
            "burst_release_tau": burst_tau,
            "burst_size": burst_size,
            "n_grid_snapshots": len(grid_keys),
            "n_excluded_no_source": n_excluded,
            "n_qualified": len(snaps),
            "snapshots": snaps,
        }
    finally:
        h.close()
        if tmp:
            tmp.unlink(missing_ok=True)


def median_or_none(vals):
    vals = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not vals:
        return None
    return float(np.median(np.asarray(vals, dtype=float)))


def jsonable(x):
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return None
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    if isinstance(x, list):
        return [jsonable(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    return x


def seed_aggregate(snaps: list[dict]) -> dict:
    keys = [
        "near_field_0_10",
        "intermediate_10_25",
        "far_field_50_100",
        "r50_um",
        "r90_um",
        "r50_shape_um",
        "r90_shape_um",
        "n_active_sources",
        "active_source_weight_sum",
    ]
    out = {"n_qualified_snapshots": len(snaps)}
    for k in keys:
        out[k] = median_or_none([s[k] for s in snaps])
    for label in BIN_LABELS:
        out[f"bin_{label}_mean"] = median_or_none([s[f"bin_{label}_mean"] for s in snaps])
    return out


def ordered_nondecreasing(vals) -> bool:
    return all(a <= b + 1e-15 for a, b in zip(vals, vals[1:]))


def ordered_nonincreasing(vals) -> bool:
    return all(a >= b - 1e-15 for a, b in zip(vals, vals[1:]))


def write_profiles_png(path: Path, seed_rows: list[dict]):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    amps = sorted({int(r["amplitude"]) for r in seed_rows})
    colors = {0: "#4C78A8", 15: "#F58518", 60: "#54A24B"}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    # Left: radial bin medians (seed-level medians averaged for display as three lines)
    ax = axes[0]
    centers = 0.5 * (BIN_EDGES_UM[:-1] + BIN_EDGES_UM[1:])
    for amp in amps:
        rows = [r for r in seed_rows if int(r["amplitude"]) == amp]
        ys = []
        for label in BIN_LABELS:
            ys.append(median_or_none([r[f"bin_{label}_mean"] for r in rows]))
        ax.plot(centers, ys, marker="o", label=f"amp {amp}", color=colors.get(amp))
    ax.set_xlabel("Source distance (µm)")
    ax.set_ylabel("Source-weight-normalized mean conc.")
    ax.set_title("Radial bacteriocin_BtuB profiles")
    ax.legend(frameon=False)

    # Right: near-field and r50 by amplitude (three seeds)
    ax = axes[1]
    for amp in amps:
        rows = sorted([r for r in seed_rows if int(r["amplitude"]) == amp], key=lambda x: x["seed"])
        xs = [amp + 0.15 * (i - 1) for i in range(len(rows))]
        ax.scatter(xs, [r["near_field_0_10"] for r in rows], color=colors.get(amp), label=f"near amp {amp}")
    ax.set_xticks(amps)
    ax.set_xlabel("Mucin-charge amplitude")
    ax.set_ylabel("Near-field 0–10 µm (normalized)")
    ax.set_title("Near-field by seed")
    ax.legend(frameon=False, fontsize=8)

    fig.savefig(path, dpi=140)
    plt.close(fig)
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-root", type=Path, default=ROOT / "generated" / "stage_C" / "results")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "analysis")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    runs = stage_c_runs(args.results_root)
    producers = [r for r in runs if r["arm"] == "producer"]
    nulls = [r for r in runs if r["arm"] == "shared_null"]

    null_reports = [analyze_null(r) for r in nulls]
    producer_reports = [analyze_producer(r) for r in producers]

    profile_rows = []
    for rep in producer_reports:
        profile_rows.extend(rep.get("snapshots") or [])

    # CSV of every qualified snapshot
    csv_path = args.output_dir / "transport_profiles.csv"
    if profile_rows:
        keys = list(profile_rows[0].keys())
        with csv_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(profile_rows)
    else:
        csv_path.write_text("")

    # Aggregate: within amplitude × seed, median across qualified snapshots
    by_amp_seed: dict[tuple[float, int], list[dict]] = defaultdict(list)
    for row in profile_rows:
        by_amp_seed[(float(row["amplitude"]), int(row["seed"]))].append(row)

    seed_level = []
    for (amp, seed), snaps in sorted(by_amp_seed.items()):
        agg = seed_aggregate(snaps)
        seed_level.append({"amplitude": amp, "seed": seed, **agg})

    # Across-seed summary per amplitude
    amp_summary = {}
    for amp in (0.0, 15.0, 60.0):
        rows = [r for r in seed_level if r["amplitude"] == amp]
        seeds = [r["seed"] for r in rows]

        def pack(key):
            vals = [r[key] for r in rows]
            finite = [v for v in vals if v is not None and not math.isnan(v)]
            return {
                "by_seed": {str(s): v for s, v in zip(seeds, vals)},
                "median": median_or_none(finite),
                "min": min(finite) if finite else None,
                "max": max(finite) if finite else None,
            }

        amp_summary[str(int(amp))] = {
            "n_seeds": len(rows),
            "qualified_snapshots_by_seed": {str(r["seed"]): r["n_qualified_snapshots"] for r in rows},
            "near_field_0_10": pack("near_field_0_10"),
            "intermediate_10_25": pack("intermediate_10_25"),
            "far_field_50_100": pack("far_field_50_100"),
            "r50_um": pack("r50_um"),
            "r90_um": pack("r90_um"),
        }

    # Burst identity across producers
    taus = {r.get("burst_release_tau") for r in producer_reports if r.get("status") == "ok"}
    burst_tau_ok = taus == {BURST_TAU_S}

    # pI-law resolution check
    transport_law = {}
    law_ok = True
    for amp in (0, 15, 60):
        r = retardation_from_pI(float(amp))
        deff = D_FREE / r
        exp = EXPECTED_TRANSPORT[amp]
        r_ok = abs(r - exp["retardation"]) <= RETARDATION_RTOL * abs(exp["retardation"])
        d_ok = abs(deff - exp["D_eff"]) <= DEFF_RTOL * abs(exp["D_eff"])
        transport_law[str(amp)] = {
            "retardation_resolved": r,
            "retardation_expected": exp["retardation"],
            "D_eff_resolved": deff,
            "D_eff_expected": exp["D_eff"],
            "match": bool(r_ok and d_ok),
        }
        law_ok = law_ok and r_ok and d_ok

    # Gate checks
    null_ok = all(n.get("toxin_free") for n in null_reports) and len(null_reports) == 3
    snapshot_ok = True
    for amp in (0.0, 15.0, 60.0):
        for seed in (20260911, 20260913, 20260917):
            nq = next(
                (r["n_qualified_snapshots"] for r in seed_level if r["amplitude"] == amp and r["seed"] == seed),
                0,
            )
            if nq < 2:
                snapshot_ok = False

    near_order_ok = True
    r50_order_ok = True
    r90_order_ok = True
    order_detail = {}
    for seed in (20260911, 20260913, 20260917):
        near = []
        r50 = []
        r90 = []
        for amp in (0.0, 15.0, 60.0):
            row = next((r for r in seed_level if r["amplitude"] == amp and r["seed"] == seed), None)
            near.append(None if row is None else row["near_field_0_10"])
            r50.append(None if row is None else row["r50_um"])
            r90.append(None if row is None else row["r90_um"])
        near_ok = all(v is not None and not math.isnan(v) for v in near) and ordered_nondecreasing(near)
        r50_ok = all(v is not None and not math.isnan(v) for v in r50) and ordered_nonincreasing(r50)
        r90_ok = all(v is not None and not math.isnan(v) for v in r90) and ordered_nonincreasing(r90)
        order_detail[str(seed)] = {
            "near_field_0_10": near,
            "near_ordered_0_le_15_le_60": near_ok,
            "r50_um": r50,
            "r50_ordered_0_ge_15_ge_60": r50_ok,
            "r90_um": r90,
            "r90_ordered_0_ge_15_ge_60": r90_ok,
        }
        near_order_ok = near_order_ok and near_ok
        r50_order_ok = r50_order_ok and r50_ok
        r90_order_ok = r90_order_ok and r90_ok

    gate_pass = bool(
        null_ok and snapshot_ok and near_order_ok and r50_order_ok and r90_order_ok and law_ok and burst_tau_ok
    )

    metrics = {
        "campaign_id": "receptor_selection_campaign_v1",
        "stage": "C",
        "results_root": str(args.results_root),
        "bin_edges_um": BIN_EDGES_UM.tolist(),
        "burst_release_tau_s": BURST_TAU_S,
        "burst_prune_s": BURST_PRUNE_S,
        "burst_release_tau_identical": burst_tau_ok,
        "observed_burst_release_tau": sorted(t for t in taus if t is not None),
        "null_checks": null_reports,
        "producer_run_status": [
            {
                "run_id": r["run_id"],
                "seed": r["seed"],
                "amplitude": r["amplitude"],
                "status": r["status"],
                "n_qualified": r.get("n_qualified"),
                "n_excluded_no_source": r.get("n_excluded_no_source"),
                "burst_release_tau": r.get("burst_release_tau"),
            }
            for r in producer_reports
        ],
        "seed_level_medians": seed_level,
        "amplitude_summary": amp_summary,
        "transport_law": transport_law,
        "order_checks_by_seed": order_detail,
        "n_profile_rows": len(profile_rows),
    }
    (args.output_dir / "transport_metrics.json").write_text(
        json.dumps(jsonable(metrics), indent=2, allow_nan=False) + "\n"
    )

    gate = {
        "status": "PASS" if gate_pass else "FAIL",
        "C_transport_gate": gate_pass,
        "checks": {
            "nulls_toxin_free": null_ok,
            "ge_2_qualified_snapshots_per_amp_seed": snapshot_ok,
            "near_field_ordered_0_le_15_le_60_all_seeds": near_order_ok,
            "r50_ordered_0_ge_15_ge_60_all_seeds": r50_order_ok,
            "r90_ordered_0_ge_15_ge_60_all_seeds": r90_order_ok,
            "transport_law_match": law_ok,
            "burst_release_tau_identical": burst_tau_ok,
        },
        "selected_mucin_charge_amplitude_if_pass": 15,
        "null_max_abs": {n["run_id"]: n.get("max_abs_bacteriocin_BtuB") for n in null_reports},
        "order_detail": order_detail,
        "notes": [
            "r50/r90 are volume-weighted concentration-mass radii vs nearest active source.",
            "Early sparse-source snapshots have large Voronoi cells and inflate radii; this can break amplitude ordering of seed medians even when late, density-matched snapshots are ordered.",
            "r50_shape_um/r90_shape_um are diagnostic shell-mean (measure dr) radii and are not gate inputs.",
        ],
    }
    (args.output_dir / "transport_gate.json").write_text(
        json.dumps(jsonable(gate), indent=2, allow_nan=False) + "\n"
    )

    png_ok = write_profiles_png(args.output_dir / "transport_profiles.png", seed_level)
    metrics["transport_profiles_png_written"] = png_ok
    (args.output_dir / "transport_metrics.json").write_text(
        json.dumps(jsonable(metrics), indent=2, allow_nan=False) + "\n"
    )

    print(json.dumps({"transport_gate": gate["status"], "C_transport_gate": gate_pass, "n_profile_rows": len(profile_rows)}, indent=2))
    return 0 if gate_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
