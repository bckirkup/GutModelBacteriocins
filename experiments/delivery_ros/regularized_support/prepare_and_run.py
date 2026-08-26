#!/usr/bin/env python3
"""Population-scale acceptance gate for the regularized delivery sink (PR #335).

Arms are identical to the #334 far-field campaign except for grid_dx and
metabolism.delivery_far_field_radius, which now governs both the far-field
concentration read and the support the prescribed mass is deposited over.

Acceptance criterion is population-scale resolution invariance:
final N, funded/demanded fraction, and delivery retry withholding must stop
tracking grid_dx with the radius on, and must still track it with it off.
"""
import json
import os
import subprocess
import time

BASE = "/home/ubuntu/gutibm-campaign/far-field-invariance/C_ff10_res2/input.json"
ROOT = "/home/ubuntu/gutibm-campaign/regularized-support"
BINARY = "/home/ubuntu/campaign-src-335/build-serial/gut_ibm"

ARMS = {
    "P_ff10_res2": (2.0e-6, 1.0e-5),
    "P_ff10_res4": (4.0e-6, 1.0e-5),
    "P_ff10_res6": (6.0e-6, 1.0e-5),
    "P_ff0_res2": (2.0e-6, 0.0),
    "P_ff0_res4": (4.0e-6, 0.0),
    "P_ff0_res6": (6.0e-6, 0.0),
}

COMMENT = [
    "Population-scale acceptance gate for PR #335 (regularized delivery sink).",
    "Nothing here is fitted. Only grid_dx and metabolism.delivery_far_field_radius vary.",
    "ROS hazard is off (oxygen.k_ROS = 0) so the measurement is of growth and supply,",
    "not of the uncited ambient-ROS mortality retired in #332.",
    "Read per arm: final N, cumulative funded vs demanded growth oxygen and carbon,",
    "delivery_reduction_cumulative and delivery_retry_events_cumulative for carbon and oxygen,",
    "reaction clips (must stay 0), funded <= realized, and /run_provenance termination cause.",
    "Radius-on arms must agree across 2/4/6 um; radius-off arms must disagree.",
]


def build_configs():
    with open(BASE) as fh:
        base = json.load(fh)
    for arm, (dx, radius) in ARMS.items():
        cfg = dict(base)
        cfg["grid_dx"] = dx
        cfg["metabolism.delivery_far_field_radius"] = radius
        cfg["_comment"] = COMMENT + [
            f"arm = {arm}. grid_dx = {dx}. delivery_far_field_radius = {radius}."
        ]
        d = os.path.join(ROOT, arm)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "input.json"), "w") as fh:
            json.dump(cfg, fh, indent=1)


def run():
    for arm in ARMS:
        d = os.path.join(ROOT, arm)
        start = time.monotonic()
        with open(os.path.join(d, "full.log"), "w") as log:
            rc = subprocess.call([BINARY, "input.json"], cwd=d,
                                 stdout=log, stderr=subprocess.STDOUT)
        wall = time.monotonic() - start
        with open(os.path.join(d, "wrapper_wall_seconds.txt"), "w") as fh:
            fh.write(f"{wall}\n")
        print(f"{arm}: rc={rc} wall={wall:.1f}s", flush=True)


if __name__ == "__main__":
    assert os.path.exists(BINARY), BINARY
    os.makedirs(ROOT, exist_ok=True)
    build_configs()
    run()
