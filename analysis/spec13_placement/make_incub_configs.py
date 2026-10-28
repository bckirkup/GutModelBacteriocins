"""Emit the incubator-variant arm configs for the spec13 placement study.

Same C/S/R live-patch arms as make_configs.py, but with every loss
channel damped (contraction, disruption, single-cell loss, advective
turnover, distal transit) and the corrinoid pool lowered to b12=1e-5 so
colicin E faces a penetrable apparent Kd (docs/SPEC14_SECTION8_VALIDATION
landmine: shipped b12=1e-3 leaves the toxin ~inert).  This is the
condition under which a within-patch selection signal can express at
all; the churn variants exercise the physiological baseline.

Outputs land flat in out/ as incub_<ret>_ly<ly>_s<seed>_*.{csv,json}.
"""

import itertools
import json
import os

from make_configs import OUT, SEEDS, arm_config

INCUB_DELTAS = {
    "layer2.contraction_rate_per_min": 0.0,
    "layer2.crypt_disruption_prob": 0.0,
    "layer2.crypt_single_cell_loss_per_h": 0.0,
    "radial_turnover": 1e12,
    "distal_transit": 1e12,
    "b12_initial_conc": 1e-5,
}


def main():
    manifest = []
    for seed, (ret, lysis) in itertools.product(
        SEEDS, [("lit", 0.0), ("lit", 0.10), ("ship", 0.0), ("ship", 0.10)]
    ):
        stem = f"incub_{ret}_ly{lysis:.3g}_s{seed}"
        cfg = arm_config("crypt", ret, lysis, seed, OUT)
        cfg.update(INCUB_DELTAS)
        cfg["layer2.audit_lineage_file"] = os.path.join(OUT, f"{stem}_lineage.csv")
        cfg["layer2.provenance_file"] = os.path.join(OUT, f"{stem}_prov.json")
        cfg["layer2.timeseries_file"] = os.path.join(OUT, f"{stem}_ts.csv")
        path = os.path.join(OUT, f"cfg_{stem}.json")
        with open(path, "w") as fh:
            json.dump(cfg, fh, indent=2)
            fh.write("\n")
        manifest.append(path)
    print(f"wrote {len(manifest)} configs")


if __name__ == "__main__":
    main()
