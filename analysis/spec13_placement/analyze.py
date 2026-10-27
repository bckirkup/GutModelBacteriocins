"""Summarize the spec13 placement runs into readout_summary.json.

Churn arms read out/<name>_s<seed>/{lineage,prov}.*; incubator arms read
out/incub_<ret>_ly<ly>_s<seed>_*.{csv,json}.  Per arm we record the
per-strain last-alive hour, final alive counts, per-strain emigrant
totals, the aggregate ledger, and whether the ledger closed.
"""

import csv
import glob
import json
import os
import re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def strain(lineage_id):
    lin = int(lineage_id)
    if lin <= 10:
        return "C"
    if lin <= 20:
        return "S"
    return "R"


def summarize_lineage(path):
    last_alive_h = defaultdict(float)
    fin = defaultdict(int)
    emig = defaultdict(int)
    t_max = 0.0
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        t = float(r["time_s"])
        t_max = max(t_max, t)
        s = strain(r["lineage_id"])
        if int(r["alive_count"]) > 0:
            last_alive_h[s] = max(last_alive_h[s], t / 3600.0)
        emig[s] = max(emig[s], int(r["emigrants_cum"]))
    for r in rows:
        if float(r["time_s"]) == t_max:
            fin[strain(r["lineage_id"])] += int(r["alive_count"])
    return last_alive_h, fin, emig


def summarize(label, lineage_csv, prov_json):
    with open(prov_json) as fh:
        p = json.load(fh)
    last_alive_h, fin, emig = summarize_lineage(lineage_csv)
    return {
        "arm": label,
        "last_alive_h": {k: round(v, 2) for k, v in last_alive_h.items()},
        "final_alive": dict(fin),
        "emigrants": dict(emig),
        "ledger": p["ledger"],
        "ledger_closed": p["ledger_closed"],
    }


def main():
    rows = []
    for pf in sorted(
        glob.glob(os.path.join(OUT, "crypt_lit_ly*_s*", "prov.json"))
        + glob.glob(os.path.join(OUT, "prox_lit_ly*_s*", "prov.json"))
    ):
        m = re.search(r"out/([\w.]+_s\d+)/prov\.json", pf)
        rows.append(
            summarize(m.group(1), pf.replace("prov.json", "lineage.csv"), pf)
        )
    for pf in sorted(glob.glob(os.path.join(OUT, "incub_*_prov.json"))):
        stem = pf[: -len("_prov.json")]
        rows.append(
            summarize(os.path.basename(stem), stem + "_lineage.csv", pf)
        )
    dest = os.path.join(HERE, "readout_summary.json")
    with open(dest, "w") as fh:
        json.dump(rows, fh, indent=2)
        fh.write("\n")
    print(f"wrote {len(rows)} rows to {dest}")


if __name__ == "__main__":
    main()
