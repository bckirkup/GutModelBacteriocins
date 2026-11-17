"""Summarize the spec13 placement runs into readout_summary.json.

Churn arms read out/<name>_s<seed>/{lineage,prov}.*; incubator arms read
out/incub_<ret>_ly<ly>_s<seed>_*.{csv,json}.  Per arm we record the
per-strain last-alive hour, final alive counts (deduplicating repeated
terminal census emissions per lineage), per-strain emigrant totals
(summing each lineage's cumulative emigrants within strain), the aggregate
ledger, and whether the ledger closed.
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


def deduplicate_lineage_rows(rows, path):
    dedup = {}
    for row in rows:
        key = (round(float(row["time_s"]) * 1000), int(row["lineage_id"]))
        if key in dedup and dedup[key] != row:
            raise ValueError(f"Conflicting duplicate row in {path}: {key}")
        dedup[key] = row
    return list(dedup.values())


def summarize_lineage(path):
    last_alive_h = defaultdict(float)
    fin = defaultdict(int)
    max_emig_by_lineage = defaultdict(int)
    with open(path, encoding="utf-8") as fh:
        rows = deduplicate_lineage_rows(list(csv.DictReader(fh)), path)
    t_max_ms = max((round(float(r["time_s"]) * 1000) for r in rows), default=0)
    for r in rows:
        t = float(r["time_s"])
        t_ms = round(t * 1000)
        lid = int(r["lineage_id"])
        s = strain(lid)
        if int(r["alive_count"]) > 0:
            last_alive_h[s] = max(last_alive_h[s], t / 3600.0)
        max_emig_by_lineage[lid] = max(
            max_emig_by_lineage[lid], int(r["emigrants_cum"])
        )
        if t_ms == t_max_ms:
            fin[s] += int(r["alive_count"])
    emig = {"C": 0, "S": 0, "R": 0}
    for lid, count in max_emig_by_lineage.items():
        emig[strain(lid)] += count
    return last_alive_h, fin, emig


def summarize(label, lineage_csv, prov_json):
    with open(prov_json, encoding="utf-8") as fh:
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
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=2)
        fh.write("\n")
    print(f"wrote {len(rows)} rows to {dest}")


if __name__ == "__main__":
    main()
