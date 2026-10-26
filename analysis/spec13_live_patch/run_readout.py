#!/usr/bin/env python3
"""Compare a segment live-patch lineage CSV against a standalone gut_ibm
HDF5 output (Layer-1 reference) on lineage-resolved alive counts.

Usage: run_readout.py iso_lineage.csv l1_ref.h5
Exits 0 when every shared step matches on every lineage.
"""
import collections
import csv
import sys

import h5py


def main() -> int:
    seg = collections.defaultdict(dict)
    with open(sys.argv[1]) as fh:
        for row in csv.DictReader(fh):
            seg[int(row["step"])][int(row["lineage_id"])] = (
                int(row["alive_count"]),
                int(row["emigrants_cum"]),
            )
    mismatches = []
    with h5py.File(sys.argv[2]) as f:
        steps = sorted(int(k.split("_")[1]) for k in f["agents"])
        for s in steps:
            counts = collections.Counter(
                int(x) for x in f[f"agents/step_{s:06d}/lineage_id"][:]
            )
            for lin, n in sorted(counts.items()):
                sa = seg.get(s, {}).get(lin, (0, 0))[0]
                if sa != n:
                    mismatches.append((s, lin, n, sa))
            for lin, (sa, _em) in sorted(seg.get(s, {}).items()):
                if sa != counts.get(lin, 0):
                    mismatches.append((s, lin, counts.get(lin, 0), sa))
    print(f"ref steps {steps}; seg steps {sorted(seg)}; "
          f"mismatches {len(mismatches)}")
    for m in mismatches[:10]:
        print("  step=%d lineage=%d ref=%d seg=%d" % m)
    return 1 if mismatches else 0


if __name__ == "__main__":
    raise SystemExit(main())
