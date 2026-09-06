set +e
job_status=0
work=/tmp/gutibm-gpux
rm -rf "${work}"
mkdir -p "${work}"
if [[ "${GPUX_GPU}" == "1" ]]; then
  nvidia-smi || job_status=1
fi
base=examples/diversity_paradox/input.json
[[ -f "${base}" ]] || base=/src/examples/diversity_paradox/input.json
label="${GPUX_ARM}_seed${GPUX_SEED}"
cfg="${work}/${label}.json"
h5="${work}/${label}.outcome.h5"
python3 - "${base}" "${cfg}" "${h5}" <<'PY'
import json, os, sys
base, cfg, h5 = sys.argv[1:4]
c = json.load(open(base))
c.pop("_comment", None)
c["domain_x"] = 0.001
c["domain_y"] = 0.001
c["total_time"] = 86400
c["output_interval"] = 3600
c["seed"] = int(os.environ["GPUX_SEED"])
c["gpu_enabled"] = os.environ["GPUX_GPU"] == "1"
c["profile_steps"] = False
scale = int(os.environ["GPUX_AGENT_SCALE"])
for s in c["initial_strains"]:
    s["count"] = s["count"] * scale
c["immigration"]["rate"] = c["immigration"]["rate"] * scale
if scale > 1:
    c["dysbiosis_threshold"] = 0
c["hdf5_file"] = h5
c["hdf5"] = {"enabled": True, "compression": "gzip", "compression_level": 4,
             "schedule": {"summary": 1, "provenance": 1, "agents": 360,
                          "grid": 0, "lineage": 0, "genome": 0, "grid_species": []}}
json.dump(c, open(cfg, "w"), indent=1)
print("===GPUX_CONFIG_BEGIN " + os.environ["GPUX_ARM"] + "_seed" + os.environ["GPUX_SEED"] + "===")
print(json.dumps({k: c[k] for k in ("domain_x","domain_y","domain_z","grid_dx","total_time","bio_dt","seed","gpu_enabled","initial_strains","immigration","use_fmm") if k in c} | {"dysbiosis_threshold": c.get("dysbiosis_threshold", "default")}))
print("===GPUX_CONFIG_END===")
PY
[[ $? -eq 0 ]] || job_status=1
start=$(date +%s)
( cd /src/build && ./gut_ibm "${cfg}" ) > "${work}/${label}.stdout" 2> "${work}/${label}.stderr"
run_status=$?
wall=$(( $(date +%s) - start ))
echo "===GPUX_RUN ${label} exit=${run_status} wall_seconds=${wall}==="
tail -n 40 "${work}/${label}.stderr"
tail -n 20 "${work}/${label}.stdout"
[[ ${run_status} -eq 0 ]] || job_status=1
if [[ -f "${h5}" ]]; then
  python3 - "${h5}" "${label}" <<'PY'
import json, sys
import h5py
import numpy as np
path, label = sys.argv[1:3]
out = {"label": label}
with h5py.File(path, "r") as f:
    if "run_provenance" in f:
        g = f["run_provenance"]
        prov = {}
        for k in g:
            if not isinstance(g[k], h5py.Dataset):
                continue
            v = g[k][()]
            prov[k] = v.decode() if isinstance(v, bytes) else (v.tolist() if hasattr(v, "tolist") else v)
        for k, v in g.attrs.items():
            prov["attr:" + k] = v.decode() if isinstance(v, bytes) else (v.tolist() if hasattr(v, "tolist") else v)
        out["run_provenance"] = prov
    counts = {}
    if "agents" in f:
        for step in sorted(f["agents"], key=lambda s: int(s.split("_")[-1]) if s.split("_")[-1].isdigit() else 0):
            t = np.asarray(f["agents"][step]["type"])
            u, n = np.unique(t, return_counts=True)
            counts[step] = {int(a): int(b) for a, b in zip(u, n)}
    out["type_counts_by_agent_dump"] = counts
print("===GPUX_META_BEGIN " + label + "===")
print(json.dumps(out, default=str))
print("===GPUX_META_END===")
PY
  python3 -m gut_ibm_tools.gpu_precision export --hdf5 "${h5}" --arm "${GPUX_ARM}" --seed "${GPUX_SEED}" --output "${work}/${label}.precision.jsonl"
  if [[ $? -eq 0 && -f "${work}/${label}.precision.jsonl" ]]; then
    echo "===BENCH_PRECISION_BEGIN ${label}==="
    cat "${work}/${label}.precision.jsonl"
    echo "===BENCH_PRECISION_END==="
  else
    echo "WARNING: precision export failed for ${h5}" >&2
    job_status=1
  fi
else
  echo "WARNING: missing outcome HDF5 ${h5}" >&2
  job_status=1
fi
echo "GPUX_JOB_STATUS=${job_status}"
exit "${job_status}"
