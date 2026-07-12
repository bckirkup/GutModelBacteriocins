# Diversity Campaign

Parameter sweep testing the Advective Double-Bind / Combinatorial Washout Trap
with all mechanisms engaged.

## Conditions

| Config | O2 | Acetate | Crypts | Adaptive dt | `kd_corrinoid_btuB` | Purpose |
|--------|----|---------|--------|-------------|----------------------|---------|
| `baseline.json` | OFF | OFF | OFF | OFF | default (`1e-9`) | Reference |
| `full_mechanisms.json` | ON | ON | ON | ON | default (`1e-9`) | Full mechanisms |
| `kd_sweep_1e-9.json` | ON | ON | ON | ON | `1e-9` | Strong block; colicin E weak |
| `kd_sweep_1e-8.json` | ON | ON | ON | ON | `1e-8` | Moderate block |
| `kd_sweep_1e-7.json` | ON | ON | ON | ON | `1e-7` | Weak block |
| `kd_sweep_1e-6.json` | ON | ON | ON | ON | `1e-6` | Minimal block; colicin E potent |

All configs use a 2 mm × 2 mm × 100 µm domain, 2 µm grid spacing, seven days,
500 residents plus 100 immigrants, FMM-accelerated QSSA, and motility with
substep integration.

## Resource warning

These are production-scale experiments, not smoke tests. The grid contains
`1000 × 1000 × 50 = 50,000,000` cells. Concentration and reaction arrays alone
require approximately 7.2 GB for the nine-species baseline and 8.0 GB when
oxygen adds a tenth species, before agents, HDF5 buffers, MPI, or GPU mirrors.
On WSL2, review `docs/WSL2_SETUP.md`, keep enough RAM available to Windows, and
start with one MPI rank before launching the full campaign.

HDF5 schedule values are measured in simulation steps. The baseline's
60-step schedule is hourly because its timestep is fixed at 60 seconds. Output
spacing in the adaptive-timestep experiments varies with the selected timestep.

## Batch runs

```bash
# Build
cmake -B build -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON
cmake --build build -j"${GUTIBM_BUILD_JOBS:-$(nproc)}"

# Inspect job expansion without running simulations
.venv/bin/python -m gut_ibm_tools.batch_runner \
  experiments/diversity_campaign/batch_kd_sweep.json --dry-run

# Baseline × 3 seeds (3 runs)
.venv/bin/python -m gut_ibm_tools.batch_runner \
  experiments/diversity_campaign/batch_baseline.json

# Full mechanisms × 3 seeds (3 runs)
.venv/bin/python -m gut_ibm_tools.batch_runner \
  experiments/diversity_campaign/batch_full_mechanisms.json

# Four Kd values × 3 seeds (12 runs)
.venv/bin/python -m gut_ibm_tools.batch_runner \
  experiments/diversity_campaign/batch_kd_sweep.json
```

Total: 18 runs. The batch manifests write isolated outputs beneath
`batch_results/diversity_campaign/` and can be selected through
`./rebuild_and_run.sh`.

## Single run

```bash
# From the repository root
mpirun --bind-to none -np 1 build/gut_ibm \
  experiments/diversity_campaign/full_mechanisms.json

# Or isolate output in a run directory
mkdir -p runs/full_mechanisms
cd runs/full_mechanisms
mpirun --bind-to none -np 1 ../../build/gut_ibm \
  ../../experiments/diversity_campaign/full_mechanisms.json
```

Increase MPI ranks only after the one-rank run initializes within the available
WSL2 memory.

## GPU

To opt into a CUDA-enabled build, add these top-level keys to a copied
simulation config:

```json
"gpu_enabled": true,
"gpu_device_id": -1
```

`gpu_device_id: -1` lets ranks select devices by local rank. The full grid is
mirrored on the GPU, so verify both host RAM and VRAM capacity first.

## Analysis

```python
from gut_ibm_tools import GutIBMData, validation

with GutIBMData("full_mechanisms.h5") as data:
    spatial = validation.validate_spatial_signatures(data, data.steps[-1])
    genomic = validation.validate_genomic_signatures(data)
    print(f"Monochromatic score: {spatial['monochromatic_score']:.3f}")
    print(f"Resident retention: {genomic['resident_retention']:.1%}")
```

## Validation targets

- Resident retention: 70–80% after seven days
- Monochromatic patchiness: greater than 0.7
- BtuB/FepA-downregulated immigrants wash out (`mu < gamma_flow`)
