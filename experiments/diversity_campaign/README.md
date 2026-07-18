# Diversity Campaign v3 — Three-Stage Experimental Plan

Configs live under `experiments/diversity_campaign/` and are discovered by
`./rebuild_and_run.sh` in stage order (stage1 → stage2 → stage3).

## Stage overview

| Stage | Folder | Scale | Purpose |
|-------|--------|-------|---------|
| 1 | `stage1_motility_validation/` | 500 µm, 1 h, ~100 agents | Validate motility / taxis modes |
| 2 | `stage2_mechanism_validation/` | 500 µm, 1 h, motility OFF | Validate O₂, acetate, crypts, CDI |
| 3 | `stage3_campaign/` | 2 mm, 7 d, full grid | Science runs + Kd sweep |

Run Stages 1 and 2 before Stage 3. Stage 3 configs are production-scale (see
resource warning below).

## Using `rebuild_and_run.sh`

```bash
# Interactive: pick a single config or batch (files grouped by stage folder)
./rebuild_and_run.sh

# Or after build: reuse and open the experiment menu
./rebuild_and_run.sh --reuse-build

# Menu option "Run a diversity-campaign stage" executes every single-run
# JSON in that stage folder in sorted order (1a→1f, 2a→2f, …).
# Each successful run's HDF5 is gzipped to *.h5.gz by default.

# Direct single / batch (non-interactive)
./rebuild_and_run.sh --reuse-build --mode single \
  --config experiments/diversity_campaign/stage1_motility_validation/1a_motility_off.json \
  --mpi-ranks 1

./rebuild_and_run.sh --reuse-build --mode batch \
  --config experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json \
  --batch-action dry-run

# Leave HDF5 uncompressed (e.g. for immediate analysis)
./rebuild_and_run.sh --reuse-build --no-gzip-hdf5 --mode stage \
  --config experiments/diversity_campaign/stage1_motility_validation
```

## Stage 1: Motility validation (1-hour runs, small domain)

| Config | Motility | O₂ | Aerotaxis | Energy taxis | Carbon chemo | AI-2 | Expected |
|--------|----------|-----|-----------|--------------|--------------|------|----------|
| `1a_motility_off` | OFF | — | — | — | — | — | Stable population (reference) |
| `1b_motility_blind` | ON blind | — | OFF | OFF | — | — | Population collapse (~2 agents) |
| `1c_aerotaxis_only` | ON | ON | ON | OFF | — | — | Directional bias, partial recovery |
| `1d_aerotaxis_energy` | ON | ON | ON | ON | — | — | Bias + speed reduction → stable |
| `1e_full_taxis` | ON | ON | ON | ON | ON | — | Strongest epithelial bias |
| `1f_full_taxis_ai2` | ON | ON | ON | ON | ON | ON | Bias + clustering |

Each run: ~500 µm domain, 100 agents, 1 hour, dt = 60 s, dx = 5 µm.
Expected wall time: ~1–5 min each on a laptop.

## Stage 2: Mechanism validation (1-hour runs, motility OFF)

| Config | O₂ | Acetate | Crypts | Adaptive dt | CDI | Expected |
|--------|-----|---------|--------|-------------|-----|----------|
| `2a_baseline` | — | — | — | — | — | Baseline reference |
| `2b_o2_only` | ON | — | — | — | — | z-dependent µ (aerobic boost) |
| `2c_o2_acetate` | ON | ON | — | — | — | MetE penalty for BtuB-low agents |
| `2d_o2_acetate_crypts` | ON | ON | ON | — | — | Crypt agents survive washout |
| `2e_full_mechanisms` | ON | ON | ON | ON | — | Adaptive dt varies; population bounded |
| `2f_full_mech_cdi` | ON | ON | ON | — | types | CDI kills at strain boundaries |

## Stage 3: Full campaign (7-day runs, full domain)

| Config | Description | Motility | Kd |
|--------|-------------|----------|-----|
| `3a_baseline` | Baseline | OFF | default (1 nM) |
| `3b_full_mechanisms` | Full mechanisms | OFF | default (1 nM) |
| `3c_full_mech_motile` | Full mechanisms + motile | ON (aerotaxis+energy) | default (1 nM) |
| `3_kd_sweep_1e-9` … `1e-6` | Kd sweep singles | OFF | 1 nM → 1 µM |

### Batch manifests (Stage 3)

```bash
# Baseline × 3 seeds
./rebuild_and_run.sh --cuda on --reuse-build --mode batch \
  --config experiments/diversity_campaign/stage3_campaign/batch_baseline.json

# Full mechanisms × 3 seeds
./rebuild_and_run.sh --cuda on --reuse-build --mode batch \
  --config experiments/diversity_campaign/stage3_campaign/batch_full_mechanisms.json

# Full mechanisms + motility × 3 seeds (after Stage 1 passes)
./rebuild_and_run.sh --cuda on --reuse-build --mode batch \
  --config experiments/diversity_campaign/stage3_campaign/batch_full_mech_motile.json

# Four Kd values × 3 seeds (12 runs)
./rebuild_and_run.sh --cuda on --reuse-build --mode batch \
  --config experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json \
  --batch-action dry-run
```

Stage 3 batches use **1 MPI rank** by default (full-grid memory per rank).
Outputs land under `batch_results/diversity_campaign/`.

## Decision points

After Stage 1:
- If `1d` sustains a viable population → consider motility ON for Stage 3 (`3c`)
- If `1d` still collapses → keep motility OFF for Stage 3; investigate parameters

After Stage 2:
- If `2e` works → proceed to Stage 3 with full mechanisms
- If any mechanism crashes or yields NaN → fix before Stage 3

After Stage 3:
- Analyze Kd sweep; compare motile vs non-motile (`3c` vs `3b`) if Stage 1 passed

## Resource warning (Stage 3)

Stage 3 uses a 2 mm × 2 mm × 100 µm domain at 2 µm grid spacing
(`1000 × 1000 × 50 = 50,000,000` cells) and ~10 chemical species
(O₂ + acetate + defaults). **Every MPI rank stores the full grid** (no
slab-local chemistry yet).

### Exact chemistry memory (dominant term)

```
bytes ≈ ncells × nspecies × 2 (conc+reac) × 8
      = 50e6 × 10 × 2 × 8
      = 8.0 GB host   (ChemicalField)
      + 8.0 GB VRAM   (ChemicalFieldGpu mirrors, when gpu_enabled)
```

Agents are cheap at campaign scale (~600 start; even 10⁵ is ≪ 1 GB). Transient
GPU sync copies allocate another ~0.4 GB host temporary per species. FMM/QSSA
scratch adds ~0.4–1 GB.

| Setup | Host chem | Device chem | Practical minimum |
|-------|-----------|-------------|-------------------|
| 1 rank, CPU | ~8 GB (+~3 GB overhead) | — | **≥16 GB** WSL RAM |
| 1 rank, GPU | ~8 GB (+ overhead/temps) | **8 GB** (+ scratch) | **≥16–24 GB** WSL RAM **and ≥12 GB** VRAM (16 GB safer) |
| 4 ranks, GPU | ~32 GB+ aggregate | 8 GB × devices | workstation / HPC only |

| GPU VRAM | Full Stage 3 (`dx=2 µm`) |
|----------|--------------------------|
| 6–8 GB | **Will not fit** (chem mirrors alone = 8 GB) |
| 10 GB | Marginal — likely crash after init/early steps |
| 12 GB | Possible if little else on the GPU |
| 16 GB+ | Recommended for GPU Stage 3 |

If a run prints `Step 1` (or later) then the terminal/WSL session dies with no
tidy traceback, treat it as **OOM** (host or VRAM), not a GutIBM logic error.
Check `dmesg`, WSL logs, and `nvidia-smi` for killed processes / Xid errors.

Batch manifests default to `"mpi_ranks": 1`. Raise the WSL cap in
`%UserProfile%\.wslconfig` (see `docs/WSL2_SETUP.md`) before trying more ranks.

### Laptop-safe grid (if full Stage 3 OOMs)

Chemistry memory scales as `1/dx³` (and as `L²` in x/y). Same 2 mm domain:

| `grid_dx` | Cells | Host chem | GPU chem |
|-----------|------:|----------:|---------:|
| 2 µm (current) | 50.0 M | 8.00 GB | 8.00 GB |
| 4 µm | 6.2 M | 1.00 GB | 1.00 GB |
| 5 µm | 3.2 M | 0.51 GB | 0.51 GB |

Or keep `dx=2 µm` but shrink to 1 mm × 1 mm → 12.5 M cells → **2 GB** chem
each on host/device.

**GPU is on by default for Stage 3** (`gpu_enabled: true`, `gpu_device_id: -1`).
Stages 1–2 stay on CPU — the 200k-cell grids are too small for GPU transfer
overhead to pay off.

Build a CUDA binary before launching Stage 3:

```bash
./rebuild_and_run.sh --cuda on --reuse-build --mode batch \
  --config experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json \
  --batch-action dry-run
```

If `nvcc` / a GPU is missing, `rebuild_and_run.sh` warns when a config requests
GPU but the binary was built without CUDA. Set `"gpu_enabled": false` only as a
CPU fallback (still needs ~8+ GB host RAM at full resolution).

## Validation targets (Stage 3)

- Resident retention: 70–80% after seven days
- Monochromatic patchiness: greater than 0.7
- BtuB/FepA-downregulated immigrants wash out (`mu < gamma_flow`)
- Kd sweep: higher Kd → colicin E more potent → higher resident retention
