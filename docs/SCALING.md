# Scaling Benchmarks and Profiling

GutIBM targets **10⁶–10⁷** discrete *Enterobacteriaceae* agents. This document
records how to measure performance, interpret results, and choose runtime
parameters for large runs.

## Quick start

```bash
# CI-scale smoke (500 / 1.5k / 3k agents, 1 step, brute-force QSSA)
cd build && ctest -R scaling_benchmark --output-on-failure

# Manual sweep (configurable agent counts and MPI ranks)
bash scripts/run_scaling_benchmark.sh

# Single large run with hotspot profiling
cd build
mpirun -np 4 ./gut_ibm ../examples/scaling_benchmark/input_1e6.json
```

Example configs live in `examples/scaling_benchmark/` (`input_1e5.json`,
`input_1e6.json`).

## Benchmark tooling

| Artifact | Purpose |
|----------|---------|
| `tests/test_scaling_benchmark.cpp` | Fast smoke test; prints `BENCHMARK ...` lines |
| `scripts/run_scaling_benchmark.sh` | Sweeps agent counts × MPI ranks → CSV |
| `profile_steps true` | Per-step hotspot breakdown in `gut_ibm` stdout |

### Environment variables (`run_scaling_benchmark.sh`)

| Variable | Default | Meaning |
|----------|---------|---------|
| `AGENT_COUNTS` | `10000 50000 100000 500000` | Space-separated population sizes |
| `MPI_RANKS` | `1 2 4` | MPI process counts |
| `STEPS` | `5` | Biological steps per run |
| `USE_FMM_THRESHOLD` | `50000` | Enable Barnes-Hut FMM at or above this count |
| `OUT_DIR` | `benchmark_results/` | CSV and log output directory |

Results CSV columns: `timestamp, agents, ranks, use_fmm, wall_s, step_ms,
chemistry_ms, biology_ms, hash_ms, rss_mb`.

### Profiling output

When `profile_steps` is enabled, rank 0 prints mean wall time per phase after
`run()` completes:

```
Step profile (mean wall time per step, s):
  ghost_exchange=...
  spatial_hash=...
  biology=...
  chemistry=...
  physics=...
  mpi_migrate=...
  cleanup=...
PROFILE_CSV steps=... ghost_s=... hash_s=... biology_s=... chemistry_s=...
```

**Typical hotspots by scale:**

| Scale | Dominant cost | Mitigation |
|-------|---------------|------------|
| &lt; 10⁵ agents | Biology (Fix compute, mechanics) | Reduce enabled fixes; coarser `grid_dx` |
| 10⁵–10⁶ | QSSA chemistry (Green's superposition) | `use_fmm true`, widen `toxin_cutoff` carefully |
| 10⁶+ | Chemistry + MPI ghost exchange | Barnes-Hut FMM + ≥4 MPI ranks; disable HDF5 during tuning |
| Multi-rank | `mpi_migrate`, ghost exchange | 1D slab along x; keep agents well mixed in x |

## Memory budget

### Chemical grid (fixed per domain)

Default domain 1 mm × 1 mm × 100 µm at `grid_dx = 2 µm`:

```
nx = 500, ny = 500, nz = 50  →  12.5 × 10⁶ cells
```

Six species × (`conc` + `reac`) × 8 bytes ≈ **1.2 GB** host memory, independent
of agent count. Doubling domain linearly increases grid memory.

### Per-agent host memory (AoS layout)

Agents are stored as an **array of structs** (`std::vector<Agent>`). Each
`Agent` holds scalars inline; `Genome::bi_loci` is a heap-allocated
`std::vector` (typically 1–2 BI clusters per resident strain).

| Component | Approximate size |
|-----------|------------------|
| `Agent` scalars + `Genome` fixed fields | ~400–500 B |
| Per BI locus (`BICluster`) | ~80 B |
| Vector / allocator overhead | ~50–100 B per agent |

**Rule of thumb:** ~0.5–1.5 KB per agent with default plasmids (`ColE1`+`ColB`).
At 10⁷ agents expect **5–15 GB** for the agent pool alone, before spatial hash,
lineage tracker, and MPI buffers.

Measure on your machine:

```bash
# VmRSS after init, from benchmark test output:
BENCHMARK ... rss_mb=... bytes_per_agent=...
```

### GPU memory (SoA layout)

When CUDA is enabled (`gpu_enabled true`), scalar agent fields and the chemical
grid are mirrored in **structure-of-arrays** device buffers for coalesced kernel
access. See `src/gpu/README.md`.

| Layout | Where | Strength | Weakness |
|--------|-------|----------|----------|
| **AoS** (host) | `AgentPool` | Simple Fix logic; variable `bi_loci` | Poor GPU coalescing; pointer chasing |
| **SoA** (device) | `AgentPoolGpu` | Fast metabolism / GF kernels | Genome still on host; sync each step |

Genomes (`bi_loci`, mutation state) remain on the host by design — full SoA for
variable-length genomes would require indirection tables or a fixed-capacity
locus array. The hybrid approach keeps HGT/mutation on CPU while offloading
Monod growth and QSSA field updates.

**VRAM estimate:** grid (~1.2 GB) + agent SoA (~32 B × N scalars) + hash buffers.
At 10⁶ agents, plan for **≥ 2 GB** device memory with default grid.

## Recommended parameters by target scale

| Target agents | `use_fmm` | `fmm_theta` | `toxin_cutoff` | MPI ranks | Notes |
|---------------|-----------|-------------|----------------|-----------|-------|
| ≤ 10⁴ | `false` | — | 200 µm | 1 | Brute-force QSSA is fine |
| 10⁵ | `true` | 0.5 | 200 µm | 1–2 | Enable FMM before chemistry dominates |
| 10⁶ | `true` | 0.5 | 200 µm | 4–8 | Disable HDF5 (`hdf5_every 0`) while tuning |
| 10⁷ | `true` | 0.3–0.5 | 150–200 µm | 16+ | May need coarser `grid_dx` or smaller domain |

Additional tuning:

- Set `hdf5_every 0` or `cfg.hdf5.enabled = false` during performance runs.
- `adaptive_dt_enabled false` gives stable step times for benchmarking.
- `profile_steps true` for hotspot identification; adds minor timer overhead.
- Barnes-Hut FMM is the supported far-field path; true higher-order FMM (#29) is
  not required for these benchmarks.

## Strong vs weak scaling

**Strong scaling:** fixed global agent count, increase MPI ranks. Expect diminishing
returns once ghost layers and migration exceed compute — the 1D x-slab decomposition
limits surface-to-volume improvement.

**Weak scaling:** grow agent count proportionally with ranks (each rank holds
≈ N/P agents). Ideal chemistry time stays flat if FMM + spatial hash remain O(N).

The benchmark script records both dimensions in CSV form for offline plotting.

## CI integration

- **Pull requests:** `scaling_benchmark` CTest (label `benchmark;slow`) runs 500 /
  1.5k / 3k agents (brute-force QSSA) on a compact domain. Barnes-Hut FMM sweeps
  are exercised via `run_scaling_benchmark.sh` (nightly workflow).
- **Nightly:** `.github/workflows/benchmark-nightly.yml` runs a larger sweep on
  `ubuntu-latest` and uploads `benchmark_results/*.csv` as an artifact.

## Known limits (v0.1.0)

- All CTest targets still default to `nprocs=1`; MPI scaling requires
  `run_scaling_benchmark.sh` or manual `mpirun`.
- GPU Barnes-Hut FMM is not offloaded — `use_fmm` far-field stays on CPU.
- Agent memory includes per-agent heap allocations; pool allocators (#future) could
  reduce `bytes_per_agent` by ~30%.
- Grid memory is the floor cost; shrinking `domain_x/y` reduces grid cells linearly.

## References

- `docs/PARAMETERS.md` — QSSA / FMM parameters
- `src/gpu/README.md` — SoA GPU buffers and VRAM
- `examples/diversity_paradox/input.json` — production-scale FMM example
- Issue [#55](https://github.com/bckirkup/GutModelBacteriocins/issues/55)
