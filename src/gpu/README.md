# GPU Acceleration

CUDA GPU acceleration for GutIBM (issue #33). OpenMP remains the CPU shared-memory path when CUDA is disabled or no device is available.

## Implemented Kernels

| File | Role |
|------|------|
| `greens_kernel.cu` | QSSA Green's function superposition (brute-force path) |
| `field_update_kernel.cu` | Reaction integration, boundaries, nutrient depletion, grid coupling |
| `agent_update_kernel.cu` | Metabolism Monod growth, biomass update, Fur modulation, acetate overflow/scavenging, delivery prepare/commit, and funded respiration (`receptor_expr` / `ligand_affinity` are **SoA**: `r * num_agents + i`) |
| `chemistry_kernel.cu` | Agent O₂ respiration + VBF nutrient coupling (implicit carbon Monod, O₂, and iron sinks) |
| `spatial_hash_kernel.cu` | Parallel cell-key assignment + CSR neighbor lists |
| `diffusion_kernel.cu` | PCR tridiagonal solver for backward-Euler directional diffusion (x/y periodic, z bounded) + Route-B device delivery |
| `fmm_far_kernel.cu` | FMM local-expansion far-field grid deposit (when M2L locals ready) |
| `receptor_kernel.cu` | Receptor competitive-binding kill probability on device |
| `mechanics_kernel.cu` | Overdamped Stokes Hertzian repulsion + EPS adhesion via CSR neighbor lists (including shared host/device CDI corpse participation predicate) |

Host-side facades: `device.cpp`, `dispatch.cpp`, `chemistry_pipeline.cpp`, `greens_function_gpu.cpp`, `chemical_field_gpu.cpp`, `agent_pool_gpu.cpp`, `spatial_hash_gpu.cpp`, `fmm_gpu.cpp`, `receptor_gpu.cpp`, `mechanics_gpu.cpp`, `qssa_gpu.cpp`, `vbf_gpu.cpp`, `diffusion_gpu.cpp`, `gpu_profile.cpp`.

### Production chemistry path (Spec 9 PR5)

`run_chemistry_pipeline()` in `chemistry_pipeline.cpp` owns the GPU-first chemistry orchestration previously inlined in `Simulation::module_chemistry()`:

1. Agent O₂ depletion on device (or CPU fallback)
2. Single host download for `MPI_Allreduce` on reaction grids (or device Allreduce when CUDA-aware MPI is active)
3. VBF coupling on device when the O₂ GPU path succeeded (or CPU fallback + upload)
4. Reaction integration, implicit diffusion (including Route-B device delivery for replicated single-rank lines ≤ 512), and boundaries on device
5. Per-species host/device residency tracking minimizes unnecessary full-grid syncs (see `docs/LANDMINES_ARCHIVE.md`)

A dedicated `cudaStream_t` from `gpu_compute_stream()` serializes chemistry kernels without per-kernel device-wide syncs. QSSA near-field superposition runs on GPU via `gpu_superpose_to_device`; FMM far-field uses `fmm_far_kernel.cu` when dense M2L locals are ready (trees with >256 nodes fall back to CPU `traverse_far`).

Chemical-field device mirrors use global y/z extents with a local x storage
stride. In slab mode, only owned x cells are updated or included in scalar
accounting; concentration halos are refreshed separately and never contribute
to boundary, reaction-clip, uptake, or VBF totals. Periodic x diffusion uses
the exact host slab transpose path before local y/z device operators. Host
slab halos are exchanged before each mirror upload so agent biology reads
current halo copies; halo cells remain excluded from scalar totals. The
host/device round trip is accumulated as `gpu_slab_x_roundtrip_s` in
`StepProfile` when profiling is enabled. `StepProfile` also records
`gpu_h2d_s`, `gpu_d2h_s`, `mpi_reaction_reduce_s`, and `hdf5_s`. Use
`scripts/run_gpu_scaling_benchmark.sh` for CPU vs GPU sweeps.

### CUDA-aware MPI reaction reduce (#156)

When `GUTIBM_CUDA_AWARE_MPI=1` **and** the linked MPI library reports CUDA buffer support (`MPIX_Query_cuda_support` on Open MPI), `ChemicalFieldGpu::try_sum_reactions_on_device()` performs `MPI_Allreduce` on device reaction buffers and skips the host round-trip. Without CUDA-aware MPI, the chemistry pipeline falls back to host `sum_reactions_across_ranks()` automatically.

| Variable | When to set |
|----------|-------------|
| `GUTIBM_CUDA_AWARE_MPI=1` | Opt in on clusters with CUDA-aware Open MPI / MVAPICH2-GDR |
| `GUTIBM_CUDA_AWARE_MPI_FORCE=1` | Manual override when auto-detection is unavailable but buffers are known CUDA-aware |

Validation: `mpirun -np 2 ./test_cuda_aware_mpi_reaction` (device parity runs only when CUDA-aware MPI is detected). Broader MPI smoke: `bash scripts/run_mpi_scaling_smoke.sh`.

Parity and production smoke: 23 `gpu`-labelled CTest targets (`ctest -L gpu`), including `gpu_smoke`, `qssa_gpu_parity`, `gpu_scaling_benchmark`, `spatial_hash_gpu_csr`, `mechanics_gpu_parity`, `gpu_feature_combinations`, `gpu_production_path`, `gpu_uptake_limit`, `gpu_residency_audit`, `gpu_reaction_residency`, `mpi_gpu_multi_rank`, `cuda_aware_mpi_reaction`, plus `scripts/compare_gpu_parity.sh`.

## Build

```bash
cmake -B build -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON -DGUTIBM_USE_CUDA=ON
cmake --build build -j$(nproc)
```

Prerequisites:
- NVIDIA GPU with Compute Capability >= 6.0 (Pascal+)
- CUDA Toolkit >= 11.0 (or `nvidia-cuda-toolkit` on Ubuntu)
- NVIDIA driver for runtime execution

## Runtime Configuration

Input file keys:
```
gpu_enabled true
gpu_device_id 0    # omit or -1 for auto (MPI rank % device count)
```

C++ API:
```cpp
cfg.gpu.enabled = true;
cfg.gpu.device_id = 0;  // -1 = auto
```

## Architecture

```
Simulation::step()
  ├─ GPU: sync agents/fields to device after MPI ghost exchange
  ├─ GPU: grid coupling (async stream)
  ├─ CPU: spatial hash rebuild (host)
  ├─ GPU: FixMetabolism compute (Monod + biomass + O₂ boost + Fur modulation +
  │      acetate overflow/scavenging + delivery prepare/commit + funded respiration);
  │      siderophore field chemistry host-side
  ├─ GPU: QSSA bacteriocin near-field (`gpu_superpose_to_device`); FMM far-field GPU when locals ready
  ├─ GPU: chemistry_pipeline (O₂, VBF, reactions, diffusion / Route-B delivery)
  ├─ GPU: receptor kill-prob kernel (Bernoulli kills on host)
  ├─ GPU: mechanics (overdamped Stokes Hertzian + adhesion via spatial_hash_gpu CSR,
  │      including CDI corpse participation via shared host/device predicate)
  └─ GPU: sync agents to device before MPI migration
```

Agent genomes (`bi_loci` variable length) remain on the host; SoA device buffers carry scalars and `bi_loci_count` for metabolism penalties. CDI corpse mechanics participation uses the shared host/device predicate in `src/fixes/mechanics_participation.h`.

## Memory

Default 1 mm × 1 mm × 100 µm grid at 2 µm resolution:
- ~12.5M cells × 6 species × (`conc` + `reac`) ≈ **1.2 GB** device memory for `ChemicalFieldGpu`
- Agent SoA scales ~O(N) with additional VRAM

Production example (`diversity_paradox` scale, 2 mm × 2 µm): nx=1000 fits the standard PCR line limit (1024); Route-B device delivery supports line lengths up to 512 (longer lines fall back to host-forced delivery recorded in `/run_provenance/chemistry_placement`).

## Tests

```bash
cd build && ctest -L gpu --output-on-failure
mpirun -np 2 ./tests/test_mpi_gpu_multi_rank   # MPI + GPU chemistry checksum
mpirun -np 2 ./tests/test_cuda_aware_mpi_reaction  # CUDA-aware MPI reaction reduce (#156)
mpirun -np 4 --oversubscribe ./tests/test_mpi_four_rank  # 4-rank MPI smoke (#154)
bash scripts/compare_gpu_parity.sh
bash scripts/run_mpi_scaling_smoke.sh
```

Tests return CTest skip code 77 when no CUDA device is present (compile-only CI still builds all kernels); set `REQUIRE_GPU=1` to fail on skip.

## Fallback

All GPU entry points return `false` or no-op when:
1. Built without `GUTIBM_USE_CUDA`
2. `gpu_enabled` is false in config
3. No CUDA device is available
4. Grid line length exceeds the PCR shared-memory limit (nx/ny/nz > 1024 for standard diffusion, or > 512 for Route-B device delivery) — diffusion falls back to CPU

The existing OpenMP/serial CPU implementations are used unchanged.

## Not Yet on GPU

- Barnes-Hut FMM octree traversal (`use_fmm=true` far-field M2L is CPU tree-walk; L-eval on GPU when locals ready; #159)
- Siderophore secretion, chelation, and ferric-enterobactin reimport remain deliberately host-side because they aggregate per-cell field chemistry (per-agent Fur receptor modulation and acetate overflow/scavenging run in the metabolism kernel)
- Conjugation, mutation, receptor kill RNG
- HDF5 checkpoint I/O
- Multi-GPU NCCL
- OpenCL / HIP portability
