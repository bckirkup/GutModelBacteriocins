# GPU Acceleration

CUDA GPU acceleration for GutIBM (issue #33). OpenMP remains the CPU shared-memory path when CUDA is disabled or no device is available.

## Implemented Kernels

| File | Role |
|------|------|
| `greens_kernel.cu` | QSSA Green's function superposition (brute-force path) |
| `field_update_kernel.cu` | Reaction integration, boundaries, nutrient depletion, grid coupling |
| `agent_update_kernel.cu` | Metabolism Monod growth + biomass update |
| `chemistry_kernel.cu` | Agent O₂ respiration + VBF nutrient coupling |
| `spatial_hash_kernel.cu` | Parallel cell-key assignment (reserved; host hash used in production path) |
| `diffusion_kernel.cu` | PCR tridiagonal solver for backward-Euler directional diffusion (x/y periodic, z bounded) |

Host-side facades: `device.cpp`, `dispatch.cpp`, `chemistry_pipeline.cpp`, `greens_function_gpu.cpp`, `chemical_field_gpu.cpp`, `agent_pool_gpu.cpp`, `spatial_hash_gpu.cpp`, `qssa_gpu.cpp`, `vbf_gpu.cpp`, `diffusion_gpu.cpp`.

### Production chemistry path (Spec 9 PR5)

`run_chemistry_pipeline()` in `chemistry_pipeline.cpp` owns the GPU-first chemistry orchestration previously inlined in `Simulation::module_chemistry()`:

1. Agent O₂ depletion on device (or CPU fallback)
2. Single host download for `MPI_Allreduce` on reaction grids
3. VBF coupling on device when the O₂ GPU path succeeded (or CPU fallback + upload)
4. Reaction integration, implicit diffusion, and boundaries on device
5. One concentration sync back to host (or CPU fallback)

A dedicated `cudaStream_t` from `gpu_compute_stream()` serializes chemistry kernels without per-kernel device-wide syncs. Barnes-Hut FMM (`use_fmm=true`) still evaluates the bacteriocin far-field on CPU while the rest of the chemistry stack runs on GPU.

Parity and production smoke: `test_gpu_smoke`, `test_gpu_feature_combinations`, `test_gpu_production_path`, `test_mpi_gpu_multi_rank`, `scripts/compare_gpu_parity.sh`.

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
  ├─ GPU: FixMetabolism compute (Monod + biomass) when fur disabled
  ├─ CPU: QSSA bacteriocin (GF brute-force on GPU; FMM far-field on CPU)
  ├─ GPU: chemistry_pipeline (O₂, VBF, reactions, diffusion)
  ├─ CPU: mechanics, conjugation, mutation, receptor RNG
  └─ GPU: sync agents to device before MPI migration
```

Agent genomes (`bi_loci` variable length) remain on the host; SoA device buffers carry scalars and `bi_loci_count` for metabolism penalties.

## Memory

Default 1 mm × 1 mm × 100 µm grid at 2 µm resolution:
- ~12.5M cells × 6 species × (`conc` + `reac`) ≈ **1.2 GB** device memory for `ChemicalFieldGpu`
- Agent SoA scales ~O(N) with additional VRAM

Production example (`diversity_paradox` scale, 2 mm × 2 µm): nx=1000 fits the PCR line limit (1024).

## Tests

```bash
cd build && ctest -R 'greens_function_gpu|gpu_diffusion|gpu_chemical_field|gpu_feature_combinations|gpu_production_path|gpu_smoke' --output-on-failure
mpirun -np 2 ./test_mpi_gpu_multi_rank   # MPI + GPU chemistry checksum
bash scripts/compare_gpu_parity.sh
```

Tests skip the GPU execution path when no CUDA device is present (compile-only CI still builds all kernels).

## Fallback

All GPU entry points return `false` or no-op when:
1. Built without `GUTIBM_USE_CUDA`
2. `gpu_enabled` is false in config
3. No CUDA device is available
4. Grid line length exceeds the PCR shared-memory limit (nx/ny/nz > 1024) — diffusion falls back to CPU

The existing OpenMP/serial CPU implementations are used unchanged.

## Not Yet on GPU

- Barnes-Hut FMM octree traversal (`use_fmm=true` far-field stays on CPU)
- Fur iron regulation (`fur.enabled` disables GPU metabolism)
- Mechanics, conjugation, mutation, receptor kill RNG
- HDF5 checkpoint I/O
- Multi-GPU NCCL / device-side MPI reaction reduction
- OpenCL / HIP portability
