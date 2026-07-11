# GPU Acceleration

CUDA GPU acceleration for GutIBM (issue #33). OpenMP remains the CPU shared-memory path when CUDA is disabled or no device is available.

## Implemented Kernels

| File | Role |
|------|------|
| `greens_kernel.cu` | QSSA Green's function superposition (brute-force path) |
| `field_update_kernel.cu` | Reaction integration, boundaries, nutrient depletion, grid coupling |
| `agent_update_kernel.cu` | Metabolism Monod growth + biomass update |
| `chemistry_kernel.cu` | Agent O₂ respiration + VBF nutrient coupling |
| `spatial_hash_kernel.cu` | Parallel cell-key assignment for agent sorting |
| `diffusion_kernel.cu` | PCR tridiagonal solver for backward-Euler directional diffusion (x/y periodic, z bounded) |

Host-side facades: `device.cpp`, `dispatch.cpp`, `greens_function_gpu.cpp`, `chemical_field_gpu.cpp`, `agent_pool_gpu.cpp`, `spatial_hash_gpu.cpp`, `qssa_gpu.cpp`, `vbf_gpu.cpp`, `diffusion_gpu.cpp`.

### GPU diffusion (Spec 9)

`gpu_apply_species_diffusion()` mirrors `ChemicalField::apply_diffusion` for a single species:
periodic x/y via Sherman–Morrison PCR, epithelial Dirichlet + luminal Neumann z boundaries,
optional z-gradient background shift, and non-negativity clamp. `ChemicalFieldGpu::apply_diffusion`
and `apply_boundaries` run inside `Simulation::module_chemistry()` when GPU reactions succeed.
Agent O₂ respiration and VBF nutrient coupling also run on device on the GPU path (PR4).
Host concentrations sync once after the GPU pass (or after CPU fallback). Parity:
`test_gpu_diffusion`, `test_gpu_chemical_field`, `test_gpu_smoke`, `test_gpu_feature_combinations`.

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
  ├─ GPU: grid coupling, spatial hash cell keys
  ├─ GPU: FixMetabolism compute (Monod + biomass)
  ├─ GPU: QSSA GF superposition + nutrient depletion + field update
  ├─ CPU: VBF coupling, mechanics, conjugation, mutation, receptor RNG
  └─ GPU: sync back to host before MPI migration
```

Agent genomes (`bi_loci` variable length) remain on the host; SoA device buffers carry scalars and `bi_loci_count` for metabolism penalties.

## Memory

Default 1 mm × 1 mm × 100 µm grid at 2 µm resolution:
- ~12.5M cells × 6 species × (`conc` + `reac`) ≈ **1.2 GB** device memory for `ChemicalFieldGpu`
- Agent SoA scales ~O(N) with additional VRAM

## Tests

```bash
cd build && ctest -R 'greens_function_gpu|gpu_diffusion|gpu_chemical_field|gpu_feature_combinations|gpu_smoke' --output-on-failure
bash scripts/compare_gpu_parity.sh
```

Tests skip the GPU execution path when no CUDA device is present (compile-only CI still builds all kernels).

## Fallback

All GPU entry points return `false` or no-op when:
1. Built without `GUTIBM_USE_CUDA`
2. `gpu_enabled` is false in config
3. No CUDA device is available

The existing OpenMP/serial CPU implementations are used unchanged.

## Not Yet on GPU

- Barnes-Hut FMM (`use_fmm=true` near/far field stays on CPU)
- Mechanics, conjugation, mutation, receptor kill RNG
- HDF5 checkpoint I/O
- OpenCL / HIP portability
