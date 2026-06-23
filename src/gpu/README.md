# GPU Acceleration — Planned Architecture

This directory is a stub for the planned CUDA/OpenCL GPU acceleration of GutIBM.
The current implementation uses OpenMP threading as an intermediate parallelism step
(see issue #18). A full GPU port is tracked as a separate derivative issue.

## Planned Kernel Targets

### 1. Green's Function Superposition (`greens_kernel.cu`)
- Each CUDA thread evaluates one source-cell pair
- Grid: `(num_sources, num_grid_cells_in_cutoff)`
- Shared memory for source positions and GF parameters
- Expected speedup: 50-100x over serial (embarrassingly parallel)

### 2. Spatial Hash Build & Query (`spatial_hash_kernel.cu`)
- Parallel radix sort on cell keys for cache-coherent queries
- Atomic insert into hash buckets
- Neighbor queries: one thread per agent, shared-memory tile for bucket contents

### 3. Agent Update (`agent_update_kernel.cu`)
- Monod kinetics, growth, advection: one thread per agent
- Structure-of-Arrays (SoA) layout for coalesced memory access
- Separate kernels for read-only (growth rate) and write (state update) phases

### 4. Chemical Field Update (`field_update_kernel.cu`)
- Grid-parallel reaction application
- Boundary condition enforcement
- Trivially parallel: one thread per grid cell per species

## Memory Layout

Current AoS (Array of Structures) agent layout should be converted to SoA
for GPU coalescing:

```
// CPU (current):  Agent[] — each Agent has x, v, biomass, mu, ...
// GPU (planned):  x[], v[], biomass[], mu[], ...  (separate arrays)
```

## Build Integration

```cmake
option(GUTIBM_USE_CUDA "Enable CUDA GPU acceleration" OFF)
if(GUTIBM_USE_CUDA)
  enable_language(CUDA)
  add_definitions(-DGUTIBM_CUDA)
  # CUDA sources would be added here
endif()
```

## Prerequisites

- NVIDIA GPU with Compute Capability >= 6.0 (Pascal or newer)
- CUDA Toolkit >= 11.0
- Or: OpenCL 1.2+ runtime for vendor-neutral GPU support

## Fallback Strategy

The GPU kernels will be guarded by `#ifdef GUTIBM_CUDA` with automatic
fallback to the OpenMP-threaded CPU path when CUDA is unavailable.
This ensures the code remains portable across all platforms.
