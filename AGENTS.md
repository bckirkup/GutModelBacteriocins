# AGENTS.md — AI Agent Guidelines for GutIBM

## Repository Purpose
Massively parallel 3D Individual-based Model (IbM) for simulating
*Enterobacteriaceae* population dynamics and genomic diversity in the colonic
mucus layer. C++17 core with Python analysis tools. Built on NUFEB-2
framework philosophy with MPI domain decomposition and HDF5 I/O.

## Setup
```bash
mkdir -p build && cd build
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON
make -j$(nproc)
```

Python analysis tools:
```bash
cd python && pip install -e ".[dev]"
```

## Validation Commands
Run these before committing:
```bash
# C++ build (with warnings)
cd build && cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON -DCMAKE_CXX_FLAGS="-Wall -Wextra" && make -j$(nproc)

# C++ tests
cd build && ctest --output-on-failure

# Python tools (if modified)
cd python && ruff check . && pytest
```

## Architecture Rules
- **NUFEB-2 Fix architecture** — biological rules are modular Fix plugins
- **QSSA for chemistry** — no explicit PDE solvers; analytical Green's function kernels
- **VBF for anaerobes** — 99% of background flora is a continuum field, not discrete agents
- **Spatial hashing** — O(N) neighbor lookups, not O(N²)
- **MPI domain decomposition** — code must be parallelizable across HPC nodes
- **Never modify tests to make them pass** — fix the implementation

## Key Files
| Path | Purpose |
|------|---------|
| `src/core/` | Simulation engine, domain management, spatial hashing |
| `src/fields/` | Spatial environment and chemical gradient management |
| `src/diffusion/` | Analytical advection-diffusion kernel implementations |
| `src/fixes/` | Biological rule modules (metabolism, lysis, conjugation) |
| `src/genome/` | Genomic state, BI locus evolution, mutation logic |
| `src/io/` | Parallel HDF5 data persistence and checkpointing |
| `python/gut_ibm_tools/` | Analysis suite for spatial and genomic signatures |
| `examples/` | Configuration templates for specific biological scenarios |
| `tests/` | C++ unit tests (CTest) |

## Key Concepts
- **Advective Double-Bind (EARI)**: Growth/resistance trade-off leading to washout
- **Combinatorial Washout Trap (VADI)**: Physical expulsion when μ < γ_flow
- **Lethal Core/Halo**: Spatial toxin zones based on isoelectric point (pI)
- **Comet-tails**: Downstream-elongated inhibitory zones from mucus flow
- **BI Locus**: Bacteriocin-Immunity genetic region subject to recombination
- **Method of Images**: Boundary condition technique for diffusion kernels
- **MPS (Mating-Pair Stabilization)**: Shear-dependent HGT probability

## Spec Documents
- `EARI.md` — Eco-Advective Receptor Interference framework
- `VADI.md` — Viscous Advective-Diffusion Interference framework

## Code Conventions
- C++17 with modern idioms (smart pointers, RAII, move semantics)
- CMake build system
- MPI for parallelism (guard all MPI calls with rank checks)
- HDF5 for I/O (parallel HDF5 when MPI is enabled)
- Python analysis code uses NumPy, SciPy, h5py

## PR Requirements
- Clean build with -Wall -Wextra (no new warnings)
- All CTest tests pass
- New features include tests
- MPI-safe (no deadlocks, proper collective operations)
- Update examples/ if adding new biological scenarios
