---
name: gut-ibm-development
description: Build, test, and develop the GutIBM C++/MPI simulation. Covers CMake configuration, running tests, adding Fix modules, and using the Python analysis tools.
---

# GutIBM Development Skill

## Build

```bash
mkdir -p build && cd build
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON
make -j$(nproc)
```

Debug build:
```bash
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
```

## Run Tests

```bash
cd build && ctest --output-on-failure
```

## Run Simulation

```bash
cd build
mpirun -np 4 ./gut_ibm ../examples/<scenario>/config.yaml
```

Single-process (debugging):
```bash
./gut_ibm ../examples/<scenario>/config.yaml
```

## Python Analysis Tools

```bash
cd python && pip install -e ".[dev]"
```

Key analysis scripts in `python/gut_ibm_tools/`:
- Spatial signature analysis (colony morphology, comet-tail visualization)
- Genomic diversity metrics (BI locus allele frequencies)
- HDF5 checkpoint readers

## Adding a New Fix Module

1. Create `src/fixes/fix_<name>.h` and `src/fixes/fix_<name>.cpp`
2. Inherit from `Fix` base class
3. Implement `init()`, `pre_force()` or `post_force()`, and `end_of_step()` as needed
4. Register in `src/fixes/fix_registry.cpp`
5. Add to `CMakeLists.txt` source list
6. Add test in `tests/test_fix_<name>.cpp`
7. Add example configuration in `examples/`

## Adding a New Diffusion Kernel

1. Create kernel in `src/diffusion/`
2. Implement the Green's function solution with Method of Images for boundary conditions
3. Ensure QSSA compatibility (no explicit time-stepping)
4. Add unit test verifying against analytical solution
5. Wire into the field update in `src/fields/`

## Key Configuration Parameters

| Parameter | Description |
|-----------|-------------|
| `domain_size` | Physical dimensions of simulation box (μm) |
| `mucus_flow_rate` | Advective velocity in mucus layer (μm/s) |
| `dt_bio` | Biological timestep (typically 60s) |
| `n_initial_cells` | Starting population per strain |
| `mpi_decomposition` | Domain split across MPI ranks |
| `hdf5_checkpoint_interval` | Steps between HDF5 snapshots |

## Architecture Notes

- Chemistry is instantaneous (QSSA) — no CFL constraint
- Bio timestep (60s) is decoupled from chemistry
- Spatial hashing gives O(N) neighbor queries for millions of agents
- VBF (Viscoelastic Background Field) represents anaerobic majority as continuum
- All toxin fields use analytical advection-diffusion kernels, not grid PDE
- HDF5 I/O is parallel when MPI is enabled

## Spec Documents

- `EARI.md` — Eco-Advective Receptor Interference (the advective double-bind)
- `VADI.md` — Viscous Advective-Diffusion Interference (combinatorial washout trap)
