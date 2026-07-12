# GutIBM — 3D Individual-based Model for Enterobacteriaceae Gut Dynamics

A massively parallel, 3D Individual-based Model (IbM) built to solve the
*Enterobacteriaceae* "diversity paradox" in the mammalian gut. Despite spatial
fragmentation that should theoretically support high diversity, empirical data
shows these populations maintain limited coexistence and monochromatic strain
clustering.

GutIBM integrates mucus fluid dynamics, bacteriocin (toxin) diffusion, and
metabolic trade-offs governed by TonB-dependent transporters (TBDTs), synthesizing
the **EARI** (Eco-Advective Receptor Interference) and **VADI** (Viscous
Advective-Diffusion Interference) theoretical frameworks.

**Version:** 0.1.0 (research prototype — see [AGENTS.md](AGENTS.md) for HPC/MPI scaling notes)

## Architecture

Built upon the **NUFEB-2** (Newcastle University Frontiers in Engineering Biology)
framework philosophy — a LAMMPS-inspired modular `Fix` architecture for biological
rule execution, with two-way coupled CFD-DEM capabilities.

### Key Design Decisions

| Challenge | NUFEB Approach | GutIBM Approach |
|---|---|---|
| Diffusion solver | Explicit FTCS grid PDE | **QSSA + Green's function kernels** (no CFL constraint) |
| Background flora | Discrete agents | **Viscoelastic Background Field** (VBF) continuum |
| Timescale coupling | Single dt | **Decoupled**: bio (60s) → chemistry (instantaneous QSSA) → physics (60s) |
| Spatial lookups | Grid-based | **Spatial hashing** (O(N) neighbor queries) |
| Toxin fields | Grid diffusion | **Analytical advection-diffusion kernels** with Method of Images |
| Large N toxin sums | O(N×M) pairwise | **Barnes–Hut FMM** optional far-field (`use_fmm`) |

### Computational Scaling

The simulation targets **10⁶–10⁷** discrete cells by:

1. **Avoiding explicit FTCS solvers** — QSSA Green's function kernels compute steady-state fields analytically.
2. **Spatial hashing** — O(N) neighbor lookups instead of O(N²) pairwise.
3. **VBF abstraction** — 99% of anaerobic microbiota are a continuous medium, not discrete agents.
4. **MPI domain decomposition** — 1D slab along x-axis with ghost exchange and agent migration.
5. **Optional CUDA** — GPU kernels for Green's function superposition and field updates ([`src/gpu/README.md`](src/gpu/README.md)).

See [docs/SCALING.md](docs/SCALING.md) for benchmark scripts, memory budgets, and recommended `use_fmm` / MPI settings.

### The Advective Double-Bind (EARI) / Combinatorial Washout Trap (VADI)

Immigrant strains entering a resident's bacteriocin "comet-tail" face an
evolutionary trap:

```
Immigrant enters comet-tail
    │
    ├─► Keep receptors (BtuB/FepA) → KILLED by colicin
    │
    └─► Downregulate receptors → SURVIVE toxin
            │
            └─► Increased Km → reduced μ_realized
                    │
                    └─► μ_realized < γ_flow → WASHED OUT
```

Crypt refugia (`crypts_enabled`) and partial receptor missense mutations provide escape routes — see [docs/MECHANISMS.md](docs/MECHANISMS.md).

## Implemented Fix Modules

| Fix | Purpose |
|---|---|
| `fix_metabolism` | Triple Monod kinetics, graded iron uptake (FepA/IroN/IutA/Fiu), Fur regulation, acetate–MetE penalty, division, death |
| `fix_bacteriocin` | SOS / phage / continuous microcin release (Spec 2), pI-dependent diffusion classes |
| `fix_receptor` | Competitive TBDT binding, immunity escape affinity, partial missense resistance |
| `fix_motility` | Mucus-adapted run-and-reverse swimming, chemotaxis (Spec 3) |
| `fix_conjugation` | F-pili HGT, shear-dependent MPS, pili length heterogeneity |
| `fix_cdi` | Contact-dependent inhibition, corpse barriers (Spec 3) |
| `fix_mutation` | BI duplication/recombination, receptor downregulation, super-killers, compensatory amelioration |
| `fix_mechanics` | Hertzian soft-sphere repulsion, optional EPS adhesion |

Select plugins with the `fixes` JSON array or enable subsystems via dot keys (`oxygen.enabled`, `fur.enabled`, …). Full reference: [docs/PARAMETERS.md](docs/PARAMETERS.md).

### Chemical environment (Spec 1)

Optional continuum fields coupled to metabolism and QSSA:

| Subsystem | Key toggle | Effect |
|-----------|------------|--------|
| Oxygen | `oxygen.enabled` | Epithelial O₂ gradient, aerobic growth boost, ROS→SOS coupling |
| Acetate | `acetate.enabled` | Dynamic fermentation/scavenging; MetE penalty amplifier |
| Mucin | `mucin.enabled` | Glycoprotein field, goblet secretion, liberation kinetics |
| Protease | `protease.enabled` | First-order colicin decay (per-plasmid half-lives) |
| z-gradients | `carbon_z_gradient`, … | Epithelium-to-lumen nutrient profiles |

### Toxin Biophysics

Bacteriocin diffusion is calibrated by isoelectric point (pI):

- **Lethal Cores** (pI > 8.5): Basic colicins bind mucin → high retardation
- **Lethal Halos** (pI < 7.0): Acidic colicins spread widely downstream
- **Comet Tails**: Mucus flow distorts radial diffusion into elongated exclusion zones
- **Peristaltic mixing**: Oscillatory flow modulation (`peristaltic_*` keys)

### Plasmid Library

| Name | pI | Class | Receptor | Release mode |
|---|---|---|---|---|
| Colicin E1 | 9.0 | Lethal Core | BtuB | SOS lysis |
| Colicin E2 | 6.5* | Lethal Halo | BtuB | SOS lysis |
| Colicin B | 5.4 | Lethal Halo | FepA | Phage lysis |
| Colicin Ia | 7.2 | Neutral | CirA | Phage lysis |
| Colicin M | 9.3 | Lethal Core | FhuA | SOS lysis |
| Microcin V | 5.0 | Lethal Halo | CirA | Continuous secretion |

*\*ColE2: secreted as equimolar complex with acidic Im2 protein*

## Building

### Prerequisites

```bash
sudo apt-get install cmake libopenmpi-dev openmpi-bin libhdf5-mpi-dev
```

For CUDA, use the toolkit version appropriate to the host. WSL2 users must
follow [`docs/WSL2_SETUP.md`](docs/WSL2_SETUP.md) and install the WSL-specific
toolkit without a Linux display driver.

### Compile

```bash
mkdir build && cd build
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON
make -j$(nproc)
```

CUDA build: add `-DGUTIBM_USE_CUDA=ON`. OpenMP: `-DGUTIBM_USE_OPENMP=ON`.

### Run Tests

```bash
cd build
ctest --output-on-failure          # full suite
ctest -L unit -LE slow             # fast CI shard
```

### Clean rebuild and interactive experiments

From the repository root:

```bash
./rebuild_and_run.sh
```

The helper safely removes `build/`, configures MPI/HDF5 with automatic CUDA
detection, creates or refreshes `.venv`, rebuilds, runs the full CTest suite,
then prompts for single or batch JSON files in
[`experiments/`](experiments/README.md). See
[`docs/WSL2_SETUP.md`](docs/WSL2_SETUP.md) for NVIDIA toolkit, MPI, HDF5, and
WSL memory guidance.

## Usage

### Single simulation

```bash
# Strict JSON config — see docs/CONFIG_FORMAT.md
./gut_ibm ../examples/diversity_paradox/input.json

# MPI parallel
mpirun -np 8 ./gut_ibm ../examples/diversity_paradox/input.json

# Checkpoint restart
./gut_ibm ../examples/my_run/input.json   # with checkpoint_file in JSON
```

Strict config parsing: set `GUTIBM_STRICT_CONFIG=1` to abort on invalid numerics.

### Batch parameter scans

Resumable multi-run CLI for sweeps and experiment grids:

```bash
cd python && pip install -e ".[dev]"

python -m gut_ibm_tools.batch_runner ../examples/batch_scan/batch.json --dry-run
python -m gut_ibm_tools.batch_runner ../examples/batch_scan/batch.json
python -m gut_ibm_tools.batch_runner ../examples/batch_scan/batch.json --resume
```

Each job writes to `{output_dir}/jobs/{job_id}/` with its own config, HDF5, and log.
Progress is tracked in `batch_manifest.json` — safe to interrupt with `Ctrl+C`.

See [docs/BATCH_RUNNER.md](docs/BATCH_RUNNER.md) for sweep vs explicit-run JSON schemas.

## Examples

| Directory | Description |
|-----------|-------------|
| `examples/single_colony/` | Comet-tail formation, peristaltic mixing |
| `examples/diversity_paradox/` | Full 7-day resident vs immigrant EARI/VADI run |
| `examples/eari_vadi_validation/` | Short CI regression scenario (15 min sim) |
| `examples/cell_biology/` | Fur, CDI, motility (Spec 3) |
| `examples/scaling_benchmark/` | Large-agent-count templates (10⁵–10⁶) |
| `examples/batch_scan/` | Batch runner sweep manifests |

## Python Analysis Toolkit

```bash
cd python
pip install -e ".[viz]"    # analysis + optional matplotlib
pip install -e ".[dev]"   # adds pytest, ruff
```

```python
from gut_ibm_tools import GutIBMData, validation, analysis

with GutIBMData("output.h5") as data:
    spatial = validation.validate_spatial_signatures(data, data.steps[-1])
    genomic = validation.validate_genomic_signatures(data)
    print(f"Monochromatic score: {spatial['monochromatic_score']:.3f}")
    print(f"Resident retention:  {genomic['resident_retention']:.1%}")
```

**CLI tools:**

```bash
# Golden regression (CI)
python -m gut_ibm_tools.validation_regression output.h5 \
  --golden python/tests/fixtures/eari_vadi_ci_golden.json --check-fish-targets

# Batch runner
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json
# or: gut-ibm-batch examples/batch_scan/batch.json
```

## Output Format

HDF5 file structure (compatible with `nufeb_tools`):

```
/step_000000/
  atoms/       id, type, x, y, z, radius, biomass, mu, state, lineage
  grid/        carbon, iron, b12, bacteriocin, …
  metadata/    time, step, num_agents, num_lineages
  lineage/     btuB_expression, fepA_expression, num_bi_loci, generation
  genome/      has_conjugative_plasmid, plasmid_cost_amelioration, …
```

Parallel HDF5 when built with MPI. Checkpoint groups support restart.

## Validation Targets

| Metric | Target | Source |
|---|---|---|
| Resident retention | 70–80% | Longitudinal metagenomics (multi-day runs) |
| Monochromatic patchiness | > 0.7 | HiPR-FISH / exclusion-radius clustering |
| Comet-tail asymmetry | > 1.5 | EARI/VADI advection models |
| Hopkins statistic | > 0.7 | Spatial clustering (VADI §75) |

Short CI runs use **golden-file regression** instead of full targets:
`bash scripts/validate_eari_vadi.sh`

## Project Structure

```
src/
  core/           Agent, Domain, SpatialHash, Simulation engine
  fields/         ChemicalField, AdvectionField, VBF
  diffusion/      Green's function kernels, QSSA solver, Barnes–Hut FMM
  fixes/          Biological rule modules (8 Fix plugins)
  genome/         Lineage tracker, plasmid library
  io/             HDF5 reader/writer, JSON input parser
  gpu/            Optional CUDA kernels
python/
  gut_ibm_tools/  HDF5 reader, analysis, validation, batch runner, FISH models
examples/         Scenario JSON + batch manifests
scripts/          CI helpers (validate, smoke, parity, scaling benchmark)
tests/            CTest unit/integration/MPI/GPU targets
docs/             Mechanisms, parameters, API, batch runner, scaling
```

## Documentation

| Document | Description |
|----------|-------------|
| [docs/CONFIG_FORMAT.md](docs/CONFIG_FORMAT.md) | Strict JSON input format, Fix selection, feature toggles |
| [docs/BATCH_RUNNER.md](docs/BATCH_RUNNER.md) | Resumable parameter-scan CLI and manifest schema |
| [docs/MECHANISMS.md](docs/MECHANISMS.md) | Biological mechanism deep-dives for each Fix module |
| [docs/PARAMETERS.md](docs/PARAMETERS.md) | All configurable parameters with defaults and units |
| [docs/API.md](docs/API.md) | C++ and Python class/function reference |
| [docs/SCALING.md](docs/SCALING.md) | Agent-count benchmarks and profiling |
| [EARI.md](EARI.md) | Eco-Advective Receptor Interference blueprint |
| [VADI.md](VADI.md) | Viscous Advective-Diffusion Interference blueprint |
| [AGENTS.md](AGENTS.md) | Developer/agent guidelines, CI map, known landmines |

Spec design documents: `GutIBM Spec 1_ Chemical Environment.md`, `GutIBM Spec 2_ Bacteriocin Induction.md`, `GutIBM Spec 3_ Cell Biology.md`.

## CI

GitHub Actions (`.github/workflows/ci.yml`) on every push/PR to `main`:

| Job | What it runs |
|-----|--------------|
| `unit-tests` | Fast CTest unit shard (Release, MPI, HDF5) |
| `integration-tests` | Smoke, config diversity, HDF5, batch runner smoke |
| `openmp-parity` | Serial vs OpenMP simulation fingerprints |
| `cuda-compile` | CUDA compile + GPU test targets (single arch) |
| `eari-vadi-validation` | Short sim → HDF5 → golden + FISH regression |
| `python-lint` | JSON syntax, ruff, pytest, batch dry-run |

Helper scripts: `scripts/validate_eari_vadi.sh`, `scripts/smoke_batch_runner.sh`, `scripts/run_scaling_benchmark.sh`.

## References

- NUFEB-2: Li et al., *NUFEB: A massively parallel simulator for individual-based modelling of microbial communities*
- EARI model: Eco-Advective Receptor Interference framework ([EARI.md](EARI.md))
- VADI model: Viscous Advective-Diffusion Interference framework ([VADI.md](VADI.md))

## License

See [LICENSE](LICENSE) for details.
