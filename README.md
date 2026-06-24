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

### Computational Scaling

The simulation targets **10⁶–10⁷** discrete cells by:

1. **Avoiding explicit FTCS solvers** — at 1 μm resolution, the CFL stability
   criterion forces sub-millisecond timesteps. Instead, QSSA Green's function
   kernels compute steady-state fields analytically in O(1) per source.
2. **Spatial hashing** — O(N) neighbor lookups instead of O(N²) pairwise.
3. **VBF abstraction** — 99% of anaerobic microbiota are a continuous medium,
   not discrete agents.
4. **MPI domain decomposition** — scales across HPC nodes.

## Biological Mechanisms

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

### Implemented Fix Modules

| Fix | Purpose |
|---|---|
| `fix_metabolism` | Triple Monod kinetics (carbon, iron, B12), receptor-dependent Km, division, death |
| `fix_receptor` | Competitive binding at TBDTs (BtuB, FepA, CirA), toxin occupancy with ligand competition |
| `fix_bacteriocin` | Stochastic SOS lysis (1%/division), pI-dependent diffusion (Lethal Core vs. Halo) |
| `fix_conjugation` | Contact-dependent HGT via F-pili, shear-dependent MPS probability |
| `fix_mutation` | BI locus duplication/recombination, receptor downregulation, super-killer emergence |

### Toxin Biophysics

Bacteriocin diffusion is calibrated by isoelectric point (pI):

- **Lethal Cores** (pI > 8.5): Basic colicins (e.g., ColE1) bind negatively charged
  mucin glycoproteins → high retardation → concentrated near producer
- **Lethal Halos** (pI < 6.0): Acidic colicins (e.g., ColB) are repelled by mucin →
  wider, diffuse halos reaching further downstream
- **Comet Tails**: Mucus flow distorts radial diffusion into elongated downstream
  exclusion zones (validated by Péclet number ≥ 1)

### Plasmid Library

Pre-characterized *E. coli* bacteriocin systems:

| Name | pI | Class | Receptor | Type |
|---|---|---|---|---|
| Colicin E1 | 9.0 | Lethal Core | BtuB | Pore-forming |
| Colicin E2 | 6.5* | Lethal Halo | BtuB | Nuclease |
| Colicin B | 5.4 | Lethal Halo | FepA | Pore-forming |
| Colicin Ia | 7.2 | Neutral | CirA | Pore-forming |
| Colicin M | 9.3 | Lethal Core | FhuA | Murein inhibitor |
| Microcin V | 5.0 | Lethal Halo | CirA | Small peptide |

*\*ColE2 corrected: secreted as equimolar complex with acidic Im2 protein*

## Building

### Prerequisites

```bash
sudo apt-get install cmake libopenmpi-dev openmpi-bin libhdf5-mpi-dev
```

### Compile

```bash
mkdir build && cd build
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON
make -j$(nproc)
```

### Run Tests

```bash
cd build
ctest --output-on-failure
```

## Usage

```bash
# Default configuration (24h, 600 agents)
./gut_ibm

# Custom configuration
./gut_ibm ../examples/diversity_paradox/input.json

# MPI parallel
mpirun -np 8 ./gut_ibm ../examples/diversity_paradox/input.json
```

## Python Analysis Toolkit

```bash
cd python
pip install -e ".[viz]"
```

```python
from gut_ibm_tools import GutIBMData, validation, analysis

with GutIBMData("output.h5") as data:
    # Spatial validation (HiPR-FISH analog)
    spatial = validation.validate_spatial_signatures(data, data.steps[-1])
    
    # Genomic validation (longitudinal metagenomics)
    genomic = validation.validate_genomic_signatures(data)
    
    print(f"Monochromatic score: {spatial['monochromatic_score']:.3f}")
    print(f"Resident retention:  {genomic['resident_retention']:.1%}")
    print(f"Comet-tail ratio:    {spatial['comet_tail_ratio']:.2f}")
```

## Output Format

HDF5 file structure (compatible with `nufeb_tools`):

```
/step_000000/
  atoms/
    id, type, x, y, z, radius, biomass, mu, state, lineage
  grid/
    carbon, iron, b12, bacteriocin
  metadata/
    time, step, num_agents, num_lineages
  lineage/
    btuB_expression, fepA_expression, num_bi_loci, generation
```

## Validation Targets

| Metric | Target | Source |
|---|---|---|
| Resident retention | 70–80% | Human longitudinal metagenomics |
| Monochromatic patchiness | > 0.7 | HiPR-FISH mouse gut imaging |
| Comet-tail asymmetry | > 1.5 | Predicted by EARI/VADI models |
| Sweep-and-stasis pattern | Present | Allele time-series analysis |

## Project Structure

```
src/
  core/           Agent, Domain, SpatialHash, Simulation engine
  fields/         ChemicalField, AdvectionField, VBF
  diffusion/      Green's function kernels, QSSA solver
  fixes/          Biological rule modules (metabolism, receptor, etc.)
  genome/         Lineage tracker, plasmid library
  io/             HDF5 writer, input parser
python/
  gut_ibm_tools/  Analysis, validation, visualization
examples/
  single_colony/  Comet-tail formation test
  diversity_paradox/  Full Advective Double-Bind simulation
  eari_vadi_validation/  Short CI regression scenario (issue #56)
tests/            Unit tests for spatial hash, Green's functions, agents
```

## Documentation

| Document | Description |
|----------|-------------|
| [docs/MECHANISMS.md](docs/MECHANISMS.md) | Biological mechanism deep-dives for each Fix module |
| [docs/API.md](docs/API.md) | Class and function reference |
| [docs/PARAMETERS.md](docs/PARAMETERS.md) | All configurable parameters with defaults and units |
| [EARI.md](EARI.md) | Eco-Advective Receptor Interference blueprint |
| [VADI.md](VADI.md) | Viscous Advective-Diffusion Interference blueprint |

## CI

GitHub Actions runs on every push/PR to `main`:
- **Build matrix**: Release + Debug, OpenMP ON/OFF, Ubuntu latest, cmake + OpenMPI + HDF5
- **Tests**: CTest unit tests including MPI multi-rank and HDF5 round-trip
- **Python**: ruff lint + pytest on `gut_ibm_tools`
- **OpenMP parity**: serial vs OpenMP simulation fingerprint comparison
- **EARI/VADI validation**: short simulation → HDF5 → golden-file metric regression (`scripts/validate_eari_vadi.sh`)

## References

- NUFEB-2: Li et al., *NUFEB: A massively parallel simulator for individual-based
  modelling of microbial communities*
- EARI model: Eco-Advective Receptor Interference framework
- VADI model: Viscous Advective-Diffusion Interference framework

## License

See [LICENSE](LICENSE) for details.
