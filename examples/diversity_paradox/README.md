# Diversity Paradox Full Simulation

Config format: strict JSON with optional `"_comment"` — see [CONFIG_FORMAT.md](../../docs/CONFIG_FORMAT.md).

Week-long simulation testing the Advective Double-Bind hypothesis.
Resident B2-phylogroup colonies (expressing Colicin E1 + Colicin B
as a "super-killer" array) compete against naive immigrant strains.

## Biological Setup

- **Residents** (Type 1, N=500): B2 phylogroup, carry ColE1+ColB plasmids,
  conjugative F-plasmid. Located near epithelium.
- **Immigrants** (Type 2, initial cohort N=100): No bacteriocin arsenal,
  plus continuous arrivals at the luminal boundary — one cell per hour on
  average, placed in the top 10 µm (`immigration` block). Must either acquire
  immunity via HGT or downregulate receptors to survive.

## Validation Targets

1. **Spatial**: Monochromatic patchiness score > 0.7 (HiPR-FISH analog)
2. **Genomic**: 70-80% resident lineage retention after 7 days
3. **Mechanism**: Immigrants that downregulate BtuB/FepA should show
   growth rates below the washout threshold and be flushed

## Running

```bash
# Full-scale (requires HPC)
mpirun -np 16 gut_ibm input.json

# Quick test (reduce domain and time in JSON, or use eari_vadi_validation)
gut_ibm input.json

# Large runs (10⁶+ agents): enable use_fmm in input.json — see docs/SCALING.md
```

## Outputs

Summary, agent, lineage, and genome data are written on their configured
schedules. Grid dumps of carbon, bacteriocin BtuB, and bacteriocin FepA are
written every 6 hours.

## Analysis

```python
from gut_ibm_tools import GutIBMData, validation

with GutIBMData("diversity_paradox_output.h5") as data:
    spatial = validation.validate_spatial_signatures(data, data.steps[-1])
    genomic = validation.validate_genomic_signatures(data)
    
    print(f"Monochromatic score: {spatial['monochromatic_score']:.3f}")
    print(f"Resident retention: {genomic['resident_retention']:.1%}")
```

## Measured behaviour at shipped parameters (12 h scaled run)

Enabling the documented immigration mechanism does **not** change this
scenario's regime. In a 400×400×100 µm scaled arm (x and y reduced 5×, all
biology, `grid_dx`, mucus thickness, turnover and transit unchanged, founders
scaled to the same areal density: 20 residents + 4 immigrants, seed 42,
720 steps of 60 s), immigration fired 14 times and every arriving lineage was
extinct within at most 4800 s of observed residence (median 600 s). Radial
advection carries a cell from z = 95 µm to the lumen boundary in a few hundred
seconds, against a resident doubling time of ln2/5.0e-4 ≈ 1390 s, so an
immigrant arriving in the luminal band is expelled before it can divide twice.
Cumulative boundary outflow was 61 with immigration against 33 without; live
residents at the horizon 121 against 195; live immigrants 0 in both arms.

Consequences for the validation targets above:

- The 70–80% resident-retention target is **not discriminating** at these
  parameters: type-based resident retention is 100% in both arms, because the
  initial immigrant cohort dies out with or without immigration.
- `monochromatic_score` was 1.000 in both arms and is degenerate once a single
  type remains, so it cannot evidence the spatial claim in this regime.
- Grid output now emits, and shows the bacteriocin fields are identically zero
  in most dumps: colicin release is tied to SOS lysis (3–4 events over 12 h)
  and the field decays within a few steps. Checking the comet-tail claim needs
  a dedicated high-cadence grid window around an induction, not the 6 h stride
  shipped here.

The full 7-day, 2×2 mm horizon has not been run: at that grid size a single
all-species grid dump is hundreds of MB per species. The targets above are
therefore unvalidated, not met — see `docs/EXTERNAL_AUDIT_2026-08.md` claim 10.
