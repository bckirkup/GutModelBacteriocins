# Diversity Paradox Full Simulation

Config format: strict JSON with optional `"_comment"` — see [CONFIG_FORMAT.md](../../docs/CONFIG_FORMAT.md).

Week-long simulation testing the Advective Double-Bind hypothesis.
Resident B2-phylogroup colonies (expressing Colicin E1 + Colicin B
as a "super-killer" array) compete against naive immigrant strains.

## Biological Setup

- **Residents** (Type 1, N=500): B2 phylogroup, carry ColE1+ColB plasmids,
  conjugative F-plasmid. Located near epithelium.
- **Immigrants** (Type 2, N=100): No bacteriocin arsenal, introduced
  periodically via the lumen boundary. Must either acquire immunity
  via HGT or downregulate receptors to survive.

## Validation Targets

1. **Spatial**: Monochromatic patchiness score > 0.7 (HiPR-FISH analog)
2. **Genomic**: 70-80% resident lineage retention after 7 days
3. **Mechanism**: Immigrants that downregulate BtuB/FepA should show
   growth rates below the washout threshold and be flushed

## Running

```bash
# Full-scale (requires HPC)
mpirun -np 16 gut_ibm input.json

# Quick test (reduce domain and time)
gut_ibm input.json
```

## Analysis

```python
from gut_ibm_tools import GutIBMData, validation

with GutIBMData("diversity_paradox_output.h5") as data:
    spatial = validation.validate_spatial_signatures(data, data.steps[-1])
    genomic = validation.validate_genomic_signatures(data)
    
    print(f"Monochromatic score: {spatial['monochromatic_score']:.3f}")
    print(f"Resident retention: {genomic['resident_retention']:.1%}")
```
