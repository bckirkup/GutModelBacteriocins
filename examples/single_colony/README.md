# Single Colony Comet-Tail Test

Tests the formation of advective "comet-tail" chemical exclusion zones
from a single resident Enterobacteriaceae colony expressing Colicin E1
(Lethal Core, pI ~9.0) and Colicin B (Lethal Halo, pI ~5.4).

## Expected Behavior

1. Resident colony maintains position near the epithelium
2. SOS-mediated lysis events (~1% per division) release toxin bursts
3. Mucus flow distorts radial diffusion into downstream comet tails
4. Colicin E1 (basic) forms a concentrated core near the producer
5. Colicin B (acidic) forms a wider, diffuse halo reaching further downstream

## Running

```bash
mpirun -np 4 gut_ibm input.json
```

## Analysis

```python
from gut_ibm_tools import GutIBMData, analysis

with GutIBMData("single_colony_output.h5") as data:
    step = data.steps[-1]
    grid = data.get_grid(step)
    comet = analysis.comet_tail_index(
        grid_positions, grid["bacteriocin"]
    )
    print(f"Comet tail ratio: {comet:.2f}")
```
