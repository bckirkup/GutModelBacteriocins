# EARI/VADI CI Validation Scenario

Short deterministic simulation used by the CI **eari-vadi-validation** job
(issue #56). Runs a resident–immigrant diversity-paradox setup for 1 simulated
hour and writes HDF5 output for Python regression checks.

## Biological focus

| Mechanism | Framework | What this run exercises |
|-----------|-----------|-------------------------|
| Advective comet-tails | EARI / VADI §35 | Mucus flow + SOS lysis → downstream toxin field |
| Spatial clustering | VADI §75 | Agent-type segregation metrics |
| Washout trap | VADI | Immigrant survival vs advective clearance |
| Lineage retention | EARI | Longitudinal resident fraction (short-run analog) |

Full validation targets (70–80% retention, monochromatic > 0.7, comet-tail > 1.5)
apply to multi-day runs such as `examples/diversity_paradox/`. This example
pins a **golden-file baseline** for regression instead.

## Run locally

```bash
cd build
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON
make -j$(nproc) gut_ibm
mpirun -np 1 ./gut_ibm ../examples/eari_vadi_validation/input.json

cd ../python
pip install -e ".[dev]"
python -m gut_ibm_tools.validation_regression \
  ../build/eari_vadi_validation.h5 \
  --golden tests/fixtures/eari_vadi_ci_golden.json \
  --fish-golden tests/fixtures/eari_vadi_ci_fish_golden.json \
  --check-fish-targets
```

Or use the CI helper script from the repo root:

```bash
bash scripts/validate_eari_vadi.sh
```

## Updating the golden baseline

If biology or numerics change intentionally, regenerate after a clean run:

```bash
python -m gut_ibm_tools.validation_regression \
  ../build/eari_vadi_validation.h5 \
  --write-golden tests/fixtures/eari_vadi_ci_golden.json \
  --scenario eari_vadi_ci

# FISH observability golden (issue #25):
python -m gut_ibm_tools.validation_regression \
  ../build/eari_vadi_validation.h5 \
  --write-fish-golden tests/fixtures/eari_vadi_ci_fish_golden.json \
  --fish-scenario eari_vadi_ci_fish
```
