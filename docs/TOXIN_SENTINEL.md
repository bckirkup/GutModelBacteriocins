# Toxin Sentinel Measurement

## Scope

This document records a mortality-attribution measurement for the **EARI/VADI
CI fixture**, not for the full-length EARI/VADI campaign.

The CI fixture runs for 900 s with `bio_dt = 60 s`, giving 15 biological
steps. It starts with 160 agents: 120 resident ColE1/ColB producers and
40 susceptible immigrant agents. The approximately 21-minute doubling time
is longer than this 15-minute run, and both measurements recorded zero
divisions.

The measurements were run with the GCC 13 rank-1 serial-compatible build
because Open MPI cannot initialize on this VM due to the known
`libhwloc` integer-divide-by-zero failure.

## Configurations

Both runs used the same EARI/VADI CI input and all default corrected-unit
parameters except for the explicitly stated override:

| Parameter | Corrected-unit run | Legacy-potency run |
|-----------|-------------------:|-------------------:|
| B12 field | `1e-3 mol/m³` | `1e-3 mol/m³` |
| `kd_b12_btuB` | `1e-6 mol/m³` | `1e-6 mol/m³` |
| `kd_colicinE_btuB` | `5e-7 mol/m³` | `5e-10 mol/m³` |
| BtuB competition factor | `1 + 1e-3/1e-6 = 1001` | `1001` |
| Apparent ColE1 Kd | `5.005e-4 mol/m³` | `5.005e-7 mol/m³` |

The legacy-potency run therefore differs by exactly one toxin-Kd override.
It restores the pre-units-migration apparent Kd while leaving the corrected
B12 field and corrinoid competition factor in place.

## Run-total event counters

Event counters were read from each HDF5 summary checkpoint and summed across
the 15-step run.

| Counter | Corrected units | Legacy toxin potency |
|---------|---------------:|---------------------:|
| `mortality_colicin` | 0 | 3 |
| `mortality_cdi` | 0 | 0 |
| `outflow_washout` | 60 | 60 |
| `outflow_boundary` | 0 | 0 |
| `mortality_lysis` | 0 | 0 |
| **Total counted deaths** | **60** | **63** |
| `divisions` | **0** | **0** |

### Death fractions

Percentages use the total counted death events in the corresponding run.

| Cause | Corrected units | Legacy toxin potency |
|-------|----------------:|---------------------:|
| Colicin | **0 / 60 = 0.00%** | **3 / 63 = 4.76%** |
| CDI | 0.00% | 0.00% |
| Washout | **60 / 60 = 100.00%** | **60 / 63 = 95.24%** |
| Boundary | 0.00% | 0.00% |
| Starvation | 0.00% | 0.00% |

## Interpretation

For this **CI fixture**, mortality is washout-dominated rather than
bacteriocin-dominated: corrected units produced zero colicin deaths, while
restoring the former toxin potency produced only 3 of 63 counted deaths.
Resident retention is therefore primarily reporting washout geometry in this
short assay.

This result says nothing measured about the full-length EARI/VADI campaign.
It characterises only the 900-second CI sentinel. In particular, the absence
of divisions is a property of this short run; it should not be extrapolated
to a longer campaign.

The microcolony arithmetic in [`UNITS_AUDIT.md`](UNITS_AUDIT.md) gives an
untested hypothesis for generalisation: toxin killing should be negligible
where co-located producer counts remain below roughly `10²`, with the relevant
threshold depending on distance and the simulated source dynamics. Checking
the campaign's realised local producer densities against that threshold is
the measurement needed to determine whether the CI-fixture conclusion
generalises.

## Testing consequence

The EARI/VADI CI fixture's population metrics are a weak sentinel for receptor
or toxin changes. They should not be cited as validation of toxin physics.
Re-anchoring those metrics after a toxin change is a consistency step, not
evidence that the changed toxin mechanism is correct.

This is the configuration-sensitivity gap described by the `ci-test-design`
skill: a golden regression can protect fixed outputs while remaining
insensitive to the parameter under test. A toxin-sensitive colony challenge
with controlled producer number and target distance is required alongside
the CI fixture.

## T2 colony-challenge assay

The companion CTest `toxin_sentinel` is the controlled T2 experiment called
for by the coherence diagnosis. It uses the corrected default
`kill_rate_colicin = 1e-3 s⁻¹`, a 5 µm grid, a 600 s run, and a compact
100-target cluster. ColE1 producers are placed at the same controlled
location and forced into `SOS_INDUCED` with the model's 300 s lysis delay.
Targets are placed either 10 µm or 50 µm away. Washout, motility, CDI,
crypt migration, and mechanics are disabled so that `mortality_colicin` is the
unambiguous outcome counter; no production defaults are changed.

The assay sweeps producer counts `{1, 10, 100, 1000, 10000}`. Fractions below
are `mortality_colicin / 100`; all runs had zero divisions.

| Distance | Producers | Killed fraction | Expected fraction* |
|---------:|----------:|----------------:|-------------------:|
| 10 µm | 1 | 0.00 | 0.001 |
| 10 µm | 10 | 0.00 | 0.01 |
| 10 µm | 100 | 0.10 | 0.08 |
| 10 µm | 1000 | 0.20 | 0.20 |
| 10 µm | 10000 | 0.20 | 0.30 |
| 50 µm | 1 | 0.00 | 0.0002 |
| 50 µm | 10 | 0.00 | 0.002 |
| 50 µm | 100 | 0.02 | 0.02 |
| 50 µm | 1000 | 0.20 | 0.10 |
| 50 µm | 10000 | 0.20 | 0.20 |

The exact deterministic run therefore brackets the switch-on between 100 and
1000 producers at 50 µm: 100 producers gives only 2% killing, while 1000
produces 20%. This is consistent with the audit's order-`10^3` estimate for
50 µm, rather than the earlier periodic-grid result that appeared to switch
on near 100 producers. The 10 µm result responds earlier, as expected.

The plateau at 0.20 for the upper producer counts is not occupancy saturation.
With unit occupancy, five 60-second hazard steps at the default
`1e-3 s⁻¹` kill rate give `H = 0.30` and
`1 - exp(-H) = 0.259`. Thus this 600-second assay has a dose-independent
measurement ceiling of about 26% deaths; the 1000- and 10000-producer rows
are hazard-limited and should not be interpreted as a saturated dose response.

### Analytical versus simulated reconciliation

For the 50 µm, 100-producer point, the toxin field sampled at the target grid
cell at the first post-lysis exposure step was
`5.60936e-5 mol/m³`. The corresponding single-source value is
`5.60936e-7 mol/m³`; the source-to-target agent distance is nominally 50 µm,
while the target grid-cell centre is offset by the 5 µm discretisation
(approximately 52.56 µm from the producer coordinate). The bounded
Green's-function calculation, using the actual source position, sampled-cell centre, configured
`D_eff = 8e-13 m²/s`, release age, decay, and z-boundary images, gives the
same single-source concentration (`~5.6e-7 mol/m³`). Thus the field itself,
not only the downstream hazard, is now reconciled with the transport kernel.
With the corrected BtuB apparent Kd,

```text
Kd_app = 5e-7 × (1 + 1e-3/1e-6) = 5.005e-4 mol/m³
occupancy(step 6) = 5.60936e-5 / (5.005e-4 + 5.60936e-5)
                  = 0.10078
```

Using the sampled field and its observed finite-burst decay, the assay
calculated five 60-second hazards beginning at that step:

```text
C(age) = C(step 6) × exp[-age/300 - (ln 2) age/1800]
H      = Σ (1e-3 s⁻¹ × occupancy(C(age)) × 60 s)
       = 0.02089
P      = 1 - exp(-H) = 0.02067
```

The simulation killed `2/100 = 0.02` targets at this point. The 0.13
percentage-point difference is within the assay's deterministic finite-target
resolution and is not a material disagreement. The test pins both the
field-to-Green's-function agreement and this analytical death probability.

An earlier version co-located targets while leaving the linear mechanics Fix
active. Repulsion moved targets between grid cells (the sampled cell changed
from the intended 10 µm target cell to a cell 40 µm away), producing an
apparent analytical mismatch. A second, independent discrepancy was then
traced to periodic-grid aliasing: with the initial 100 × 50 µm domain and a
60 µm cutoff, the near-field loops wrapped across the 20-cell x and
10-cell y dimensions and visited some cells multiple times. That inflated the
sampled field by several-fold, including the roughly 10.8-fold excess at the
50 µm point. The final assay disables x/y periodicity so each wrapped cell is
visited once; its sampled field now agrees with the bounded Green's function.
Both target-cell quantisation and the toxin cutoff were checked: the assay
uses a 60 µm cutoff for speed, while the production default is 200 µm; both
exceed the maximum requested 50 µm distance. The corrected field scales
linearly with producer count, confirming burst superposition.

The result supports the microcolony hypothesis without claiming that the
single-lysis steady-state estimates are a dynamic CTest threshold: finite
burst release, the 300 s lysis delay, grid sampling, and the compact-cluster
geometry all matter.
