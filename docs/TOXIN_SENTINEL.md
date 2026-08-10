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
| `colicin_kills` | 0 | 3 |
| `cdi_kills` | 0 | 0 |
| `washout_deaths` | 60 | 60 |
| `boundary_deaths` | 0 | 0 |
| `starvation_deaths` | 0 | 0 |
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
