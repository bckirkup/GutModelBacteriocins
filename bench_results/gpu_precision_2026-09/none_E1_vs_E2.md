# GPU-vs-CPU outcome comparison: E1 vs E2

Paired seeds: [55, 57, 58, 59]
Comparison step: 360 (time=21600.0)

## Reproducibility control

Device repeats bit-identical: False
Dynamical reproducibility: True
Accounting reproducibility: False

### Repeat mismatch details

| estimand | class | mismatch steps | max ULP | max relative difference |
|---|---|---:|---:|---:|
| `reaction_clip_cumulative_total` | accounting | 69 | 1 | 1.26218e-29 |
| `reaction_clip_cumulative_total` | accounting | 110 | 1 | 1.26218e-29 |

### Accounting repeat-noise caveat

Accounting repeat noise was observed in: reaction_clip_cumulative_total (max ULP 1, max relative difference 1.26218e-29); reaction_clip_cumulative_total (max ULP 1, max relative difference 1.26218e-29). Verdicts stand; these accounting counters are judged against their measured repeat noise, and a host/device difference below that floor is not attributable to the backend.

## Distributional verdict

| estimand | kind | host median | device median | shift | seed spread | shift/spread | verdict |
|---|---|---|---|---|---|---|---|
| `bacteriostatic_live_agents` | stochastic | 3.5 | 3 | 0.5 | 5 | 0.1 | interchangeable |
| `cumulative_conjugation_transfers` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_divisions` | stochastic | 109 | 102 | 7 | 38.5 | 0.182 | interchangeable |
| `cumulative_mortality_cdi` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_mortality_colicin` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_mortality_lysis` | stochastic | 4 | 5 | 1 | 6 | 0.167 | interchangeable |
| `cumulative_outflow_boundary` | stochastic | 530 | 526 | 4 | 11 | 0.364 | interchangeable |
| `cumulative_outflow_washout` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_sos_inductions` | stochastic | 4 | 5 | 1 | 5 | 0.2 | interchangeable |
| `delivery_infeasible_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `delivery_reduction_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `delivery_retry_events_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `halt_reason_code` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `maintenance_shortfall_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `max_toxin_BtuB` | stochastic | 0 | 0 | 0 | 3.1253e-05 | 0 | interchangeable |
| `max_toxin_CirA` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `max_toxin_FepA` | stochastic | 0 | 0 | 0 | 8.0761e-08 | 0 | interchangeable |
| `max_toxin_FhuA` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `mean_carbon` | stochastic | 0.00120974 | 0.00120974 | 1.12059e-11 | 1.60705e-10 | 0.0697 | interchangeable |
| `mean_iron` | stochastic | 9.99496e-05 | 9.99494e-05 | 2.24116e-10 | 1.44594e-09 | 0.155 | interchangeable |
| `mean_oxygen` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `mean_realized_fermentation_fraction` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `n_total` | stochastic | 76 | 77 | 1 | 42 | 0.0238 | interchangeable |
| `num_lineages` | stochastic | 16 | 17 | 1 | 4.5 | 0.222 | interchangeable |
| `reaction_clip_cumulative_total` | invariant | 7.50771e-14 | 7.40104e-14 | 1.06673e-15 | 1.24463e-14 | 0.0857 | identical |
| `uptake_shortfall_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `washout_trapped_live_agents` | stochastic | 28.5 | 29 | 0.5 | 4.5 | 0.111 | interchangeable |

## Interpretation

Accounting repeat noise was observed in: reaction_clip_cumulative_total (max ULP 1, max relative difference 1.26218e-29); reaction_clip_cumulative_total (max ULP 1, max relative difference 1.26218e-29). Verdicts stand; these accounting counters are judged against their measured repeat noise, and a host/device difference below that floor is not attributable to the backend.
