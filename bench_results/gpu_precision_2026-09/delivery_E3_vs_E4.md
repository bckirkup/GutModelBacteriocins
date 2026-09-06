# GPU-vs-CPU outcome comparison: E3 vs E4

Paired seeds: [55, 56, 57, 58, 59]
Comparison step: 360 (time=21600.0)

## Reproducibility control

Device repeats bit-identical: True
Dynamical reproducibility: True
Accounting reproducibility: True

### Repeat mismatch details

No repeat mismatches recorded.

## Distributional verdict

| estimand | kind | host median | device median | shift | seed spread | shift/spread | verdict |
|---|---|---|---|---|---|---|---|
| `bacteriostatic_live_agents` | stochastic | 6 | 6 | 0 | 4 | 0 | interchangeable |
| `cumulative_conjugation_transfers` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_divisions` | stochastic | 100 | 100 | 0 | 37 | 0 | interchangeable |
| `cumulative_mortality_cdi` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_mortality_colicin` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_mortality_lysis` | stochastic | 5 | 5 | 0 | 6 | 0 | interchangeable |
| `cumulative_outflow_boundary` | stochastic | 522 | 520 | 2 | 16 | 0.125 | interchangeable |
| `cumulative_outflow_washout` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_sos_inductions` | stochastic | 5 | 5 | 0 | 5.5 | 0 | interchangeable |
| `delivery_infeasible_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `delivery_reduction_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `delivery_retry_events_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `halt_reason_code` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `maintenance_shortfall_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `max_toxin_BtuB` | stochastic | 0 | 0 | 0 | 7.17978e-06 | 0 | interchangeable |
| `max_toxin_CirA` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `max_toxin_FepA` | stochastic | 0 | 0 | 0 | 1.89977e-08 | 0 | interchangeable |
| `max_toxin_FhuA` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `mean_carbon` | stochastic | 0.00120969 | 0.00120969 | 1.25431e-12 | 1.81097e-08 | 6.93e-05 | interchangeable |
| `mean_iron` | stochastic | 9.9949e-05 | 9.99494e-05 | 3.78873e-10 | 1.36792e-09 | 0.277 | interchangeable |
| `mean_oxygen` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `mean_realized_fermentation_fraction` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `n_total` | stochastic | 75 | 75 | 0 | 29.5 | 0 | interchangeable |
| `num_lineages` | stochastic | 17 | 17 | 0 | 7 | 0 | interchangeable |
| `reaction_clip_cumulative_total` | invariant | 3.25843e-15 | 3.24936e-15 | 9.072e-18 | 4.07078e-16 | 0.0223 | identical |
| `uptake_shortfall_cumulative_total` | invariant | 3.88826e-16 | 3.88826e-16 | 1.3256e-23 | 1.51545e-16 | 8.75e-08 | identical |
| `washout_trapped_live_agents` | stochastic | 31 | 30 | 1 | 12.5 | 0.08 | interchangeable |

## Interpretation

Paired divergence is descriptive: a stochastic nonlinear model separates once a field difference flips one draw. The verdict is distributional; 'interchangeable' means the backend shift in the median is within the within-backend seed spread, which is reported alongside so the resolution of the measurement is visible.
