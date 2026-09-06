# GPU-vs-CPU outcome comparison: GSH vs GSD

Paired seeds: [62, 63]
Comparison step: 1345 (time=80700.0)


## WARNING: TRUNCATED RUNS PRESENT

**WARNING: supplied runs have different trajectory depths; all verdicts use common step 1345, not each run's own final step.**

Run last steps: `{'GSH:62': 1440, 'GSH:63': 1425, 'GSD:62': 1440, 'GSD:63': 1365, 'GSD:61': 1345}`

## Reproducibility control

Device repeats bit-identical: True
Dynamical reproducibility: True
Accounting reproducibility: True

### Repeat mismatch details

No repeat mismatches recorded.

## Distributional verdict

| estimand | kind | host median | device median | shift | seed spread | shift/spread | verdict |
|---|---|---|---|---|---|---|---|
| `bacteriostatic_live_agents` | stochastic | 75.5 | 198 | 122.5 | 310 | 0.395 | interchangeable |
| `cumulative_conjugation_transfers` | stochastic | 243.5 | 183 | 60.5 | 944 | 0.0641 | interchangeable |
| `cumulative_divisions` | stochastic | 4765 | 13008 | 8243 | 15499 | 0.532 | interchangeable |
| `cumulative_mortality_cdi` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_mortality_colicin` | stochastic | 10 | 24 | 14 | 34 | 0.412 | interchangeable |
| `cumulative_mortality_lysis` | stochastic | 54.5 | 193 | 138.5 | 376 | 0.368 | interchangeable |
| `cumulative_outflow_boundary` | stochastic | 1532.5 | 2715 | 1182.5 | 3558 | 0.332 | interchangeable |
| `cumulative_outflow_washout` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `cumulative_sos_inductions` | stochastic | 54.5 | 194 | 139.5 | 377 | 0.37 | interchangeable |
| `delivery_infeasible_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `delivery_reduction_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `delivery_retry_events_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `halt_reason_code` | invariant | 0 | 0 | 0 | 1 | 0 | identical |
| `maintenance_shortfall_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `max_toxin_BtuB` | stochastic | 3.05657e-05 | 0.000104342 | 7.37764e-05 | 0.000118705 | 0.622 | interchangeable |
| `max_toxin_CirA` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `max_toxin_FepA` | stochastic | 7.91006e-08 | 2.69425e-07 | 1.90324e-07 | 3.10101e-07 | 0.614 | interchangeable |
| `max_toxin_FhuA` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `mean_carbon` | stochastic | 0.00120973 | 0.0012097 | 2.69845e-08 | 4.70673e-08 | 0.573 | interchangeable |
| `mean_iron` | stochastic | 9.98278e-05 | 9.96096e-05 | 2.18197e-07 | 4.49814e-07 | 0.485 | interchangeable |
| `mean_oxygen` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `mean_realized_fermentation_fraction` | stochastic | 0 | 0 | 0 | 0 | — | identical |
| `n_total` | stochastic | 3794 | 10691 | 6897 | 11540 | 0.598 | interchangeable |
| `num_lineages` | stochastic | 3.5 | 5 | 1.5 | 6 | 0.25 | interchangeable |
| `reaction_clip_cumulative_total` | invariant | 1.98323e-12 | 5.30776e-12 | 3.32453e-12 | 6.22335e-12 | 0.534 | invariant_differs |
| `uptake_shortfall_cumulative_total` | invariant | 0 | 0 | 0 | 0 | — | identical |
| `washout_trapped_live_agents` | stochastic | 536.5 | 1388 | 851.5 | 2035 | 0.418 | interchangeable |

## Interpretation

Paired divergence is descriptive: a stochastic nonlinear model separates once a field difference flips one draw. The verdict is distributional; 'interchangeable' means the backend shift in the median is within the within-backend seed spread, which is reported alongside so the resolution of the measurement is visible.

## Flagged

These estimands moved by more than the within-backend seed spread, or differ on an invariant channel:

- `reaction_clip_cumulative_total`
