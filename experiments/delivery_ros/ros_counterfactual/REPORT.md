# ROS counterfactual campaign report

Valid arms are the two Group A arms and four corrected Group B arms. Full per-step trajectories, SOS component ledgers, and all nonzero clip intervals are in `metrics.json`; `metrics.csv` contains the per-arm summary rows.

## Build and policy

- Commit: `abf393750453673ebb1f4c4ab54f447797a647fd`
- Binary SHA-256: `859fc77a8e5280b2eb0271ae0603cbf4965adf891ff5edc2ba4fd2e3817221d6`
- Stop policy: relative clip threshold `1e-6` against cumulative agent growth uptake; funded-minus-realized checked at every ledger record.
- `delivery_reduction_*`: absent from every emitted HDF5 output.

## Run summary

| Arm | dx (µm) | k_ROS | founders | final N | divisions | termination | time (s) | wall (s) |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| A_ros0_res2 | 2 | 0 | 4 | 5 | 5 | horizon_reached | 21600 | 63.37 |
| A_ctrl_res2 | 2 | 100 | 4 | 1 | 0 | population_stop | 3540 | 12.09 |
| B_ros0_res2 | 2 | 0 | 80 | 76 | 103 | horizon_reached | 21600 | 83.77 |
| B_ctrl_res2 | 2 | 100 | 80 | 21 | 44 | horizon_reached | 21600 | 209.18 |
| B_ros0_res6 | 6 | 0 | 80 | 171 | 263 | horizon_reached | 21600 | 2.83 |
| B_ctrl_res6 | 6 | 100 | 80 | 68 | 99 | horizon_reached | 21600 | 7.04 |

## Population loss closure

| Arm | lysis | colicin | CDI | washout-trapped | boundary outflow | channel sum | observed N decrement | divisions − losses | closure |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| A_ros0_res2 | 0 | 0 | 0 | 0 | 4 | 4 | -1 | 4 | True |
| A_ctrl_res2 | 3 | 0 | 0 | 0 | 0 | 3 | 3 | 3 | True |
| B_ros0_res2 | 5 | 0 | 0 | 0 | 102 | 107 | 4 | 107 | True |
| B_ctrl_res2 | 75 | 1 | 0 | 0 | 27 | 103 | 59 | 103 | True |
| B_ros0_res6 | 6 | 0 | 0 | 0 | 166 | 172 | -91 | 172 | True |
| B_ctrl_res6 | 72 | 2 | 0 | 0 | 37 | 111 | 12 | 111 | True |

Loss channels include mortality and transport. The channel sum closes against `initial N + divisions + immigrations − final N`; therefore boundary/washout entries are transported out, while lysis/colicin/CDI entries are deaths.

## Oxygen and carbon ledgers

| Arm | O₂ growth funded/demanded | O₂ growth fraction | O₂ maintenance funded/demanded | O₂ maint shortfall | O₂ agent/flora removal | agent share | C growth funded/demanded | C maintenance funded/demanded | C maint shortfall | C clips (relative) |
|---|---|---:|---|---:|---|---:|---|---|---:|---|
| A_ros0_res2 | 4.891759e-14 / 6.6242768e-14 | 0.73845932 | 1.0246176e-13 / 1.0818e-13 | 5.7182442e-15 | 1.5137935e-13 / 0 | 1 | 3.2906054e-15 / 4.660169e-15 | 1.8598495e-15 / 1.8687605e-15 | 8.9109916e-18 | 0 (0 final; max 0) |
| A_ctrl_res2 | 9.8255665e-15 / 1.1427706e-14 | 0.85980219 | 1.0585099e-14 / 1.098e-14 | 3.9490064e-16 | 2.0410666e-14 / 0 | 1 | 4.7899141e-16 / 4.7899141e-16 | 1.5797625e-16 / 1.5797625e-16 | 0 | 0 (0 final; max 0) |
| B_ros0_res2 | 6.696129e-13 / 1.0799575e-12 | 0.62003632 | 2.0314399e-12 / 2.47878e-12 | 4.4734013e-13 | 2.7010528e-12 / 3.5807843e-12 | 0.42997816 | 7.5149352e-14 / 1.0947009e-13 | 3.9499169e-14 / 3.9613475e-14 | 1.1430617e-16 | 1.1507041e-20 (1.5312229e-07 final; max 1.640424e-07) |
| B_ctrl_res2 | 4.355552e-13 / 6.2831981e-13 | 0.69320622 | 9.3876742e-13 / 1.03314e-12 | 9.4372577e-14 | 1.3743226e-12 / 3.8972422e-12 | 0.26070487 | 3.6640982e-14 / 4.7232205e-14 | 1.6691953e-14 / 1.6746181e-14 | 5.4227894e-17 | 0 (0 final; max 0) |
| B_ros0_res6 | 2.590145e-12 / 2.7326248e-12 | 0.94785973 | 3.5429211e-12 / 3.63468e-12 | 9.175893e-14 | 6.133066e-12 / 1.3629827e-12 | 0.81817318 | 1.2775621e-13 / 1.2811721e-13 | 6.1295834e-14 / 6.1339543e-14 | 4.3708774e-17 | 0 (0 final; max 0) |
| B_ctrl_res6 | 1.0912699e-12 / 1.1043346e-12 | 0.98816962 | 1.1873775e-12 / 1.19334e-12 | 5.9625492e-15 | 2.2786474e-12 / 1.8048421e-12 | 0.55801475 | 4.4916233e-14 / 4.5045415e-14 | 1.9928115e-14 / 1.9952085e-14 | 2.3969636e-17 | 0 (0 final; max 0) |

All funded-minus-realized final residuals are zero at reported precision. The largest positive cumulative floating residual seen during the stepwise check was `4.04e-28 mol` (round-off).

## Clips

| Arm | O₂ cumulative clip | O₂ relative | C cumulative clip | C relative | nonzero interval records |
|---|---:|---:|---:|---:|---:|
| A_ros0_res2 | 0 | 0 | 0 | 0 | 0 |
| A_ctrl_res2 | 0 | 0 | 0 | 0 | 0 |
| B_ros0_res2 | 0 | 0 | 1.1507041e-20 | 1.640424e-07 | 241 |
| B_ctrl_res2 | 0 | 0 | 0 | 0 | 0 |
| B_ros0_res6 | 0 | 0 | 0 | 0 | 0 |
| B_ctrl_res6 | 0 | 0 | 0 | 0 | 0 |

For `B_ros0_res2`, the first nonzero carbon interval is recorded in `metrics.json` at step 120 / 7200 s; maximum interval-relative carbon clip is below `1e-6`, so the run continued. Other arms have no nonzero clips.

## SOS hazard components

| Arm | cumulative basal | post-division | nuclease cross-induction | ROS | total | SOS inductions | final interval components (basal, post, nuclease, ROS) | final per-live-agent interval components |
|---|---:|---:|---:|---:|---:|---:|---|---|
| A_ros0_res2 | 0.001803 | 0.0016666667 | 0 | 0 | 0.0034696667 | 0 | (4.9999999999999996e-06, 0.0, 0.0, 0.0) | (1e-06, 0.0, 0.0, 0.0) |
| A_ctrl_res2 | 0.000171 | 0 | 0 | 0.049943297 | 0.050114297 | 3 | (1e-06, 0.0, 0.0, 0.0004888456129685165) | (1e-06, 0.0, 0.0, 0.0004888456129685165) |
| B_ros0_res2 | 0.030334 | 0.026 | 0 | 0 | 0.056334 | 5 | (5.099999999999995e-05, 0.0, 0.0, 0.0) | (6.710526315789467e-07, 0.0, 0.0, 0.0) |
| B_ctrl_res2 | 0.006327 | 0.0063333333 | 0 | 1.3075139 | 1.3201742 | 74 | (0.0, 0.0, 0.0, 0.0) | (0.0, 0.0, 0.0, 0.0) |
| B_ros0_res6 | 0.044664 | 0.066666667 | 0 | 0 | 0.11133067 | 6 | (0.00012399999999999976, 0.0, 0.0, 0.0) | (7.251461988304079e-07, 0.0, 0.0, 0.0) |
| B_ctrl_res6 | 0.003509 | 0.004 | 0 | 1.1117678 | 1.1192768 | 71 | (0.0, 0.0, 0.0, 0.0) | (0.0, 0.0, 0.0, 0.0) |

Each SOS component sum equals the emitted total (`cumulative_component_sum_matches_total = true` for all arms). Full interval and cumulative records include per-live-agent means in JSON.

## Fermentation, acetate, grids

| Arm | mean fermentation | final fermentation | final acetate mean | acid inhibition | final O₂ mean/min | final C mean/min |
|---|---:|---:|---:|---|---|---|
| A_ros0_res2 | 0.17132187 | 0.31899901 | 7.2733538e-05 | True | 0.049475856209215 / 0.00023225590366592551 | 0.002885852358427332 / 0.00033207867238222995 |
| A_ctrl_res2 | 0.051038816 | 0.024616678 | 7.2737968e-05 | True | 0.04960218668939039 / 0.018441696894782988 | 0.0028862802581340243 / 0.0011634530688178222 |
| B_ros0_res2 | 0.34609819 | 0.58930594 | 7.5780135e-05 | True | 0.035058103088567306 / 0.0037828851652289305 | 0.002826983716323394 / 0.00013236787469423126 |
| B_ctrl_res2 | 0.26616127 | 0.40640831 | 7.3272006e-05 | True | 0.03806984032972317 / 0.00015311606516340315 | 0.0028686188095472295 / 4.210024795162863e-05 |
| B_ros0_res6 | 0.052357127 | 0.095782808 | 5.9396493e-05 | True | 0.021586165628110374 / 0.000736874757927254 | 0.0026662708847324736 / 0.000722607838676039 |
| B_ctrl_res6 | 0.009247819 | 0.016298485 | 5.9396493e-05 | True | 0.03178281018733855 / 0.004189750075658367 | 0.002813662743109695 / 0.0008220225932864078 |

Final-grid agent-cell concentrations and cell/domain-mean ratios are retained per species in `metrics.json` under `final_grid`.

## Sampled trajectories

Samples use exact saved agent/summary times where available. `—` means the run ended before that sample.

### A_ros0_res2

| time (s) | N | divisions interval | mean mu_realized (s⁻¹) | mean biomass (mol) | bacteriostatic | washout-trapped |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 4 | 0 | 0.00055 | 5.7595865e-16 | 0 | 0 |
| 600 | 4 | 0 | 0.00011295875 | 6.2076518e-16 | 0 | 2 |
| 1200 | 4 | 0 | 0.00011163054 | 6.6216141e-16 | 0 | 2 |
| 3600 | 4 | 0 | 8.6237632e-05 | 8.472627e-16 | 0 | 2 |
| 7200 | 6 | 2 | 5.3602485e-05 | 7.7782458e-16 | 0 | 6 |
| 10800 | 6 | 4 | 5.2275519e-05 | 7.5640211e-16 | 0 | 6 |
| 14400 | 5 | 5 | 4.6035011e-05 | 8.0418115e-16 | 0 | 5 |
| 18000 | 5 | 5 | 2.9285114e-06 | 8.8300916e-16 | 3 | 5 |
| 21600 | 5 | 5 | 1.0998992e-05 | 9.3662128e-16 | 0 | 5 |

### A_ctrl_res2

| time (s) | N | divisions interval | mean mu_realized (s⁻¹) | mean biomass (mol) | bacteriostatic | washout-trapped |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 4 | 0 | 0.00055 | 5.7595865e-16 | 0 | 0 |
| 600 | 4 | 0 | 0.00011295875 | 6.2076518e-16 | 0 | 2 |
| 1200 | 4 | 0 | 0.00011163054 | 6.6216141e-16 | 0 | 2 |

### B_ros0_res2

| time (s) | N | divisions interval | mean mu_realized (s⁻¹) | mean biomass (mol) | bacteriostatic | washout-trapped |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 80 | 0 | 0.0005375 | 5.7595865e-16 | 0 | 0 |
| 600 | 80 | 0 | 0.00011208383 | 6.2089673e-16 | 0 | 28 |
| 1200 | 80 | 0 | 0.00010738353 | 6.6186123e-16 | 0 | 33 |
| 3600 | 80 | 0 | 9.7198396e-05 | 8.3734155e-16 | 0 | 51 |
| 7200 | 94 | 14 | 7.6634017e-05 | 9.6355858e-16 | 9 | 85 |
| 10800 | 156 | 79 | 3.9936166e-05 | 6.9049569e-16 | 0 | 156 |
| 14400 | 146 | 79 | 2.46498e-05 | 7.9253913e-16 | 17 | 146 |
| 18000 | 116 | 98 | 1.7331948e-05 | 7.5999555e-16 | 3 | 116 |
| 21600 | 76 | 103 | 8.8047134e-06 | 7.6490078e-16 | 8 | 76 |

### B_ctrl_res2

| time (s) | N | divisions interval | mean mu_realized (s⁻¹) | mean biomass (mol) | bacteriostatic | washout-trapped |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 80 | 0 | 0.0005375 | 5.7595865e-16 | 0 | 0 |
| 600 | 79 | 0 | 0.00011249928 | 6.208924e-16 | 0 | 27 |
| 1200 | 73 | 0 | 0.00010810603 | 6.6214559e-16 | 0 | 28 |
| 3600 | 46 | 0 | 0.00010399209 | 8.4263996e-16 | 0 | 23 |
| 7200 | 61 | 27 | 9.0763704e-05 | 6.7218215e-16 | 3 | 52 |
| 10800 | 49 | 36 | 5.530571e-05 | 7.3017886e-16 | 0 | 49 |
| 14400 | 48 | 39 | 2.96016e-05 | 8.0823583e-16 | 5 | 48 |
| 18000 | 37 | 44 | 1.8305172e-05 | 7.7895883e-16 | 2 | 37 |
| 21600 | 21 | 44 | 1.4810451e-05 | 8.8770556e-16 | 5 | 21 |

### B_ros0_res6

| time (s) | N | divisions interval | mean mu_realized (s⁻¹) | mean biomass (mol) | bacteriostatic | washout-trapped |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 80 | 0 | 0.0005375 | 5.7595865e-16 | 0 | 0 |
| 600 | 80 | 0 | 0.00013808361 | 6.2904685e-16 | 0 | 22 |
| 1200 | 80 | 0 | 0.00013600127 | 6.8280075e-16 | 0 | 26 |
| 3600 | 80 | 0 | 0.00012524513 | 9.335752e-16 | 0 | 42 |
| 7200 | 157 | 77 | 0.00010549335 | 7.20411e-16 | 0 | 137 |
| 10800 | 187 | 110 | 8.6369188e-05 | 8.4186486e-16 | 0 | 187 |
| 14400 | 253 | 186 | 6.2842171e-05 | 7.7273013e-16 | 0 | 253 |
| 18000 | 204 | 214 | 5.2479847e-05 | 8.520946e-16 | 0 | 204 |
| 21600 | 171 | 263 | 4.9361785e-05 | 7.7352937e-16 | 0 | 171 |

### B_ctrl_res6

| time (s) | N | divisions interval | mean mu_realized (s⁻¹) | mean biomass (mol) | bacteriostatic | washout-trapped |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 80 | 0 | 0.0005375 | 5.7595865e-16 | 0 | 0 |
| 600 | 75 | 0 | 0.00013814803 | 6.2904181e-16 | 0 | 20 |
| 1200 | 67 | 0 | 0.00013743237 | 6.8355961e-16 | 0 | 20 |
| 3600 | 37 | 0 | 0.00012748192 | 9.3512006e-16 | 0 | 20 |
| 7200 | 48 | 29 | 0.00011161094 | 7.1271373e-16 | 0 | 45 |
| 10800 | 42 | 33 | 9.1282102e-05 | 9.5828875e-16 | 1 | 42 |
| 14400 | 71 | 67 | 7.5128512e-05 | 7.0903575e-16 | 0 | 71 |
| 18000 | 56 | 71 | 6.2764583e-05 | 9.0073932e-16 | 0 | 56 |
| 21600 | 68 | 99 | 5.6766137e-05 | 6.7034258e-16 | 0 | 68 |

## Invalid superseded arm

- `B_ros0_res2` 100-founder output is excluded from `arms` metrics.
- Preserved at `invalid_B_ros0_res2_100_founders/` and recorded in `metrics.json` under `invalid_arms`.
- Reason: first dispatch used 80 strain-1 founders plus 20 strain-2 founders instead of 60 + 20 = 80 total.
