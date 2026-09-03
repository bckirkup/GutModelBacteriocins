# GPU cost/benefit measurement: A2/A3/A5/A6 s1

- Source/image git SHA: `310e97b04f386ff434c724218643d5602558b41f`
- Image tag: `gpubench-310e97b04f38`
- Image digest: `sha256:f5a03461190a111ea1e6594913c146c484bfc8c140375f360163811a6a5b4ce4`
- AWS rate used for estimate: `$0.526/instance-hour`
- Seeds: `55`, `56`, `57`

## A2/s1

- Job ID: `b2ffb5d0-87d7-49a7-9b10-38a341a4532d`; Batch status: **SUCCEEDED**; exit code: `0`; status reason: `Essential container in task exited`
- Container duration: `2653.457 s`; estimated container cost: `$0.388`
- Expected chemistry placement: `host`; reportable as expected placement: **True**
- Failure reasons: `[]`; expected/found blocks: `3/3` for each seed

- `openmp_compiled`: `1`; `mpi_rank_count`: `1`; GPU: `Tesla T4` driver `580.159.03`; recorded git_sha: `310e97b04f386ff434c724218643d5602558b41f`

### Phase/timing and scalar transfer metrics

| seed | cost wall (s) | chemistry (s) | biology (s) | physics (s) | spatial hash (s) | ghost exchange (s) | H2D bytes | D2H bytes | H2D calls | D2H calls | H2D s | D2H s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 55 | 370.964214 | 21.853676 | 452.391230 | 2.907400 | 0.244143 | 10.669697 | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| 56 | 411.423984 | 20.673585 | 467.452157 | 2.817311 | 0.238504 | 10.678035 | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| 57 | 363.490018 | 21.597635 | 452.192528 | 2.824691 | 0.233697 | 10.678399 | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| **median** | **370.964214** | **21.597635** | **452.391230** | **2.824691** | **0.238504** | **10.678035** | **0** | **0** | **0** | **0** | **0.000000** | **0.000000** |

### Complete per-call-site transfer table

Each site has `h2d_bytes`, `d2h_bytes`, `h2d_calls`, `d2h_calls`, `h2d_s`, `d2h_s`; zero rows below mean the label was absent or had zero traffic in the saved profile.

#### seed 55

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| step_agents | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_reactions | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_diffusion_result | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_fmm | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_near_field | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| mechanics | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| unattributed | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |

#### seed 56

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| step_agents | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_reactions | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_diffusion_result | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_fmm | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_near_field | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| mechanics | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| unattributed | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |

#### seed 57

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| step_agents | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_reactions | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_diffusion_result | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_fmm | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_near_field | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| mechanics | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| unattributed | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |

#### median

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| step_agents | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_reactions | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_diffusion_result | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_fmm | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_near_field | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| mechanics | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| unattributed | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |

- Site/scalar invariant differences (all four fields) by seed: `[{'seed': 55, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}, {'seed': 56, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}, {'seed': 57, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}]`

## A3/s1

- Job ID: `a7a82649-07b4-4419-b050-0ea5c5d44381`; Batch status: **SUCCEEDED**; exit code: `0`; status reason: `Essential container in task exited`
- Container duration: `3445.102 s`; estimated container cost: `$0.503`
- Expected chemistry placement: `host_forced_delivery`; reportable as expected placement: **False**
- Failure reasons: `[]`; expected/found blocks: `3/3` for each seed

- **Placement gate mismatch:** recorded pass placements were seed 55: cost=host, profile=host, seed 56: cost=host, profile=host, seed 57: cost=host, profile=host. The values are retained but are not reportable under the expected placement.

- `openmp_compiled`: `1`; `mpi_rank_count`: `1`; GPU: `Tesla T4` driver `580.159.03`; recorded git_sha: `310e97b04f386ff434c724218643d5602558b41f`

### Phase/timing and scalar transfer metrics

| seed | cost wall (s) | chemistry (s) | biology (s) | physics (s) | spatial hash (s) | ghost exchange (s) | H2D bytes | D2H bytes | H2D calls | D2H calls | H2D s | D2H s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 55 | 428.719470 | 75.375626 | 406.254611 | 10.756128 | 0.244059 | 10.686532 | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| 56 | 523.567777 | 75.844857 | 545.994500 | 10.796944 | 0.259987 | 10.681945 | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| 57 | 581.270308 | 75.593973 | 636.500428 | 11.273477 | 0.266482 | 10.683934 | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| **median** | **523.567777** | **75.593973** | **545.994500** | **10.796944** | **0.259987** | **10.683934** | **0** | **0** | **0** | **0** | **0.000000** | **0.000000** |

### Complete per-call-site transfer table

Each site has `h2d_bytes`, `d2h_bytes`, `h2d_calls`, `d2h_calls`, `h2d_s`, `d2h_s`; zero rows below mean the label was absent or had zero traffic in the saved profile.

#### seed 55

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| step_agents | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_reactions | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_diffusion_result | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_fmm | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_near_field | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| mechanics | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| unattributed | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |

#### seed 56

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| step_agents | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_reactions | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_diffusion_result | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_fmm | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_near_field | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| mechanics | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| unattributed | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |

#### seed 57

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| step_agents | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_reactions | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_diffusion_result | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_fmm | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_near_field | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| mechanics | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| unattributed | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |

#### median

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| step_agents | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_reactions | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_diffusion_result | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_fmm | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_near_field | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| mechanics | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| unattributed | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |

- Site/scalar invariant differences (all four fields) by seed: `[{'seed': 55, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}, {'seed': 56, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}, {'seed': 57, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}]`

## A5/s1

- Job ID: `19514d8d-a698-4e74-bf00-93711d107fa2`; Batch status: **SUCCEEDED**; exit code: `0`; status reason: `Essential container in task exited`
- Container duration: `775.282 s`; estimated container cost: `$0.113`
- Expected chemistry placement: `device`; reportable as expected placement: **True**
- Failure reasons: `[]`; expected/found blocks: `3/3` for each seed

- `openmp_compiled`: `1`; `mpi_rank_count`: `1`; GPU: `Tesla T4` driver `580.159.03`; recorded git_sha: `310e97b04f386ff434c724218643d5602558b41f`

### Phase/timing and scalar transfer metrics

| seed | cost wall (s) | chemistry (s) | biology (s) | physics (s) | spatial hash (s) | ghost exchange (s) | H2D bytes | D2H bytes | H2D calls | D2H calls | H2D s | D2H s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 55 | 128.502745 | 11.148280 | 92.730560 | 0.341665 | 2.796374 | 10.685200 | 31,470,969,816 | 29,366,712,872 | 1,994 | 692 | 6.671676 | 6.052549 |
| 56 | 128.509873 | 11.160057 | 93.895021 | 0.360729 | 2.831464 | 10.675912 | 29,969,728,056 | 28,166,569,720 | 1,940 | 650 | 6.381078 | 5.841900 |
| 57 | 129.509308 | 11.153103 | 94.463365 | 0.328503 | 2.834161 | 10.679637 | 30,967,776,632 | 28,966,346,232 | 1,976 | 678 | 6.571494 | 5.965328 |
| **median** | **128.509873** | **11.153103** | **93.895021** | **0.341665** | **2.831464** | **10.679637** | **30,967,776,632** | **28,966,346,232** | **1,976** | **678** | **6.571494** | **5.965328** |

### Complete per-call-site transfer table

Each site has `h2d_bytes`, `d2h_bytes`, `h2d_calls`, `d2h_calls`, `h2d_s`, `d2h_s`; zero rows below mean the label was absent or had zero traffic in the saved profile.

#### seed 55

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 11,000,000,000 | 0 | 110 | 0 | 2.282015 | 0.000000 |
| step_agents | 229,943,040 | 37,365,692 | 220 | 70 | 0.065581 | 0.012726 |
| chem_reactions | 11,000,000,000 | 11,000,000,000 | 110 | 110 | 2.292073 | 2.245064 |
| chem_diffusion_result | 0 | 11,000,000,000 | 0 | 110 | 0.000000 | 2.281113 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 4,000,000,000 | 0 | 40 | 0.000000 | 0.820135 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 2,400,000,000 | 0 | 24 | 0 | 0.499328 | 0.000000 |
| qssa_fmm | 5,600,005,408 | 3,200,000,000 | 96 | 32 | 1.172903 | 0.653239 |
| greens_near_field | 11,904 | 640 | 64 | 80 | 0.000513 | 0.000821 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 68,982,816 | 5,748,568 | 20 | 10 | 0.016628 | 0.001976 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 21,748,648 | 5,748,568 | 60 | 20 | 0.008006 | 0.002108 |
| mechanics | 0 | 80 | 0 | 20 | 0.000000 | 0.000269 |
| unattributed | 1,150,278,000 | 117,849,324 | 1,290 | 200 | 0.334630 | 0.035097 |

#### seed 56

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 11,000,000,000 | 0 | 110 | 0 | 2.287900 | 0.000000 |
| step_agents | 229,745,920 | 37,333,660 | 220 | 70 | 0.066908 | 0.012859 |
| chem_reactions | 11,000,000,000 | 11,000,000,000 | 110 | 110 | 2.298530 | 2.263356 |
| chem_diffusion_result | 0 | 11,000,000,000 | 0 | 110 | 0.000000 | 2.285887 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 4,000,000,000 | 0 | 40 | 0.000000 | 0.824795 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 3,000,000,000 | 0 | 30 | 0 | 0.631316 | 0.000000 |
| qssa_fmm | 3,500,005,824 | 2,000,000,000 | 60 | 20 | 0.733280 | 0.414868 |
| greens_near_field | 16,512 | 400 | 40 | 50 | 0.000329 | 0.000563 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 68,923,680 | 5,743,640 | 20 | 10 | 0.016812 | 0.001998 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 21,743,720 | 5,743,640 | 60 | 20 | 0.008055 | 0.002123 |
| mechanics | 0 | 80 | 0 | 20 | 0.000000 | 0.000274 |
| unattributed | 1,149,292,400 | 117,748,300 | 1,290 | 200 | 0.337948 | 0.035176 |

#### seed 57

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 11,000,000,000 | 0 | 110 | 0 | 2.293340 | 0.000000 |
| step_agents | 229,437,440 | 37,283,532 | 220 | 70 | 0.067709 | 0.012662 |
| chem_reactions | 11,000,000,000 | 11,000,000,000 | 110 | 110 | 2.287127 | 2.242890 |
| chem_diffusion_result | 0 | 11,000,000,000 | 0 | 110 | 0.000000 | 2.276079 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 4,000,000,000 | 0 | 40 | 0.000000 | 0.820828 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 2,600,000,000 | 0 | 26 | 0 | 0.541413 | 0.000000 |
| qssa_fmm | 4,900,007,072 | 2,800,000,000 | 84 | 28 | 1.022407 | 0.572853 |
| greens_near_field | 14,976 | 560 | 56 | 70 | 0.000443 | 0.000793 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 68,831,136 | 5,735,928 | 20 | 10 | 0.016580 | 0.002000 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 21,736,008 | 5,735,928 | 60 | 20 | 0.008003 | 0.002086 |
| mechanics | 0 | 80 | 0 | 20 | 0.000000 | 0.000269 |
| unattributed | 1,147,750,000 | 117,590,204 | 1,290 | 200 | 0.334471 | 0.034867 |

#### median

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 11,000,000,000 | 0 | 110 | 0 | 2.287900 | 0.000000 |
| step_agents | 229,745,920 | 37,333,660 | 220 | 70 | 0.066908 | 0.012726 |
| chem_reactions | 11,000,000,000 | 11,000,000,000 | 110 | 110 | 2.292073 | 2.245064 |
| chem_diffusion_result | 0 | 11,000,000,000 | 0 | 110 | 0.000000 | 2.281113 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 4,000,000,000 | 0 | 40 | 0.000000 | 0.820828 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 2,600,000,000 | 0 | 26 | 0 | 0.541413 | 0.000000 |
| qssa_fmm | 4,900,007,072 | 2,800,000,000 | 84 | 28 | 1.022407 | 0.572853 |
| greens_near_field | 14,976 | 560 | 56 | 70 | 0.000443 | 0.000793 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 68,923,680 | 5,743,640 | 20 | 10 | 0.016628 | 0.001998 |
| delivery | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| spatial_hash | 21,743,720 | 5,743,640 | 60 | 20 | 0.008006 | 0.002108 |
| mechanics | 0 | 80 | 0 | 20 | 0.000000 | 0.000269 |
| unattributed | 1,149,292,400 | 117,748,300 | 1,290 | 200 | 0.334630 | 0.035097 |

- Site/scalar invariant differences (all four fields) by seed: `[{'seed': 55, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}, {'seed': 56, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}, {'seed': 57, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}]`

## A6/s1

- Job ID: `7d0b7db8-eba0-4962-bc29-cbd96ce7a6e6`; Batch status: **SUCCEEDED**; exit code: `0`; status reason: `Essential container in task exited`
- Container duration: `1290.199 s`; estimated container cost: `$0.189`
- Expected chemistry placement: `device_delivery`; reportable as expected placement: **True**
- Failure reasons: `[]`; expected/found blocks: `3/3` for each seed

- `openmp_compiled`: `1`; `mpi_rank_count`: `1`; GPU: `Tesla T4` driver `580.159.03`; recorded git_sha: `310e97b04f386ff434c724218643d5602558b41f`

### Phase/timing and scalar transfer metrics

| seed | cost wall (s) | chemistry (s) | biology (s) | physics (s) | spatial hash (s) | ghost exchange (s) | H2D bytes | D2H bytes | H2D calls | D2H calls | H2D s | D2H s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 55 | 213.882788 | 25.706945 | 154.229047 | 8.234957 | 2.873914 | 10.684417 | 35,551,506,696 | 19,284,079,452 | 1,738 | 554 | 7.523885 | 3.960502 |
| 56 | 213.833239 | 25.731305 | 158.488901 | 8.259986 | 2.868648 | 10.682613 | 35,550,236,376 | 19,283,992,452 | 1,738 | 554 | 7.517342 | 3.987625 |
| 57 | 216.201254 | 25.767554 | 158.394743 | 8.091667 | 2.842206 | 10.685353 | 35,546,809,200 | 19,283,758,944 | 1,738 | 554 | 7.537364 | 4.006021 |
| **median** | **213.882788** | **25.731305** | **158.394743** | **8.234957** | **2.868648** | **10.684417** | **35,550,236,376** | **19,283,992,452** | **1,738** | **554** | **7.523885** | **3.987625** |

### Complete per-call-site transfer table

Each site has `h2d_bytes`, `d2h_bytes`, `h2d_calls`, `d2h_calls`, `h2d_s`, `d2h_s`; zero rows below mean the label was absent or had zero traffic in the saved profile.

#### seed 55

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 11,000,000,000 | 0 | 110 | 0 | 2.309804 | 0.000000 |
| step_agents | 231,927,040 | 37,688,092 | 220 | 70 | 0.067052 | 0.013134 |
| chem_reactions | 11,000,000,000 | 0 | 110 | 0 | 2.286833 | 0.000000 |
| chem_diffusion_result | 0 | 11,000,000,000 | 0 | 110 | 0.000000 | 2.275313 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 4,000,000,000 | 0 | 40 | 0.000000 | 0.797388 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 3,000,000,000 | 0 | 30 | 0 | 0.626504 | 0.000000 |
| qssa_fmm | 3,500,003,744 | 2,000,000,000 | 60 | 20 | 0.729692 | 0.408354 |
| greens_near_field | 9,408 | 400 | 40 | 50 | 0.000328 | 0.000607 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 69,578,016 | 5,798,168 | 20 | 10 | 0.016956 | 0.002035 |
| delivery | 5,800,000,000 | 2,200,002,576 | 58 | 114 | 1.207674 | 0.448976 |
| spatial_hash | 21,798,248 | 5,798,168 | 60 | 20 | 0.008060 | 0.002122 |
| mechanics | 0 | 80 | 0 | 20 | 0.000000 | 0.000299 |
| unattributed | 928,190,240 | 34,791,968 | 1,030 | 100 | 0.270982 | 0.012275 |

#### seed 56

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 11,000,000,000 | 0 | 110 | 0 | 2.293392 | 0.000000 |
| step_agents | 231,687,040 | 37,649,092 | 220 | 70 | 0.067405 | 0.012932 |
| chem_reactions | 11,000,000,000 | 0 | 110 | 0 | 2.297631 | 0.000000 |
| chem_diffusion_result | 0 | 11,000,000,000 | 0 | 110 | 0.000000 | 2.280489 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 4,000,000,000 | 0 | 40 | 0.000000 | 0.812510 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 3,000,000,000 | 0 | 30 | 0 | 0.625491 | 0.000000 |
| qssa_fmm | 3,500,006,240 | 2,000,000,000 | 60 | 20 | 0.729013 | 0.409415 |
| greens_near_field | 14,592 | 400 | 40 | 50 | 0.000322 | 0.000591 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 69,506,016 | 5,792,168 | 20 | 10 | 0.016922 | 0.002057 |
| delivery | 5,800,000,000 | 2,200,002,576 | 58 | 114 | 1.208100 | 0.454947 |
| spatial_hash | 21,792,248 | 5,792,168 | 60 | 20 | 0.008070 | 0.002118 |
| mechanics | 0 | 80 | 0 | 20 | 0.000000 | 0.000291 |
| unattributed | 927,230,240 | 34,755,968 | 1,030 | 100 | 0.270995 | 0.012276 |

#### seed 57

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 11,000,000,000 | 0 | 110 | 0 | 2.291659 | 0.000000 |
| step_agents | 231,042,880 | 37,544,416 | 220 | 70 | 0.066636 | 0.012908 |
| chem_reactions | 11,000,000,000 | 0 | 110 | 0 | 2.311752 | 0.000000 |
| chem_diffusion_result | 0 | 11,000,000,000 | 0 | 110 | 0.000000 | 2.285665 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 4,000,000,000 | 0 | 40 | 0.000000 | 0.815791 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 3,000,000,000 | 0 | 30 | 0 | 0.627682 | 0.000000 |
| qssa_fmm | 3,500,007,488 | 2,000,000,000 | 60 | 20 | 0.734150 | 0.412982 |
| greens_near_field | 16,320 | 400 | 40 | 50 | 0.000332 | 0.000579 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 69,312,768 | 5,776,064 | 20 | 10 | 0.016815 | 0.002040 |
| delivery | 5,800,000,000 | 2,200,002,576 | 58 | 114 | 1.210135 | 0.461260 |
| spatial_hash | 21,776,144 | 5,776,064 | 60 | 20 | 0.008106 | 0.002129 |
| mechanics | 0 | 80 | 0 | 20 | 0.000000 | 0.000288 |
| unattributed | 924,653,600 | 34,659,344 | 1,030 | 100 | 0.270096 | 0.012381 |

#### median

| label | h2d_bytes | d2h_bytes | h2d_calls | d2h_calls | h2d_s | d2h_s |
|---|---:|---:|---:|---:|---:|---:|
| step_conc_upload | 11,000,000,000 | 0 | 110 | 0 | 2.293392 | 0.000000 |
| step_agents | 231,687,040 | 37,649,092 | 220 | 70 | 0.067052 | 0.012932 |
| chem_reactions | 11,000,000,000 | 0 | 110 | 0 | 2.297631 | 0.000000 |
| chem_diffusion_result | 0 | 11,000,000,000 | 0 | 110 | 0.000000 | 2.280489 |
| chem_host_fallback | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_host_sync | 0 | 4,000,000,000 | 0 | 40 | 0.000000 | 0.812510 |
| qssa_toxin_robin_mix | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| qssa_toxin_zero | 3,000,000,000 | 0 | 30 | 0 | 0.626504 | 0.000000 |
| qssa_fmm | 3,500,006,240 | 2,000,000,000 | 60 | 20 | 0.729692 | 0.409415 |
| greens_near_field | 14,592 | 400 | 40 | 50 | 0.000328 | 0.000591 |
| greens_grid_download | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 |
| receptor | 69,506,016 | 5,792,168 | 20 | 10 | 0.016922 | 0.002040 |
| delivery | 5,800,000,000 | 2,200,002,576 | 58 | 114 | 1.208100 | 0.454947 |
| spatial_hash | 21,792,248 | 5,792,168 | 60 | 20 | 0.008070 | 0.002122 |
| mechanics | 0 | 80 | 0 | 20 | 0.000000 | 0.000291 |
| unattributed | 927,230,240 | 34,755,968 | 1,030 | 100 | 0.270982 | 0.012276 |

- Site/scalar invariant differences (all four fields) by seed: `[{'seed': 55, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}, {'seed': 56, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}, {'seed': 57, 'h2d_bytes': 0, 'd2h_bytes': 0, 'h2d_calls': 0, 'd2h_calls': 0}]`
