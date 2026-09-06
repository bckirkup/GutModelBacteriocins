# What GPU acceleration buys GutIBM — campaign summary (2026-09-04..06)

Queue `gutibm-gpu-practice`, g4dn.xlarge (T4), image `gpubench-7f6f84521f25`.
Total spend across both campaigns: ~$20 of the $30 authorization.

## 1. Precision campaign (p1: 500-agent collapse, 360 steps, seeds 55-59)

| Comparison | Paired seeds | Verdict | Wall host | Wall device | Speedup |
|---|---|---|---|---|---|
| E1 host vs E2 device (uptake=none) | 55,57,58,59 | 27/27 estimands interchangeable | 1716-2422 s | 996-1044 s | 1.7-2.4x |
| E3 host vs E4 device_delivery (uptake=delivery) | 55-59 | 27/27 interchangeable, device repeats bit-identical | 2735-4232 s | 1303-1345 s | 2.1-3.2x |

Placement provenance: `device` (E2), `device_delivery` (E4) on every run. Only
repeat noise: 1 ULP in `reaction_clip_cumulative_total` (accounting, not dynamics).

## 2. Science campaign (GS: diversity_paradox, 1 mm x 1 mm x 0.1 mm, 1440 steps, seeds 61-63)

Residents (type 1, 500) + immigrants (type 2, 100, continuous 1/h). The
scenario is bimodal: seed 62 collapses (resident-only remnant), seeds 61/63
bloom to ~11-12k agents and trip the dysbiosis guard near the end.

| Seed | Regime | Host outcome | Device outcome | Wall host | Wall device | Speedup |
|---|---|---|---|---|---|---|
| 61 | bloom | timed out at 8 h (no data) | 11,582 agents, guard halt step 1345 | >28,800 s | 5,684 s | >5.1x |
| 62 | collapse | 137 agents, full 1440 steps | 51 agents, full 1440 steps | 8,388 s | 4,039 s | 2.1x |
| 63 | bloom | 11,651 agents, guard halt step 1425 | 12,123 agents, guard halt step 1365 | 19,632 s | 5,074 s | 3.9x |

Comparator (`gpu_precision`, host 62/63 vs device 61/62/63, compared at step
1345): both paired seeds bit-identical until step ~100-180 (mean_carbon onset
101 / 184), then stochastic divergence with the same regime on both backends;
all dynamical estimands `interchangeable`; the only flagged estimand is the
accounting-only `reaction_clip_cumulative_total`. Placement provenance
`host`/`device` as intended on all five runs.

## 3. Reading

- Scientific: no evidence of a backend bias in either the grid-dominated or the
  agent-dense regime; device runs are drop-in replicates.
- Operational: speedup scales with agent load — ~2x when the chemistry grid
  dominates (p1, collapse seed), ~4x at ~10k agents, >5x on the densest seed.
  At the 1440-step horizon the host arm needs 2.3-8+ h per replicate against
  1.1-1.6 h on the T4, i.e. the GPU is what makes production-horizon replicate
  sweeps feasible on one instance-day.
- Not answered: the 10x-dense arm (GDD/GDH, guard off) exceeded 8 h even on the
  device; a dense benchmark needs the guard on or a shorter horizon.

Artifacts: `bench_results/gpu_precision_2026-09/` (comparator JSON + markdown, campaign design,
job script, job IDs). Raw CloudWatch logs and per-run precision exports are
kept out of the repo (3.4 MB); job IDs in `gpux_submitted_jobs.txt` locate them.
