# GPU cost/benefit measurement, 2026-08

Measured on the `gutibm-gpu-practice` AWS Batch queue, `g4dn.xlarge` On-Demand
(1× Tesla T4, driver 580.159.03, CUDA 13.0 runtime reported by `nvidia-smi`),
4 vCPU, one MPI rank, replicated chemical decomposition. Every record carries
`/run_provenance/chemistry_placement` and `gpu_compiled`, and no record outside
its arm's accepted placement set is included in any ratio below.

- Raw merged records: `bench_results/gpu_cost_benefit_2026-08/merged.json`
- Generated report: `bench_results/gpu_cost_benefit_2026-08/report.md`
- Image: `gutibm:gpubench-cfa62d65b5bb`,
  digest `sha256:3251eabcc80b2a87d61f79985cf089fca28c27021e271b88b4509c9058864c8e`
- Total measured container time 20,106 s ≈ 5.59 instance-hours ≈ $2.94 at
  $0.526/hour, inside the authorized $5-10.

The report is rendered from the merged JSON only; nothing in it is hand-edited.

## Scales

| Scale | Grid | Agents | Steps | Seeds |
|---|---|---:|---:|---|
| `s1` | 500 × 500 × 50 | 1e5 | 10 | 55, 56, 57 |
| `s2` | 500 × 500 × 50 | 1e6 | 5 | 55 |

`s2` is a cost scale at ten times the agent density of `s1`, not a
scientifically equivalent stand-in for the shipped million-agent grid: that
configuration needs ≈37.9 GB of host memory by the measured law
`RSS = 1.60e7 + 731.3·cells + 1352.8·agents` bytes, above the 32 GB of every
GPU instance available here, and outside Route B's ≤512-line device-delivery
scope. Both scales here sit inside that scope.

## What the GPU buys, by uptake mode (`s1`, 3 seeds, one instance)

Medians over three seeds. Chemistry is the phase the device chemistry port
touches; total wall is the whole run.

| Uptake mode | Host chemistry (s) | Device chemistry (s) | Host wall (s) | Device wall (s) |
|---|---:|---:|---:|---:|
| `none` (shipped) | 38.71 | 49.87 | 674.9 | 178.4 |
| `sherwood` | 38.60 | 49.42 | 666.9 | 176.0 |
| `delivery` | 88.07 | 44.65 | 859.3 | 237.5 |

Placements observed: `host` for A1/A2/A3, `device` for A4/A5,
`device_delivery` for A6 — Route B ran at `s1`, it was not host-forced.

Three results, and the first contradicts the premise the benchmark was designed
around:

1. **Chemistry on device is a regression for the two non-delivery modes** —
   1.29× slower for `none`, 1.28× for `sherwood`. H2D plus D2H transfer is
   19.3 s of the device's 49.9 s, so the copies alone exceed a third of the
   entire host chemistry phase. The device chemistry kernel does not pay for
   its own data movement at this grid size.
2. **For `delivery`, chemistry on device does pay: 1.97× faster** (88.07 →
   44.65 s). Delivery is the mode with substantially more host chemistry work
   to displace — it costs 2.28× the host chemistry of `none` — and Route B
   carries both sink terms on device rather than round-tripping per species.
   So the device delivery port, the narrowest and most recently landed piece of
   GPU work, is the only part of the chemistry story that is a win.
3. **Total wall improves 3.6-3.8× in every mode, and almost none of that is
   chemistry.** It is the QSSA Green's-function toxin sampling inside biology:
   602.99 → 92.49 s at seed 55, with kernel evaluation counts of the same order
   on both backends (1.53e10 host, 1.50e10 device), so this is throughput on
   identical work rather than work skipped — ≈39 ns/evaluation on 4 vCPUs
   against ≈6 ns on the T4.

## The speedup does not survive ten times the agent density (`s2`, 1 seed)

| Phase (s) | Host (A1) | Device (A4) | Ratio |
|---|---:|---:|---:|
| `biology_s` | 719.27 | 559.06 | 1.29× |
| `chemistry_s` | 20.92 | 26.33 | 0.79× |
| `physics_s` | 83.27 | 3.99 | 20.9× |
| `total_s` (profiled) | 833.74 | 626.86 | 1.33× |
| cost-pass wall | 843.93 | 621.57 | **1.36×** |

End-to-end benefit falls from 3.78× at 1e5 agents to 1.36× at 1e6 agents on the
same grid. The mechanism is visible in the profile rather than inferred:
Green's-function kernel evaluations drop to 2.10e9 (host) / 1.86e9 (device) —
7.3× fewer than `s1` despite ten times the agents — so the phase the device
actually accelerates is a much smaller share of biology at high density, while
the per-agent host work that dominates the rest is unchanged. Ordinary Amdahl
behaviour, with the accelerated fraction shrinking as density grows. *Why* the
evaluation count falls that steeply is not established here; the octree/FMM
approximation path becoming active at higher source density is the obvious
candidate and remains a hypothesis.

Mechanics is the second real device win and it was not part of the stated
rationale: `physics_s` 83.27 → 3.99 s, which only becomes material at `s2`
where mechanics is 10% of the host run rather than 0.5%.

Cost per completed replicate at $0.526/hour: `s1` $0.099 host vs $0.026 device;
`s2` $0.123 host vs $0.091 device.

## Unexplained, and deliberately not smoothed over

- `spatial_hash_s` is 0.27 s host / 13.67 s device at `s1` and 1.46 / 10.17 at
  `s2`, and `mpi_reaction_reduce_s` is ~3e-5 s host / 10.84 s device at `s1`,
  in single-rank runs where there is no peer to reduce with. Together that is
  ≈13% of the device's `s1` wall time. Neither phase runs on the device, so the
  likely explanation is that these host-side timers are absorbing
  synchronization stalls for device work queued earlier in the step; that is a
  measurement-attribution hypothesis, not a verified cause, and the raw values
  are retained rather than netted out.
- `neumann_drift_envelope_evaluations` reads 0 on every device record while the
  host records 1.24e8 sources outside the validated envelope for the same
  configuration. The counter is incremented in the host sealed concentration
  path only, so the diagnostic is blind on device — it is not evidence that
  device runs stayed inside the envelope. The pre-run `drift_envelope_policy`
  gate still applies to both backends.

## What this does and does not decide

It does not move any scientific default. The `delivery` uptake mode costs 1.27×
the host wall time of shipped `none` at `s1` (859.3 vs 674.9 s) and 2.28× its
host chemistry — but on device its chemistry is *cheaper* than shipped `none`'s
(44.65 vs 49.87 s), so device execution removes the chemistry-cost argument
against funded uptake at this scale. Whether funded uptake becomes a default
still rests on outcome movement, not on these timings.

For anyone sizing a run: the T4 is worth its price for toxin-sampling-dominated
configurations (low agent count, many Green's-function sources) and for
mechanics-dominated ones, roughly breaks even on non-delivery chemistry, and
returns 1.36× at a million agents — not the order-of-magnitude the phrase "GPU
port" tends to imply.

## Arms not measured

The B (image-series precision), C (Robin transfer length) and D (drift
correction) axes are declared at `s0` and were not run; the report carries them
as `missing` rather than dropping them. D3 stays `blocked` because the
corrected sealed kernel routes to host by design (#394), so there is no device
arm to measure. No precision or scientific-outcome column is populated by this
campaign — `precision_summary` reads "not recorded" on every record, and any
GPU-vs-CPU numerical comparison remains unmeasured.
