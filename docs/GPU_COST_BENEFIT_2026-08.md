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

## Two profile fields were misattributed, and both are host↔device copies

The first version of this document listed the following as unexplained and
hypothesized synchronization stalls. That hypothesis was wrong; both fields are
full-field marshalling counted under a phase name that does not describe it.
The raw values are retained rather than netted out.

- `spatial_hash_s` is 0.27 s host / 13.67 s device at `s1` and 1.46 / 10.17 at
  `s2`. The step-start device mirror upload — every species' `conc` *and*
  `reac` — happens inside that profiler lap, so on device the field times a
  whole-field H2D pass rather than the hash rebuild.
- `mpi_reaction_reduce_s` is ~3e-5 s host / 10.84 s device at `s1`, in
  single-rank runs where there is no peer to reduce with. On a single rank
  `try_sum_reactions_on_device` declines at `ranks <= 1`, so the device path
  pulled every reaction field back to host, copied it cell-by-cell into the
  field, called a rank sum that returns immediately for one rank, and the
  pipeline re-uploaded all of it.

Together that is ≈13% of the device's `s1` wall time, and with the counted
19.3 s of `gpu_h2d_s` + `gpu_d2h_s` it puts ≈43 s of the device run's 199 s
into host↔device marshalling — several times the 11 s chemistry regression the
campaign attributed to device execution. The copies these two fields measure
are removed in the follow-up work; the timings above are pre-removal.

- Every pass in this campaign reports `openmp_compiled: 0`. The benchmark image
  did not pass `-DGUTIBM_USE_OPENMP=ON` and the CMake option defaults off, so
  each host arm ran single-threaded on a 4-vCPU instance, including the QSSA
  superposition that has an OpenMP path. The host/device comparison is
  same-image and therefore internally valid, but every ratio here is against a
  serial CPU baseline, not against this machine's CPU capability.
- `neumann_drift_envelope_evaluations` reads 0 on every device record while the
  host records 1.24e8 sources outside the validated envelope for the same
  configuration. The counter is incremented in the host sealed concentration
  path only, so the diagnostic is blind on device — it is not evidence that
  device runs stayed inside the envelope. The pre-run `drift_envelope_policy`
  gate still applies to both backends.

## Remeasurement after #400 and #401: the chemistry regression was placement

The `A1`/`A4` `s1` pair was rerun on the same practice queue and instance type
with both mechanistic fixes in the image (`gutibm:gpubench-6bd6a522e3ad`,
digest `sha256:0dce54243b3503b5d1e4d9c4b3ece3b7f209d9f3116e9103eadf52f7a1855645`,
`git_sha 6bd6a522e3ad68ab9438a9527d32b47728451fb8`). Raw records:
`bench_results/gpu_cost_benefit_2026-08_postfix/`. Container time 3,484 s
≈ 0.97 instance-hours ≈ $0.51. Medians over seeds 55/56/57; both arms
`SUCCEEDED` with `openmp_compiled: 1` on every pass and
`chemistry_placement: host`/`device` as declared, so neither arm was silently
host-forced.

| Phase (s), `s1` | Host A1 pre-fix | Host A1 now | Device A4 pre-fix | Device A4 now |
|---|---:|---:|---:|---:|
| `chemistry_s` | 38.71 | 21.41 | 49.87 | **11.18** |
| `biology_s` | 602.99 | 445.16 | 92.49 | 94.46 |
| `spatial_hash_s` | 0.27 | 0.25 | 13.67 | 2.85 |
| `mpi_reaction_reduce_s` | ~3e-5 | ~3e-5 | 10.84 | 2.7e-5 |
| cost-pass wall | 674.9 | 369.4 | 178.4 | **129.0** |

Three corrections to the campaign above follow from this, and the first
retracts its headline finding:

1. **Device chemistry is not a regression; it is 1.91× faster than the host**
   (21.41 vs 11.18 s). The earlier 1.29× regression was host↔device
   marshalling charged to the chemistry phase, not device arithmetic: removing
   the dead whole-field uploads and the single-rank reaction round trip (#400)
   cut device chemistry 4.46× on unchanged numerics.
2. **The host baseline was serial.** With OpenMP compiled in, host chemistry
   improves 1.81× and host wall 1.83× on the same 4 vCPUs, so the campaign's
   3.78× end-to-end figure was T4-against-one-core. Against the CPU this
   machine actually has, the device wins **2.86×** end-to-end at `s1`
   (biology 4.71×, mechanics 8.8×).
3. **The two misattributed profile fields resolve as predicted.** The
   single-rank reaction round trip is gone entirely (10.84 s → 2.7e-5 s);
   `spatial_hash_s` falls 13.67 → 2.85 s, the residue being the concentration
   upload that still happens inside that lap after the concentration halo
   exchange.

### Transfers, now measured rather than inferred

With #401 the transfer profile reports bytes and call counts, and the timed
interval excludes queued kernel drain. Device `s1`, medians:

| Direction | Seconds | Bytes | Calls | Effective bandwidth |
|---|---:|---:|---:|---:|
| H2D | 6.07 | 28.37 GB | 1,950 | 4.67 GB/s |
| D2H | 6.03 | 28.97 GB | 678 | 4.80 GB/s |

Transfers are 12.1 s of the 133 s profiled run (9.1%), at ≈22 MB per call, so
per-call launch overhead is not the cost — bandwidth is, and 4.7-4.8 GB/s is
the expected pageable-memory rate on this link rather than an anomaly. Pinned
host staging would raise that toward the ≈10 GB/s pinned rate, saving ≈6 s of a
129 s run (≈5%), and it is not the largest remaining lever: the run still moves
57 GB in ten steps, ≈5.7 GB/step against a 100 MB species field, so eliminating
copies buys more than making the same copies faster. Neither is implemented,
and no bandwidth claim here is extrapolated — every number is from the
committed records.

## What this does and does not decide

It does not move any scientific default. The `delivery` uptake mode costs 1.27×
the host wall time of shipped `none` at `s1` (859.3 vs 674.9 s) and 2.28× its
host chemistry — but on device its chemistry is *cheaper* than shipped `none`'s
(44.65 vs 49.87 s), so device execution removes the chemistry-cost argument
against funded uptake at this scale. Whether funded uptake becomes a default
still rests on outcome movement, not on these timings. Those two chemistry
figures are pre-#400 and pre-OpenMP; the A2/A3/A5/A6 arms have not been rerun,
so the delivery-vs-`none` chemistry comparison is not current, while the
conclusion it supports — that device execution removes the chemistry-cost
argument — only strengthens now that device chemistry is faster in every mode
measured.

For anyone sizing a run: the T4 is worth its price for toxin-sampling-dominated
configurations (low agent count, many Green's-function sources) and for
mechanics-dominated ones, wins ≈1.9× on non-delivery chemistry after #400, and
returns 1.36× at a million agents against a serial host baseline — not the
order-of-magnitude the phrase "GPU port" tends to imply. The `s2` ratios and
the A2/A3/A5/A6 arms were measured before both fixes and against a
single-threaded CPU; only the `A1`/`A4` `s1` pair above is current.

## Arms not measured

The B (image-series precision), C (Robin transfer length) and D (drift
correction) axes are declared at `s0` and were not run; the report carries them
as `missing` rather than dropping them. D3 stays `blocked` because the
corrected sealed kernel routes to host by design (#394), so there is no device
arm to measure. No precision or scientific-outcome column is populated by this
campaign — `precision_summary` reads "not recorded" on every record, and any
GPU-vs-CPU numerical comparison remains unmeasured.
