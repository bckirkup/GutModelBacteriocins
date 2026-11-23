# GutIBM Project Board

Living kanban for open work on [GutModelBacteriocins](https://github.com/bckirkup/GutModelBacteriocins).

**Last updated:** 2026-11-22

## One-click setup (local)

With `gh` authenticated and `project` scope:

```bash
./scripts/setup_project_board.sh
```

This creates GitHub labels, milestones, a Projects v2 board, and links open issues. See [scripts/setup_project_board.sh](../scripts/setup_project_board.sh).

**Manual UI:** Repository → **Projects** → **New project** → **Board** → name it `GutIBM Roadmap`, then drag issues from the tables below.

---

## Board columns

| Column | Meaning |
|--------|---------|
| **Done** | Merged to `main`; close the issue |
| **Backlog** | Scoped but not started; no active branch |
| **Ready** | PR bundle defined; pick up next |
| **In Progress** | Branch open / agent or human working |
| **In Review** | PR open, awaiting merge |

---

## Current board state (Nov 2026 — Post–Cycle 2 Hardening)

The Sep–Nov 2026 campaigns, Layer 2/3 multiscale architecture waves (#307–#473), and all Cycle 2 implementation PRs (**#484–#495**) are **merged**, closing the entire `#464–#483` Cycle 2 issue queue. Two issues remain open on GitHub:

### Open Issues (2 open as of 2026-11-22)

| Issue | Track | Priority | Scope |
| :--- | :--- | :--- | :--- |
| #446 | Campaign Ops | `maintenance` | Campaign storage index (automated S3 collector) |
| #159 | GPU FMM Far-Field | `backlog` | Sub-quadratic CPU M2L delivered in #164; full device-side octree traversal remains open |

### Done (Nov 2026 Cycle 2 wave — PRs #484–#495)

| Issue(s) | Track | PR | Scope |
| :--- | :--- | :--- | :--- |
| #483 | Docs / Tracker Hygiene | #486 | Reconciled Spec 12/13 prose, project board, and legacy P6 issues (#155–#161) |
| #482 | Campaign Governance | #485 | Stage D decision-record axis enforcement in `preflight.py` and HDF5 `container_image_digest` auth in `analyze.py` |
| #471, #474, #481 | Observation & Analysis | #484 | Optical-slice `depth_window_m`, abundance-weighted retention, pellet/cage LOD censoring, and canonical division-wave slope/tail-median estimators |
| #465, #476 | Spec 13 (`GAP-02` + I/O) | #487 | Multi-slot live patches, active `try_deposit`, split `sos_lysis`/`toxin_kill` ledger, and `segment_io.cpp` fixes |
| #464, #468, #477 | Spec 14 & Bacteriocin | #488 | Per-division pair SOS normalization, `D_free_colicin` wiring, BtuB/B12 decoupling (`metabolism.b12_growth_essential`), and carriage-matched controls |
| #467, #470 | Spec 13 & Oxygen BCs | #489 | Post-48h occupancy readout, 2×3 `reattach_p0` sensitivity, and Dirichlet/Robin O₂ BC presets |
| #466, #469, #474 | Spec 13 Layer 2 | #491 | Configurable `live_physics_mode`, delivery alignment & wet-mass docs (`1,100 kg/m³`), and perturbation pulse schedules |
| #470, #472, #479 | Metabolism, GPU & HDF5 | #490 | Implicit iron VBF sink, CUDA delivery/funded metabolism, bitwise checkpoint RNG replay (`hdf5_rng_replay`), and QSSA transient bounds |
| Wave 2 hardening | Checkpoint & Layer 2/3 | #492 | Seed-fork RNG/tag hardening, live O₂/stride inheritance, and Layer-3 pulse/death checkpoint restore |
| #478 | Numerics / QSSA | #495 | Hybrid Sturm-Liouville / Hankel wall-normal Peclet sealed kernel on CPU + CUDA GPU and `qssa.drift_envelope_policy=auto` |
| #480 | HPC / GPU & MPI Memory | #493 | Pre-allocation `ChemicalField` multi-rank OOM guard (`auto_slab`/`error`/`warn`), lazy delivery vectors, device-resident GPU stepping, and `PinnedHostBuffer<T>` staging |
| #475 | Spec 14 (`GAP-10`) | #494 | Default-off `FixProphage` temperate prophage carriage, DNA-damage Hill induction, virion PDE/adsorption ledger, CUDA parity, and human/murine presets |

### Historical Jul/Aug 2026 summary

The Jul/Aug 2026 implementation wave (#256–#306, with #303 superseded before
merge) is **closed**. It completed the MPI ownership and ledger work, the
physical-T4 GPU gate and the device defects it exposed, provenance and halt
reason output, finite-rate epithelial delivery, carbon accounting, the RPS
configuration surface, pI-driven transport, and the Pirt/metabolic-mode
maintenance mechanisms.

### In Review

_None_

### Ready (recommended next bundles)

| Issue | Title | Track | Priority |
|-------|-------|-------|----------|
| — | Delivery re-bracket after the maintenance-budget fix; see [density-limitation README §6](../experiments/density_limitation/README.md#6-next-experiment) | `track:campaign` | high |

### Backlog

| Issue | Title | Track | Priority |
|-------|-------|-------|----------|
| — | AWS Batch Spot + CUDA (`us-east-1`): smoke on `g4dn.xlarge` then Stage 3 ([docs/AWS_BATCH.md](AWS_BATCH.md)) | `track:infra` | high (practice CUDA smoke first) |
| — | MPI/HPC validation `mpirun -np 8+` | `track:mpi` | low (manual HPC) |
| — | GPU FMM octree traversal on device | `track:gpu` | low |

### Done (Jul 2026 wave)

| Issue | Title | PR |
|-------|-------|-----|
| #152 | GPU ROI (profiling, QSSA, Fur/O₂, FMM far-field, receptor CSR) | #152 |
| #155 | Python integration pytest in CI | #162 |
| #158 | GPU CPU/GPU parity CI | #162 |
| #154 | MPI four-rank validation | #163 |
| #156 | CUDA-aware MPI reaction reduce | #163 |
| #157 | GPU mechanics force kernel | #164 |
| #159 | Sub-quadratic FMM M2L | #164 |
| #160 | Metabolic washout trap long-horizon regression | #167 |
| #161 | OpenMP stochastic toxin-kill parity | #168 |

### Done (Aug 2026 wave)

| Issue | Title | PR |
|-------|-------|-----|
| #256 | Remove nightly benchmark | #256 |
| #257 | Boundary-flux accounting | #257 |
| #258 | Implicit carbon sink and reaction-clip accounting | #258 |
| #259 | Remove trajectory goldens | #259 |
| #260–#268 | Chemistry decomposition, cross-rank toxin sources, slab invariants/storage, shared MPI paths, and per-axis chemistry stride | #260–#268 |
| #269–#273 | GPU slab chemistry, agent-sampled QSSA, toxin lumping/species subset, and grid-halo configuration | #269–#273 |
| #274–#282 | Nested batch overrides, MPI event counters, process-group smoke, counted lysis, founder placement, population stocks, overdamped mechanics, and population-ledger semantics | #274–#282 |
| #283 | Reachable bacteriostasis threshold | #283 |
| #284–#287 | Runtime entrypoint packaging, run provenance, AWS quota/log reads, and noise-tolerant dysbiosis guard | #284–#287 |
| #288–#290 | GPU toxin-burst parity, halt provenance, and host reaction-loss accounting | #288–#290 |
| #291–#294 | Direct GPU kernel tests, device diagnostics, GPU Fur metabolism, and the physical T4 device gate | #291–#294 |
| #295–#298 | Finite-rate epithelial boundaries, Sonar cleanup, transport-limited uptake, and ECR image retention | #295–#298 |
| #299 | Nuclease-only SOS induction and immunity | #299 |
| #300 | RPS strain and plasmid configuration surface | #300 |
| #301 | pI-derived mucin retardation | #301 |
| #302 | Pirt carbon maintenance sink | #302 |
| #304–#305 | Oxygen-dependent metabolic modes, acid inhibition, CUDA signatures, and negative-growth host/device parity | #304–#305 |
| #306 | Maintenance routing through uptake limits and corrected carbon-ledger semantics | #306 |

### Done (Jun 2026 wave — reference)

| Issue | Title | Notes |
|-------|-------|-------|
| #40–#43 | Washout, MPI transfer, plasmids, multi-rank tests | |
| #44, #52, #59, #80, #81 | HDF5 checkpoint + parallel I/O | |
| #56 | EARI/VADI CI validation | |
| #55 | Scaling benchmark smoke + driver | |
| #25 | HCR-FISH / DNA-FISH models | |
| #29 | Higher-order FMM (CPU) | |
| #33 | GPU acceleration (phase 1) | Extended by #152 |
| #75–#79, #76 | MPI/parser/config/docs hygiene | |

---

## Merge order (remaining)

_No open P6 merge queue._

---

## What not to bundle

| Don't combine | Reason |
|---------------|--------|
| #160 + #161 | Different subsystems (biology regression vs OpenMP RNG) | Resolved (#160 washout_trap, #161 openmp parity) |
| Functional PRs + doc-only sweeps | Keep review scope narrow |

---

## Milestones

| Milestone | Issues |
|-----------|--------|
| P6 — HPC & GPU phase 2 | #154–#161 (done) |

Create via `./scripts/setup_project_board.sh`.

---

## Maintaining this board

After each merged PR:

1. Move issue(s) to **Done** and close on GitHub
2. Update the tables above
3. Link new PRs in the **In Review** row

When opening a new PR, reference the issue (`#160`) in the title or body.
