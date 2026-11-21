# GutIBM Project Board

Living kanban for open work on [GutModelBacteriocins](https://github.com/bckirkup/GutModelBacteriocins).

**Last updated:** 2026-11-21

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

## Current board state (Nov 2026 — Cycle 2 Empirical Calibration & Multiscale Dispatch)

The Sep–Nov 2026 campaigns and Layer 2/3 multiscale architecture waves (#307–#473) are **merged**, delivering Spec 12 delivery/ROS campaigns, Gate C single-source & ecological refinement adjudication (#429), Layer 2/3 `MucusSegment` + `CompartmentChain` + live audit patches (#453–#463), and the Cycle 2 210-DOI empirical calibration synthesis & Spec 13–14 audit (#456, #457, #473). Active Cycle 2 work is tracked across **#464–#483** (with automated S3 collector index maintained in #446).

### Active Cycle 2 Queue (#464–#483)

| Issue | Track | Priority | Scope |
| :--- | :--- | :--- | :--- |
| #464 | Spec 14 (`GAP-01`) | `critical` | Decouple BtuB receptor corrinoid shielding from obligate B12 Monod growth arrest + carriage-matched ColE2/ColIa/ColB controls |
| #465, #476 | Spec 13 (`GAP-02` + I/O) | `critical` | Multi-slot (`0..N`) live patches, Active `try_deposit`, split `sos_lysis` vs `toxin_kill` ledger, and `segment_io.cpp` fixes |
| #466 | Spec 13 (`GAP-03`) | `important` | Configurable `live_config_for` per-PatchType crypt flow refugia and oxygen/ROS inheritance |
| #467 | Spec 13 (`GAP-05`) | `important` | Layer 3 mucosal occupancy sensitivity (`reattach_p0` vs luminal growth) and post-48h window vs 72h endpoint separation |
| #468, #477 | Bacteriocin (`GAP-06`) | `important` | Per-division pair SOS normalization (`sos_lysis_prob`) and `D_free_colicin` vs per-plasmid `diff_coeff` contract |
| #469, #470 | Spec 12 (`GAP-07`/`08`) | `important` | Wet-mass (`1,100 kg/m^3`) documentation, live-slot `delivery` sync, and Dirichlet vs Robin O2 presets + GPU maintenance guard |
| #471, #474, #481 | Observation & Analysis | `important` | Optical-slice `depth_window_m`, abundance-weighted retention, Kirkup & Riley pellet/cage LOD censoring, and division-wave slope guards |
| #472, #478, #479, #480 | Numerics, Metabolism & HPC | `important` | Layer-1 checkpoint RNG serialization, QSSA transient/wall-normal drift bounds, implicit iron VBF sink, GPU funded respiration, and MPI/GPU memory |
| #475 | Spec 14 (`GAP-10`) | `gated (#464)` | Temperate prophage carriage, DNA-damage Hill induction, superinfection immunity, and lysogen conversion |
| #482 | Campaign Governance | `important` | Stage D decision-record axis enforcement in `preflight.py` and HDF5 `container_image_digest` authentication in `analyze.py` |
| #159 | GPU FMM Far-Field | `backlog` | Sub-quadratic CPU M2L delivered in #164; full device-side octree traversal remains low-priority backlog |

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
