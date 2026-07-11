# GutIBM Project Board

Living kanban for open work on [GutModelBacteriocins](https://github.com/bckirkup/GutModelBacteriocins).

**Last updated:** 2026-07-11

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

## Current board state (Jul 2026)

The Jun 2026 queue (#40–#81, #25, #29, #33, #55) is **closed**. Post–GPU ROI work (#154–#161) is largely complete through PRs #162–#164.

### In Review

_None_

### Ready (recommended next bundles)

| Bundle | Issues | Track | Notes |
|--------|--------|-------|-------|
| **Science regressions** | #161 | `track:ci` | OpenMP stochastic parity |

### Backlog

| Issue | Title | Track | Priority |
|-------|-------|-------|----------|
| #161 | OpenMP parity on stochastic toxin-kill | `track:ci` | low |
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
| #160 | Metabolic washout trap long-horizon regression | (PR pending) |

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

1. **#161** — OpenMP stochastic parity

---

## What not to bundle

| Don't combine | Reason |
|---------------|--------|
| #160 + #161 | Different subsystems (biology regression vs OpenMP RNG) |
| Functional PRs + doc-only sweeps | Keep review scope narrow |

---

## Milestones

| Milestone | Issues |
|-----------|--------|
| P6 — HPC & GPU phase 2 | #154–#160 (done); #161 (backlog) |

Create via `./scripts/setup_project_board.sh`.

---

## Maintaining this board

After each merged PR:

1. Move issue(s) to **Done** and close on GitHub
2. Update the tables above
3. Link new PRs in the **In Review** row

When opening a new PR, reference the issue (`#160`) in the title or body.
