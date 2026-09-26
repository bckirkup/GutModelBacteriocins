# Checkpoint Analysis Tool

Command-line analysis for Spec-4 GutIBM HDF5 closed-restart checkpoints
(`restart.directory/step_NNNNNN.h5`). Designed for large runs: stream one file at
a time, never load the full campaign into memory.

## Install

```bash
cd python
pip install -e ".[viz]"   # matplotlib needed for --snapshot / --timeseries
```

Entry point: `gut-ibm-analyze` (also `python -m gut_ibm_tools.checkpoint_analyze`).

## Workflow

1. **Summarize** all checkpoints → small CSV (~1 MB for dozens of files)
2. Inspect trends (upload CSV, or run `--timeseries`)
3. **Snapshot** only the interesting steps for spatial plots
4. Iterate

### Mode 1 — Summarize

```bash
gut-ibm-analyze --summarize \
  --checkpoint-dir practice_outputs/baseline_seed1001/ckpt \
  --output summary.csv
```

Writes:

- `summary.csv` — one row per readable checkpoint
- `summary_meta.json` — column list, skip count, source directory

Options:

| Flag | Meaning |
|------|---------|
| `--stride N` | Keep every Nth file |
| `--max-checkpoints N` | Cap how many files to read |
| `--quiet` | Hide per-file progress |

Progress goes to stderr. Corrupt or truncated files (and AWS multipart leftovers
like `step_000030.h5.6cF6bCBE`) are skipped with a warning.

### Mode 2 — Snapshot

```bash
gut-ibm-analyze --snapshot \
  --checkpoint practice_outputs/baseline_seed1001/ckpt/step_002520.h5 \
  --output snap_2520
```

`snap_2520/` must resolve under the current working directory (path-safety rule).
Typical figures:

- `agents_xz.png`, `agents_xy.png` — cells by strain type
- `density_vs_depth.png` — z histogram by type
- `heatmap_<species>.png` — mid-y chemical slices (when present)
- `kill_zone.png` — bacteriocin field + agent overlay

### Mode 3 — Timeseries

```bash
gut-ibm-analyze --timeseries \
  --summary summary.csv \
  --output timeseries
```

Produces population, kill-event, chemistry, spatial-metric, and optional
population-vs-toxin phase plots under `timeseries/`.

## Schema mapping (aspirational → Spec 4)

Closed restarts use the same layered layout as analysis trails:

```
/agents/step_NNNNNN/{id,type,state,x,y,z,mu_realized,mu_max,...}
/summary/step_NNNNNN/{time,step,n_total,events/*,chem/*,...}
/grid/step_NNNNNN/{carbon,iron,acetate,bacteriocin_*,...}
```

| Requested concept | Actual source |
|-------------------|---------------|
| positions / strain | `x,y,z` / `type` |
| growth rate | `mu_realized` |
| alive | `state != 3` (DEAD) |
| glucose | `carbon` |
| colicin field | `bacteriocin_BtuB` (+ FepA/CirA/FhuA) |
| kills / divisions | `summary/.../events/*` |
| O₂ / AI-2 grids | only if those species were enabled |

Missing datasets are omitted from the CSV / skipped in plots — the tool
auto-discovers what each file contains.

Not persisted today (skipped in v1): per-agent energy, motility mode,
bacteriocin induction flag, CDI concentration grid.

## Library API

```python
from gut_ibm_tools.checkpoint_scan import discover_checkpoints
from gut_ibm_tools.checkpoint_summary import extract_checkpoint_row, extract_many

paths = discover_checkpoints("path/to/ckpt")
row = extract_checkpoint_row(paths[0])       # opens, extracts, closes
rows, skips = extract_many(paths)            # stream all
```

Spatial metrics reuse `gut_ibm_tools.analysis` (Hopkins, monochromatic score).

## Memory / performance notes

- One `h5py.File` open at a time; agents + needed grids for that step only.
- A ~500 MB Tier-2 restart typically yields a ~100–150 column summary row in
  tens of seconds on a laptop (grid means dominate).
- Prefer `--summarize` with `--stride` for first-look on 80+ file campaigns.
