# GPUX campaign: what does the GPU buy on the science scenario?

Submitted 2026-09-04 to `gutibm-gpu-practice` (job def `gutibm-gpubench:11`, same
image digest as the p1 precision campaign). Command override in `job_command.sh`.

## Why
The p1 precision campaign (E1-E4) showed the device is outcome-interchangeable
with the host, but at 500 collapsing agents the whole cost is the 12.5M-cell
chemistry grid and agent-side kernels are idle. Question left open: on the
flagship science configuration, at the 1440-step horizon, does the device
(a) still produce interchangeable retention/diversity outcomes and (b) how much
wall time does it save when agents are dense enough to matter?

## Design
Base: `examples/diversity_paradox/input.json` (residents type 1 vs immigrants
type 2, continuous luminal immigration, FMM on), with:
- domain 1 mm x 1 mm x 0.1 mm (same grid as p1; the 2 mm example grid is 50M
  cells and exceeds the 13 GB job memory envelope / 512-line Route B scope)
- total_time 86400 s = 1440 steps at bio_dt 60
- HDF5 summary+provenance every step, agents dump every 360 steps (per-type
  counts = resident-retention proxy), no grid/lineage/genome dumps
- outputs exported to CloudWatch via `===BENCH_PRECISION_BEGIN` blocks (S3
  PutObject is denied to the instance role; HDF5 is not durable)

| Arm | Placement | Agents (t=0) | Immigration rate | Seeds | Timeout |
|---|---|---|---|---|---|
| GSH | host | 500 + 100 | 1/h | 61,62,63 | 4 h |
| GSD | device | 500 + 100 | 1/h | 61,62,63 | 4 h |
| GDH | host | 10000 + 2000 | 20/h | 61,62 | 6 h |
| GDD | device | 10000 + 2000 | 20/h | 61,62 | 6 h |

Estimates from p1 per-step costs (host 4.8 s/step, device 2.8 s/step chemistry;
agent side ~3 s/step host per 1e4 agents from s1): GS host ~1.9 h, GS device
~1.1 h, GD host ~3-4 h, GD device ~1.5-2 h. Cost at $0.526/h: ~$10 total.
Queue is one g4dn.xlarge, so serial wall ~12 h after the 9 queued p1 jobs.

## Analysis plan
- `gpu_precision --host GSH:.. --device GSD:..` (and GDH/GDD) on harvested
  exports; distributional verdicts + divergence onset.
- Resident retention = type-1 count at step 1440 / 500 (or 10000), from
  `===GPUX_META` blocks; compare host vs device medians against seed spread.
- Speedup = host wall / device wall per tier; compare against the 2.1x seen at
  p1 to see whether agent density moves it.
