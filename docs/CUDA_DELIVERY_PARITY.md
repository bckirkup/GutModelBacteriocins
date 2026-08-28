# CUDA parity for delivery-limited uptake

Scope: what has to exist before `metabolism.uptake_limit` can ship as
`delivery`, measured against the code as of #352 rather than argued. The parser
currently refuses `uptake_limit="delivery"` with `gpu_enabled=true`
(`src/io/input_parser.cpp`), and that refusal is the only thing standing between
a GPU configuration and a silently unfunded run. This document is the inventory
behind it and the two routes out; no route is implemented yet.

## What parity has to mean here

Not "the fields look similar". Delivery mode's whole purpose is that growth is
funded only from mass the field actually paid, so parity is an accounting claim:
for the same configuration and seed, host and device must agree on realized
removal per agent, on the rationing factor, and hence on population. A device
path that reproduces the concentration field to plotting accuracy while
mis-attributing realized removal is exactly the defect class of #338 (mass-exact
and negative) and #344 (priced at `P + G(z)` twice) — it looks like biology.

## Inventory

| Piece | Host | Device | Gap |
|---|---|---|---|
| Agent-side delivery prep and funded growth | `FixMetabolism::prepare_delivery_uptake`, `commit_delivery_uptake`, `post_chemistry` | `try_gpu_metabolism` returns `false` in delivery mode | None. Agent metabolism already falls back to host; this is a throughput cost, not a correctness gap. |
| Prescribed z-gradient reference profile | `z_gradient_reference` + `shift_z_gradient` | `shift_z_gradient_kernel` | None. Verified identical: `boundary_conc` at `iz == 0`, `initial_conc · exp(-z_rel/λ)` above, with the same `nz-1 → nz-2` profile clamp. |
| Epithelial boundary delivery (Dirichlet and flux) | `diffuse_bounded_z_delivery` | `solve_delivery_line` (`beta`, `flux_source`) | None. |
| **Per-cell first-order agent sink inside the implicit z solve** | `diffuse_bounded_z_delivery_with_sink` | absent — `solve_delivery_line` builds `1 + 2α` diagonals with no `s·dt` term | **Real.** This is the funded-uptake mechanism itself. |
| **Realized removal per cell, read back to fund biomass** | `split_delivery_sink_realized`, `add_delivery_reduction` | absent | **Real.** Without it there is nothing to fund growth from. |
| **Positivity rationing** | `run_delivery_rationing` — up to 4 local dilation retries, then a 12-iteration bisection, each iteration a full re-solve, with a collective negativity test (`collective_negative`) | absent | **Real**, and the expensive one. |

## The failure mode if the refusal is simply deleted

`ChemicalFieldGpu::apply_diffusion` is driven per species off `ChemicalSpec`; it
has no sink array and returns no realized removal. In
`run_chemistry_pipeline`, device diffusion is taken whenever reactions ran on
device, and the host delivery solve is reached only through the
`!result.diffusion_on_gpu` fallback. So permitting `gpu_enabled=true` with
delivery, and changing nothing else, would run carbon through a device solve
with the agent sink dropped: agents would be funded from mass the field never
removed, with no error and no clip. That is the pre-#332 unfunded behaviour
wearing delivery mode's provenance.

## Route A — permit GPU with host-forced delivery chemistry

Keep the numerics on the host, and make that a declared property of the run
rather than an accident of a fallback branch:

- allow `gpu_enabled=true` with `uptake_limit="delivery"`;
- force the delivery chemistry path to host, unconditionally, rather than
  relying on `reactions_on_gpu` being false;
- record it in `/run_provenance/` (e.g. a chemistry-placement field reading
  `host_forced_delivery`) so no archived run can be mistaken for a device solve;
- gate it with a test asserting that with `gpu_enabled=true` and delivery, the
  device diffusion path is *not* taken and realized removal is nonzero.

Cost: delivery runs get no chemistry speedup from the GPU; agent kernels and
non-delivery species keep theirs. Benefit: the shipped default can move to
`delivery` without breaking any GPU configuration, which is the decision this
blocks, and nothing scientific is left to a silent branch.

## Route B — port the delivery solve to the device

Add the `s·dt` diagonal term and per-cell realized-removal output to
`solve_delivery_line`, download realized removal for agent funding, and run the
rationing loop with the negativity test as a device reduction.

The cost structure is the thing to note before starting: rationing can require
up to 18 full three-axis solves in a single biology step, each followed by a
collective negativity decision. A port therefore cannot avoid a per-iteration
reduction — but that same arithmetic is why delivery mode is the most expensive
thing in the model on CPU, and why the port is worth doing eventually rather
than never.

## The parity gate, either way

Per the #334 lesson, a per-agent invariance test is not evidence. The gate is a
population-scale paired comparison: same config, same seed, host versus device,
asserting agreement on funded mass, rationing factor, and final population, and
asserting a *contrast* against a `uptake_limit="none"` arm so that a gate which
has stopped exercising delivery fails loudly instead of passing vacuously.

## Related default still open

`closure.enforce_reaction_clip` ships `false` while
`closure.enforce_delivery_realization` ships `true`, and a violation of either
already terminates with `/run_provenance/termination_cause_code =
closure_violation`. So the clip-halt provenance mechanism exists; what is
undecided is whether a clipped reaction should halt a shipped run by default.
Under delivery mode a clip means the field could not pay what was prescribed,
which is precisely the condition the funded accounting exists to expose.
