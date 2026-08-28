# Devin Implementation Brief — GutIBM Spec 12 Addendum A + Spec 14

**Date**: 2026-08-28
**Branch from**: `main`
**Suggested branches**: `spec12-addendumA-damkohler`, `spec14-prophage`

Two independent workstreams. Neither changes existing behaviour; both are
default-off. They can be implemented in parallel by separate PRs.

---

## Read first

| Document | Purpose |
|---|---|
| `/workspace/gutibm_spec12_addendumA_scale_gating.md` | Full Addendum A spec |
| `/workspace/gutibm_spec14_prophage_induction.md` | Full Spec 14 spec |
| `/workspace/spec12_da_report.md` | Reaction-diffusion analysis backing Addendum A |
| `/workspace/spec12_da_sweep.csv` | Sensitivity results |
| `/workspace/spec12_da_analysis.py` | Reference implementation of the Da solve |

---

## Ordering — please respect this

1. **Spec 14 §8 validation run.** No code changes. Runs on current `main`.
   Two strains, well-mixed vs full mucus. This may show that spatial structure
   alone explains the target phenomenon, which would change the priority of
   everything below it.
2. **Addendum A Change A1** — per-voxel Damköhler diagnostic (scalars only).
3. **Addendum A Change A2** — aggregate-radius sweep harness.
4. **Spec 14 Change 1** — accumulated DNA damage state.
5. **Spec 14 Change 2** — lysogeny state + Hill-function induction.
6. **Spec 14 Change 3** — free phage field. *Defer; may prove unnecessary.*

---

## Workstream 1: Spec 12 Addendum A

### The one thing to understand before coding

Spec 12 recorded that `agent_carbon_coupling` was a "no-op". We now know **why**,
and the reason is not what we assumed.

We had assumed the populations were too dilute. That is wrong. At `grid_dx = 2 µm`
a *single* agent in one voxel already represents a local density of
**1.25e11 cells/mL** — squarely in the dense-microcolony regime where the
literature reports oxygen limitation. Density is not the problem.

The problem is that `Da` scales as `R^2`, and one voxel gives `R = 1 µm`:

```
Da(single agent, dx=2µm) = 0.003
R_50 at that same density = 19.6 µm  = ~10 voxels radius = ~3,900 agents
```

**The gate is spatial coherence, not density.**

Consequence, and the reason this brief exists:

> **Do not "fix" the coupling by increasing `agent_carbon_coupling`, raising
> per-cell uptake, or removing the continuum VBF sink.** The coefficient is the
> right magnitude. Compensating for a geometry deficit with inflated per-cell
> metabolism would produce unphysical parameters that then contaminate every
> other result. If you find yourself tuning a rate constant to make the coupling
> visible, stop and flag it.

### Change A1 — per-voxel Damköhler diagnostic

Spec §5.1. Key points:

- Config struct `DamkohlerDiagnosticConfig`, `enabled = false` by default.
- **Use a coherence-aware length scale, not `dx`.** `R = dx/2` always reports
  `Da << 1` by construction and is useless. Compute the equivalent-sphere radius
  of the *occupied* volume in a neighbourhood of half-width `W = 5`.
- **Use the realized per-agent uptake, not the config max.** For O2 that is
  `q_consumption * max(mu_realized,0) + q_maintenance`, matching
  `src/diffusion/qssa_solver.cpp:496-497`. Reporting `Da` from `q_max` while the
  model applies a Monod-limited rate will overstate coupling exactly where it matters.
- **Implement the neighbourhood sum as three separable 1-D box-filter passes.**
  A naive triple loop at `W=5` is 1,331 reads/voxel and will dominate step time.
- Extend the existing agent loop in `qssa_solver.cpp` (~line 480–503); do not add
  a second pass over agents.
- Default to scalar summaries only (`da_max`, `da_p99`, `da_mean_occupied`,
  `frac_voxels_da_gt_0p1`). The per-voxel field is as large as a chemical field.

**The regression test that matters** (spec §5.1):
- 4-voxel-radius blob at 1 agent/voxel → `da_max < 0.05`
- 15-voxel-radius blob at the same density → `da_max > 0.5`

Same density, different extent, order-of-magnitude different `Da`. If both report
similar values, the coherence length is wrong and the diagnostic is worthless.
Please treat this test as the acceptance criterion for A1.

Also: with `damkohler.enabled = false`, output must be byte-identical to `main`.

### Change A2 — aggregate-radius sweep harness

Spec §5.2. Spherical clump of specified radius and density, fields on,
**growth/division/motility OFF** (if agents divide, clump geometry changes and the
measurement is meaningless). Sweep radius {2,5,10,20,40,80} µm × density
{1e10,1e11,1e12} cells/mL. Metric: `C(centre)/C(edge)` vs radius.

Expected crossings: ~22 µm at 1e11 cells/mL, ~7 µm at 1e12.

If the model's crossing is far from prediction, that is a **real finding** about
the field solver — likely related to the explicit-vs-implicit ordering issue we
already hit with the VBF O2 sink. Report it; do not tune it away.

---

## Workstream 2: Spec 14 — Prophage induction

### Start with §8. It needs no code.

Reproduce the Henrot et al. 2026 in vitro / in vivo reversal on current `main`:
a colicin producer excludes a sensitive competitor within ~6 h well-mixed, but
coexists for 10 days in the full spatial model. Two strains, two configs.

If spatial refugia plus density limitation reproduce this, it is a strong result
on its own and Changes 1–3 become refinement rather than necessity.

If it does **not** reproduce, check before concluding the model lacks a mechanism:
- is the well-mixed colicin actually reaching all agents, or is `toxin_cutoff` in
  the QSSA artificially sparing them?
- is colicin potency high enough? (`kd_colicinE_btuB` vs the corrinoid pool — this
  connects to the open Kd sweep.)

### Most of the machinery already exists

See spec §2 for the full table. Notably `BICluster::is_nuclease` is already the
DNA-damaging flag and is already `true` for the nuclease colicin in
`src/genome/plasmid.cpp:60`. **Gate induction on it; do not add a parallel flag.**
`BICluster::phage_lysogeny_rate` is commented "reserved for future HGT pathway" —
this spec is that pathway.

Three things are genuinely missing: (a) `is_lysogen` as state independent of
`bi_loci` (a lysogen victim carries no colicin), (b) *accumulated* rather than
instantaneous DNA damage, (c) a free-phage pool.

### Two traps

**1. MPI transfer.** New per-agent state (`dna_damage`, `is_lysogen`) **must** be
added to the transfer struct in `src/core/agent_transfer.cpp`. Follow the
`is_nuclease` pattern at lines 98 / 170 / 240. Omitting this silently zeros the
state for agents crossing rank boundaries and produces rank-count-dependent
results — a genuinely nasty bug to chase later.

**2. `check_phage_induction` signature.** It currently takes a `const BICluster&`
and is called from the per-locus loop. A lysogen may have **no** `bi_loci` at all.
Add a separate `check_prophage_induction(Agent&, Real dt)` called for any agent
with `is_lysogen`; do not force it through the locus loop.

Also: `is_lysogen` must be inherited on division.

### Parameter warning

Use the **in vivo** spontaneous induction rate, not the in vitro one. The gap is
six orders of magnitude (1–2% of population per day vs 1e-8 per cell/generation).
De Paepe et al. 2016 found λ carriage is costly in the mouse gut *precisely
because* reactivation is frequent. Using 1e-8 would delete the phenomenon.
Default: **1.75e-7 /s**.

Likewise use in vivo burst size (12.1 PFU/cell) not the classical in vitro ~60.

The Hill function for damage → induction is **invented**. No dose-response curve
exists in the literature. `hill_n` and `damage_half_max` are calibration
parameters; please comment them as such in the code so nobody later mistakes them
for measured values. Sweep `hill_n` in {1, 2.5, 5}.

### Change 3 — defer

The free phage field (spec §7) is the largest piece and may prove unnecessary
depending on the §8 result. When it happens, reuse the `ToxinBurstSource`
machinery rather than building new infrastructure — a phage burst is structurally
identical to a colicin burst. Note in the code that treating ~12 virions as a
continuum concentration is a real approximation at low copy number.

---

## Backward compatibility (both workstreams)

- All new config flags default off/zero.
- Existing code paths preserved behind if/else, following the pattern Spec 12
  Change 2 used for the legacy O2 boost.
- Both CPU and GPU paths must be updated together; GPU per-agent buffers follow
  the `d_f_ferm_realized_` pattern from Spec 12.
- Byte-identical output with all new features disabled — please verify explicitly.

## Questions to raise rather than guess

- If the A1 diagnostic shows `Da > 0.1` in a meaningful fraction of occupied
  voxels in *existing* campaign runs, tell us — that would mean coupling is
  already active in runs we interpreted as uncoupled, and several conclusions
  would need revisiting.
- If the A2 sweep crossing disagrees with prediction by more than ~2x, stop and
  report before adjusting parameters.
- If §8 shows coexistence in the well-mixed condition too, that is a colicin
  potency or QSSA cutoff issue, not a missing mechanism. Flag it.
