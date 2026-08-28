# GutIBM Spec 14 — Prophage Induction as a Second Bacteriocin Killing Channel

**Status**: Proposed. Ready for implementation.
**Date**: 2026-08-28
**Branch from**: `main`
**Default state**: OFF. All behaviour gated behind `prophage.enabled = false`.

---

## 1. Motivation

Henrot, Debarbieux & Petit (bioRxiv 2026, doi:10.64898/2026.05.26.727859) screened
1,768 human fecal *E. coli* isolates, selected 30 bacteriocin producers, identified
74 bacteriocin genes, and tested induction of a λ-related coliphage. Two results
matter to us:

1. **Only DNA-damaging bacteriocins induce prophage.** E-type endonuclease colicins
   (E2/E7/E8/E9) and the microcin Mcc1229 induced a broad panel of temperate phage
   genera (5–120x and 21–89x fold increases in titer respectively). MccB17 did not
   induce the tested prophage. RNase colicins (E3, D) would not be expected to.
2. **The competitive outcome reverses between *in vitro* and *in vivo*.** An E-type
   endonuclease colicin producer outcompeted a λ-lysogen within **6 hours** in
   liquid culture, yet both populations persisted at **comparable levels over 10
   days in dixenic mice**.

Result 2 is the single cleanest validation target we have found for GutIBM: the same
two strains, two environments, opposite outcomes. Our spatial refugia plus
density limitation should reproduce the reversal **with no new mechanism at all**.

Result 1 motivates a genuinely new mechanism: prophage induction is a *second
killing channel* driven by the same DNA-damage variable we already track, plus a
horizontal-transfer channel (new lysogens among the producer's kin).

**Sequencing note.** Section 8 (the validation target) requires **no code changes**
and should be run **first**, on current `main`. If spatial structure alone
reproduces the reversal, that is a publishable result and prophage induction
becomes a refinement rather than a necessity. Do not implement Sections 5–7 before
running Section 8.

---

## 2. What the repo already has

Substantial machinery exists. This spec is mostly a wiring and state-tracking job,
not new subsystems.

| Existing | Location | Reuse |
|---|---|---|
| `PhenoState::SOS_INDUCED` | `src/core/types.h:57` | Target state on induction |
| `ReleaseMode::PHAGE_LYSIS` | `src/core/types.h:71` | Already distinguishes phage-mediated release |
| `BICluster::phage_induction_rate` | `src/core/types.h` | Per-generation spontaneous rate |
| `BICluster::phage_burst_size` | `src/core/types.h` | Currently unused |
| `BICluster::phage_lysogeny_rate` | `src/core/types.h` | Comment says "reserved for future HGT pathway" — **this spec is that pathway** |
| `BICluster::is_nuclease` | `src/core/types.h:89` | **This is the DNA-damaging flag.** Already set true for the nuclease colicin in `src/genome/plasmid.cpp:60` |
| `FixBacteriocin::check_phage_induction` | `src/fixes/fix_bacteriocin.cpp:102` | Extend with damage-driven term |
| `sos_cross_induction_rate` | `src/fixes/fix_bacteriocin.h` | Already couples local nuclease toxin to SOS |
| `Simulation::local_nuclease_toxin` | `src/core/simulation.cpp:1296` | Returns local BtuB-bacteriocin conc |
| `AgentTimers::sos_timer` | `src/core/agent.h:29` | Lysis delay countdown |
| `StepEvents::phage_inductions` | `src/core/step_events.h` | Counter already written to HDF5 |

**Three things are genuinely missing:**

**(a) Lysogeny is conflated with bacteriocin carriage.** Right now the only way an
agent can undergo phage-mediated lysis is to carry a `BICluster` with
`release_mode == PHAGE_LYSIS`. But a lysogen is not necessarily a bacteriocin
producer — in the target experiment the λ-lysogen is the *victim*, carrying a
prophage and no colicin. We need `is_lysogen` as agent state independent of
`bi_loci`.

**(b) DNA damage is instantaneous, not accumulated.** `check_sos_induction` uses
`sos_cross_induction_rate * nuclease_conc` as an instantaneous hazard. Prophage
induction is better described as a threshold on *accumulated* damage with repair,
which is what produces the switch-like λ response. We need a per-agent damage state
variable.

**(c) There is no free-phage pool.** `lyse_agent` releases toxin only. Without free
phage there is no lysogenization channel and no second killing channel — only
the producer's own suicide.

---

## 3. Biology being modelled

```
  colicin E (nuclease) enters target via BtuB
            |
            v
  DNA double-strand breaks accumulate  <---- repair (RecBCD/RecA)
            |
            v  (Hill threshold)
  SOS derepression --> if is_lysogen: prophage excision
            |                              |
            v                              v
  [existing] target dies            lysis + burst of free phage
                                           |
                        +------------------+------------------+
                        v                                     v
              infect nearby susceptible              decay / washout
                        |
              +---------+---------+
              v                   v
        lysogeny (~19%)      lytic death (~81%)
              |
              v
     new lysogen; gains superinfection immunity
```

The ecologically interesting consequence: a colicin producer that kills lysogens
**manufactures phage that can lysogenize its own kin**. Whether this helps or hurts
the producer at gut densities is a genuinely open question and the main scientific
payoff of the spec.

---

## 4. Parameters

Anchored where measurements exist; assumptions flagged explicitly.
task:44e9e0b4-233f-4123-8301-93399cb96056

### 4.1 Measured

| Parameter | Value | Units | Source | Confidence |
|---|---|---|---|---|
| λ latent period, in vitro | 40–80 | min | classical; growth-rate dependent | High |
| λ lysis time after induction | ~90 | min | heat-induction | Medium |
| Burst size, coliphage in vitro | ~60 (2–200) | PFU/cell | Ellis & Delbrück 1939 | High |
| **Burst size, λ in vivo (mouse gut)** | **12.1** | PFU/cell | De Paepe et al. 2016 | Medium |
| **Lysogenization fraction, in vivo** | **0.19** | — | De Paepe et al. 2016 | Medium |
| Lysogenization threshold | >= 2 | genomes/cell | Kourilsky 1973 | High |
| λ adsorption const. (no side fibers) | 1.3e-9 | mL/min | Shao & Wang 2008 | High |
| λ adsorption const. (with side fibers) | 9.9e-9 | mL/min | Shao & Wang 2008 | High |
| Spontaneous induction, in vitro | ~1e-8 | per cell/gen | review | Low |
| **Spontaneous induction, in vivo** | **1–2** | % of population | De Paepe et al. 2016 | High |
| E-type colicin induction | 5–120x | fold titer | Henrot et al. 2026 | Medium |
| Mcc1229 induction | 21–89x | fold titer | Henrot et al. 2026 | Medium |

**Use the in vivo values.** The in vitro / in vivo gap in spontaneous induction is
six orders of magnitude (1e-8 per cell/gen vs 1–2% of the population). De Paepe's
central finding is that λ carriage is *costly in the mouse gut because reactivation
is frequent*. Using the in vitro 1e-8 would eliminate the phenomenon we are
modelling. Set `spontaneous_induction_rate` from the in vivo figure.

Converting 1.5% of population per day to a per-second hazard:
`-ln(1 - 0.015) / 86400 = 1.75e-7 /s`. Use **1.75e-7 /s** as default.

Adsorption: `9.9e-9 mL/min = 1.65e-16 m^3/s`. Mucus will reduce this; see 4.2.

### 4.2 Assumptions — no measurement available

Flag these clearly in code comments. They are calibration targets, not knowns.

| Parameter | Proposed | Units | Justification |
|---|---|---|---|
| `hill_n` (damage → induction) | 2.5 | — | **No measured Hill coefficient exists.** λ's Cro/CI switch is known bistable, so `n > 1`. Sweep {1, 2.5, 5}. |
| `damage_half_max` | calibrate | damage units | No dose-response curve published. Calibrate so local nuclease colicin at killing concentrations gives induction comparable to the 5–120x observed titer increase. |
| `damage_repair_rate` | 1e-3 | 1/s | ~17 min repair timescale. Unmeasured in this context. |
| `phage_D_free` | 1e-12 | m^2/s | λ virion ~60 nm. Stokes–Einstein gives ~7e-12 in water; mucus retards **virions strongly** (unlike O2 — see Spec 12 Addendum A §4.1). Assume ~10x retardation. Sweep. |
| `phage_decay_rate` | 1e-4 | 1/s | ~2 h half-life. **Not measured in intestinal contents.** |
| `adsorption_rate_mucus` | 1.65e-17 | m^3/s | In vitro value reduced 10x for hindered diffusion. Unmeasured. |
| `prophage_carriage_cost` | 0.02 | frac. of mu_max | **No measurement found.** De Paepe shows carriage *is* costly, mechanism is frequent reactivation — which this model produces endogenously. Keep the direct cost small to avoid double-counting. Sweep {0, 0.02, 0.05}. |
| `lysogen_prevalence_init` | 0.3 | — | Fraction of founders that are lysogens. Poorly constrained for commensal *E. coli*. Sweep {0, 0.3, 0.7}. |
| `superinfection_immunity` | 1.0 | — | Full immunity of lysogens to the same phage. Well-established qualitatively. |

### 4.3 Which bacteriocins induce

Gate on the **existing** `BICluster::is_nuclease` flag; do not add a parallel flag.

| Bacteriocin | DNA-damaging | Induces? |
|---|---|---|
| Colicins E2, E7, E8, E9 | Yes (DNase) | **Yes** |
| Mcc1229 | Yes | **Yes** |
| Colicin E3, colicin D | No (RNase) | No |
| MccB17 | Gyrase inhibitor | **No** (measured, one prophage) |
| MccE492 | Pore-forming | No (assumed; untested) |
| Pore formers (E1, Ia, A) | No | No |

MccB17 inhibits DNA gyrase and might be expected to induce SOS, but was measured
*not* to induce the tested prophage. Follow the measurement. Note in code that this
is a single-prophage result and may not generalize.

---

## 5. Change 1 — Accumulated DNA damage state

**What**: Replace instantaneous nuclease hazard with an accumulating, repairable
damage variable.

```
// per agent, per step:
damage += k_damage_per_conc * local_nuclease_conc * dt      // accrual
damage -= damage_repair_rate * damage * dt                  // first-order repair
damage  = max(damage, 0)
```

**Files**:
- `src/core/agent.h` — add `Real dna_damage = 0.0;` to `AgentTimers`, or better to a
  new `AgentDamage` sub-struct to keep `AgentTimers` semantically clean.
- `src/core/agent_transfer.cpp` — **add `dna_damage` to the MPI transfer struct.**
  Follow the `is_nuclease` pattern at lines 98/170/240. Missing this will silently
  zero damage for agents crossing rank boundaries, producing rank-count-dependent
  results. This is the easiest bug to introduce here.
- `src/fixes/fix_bacteriocin.cpp` — update damage in `compute()`.
- `src/gpu/agent_update_kernel.cu`, `src/gpu/agent_pool_gpu.h/.cpp` — mirror on GPU
  with a `DeviceBuffer<double> d_dna_damage_`, following the `d_f_ferm_realized_`
  pattern from Spec 12 Change 2.

**Backward compatibility**: when `prophage.enabled = false`, keep the existing
instantaneous `sos_cross_induction_rate * nuclease_conc` path in
`check_sos_induction` untouched. Use if/else, exactly as Spec 12 Change 2 did for
the O2 boost.

**Test**: agent in constant nuclease field → damage rises to plateau
`k_damage_per_conc * conc / damage_repair_rate`; remove field → exponential decay
with time constant `1/damage_repair_rate`.

---

## 6. Change 2 — Lysogeny state and damage-driven induction

**What**: `is_lysogen` as agent state independent of `bi_loci`; Hill-function
induction on accumulated damage.

```
p_induce_per_s = spontaneous_induction_rate
               + max_induction_rate * damage^n / (damage_half_max^n + damage^n)
p = 1 - exp(-p_induce_per_s * dt)
```

**Files**:
- `src/core/agent.h` — add `bool is_lysogen = false;` and
  `Real prophage_cost_applied = false;` to `AgentFlags`.
- `src/core/agent_transfer.cpp` — transfer `is_lysogen` too (same warning as above).
- `src/fixes/fix_bacteriocin.cpp` — in `check_phage_induction`, add the
  damage-driven term. Note the current signature takes a `const BICluster&`; a
  lysogen may have **no** `bi_loci` at all, so add a separate
  `check_prophage_induction(Agent&, Real dt)` called for any agent with
  `is_lysogen`, rather than trying to force it through the per-locus loop.
- `src/fixes/fix_bacteriocin.h` — declare it; add config fields.
- Inheritance: `is_lysogen` must be inherited on division. Check wherever
  `just_divided` is handled / the daughter agent is constructed.

**Carriage cost**: apply once, like `microcin_penalty_applied`:
`if (is_lysogen && !prophage_cost_applied) { mu_max *= (1 - cost); ... }`

**Superinfection immunity**: a lysogen is immune to infection by the same phage.
Since we model one phage type, `is_lysogen == true` → not a valid infection target.

**Test**: lysogen in high nuclease field induces at high rate; non-lysogen in the
same field does not induce (but may still be killed by the existing colicin path).
Zero damage → induction at exactly `spontaneous_induction_rate`.

---

## 7. Change 3 — Free phage field and lysogenization

**What**: A diffusible phage field; induced lysogens release a burst; free phage
adsorb to susceptible agents and either kill or lysogenize.

This is the largest change. **Recommend deferring until Changes 1–2 and the
Section 8 validation are done** — the second killing channel may prove unnecessary.

**Representation**: reuse the existing chemical field / toxin-burst machinery
rather than building new infrastructure. Phage is a diffusing, decaying,
burst-released species — structurally identical to a colicin burst
(`ToxinBurstSource` in `qssa_solver.h`). Add species `PHAGE_LAMBDA`.

Caveat to record: phage are *particles*, not a concentration, and treating a burst
of 12 as a continuum field is a real approximation at low copy number. At the
densities where this matters (a handful of virions per voxel) a stochastic
particle representation would be more correct. Note it as a limitation; do not
build particles in the first pass.

**On lysis** (`lyse_agent`), if `is_lysogen`:
```
n_phage = phage_burst_size            // default 12.1 (in vivo)
```
released as a `ToxinBurstSource`-analogue with `phage_D_free`, `phage_decay_rate`.

**Infection** (new `FixPhage`, or extend `FixReceptor`):
```
// for each susceptible (non-lysogen, non-dead) agent:
p_adsorb = 1 - exp(-adsorption_rate * phage_conc * N_A * dt)
if (adsorbed):
    if (rng.bernoulli(lysogenization_fraction))   // 0.19
        is_lysogen = true                          // HGT channel
    else
        state = DEAD; step_events().phage_kills++  // second killing channel
```

Watch the unit conversion: `phage_conc` is mol/m^3 in the field but adsorption
constants are per-particle volumes. Multiply by Avogadro carefully, and sanity-check
that a burst of 12 particles in an 8e-18 m^3 voxel gives a plausible
concentration (12 / 6.022e23 / 8e-18 = 2.5e-6 mol/m^3 — same order as the
colicin concentrations already in the model, which is reassuring).

**New counters** in `StepEvents`: `phage_kills`, `new_lysogens`,
`prophage_inductions` (keep the existing `phage_inductions` meaning
bacteriocin-locus-driven release to avoid breaking existing output).

**Files**: `src/core/step_events.h`, `src/io/hdf5_writer.cpp`,
`src/fixes/fix_phage.{h,cpp}` (new), `src/fixes/fix_registry.cpp`,
`src/diffusion/qssa_solver.cpp`, `src/io/input_parser.cpp`,
`src/core/species` name table.

---

## 8. Validation target — RUN THIS FIRST, NO CODE CHANGES NEEDED

Reproduce the Henrot et al. in vitro / in vivo reversal.

| | In vitro analogue | In vivo analogue |
|---|---|---|
| Spatial structure | well-mixed, single voxel or tiny domain | full mucus domain |
| Density limitation | off | on (Spec 12) |
| Refugia | none | crypts on |
| Mucus retardation of colicin | minimal | pI-dependent (existing) |
| Washout | none | on |
| Expected | producer excludes lysogen within ~6 h | coexistence at comparable levels to 10 d |

Two strains: (1) E-type endonuclease colicin producer, `is_nuclease = true`;
(2) λ-lysogen, no bacteriocin. In the pre-Change-1 version strain 2 is simply a
sensitive non-producer, which is sufficient to test the reversal — the reversal is
claimed to be about *spatial structure*, not about phage.

**Metric**: ratio of strain abundances vs time; time to 100-fold exclusion.

**Interpretation guide**:
- Reversal reproduced without prophage → strong result. Spatial refugia and density
  limitation suffice; Changes 1–3 become a mechanistic refinement.
- Reversal *not* reproduced → check whether colicin diffusion in the well-mixed
  case is actually reaching all agents (the QSSA cutoff `toxin_cutoff` could
  artificially spare agents), before concluding the model lacks a mechanism.
- Coexistence in *both* conditions → the in vitro case is not well-mixed enough,
  or colicin potency is too low. Check `kd_colicinE_btuB` and the corrinoid pool;
  this connects to the open Kd sweep.

---

## 9. Suggested campaign

Once Changes 1–2 land:

| Axis | Values |
|---|---|
| `lysogen_prevalence_init` | 0, 0.3, 0.7 |
| `hill_n` | 1, 2.5, 5 |
| `prophage_carriage_cost` | 0, 0.02, 0.05 |
| spatial structure | well-mixed vs full mucus |

Primary question: **does prophage induction help or hurt the colicin producer?**
It kills lysogen competitors (helps) but manufactures phage that can lysogenize
the producer's kin, imposing carriage cost and creating a reservoir (hurts). The
sign of the net effect as a function of `lysogen_prevalence_init` is the result
worth publishing.

---

## 10. Priority and effort

| Change | Effort | Value | Priority |
|---|---|---|---|
| §8 validation on current `main` | **Low — no code** | **High** | **1** |
| Change 1: accumulated damage | Low | Medium | 2 |
| Change 2: lysogeny + Hill induction | Medium | High | 3 |
| Change 3: free phage field | High | Medium | 4 (defer) |

---

## 11. Limitations

- **The Hill function is invented.** No dose-response curve for prophage induction
  by colicins exists in the literature. `hill_n` and `damage_half_max` are
  calibration parameters and must be reported as such. Do not present the induction
  curve as literature-derived.
- Single phage type; real gut *E. coli* carry multiple prophages with distinct
  induction thresholds.
- Continuum treatment of phage particles is questionable at low copy number (§7).
- Adsorption and decay rates in mucus are unmeasured; both were scaled from broth
  values by an assumed factor.
- Prophage prevalence in commensal *E. coli* is not well constrained, so
  `lysogen_prevalence_init` is effectively a free parameter.
- Henrot et al. is a **preprint** (May 2026, 0 citations at time of writing). The
  key quantitative claims should be re-checked against the peer-reviewed version.
