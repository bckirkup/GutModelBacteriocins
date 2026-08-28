# Selecting the lysis prior against in-vivo displacement, not fitting it

Design document for work order item 3. Companion to `SPEC13_MULTISCALE.md`
(architecture), `SPEC13_IMPLEMENTATION_REVIEW.md` (seams), and
`DELIVERY_ROS_CAMPAIGN.md` (where the candidate priors came from). Nothing here
is implemented. Every code claim below is a file/line reference or an arithmetic
consequence of a shipped default, not a statement of intent.

The premise of the work order is that the per-generation lysis probability is
not separately observable in vivo, so the 1% / 2% / 5% culture-derived priors
must be *selected* by which one reproduces the Kirkup & Riley 2004 displacement
kinetics, rather than fitted as a gut lysis rate. This document is about what
has to be true of the model before that selection means anything. Three of the
findings below change the campaign design, and one of them changes what the
`k_ROS_funded` coefficients mean.

Two prerequisites have since been resolved and are folded in below. The
comparison observable is **per-animal stool**, sampled from individual faecal
pellets that were subdivided for plating (§4d) — which forces Layer 3 and, more
consequentially, makes the target a *categorical dominance sequence* with a 1%
detection limit rather than a continuous displacement curve. And the resistant
strain is **already configurable** through per-strain `receptor_expression`
(§4a); the earlier audit finding that it was not is withdrawn.

## 1. The release chain, as implemented

Read as of `src/fixes/fix_bacteriocin.cpp` on `main`:

| Step | Where | Behaviour |
|---|---|---|
| Induction hazard | `check_sos_induction()` | four additive rates: basal, post-division, nuclease cross-induction, ROS; `p = 1 - exp(-rate_total·dt)` |
| Gating | `has_release_mode(agent, SOS_LYSIS)` | only carriers of an SOS-lysis locus are ever evaluated — a plasmid-free strain has zero lysis hazard of any kind |
| Delay | `k_sos_lysis_delay_s = 300.0` | induced cells die 300 s later; phage route 60 s |
| Death + release | `post_step()` → `lyse_agent()` | marks `DEAD`, books `mortality_lysis`, records `ProvenanceCause::LYSIS`, and *only here* emits `ToxinBurstSource` |
| Kill of sensitives | `FixReceptor::compute_kill_prob()` | receptor-bound toxin hazard, immunity factor, provenance-recorded |

Two consequences matter for the campaign. First, **lysis is the sole release
path for every colicin in the library except MccV** (`plasmid.cpp`: ColE1,
ColE2, ColM are `SOS_LYSIS`; ColB and ColIa are `PHAGE_LYSIS`; only microcin V
is `CONTINUOUS`). So the lysis prior is not a mortality knob with a side effect
on toxin — it *is* the dose schedule. Displacement is monotone in it only
because two opposing effects are bundled: more lysis is more producer death and
more toxin. That is exactly why the prior is selectable from displacement at
all, and also why the selection is only meaningful if the release parameters it
multiplies are right (§4).

Second, the nuclease cross-induction route (`sos_cross_induction_rate = 1e3` per
mol/m³, ColE2 is `is_nuclease`) makes the hazard *state-dependent*: lysis
releases nuclease colicin, which induces SOS in neighbours, which lyse. Within
a producer cluster the realized per-generation lysis is therefore not the
configured prior; it is the prior plus a density-dependent term. The prior must
be reported as measured (`mortality_lysis / divisions`), never as configured.

## 2. The model already lyses at roughly the top of the in-vitro band, uncited

Shipped defaults (`src/fixes/fix_bacteriocin.h`):

```
sos_lysis_prob = 0.01      // "1% per division"
sos_basal_rate = 1.0e-6    // 1/s, spontaneous
```

Both are charged to producers *before any oxygen is funded and with
`oxygen.k_ROS = 0.0`*, i.e. in the configuration we ship today:

- **Post-division.** `fix_metabolism.cpp:893-894` sets `just_divided` on *both*
  the mother and the daughter, and `simulation.cpp` clears it at the start of
  each step, so each division draws the 1% hazard twice — once per product.
  Over a cell's life it is marked once when born and once when it divides, so
  the per-generation contribution is **≈2%**, not the 1% the comment claims.
- **Basal.** 1e-6/s is a per-*time* hazard. At the measured `T_gen ≈ 2.3 h`
  uncrowded that is `1 - exp(-1e-6 × 8280) = 0.82%` per generation, rising to
  1.2% at the `T_gen ≈ 3.3 h` measured at 80 founders.

Total shipped lysis is therefore **≈2.8% per generation** for a producer,
before cross-induction and before any ROS term. This has three consequences:

1. The `k_ROS_funded = 6.2e11 / 1.24e12 / 3.1e12 mol⁻¹` coefficients from #339
   are *increments*, not the model's per-generation lysis probability. Setting
   the "1%" coefficient produces ≈3.8%/generation, not 1%. Nothing in
   `DELIVERY_ROS_CAMPAIGN.md` is wrong — it says "1% per generation against the
   measured Q", which is what the coefficient buys — but a campaign that labels
   its arms 1/2/5% while running those coefficients would mislabel every arm.
2. The 1–5% in-vitro band is *already spanned by the shipped defaults*, through
   two routes that have no citation in the tree and no oxygen dependence. The
   funded-ROS coefficient is not the only, or even the dominant, thing the
   selection would be probing.
3. `sos_lysis_prob` is the one route that is per-generation *by construction*
   and needs no oxygen at all — which is what makes the campaign in §5 runnable
   before CUDA parity exists.

## 3. Sweep the per-division route, convert afterwards

Recommended: run the selection on `bacteriocin.sos_lysis_prob` with
`sos_basal_rate = 0` and `oxygen.k_ROS = 0`, so the swept quantity is the total
per-generation lysis probability and nothing else contributes an unlabelled
share. Then convert the selected value to `k_ROS_funded` post hoc via the
measured `Q_O2_per_gen`, and confirm with a single funded-mode arm on CPU.

Why not sweep the funded-ROS coefficient directly:

- `metabolism.uptake_limit=delivery` is required for a funded flux to exist at
  all, and the parser refuses it with `gpu_enabled=true` (no CUDA parity), so
  the whole campaign would be CPU-only *and* would be selecting a coefficient
  whose shipped default we have not yet decided (`AGENTS.md`, Partially
  decided).
- The conversion is only valid where funded O₂ per generation equals the
  uncrowded `1.63e-14 mol`. It does not: `Q` rises to `2.31e-14` by 20 founders
  and bloom-density funding is **1.5% of the analytic draw, all of it
  maintenance** (`DELIVERY_ROS_CAMPAIGN.md` §3). A funded-ROS run at the
  densities a displacement experiment produces therefore realizes *far less*
  lysis per generation than its coefficient nominally encodes. Report that gap;
  do not correct for it, and do not select a coefficient whose realized effect
  is set by crowding we have not characterized.

So: select a per-generation probability, publish the implied coefficient as an
implication with its validity range attached, and keep the funded route as a
confirmation arm rather than the sweep axis.

## 4. Prerequisites — each of these would be absorbed into the selected prior

A prior selected against displacement absorbs every error in the machinery
between lysis and the death of a sensitive cell. Four were identified; (a) and
(d) are now resolved, (b) and (c) remain open blockers.

**(a) The resistant strain is configurable — the audit was wrong, and the real
gap is the cost, not the genotype.** `experiments/rps_campaign/rps_spec_audit.md`
§C states that `InitialStrain` has no receptor-genotype field and that resistance
can arise only from `FixMutation`. That is false, and the AGENTS.md landmine row
repeating it has been corrected. A per-strain receptor genotype is parsed,
validated, tested, and applied on the shipped path:

- `src/io/config_json.cpp::parse_receptor_expression_object` accepts a
  `receptor_expression` object on each entry of `initial_strains` (aliases
  `receptor_genotype`, `receptors`), rejects unknown receptor names, and
  requires each value finite in `[0, 1]`.
- `Simulation::create_strain_agent` writes it into `receptor_expr_base`,
  `receptor_expr`, and `genome.receptor_expression`, and tags the founder
  `PhenoState::RESISTANT` when `receptor_expression_is_resistant()` (`< 0.2`).
- `tests/fixtures/parser_strains.json` already carries `{"BtuB": 0.0,
  "FepA": 0.35}` and `tests/test_config_ingestion.cpp` asserts both land.
- Division copies `receptor_expr_base` to the daughter (`fix_metabolism.cpp`),
  `agent_transfer.cpp` ships it across ranks, and `restore_receptor_fields`
  restores it from checkpoint. The genotype is not lost by MPI or restart.

ColE1 and ColE2 both enter through BtuB, and the kill hazard is
`kill_rate × occupancy × immunity` with `occupancy` *linear* in
`agent.receptor_expr[BtuB]` (`fix_receptor.cpp::toxin_occupancy`). So
`{"BtuB": 0.0}` is a founder that is completely resistant to both colicins by
the same mechanism as the paper's R strain — Kirkup & Riley's R is a
spontaneous `btuB` mutant, not an immunity carrier. **The C/S/R system is
constructible today.**

What is *not* settled is the cost, and the RPS cycle stands or falls on it: S
beats R only if losing BtuB costs something. The cost route is real but weak at
shipped defaults. `fix_metabolism.cpp` computes
`Km_b12 = km_b12 / (expr_btuB · lig_aff)` with `expr_btuB` clamped to a floor of
`0.01`, so a knockout inflates `Km_b12` by 100× — but only that, no matter how
far below `0.01` the genotype goes. Against the shipped `km_b12 = 1e-6` and
`b12_initial_conc = 1e-3 mol/m³` the knockout is still corrinoid-saturated:
`monod_b12` falls only from 0.999 to 0.909, a ~9% growth penalty. Three
consequences for the campaign:

1. The resistance cost must be *set by corrinoid supply*, not by an added
   `mu_max` penalty — an extra `mu_max` cost charges resistance twice, which is
   the double-counting the audit correctly warned about.
2. `b12_initial_conc` therefore becomes a campaign axis, and the S-beats-R leg
   has to be demonstrated at population scale before any RPS claim, not
   assumed from the Km algebra. That is what the probe in §4a-probe set out to
   measure, and §4a-result records that it could not: the S-beats-R leg remains
   undemonstrated at population scale, so no RPS claim is available yet.
3. Two bookkeeping caveats for strain accounting: `PhenoState::RESISTANT` is
   the same field that later holds `SOS_INDUCED`/`DEAD`, so classify strains by
   `receptor_expr[BtuB]` or agent `type`, never by `state`; and `FixMutation`
   writes `receptor_expression` at 1e-7/step, so a selection arm that wants a
   fixed three-strain composition must say what it does with de-novo
   resistance rather than let it accumulate silently.

**(a-probe) Population-scale check of the above.** The claims in (a) are read
from source; per-agent reasoning about a kill hazard is exactly the class of
claim that #334 established does not transfer to a population. The probe under
`experiments/rps_probe/` (config generator and analysis committed; outputs are
not) asserts three contrasts on 12 local serial arms at 2× the shipped carbon
amplitude, chosen because the corrected ladder puts patch capacity at 122 agents
from 100 founders at 1× and a fitness difference has no headroom to express
itself there:

- **A1** with a ColE1 producer present, the `{"BtuB": 0.0}` founder ends above
  the otherwise-identical sensitive founder;
- **A2** in a null arm with the identical three founders and no ColE1 anywhere,
  that contrast disappears — otherwise the separation is not colicin;
- **B** the no-producer R/S ratio is non-increasing as `b12_initial_conc` falls
  through `1e-3 → 1e-4 → 3e-5 → 1e-5`, which is the resistance cost the RPS
  cycle needs and locates the corrinoid level at which it becomes real.

A FAIL on any of these is a finding about the model, not an arm to retune.

**(a-result) Ran at `3c166c1`, serial, 12 arms. A1 FAIL, A2 FAIL, B PASS but
the PASS is an artifact. The probe did not test what it was built to test, and
why it could not is the finding.**

Two regimes, split by seed and not by treatment. Every arm loses most of its
founders in the first ~2 h — 120–180 founders down to ~30 live by step 120,
against 134–745 cumulative `outflow_boundary` — and then the patch either
escapes or does not. Cumulative divisions are bimodal with nothing in between:
22–297 in the collapse regime, 2383–7629 in the escape regime. Whatever few
lineages survive the crash own the patch, so the final composition is a founder
lottery.

| Arm | C | S | R | divisions | colicin kills |
|---|---:|---:|---:|---:|---:|
| `A_three_strain_s20260901` | 0 | 0 | 2 | 84 | 0 |
| `A_three_strain_s20260902` | 0 | 3455 | 0 | 3851 | 0 |
| `A_three_strain_s20260903` | 713 | 6252 | 2 | 7629 | 11 |
| `A_null_no_producer_s20260901` | 7 | 0 | 395 | 553 | — (no producer) |

Three identical-treatment seeds produce R-only, S-only, and C+S. **The null arm
— no ColE1 anywhere in the config — produces the single largest R/S separation
in the whole probe, 395 against the producer arms' 0.67.** That is A2 failing,
and it is unambiguous: at this scale, type composition is set by which founders
happened to survive, not by colicin. Nothing about the `{"BtuB": 0.0}` genotype
can be concluded from A1's FAIL, because the treatment it was contrasting
against barely exists (below).

**The producer treatment is nearly inert: 11 receptor-mediated kills against
7629 divisions**, in the one arm where the producer grew at all, and zero
everywhere else. This is (b) measured rather than argued — `retardation = 50`
puts `D_eff` 72× below the literature analogue, so the kill zone is ~1/72 the
area, and a producer at these densities almost never has a sensitive cell
inside it. Fixing (b) is therefore a prerequisite for *any* colicin efficacy
arm, not just for the lysis prior's accuracy.

**B's PASS must not be quoted; the assertion was wrong.** The monotone R/S
sequence (103.0, 7.0, 4.0, 1.0) is the same seed split: the `s20260901` arms are
R-only, the `s20260902` arms are S-only, and one strain is extinct in 8 of the
12 arms, so most of those "ratios" were `n3 / max(0, 1)` — an arbitrary number
divided by a strain that is not there. This is the failure mode the ladder's λ
arms already showed once, a monotone single-seed trend that is not a treatment
effect, and here it passed the gate. `analyze.py` now refuses a ratio against
an extinct denominator and reports `INSUFFICIENT` with the usable-arm count
instead of averaging it away; re-run against the same outputs, all three
contrasts are FAIL or INSUFFICIENT, which is the honest reading.

There is, however, a **controlled** signal hiding inside the B group, because
seed fixes founder placement across the corrinoid axis. Within `s20260901`
(R-dominated) the final population falls `206 → 14 → 8 → 2` as
`b12_initial_conc` falls `1e-3 → 1e-4 → 3e-5 → 1e-5`; within `s20260902`
(S-dominated) it does not move at all: `2360 → 2739 → 2144 → 3463`. A
BtuB-null population is corrinoid-limited over this range and a BtuB-intact one
is not, which is the direction the `Km_b12` clamp predicts. Treat it as a
hypothesis with n=1 seed per regime, not as the measured cost: it is confounded
with the escape/collapse split, and the two regimes are not the same population
being compared.

**What a probe with power has to change.** Not the arm labels — the scale. The
patch is at the self-sustainment threshold by construction
(`CARBON_LADDER_CAMPAIGN.md`: 194 divisions against 190 exports at the shipped
default), and a marginal patch destroys composition observables before it
destroys population ones, because it hands the whole patch to a handful of
lineages. Concretely: (i) the founder crash has to be removed or measured —
either a domain with enough capacity that 120 founders are not competing for a
few escape sites, or explicit reseeding, which is Layer 2's job; (ii) seeds must
be ≥5 per arm with the paired within-seed contrast as the statistic, never a
ratio of between-arm means; (iii) (b) must be fixed first, or the producer arm
is a null arm with extra steps. **This is a second, independent argument for the
same conclusion §4d reaches from the data side: strain composition is not
observable at patch scale, and the RPS comparison needs Layers 2 and 3.**

**(b) Colicin transport is 72× too slow, and not configurable.** ColE1/E2
`retardation = 50` is hardcoded in `plasmid.cpp` and copied into the burst;
`D_eff = 4e-11/50 = 8e-13 m²/s` against a literature analogue of
`7e-11/1.2 = 5.8e-11` (audit §B). Kill-zone radius scales ~8.5×, area ~72×. The
`retardation_basic/neutral/acidic` config keys have no live consumer. `burst_size`
(5e4 molecules for ColE2) is hardcoded the same way. **These are larger levers
on colicin efficacy than the lysis prior itself.** Selecting a lysis prior at
`retardation = 50` would compensate a transport error with a biological
coefficient — the same error class as calibrating `Q_O2_per_gen` against a
grid-dependent divisor (`DELIVERY_ROS_CAMPAIGN.md` §5). Requirement: make
`retardation` and `burst_size` configurable per BI locus, source both, and run
the selection at *both* the current and the literature value. If the admissible
prior moves between the two, the campaign has measured the transport parameter,
not the lysis rate, and must say so.

**(c) A single patch cannot host a weeks-long displacement experiment.** With
ROS off, every loss in the delivery/ROS campaign was boundary export, and the
corrected carbon ladder puts the *shipped* carbon default at the threshold of
net growth (`CARBON_LADDER_CAMPAIGN.md`: peak 122 from 100 founders, 194
divisions against 190 exports). A patch that is marginal over 12 h cannot be
run for the 2–4 week displacement horizon; the population is not persistent
enough for a ratio between two strains to mean anything at the end. **Layer 2
reseeding (Phase 1 of `SPEC13_MULTISCALE.md`) is a hard prerequisite**, and it
is the reason this campaign is a Spec 13 item rather than a patch-scale one.

**(d) The comparison observable is settled: per-animal stool, and it is
categorical.** Kirkup & Riley sampled faecal pellets taken directly from each
mouse — not mucosa — and individual pellets were subdivided for plating. Two
things follow, and the second is the one that reshapes the campaign.

*It forces Layer 3.* A mucosal ratio would have been satisfiable at Layer 2. A
stool count is not: the patch's only export channel books `outflow_boundary`
and discards the cell, so today the model has no stool observable at all. The
luminal compartment of `SPEC13_MULTISCALE.md` is a hard prerequisite, and the
observation model in that document (§ *Stool observation model*) — animal, cage,
pellet, subdivision, detection limit — is the thing that has to exist before a
selection statistic can be computed.

*It is not a time series of ratios.* The paper is explicit: "Only rarely was a
mixed population recovered from individual faecal pellets; at any given time
each mouse was dominated by a single strain. Detection limits for a second
strain were 1% of the population." Figures 1–3 are grids of mouse × half-week
coloured by *dominant strain*. So there is no continuous displacement curve to
fit and nothing to digitize off a y-axis; the data are a censored categorical
state sequence per animal. What is available to select against, with the
numbers as published:

| Statistic | Colicin E1 | Colicin E2 |
|---|---|---|
| Cages × mice | 12 × 3 = 36 | 12 × 3 = 36 |
| Sampling | half-weekly, 12 weeks | half-weekly, 12 weeks |
| Dominance transitions | 123 (">98 putative") | 43 |
| Mixed pellets | rare | rare |
| Second-strain detection limit | 1% | 1% |
| C fate | no cage fixed on C; C eliminated in all but two cages by week 7 | — |
| Early dynamics | — | 4 of 36 mice displaced within one week |
| Directions resolved (one-tailed Wald, 95%) | C→S, R→C, S→R all rejected as null | C→S at 98.7%, R→C at 97.7%; S-vs-R insufficient |
| Reappearance after loss | R 9×, S 9×, C 14× | S never reappeared once eliminated |
| Steady density | ~1e6 CFU/g faeces, stable >4 weeks | same |
| Reintroduction assay | — | S invaded and displaced R in 11 of 16 mice, Fisher P = 0.055 |

The transition-count contrast between the two producers (123 vs 43, same
hosts, same protocol, different colicin) is the single most useful number here:
it is a *rate* difference that a lysis prior can plausibly move, measured
within one experimental system, and it does not require any absolute CFU
calibration to compare against.

The binding caveat: at 12 cages the transition count is the aggregate of a
strongly non-independent design. The paper's own analysis presumes exponential
cooperativity among co-caged mice and homogenizes mouse identity into a 10×10
matrix for the Wald test. Any model statistic must be computed with the same
nesting (pellet within animal within cage) or it will overstate its own
precision — and the reappearance counts show that strains persist below the 1%
detection limit, so a model that lets a strain go extinct at zero is not
measuring the same process.

**Where a reappearing strain comes from, and why it is a confound.** Per the
experimentalist: a mouse could carry a small percentage of a second or third
strain, but the refugium for reinoculation was *the other mice in the cage*,
not another site within the same mouse. That makes the coprophagic exchange
kernel the mechanism that generates the observable rather than a nuisance term
to be kept small, and it creates a degeneracy: transfer rate and lysis prior
both move the transition count, in the same direction, and cannot be separated
by that statistic alone. The exit is to pin them on different observables —
reappearance counts (R 9×, S 9×, C 14×) and C's elimination from 10 of 12 cages
by week 7 constrain transfer and producer mortality; the transition count is
then left to the prior. Fix and report the transfer rate against the
reappearance counts *before* the prior sweep starts, not alongside it.

## 5. The selection procedure

Deliberately structured so that the answer can be "none of them", and so that a
null arm can falsify the observable itself.

**Arms.** `sos_lysis_prob ∈ {0, 0.005, 0.01, 0.02, 0.05, 0.10}` (per generation,
with basal and ROS routes zeroed) crossed with `retardation ∈ {50, 1.2}` = 12
arms. The 0 arm is the null: with release lysis-gated, it emits no toxin at all.
The 0.10 arm brackets above the in-vitro band so that a monotone-but-unbounded
response is visible as such.

**Founders.** Three strains per animal-equivalent, matching the paper's design:
C carrying ColE1 (or ColE2), S with default receptors, and R declared
`receptor_expression: {"BtuB": 0.0}` per §4a. Each animal is inoculated with a
*single* strain, as in the paper — one C, one S, one R mouse per cage — so
displacement can only follow inter-host transfer. That makes the between-host
exchange rate (coprophagy; the model's `immigration` path already injects a
chosen `initial_strains` index) a structural parameter of the design and not a
free knob: it must be fixed and reported before the sweep, or it will be tuned
against the transition count along with everything else.

**Replication.** ≥5 seeds per arm. Justification: the λ arms of the carbon
ladder showed within-λ peak-N spread of 1.09 at 100 founders and a single-seed
monotonic result that did not survive replication; a displacement statistic is a
*ratio* of two stochastic populations and is noisier than either. Report the
within-arm band alongside every between-arm difference — the gate is that
between-prior separation exceeds within-prior spread, and if it does not, the
campaign has no resolving power and no prior is selected.

**Statistic.** Categorical, because the data are (§4d). Each simulated animal
is sampled on the paper's cadence — half-weekly for 12 weeks — through the
stool observation model: one pellet, subdivided, plated, with any strain below
1% of that pellet's recovered population recorded as *not detected*. Each
sample yields a dominant-strain call, so a run produces the same object the
paper publishes: a mouse × half-week grid of dominance states. Five numbers per
replicate, all computed on that grid and never on the underlying uncensored
counts:

1. **Transition count** — dominance changes summed over animals. This is the
   primary target: 123 for E1, 43 for E2.
2. **Time to C elimination** — fraction of cages with no C by week 7 (in vivo:
   10 of 12), and whether any cage fixes on C (in vivo: none).
3. **Directional asymmetry** — the transition matrix, tested with the same
   one-tailed Wald test the paper used, against the same three null hypotheses.
4. **Mixed-pellet frequency** — the paper's pellets were rarely mixed. A model
   that produces routinely mixed pellets has the wrong within-host competition,
   whatever its transition count.
5. **Transition sharpness** — the dwell time in a mixed or ambiguous dominance
   state, in half-week samples. In vivo the switches are sharp: a mouse flips
   from one dominant strain to another rather than drifting through a mixed
   interval. This is (4) read along the time axis, and it is a separate
   falsifier: an arm can hit the transition count with slow ratio drift that
   happens to cross the 50% line the right number of times, and that arm is
   wrong.

**Selection rule.** A prior is *admissible* if its replicate band contains the
in-vivo value of (1) and does not contradict (2)–(5). Report the admissible
**set**, not a point estimate; with five priors and this much stochasticity the
honest answer is likely an interval.

**Two-sided, and that is what makes it selective.** A displacement rate alone
would be satisfied by an arbitrarily strong producer. The paper's C strain is
also the *loser*: it fixed in no cage and was gone from 10 of 12 by week 7. A
lysis prior high enough to displace S quickly also kills producers faster than
they divide, so (1) and (2) pull the coefficient in opposite directions and
bound it from both sides. Any arm that reproduces the transition count while
leaving C established at week 7 has failed, not passed.

Two informative failure modes, both of which must be reported rather than
tuned away: if *every* prior is admissible the observable does not constrain
lysis; if *none* is, the release/transport parameters of §4b (or the missing
routes of §6) are wrong, and no value of the lysis coefficient will fix that.

**Falsifier.** The null arm must show no displacement. If the sensitive strain
declines at `sos_lysis_prob = 0` — through carbon competition, washout
asymmetry, or a cost difference — then the observable is not measuring colicin
antagonism and every selection made with it is void. This is the population-scale
asserted contrast the delivery campaign learned to demand: per-agent kill tests
passing does not establish that the population-scale observable responds.

**Forbidden.** Do not tune `retardation`, `burst_size`, resistance cost,
inter-host transfer rate, or inoculum against the same statistic
simultaneously with the prior. Four
parameters against one displacement curve is degenerate, and the model has
already shown that a free coefficient times an unaudited quantity looks like
biology.

## 6. What a selected prior does and does not mean

It is a model-selection output conditional on the architecture, not a measured
gut lysis rate — the work order's framing, and the code supports it: SOS
induction here has exactly four routes, and the in-vivo cell has more. Bile-salt
induction, nutrient-stress induction, and SOS-independent operon expression are
not represented. A prior selected against in-vivo displacement therefore
*absorbs* those routes: it is the per-generation lysis the model needs to
reproduce the observed kinetics through its one respiratory route, and it should
be expected to exceed a true measured respiratory lysis rate. Report it as an
upper bound on the respiratory contribution, with the unmodelled routes named.

## 7. Ordering

1. Fix the microcin penalty compounding defect (§8) — cheap, and it corrupts any
   continuous-release comparison arm.
2. ~~Settle the comparison observable~~ — done (§4d): per-animal stool,
   categorical dominance, 1% detection limit. This forces Layer 3.
3. Run the btuB-null probe (§4a-probe) to establish at population scale that
   the resistant genotype is both effective and costly, and to locate the
   corrinoid level at which the S-beats-R leg exists. Cheap, local, serial, and
   it gates the RPS claim.
4. Make `retardation` and `burst_size` configurable per BI locus, sourced, with
   parser + config-ingestion coverage (§4b).
5. Layer 2 patch network with reseeding, CPU-only, lookup patches, per the Phase
   1 plan and gates in `SPEC13_IMPLEMENTATION_REVIEW.md` (§4c).
6. Layer 3 luminal compartment and the stool observation model — animal, cage,
   pellet, subdivision, detection limit — per `SPEC13_MULTISCALE.md`. Without
   it there is no observable to select against.
7. Run the 12-arm × 5-seed selection with the null-arm falsifier (§5), locally
   and serial; no AWS or GPU without explicit authorization.
8. Publish the admissible set, the implied `k_ROS_funded` range with its
   uncrowded-`Q` validity caveat (§3), and the absorbed-routes caveat (§6).

## 8. Defect found while auditing the release chain

`Agent::flags.microcin_penalty_applied` is reset every step in
`simulation.cpp` alongside `just_divided`, but it is the only guard in
`FixBacteriocin::apply_microcin_secretion()` before
`agent.mu_max *= (1 - microcin_mu_penalty)`. Nothing restores `mu_max`, so the
documented "static 2–5% penalty" is re-charged every biology step: at the
shipped `0.03` a microcin producer retains `0.97^720 ≈ 3e-10` of its `mu_max`
over a 12 h run. `tests/test_bacteriocin.cpp` asserted a single application,
which is why it passed. No shipped example carries MccV, so no shipped result is
affected — but any Spec 13 arm comparing a continuous-release producer against a
lysis-release producer would have been comparing against a dead strain.
