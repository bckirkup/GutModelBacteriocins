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
between lysis and the death of a sensitive cell. Four are known and open.

**(a) There is no resistant strain.** `experiments/rps_campaign/rps_spec_audit.md`
§C: `SimulationConfig::InitialStrain` has no receptor-genotype field and
`parse_strain_object` ignores unknown keys, so `receptor_mutations` parses
silently and does nothing; every founder gets `receptor_expr_base.fill(1.0)`.
Resistance arises only from `FixMutation` at 1e-7/step. **Kirkup & Riley 2004 is
a three-strain non-transitive result; the model can currently express only C vs
S.** A C-vs-S displacement selection is still well posed (and is the pairwise
leg the paper reports as resolved), but it cannot be presented as reproducing
the RPS result, and the R-strain gap is a prerequisite for the coexistence
claim. Note also the audit's double-counting warning: a receptor knockout
already costs growth through `Km_b12 = km_b12/(expr·affinity)` (9% at 1 µM
corrinoid, 25% at 0.3 µM), so an *additional* `mu_max` cost charges resistance
twice.

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

**(d) The comparison observable is not yet decided, and it decides the layer.**
Kirkup & Riley recovered strain frequencies from animals over weeks. If the
target series is *fecal/stool* counts, the model has no such observable until
the luminal compartment exists (Layer 3), because the patch's only export
channel books `outflow_boundary` and discards the cell. If a *mucosal* ratio is
the fair comparison, Layer 2 suffices. This is a question about the source data,
not about the code, and it is the one input I cannot derive: **the digitized
time series, its sampling cadence, replicate count, detection limit, and
inoculum ratio have to come from the paper's figures before the selection
statistic can be defined.**

## 5. The selection procedure

Deliberately structured so that the answer can be "none of them", and so that a
null arm can falsify the observable itself.

**Arms.** `sos_lysis_prob ∈ {0, 0.005, 0.01, 0.02, 0.05, 0.10}` (per generation,
with basal and ROS routes zeroed) crossed with `retardation ∈ {50, 1.2}` = 12
arms. The 0 arm is the null: with release lysis-gated, it emits no toxin at all.
The 0.10 arm brackets above the in-vitro band so that a monotone-but-unbounded
response is visible as such.

**Replication.** ≥5 seeds per arm. Justification: the λ arms of the carbon
ladder showed within-λ peak-N spread of 1.09 at 100 founders and a single-seed
monotonic result that did not survive replication; a displacement statistic is a
*ratio* of two stochastic populations and is noisier than either. Report the
within-arm band alongside every between-arm difference — the gate is that
between-prior separation exceeds within-prior spread, and if it does not, the
campaign has no resolving power and no prior is selected.

**Statistic.** Computed on whatever quantity the source data reports (§4d), with
the model censored to match: the in-vivo sampling cadence, the detection limit,
and the same inoculum ratio. Two numbers per replicate — displacement half-time
(time for the sensitive fraction to fall to half its initial value) and the
strain ratio at the final in-vivo sampling time.

**Selection rule.** A prior is *admissible* if its replicate band contains both
in-vivo statistics. Report the admissible **set**, not a point estimate; with
five priors and this much stochasticity the honest answer is likely an interval.
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

**Forbidden.** Do not tune `retardation`, `burst_size`, resistance cost, or
inoculum against the same time series simultaneously with the prior. Four
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
2. Obtain the source time series and settle the comparison observable (§4d).
   This determines whether Layer 2 or Layer 3 is required and is a blocking
   input, not an implementation task.
3. Make `retardation` and `burst_size` configurable per BI locus, sourced, with
   parser + config-ingestion coverage (§4b).
4. Add the configurable resistant strain — receptor genotype on
   `InitialStrain`, cost taken mechanistically and reported, not double-charged
   (§4a). Required for the RPS claim, not for the C-vs-S selection.
5. Layer 2 patch network with reseeding, CPU-only, lookup patches, per the Phase
   1 plan and gates in `SPEC13_IMPLEMENTATION_REVIEW.md` (§4c).
6. Run the 12-arm × 5-seed selection with the null-arm falsifier (§5), locally
   and serial; no AWS or GPU without explicit authorization.
7. Publish the admissible set, the implied `k_ROS_funded` range with its
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
