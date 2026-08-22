# RPS campaign spec: what the code can ingest today, and what the numbers say

Audit of `rps_campaign_config_spec.md` against the tree at current `main`.
Every claim below is a code reference or a computed value, not an assertion
about intent.

## Unit convention (used throughout)

`mol/m^3 = mmol/L`, so **1 nM = 1e-6 mol/m^3** and **1 µM = 1e-3 mol/m^3**.
The defaults follow this: B12 field 1e-3 mol/m^3 = 1 µM (matching the 0.97 µM
fecal figure), `kd_b12_btuB` 1e-6 = 1 nM.

## A. Three translation errors in the proposed JSON

1. **Kd values are 1000x too tight.** The literature table gives BtuB–colicin E
   Kd = 1 nM and ColB–FepA = 2 nM, which are `1e-6` and `2e-6` mol/m^3. The JSON
   writes `kd_colicinE_btuB: 1e-09` and `kd_colicinB_fepA: 2e-09`, i.e. 1 pM and
   2 pM. Note the current defaults (5e-7, 2e-6) are already at the literature
   values, so no change is actually indicated here.
2. **The competition factor doesn't come out at 11.** The factor is
   `1 + [corrinoid]/kd_b12_btuB`. The JSON sets `kd_b12_btuB = 1e-4` *and*
   reduces corrinoid to 3e-4, which gives `1 + 3 = 4`, not 11. The 11 came from
   my sweep at corrinoid 1e-3. The two changes were derived independently and
   compound.
3. **Literature Kd + literature corrinoid lands back in the dead arm.** Taking
   the file's own BtuB–cobalamin Kd of 0.3 nM (3e-7) with 0.3 µM mucus corrinoid
   gives a factor of 1001 — exactly the arm where I measured 0.66 expected kills
   and zero realized. The escape is the file's own analogue argument: 98% of the
   pool is non-Cbl analogue at Kd ~10 nM. Using the two-ligand form the file
   recommends,
   `1 + [Cbl]/Kd_Cbl + [analogue]/Kd_an = 1 + 10/0.3 + 290/10 ~= 63`.
   That is a genuinely literature-grounded competition factor, and it sits
   between my 101x and 11x arms — i.e. **colicin weak but non-zero**, roughly
   4–6 expected kills per 17 h rather than 0.7 or 16.
   The model has one competing ligand per receptor, so the two-species pool must
   be collapsed into an effective Kd: `[L]/Kd_eff = 62` at `[L] = 3e-4` gives
   **`kd_b12_btuB ~= 5e-6` (5 nM)** — which is a *different* number from both the
   current default and the spec's proposal, and is the one I'd actually run.

## B. Retardation: the spec is aimed at the wrong number

There are three retardation values in the tree, and the one the spec proposes to
change is the one that does nothing.

| Where | Value for ColE1 | Actually used? |
|---|---|---|
| `plasmid.cpp` `colicin_E1().retardation` | **50.0** | **Yes** — `fix_bacteriocin.cpp:152` copies it into the burst, and `greens_function.cpp:304` computes `D_eff = diff_coeff / retardation`. Hardcoded, no config key. |
| `chemicals` table, `BACTERIOCIN_*` | 10.0 | No — toxin species have `diffusion_enabled = false`; QSSA owns them. This is the "10" the spec cites. |
| `retardation_basic/neutral/acidic` (50/5/1.5) | 50.0 | **No** — their only consumer is `FixBacteriocin::retardation_for_pI`, which is **dead code, called from nowhere**. Sweeping these keys would produce identical runs. |

So the effective ColE1 transport is `D_eff = 4e-11/50 = 8e-13 m^2/s`, against the
literature analogue `7e-11/1.2 = 5.8e-11` — **72x slower**. Per 60 s step the rms
displacement is 9.8 µm versus 83 µm; kill-zone radius ~8.5x, area ~72x. On the
numbers this is a larger lever on colicin efficacy than the corrinoid knob, and
it is the one thing in the spec that **cannot be swept without a code change**.
`D_free` and `burst_size` are hardcoded per plasmid the same way.

## C. The Resistant strain cannot be configured at all

`SimulationConfig::InitialStrain` (`src/io/input_parser.h:118`) is
`{type, count, mu_max, plasmids, conjugative, cdi_type, cdi_immunity}`. There is
no receptor-genotype field, and `parse_strain_object` ignores unknown keys, so
`receptor_mutations: {BtuB: null, FepA: null}` and `label` parse silently and do
nothing. `receptor_expr_base.fill(1.0)` for every founder
(`src/core/agent.cpp:48`); receptor knockouts arise only stochastically from
`FixMutation` at 1e-7/step. **This is the blocking gap: RPS has no R strain.**

Two related consequences once it exists:

- **The cost is emergent, and the spec would double-count it.** BtuB expression
  feeds `Km_b12 = km_b12 / (expr * ligand_affinity)` with `expr` floored at 0.01
  (`fix_metabolism.cpp:397`, deliberately, to avoid a divide-by-zero). At the
  hardcoded `km_b12 = 1e-6` a BtuB-null founder gets `monod_b12` = 0.91 at 1 µM
  corrinoid (**~9% cost**) and 0.75 at 0.3 µM (**~25% cost**) — not 2%, and it
  moves with the same corrinoid knob that sets colicin competition, so the
  "resistance cost" and "corrinoid" factors in the core sweep are not
  independent. Setting a 2% cost *via `mu_max`* on top of a knocked-out receptor
  charges the cell twice. I'd take the cost mechanistically and report it, or set
  `mu_max` costs with receptors intact — not both.
- FepA-null similarly reduces the iron Monod term rather than costing a fixed 5%.

## D. Other keys, and whether they exist

| Spec key | Status |
|---|---|
| `kill_rate_colicin: 0.001` | Exists; unchanged from default. Independently supported: 1e-3/s at occupancy 1 gives mean time-to-death ~17 min against the file's 15 min (2–60) commitment-to-death. |
| `kd_b12_btuB`, `kd_colicinE_btuB`, `kd_colicinB_fepA`, `kd_enterobactin` | All real keys. |
| `corrinoid_field.concentration_mol_m3: 3e-4` | **No key.** B12 initial/boundary conc is hardcoded in `input_parser.cpp:125`. Only carbon, mucin and ferrichrome have concentration keys. Needs code. |
| `sos_lysis_prob`, `burst_release_tau`, `toxin_cutoff` | Real keys. |
| `burst_size_molecules`, `D_free_colicin_m2s`, `retardation_mucus` | **No keys** — hardcoded per plasmid (§B). Note `D_free_colicin` in the resolved config is a *different* setting from `colicin_E1().diff_coeff`, which is what the burst uses. |
| `delivery.mode: flux`, `epithelial_flux` | Real keys; this is the flux-0.3x arm. |
| `vbf.carbon_sink_vmax` | Real key. |
| `immigration.rate_per_day`, `types_introduced: [1,2,3]` | Partly: `immigration.{enabled,rate,count,strain_index,schedule,placement}` exist, but `strain_index` is a **single** type. Random-type immigration needs code. Also the units are per-step rate, not per-day. |
| Immunity cost 0.2% | No immunity-specific cost exists. `plasmid_copy_cost` = 0.02 is a 2% penalty per *transferred* plasmid via conjugation, not a founder cost. |
| Lysis fraction 0.2 of induced, 90 min induction-to-lysis | Not representable: induction leads to lysis with a fixed delay constant, no survive-and-repair branch. Would need code; today every induced cell lyses. |

## E. Two campaign-design blockers unrelated to parameters

1. **The 7-day horizon does not survive the dysbiosis guard.** In the flux-0.3x
   arm the guard halted 8 of my 12 corrinoid runs at 17–21 h, and the four that
   reached 24 h did so only because the guard's trajectory test didn't trip.
   Nothing in this configuration stays inside the envelope for 168 h; the only
   arm that never bloomed was flux-0.1x, which decays. A 7-day RPS run needs
   either that arm, an explicit carrying-capacity mechanism, or the guard off
   with the bloom accepted as a known artifact.
2. **205 runs is not an overnight job.** The queue runs strictly serially (4
   vCPU cap, one g4dn job takes all 4), and 24 h of simulated time cost ~200–275
   s wall. A 7-day horizon is 7x the steps at larger populations, so ≥30–60
   min/run, i.e. **4–8 days of continuous queue time** for 205 runs. Raising the
   compute-environment vCPU cap needs `batch:UpdateComputeEnvironment`, which the
   Devin role does not have.

## F. What I'd do, in order

1. Code: per-strain receptor genotype in `InitialStrain` (unblocks R), a
   corrinoid concentration key, per-plasmid transport/burst overrides, and either
   wire or delete the dead `retardation_*` keys. Small, independent PRs.
2. Then a **retardation** arm before the full matrix — on these numbers it
   dominates, and it has never been varied.
3. Run the corrinoid axis at the analogue-corrected `kd_b12_btuB ~= 5e-6` with
   corrinoid 3e-4, not at either the current default or the spec's 1e-4.
4. Report the R-strain cost as measured rather than imposing 2% on top of it.
