# Devin Task: Rework Bacteriocin Induction Model + Add Phage-Mediated Lysis

## Context

GutIBM is a 3D agent-based model of E. coli bacteriocin ecology in the gut.
The current `fix_bacteriocin.cpp` has a single SOS induction pathway with a
constant basal rate (`sos_basal_rate = 1e-6 /s`). This is biologically
incomplete: it misses per-division induction, toxin-triggered cross-induction,
and the entirely separate phage-mediated lysis pathway used by Group B colicins.

The config parameter `sos_lysis_prob` (per-division SOS probability) is parsed
and stored in `BacteriocinConfig` but never used — this is a bug.

## What Needs to Change

### 1. Rework SOS Induction in `fix_bacteriocin.cpp`

Replace the current single-rate `check_sos_induction()` with a multi-trigger model.
The SOS induction probability per timestep should be:

```
p_sos = 1 - exp(-(rate_basal + rate_replication + rate_toxin) * dt)
```

Where:

**a) Basal rate** (existing, keep as-is):
```
rate_basal = cfg_.sos_basal_rate   // default 1e-6 /s
```

**b) Replication-coupled rate** (NEW — uses the currently unused `sos_lysis_prob`):
```
rate_replication = agent.just_divided ? cfg_.sos_lysis_prob / cfg_.bio_dt : 0
```

This requires a `bool just_divided` flag on Agent, set by `fix_metabolism.cpp`
during division and cleared at the start of each timestep. The flag converts
the per-division probability into a per-timestep rate only during the step
in which division occurred.

Biological basis: >75% of enteric bacteriocin promoters are SOS-regulated via
LexA cleavage. Replication stress during division is a primary trigger. The
knowledge base reports ~3% of colicin-producing cells undergo lysis to release
ColE2.

**c) Toxin-induced cross-induction rate** (NEW):
```
rate_toxin = k_cross * [DNA_damaging_toxin_at_agent_position]
```

Where `k_cross` is a new config parameter (`sos_cross_induction_rate`, default
1e3 /s per mol/m³). This applies ONLY when the local toxin concentration at
the agent's grid cell includes DNA-damaging colicins (nuclease types: ColE2,
E7, E8, E9). Pore-forming colicins (ColE1, ColB, ColIa) do NOT trigger SOS.

Implementation: Add a `bool is_nuclease` field to `BICluster`. Set it true for
ColE2 (toxin_id 2). In `check_sos_induction()`, read the local bacteriocin
concentration from the chemical field and multiply by `k_cross` only if the
dominant local toxin is nuclease-type.

Biological basis: The "provoker" mechanism — sublethal doses of nuclease
colicins (E2, E7) trigger SOS in neighboring cells, inducing them to lyse and
release their own colicins. This creates autocatalytic cascading.

### 2. Add Phage-Mediated Lysis Pathway

Group B colicins (ColB, ColIa, and any on large conjugative plasmids) are
released via temperate phage-mediated bacteriolysis, NOT SOS-mediated suicide.
This is a fundamentally different pathway:

- The producing cell carries a lysogenic prophage
- Spontaneous prophage induction occurs at a low rate (~1e-4 to 1e-3 per cell
  per generation)
- Upon induction, the phage lytic cycle produces virions AND releases colicin
- The phage can also lysogenize nearby susceptible cells (spreading the colicin
  gene)

**Implementation:**

Add new fields to `BICluster`:
```cpp
enum class ReleaseMode : uint8_t {
    SOS_LYSIS = 0,      // Group A: SOS-triggered suicide lysis
    PHAGE_LYSIS = 1,    // Group B: temperate phage-mediated
    CONTINUOUS = 2       // Microcins: secreted without lysis
};

// In BICluster:
ReleaseMode release_mode = ReleaseMode::SOS_LYSIS;
Real phage_induction_rate = 0.0;  // only for PHAGE_LYSIS mode
Real phage_burst_size = 0;        // virions produced per lysis
Real phage_lysogeny_rate = 0.0;   // per virion per susceptible per s
```

Update the plasmid library:
- ColE1, ColE2, ColM: `release_mode = SOS_LYSIS`
- ColB, ColIa: `release_mode = PHAGE_LYSIS`, `phage_induction_rate = 1e-4 /generation`
- MccV: `release_mode = CONTINUOUS`

In `fix_bacteriocin.cpp`, the `compute()` loop should branch on release mode:

```cpp
for (Agent& a : agents) {
    for (auto& bi : a.genome.bi_loci) {
        switch (bi.release_mode) {
        case ReleaseMode::SOS_LYSIS:
            check_sos_induction(a, dt);  // existing (reworked per §1)
            break;
        case ReleaseMode::PHAGE_LYSIS:
            check_phage_induction(a, bi, dt);  // NEW
            break;
        case ReleaseMode::CONTINUOUS:
            // Already handled: microcin secretion applies mu_penalty
            break;
        }
    }
}
```

**`check_phage_induction()`:**
```cpp
void FixBacteriocin::check_phage_induction(Agent& agent, BICluster& bi, Real dt) {
    if (agent.state == PhenoState::SOS_INDUCED) return;

    // Convert per-generation rate to per-second
    Real gen_time = (agent.mu_realized > 0) ? 
        std::log(2.0) / agent.mu_realized : 1e6;
    Real rate_per_s = bi.phage_induction_rate / gen_time;

    Real p_induction = 1.0 - std::exp(-rate_per_s * dt);

    if (sim_.rng().bernoulli(p_induction)) {
        agent.state = PhenoState::SOS_INDUCED;  // reuse state for phage lysis
        agent.sos_timer = 60.0;  // phage lytic cycle ~1 hr, but late genes ~1 min
    }
}
```

**Phage lysogeny (optional, lower priority):**
When a phage-induced cell lyses, in addition to the colicin burst, it produces
`phage_burst_size` virions. These diffuse locally (use the existing toxin
diffusion machinery with a separate chemical species or just immediate
neighbor check). Susceptible cells within contact range have a probability
`phage_lysogeny_rate * dt` of acquiring the BI cluster — effectively another
HGT pathway alongside conjugation.

### 3. Per-Colicin Burst Size

Currently `burst_molecules = 1e4` is a single constant. It should be per-colicin:

- ColE1 plasmid: ~20 copies, high expression from SOS promoter → ~1e5 molecules
- ColE2: ~20 copies, but secreted as complex → ~5e4 molecules  
- ColB: monocopy plasmid, phage-mediated → ~1e4 molecules
- ColIa: monocopy, phage-mediated → ~1e4 molecules
- ColM: small protein (29.5 kDa), high copy → ~2e5 molecules
- MccV: continuous secretion, not burst (already handled separately)

Add a `burst_size` field to `BICluster` and set it in the plasmid library.
Remove the global `burst_molecules` from `BacteriocinConfig` (or keep as
fallback default).

### 4. Agent Flag for Division

Add to the Agent struct:
```cpp
bool just_divided = false;
```

In `fix_metabolism.cpp`, set `a.just_divided = true` and `daughter.just_divided = true`
during division. At the top of `Simulation::step()`, clear all `just_divided` flags.

### 5. New Config Parameters

Add to `BacteriocinConfig`:
```cpp
Real sos_cross_induction_rate = 1e3;  // 1/s per mol/m³ of nuclease toxin
```

Add to `InputParser`:
```cpp
if (key == "sos_cross_induction_rate") { cfg.bacteriocin.sos_cross_induction_rate = parse_config_real(key, val); return true; }
```

### 6. Tests to Add/Modify

**Modify `test_bacteriocin.cpp`:**
- `test_sos_induction_on_division`: Create agent with BI locus, set 
  `just_divided = true`, high `sos_lysis_prob`. Verify SOS is triggered.
- `test_no_sos_without_division`: `just_divided = false`, only basal rate.
  Verify SOS probability matches basal expectation.
- `test_phage_induction`: Agent with ColB (PHAGE_LYSIS mode). Verify phage
  induction occurs at correct rate.
- `test_phage_does_not_trigger_sos_path`: Agent with ColB should not use
  the SOS induction pathway.
- `test_microcin_no_lysis`: Agent with MccV should never enter SOS state.
- `test_cross_induction`: Place nuclease toxin in agent's grid cell. Verify
  elevated SOS rate.

**Modify `test_config_diversity.cpp`:**
- Fix the assertion failure at line 215 (`test_parsed_fix_list_is_honored`).
  This is an existing bug where the Fix subset parsing doesn't produce
  different fingerprints. Investigate whether the `fixes` array is actually
  being honored in `Simulation::init()`.

## Files to Modify

1. `src/core/types.h` — Add `ReleaseMode` enum, add `release_mode`, 
   `phage_induction_rate`, `phage_burst_size`, `phage_lysogeny_rate`, 
   `is_nuclease`, `burst_size` to BICluster. Add `just_divided` to Agent.
2. `src/fixes/fix_bacteriocin.h` — Add `sos_cross_induction_rate` to config.
   Add `check_phage_induction()` declaration.
3. `src/fixes/fix_bacteriocin.cpp` — Rework `check_sos_induction()` per §1.
   Add `check_phage_induction()` per §2. Branch `compute()` on release_mode.
   Use per-colicin burst_size in `lyse_agent()`.
4. `src/fixes/fix_metabolism.cpp` — Set `just_divided = true` during division.
5. `src/core/simulation.cpp` — Clear `just_divided` flags at top of `step()`.
6. `src/genome/plasmid.cpp` — Set per-colicin `release_mode`, `burst_size`,
   `is_nuclease`, and phage parameters.
7. `src/io/input_parser.cpp` — Parse `sos_cross_induction_rate`.
8. `tests/test_bacteriocin.cpp` — Add/modify tests per §6.
9. `tests/test_config_diversity.cpp` — Fix line 215 assertion.
10. `docs/PARAMETERS.md` — Document new parameters.
11. `docs/MECHANISMS.md` — Update fix_bacteriocin section.

## Biological References

- >75% of enteric bacteriocin promoters are SOS-regulated (LexA cleavage)
- ~3% of colicin-producing cells undergo lysis to release ColE2 (Götz 2020)
- Group B colicin release relies entirely on temperate phage-mediated lysis
  (knowledge_base_report.md §Regulatory Triggers)
- ColE1 plasmid: ~20 copies/cell; Group B: monocopy (Schmidt & Inselburg 1982)
- Nuclease colicins (E2, E7) trigger SOS in neighboring cells via DNA damage
  ("provoker" mechanism, knowledge_base_report.md §Bacteriocin Mediated Interference)
- Prophage lysis genes required for ColIb release in S. Typhimurium
- Temperate bacteriophages = 20-50% of total gut phage population

## Priority

1. §1 (SOS rework) + §4 (just_divided flag) — highest, fixes the bug
2. §3 (per-colicin burst size) — straightforward
3. §2 (phage-mediated lysis) — important for Group B biology
4. §6 (tests) — validate everything
5. Phage lysogeny (§2 optional) — lower priority, can be a follow-up
