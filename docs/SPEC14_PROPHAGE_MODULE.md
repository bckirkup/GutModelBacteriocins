# Spec 14 Temperate Prophage Module (`GAP-10`, Issue #475)

## 1. Purpose and Scientific Status

This document specifies the temperate prophage dynamics module (`FixProphage`, registered as `"prophage"` in `FixRegistry`) implemented in GutIBM following the resolution of the `GAP-01` BtuB competitive-shielding prerequisite (`#464` / PR `#488`).

> **Scientific caveat (*in vivo* lysogen displacement):** Consistent with `Henrot et al. (2026 preprint)` (`10.64898/2026.01.07.698133`), which reported null *in vivo* bacteriocin enhancement of temperate prophage transfer in the murine gut despite positive *in vitro* stimulation, *in vivo* bacteriocin-driven lysogen displacement remains **hypothesis / regime-dependent** and is **not** treated as a demonstrated default gut outcome. The module ships **default-off** (`prophage.enabled = false`) so existing trajectories remain bitwise unchanged unless explicitly enabled.

## 2. Separate Prophage Carriage vs. Active Inducibility

Following `docs/SPEC14_PRIOR_REVIEW.md` §3 and the repository calibration guardrail (*a sequence census is an upper bound on a functionally defined parameter*), genomic prophage carriage is decoupled from active excision/lysis capability:

- **Genomic carriage (`lysogen_prevalence_init` in `[0, 1]`):** fraction of founders carrying integrated prophage sequence (e.g., Dikareva et al. 2023 upper-bound prior `0.80` for human fecal *Enterobacteriaceae*). Every carrier receives `agent.flags.has_prophage_immunity = true`, conferring repressor-mediated superinfection immunity even when the prophage is cryptic or defective.
- **Active inducible subset (`active_inducible_fraction` in `[0, 1]`):** conditional fraction of carriers whose prophage can excise, lyse the host, and produce free virions. These cells additionally receive `agent.flags.has_inducible_prophage = true` (effective population prevalence `lysogen_prevalence_init * active_inducible_fraction`).
- **Cryptic carriers (`has_prophage_immunity == true` and `has_inducible_prophage == false`):** never undergo spontaneous or DNA-damage induction, yet still adsorb and neutralize free virions with zero lytic kill and zero state transition.
- **Per-strain founder overrides:** `initial_strains[].lysogen_prevalence` and `initial_strains[].active_inducible_fraction` override the global defaults per founder population and propagate across division inheritance, `agent_transfer.cpp` MPI serialization, Layer-2/3 `founder_genotype.cpp`, HDF5 checkpoints (backward-compatible `false` defaults when absent), and `AgentPoolGpu` SoA buffers.

## 3. DNA-Damage Hill Induction and Dual-Regime Presets

Each live, non-SOS-induced, actively inducible lysogen experiences the hazard (1/s):

$$h_{\mathrm{ind}}(D) = k_{\mathrm{spont}} + k_{\max,\mathrm{damage}} \cdot \frac{D^{n}}{K_{\mathrm{damage}}^{n} + D^{n}}$$

with local DNA-damage driver signal:

$$D = D_{\mathrm{stimulus}} + w_{\mathrm{nuc}} \cdot C_{\mathrm{nuclease,eff}} + w_{\mathrm{ROS}} \cdot r_{\mathrm{ROS}}$$

where:

- $C_{\mathrm{nuclease,eff}}$ sums receptor-expression- and cognate-immunity-weighted local nuclease colicin exposure (e.g., `ColE2` DNase via `QSSASolver::sampled_nuclease_conc`), protecting cognate `ColE2` producers via `immunity_factor` and gating entry by receptor expression (`BtuB`).
- $r_{\mathrm{ROS}}$ reads `Simulation::ros_induction_rate(agent)` (funded respiratory or ambient ROS when oxygen/ROS is active).
- $D_{\mathrm{stimulus}}$ (`prophage.damage_stimulus`) provides an explicit exogenous SOS trigger (e.g., mitomycin-C-like stress or deterministic Hill curve probes).
- Agents already in `PhenoState::SOS_INDUCED` (committed to `FixBacteriocin` SOS lysis) are skipped by `FixProphage` to avoid double-counting colicin SOS lysis.

### Empirical Hazard Presets (`examples/prophage_presets/`)

| Preset | `lysogen_prevalence_init` | `active_inducible_fraction` | Effective active fraction | `k_spont` (1/s) | Implied baseline rate (/bacterium/day) | Empirical reference band |
|---|---|---|---|---|---|---|
| `human_commensal_community.json` | `0.8` | `0.375` | `0.30` | `1.2e-7` | `~0.0031` | Lopez et al. (2025) / Current Biology 2025 human-gut community band (`0.001`–`0.01` /bacterium/day) |
| `murine_monoxenic_sos.json` | `1.0` | `1.0` | `1.00` | `2.0e-6` | `~0.173` | De Paepe et al. (2016) / Frazao et al. (2022) monoxenic mouse *E. coli* (`1.70e-6`–`4.44e-6` 1/s) |

Realized community daily induction rate (`prophage_induction_rate_per_bacterium_day = (prophage_inductions / prophage_bacterium_seconds) * 86400`) is emitted alongside interval and cumulative `StepEvents` counters (`prophage_inductions`, `prophage_lytic_kills`, `prophage_lysogen_conversions`, `virions_released`, `virions_adsorbed`, `virions_decayed`) in both in-memory summaries and HDF5 `summary/*/events/`.

## 4. Virion Transport, Adsorption, Superinfection Immunity, and Conservation Ledger

When `prophage.enabled = true`, `InputParser::finalize_config()` registers `species::VIRION` (`"virion"`, stored in virions/m³) with:

- Implicit L-stable backward-Euler 3D diffusion (`D_virion`, default `4.0e-12` m²/s for a ~50–100 nm phage particle in mucus; optional additional `virion_retardation >= 1.0`), periodic in `x`/`y` and strictly zero-flux (`EpithelialBoundaryMode::Flux`, `epithelial_flux = 0.0`) at both `z = 0` and `z = H`.
- Exact per-voxel competing-rate removal over timestep `dt`:
  - Total clearance rate in voxel $c$: $\Lambda_c = K_c + A_c$, where $K_c = k_{\mathrm{virion\_decay}} + s_{\mathrm{washout}} \cdot \gamma_{\mathrm{flow}}(z_c)$ couples first-order inactivation with `AdvectionField::washout_rate(z)`, and $A_c = \sum_{i \in c} a_i$ with $a_i = k_{\mathrm{adsorb}} r_i / \Delta V$.
  - Removed virions in voxel $c$: $V_{\mathrm{rem}}(c) = V_0(c)\,(1 - e^{-\Lambda_c \Delta t}) \le V_0(c)$, partitioned analytically into $V_{\mathrm{decayed}}(c) = V_{\mathrm{rem}}(c) K_c / \Lambda_c$ and per-agent adsorbed virions $v_i = V_{\mathrm{rem}}(c) a_i / \Lambda_c$.
  - Because $V_{\mathrm{rem}}(c) \le V_0(c)$ identically, free virion concentrations remain strictly non-negative with **zero reaction clipping**.
- Adsorption outcomes for live owned cells:
  - **Existing lysogens (`has_prophage_immunity == true`, active or cryptic):** adsorbed virions $v_i$ are neutralized and booked into `virions_adsorbed` with zero lytic kill and zero state change.
  - **Susceptible non-lysogens (`has_prophage_immunity == false`):** infected with Poisson hit probability $1 - e^{-v_i}$; upon infection, convert into an active lysogen (`has_prophage_immunity = true`, `has_inducible_prophage = true`, `prophage_lysogen_conversions++`) with probability `p_lysogeny` (default `0.19`, De Paepe et al. 2016), or undergo lytic death (`prophage_lytic_kills++`, releasing `virion_burst_size`, default `12.1`) with probability `1 - p_lysogeny`.
- **Exact conservation closure (CPU and GPU):**
  - Virions: $V_{\mathrm{grid}}(T) + V_{\mathrm{adsorbed,cum}}(T) + V_{\mathrm{decayed,cum}}(T) = V_{\mathrm{grid}}(0) + V_{\mathrm{released,cum}}(T)$ to double-precision rounding (`< 1e-9` relative).
  - Cells (closed patch): $N_{\mathrm{live}}(T) = N_{\mathrm{live}}(0) + N_{\mathrm{divisions}} - N_{\mathrm{mortality\_lysis}} - N_{\mathrm{prophage\_lytic\_kills}}$.
- **GPU parity:** `src/gpu/prophage_gpu.{h,cpp}` and `src/gpu/prophage_kernel.cu` evaluate Hill induction probabilities, per-cell adsorption accumulation, competing-rate decay/adsorption splitting, and susceptible infection probabilities on CUDA devices using `AgentPoolGpu` SoA prophage buffers, while preserving host-serial `RNG` draws for exact CPU/GPU parity (`tests/test_gpu_prophage_parity.cpp`).

## 5. Configuration Reference

| Key | Type | Default | Description |
|---|---|---|---|
| `prophage.enabled` | bool | `false` | Master switch; registers `species::VIRION` when `true` |
| `prophage.lysogen_prevalence_init` | real `[0, 1]` | `0.8` | Genomic prophage carriage fraction (`has_prophage_immunity`) |
| `prophage.active_inducible_fraction` | real `[0, 1]` | `0.375` | Active inducible subset among carriers (`has_inducible_prophage`) |
| `prophage.k_spont` | real `>= 0` | `1.2e-7` | Spontaneous induction hazard per active lysogen (1/s) |
| `prophage.k_max_damage` | real `>= 0` | `1.0e-3` | Max DNA-damage Hill induction hazard (1/s) |
| `prophage.K_damage` | real `> 0` | `1.0e-6` | Hill half-saturation DNA-damage signal |
| `prophage.hill_n` | real `> 0` | `2.0` | Hill cooperativity exponent |
| `prophage.damage_stimulus` | real `>= 0` | `0.0` | Exogenous DNA-damage trigger added to signal `D` |
| `prophage.nuclease_damage_weight` | real `>= 0` | `1.0` | Weight on receptor- and immunity-scaled local nuclease colicin |
| `prophage.ros_damage_weight` | real `>= 0` | `1.0` | Weight on local ROS induction signal |
| `prophage.lysis_delay_s` | real `>= 0` | `0.0` | Latent lysis timer after induction (`0.0` = immediate lysis) |
| `prophage.virion_burst_size` | real `>= 0` | `12.1` | Free virions released per prophage induction or lytic kill |
| `prophage.D_virion` | real `>= 0` | `4.0e-12` | Free/mucus virion diffusivity (m²/s) |
| `prophage.virion_retardation` | real `>= 1` | `1.0` | Additional mucin retardation factor (`D_eff = D_virion / virion_retardation`) |
| `prophage.k_virion_decay` | real `>= 0` | `1.0e-4` | First-order free virion decay rate (1/s) |
| `prophage.virion_washout_scale` | real `>= 0` | `1.0` | Coupling scale to `AdvectionField::washout_rate(z)` |
| `prophage.k_adsorb` | real `>= 0` | `1.0e-15` | Receptor-mediated adsorption clearance constant (m³/s/cell) |
| `prophage.p_lysogeny` | real `[0, 1]` | `0.19` | Lysogenization probability on susceptible infection |
| `prophage.adsorption_receptor` | string | `"none"` | `"none"` (constitutive receptor `= 1`) or a `ReceptorType` name (e.g., `"FhuA"`, `"BtuB"`) |
| `initial_strains[].lysogen_prevalence` | optional real `[0, 1]` | unset | Per-strain founder override for `has_prophage_immunity` |
| `initial_strains[].active_inducible_fraction` | optional real `[0, 1]` | unset | Per-strain founder override for active subset of carriers |

## 6. Validation Coverage (`tests/test_prophage.cpp`, `tests/test_gpu_prophage_parity.cpp`)

1. **Default-off backward parity:** bitwise identical simulation fingerprint and zero `species::VIRION` registration when `prophage.enabled = false`.
2. **Cryptic vs. active lysogen separation:** cryptic carriers (`has_prophage_immunity = true`, `has_inducible_prophage = false`) suffer zero lysis under `D = 50 * K_damage` and remain 100% superinfection-immune under dense virion exposure while still adsorbing virions.
3. **Hill monotonicity:** exact analytic and live stochastic sweep across `D in {0, 0.5 * K_damage, K_damage, 5 * K_damage}`.
4. **Cell + virion conservation closure and inheritance:** exact closed-patch cell balance, `< 1e-9` relative virion mass/count closure with zero reaction clipping, latent lysis timer progression, and division inheritance of cryptic/active flags.
5. **Henrot et al. (2026) positive *in vitro* vs. suppressed *in vivo* contrast:** paired well-mixed *in vitro* control (`high D_virion`, zero mucus washout, low cryptic blockade, unshielded `ColE2` vs. plasmid-free null) produces measurable `ColE2`-provoked lysogen conversions, whereas mucosal retardation + washout/decay + high background cryptic immunity suppresses transfer by more than 5×.
6. **Persistence and capability guardrails:** MPI `agent_transfer` pack/unpack, HDF5 checkpoint round-trip + legacy backward-compatible defaults, exhaustive `test_config_ingestion.cpp` tracking, and `mechanism_capability.cpp` refusal on tabulated-only Layer-2/3 configurations.
