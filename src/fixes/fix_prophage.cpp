/* -----------------------------------------------------------------------
   GutIBM – Temperate prophage Fix implementation (Spec 14, GAP-10)
   ----------------------------------------------------------------------- */

#include "fix_prophage.h"
#include "prophage_gpu.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

namespace gutibm {

namespace {

Int map_agent_cell_to_storage(const ChemicalField& chem, const Domain& domain,
                              Int global_cell) {
  if (global_cell < 0) return -1;
  if (chem.slab_mode()) {
    return chem.global_to_storage_cell(global_cell);
  }
  const Int global_plane = domain.nx() * domain.ny();
  const Int rem = global_cell % global_plane;
  const Int ix = rem % domain.nx();
  if (ix < domain.local_grid_x_begin() || ix >= domain.local_grid_x_end()) {
    return -1;
  }
  return global_cell < chem.ncells() ? global_cell : -1;
}

Real agent_adsorption_expression(const Agent& agent, Int receptor_idx) {
  if (receptor_idx < 0 || receptor_idx >= NUM_RECEPTORS) return 1.0;
  return std::max(agent.receptor_expr[static_cast<size_t>(receptor_idx)], 0.0);
}

bool agent_participates_in_adsorption(const Agent& agent) {
  return !agent.flags.is_ghost
      && agent.state != PhenoState::DEAD
      && agent.state != PhenoState::SOS_INDUCED;
}

bool agent_eligible_for_induction(const Agent& agent) {
  return !agent.flags.is_ghost
      && agent.state != PhenoState::DEAD
      && agent.state != PhenoState::SOS_INDUCED
      && agent.flags.has_inducible_prophage
      && !agent.flags.prophage_induced;
}

}  // namespace

Real prophage_hill_induction_hazard(Real damage_signal,
                                    const ProphageConfig& cfg) {
  Real hazard = std::max(cfg.k_spont, 0.0);
  if (damage_signal > 0.0 && cfg.k_max_damage > 0.0 && cfg.hill_n > 0.0) {
    if (cfg.K_damage <= 0.0) {
      hazard += cfg.k_max_damage;
    } else {
      const Real ratio = cfg.K_damage / damage_signal;
      const Real denom = 1.0 + std::pow(ratio, cfg.hill_n);
      if (std::isfinite(denom) && denom > 0.0) {
        hazard += cfg.k_max_damage / denom;
      }
    }
  }
  return hazard;
}

Int prophage_adsorption_receptor_index(const ProphageConfig& cfg) {
  if (cfg.adsorption_receptor.empty() || cfg.adsorption_receptor == "none") {
    return -1;
  }
  const auto receptor = receptor_type_from_name(cfg.adsorption_receptor);
  if (!receptor.has_value()) return -1;
  return static_cast<Int>(to_underlying(*receptor));
}

FixProphage::FixProphage(Simulation& sim, const ProphageConfig& cfg)
    : Fix("prophage", sim), cfg_(cfg) {}

void FixProphage::init() { /* no-op: configuration is read at construction */ }

Real FixProphage::compute_damage_signal(const Agent& agent,
                                        Int agent_index) const {
  Real nuclease_signal = 0.0;
  if (cfg_.nuclease_damage_weight > 0.0) {
    for (const ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
      const Real conc = sim_.qssa().sampled_nuclease_conc(agent_index, target);
      if (conc <= 0.0) continue;
      Real immunity_factor = 1.0;
      for (const auto& bi : agent.genome.bi_loci) {
        if (bi.target != target || !bi.is_nuclease) continue;
        const Real candidate =
            sim_.config().fixes.receptor.immunity_factor
            * bi.immunity_binding_affinity;
        immunity_factor = std::min(immunity_factor, candidate);
      }
      const Real expr =
          std::max(agent.receptor_expr[to_underlying(target)], 0.0);
      nuclease_signal += conc * expr * immunity_factor;
    }
  }

  const Real ros_signal = cfg_.ros_damage_weight > 0.0
      ? std::max(sim_.ros_induction_rate(agent), 0.0)
      : 0.0;

  return std::max(cfg_.damage_stimulus, 0.0)
      + cfg_.nuclease_damage_weight * nuclease_signal
      + cfg_.ros_damage_weight * ros_signal;
}

void FixProphage::compute(Real dt) {
  if (!cfg_.enabled || dt <= 0.0) return;

  auto& chem = sim_.chemical_field();
  const Int i_virion = chem.find(species::VIRION);
  if (i_virion < 0) return;

  const Real cell_vol = sim_.domain().cell_volume();
  if (cell_vol <= 0.0) return;

  Int live_owned = 0;
  for (Agent& a : sim_.agents()) {
    if (a.state != PhenoState::DEAD && !a.flags.is_ghost) {
      Int ix = 0;
      Int iy = 0;
      Int iz = 0;
      sim_.domain().pos_to_grid(a.x, ix, iy, iz);
      a.grid_cell = sim_.domain().cell_index(ix, iy, iz);
      ++live_owned;
    }
  }
  sim_.step_events().prophage_bacterium_seconds +=
      static_cast<Real>(live_owned) * dt;

  // Phase 1: Complete pending latent lysis for previously induced lysogens.
  resolve_latent_lysis(i_virion, cell_vol, dt);

  // Phase 2: Spontaneous + DNA-damage Hill prophage induction.
  run_hill_induction(i_virion, cell_vol, dt);

  // Phase 3: Virion decay/washout, receptor-mediated adsorption, superinfection
  // immunity, lysogenic conversion, and lytic infection.
  run_virion_adsorption_and_decay(i_virion, cell_vol, dt);

  chem.mark_host_reac_dirty(i_virion);
  sim_.step_events().refresh_prophage_daily_rate();
}

void FixProphage::resolve_latent_lysis(Int i_virion, Real cell_vol, Real dt) {
  for (Agent& a : sim_.agents()) {
    if (a.flags.is_ghost || a.state == PhenoState::DEAD) continue;
    if (!a.flags.prophage_induced) continue;

    a.timers.prophage_lysis_timer -= dt;
    if (a.timers.prophage_lysis_timer <= 0.0) {
      lyse_induced_agent(a, i_virion, cell_vol, dt);
    }
  }
}

void FixProphage::run_hill_induction(Int i_virion, Real cell_vol, Real dt) {
  auto& agents = sim_.agents();
  const Int n = agents.size();
  if (n <= 0) return;

  std::vector<double> damage_signals(static_cast<size_t>(n), 0.0);
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int i = 0; i < n; ++i) {
    const Agent& a = agents[i];
    if (!agent_eligible_for_induction(a)) continue;
    damage_signals[static_cast<size_t>(i)] = compute_damage_signal(a, i);
  }

  std::vector<double> induction_probs(static_cast<size_t>(n), 0.0);
  bool gpu_ok = false;
  if (sim_.gpu_active()) {
    sim_.agents_gpu().sync_from_host(agents);
    gpu_ok = gpu_compute_prophage_induction_probs(
        sim_.agents_gpu(), damage_signals, cfg_, dt, induction_probs);
  }

  if (!gpu_ok) {
    #ifdef GUTIBM_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (Int i = 0; i < n; ++i) {
      const Agent& a = agents[i];
      if (!agent_eligible_for_induction(a)) continue;
      const Real h = prophage_hill_induction_hazard(
          damage_signals[static_cast<size_t>(i)], cfg_);
      const Real p = h > 0.0 ? -std::expm1(-h * dt) : 0.0;
      induction_probs[static_cast<size_t>(i)] = std::clamp(p, 0.0, 1.0);
    }
  }

  auto& rng = sim_.rng();
  for (Int i = 0; i < n; ++i) {
    Agent& a = agents[i];
    if (!agent_eligible_for_induction(a)) continue;
    const Real p_ind = induction_probs[static_cast<size_t>(i)];
    if (p_ind > 0.0 && rng.bernoulli(p_ind)) {
      sim_.step_events().prophage_inductions++;
      if (cfg_.lysis_delay_s <= 0.0) {
        lyse_induced_agent(a, i_virion, cell_vol, dt);
      } else {
        a.flags.prophage_induced = true;
        a.timers.prophage_lysis_timer = cfg_.lysis_delay_s;
      }
    }
  }
}

void FixProphage::run_virion_adsorption_and_decay(Int i_virion, Real cell_vol,
                                                  Real dt) {
  auto& agents = sim_.agents();
  auto& chem = sim_.chemical_field();
  const auto& domain = sim_.domain();
  const auto& advection = sim_.advection();
  const Int n = agents.size();
  const Int ncells = chem.ncells();

  std::vector<double> infect_probs(static_cast<size_t>(n), 0.0);
  double total_decayed = 0.0;
  double total_adsorbed = 0.0;

  bool gpu_ok = false;
  if (sim_.gpu_active()) {
    sim_.agents_gpu().sync_from_host(agents);
    std::vector<double> reac_delta;
    gpu_ok = gpu_compute_prophage_adsorption(
        sim_.agents_gpu(), agents, sim_.chem_gpu(), chem, domain, advection,
        cfg_, i_virion, dt, infect_probs, reac_delta,
        total_decayed, total_adsorbed);
    if (gpu_ok) {
      for (Int c = 0; c < ncells; ++c) {
        const Real delta = reac_delta[static_cast<size_t>(c)];
        if (delta < 0.0 || delta > 0.0) {
          chem.reac(i_virion, c) += delta;
        }
      }
    }
  }

  if (!gpu_ok) {
    const Int receptor_idx = prophage_adsorption_receptor_index(cfg_);
    std::vector<Real> cell_adsorption_rate(static_cast<size_t>(ncells), 0.0);
    std::vector<Real> cell_adsorb_factor(static_cast<size_t>(ncells), 0.0);

    if (cfg_.k_adsorb > 0.0) {
      for (Int i = 0; i < n; ++i) {
        const Agent& a = agents[i];
        if (!agent_participates_in_adsorption(a)) continue;
        const Int cell = map_agent_cell_to_storage(chem, domain, a.grid_cell);
        if (cell < 0) continue;
        const Real expr = agent_adsorption_expression(a, receptor_idx);
        cell_adsorption_rate[static_cast<size_t>(cell)] +=
            cfg_.k_adsorb * expr / cell_vol;
      }
    }

    const Int x_begin = chem.slab_mode()
        ? chem.owned_storage_x_begin() : domain.local_grid_x_begin();
    const Int x_end = chem.slab_mode()
        ? chem.owned_storage_x_end() : domain.local_grid_x_end();
    const Int storage_nx = chem.slab_mode() ? chem.storage_nx() : domain.nx();
    const Int ny = chem.global_ny();
    const Int nz = chem.global_nz();

    for (Int iz = 0; iz < nz; ++iz) {
      const Real z_center =
          domain.lo()[2] + (static_cast<Real>(iz) + 0.5) * domain.dx_z();
      const Real washout = cfg_.virion_washout_scale > 0.0
          ? cfg_.virion_washout_scale * advection.washout_rate(z_center)
          : 0.0;
      const Real k_decay = std::max(cfg_.k_virion_decay, 0.0) + washout;

      for (Int iy = 0; iy < ny; ++iy) {
        for (Int ix = x_begin; ix < x_end; ++ix) {
          const Int cell = iz * (storage_nx * ny) + iy * storage_nx + ix;
          const Real conc = chem.conc(i_virion, cell);
          if (conc <= 0.0) continue;

          const Real a_total = cell_adsorption_rate[static_cast<size_t>(cell)];
          const Real lambda = k_decay + a_total;
          if (lambda <= 0.0) continue;

          const Real frac_removed = -std::expm1(-lambda * dt);
          const Real removed_conc = conc * frac_removed;
          const Real removed_virions = removed_conc * cell_vol;
          const Real decayed = removed_virions * (k_decay / lambda);
          const Real adsorbed = removed_virions * (a_total / lambda);

          chem.reac(i_virion, cell) -= removed_conc / dt;
          cell_adsorb_factor[static_cast<size_t>(cell)] =
              removed_virions / lambda;
          total_decayed += decayed;
          total_adsorbed += adsorbed;
        }
      }
    }

    if (cfg_.k_adsorb > 0.0) {
      #ifdef GUTIBM_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (Int i = 0; i < n; ++i) {
        const Agent& a = agents[i];
        if (!agent_participates_in_adsorption(a)
            || a.flags.has_prophage_immunity) {
          continue;
        }
        const Int cell = map_agent_cell_to_storage(chem, domain, a.grid_cell);
        if (cell < 0) continue;
        const Real factor = cell_adsorb_factor[static_cast<size_t>(cell)];
        if (factor <= 0.0) continue;
        const Real expr = agent_adsorption_expression(a, receptor_idx);
        const Real v_adsorbed = (cfg_.k_adsorb * expr / cell_vol) * factor;
        const Real p = v_adsorbed > 0.0 ? -std::expm1(-v_adsorbed) : 0.0;
        infect_probs[static_cast<size_t>(i)] = std::clamp(p, 0.0, 1.0);
      }
    }
  }

  sim_.step_events().virions_decayed += total_decayed;
  sim_.step_events().virions_adsorbed += total_adsorbed;

  auto& rng = sim_.rng();
  for (Int i = 0; i < n; ++i) {
    Agent& a = agents[i];
    if (!agent_participates_in_adsorption(a) || a.flags.has_prophage_immunity) {
      continue;
    }
    const Real p_inf = infect_probs[static_cast<size_t>(i)];
    if (p_inf > 0.0 && rng.bernoulli(p_inf)) {
      const bool lysogenize = cfg_.p_lysogeny > 0.0
          && (cfg_.p_lysogeny >= 1.0 || rng.bernoulli(cfg_.p_lysogeny));
      if (lysogenize) {
        a.flags.has_prophage_immunity = true;
        a.flags.has_inducible_prophage = true;
        sim_.step_events().prophage_lysogen_conversions++;
      } else {
        record_lytic_kill(a, i_virion, cell_vol, dt);
      }
    }
  }
}

void FixProphage::lyse_induced_agent(Agent& agent, Int i_virion, Real cell_vol,
                                     Real dt) {
  auto& chem = sim_.chemical_field();
  if (cfg_.virion_burst_size > 0.0) {
    const Int cell =
        map_agent_cell_to_storage(chem, sim_.domain(), agent.grid_cell);
    if (cell >= 0) {
      chem.reac(i_virion, cell) += cfg_.virion_burst_size / (cell_vol * dt);
      sim_.step_events().virions_released += cfg_.virion_burst_size;
    }
  }

  agent.state = PhenoState::DEAD;
  agent.flags.prophage_induced = false;
  agent.timers.prophage_lysis_timer = -1.0;
  sim_.step_events().mortality_lysis++;

  if (sim_.provenance_enabled()) {
    KillProvenanceEvent event;
    event.victim_id = agent.identity.tag;
    event.position = agent.x;
    event.strain = agent.identity.type;
    event.cause = ProvenanceCause::LYSIS;
    sim_.record_kill_provenance(event);
  }
  sim_.lineage_tracker().record_lysis(agent.identity.tag, agent.x,
                                      agent.genome.lineage_id);
}

void FixProphage::record_lytic_kill(Agent& agent, Int i_virion, Real cell_vol,
                                    Real dt) {
  auto& chem = sim_.chemical_field();
  if (cfg_.virion_burst_size > 0.0) {
    const Int cell =
        map_agent_cell_to_storage(chem, sim_.domain(), agent.grid_cell);
    if (cell >= 0) {
      chem.reac(i_virion, cell) += cfg_.virion_burst_size / (cell_vol * dt);
      sim_.step_events().virions_released += cfg_.virion_burst_size;
    }
  }

  agent.state = PhenoState::DEAD;
  sim_.step_events().prophage_lytic_kills++;

  if (sim_.provenance_enabled()) {
    KillProvenanceEvent event;
    event.victim_id = agent.identity.tag;
    event.position = agent.x;
    event.strain = agent.identity.type;
    event.cause = ProvenanceCause::PROPHAGE_LYTIC;
    sim_.record_kill_provenance(event);
  }
}

Real FixProphage::total_virions_in_field(const Simulation& sim) {
  const auto& chem = sim.chemical_field();
  const Int i_virion = chem.find(species::VIRION);
  if (i_virion < 0) return 0.0;

  const auto& domain = sim.domain();
  const Real cell_vol = domain.cell_volume();
  const Int x_begin = chem.slab_mode()
      ? chem.owned_storage_x_begin() : domain.local_grid_x_begin();
  const Int x_end = chem.slab_mode()
      ? chem.owned_storage_x_end() : domain.local_grid_x_end();
  const Int storage_nx = chem.slab_mode() ? chem.storage_nx() : domain.nx();
  const Int ny = chem.global_ny();
  const Int nz = chem.global_nz();

  Real local_sum = 0.0;
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = x_begin; ix < x_end; ++ix) {
        const Int cell = iz * (storage_nx * ny) + iy * storage_nx + ix;
        local_sum += chem.conc(i_virion, cell) * cell_vol;
      }
    }
  }

  Real global_sum = local_sum;
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized != 0 && finalized == 0 && domain.nprocs() > 1) {
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  }
#endif
  return global_sum;
}

}  // namespace gutibm
