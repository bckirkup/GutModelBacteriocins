#include "mechanism_capability.h"

namespace gutibm {

namespace {

void flag_if(bool differs, const char* key, std::vector<std::string>& out) {
  if (differs) {
    out.emplace_back(key);
  }
}

}  // namespace

// Every entry is "configured value != shipped default" — a sweep axis.
// When adding a mechanism-dynamics field consumed only by the live
// engine, flag it here or the capability limit silently widens.
std::vector<std::string> unexecuted_mechanism_params(
    const SimulationConfig& cfg) {
  const SimulationConfig def = InputParser::default_config();
  std::vector<std::string> out;

  // Fix plugin execution itself is unsupported — both lists only steer
  // the live FixRegistry.
  flag_if(!cfg.enabled_fixes.empty(), "enabled_fixes", out);
  flag_if(!cfg.disabled_fixes.empty(), "disabled_fixes", out);

  // Bacteriocin release / SOS induction.
  {
    const auto& v = cfg.fixes.bacteriocin;
    const auto& d = def.fixes.bacteriocin;
    flag_if(v.sos_lysis_prob != d.sos_lysis_prob,
            "fixes.bacteriocin.sos_lysis_prob", out);
    flag_if(v.sos_basal_rate != d.sos_basal_rate,
            "fixes.bacteriocin.sos_basal_rate", out);
    flag_if(v.sos_cross_induction_rate != d.sos_cross_induction_rate,
            "fixes.bacteriocin.sos_cross_induction_rate", out);
    flag_if(v.burst_release_tau != d.burst_release_tau,
            "fixes.bacteriocin.burst_release_tau", out);
    flag_if(v.microcin_mu_penalty != d.microcin_mu_penalty,
            "fixes.bacteriocin.microcin_mu_penalty", out);
    flag_if(v.D_free_colicin != d.D_free_colicin,
            "fixes.bacteriocin.D_free_colicin", out);
    flag_if(v.mucin_charge.r_min != d.mucin_charge.r_min,
            "fixes.bacteriocin.mucin_charge.r_min", out);
    flag_if(v.mucin_charge.amplitude != d.mucin_charge.amplitude,
            "fixes.bacteriocin.mucin_charge.amplitude", out);
    flag_if(v.mucin_charge.dz_half != d.mucin_charge.dz_half,
            "fixes.bacteriocin.mucin_charge.dz_half", out);
    flag_if(v.mucin_charge.width != d.mucin_charge.width,
            "fixes.bacteriocin.mucin_charge.width", out);
    flag_if(v.mucin_charge.ph != d.mucin_charge.ph,
            "fixes.bacteriocin.mucin_charge.ph", out);
  }

  // Receptor affinities and kill rates.
  {
    const auto& v = cfg.fixes.receptor;
    const auto& d = def.fixes.receptor;
    flag_if(v.kd_b12_btuB != d.kd_b12_btuB,
            "fixes.receptor.kd_b12_btuB", out);
    flag_if(v.kd_colicinE_btuB != d.kd_colicinE_btuB,
            "fixes.receptor.kd_colicinE_btuB", out);
    flag_if(v.kd_enterobactin != d.kd_enterobactin,
            "fixes.receptor.kd_enterobactin", out);
    flag_if(v.kd_colicinB_fepA != d.kd_colicinB_fepA,
            "fixes.receptor.kd_colicinB_fepA", out);
    flag_if(v.kd_lin_enterobactin != d.kd_lin_enterobactin,
            "fixes.receptor.kd_lin_enterobactin", out);
    flag_if(v.kd_colicinIa_cirA != d.kd_colicinIa_cirA,
            "fixes.receptor.kd_colicinIa_cirA", out);
    flag_if(v.kd_colicinM_fhuA != d.kd_colicinM_fhuA,
            "fixes.receptor.kd_colicinM_fhuA", out);
    flag_if(v.kd_ferrichrome != d.kd_ferrichrome,
            "fixes.receptor.kd_ferrichrome", out);
    flag_if(v.cirA_linearized_fraction != d.cirA_linearized_fraction,
            "fixes.receptor.cirA_linearized_fraction", out);
    flag_if(v.kill_rate_colicin != d.kill_rate_colicin,
            "fixes.receptor.kill_rate_colicin", out);
    flag_if(v.kill_rate_microcin != d.kill_rate_microcin,
            "fixes.receptor.kill_rate_microcin", out);
    flag_if(v.immunity_factor != d.immunity_factor,
            "fixes.receptor.immunity_factor", out);
  }

  // HGT / genome dynamics.
  {
    const auto& v = cfg.fixes.conjugation;
    const auto& d = def.fixes.conjugation;
    flag_if(v.pili_length != d.pili_length,
            "fixes.conjugation.pili_length", out);
    flag_if(v.base_transfer_rate != d.base_transfer_rate,
            "fixes.conjugation.base_transfer_rate", out);
    flag_if(v.shear_critical != d.shear_critical,
            "fixes.conjugation.shear_critical", out);
    flag_if(v.plasmid_copy_cost != d.plasmid_copy_cost,
            "fixes.conjugation.plasmid_copy_cost", out);
    flag_if(v.pili_heterogeneity != d.pili_heterogeneity,
            "fixes.conjugation.pili_heterogeneity", out);
    flag_if(v.pili_length_min != d.pili_length_min,
            "fixes.conjugation.pili_length_min", out);
    flag_if(v.pili_length_max != d.pili_length_max,
            "fixes.conjugation.pili_length_max", out);
  }
  {
    const auto& v = cfg.fixes.mutation;
    const auto& d = def.fixes.mutation;
    flag_if(v.bi_duplication_rate != d.bi_duplication_rate,
            "fixes.mutation.bi_duplication_rate", out);
    flag_if(v.bi_recombination_rate != d.bi_recombination_rate,
            "fixes.mutation.bi_recombination_rate", out);
    flag_if(v.receptor_mutation_rate != d.receptor_mutation_rate,
            "fixes.mutation.receptor_mutation_rate", out);
    flag_if(v.super_killer_rate != d.super_killer_rate,
            "fixes.mutation.super_killer_rate", out);
    flag_if(v.partial_resistance_rate != d.partial_resistance_rate,
            "fixes.mutation.partial_resistance_rate", out);
    flag_if(v.receptor_reduction != d.receptor_reduction,
            "fixes.mutation.receptor_reduction", out);
    flag_if(v.max_bi_loci != d.max_bi_loci,
            "fixes.mutation.max_bi_loci", out);
    flag_if(v.immunity_escape_prob != d.immunity_escape_prob,
            "fixes.mutation.immunity_escape_prob", out);
    flag_if(v.escape_affinity_lo != d.escape_affinity_lo,
            "fixes.mutation.escape_affinity_lo", out);
    flag_if(v.escape_affinity_hi != d.escape_affinity_hi,
            "fixes.mutation.escape_affinity_hi", out);
  }

  // Contact mechanics.
  {
    const auto& v = cfg.fixes.mechanics;
    const auto& d = def.fixes.mechanics;
    flag_if(v.hertz_k != d.hertz_k, "fixes.mechanics.hertz_k", out);
    flag_if(v.hertzian_enabled != d.hertzian_enabled,
            "fixes.mechanics.hertzian_enabled", out);
    flag_if(v.adhesion_enabled != d.adhesion_enabled,
            "fixes.mechanics.adhesion_enabled", out);
    flag_if(v.adhesion_strength != d.adhesion_strength,
            "fixes.mechanics.adhesion_strength", out);
    flag_if(v.adhesion_range != d.adhesion_range,
            "fixes.mechanics.adhesion_range", out);
  }

  // Toxin transport override axes (retardation/diffusion/burst): carried
  // on the cluster as state but never read at these layers.
  flag_if(!cfg.plasmid_overrides.empty(), "plasmid_overrides", out);

  // Oxygen / ROS mechanisms. The three anaerobic factors the lumen
  // physio reads are excluded — they are consumed.
  {
    const auto& v = cfg.chem_env.oxygen;
    const auto& d = def.chem_env.oxygen;
    flag_if(v.enabled != d.enabled, "chem_env.oxygen.enabled", out);
    flag_if(v.epithelial_conc != d.epithelial_conc,
            "chem_env.oxygen.epithelial_conc", out);
    flag_if(v.delivery_uptake_enabled != d.delivery_uptake_enabled,
            "chem_env.oxygen.delivery_uptake_enabled", out);
    flag_if(v.respiration_driver != d.respiration_driver,
            "chem_env.oxygen.respiration_driver", out);
    flag_if(v.ros_driver != d.ros_driver, "chem_env.oxygen.ros_driver", out);
    flag_if(v.D_free != d.D_free, "chem_env.oxygen.D_free", out);
    flag_if(v.Km != d.Km, "chem_env.oxygen.Km", out);
    flag_if(v.boost_max != d.boost_max, "chem_env.oxygen.boost_max", out);
    flag_if(v.metabolic_switch_enabled != d.metabolic_switch_enabled,
            "chem_env.oxygen.metabolic_switch_enabled", out);
    flag_if(v.mu_crit != d.mu_crit, "chem_env.oxygen.mu_crit", out);
    flag_if(v.aerobic_mu_factor != d.aerobic_mu_factor,
            "chem_env.oxygen.aerobic_mu_factor", out);
    flag_if(v.aerobic_carbon_cost_factor != d.aerobic_carbon_cost_factor,
            "chem_env.oxygen.aerobic_carbon_cost_factor", out);
    flag_if(v.tau_metabolic_switch != d.tau_metabolic_switch,
            "chem_env.oxygen.tau_metabolic_switch", out);
    flag_if(v.anaerobic_maintenance_factor != d.anaerobic_maintenance_factor,
            "chem_env.oxygen.anaerobic_maintenance_factor", out);
    flag_if(v.q_consumption != d.q_consumption,
            "chem_env.oxygen.q_consumption", out);
    flag_if(v.q_maintenance != d.q_maintenance,
            "chem_env.oxygen.q_maintenance", out);
    flag_if(v.vbf_sink != d.vbf_sink, "chem_env.oxygen.vbf_sink", out);
    flag_if(v.k_ROS != d.k_ROS, "chem_env.oxygen.k_ROS", out);
    flag_if(v.k_ROS_respiratory != d.k_ROS_respiratory,
            "chem_env.oxygen.k_ROS_respiratory", out);
    flag_if(v.k_ROS_funded != d.k_ROS_funded,
            "chem_env.oxygen.k_ROS_funded", out);
  }

  // Receptor-kill substrate and B12 competition pool.
  flag_if(cfg.b12_initial_conc != def.b12_initial_conc, "b12_initial_conc",
          out);

  // Agent behaviors that only exist in the live engine.
  {
    const auto& v = cfg.cell_bio.motility;
    const auto& d = def.cell_bio.motility;
    flag_if(v.enabled != d.enabled, "cell_bio.motility.enabled", out);
    flag_if(v.swim_speed != d.swim_speed, "cell_bio.motility.swim_speed",
            out);
    flag_if(v.run_mean_duration != d.run_mean_duration,
            "cell_bio.motility.run_mean_duration", out);
    flag_if(v.stop_probability != d.stop_probability,
            "cell_bio.motility.stop_probability", out);
    flag_if(v.stop_duration != d.stop_duration,
            "cell_bio.motility.stop_duration", out);
    flag_if(v.chemotaxis_enabled != d.chemotaxis_enabled,
            "cell_bio.motility.chemotaxis_enabled", out);
    flag_if(v.chi_carbon != d.chi_carbon, "cell_bio.motility.chi_carbon",
            out);
    flag_if(v.chemotaxis_threshold != d.chemotaxis_threshold,
            "cell_bio.motility.chemotaxis_threshold", out);
    flag_if(v.aerotaxis_enabled != d.aerotaxis_enabled,
            "cell_bio.motility.aerotaxis_enabled", out);
    flag_if(v.aerotaxis_sensitivity != d.aerotaxis_sensitivity,
            "cell_bio.motility.aerotaxis_sensitivity", out);
    flag_if(v.energy_taxis_enabled != d.energy_taxis_enabled,
            "cell_bio.motility.energy_taxis_enabled", out);
    flag_if(v.energy_taxis_floor != d.energy_taxis_floor,
            "cell_bio.motility.energy_taxis_floor", out);
    flag_if(v.surface_sensing_enabled != d.surface_sensing_enabled,
            "cell_bio.motility.surface_sensing_enabled", out);
    flag_if(v.surface_sensing_depth != d.surface_sensing_depth,
            "cell_bio.motility.surface_sensing_depth", out);
  }
  {
    const auto& v = cfg.cell_bio.fur;
    const auto& d = def.cell_bio.fur;
    flag_if(v.enabled != d.enabled, "cell_bio.fur.enabled", out);
    flag_if(v.Km != d.Km, "cell_bio.fur.Km", out);
    flag_if(v.upregulation_max != d.upregulation_max,
            "cell_bio.fur.upregulation_max", out);
    flag_if(v.receptor_max != d.receptor_max, "cell_bio.fur.receptor_max",
            out);
  }
  {
    const auto& v = cfg.cell_bio.cdi;
    const auto& d = def.cell_bio.cdi;
    flag_if(v.enabled != d.enabled, "cell_bio.cdi.enabled", out);
    flag_if(v.kill_rate != d.kill_rate, "cell_bio.cdi.kill_rate", out);
    flag_if(v.contact_radius != d.contact_radius,
            "cell_bio.cdi.contact_radius", out);
    flag_if(v.corpse_persistence != d.corpse_persistence,
            "cell_bio.cdi.corpse_persistence", out);
  }
  {
    const auto& v = cfg.quorum_sensing;
    const auto& d = def.quorum_sensing;
    flag_if(v.enabled != d.enabled, "quorum_sensing.enabled", out);
    flag_if(v.ai2_basal_rate != d.ai2_basal_rate,
            "quorum_sensing.ai2_basal_rate", out);
    flag_if(v.ai2_growth_coupled != d.ai2_growth_coupled,
            "quorum_sensing.ai2_growth_coupled", out);
    flag_if(v.lsr_vmax != d.lsr_vmax, "quorum_sensing.lsr_vmax", out);
    flag_if(v.lsr_km != d.lsr_km, "quorum_sensing.lsr_km", out);
    flag_if(v.ai2_D_free != d.ai2_D_free, "quorum_sensing.ai2_D_free", out);
    flag_if(v.ai2_decay_rate != d.ai2_decay_rate,
            "quorum_sensing.ai2_decay_rate", out);
    flag_if(v.ai2_chemotaxis_enabled != d.ai2_chemotaxis_enabled,
            "quorum_sensing.ai2_chemotaxis_enabled", out);
    flag_if(v.chi_ai2 != d.chi_ai2, "quorum_sensing.chi_ai2", out);
  }
  flag_if(cfg.immigration.enabled, "immigration.enabled", out);

  // QSSA toxin transport solver and its policies.
  {
    const auto& v = cfg.qssa;
    const auto& d = def.qssa;
    flag_if(v.toxin_cutoff != d.toxin_cutoff, "qssa.toxin_cutoff", out);
    flag_if(v.lumen_transfer_length != d.lumen_transfer_length,
            "qssa.lumen_transfer_length", out);
    flag_if(v.lumen_transfer_basis != d.lumen_transfer_basis,
            "qssa.lumen_transfer_basis", out);
    flag_if(v.image_series_relative_tolerance !=
                d.image_series_relative_tolerance,
            "qssa.image_series_relative_tolerance", out);
    flag_if(v.image_series_max_shells != d.image_series_max_shells,
            "qssa.image_series_max_shells", out);
    flag_if(v.image_series_mode != d.image_series_mode,
            "qssa.image_series_mode", out);
    flag_if(v.low_screening_policy != d.low_screening_policy,
            "qssa.low_screening_policy", out);
    flag_if(v.drift_envelope_policy != d.drift_envelope_policy,
            "qssa.drift_envelope_policy", out);
    flag_if(v.drift_correction != d.drift_correction,
            "qssa.drift_correction", out);
    flag_if(v.nutrient_cutoff != d.nutrient_cutoff, "qssa.nutrient_cutoff",
            out);
    flag_if(v.colicin_release_rate != d.colicin_release_rate,
            "qssa.colicin_release_rate", out);
    flag_if(v.microcin_secretion != d.microcin_secretion,
            "qssa.microcin_secretion", out);
    flag_if(v.use_fmm != d.use_fmm, "qssa.use_fmm", out);
    flag_if(v.toxin_evaluation != d.toxin_evaluation,
            "qssa.toxin_evaluation", out);
    flag_if(v.toxin_lumping != d.toxin_lumping, "qssa.toxin_lumping", out);
    flag_if(v.fmm_theta != d.fmm_theta, "qssa.fmm_theta", out);
    flag_if(v.fmm_expansion_order != d.fmm_expansion_order,
            "qssa.fmm_expansion_order", out);
  }

  return out;
}

}  // namespace gutibm
