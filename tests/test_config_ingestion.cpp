/* -----------------------------------------------------------------------
   GutIBM – Exhaustive config-ingestion tracking test

   Guarantees that EVERY input variable accepted by the config parser is
   actually plumbed into SimulationConfig (the struct Simulation::init
   consumes). For each key we set a non-default sentinel value and assert:
     1. the value lands in the expected SimulationConfig field, and
     2. the field differs from its default (i.e. the key is not dead-wired
        to a constant / silently ignored).

   Both parser entry points are exercised:
     - InputParser::apply_flat_key (used by the legacy line parser and by
       the JSON scalar path), and
     - the full JSON document path (InputParser::parse), which routes
       scalars through ConfigJson::apply_json_scalar -> apply_flat_key.

   A completeness guard scans the parser sources for every `key == "..."`
   literal and asserts the set of parsed keys is exactly the set covered
   here — so adding a new config key without a probe fails CI, and removing
   a key's wiring (regression) also fails CI.
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "config_json.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifndef GUTIBM_SOURCE_DIR
#define GUTIBM_SOURCE_DIR "."
#endif

using namespace gutibm;

namespace {

struct FailureCounter {
  int value = 0;
};

FailureCounter& failure_counter() {
  static FailureCounter counter;
  return counter;
}

void record_failure(const std::string& msg) {
  std::cerr << "FAIL: " << msg << "\n";
  ++failure_counter().value;
}

void expect(bool cond, const std::string& msg) {
  if (!cond) record_failure(msg);
}

enum class Kind { Real, Int, Bool, Str };

using RealGet = std::function<double(const SimulationConfig&)>;
using IntGet = std::function<long long(const SimulationConfig&)>;
using BoolGet = std::function<bool(const SimulationConfig&)>;
using StrGet = std::function<std::string(const SimulationConfig&)>;

// A single tracked input variable and how to read it back out of the config.
struct Probe {
  std::string key;
  Kind kind = Kind::Real;
  bool primary = true;  // include in the combined full-document JSON test
  RealGet gr;
  IntGet gi;
  BoolGet gb;
  StrGet gs;
  std::string sval;  // sentinel for string keys
};

Probe R(std::string key, RealGet g, bool primary = true) {
  Probe p;
  p.key = std::move(key);
  p.kind = Kind::Real;
  p.primary = primary;
  p.gr = std::move(g);
  return p;
}

Probe I(std::string key, IntGet g, bool primary = true) {
  Probe p;
  p.key = std::move(key);
  p.kind = Kind::Int;
  p.primary = primary;
  p.gi = std::move(g);
  return p;
}

Probe B(std::string key, BoolGet g, bool primary = true) {
  Probe p;
  p.key = std::move(key);
  p.kind = Kind::Bool;
  p.primary = primary;
  p.gb = std::move(g);
  return p;
}

Probe S(std::string key, StrGet g, std::string sentinel, bool primary = true) {
  Probe p;
  p.key = std::move(key);
  p.kind = Kind::Str;
  p.primary = primary;
  p.gs = std::move(g);
  p.sval = std::move(sentinel);
  return p;
}

// Namespaced keys (Spec 1 / Spec 3) accept both a dotted and an underscore
// spelling that map to the same field. Track both; only the dotted form is
// "primary" for the combined-document test to avoid writing a field twice.
void add_ns_real(std::vector<Probe>& v, const char* dotted, const char* under, const RealGet& g) {
  v.push_back(R(dotted, g, true));
  v.push_back(R(under, g, false));
}

void add_ns_bool(std::vector<Probe>& v, const char* dotted, const char* under, const BoolGet& g) {
  v.push_back(B(dotted, g, true));
  v.push_back(B(under, g, false));
}

void add_ns_int(std::vector<Probe>& v, const char* dotted, const char* under, const IntGet& g) {
  v.push_back(I(dotted, g, true));
  v.push_back(I(under, g, false));
}

const ChemicalSpec& carbon_spec(const SimulationConfig& c) {
  for (const auto& s : c.chemicals) {
    if (s.name == "carbon") return s;
  }
  static const ChemicalSpec kMissing{};
  return kMissing;
}

std::vector<Probe> build_probes() {
  std::vector<Probe> v;

  // ── Time control ───────────────────────────────────────────────────────
  v.push_back(R("total_time", [](const SimulationConfig& c) { return c.time.total_time; }));
  v.push_back(R("bio_dt", [](const SimulationConfig& c) { return c.time.bio_dt; }));
  v.push_back(R("output_interval", [](const SimulationConfig& c) { return c.time.output_interval; }));
  v.push_back(I("seed", [](const SimulationConfig& c) { return static_cast<long long>(c.seed); }));

  // ── Domain ───────────────────────────────────────────────────────────────
  v.push_back(R("grid_dx", [](const SimulationConfig& c) { return c.domain.grid_dx; }));
  v.push_back(R("domain_x", [](const SimulationConfig& c) { return c.domain.hi[0]; }));
  v.push_back(R("domain_y", [](const SimulationConfig& c) { return c.domain.hi[1]; }));
  v.push_back(R("domain_z", [](const SimulationConfig& c) { return c.domain.hi[2]; }));
  v.push_back(R("hash_cell_size", [](const SimulationConfig& c) { return c.domain.hash_cell_size; }));
  v.push_back(R("ghost_width", [](const SimulationConfig& c) { return c.domain.ghost_width; }));

  // ── Advection / crypts / peristalsis ──────────────────────────────────────
  v.push_back(R("mucus_thickness", [](const SimulationConfig& c) { return c.advection.mucus_thickness; }));
  v.push_back(R("radial_turnover", [](const SimulationConfig& c) { return c.advection.radial_turnover; }));
  v.push_back(R("distal_transit", [](const SimulationConfig& c) { return c.advection.distal_transit_time; }));
  v.push_back(R("distal_length", [](const SimulationConfig& c) { return c.advection.distal_length; }));
  v.push_back(R("profile_alpha", [](const SimulationConfig& c) { return c.advection.profile_alpha; }));
  v.push_back(B("taylor_aris_enabled", [](const SimulationConfig& c) { return c.advection.taylor_aris_enabled; }));
  v.push_back(B("peristaltic_enabled", [](const SimulationConfig& c) { return c.advection.peristaltic_enabled; }));
  v.push_back(R("peristaltic_period", [](const SimulationConfig& c) { return c.advection.peristaltic_period; }));
  v.push_back(R("peristaltic_amplitude", [](const SimulationConfig& c) { return c.advection.peristaltic_amplitude; }));
  v.push_back(R("peristaltic_wavelength", [](const SimulationConfig& c) { return c.advection.peristaltic_wavelength; }));
  v.push_back(B("crypts_enabled", [](const SimulationConfig& c) { return c.advection.crypts_enabled; }));
  v.push_back(R("crypt_depth", [](const SimulationConfig& c) { return c.advection.crypt_depth; }));
  v.push_back(R("crypt_exit_rate", [](const SimulationConfig& c) { return c.advection.crypt_exit_rate; }));
  v.push_back(R("crypt_entry_rate", [](const SimulationConfig& c) { return c.advection.crypt_entry_rate; }));
  v.push_back(I("crypt_carrying_capacity", [](const SimulationConfig& c) { return static_cast<long long>(c.advection.crypt_carrying_capacity); }));

  // ── QSSA / FMM ────────────────────────────────────────────────────────────
  v.push_back(R("toxin_cutoff", [](const SimulationConfig& c) { return c.qssa.toxin_cutoff; }));
  v.push_back(R("nutrient_cutoff", [](const SimulationConfig& c) { return c.qssa.nutrient_cutoff; }));
  v.push_back(R("colicin_release_rate", [](const SimulationConfig& c) { return c.qssa.colicin_release_rate; }));
  v.push_back(R("microcin_secretion", [](const SimulationConfig& c) { return c.qssa.microcin_secretion; }));
  v.push_back(B("use_fmm", [](const SimulationConfig& c) { return c.qssa.use_fmm; }));
  v.push_back(R("fmm_theta", [](const SimulationConfig& c) { return c.qssa.fmm_theta; }));
  v.push_back(I("fmm_expansion_order", [](const SimulationConfig& c) { return static_cast<long long>(c.qssa.fmm_expansion_order); }));

  // ── VBF ───────────────────────────────────────────────────────────────────
  v.push_back(R("vbf_density", [](const SimulationConfig& c) { return c.vbf.density; }));
  v.push_back(R("vbf_viscosity", [](const SimulationConfig& c) { return c.vbf.viscosity; }));
  v.push_back(R("vbf_drag_coeff", [](const SimulationConfig& c) { return c.vbf.drag_coeff; }));
  v.push_back(R("vbf_nutrient_sink", [](const SimulationConfig& c) { return c.vbf.nutrient_sink; }));
  v.push_back(R("vbf_mucin_liberation", [](const SimulationConfig& c) { return c.vbf.mucin_liberation; }));
  v.push_back(R("vbf_carrying_cap", [](const SimulationConfig& c) { return c.vbf.carrying_cap; }));
  v.push_back(B("vbf_mucin_z_gradient", [](const SimulationConfig& c) { return c.vbf.mucin_z_gradient_enabled; }));
  v.push_back(R("vbf_mucin_z_lambda", [](const SimulationConfig& c) { return c.vbf.mucin_z_gradient_lambda; }));
  v.push_back(R("vbf_carbon_sink_vmax", [](const SimulationConfig& c) { return c.vbf.carbon_sink_vmax; }));
  v.push_back(R("vbf_carbon_sink_km", [](const SimulationConfig& c) { return c.vbf.carbon_sink_km; }));

  // ── Carbon z-gradient + bacteriocin SOS ───────────────────────────────────
  v.push_back(B("carbon_z_gradient", [](const SimulationConfig& c) { return carbon_spec(c).z_gradient_enabled; }));
  v.push_back(R("carbon_z_lambda", [](const SimulationConfig& c) { return carbon_spec(c).z_gradient_lambda; }));
  v.push_back(R("sos_lysis_prob", [](const SimulationConfig& c) { return c.fixes.bacteriocin.sos_lysis_prob; }));
  v.push_back(R("sos_basal_rate", [](const SimulationConfig& c) { return c.fixes.bacteriocin.sos_basal_rate; }));
  v.push_back(R("sos_cross_induction_rate", [](const SimulationConfig& c) { return c.fixes.bacteriocin.sos_cross_induction_rate; }));
  v.push_back(R("retardation_basic", [](const SimulationConfig& c) { return c.fixes.bacteriocin.retardation_basic; }));
  v.push_back(R("retardation_acidic", [](const SimulationConfig& c) { return c.fixes.bacteriocin.retardation_acidic; }));
  v.push_back(R("retardation_neutral", [](const SimulationConfig& c) { return c.fixes.bacteriocin.retardation_neutral; }));
  v.push_back(R("D_free_colicin", [](const SimulationConfig& c) { return c.fixes.bacteriocin.D_free_colicin; }));
  v.push_back(R("burst_molecules", [](const SimulationConfig& c) { return c.fixes.bacteriocin.burst_molecules; }));
  v.push_back(R("microcin_mu_penalty", [](const SimulationConfig& c) { return c.fixes.bacteriocin.microcin_mu_penalty; }));

  // ── Metabolism Fix tunables ───────────────────────────────────────────────
  v.push_back(R("division_threshold", [](const SimulationConfig& c) { return c.fixes.metabolism.division_threshold; }));
  v.push_back(R("maintenance_rate", [](const SimulationConfig& c) { return c.fixes.metabolism.maintenance_rate; }));
  v.push_back(R("metE_penalty", [](const SimulationConfig& c) { return c.fixes.metabolism.metE_penalty; }));
  v.push_back(R("metE_acetate_km", [](const SimulationConfig& c) { return c.fixes.metabolism.metE_acetate_km; }));
  v.push_back(R("metE_acetate_max_factor", [](const SimulationConfig& c) { return c.fixes.metabolism.metE_acetate_max_factor; }));
  v.push_back(R("eut_km", [](const SimulationConfig& c) { return c.fixes.metabolism.eut_km; }));
  v.push_back(R("eut_max_penalty", [](const SimulationConfig& c) { return c.fixes.metabolism.eut_max_penalty; }));
  v.push_back(R("km_iron_primary", [](const SimulationConfig& c) { return c.fixes.metabolism.km_iron_primary; }));
  v.push_back(R("km_iron_iroN", [](const SimulationConfig& c) { return c.fixes.metabolism.km_iron_iroN; }));
  v.push_back(R("km_iron_iutA", [](const SimulationConfig& c) { return c.fixes.metabolism.km_iron_iutA; }));
  v.push_back(R("km_iron_fiu", [](const SimulationConfig& c) { return c.fixes.metabolism.km_iron_fiu; }));

  // ── Mechanics Fix tunables ────────────────────────────────────────────────
  v.push_back(R("hertz_k", [](const SimulationConfig& c) { return c.fixes.mechanics.hertz_k; }));
  v.push_back(B("hertzian_enabled", [](const SimulationConfig& c) { return c.fixes.mechanics.hertzian_enabled; }));
  v.push_back(B("adhesion_enabled", [](const SimulationConfig& c) { return c.fixes.mechanics.adhesion_enabled; }));
  v.push_back(R("adhesion_strength", [](const SimulationConfig& c) { return c.fixes.mechanics.adhesion_strength; }));
  v.push_back(R("adhesion_range", [](const SimulationConfig& c) { return c.fixes.mechanics.adhesion_range; }));

  // ── Receptor Fix tunables ─────────────────────────────────────────────────
  v.push_back(R("kd_b12_btuB", [](const SimulationConfig& c) { return c.fixes.receptor.kd_b12_btuB; }));
  // Alias for kd_b12_btuB (writes the same field); non-primary so it does not
  // collide with kd_b12_btuB in the combined full-document JSON test.
  v.push_back(R("kd_corrinoid_btuB", [](const SimulationConfig& c) { return c.fixes.receptor.kd_b12_btuB; }, false));
  v.push_back(R("kd_colicinE_btuB", [](const SimulationConfig& c) { return c.fixes.receptor.kd_colicinE_btuB; }));
  v.push_back(R("kd_enterobactin", [](const SimulationConfig& c) { return c.fixes.receptor.kd_enterobactin; }));
  v.push_back(R("kd_colicinB_fepA", [](const SimulationConfig& c) { return c.fixes.receptor.kd_colicinB_fepA; }));
  v.push_back(R("kd_lin_enterobactin", [](const SimulationConfig& c) { return c.fixes.receptor.kd_lin_enterobactin; }));
  v.push_back(R("kd_colicinIa_cirA", [](const SimulationConfig& c) { return c.fixes.receptor.kd_colicinIa_cirA; }));
  v.push_back(R("kill_rate_colicin", [](const SimulationConfig& c) { return c.fixes.receptor.kill_rate_colicin; }));
  v.push_back(R("kill_rate_microcin", [](const SimulationConfig& c) { return c.fixes.receptor.kill_rate_microcin; }));
  v.push_back(R("immunity_factor", [](const SimulationConfig& c) { return c.fixes.receptor.immunity_factor; }));

  // ── Conjugation Fix tunables ──────────────────────────────────────────────
  v.push_back(R("pili_length", [](const SimulationConfig& c) { return c.fixes.conjugation.pili_length; }));
  v.push_back(R("base_transfer_rate", [](const SimulationConfig& c) { return c.fixes.conjugation.base_transfer_rate; }));
  v.push_back(R("shear_critical", [](const SimulationConfig& c) { return c.fixes.conjugation.shear_critical; }));
  v.push_back(R("plasmid_copy_cost", [](const SimulationConfig& c) { return c.fixes.conjugation.plasmid_copy_cost; }));
  v.push_back(B("pili_heterogeneity", [](const SimulationConfig& c) { return c.fixes.conjugation.pili_heterogeneity; }));
  v.push_back(R("pili_length_min", [](const SimulationConfig& c) { return c.fixes.conjugation.pili_length_min; }));
  v.push_back(R("pili_length_max", [](const SimulationConfig& c) { return c.fixes.conjugation.pili_length_max; }));

  // ── Mutation Fix tunables ─────────────────────────────────────────────────
  v.push_back(R("bi_duplication_rate", [](const SimulationConfig& c) { return c.fixes.mutation.bi_duplication_rate; }));
  v.push_back(R("bi_recombination_rate", [](const SimulationConfig& c) { return c.fixes.mutation.bi_recombination_rate; }));
  v.push_back(R("receptor_mutation_rate", [](const SimulationConfig& c) { return c.fixes.mutation.receptor_mutation_rate; }));
  v.push_back(R("super_killer_rate", [](const SimulationConfig& c) { return c.fixes.mutation.super_killer_rate; }));
  v.push_back(R("partial_resistance_rate", [](const SimulationConfig& c) { return c.fixes.mutation.partial_resistance_rate; }));
  v.push_back(R("receptor_reduction", [](const SimulationConfig& c) { return c.fixes.mutation.receptor_reduction; }));
  v.push_back(I("max_bi_loci", [](const SimulationConfig& c) { return static_cast<long long>(c.fixes.mutation.max_bi_loci); }));
  v.push_back(R("immunity_escape_prob", [](const SimulationConfig& c) { return c.fixes.mutation.immunity_escape_prob; }));
  v.push_back(R("escape_affinity_lo", [](const SimulationConfig& c) { return c.fixes.mutation.escape_affinity_lo; }));
  v.push_back(R("escape_affinity_hi", [](const SimulationConfig& c) { return c.fixes.mutation.escape_affinity_hi; }));
  v.push_back(R("compensatory_rate", [](const SimulationConfig& c) { return c.fixes.mutation.compensatory_rate; }));
  v.push_back(R("compensatory_reduction", [](const SimulationConfig& c) { return c.fixes.mutation.compensatory_reduction; }));

  // ── Output / checkpoint ───────────────────────────────────────────────────
  v.push_back(S("hdf5_file", [](const SimulationConfig& c) { return c.hdf5.filename; }, "probe_output.h5"));
  v.push_back(S("hdf5.file", [](const SimulationConfig& c) { return c.hdf5.filename; }, "probe_output.h5"));
  add_ns_int(v, "hdf5.schedule.summary", "hdf5_schedule_summary",
             [](const SimulationConfig& c) { return static_cast<long long>(c.hdf5.schedule.summary); });
  add_ns_int(v, "hdf5.schedule.agents", "hdf5_schedule_agents",
             [](const SimulationConfig& c) { return static_cast<long long>(c.hdf5.schedule.agents); });
  add_ns_int(v, "hdf5.schedule.grid", "hdf5_schedule_grid",
             [](const SimulationConfig& c) { return static_cast<long long>(c.hdf5.schedule.grid); });
  add_ns_int(v, "hdf5.schedule.lineage", "hdf5_schedule_lineage",
             [](const SimulationConfig& c) { return static_cast<long long>(c.hdf5.schedule.lineage); });
  add_ns_int(v, "hdf5.schedule.genome", "hdf5_schedule_genome",
             [](const SimulationConfig& c) { return static_cast<long long>(c.hdf5.schedule.genome); });
  add_ns_bool(v, "hdf5.enabled", "hdf5_enabled", [](const SimulationConfig& c) { return c.hdf5.enabled; });
  v.push_back(S("hdf5.compression", [](const SimulationConfig& c) { return c.hdf5.compression; }, "gzip"));
  add_ns_int(v, "hdf5.compression_level", "hdf5_compression_level",
             [](const SimulationConfig& c) { return static_cast<long long>(c.hdf5.compression_level); });
  v.push_back(S("checkpoint_file", [](const SimulationConfig& c) { return c.checkpoint.file; }, "probe_checkpoint.h5"));
  v.push_back(S("checkpoint_step", [](const SimulationConfig& c) { return c.checkpoint.step; }, "step_000007"));

  // ── Adaptive timestep ─────────────────────────────────────────────────────
  v.push_back(B("adaptive_dt_enabled", [](const SimulationConfig& c) { return c.adaptive_dt.enabled; }));
  v.push_back(R("dt_min", [](const SimulationConfig& c) { return c.adaptive_dt.min; }));
  v.push_back(R("dt_max", [](const SimulationConfig& c) { return c.adaptive_dt.max; }));
  v.push_back(R("dt_safety", [](const SimulationConfig& c) { return c.adaptive_dt.safety; }));
  v.push_back(R("dt_growth_limit", [](const SimulationConfig& c) { return c.adaptive_dt.growth_limit; }));

  // ── Misc / GPU / profiling ────────────────────────────────────────────────
  v.push_back(B("gpu_enabled", [](const SimulationConfig& c) { return c.gpu.enabled; }));
  v.push_back(I("gpu_device_id", [](const SimulationConfig& c) { return static_cast<long long>(c.gpu.device_id); }));
  v.push_back(B("profile_steps", [](const SimulationConfig& c) { return c.profile_steps; }));
  v.push_back(R("dysbiosis_threshold", [](const SimulationConfig& c) { return c.dysbiosis_threshold; }));

  // ── Oxygen (Spec 1) ───────────────────────────────────────────────────────
  add_ns_bool(v, "oxygen.enabled", "oxygen_enabled", [](const SimulationConfig& c) { return c.chem_env.oxygen.enabled; });
  add_ns_real(v, "oxygen.epithelial_conc", "oxygen_epithelial_conc", [](const SimulationConfig& c) { return c.chem_env.oxygen.epithelial_conc; });
  add_ns_real(v, "oxygen.D_free", "oxygen_D_free", [](const SimulationConfig& c) { return c.chem_env.oxygen.D_free; });
  add_ns_real(v, "oxygen.Km", "oxygen_Km", [](const SimulationConfig& c) { return c.chem_env.oxygen.Km; });
  add_ns_real(v, "oxygen.boost_max", "oxygen_boost_max", [](const SimulationConfig& c) { return c.chem_env.oxygen.boost_max; });
  add_ns_real(v, "oxygen.q_consumption", "oxygen_q_consumption", [](const SimulationConfig& c) { return c.chem_env.oxygen.q_consumption; });
  add_ns_real(v, "oxygen.q_maintenance", "oxygen_q_maintenance", [](const SimulationConfig& c) { return c.chem_env.oxygen.q_maintenance; });
  add_ns_real(v, "oxygen.vbf_sink", "oxygen_vbf_sink", [](const SimulationConfig& c) { return c.chem_env.oxygen.vbf_sink; });
  add_ns_real(v, "oxygen.k_ROS", "oxygen_k_ROS", [](const SimulationConfig& c) { return c.chem_env.oxygen.k_ROS; });

  // ── Acetate (Spec 1) ──────────────────────────────────────────────────────
  add_ns_bool(v, "acetate.enabled", "acetate_enabled", [](const SimulationConfig& c) { return c.chem_env.acetate.enabled; });
  add_ns_real(v, "acetate.D_free", "acetate_D_free", [](const SimulationConfig& c) { return c.chem_env.acetate.D_free; });
  add_ns_real(v, "acetate.vbf_production", "acetate_vbf_production", [](const SimulationConfig& c) { return c.chem_env.acetate.vbf_production; });
  add_ns_real(v, "acetate.vbf_consumption", "acetate_vbf_consumption", [](const SimulationConfig& c) { return c.chem_env.acetate.vbf_consumption; });
  add_ns_real(v, "acetate.overflow_threshold", "acetate_overflow_threshold", [](const SimulationConfig& c) { return c.chem_env.acetate.overflow_threshold; });
  add_ns_real(v, "acetate.overflow_rate", "acetate_overflow_rate", [](const SimulationConfig& c) { return c.chem_env.acetate.overflow_rate; });
  add_ns_real(v, "acetate.scavenge_Km", "acetate_scavenge_Km", [](const SimulationConfig& c) { return c.chem_env.acetate.scavenge_Km; });
  add_ns_real(v, "acetate.scavenge_rate", "acetate_scavenge_rate", [](const SimulationConfig& c) { return c.chem_env.acetate.scavenge_rate; });
  add_ns_real(v, "acetate.epithelial_uptake", "acetate_epithelial_uptake", [](const SimulationConfig& c) { return c.chem_env.acetate.epithelial_uptake; });

  // ── Mucin (Spec 1) ────────────────────────────────────────────────────────
  add_ns_bool(v, "mucin.enabled", "mucin_enabled", [](const SimulationConfig& c) { return c.chem_env.mucin.enabled; });
  add_ns_real(v, "mucin.secretion_rate", "mucin_secretion_rate", [](const SimulationConfig& c) { return c.chem_env.mucin.secretion_rate; });
  add_ns_real(v, "mucin.Km_degradation", "mucin_Km_degradation", [](const SimulationConfig& c) { return c.chem_env.mucin.Km_degradation; });
  add_ns_real(v, "mucin.k_liberation", "mucin_k_liberation", [](const SimulationConfig& c) { return c.chem_env.mucin.k_liberation; });
  add_ns_real(v, "mucin.initial_conc", "mucin_initial_conc", [](const SimulationConfig& c) { return c.chem_env.mucin.initial_conc; });

  // ── Protease (Spec 1) ─────────────────────────────────────────────────────
  add_ns_bool(v, "protease.enabled", "protease_enabled", [](const SimulationConfig& c) { return c.chem_env.protease.enabled; });
  add_ns_real(v, "protease.default_half_life", "protease_default_half_life", [](const SimulationConfig& c) { return c.chem_env.protease.default_half_life; });
  add_ns_real(v, "protease.dilution_rate", "protease_dilution_rate", [](const SimulationConfig& c) { return c.chem_env.protease.dilution_rate; });

  // ── Siderophore (Spec 4) ───────────────────────────────────────────────────
  add_ns_bool(v, "siderophore.enabled", "siderophore_enabled", [](const SimulationConfig& c) { return c.chem_env.siderophore.enabled; });
  add_ns_real(v, "siderophore.secretion_rate", "siderophore_secretion_rate", [](const SimulationConfig& c) { return c.chem_env.siderophore.secretion_rate; });
  add_ns_real(v, "siderophore.D_free", "siderophore_D_free", [](const SimulationConfig& c) { return c.chem_env.siderophore.D_free; });
  add_ns_real(v, "siderophore.chelation_rate", "siderophore_chelation_rate", [](const SimulationConfig& c) { return c.chem_env.siderophore.chelation_rate; });
  add_ns_real(v, "siderophore.Km_reimport", "siderophore_Km_reimport", [](const SimulationConfig& c) { return c.chem_env.siderophore.Km_reimport; });
  add_ns_real(v, "siderophore.recapture_fraction", "siderophore_recapture_fraction", [](const SimulationConfig& c) { return c.chem_env.siderophore.recapture_fraction; });

  // ── Fur (Spec 3) ──────────────────────────────────────────────────────────
  add_ns_bool(v, "fur.enabled", "fur_enabled", [](const SimulationConfig& c) { return c.cell_bio.fur.enabled; });
  add_ns_real(v, "fur.Km", "fur_Km", [](const SimulationConfig& c) { return c.cell_bio.fur.Km; });
  add_ns_real(v, "fur.upregulation_max", "fur_upregulation_max", [](const SimulationConfig& c) { return c.cell_bio.fur.upregulation_max; });
  add_ns_real(v, "fur.receptor_max", "fur_receptor_max", [](const SimulationConfig& c) { return c.cell_bio.fur.receptor_max; });

  // ── CDI (Spec 3) ──────────────────────────────────────────────────────────
  add_ns_bool(v, "cdi.enabled", "cdi_enabled", [](const SimulationConfig& c) { return c.cell_bio.cdi.enabled; });
  add_ns_real(v, "cdi.kill_rate", "cdi_kill_rate", [](const SimulationConfig& c) { return c.cell_bio.cdi.kill_rate; });
  add_ns_real(v, "cdi.contact_radius", "cdi_contact_radius", [](const SimulationConfig& c) { return c.cell_bio.cdi.contact_radius; });
  add_ns_real(v, "cdi.corpse_persistence", "cdi_corpse_persistence", [](const SimulationConfig& c) { return c.cell_bio.cdi.corpse_persistence; });

  // ── Motility (Spec 3) ─────────────────────────────────────────────────────
  add_ns_bool(v, "motility.enabled", "motility_enabled", [](const SimulationConfig& c) { return c.cell_bio.motility.enabled; });
  add_ns_real(v, "motility.swim_speed", "motility_swim_speed", [](const SimulationConfig& c) { return c.cell_bio.motility.swim_speed; });
  add_ns_real(v, "motility.run_mean_duration", "motility_run_mean_duration", [](const SimulationConfig& c) { return c.cell_bio.motility.run_mean_duration; });
  add_ns_real(v, "motility.stop_probability", "motility_stop_probability", [](const SimulationConfig& c) { return c.cell_bio.motility.stop_probability; });
  add_ns_real(v, "motility.stop_duration", "motility_stop_duration", [](const SimulationConfig& c) { return c.cell_bio.motility.stop_duration; });
  add_ns_bool(v, "motility.chemotaxis_enabled", "motility_chemotaxis_enabled", [](const SimulationConfig& c) { return c.cell_bio.motility.chemotaxis_enabled; });
  add_ns_real(v, "motility.chi_carbon", "motility_chi_carbon", [](const SimulationConfig& c) { return c.cell_bio.motility.chi_carbon; });
  add_ns_real(v, "motility.chi_oxygen", "motility_chi_oxygen", [](const SimulationConfig& c) { return c.cell_bio.motility.chi_oxygen; });
  add_ns_real(v, "motility.cluster_suppress_radius", "motility_cluster_suppress_radius", [](const SimulationConfig& c) { return c.cell_bio.motility.cluster_suppress_radius; });
  add_ns_int(v, "motility.cluster_suppress_threshold", "motility_cluster_suppress_threshold", [](const SimulationConfig& c) { return static_cast<long long>(c.cell_bio.motility.cluster_suppress_threshold); });
  add_ns_real(v, "motility.cluster_tumble_factor", "motility_cluster_tumble_factor", [](const SimulationConfig& c) { return c.cell_bio.motility.cluster_tumble_factor; });

  return v;
}

// Strain-object and top-level array keys handled in config_json.cpp. These are
// verified end-to-end in test_strain_and_array_keys(); listed here so the
// completeness guard treats them as covered.
const std::set<std::string, std::less<>>& array_and_strain_keys() {
  static const std::set<std::string, std::less<>> keys = {
      "initial_strains", "fixes", "hdf5", "schedule", "grid_species",
      "type",         "count",
      "mu_max",          "plasmids", "conjugative", "cdi_type",
      "cdi_immunity"};
  return keys;
}

// Relative tolerance is 1e-4 (0.01%): tight enough to catch a value routed to
// the wrong field, but loose enough to absorb the JSON scalar path's ~6
// significant-digit string round-trip (ConfigJson::apply_json_scalar reformats
// numbers with default stream precision). The distinctness check still holds
// because sentinels differ from defaults by a factor of ~3.
bool close(double a, double b) {
  const double scale = std::max(1.0, std::fabs(b));
  return std::fabs(a - b) <= 1e-4 * scale;
}

std::string fmt_real(double v) {
  std::ostringstream o;
  o << std::setprecision(17) << v;
  return o.str();
}

// Scale-robust sentinel: distinct from the default even for large magnitudes
// (default + 1.0 would be lost to rounding for values like 1e11).
double real_sentinel(double def) { return def * 3.0 + 7.0; }

// Scalar string fed to the parser for a probe (unquoted).
std::string sentinel_scalar(const Probe& p, const SimulationConfig& def) {
  switch (p.kind) {
    case Kind::Real:
      return fmt_real(real_sentinel(p.gr(def)));
    case Kind::Int:
      return std::to_string(p.gi(def) + 7);
    case Kind::Bool:
      return (!p.gb(def)) ? "true" : "false";
    case Kind::Str:
      return p.sval;
  }
  return "";
}

// Assert the parsed config `got` holds the sentinel for probe `p`, and that the
// sentinel actually differs from the default (guards against dead wiring).
void check_ingested(const Probe& p, const SimulationConfig& def, const SimulationConfig& got) {
  switch (p.kind) {
    case Kind::Real: {
      if (const double s = real_sentinel(p.gr(def));
          !close(p.gr(got), s) || close(p.gr(got), p.gr(def))) {
        std::ostringstream m;
        m << "real key '" << p.key << "' got=" << p.gr(got) << " expected=" << s
          << " default=" << p.gr(def) << " (did not ingest into SimulationConfig)";
        record_failure(m.str());
      }
      break;
    }
    case Kind::Int: {
      if (const long long s = p.gi(def) + 7; p.gi(got) != s || p.gi(got) == p.gi(def)) {
        std::ostringstream m;
        m << "int key '" << p.key << "' got=" << p.gi(got) << " expected=" << s
          << " default=" << p.gi(def) << " (did not ingest into SimulationConfig)";
        record_failure(m.str());
      }
      break;
    }
    case Kind::Bool: {
      if (const bool s = !p.gb(def); p.gb(got) != s) {
        std::ostringstream m;
        m << "bool key '" << p.key << "' got=" << p.gb(got) << " expected=" << s
          << " (did not ingest into SimulationConfig)";
        record_failure(m.str());
      }
      break;
    }
    case Kind::Str: {
      if (p.gs(got) != p.sval || p.gs(got) == p.gs(def)) {
        std::ostringstream m;
        m << "string key '" << p.key << "' got='" << p.gs(got) << "' expected='"
          << p.sval << "' (did not ingest into SimulationConfig)";
        record_failure(m.str());
      }
      break;
    }
  }
}

std::string read_file(const std::string& path) {
  std::ifstream in(path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// Extract every `key == "LITERAL"` string literal from parser source. Value
// comparisons in the parser use `val == ...`, so this pattern captures only
// config keys.
void scan_source_keys(const std::string& path, std::set<std::string, std::less<>>& out) {
  const std::string src = read_file(path);
  const std::string tok = "key == \"";
  size_t pos = 0;
  while ((pos = src.find(tok, pos)) != std::string::npos) {
    const size_t start = pos + tok.size();
    const size_t end = src.find('"', start);
    if (end == std::string::npos) break;
    out.insert(src.substr(start, end - start));
    pos = end + 1;
  }
}

std::string source_path(const char* rel) {
  return std::string(GUTIBM_SOURCE_DIR) + "/" + rel;
}

}  // namespace

// Every flat key ingests via InputParser::apply_flat_key (legacy + JSON scalar
// entry point), landing in the right SimulationConfig field.
void test_every_flat_key_ingests_via_apply() {
  const std::vector<Probe> probes = build_probes();
  const SimulationConfig def = InputParser::default_config();

  for (const Probe& p : probes) {
    SimulationConfig cfg = def;
    InputParser::apply_flat_key(cfg, p.key, sentinel_scalar(p, def));
    check_ingested(p, def, cfg);
  }

  std::cout << "  test_every_flat_key_ingests_via_apply: PASSED (" << probes.size()
            << " keys)\n";
}

// Every flat key survives the full JSON document path (InputParser::parse ->
// ConfigJson::parse_document -> apply_json_scalar -> apply_flat_key).
void test_every_flat_key_ingests_via_json_document() {
  const std::vector<Probe> probes = build_probes();
  const SimulationConfig def = InputParser::default_config();

  std::ostringstream doc;
  doc << "{\n  \"_comment\": \"exhaustive config ingestion probe\"";
  for (const Probe& p : probes) {
    if (!p.primary) continue;
    doc << ",\n  \"" << p.key << "\": ";
    if (p.kind == Kind::Str) {
      doc << '"' << sentinel_scalar(p, def) << '"';
    } else {
      doc << sentinel_scalar(p, def);
    }
  }
  doc << "\n}\n";

  const std::string path = source_path("tests/fixtures/_config_ingestion_probe.json");
  {
    std::ofstream out(path);
    out << doc.str();
  }

  const SimulationConfig got = InputParser::parse(path);
  std::remove(path.c_str());

  int checked = 0;
  for (const Probe& p : probes) {
    if (!p.primary) continue;
    check_ingested(p, def, got);
    ++checked;
  }

  std::cout << "  test_every_flat_key_ingests_via_json_document: PASSED (" << checked
            << " keys)\n";
}

// Strain-object fields and the top-level initial_strains / fixes arrays ingest.
void test_strain_and_array_keys() {
  const std::string json = R"({
    "_comment": "strain + array key ingestion probe",
    "initial_strains": [
      {"type": 3, "count": 9, "mu_max": 7e-4, "plasmids": ["ColE1", "ColB"],
       "conjugative": true, "cdi_type": 5, "cdi_immunity": 6}
    ],
    "fixes": ["metabolism", "mechanics"]
  })";

  const std::string path = source_path("tests/fixtures/_config_ingestion_strains.json");
  {
    std::ofstream out(path);
    out << json;
  }

  const SimulationConfig cfg = InputParser::parse(path);
  std::remove(path.c_str());

  expect(cfg.initial_strains.size() == 1, "initial_strains should parse one strain");
  if (cfg.initial_strains.size() == 1) {
    const auto& s = cfg.initial_strains[0];
    expect(s.type == 3, "strain 'type' not ingested");
    expect(s.count == 9, "strain 'count' not ingested");
    expect(std::abs(s.mu_max - 7e-4) < 1e-12, "strain 'mu_max' not ingested");
    expect(s.plasmids.size() == 2 && s.plasmids[0] == "ColE1" && s.plasmids[1] == "ColB",
           "strain 'plasmids' not ingested");
    expect(s.conjugative, "strain 'conjugative' not ingested");
    expect(s.cdi_type == 5, "strain 'cdi_type' not ingested");
    expect(s.cdi_immunity == 6, "strain 'cdi_immunity' not ingested");
  }

  expect(cfg.enabled_fixes.size() == 2, "'fixes' array not ingested");
  if (cfg.enabled_fixes.size() == 2) {
    expect(cfg.enabled_fixes[0] == "metabolism" && cfg.enabled_fixes[1] == "mechanics",
           "'fixes' array contents not ingested");
  }

  std::cout << "  test_strain_and_array_keys: checked strain + array keys\n";
}

// Completeness guard: the set of keys parsed by the sources must equal the set
// of keys tracked here. Adding a parser key without a probe (or removing a
// probe for a still-parsed key) fails CI.
void test_all_parser_keys_are_tracked() {
  std::set<std::string, std::less<>> source_keys;
  scan_source_keys(source_path("src/io/input_parser.cpp"), source_keys);
  scan_source_keys(source_path("src/io/config_json.cpp"), source_keys);
  expect(!source_keys.empty(), "failed to scan parser sources for config keys");

  std::set<std::string, std::less<>> covered;
  for (const Probe& p : build_probes()) covered.insert(p.key);
  for (const std::string& k : array_and_strain_keys()) covered.insert(k);

  std::vector<std::string> untracked;   // parsed but no probe
  std::vector<std::string> stale;       // probe but no longer parsed
  for (const std::string& k : source_keys) {
    if (!covered.contains(k)) untracked.push_back(k);
  }
  for (const std::string& k : covered) {
    if (!source_keys.contains(k)) stale.push_back(k);
  }

  if (!untracked.empty()) {
    std::cerr << "ERROR: config keys parsed but NOT tracked by an ingestion probe:\n";
    for (const std::string& k : untracked) std::cerr << "  - " << k << "\n";
    std::cerr << "Add a probe in build_probes() (or array_and_strain_keys()) in "
                 "tests/test_config_ingestion.cpp.\n";
  }
  if (!stale.empty()) {
    std::cerr << "ERROR: ingestion probes for keys the parser no longer recognizes:\n";
    for (const std::string& k : stale) std::cerr << "  - " << k << "\n";
  }
  expect(untracked.empty(), "every parsed config key must have an ingestion probe");
  expect(stale.empty(), "remove ingestion probes for keys no longer parsed");

  std::cout << "  test_all_parser_keys_are_tracked: " << source_keys.size()
            << " parser keys checked against tracked probes\n";
}

int main() {
  std::cout << "=== Config Ingestion Tracking Tests ===\n";
  test_every_flat_key_ingests_via_apply();
  test_every_flat_key_ingests_via_json_document();
  test_strain_and_array_keys();
  test_all_parser_keys_are_tracked();
  if (failure_counter().value != 0) {
    std::cerr << failure_counter().value << " config ingestion check(s) FAILED.\n";
    return 1;
  }
  std::cout << "All config ingestion tracking tests passed.\n";
  return 0;
}
