/* -----------------------------------------------------------------------
   GutIBM – Shared founder-genotype test helpers
   ----------------------------------------------------------------------- */

#pragma once

#include "agent.h"
#include "error.h"
#include "input_parser.h"
#include "plasmid.h"

#include <cassert>
#include <cmath>
#include <initializer_list>
#include <string>

namespace gutibm::test {

// The declared resistant founder used across the layer tests:
// receptor_expression {"BtuB": 0.0} -> RESISTANT, cdi 7/3, ColE1.
inline void declare_resistant_founder(SimulationConfig::InitialStrain& s) {
  s.receptor_expression = {{"BtuB", 0.0}};
  s.cdi_type = 7;
  s.cdi_immunity = 3;
}

// The retardation the configured_plasmid_cluster path derives for the
// founder's ColE1 cluster under this config.
inline Real declared_cluster_retardation(const SimulationConfig& cfg) {
  return retardation_from_pI(PlasmidLibrary::find("ColE1")->cluster.pI,
                             cfg.fixes.bacteriocin.mucin_charge);
}

inline void expect_declared_r_phenotype(const Agent& a, Int btub) {
  assert(a.state == PhenoState::RESISTANT);
  assert(a.receptor_expr[btub] == 0.0);
  assert(a.receptor_expr_base[btub] == 0.0);
  assert(a.genome.receptor_expression[btub] == 0.0);
  assert(a.genome.cdi_type == 7);
  assert(a.genome.cdi_immunity == 3);
}

inline void expect_configured_cluster(const Agent& a, Real expected_ret) {
  assert(a.genome.bi_loci.size() == 1);
  assert(std::abs(a.genome.bi_loci[0].retardation - expected_ret) < 1e-9);
}

inline bool carries_declared_r_genotype(const Agent& a, Int btub) {
  return a.genome.bi_loci.size() == 1 &&
         a.state == PhenoState::RESISTANT && a.receptor_expr[btub] == 0.0 &&
         a.genome.receptor_expression[btub] == 0.0 &&
         a.genome.cdi_type == 7 && a.genome.cdi_immunity == 3;
}

// ── Capability-gate sweep machinery ───────────────────────────────────
// init_with(cfg) must throw ConfigError naming the dead key.
struct DeadKnobCase {
  const char* name;
  const char* expected_key;
  void (*set)(SimulationConfig&);
};

// Shared dead-knob sweep axes refused by every tabulated layer.
inline void knob_sos_lysis(SimulationConfig& c) {
  c.fixes.bacteriocin.sos_lysis_prob = 0.05;
}
inline void knob_retardation_override(SimulationConfig& c) {
  c.plasmid_overrides["ColE1"].retardation = 1.2;
}
inline void knob_fix_list(SimulationConfig& c) {
  c.enabled_fixes = {"bacteriocin"};
}

template <typename InitFn>
inline int expect_dead_knobs_refused(
    std::initializer_list<DeadKnobCase> cases, SimulationConfig base,
    InitFn init_with) {
  int refused = 0;
  for (const DeadKnobCase& tc : cases) {
    SimulationConfig cfg = base;
    tc.set(cfg);
    bool threw = false;
    try {
      init_with(cfg);
    } catch (const ConfigError& e) {
      threw = true;
      assert(std::string(e.what()).find(tc.expected_key) !=
             std::string::npos);
    }
    assert(threw);
    ++refused;
  }
  return refused;
}

}  // namespace gutibm::test
