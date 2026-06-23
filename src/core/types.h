/* -----------------------------------------------------------------------
   GutIBM - 3D Individual-based Model for Enterobacteriaceae Gut Dynamics
   Inspired by NUFEB-2 / LAMMPS architecture
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_TYPES_H
#define GUTIBM_TYPES_H

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace gutibm {

using Real    = double;
using Int     = int;
using TagID   = int64_t;
using Vec3    = std::array<Real, 3>;

static constexpr Real PI        = 3.14159265358979323846;
static constexpr Real BOLTZMANN = 1.380649e-23;    // J/K
static constexpr Real AVOGADRO  = 6.02214076e23;

// Biological constants (SI)
static constexpr Real CELL_RADIUS_DEFAULT  = 0.5e-6;   // 0.5 um
static constexpr Real CELL_DENSITY_DEFAULT = 1100.0;    // kg/m^3 (slightly > water)
static constexpr Real GUT_TEMPERATURE      = 310.15;    // 37 C in Kelvin
static constexpr Real GUT_PH              = 6.8;        // colonic pH

// Receptor types (TonB-Dependent Transporters)
enum class ReceptorType : uint8_t {
  BtuB  = 0,   // vitamin B12 / colicin E1,E2,E5,E7,E8,E9,K,Ia,Ib
  FepA  = 1,   // enterobactin-Fe / colicin B,D
  Tsx   = 2,   // nucleosides / colicin K  (secondary)
  FhuA  = 3,   // ferrichrome / colicin M, phage T5
  IroN  = 4,   // salmochelin (secondary enterobactin)
  Fiu   = 5,   // catechol siderophores
  CirA  = 6,   // linearized enterobactin / colicin Ia
  IutA  = 7,   // aerobactin (secondary iron uptake)
  NUM_RECEPTOR_TYPES
};

static constexpr int NUM_RECEPTORS = static_cast<int>(ReceptorType::NUM_RECEPTOR_TYPES);

// Bacteriocin classification by isoelectric point
enum class BacteriocinClass : uint8_t {
  LETHAL_CORE = 0,   // pI > 8.5 → binds mucin → concentrated near producer
  LETHAL_HALO = 1,   // pI < 7.0 → repelled by mucin → wider diffuse halo
  NEUTRAL     = 2    // 7.0 <= pI <= 8.5
};

// Agent phenotype states
enum class PhenoState : uint8_t {
  NORMAL       = 0,
  RESISTANT    = 1,   // downregulated receptor(s) → metabolic penalty
  SOS_INDUCED  = 2,   // about to lyse
  DEAD         = 3
};

// Plasmid / BI-locus representation
struct BICluster {
  uint16_t toxin_id;        // bacteriocin identity
  uint16_t immunity_id;     // cognate immunity protein
  ReceptorType target;      // which receptor the toxin hijacks
  BacteriocinClass bclass;  // lethal core vs halo
  Real pI;                  // isoelectric point
  Real diff_coeff;          // free diffusion coefficient (m^2/s)
  Real retardation;         // mucin retardation factor (1 = no retardation)
  Real molecular_weight;    // Da
  Real immunity_binding_affinity = 1.0;  // 1.0 = full cognate protection, 0.0 = none
};

// Agent genome (compact representation for 10^7 agents)
struct Genome {
  TagID lineage_id;                  // ancestral lineage tag
  TagID parent_id;                   // immediate parent
  uint32_t generation;               // division count from ancestor
  std::vector<BICluster> bi_loci;    // bacteriocin-immunity clusters
  std::array<Real, NUM_RECEPTORS> receptor_expression;  // 0.0–1.0 expression
  bool has_conjugative_plasmid;      // can initiate HGT
  uint32_t mutations;                // accumulated mutation count
  Real plasmid_cost_amelioration;    // compensatory mutations reduce per-locus cost
};

inline Real sphere_volume(Real radius) {
  return (4.0 / 3.0) * PI * radius * radius * radius;
}

inline Real sphere_mass(Real radius, Real density) {
  return sphere_volume(radius) * density;
}

}  // namespace gutibm

#endif  // GUTIBM_TYPES_H
