/* -----------------------------------------------------------------------
   GutIBM – Plasmid library implementation
   ----------------------------------------------------------------------- */

#include "plasmid.h"
#include <unordered_map>

namespace gutibm {

BacteriocinClass classify_by_pI(Real pI) {
  if (pI > 8.5) return BacteriocinClass::LETHAL_CORE;
  if (pI < 7.0) return BacteriocinClass::LETHAL_HALO;
  return BacteriocinClass::NEUTRAL;
}

BICluster PlasmidLibrary::colicin_E1() {
  return BICluster{
    .toxin_id        = 1,
    .immunity_id     = 1,
    .target          = ReceptorType::BtuB,
    .bclass          = classify_by_pI(9.0),
    .pI              = 9.0,
    .diff_coeff      = 4.0e-11,    // ~50 kDa protein
    .retardation     = 50.0,       // basic → binds mucin strongly
    .molecular_weight = 57000.0,
    .protease_half_life = 1800.0
  };
}

BICluster PlasmidLibrary::colicin_E2() {
  // Note: secreted as equimolar complex with Im2 (acidic)
  // Net complex pI < 7.0 → actually a Lethal Halo
  return BICluster{
    .toxin_id        = 2,
    .immunity_id     = 2,
    .target          = ReceptorType::BtuB,
    .bclass          = classify_by_pI(6.5),
    .pI              = 6.5,        // complex pI (corrected per VADI review)
    .diff_coeff      = 3.5e-11,
    .retardation     = 3.0,        // modest retardation
    .molecular_weight = 61500.0,    // toxin + immunity complex
    .protease_half_life = 1800.0
  };
}

BICluster PlasmidLibrary::colicin_B() {
  return BICluster{
    .toxin_id        = 3,
    .immunity_id     = 3,
    .target          = ReceptorType::FepA,
    .bclass          = classify_by_pI(5.4),
    .pI              = 5.4,
    .diff_coeff      = 4.0e-11,
    .retardation     = 1.5,        // acidic → repelled by mucin
    .molecular_weight = 54800.0,
    .protease_half_life = 900.0
  };
}

BICluster PlasmidLibrary::colicin_Ia() {
  return BICluster{
    .toxin_id        = 4,
    .immunity_id     = 4,
    .target          = ReceptorType::CirA,
    .bclass          = classify_by_pI(7.2),
    .pI              = 7.2,
    .diff_coeff      = 4.0e-11,
    .retardation     = 5.0,
    .molecular_weight = 69400.0,
    .protease_half_life = 2400.0
  };
}

BICluster PlasmidLibrary::colicin_M() {
  return BICluster{
    .toxin_id        = 5,
    .immunity_id     = 5,
    .target          = ReceptorType::FhuA,
    .bclass          = classify_by_pI(9.3),
    .pI              = 9.3,
    .diff_coeff      = 5.0e-11,
    .retardation     = 60.0,
    .molecular_weight = 29500.0,
    .protease_half_life = 900.0
  };
}

BICluster PlasmidLibrary::microcin_V() {
  return BICluster{
    .toxin_id        = 10,
    .immunity_id     = 10,
    .target          = ReceptorType::CirA,
    .bclass          = classify_by_pI(5.0),
    .pI              = 5.0,
    .diff_coeff      = 1.0e-10,    // small peptide → faster diffusion
    .retardation     = 1.2,
    .molecular_weight = 8900.0,
    .protease_half_life = 7200.0
  };
}

const std::vector<PlasmidEntry>& PlasmidLibrary::entries() {
  static std::vector<PlasmidEntry> lib = {
    {"ColE1", colicin_E1(), false,
     "Pore-forming colicin, BtuB receptor, basic pI → Lethal Core"},
    {"ColE2", colicin_E2(), false,
     "Nuclease colicin, BtuB receptor, secreted with Im2 → Lethal Halo (corrected)"},
    {"ColB",  colicin_B(), true,
     "Pore-forming colicin, FepA receptor, acidic pI → Lethal Halo, conjugative"},
    {"ColIa", colicin_Ia(), true,
     "Pore-forming colicin, CirA receptor, neutral pI, conjugative"},
    {"ColM",  colicin_M(), false,
     "Murein inhibitor, FhuA receptor, basic pI → Lethal Core"},
    {"MccV",  microcin_V(), true,
     "Small peptide microcin, continuous secretion, CirA receptor"},
  };
  return lib;
}

const PlasmidEntry* PlasmidLibrary::find(const std::string& name) {
  static const std::unordered_map<std::string, std::string> aliases = {
    {"colicin_E1", "ColE1"},
    {"colicin_E2", "ColE2"},
    {"colicin_B",  "ColB"},
    {"colicin_Ia", "ColIa"},
    {"colicin_M",  "ColM"},
    {"microcin_V", "MccV"},
  };

  std::string key = name;
  if (auto it = aliases.find(name); it != aliases.end()) key = it->second;

  for (const auto& entry : entries()) {
    if (entry.name == key) return &entry;
  }
  return nullptr;
}

}  // namespace gutibm
