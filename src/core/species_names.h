/* -----------------------------------------------------------------------
   GutIBM – Canonical chemical species names

   Single source of truth for the string identifiers used to register and
   look up chemical species in the ChemicalField. Using these constants
   instead of scattered string literals prevents silent typo bugs (a
   mistyped name returns index -1 and the coupled physics is silently
   skipped) and makes it obvious where a new species must be wired in.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_SPECIES_NAMES_H
#define GUTIBM_SPECIES_NAMES_H

#include "types.h"
#include <array>

namespace gutibm::species {

inline constexpr const char* CARBON            = "carbon";
inline constexpr const char* IRON              = "iron";
inline constexpr const char* B12               = "b12";
inline constexpr const char* BACTERIOCIN_BTUB  = "bacteriocin_BtuB";
inline constexpr const char* BACTERIOCIN_FEPA  = "bacteriocin_FepA";
inline constexpr const char* BACTERIOCIN_CIRA  = "bacteriocin_CirA";
inline constexpr const char* BACTERIOCIN_FHUA  = "bacteriocin_FhuA";
inline constexpr const char* ACETATE            = "acetate";
inline constexpr const char* ETHANOLAMINE       = "ethanolamine";
inline constexpr const char* OXYGEN             = "oxygen";
inline constexpr const char* MUCIN              = "mucin";
inline constexpr const char* SIDEROPHORE         = "siderophore";
inline constexpr const char* AI2                 = "ai2";

// Per-receptor toxin field for a TonB-dependent transporter target.
inline const char* bacteriocin_species_for(ReceptorType target) {
  using enum ReceptorType;
  switch (target) {
    case BtuB: return BACTERIOCIN_BTUB;
    case FepA: return BACTERIOCIN_FEPA;
    case CirA: return BACTERIOCIN_CIRA;
    case FhuA: return BACTERIOCIN_FHUA;
    default: return nullptr;
  }
}

inline constexpr std::array<ReceptorType, 4> BACTERIOCIN_RECEPTOR_TARGETS = {
    ReceptorType::BtuB, ReceptorType::FepA, ReceptorType::CirA, ReceptorType::FhuA};

}  // namespace gutibm::species

#endif  // GUTIBM_SPECIES_NAMES_H
