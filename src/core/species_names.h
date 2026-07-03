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

namespace gutibm::species {

inline constexpr const char* CARBON         = "carbon";
inline constexpr const char* IRON           = "iron";
inline constexpr const char* B12            = "b12";
inline constexpr const char* BACTERIOCIN    = "bacteriocin";
inline constexpr const char* NUCLEASE_TOXIN = "nuclease_toxin";
inline constexpr const char* ACETATE        = "acetate";
inline constexpr const char* ETHANOLAMINE   = "ethanolamine";
inline constexpr const char* OXYGEN         = "oxygen";
inline constexpr const char* MUCIN          = "mucin";

}  // namespace gutibm::species

#endif  // GUTIBM_SPECIES_NAMES_H
