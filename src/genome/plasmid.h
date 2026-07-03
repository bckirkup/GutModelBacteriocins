/* -----------------------------------------------------------------------
   GutIBM – Plasmid and BI locus database
   Predefined bacteriocin-immunity clusters based on E. coli genetics.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_PLASMID_H
#define GUTIBM_PLASMID_H

#include "types.h"
#include <vector>
#include <string>

namespace gutibm {

struct PlasmidEntry {
  std::string name;
  BICluster cluster;
  bool conjugative;
  std::string description;
};

// Classify bacteriocin by isoelectric point.
// Thresholds: pI > 8.5 → CORE, pI < 7.0 → HALO, else NEUTRAL (VADI §74).
// SINGLE SOURCE OF TRUTH for pI → BacteriocinClass. All other code must use
// the pre-computed BICluster.bclass rather than re-classifying from pI.
BacteriocinClass classify_by_pI(Real pI);

// Library of well-characterized E. coli bacteriocin systems
class PlasmidLibrary {
 public:
  static const std::vector<PlasmidEntry>& entries();

  // Resolve a plasmid by canonical name or legacy alias (e.g. colicin_E1 → ColE1).
  // Returns nullptr if not found.
  static const PlasmidEntry* find(const std::string& name);

  // Create a default BI cluster for a named colicin
  static BICluster colicin_E1();   // pI ~9.0, pore-forming, BtuB
  static BICluster colicin_E2();   // nuclease, BtuB (secreted as acidic complex)
  static BICluster colicin_B();    // pI ~5.4, pore-forming, FepA
  static BICluster colicin_Ia();   // pore-forming, CirA
  static BICluster colicin_M();    // murein synthesis inhibitor, FhuA
  static BICluster microcin_V();   // small peptide, continuous secretion
};

}  // namespace gutibm

#endif  // GUTIBM_PLASMID_H
