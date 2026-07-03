/* -----------------------------------------------------------------------
   GutIBM – Receptor classification helpers
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_RECEPTOR_UTILS_H
#define GUTIBM_RECEPTOR_UTILS_H

#include "types.h"

namespace gutibm {

inline bool is_iron_receptor(int r) {
  return r == to_underlying(ReceptorType::FepA)
      || r == to_underlying(ReceptorType::FhuA)
      || r == to_underlying(ReceptorType::IroN)
      || r == to_underlying(ReceptorType::IutA)
      || r == to_underlying(ReceptorType::Fiu)
      || r == to_underlying(ReceptorType::CirA);
}

}  // namespace gutibm

#endif  // GUTIBM_RECEPTOR_UTILS_H
