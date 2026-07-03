/* -----------------------------------------------------------------------
   GutIBM – Receptor classification helpers
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_RECEPTOR_UTILS_H
#define GUTIBM_RECEPTOR_UTILS_H

#include "types.h"

namespace gutibm {

inline bool is_iron_receptor(int r) {
  using enum ReceptorType;
  return r == to_underlying(FepA)
      || r == to_underlying(FhuA)
      || r == to_underlying(IroN)
      || r == to_underlying(IutA)
      || r == to_underlying(Fiu)
      || r == to_underlying(CirA);
}

}  // namespace gutibm

#endif  // GUTIBM_RECEPTOR_UTILS_H
