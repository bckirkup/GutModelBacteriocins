/* -----------------------------------------------------------------------
   GutIBM – Contact-dependent inhibition (CDI) fix (Spec 3)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_CDI_H
#define GUTIBM_FIX_CDI_H

#include "fix.h"
#include "cdi_config.h"

namespace gutibm {

class FixCdi : public Fix {
 public:
  FixCdi(Simulation& sim, const CdiConfig& cfg);

  void compute(Real dt) override;

 private:
  CdiConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_CDI_H
