/* -----------------------------------------------------------------------
   GutIBM – Base fix interface
   Follows NUFEB/LAMMPS pattern: each Fix governs one aspect of
   agent behavior and is called at specific points in the timestep.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_H
#define GUTIBM_FIX_H

#include "types.h"
#include <string>

namespace gutibm {

class Simulation;

class Fix {
 public:
  Fix(const std::string& name, Simulation& sim) : name_(name), sim_(sim) {}
  virtual ~Fix() = default;

  // Called once at simulation start
  virtual void init() {}

  // Called at the start of each biological timestep
  virtual void pre_step(Real dt) {}

  // Main computation
  virtual void compute(Real dt) = 0;

  // Called at end of each biological timestep (post-processing)
  virtual void post_step(Real dt) {}

  const std::string& name() const { return name_; }

 protected:
  std::string name_;
  Simulation& sim_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_H
