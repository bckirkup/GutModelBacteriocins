/* -----------------------------------------------------------------------
   GutIBM – HDF5 output writer
   Compatible with nufeb_tools analysis framework.
   
   Output structure:
     /step_NNNN/
       atoms/          (agent data)
         id, type, x, y, z, radius, biomass, mu, state, lineage
       grid/           (chemical fields)
         carbon, iron, b12, bacteriocin
       metadata/
         time, num_agents, num_lineages
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_HDF5_WRITER_H
#define GUTIBM_HDF5_WRITER_H

#include "types.h"
#include <string>

namespace gutibm {

class Simulation;

struct HDF5Config {
  std::string filename   = "gut_ibm_output.h5";
  Int         dump_every = 100;   // steps between dumps
  bool        parallel   = false; // MPI-parallel I/O
};

class HDF5Writer {
 public:
  HDF5Writer() = default;

  void init(const HDF5Config& cfg);
  void write_step(const Simulation& sim, Int step, Real time);
  void finalize();

  bool is_enabled() const { return enabled_; }

 private:
  void write_agents(const Simulation& sim, const std::string& group);
  void write_grid(const Simulation& sim, const std::string& group);
  void write_metadata(const Simulation& sim, const std::string& group,
                       Int step, Real time);
  void write_lineage(const Simulation& sim, const std::string& group);

  HDF5Config cfg_;
  bool enabled_ = false;

#ifdef GUTIBM_HDF5
  // HDF5 file handle stored as int64 to avoid including hdf5.h in header
  int64_t file_id_ = -1;
#endif
};

}  // namespace gutibm

#endif  // GUTIBM_HDF5_WRITER_H
