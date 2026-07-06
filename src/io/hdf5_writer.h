/* -----------------------------------------------------------------------
   GutIBM – HDF5 output writer (Spec 4 layered schema)
   
   Output structure:
     /file attrs: gutibm_version, nx, ny, nz, grid_dx, domain_lo/hi
     /summary/step_NNNNNN/   — scalars + event counters + chem summaries
     /agents/step_NNNNNN/    — lightweight per-agent arrays
     /grid/step_NNNNNN/      — 3D species datasets (optional gzip)
     /lineage/step_NNNNNN/   — lineage tracker arrays
     /genome/step_NNNNNN/    — full genome / BI locus tables
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_HDF5_WRITER_H
#define GUTIBM_HDF5_WRITER_H

#include "types.h"
#include <string>
#include <vector>

namespace gutibm {

class Simulation;

struct HDF5Schedule {
  Int summary = 1;
  Int agents  = 5;
  Int grid    = 0;   // 0 = disabled
  Int lineage = 100;
  Int genome  = 100;
  std::vector<std::string> grid_species;  // empty = none; ["all"] = all species
};

struct HDF5Config {
  std::string filename      = "gut_ibm_output.h5";
  bool        enabled       = true;
  bool        parallel      = false;  // reserved
  HDF5Schedule schedule;
  std::string compression   = "none";  // "none" | "gzip"
  Int         compression_level = 4;
};

class HDF5Writer {
 public:
  HDF5Writer() = default;

  void init(const HDF5Config& cfg, const class Domain& domain);
  void write_step(Simulation& sim, Int step, Real time, Real dt) const;
  void finalize();

  bool is_enabled() const { return enabled_; }

 private:
  bool layer_due(Int interval, Int step) const;
  bool should_write_species(const std::string& name) const;

  void write_summary(Simulation& sim, const std::string& step_group,
                     Int step, Real time, Real dt) const;
  void write_agents_layer(const Simulation& sim, const std::string& step_group) const;
  void write_grid_layer(const Simulation& sim, const std::string& step_group) const;
  void write_lineage_layer(const Simulation& sim, const std::string& step_group) const;
  void write_genome_layer(const Simulation& sim, const std::string& step_group) const;

  HDF5Config cfg_;
  bool enabled_ = false;
  Int nx_ = 0;
  Int ny_ = 0;
  Int nz_ = 0;
  Real grid_dx_ = 0.0;
  Vec3 domain_lo_{};
  Vec3 domain_hi_{};

#ifdef GUTIBM_HDF5
  int64_t file_id_ = -1;
#endif
};

}  // namespace gutibm

#endif  // GUTIBM_HDF5_WRITER_H
