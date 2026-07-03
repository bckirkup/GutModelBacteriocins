/* -----------------------------------------------------------------------
   GutIBM – HDF5 snapshot reader
   Mirrors the dataset layout produced by HDF5Writer.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_HDF5_READER_H
#define GUTIBM_HDF5_READER_H

#include "types.h"
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace gutibm {

struct HDF5CheckpointMetadata {
  Int   step         = 0;
  Real  time         = 0.0;
  Int   num_agents   = 0;
  Int   num_lineages = 0;
};

struct HDF5CheckpointAgents {
  std::vector<int64_t> id;
  std::vector<int32_t> type;
  std::vector<int32_t> state;
  std::vector<double>  x;
  std::vector<double>  y;
  std::vector<double>  z;
  std::vector<double>  radius;
  std::vector<double>  biomass;
  std::vector<double>  mu;
  std::vector<int64_t> lineage;
};

struct HDF5CheckpointLineage {
  std::vector<double>  btuB_expression;
  std::vector<double>  fepA_expression;
  std::vector<int32_t> num_bi_loci;
  std::vector<int32_t> generation;
};

struct HDF5CheckpointGenome {
  bool present = false;

  std::vector<int64_t> parent_id;
  std::vector<int32_t> mutations;
  std::vector<int32_t> has_conjugative_plasmid;
  std::vector<double>  plasmid_cost_amelioration;
  std::vector<int32_t> cdi_type;
  std::vector<int32_t> cdi_immunity;
  std::vector<double>  receptor_expression;  // N * NUM_RECEPTORS, row-major
  std::vector<double>  toxin_affinity;
  std::vector<double>  ligand_affinity;

  // Flat BI locus arrays (length = sum of num_bi_loci)
  std::vector<int32_t> bi_toxin_id;
  std::vector<int32_t> bi_immunity_id;
  std::vector<int32_t> bi_target;
  std::vector<int32_t> bi_bclass;
  std::vector<double>  bi_pI;
  std::vector<double>  bi_diff_coeff;
  std::vector<double>  bi_retardation;
  std::vector<double>  bi_molecular_weight;
  std::vector<double>  bi_immunity_binding_affinity;
};

struct HDF5CheckpointGrid {
  std::map<std::string, std::vector<double>, std::less<>> species;
};

struct HDF5CheckpointSnapshot {
  std::string              step_name;
  HDF5CheckpointMetadata   metadata;
  HDF5CheckpointAgents     agents;
  HDF5CheckpointLineage    lineage;
  HDF5CheckpointGenome     genome;
  HDF5CheckpointGrid       grid;
};

class HDF5Reader {
 public:
  HDF5Reader() = default;
  ~HDF5Reader();

  HDF5Reader(const HDF5Reader&) = delete;
  HDF5Reader& operator=(const HDF5Reader&) = delete;

  bool open(const std::string& filename);
  void close();

  bool is_open() const { return open_; }

  std::vector<std::string> list_steps() const;
  std::string latest_step() const;

  HDF5CheckpointSnapshot load_step(const std::string& step_name) const;

  // Convenience: open, load latest (or named) step, close.
  static HDF5CheckpointSnapshot load_snapshot(
      const std::string& filename,
      const std::string& step_name = "");

 private:
  bool open_ = false;
  std::string filename_;

#ifdef GUTIBM_HDF5
  int64_t file_id_ = -1;
#endif
};

}  // namespace gutibm

#endif  // GUTIBM_HDF5_READER_H
