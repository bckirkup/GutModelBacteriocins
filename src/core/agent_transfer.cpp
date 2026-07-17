/* -----------------------------------------------------------------------
   GutIBM – MPI agent serialization implementation
   ----------------------------------------------------------------------- */

#include "agent_transfer.h"
#include "types.h"
#include <cstddef>
#include <cstring>
#include <span>
#include <vector>

namespace gutibm {

namespace {

struct AgentIdentityXfer {
  int64_t tag;
  int32_t type;
  int32_t owner_rank;
};

struct AgentSpatialXfer {
  std::array<double, 3> x;
  std::array<double, 3> v;
  double radius;
  double mass;
  double outer_radius;
};

struct AgentMetabolismXfer {
  double mu_max;
  double mu_realized;
  double biomass;
  double maintenance;
  std::array<double, NUM_RECEPTORS> receptor_expr;
  std::array<double, NUM_RECEPTORS> receptor_expr_base;
  double km_iron;
  double km_b12;
  double km_carbon;
};

struct AgentLifecycleXfer {
  int32_t state;
  double age;
  double sos_timer;
  double death_time;
  int32_t grid_cell;
  int32_t in_crypt;
};

struct AgentMotilityXfer {
  std::array<double, 3> swim_direction;
  double swim_speed;
  double run_timer;
  int32_t is_stopped;
  double stop_timer;
  double prev_carbon;
  double prev_oxygen;
  double prev_ai2;
};

struct AgentGenomeXfer {
  int64_t lineage_id;
  int64_t parent_id;
  uint32_t generation;
  uint32_t mutations;
  int32_t has_conjugative_plasmid;
  double plasmid_cost_amelioration;
  uint16_t cdi_type;
  uint16_t cdi_immunity;
  std::array<double, NUM_RECEPTORS> genome_receptor_expr;
  std::array<double, NUM_RECEPTORS> toxin_affinity;
  std::array<double, NUM_RECEPTORS> ligand_affinity;
  int32_t num_bi_loci;
};

struct AgentTransferData {
  AgentIdentityXfer identity;
  AgentSpatialXfer spatial;
  AgentMetabolismXfer metabolism;
  AgentLifecycleXfer lifecycle;
  AgentMotilityXfer motility;
  AgentGenomeXfer genome;
};

struct BIClusterTransferData {
  uint16_t toxin_id;
  uint16_t immunity_id;
  int32_t  target;
  int32_t  bclass;
  double   pI;
  double   diff_coeff;
  double   retardation;
  double   molecular_weight;
  double   immunity_binding_affinity;
  double   protease_half_life;
  int32_t  release_mode;
  int32_t  is_nuclease;
  double   burst_size;
  double   phage_induction_rate;
  double   phage_burst_size;
  double   phage_lysogeny_rate;
};

void pack_agent(const Agent& a, AgentTransferData& d) {
  d.identity.tag = a.identity.tag;
  d.identity.type = a.identity.type;
  d.identity.owner_rank = a.identity.owner_rank;
  for (int k = 0; k < 3; ++k) {
    d.spatial.x[k] = a.x[k];
    d.spatial.v[k] = a.v[k];
  }
  d.spatial.radius = a.radius;
  d.spatial.mass = a.mass;
  d.spatial.outer_radius = a.outer_radius;
  d.metabolism.mu_max = a.mu_max;
  d.metabolism.mu_realized = a.mu_realized;
  d.metabolism.biomass = a.biomass;
  d.metabolism.maintenance = a.maintenance;
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    d.metabolism.receptor_expr[k] = a.receptor_expr[k];
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    d.metabolism.receptor_expr_base[k] = a.receptor_expr_base[k];
  d.metabolism.km_iron = a.km.km_iron;
  d.metabolism.km_b12 = a.km.km_b12;
  d.metabolism.km_carbon = a.km.km_carbon;
  d.lifecycle.state = static_cast<int32_t>(to_underlying(a.state));
  d.lifecycle.age = a.timers.age;
  d.lifecycle.sos_timer = a.timers.sos_timer;
  d.lifecycle.death_time = a.timers.death_time;
  d.lifecycle.grid_cell = a.grid_cell;
  d.lifecycle.in_crypt = a.flags.in_crypt ? 1 : 0;
  for (int k = 0; k < 3; ++k)
    d.motility.swim_direction[k] = a.motility.swim_direction[k];
  d.motility.swim_speed = a.motility.swim_speed;
  d.motility.run_timer = a.motility.run_timer;
  d.motility.is_stopped = a.motility.is_stopped ? 1 : 0;
  d.motility.stop_timer = a.motility.stop_timer;
  d.motility.prev_carbon = a.motility.prev_carbon;
  d.motility.prev_oxygen = a.motility.prev_oxygen;
  d.motility.prev_ai2 = a.motility.prev_ai2;
  d.genome.lineage_id = a.genome.lineage_id;
  d.genome.parent_id = a.genome.parent_id;
  d.genome.generation = a.genome.generation;
  d.genome.mutations = a.genome.mutations;
  d.genome.has_conjugative_plasmid = a.genome.has_conjugative_plasmid ? 1 : 0;
  d.genome.plasmid_cost_amelioration = a.genome.plasmid_cost_amelioration;
  d.genome.cdi_type = a.genome.cdi_type;
  d.genome.cdi_immunity = a.genome.cdi_immunity;
  for (int k = 0; k < NUM_RECEPTORS; ++k) {
    d.genome.genome_receptor_expr[k] = a.genome.receptor_expression[k];
    d.genome.toxin_affinity[k] = a.genome.toxin_affinity[k];
    d.genome.ligand_affinity[k] = a.genome.ligand_affinity[k];
  }
  d.genome.num_bi_loci = static_cast<int32_t>(a.genome.bi_loci.size());
}

void pack_bi_cluster(const BICluster& c, BIClusterTransferData& d) {
  d.toxin_id = c.toxin_id;
  d.immunity_id = c.immunity_id;
  d.target = static_cast<int32_t>(to_underlying(c.target));
  d.bclass = static_cast<int32_t>(to_underlying(c.bclass));
  d.pI = c.pI;
  d.diff_coeff = c.diff_coeff;
  d.retardation = c.retardation;
  d.molecular_weight = c.molecular_weight;
  d.immunity_binding_affinity = c.immunity_binding_affinity;
  d.protease_half_life = c.protease_half_life;
  d.release_mode = static_cast<int32_t>(to_underlying(c.release_mode));
  d.is_nuclease = c.is_nuclease ? 1 : 0;
  d.burst_size = c.burst_size;
  d.phage_induction_rate = c.phage_induction_rate;
  d.phage_burst_size = c.phage_burst_size;
  d.phage_lysogeny_rate = c.phage_lysogeny_rate;
}

Agent unpack_agent(const AgentTransferData& d, const BIClusterTransferData* bi_data) {
  Agent a = Agent::create_default(d.identity.tag, d.identity.type,
                                  {d.spatial.x[0], d.spatial.x[1], d.spatial.x[2]},
                                  d.metabolism.mu_max);
  a.identity.owner_rank = d.identity.owner_rank;
  for (int k = 0; k < 3; ++k) a.v[k] = d.spatial.v[k];
  a.radius = d.spatial.radius;
  a.mass = d.spatial.mass;
  a.outer_radius = d.spatial.outer_radius;
  a.mu_realized = d.metabolism.mu_realized;
  a.biomass = d.metabolism.biomass;
  a.maintenance = d.metabolism.maintenance;
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    a.receptor_expr[k] = d.metabolism.receptor_expr[k];
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    a.receptor_expr_base[k] = d.metabolism.receptor_expr_base[k];
  a.km.km_iron = d.metabolism.km_iron;
  a.km.km_b12 = d.metabolism.km_b12;
  a.km.km_carbon = d.metabolism.km_carbon;
  a.state = static_cast<PhenoState>(d.lifecycle.state);
  a.timers.age = d.lifecycle.age;
  a.timers.sos_timer = d.lifecycle.sos_timer;
  a.timers.death_time = d.lifecycle.death_time;
  a.grid_cell = d.lifecycle.grid_cell;
  a.flags.in_crypt = (d.lifecycle.in_crypt != 0);
  for (int k = 0; k < 3; ++k)
    a.motility.swim_direction[k] = d.motility.swim_direction[k];
  a.motility.swim_speed = d.motility.swim_speed;
  a.motility.run_timer = d.motility.run_timer;
  a.motility.is_stopped = (d.motility.is_stopped != 0);
  a.motility.stop_timer = d.motility.stop_timer;
  a.motility.prev_carbon = d.motility.prev_carbon;
  a.motility.prev_oxygen = d.motility.prev_oxygen;
  a.motility.prev_ai2 = d.motility.prev_ai2;
  a.genome.lineage_id = d.genome.lineage_id;
  a.genome.parent_id = d.genome.parent_id;
  a.genome.generation = d.genome.generation;
  a.genome.mutations = d.genome.mutations;
  a.genome.has_conjugative_plasmid = (d.genome.has_conjugative_plasmid != 0);
  a.genome.plasmid_cost_amelioration = d.genome.plasmid_cost_amelioration;
  a.genome.cdi_type = d.genome.cdi_type;
  a.genome.cdi_immunity = d.genome.cdi_immunity;
  for (int k = 0; k < NUM_RECEPTORS; ++k) {
    a.genome.receptor_expression[k] = d.genome.genome_receptor_expr[k];
    a.genome.toxin_affinity[k] = d.genome.toxin_affinity[k];
    a.genome.ligand_affinity[k] = d.genome.ligand_affinity[k];
  }

  a.genome.bi_loci.resize(d.genome.num_bi_loci);
  for (int k = 0; k < d.genome.num_bi_loci; ++k) {
    const auto& bd = bi_data[k];
    a.genome.bi_loci[k].toxin_id = bd.toxin_id;
    a.genome.bi_loci[k].immunity_id = bd.immunity_id;
    a.genome.bi_loci[k].target = static_cast<ReceptorType>(bd.target);
    a.genome.bi_loci[k].bclass = static_cast<BacteriocinClass>(bd.bclass);
    a.genome.bi_loci[k].pI = bd.pI;
    a.genome.bi_loci[k].diff_coeff = bd.diff_coeff;
    a.genome.bi_loci[k].retardation = bd.retardation;
    a.genome.bi_loci[k].molecular_weight = bd.molecular_weight;
    a.genome.bi_loci[k].immunity_binding_affinity = bd.immunity_binding_affinity;
    a.genome.bi_loci[k].protease_half_life = bd.protease_half_life;
    a.genome.bi_loci[k].release_mode =
        static_cast<ReleaseMode>(bd.release_mode);
    a.genome.bi_loci[k].is_nuclease = (bd.is_nuclease != 0);
    a.genome.bi_loci[k].burst_size = bd.burst_size;
    a.genome.bi_loci[k].phage_induction_rate = bd.phage_induction_rate;
    a.genome.bi_loci[k].phage_burst_size = bd.phage_burst_size;
    a.genome.bi_loci[k].phage_lysogeny_rate = bd.phage_lysogeny_rate;
  }
  return a;
}

}  // namespace

namespace {

template <typename T>
void append_bytes(std::vector<char>& buf, const T& value) {
  const auto bytes = std::as_bytes(std::span<const T, 1>(&value, 1));
  buf.insert(buf.end(), reinterpret_cast<const char*>(bytes.data()),
             reinterpret_cast<const char*>(bytes.data() + bytes.size()));
}

}  // namespace

void agent_transfer_serialize(const std::vector<Agent>& agents,
                              std::vector<char>& buf) {
  buf.clear();
  auto count = static_cast<int32_t>(agents.size());
  append_bytes(buf, count);
  for (const auto& a : agents) {
    AgentTransferData d;
    pack_agent(a, d);
    append_bytes(buf, d);
    for (const auto& c : a.genome.bi_loci) {
      BIClusterTransferData bd;
      pack_bi_cluster(c, bd);
      append_bytes(buf, bd);
    }
  }
}

std::vector<Agent> agent_transfer_deserialize(const std::vector<char>& buf) {
  std::vector<Agent> result;
  if (buf.size() < sizeof(int32_t)) return result;
  size_t offset = 0;
  int32_t count;
  std::memcpy(&count, buf.data() + offset, sizeof(count));
  offset += sizeof(count);

  result.reserve(count);
  for (int32_t i = 0; i < count; ++i) {
    if (offset + sizeof(AgentTransferData) > buf.size()) break;
    AgentTransferData d;
    std::memcpy(&d, buf.data() + offset, sizeof(d));
    offset += sizeof(d);

    std::vector<BIClusterTransferData> bi_data(d.genome.num_bi_loci);
    for (int32_t k = 0; k < d.genome.num_bi_loci; ++k) {
      if (offset + sizeof(BIClusterTransferData) > buf.size()) break;
      std::memcpy(&bi_data[k], buf.data() + offset, sizeof(bi_data[k]));
      offset += sizeof(bi_data[k]);
    }
    result.push_back(unpack_agent(d, bi_data.data()));
  }
  return result;
}

}  // namespace gutibm
