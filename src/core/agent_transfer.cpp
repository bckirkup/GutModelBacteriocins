/* -----------------------------------------------------------------------
   GutIBM – MPI agent serialization implementation
   ----------------------------------------------------------------------- */

#include "agent_transfer.h"
#include "types.h"
#include <cstring>

namespace gutibm {

namespace {

struct AgentTransferData {
  int64_t  tag;
  int32_t  type;
  int32_t  owner_rank;
  double   x[3];
  double   v[3];
  double   radius;
  double   mass;
  double   outer_radius;
  double   mu_max;
  double   mu_realized;
  double   biomass;
  double   maintenance;
  double   receptor_expr[NUM_RECEPTORS];
  double   km_iron;
  double   km_b12;
  double   km_carbon;
  int32_t  state;
  double   age;
  double   sos_timer;
  int32_t  grid_cell;
  int32_t  in_crypt;
  int64_t  lineage_id;
  int64_t  parent_id;
  uint32_t generation;
  uint32_t mutations;
  int32_t  has_conjugative_plasmid;
  double   plasmid_cost_amelioration;
  double   genome_receptor_expr[NUM_RECEPTORS];
  double   toxin_affinity[NUM_RECEPTORS];
  double   ligand_affinity[NUM_RECEPTORS];
  int32_t  num_bi_loci;
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
};

void pack_agent(const Agent& a, AgentTransferData& d) {
  d.tag = a.tag;
  d.type = a.type;
  d.owner_rank = a.owner_rank;
  for (int k = 0; k < 3; ++k) {
    d.x[k] = a.x[k];
    d.v[k] = a.v[k];
  }
  d.radius = a.radius;
  d.mass = a.mass;
  d.outer_radius = a.outer_radius;
  d.mu_max = a.mu_max;
  d.mu_realized = a.mu_realized;
  d.biomass = a.biomass;
  d.maintenance = a.maintenance;
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    d.receptor_expr[k] = a.receptor_expr[k];
  d.km_iron = a.km_iron;
  d.km_b12 = a.km_b12;
  d.km_carbon = a.km_carbon;
  d.state = static_cast<int32_t>(a.state);
  d.age = a.age;
  d.sos_timer = a.sos_timer;
  d.grid_cell = a.grid_cell;
  d.in_crypt = a.in_crypt ? 1 : 0;
  d.lineage_id = a.genome.lineage_id;
  d.parent_id = a.genome.parent_id;
  d.generation = a.genome.generation;
  d.mutations = a.genome.mutations;
  d.has_conjugative_plasmid = a.genome.has_conjugative_plasmid ? 1 : 0;
  d.plasmid_cost_amelioration = a.genome.plasmid_cost_amelioration;
  for (int k = 0; k < NUM_RECEPTORS; ++k) {
    d.genome_receptor_expr[k] = a.genome.receptor_expression[k];
    d.toxin_affinity[k] = a.genome.toxin_affinity[k];
    d.ligand_affinity[k] = a.genome.ligand_affinity[k];
  }
  d.num_bi_loci = static_cast<int32_t>(a.genome.bi_loci.size());
}

void pack_bi_cluster(const BICluster& c, BIClusterTransferData& d) {
  d.toxin_id = c.toxin_id;
  d.immunity_id = c.immunity_id;
  d.target = static_cast<int32_t>(c.target);
  d.bclass = static_cast<int32_t>(c.bclass);
  d.pI = c.pI;
  d.diff_coeff = c.diff_coeff;
  d.retardation = c.retardation;
  d.molecular_weight = c.molecular_weight;
  d.immunity_binding_affinity = c.immunity_binding_affinity;
}

Agent unpack_agent(const AgentTransferData& d, const BIClusterTransferData* bi_data) {
  Agent a = Agent::create_default(d.tag, d.type, {d.x[0], d.x[1], d.x[2]}, d.mu_max);
  a.owner_rank = d.owner_rank;
  for (int k = 0; k < 3; ++k) a.v[k] = d.v[k];
  a.radius = d.radius;
  a.mass = d.mass;
  a.outer_radius = d.outer_radius;
  a.mu_realized = d.mu_realized;
  a.biomass = d.biomass;
  a.maintenance = d.maintenance;
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    a.receptor_expr[k] = d.receptor_expr[k];
  a.km_iron = d.km_iron;
  a.km_b12 = d.km_b12;
  a.km_carbon = d.km_carbon;
  a.state = static_cast<PhenoState>(d.state);
  a.age = d.age;
  a.sos_timer = d.sos_timer;
  a.grid_cell = d.grid_cell;
  a.in_crypt = (d.in_crypt != 0);
  a.genome.lineage_id = d.lineage_id;
  a.genome.parent_id = d.parent_id;
  a.genome.generation = d.generation;
  a.genome.mutations = d.mutations;
  a.genome.has_conjugative_plasmid = (d.has_conjugative_plasmid != 0);
  a.genome.plasmid_cost_amelioration = d.plasmid_cost_amelioration;
  for (int k = 0; k < NUM_RECEPTORS; ++k) {
    a.genome.receptor_expression[k] = d.genome_receptor_expr[k];
    a.genome.toxin_affinity[k] = d.toxin_affinity[k];
    a.genome.ligand_affinity[k] = d.ligand_affinity[k];
  }

  a.genome.bi_loci.resize(d.num_bi_loci);
  for (int k = 0; k < d.num_bi_loci; ++k) {
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
  }
  return a;
}

}  // namespace

void agent_transfer_serialize(const std::vector<Agent>& agents,
                              std::vector<char>& buf) {
  buf.clear();
  int32_t count = static_cast<int32_t>(agents.size());
  buf.insert(buf.end(), reinterpret_cast<const char*>(&count),
             reinterpret_cast<const char*>(&count) + sizeof(count));
  for (const auto& a : agents) {
    AgentTransferData d;
    pack_agent(a, d);
    buf.insert(buf.end(), reinterpret_cast<const char*>(&d),
               reinterpret_cast<const char*>(&d) + sizeof(d));
    for (const auto& c : a.genome.bi_loci) {
      BIClusterTransferData bd;
      pack_bi_cluster(c, bd);
      buf.insert(buf.end(), reinterpret_cast<const char*>(&bd),
                 reinterpret_cast<const char*>(&bd) + sizeof(bd));
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

    std::vector<BIClusterTransferData> bi_data(d.num_bi_loci);
    for (int32_t k = 0; k < d.num_bi_loci; ++k) {
      if (offset + sizeof(BIClusterTransferData) > buf.size()) break;
      std::memcpy(&bi_data[k], buf.data() + offset, sizeof(bi_data[k]));
      offset += sizeof(bi_data[k]);
    }
    result.push_back(unpack_agent(d, bi_data.data()));
  }
  return result;
}

}  // namespace gutibm
