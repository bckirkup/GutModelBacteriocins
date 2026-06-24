/* -----------------------------------------------------------------------
   GutIBM – MPI agent serialization (pack/unpack for rank migration)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_AGENT_TRANSFER_H
#define GUTIBM_AGENT_TRANSFER_H

#include "agent.h"
#include <vector>

namespace gutibm {

// Serialize agents to a byte buffer (count prefix + packed records).
void agent_transfer_serialize(const std::vector<Agent>& agents,
                              std::vector<char>& buf);

// Deserialize agents from a byte buffer produced by agent_transfer_serialize.
std::vector<Agent> agent_transfer_deserialize(const std::vector<char>& buf);

}  // namespace gutibm

#endif  // GUTIBM_AGENT_TRANSFER_H
