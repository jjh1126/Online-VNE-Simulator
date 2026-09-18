#pragma once

#include <optional>
#include <vector>

#include "core/request.hpp"

namespace vne {

struct Embedding {
  int departure = 0;
  std::vector<NodeId> virtualToSubstrate;
  std::vector<double> cpuAllocations;
  std::vector<std::vector<EdgeId>> linkPaths;
  // D-/R-ViNE-LB use splittable flows.  Each virtual link owns one or more
  // (path, bandwidth) allocations; linkPaths remains populated for the
  // unsplittable algorithms.
  std::vector<std::vector<std::pair<std::vector<EdgeId>, double>>> linkFlowPaths;
  double revenue = 0.0;
  double cost = 0.0;
  bool boundaryLossComponentsValid = false;
  double nodeBoundaryLossCost = 0.0;
  double linkBoundaryLossCost = 0.0;
  double boundaryLossObjective = 0.0;
};

class VNEAlgorithm {
 public:
  virtual ~VNEAlgorithm() = default;
  virtual std::optional<Embedding> embed(
      SubstrateGraph& substrate,
      const VirtualNetworkRequest& request,
      bool firstRequest,
      RejectReason* reason) = 0;
};

}  // namespace vne
