#pragma once

#include <utility>

#include "algorithms/vne_algorithm.hpp"

namespace vne {

class RFDAlgorithm final : public VNEAlgorithm {
 public:
  RFDAlgorithm(RFDVariant variant, AlgorithmConfig config);
  void initialize(const SubstrateGraph& substrate);
  RFDMetric measure(const SubstrateGraph& substrate) const;
  std::optional<Embedding> embed(
      SubstrateGraph& substrate,
      const VirtualNetworkRequest& request,
      bool firstRequest,
      RejectReason* reason) override;

 private:
  struct StaticData {
    int nodes = 0;
    std::vector<int> nodeDegree;
    std::vector<std::vector<std::vector<std::vector<EdgeId>>>> paths;
    std::vector<std::vector<std::pair<EdgeId, NodeId>>> adjacentLinks;
  };

  RFDVariant variant_;
  AlgorithmConfig config_;
  StaticData static_;
};

}  // namespace vne
