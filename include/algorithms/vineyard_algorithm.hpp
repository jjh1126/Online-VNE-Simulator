#pragma once

#include "algorithms/vne_algorithm.hpp"

namespace vne {

// ViNEYard's location-constrained augmented-LP embedding.  Deterministic is
// D-ViNE-LB; randomized rounding is R-ViNE-LB.
enum class VineyardMode { Deterministic, Randomized };

class VineyardAlgorithm final : public VNEAlgorithm {
 public:
  VineyardAlgorithm(VineyardMode mode, AlgorithmConfig config, int substrateSeed, int requestSeed);
  std::optional<Embedding> embed(
      SubstrateGraph& substrate, const VirtualNetworkRequest& request,
      bool firstRequest, RejectReason* reason) override;
  const VineyardDiagnostic& lastDiagnostic() const { return lastDiagnostic_; }

 private:
  VineyardMode mode_;
  AlgorithmConfig config_;
  int substrateSeed_;
  int requestSeed_;
  VineyardDiagnostic lastDiagnostic_;
};

}  // namespace vne
