#pragma once

#include "algorithms/rfd_algorithm.hpp"
#include "simulator/metrics.hpp"

namespace vne {

struct SimulatorConfig {
  AlgorithmConfig algorithm;
  int horizon = -1;
  int warmup = 0;
  int substrateSeed = 0;
  int requestSeed = 0;
  bool useBcp = false;
  bool useVineyard = false;
  bool randomizedVineyard = false;
};

class Simulator {
 public:
  Simulator(
      const SubstrateGraph& substrate,
      const std::vector<VirtualNetworkRequest>& requests,
      RFDVariant variant,
      SimulatorConfig config);
  SimulationResult run();

 private:
  SubstrateGraph substrate_;
  std::vector<VirtualNetworkRequest> requests_;
  RFDVariant variant_;
  SimulatorConfig config_;
};

}  // namespace vne
