#include "simulator/simulator.hpp"

#include "algorithms/vineyard_algorithm.hpp"

#include <algorithm>
#include <chrono>
#include <optional>
#include <queue>

namespace vne {
std::optional<Embedding> embedBcp(
    SubstrateGraph& graph, const VirtualNetworkRequest& request,
    bool firstRequest, AlgorithmConfig config, RejectReason* reason);

namespace {

double cpuUtilization(const SubstrateGraph& graph) {
  double capacity = 0.0, residual = 0.0;
  for (int node = 0; node < graph.nodeCount; ++node) { capacity += graph.cpuCapacity[node]; residual += graph.cpuResidual[node]; }
  return capacity > 0.0 ? 1.0 - residual / capacity : 0.0;
}
double linkUtilization(const SubstrateGraph& graph) {
  double capacity = 0.0, residual = 0.0;
  for (const auto& edge : graph.edges) { capacity += edge.capacity; residual += edge.residual; }
  return capacity > 0.0 ? 1.0 - residual / capacity : 0.0;
}
void releaseEmbedding(SubstrateGraph& graph, const Embedding& embedding, const VirtualNetworkRequest& request) {
  for (int node = 0; node < static_cast<int>(embedding.virtualToSubstrate.size()); ++node) graph.releaseCpu(embedding.virtualToSubstrate[node], embedding.cpuAllocations[node]);
  if (!embedding.linkFlowPaths.empty()) {
    for (int edge = 0; edge < static_cast<int>(embedding.linkFlowPaths.size()); ++edge)
      for (const auto& flow : embedding.linkFlowPaths[edge]) graph.releaseBandwidth(flow.first, flow.second);
  } else {
    for (int edge = 0; edge < static_cast<int>(embedding.linkPaths.size()); ++edge) graph.releaseBandwidth(embedding.linkPaths[edge], request.edges[edge].bandwidth);
  }
}

}  // namespace

Simulator::Simulator(
    const SubstrateGraph& substrate,
    const std::vector<VirtualNetworkRequest>& requests,
    RFDVariant variant,
    SimulatorConfig config)
    : substrate_(substrate), requests_(requests), variant_(variant),
      config_(config) {}

SimulationResult Simulator::run() {
  const auto started = std::chrono::steady_clock::now();
  RFDAlgorithm algorithm(variant_, config_.algorithm);
  VineyardAlgorithm vineyard(
      config_.randomizedVineyard ? VineyardMode::Randomized : VineyardMode::Deterministic,
      config_.algorithm, config_.substrateSeed, config_.requestSeed);
  if (!config_.useBcp && !config_.useVineyard)
    algorithm.initialize(substrate_);
  int horizon = 1;
  for (const auto& request : requests_) horizon = std::max(horizon, request.departure + 1);
  if (config_.horizon > 0) horizon = config_.horizon;
  if (config_.warmup < 0 || config_.warmup >= horizon) throw std::runtime_error("warmup must be within the simulation horizon");
  std::vector<std::vector<int>> arrivals(horizon);
  for (int index = 0; index < static_cast<int>(requests_.size()); ++index) if (requests_[index].arrival < horizon) arrivals[requests_[index].arrival].push_back(index);

  struct Active { Embedding embedding; int request = -1; };
  std::vector<std::optional<Active>> active(requests_.size());
  std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> departures;
  SimulationResult result;
  if (config_.useBcp) {
    result.summary.algorithm = "VNE-BCP";
  } else if (config_.useVineyard) {
    result.summary.algorithm = config_.randomizedVineyard ? "R-ViNE-LB" : "D-ViNE-LB";
  } else {
    result.summary.algorithm = algorithmName(variant_);
  }
  result.summary.substrateSeed = config_.substrateSeed;
  result.summary.requestSeed = config_.requestSeed;
  result.summary.warmupTime = config_.warmup;
  result.summary.measurementEnd = horizon;
  bool firstRequest = true;
  double cpuSum = 0.0, linkSum = 0.0;

  for (int time = 0; time < horizon; ++time) {
    while (!departures.empty() && departures.top().first <= time) {
      const int requestIndex = departures.top().second; departures.pop();
      if (active[requestIndex]) { releaseEmbedding(substrate_, active[requestIndex]->embedding, requests_[requestIndex]); active[requestIndex].reset(); }
    }
    for (int requestIndex : arrivals[time]) {
      const auto& request = requests_[requestIndex];
      RejectReason reason = RejectReason::Route;
      auto embedding = config_.useBcp
          ? embedBcp(substrate_, request, firstRequest,
                     config_.algorithm, &reason)
          : config_.useVineyard
              ? vineyard.embed(substrate_, request, firstRequest, &reason)
              : algorithm.embed(substrate_, request, firstRequest, &reason);
      // Algorithm 3 applies to the first VNR attempt, not the first accepted VNR.
      firstRequest = false;
      const bool measured = request.arrival >= config_.warmup;
      RequestRecord record;
      record.algorithm = result.summary.algorithm;
      record.substrateSeed = config_.substrateSeed;
      record.requestSeed = config_.requestSeed;
      record.requestId = request.id;
      record.arrival = request.arrival;
      record.departure = request.departure;
      record.virtualNodes = static_cast<int>(request.nodes.size());
      record.virtualEdges = static_cast<int>(request.edges.size());
      record.accepted = embedding.has_value();
      record.reason = reason;
      record.revenue = embedding ? embedding->revenue : 0.0;
      record.cost = embedding ? embedding->cost : 0.0;
      record.cpuUtilization = cpuUtilization(substrate_);
      record.linkUtilization = linkUtilization(substrate_);
      if (embedding && embedding->boundaryLossComponentsValid) {
        record.boundaryLossComponentsValid = true;
        record.nodeBoundaryLossCost = embedding->nodeBoundaryLossCost;
        record.linkBoundaryLossCost = embedding->linkBoundaryLossCost;
        record.boundaryLossObjective = embedding->boundaryLossObjective;
      }
      if (config_.useVineyard) record.vineyard = vineyard.lastDiagnostic();
      if (measured) {
        result.records.push_back(record);
        ++result.summary.totalRequests;
      } else {
        ++result.summary.warmupRequests;
      }
      if (embedding) {
        if (measured) {
          ++result.summary.accepted;
          result.summary.revenue += embedding->revenue;
          result.summary.cost += embedding->cost;
        } else {
          ++result.summary.warmupAccepted;
        }
        active[requestIndex] = Active{std::move(*embedding), requestIndex};
        departures.push({request.departure, requestIndex});
      }
    }
    if (time >= config_.warmup) {
      cpuSum += cpuUtilization(substrate_);
      linkSum += linkUtilization(substrate_);
    }
  }
  const int measurementDuration = horizon - config_.warmup;
  result.summary.averageCpuUtilization = cpuSum / static_cast<double>(measurementDuration);
  result.summary.averageLinkUtilization = linkSum / static_cast<double>(measurementDuration);
  result.summary.runtimeMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started).count();
  return result;
}

}  // namespace vne
