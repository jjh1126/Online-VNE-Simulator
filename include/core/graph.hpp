#pragma once

#include <filesystem>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace vne {

struct SubstrateEdge {
  NodeId u = -1;
  NodeId v = -1;
  double capacity = 0.0;
  double residual = 0.0;
};

struct SubstrateSnapshot {
  std::vector<double> cpu;
  std::vector<double> edge;
};

struct NodeLocation {
  bool present = false;
  double x = 0.0;
  double y = 0.0;
};

struct VirtualNode {
  double cpu = 0.0;
  bool hasLocation = false;
  double x = 0.0;
  double y = 0.0;
  double radius = 0.0;
};
struct VirtualEdge { NodeId a = -1; NodeId b = -1; double bandwidth = 0.0; };

struct SubstrateGraph {
  int nodeCount = 0;
  std::vector<double> cpuCapacity;
  std::vector<double> cpuResidual;
  std::vector<NodeLocation> nodeLocations;
  std::vector<SubstrateEdge> edges;
  std::vector<std::vector<EdgeId>> adjacency;
  std::vector<std::vector<EdgeId>> edgeIndex;

  void initialize(int nodes);
  void addEdge(NodeId u, NodeId v, double capacity);
  NodeId other(EdgeId edge, NodeId node) const;
  void sortAdjacency();
  SubstrateSnapshot snapshot() const;
  void restore(const SubstrateSnapshot& snapshot);
  bool allocateCpu(NodeId node, double amount);
  void releaseCpu(NodeId node, double amount);
  void allocateBandwidth(const std::vector<EdgeId>& path, double amount);
  bool tryAllocateBandwidth(
      const std::vector<std::pair<std::vector<EdgeId>, double>>& allocations,
      double relativeTolerance);
  void releaseBandwidth(const std::vector<EdgeId>& path, double amount);
};

SubstrateGraph loadSubstrate(const std::filesystem::path& path);

}  // namespace vne
