#include "algorithms/vne_algorithm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

namespace vne {
namespace {

struct BcpPath {
  std::vector<NodeId> nodes;
  std::vector<EdgeId> edges;
};

struct BcpMappingState {
  std::vector<NodeId> nodes;
  std::vector<std::vector<EdgeId>> paths;
  std::vector<char> usedNodes;
  std::vector<char> routedEdges;
};

struct BcpCandidate {
  NodeId substrate = -1;
  double score = 0.0;
  BcpMappingState mapping;
  SubstrateSnapshot snapshot;
};

struct BoundaryLossComponents {
  double node = 0.0;
  double link = 0.0;
};

std::vector<std::vector<int>> bcpIncidence(
    const VirtualNetworkRequest& request) {
  std::vector<std::vector<int>> output(request.nodes.size());
  for (int edge = 0; edge < static_cast<int>(request.edges.size()); ++edge) {
    output[request.edges[edge].a].push_back(edge);
    output[request.edges[edge].b].push_back(edge);
  }
  return output;
}

double nodeBoundaryLoss(
    const SubstrateGraph& graph, const SubstrateSnapshot& state, NodeId node) {
  double capacity = 0.0;
  double residual = 0.0;
  for (EdgeId edge : graph.adjacency[node]) {
    capacity += graph.edges[edge].capacity;
    residual += state.edge[edge];
  }
  if (capacity <= kEps) return 1.0;
  return 1.0 - std::clamp(residual / capacity, 0.0, 1.0);
}

double linkBoundaryLoss(
    const SubstrateGraph& graph, const SubstrateSnapshot& state, EdgeId edge) {
  const auto endpointLoss = [&](NodeId node) -> std::optional<double> {
    double capacity = 0.0;
    double residual = 0.0;
    for (EdgeId other : graph.adjacency[node]) {
      if (other == edge) continue;
      capacity += graph.edges[other].capacity;
      residual += state.edge[other];
    }
    if (capacity <= kEps) return std::nullopt;
    return 1.0 - std::clamp(residual / capacity, 0.0, 1.0);
  };
  const auto uLoss = endpointLoss(graph.edges[edge].u);
  const auto vLoss = endpointLoss(graph.edges[edge].v);
  if (uLoss && vLoss) return std::max(*uLoss, *vLoss);
  if (uLoss) return *uLoss;
  if (vLoss) return *vLoss;
  return 1.0;
}

BoundaryLossComponents boundaryLossComponents(
    const SubstrateGraph& graph, const SubstrateSnapshot& requestStart,
    const SubstrateSnapshot& after) {
  BoundaryLossComponents output;
  for (NodeId node = 0; node < graph.nodeCount; ++node) {
    output.node += std::max(0.0, after.cpu[node]) *
        (nodeBoundaryLoss(graph, after, node) -
         nodeBoundaryLoss(graph, requestStart, node));
  }
  for (EdgeId edge = 0; edge < static_cast<int>(graph.edges.size()); ++edge) {
    output.link += std::max(0.0, after.edge[edge]) *
        (linkBoundaryLoss(graph, after, edge) -
         linkBoundaryLoss(graph, requestStart, edge));
  }
  return output;
}

double boundaryLossObjective(
    const SubstrateGraph& graph, const SubstrateSnapshot& requestStart,
    const SubstrateSnapshot& after, const AlgorithmConfig& config) {
  const auto components = boundaryLossComponents(
      graph, requestStart, after);
  return config.beta * components.link + config.gamma * components.node;
}

std::optional<BcpPath> bcpMinimumHopPath(
    const SubstrateGraph& graph, NodeId source, NodeId target, double demand) {
  if (source == target) return BcpPath{{source}, {}};
  std::vector<NodeId> parent(graph.nodeCount, -1);
  std::vector<EdgeId> parentEdge(graph.nodeCount, -1);
  std::queue<NodeId> queue;
  parent[source] = source;
  queue.push(source);
  while (!queue.empty()) {
    const NodeId current = queue.front();
    queue.pop();
    std::vector<EdgeId> edges = graph.adjacency[current];
    std::sort(edges.begin(), edges.end(), [&](EdgeId left, EdgeId right) {
      const NodeId leftNode = graph.other(left, current);
      const NodeId rightNode = graph.other(right, current);
      return leftNode != rightNode ? leftNode < rightNode : left < right;
    });
    for (EdgeId edge : edges) {
      if (graph.edges[edge].residual + kEps < demand) continue;
      const NodeId next = graph.other(edge, current);
      if (parent[next] >= 0) continue;
      parent[next] = current;
      parentEdge[next] = edge;
      if (next == target) {
        while (!queue.empty()) queue.pop();
        break;
      }
      queue.push(next);
    }
  }
  if (parent[target] < 0) return std::nullopt;
  BcpPath output;
  for (NodeId node = target;; node = parent[node]) {
    output.nodes.push_back(node);
    if (node == source) break;
    output.edges.push_back(parentEdge[node]);
  }
  std::reverse(output.nodes.begin(), output.nodes.end());
  std::reverse(output.edges.begin(), output.edges.end());
  return output;
}

int bcpRootNode(
    const VirtualNetworkRequest& request,
    const std::vector<std::vector<int>>& incident) {
  int best = 0;
  int bestDegree = -1;
  double bestCpu = -1.0;
  double bestBandwidth = -1.0;
  for (int node = 0; node < static_cast<int>(request.nodes.size()); ++node) {
    double bandwidth = 0.0;
    for (int edge : incident[node]) bandwidth += request.edges[edge].bandwidth;
    const int degree = static_cast<int>(incident[node].size());
    if (degree > bestDegree ||
        (degree == bestDegree &&
         (request.nodes[node].cpu > bestCpu + kEps ||
          (std::abs(request.nodes[node].cpu - bestCpu) <= kEps &&
           (bandwidth > bestBandwidth + kEps ||
            (std::abs(bandwidth - bestBandwidth) <= kEps && node < best)))))) {
      best = node;
      bestDegree = degree;
      bestCpu = request.nodes[node].cpu;
      bestBandwidth = bandwidth;
    }
  }
  return best;
}

std::vector<int> bcpBfsOrder(
    const VirtualNetworkRequest& request,
    const std::vector<std::vector<int>>& incident, int root) {
  std::vector<int> output;
  std::vector<char> seen(request.nodes.size(), 0);
  const auto bfs = [&](int start) {
    std::queue<int> queue;
    seen[start] = 1;
    queue.push(start);
    while (!queue.empty()) {
      const int node = queue.front();
      queue.pop();
      output.push_back(node);
      std::vector<int> neighbors;
      for (int edge : incident[node]) {
        neighbors.push_back(request.edges[edge].a == node
            ? request.edges[edge].b : request.edges[edge].a);
      }
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()),
                      neighbors.end());
      for (int next : neighbors) {
        if (seen[next]) continue;
        seen[next] = 1;
        queue.push(next);
      }
    }
  };
  bfs(root);
  for (int node = 0; node < static_cast<int>(request.nodes.size()); ++node)
    if (!seen[node]) bfs(node);
  return output;
}

std::optional<NodeId> bcpFirstNodeMapping(
    const SubstrateGraph& graph, const VirtualNetworkRequest& request,
    int root) {
  const SubstrateSnapshot state = graph.snapshot();
  NodeId best = -1;
  double bestBoundaryLoss = std::numeric_limits<double>::infinity();
  double bestResidualCpu = -1.0;
  for (NodeId candidate = 0; candidate < graph.nodeCount; ++candidate) {
    if (graph.cpuResidual[candidate] + kEps < request.nodes[root].cpu)
      continue;
    const double loss = nodeBoundaryLoss(graph, state, candidate);
    const double residualCpu = graph.cpuResidual[candidate];
    const bool better = best < 0 || loss < bestBoundaryLoss - kEps ||
        (std::abs(loss - bestBoundaryLoss) <= kEps &&
         (residualCpu > bestResidualCpu + kEps ||
          (std::abs(residualCpu - bestResidualCpu) <= kEps &&
           candidate < best)));
    if (better) {
      best = candidate;
      bestBoundaryLoss = loss;
      bestResidualCpu = residualCpu;
    }
  }
  return best < 0 ? std::nullopt : std::optional<NodeId>(best);
}

class BcpAlgorithm final : public VNEAlgorithm {
 public:
  explicit BcpAlgorithm(AlgorithmConfig config) : config_(config) {}

  std::optional<Embedding> embed(
      SubstrateGraph& graph, const VirtualNetworkRequest& request,
      bool firstRequest, RejectReason* reason) override {
    (void)firstRequest;
    if (reason) *reason = RejectReason::Route;
    if (request.nodes.size() > static_cast<size_t>(graph.nodeCount)) {
      if (reason) *reason = RejectReason::Cpu;
      return std::nullopt;
    }

    const SubstrateSnapshot requestStart = graph.snapshot();
    const auto incident = bcpIncidence(request);
    const int root = bcpRootNode(request, incident);
    const auto order = bcpBfsOrder(request, incident, root);
    const auto selectedRoot = bcpFirstNodeMapping(graph, request, root);
    if (!selectedRoot) {
      if (reason) *reason = RejectReason::Cpu;
      return std::nullopt;
    }
    BcpMappingState mapping{
        std::vector<NodeId>(request.nodes.size(), -1),
        std::vector<std::vector<EdgeId>>(request.edges.size()),
        std::vector<char>(graph.nodeCount, 0),
        std::vector<char>(request.edges.size(), 0)};

    const auto evaluate = [&](int virtualNode, NodeId substrate)
        -> std::optional<BcpCandidate> {
      if (mapping.usedNodes[substrate] ||
          graph.cpuResidual[substrate] + kEps <
              request.nodes[virtualNode].cpu) return std::nullopt;
      const SubstrateSnapshot committed = graph.snapshot();
      const BcpMappingState backup = mapping;
      graph.allocateCpu(substrate, request.nodes[virtualNode].cpu);
      mapping.nodes[virtualNode] = substrate;
      mapping.usedNodes[substrate] = 1;

      std::vector<int> pending;
      for (int edge : incident[virtualNode]) {
        const int neighbor = request.edges[edge].a == virtualNode
            ? request.edges[edge].b : request.edges[edge].a;
        if (!mapping.routedEdges[edge] && mapping.nodes[neighbor] >= 0)
          pending.push_back(edge);
      }
      std::sort(pending.begin(), pending.end());

      bool feasible = true;
      for (int edge : pending) {
        const auto& virtualEdge = request.edges[edge];
        const int neighbor = virtualEdge.a == virtualNode
            ? virtualEdge.b : virtualEdge.a;
        const auto path = bcpMinimumHopPath(
            graph, substrate, mapping.nodes[neighbor], virtualEdge.bandwidth);
        if (!path) {
          feasible = false;
          break;
        }
        graph.allocateBandwidth(path->edges, virtualEdge.bandwidth);
        mapping.paths[edge] = path->edges;
        mapping.routedEdges[edge] = 1;
      }

      std::optional<BcpCandidate> output;
      if (feasible) {
        output = BcpCandidate{
            substrate,
            boundaryLossObjective(graph, requestStart, graph.snapshot(), config_),
            mapping, graph.snapshot()};
      }
      graph.restore(committed);
      mapping = backup;
      return output;
    };

    for (int index = 0; index < static_cast<int>(order.size()); ++index) {
      const int virtualNode = order[index];
      std::vector<BcpCandidate> candidates;
      bool hasCpuCandidate = false;
      if (index == 0) {
        hasCpuCandidate = true;
        if (const auto candidate = evaluate(virtualNode, *selectedRoot))
          candidates.push_back(*candidate);
      } else {
        for (NodeId substrate = 0; substrate < graph.nodeCount; ++substrate) {
          if (mapping.usedNodes[substrate] ||
              graph.cpuResidual[substrate] + kEps <
                  request.nodes[virtualNode].cpu) continue;
          hasCpuCandidate = true;
          if (const auto candidate = evaluate(virtualNode, substrate))
            candidates.push_back(*candidate);
        }
      }
      if (candidates.empty()) {
        graph.restore(requestStart);
        if (reason)
          *reason = hasCpuCandidate ? RejectReason::Route : RejectReason::Cpu;
        return std::nullopt;
      }
      const auto best = std::min_element(
          candidates.begin(), candidates.end(),
          [](const BcpCandidate& left, const BcpCandidate& right) {
            return left.score < right.score - kEps ||
                (std::abs(left.score - right.score) <= kEps &&
                 left.substrate < right.substrate);
          });
      graph.restore(best->snapshot);
      mapping = best->mapping;
    }

    Embedding output;
    output.departure = request.departure;
    output.virtualToSubstrate = mapping.nodes;
    output.linkPaths = mapping.paths;
    for (const auto& node : request.nodes) {
      output.cpuAllocations.push_back(node.cpu);
      output.revenue += node.cpu;
      output.cost += node.cpu;
    }
    for (int edge = 0; edge < static_cast<int>(request.edges.size()); ++edge) {
      output.revenue += request.edges[edge].bandwidth;
      output.cost += request.edges[edge].bandwidth *
          static_cast<double>(mapping.paths[edge].size());
    }
    const auto components = boundaryLossComponents(
        graph, requestStart, graph.snapshot());
    output.boundaryLossComponentsValid = true;
    output.nodeBoundaryLossCost = components.node;
    output.linkBoundaryLossCost = components.link;
    output.boundaryLossObjective = config_.beta * components.link +
        config_.gamma * components.node;
    if (reason) *reason = RejectReason::None;
    return output;
  }

 private:
  AlgorithmConfig config_;
};

}  // namespace

std::optional<Embedding> embedBcp(
    SubstrateGraph& graph, const VirtualNetworkRequest& request,
    bool firstRequest, AlgorithmConfig config, RejectReason* reason) {
  return BcpAlgorithm(config).embed(graph, request, firstRequest, reason);
}

}  // namespace vne
