#include "algorithms/rfd_algorithm.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <random>
#include <stdexcept>

namespace vne {
namespace {

struct Path { std::vector<NodeId> nodes; std::vector<EdgeId> edges; };
struct MappingState { std::vector<NodeId> nodes; std::vector<std::vector<EdgeId>> paths; std::vector<char> used; };
struct Candidate { NodeId substrate = -1; double score = 0.0; };

std::vector<std::vector<int>> incidence(const VirtualNetworkRequest& request) {
  std::vector<std::vector<int>> output(request.nodes.size());
  for (int edge = 0; edge < static_cast<int>(request.edges.size()); ++edge) {
    output[request.edges[edge].a].push_back(edge);
    output[request.edges[edge].b].push_back(edge);
  }
  return output;
}

double bottleneckRatio(const SubstrateGraph& graph, const std::vector<EdgeId>& path) {
  if (path.empty()) return 0.0;
  EdgeId bottleneck = path.front();
  for (EdgeId edge : path) {
    if (graph.edges[edge].residual < graph.edges[bottleneck].residual - kEps ||
        (std::abs(graph.edges[edge].residual - graph.edges[bottleneck].residual) <= kEps && edge < bottleneck)) bottleneck = edge;
  }
  return graph.edges[bottleneck].capacity > 0.0 ? std::clamp(graph.edges[bottleneck].residual / graph.edges[bottleneck].capacity, 0.0, 1.0) : 0.0;
}

bool pathLess(const Path& a, const Path& b) {
  return a.edges.size() != b.edges.size() ? a.edges.size() < b.edges.size() : a.edges < b.edges;
}

double linkSelectionScore(
    const AlgorithmConfig& config, double bandwidth, const Path& path,
    double before, double after) {
  return config.beta * bandwidth * static_cast<double>(path.edges.size()) +
      config.lambda * (after - before);
}

double candidateCost(
    const AlgorithmConfig& config, double cpuCost, double pathCost,
    double before, double after) {
  return config.alpha * cpuCost + config.beta * pathCost +
      config.lambda * (after - before);
}

std::optional<Path> shortestFeasiblePath(
    const SubstrateGraph& graph, NodeId source, NodeId target, double demand,
    int maxHops, const std::vector<char>& bannedNodes,
    const std::vector<char>& bannedEdges) {
  if (source == target) return Path{{source}, {}};
  std::vector<NodeId> parent(graph.nodeCount, -1);
  std::vector<EdgeId> parentEdge(graph.nodeCount, -1);
  std::vector<int> depth(graph.nodeCount, -1);
  std::queue<NodeId> queue;
  parent[source] = source; depth[source] = 0; queue.push(source);
  while (!queue.empty()) {
    const NodeId current = queue.front(); queue.pop();
    if (depth[current] >= maxHops) continue;
    std::vector<EdgeId> edges = graph.adjacency[current];
    std::sort(edges.begin(), edges.end(), [&](EdgeId a, EdgeId b) {
      const NodeId va = graph.other(a, current), vb = graph.other(b, current);
      return va != vb ? va < vb : a < b;
    });
    for (EdgeId edge : edges) {
      if (bannedEdges[edge] || graph.edges[edge].residual + kEps < demand) continue;
      const NodeId next = graph.other(edge, current);
      if ((bannedNodes[next] && next != target) || parent[next] >= 0) continue;
      parent[next] = current; parentEdge[next] = edge; depth[next] = depth[current] + 1;
      if (next == target) {
        while (!queue.empty()) queue.pop();
        break;
      }
      queue.push(next);
    }
  }
  if (parent[target] < 0) return std::nullopt;
  Path output;
  for (NodeId node = target;; node = parent[node]) {
    output.nodes.push_back(node);
    if (node == source) break;
    output.edges.push_back(parentEdge[node]);
  }
  std::reverse(output.nodes.begin(), output.nodes.end());
  std::reverse(output.edges.begin(), output.edges.end());
  return output;
}

std::vector<Path> yenPaths(
    const SubstrateGraph& graph, NodeId source, NodeId target, double demand,
    int maxHops, int count) {
  if (count <= 0 || maxHops < 0) return {};
  std::vector<char> noNodes(graph.nodeCount, 0), noEdges(graph.edges.size(), 0);
  auto first = shortestFeasiblePath(graph, source, target, demand, maxHops, noNodes, noEdges);
  if (!first) return {};
  std::vector<Path> accepted{*first};
  std::vector<Path> candidates;
  const auto contains = [](const std::vector<Path>& paths, const Path& value) {
    return std::any_of(paths.begin(), paths.end(), [&](const Path& candidate) { return candidate.edges == value.edges; });
  };
  for (int k = 1; k < count; ++k) {
    const Path previous = accepted.back();
    for (int spur = 0; spur + 1 < static_cast<int>(previous.nodes.size()); ++spur) {
      std::vector<char> bannedNodes(graph.nodeCount, 0), bannedEdges(graph.edges.size(), 0);
      for (int i = 0; i < spur; ++i) bannedNodes[previous.nodes[i]] = 1;
      for (const Path& path : accepted) {
        if (static_cast<int>(path.nodes.size()) <= spur + 1 ||
            !std::equal(previous.nodes.begin(), previous.nodes.begin() + spur + 1, path.nodes.begin())) continue;
        bannedEdges[path.edges[spur]] = 1;
      }
      auto tail = shortestFeasiblePath(
          graph, previous.nodes[spur], target, demand, maxHops - spur,
          bannedNodes, bannedEdges);
      if (!tail) continue;
      Path joined;
      joined.nodes.insert(joined.nodes.end(), previous.nodes.begin(), previous.nodes.begin() + spur);
      joined.nodes.insert(joined.nodes.end(), tail->nodes.begin(), tail->nodes.end());
      joined.edges.insert(joined.edges.end(), previous.edges.begin(), previous.edges.begin() + spur);
      joined.edges.insert(joined.edges.end(), tail->edges.begin(), tail->edges.end());
      if (static_cast<int>(joined.edges.size()) > maxHops) continue;
      if (!contains(accepted, joined) && !contains(candidates, joined)) candidates.push_back(std::move(joined));
    }
    if (candidates.empty()) break;
    const auto best = std::min_element(candidates.begin(), candidates.end(), pathLess);
    accepted.push_back(*best);
    candidates.erase(best);
  }
  return accepted;
}

int hopDistance(const SubstrateGraph& graph, NodeId source, NodeId target) {
  std::vector<int> distance(graph.nodeCount, -1); std::queue<NodeId> queue;
  distance[source] = 0; queue.push(source);
  while (!queue.empty()) {
    const NodeId current = queue.front(); queue.pop();
    if (current == target) return distance[current];
    for (EdgeId edge : graph.adjacency[current]) {
      const NodeId next = graph.other(edge, current);
      if (distance[next] < 0) { distance[next] = distance[current] + 1; queue.push(next); }
    }
  }
  return std::numeric_limits<int>::max() / 4;
}

int rootNode(const VirtualNetworkRequest& request, const std::vector<std::vector<int>>& incident) {
  int best = 0; int bestDegree = -1; double bestCpu = -1.0, bestBandwidth = -1.0;
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
      best = node; bestDegree = degree; bestCpu = request.nodes[node].cpu; bestBandwidth = bandwidth;
    }
  }
  return best;
}

std::vector<int> visitOrder(
    const VirtualNetworkRequest& request, const std::vector<std::vector<int>>& incident,
    int root, RFDVariant variant) {
  std::vector<int> output; std::vector<char> seen(request.nodes.size(), 0);
  const auto neighbors = [&](int node) {
    std::vector<int> result;
    for (int edge : incident[node]) result.push_back(request.edges[edge].a == node ? request.edges[edge].b : request.edges[edge].a);
    std::sort(result.begin(), result.end()); result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
  };
  std::function<void(int)> dfs = [&](int node) { seen[node] = 1; output.push_back(node); for (int next : neighbors(node)) if (!seen[next]) dfs(next); };
  const auto bfs = [&](int start) { std::queue<int> queue; seen[start] = 1; queue.push(start); while (!queue.empty()) { const int node = queue.front(); queue.pop(); output.push_back(node); for (int next : neighbors(node)) if (!seen[next]) { seen[next] = 1; queue.push(next); } } };
  if (variant == RFDVariant::BreadthFirst) bfs(root); else dfs(root);
  for (int node = 0; node < static_cast<int>(request.nodes.size()); ++node) if (!seen[node]) { if (variant == RFDVariant::BreadthFirst) bfs(node); else dfs(node); }
  return output;
}

}  // namespace

RFDAlgorithm::RFDAlgorithm(RFDVariant variant, AlgorithmConfig config)
    : variant_(variant), config_(config) {}

void RFDAlgorithm::initialize(const SubstrateGraph& graph) {
  static_ = {};
  static_.nodes = graph.nodeCount;
  static_.nodeDegree.assign(graph.nodeCount, 0);
  static_.paths.assign(graph.nodeCount, std::vector<std::vector<std::vector<EdgeId>>>(graph.nodeCount));
  for (NodeId source = 0; source < graph.nodeCount; ++source) {
    // CPU-less nodes are forwarding-only switches.  They remain available to
    // the path enumerator as transit nodes, but are not meaningful RFD path
    // endpoints.
    if (graph.cpuCapacity[source] <= kEps) continue;
    for (NodeId target = source + 1; target < graph.nodeCount; ++target) {
      if (graph.cpuCapacity[target] <= kEps) continue;
      std::vector<char> used(graph.nodeCount, 0); std::vector<EdgeId> path; std::vector<std::vector<EdgeId>> all;
      used[source] = 1;
      std::function<void(NodeId)> enumerate = [&](NodeId current) {
        if (current == target) { all.push_back(path); return; }
        if (static_cast<int>(path.size()) >= config_.tau) return;
        for (EdgeId edge : graph.adjacency[current]) {
          const NodeId next = graph.other(edge, current);
          if (!used[next]) { used[next] = 1; path.push_back(edge); enumerate(next); path.pop_back(); used[next] = 0; }
        }
      };
      enumerate(source);
      static_.paths[source][target] = static_.paths[target][source] = std::move(all);
      if (!static_.paths[source][target].empty()) ++static_.nodeDegree[source], ++static_.nodeDegree[target];
    }
  }
  static_.adjacentLinks.assign(graph.edges.size(), {});
  for (NodeId node = 0; node < graph.nodeCount; ++node) {
    for (EdgeId from : graph.adjacency[node]) for (EdgeId to : graph.adjacency[node]) {
      if (from != to) static_.adjacentLinks[from].push_back({to, node});
    }
  }
}

RFDMetric RFDAlgorithm::measure(const SubstrateGraph& graph) const {
  if (static_.nodes != graph.nodeCount) throw std::runtime_error("RFD algorithm was not initialized for this substrate");
  const int nodes = graph.nodeCount, edges = static_cast<int>(graph.edges.size());
  std::vector<double> nodeRatio(nodes), edgeRatio(edges), nodeConnectivity(nodes, 0.0), edgeConnectivity(edges, 0.0);
  for (int node = 0; node < nodes; ++node) nodeRatio[node] = graph.cpuCapacity[node] > 0.0 ? std::clamp(graph.cpuResidual[node] / graph.cpuCapacity[node], 0.0, 1.0) : 0.0;
  for (int edge = 0; edge < edges; ++edge) edgeRatio[edge] = graph.edges[edge].capacity > 0.0 ? std::clamp(graph.edges[edge].residual / graph.edges[edge].capacity, 0.0, 1.0) : 0.0;
  for (int source = 0; source < nodes; ++source) {
    if (graph.cpuCapacity[source] <= kEps) continue;
    for (int target = 0; target < nodes; ++target) {
      if (graph.cpuCapacity[target] <= kEps) continue;
      if (static_.nodeDegree[target] == 0) continue;
      const auto& paths = static_.paths[source][target];
      if (source == target || paths.empty()) continue;
      double average = 0.0;
      for (const auto& path : paths) average += bottleneckRatio(graph, path);
      average /= static_cast<double>(paths.size());
      // Eq. (2)--(4): d_j^tau normalizes the connectivity being computed for j.
      nodeConnectivity[target] += nodeRatio[source] * average / static_cast<double>(static_.nodeDegree[target]);
    }
  }
  for (int source = 0; source < edges; ++source) {
    for (const auto& [target, common] : static_.adjacentLinks[source]) {
      const size_t degree = static_.adjacentLinks[target].size();
      if (degree == 0) continue;
      // Eq. (6)--(8): d_j^l belongs to the target link, not its neighbor.
      // A zero-capacity CPU node is a forwarding-only switch, so CPU must not
      // constrain link connectivity at that common endpoint.
      const double commonRatio = graph.cpuCapacity[common] <= kEps
          ? 1.0 : nodeRatio[common];
      edgeConnectivity[target] += edgeRatio[source] * commonRatio / static_cast<double>(degree);
    }
  }
  RFDMetric output;
  for (int node = 0; node < nodes; ++node) {
    if (graph.cpuCapacity[node] <= kEps) continue;
    const double rfd = 1.0 - std::clamp(nodeConnectivity[node], 0.0, 1.0);
    output.node += rfd;
    output.fragmentationCost += std::max(0.0, graph.cpuResidual[node]) * rfd;
  }
  for (int edge = 0; edge < edges; ++edge) { const double rfd = 1.0 - std::clamp(edgeConnectivity[edge], 0.0, 1.0); output.link += rfd; output.fragmentationCost += std::max(0.0, graph.edges[edge].residual) * rfd; }
  output.total = output.node + output.link;
  return output;
}

std::optional<Embedding> RFDAlgorithm::embed(
    SubstrateGraph& graph, const VirtualNetworkRequest& request,
    bool firstRequest, RejectReason* reason) {
  if (reason) *reason = RejectReason::Route;
  if (request.nodes.size() > static_cast<size_t>(graph.nodeCount)) { if (reason) *reason = RejectReason::Cpu; return std::nullopt; }
  const auto start = graph.snapshot();
  const auto incident = incidence(request);
  const int root = rootNode(request, incident);
  const auto order = visitOrder(request, incident, root, variant_);
  MappingState state{std::vector<NodeId>(request.nodes.size(), -1), std::vector<std::vector<EdgeId>>(request.edges.size()), std::vector<char>(graph.nodeCount, 0)};

  const auto refine = [&]() -> std::optional<NodeId> {
    std::vector<std::vector<NodeId>> feasible(request.nodes.size());
    for (int node = 0; node < static_cast<int>(request.nodes.size()); ++node) for (NodeId candidate = 0; candidate < graph.nodeCount; ++candidate) if (graph.cpuResidual[candidate] + kEps >= request.nodes[node].cpu) feasible[node].push_back(candidate);
    for (const auto& candidates : feasible) if (candidates.empty()) return std::nullopt;
    std::mt19937_64 random(static_cast<std::uint64_t>(request.id + 1) * 0x9e3779b97f4a7c15ULL);
    std::uniform_int_distribution<int> choose(0, static_cast<int>(feasible[root].size()) - 1);
    NodeId current = feasible[root][choose(random)];
    std::vector<int> neighbors; for (int edge : incident[root]) neighbors.push_back(request.edges[edge].a == root ? request.edges[edge].b : request.edges[edge].a);
    std::sort(neighbors.begin(), neighbors.end()); neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    for (int iteration = 0; iteration < 2048; ++iteration) {
      const NodeId previous = current; std::vector<NodeId> closest;
      for (int virtualNode : neighbors) { NodeId best = -1; int distance = std::numeric_limits<int>::max(); for (NodeId candidate : feasible[virtualNode]) { const int nextDistance = hopDistance(graph, current, candidate); if (nextDistance < distance || (nextDistance == distance && (best < 0 || candidate < best))) best = candidate, distance = nextDistance; } closest.push_back(best); }
      NodeId best = -1; long long sumBest = std::numeric_limits<long long>::max();
      for (NodeId candidate : feasible[root]) { long long sum = 0; for (NodeId neighbor : closest) sum += hopDistance(graph, candidate, neighbor); if (sum < sumBest || (sum == sumBest && (best < 0 || candidate < best))) best = candidate, sumBest = sum; }
      current = best;
      if (hopDistance(graph, previous, current) < config_.maxD) break;
    }
    return current;
  };
  // Algorithm 3 is defined for the first VNR while the substrate is unused.
  const std::optional<NodeId> forcedRoot = firstRequest ? refine() : std::nullopt;
  if (firstRequest && !forcedRoot) { if (reason) *reason = RejectReason::Cpu; return std::nullopt; }

  const auto routePending = [&](int virtualNode, double* added) {
    std::vector<int> pending;
    for (int edge : incident[virtualNode]) { const int neighbor = request.edges[edge].a == virtualNode ? request.edges[edge].b : request.edges[edge].a; if (state.paths[edge].empty() && state.nodes[neighbor] >= 0) pending.push_back(edge); }
    std::sort(pending.begin(), pending.end());
    double total = 0.0, baseline = measure(graph).fragmentationCost;
    for (int edge : pending) {
      const auto& virtualEdge = request.edges[edge]; const int neighbor = virtualEdge.a == virtualNode ? virtualEdge.b : virtualEdge.a;
      const auto paths = yenPaths(
          graph, state.nodes[neighbor], state.nodes[virtualNode],
          virtualEdge.bandwidth, std::max(0, graph.nodeCount - 1), config_.ksp);
      if (paths.empty()) return false;
      Path selected; double best = std::numeric_limits<double>::infinity();
      for (const auto& path : paths) {
        const auto snapshot = graph.snapshot(); graph.allocateBandwidth(path.edges, virtualEdge.bandwidth);
        const double score = linkSelectionScore(
            config_, virtualEdge.bandwidth, path, baseline,
            measure(graph).fragmentationCost);
        graph.restore(snapshot);
        if (score < best - kEps || (std::abs(score - best) <= kEps && (selected.edges.empty() || path.edges < selected.edges))) selected = path, best = score;
      }
      graph.allocateBandwidth(selected.edges, virtualEdge.bandwidth); state.paths[edge] = selected.edges;
      baseline = measure(graph).fragmentationCost;
      // The path's RFD delta chose this path; final candidate accounting adds
      // only its resource cost because the complete RFD delta is added once.
      total += virtualEdge.bandwidth * static_cast<double>(selected.edges.size());
    }
    if (added) *added = total;
    return true;
  };
  const auto evaluate = [&](int virtualNode, NodeId substrate) -> std::optional<Candidate> {
    if (state.used[substrate] || graph.cpuResidual[substrate] + kEps < request.nodes[virtualNode].cpu) return std::nullopt;
    const auto snapshot = graph.snapshot(); const MappingState backup = state; const double before = measure(graph).fragmentationCost;
    graph.allocateCpu(substrate, request.nodes[virtualNode].cpu); state.nodes[virtualNode] = substrate; state.used[substrate] = 1;
    double pathCost = 0.0;
    const bool ok = routePending(virtualNode, &pathCost);
    const double score = candidateCost(
        config_, request.nodes[virtualNode].cpu, pathCost, before,
        measure(graph).fragmentationCost);
    graph.restore(snapshot); state = backup;
    return ok ? std::optional<Candidate>({substrate, score}) : std::nullopt;
  };
  const auto commit = [&](int virtualNode, NodeId substrate) {
    if (state.used[substrate] || !graph.allocateCpu(substrate, request.nodes[virtualNode].cpu)) return false;
    state.nodes[virtualNode] = substrate; state.used[substrate] = 1;
    return routePending(virtualNode, nullptr);
  };
  int backtracks = 0; bool limited = false;
  std::function<bool(int)> place = [&](int index) {
    if (index == static_cast<int>(order.size())) return true;
    const int virtualNode = order[index]; std::vector<Candidate> candidates;
    if (firstRequest && index == 0) { if (const auto candidate = evaluate(virtualNode, *forcedRoot)) candidates.push_back(*candidate); }
    else for (NodeId substrate = 0; substrate < graph.nodeCount; ++substrate) if (const auto candidate = evaluate(virtualNode, substrate)) candidates.push_back(*candidate);
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) { return a.score < b.score - kEps || (std::abs(a.score - b.score) <= kEps && a.substrate < b.substrate); });
    for (const auto& candidate : candidates) {
      const auto snapshot = graph.snapshot(); const MappingState backup = state;
      if (commit(virtualNode, candidate.substrate) && place(index + 1)) return true;
      graph.restore(snapshot); state = backup;
      if (++backtracks > config_.delta) { limited = true; return false; }
    }
    return false;
  };
  if (!place(0)) { graph.restore(start); if (reason) *reason = limited ? RejectReason::Backtrack : RejectReason::Route; return std::nullopt; }
  Embedding output; output.departure = request.departure; output.virtualToSubstrate = state.nodes; output.linkPaths = state.paths;
  for (const auto& node : request.nodes) { output.cpuAllocations.push_back(node.cpu); output.revenue += node.cpu; output.cost += node.cpu; }
  for (int edge = 0; edge < static_cast<int>(request.edges.size()); ++edge) { output.revenue += request.edges[edge].bandwidth; output.cost += request.edges[edge].bandwidth * static_cast<double>(state.paths[edge].size()); }
  if (reason) *reason = RejectReason::None;
  return output;
}

}  // namespace vne
