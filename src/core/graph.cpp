#include "core/graph.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>

namespace vne {
namespace {

std::vector<std::string> lines(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open " + path.string());
  std::vector<std::string> output;
  std::string line;
  while (std::getline(input, line)) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    const auto begin = line.find_first_not_of(" \t\r\n");
    if (begin != std::string::npos) output.push_back(line.substr(begin, line.find_last_not_of(" \t\r\n") - begin + 1));
  }
  return output;
}
std::vector<std::string> split(const std::string& line) {
  std::istringstream input(line);
  std::vector<std::string> output;
  std::string value;
  while (input >> value) output.push_back(value);
  return output;
}
int toInt(const std::string& value, const char* name) {
  size_t end = 0;
  try { const int result = std::stoi(value, &end); if (end == value.size()) return result; } catch (...) {}
  throw std::runtime_error(std::string("invalid ") + name + ": " + value);
}
double toDouble(const std::string& value, const char* name) {
  size_t end = 0;
  try { const double result = std::stod(value, &end); if (end == value.size() && std::isfinite(result)) return result; } catch (...) {}
  throw std::runtime_error(std::string("invalid ") + name + ": " + value);
}

}  // namespace

void SubstrateGraph::initialize(int nodes) {
  if (nodes <= 0) throw std::runtime_error("substrate node count must be positive");
  nodeCount = nodes;
  cpuCapacity.assign(nodes, 0.0);
  cpuResidual.assign(nodes, 0.0);
  nodeLocations.assign(nodes, {});
  adjacency.assign(nodes, {});
  edgeIndex.assign(nodes, std::vector<EdgeId>(nodes, -1));
  edges.clear();
}
void SubstrateGraph::addEdge(NodeId u, NodeId v, double capacity) {
  if (u < 0 || v < 0 || u >= nodeCount || v >= nodeCount || u == v || capacity < 0.0 || edgeIndex[u][v] >= 0) throw std::runtime_error("invalid substrate edge");
  const EdgeId id = static_cast<EdgeId>(edges.size());
  edges.push_back({u, v, capacity, capacity});
  adjacency[u].push_back(id);
  adjacency[v].push_back(id);
  edgeIndex[u][v] = edgeIndex[v][u] = id;
}
NodeId SubstrateGraph::other(EdgeId edge, NodeId node) const { return edges.at(edge).u == node ? edges.at(edge).v : edges.at(edge).u; }
void SubstrateGraph::sortAdjacency() {
  for (NodeId node = 0; node < nodeCount; ++node) std::sort(adjacency[node].begin(), adjacency[node].end(), [&](EdgeId a, EdgeId b) {
    const NodeId va = other(a, node), vb = other(b, node);
    return va != vb ? va < vb : a < b;
  });
}
SubstrateSnapshot SubstrateGraph::snapshot() const {
  SubstrateSnapshot value{cpuResidual, {}};
  value.edge.reserve(edges.size());
  for (const auto& edge : edges) value.edge.push_back(edge.residual);
  return value;
}
void SubstrateGraph::restore(const SubstrateSnapshot& value) {
  cpuResidual = value.cpu;
  for (size_t i = 0; i < edges.size(); ++i) edges[i].residual = value.edge.at(i);
}
bool SubstrateGraph::allocateCpu(NodeId node, double amount) {
  if (!std::isfinite(amount) || amount < 0.0) throw std::runtime_error("invalid CPU allocation");
  if (cpuResidual.at(node) + kEps < amount) return false;
  cpuResidual[node] = std::max(0.0, cpuResidual[node] - amount);
  return true;
}
void SubstrateGraph::releaseCpu(NodeId node, double amount) {
  if (!std::isfinite(amount) || amount < 0.0) throw std::runtime_error("invalid CPU release");
  cpuResidual.at(node) = std::min(cpuCapacity.at(node), cpuResidual.at(node) + amount);
}
void SubstrateGraph::allocateBandwidth(const std::vector<EdgeId>& path, double amount) {
  if (!std::isfinite(amount) || amount < 0.0) throw std::runtime_error("invalid bandwidth allocation");
  for (EdgeId edge : path) if (edges.at(edge).residual + kEps < amount) throw std::runtime_error("insufficient bandwidth");
  for (EdgeId edge : path) edges.at(edge).residual = std::max(0.0, edges.at(edge).residual - amount);
}
bool SubstrateGraph::tryAllocateBandwidth(
    const std::vector<std::pair<std::vector<EdgeId>, double>>& allocations,
    double relativeTolerance) {
  if (!std::isfinite(relativeTolerance) || relativeTolerance < 0.0)
    throw std::runtime_error("invalid bandwidth allocation tolerance");
  std::vector<double> totals(edges.size(), 0.0);
  for (const auto& allocation : allocations) {
    const double amount = allocation.second;
    if (!std::isfinite(amount) || amount < 0.0)
      throw std::runtime_error("invalid bandwidth allocation");
    for (EdgeId edge : allocation.first) totals.at(edge) += amount;
  }
  for (size_t edge = 0; edge < edges.size(); ++edge) {
    if (!std::isfinite(totals[edge])) throw std::runtime_error("invalid bandwidth allocation");
    const double tolerance = relativeTolerance * std::max(1.0, edges[edge].capacity);
    if (totals[edge] > edges[edge].residual + tolerance) return false;
  }
  for (size_t edge = 0; edge < edges.size(); ++edge)
    edges[edge].residual = std::max(0.0, edges[edge].residual - totals[edge]);
  return true;
}
void SubstrateGraph::releaseBandwidth(const std::vector<EdgeId>& path, double amount) {
  if (!std::isfinite(amount) || amount < 0.0) throw std::runtime_error("invalid bandwidth release");
  for (EdgeId edge : path) edges.at(edge).residual = std::min(edges.at(edge).capacity, edges.at(edge).residual + amount);
}

SubstrateGraph loadSubstrate(const std::filesystem::path& path) {
  const auto input = lines(path);
  if (input.empty()) throw std::runtime_error("empty substrate file");
  const auto header = split(input.front());
  if (header.size() != 3 || header[0] != "SUBSTRATE") throw std::runtime_error("invalid substrate header");
  const int nodes = toInt(header[1], "node count"), edgeCount = toInt(header[2], "edge count");
  SubstrateGraph graph;
  graph.initialize(nodes);
  int position = 1;
  for (int node = 0; node < nodes; ++node, ++position) {
    if (position >= static_cast<int>(input.size())) throw std::runtime_error("missing NODE record");
    const auto item = split(input[position]);
    if ((item.size() != 3 && item.size() != 5) || item[0] != "NODE" || toInt(item[1], "node id") != node) throw std::runtime_error("invalid NODE record");
    const double cpu = toDouble(item[2], "node CPU");
    if (cpu < 0.0) throw std::runtime_error("negative node CPU");
    graph.cpuCapacity[node] = graph.cpuResidual[node] = cpu;
    if (item.size() == 5) {
      graph.nodeLocations[node] = {true, toDouble(item[3], "node x"), toDouble(item[4], "node y")};
    }
  }
  for (; position < static_cast<int>(input.size()); ++position) {
    const auto item = split(input[position]);
    if (item.size() != 4 || item[0] != "EDGE") throw std::runtime_error("invalid EDGE record");
    graph.addEdge(toInt(item[1], "edge endpoint"), toInt(item[2], "edge endpoint"), toDouble(item[3], "edge bandwidth"));
  }
  if (static_cast<int>(graph.edges.size()) != edgeCount) throw std::runtime_error("substrate edge count mismatch");
  graph.sortAdjacency();
  std::vector<char> seen(nodes, 0); std::queue<NodeId> queue; queue.push(0); seen[0] = 1;
  while (!queue.empty()) { const NodeId u = queue.front(); queue.pop(); for (EdgeId edge : graph.adjacency[u]) { const NodeId v = graph.other(edge, u); if (!seen[v]) { seen[v] = 1; queue.push(v); } } }
  if (std::find(seen.begin(), seen.end(), 0) != seen.end()) throw std::runtime_error("substrate graph must be connected");
  return graph;
}

}  // namespace vne
