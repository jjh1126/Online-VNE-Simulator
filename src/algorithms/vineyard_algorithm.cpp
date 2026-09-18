#include "algorithms/vineyard_algorithm.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include <glpk.h>

namespace vne {
namespace {

// Match GLPK's default bound feasibility tolerance.  Flow allocations are
// checked per physical edge in one transaction before the graph is mutated.
constexpr double kLpBoundTolerance = 1e-7;

struct AugEdge { int a; int b; int substrateEdge; int virtualNode; int substrateNode; double capacity; };
struct Scores { std::vector<std::vector<double>> score, x, flow; };
struct Arc { int to; EdgeId edge; int direction; double flow; };

double distance(double ax, double ay, double bx, double by) {
  const double dx = ax - bx, dy = ay - by;
  return std::sqrt(dx * dx + dy * dy);
}

void requireLocations(const SubstrateGraph& graph, const VirtualNetworkRequest& request) {
  for (const auto& location : graph.nodeLocations)
    if (!location.present) throw std::runtime_error("D-ViNE-LB/R-ViNE-LB require NODE <id> <cpu> <x> <y> records");
  for (const auto& node : request.nodes)
    if (!node.hasLocation) throw std::runtime_error("D-ViNE-LB/R-ViNE-LB require V <id> <cpu> <x> <y> <radius> records");
}

std::vector<std::vector<NodeId>> candidates(const SubstrateGraph& graph, const VirtualNetworkRequest& request) {
  std::vector<std::vector<NodeId>> out(request.nodes.size());
  for (int v = 0; v < static_cast<int>(request.nodes.size()); ++v) {
    for (NodeId s = 0; s < graph.nodeCount; ++s) {
      if (graph.cpuResidual[s] + kEps >= request.nodes[v].cpu &&
          distance(request.nodes[v].x, request.nodes[v].y, graph.nodeLocations[s].x, graph.nodeLocations[s].y) <= request.nodes[v].radius + kEps)
        out[v].push_back(s);
    }
  }
  return out;
}

void addRow(glp_prob* lp, int type, double lower, double upper, const std::vector<std::pair<int, double>>& values) {
  const int row = glp_add_rows(lp, 1);
  glp_set_row_bnds(lp, row, type, lower, upper);
  std::vector<int> indices(1); std::vector<double> coefficients(1);
  for (const auto& value : values) { indices.push_back(value.first); coefficients.push_back(value.second); }
  glp_set_mat_row(lp, row, static_cast<int>(values.size()), indices.data(), coefficients.data());
}

struct GlpkDeleter { void operator()(glp_prob* value) const { glp_delete_prob(value); } };
using GlpkProblem = std::unique_ptr<glp_prob, GlpkDeleter>;

bool solve(glp_prob* lp, int timeoutMs, bool presolve, VineSolverDiagnostic* diagnostic) {
  glp_smcp options; glp_init_smcp(&options); options.msg_lev = GLP_MSG_OFF; options.tm_lim = timeoutMs;
  options.presolve = presolve ? GLP_ON : GLP_OFF;
  const auto started = std::chrono::steady_clock::now();
  const int result = glp_simplex(lp, &options);
  const int status = glp_get_status(lp);
  if (diagnostic) {
    diagnostic->returnCode = result;
    diagnostic->status = status;
    diagnostic->iterations = glp_get_it_cnt(lp);
    diagnostic->runtimeMs = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
  }
  return result == 0 && (status == GLP_OPT || status == GLP_FEAS);
}

std::optional<Scores> relaxation(
    const SubstrateGraph& graph, const VirtualNetworkRequest& request,
    const std::vector<std::vector<NodeId>>& candidateSets, int timeoutMs,
    bool presolve, VineSolverDiagnostic* diagnostic) {
  const int n = graph.nodeCount, nv = static_cast<int>(request.nodes.size()), nk = static_cast<int>(request.edges.size());
  std::vector<AugEdge> edges;
  for (EdgeId e = 0; e < static_cast<int>(graph.edges.size()); ++e)
    if (graph.edges[e].residual > kEps)
      edges.push_back({graph.edges[e].u, graph.edges[e].v, e, -1, -1, graph.edges[e].residual});
  for (int v = 0; v < nv; ++v) for (NodeId s : candidateSets[v]) edges.push_back({n + v, s, -1, v, s, 9999.0});
  const int ne = static_cast<int>(edges.size()), na = n + nv;
  GlpkProblem lp(glp_create_prob()); glp_set_obj_dir(lp.get(), GLP_MIN);
  const int xOffset = 1;
  const int fOffset = xOffset + ne;
  const int totalColumns = ne + ne * 2 * nk;
  glp_add_cols(lp.get(), totalColumns);
  auto f = [=](int e, int dir, int k) { return fOffset + (e * 2 + dir) * nk + k; };
  for (int e = 0; e < ne; ++e) {
    glp_set_col_bnds(lp.get(), xOffset + e, GLP_DB, 0.0, 1.0);
    if (edges[e].substrateEdge < 0) glp_set_obj_coef(lp.get(), xOffset + e, request.nodes[edges[e].virtualNode].cpu / (graph.cpuResidual[edges[e].substrateNode] + 1e-6));
    for (int dir = 0; dir < 2; ++dir) for (int k = 0; k < nk; ++k) {
      glp_set_col_bnds(lp.get(), f(e, dir, k), GLP_DB, 0.0, edges[e].capacity);
      if (edges[e].substrateEdge >= 0) glp_set_obj_coef(lp.get(), f(e, dir, k), 1.0 / (edges[e].capacity + 1e-6));
    }
  }
  for (int e = 0; e < ne; ++e) {
    std::vector<std::pair<int, double>> row{{xOffset + e, -edges[e].capacity}};
    for (int d = 0; d < 2; ++d) for (int k = 0; k < nk; ++k) row.push_back({f(e, d, k), 1.0});
    addRow(lp.get(), GLP_UP, 0.0, 0.0, row);
  }
  for (int k = 0; k < nk; ++k) for (int node = 0; node < na; ++node) {
    std::vector<std::pair<int, double>> row;
    for (int e = 0; e < ne; ++e) {
      if (edges[e].a == node) { row.push_back({f(e, 0, k), 1.0}); row.push_back({f(e, 1, k), -1.0}); }
      if (edges[e].b == node) { row.push_back({f(e, 1, k), 1.0}); row.push_back({f(e, 0, k), -1.0}); }
    }
    double rhs = node == n + request.edges[k].a ? request.edges[k].bandwidth : node == n + request.edges[k].b ? -request.edges[k].bandwidth : 0.0;
    addRow(lp.get(), GLP_FX, rhs, rhs, row);
  }
  for (int v = 0; v < nv; ++v) { std::vector<std::pair<int, double>> row; for (int e = 0; e < ne; ++e) if (edges[e].virtualNode == v) row.push_back({xOffset + e, 1.0}); addRow(lp.get(), GLP_FX, 1.0, 1.0, row); }
  for (NodeId s = 0; s < n; ++s) {
    std::vector<std::pair<int, double>> one, cpu;
    for (int e = 0; e < ne; ++e) if (edges[e].substrateNode == s) { one.push_back({xOffset + e, 1.0}); cpu.push_back({xOffset + e, request.nodes[edges[e].virtualNode].cpu}); }
    addRow(lp.get(), GLP_UP, 0.0, 1.0, one); addRow(lp.get(), GLP_UP, 0.0, graph.cpuResidual[s], cpu);
  }
  if (!solve(lp.get(), timeoutMs, presolve, diagnostic)) return std::nullopt;
  Scores out; out.score.assign(nv, std::vector<double>(n)); out.x.assign(nv, std::vector<double>(n)); out.flow.assign(nv, std::vector<double>(n));
  for (int e = 0; e < ne; ++e) if (edges[e].virtualNode >= 0) {
    double flow = 0.0; for (int d = 0; d < 2; ++d) for (int k = 0; k < nk; ++k) flow += glp_get_col_prim(lp.get(), f(e, d, k));
    const double x = glp_get_col_prim(lp.get(), xOffset + e); out.x[edges[e].virtualNode][edges[e].substrateNode] = x; out.flow[edges[e].virtualNode][edges[e].substrateNode] = flow; out.score[edges[e].virtualNode][edges[e].substrateNode] = x * flow;
  }
  return out;
}

std::optional<std::vector<std::vector<std::pair<std::vector<EdgeId>, double>>>> mcf(
    const SubstrateGraph& graph, const VirtualNetworkRequest& request, const std::vector<NodeId>& mapped,
    int timeoutMs, bool presolve, VineSolverDiagnostic* diagnostic) {
  const int nk = static_cast<int>(request.edges.size()), n = graph.nodeCount;
  std::vector<EdgeId> physicalEdges;
  for (EdgeId edge = 0; edge < static_cast<EdgeId>(graph.edges.size()); ++edge)
    if (graph.edges[edge].residual > kEps) physicalEdges.push_back(edge);
  const int ne = static_cast<int>(physicalEdges.size());
  GlpkProblem lp(glp_create_prob()); glp_set_obj_dir(lp.get(), GLP_MIN); glp_add_cols(lp.get(), ne * 2 * nk);
  auto f = [=](int e, int dir, int k) { return 1 + (e * 2 + dir) * nk + k; };
  for (int e = 0; e < ne; ++e) for (int d = 0; d < 2; ++d) for (int k = 0; k < nk; ++k) { glp_set_col_bnds(lp.get(), f(e,d,k), GLP_DB, 0.0, graph.edges[physicalEdges[e]].residual); glp_set_obj_coef(lp.get(), f(e,d,k), 1.0); }
  for (int e = 0; e < ne; ++e) { std::vector<std::pair<int,double>> row; for (int d = 0; d < 2; ++d) for (int k = 0; k < nk; ++k) row.push_back({f(e,d,k),1.0}); addRow(lp.get(), GLP_UP, 0.0, graph.edges[physicalEdges[e]].residual, row); }
  for (int k = 0; k < nk; ++k) for (NodeId node = 0; node < n; ++node) {
    std::vector<std::pair<int,double>> row;
    for (int e = 0; e < ne; ++e) { const EdgeId edgeId = physicalEdges[e]; if (graph.edges[edgeId].u != node && graph.edges[edgeId].v != node) continue; const bool forward = graph.edges[edgeId].u == node; row.push_back({f(e, forward ? 0 : 1, k), 1.0}); row.push_back({f(e, forward ? 1 : 0, k), -1.0}); }
    const auto& edge = request.edges[k]; const double rhs = node == mapped[edge.a] ? edge.bandwidth : node == mapped[edge.b] ? -edge.bandwidth : 0.0; addRow(lp.get(), GLP_FX, rhs, rhs, row);
  }
  if (!solve(lp.get(), timeoutMs, presolve, diagnostic)) return std::nullopt;
  std::vector<std::vector<std::pair<std::vector<EdgeId>, double>>> result(nk);
  for (int k = 0; k < nk; ++k) {
    std::vector<std::vector<Arc>> flow(n);
    for (int e = 0; e < ne; ++e) for (int d = 0; d < 2; ++d) { const double value = glp_get_col_prim(lp.get(), f(e,d,k)); if (value > 1e-7) { const EdgeId edgeId = physicalEdges[e]; const NodeId from = d == 0 ? graph.edges[edgeId].u : graph.edges[edgeId].v; const NodeId to = graph.other(edgeId, from); flow[from].push_back({to,edgeId,d,value}); } }
    double remaining = request.edges[k].bandwidth; const NodeId src = mapped[request.edges[k].a], dst = mapped[request.edges[k].b];
    while (remaining > 1e-6) {
      std::vector<NodeId> parent(n, -1); std::vector<int> parentArc(n, -1); std::queue<NodeId> q; parent[src] = src; q.push(src);
      while (!q.empty() && parent[dst] < 0) { NodeId u = q.front(); q.pop(); for (int i = 0; i < static_cast<int>(flow[u].size()); ++i) if (flow[u][i].flow > 1e-7 && parent[flow[u][i].to] < 0) { parent[flow[u][i].to] = u; parentArc[flow[u][i].to] = i; q.push(flow[u][i].to); } }
      if (parent[dst] < 0) return std::nullopt;
      double amount = remaining; std::vector<EdgeId> path; for (NodeId cur = dst; cur != src; cur = parent[cur]) { const Arc& arc = flow[parent[cur]][parentArc[cur]]; amount = std::min(amount, arc.flow); path.push_back(arc.edge); }
      if (amount <= 1e-7) return std::nullopt;
      std::reverse(path.begin(), path.end());
      for (NodeId cur = dst; cur != src; cur = parent[cur]) flow[parent[cur]][parentArc[cur]].flow -= amount;
      result[k].push_back({std::move(path), amount}); remaining -= amount;
    }
  }
  return result;
}

uint32_t seed(unsigned int base, int substrateSeed, int requestSeed, int requestId) {
  uint64_t value = base; value = value * 1315423911ULL + static_cast<uint32_t>(substrateSeed); value = value * 2654435761ULL + static_cast<uint32_t>(requestSeed); value = value * 2246822519ULL + static_cast<uint32_t>(requestId); return static_cast<uint32_t>(value);
}

std::string failureStage(const char* prefix, const VineSolverDiagnostic& solver) {
  if (solver.returnCode == GLP_ETMLIM) return std::string(prefix) + "_timeout";
  if (solver.returnCode == GLP_EBOUND) return std::string(prefix) + "_invalid_bounds";
  if (solver.status == GLP_NOFEAS || solver.status == GLP_INFEAS)
    return std::string(prefix) + "_infeasible";
  return std::string(prefix) + "_solver_failure";
}

}  // namespace

VineyardAlgorithm::VineyardAlgorithm(VineyardMode mode, AlgorithmConfig config, int substrateSeed, int requestSeed)
    : mode_(mode), config_(config), substrateSeed_(substrateSeed), requestSeed_(requestSeed) {}

std::optional<Embedding> VineyardAlgorithm::embed(SubstrateGraph& graph, const VirtualNetworkRequest& request, bool, RejectReason* reason) {
  lastDiagnostic_ = {};
  lastDiagnostic_.enabled = config_.vineDiagnostics;
  lastDiagnostic_.stage = config_.vineDiagnostics ? "started" : "disabled";
  if (reason) *reason = RejectReason::Route;
  requireLocations(graph, request);
  if (request.nodes.size() > static_cast<size_t>(graph.nodeCount)) {
    if (reason) *reason = RejectReason::Cpu;
    if (config_.vineDiagnostics) lastDiagnostic_.stage = "too_many_virtual_nodes";
    return std::nullopt;
  }
  const auto snapshot = graph.snapshot();
  const auto candidateSets = candidates(graph, request);
  if (config_.vineDiagnostics && !candidateSets.empty()) {
    lastDiagnostic_.candidateMin = static_cast<int>(candidateSets.front().size());
    for (const auto& set : candidateSets) {
      lastDiagnostic_.candidateMin = std::min(lastDiagnostic_.candidateMin, static_cast<int>(set.size()));
      lastDiagnostic_.candidateMax = std::max(lastDiagnostic_.candidateMax, static_cast<int>(set.size()));
    }
  }
  for (const auto& set : candidateSets) if (set.empty()) {
    if (reason) *reason = RejectReason::Cpu;
    if (config_.vineDiagnostics) lastDiagnostic_.stage = "candidate_empty";
    return std::nullopt;
  }
  const auto score = relaxation(graph, request, candidateSets, config_.vineTimeoutMs, false,
      config_.vineDiagnostics ? &lastDiagnostic_.cnlm : nullptr);
  if (!score) {
    if (config_.vineDiagnostics) {
      lastDiagnostic_.stage = failureStage("cnlm", lastDiagnostic_.cnlm);
      (void)relaxation(graph, request, candidateSets, config_.vineTimeoutMs, true,
          &lastDiagnostic_.cnlmPresolve);
    }
    return std::nullopt;
  }
  std::mt19937 random(seed(config_.vineRandomSeed, substrateSeed_, requestSeed_, request.id));
  std::vector<NodeId> mapped(request.nodes.size(), -1); std::vector<char> used(graph.nodeCount);
  for (int v = 0; v < static_cast<int>(request.nodes.size()); ++v) {
    std::vector<NodeId> available;
    for (NodeId s : candidateSets[v]) {
      if (!used[s]) available.push_back(s);
      else if (config_.vineDiagnostics) ++lastDiagnostic_.unavailableCandidates;
    }
    if (available.empty()) {
      if (reason) *reason = RejectReason::Cpu;
      if (config_.vineDiagnostics) lastDiagnostic_.stage = "rounding_no_unused_candidate";
      return std::nullopt;
    }
    NodeId selected = available.front();
    if (mode_ == VineyardMode::Deterministic) {
      for (NodeId s : available) {
        const bool betterScore = score->score[v][s] > score->score[v][selected] + kEps;
        const bool tiedScore = std::abs(score->score[v][s] - score->score[v][selected]) <= kEps;
        const bool betterX = score->x[v][s] > score->x[v][selected] + kEps;
        const bool tiedX = std::abs(score->x[v][s] - score->x[v][selected]) <= kEps;
        if (betterScore || (tiedScore && (betterX || (tiedX && s < selected)))) selected = s;
      }
    } else {
      std::vector<double> weights; double total = 0.0;
      for (NodeId s : available) { weights.push_back(std::max(0.0, score->score[v][s])); total += weights.back(); }
      if (total <= kEps) { total = 0.0; weights.clear(); for (NodeId s : available) { weights.push_back(std::max(0.0, score->x[v][s])); total += weights.back(); } }
      if (total <= kEps) { std::uniform_int_distribution<int> distribution(0, static_cast<int>(available.size()) - 1); selected = available[distribution(random)]; }
      else selected = available[std::discrete_distribution<int>(weights.begin(), weights.end())(random)];
    }
    mapped[v] = selected; used[selected] = 1;
    if (config_.vineDiagnostics) {
      for (NodeId s : available) if (score->score[v][s] > 1e-12) ++lastDiagnostic_.positiveScoreCandidates;
      lastDiagnostic_.selectedXSum += score->x[v][selected];
      lastDiagnostic_.selectedFlowSum += score->flow[v][selected];
      lastDiagnostic_.selectedScoreSum += score->score[v][selected];
      int rank = 1; for (NodeId s : available) if (score->x[v][s] > score->x[v][selected] + kEps) ++rank;
      lastDiagnostic_.selectedXRankSum += rank;
    }
  }
  if (config_.vineDiagnostics) {
    std::vector<NodeId> shadow(request.nodes.size(), -1); std::vector<char> shadowUsed(graph.nodeCount);
    bool complete = true;
    for (int v = 0; v < static_cast<int>(request.nodes.size()); ++v) {
      NodeId selected = -1;
      for (NodeId s : candidateSets[v]) if (!shadowUsed[s] &&
          (selected < 0 || score->x[v][s] > score->x[v][selected] + kEps ||
           (std::abs(score->x[v][s] - score->x[v][selected]) <= kEps && s < selected))) selected = s;
      if (selected < 0) { complete = false; break; }
      shadow[v] = selected; shadowUsed[selected] = 1;
    }
    VineSolverDiagnostic ignored;
    lastDiagnostic_.shadowXFeasible = complete && mcf(graph, request, shadow,
        config_.vineTimeoutMs, true, &ignored).has_value() ? 1 : 0;
  }
  for (int v = 0; v < static_cast<int>(request.nodes.size()); ++v) if (!graph.allocateCpu(mapped[v], request.nodes[v].cpu)) {
    graph.restore(snapshot); if (reason) *reason = RejectReason::Cpu;
    if (config_.vineDiagnostics) lastDiagnostic_.stage = "rounding_cpu_allocation";
    return std::nullopt;
  }
  const auto flows = mcf(graph, request, mapped, config_.vineTimeoutMs, false,
      config_.vineDiagnostics ? &lastDiagnostic_.mcf : nullptr);
  if (!flows) {
    if (config_.vineDiagnostics) {
      lastDiagnostic_.stage = failureStage("mcf", lastDiagnostic_.mcf);
      (void)mcf(graph, request, mapped, config_.vineTimeoutMs, true, &lastDiagnostic_.mcfPresolve);
    }
    graph.restore(snapshot); return std::nullopt;
  }
  std::vector<std::pair<std::vector<EdgeId>, double>> allocations;
  for (const auto& edgeFlows : *flows)
    allocations.insert(allocations.end(), edgeFlows.begin(), edgeFlows.end());
  if (!graph.tryAllocateBandwidth(allocations, kLpBoundTolerance)) {
    graph.restore(snapshot);
    if (config_.vineDiagnostics) lastDiagnostic_.stage = "mcf_capacity_violation";
    return std::nullopt;
  }
  Embedding embedding; embedding.departure = request.departure; embedding.virtualToSubstrate = mapped; embedding.linkFlowPaths = *flows; embedding.cpuAllocations.reserve(request.nodes.size());
  for (const auto& node : request.nodes) { embedding.cpuAllocations.push_back(node.cpu); embedding.revenue += node.cpu; embedding.cost += node.cpu; }
  for (int edge = 0; edge < static_cast<int>(request.edges.size()); ++edge) { embedding.revenue += request.edges[edge].bandwidth; for (const auto& flow : (*flows)[edge]) embedding.cost += flow.second * static_cast<double>(flow.first.size()); }
  if (reason) *reason = RejectReason::None;
  if (config_.vineDiagnostics) lastDiagnostic_.stage = "accepted";
  return embedding;
}

}  // namespace vne
