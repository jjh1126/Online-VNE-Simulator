#pragma once

#include <filesystem>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "core/types.hpp"

namespace vne {

struct RequestRecord {
  std::string algorithm;
  int substrateSeed = 0;
  int requestSeed = 0;
  int requestId = -1;
  int arrival = 0;
  int departure = 0;
  int virtualNodes = 0;
  int virtualEdges = 0;
  bool accepted = false;
  RejectReason reason = RejectReason::Route;
  double revenue = 0.0;
  double cost = 0.0;
  double cpuUtilization = 0.0;
  double linkUtilization = 0.0;
  bool boundaryLossComponentsValid = false;
  double nodeBoundaryLossCost = 0.0;
  double linkBoundaryLossCost = 0.0;
  double boundaryLossObjective = 0.0;
  VineyardDiagnostic vineyard;
};

struct RunSummary {
  std::string algorithm;
  int substrateSeed = 0;
  int requestSeed = 0;
  int totalRequests = 0;
  int accepted = 0;
  int warmupRequests = 0;
  int warmupAccepted = 0;
  int warmupTime = 0;
  int measurementEnd = 0;
  double revenue = 0.0;
  double cost = 0.0;
  double averageCpuUtilization = 0.0;
  double averageLinkUtilization = 0.0;
  double runtimeMs = 0.0;
};

struct SimulationResult {
  RunSummary summary;
  std::vector<RequestRecord> records;
};

inline std::string compactDouble(double value) {
  std::ostringstream output;
  output << value;
  return output.str();
}

inline std::string fileSafeName(std::string value) {
  for (char& item : value) {
    if (item == '/' || item == '\\') item = '_';
  }
  return value;
}

inline std::string algorithmParameterSuffix(
    const std::string& algorithm, const AlgorithmConfig& config, int seedCount,
    int warmup, int measurementWindow, const std::string& sampleMode) {
  std::ostringstream output;
  if (algorithm == "VNE-RFD-B" || algorithm == "VNE-RFD-D") {
    output << "_tau" << config.tau << "_ksp" << config.ksp
           << "_delta" << config.delta << "_maxd" << config.maxD
           << "_a" << compactDouble(config.alpha)
           << "_b" << compactDouble(config.beta)
           << "_l" << compactDouble(config.lambda);
  } else if (algorithm == "VNE-BCP") {
    output << "_b" << compactDouble(config.beta)
           << "_g" << compactDouble(config.gamma);
  } else if (algorithm == "D-ViNE-LB") {
    output << "_timeout" << config.vineTimeoutMs;
  } else if (algorithm == "R-ViNE-LB") {
    output << "_timeout" << config.vineTimeoutMs
           << "_vseed" << config.vineRandomSeed;
  }
  output << "_" << sampleMode << "_samples" << seedCount;
  if (warmup > 0 || measurementWindow > 0)
    output << "_warmup" << warmup << "_window" << measurementWindow;
  return output.str();
}

inline void writeCompactCsv(
    const std::filesystem::path& resultsDirectory,
    std::vector<SimulationResult> results,
    const std::string& caseName,
    int seedCount,
    const AlgorithmConfig& config,
    int threads,
    int warmup,
    int measurementWindow,
    const std::string& sampleMode) {
  std::sort(results.begin(), results.end(), [](const SimulationResult& a, const SimulationResult& b) {
    const auto& x = a.summary; const auto& y = b.summary;
    return std::tie(x.algorithm, x.substrateSeed, x.requestSeed) < std::tie(y.algorithm, y.substrateSeed, y.requestSeed);
  });
  std::filesystem::create_directories(resultsDirectory);
  struct Aggregate { int runs = 0, total = 0, accepted = 0; double revenue = 0, cost = 0, cpu = 0, link = 0, runtime = 0; };
  std::unordered_map<std::string, Aggregate> aggregate;
  for (const auto& result : results) {
    auto& item = aggregate[result.summary.algorithm]; const auto& row = result.summary;
    ++item.runs; item.total += row.totalRequests; item.accepted += row.accepted; item.revenue += row.revenue; item.cost += row.cost; item.cpu += row.averageCpuUtilization; item.link += row.averageLinkUtilization; item.runtime += row.runtimeMs;
  }
  std::vector<std::string> names;
  for (const auto& [name, ignored] : aggregate) { (void)ignored; names.push_back(name); }
  std::sort(names.begin(), names.end());
  // Result directory names use Korean Standard Time regardless of the host's
  // configured timezone. KST is UTC+09:00 and does not observe DST.
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  now += 9 * 60 * 60;
  std::tm timestamp{};
#if defined(_WIN32)
  gmtime_s(&timestamp, &now);
#else
  gmtime_r(&now, &timestamp);
#endif
  std::ostringstream timestampText;
  timestampText << std::put_time(&timestamp, "%Y%m%d-%H%M%S");
  for (const auto& name : names) {
    const auto& item = aggregate.at(name); const double runs = static_cast<double>(item.runs);
    const std::string base = timestampText.str() + "_" + fileSafeName(caseName) + "_" + name +
        algorithmParameterSuffix(name, config, seedCount, warmup, measurementWindow, sampleMode);
    std::filesystem::path runDirectory = resultsDirectory / base;
    for (int suffix = 1; std::filesystem::exists(runDirectory); ++suffix) {
      runDirectory = resultsDirectory / (base + "_" + std::to_string(suffix));
    }
    std::filesystem::create_directories(runDirectory);
    std::cerr << "[RESULT] " << runDirectory.string() << '\n';
    std::ofstream summary(runDirectory / "summary.csv");
    std::ofstream requests(runDirectory / "requests.csv");
    std::ofstream runsOutput(runDirectory / "runs.csv");
    if (!summary || !requests || !runsOutput) throw std::runtime_error("cannot create compact CSV output");
    summary << std::setprecision(12);
    requests << std::setprecision(12);
    int warmupRequests = 0, warmupAccepted = 0;
    for (const auto& result : results) if (result.summary.algorithm == name) {
      warmupRequests += result.summary.warmupRequests;
      warmupAccepted += result.summary.warmupAccepted;
    }
    summary << "algo,case,seed_count,sample_mode,sample_count,threads,tau,ksp,delta,max_d,runs,warmup_time,measurement_window,measurement_end,warmup_requests,warmup_accepted,total_requests,accepted,rejected,total_revenue,total_cost,avg_cpu_utilization,avg_link_utilization,total_runtime_ms\n";
    summary << name << ',' << caseName << ',' << seedCount << ',' << sampleMode << ',' << seedCount << ',' << threads << ',' << config.tau << ',' << config.ksp << ',' << config.delta << ',' << config.maxD << ',' << item.runs << ',' << warmup << ',' << measurementWindow << ',' << warmup + measurementWindow << ',' << warmupRequests << ',' << warmupAccepted << ',' << item.total << ',' << item.accepted << ',' << item.total - item.accepted << ',' << item.revenue << ',' << item.cost << ',' << item.cpu / runs << ',' << item.link / runs << ',' << item.runtime << '\n';
    requests << "algo,substrate_seed,request_seed,request_id,arrival,departure,virtual_nodes,virtual_edges,accepted,reason,revenue,cost,cpu_utilization,link_utilization,boundary_loss_components_valid,node_boundary_loss_cost,link_boundary_loss_cost,boundary_loss_objective\n";
    runsOutput << std::setprecision(12);
    runsOutput << "algo,sample_id,substrate_seed,request_seed,total_requests,accepted,rejected,acceptance_ratio,total_revenue,total_cost,avg_cpu_utilization,avg_link_utilization,runtime_ms\n";
    for (const auto& result : results) {
      if (result.summary.algorithm != name) continue;
      const auto& run = result.summary;
      runsOutput << run.algorithm << ',' << run.substrateSeed << ',' << run.substrateSeed << ',' << run.requestSeed << ',' << run.totalRequests << ',' << run.accepted << ',' << run.totalRequests - run.accepted << ',' << (run.totalRequests ? static_cast<double>(run.accepted) / run.totalRequests : 0.0) << ',' << run.revenue << ',' << run.cost << ',' << run.averageCpuUtilization << ',' << run.averageLinkUtilization << ',' << run.runtimeMs << '\n';
      for (const auto& row : result.records) requests << row.algorithm << ',' << row.substrateSeed << ',' << row.requestSeed << ',' << row.requestId << ',' << row.arrival << ',' << row.departure << ',' << row.virtualNodes << ',' << row.virtualEdges << ',' << (row.accepted ? 1 : 0) << ',' << rejectReasonName(row.reason) << ',' << row.revenue << ',' << row.cost << ',' << row.cpuUtilization << ',' << row.linkUtilization << ',' << (row.boundaryLossComponentsValid ? 1 : 0) << ',' << row.nodeBoundaryLossCost << ',' << row.linkBoundaryLossCost << ',' << row.boundaryLossObjective << '\n';
    }
    bool hasVineyardDiagnostics = false;
    for (const auto& result : results) if (result.summary.algorithm == name)
      for (const auto& row : result.records) hasVineyardDiagnostics |= row.vineyard.enabled;
    if (hasVineyardDiagnostics) {
      std::ofstream diagnostics(runDirectory / "vineyard_diagnostics.csv");
      diagnostics << std::setprecision(12);
      diagnostics << "algo,substrate_seed,request_seed,request_id,arrival,accepted,stage,candidate_min,candidate_max,unavailable_candidates,positive_score_candidates,selected_x_sum,selected_flow_sum,selected_score_sum,selected_x_rank_sum,cnlm_rc,cnlm_status,cnlm_iterations,cnlm_ms,cnlm_presolve_rc,cnlm_presolve_status,cnlm_presolve_ms,mcf_rc,mcf_status,mcf_iterations,mcf_ms,mcf_presolve_rc,mcf_presolve_status,mcf_presolve_ms,shadow_x_feasible\n";
      for (const auto& result : results) {
        if (result.summary.algorithm != name) continue;
        for (const auto& row : result.records) if (row.vineyard.enabled) {
          const auto& d = row.vineyard;
          diagnostics << row.algorithm << ',' << row.substrateSeed << ',' << row.requestSeed << ',' << row.requestId << ',' << row.arrival << ',' << (row.accepted ? 1 : 0) << ',' << d.stage << ',' << d.candidateMin << ',' << d.candidateMax << ',' << d.unavailableCandidates << ',' << d.positiveScoreCandidates << ',' << d.selectedXSum << ',' << d.selectedFlowSum << ',' << d.selectedScoreSum << ',' << d.selectedXRankSum << ',' << d.cnlm.returnCode << ',' << d.cnlm.status << ',' << d.cnlm.iterations << ',' << d.cnlm.runtimeMs << ',' << d.cnlmPresolve.returnCode << ',' << d.cnlmPresolve.status << ',' << d.cnlmPresolve.runtimeMs << ',' << d.mcf.returnCode << ',' << d.mcf.status << ',' << d.mcf.iterations << ',' << d.mcf.runtimeMs << ',' << d.mcfPresolve.returnCode << ',' << d.mcfPresolve.status << ',' << d.mcfPresolve.runtimeMs << ',' << d.shadowXFeasible << '\n';
        }
      }
    }
  }
}

}  // namespace vne
