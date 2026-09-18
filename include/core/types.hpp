#pragma once

#include <string>

namespace vne {

using NodeId = int;
using EdgeId = int;
inline constexpr double kEps = 1e-9;

enum class RFDVariant { BreadthFirst, DepthFirst };
enum class RejectReason { None, Cpu, Route, Backtrack };

inline const char* algorithmName(RFDVariant value) {
  return value == RFDVariant::BreadthFirst ? "VNE-RFD-B" : "VNE-RFD-D";
}
inline const char* rejectReasonName(RejectReason value) {
  switch (value) {
    case RejectReason::None: return "accepted";
    case RejectReason::Cpu: return "cpu";
    case RejectReason::Route: return "route";
    case RejectReason::Backtrack: return "backtrack";
  }
  return "route";
}

struct AlgorithmConfig {
  int tau = 3;
  int ksp = 10;
  int delta = 100;
  int maxD = 1;
  double alpha = 1.0;
  double beta = 1.0;
  double lambda = 1.0;
  double gamma = 1.0;
  int vineTimeoutMs = 10000;
  unsigned int vineRandomSeed = 1;
  bool vineDiagnostics = false;
};

struct VineSolverDiagnostic {
  int returnCode = -1;
  int status = 0;
  int iterations = 0;
  double runtimeMs = 0.0;
};

struct VineyardDiagnostic {
  bool enabled = false;
  std::string stage = "disabled";
  int candidateMin = 0;
  int candidateMax = 0;
  int unavailableCandidates = 0;
  int positiveScoreCandidates = 0;
  double selectedXSum = 0.0;
  double selectedFlowSum = 0.0;
  double selectedScoreSum = 0.0;
  int selectedXRankSum = 0;
  VineSolverDiagnostic cnlm;
  VineSolverDiagnostic cnlmPresolve;
  VineSolverDiagnostic mcf;
  VineSolverDiagnostic mcfPresolve;
  int shadowXFeasible = -1;
};

struct RFDMetric {
  double node = 0.0;
  double link = 0.0;
  double total = 0.0;
  double fragmentationCost = 0.0;
};

}  // namespace vne
