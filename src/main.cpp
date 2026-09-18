#include <algorithm>
#include <atomic>
#include <cctype>
#include <exception>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "core/graph.hpp"
#include "core/request.hpp"
#include "simulator/simulator.hpp"

namespace fs = std::filesystem;
namespace {

struct Options {
  std::string caseName;
  int seedCount = 0;
  int seedOffset = 0;
  int horizon = -1;
  int warmup = 0;
  int measurementWindow = -1;
  bool horizonProvided = false;
  int threads = 0;
  bool threadsProvided = false;
  std::string algorithm = "VNE-RFD-B";
  std::string sampleMode = "paired";
  fs::path resultsDirectory = "results";
  vne::AlgorithmConfig rfd;
};

int parseInt(const std::string& value, const char* name) {
  size_t end = 0;
  try { const int result = std::stoi(value, &end); if (end == value.size()) return result; } catch (...) {}
  throw std::runtime_error(std::string("invalid ") + name + ": " + value);
}
double parseDouble(const std::string& value, const char* name) {
  size_t end = 0;
  try { const double result = std::stod(value, &end); if (end == value.size()) return result; } catch (...) {}
  throw std::runtime_error(std::string("invalid ") + name + ": " + value);
}
void usage() {
  std::cout << "Usage: vne_sim --case CASE --seed N [--seed-offset N] [--sample-mode paired|cartesian] [--algo VNE-RFD-B|VNE-RFD-D|VNE-BCP|D-ViNE-LB|R-ViNE-LB] [--results-dir DIR] ...\n";
}
Options parse(int argc, char** argv) {
  Options output;
  for (int index = 1; index < argc; ++index) {
    const std::string key = argv[index];
    if (key == "--help") { usage(); std::exit(0); }
    if (key.rfind("--", 0) != 0 || index + 1 == argc) throw std::runtime_error("invalid option: " + key);
    const std::string value = argv[++index];
    if (key == "--case") output.caseName = value;
    else if (key == "--seed") output.seedCount = parseInt(value, "seed count");
    else if (key == "--seed-offset") output.seedOffset = parseInt(value, "seed offset");
    else if (key == "--algo") output.algorithm = value;
    else if (key == "--sample-mode") output.sampleMode = value;
    else if (key == "--results-dir") output.resultsDirectory = value;
    else if (key == "--tau") output.rfd.tau = parseInt(value, "tau");
    else if (key == "--ksp") output.rfd.ksp = parseInt(value, "ksp");
    else if (key == "--delta") output.rfd.delta = parseInt(value, "delta");
    else if (key == "--max-d") output.rfd.maxD = parseInt(value, "max-d");
    else if (key == "--alpha") output.rfd.alpha = parseDouble(value, "alpha");
    else if (key == "--beta") output.rfd.beta = parseDouble(value, "beta");
    else if (key == "--lambda") output.rfd.lambda = parseDouble(value, "lambda");
    else if (key == "--gamma") output.rfd.gamma = parseDouble(value, "gamma");
    else if (key == "--vine-timeout") output.rfd.vineTimeoutMs = parseInt(value, "vine timeout");
    else if (key == "--vine-seed") output.rfd.vineRandomSeed = static_cast<unsigned int>(parseInt(value, "vine seed"));
    else if (key == "--vine-diagnostics") output.rfd.vineDiagnostics = parseInt(value, "vine diagnostics") != 0;
    else if (key == "--horizon") { output.horizon = parseInt(value, "horizon"); output.horizonProvided = true; }
    else if (key == "--warmup") output.warmup = parseInt(value, "warmup");
    else if (key == "--measurement-window") output.measurementWindow = parseInt(value, "measurement window");
    else if (key == "--threads") { output.threads = parseInt(value, "threads"); output.threadsProvided = true; }
    else throw std::runtime_error("unknown option: " + key);
  }
  if (output.caseName.empty()) throw std::runtime_error("--case is required");
  if (fs::path(output.caseName).is_absolute() || output.caseName.find("..") != std::string::npos) throw std::runtime_error("case must be a directory below examples");
  if (output.resultsDirectory.empty()) throw std::runtime_error("--results-dir must not be empty");
  if (output.seedCount <= 0) throw std::runtime_error("--seed must be positive; it is the number of substrate/request seeds");
  if (output.seedOffset < 0) throw std::runtime_error("--seed-offset must not be negative");
  if (output.rfd.tau <= 0 || output.rfd.ksp <= 0 || output.rfd.delta < 0 || output.rfd.maxD <= 0 || output.rfd.alpha < 0.0 || output.rfd.beta < 0.0 || output.rfd.lambda < 0.0 || output.rfd.gamma < 0.0 || output.rfd.vineTimeoutMs <= 0 || output.horizon == 0 || output.warmup < 0 || output.measurementWindow == 0 || (output.threadsProvided && output.threads <= 0)) throw std::runtime_error("invalid numeric option");
  if (output.horizonProvided && output.measurementWindow > 0) throw std::runtime_error("--horizon and --measurement-window are mutually exclusive");
  if (output.measurementWindow > 0) output.horizon = output.warmup + output.measurementWindow;
  if (output.horizon > 0 && output.warmup >= output.horizon) throw std::runtime_error("warmup must be smaller than the simulation horizon");
  if (output.algorithm != "VNE-RFD-B" && output.algorithm != "VNE-RFD-D" && output.algorithm != "VNE-BCP" && output.algorithm != "D-ViNE-LB" && output.algorithm != "R-ViNE-LB") throw std::runtime_error("unknown --algo value");
  if (output.sampleMode != "paired" && output.sampleMode != "cartesian") throw std::runtime_error("--sample-mode must be paired or cartesian");
  if (!output.threadsProvided) {
    output.threads = static_cast<int>(std::thread::hardware_concurrency());
    if (output.threads <= 0) output.threads = 1;
  }
  return output;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse(argc, argv);
    const fs::path input = fs::path("examples") / options.caseName;
    std::vector<vne::SubstrateGraph> substrates;
    std::vector<std::vector<vne::VirtualNetworkRequest>> requestSets;
    for (int seed = 0; seed < options.seedCount; ++seed) {
      const int inputSeed = options.seedOffset + seed;
      substrates.push_back(vne::loadSubstrate(input / ("substrate_" + std::to_string(inputSeed) + ".txt")));
      requestSets.push_back(vne::loadRequests(input / ("requests_" + std::to_string(inputSeed) + ".txt")));
    }
    struct Job {
      vne::RFDVariant variant;
      bool useBcp = false;
      bool useVineyard = false;
      bool randomizedVineyard = false;
      int substrateIndex;
      int requestIndex;
      int substrateSeed;
      int requestSeed;
    };
    std::vector<Job> jobs;
    const auto addJobs = [&](vne::RFDVariant variant, bool useBcp, bool useVineyard, bool randomizedVineyard) {
      for (int substrateSeed = 0; substrateSeed < options.seedCount; ++substrateSeed) {
        const int requestBegin = options.sampleMode == "paired" ? substrateSeed : 0;
        const int requestEnd = options.sampleMode == "paired" ? substrateSeed + 1 : options.seedCount;
        for (int requestSeed = requestBegin; requestSeed < requestEnd; ++requestSeed) {
          jobs.push_back({variant, useBcp, useVineyard, randomizedVineyard,
                          substrateSeed, requestSeed,
                          options.seedOffset + substrateSeed, options.seedOffset + requestSeed});
        }
      }
    };
    if (options.algorithm == "VNE-RFD-B") addJobs(vne::RFDVariant::BreadthFirst, false, false, false);
    if (options.algorithm == "VNE-RFD-D") addJobs(vne::RFDVariant::DepthFirst, false, false, false);
    if (options.algorithm == "VNE-BCP") addJobs(vne::RFDVariant::BreadthFirst, true, false, false);
    if (options.algorithm == "D-ViNE-LB") addJobs(vne::RFDVariant::BreadthFirst, false, true, false);
    if (options.algorithm == "R-ViNE-LB") addJobs(vne::RFDVariant::BreadthFirst, false, true, true);
    std::atomic<size_t> next{0};
    std::atomic<size_t> completed{0};
    std::mutex resultMutex, errorMutex;
    std::vector<vne::SimulationResult> results;
    std::exception_ptr error;
    const int workers = std::min<int>(options.threads, jobs.size());
    std::vector<std::thread> threads;
    threads.reserve(workers);
    for (int worker = 0; worker < workers; ++worker) {
      threads.emplace_back([&] {
        while (true) {
          const size_t index = next.fetch_add(1);
          if (index >= jobs.size()) return;
          try {
            const Job job = jobs[index];
            vne::SimulatorConfig config;
            config.algorithm = options.rfd;
            config.useBcp = job.useBcp;
            config.useVineyard = job.useVineyard;
            config.randomizedVineyard = job.randomizedVineyard;
            config.horizon = options.horizon;
            config.warmup = options.warmup;
            config.substrateSeed = job.substrateSeed;
            config.requestSeed = job.requestSeed;
            vne::Simulator simulator(substrates[job.substrateIndex], requestSets[job.requestIndex], job.variant, config);
            auto result = simulator.run();
            std::lock_guard<std::mutex> lock(resultMutex);
            results.push_back(std::move(result));
          } catch (...) {
            std::lock_guard<std::mutex> lock(errorMutex);
            if (!error) error = std::current_exception();
          }
          const size_t done = completed.fetch_add(1) + 1;
          {
            std::lock_guard<std::mutex> lock(resultMutex);
            std::cerr << "\r[" << done << "/" << jobs.size() << "]" << std::flush;
          }
        }
      });
    }
    for (auto& thread : threads) thread.join();
    std::cerr << '\n';
    if (error) std::rethrow_exception(error);
    const int measurementWindow = options.measurementWindow > 0 ? options.measurementWindow :
        (options.horizon > 0 ? options.horizon - options.warmup : -1);
    vne::writeCompactCsv(options.resultsDirectory, std::move(results), options.caseName, options.seedCount, options.rfd, workers, options.warmup, measurementWindow, options.sampleMode);
    std::cerr << "[INFO] wrote timestamped result directories under " << options.resultsDirectory.string() << "\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "[ERROR] " << error.what() << '\n';
    return 1;
  }
}
