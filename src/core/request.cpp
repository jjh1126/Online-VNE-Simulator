#include "core/request.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace vne {
namespace {
std::vector<std::string> lines(const std::filesystem::path& path) {
  std::ifstream input(path); if (!input) throw std::runtime_error("cannot open " + path.string());
  std::vector<std::string> output; std::string line;
  while (std::getline(input, line)) { const auto comment = line.find('#'); if (comment != std::string::npos) line.erase(comment); const auto begin = line.find_first_not_of(" \t\r\n"); if (begin != std::string::npos) output.push_back(line.substr(begin, line.find_last_not_of(" \t\r\n") - begin + 1)); }
  return output;
}
std::vector<std::string> split(const std::string& line) { std::istringstream input(line); std::vector<std::string> output; std::string value; while (input >> value) output.push_back(value); return output; }
int toInt(const std::string& value, const char* name) { size_t end = 0; try { const int result = std::stoi(value, &end); if (end == value.size()) return result; } catch (...) {} throw std::runtime_error(std::string("invalid ") + name + ": " + value); }
double toDouble(const std::string& value, const char* name) { size_t end = 0; try { const double result = std::stod(value, &end); if (end == value.size() && std::isfinite(result)) return result; } catch (...) {} throw std::runtime_error(std::string("invalid ") + name + ": " + value); }
}  // namespace

std::vector<VirtualNetworkRequest> loadRequests(const std::filesystem::path& path) {
  const auto input = lines(path); if (input.empty()) throw std::runtime_error("empty requests file");
  const auto header = split(input.front()); if (header.size() != 2 || header[0] != "REQUESTS") throw std::runtime_error("invalid requests header");
  const int count = toInt(header[1], "request count"); if (count < 0) throw std::runtime_error("negative request count");
  std::vector<VirtualNetworkRequest> requests; int position = 1;
  while (position < static_cast<int>(input.size())) {
    const auto requestHeader = split(input[position++]); if (requestHeader.size() != 4 || requestHeader[0] != "REQ") throw std::runtime_error("invalid REQ record");
    VirtualNetworkRequest request;
    request.id = toInt(requestHeader[1], "request id");
    request.arrival = toInt(requestHeader[2], "arrival");
    request.duration = toInt(requestHeader[3], "duration");
    if (request.arrival < 0 || request.duration < 0) throw std::runtime_error("negative arrival or duration");
    request.departure = request.arrival + request.duration;
    if (position >= static_cast<int>(input.size())) throw std::runtime_error("missing NV record");
    auto item = split(input[position++]);
    if (item.size() != 2 || item[0] != "NV") throw std::runtime_error("invalid NV record");
    const int nodes = toInt(item[1], "virtual node count");
    if (nodes <= 0) throw std::runtime_error("virtual node count must be positive");
    if (position >= static_cast<int>(input.size())) throw std::runtime_error("missing VE record");
    item = split(input[position++]);
    if (item.size() != 2 || item[0] != "VE") throw std::runtime_error("invalid VE record");
    const int edges = toInt(item[1], "virtual edge count");
    if (edges < 0) throw std::runtime_error("negative virtual edge count");
    request.nodes.resize(nodes);
    for (int node = 0; node < nodes; ++node) { if (position >= static_cast<int>(input.size())) throw std::runtime_error("missing V record"); item = split(input[position++]); if ((item.size() != 3 && item.size() != 6) || item[0] != "V" || toInt(item[1], "virtual node id") != node) throw std::runtime_error("invalid V record"); request.nodes[node].cpu = toDouble(item[2], "virtual CPU"); if (request.nodes[node].cpu < 0.0) throw std::runtime_error("negative virtual CPU"); if (item.size() == 6) { request.nodes[node].hasLocation = true; request.nodes[node].x = toDouble(item[3], "virtual node x"); request.nodes[node].y = toDouble(item[4], "virtual node y"); request.nodes[node].radius = toDouble(item[5], "virtual node radius"); if (request.nodes[node].radius < 0.0) throw std::runtime_error("negative virtual node radius"); } }
    for (int edge = 0; edge < edges; ++edge) { if (position >= static_cast<int>(input.size())) throw std::runtime_error("missing E record"); item = split(input[position++]); if (item.size() != 4 || item[0] != "E") throw std::runtime_error("invalid E record"); VirtualEdge virtualEdge{toInt(item[1], "virtual edge endpoint"), toInt(item[2], "virtual edge endpoint"), toDouble(item[3], "virtual bandwidth")}; if (virtualEdge.a < 0 || virtualEdge.b < 0 || virtualEdge.a >= nodes || virtualEdge.b >= nodes || virtualEdge.a == virtualEdge.b || virtualEdge.bandwidth < 0.0) throw std::runtime_error("invalid virtual edge"); request.edges.push_back(virtualEdge); }
    if (position >= static_cast<int>(input.size()) || split(input[position++]) != std::vector<std::string>{"END"}) throw std::runtime_error("missing END record");
    requests.push_back(std::move(request));
  }
  if (static_cast<int>(requests.size()) != count) throw std::runtime_error("request count mismatch");
  std::sort(requests.begin(), requests.end(), [](const auto& a, const auto& b) { return a.arrival != b.arrival ? a.arrival < b.arrival : a.id < b.id; });
  return requests;
}
}  // namespace vne
