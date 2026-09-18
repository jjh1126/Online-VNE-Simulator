#pragma once

#include <filesystem>
#include <vector>

#include "core/graph.hpp"

namespace vne {

struct VirtualNetworkRequest {
  int id = -1;
  int arrival = 0;
  int duration = 0;
  int departure = 0;
  std::vector<VirtualNode> nodes;
  std::vector<VirtualEdge> edges;
};

std::vector<VirtualNetworkRequest> loadRequests(const std::filesystem::path& path);

}  // namespace vne
