cmake_minimum_required(VERSION 3.16)
project(vne_framework LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_path(GLPK_INCLUDE_DIR glpk.h REQUIRED)
find_library(GLPK_LIBRARY NAMES glpk REQUIRED)

add_executable(vne_sim
  src/main.cpp
  src/core/graph.cpp
  src/core/request.cpp
  src/algorithms/rfd_algorithm.cpp
  src/algorithms/bcp_algorithm.cpp
  src/algorithms/vineyard_algorithm.cpp
  src/simulator/simulator.cpp
)
target_include_directories(vne_sim PRIVATE include)
target_include_directories(vne_sim PRIVATE ${GLPK_INCLUDE_DIR})
target_link_libraries(vne_sim PRIVATE ${GLPK_LIBRARY})
target_compile_options(vne_sim PRIVATE -O3 -Wall -Wextra -Wpedantic)
