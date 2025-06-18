cmake_minimum_required(VERSION 3.10)
project(GraphLibrary)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)

add_executable(graph_library src/main.cpp src/Graph.cpp src/Dinic.cpp) 