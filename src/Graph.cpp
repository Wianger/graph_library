#include "Graph.h"
#include <fstream>

Graph::Graph(const std::string &filename) { loadGraphFromFile(filename); }

double Graph::averageDegree() const {
  if (adj.empty())
    return 0.0;
  return static_cast<double>(edges_count * 2) / adj.size();
}

double Graph::density() const {
  if (adj.empty())
    return 0.0;
  int maxEdges = adj.size() * (adj.size() - 1) / 2;
  return static_cast<double>(edges_count) / maxEdges;
}

void Graph::loadGraphFromFile(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open())
    throw std::runtime_error("Could not open file: " + filename);
  int numNodes, numEdges;
  file >> numNodes >> numEdges;
  int u, v;
  while (file >> u >> v) {
    if (u == v)
      continue; // Skip self-loops if not desired
    if (adj[u].insert(v).second) {
      adj[v].insert(u); // Ensure undirected graph
      edges_count++;
    }
  }
  file.close();
}