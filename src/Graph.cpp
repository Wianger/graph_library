#include "Graph.h"
#include <chrono>
#include <fstream>

Graph::Graph(const std::string &filename) { load_graph_from_file(filename); }

int Graph::get_or_create_internal_id(int original_id) {
  if (original_to_internal_id.find(original_id) ==
      original_to_internal_id.end()) {
    int internal_id = internal_to_original_id.size();
    original_to_internal_id[original_id] = internal_id;
    internal_to_original_id.push_back(original_id);
    adj.emplace_back();
    return internal_id;
  }
  return original_to_internal_id.at(original_id);
}

void Graph::add_edge(int u_orig, int v_orig) {
  if (u_orig == v_orig)
    return;
  int u_internal = get_or_create_internal_id(u_orig);
  int v_internal = get_or_create_internal_id(v_orig);
  if (adj[u_internal].insert(v_internal).second) {
    adj[v_internal].insert(u_internal);
    edge_count++;
  }
}

int Graph::get_node_count() const { return internal_to_original_id.size(); }

long long Graph::get_edge_count() const { return edge_count; }

bool Graph::has_node(int original_id) const {
  return original_to_internal_id.count(original_id);
}

int Graph::get_degree(int original_id) const {
  if (!has_node(original_id))
    return 0;
  return adj[original_to_internal_id.at(original_id)].size();
}

int Graph::get_degree_by_internal_id(int internal_id) const {
  return adj[internal_id].size();
}
const std::unordered_set<int> &
Graph::get_neighbors_by_internal_id(int internal_id) const {
  return adj[internal_id];
}
int Graph::get_original_id(int internal_id) const {
  return internal_to_original_id[internal_id];
}

void Graph::load_graph_from_file(const std::string &filename) {
  std::ifstream inputFile(filename);
  if (!inputFile.is_open())
    throw std::runtime_error("Error: Could not open file " + filename);
  long long numNodes_header, numEdges_header;
  inputFile >> numNodes_header >> numEdges_header;
  int u, v;
  while (inputFile >> u >> v)
    add_edge(u, v);
  inputFile.close();
}

void Graph::k_core_decomposition(const std::string &filename) {
  const int num_nodes = get_node_count();
  if (num_nodes == 0)
    return;
  std::ofstream outputFile(filename);
  auto start_load = std::chrono::high_resolution_clock::now();
  std::vector<int> degrees(num_nodes);
  int max_degree = 0;
  for (int i = 0; i < num_nodes; ++i) {
    degrees[i] = get_degree_by_internal_id(i);
    if (degrees[i] > max_degree)
      max_degree = degrees[i];
  }
  // --- 桶排序 ---
  std::vector<int> vert(num_nodes);
  std::vector<int> pos(num_nodes);
  std::vector<int> bin(max_degree + 1, 0);
  for (int i = 0; i < num_nodes; ++i)
    bin[degrees[i]]++;
  int start = 0;
  for (int d = 0; d <= max_degree; ++d) {
    int num = bin[d];
    bin[d] = start;
    start += num;
  }
  for (int i = 0; i < num_nodes; ++i) {
    pos[i] = bin[degrees[i]];
    vert[pos[i]] = i;
    bin[degrees[i]]++;
  }
  for (int d = max_degree; d > 0; --d)
    bin[d] = bin[d - 1];
  bin[0] = 0;
  // --- 迭代剥离 ---
  for (int i = 0; i < num_nodes; ++i) {
    int u_internal = vert[i];
    const auto &neighbors_internal = get_neighbors_by_internal_id(u_internal);
    for (int v_internal : neighbors_internal) {
      if (degrees[v_internal] > degrees[u_internal]) {
        int v_pos = pos[v_internal];
        int v_degree = degrees[v_internal];
        int start_of_v_bucket = bin[v_degree];
        int w_internal = vert[start_of_v_bucket];
        if (v_internal != w_internal) {
          vert[v_pos] = w_internal;
          vert[start_of_v_bucket] = v_internal;
          pos[v_internal] = start_of_v_bucket;
          pos[w_internal] = v_pos;
        }
        bin[v_degree]++;
        degrees[v_internal]--;
      }
    }
  }
  auto end_load = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> load_duration =
      end_load - start_load;
  outputFile << load_duration.count() << " ms\n";
  for (int i = 0; i < num_nodes; ++i)
    outputFile << get_original_id(i) << " " << degrees[i] << "\n";
  outputFile.close();
}
