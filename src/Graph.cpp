#include "Graph.h"
#include "Dinic.h"
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

void Graph::find_densest_subgraph_exact(const std::string &filename) {
  int num_graph_nodes = get_node_count();
  if (num_graph_nodes == 0)
    return;
  std::ofstream outputFile(filename);
  auto start_load = std::chrono::high_resolution_clock::now();
  long long num_graph_edges = get_edge_count();
  double low = 0.0, high = num_graph_edges;
  std::vector<std::pair<int, int>> all_edges;
  all_edges.reserve(num_graph_edges);
  for (int i = 0; i < num_graph_nodes; ++i) {
    int u_orig = get_original_id(i);
    for (int v_internal : get_neighbors_by_internal_id(i))
      if (i < v_internal)
        all_edges.push_back({i, v_internal});
  }
  std::vector<int> best_subgraph_nodes;
  // Binary search for 100 iterations for high precision
  for (int iter = 0; iter < 100; ++iter) {
    double g = (low + high) / 2.0;
    if (g < 1e-9)
      break; // Avoid issues with g being too small
    // Build the flow network for density 'g'
    int num_flow_nodes = num_graph_nodes + num_graph_edges + 2;
    int s = num_flow_nodes - 2;
    int t = num_flow_nodes - 1;
    Dinic dinic(num_flow_nodes);
    // Map graph nodes and edges to flow network nodes
    // Graph nodes: 0 to num_graph_nodes - 1
    // Edge nodes: num_graph_nodes to num_graph_nodes + num_graph_edges - 1
    // Add edges from graph nodes to sink T
    for (int i = 0; i < num_graph_nodes; ++i)
      dinic.add_edge(i, t, g);
    // Add edges involving graph edges
    for (int i = 0; i < num_graph_edges; ++i) {
      int edge_node_idx = num_graph_nodes + i;
      int u_internal = all_edges[i].first;
      int v_internal = all_edges[i].second;
      dinic.add_edge(s, edge_node_idx, 1.0);
      dinic.add_edge(edge_node_idx, u_internal, 1e18); // "Infinity"
      dinic.add_edge(edge_node_idx, v_internal, 1e18); // "Infinity"
    }
    long long flow = dinic.max_flow(s, t);
    // This condition is based on the formula C(g) = |E| - max_flow(g)
    // If C(g) > 0, a denser subgraph exists.
    if (static_cast<double>(num_graph_edges) - flow > 1e-9) {
      low = g;
      // Recover the subgraph from the min-cut
      auto s_side_nodes_internal = dinic.get_min_cut_nodes(s);
      best_subgraph_nodes.clear();
      for (int internal_id : s_side_nodes_internal)
        // We only care about the original graph nodes
        if (internal_id < num_graph_nodes)
          best_subgraph_nodes.push_back(get_original_id(internal_id));
    } else
      high = g;
  }
  // Important: if the graph is empty, the returned node list can be non-empty
  // but the density is 0. We need to handle this.
  if (best_subgraph_nodes.empty())
    // Find a single node with the highest degree if nothing else is found.
    // This is a fallback for disconnected graphs or trivial cases.
    if (num_graph_nodes > 0)
      best_subgraph_nodes.push_back(get_original_id(0));
  auto end_load = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> load_duration =
      end_load - start_load;
  outputFile << load_duration.count() << " ms\n";
  outputFile << low << "\n";
  for (int node : best_subgraph_nodes)
    outputFile << node << " ";
}

void Graph::find_densest_subgraph_approx(const std::string &filename) {
  const int num_nodes = get_node_count();
  if (num_nodes == 0)
    return;
  std::ofstream outputFile(filename);
  auto start_load = std::chrono::high_resolution_clock::now();
  std::vector<int> degrees(num_nodes);
  long long current_edges = get_edge_count();
  int max_degree = 0;
  for (int i = 0; i < num_nodes; ++i) {
    degrees[i] = get_degree_by_internal_id(i);
    if (degrees[i] > max_degree)
      max_degree = degrees[i];
  }
  // --- Step 2: Bucket Sort (same as k-core) ---
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

  // --- Greedy Peeling ---
  double max_density = 0.0;
  int best_subgraph_size = 0;
  // Iterate through nodes in increasing order of degree, "removing" them
  for (int i = 0; i < num_nodes; ++i) {
    // Calculate density of the *current* graph (before removing the node)
    int current_nodes = num_nodes - i;
    if (current_nodes > 0) {
      double current_density =
          static_cast<double>(current_edges) / current_nodes;
      if (current_density > max_density) {
        max_density = current_density;
        best_subgraph_size = current_nodes;
      }
    }
    // "Remove" the node with the lowest current degree
    int u_internal = vert[i];
    // Update total edge count
    current_edges -= degrees[u_internal];
    // Update its neighbors' degrees (same logic as k-core)
    const auto &neighbors_internal = get_neighbors_by_internal_id(u_internal);
    for (int v_internal : neighbors_internal) {
      // No need to check degree ordering here, just update everyone
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
  auto end_load = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> load_duration =
      end_load - start_load;
  // --- Reconstruct the Subgraph ---
  std::vector<int> densest_subgraph_nodes;
  // The best subgraph consists of the `best_subgraph_size` nodes
  // that were *last* to be removed. These are at the end of the `vert` array.
  densest_subgraph_nodes.reserve(best_subgraph_size);
  for (int i = num_nodes - best_subgraph_size; i < num_nodes; ++i)
    densest_subgraph_nodes.push_back(get_original_id(vert[i]));
  outputFile << load_duration.count() << " ms\n";
  outputFile << max_density << "\n";
  for (int node : densest_subgraph_nodes)
    outputFile << node << " ";
}