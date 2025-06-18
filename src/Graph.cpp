#include "Graph.h"
#include "Dinic.h"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>

void ensure_directory_exists(const std::string &path) {
  std::filesystem::path p(path);
  if (p.has_parent_path()) {
    std::filesystem::create_directories(p.parent_path());
  }
}

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

double Graph::averageDegree() const {
  if (get_node_count() == 0) {
    return 0.0;
  }
  return static_cast<double>(2 * get_edge_count()) / get_node_count();
}

double Graph::density() const {
  int n = get_node_count();
  if (n <= 1) {
    return 0.0;
  }
  return static_cast<double>(2 * get_edge_count()) /
         (static_cast<long long>(n) * (n - 1));
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

  std::string line;
  // Read and discard header line
  if (!std::getline(inputFile, line)) {
    return;
  }

  while (std::getline(inputFile, line)) {
    line.erase(0, line.find_first_not_of(" \t\n\r"));
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    int u, v;
    if (iss >> u >> v) {
      add_edge(u, v);
    }
  }
}

void Graph::bron_kerbosch_pivot(std::vector<int> R, std::unordered_set<int> P,
                                std::unordered_set<int> X,
                                std::vector<std::vector<int>> &cliques) {
  if (P.empty() && X.empty()) {
    cliques.push_back(R);
    return;
  }

  if (P.empty()) {
    return;
  }

  // Choose a pivot from P U X
  int pivot = -1;
  {
    std::unordered_set<int> P_union_X = P;
    P_union_X.insert(X.begin(), X.end());
    pivot = *P_union_X.begin();
  }

  const auto &pivot_neighbors = get_neighbors_by_internal_id(pivot);
  std::unordered_set<int> P_without_pivot_neighbors;
  for (int node : P) {
    if (pivot_neighbors.find(node) == pivot_neighbors.end()) {
      P_without_pivot_neighbors.insert(node);
    }
  }

  for (int v : P_without_pivot_neighbors) {
    const auto &v_neighbors = get_neighbors_by_internal_id(v);
    std::vector<int> next_R = R;
    next_R.push_back(v);

    std::unordered_set<int> next_P;
    for (int node : P) {
      if (v_neighbors.count(node)) {
        next_P.insert(node);
      }
    }

    std::unordered_set<int> next_X;
    for (int node : X) {
      if (v_neighbors.count(node)) {
        next_X.insert(node);
      }
    }

    bron_kerbosch_pivot(next_R, next_P, next_X, cliques);

    P.erase(v);
    X.insert(v);
  }
}

void Graph::find_k_cliques(int k, const std::string &filename) {
  ensure_directory_exists(filename);
  std::ofstream outputFile(filename);
  auto start_load = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<int>> all_maximal_cliques;
  std::unordered_set<int> P, X;
  for (int i = 0; i < get_node_count(); ++i) {
    P.insert(i);
  }
  bron_kerbosch_pivot({}, P, X, all_maximal_cliques);

  std::set<std::vector<int>> unique_k_cliques;

  std::function<void(const std::vector<int> &, int, int, std::vector<int> &)>
      find_combinations = [&](const std::vector<int> &source, int n, int r,
                              std::vector<int> &current_combination) {
        if (current_combination.size() == r) {
          std::vector<int> sorted_combo = current_combination;
          std::sort(sorted_combo.begin(), sorted_combo.end());
          unique_k_cliques.insert(sorted_combo);
          return;
        }
        if (n < 0)
          return;

        // Include current element
        current_combination.push_back(source[n]);
        find_combinations(source, n - 1, r, current_combination);
        current_combination.pop_back();

        // Exclude current element
        find_combinations(source, n - 1, r, current_combination);
      };

  for (const auto &m_clique_internal : all_maximal_cliques) {
    if (m_clique_internal.size() >= k) {
      std::vector<int> m_clique_original;
      for (int node : m_clique_internal)
        m_clique_original.push_back(get_original_id(node));

      std::vector<int> current_combination;
      find_combinations(m_clique_original, m_clique_original.size() - 1, k,
                        current_combination);
    }
  }

  auto end_load = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> load_duration =
      end_load - start_load;
  outputFile << load_duration.count() << " ms\n";

  for (const auto &clique : unique_k_cliques) {
    for (size_t i = 0; i < clique.size(); ++i) {
      outputFile << clique[i] << (i == clique.size() - 1 ? "" : " ");
    }
    outputFile << "\n";
  }
  outputFile.close();
}

void Graph::find_k_clique_decomposition(int k, const std::string &output_path) {
  ensure_directory_exists(output_path);
  std::ofstream outputFile(output_path);
  auto start_time = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<int>> all_cliques_internal;
  std::unordered_set<int> P, X;
  for (int i = 0; i < get_node_count(); ++i) {
    P.insert(i);
  }
  bron_kerbosch_pivot({}, P, X, all_cliques_internal);

  std::vector<std::vector<int>> k_cliques;
  for (const auto &clique_internal : all_cliques_internal) {
    if (clique_internal.size() >= k) {
      std::vector<int> clique_original;
      for (int internal_id : clique_internal) {
        clique_original.push_back(get_original_id(internal_id));
      }
      std::sort(clique_original.begin(), clique_original.end());
      k_cliques.push_back(clique_original);
    }
  }

  std::sort(k_cliques.begin(), k_cliques.end(),
            [](const auto &a, const auto &b) { return a.size() > b.size(); });

  std::vector<std::vector<int>> maximal_cliques;
  if (!k_cliques.empty()) {
    maximal_cliques.push_back(k_cliques[0]);
    for (size_t i = 1; i < k_cliques.size(); ++i) {
      bool is_subset = false;
      for (const auto &final_clique : maximal_cliques) {
        if (std::includes(final_clique.begin(), final_clique.end(),
                          k_cliques[i].begin(), k_cliques[i].end())) {
          is_subset = true;
          break;
        }
      }
      if (!is_subset) {
        maximal_cliques.push_back(k_cliques[i]);
      }
    }
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = end_time - start_time;
  outputFile << duration.count() << " ms\n";

  for (const auto &clique : maximal_cliques) {
    for (size_t i = 0; i < clique.size(); ++i) {
      outputFile << clique[i] << (i == clique.size() - 1 ? "" : " ");
    }
    outputFile << "\n";
  }
  outputFile.close();
}

void Graph::k_core_decomposition(const std::string &filename) {
  ensure_directory_exists(filename);
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
  ensure_directory_exists(filename);
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
    for (int v_internal : get_neighbors_by_internal_id(i))
      if (i < v_internal)
        all_edges.push_back({i, v_internal});
  }
  // Binary search for 100 iterations for high precision
  for (int iter = 0; iter < 100; ++iter) {
    double g = (low + high) / 2.0;
    if (g < 1e-9)
      break; // Avoid issues with g being too small
    // Build the flow network for density 'g'
    int num_flow_nodes = num_graph_nodes + num_graph_edges + 2;
    int s = num_flow_nodes - 2;
    int t = num_flow_nodes - 1;
    Dinic<double> dinic(num_flow_nodes);
    // Add edges from graph nodes to sink T
    for (int i = 0; i < num_graph_nodes; ++i)
      dinic.add_edge(i, t, g);
    // Add edges involving graph edges
    for (int i = 0; i < num_graph_edges; ++i) {
      int edge_node_idx = num_graph_nodes + i;
      int u_internal = all_edges[i].first;
      int v_internal = all_edges[i].second;
      dinic.add_edge(s, edge_node_idx, 1.0);
      dinic.add_edge(edge_node_idx, u_internal,
                     std::numeric_limits<double>::infinity());
      dinic.add_edge(edge_node_idx, v_internal,
                     std::numeric_limits<double>::infinity());
    }
    double flow = dinic.max_flow(s, t);
    if (static_cast<double>(num_graph_edges) - flow > 1e-9) {
      low = g;
    } else {
      high = g;
    }
  }

  std::vector<int> best_subgraph_nodes;
  {
    double g = low;
    int num_flow_nodes = num_graph_nodes + num_graph_edges + 2;
    int s = num_flow_nodes - 2;
    int t = num_flow_nodes - 1;
    Dinic<double> dinic(num_flow_nodes);
    for (int i = 0; i < num_graph_nodes; ++i)
      dinic.add_edge(i, t, g);
    for (int i = 0; i < num_graph_edges; ++i) {
      int edge_node_idx = num_graph_nodes + i;
      int u_internal = all_edges[i].first;
      int v_internal = all_edges[i].second;
      dinic.add_edge(s, edge_node_idx, 1.0);
      dinic.add_edge(edge_node_idx, u_internal,
                     std::numeric_limits<double>::infinity());
      dinic.add_edge(edge_node_idx, v_internal,
                     std::numeric_limits<double>::infinity());
    }
    dinic.max_flow(s, t);

    auto s_side_all_nodes = dinic.get_min_cut_nodes(s);
    std::vector<int> v1_internal;
    for (int internal_id : s_side_all_nodes)
      if (internal_id < num_graph_nodes)
        v1_internal.push_back(internal_id);

    auto t_side_all_nodes = dinic.get_t_side_nodes(t);
    std::vector<int> v2_internal;
    std::unordered_set<int> t_side_internal_set(t_side_all_nodes.begin(),
                                                t_side_all_nodes.end());
    for (int i = 0; i < num_graph_nodes; ++i) {
      if (t_side_internal_set.find(i) == t_side_internal_set.end()) {
        v2_internal.push_back(i);
      }
    }

    auto calculate_score = [&](const std::vector<int> &nodes_internal) {
      if (nodes_internal.empty())
        return -std::numeric_limits<double>::infinity();
      long long sub_edges = 0;
      std::unordered_set<int> node_set(nodes_internal.begin(),
                                       nodes_internal.end());
      for (int u : nodes_internal) {
        for (int v : get_neighbors_by_internal_id(u)) {
          if (u < v && node_set.count(v)) {
            sub_edges++;
          }
        }
      }
      return static_cast<double>(sub_edges) - g * nodes_internal.size();
    };

    double score1 = calculate_score(v1_internal);
    double score2 = calculate_score(v2_internal);

    std::vector<int> best_subgraph_nodes_internal;
    if (score1 > score2) {
      best_subgraph_nodes_internal = v1_internal;
    } else {
      best_subgraph_nodes_internal = v2_internal;
    }

    for (int internal_id : best_subgraph_nodes_internal) {
      best_subgraph_nodes.push_back(get_original_id(internal_id));
    }
  }

  if (best_subgraph_nodes.empty() && num_graph_nodes > 0) {
    int max_degree = -1;
    int best_node = -1;
    for (int i = 0; i < num_graph_nodes; ++i) {
      if (get_degree_by_internal_id(i) > max_degree) {
        max_degree = get_degree_by_internal_id(i);
        best_node = get_original_id(i);
      }
    }
    if (best_node != -1) {
      best_subgraph_nodes.push_back(best_node);
    }
  }
  auto end_load = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> load_duration =
      end_load - start_load;
  outputFile << load_duration.count() << " ms\n";
  outputFile << low << "\n";
  for (size_t i = 0; i < best_subgraph_nodes.size(); ++i) {
    outputFile << best_subgraph_nodes[i]
               << (i == best_subgraph_nodes.size() - 1 ? "" : " ");
  }
  outputFile << "\n";
}

void Graph::find_densest_subgraph_approx(const std::string &filename) {
  ensure_directory_exists(filename);
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
      if (pos[v_internal] > pos[u_internal]) {
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
  outputFile << "\n";
}

// k-VCC Implementation
std::pair<int, std::vector<int>>
Graph::get_vertex_connectivity(int u_orig, int v_orig,
                               const std::set<int> &nodes_orig) {
  const int INF = std::numeric_limits<int>::max();
  Dinic<int> flow_network(2 * nodes_orig.size());

  std::map<int, int> node_to_flow_idx;
  int current_flow_idx = 0;
  for (int node : nodes_orig) {
    node_to_flow_idx[node] = current_flow_idx++;
  }

  int u_flow_idx = node_to_flow_idx[u_orig];
  int v_flow_idx = node_to_flow_idx[v_orig];

  for (int node : nodes_orig) {
    int flow_idx = node_to_flow_idx[node];
    int in_node = 2 * flow_idx;
    int out_node = 2 * flow_idx + 1;

    int capacity = 1;
    if (node == u_orig || node == v_orig) {
      capacity = INF;
    }
    flow_network.add_edge(in_node, out_node, capacity);
  }

  for (int node1_orig : nodes_orig) {
    if (original_to_internal_id.count(node1_orig)) {
      int node1_internal = original_to_internal_id.at(node1_orig);
      for (int node2_internal : adj[node1_internal]) {
        int node2_orig = internal_to_original_id[node2_internal];
        if (nodes_orig.count(node2_orig)) {
          int idx1_out = 2 * node_to_flow_idx[node1_orig] + 1;
          int idx2_in = 2 * node_to_flow_idx[node2_orig];
          flow_network.add_edge(idx1_out, idx2_in, INF);
        }
      }
    }
  }

  int flow = flow_network.max_flow(2 * u_flow_idx + 1, 2 * v_flow_idx);

  std::vector<int> separator;
  if (flow < INF) {
    std::vector<bool> reachable =
        flow_network.get_reachable_nodes(2 * u_flow_idx + 1);
    for (int node : nodes_orig) {
      if (node == u_orig || node == v_orig)
        continue;
      int flow_idx = node_to_flow_idx[node];
      int in_node = 2 * flow_idx;
      int out_node = 2 * flow_idx + 1;
      if (reachable[in_node] && !reachable[out_node]) {
        separator.push_back(node);
      }
    }
  }

  return {flow, separator};
}

void Graph::find_k_vcc_recursive(const std::set<int> &current_nodes, int k,
                                 std::vector<std::vector<int>> &results) {
  if (current_nodes.size() < k) {
    return;
  }

  std::vector<int> node_vec(current_nodes.begin(), current_nodes.end());

  for (size_t i = 0; i < node_vec.size(); ++i) {
    for (size_t j = i + 1; j < node_vec.size(); ++j) {
      int u = node_vec[i];
      int v = node_vec[j];

      int u_internal = original_to_internal_id.at(u);
      int v_internal = original_to_internal_id.at(v);
      bool is_adjacent = adj[u_internal].count(v_internal);

      if (is_adjacent) {
        continue;
      }

      auto [conn, separator] = get_vertex_connectivity(u, v, current_nodes);

      if (conn < k) {
        std::set<int> nodes_to_partition = current_nodes;
        for (int sep_node : separator) {
          nodes_to_partition.erase(sep_node);
        }

        std::set<int> remaining_nodes = nodes_to_partition;
        while (!remaining_nodes.empty()) {
          std::set<int> component;
          std::queue<int> q;
          int start_node_orig = *remaining_nodes.begin();

          q.push(start_node_orig);
          component.insert(start_node_orig);
          remaining_nodes.erase(start_node_orig);

          while (!q.empty()) {
            int u_bfs_orig = q.front();
            q.pop();
            int u_bfs_internal = original_to_internal_id.at(u_bfs_orig);

            for (int v_bfs_internal : adj[u_bfs_internal]) {
              int v_bfs_orig = internal_to_original_id[v_bfs_internal];
              if (remaining_nodes.count(v_bfs_orig)) {
                component.insert(v_bfs_orig);
                remaining_nodes.erase(v_bfs_orig);
                q.push(v_bfs_orig);
              }
            }
          }

          std::set<int> next_subgraph = component;
          next_subgraph.insert(separator.begin(), separator.end());
          find_k_vcc_recursive(next_subgraph, k, results);
        }
        return;
      }
    }
  }

  results.emplace_back(node_vec);
}

void Graph::find_k_vcc(int k, const std::string &output_path) {
  auto start = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<int>> results;
  std::set<int> all_nodes;
  for (const auto &orig_id : internal_to_original_id) {
    all_nodes.insert(orig_id);
  }

  if (all_nodes.size() >= k) {
    find_k_vcc_recursive(all_nodes, k, results);
  }

  for (auto &res : results) {
    std::sort(res.begin(), res.end());
  }

  std::sort(results.begin(), results.end(),
            [](const auto &a, const auto &b) { return a.size() > b.size(); });

  std::vector<std::vector<int>> final_results;
  if (!results.empty()) {
    final_results.push_back(results[0]);
    for (size_t i = 1; i < results.size(); ++i) {
      bool is_subset = false;
      for (const auto &final_res : final_results) {
        if (std::includes(final_res.begin(), final_res.end(),
                          results[i].begin(), results[i].end())) {
          is_subset = true;
          break;
        }
      }
      if (!is_subset) {
        final_results.push_back(results[i]);
      }
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "k-VCC decomposition found " << final_results.size()
            << " components in " << diff.count() << "s." << std::endl;

  ensure_directory_exists(output_path);
  std::ofstream outfile(output_path);
  if (!outfile.is_open()) {
    std::cerr << "Error opening k-vcc output file: " << output_path
              << std::endl;
    return;
  }

  outfile << diff.count() << "s" << std::endl;
  for (const auto &vcc : final_results) {
    for (size_t i = 0; i < vcc.size(); ++i) {
      outfile << vcc[i] << (i == vcc.size() - 1 ? "" : " ");
    }
    outfile << std::endl;
  }
}