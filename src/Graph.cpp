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

    if (!coreness.empty()) {
      coreness.push_back(0);
      deg.push_back(0);
      pos.push_back(vert.size());
      vert.push_back(internal_id);
    }

    return internal_id;
  }
  return original_to_internal_id.at(original_id);
}

void Graph::add_edge(int u_orig, int v_orig) {
  if (u_orig == v_orig)
    return;

  bool is_new_edge = false;

  int u_internal = get_or_create_internal_id(u_orig);
  int v_internal = get_or_create_internal_id(v_orig);

  if (adj[u_internal].find(v_internal) == adj[u_internal].end()) {
    is_new_edge = true;
    adj[u_internal].insert(v_internal);
    adj[v_internal].insert(u_internal);
    edge_count++;
  }

  if (is_new_edge && !coreness.empty()) {
    deg[u_internal]++;
    deg[v_internal]++;

    if (coreness[u_internal] > coreness[v_internal])
      std::swap(u_internal, v_internal);

    reorder_vertices(pos[u_internal]);
    reorder_vertices(pos[v_internal]);

    int k = coreness[u_internal];
    std::vector<int> S;
    std::vector<bool> in_S(get_node_count(), false);

    std::queue<int> q;
    q.push(u_internal);
    in_S[u_internal] = true;

    while (!q.empty()) {
      int curr = q.front();
      q.pop();
      S.push_back(curr);
      for (int neighbor : get_neighbors_by_internal_id(curr)) {
        if (coreness[neighbor] == k && !in_S[neighbor]) {
          q.push(neighbor);
          in_S[neighbor] = true;
        }
      }
    }

    int min_pos = get_node_count();
    for (int node : S) {
      min_pos = std::min(min_pos, pos[node]);
    }

    for (int i = min_pos; i < get_node_count(); ++i) {
      int p = vert[i];
      if (coreness[p] > k)
        break;

      int count = 0;
      for (int neighbor : get_neighbors_by_internal_id(p)) {
        if (coreness[neighbor] > k ||
            (coreness[neighbor] == k && pos[neighbor] > i)) {
          count++;
        }
      }
      if (count > k) {
        coreness[p]++;
      }
    }
  }
}

void Graph::remove_edge(int u_orig, int v_orig) {
  if (!has_node(u_orig) || !has_node(v_orig))
    return;

  int u = original_to_internal_id.at(u_orig);
  int v = original_to_internal_id.at(v_orig);

  if (adj[u].find(v) == adj[u].end())
    return;

  adj[u].erase(v);
  adj[v].erase(u);
  edge_count--;

  if (coreness.empty())
    return;

  deg[u]--;
  deg[v]--;

  if (coreness[u] > coreness[v])
    std::swap(u, v);
  reorder_vertices(pos[u]);
  reorder_vertices(pos[v]);

  int k = coreness[u];
  std::queue<int> q;
  std::vector<bool> in_queue(get_node_count(), false);
  std::vector<int> temp_deg = deg;

  if (deg[u] < k) {
    q.push(u);
    in_queue[u] = true;
  }
  if (deg[v] < k) {
    q.push(v);
    in_queue[v] = true;
  }

  while (!q.empty()) {
    int curr = q.front();
    q.pop();

    coreness[curr]--;
    reorder_vertices(pos[curr]);

    for (int neighbor : get_neighbors_by_internal_id(curr)) {
      if (coreness[neighbor] == coreness[curr] + 1) {
        temp_deg[neighbor]--;
        if (temp_deg[neighbor] < coreness[neighbor] && !in_queue[neighbor]) {
          q.push(neighbor);
          in_queue[neighbor] = true;
        }
      }
    }
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

int Graph::get_coreness(int original_id) const {
  if (!has_node(original_id)) {
    return 0;
  }
  int internal_id = original_to_internal_id.at(original_id);
  if (internal_id >= coreness.size()) {
    return 0;
  }
  return coreness[internal_id];
}

std::vector<int> Graph::get_all_node_ids() const {
  return internal_to_original_id;
}

void Graph::load_graph_from_file(const std::string &filename) {
  std::ifstream inputFile(filename);
  if (!inputFile.is_open())
    throw std::runtime_error("Error: Could not open file " + filename);

  std::string line;
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

        current_combination.push_back(source[n]);
        find_combinations(source, n - 1, r, current_combination);
        current_combination.pop_back();

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

void Graph::run_k_core_decomposition() {
  const int num_nodes = get_node_count();
  if (num_nodes == 0)
    return;

  deg.assign(num_nodes, 0);
  coreness.assign(num_nodes, 0);
  vert.assign(num_nodes, 0);
  pos.assign(num_nodes, 0);

  int max_degree = 0;
  for (int i = 0; i < num_nodes; ++i) {
    deg[i] = get_degree_by_internal_id(i);
    if (deg[i] > max_degree)
      max_degree = deg[i];
  }

  std::vector<int> bin(max_degree + 1, 0);
  for (int i = 0; i < num_nodes; ++i)
    bin[deg[i]]++;

  int start = 0;
  for (int d = 0; d <= max_degree; ++d) {
    int num = bin[d];
    bin[d] = start;
    start += num;
  }

  for (int i = 0; i < num_nodes; ++i) {
    pos[i] = bin[deg[i]];
    vert[pos[i]] = i;
    bin[deg[i]]++;
  }

  for (int d = max_degree; d > 0; --d)
    bin[d] = bin[d - 1];
  bin[0] = 0;

  for (int i = 0; i < num_nodes; ++i) {
    int u_internal = vert[i];
    coreness[u_internal] = deg[u_internal];

    const auto &neighbors_internal = get_neighbors_by_internal_id(u_internal);
    for (int v_internal : neighbors_internal) {
      if (deg[v_internal] > coreness[u_internal]) {
        int dv = deg[v_internal];
        int pv = pos[v_internal];
        int pw = bin[dv];
        int w_internal = vert[pw];

        if (v_internal != w_internal) {
          vert[pv] = w_internal;
          vert[pw] = v_internal;
          pos[v_internal] = pw;
          pos[w_internal] = pv;
        }
        bin[dv]++;
        deg[v_internal]--;
      }
    }
  }
}

void Graph::k_core_decomposition(const std::string &filename) {
  ensure_directory_exists(filename);
  std::ofstream outputFile(filename);
  auto start_load = std::chrono::high_resolution_clock::now();

  run_k_core_decomposition();

  auto end_load = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> load_duration =
      end_load - start_load;
  outputFile << load_duration.count() << " ms\n";
  const int num_nodes = get_node_count();
  for (int i = 0; i < num_nodes; ++i)
    outputFile << get_original_id(i) << " " << coreness[i] << "\n";
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
  for (int iter = 0; iter < 100; ++iter) {
    double g = (low + high) / 2.0;
    if (g < 1e-9)
      break;
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
  std::vector<int> current_degrees(num_nodes);
  long long current_edges = get_edge_count();
  for (int i = 0; i < num_nodes; ++i) {
    current_degrees[i] = get_degree_by_internal_id(i);
  }

  std::vector<bool> removed(num_nodes, false);
  std::vector<int> degree_list = current_degrees;
  int num_remaining_nodes = num_nodes;

  double max_density = 0.0;
  std::vector<int> best_subgraph;

  for (int i = 0; i < num_nodes; ++i) {
    if (num_remaining_nodes > 0) {
      double current_density =
          static_cast<double>(current_edges) / num_remaining_nodes;
      if (current_density > max_density) {
        max_density = current_density;
        best_subgraph.clear();
        for (int j = 0; j < num_nodes; ++j) {
          if (!removed[j]) {
            best_subgraph.push_back(get_original_id(j));
          }
        }
      }
    }

    int min_degree_node = -1;
    int min_degree = std::numeric_limits<int>::max();
    for (int j = 0; j < num_nodes; ++j) {
      if (!removed[j] && degree_list[j] < min_degree) {
        min_degree = degree_list[j];
        min_degree_node = j;
      }
    }

    if (min_degree_node == -1)
      break;

    removed[min_degree_node] = true;
    num_remaining_nodes--;
    current_edges -= degree_list[min_degree_node];

    const auto &neighbors = get_neighbors_by_internal_id(min_degree_node);
    for (int neighbor : neighbors) {
      if (!removed[neighbor]) {
        degree_list[neighbor]--;
      }
    }
  }

  auto end_load = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> load_duration =
      end_load - start_load;
  outputFile << load_duration.count() << " ms\n";
  outputFile << max_density << "\n";

  for (size_t i = 0; i < best_subgraph.size(); ++i) {
    outputFile << best_subgraph[i]
               << (i == best_subgraph.size() - 1 ? "" : " ");
  }
  outputFile << "\n";
}

void Graph::reorder_vertices(int start_idx) {
  for (int i = start_idx; i < get_node_count() - 1; ++i) {
    int u = vert[i];
    int v = vert[i + 1];
    if (coreness[u] > coreness[v] ||
        (coreness[u] == coreness[v] && deg[u] > deg[v])) {
      std::swap(vert[i], vert[i + 1]);
      pos[u] = i + 1;
      pos[v] = i;
    } else {
      break;
    }
  }
  for (int i = start_idx; i > 0; --i) {
    int u = vert[i];
    int v = vert[i - 1];
    if (coreness[u] < coreness[v] ||
        (coreness[u] == coreness[v] && deg[u] < deg[v])) {
      std::swap(vert[i], vert[i - 1]);
      pos[u] = i - 1;
      pos[v] = i;
    } else {
      break;
    }
  }
}

void Graph::export_to_dot(
    const std::string &filename, const std::map<int, std::string> &node_styles,
    const std::map<std::pair<int, int>, std::string> &edge_styles) {

  ensure_directory_exists(filename);
  std::ofstream dot_file(filename);
  if (!dot_file.is_open()) {
    std::cerr << "Error: Could not open file for dot export: " << filename
              << std::endl;
    return;
  }

  dot_file << "graph G {\n";
  dot_file << "  node [shape=circle, fontsize=10, width=0.5];\n";
  dot_file << "  edge [len=1.5];\n\n";

  // Write nodes
  for (int i = 0; i < get_node_count(); ++i) {
    int original_id = get_original_id(i);
    dot_file << "  " << original_id;
    auto it = node_styles.find(original_id);
    if (it != node_styles.end()) {
      dot_file << " [" << it->second << "]";
    }
    dot_file << ";\n";
  }

  dot_file << "\n";

  // Write edges
  for (int u_internal = 0; u_internal < get_node_count(); ++u_internal) {
    int u_orig = get_original_id(u_internal);
    for (int v_internal : get_neighbors_by_internal_id(u_internal)) {
      if (u_internal < v_internal) { // To avoid duplicate edges
        int v_orig = get_original_id(v_internal);
        dot_file << "  " << u_orig << " -- " << v_orig;

        std::pair<int, int> edge = {std::min(u_orig, v_orig),
                                    std::max(u_orig, v_orig)};
        auto it = edge_styles.find(edge);
        if (it != edge_styles.end()) {
          dot_file << " [" << it->second << "]";
        }
        dot_file << ";\n";
      }
    }
  }

  dot_file << "}\n";
  dot_file.close();
}