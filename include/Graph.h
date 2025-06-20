#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

class Graph {
public:
  Graph(const std::string &filename);
  double averageDegree() const;
  double density() const;
  void add_edge(int u_orig, int v_orig);
  void remove_edge(int u_orig, int v_orig);
  int get_node_count() const;
  long long get_edge_count() const;
  bool has_node(int original_id) const;
  int get_degree(int original_id) const;
  int get_coreness(int original_id) const;
  void run_k_core_decomposition();
  std::vector<int> get_all_node_ids() const;
  void k_core_decomposition(const std::string &filename);
  void find_densest_subgraph_exact(const std::string &filename);
  void find_densest_subgraph_approx(const std::string &filename);
  void find_k_cliques(int k, const std::string &filename);
  void find_k_clique_decomposition(int k, const std::string &output_path);
  void export_to_dot(
      const std::string &filename,
      const std::map<int, std::string> &node_styles = {},
      const std::map<std::pair<int, int>, std::string> &edge_styles = {});

private:
  struct Separation {
    std::vector<int> separator;
    std::vector<std::unordered_set<int>> components;
  };

  std::vector<std::unordered_set<int>> adj;
  std::unordered_map<int, int> original_to_internal_id;
  std::vector<int> internal_to_original_id;
  long long edge_count = 0;

  std::vector<int> coreness;
  std::vector<int> vert;
  std::vector<int> pos;
  std::vector<int> deg;

  int get_or_create_internal_id(int original_id);
  int get_degree_by_internal_id(int internal_id) const;
  const std::unordered_set<int> &
  get_neighbors_by_internal_id(int internal_id) const;
  int get_original_id(int internal_id) const;
  void load_graph_from_file(const std::string &filename);
  void bron_kerbosch_pivot(std::vector<int> R, std::unordered_set<int> P,
                           std::unordered_set<int> X,
                           std::vector<std::vector<int>> &cliques);
  void reorder_vertices(int start_idx);
};
