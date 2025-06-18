#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Graph {
public:
  Graph(const std::string &filename);
  double averageDegree() const;
  double density() const;
  void add_edge(int u_orig, int v_orig);
  int get_node_count() const;
  long long get_edge_count() const;
  bool has_node(int original_id) const;
  int get_degree(int original_id) const;
  void k_core_decomposition(const std::string &filename);
  void find_densest_subgraph_exact(const std::string &filename);
  void find_densest_subgraph_approx(const std::string &filename);
  void find_k_cliques(int k, const std::string &filename);

private:
  std::vector<std::unordered_set<int>> adj;
  std::unordered_map<int, int> original_to_internal_id;
  std::vector<int> internal_to_original_id;
  long long edge_count = 0;
  int get_or_create_internal_id(int original_id);
  int get_degree_by_internal_id(int internal_id) const;
  const std::unordered_set<int> &
  get_neighbors_by_internal_id(int internal_id) const;
  int get_original_id(int internal_id) const;
  void load_graph_from_file(const std::string &filename);
  void bron_kerbosch_pivot(std::vector<int> R, std::unordered_set<int> P,
                           std::unordered_set<int> X,
                           std::vector<std::vector<int>> &cliques);
};
