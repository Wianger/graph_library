#include "Graph.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>

int main(int argc, char *argv[]) {
  std::cout << "Loading complex graph from ./data/complex_test.txt"
            << std::endl;
  Graph g("./data/complex_test.txt");

  std::cout << "\n--- Running k-core decomposition ---" << std::endl;
  g.k_core_decomposition("./results/complex_k_core_output.txt");

  std::cout << "\n--- Running densest subgraph (exact) ---" << std::endl;
  g.find_densest_subgraph_exact(
      "./results/complex_densest_subgraph_exact_output.txt");

  std::cout << "\n--- Running densest subgraph (approx) ---" << std::endl;
  g.find_densest_subgraph_approx(
      "./results/complex_densest_subgraph_approx_output.txt");

  int k_clique_val = 4;
  std::cout << "\n--- Running k-clique finding for k=" << k_clique_val << " ---"
            << std::endl;
  g.find_k_cliques(k_clique_val, "./results/complex_k_clique_output.txt");

  std::cout << "\n--- Running k-clique decomposition for k=" << k_clique_val
            << " ---" << std::endl;
  g.find_k_clique_decomposition(
      k_clique_val, "./results/complex_k_clique_decomposition_output.txt");

  std::cout << "\nAll final tests completed." << std::endl;

  return 0;
}
