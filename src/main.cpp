#include "Graph.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>

int main(int argc, char *argv[]) {
  std::cout << "Loading hyper-complex graph from ./data/hyper_complex_test.txt"
            << std::endl;
  Graph g("./data/hyper_complex_test.txt");

  std::cout << "\n--- Running k-core decomposition ---" << std::endl;
  g.k_core_decomposition("./results/hyper_k_core_output.txt");

  std::cout << "\n--- Running densest subgraph (exact) ---" << std::endl;
  g.find_densest_subgraph_exact(
      "./results/hyper_densest_subgraph_exact_output.txt");

  std::cout << "\n--- Running densest subgraph (approx) ---" << std::endl;
  g.find_densest_subgraph_approx(
      "./results/hyper_densest_subgraph_approx_output.txt");

  int k_clique_val = 5;
  std::cout << "\n--- Running k-clique finding for k=" << k_clique_val << " ---"
            << std::endl;
  g.find_k_cliques(k_clique_val, "./results/hyper_k_clique_output.txt");

  std::cout << "\n--- Running k-clique decomposition for k=" << k_clique_val
            << " ---" << std::endl;
  g.find_k_clique_decomposition(
      k_clique_val, "./results/hyper_k_clique_decomposition_output.txt");

  std::cout << "\n--- Running k-VCC decomposition for k=2 ---" << std::endl;
  g.find_k_vcc(2, "./results/hyper_k_vcc_2_output.txt");

  std::cout << "\n--- Running k-VCC decomposition for k=3 ---" << std::endl;
  g.find_k_vcc(3, "./results/hyper_k_vcc_3_output.txt");

  std::cout << "\n--- Running k-VCC decomposition for k=4 ---" << std::endl;
  g.find_k_vcc(4, "./results/hyper_k_vcc_4_output.txt");

  std::cout << "\nAll tests on hyper-complex graph completed." << std::endl;

  return 0;
}
