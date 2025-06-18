#include "Graph.h"
#include <iostream>

int main() {
  Graph graph("../data/test.txt");
  graph.k_core_decomposition("../results/k_core_output.txt");
  graph.find_densest_subgraph_exact("../results/densest_subgraph_output.txt");
  graph.find_densest_subgraph_approx(
      "../results/densest_subgraph_approx_output.txt");
  graph.find_k_cliques(4, "../results/k_clique_output.txt");
  return 0;
}
