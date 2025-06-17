#include "Graph.h"
#include <iostream>

int main() {
  Graph graph("../data/Amazon.txt");
  graph.k_core_decomposition("../results/k_core_output.txt");
  graph.find_densest_subgraph_exact(
      "../results/densest_subgraph_exact_output.txt");
  graph.find_densest_subgraph_approx(
      "../results/densest_subgraph_approx_output.txt");
  return 0;
}
