#include "Graph.h"
#include <iostream>

int main() {
  Graph graph("../data/Amazon.txt");
  graph.k_core_decomposition("../results/k_core_output.txt");
  return 0;
}
