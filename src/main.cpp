#include "Graph.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <vector>

// Helper to get a color based on a value (e.g., coreness)
std::string get_color_for_value(int value, int max_value) {
  if (max_value <= 0)
    return "\"#f0f0f0\""; // light grey for 0-core
  // Simple gradient from yellow to red
  double ratio = static_cast<double>(value) / max_value;
  int r = 255;
  int g = static_cast<int>(255 * (1 - ratio));
  int b = 0;
  std::stringstream ss;
  ss << "\"#" << std::hex << std::setw(2) << std::setfill('0') << r
     << std::setw(2) << std::setfill('0') << g << std::setw(2)
     << std::setfill('0') << b << "\"";
  return ss.str();
}

void process_graph(const std::string &filename) {
  std::string prefix = filename.substr(0, filename.find('.'));
  std::cout << "\n\n=============================================\n";
  std::cout << "===  PROCESSING: " << filename << "  ===\n";
  std::cout << "=============================================\n";

  try {
    Graph g("../data/" + filename);
    std::cout << "Graph loaded. Nodes: " << g.get_node_count()
              << ", Edges: " << g.get_edge_count() << std::endl;

    // --- 1. K-Core Decomposition & Visualization ---
    std::cout << "\n[1] Running k-core decomposition...\n";
    g.run_k_core_decomposition();
    std::map<int, std::string> kcore_styles;
    int max_core = 0;
    const auto all_nodes = g.get_all_node_ids();
    for (int original_id : all_nodes) {
      max_core = std::max(max_core, g.get_coreness(original_id));
    }
    for (int original_id : all_nodes) {
      int coreness = g.get_coreness(original_id);
      if (coreness > max_core / 2) { // Style nodes in the upper half of cores
        kcore_styles[original_id] = "style=filled, fillcolor=" +
                                    get_color_for_value(coreness, max_core) +
                                    ", label=\"" + std::to_string(original_id) +
                                    " (k=" + std::to_string(coreness) + ")\"";
      }
    }
    g.export_to_dot("../results/" + prefix + "_kcore.dot", kcore_styles);
    std::cout << "  -> Generated " << prefix
              << "_kcore.dot (nodes with coreness > " << max_core / 2
              << " are highlighted)\n";

    // --- 2. Densest Subgraph (Approximate) & Visualization ---
    std::cout << "\n[2] Finding approximate densest subgraph...\n";
    g.find_densest_subgraph_approx("../results/" + prefix +
                                   "_densest_approx.txt");
    std::ifstream densest_file("../results/" + prefix + "_densest_approx.txt");
    if (densest_file.is_open()) {
      std::string line;
      std::getline(densest_file, line); // Skip time
      std::getline(densest_file, line); // Skip density
      std::getline(densest_file, line); // Get nodes
      std::stringstream ss(line);
      int node_id;
      std::map<int, std::string> densest_styles;
      while (ss >> node_id) {
        densest_styles[node_id] = "style=filled, fillcolor=red";
      }
      g.export_to_dot("../results/" + prefix + "_densest_approx.dot",
                      densest_styles);
      std::cout << "  -> Generated " << prefix
                << "_densest_approx.dot with densest subgraph highlighted\n";
    }

    // --- 3. K-Clique Decomposition (only on CondMat) ---
    if (prefix == "CondMat") {
      int k = 10;
      std::cout << "\n[3] Finding " << k << "-clique decomposition...\n";
      g.find_k_clique_decomposition(k, "../results/" + prefix +
                                           "_clique_decomp.txt");
      std::ifstream clique_file("../results/" + prefix + "_clique_decomp.txt");
      if (clique_file.is_open()) {
        std::string line;
        std::getline(clique_file, line); // skip time
        std::map<int, std::string> clique_styles;
        while (std::getline(clique_file, line)) {
          std::stringstream ss(line);
          int node_id;
          while (ss >> node_id) {
            clique_styles[node_id] = "style=filled, fillcolor=blue";
          }
        }
        g.export_to_dot("../results/" + prefix + "_clique_decomp.dot",
                        clique_styles);
        std::cout << "  -> Generated " << prefix << "_clique_decomp.dot with "
                  << k << "-cliques highlighted\n";
      }
    }

  } catch (const std::exception &e) {
    std::cerr << "Error processing " << filename << ": " << e.what()
              << std::endl;
  }
}

int main() {
  try {
    std::vector<std::string> datasets = {"CondMat.txt", "Amazon.txt",
                                         "Gowalla.txt"};

    for (const auto &ds_file : datasets) {
      process_graph(ds_file);
    }

  } catch (const std::exception &e) {
    std::cerr << "An error occurred: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
