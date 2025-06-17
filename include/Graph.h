#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Graph {
public:
  Graph(const std::string &filename);
  double averageDegree() const;
  double density() const;

private:
  int edges_count = 0;
  std::unordered_map<int, std::unordered_set<int>> adj;
  void loadGraphFromFile(const std::string &filename);
};
