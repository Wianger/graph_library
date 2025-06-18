// --- START OF FILE Dinic.h ---

#pragma once

#include <algorithm>
#include <limits>
#include <queue>
#include <vector>

// 模板化的 Edge 結構
template <typename T> struct Edge {
  int to;
  T capacity;
  int rev; // 反向邊在鄰接列表中的索引
};

// 模板化的 Dinic 類
template <typename T> class Dinic {
public:
  Dinic(int v);
  void add_edge(int from, int to, T capacity);
  T max_flow(int s, int t);
  std::vector<int> get_min_cut_nodes(int s);
  std::vector<int> get_t_side_nodes(int t);

private:
  const int V; // 頂點數量
  std::vector<std::vector<Edge<T>>> adj;
  std::vector<int> level;
  std::vector<int> iter; // 用於 DFS 的迭代器，防止重複遍歷

  bool bfs(int s, int t);
  T dfs(int v, int t, T f);
};
