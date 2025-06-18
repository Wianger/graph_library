#pragma once

#include <algorithm>
#include <limits>
#include <queue>
#include <vector>

// 模板化的 Dinic 類
template <typename T> class Dinic {
private:
  struct Edge {
    int to;
    T capacity;
    int rev; // 反向邊在鄰接列表中的索引
  };

public:
  Dinic(int v);
  void add_edge(int from, int to, T capacity);
  T max_flow(int s, int t);
  std::vector<int> get_min_cut_nodes(int s);
  std::vector<int> get_t_side_nodes(int t);
  std::vector<bool> get_reachable_nodes(int s);

private:
  int V;
  std::vector<std::vector<Edge>> adj;
  std::vector<int> level;
  std::vector<int> iter;
  bool bfs(int s, int t);
  T dfs(int u, int t, T f);
};

// Implementation below

template <typename T> inline T get_epsilon() { return 0; }
template <> inline double get_epsilon<double>() { return 1e-9; }
template <typename T> inline T get_infinity() {
  return std::numeric_limits<T>::max();
}
template <> inline double get_infinity<double>() {
  return std::numeric_limits<double>::infinity();
}

template <typename T> Dinic<T>::Dinic(int v) : V(v), adj(v) {}

template <typename T> void Dinic<T>::add_edge(int from, int to, T capacity) {
  adj[from].push_back({to, capacity, (int)adj[to].size()});
  adj[to].push_back({from, 0, (int)adj[from].size() - 1});
}

template <typename T> bool Dinic<T>::bfs(int s, int t) {
  level.assign(V, -1);
  std::queue<int> q;
  level[s] = 0;
  q.push(s);
  while (!q.empty()) {
    int u = q.front();
    q.pop();
    for (const auto &edge : adj[u]) {
      if (edge.capacity > get_epsilon<T>() && level[edge.to] < 0) {
        level[edge.to] = level[u] + 1;
        q.push(edge.to);
      }
    }
  }
  return level[t] != -1;
}

template <typename T> T Dinic<T>::dfs(int u, int t, T f) {
  if (u == t)
    return f;
  for (int &i = iter[u]; i < (int)adj[u].size(); ++i) {
    Edge &e = adj[u][i];
    if (e.capacity > get_epsilon<T>() && level[u] < level[e.to]) {
      T d = dfs(e.to, t, std::min(f, e.capacity));
      if (d > get_epsilon<T>()) {
        e.capacity -= d;
        adj[e.to][e.rev].capacity += d;
        return d;
      }
    }
  }
  return 0;
}

template <typename T> T Dinic<T>::max_flow(int s, int t) {
  T flow = 0;
  const T INF = get_infinity<T>();
  while (bfs(s, t)) {
    iter.assign(V, 0);
    T f;
    while ((f = dfs(s, t, INF)) > get_epsilon<T>()) {
      flow += f;
    }
  }
  return flow;
}

template <typename T> std::vector<int> Dinic<T>::get_min_cut_nodes(int s) {
  std::vector<int> s_nodes;
  std::vector<bool> visited(V, false);
  std::queue<int> q;
  q.push(s);
  visited[s] = true;
  while (!q.empty()) {
    int u = q.front();
    q.pop();
    s_nodes.push_back(u);
    for (const auto &edge : adj[u]) {
      if (edge.capacity > get_epsilon<T>() && !visited[edge.to]) {
        visited[edge.to] = true;
        q.push(edge.to);
      }
    }
  }
  return s_nodes;
}

template <typename T> std::vector<int> Dinic<T>::get_t_side_nodes(int t) {
  std::vector<std::vector<int>> rev_adj(V);
  for (int u = 0; u < V; ++u) {
    for (const auto &edge : adj[u]) {
      if (edge.capacity > get_epsilon<T>()) {
        rev_adj[edge.to].push_back(u);
      }
    }
  }
  std::vector<int> t_side_nodes;
  std::vector<bool> visited(V, false);
  std::queue<int> q;
  q.push(t);
  visited[t] = true;
  while (!q.empty()) {
    int u = q.front();
    q.pop();
    t_side_nodes.push_back(u);
    for (int v : rev_adj[u]) {
      if (!visited[v]) {
        visited[v] = true;
        q.push(v);
      }
    }
  }
  return t_side_nodes;
}

template <typename T> std::vector<bool> Dinic<T>::get_reachable_nodes(int s) {
  std::vector<bool> reachable(V, false);
  std::queue<int> q;
  q.push(s);
  reachable[s] = true;
  while (!q.empty()) {
    int u = q.front();
    q.pop();
    for (const auto &edge : adj[u]) {
      if (!reachable[edge.to] && edge.capacity > 0) {
        reachable[edge.to] = true;
        q.push(edge.to);
      }
    }
  }
  return reachable;
}
