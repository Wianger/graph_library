// --- START OF FILE Dinic.cpp ---

#include "Dinic.h"

template <typename T> T get_epsilon() { return 0; }
template <> double get_epsilon<double>() { return 1e-9; }
template <typename T> T get_infinity() { return std::numeric_limits<T>::max(); }
template <> double get_infinity<double>() {
  return std::numeric_limits<double>::infinity();
}

// --- 構造函數 ---
template <typename T> Dinic<T>::Dinic(int v) : V(v), adj(v) {}

// --- 添加邊 ---
template <typename T> void Dinic<T>::add_edge(int from, int to, T capacity) {
  // 正向邊
  adj[from].push_back({to, capacity, (int)adj[to].size()});
  // 反向邊（殘餘容量）
  adj[to].push_back({from, 0, (int)adj[from].size() - 1});
}

// --- BFS: 構建層級圖 ---
template <typename T> bool Dinic<T>::bfs(int s, int t) {
  level.assign(V, -1);
  std::queue<int> q;

  level[s] = 0;
  q.push(s);

  while (!q.empty()) {
    int u = q.front();
    q.pop();

    for (const auto &edge : adj[u]) {
      // 如果邊有剩餘容量且目標節點未被訪問
      if (edge.capacity > get_epsilon<T>() && level[edge.to] < 0) {
        level[edge.to] = level[u] + 1;
        q.push(edge.to);
      }
    }
  }
  // 如果能到達匯點 t，則返回 true
  return level[t] != -1;
}

// --- DFS: 在層級圖上尋找增廣路徑 ---
template <typename T> T Dinic<T>::dfs(int u, int t, T f) {
  // 到達匯點，返回流量
  if (u == t)
    return f;
  // 使用 iter 引用來避免重複訪問無效路徑（關鍵優化）
  for (int &i = iter[u]; i < adj[u].size(); ++i) {
    Edge<T> &e = adj[u][i];
    // 如果容量足夠且在下一層
    if (e.capacity > get_epsilon<T>() && level[u] < level[e.to]) {
      // 遞歸尋找路徑，流量上限是當前流量 f 和邊容量 e.capacity 的較小者
      T d = dfs(e.to, t, std::min(f, e.capacity));

      // 如果找到了增廣流量
      if (d > get_epsilon<T>()) {
        e.capacity -= d;                // 減去正向邊容量
        adj[e.to][e.rev].capacity += d; // 增加反向邊容量
        return d;
      }
    }
  }
  return 0;
}

// --- 主函數: 計算最大流 ---
template <typename T> T Dinic<T>::max_flow(int s, int t) {
  T flow = 0;
  const T INF = get_infinity<T>();

  // 當 BFS 仍然能找到從 s 到 t 的路徑時
  while (bfs(s, t)) {
    iter.assign(V, 0); // 重置 DFS 迭代器
    T f;
    // 在當前層級圖上持續尋找增廣路徑，直到找不到為止
    while ((f = dfs(s, t, INF)) > get_epsilon<T>()) {
      flow += f;
    }
  }
  return flow;
}

// --- 獲取最小割中 S 側的節點 ---
template <typename T> std::vector<int> Dinic<T>::get_min_cut_nodes(int s) {
  std::vector<int> s_nodes;
  // 在最終的殘餘圖上從 s 開始進行一次遍歷（這裡用BFS）
  // 找到所有 s 可達的節點
  std::vector<bool> visited(V, false);
  std::queue<int> q;

  q.push(s);
  visited[s] = true;

  while (!q.empty()) {
    int u = q.front();
    q.pop();
    s_nodes.push_back(u);

    for (const auto &edge : adj[u]) {
      // 如果邊有殘餘容量且目標節點未被訪問
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

// --- 顯式模板實例化 ---
// 這告訴編譯器生成 Dinic<long long> 和 Dinic<double> 的代碼。
// 如果您需要其他類型，可以在這裡添加。
template class Dinic<long long>;
template class Dinic<double>;
template class Dinic<int>;