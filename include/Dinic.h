#include <vector>

// --- Dinic's Max-Flow Algorithm Implementation ---
struct Edge {
  int to;
  long long capacity;
  int rev; // Reverse edge index in the adjacency list of 'to'
};

class Dinic {
private:
  int V;
  std::vector<std::vector<Edge>> adj;
  std::vector<int> level;
  std::vector<int> iter;

  bool bfs(int s, int t) {
    level.assign(V, -1);
    std::vector<int> q;
    q.push_back(s);
    level[s] = 0;
    int head = 0;
    while (head < q.size()) {
      int u = q[head++];
      for (const auto &edge : adj[u]) {
        if (edge.capacity > 0 && level[edge.to] < 0) {
          level[edge.to] = level[u] + 1;
          q.push_back(edge.to);
        }
      }
    }
    return level[t] != -1;
  }

  long long dfs(int u, int t, long long f) {
    if (u == t)
      return f;
    for (int &i = iter[u]; i < adj[u].size(); ++i) {
      Edge &e = adj[u][i];
      if (e.capacity > 0 && level[u] < level[e.to]) {
        long long d = dfs(e.to, t, std::min(f, e.capacity));
        if (d > 0) {
          e.capacity -= d;
          adj[e.to][e.rev].capacity += d;
          return d;
        }
      }
    }
    return 0;
  }

public:
  Dinic(int v) : V(v), adj(v) {}

  void add_edge(int from, int to, long long capacity) {
    adj[from].push_back({to, capacity, (int)adj[to].size()});
    adj[to].push_back({from, 0, (int)adj[from].size() - 1}); // Residual edge
  }

  long long max_flow(int s, int t) {
    long long flow = 0;
    while (bfs(s, t)) {
      iter.assign(V, 0);
      long long f;
      while ((f = dfs(s, t, -1)) > 0) // -1 represents infinity
        flow += f;
    }
    return flow;
  }

  // Get the nodes on the S-side of the min-cut after running max_flow
  std::vector<int> get_min_cut_nodes(int s) {
    std::vector<int> s_nodes;
    // level will be populated by the last bfs call in max_flow
    for (int i = 0; i < V; ++i)
      if (level[i] != -1)
        s_nodes.push_back(i);
    return s_nodes;
  }
};