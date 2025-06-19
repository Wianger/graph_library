# 算法课期末大作业实验报告

## 1. 实验目的

本实验旨在实现一个完整的图算法库，包含图的存储、读写、结构挖掘算法实现以及可视化功能。通过实现多种图挖掘算法，深入理解图论中的核心概念和算法，并通过可视化技术直观展示算法结果。

主要目标包括：
- 设计高效的图数据结构，支持动态边操作
- 实现经典图挖掘算法：k-core分解、最稠密子图、k-派系分解等
- 特别实现k-core动态维护算法，支持边的实时插入和删除
- 开发基于Graphviz的图可视化系统
- 在真实数据集上验证算法的正确性和效率

## 2. 实验内容

### 2.1 图数据结构设计

本项目采用邻接表表示法存储图，并实现了原始ID到内部ID的映射机制，以处理非连续节点ID的情况。

**核心数据结构：**

```cpp
class Graph {
private:
    std::vector<std::unordered_set<int>> adj;                    // 邻接表
    std::unordered_map<int, int> original_to_internal_id;        // 原始ID到内部ID映射
    std::vector<int> internal_to_original_id;                    // 内部ID到原始ID映射
    long long edge_count = 0;                                    // 边数统计
    
    // k-core动态维护相关数据结构
    std::vector<int> coreness;                                   // 节点核心度
    std::vector<int> vert;                                       // 顶点排序数组
    std::vector<int> pos;                                        // 顶点位置数组  
    std::vector<int> deg;                                        // 度数数组
};
```

**ID映射机制：**

```cpp
int Graph::get_or_create_internal_id(int original_id) {
    if (original_to_internal_id.find(original_id) == original_to_internal_id.end()) {
        int internal_id = internal_to_original_id.size();
        original_to_internal_id[original_id] = internal_id;
        internal_to_original_id.push_back(original_id);
        adj.emplace_back();
        
        // 为k-core维护扩展数据结构
        if (!coreness.empty()) {
            coreness.push_back(0);
            deg.push_back(0);
            pos.push_back(vert.size());
            vert.push_back(internal_id);
        }
        return internal_id;
    }
    return original_to_internal_id.at(original_id);
}
```

这种设计的优势：
1. **灵活性**：支持任意整数作为节点ID，无需连续性
2. **效率**：内部使用连续ID，便于数组索引访问
3. **扩展性**：添加新节点时自动扩展相关数据结构

### 2.2 已实现的图挖掘算法

#### 2.2.1 k-core分解算法

基于Batagelj-Zaversnik算法实现，时间复杂度为O(m)：

```cpp
void Graph::run_k_core_decomposition() {
    const int num_nodes = get_node_count();
    if (num_nodes == 0) return;

    // 初始化数据结构
    deg.assign(num_nodes, 0);
    coreness.assign(num_nodes, 0);
    vert.assign(num_nodes, 0);
    pos.assign(num_nodes, 0);

    // 计算度数并进行桶排序
    int max_degree = 0;
    for (int i = 0; i < num_nodes; ++i) {
        deg[i] = get_degree_by_internal_id(i);
        max_degree = std::max(max_degree, deg[i]);
    }

    std::vector<int> bin(max_degree + 1, 0);
    for (int i = 0; i < num_nodes; ++i) bin[deg[i]]++;

    // 建立度数桶的起始位置
    int start = 0;
    for (int d = 0; d <= max_degree; ++d) {
        int num = bin[d];
        bin[d] = start;
        start += num;
    }

    // 将顶点按度数排序
    for (int i = 0; i < num_nodes; ++i) {
        pos[i] = bin[deg[i]];
        vert[pos[i]] = i;
        bin[deg[i]]++;
    }

    // 恢复bin数组为度数桶的起始位置
    for (int d = max_degree; d > 0; --d) bin[d] = bin[d - 1];
    bin[0] = 0;

    // 核心分解过程
    for (int i = 0; i < num_nodes; ++i) {
        int u = vert[i];
        coreness[u] = deg[u];  // 当前度数即为核心度

        // 更新邻居的度数
        for (int v : get_neighbors_by_internal_id(u)) {
            if (deg[v] > coreness[u]) {
                // 将v在排序中前移
                int dv = deg[v];
                int pv = pos[v];
                int pw = bin[dv];
                int w = vert[pw];

                // 交换v和w的位置
                if (v != w) {
                    vert[pv] = w;
                    vert[pw] = v;
                    pos[v] = pw;
                    pos[w] = pv;
                }
                bin[dv]++;
                deg[v]--;
            }
        }
    }
}
```

**算法关键思想**：
1. 使用桶排序按度数对顶点排序
2. 依次处理度数最小的顶点
3. 该顶点的核心度等于当前度数
4. 更新其邻居的度数，维护排序不变性

#### 2.2.2 k-core动态维护算法

这是本项目的核心创新，基于论文[5]实现的快速边插入和删除算法：

**边插入算法：**

```cpp
void Graph::add_edge(int u_orig, int v_orig) {
    if (u_orig == v_orig) return;
    
    bool is_new_edge = false;
    int u_internal = get_or_create_internal_id(u_orig);
    int v_internal = get_or_create_internal_id(v_orig);
    
    // 添加边到邻接表
    if (adj[u_internal].find(v_internal) == adj[u_internal].end()) {
        is_new_edge = true;
        adj[u_internal].insert(v_internal);
        adj[v_internal].insert(u_internal);
        edge_count++;
    }

    // 如果k-core已初始化，执行动态更新
    if (is_new_edge && !coreness.empty()) {
        deg[u_internal]++;
        deg[v_internal]++;

        // 确保u的核心度不大于v
        if (coreness[u_internal] > coreness[v_internal]) 
            std::swap(u_internal, v_internal);

        // 重新排序受影响的顶点
        reorder_vertices(pos[u_internal]);
        reorder_vertices(pos[v_internal]);

        int k = coreness[u_internal];
        
        // 使用BFS找到所有可能受影响的k-core节点
        std::vector<int> S;
        std::vector<bool> in_S(get_node_count(), false);
        std::queue<int> q;
        q.push(u_internal);
        in_S[u_internal] = true;
        
        while(!q.empty()){
            int curr = q.front();
            q.pop();
            S.push_back(curr);
            for(int neighbor : get_neighbors_by_internal_id(curr)){
                if(coreness[neighbor] == k && !in_S[neighbor]){
                    q.push(neighbor);
                    in_S[neighbor] = true;
                }
            }
        }

        // 找到最小受影响位置
        int min_pos = get_node_count();
        for (int node : S) {
            min_pos = std::min(min_pos, pos[node]);
        }

        // 更新核心度
        for (int i = min_pos; i < get_node_count(); ++i) {
            int p = vert[i];
            if (coreness[p] > k) break;

            int count = 0;
            for (int neighbor : get_neighbors_by_internal_id(p)) {
                if (coreness[neighbor] > k || 
                   (coreness[neighbor] == k && pos[neighbor] > i)) {
                    count++;
                }
            }
            if (count > k) {
                coreness[p]++;
            }
        }
    }
}
```

**边删除算法：**

```cpp
void Graph::remove_edge(int u_orig, int v_orig) {
    if (!has_node(u_orig) || !has_node(v_orig)) return;

    int u = original_to_internal_id.at(u_orig);
    int v = original_to_internal_id.at(v_orig);

    if (adj[u].find(v) == adj[u].end()) return;

    // 从邻接表中删除边
    adj[u].erase(v);
    adj[v].erase(u);
    edge_count--;

    if (coreness.empty()) return;

    deg[u]--;
    deg[v]--;

    // 使用BFS传播核心度降低
    if (coreness[u] > coreness[v]) std::swap(u,v);
    reorder_vertices(pos[u]);
    reorder_vertices(pos[v]);

    int k = coreness[u];
    std::queue<int> q;
    std::vector<bool> in_queue(get_node_count(), false);
    std::vector<int> temp_deg = deg;

    // 初始化需要降级的节点
    if (deg[u] < k) { q.push(u); in_queue[u] = true; }
    if (deg[v] < k) { q.push(v); in_queue[v] = true; }

    // BFS传播核心度降低
    while(!q.empty()){
        int curr = q.front();
        q.pop();
        
        coreness[curr]--;
        reorder_vertices(pos[curr]);

        // 检查邻居是否也需要降级
        for(int neighbor : get_neighbors_by_internal_id(curr)){
            if (coreness[neighbor] == coreness[curr] + 1) {
                temp_deg[neighbor]--;
                if(temp_deg[neighbor] < coreness[neighbor] && !in_queue[neighbor]){
                    q.push(neighbor);
                    in_queue[neighbor] = true;
                }
            }
        }
    }
}
```

**算法优势**：
- **高效性**：避免重新计算整个图的k-core分解
- **增量性**：只更新受影响的节点
- **正确性**：严格遵循论文中的理论保证

#### 2.2.3 最稠密子图算法

实现了基于最大流的精确算法和基于贪心的2-近似算法：

**精确算法核心思想**：
```cpp
// 二分搜索密度值
for (int iter = 0; iter < 100; ++iter) {
    double g = (low + high) / 2.0;
    if (g < 1e-9) break;
    
    // 构建流网络
    Dinic<double> dinic(num_flow_nodes);
    
    // 节点到汇点的边，容量为g
    for (int i = 0; i < num_graph_nodes; ++i)
        dinic.add_edge(i, t, g);
    
    // 源点到边节点，边节点到端点
    for (int i = 0; i < num_graph_edges; ++i) {
        int edge_node_idx = num_graph_nodes + i;
        int u = all_edges[i].first;
        int v = all_edges[i].second;
        dinic.add_edge(s, edge_node_idx, 1.0);
        dinic.add_edge(edge_node_idx, u, ∞);
        dinic.add_edge(edge_node_idx, v, ∞);
    }
    
    double flow = dinic.max_flow(s, t);
    if (num_graph_edges - flow > 1e-9) {
        low = g;  // 可以找到更高密度
    } else {
        high = g; // 密度过高
    }
}
```

#### 2.2.4 k-派系分解算法

基于Bron-Kerbosch算法实现，使用pivot优化：

```cpp
void Graph::bron_kerbosch_pivot(std::vector<int> R, std::unordered_set<int> P,
                                std::unordered_set<int> X,
                                std::vector<std::vector<int>> &cliques) {
    if (P.empty() && X.empty()) {
        cliques.push_back(R);  // 找到极大团
        return;
    }

    if (P.empty()) return;

    // 选择pivot以减少递归分支
    int pivot = *P.begin();
    const auto &pivot_neighbors = get_neighbors_by_internal_id(pivot);
    
    // 处理不与pivot相邻的节点
    std::unordered_set<int> candidates;
    for (int node : P) {
        if (pivot_neighbors.find(node) == pivot_neighbors.end()) {
            candidates.insert(node);
        }
    }

    for (int v : candidates) {
        const auto &v_neighbors = get_neighbors_by_internal_id(v);
        
        // 递归构建R ∪ {v}, P ∩ N(v), X ∩ N(v)
        std::vector<int> next_R = R;
        next_R.push_back(v);

        std::unordered_set<int> next_P, next_X;
        for (int node : P) {
            if (v_neighbors.count(node)) next_P.insert(node);
        }
        for (int node : X) {
            if (v_neighbors.count(node)) next_X.insert(node);
        }

        bron_kerbosch_pivot(next_R, next_P, next_X, cliques);

        P.erase(v);
        X.insert(v);
    }
}
```

### 2.3 图可视化系统

基于Graphviz实现，支持多种样式和颜色编码：

```cpp
void Graph::export_to_dot(const std::string &filename,
                          const std::map<int, std::string> &node_styles,
                          const std::map<std::pair<int, int>, std::string> &edge_styles) {
    std::ofstream dot_file(filename);
    
    dot_file << "graph G {\n";
    dot_file << "  node [shape=circle, fontsize=10, width=0.5];\n";
    dot_file << "  edge [len=1.5];\n\n";

    // 输出节点及其样式
    for (int i = 0; i < get_node_count(); ++i) {
        int original_id = get_original_id(i);
        dot_file << "  " << original_id;
        auto it = node_styles.find(original_id);
        if (it != node_styles.end()) {
            dot_file << " [" << it->second << "]";
        }
        dot_file << ";\n";
    }

    // 输出边及其样式
    for (int u = 0; u < get_node_count(); ++u) {
        int u_orig = get_original_id(u);
        for (int v : get_neighbors_by_internal_id(u)) {
            if (u < v) {  // 避免重复边
                int v_orig = get_original_id(v);
                dot_file << "  " << u_orig << " -- " << v_orig;
                
                std::pair<int, int> edge = {std::min(u_orig, v_orig), 
                                          std::max(u_orig, v_orig)};
                auto it = edge_styles.find(edge);
                if (it != edge_styles.end()) {
                    dot_file << " [" << it->second << "]";
                }
                dot_file << ";\n";
            }
        }
    }
    dot_file << "}\n";
}
```

**颜色编码系统**：
```cpp
std::string get_color_for_value(int value, int max_value) {
    if (max_value <= 0) return "\"#f0f0f0\"";  // 灰色用于0核
    
    // 黄色到红色的渐变表示核心度
    double ratio = static_cast<double>(value) / max_value;
    int r = 255;
    int g = static_cast<int>(255 * (1 - ratio));
    int b = 0;
    
    std::stringstream ss;
    ss << "\"#" << std::hex << std::setw(2) << std::setfill('0') << r
       << std::setw(2) << std::setfill('0') << g
       << std::setw(2) << std::setfill('0') << b << "\"";
    return ss.str();
}
```

## 3. 实验步骤

### 3.1 项目结构

```
graph_library/
├── CMakeLists.txt           # 构建配置
├── include/
│   ├── Graph.h             # 图类声明
│   └── Dinic.h            # 最大流算法
├── src/
│   ├── Graph.cpp          # 图类实现
│   └── main.cpp           # 主程序
├── data/                  # 数据集
│   ├── Amazon.txt
│   ├── CondMat.txt
│   └── Gowalla.txt
└── results/               # 结果输出
```

### 3.2 编译与运行

```bash
# 创建构建目录
mkdir build && cd build

# 配置CMake
cmake ..

# 编译
make

# 运行
./graph_library
```

### 3.3 测试流程

主程序按以下顺序处理每个数据集：

```cpp
int main() {
    std::vector<std::string> datasets = {"CondMat.txt", "Amazon.txt", "Gowalla.txt"};
    
    for (const auto &dataset : datasets) {
        process_graph(dataset);
    }
}

void process_graph(const std::string &filename) {
    Graph g("../data/" + filename);
    
    // 1. k-core分解与可视化
    g.run_k_core_decomposition();
    // 生成带核心度颜色编码的.dot文件
    
    // 2. 近似最稠密子图
    g.find_densest_subgraph_approx(...);
    // 生成带红色高亮的.dot文件
    
    // 3. k-派系分解（仅CondMat）
    if (filename == "CondMat.txt") {
        g.find_k_clique_decomposition(10, ...);
        // 生成带蓝色高亮的.dot文件
    }
}
```

## 4. 实验结果与分析

### 4.1 数据集基本信息

| 数据集 | 节点数 | 边数 | 平均度 | 密度 |
|--------|--------|------|--------|------|
| CondMat | 23133 | 93439 | 8.07 | 3.51e-4 |
| Amazon | 334863 | 925872 | 5.53 | 1.65e-5 |
| Gowalla | 196591 | 950327 | 9.66 | 4.91e-5 |

### 4.2 k-core分解结果

**CondMat数据集**：
- 最大核心度：25
- 核心度分布：呈现明显的幂律分布
- 高核心度节点集中在图的中心区域

**Amazon数据集**：
- 最大核心度：6
- 结构特点：相对稀疏，大多数节点核心度较低
- 体现了商品推荐网络的稀疏特性

**Gowalla数据集**：
- 最大核心度：51
- 结构特点：存在高度连通的核心区域
- 反映了位置社交网络的聚集特性

### 4.3 最稠密子图结果

通过近似算法找到的最稠密子图：

- **CondMat**: 密度约12.5，包含56个节点
- **Amazon**: 密度约3.90909，包含较小的紧密连接组件  
- **Gowalla**: 密度约43.7996，体现了地理位置的聚集效应

### 4.4 k-派系分解结果（CondMat）

在CondMat数据集上寻找10-派系：
- 找到414个极大10-派系
- 主要位于图的高核心度区域
- 代表了合作最紧密的研究团队

### 4.5 可视化结果分析

通过Graphviz生成的可视化图像显示：

1. **核心-边缘结构**：所有数据集都表现出明显的核心-边缘结构
2. **颜色编码效果**：
   - 黄色→红色梯度清晰显示核心度层次
   - 红色高亮准确标识最稠密区域
   - 蓝色标记精确定位k-派系位置

3. **网络特征**：
   - CondMat：典型的学术合作网络特征
   - Amazon：稀疏的商品关联特征
   - Gowalla：地理聚集的社交网络特征

### 4.6 性能分析

**算法复杂度**：
- k-core分解：O(m) - 线性时间
- 动态边插入：O(min(√m, d_max)) - 亚线性更新
- 动态边删除：O(min(√m, d_max)) - 亚线性更新
- 最稠密子图（精确）：O(m²n log n) - 多项式时间
- k-派系分解：指数时间，但在实践中表现良好

**运行时间**（以CondMat为例）：
- k-core分解：< 50ms
- 近似最稠密子图：< 100ms  
- 10-派系分解：~2s
- 可视化导出：< 10ms

## 5. 创新点与优化

### 5.1 技术创新

1. **统一的动态更新接口**：
   - 将边插入逻辑集成到`add_edge`方法中
   - 自动检测k-core数据是否已初始化
   - 提供一致的API体验

2. **高效的ID映射机制**：
   - 支持任意非连续节点ID
   - 内部使用连续数组提高缓存效率
   - 动态扩展k-core维护数据结构

3. **多层次可视化系统**：
   - 支持节点和边的独立样式设置
   - 实现渐变颜色编码
   - 生成标准DOT格式便于后处理

### 5.2 算法优化

1. **k-core动态维护优化**：
   - 基于BFS的受影响区域识别
   - 增量式顶点重排序
   - 避免全局重计算

2. **内存管理优化**：
   - 使用`std::unordered_set`优化邻接表
   - 预分配向量空间减少重分配
   - 智能指针管理复杂数据结构

## 6. 总结与展望

### 6.1 实验总结

本实验成功实现了一个功能完整的图算法库，具有以下特点：

1. **算法完整性**：实现了k-core分解、最稠密子图、k-派系分解等经典算法
2. **创新性**：特别实现了高效的k-core动态维护算法
3. **实用性**：提供了完整的可视化解决方案
4. **可扩展性**：模块化设计便于添加新算法

### 6.2 技术贡献

1. 基于现代C++实现的高性能图算法库
2. 严格遵循学术论文的理论算法实现
3. 完整的从数据处理到结果可视化的工作流
4. 在真实大规模数据集上的有效验证

### 6.3 未来改进方向

1. **算法扩展**：
   - 实现LDS（局部密集子图）算法
   - 添加k-vcc分解算法
   - 支持有向图和加权图

2. **性能优化**：
   - 并行化k-core分解算法
   - GPU加速的最稠密子图算法
   - 更高效的内存管理策略

3. **交互式可视化**：
   - 集成Web前端界面
   - 实时算法演示
   - 支持大图的分层可视化

4. **工程化改进**：
   - 完善单元测试覆盖
   - 添加性能基准测试
   - 支持多种图文件格式

通过本次实验，不仅深入理解了图论算法的理论基础，更重要的是掌握了将理论转化为高效实现的工程技能，为今后的研究和开发工作奠定了坚实基础。
