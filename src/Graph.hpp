/**
 ********************************************
 * @file    :Graph.hpp
 * @author  :XXY
 * @brief   :Graph Structure And Operations
 * @date    :2025/11/30
 ********************************************
 */

#ifndef POP_GRAPH_HPP
#define POP_GRAPH_HPP

#include "CH.hpp"
#include "PL.hpp"
#include "common.hpp"
#include <boost/heap/fibonacci_heap.hpp>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullError.h>

/**
 * @brief 距离比较器，用于优先队列
 */
struct DistanceComparator {
    bool operator()(const std::tuple<unsigned, unsigned, unsigned>& lhs, const std::tuple<unsigned, unsigned, unsigned>& rhs) const {
        return std::get<0>(lhs) > std::get<0>(rhs);
    }
};

/**
 * @brief 图中的边
 */
typedef CH* Edge;

/**
 * @brief 图中的顶点
 */
struct Vertex {
    bool boundary = false;
    unsigned valid{};
    unsigned partition{};
    map<unsigned, Edge> edges{};
    vector<vector<unsigned>> dis2Par{};

    Vertex() = default;

    explicit Vertex(istringstream &in) {
        in >> boundary >> valid >> partition;
        unsigned edgesSize;
        in >> edgesSize;
        for (unsigned i = 0; i != edgesSize; ++i) {
            unsigned t;
            in >> t;
            edges[t] = new CH(in);
        }
    }

    void Print(ofstream &out) {
        out << boundary << ' ' << valid << ' ' << partition << ' ';
        out << edges.size() << ' ';
        for (const auto &d : edges) {
            out << d.first << ' ';
            d.second->Print(out);
        }
        out << endl;
    }
};

/**
 * @brief 图结构
 */
class Graph {
    friend class TD;

public:
    unsigned n{}, m{}, w{D};
    vector<Vertex> vertexList{};
    vector<vector<unsigned>> partitions;
    vector<set<unsigned>> borders;
private:
//    void SerialCliqueConstruction(unsigned v, const vector<unsigned>& neighbors);
//    void ParallelCliqueConstruction(unsigned v, const vector<unsigned>& neighbors);
//    void ConvertToMetisFormat(const string &dataPath,
//                              vector<idx_t> &xadj,
//                              vector<idx_t> &adjncy,
//                              vector<idx_t> &adjwgt);
//    void OptimizePartitionBoundaries();
//    void RecalculateBoundaries();
//    void GreedyRebalance(vector<idx_t>& part, unsigned nparts);

public:
    Graph(const string &dataPath, unsigned costNum, unsigned randomNum);

    explicit Graph(const string &borderGraph);

    void VertexContract(unsigned v, set<OrderedNode> &degreeOrder);

    [[nodiscard]] std::tuple<std::vector<double>, double, unsigned> Dijkstra(unsigned s, unsigned t, const vector<double> &weight) const;

//    [[nodiscard]] Path AStar(unsigned s,
//                             unsigned t,
//                             const vector<double> &weight,
//                             map<unsigned, CH*> &s2b,
//                             map<unsigned, CH*> &b2t) const;

//    std::tuple<std::vector<double>, double, double>
//    AStar(unsigned s, unsigned t, const std::vector<double> &weight,
//          std::map<unsigned, CH *> &s2b, std::map<unsigned, CH *> &b2t) const;

//    [[nodiscard]] Path AStar_pl(unsigned s,
//                             unsigned t,
//                             const vector<double> &weight,
//                             map<unsigned, CH*> &s2b,
//                             map<unsigned, CH*> &b2t) const;

    void SetParNum(const string &path);

    void PrepareAStar();

    void PrintBorder(const string &path);

    void PrintBorderGraphDot(const string &path);

    void PrintBorderGraph(const string &path);

    int ComputePointLocation(unsigned int p);

    // 返回 <naive_ms, pl_ms>，同一缓存热度
    std::pair<double,double> pairedMs(unsigned s, unsigned t,
                                      const std::vector<double>& w,
                                      std::map<unsigned,CH*>& s2b,
                                      std::map<unsigned,CH*>& b2t) const;

    Path AStarTest(unsigned s, unsigned t, const vector<double>& weight,
                   map<unsigned, CH*>& s2b, map<unsigned, CH*>& b2t, bool usePL) const;
    void SetParNumWithMetis(const string &dataPath, unsigned nodePerPartition);
};

/**
 * @brief graph constructor
 * @param dataPath the path to the dataset
 * @param costNum the number of cost to get for the dataset
 * @param randomNum the number of random cost to generate
 */
Graph::Graph(const string &dataPath, unsigned costNum, unsigned randomNum) {
    w = costNum + randomNum;
    ifstream data(dataPath);
    string line;
    unsigned s, t, c;

    // 读取顶点和边信息
    getline(data, line);
    istringstream row(line);
    row >> n >> m >> c;
    assert(costNum <= c);
    vertexList.resize(n);

    while (getline(data, line)) {
        row.clear();
        row.str(line);
        row >> s >> t;
        vector<unsigned> cost;
        for (unsigned i = 0; i != costNum; ++i) {
            row >> c;
            cost.push_back(c);
//            cost.push_back(1);
        }
        for (unsigned i = 0; i != randomNum; ++i) {
            cost.push_back(DISTRIBUTION(DRE));
//            cost.push_back(1);
        }

        if (vertexList[s].edges.find(t) == vertexList[s].edges.end()) {
            Edge newEdge = new CH(cost);
            vertexList[s].edges[t] = newEdge;
            vertexList[t].edges[s] = newEdge;
        }
    }

    data.close();
    Now("Graph Loaded");
}


Graph::Graph(const string &borderGraph) {
    ifstream in(borderGraph);
    string line;
    getline(in, line);
    istringstream row(line);
    row >> n >> m >> w;
    vertexList.resize(n);
    for (unsigned i = 0; i != n; ++i) {
        getline(in, line);
        row.clear();
        row.str(line);
        vertexList[i] = Vertex(row);
    }
    getline(in, line);
    row.clear();
    row.str(line);
    unsigned p, temp;
    row >> p;
    partitions.resize(p);
    for (unsigned i = 0; i != p; ++i) {
        getline(in, line);
        row.clear();
        row.str(line);
        while (row >> temp) {
            partitions[i].emplace_back(temp);
        }
    }
    getline(in, line);
    row.clear();
    row.str(line);
    row >> p;
    borders.resize(p);
    for (unsigned i = 0; i != p; ++i) {
        getline(in, line);
        row.clear();
        row.str(line);
        while (row >> temp) {
            borders[i].insert(temp);
        }
    }

    for (unsigned i = 0; i != n; ++i) {
        for (auto v : vertexList[i].edges) {
            if (vertexList[v.first].edges[i] != v.second) {
                delete vertexList[v.first].edges[i];
                vertexList[v.first].edges[i] = v.second;
            }
        }
    }
    in.close();
}


void Graph::VertexContract(unsigned v, set<OrderedNode> &degreeOrder) {
    // get neighbors
    vector<unsigned> neighbors, originDegrees;
    for (auto e : vertexList[v].edges) {
        neighbors.push_back(e.first);
    }
    for (auto neighbor : neighbors) {
        originDegrees.push_back(vertexList[neighbor].edges.size());
    }

    // form a clique of neighbors
    for (unsigned i = 0; i != neighbors.size(); ++i) {
        unsigned n1 = neighbors[i];
        for (unsigned j = i + 1; j != neighbors.size(); ++j) {
            unsigned n2 = neighbors[j];
            vector<Point> newPoints;
            LSCatenation(vertexList[v].edges[n1]->points, vertexList[v].edges[n2]->points, newPoints);
            if (vertexList[n1].edges.find(n2) != vertexList[n1].edges.end()) {
                vertexList[n1].edges[n2]->Insert(newPoints);
            } else {
                Edge newEdge = new CH(newPoints);
                vertexList[n1].edges[n2] = newEdge;
                vertexList[n2].edges[n1] = newEdge;
            }
        }
    }

    // remove edges & update degree order
    unsigned i = 0;
    for (auto neighbor : neighbors) {
        vertexList[neighbor].edges.erase(v);

        // Update degree order if the neighbor is not a boundary vertex
        if (!vertexList[neighbor].boundary) {
            degreeOrder.erase(degreeOrder.find({originDegrees[i], neighbor}));
            degreeOrder.insert({static_cast<unsigned>(vertexList[neighbor].edges.size()), neighbor});
        }
        ++i;
    }
    vertexList[v].edges.clear();
}

/**
 * @brief Dijkstra algorithm
 * @param s the source vertex
 * @param t the target vertex
 * @param weight the weight of each cost
 * @return the cost vector and the distance of the shortest path
 */
std::tuple<std::vector<double>, double, unsigned> Graph::Dijkstra(unsigned s, unsigned t, const std::vector<double> &weight) const {
    vector<bool> visitFlag(n, false);
    vector<double> distance(n, DBL_MAX);
    vector<unsigned> pathLog(n, UINT_MAX);
    distance[s] = 0;

    // 优先队列使用总代价排序（与A*一致）
    set<pair<double, unsigned>> distanceOrder;
    distanceOrder.insert({0, s});

    while (!distanceOrder.empty()) {
        auto v = distanceOrder.begin()->second;
        distanceOrder.erase(distanceOrder.begin());

        if (v == t) break;
        if (visitFlag[v]) continue;
        visitFlag[v] = true;

        for (auto e : vertexList[v].edges) {
            if (visitFlag[e.first]) continue;

            // 修复1：使用NaiveQuery动态获取当前权重下的最优代价 ------------
            double edgeCost = e.second->NaiveQuery(weight).second; // 动态选择最优点的总代价
            auto d = distance[v] + edgeCost; // ----------------------------

            if (d < distance[e.first]) {
                distanceOrder.erase({distance[e.first], e.first});
                distance[e.first] = d;
                distanceOrder.insert({d, e.first});
                pathLog[e.first] = v;
            }
        }
    }

    // 修复2：路径累加时使用最优点的各维度值 --------------------------
    vector<double> cost(D, 0.0);
    vector<unsigned> pathVertex;
    unsigned back = t;
    unsigned nodeCount = 0;

    while (back != s) {
        pathVertex.push_back(back);
        auto prev = pathLog[back];

        // 获取边的最优点的各维度值
        auto& edge = vertexList[back].edges.at(prev);
        auto optimalPoint = edge->NaiveQuery(weight).first; // 动态选择最优点

        for (unsigned i = 0; i < D; i++) {
            cost[i] += optimalPoint[i]; // 累加实际最优点的各维度值
        }

        back = prev;
        nodeCount++;
    }

    pathVertex.push_back(s);
    reverse(pathVertex.begin(), pathVertex.end());

    return make_tuple(cost, distance[t], nodeCount + 1);
}



//Path Graph::BorderDijkstra(unsigned s,
//                           unsigned t,
//                           const vector<double> &weight,
//                           map<unsigned, CH *> &s2b,
//                           map<unsigned, CH *> &b2t) const {
//    // initialization
//    vector<bool> visitFlag(n, false);
//    vector<double> distance(n, DBL_MAX);
//    vector<unsigned> pathLog(n, UINT_MAX);
//    distance[s] = 0;
//
//    // Dijkstra with priority queue(pair<distance, vertex>)
//    set<pair<double, unsigned>> distanceOrder;
//    distanceOrder.insert({0, s});
//    while (!distanceOrder.empty()) {
//        auto v = distanceOrder.begin()->second;
//        distanceOrder.erase(distanceOrder.begin());
//        if (v == t) {
//            break;
//        }
//        if (visitFlag[v]) {
//            continue;
//        }
//        visitFlag[v] = true;
//        auto neighbors = vertexList[v].edges;
//        if (v == s) {
//            neighbors = s2b;
//        } else if (b2t.find(v) != b2t.end()) {
//            neighbors[t] = b2t[v];
//        }
//        for (auto e : neighbors) {
//            if (visitFlag[e.first]) {
//                continue;
//            }
//            auto d = distance[v] + e.second->NaiveQuery(weight).second;
//            if (d < distance[e.first]) {
//                if (distanceOrder.find({distance[e.first], e.first}) != distanceOrder.end()) {
//                    distanceOrder.erase({distance[e.first], e.first});
//                }
//                distance[e.first] = d;
//                distanceOrder.insert({d, e.first});
//                pathLog[e.first] = v;
//            }
//        }
//    }
//
//    // get path
//    //TODO: return cost
//    vector<double> cost(D, 0.0); // 使用 double 类型
///*    vector<unsigned> pathVertex;
//    pathVertex.emplace_back(t);
//    unsigned back = t;
//    while(back != s){
//        cost = VectorAdd(cost, (vertexList[back].edges.find(pathLog[back]))->second->points[0]);
//        back = pathLog[back];
//        pathVertex.emplace_back(back);
//    }
//    pathVertex.emplace_back(s);
//    reverse(pathVertex.begin(), pathVertex.end());*/
//    return make_pair(cost, distance[t]);
//}

//Path Graph::AStar(unsigned s,
//                  unsigned t,
//                  const vector<double> &weight,
//                  map<unsigned, CH *> &s2b,
//                  map<unsigned, CH *> &b2t) const {
//    // 获取目标点的分区号
//    auto pt = vertexList[t].partition;
//
//    // 初始化 OpenSet 优先队列 (用于存储 fScore 和对应的节点)
//    set<pair<double, unsigned>> openSet;
//    openSet.insert(make_pair(0, s));
//
//    // 路径记录数组
//    vector<unsigned> pathLog(n, UINT_MAX);
//
//    // 初始化 gScore 和 fScore
//    vector<double> gScore(n, DBL_MAX);
//    gScore[s] = 0;
//
//    vector<double> fScore(n, DBL_MAX);
//    fScore[s] = 0;
//
//    // A* 搜索主循环
//    while (!openSet.empty()) {
//        auto v = openSet.begin()->second;  // 当前节点
//        openSet.erase(openSet.begin());
//
//        if (v == t) {
//            break; // 找到目标节点，退出
//        }
//
//        // 获取邻居节点
//        auto neighbors = vertexList[v].edges;
//        if (v == s) {
//            neighbors = s2b; // 起始点的邻居
//        } else if (b2t.find(v) != b2t.end()) {
//            neighbors[t] = b2t[v]; // 终点的邻居
//        }
//
//        for (auto e : neighbors) {
//            // 使用点定位
////            double edgeCost = (e.second->pl != nullptr)
////                              ? e.second->pl->PLQuery(weight).second // 如果存在 pl 结构，使用 PLQuery
////                              : e.second->NaiveQuery(weight).second; // 如果没有 pl 结构，使用 NaiveQuery
////            auto tentativeGScore = gScore[v] + edgeCost; // 计算 gScore
//
//            //不使用点定位
//            auto tentativeGScore = gScore[v] + e.second->NaiveQuery(weight).second;
//
//            if (tentativeGScore < gScore[e.first]) {
//                // 更新 OpenSet 中的节点
//                if (openSet.find({fScore[e.first], e.first}) != openSet.end()) {
//                    openSet.erase({fScore[e.first], e.first});
//                }
//
//                // 更新 gScore 和 fScore
//                gScore[e.first] = tentativeGScore;
//                fScore[e.first] = tentativeGScore + DotProduct(weight,
//                                                               std::vector<double>(vertexList[e.first].dis2Par[pt].begin(), vertexList[e.first].dis2Par[pt].end()));
//
//                pathLog[e.first] = v;
//
//                // 插入新的 fScore 和对应节点
//                openSet.insert({fScore[e.first], e.first});
//            }
//        }
//    }
//    vector<double> cost(D, 0); // 假设返回的 cost 矢量
//    return make_pair(cost, fScore[t]); // 返回路径代价和最终 fScore
//}

//std::tuple<std::vector<double>, double, double>
//Graph::AStar(unsigned s, unsigned t, const vector<double> &weight,
//             map<unsigned, CH *> &s2b, map<unsigned, CH *> &b2t) const {
//    // 获取目标点的分区号
//    auto pt = vertexList[t].partition;
//
//    // 初始化 OpenSet 优先队列 (用于存储 fScore 和对应的节点)
//    set<pair<double, unsigned>> openSet;
//    openSet.insert(make_pair(0, s));
//
//    // 路径记录数组
//    vector<unsigned> pathLog(n, UINT_MAX);
//
//    // 初始化 gScore 和 fScore
//    vector<double> gScore(n, DBL_MAX);
//    gScore[s] = 0;
//
//    vector<double> fScore(n, DBL_MAX);
//    fScore[s] = 0;
//
//    // A* 搜索主循环
//    while (!openSet.empty()) {
//        auto v = openSet.begin()->second;  // 当前节点
//        openSet.erase(openSet.begin());
//
//        if (v == t) {
//            break; // 找到目标节点，退出
//        }
//
//        // 获取邻居节点
//        auto neighbors = vertexList[v].edges;
//        if (v == s) {
//            neighbors = s2b; // 起始点的邻居
//        } else if (b2t.find(v) != b2t.end()) {
//            neighbors[t] = b2t[v]; // 终点的邻居
//        }
//
//        for (auto e : neighbors) {
//            // 不使用点定位
//            auto tentativeGScore = gScore[v] + e.second->NaiveQuery(weight).second;
//
//            if (tentativeGScore < gScore[e.first]) {
//                // 更新 OpenSet 中的节点
//                if (openSet.find({fScore[e.first], e.first}) != openSet.end()) {
//                    openSet.erase({fScore[e.first], e.first});
//                }
//
//                // 更新 gScore 和 fScore
//                gScore[e.first] = tentativeGScore;
//                fScore[e.first] = tentativeGScore + DotProduct(weight,
//                                                               std::vector<double>(vertexList[e.first].dis2Par[pt].begin(), vertexList[e.first].dis2Par[pt].end()));
//
//                pathLog[e.first] = v;
//
//                // 插入新的 fScore 和对应节点
//                openSet.insert({fScore[e.first], e.first});
//            }
//        }
//    }
//
//    // 重建路径并计算平均向量数量
//    vector<double> cost(D, 0.0);
//    unsigned back = t;
//    unsigned totalVectors = 0;
//    unsigned edgeCount = 0;
//
//    // 创建一个映射来跟踪特殊边
//    map<pair<unsigned, unsigned>, CH*> specialEdges;
//
//    // 填充特殊边映射
//    for (const auto& edge : s2b) {
//        specialEdges[{s, edge.first}] = edge.second;
//    }
//    for (const auto& edge : b2t) {
//        specialEdges[{edge.first, t}] = edge.second;
//    }
//
//    while (back != s) {
//        auto prev = pathLog[back];
//
//        CH* edge = nullptr;
//
//        // 首先检查是否是特殊边
//        auto specialEdgeKey = make_pair(prev, back);
//        if (specialEdges.find(specialEdgeKey) != specialEdges.end()) {
//            edge = specialEdges[specialEdgeKey];
//        }
//            // 如果不是特殊边，检查常规边
//        else if (vertexList[back].edges.find(prev) != vertexList[back].edges.end()) {
//            edge = vertexList[back].edges.at(prev);
//        }
//        else {
//            // 如果既不是特殊边也不是常规边，尝试反向查找
//            specialEdgeKey = make_pair(back, prev);
//            if (specialEdges.find(specialEdgeKey) != specialEdges.end()) {
//                edge = specialEdges[specialEdgeKey];
//            }
//            else {
//                edge = vertexList[prev].edges.at(back);
//            }
//        }
//
//        // 统计向量数量
//        totalVectors += edge->points.size();
//        edgeCount++;
//
//        // 获取最优点的各维度值
//        auto optimalPoint = edge->NaiveQuery(weight).first;
//        for (unsigned i = 0; i < D; i++) {
//            cost[i] += optimalPoint[i];
//        }
//
//        back = prev;
//    }
//
//    // 计算平均向量数量
//    double avgVectors = edgeCount > 0 ? static_cast<double>(totalVectors) / edgeCount : 0.0;
//
//    return make_tuple(cost, fScore[t], avgVectors);
//}

//Path Graph::AStar_pl(unsigned s,
//                  unsigned t,
//                  const vector<double> &weight,
//                  map<unsigned, CH *> &s2b,
//                  map<unsigned, CH *> &b2t) const {
//    // 获取目标点的分区号
//    auto pt = vertexList[t].partition;
//
//    // 初始化 OpenSet 优先队列 (用于存储 fScore 和对应的节点)
//    set<pair<double, unsigned>> openSet;
//    openSet.insert(make_pair(0, s));
//
//    // 路径记录数组
//    vector<unsigned> pathLog(n, UINT_MAX);
//
//    // 初始化 gScore 和 fScore
//    vector<double> gScore(n, DBL_MAX);
//    gScore[s] = 0;
//
//    vector<double> fScore(n, DBL_MAX);
//    fScore[s] = 0;
//
//    // A* 搜索主循环
//    while (!openSet.empty()) {
//        auto v = openSet.begin()->second;  // 当前节点
//        openSet.erase(openSet.begin());
//
//        if (v == t) {
//            break; // 找到目标节点，退出
//        }
//
//        // 获取邻居节点
//        auto neighbors = vertexList[v].edges;
//        if (v == s) {
//            neighbors = s2b; // 起始点的邻居
//        } else if (b2t.find(v) != b2t.end()) {
//            neighbors[t] = b2t[v]; // 终点的邻居
//        }
//
//        for (auto e : neighbors) {
//            // 使用点定位
////            double edgeCost = (e.second->pl != nullptr)
////                              ? e.second->pl->PLQuery(weight).second // 如果存在 pl 结构，使用 PLQuery
////                              : e.second->NaiveQuery(weight).second; // 如果没有 pl 结构，使用 NaiveQuery
////            auto tentativeGScore = gScore[v] + edgeCost; // 计算 gScore
//
//            //不使用点定位
//            auto tentativeGScore = gScore[v] + e.second->NaiveQuery(weight).second;
//
//            if (tentativeGScore < gScore[e.first]) {
//                // 更新 OpenSet 中的节点
//                if (openSet.find({fScore[e.first], e.first}) != openSet.end()) {
//                    openSet.erase({fScore[e.first], e.first});
//                }
//
//                // 更新 gScore 和 fScore
//                gScore[e.first] = tentativeGScore;
//                fScore[e.first] = tentativeGScore + DotProduct(weight,
//                                                               std::vector<double>(vertexList[e.first].dis2Par[pt].begin(), vertexList[e.first].dis2Par[pt].end()));
//
//                pathLog[e.first] = v;
//
//                // 插入新的 fScore 和对应节点
//                openSet.insert({fScore[e.first], e.first});
//            }
//        }
//    }
//    vector<double> cost(D, 0); // 假设返回的 cost 矢量
//    return make_pair(cost, fScore[t]); // 返回路径代价和最终 fScore
//}

///**
// * @brief print the graph to file
// * @param path the path of the output file
// */
//void Graph::PrintGraph(const string &path){
//    ofstream out(path);
//    out << n << ' ' << m << ' ' << w << endl;
//    for(const auto& v : vertexList){
//        for(const auto& e : v.edges){
//            out << "[" << e.first << ":";
//            for(const auto& p : e.second->points){
//                for(const auto& c : p){
//                    out << ' ' << c;
//                }
//                out << '|';
//            }
//            out << "] ";
//        }
//        out << endl;
//    }
//}

/**
 * @brief set the partition number of each vertex and record borders
 * @param path the path of the partition result
 */
void Graph::SetParNum(const string &path) {
    // read partition result
    ifstream par(path);
    string line, temp;
    unsigned parNum, vid;

    getline(par, line);
    istringstream ss(line);
    ss >> parNum;
    partitions.resize(parNum);
    borders.resize(parNum);

    while (getline(par, line)) {
        ss.clear();
        ss.str(line);
        ss >> temp;
        assert(temp == "partition");
        ss >> parNum;
        getline(par, line);
        istringstream ss2(line);
        while (ss2 >> vid) {
            vertexList[vid - 1].partition = parNum;
            partitions[parNum].push_back(vid - 1);
        }
    }
    par.close();

    // record borders
    for (unsigned i = 0; i != n; ++i) {
        for (auto e : vertexList[i].edges) {
            if (vertexList[i].partition != vertexList[e.first].partition) {
                borders[vertexList[i].partition].insert(i);
                vertexList[i].boundary = true;
                break;
            }
        }
    }
}

//void Graph::PrepareAStar() {
//    Now("Preparing For A*...");
//    for (const auto &v : vertexList) {
//        if (v.boundary) {
//            for (auto e : v.edges) {
//                e.second->ComputeExPoints();
//            }
//        }
//    }
//
//    unsigned num = 0;
//    unsigned total = 0;
//    for (const auto &b : borders) {
//        total += b.size();
//    }
//
//    // 使用 Fibonacci 堆优化
//    using DistanceState = std::tuple<unsigned, unsigned, unsigned>;
//    boost::heap::fibonacci_heap<DistanceState, boost::heap::compare<DistanceComparator>> pq;
//
//    for (unsigned ver = 0; ver != n; ++ver) {
//        if (vertexList[ver].boundary) {
//            // 初始化每个边界节点到各个分区的距离
//            for (unsigned i = 0; i != partitions.size(); ++i) {
//                vertexList[ver].dis2Par.emplace_back(D, UINT_MAX);
//            }
//
//            for (unsigned i = 0; i != D; ++i) {
//                std::vector<bool> visitFlag(n, false);
//                std::vector<unsigned> distance(n, UINT_MAX);
//                distance[ver] = 0;
//                vertexList[ver].dis2Par[vertexList[ver].partition][i] = 0;
//
//                // 使用 Fibonacci 堆初始化优先队列
//                pq.emplace(std::make_tuple(0, ver, i));
//                unsigned count = 1;
//
//                // 广度优先搜索（BFS）
//                while (!pq.empty() && count != partitions.size()) {
//                    auto [currentDist, currentNode, dim] = pq.top();
//                    pq.pop();
//
//                    if (visitFlag[currentNode]) {
//                        continue;
//                    }
//                    visitFlag[currentNode] = true;
//
//                    // 更新 dis2Par 中的最短距离
//                    if (currentDist < vertexList[ver].dis2Par[vertexList[currentNode].partition][dim]) {
//                        ++count;
//                        vertexList[ver].dis2Par[vertexList[currentNode].partition][dim] = currentDist;
//                    }
//
//                    // 遍历邻接节点
//                    for (auto &e : vertexList[currentNode].edges) {
//                        unsigned neighborID = e.first;
//                        if (visitFlag[neighborID]) {
//                            continue;
//                        }
//                        unsigned newDist = currentDist + e.second->exPoints[dim].second;
//                        if (newDist < distance[neighborID]) {
//                            distance[neighborID] = newDist;
//                            pq.emplace(std::make_tuple(newDist, neighborID, dim));
//                        }
//                    }
//                }
//            }
//            if (num % 100 == 0) { // 每处理100个点输出一次
//                Now("Vertex " + std::to_string(ver) + " is prepared, " + std::to_string(++num) + '/' + std::to_string(total), true);
//            }
//            else {
//                ++num; // 其他情况下，正常增加 num
//            }
//        }
//    }
//}

void Graph::PrepareAStar() {
    Now("Preparing For A*...");

    // 阶段1: 串行计算所有边界节点的 exPoints（避免并发问题）
    for (size_t v = 0; v < vertexList.size(); ++v) {
        if (vertexList[v].boundary) {
            for (auto& e : vertexList[v].edges) {
                e.second->ComputeExPoints();
            }
        }
    }

    // 阶段2: 并行预计算 dis2Par
    unsigned total = 0;
    for (const auto& b : borders) total += b.size();

#pragma omp parallel for schedule(dynamic)
    for (unsigned ver = 0; ver < n; ++ver) {
        if (!vertexList[ver].boundary) continue;

        // 每个线程维护独立的优先队列和临时数据
        using DistanceState = std::tuple<unsigned, unsigned, unsigned>;
        boost::heap::fibonacci_heap<DistanceState, boost::heap::compare<DistanceComparator>> pq;

        // 初始化 dis2Par（无需锁，每个 ver 独立）
        vertexList[ver].dis2Par.resize(partitions.size());
        for (auto& vec : vertexList[ver].dis2Par) {
            vec.resize(D, UINT_MAX);
        }

        // 遍历每个维度
        for (unsigned dim = 0; dim < D; ++dim) {
            std::vector<bool> visitFlag(n, false);
            std::vector<unsigned> distance(n, UINT_MAX);
            distance[ver] = 0;
            vertexList[ver].dis2Par[vertexList[ver].partition][dim] = 0;

            pq.emplace(std::make_tuple(0, ver, dim));
            unsigned count = 1;

            // BFS 主循环
            while (!pq.empty() && count < partitions.size()) {
                auto [currentDist, currentNode, currentDim] = pq.top();
                pq.pop();

                if (visitFlag[currentNode]) continue;
                visitFlag[currentNode] = true;

                // 更新 dis2Par
                if (currentDist < vertexList[ver].dis2Par[vertexList[currentNode].partition][currentDim]) {
                    vertexList[ver].dis2Par[vertexList[currentNode].partition][currentDim] = currentDist;
                    ++count;
                }

                // 处理邻居
                for (auto& e : vertexList[currentNode].edges) {
                    unsigned neighborID = e.first;
                    if (visitFlag[neighborID]) continue;

                    unsigned newDist = currentDist + e.second->exPoints[dim].second;
                    if (newDist < distance[neighborID]) {
                        distance[neighborID] = newDist;
                        pq.emplace(std::make_tuple(newDist, neighborID, dim));
                    }
                }
            }
            pq.clear(); // 清空队列以复用
        }

        // 进度输出（线程安全）
#pragma omp critical
        {
            static unsigned num = 0;
            if (++num % 100 == 0) {
                Now("Processed " + std::to_string(num) + '/' + std::to_string(total) + " boundary vertices", true);
            }
        }
    }
}

/**
 * @brief print the border of each partition to file
 * @param path the path of the output file
 */
void Graph::PrintBorder(const string &path){
    ofstream out(path);
    out << "There are " << borders.size() << " partitions" << endl;
    unsigned total = 0;
    for(unsigned i = 0; i != borders.size(); ++i){
        total += borders[i].size();
        out << "partition " << i << " has " << borders[i].size() << " borders with " << partitions[i].size() << " vertices" << endl;
    }
    out << total << " borders in total" << endl;
}

void Graph::PrintBorderGraphDot(const string &path) {
    ofstream out(path);
    vector<set<unsigned>> connected(partitions.size());
    out << "graph BG{" << endl << '\t' << "node[shape = circle];" << endl;
    for (unsigned i = 0; i != n; ++i) {
        if (vertexList[i].boundary) {
            for (auto e : vertexList[i].edges) {
                if (vertexList[e.first].boundary) {
                    connected[vertexList[i].partition].insert(vertexList[e.first].partition);
                }
            }
        }
    }
    for (auto i = 0; i != connected.size(); ++i) {
        for (auto j : connected[i]) {
            if (i < j) {
                out << '\t' << i << "--" << j << ';' << endl;
            }
        }
    }
    out << '}';
}

void Graph::PrintBorderGraph(const string &path) {
    ofstream out(path);
    out << n << ' ' << m << ' ' << w << endl;
    for (unsigned i = 0; i != n; ++i) {
        vertexList[i].Print(out);
    }
    out << partitions.size() << endl;
    for (auto &partition : partitions) {
        for (auto v : partition) {
            out << v << ' ';
        }
        out << endl;
    }
    out << borders.size() << endl;
    for (auto &border : borders) {
        for (auto v : border) {
            out << v << ' ';
        }
        out << endl;
    }
    out.close();
}

int Graph::ComputePointLocation(unsigned int p) {
    std::vector<CH*> edgesToProcess;

    // 阶段1：收集需要构建 PL 的边（串行）
    for (auto v : partitions[p]) {
        if (vertexList[v].boundary) {
            for (auto& neighbor : vertexList[v].edges) {
                if (neighbor.second->points.size() > 50 && !neighbor.second->pl) {
                    edgesToProcess.push_back(neighbor.second);
                }
            }
        }
    }

    // 阶段2：并行构建 PL 对象，减少锁竞争
    std::vector<std::unique_ptr<PL>> plObjects(edgesToProcess.size());

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < edgesToProcess.size(); ++i) {
        // 不持有任何锁的情况下创建 PL 对象
        plObjects[i] = std::make_unique<PL>(edgesToProcess[i]->points);
    }

    // 阶段3：分配 PL 对象（串行，避免竞争条件）
    int plCount = 0;
    for (size_t i = 0; i < edgesToProcess.size(); ++i) {
        if (!edgesToProcess[i]->pl) {  // 双重检查，防止并发调用
            edgesToProcess[i]->pl = std::move(plObjects[i]);
            plCount++;
        }
    }

    return plCount;
}


std::pair<double, double> Graph::pairedMs(unsigned s, unsigned t,
                                          const std::vector<double>& w,
                                          std::map<unsigned, CH*>& s2b,
                                          std::map<unsigned, CH*>& b2t) const {
    using clk = std::chrono::steady_clock;

    // 预热运行，确保缓存状态一致
    AStarTest(s, t, w, s2b, b2t, false);
    AStarTest(s, t, w, s2b, b2t, true);

    // 多次测量取平均值
    const int numRuns = 3;
    double totalNaive = 0.0;
    double totalPL = 0.0;

    for (int i = 0; i < numRuns; i++) {
        // 交替顺序以消除缓存效应
        if (i % 2 == 0) {
            // 先测 Naive
            auto t0 = clk::now();
            AStarTest(s, t, w, s2b, b2t, false);
            auto t1 = clk::now();
            totalNaive += std::chrono::duration<double, std::milli>(t1 - t0).count();

            // 再测 PL
            auto t2 = clk::now();
            AStarTest(s, t, w, s2b, b2t, true);
            auto t3 = clk::now();
            totalPL += std::chrono::duration<double, std::milli>(t3 - t2).count();
        } else {
            // 先测 PL
            auto t0 = clk::now();
            AStarTest(s, t, w, s2b, b2t, true);
            auto t1 = clk::now();
            totalPL += std::chrono::duration<double, std::milli>(t1 - t0).count();

            // 再测 Naive
            auto t2 = clk::now();
            AStarTest(s, t, w, s2b, b2t, false);
            auto t3 = clk::now();
            totalNaive += std::chrono::duration<double, std::milli>(t3 - t2).count();
        }
    }

    return { totalNaive / numRuns, totalPL / numRuns };
}

Path Graph::AStarTest(unsigned s, unsigned t, const vector<double>& weight,
                      map<unsigned, CH*>& s2b, map<unsigned, CH*>& b2t, bool usePL) const {
    // 获取目标点的分区号
    auto pt = vertexList[t].partition;

    // 初始化 OpenSet 优先队列
    set<pair<double, unsigned>> openSet;
    openSet.insert(make_pair(0, s));

    // 路径记录数组
    vector<unsigned> pathLog(n, UINT_MAX);

    // 初始化 gScore 和 fScore
    vector<double> gScore(n, DBL_MAX);
    gScore[s] = 0;

    vector<double> fScore(n, DBL_MAX);
    fScore[s] = 0;

    // A* 搜索主循环
    while (!openSet.empty()) {
        auto v = openSet.begin()->second;
        openSet.erase(openSet.begin());

        if (v == t) {
            break;
        }

        // 获取邻居节点
        auto neighbors = vertexList[v].edges;
        if (v == s) {
            neighbors = s2b;
        } else if (b2t.find(v) != b2t.end()) {
            neighbors[t] = b2t[v];
        }

        for (auto e : neighbors) {
            double edgeCost;

            // 根据 usePL 参数选择查询方式
            if (usePL && e.second->pl != nullptr) {
//                edgeCost = e.second->NaiveQuery(weight).second;
                edgeCost = e.second->pl->PLQuery(weight).second;
            } else {
                edgeCost = e.second->NaiveQuery(weight).second;
            }

            auto tentativeGScore = gScore[v] + edgeCost;

            if (tentativeGScore < gScore[e.first]) {
                if (openSet.find({fScore[e.first], e.first}) != openSet.end()) {
                    openSet.erase({fScore[e.first], e.first});
                }

                gScore[e.first] = tentativeGScore;
                fScore[e.first] = tentativeGScore + DotProduct(weight,
                                                               std::vector<double>(vertexList[e.first].dis2Par[pt].begin(),
                                                                                   vertexList[e.first].dis2Par[pt].end()));

                pathLog[e.first] = v;
                openSet.insert({fScore[e.first], e.first});
            }
        }
    }

    vector<double> cost(D, 0);
    return make_pair(cost, fScore[t]);
}
//
//void Graph::SerialCliqueConstruction(unsigned v, const vector<unsigned>& neighbors) {
//    for (unsigned i = 0; i != neighbors.size(); ++i) {
//        unsigned n1 = neighbors[i];
//        for (unsigned j = i + 1; j != neighbors.size(); ++j) {
//            unsigned n2 = neighbors[j];
//            vector<Point> newPoints;
//            LSCatenation(vertexList[v].edges[n1]->points, vertexList[v].edges[n2]->points, newPoints);
//            if (vertexList[n1].edges.find(n2) != vertexList[n1].edges.end()) {
//                vertexList[n1].edges[n2]->Insert(newPoints);
//            } else {
//                Edge newEdge = new CH(newPoints);
//                vertexList[n1].edges[n2] = newEdge;
//                vertexList[n2].edges[n1] = newEdge;
//            }
//        }
//    }
//}
//
//void Graph::ParallelCliqueConstruction(unsigned v, const vector<unsigned>& neighbors) {
//    const size_t totalPairs = (neighbors.size() * (neighbors.size() - 1)) / 2;
//
//    // 使用 tuple 减少内存分配，直接存储结果
//    vector<tuple<unsigned, unsigned, vector<Point>>> edgeData;
//    edgeData.reserve(totalPairs);
//
//    // 预填充所有边对信息
//    for (unsigned i = 0; i < neighbors.size(); ++i) {
//        for (unsigned j = i + 1; j < neighbors.size(); ++j) {
//            edgeData.emplace_back(neighbors[i], neighbors[j], vector<Point>());
//        }
//    }
//
//    // 并行计算 LSCatenation，使用静态调度减少调度开销
//#pragma omp parallel for schedule(static)
//    for (size_t idx = 0; idx < edgeData.size(); ++idx) {
//        auto& [n1, n2, newPoints] = edgeData[idx];
//        LSCatenation(vertexList[v].edges[n1]->points,
//                     vertexList[v].edges[n2]->points,
//                     newPoints);
//    }
//
//    // 串行部分：更新图结构
//    for (const auto& [n1, n2, newPoints] : edgeData) {
//        if (vertexList[n1].edges.find(n2) != vertexList[n1].edges.end()) {
//            vertexList[n1].edges[n2]->Insert(const_cast<vector<Point>&>(newPoints));
//        } else {
//            Edge newEdge = new CH(newPoints);
//            vertexList[n1].edges[n2] = newEdge;
//            vertexList[n2].edges[n1] = newEdge;
//        }
//    }
//}
//
///**
// * @brief 使用METIS递归二分算法进行图分区，优化边界节点数量
// * @param dataPath 图数据文件路径
// * @param nodePerPartition 每个分区期望的节点数
// */
//void Graph::SetParNumWithMetis(const string &dataPath, unsigned nodePerPartition) {
//    Now("开始使用METIS图分区（针对困难图结构优化）...");
//
//    // 第一步：读取图并转换为METIS格式
//    vector<idx_t> xadj, adjncy, adjwgt;
//    ConvertToMetisFormat(dataPath, xadj, adjncy, adjwgt);
//
//    // 修正：计算合理的分区数
//    unsigned nparts = n / nodePerPartition;
//    nparts = max(2u, nparts);
//    nparts = min(nparts, n / 50);
//    nparts = min(nparts, 1000u);
//
//    Now("图转换完成，节点数: " + to_string(n) + ", 边数: " + to_string(adjncy.size()/2) +
//        ", 目标分区数: " + to_string(nparts) + ", 平均每分区节点数: " + to_string(n/nparts));
//
//    // **预处理：分析图结构并应用启发式改进**
//    Now("分析图结构特性...");
//
//    // 计算度数分布
//    vector<unsigned> degrees(n);
//    for (unsigned i = 0; i < n; ++i) {
//        degrees[i] = xadj[i + 1] - xadj[i];
//    }
//
//    // 找出高度数节点（可能是造成不平衡的原因）
//    sort(degrees.rbegin(), degrees.rend());
//    unsigned maxDegree = degrees[0];
//    unsigned avgDegree = accumulate(degrees.begin(), degrees.end(), 0u) / n;
//
//    cout << "度数统计：最大度数=" << maxDegree << ", 平均度数=" << avgDegree
//         << ", 度数比值=" << (double)maxDegree / avgDegree << endl;
//
//    // **多策略尝试**
//    vector<idx_t> bestPart;
//    idx_t bestObjval = UINT_MAX;
//    double bestBalance = DBL_MAX;
//
//    for (int strategy = 0; strategy < 4; ++strategy) {
//        Now("尝试策略 " + to_string(strategy + 1) + "/4...");
//
//        idx_t options[METIS_NOPTIONS];
//        METIS_SetDefaultOptions(options);
//
//        switch (strategy) {
//            case 0: // 保守平衡策略 - 使用最稳定的参数
//                options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
//                options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
//                options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
//                options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
//                options[METIS_OPTION_UFACTOR] = 1050;
//                options[METIS_OPTION_NCUTS] = 100;
//                options[METIS_OPTION_NITER] = 100;
//                break;
//
//            case 1: // 经过验证的成功策略
//                options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
//                options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
//                options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
//                options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
//                options[METIS_OPTION_UFACTOR] = 1100;
//                options[METIS_OPTION_NCUTS] = 150;
//                options[METIS_OPTION_NITER] = 150;
//                break;
//
//            case 2: // 重边匹配策略
//                options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
//                options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
//                options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
//                options[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;
//                options[METIS_OPTION_UFACTOR] = 1200;
//                options[METIS_OPTION_NCUTS] = 100;
//                options[METIS_OPTION_NITER] = 100;
//                break;
//
//            case 3: // 宽松策略
//                options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
//                options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
//                options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
//                options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
//                options[METIS_OPTION_UFACTOR] = 1500;
//                options[METIS_OPTION_NCUTS] = 200;
//                options[METIS_OPTION_NITER] = 50;
//                break;
//        }
//
//        options[METIS_OPTION_NUMBERING] = 0;
//        options[METIS_OPTION_SEED] = 42 + strategy * 17; // 不同种子
//
//        idx_t nvtxs = n;
//        idx_t ncon = 1;
//        idx_t numParts = nparts;
//        vector<idx_t> part(nvtxs);
//        idx_t objval;
//
//        // 尝试k-way分区
//        int result = METIS_PartGraphKway(
//                &nvtxs, &ncon, xadj.data(), adjncy.data(),
//                nullptr, nullptr, adjwgt.data(), &numParts,
//                nullptr, nullptr, options, &objval, part.data()
//        );
//
//        if (result == METIS_OK) {
//            // 计算平衡度
//            vector<unsigned> partSizes(numParts, 0);
//            for (unsigned i = 0; i < n; ++i) {
//                partSizes[part[i]]++;
//            }
//
//            auto [minIt, maxIt] = minmax_element(partSizes.begin(), partSizes.end());
//            double balance = (double)*maxIt / *minIt;
//
//            cout << "策略 " << strategy + 1 << " 结果：平衡度=" << balance
//                 << ", 边切割=" << objval << endl;
//
//            // 选择最平衡的结果
//            if (balance < bestBalance) {
//                bestBalance = balance;
//                bestObjval = objval;
//                bestPart = part;
//            }
//        } else {
//            cout << "策略 " << strategy + 1 << " 失败，错误码: " << result << endl;
//        }
//    }
//
//    if (bestPart.empty()) {
//        throw runtime_error("所有METIS策略都失败了");
//    }
//
//    Now("选择最佳策略，平衡度: " + to_string(bestBalance) + ", 边切割数: " + to_string(bestObjval));
//
//    // **如果仍然极度不平衡，使用贪心重分配**
//    if (bestBalance > 100.0) {
//        Now("METIS结果仍然不平衡，启动贪心重分配算法...");
//        GreedyRebalance(bestPart, nparts);
//    }
//
//    // 应用最佳分区结果
//    partitions.clear();
//    borders.clear();
//    partitions.resize(nparts);
//    borders.resize(nparts);
//
//    for (unsigned i = 0; i < n; ++i) {
//        vertexList[i].partition = bestPart[i];
//        partitions[bestPart[i]].push_back(i);
//    }
//
//    // 计算边界
//    RecalculateBoundaries();
//
//    // 最终统计
//    unsigned totalBorders = 0;
//    for (const auto& border : borders) {
//        totalBorders += border.size();
//    }
//
//    cout << "最终分区统计：" << endl;
//    cout << "实际分区数: " << nparts << endl;
//
//    auto sizes = vector<unsigned>(nparts);
//    for (unsigned i = 0; i < nparts; ++i) {
//        sizes[i] = partitions[i].size();
//    }
//    auto [minIt, maxIt] = minmax_element(sizes.begin(), sizes.end());
//    double finalBalance = (double)*maxIt / *minIt;
//    cout << "最终平衡度：" << finalBalance << endl;
//
//    for (unsigned i = 0; i < nparts; ++i) {
//        cout << "分区 " << i << ": " << partitions[i].size() << " 个节点, "
//             << borders[i].size() << " 个边界节点" << endl;
//    }
//    cout << "总边界节点数: " << totalBorders << " ("
//         << (totalBorders * 100.0 / n) << "% of total nodes)" << endl;
//
//    Now("METIS图分区完成");
//}
//
//
//// 新增：贪心重分配函数
///**
// * @brief 增强的贪心重分配函数，专注于减少边界节点和跨分区边
// * @param part 分区分配数组
// * @param nparts 分区数量
// */
//void Graph::GreedyRebalance(vector<idx_t>& part, unsigned nparts) {
//    unsigned targetSize = n / nparts;
//    vector<unsigned> partSizes(nparts, 0);
//
//    // 统计初始分区大小
//    for (unsigned i = 0; i < n; ++i) {
//        partSizes[part[i]]++;
//    }
//
//    cout << "开始平衡优先的重分配..." << endl;
//
//    // **阶段1：严格平衡约束下的优化**
//    for (int iter = 0; iter < 10; ++iter) {
//        bool improved = false;
//
//        // 找出最不平衡的分区对
//        auto [minIt, maxIt] = minmax_element(partSizes.begin(), partSizes.end());
//        unsigned minPart = distance(partSizes.begin(), minIt);
//        unsigned maxPart = distance(partSizes.begin(), maxIt);
//
//        if (*maxIt <= targetSize * 1.3) break; // 已经足够平衡
//
//        // 从过大分区向过小分区移动节点
//        vector<pair<int, unsigned>> candidates; // <收益, 节点ID>
//
//        for (unsigned nodeId = 0; nodeId < n; ++nodeId) {
//            if (part[nodeId] != maxPart) continue;
//
//            // 计算移动到目标分区的收益
//            int benefit = 0;
//            bool hasConnectionToMin = false;
//
//            for (const auto& edge : vertexList[nodeId].edges) {
//                if (part[edge.first] == minPart) {
//                    hasConnectionToMin = true;
//                    benefit += 2; // 内部边奖励
//                } else if (part[edge.first] != maxPart) {
//                    benefit -= 1; // 跨分区边惩罚
//                }
//            }
//
//            // 优先选择有连接或边界节点
//            if (hasConnectionToMin || benefit >= 0) {
//                candidates.push_back({benefit, nodeId});
//            }
//        }
//
//        // 按收益排序并执行移动
//        sort(candidates.rbegin(), candidates.rend());
//
//        unsigned moveCount = min(
//                static_cast<unsigned>(candidates.size()),
//                (*maxIt - *minIt) / 2  // 不超过平衡所需数量
//        );
//
//        for (unsigned i = 0; i < moveCount && i < 50; ++i) {
//            unsigned nodeId = candidates[i].second;
//            part[nodeId] = minPart;
//            partSizes[maxPart]--;
//            partSizes[minPart]++;
//            improved = true;
//        }
//
//        if (improved) {
//            cout << "第 " << iter << " 轮：分区 " << maxPart << " -> "
//                 << minPart << "，移动 " << moveCount << " 个节点" << endl;
//        } else {
//            break;
//        }
//    }
//
//    // **阶段2：渐进式边界优化**（在保持平衡的前提下）
//    cout << "开始边界优化..." << endl;
//
//    for (int iter = 0; iter < 8; ++iter) {
//        bool improved = false;
//
//        // 更严格的平衡约束
//        auto [minIt, maxIt] = minmax_element(partSizes.begin(), partSizes.end());
//        if ((double)*maxIt / *minIt > 2.0) {
//            cout << "平衡度过差，跳过边界优化" << endl;
//            break;
//        }
//
//        for (unsigned nodeId = 0; nodeId < n; ++nodeId) {
//            unsigned currentPart = part[nodeId];
//
//            // 计算各分区的连接强度
//            map<unsigned, unsigned> connections;
//            for (const auto& edge : vertexList[nodeId].edges) {
//                connections[part[edge.first]]++;
//            }
//
//            // 寻找更好的分区
//            for (const auto& conn : connections) {
//                unsigned targetPart = conn.first;
//                if (targetPart == currentPart) continue;
//
//                // 严格平衡约束检查
//                if (partSizes[targetPart] >= targetSize * 1.4 ||
//                    partSizes[currentPart] <= targetSize * 0.6) continue;
//
//                // 边界改善检查
//                unsigned currentCrossEdges = 0;
//                unsigned targetCrossEdges = 0;
//
//                for (const auto& edge : vertexList[nodeId].edges) {
//                    if (part[edge.first] != currentPart) currentCrossEdges++;
//                    if (part[edge.first] != targetPart) targetCrossEdges++;
//                }
//
//                // 只有在显著改善边界时才移动
//                if (targetCrossEdges < currentCrossEdges - 1 &&
//                    conn.second > connections[currentPart]) {
//
//                    part[nodeId] = targetPart;
//                    partSizes[currentPart]--;
//                    partSizes[targetPart]++;
//                    improved = true;
//                    break;
//                }
//            }
//        }
//
//        if (!improved) break;
//    }
//
//    // **阶段3：最终统计和报告**
//    unsigned totalBoundaryNodes = 0;
//    unsigned totalCrossEdges = 0;
//
//    for (unsigned i = 0; i < n; ++i) {
//        bool isBoundary = false;
//        for (const auto& e : vertexList[i].edges) {
//            if (part[e.first] != part[i]) {
//                totalCrossEdges++;
//                isBoundary = true;
//            }
//        }
//        if (isBoundary) totalBoundaryNodes++;
//    }
//
//    totalCrossEdges /= 2;
//
//    // 重新计算平衡度
//    fill(partSizes.begin(), partSizes.end(), 0);
//    for (unsigned i = 0; i < n; ++i) {
//        partSizes[part[i]]++;
//    }
//
//    auto [finalMinIt, finalMaxIt] = minmax_element(partSizes.begin(), partSizes.end());
//    double finalBalance = (double)*finalMaxIt / *finalMinIt;
//
//    cout << "平衡优先重分配完成：" << endl;
//    cout << "  最终平衡度: " << finalBalance << endl;
//    cout << "  边界节点数: " << totalBoundaryNodes << " ("
//         << (100.0 * totalBoundaryNodes / n) << "%)" << endl;
//    cout << "  跨分区边数: " << totalCrossEdges << endl;
//
//    for (unsigned i = 0; i < nparts; ++i) {
//        cout << "  分区 " << i << ": " << partSizes[i] << " 个节点" << endl;
//    }
//}
//
//
///**
// * @brief 将图数据转换为METIS格式
// * @param dataPath 输入图文件路径
// * @param xadj 输出：邻接表索引数组
// * @param adjncy 输出：邻接表数组
// * @param adjwgt 输出：边权重数组
// */
//void Graph::ConvertToMetisFormat(const string &dataPath,
//                                 vector<idx_t> &xadj,
//                                 vector<idx_t> &adjncy,
//                                 vector<idx_t> &adjwgt) {
//    Now("读取并转换图数据为METIS格式...");
//
//    ifstream data(dataPath);
//    if (!data.is_open()) {
//        throw runtime_error("无法打开文件: " + dataPath);
//    }
//
//    string line;
//    unsigned s, t;
//
//    // 读取第一行
//    getline(data, line);
//    istringstream row(line);
//    unsigned temp_c;
//    row >> n >> m >> temp_c;
//
//    Now("图信息：节点数=" + to_string(n) + ", 边数=" + to_string(m) + ", 维度数=" + to_string(temp_c));
//
//    // 构建邻接表
//    vector<set<unsigned>> tempAdj(n); // 使用set自动去重
//    unsigned actualEdges = 0;
//
//    // 读取边数据
//    while (getline(data, line)) {
//        if (line.empty()) continue;
//
//        row.clear();
//        row.str(line);
//
//        if (!(row >> s >> t)) continue;
//
//        // 验证节点ID
//        if (s >= n || t >= n) {
//            cerr << "警告：发现无效的节点ID: " << s << " 或 " << t << endl;
//            continue;
//        }
//
//        // 跳过权重数据
//        for (unsigned i = 0; i < temp_c; ++i) {
//            unsigned dummy;
//            if (!(row >> dummy)) break;
//        }
//
//        // 添加边（set自动处理重复）
//        if (s != t) {
//            bool isNewEdge = (tempAdj[s].find(t) == tempAdj[s].end());
//            tempAdj[s].insert(t);
//            tempAdj[t].insert(s);
//            if (isNewEdge) actualEdges++;
//        }
//    }
//    data.close();
//
//    cout << "实际读取的唯一边数: " << actualEdges << endl;
//
//    // **关键诊断：连通组件分析**
//    vector<bool> visited(n, false);
//    vector<unsigned> componentSizes;
//
//    for (unsigned i = 0; i < n; ++i) {
//        if (!visited[i]) {
//            // BFS找连通组件
//            queue<unsigned> q;
//            q.push(i);
//            visited[i] = true;
//            unsigned componentSize = 0;
//
//            while (!q.empty()) {
//                unsigned v = q.front();
//                q.pop();
//                componentSize++;
//
//                for (unsigned neighbor : tempAdj[v]) {
//                    if (!visited[neighbor]) {
//                        visited[neighbor] = true;
//                        q.push(neighbor);
//                    }
//                }
//            }
//            componentSizes.push_back(componentSize);
//        }
//    }
//
//    // 输出连通性分析
//    sort(componentSizes.rbegin(), componentSizes.rend());
//    cout << "连通组件分析：" << endl;
//    cout << "  连通组件数量: " << componentSizes.size() << endl;
//    for (size_t i = 0; i < min(componentSizes.size(), size_t(10)); ++i) {
//        cout << "  组件 " << i << ": " << componentSizes[i] << " 个节点 ("
//             << (100.0 * componentSizes[i] / n) << "%)" << endl;
//    }
//
//    // **问题检测**
//    if (componentSizes.size() > 1 && componentSizes[0] > 0.9 * n) {
//        cout << "警告：检测到严重不平衡的连通结构！" << endl;
//        cout << "主连通组件占 " << (100.0 * componentSizes[0] / n) << "% 的节点" << endl;
//    }
//
//    // 转换为CSR格式
//    xadj.resize(n + 1);
//    adjncy.clear();
//    adjwgt.clear();
//
//    xadj[0] = 0;
//    for (unsigned i = 0; i < n; ++i) {
//        for (unsigned neighbor : tempAdj[i]) {
//            adjncy.push_back(neighbor);
//            adjwgt.push_back(1); // 统一权重
//        }
//        xadj[i + 1] = adjncy.size();
//    }
//
//    cout << "METIS CSR格式：" << endl;
//    cout << "  邻接表大小: " << adjncy.size() << endl;
//    cout << "  平均节点度数: " << (n > 0 ? (double)adjncy.size() / n : 0) << endl;
//
//    Now("图转换完成，邻接表大小: " + to_string(adjncy.size()));
//}
//
///**
// * @brief 优化分区边界，通过节点迁移减少边界节点数量
// */
//void Graph::OptimizePartitionBoundaries() {
//    bool improved = true;
//    int iterations = 0;
//    const int maxIterations = 5;
//
//    while (improved && iterations < maxIterations) {
//        improved = false;
//        iterations++;
//
//        vector<unsigned> candidateNodes;
//        for (unsigned i = 0; i < n; ++i) {
//            if (vertexList[i].boundary) {
//                candidateNodes.push_back(i);
//            }
//        }
//
//        for (unsigned nodeId : candidateNodes) {
//            if (!vertexList[nodeId].boundary) continue; // 可能已被优化
//
//            // 统计邻居分区
//            map<unsigned, unsigned> neighborPartitions;
//            for (const auto& edge : vertexList[nodeId].edges) {
//                unsigned neighborPart = vertexList[edge.first].partition;
//                neighborPartitions[neighborPart]++;
//            }
//
//            // 找到最佳迁移目标分区
//            unsigned currentPart = vertexList[nodeId].partition;
//            unsigned bestPart = currentPart;
//            unsigned maxConnections = neighborPartitions[currentPart];
//
//            for (const auto& np : neighborPartitions) {
//                if (np.first != currentPart && np.second > maxConnections) {
//                    maxConnections = np.second;
//                    bestPart = np.first;
//                }
//            }
//
//            // 如果迁移有利，则执行
//            if (bestPart != currentPart &&
//                maxConnections > neighborPartitions[currentPart] + 1) {
//
//                // 检查分区大小平衡
//                if (partitions[bestPart].size() < partitions[currentPart].size()) {
//                    // 执行节点迁移
//                    vertexList[nodeId].partition = bestPart;
//
//                    // 更新分区列表
//                    auto it = find(partitions[currentPart].begin(),
//                                   partitions[currentPart].end(), nodeId);
//                    partitions[currentPart].erase(it);
//                    partitions[bestPart].push_back(nodeId);
//
//                    improved = true;
//                }
//            }
//        }
//
//        if (improved) {
//            RecalculateBoundaries();
//        }
//    }
//}
//
///**
// * @brief 重新计算所有边界节点
// */
//void Graph::RecalculateBoundaries() {
//    // 清除现有边界信息
//    for (auto& border : borders) border.clear();
//
//    for (unsigned i = 0; i < n; ++i) {
//        vertexList[i].boundary = false;
//
//        for (const auto& e : vertexList[i].edges) {
//            if (vertexList[i].partition != vertexList[e.first].partition) {
//                if (!vertexList[i].boundary) {
//                    borders[vertexList[i].partition].insert(i);
//                    vertexList[i].boundary = true;
//                }
//                break;
//            }
//        }
//    }
//}



#endif //POP_GRAPH_HPP
