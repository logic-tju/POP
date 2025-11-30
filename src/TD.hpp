/**
 ********************************************
 * @file    :TD.hpp
 * @author  :XXY
 * @brief   :Tree Decomposition
 * @date    :2025/11/30
 ********************************************
 */

#ifndef POP_TD_HPP
#define POP_TD_HPP

#include "Graph.hpp"
#include "CH.hpp"
#include "common.hpp"

/**
 * @brief Tree node structure for tree decomposition
 */
struct TreeNode {
    unsigned nid{};        // Node ID
    unsigned order{};      // Order of the node in the elimination process
    unsigned layer{};      // Layer of the node in the tree

    TreeNode* parent{nullptr};                 // Parent node
    vector<TreeNode*> children{};              // Children nodes

    map<unsigned, CH*> neighbors;              // Neighbors of the node
    map<unsigned, CH*> dis2Borders;            // Distance to border nodes

    void Print(ofstream& out) {
        out << nid << ' ' << order << ' ' << layer << ' '
            << (parent == nullptr ? UINT_MAX : parent->nid) << ' '
            << children.size() << ' ';
        for (const auto& c : children) {
            out << c->nid << ' ';
        }
        out << dis2Borders.size() << ' ';
        for (const auto& d : dis2Borders) {
            out << d.first << ' ';
            d.second->Print(out);
        }
        out << endl;
    }

    void Read(istringstream& in, map<unsigned, TreeNode*>& treeNodes) {
        in >> nid >> order >> layer;
        unsigned parentID;
        in >> parentID;
        if (parentID != UINT_MAX) {
            parent = treeNodes[parentID];
        }
        unsigned childrenSize;
        in >> childrenSize;
        for (unsigned i = 0; i < childrenSize; ++i) {
            unsigned childID;
            in >> childID;
            children.push_back(treeNodes[childID]);
        }
        unsigned dis2BordersSize;
        in >> dis2BordersSize;
        for (unsigned i = 0; i < dis2BordersSize; ++i) {
            unsigned border;
            in >> border;
            dis2Borders[border] = new CH(in);
        }
    }
};

/**
 * @brief Tree decomposition class
 */
class TD {
public:
    // there may be many roots if partition is not connected
    vector<TreeNode*> root{};        // Roots of the tree
    map<unsigned, TreeNode*> treeNodes;  // Map to store tree nodes by ID
    set<unsigned> borders;  // Set to store border nodes

    // Constructor
    TD();

    // Constructor for Tree Decomposition based on Graph and partition number
    TD(Graph& G, unsigned p);

    // Function to set label for nodes in the tree
    void SetLabel();

    // Function to print the tree decomposition to a file
    void Print(const string& path);

    // Function to read the tree decomposition from a file
    void Read(const string& path);
};

TD::TD() = default;

TD::TD(Graph& G, unsigned p) {
    borders = G.borders[p];
    set<OrderedNode> degreeOrder;

    for (auto v : G.partitions[p]) {
        if (!G.vertexList[v].boundary) {
            degreeOrder.insert({static_cast<unsigned>(G.vertexList[v].edges.size()), v});
        }
    }

    unsigned count = 0;
    while (!degreeOrder.empty()) {
        unsigned minDegreeNode = degreeOrder.begin()->id;
        degreeOrder.erase(degreeOrder.begin());

        auto* newNode = new TreeNode;
        newNode->nid = minDegreeNode;
        newNode->order = count++;

        for (const auto& e : G.vertexList[minDegreeNode].edges) {
            newNode->neighbors[e.first] = e.second;
        }
        treeNodes[minDegreeNode] = newNode;

        G.vertexList[minDegreeNode].valid = count;
        G.VertexContract(minDegreeNode, degreeOrder);
        // 每处理 100 个节点输出一次进度，减少输出频率
//        if (count % 100 == 0) {
//            std::cout << "MDE Step: " << count << "/" << G.partitions[p].size() << " Nodes processed." << std::endl;
//        }
    }

    for (auto t : treeNodes) {
        unsigned parent = UINT_MAX;
        auto it = t.second->neighbors.begin();
        while (it != t.second->neighbors.end() && G.vertexList[it->first].boundary) {
            ++it;
        }

        if (it == t.second->neighbors.end()) {
            root.emplace_back(t.second);
            continue;
        } else {
            parent = it->first;
            for (; it != t.second->neighbors.end(); ++it) {
                if (!G.vertexList[it->first].boundary &&
                    treeNodes[it->first]->order < treeNodes[parent]->order) {
                    parent = it->first;
                }
            }
        }

        t.second->parent = treeNodes[parent];
        treeNodes[parent]->children.push_back(t.second);
    }
//    Now("Tree Decomposition Constructed");
}

void TD::SetLabel() {
    unsigned count = 0;
    for (auto r : root) {
        queue<TreeNode*> nodeQueue;
        nodeQueue.push(r);

        while (!nodeQueue.empty()) {
            auto node = nodeQueue.front();
            nodeQueue.pop();

            for (auto child : node->children) {
                nodeQueue.push(child);
            }

            for (auto b : r->neighbors) {
                vector<Point> newPoints;
                map<unsigned, CH*> innerNeighbor, borderNeighbor;

                for (auto n : node->neighbors) {
                    if (borders.find(n.first) != borders.end()) {
                        borderNeighbor.insert(n);
                    } else {
                        innerNeighbor.insert(n);
                    }
                }

                for (auto neighbor : innerNeighbor) {
                    vector<vector<double>> l1 = neighbor.second->points;
                    vector<vector<double>> l2 = treeNodes[neighbor.first]->dis2Borders[b.first]->points;
                    LSCatenation(l1, l2, newPoints);
                }

                if (borderNeighbor.find(b.first) != borderNeighbor.end()) {
                    for (const auto& p : borderNeighbor[b.first]->points) {
                        Point point;
                        std::fill(point.begin(), point.end(), 0.0); // 初始化为 0
                        size_t i = 0;
                        for (auto val : p) {
                            if (i < D) {
                                point[i++] = static_cast<double>(val); // 将 unsigned 转换为 double
                            } else {
                                break; // 确保不超过 D 的维度
                            }
                        }
                        newPoints.emplace_back(point);
                    }
                }

                node->dis2Borders[b.first] = new CH(newPoints);
            }

            for (auto& n : node->neighbors) {
                delete n.second;
            }
            if (node != r) {
                node->neighbors.clear();
            }
        }
    }
//    Now("label set");
}

void TD::Print(const string &path){
    ofstream out(path);
    out << treeNodes.size() << ' ';
    for(auto tn : treeNodes){
        out << tn.first << ' ';
    }
    out << endl;
    for(auto tn : treeNodes){
        tn.second->Print(out);
    }
    for(auto r : root){
        out << r->nid << " ";
    }
    out << endl;
    for(auto b : borders){
        out << b << " ";
    }
    out.close();
}

void TD::Read(const string &path){
    ifstream in(path);
    string line;
    getline(in, line);
    istringstream ss(line);
    unsigned n;
    ss >> n;
    for(unsigned i = 0; i < n; i++){
        unsigned nid;
        ss >> nid;
        treeNodes[nid] = new TreeNode();
    }

    // 为每行创建新的istringstream对象
    for(auto & treeNode : treeNodes){
        getline(in, line);
        istringstream ss2(line);  // 每行使用新的istringstream
        treeNode.second->Read(ss2, treeNodes);
    }

    getline(in, line);
    istringstream ss3(line);  // 新实例
    while(ss3 >> n){
        root.emplace_back(treeNodes[n]);
    }

    getline(in, line);
    istringstream ss4(line);  // 新实例
    while(ss4 >> n){
        borders.insert(n);
    }
    in.close();
}


#endif // POP_TD_HPP
