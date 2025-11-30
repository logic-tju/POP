/**
 ********************************************
 * @file    :main.cpp
 * @author  :XXY
 * @brief   :Main entry point for POP
 * @date    :2025/11/30
 ********************************************
 */

#include "TD.hpp"

using namespace std;

[[maybe_unused]]void Read(const string &dataPath,
                          const string &borderPath,
                          const string &indexPath,
                          const string &queryPath,
                          const vector<double> &weight,
                          const unsigned &costNum,
                          const unsigned &randomNum,
                          const unsigned &partitionNum);
void Write(const string &dataPath,
           const string &clusterPath,
           const string &dotPath,
           const string &indexPath,
           const string &borderPath,
           const unsigned &costNum,
           const unsigned &randomNum,
           const unsigned &partitionNum,
           const vector<double> &weight,
           const string &queryPath,
           const string &nodePerPartition);

vector<double> DijkQ(const vector<pair<unsigned, unsigned>> &query,
                     const vector<double> &weight,
                     const string &queryPath,
                     const Graph &g);

vector<double> ASQ(const vector<pair<unsigned, unsigned>> &query,
                   const vector<double> &weight,
                   const string &queryPath,
                   Graph &g,
                   const vector<TD *> &trees);
vector<double> ASQ_pl(const vector<pair<unsigned, unsigned>> &query,
                      const vector<double> &weight,
                      const string &queryPath,
                      Graph &g,
                      const vector<TD *> &trees);
vector<pair<unsigned, unsigned>> queryGenerator(Graph &g, unsigned num);

int main(int argc, char *argv[]) {
    // usage: POP DATASET NODEPERPARTITION PATH2ROOTDIR DATASET_ COSTNUM RANDOMNUM WEIGHT
    assert(argc >= 6);
    string DATASET = argv[1];
    string NODEPERPARTITION = argv[2];
    string PATH2ROOTDIR = argv[3];
    string DATASET_ = argv[4];

    auto command = "cd " + PATH2ROOTDIR + "/natural_cut && bash run.sh " + DATASET + " " + NODEPERPARTITION;
    auto dataPath = PATH2ROOTDIR + "/datasets/" + DATASET_;
    auto clusterPath = PATH2ROOTDIR + "/natural_cut/data/" + DATASET + "/node_clusters.txt";
    auto dotPath = PATH2ROOTDIR + "/output/" + DATASET + "/border_graph.dot";
    auto indexPath = PATH2ROOTDIR + "/index_storage/" + DATASET + "/td";
    auto borderPath = PATH2ROOTDIR + "/index_storage/" + DATASET + "/border_graph";
    auto queryPath = PATH2ROOTDIR + "/output/" + DATASET + "/";
    auto costNum = stoi(argv[5]);
    auto randomNum = stoi(argv[6]);

    cout << "Data path: " << dataPath << endl;
    cout << "Cluster path: " << clusterPath << endl;
    cout << "Dot path: " << dotPath << endl;
    cout << "Index path: " << indexPath << endl;
    cout << "Border path: " << borderPath << endl;
    cout << "Query path: " << queryPath << endl;
    cout << "costNum: " << costNum << endl;
    cout << "randomNum: " << randomNum << endl;
    cout << "dimension: " << D << endl;

    vector<double> weight;
    for (int i = 7; i != argc; ++i) {
        weight.emplace_back(stod(argv[i]));
    }
    auto partitionNum = 0;

    int nc = system(command.c_str());
    assert(nc == 0);

    ifstream in(clusterPath);
    in >> partitionNum;
    in.close();

    InfPointsInit();
    Write(dataPath, clusterPath, dotPath, indexPath, borderPath, costNum, randomNum, partitionNum, weight, queryPath, NODEPERPARTITION);
    Read(dataPath, borderPath, indexPath, queryPath, weight, costNum, randomNum, partitionNum);

    cout << "FINISHED" << endl;
    return 0;
}

void Write(const string &dataPath,
           const string &clusterPath,
           const string &dotPath,
           const string &indexPath,
           const string &borderPath,
           const unsigned &costNum,
           const unsigned &randomNum,
           const unsigned &partitionNum,
           const vector<double> &weight,
           const string &queryPath,
           const string &nodePerPartition) {
    Graph g(dataPath, costNum, randomNum);
    g.SetParNum(clusterPath);

    Now("begin");

    int progress_interval = max(1, static_cast<int>(partitionNum * 0.1));
    int completed = 0;

#pragma omp parallel for schedule(dynamic) shared(completed, progress_interval, cout)
    for (unsigned i = 0; i < partitionNum; ++i) {
        TD *newTree = new TD(g, i);
        newTree->SetLabel();

        newTree->Print(indexPath + to_string(i));

        delete newTree;

#pragma omp atomic
        completed++;

        if (completed % progress_interval == 0) {
#pragma omp critical
            {
                double progress = static_cast<double>(completed) / partitionNum * 100.0;
                cout << "Progress: " << fixed << setprecision(1) << progress
                     << "% (" << completed << "/" << partitionNum << " partitions processed)" << endl;
            }
        }
    }

#pragma omp barrier
    Now("All partitions have been processed.");
    g.PrintBorderGraph(borderPath);
    Now("Print Border Graph.");
}

[[maybe_unused]] void Read(const std::string &dataPath,
                           const std::string &borderPath,
                           const std::string &indexPath,
                           const std::string &queryPath,
                           const std::vector<double> &weight,
                           const unsigned &costNum,
                           const unsigned &randomNum,
                           const unsigned &partitionNum) {
    Graph g_(dataPath, costNum, randomNum);
    Graph g(borderPath);

    Now("Building PL structures...");
    int totalPLCount = 0;
    for (unsigned i = 0; i < partitionNum; ++i) {
        totalPLCount += g.ComputePointLocation(i);
    }
    std::cout << "Total PL structures built: " << totalPLCount << std::endl;
    Now("PL structures completed");

    Now("Calculating PL index memory usage...");
    calculatePLIndexSizeToFile(g, queryPath);
    Now("PL index memory calculation completed");

    g.PrepareAStar();

    std::vector<TD *> trees(partitionNum, nullptr);

    Now("Loading tree structures in parallel...");

#pragma omp parallel for schedule(dynamic) num_threads(std::min(partitionNum, 8u))
    for (unsigned i = 0; i < partitionNum; ++i) {
        TD *newTree = new TD();
        try {
            newTree->Read(indexPath + std::to_string(i));
            trees[i] = newTree;
        } catch (const std::exception& e) {
#pragma omp critical
            {
                std::cerr << "Error reading tree " << i << ": " << e.what() << std::endl;
            }
            delete newTree;
            trees[i] = nullptr;
        }
    }

    Now("Tree structures loaded");

    bool all_loaded = true;
    for (unsigned i = 0; i < partitionNum; ++i) {
        if (trees[i] == nullptr) {
            std::cerr << "Failed to load tree " << i << std::endl;
            all_loaded = false;
        }
    }

    if (!all_loaded) {
        for (auto tree : trees) {
            if (tree != nullptr) {
                delete tree;
            }
        }
        throw std::runtime_error("Failed to load all trees");
    }

    auto query = queryGenerator(g, 100);

    auto dij = DijkQ(query, weight, queryPath, g_);
    auto astar = ASQ(query, weight, queryPath, g, trees);
    auto astar_pl = ASQ_pl(query, weight, queryPath, g, trees);

    std::cout << "Dijkstra: " << std::accumulate(dij.begin(), dij.end(), 0.0) / (double) dij.size() << std::endl;
    std::cout << "A*: " << std::accumulate(astar.begin(), astar.end(), 0.0) / (double) astar.size() << std::endl;
    std::cout << "A*_pl: " << std::accumulate(astar_pl.begin(), astar_pl.end(), 0.0) / (double) astar_pl.size() << std::endl;

    for (unsigned i = 0; i != dij.size(); ++i) {
        astar[i] = dij[i] / astar[i];
        astar_pl[i] = dij[i] / astar_pl[i];
    }
    auto sp1 = std::accumulate(astar.begin(), astar.end(), 0.0) / (double) astar.size();
    auto sp2 = std::accumulate(astar_pl.begin(), astar_pl.end(), 0.0) / (double) astar_pl.size();
    std::cout << "A* is " << sp1 << " times faster than Dijkstra." << std::endl;
    std::cout << "A*_pl is " << sp2 << " times faster than Dijkstra." << std::endl;

    for (auto tree : trees) {
        if (tree != nullptr) {
            delete tree;
        }
    }
}

vector<pair<unsigned, unsigned>> queryGenerator(Graph &g, unsigned num) {
    vector<pair<unsigned, unsigned>> query;
    uniform_int_distribution<unsigned> Q(1, g.n);
    while (query.size() != num) {
        auto s = Q(DRE), t = Q(DRE);
        if ((!g.vertexList[s].boundary) && (!g.vertexList[t].boundary)
            && (g.vertexList[s].partition != g.vertexList[t].partition)) {
            query.emplace_back(s, t);
        }
    }
    return query;
}

vector<double> DijkQ(const vector<pair<unsigned, unsigned>> &query,
                     const vector<double> &weight,
                     const string &queryPath,
                     const Graph &g) {
    cout << "Dijkstra is running..." << endl;

    struct QueryResult {
        double time_ms;
        double dist;
        unsigned int numNodes;
    };

    vector<QueryResult> results(query.size());

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < query.size(); ++i) {
        const auto& p = query[i];

        auto start = std::chrono::steady_clock::now();
        auto [cost, dist, numNodes] = g.Dijkstra(p.first, p.second, weight);
        auto end = std::chrono::steady_clock::now();
        auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

        results[i] = {time_ms, dist, numNodes};
    }

    ofstream out(queryPath + "dijkstra_query");
    vector<double> dij;
    dij.reserve(query.size());

    for (size_t i = 0; i < query.size(); ++i) {
        const auto& p = query[i];
        const auto& result = results[i];

        dij.push_back(result.time_ms);
        out << '<' << p.first << ", " << p.second << ">: " << setprecision(10)
            << result.dist << " ," << result.time_ms << " ," << result.numNodes << endl;
    }

    out.close();
    return dij;
}

vector<double> ASQ(const vector<pair<unsigned, unsigned>> &query,
                   const vector<double> &weight,
                   const string &queryPath,
                   Graph &g,
                   const vector<TD *> &trees) {
    cout << "A* is running..." << endl;
    vector<double> as(query.size());
    ofstream out(queryPath + "astar_query");

    for (size_t i = 0; i < query.size(); ++i) {
        const auto& p = query[i];

        g.vertexList[p.second].dis2Par.clear();
        g.vertexList[p.second].dis2Par.resize(g.partitions.size());
        g.vertexList[p.second].dis2Par[g.vertexList[p.second].partition] = vector<unsigned>(D, 0);

        auto db1 = trees[g.vertexList[p.first].partition]->treeNodes[p.first]->dis2Borders;
        auto db2 = trees[g.vertexList[p.second].partition]->treeNodes[p.second]->dis2Borders;

        auto start = std::chrono::steady_clock::now();
        auto [naiveMs, plMs] = g.pairedMs(p.first, p.second, weight, db1, db2);
        auto path = g.AStarTest(p.first, p.second, weight, db1, db2, false);
        auto end = std::chrono::steady_clock::now();
        auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

        as[i] = naiveMs;
        out << '<' << p.first << ", " << p.second << ">: " << setprecision(10)
            << path.second << " ," << naiveMs << endl;
    }

    out.close();
    return as;
}

vector<double> ASQ_pl(const vector<pair<unsigned, unsigned>> &query,
                      const vector<double> &weight,
                      const string &queryPath,
                      Graph &g,
                      const vector<TD *> &trees) {
    cout << "A*_pl is running..." << endl;
    vector<double> as_pl(query.size());
    ofstream out(queryPath + "astar_pl_query");

    for (size_t i = 0; i < query.size(); ++i) {
        const auto& p = query[i];

        g.vertexList[p.second].dis2Par.clear();
        g.vertexList[p.second].dis2Par.resize(g.partitions.size());
        g.vertexList[p.second].dis2Par[g.vertexList[p.second].partition] = vector<unsigned>(D, 0);

        auto db1 = trees[g.vertexList[p.first].partition]->treeNodes[p.first]->dis2Borders;
        auto db2 = trees[g.vertexList[p.second].partition]->treeNodes[p.second]->dis2Borders;

        auto start = std::chrono::steady_clock::now();
        auto [naiveMs, plMs] = g.pairedMs(p.first, p.second, weight, db1, db2);
        auto path = g.AStarTest(p.first, p.second, weight, db1, db2, true);
        auto end = std::chrono::steady_clock::now();
        auto time_ms = std::chrono::duration<double, std::milli>(end - start).count();

        as_pl[i] = plMs;
        out << '<' << p.first << ", " << p.second << ">: " << setprecision(10)
            << path.second << " ," << plMs << endl;
    }

    out.close();
    return as_pl;
}