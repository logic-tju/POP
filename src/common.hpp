/**
 ********************************************
 * @file    :common.hpp
 * @author  :XXY
 * @brief   :Basic Structs And Headers
 * @date    :2025/11/30
 ********************************************
 */

#ifndef POP_COMMON_HPP
#define POP_COMMON_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <array>
#include <stack>
#include <queue>
#include <map>
#include <set>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <random>
#include <numeric>
#include <memory>
#include <omp.h>
#include <mutex>
#include <string>
#include <chrono>
#include <limits>
#include <tuple>
#include <thread>
#include <stdexcept>
#include <unordered_set>
#include <fstream>
#include <climits>
//#include <metis.h>


// 引入 Qhull 的 C++ 接口头文件
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullError.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullVertexSet.h>

using namespace std;

// 定义全局维度 D，修改此值即可控制代码运行的维度
const int D = 3;

// 定义 Qhull 使用的数据类型
typedef std::array<double, D> Point;  // 使用 std::array<double, D> 表示点坐标
typedef pair<vector<double>, double> Path;

// 定义无穷点集合（仅正向）
inline std::vector<Point> INFINITYPOINTS;

default_random_engine DRE;
uniform_int_distribution<unsigned> DISTRIBUTION(1, 100000);

/**
 * @brief sort by degree and find by vertex id
 */
struct OrderedNode{
    unsigned degree;
    unsigned id;

    bool operator< (const OrderedNode& n) const{
        return degree < n.degree || (degree == n.degree && id < n.id);
    }
//    bool operator== (const OrderedNode& n) const{
//        return id == n.id;
//    }
};

/**
 * @brief 按维度初始化无限点
 */
inline void InfPointsInit() {
    INFINITYPOINTS.clear();
    for (unsigned i = 0; i < D; ++i) {
        Point temp;
        temp.fill(0.0);
        temp[i] = 1e9;  // 使用一个非常大的值表示无穷点
        INFINITYPOINTS.push_back(temp);
    }
}

/**
 * @brief 计算权重和成本向量的点积
 * @param weight 权重向量
 * @param cost 成本向量
 * @return 点积结果
 */
inline double DotProduct(const std::vector<double> &weight, const std::vector<double> &cost) {
    double result = 0.0;
    for (size_t i = 0; i < weight.size(); ++i) {
        result += weight[i] * cost[i];
    }
    return result;
}

/**
 * @brief 向量加法
 * @param v1 第一个向量
 * @param v2 第二个向量
 * @return 两个向量相加的结果
 */
inline vector<double> VectorAdd(const vector<double>& v1, const vector<double>& v2) {
    vector<double> result(D);
    for (size_t i = 0; i < D; ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

/**
 * @brief 列表连接函数，用于合并两个邻接点列表
 * @param l1 邻居1的列表
 * @param l2 邻居2的列表
 * @param result 合并后的结果列表
 */
inline void LSCatenation(const vector<vector<double>>& l1, const vector<vector<double>>& l2, vector<Point>& result) {
    for (const auto& v1 : l1) {
        for (const auto& v2 : l2) {
            vector<double> sum = VectorAdd(v1, v2);

            // 转换为 Point 类型
            Point point;
            for (size_t i = 0; i < D; ++i) {
                point[i] = sum[i];
            }
            result.emplace_back(point);
        }
    }
}

/**
 * @brief 数据集合并函数
 * @param input1 第一个数据集路径
 * @param input2 第二个数据集路径
 * @param output 合并后数据集路径
 */
[[maybe_unused]] inline void DatasetMerge(const string& input1, const string& input2, const string& output) {
    ifstream in1(input1);
    ifstream in2(input2);
    ofstream out(output);
    string line, temp;
    unsigned s, t, c;

    queue<unsigned> cost;
    while (getline(in1, line)) {
        istringstream row(line);
        row >> temp;
        if (temp == "a") {
            row >> s >> t >> c;
            cost.push(c);
        }
    }

    while (getline(in2, line)) {
        istringstream row(line);
        row >> temp;
        if (temp == "p") {
            row >> temp >> s >> t;
            out << s << ' ' << t << ' ' << 2 << endl;
        } else if (temp == "a") {
            row >> s >> t >> c;
            if (!cost.empty()) {
                out << (s - 1) << ' ' << (t - 1) << ' ' << c << ' ' << cost.front() << endl;
                cost.pop();
            } else {
                out << (s - 1) << ' ' << (t - 1) << ' ' << c << " 0" << endl; // 默认成本为0
            }
        }
    }

    in1.close();
    in2.close();
    out.close();
}

/**
 * @brief 显示进度的辅助函数
 * @param output 输出的字符串
 * @param up 是否将光标上移
 */
inline void Now(const string& output, bool up = false) {
    time_t now = time(nullptr);
    cout << "\033[K\033[32;40m" << output << " At " << ctime(&now) << (up ? "\033[1A" : "") << flush;
}

/**
 * @brief 获取当前时间字符串
 * @return 格式化的时间字符串
 */
inline std::string getCurrentTimeString() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

/**
 * @brief 简化的PL索引内存统计函数
 * @param g 图对象引用
 * @param queryPath 查询路径
 */
template<typename GraphType>
void calculatePLIndexSizeToFile(const GraphType& g, const std::string& queryPath) {
    std::string outputFile = queryPath + "/pl_memory_statistics.txt";
    std::ofstream outFile(outputFile);

    if (!outFile.is_open()) {
        std::cerr << "Error: Cannot open file " << outputFile << " for writing." << std::endl;
        return;
    }

    size_t total_pl_size = 0;
    int pl_count = 0;

    // 遍历所有顶点和边，收集统计信息
    for (unsigned i = 0; i < g.vertexList.size(); ++i) {
        for (const auto& edge_pair : g.vertexList[i].edges) {
            if (edge_pair.second->pl != nullptr) {
                total_pl_size += edge_pair.second->pl->getIndexSize();
                pl_count++;
            }
        }
    }

    // 输出简单的统计信息
    outFile << "Total PL structures: " << pl_count << std::endl;
    outFile << "Total PL memory: " << total_pl_size << " bytes ("
            << std::fixed << std::setprecision(2) << total_pl_size / (1024.0 * 1024.0) << " MB)" << std::endl;

    outFile.close();

    // 控制台输出
    std::cout << "PL structures: " << pl_count << ", Total memory: "
              << std::fixed << std::setprecision(2) << total_pl_size / (1024.0 * 1024.0) << " MB" << std::endl;
}

#endif // POP_COMMON_HPP
