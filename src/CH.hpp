/**
 ********************************************
 * @file    :CH.hpp
 * @author  :XXY
 * @brief   :Convex Hull & Linear Skyline
 * @date    :2025/11/30
 ********************************************
 */

#ifndef POP_CH_HPP
#define POP_CH_HPP

#include "common.hpp"
#include "PL.hpp"

// 使用 Qhull C++ 命名空间
using namespace orgQhull;
using namespace std;

/**
 * @brief Convex Hull Class using Qhull
 */
class CH {
    friend class Graph;
    friend class TD;

public:
    // 凸包顶点集合，使用 double
    std::vector<std::vector<double>> points;
    // 扩展点集合
    vector<pair<vector<double>, double>> exPoints;
    // 使用 std::unique_ptr 管理 PL 对象
    std::unique_ptr<PL> pl{nullptr};  // 点定位结构
    // 标记是否已初始化 Qhull
    bool initialized{false};

public:
    // 构造函数
    CH();
    explicit CH(std::istringstream& in);
    explicit CH(const std::vector<unsigned>& cost);
    explicit CH(const std::vector<std::vector<unsigned>>& inputPoints);
    explicit CH(const std::vector<Point>& qhullPoints);

    // 删除拷贝构造函数和拷贝赋值运算符
    CH(const CH&) = delete;
    CH& operator=(const CH&) = delete;

    // 移动构造函数和移动赋值运算符
    CH(CH&& other) noexcept;
    CH& operator=(CH&& other) noexcept;

    // 析构函数
    ~CH();

    // 成员函数
    void Insert(std::vector<Point>& newPoints);
    void ComputeCHPoints();
    Path NaiveQuery(const std::vector<double>& weight);
    void ComputeExPoints();
    void Print(std::ofstream& out);
};

// 构造函数实现
inline CH::CH() {
    ComputeCHPoints();
}

// 从输入流构造
inline CH::CH(std::istringstream& in) {
    unsigned pointsSize;
    in >> pointsSize;
    points.reserve(pointsSize);
    for (unsigned j = 0; j < pointsSize; ++j) {
        std::vector<double> point(D);
        for (unsigned k = 0; k < D; ++k) {
            in >> point[k];
        }
        points.push_back(point);
    }
}

// 从成本向量构造
inline CH::CH(const std::vector<unsigned>& cost) {
    points.reserve(1);  // 仅存储传入的点
    std::vector<double> doubleCost(cost.begin(), cost.end());
    points.emplace_back(doubleCost);
}

// 从点集合构造
inline CH::CH(const std::vector<std::vector<unsigned>>& inputPoints) {
    points.reserve(inputPoints.size());

    // 逐点插入，使用 emplace_back 以提高效率
    for (const auto& point : inputPoints) {
        std::vector<double> doublePoint(point.begin(), point.end());
        points.emplace_back(doublePoint);
    }
}

// 从 Qhull 点集合构造
inline CH::CH(const std::vector<Point>& qhullPoints) {
    points.reserve(qhullPoints.size());

    // 插入点集
    for (const auto& p : qhullPoints) {
        std::vector<double> doublePoint(p.begin(), p.end());
        points.emplace_back(std::move(doublePoint));
    }
    ComputeCHPoints();
}

// 移动构造函数
inline CH::CH(CH&& other) noexcept
        : pl(std::move(other.pl)),
          initialized(other.initialized),
          points(std::move(other.points)),
          exPoints(std::move(other.exPoints))
{
    other.initialized = false;
}

// 移动赋值运算符
inline CH& CH::operator=(CH&& other) noexcept {
    if (this != &other) {
// 释放当前资源（Qhull C++ 接口的资源由 Qhull 类自动管理，无需手动释放）
        pl = std::move(other.pl);
        initialized = other.initialized;
        points = std::move(other.points);
        exPoints = std::move(other.exPoints);

        other.initialized = false;
    }
    return *this;
}

// 析构函数，Qhull C++ 接口自动管理资源
inline CH::~CH() {
    // std::unique_ptr 会自动释放 pl
    // Qhull C++ 接口无需手动释放
}

/**
 * @brief 插入新的点集并重新计算凸包
 */
inline void CH::Insert(std::vector<Point>& newPoints) {
    // 批量插入点集
    points.reserve(points.size() + newPoints.size());

    for (const auto& p : newPoints) {
        std::vector<double> doublePoint(p.begin(), p.end());
        points.emplace_back(std::move(doublePoint));
    }

    // 重新计算凸包
    ComputeCHPoints();
}

/**
 * @brief 使用 Qhull 计算凸包点集
 */
inline void CH::ComputeCHPoints() {
    try {
        // 将 points (std::vector<std::vector<double>>) 转换为 allPoints (std::vector<Point>)
        std::vector<Point> allPoints;
        allPoints.reserve(points.size() + INFINITYPOINTS.size());

        for (const auto& vec : points) {
            if (vec.size() != D) {
                throw std::runtime_error("Point dimension mismatch");
            }
            Point arrPoint;
            for (int i = 0; i < D; ++i) {
                arrPoint[i] = vec[i];
            }
            allPoints.push_back(arrPoint);
        }

        // 插入无穷远点
        allPoints.insert(allPoints.end(), INFINITYPOINTS.begin(), INFINITYPOINTS.end());

        // 转换为 Qhull C++ 接口所需的 double 数组
        std::vector<double> qhullPoints;
        qhullPoints.reserve(allPoints.size() * D);
        for (const auto& p : allPoints) {
            qhullPoints.insert(qhullPoints.end(), p.begin(), p.end());
        }

        // 创建 Qhull 对象并计算凸包
        orgQhull::Qhull q;
//        q.runQhull("", D, static_cast<int>(allPoints.size()), qhullPoints.data(), "Qt Q12 Pp");
        q.runQhull("", D, static_cast<int>(allPoints.size()), qhullPoints.data(), "Qt");

        // 使用 std::set 去重
        std::set<std::vector<double>> uniquePoints; // 用于存储唯一的点

        // 提取凸包顶点
        for (const auto& facet : q.facetList()) {
            if (!facet.isGood()) {
                continue;
            }

            for (const auto& vertex : facet.vertices()) {
                std::vector<double> doublePoint(D);
                for (int i = 0; i < D; ++i) {
                    doublePoint[i] = vertex.point()[i];
                }

                Point checkPoint;
                for (int i = 0; i < D; ++i) {
                    checkPoint[i] = doublePoint[i];
                }

                // 检查是否为无穷点
                if (std::find(INFINITYPOINTS.begin(), INFINITYPOINTS.end(), checkPoint) == INFINITYPOINTS.end()) {
                    uniquePoints.insert(std::move(doublePoint)); // 自动去重
                }
            }
        }

        // 转换回 std::vector 并替换 points
        points.clear();
        points.assign(uniquePoints.begin(), uniquePoints.end());

        // 标记为已初始化
        initialized = true;
    }
    catch (const orgQhull::QhullError& e) {
        std::cerr << "Qhull Error: " << e.what() << std::endl;
        throw;
    }
}

/**
 * @brief 查询最短路径
 */
inline Path CH::NaiveQuery(const std::vector<double>& weight) {
    if (points.empty()) {
        return {{}, DBL_MAX};
    }

    double minDist = DotProduct(weight, points[0]);
    const std::vector<double>* bestPoint = &points[0];

    // 使用指针避免vector拷贝，展开循环减少函数调用
    for (size_t i = 1; i < points.size(); ++i) {
        const auto& p = points[i];
        double dist = 0.0;

        // 内联点积计算避免函数调用开销
        for (size_t j = 0; j < D; ++j) {
            dist += weight[j] * p[j];
        }

        if (dist < minDist) {
            minDist = dist;
            bestPoint = &p;
        }
    }
    return {*bestPoint, minDist};
}

/**
 * @brief 计算扩展点
 */
void CH::ComputeExPoints() {
    if (exPoints.empty()) {
        for (unsigned i = 0; i < D; ++i) {
            vector<double> weight(D, 0.0);
            weight[i] = 1.0;
            double minDist = DBL_MAX;
            vector<double> cost;

            for (auto &p : points) {
                // 假设 points 是 vector<vector<double>>
                if (p[i] < minDist) {
                    minDist = p[i];
                    cost = p; // 现在 cost 和 p 的类型相同，都是 vector<double>
                }
            }
            exPoints.emplace_back(cost, minDist);
        }
    }
}


/**
 * @brief 输出点集
 */
void CH::Print(std::ofstream& out) {
    out << points.size() << ' ';
    for (const auto& p : points) {
        for (const auto& c : p) {
            out << c << ' ';
        }
    }
}

#endif // POP_CH_HPP