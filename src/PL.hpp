/**
 ********************************************
 * @file    :PL.hpp
 * @author  :XXY
 * @brief   :Point Location
 * @date    :2025/11/30
 ********************************************
 */

#ifndef POP_POINTLOCATION_HPP
#define POP_POINTLOCATION_HPP

#include "common.hpp"
#include <Eigen/Dense>
#include <spatialindex/SpatialIndex.h>
#include <unordered_map>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <memory>
#include <random>
#include <chrono>

const auto ep = 10 * std::numeric_limits<double>::epsilon();

// 优化的预计算结果结构
struct PrecomputedResult {
    std::vector<double> direction;
    std::vector<double> optimalPoint;
    double cost;
    unsigned vertexID;
    double confidence;

    // 添加压缩存储选项
    PrecomputedResult() = default;

    // 减少内存占用的构造函数
    PrecomputedResult(const std::vector<double>& dir, const std::vector<double>& point,
                      double c, unsigned vID, double conf)
            : cost(c), vertexID(vID), confidence(conf) {
        // 只存储归一化后的方向，减少精度损失
        direction = dir;
        optimalPoint = point;
    }
};


// 优化的查询缓存结构
struct QueryCacheEntry {
    std::vector<double> result_point;
    double cost;
    std::chrono::steady_clock::time_point timestamp;

    // 添加内存优化的构造函数
    QueryCacheEntry() = default;

    QueryCacheEntry(const std::vector<double>& point, double c)
            : result_point(point), cost(c),
              timestamp(std::chrono::steady_clock::now()) {}

    // 检查是否过期
    bool isExpired() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::minutes>(now - timestamp).count() > 5;
    }
};

// 优化的 Visitor 类 - 预分配内存
class OptimizedVisitor : public SpatialIndex::IVisitor {
private:
    static constexpr size_t INITIAL_CAPACITY = 64;

public:
    std::vector<SpatialIndex::id_type> result_ids;

    OptimizedVisitor() {
        result_ids.reserve(INITIAL_CAPACITY);
    }

    void visitNode(const SpatialIndex::INode& /*n*/) override {}

    void visitData(const SpatialIndex::IData& d) override {
        result_ids.emplace_back(d.getIdentifier());
    }

    void visitData(std::vector<const SpatialIndex::IData*>& /*v*/) override {}

    inline void clear() {
        result_ids.clear();
        if (result_ids.capacity() > INITIAL_CAPACITY * 2) {
            result_ids.shrink_to_fit();
            result_ids.reserve(INITIAL_CAPACITY);
        }
    }
};

class PL {
private:
    // 核心数据结构
    SpatialIndex::ISpatialIndex* rtree = nullptr;
    SpatialIndex::IStorageManager* storage = nullptr;
    std::vector<std::pair<double, std::vector<unsigned>>> segment_data;
    std::vector<unsigned> boxid2vertexid;
    std::vector<std::vector<double>> vertexPoints;

    // 预存储优化相关
    std::vector<PrecomputedResult> precomputedResults;
    std::unordered_map<size_t, std::vector<unsigned>> directionToVertices;
    std::vector<std::vector<double>> representativeDirections;

    // 查询缓存
    mutable std::unordered_map<size_t, QueryCacheEntry> queryCache;
    static constexpr size_t MAX_CACHE_SIZE = 128;
    static constexpr std::chrono::minutes CACHE_EXPIRE_TIME{5};

    // 快速查找表
    std::vector<std::pair<double, unsigned>> sortedVerticesByCost;
    std::unordered_map<unsigned, std::vector<unsigned>> vertexNeighbors;

    // 性能优化相关
    mutable OptimizedVisitor cached_visitor;
    mutable std::vector<double> temp_buffer;
    mutable std::vector<double> temp_reduced_weight;
    mutable std::vector<double> temp_mins, temp_maxs;
    bool is_2d;
    size_t vertex_count;

    // 统计信息
    mutable size_t cache_hits = 0;
    mutable size_t total_queries = 0;

    // =========================================================================
    // L1归一化函数 - 分量和归一化
    // =========================================================================
    inline bool normalizeVector_L1(std::vector<double>& v, double& out_sum_val) const {
        double sum = 0.0;
        for (double val : v) {
            sum += val;  // 计算分量和
        }

        if (std::abs(sum) < 1e-12) { // 分量和接近零
            out_sum_val = 0.0;
            return false; // 无法归一化
        }

        out_sum_val = sum;
        for (double& val : v) {
            val /= sum;  // 每个分量除以分量和
        }
        return true;
    }


    // =========================================================================

public:
    explicit PL(const std::vector<std::vector<double>>& points);

    ~PL() {
        if (rtree) {
            delete rtree;
            rtree = nullptr;
        }
        if (storage) {
            delete storage;
            storage = nullptr;
        }
    }

    // 禁用拷贝构造和赋值
    PL(const PL&) = delete;
    PL& operator=(const PL&) = delete;

    // 允许移动构造和赋值
    PL(PL&& other) noexcept
            : rtree(other.rtree), storage(other.storage),
              segment_data(std::move(other.segment_data)),
              boxid2vertexid(std::move(other.boxid2vertexid)),
              vertexPoints(std::move(other.vertexPoints)),
              precomputedResults(std::move(other.precomputedResults)),
              directionToVertices(std::move(other.directionToVertices)),
              representativeDirections(std::move(other.representativeDirections)),
              queryCache(std::move(other.queryCache)),
              sortedVerticesByCost(std::move(other.sortedVerticesByCost)),
              vertexNeighbors(std::move(other.vertexNeighbors)),
              cached_visitor(std::move(other.cached_visitor)),
              temp_buffer(std::move(other.temp_buffer)),
              temp_reduced_weight(std::move(other.temp_reduced_weight)),
              temp_mins(std::move(other.temp_mins)),
              temp_maxs(std::move(other.temp_maxs)),
              is_2d(other.is_2d), vertex_count(other.vertex_count),
              cache_hits(other.cache_hits), total_queries(other.total_queries) {
        // 清空原对象的指针，避免双重释放
        other.rtree = nullptr;
        other.storage = nullptr;
    }

    PL& operator=(PL&& other) noexcept {
        if (this != &other) {
            // 先清理当前对象的资源
            if (rtree) delete rtree;
            if (storage) delete storage;

            // 移动资源
            rtree = other.rtree;
            storage = other.storage;
            segment_data = std::move(other.segment_data);
            boxid2vertexid = std::move(other.boxid2vertexid);
            vertexPoints = std::move(other.vertexPoints);
            precomputedResults = std::move(other.precomputedResults);
            directionToVertices = std::move(other.directionToVertices);
            representativeDirections = std::move(other.representativeDirections);
            queryCache = std::move(other.queryCache);
            sortedVerticesByCost = std::move(other.sortedVerticesByCost);
            vertexNeighbors = std::move(other.vertexNeighbors);
            cached_visitor = std::move(other.cached_visitor);
            temp_buffer = std::move(other.temp_buffer);
            temp_reduced_weight = std::move(other.temp_reduced_weight);
            temp_mins = std::move(other.temp_mins);
            temp_maxs = std::move(other.temp_maxs);
            is_2d = other.is_2d;
            vertex_count = other.vertex_count;
            cache_hits = other.cache_hits;
            total_queries = other.total_queries;

            // 清空原对象的指针
            other.rtree = nullptr;
            other.storage = nullptr;
        }
        return *this;
    }

    [[nodiscard]] Path PLQuery(const std::vector<double>& w) const;

    // 性能监控
    double getCacheHitRate() const {
        return total_queries > 0 ? double(cache_hits) / total_queries : 0.0;
    }

    size_t getIndexSize() const;

private:
    // 预存储构建函数
    void buildPrecomputedResults();
    void buildRepresentativeDirections();
    void buildVertexNeighbors();
    void buildFastLookupTables();

    // 查询优化函数 - 内联优化
    inline Path fastPrecomputedQuery(const std::vector<double>& w) const;
    inline Path cachedQuery(const std::vector<double>& w) const;
    inline Path neighborhoodQuery(const std::vector<double>& w, unsigned startVertex) const;

    // 工具函数 - 内联优化
    inline size_t hashVector(const std::vector<double>& v) const;
    inline double fastVectorSimilarity(const std::vector<double>& a, const std::vector<double>& b) const;
    inline unsigned findNearestVertex(const std::vector<double>& w) const;
    inline void cleanCache() const;

    inline bool isInfinityPoint(const std::vector<double>& p) const {
        Point checkP;
        std::copy(p.begin(), p.begin() + D, checkP.begin());
        return std::find(INFINITYPOINTS.begin(), INFINITYPOINTS.end(), checkP) != INFINITYPOINTS.end();
    }

    inline void computeBoundingBox(const std::vector<std::array<double, D-1>>& reducedNormals,
                                   std::vector<double>& mins, std::vector<double>& maxs) const {
        std::fill(mins.begin(), mins.end(), DBL_MAX);
        std::fill(maxs.begin(), maxs.end(), -DBL_MAX);

        for (const auto& rn : reducedNormals) {
            for (unsigned d = 0; d < D - 1; ++d) {
                if (rn[d] < mins[d]) mins[d] = rn[d];
                if (rn[d] > maxs[d]) maxs[d] = rn[d];
            }
        }
    }
};

// =========================================================================
// PL Constructor
// =========================================================================
PL::PL(const std::vector<std::vector<double>>& points)
        : is_2d(D == 2), temp_buffer(D), temp_reduced_weight(D-1), temp_mins(D-1), temp_maxs(D-1) {

    // 预分配内存
    std::vector<Point> allPoints;
    allPoints.reserve(points.size() + INFINITYPOINTS.size());

    // 添加输入点
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        if (static_cast<int>(p.size()) != D) {
            throw std::runtime_error("Point dimension mismatch");
        }
        Point arrPoint;
        std::copy(p.begin(), p.end(), arrPoint.begin());
        allPoints.emplace_back(std::move(arrPoint));
    }

    allPoints.insert(allPoints.end(), INFINITYPOINTS.begin(), INFINITYPOINTS.end());

    std::vector<double> coords;
    coords.reserve(allPoints.size() * D);
    for (const auto& ap : allPoints) {
        coords.insert(coords.end(), ap.begin(), ap.end());
    }

    // 使用Qhull计算凸包
    orgQhull::Qhull q;
    try {
        q.runQhull("", D, static_cast<int>(allPoints.size()), coords.data(), "Qt Qx");
        if (q.qhullStatus() != 0) {
            throw std::runtime_error("Qhull failed to compute the convex hull");
        }
    } catch (const std::exception& e) {
        std::cerr << "Qhull失败: " << e.what() << std::endl;
        throw;
    }

    // 计算有效顶点数
    size_t valid_vertices = 0;
    for (const auto& v : q.vertexList()) {
        for (int i = 0; i < D; ++i) {
            temp_buffer[i] = v.point().coordinates()[i];
        }
        if (!isInfinityPoint(temp_buffer)) {
            valid_vertices++;
        }
    }

    // 预分配内存
    vertexPoints.reserve(valid_vertices);
    vertex_count = valid_vertices;

    // R*-Tree初始化
    if (!is_2d) {
        try {
            storage = SpatialIndex::StorageManager::createNewMemoryStorageManager();
            if (!storage) {
                throw std::runtime_error("Failed to create storage manager");
            }

            SpatialIndex::id_type indexIdentifier = 1;
            rtree = SpatialIndex::RTree::createNewRTree(
                    *storage, 0.7, 32, 32, D - 1,
                    SpatialIndex::RTree::RV_RSTAR, indexIdentifier
            );
            if (!rtree) {
                delete storage;
                storage = nullptr;
                throw std::runtime_error("Failed to create R*-Tree");
            }

            boxid2vertexid.reserve(valid_vertices);
        } catch (const std::exception& e) {
            if (rtree) {
                delete rtree;
                rtree = nullptr;
            }
            if (storage) {
                delete storage;
                storage = nullptr;
            }
            throw std::runtime_error(std::string("R*-Tree initialization failed: ") + e.what());
        }
    } else {
        segment_data.reserve(valid_vertices);
    }

    // 分析面片并建立映射
    std::unordered_map<const double*, std::vector<orgQhull::QhullFacet>> vertexFacets;
    vertexFacets.reserve(valid_vertices * 2);

    for (const auto& facet : q.facetList()) {
        if (!facet.isGood()) continue;

        try {
            auto hp = facet.hyperplane();
            const coordT* coords = hp.coordinates();
            if (!coords) continue;

            // 检查法向量有效性
            bool valid_normal = false;
            for (int i = 0; i < D; ++i) {
                if (std::abs(coords[i]) > 1e-12) {
                    valid_normal = true;
                    break;
                }
            }
            if (!valid_normal) continue;

            // 记录顶点-面片关系
            for (const auto& vertex : facet.vertices()) {
                const double* vcoords = vertex.point().coordinates();
                vertexFacets[vcoords].emplace_back(facet);
            }
        } catch (const std::exception&) {
            continue;
        }
    }

    // 处理顶点并构建降维包围盒
    SpatialIndex::id_type box_id = 0;
    unsigned vertexID = 0;
    std::vector<std::array<double, D-1>> reducedNormals;
    reducedNormals.reserve(8);

    for (const auto& v : q.vertexList()) {
        for (int i = 0; i < D; ++i) {
            temp_buffer[i] = v.point().coordinates()[i];
        }

        if (isInfinityPoint(temp_buffer)) continue;

        const double* vpoint = v.point().coordinates();
        auto it = vertexFacets.find(vpoint);
        if (it == vertexFacets.end() || it->second.empty()) {
            continue;
        }

        reducedNormals.clear();

        for (const auto& f : it->second) {
            try {
                auto hp = f.hyperplane();
                const coordT* coords = hp.coordinates();
                if (!coords) continue;

                std::vector<double> current_normal(D);
                for(int i = 0; i < D; ++i) {
                    current_normal[i] = static_cast<double>(coords[i]);
                }

                double sum_val;
                if (!normalizeVector_L1(current_normal, sum_val)) {
                    continue;
                }

                // 确保最后一个分量为非负
                if (current_normal[D-1] < 0) {
                    for(unsigned d = 0; d < D; ++d) {
                        current_normal[d] = -current_normal[d];
                    }
                }

                // 降维：去掉最后一维
                std::array<double, D-1> reduced;
                for (unsigned d = 0; d < D-1; ++d) {
                    reduced[d] = current_normal[d];
                }

                reducedNormals.emplace_back(reduced);

            } catch (const std::exception&) {
                continue;
            }
        }

        if (reducedNormals.empty()) {
            continue;
        }

        // 计算包围盒
        computeBoundingBox(reducedNormals, temp_mins, temp_maxs);

        vertexPoints.emplace_back(temp_buffer);

        if (!is_2d && rtree) {
            try {
                SpatialIndex::Region box(temp_mins.data(), temp_maxs.data(), D - 1);
                rtree->insertData(0, nullptr, box, box_id);

                if (boxid2vertexid.size() <= box_id) {
                    boxid2vertexid.resize(box_id + 1);
                }
                boxid2vertexid[box_id] = vertexID;
                box_id++;
            } catch (const std::exception&) {
                // 忽略插入失败
            }
        } else if (is_2d) {
            double key = (temp_mins[0] + temp_maxs[0]) * 0.5;
            segment_data.emplace_back(key, std::vector<unsigned>{vertexID});
        }

        vertexID++;
    }

    // 2D后处理优化
    if (is_2d && !segment_data.empty()) {
        std::sort(segment_data.begin(), segment_data.end());

        auto write_it = segment_data.begin();
        for (auto read_it = segment_data.begin(); read_it != segment_data.end(); ++read_it) {
            if (write_it != read_it && std::abs(write_it->first - read_it->first) < ep) {
                write_it->second.insert(write_it->second.end(),
                                        read_it->second.begin(), read_it->second.end());
            } else if (write_it != read_it) {
                ++write_it;
                *write_it = std::move(*read_it);
            }
        }
        segment_data.erase(++write_it, segment_data.end());
        segment_data.shrink_to_fit();
    }

    // 构建预存储优化
    try {
        buildRepresentativeDirections();
        buildPrecomputedResults();
        buildVertexNeighbors();
        buildFastLookupTables();
    } catch (const std::exception&) {
        // 优化构建失败时继续运行，不抛出异常
    }
}



// =========================================================================
// buildRepresentativeDirections
// =========================================================================
void PL::buildRepresentativeDirections() {
    // 大幅减少预计算方向数量
    const size_t base_directions = std::min(static_cast<size_t>(D * 2), static_cast<size_t>(15)); // 从D*3和30减少
    const size_t max_additional = std::min(static_cast<size_t>(15), vertex_count / 20); // 从40和vertex_count/15减少
    const size_t total_directions = base_directions + max_additional;

    representativeDirections.clear();
    representativeDirections.reserve(total_directions);

    // 1. 添加坐标轴主导方向（减少每维的方向数）
    for (int i = 0; i < D && representativeDirections.size() < total_directions; ++i) {
        // 只保留正向主导和一个对角方向
        std::vector<double> dir_pos(D, 0.05);
        dir_pos[i] = 0.9;
        double sum_val;
        if (normalizeVector_L1(dir_pos, sum_val)) {
            representativeDirections.emplace_back(std::move(dir_pos));
        }

        if (representativeDirections.size() < total_directions) {
            std::vector<double> dir_diag(D, 0.15);
            dir_diag[i] = 0.3;
            if (normalizeVector_L1(dir_diag, sum_val)) {
                representativeDirections.emplace_back(std::move(dir_diag));
            }
        }
    }

    // 2. 填充剩余位置为均匀分布方向
    while (representativeDirections.size() < total_directions) {
        std::vector<double> dir(D, 1.0 / D);
        representativeDirections.emplace_back(std::move(dir));
    }

    // 立即收缩容器
    representativeDirections.shrink_to_fit();
}


// =========================================================================
// buildPrecomputedResults (No change needed, it uses already normalized directions)
// =========================================================================
void PL::buildPrecomputedResults() {
    precomputedResults.reserve(representativeDirections.size());

    for (const auto& direction : representativeDirections) {
        double minCost = DBL_MAX;
        std::vector<double> bestPoint;
        unsigned bestVertexID = 0;

        // 优化：内联计算避免函数调用
        for (size_t vID = 0; vID < vertexPoints.size(); ++vID) {
            const auto& vPoint = vertexPoints[vID];
            double cost = 0.0;

            // 展开循环提高性能
            if constexpr (D <= 4) {
                for (size_t d = 0; d < D; ++d) {
                    cost += direction[d] * vPoint[d];
                }
            } else {
                for (size_t d = 0; d < D; ++d) {
                    cost += direction[d] * vPoint[d];
                }
            }

            if (cost < minCost) {
                minCost = cost;
                bestPoint = vPoint;
                bestVertexID = static_cast<unsigned>(vID);
            }
        }

        PrecomputedResult result;
        result.direction = direction;
        result.optimalPoint = bestPoint;
        result.cost = minCost;
        result.vertexID = bestVertexID;
        result.confidence = 1.0;

        precomputedResults.emplace_back(std::move(result));

        size_t hash = hashVector(direction);
        directionToVertices[hash].push_back(bestVertexID);
    }
}

// =========================================================================
// buildVertexNeighbors (No changes needed, as it depends on point coordinates directly)
// =========================================================================
void PL::buildVertexNeighbors() {
    if (vertexPoints.empty()) return;

    vertexNeighbors.clear();

    // 更严格的邻居参数
    double neighborRadiusSq;
    size_t max_neighbors_per_vertex;

    if (vertex_count > 2000) {
        neighborRadiusSq = 0.0005;  // 更严格
        max_neighbors_per_vertex = 3;  // 进一步减少
    } else if (vertex_count > 1000) {
        neighborRadiusSq = 0.001;   // 更严格
        max_neighbors_per_vertex = 4;  // 减少
    } else if (vertex_count > 500) {
        neighborRadiusSq = 0.003;
        max_neighbors_per_vertex = 6;
    } else {
        neighborRadiusSq = 0.008;   // 稍微减少
        max_neighbors_per_vertex = 10;  // 减少
    }

    for (size_t i = 0; i < vertexPoints.size(); ++i) {
        std::vector<std::pair<double, unsigned>> candidates;
        candidates.reserve(max_neighbors_per_vertex);  // 精确预留

        for (size_t j = 0; j < vertexPoints.size(); ++j) {
            if (i == j) continue;

            double distanceSq = 0.0;
            for (size_t d = 0; d < D; ++d) {
                double diff = vertexPoints[i][d] - vertexPoints[j][d];
                distanceSq += diff * diff;
                if (distanceSq > neighborRadiusSq) break;
            }

            if (distanceSq <= neighborRadiusSq) {
                candidates.emplace_back(distanceSq, static_cast<unsigned>(j));
                // 提前截断避免过多候选
                if (candidates.size() >= max_neighbors_per_vertex * 2) {
                    break;
                }
            }
        }

        // 只保留最近的邻居
        if (candidates.size() > max_neighbors_per_vertex) {
            std::partial_sort(candidates.begin(),
                              candidates.begin() + max_neighbors_per_vertex,
                              candidates.end());
            candidates.resize(max_neighbors_per_vertex);
        }

        if (!candidates.empty()) {
            std::vector<unsigned> neighbors;
            neighbors.reserve(candidates.size());
            for (const auto& candidate : candidates) {
                neighbors.push_back(candidate.second);
            }
            neighbors.shrink_to_fit();  // 立即收缩
            vertexNeighbors[static_cast<unsigned>(i)] = std::move(neighbors);
        }
    }
}



// =========================================================================
// buildFastLookupTables (Modified to align with L2 normalization)
// =========================================================================
void PL::buildFastLookupTables() {
    // 只在小数据集上构建查找表
    if (vertex_count > 1000) {
        sortedVerticesByCost.clear();
        return;
    }

    std::vector<double> ref_direction(D, 1.0 / std::sqrt(static_cast<double>(D)));

    sortedVerticesByCost.clear();
    sortedVerticesByCost.reserve(vertexPoints.size());

    for (size_t i = 0; i < vertexPoints.size(); ++i) {
        double cost = 0.0;
        for (size_t d = 0; d < D; ++d) {
            cost += ref_direction[d] * vertexPoints[i][d];
        }
        sortedVerticesByCost.emplace_back(cost, static_cast<unsigned>(i));
    }

    std::sort(sortedVerticesByCost.begin(), sortedVerticesByCost.end());
    sortedVerticesByCost.shrink_to_fit();
}


// =========================================================================
// PLQuery
// =========================================================================
Path PL::PLQuery(const std::vector<double>& w_raw) const {
    total_queries++;

    // 1. 缓存查询
    {
        size_t hash = hashVector(w_raw);
        auto it = queryCache.find(hash);
        if (it != queryCache.end()) {
            auto now = std::chrono::steady_clock::now();
            if (now - it->second.timestamp < CACHE_EXPIRE_TIME) {
                cache_hits++;
                return {it->second.result_point, it->second.cost};
            } else {
                queryCache.erase(it);
            }
        }
    }

    // 2. 预计算结果查询
    Path precomputedResult = fastPrecomputedQuery(w_raw);
    if (precomputedResult.second != DBL_MAX) {
        size_t hash = hashVector(w_raw);
        QueryCacheEntry entry;
        entry.result_point = precomputedResult.first;
        entry.cost = precomputedResult.second;
        entry.timestamp = std::chrono::steady_clock::now();
        cleanCache();
        queryCache[hash] = entry;
        return precomputedResult;
    }

    std::vector<double> normalized_w = w_raw;
    double sum_val;

    if (!normalizeVector_L1(normalized_w, sum_val)) {
        std::cerr << "WARNING: L1归一化失败（分量和为零），回退到暴力搜索" << std::endl;
        goto brute_force_fallback;
    }

    // 3. 使用归一化后的向量进行查询
    if (!is_2d && rtree) {
        // 降维：去掉最后一维，使用归一化后的向量
        for (size_t i = 0; i < D - 1; ++i) {
            temp_reduced_weight[i] = normalized_w[i];
            temp_mins[i] = temp_reduced_weight[i] - ep;
            temp_maxs[i] = temp_reduced_weight[i] + ep;
        }

        try {
            SpatialIndex::Region queryBox(temp_mins.data(), temp_maxs.data(), D - 1);
            cached_visitor.clear();
            rtree->intersectsWithQuery(queryBox, cached_visitor);

            double min_cost = DBL_MAX;
            const std::vector<double>* best_point = nullptr;

            const auto& result_ids = cached_visitor.result_ids;
            for (size_t idx = 0; idx < result_ids.size(); ++idx) {
                SpatialIndex::id_type id = result_ids[idx];
                if (id >= boxid2vertexid.size()) continue;

                unsigned vID = boxid2vertexid[id];
                if (vID >= vertexPoints.size()) continue;

                const std::vector<double>& vPoint = vertexPoints[vID];

                // 使用原始权重计算代价
                double cost = 0.0;
                for (size_t i = 0; i < D; ++i) {
                    cost += w_raw[i] * vPoint[i];
                }

                if (cost < min_cost) {
                    min_cost = cost;
                    best_point = &vPoint;
                }
            }

            if (best_point) {
                unsigned startVertex = findNearestVertex(w_raw);
                Path neighborResult = neighborhoodQuery(w_raw, startVertex);
                if (neighborResult.second < min_cost) {
                    size_t hash = hashVector(w_raw);
                    QueryCacheEntry entry;
                    entry.result_point = neighborResult.first;
                    entry.cost = neighborResult.second;
                    entry.timestamp = std::chrono::steady_clock::now();
                    cleanCache();
                    queryCache[hash] = entry;
                    return neighborResult;
                }

                size_t hash = hashVector(w_raw);
                QueryCacheEntry entry;
                entry.result_point = *best_point;
                entry.cost = min_cost;
                entry.timestamp = std::chrono::steady_clock::now();
                cleanCache();
                queryCache[hash] = entry;

                return {*best_point, min_cost};
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: R*-Tree查询失败，回退到暴力搜索: " << e.what() << std::endl;
            goto brute_force_fallback;
        }

    } else if (is_2d) {
        // 2D情况：使用归一化后向量的第一个分量
        double target = normalized_w[0] - ep;
        double target_upper = normalized_w[0] + ep;

        auto lower = std::lower_bound(segment_data.begin(), segment_data.end(),
                                      std::make_pair(target, std::vector<unsigned>()));
        auto upper = std::upper_bound(segment_data.begin(), segment_data.end(),
                                      std::make_pair(target_upper, std::vector<unsigned>()));

        double min_cost = DBL_MAX;
        const std::vector<double>* best_point = nullptr;

        for (auto it = lower; it != upper; ++it) {
            for (unsigned vID : it->second) {
                if (vID >= vertexPoints.size()) continue;

                const std::vector<double>& vPoint = vertexPoints[vID];
                // 使用原始权重计算代价
                double cost = w_raw[0] * vPoint[0] + w_raw[1] * vPoint[1];

                if (cost < min_cost) {
                    min_cost = cost;
                    best_point = &vPoint;
                }
            }
        }

        if (best_point) {
            unsigned startVertex = findNearestVertex(w_raw);
            Path neighborResult = neighborhoodQuery(w_raw, startVertex);
            if (neighborResult.second < min_cost) {
                size_t hash = hashVector(w_raw);
                QueryCacheEntry entry;
                entry.result_point = neighborResult.first;
                entry.cost = neighborResult.second;
                entry.timestamp = std::chrono::steady_clock::now();
                cleanCache();
                queryCache[hash] = entry;
                return neighborResult;
            }

            size_t hash = hashVector(w_raw);
            QueryCacheEntry entry;
            entry.result_point = *best_point;
            entry.cost = min_cost;
            entry.timestamp = std::chrono::steady_clock::now();
            cleanCache();
            queryCache[hash] = entry;

            return {*best_point, min_cost};
        }
    }

    brute_force_fallback:
    // 暴力搜索回退
    double min_cost = DBL_MAX;
    std::vector<double> best_point;

    for (const auto& p : vertexPoints) {
        double cost = 0.0;
        if constexpr (D <= 4) {
            for (size_t i = 0; i < D; ++i) {
                cost += w_raw[i] * p[i];
            }
        } else {
            for (size_t i = 0; i < D; ++i) {
                cost += w_raw[i] * p[i];
            }
        }

        if (cost < min_cost) {
            min_cost = cost;
            best_point = p;
        }
    }

    if (!best_point.empty()) {
        size_t hash = hashVector(w_raw);
        QueryCacheEntry entry;
        entry.result_point = best_point;
        entry.cost = min_cost;
        entry.timestamp = std::chrono::steady_clock::now();
        cleanCache();
        queryCache[hash] = entry;

        // 定期清理缓冲区
        if (total_queries % 100 == 0) {
            if (cached_visitor.result_ids.capacity() > 256) {
                const_cast<OptimizedVisitor&>(cached_visitor).result_ids.shrink_to_fit();
                const_cast<OptimizedVisitor&>(cached_visitor).result_ids.reserve(64);
            }
        }

        return {best_point, min_cost};
    }

    throw std::runtime_error("No valid solution found");
}


// =========================================================================
// fastPrecomputedQuery (Modified for L2 normalization)
// =========================================================================
inline Path PL::fastPrecomputedQuery(const std::vector<double>& w_raw) const {
    if (precomputedResults.empty()) {
        return {{}, DBL_MAX};
    }

    // 对查询向量进行L1归一化
    std::vector<double> normalized_w = w_raw;
    double sum_val;

    if (!normalizeVector_L1(normalized_w, sum_val)) {
        return {{}, DBL_MAX};
    }

    double maxSimilarity = -1.0;
    const PrecomputedResult* bestResult = nullptr;

    // 计算与预存方向的相似度
    for (const auto& result : precomputedResults) {
        double similarity = 0.0;

        // 点积计算（两个L1归一化向量的点积）
        for (size_t i = 0; i < D; ++i) {
            similarity += normalized_w[i] * result.direction[i];
        }

        if (similarity > maxSimilarity) {
            maxSimilarity = similarity;
            bestResult = &result;
        }
    }

    if (bestResult && maxSimilarity > 0.82) {
        // 使用原始权重计算实际代价
        double actualCost = 0.0;
        for (size_t i = 0; i < D; ++i) {
            actualCost += w_raw[i] * bestResult->optimalPoint[i];
        }
        return {bestResult->optimalPoint, actualCost};
    }

    return {{}, DBL_MAX};
}


// =========================================================================
// cachedQuery (Uses hashVector, which takes original w. No changes needed directly)
// =========================================================================
inline Path PL::cachedQuery(const std::vector<double>& w) const {
    size_t hash = hashVector(w);
    auto it = queryCache.find(hash);

    if (it != queryCache.end()) {
        auto now = std::chrono::steady_clock::now();
        if (now - it->second.timestamp < CACHE_EXPIRE_TIME) {
            return {it->second.result_point, it->second.cost};
        } else {
            queryCache.erase(it);
        }
    }

    return {{}, DBL_MAX};
}

// =========================================================================
// neighborhoodQuery (Uses original w for cost calculation. No changes needed directly)
// =========================================================================
inline Path PL::neighborhoodQuery(const std::vector<double>& w, unsigned startVertex) const {
    if (startVertex >= vertexPoints.size()) {
        return {{}, DBL_MAX};
    }

    double bestCost = DBL_MAX;
    std::vector<double> bestPoint;

    const auto& startPoint = vertexPoints[startVertex];
    double startCost = 0.0;
    for (size_t i = 0; i < D; ++i) {
        startCost += w[i] * startPoint[i];
    }
    bestCost = startCost;
    bestPoint = startPoint;

    auto it = vertexNeighbors.find(startVertex);
    if (it != vertexNeighbors.end()) {
        for (unsigned neighborID : it->second) {
            if (neighborID >= vertexPoints.size()) continue;

            const auto& neighborPoint = vertexPoints[neighborID];

            double neighborCost = 0.0;
            for (size_t i = 0; i < D; ++i) {
                neighborCost += w[i] * neighborPoint[i];
            }

            if (neighborCost < bestCost) {
                bestCost = neighborCost;
                bestPoint = neighborPoint;
            }
        }
    }

    return {bestPoint, bestCost};
}

// =========================================================================
// hashVector (Uses raw vector elements. No changes needed directly)
// =========================================================================
inline size_t PL::hashVector(const std::vector<double>& v) const {
    size_t hash = 0;
    const size_t prime = 1099511628211ULL;
    const size_t offset_basis = 14695981039346656037ULL;

    hash = offset_basis;
    for (size_t i = 0; i < v.size(); ++i) {
        // 只处理极小值，保持原有量化精度
        // Note: This quantization for hashing might be problematic for floats.
        // Consider converting double to its bit pattern (uint64_t) for hashing.
        double val = (std::abs(v[i]) < 1e-12) ? 0.0 : v[i];
        int64_t quantized = static_cast<int64_t>(val * 10000.0); // 保持原精度
        hash ^= static_cast<size_t>(quantized);
        hash *= prime;
    }
    return hash;
}

// =========================================================================
// fastVectorSimilarity (Assumes inputs are normalized with L2. No changes needed directly)
// =========================================================================
inline double PL::fastVectorSimilarity(const std::vector<double>& a, const std::vector<double>& b) const {
    if (a.size() != b.size()) return 0.0;

    // 优化：假设向量已L2归一化，直接计算点积即为余弦相似度
    double dot_product = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        dot_product += a[i] * b[i];
    }

    return dot_product;
}

// =========================================================================
// findNearestVertex (Modified to align with L2 normalization)
// =========================================================================
inline unsigned PL::findNearestVertex(const std::vector<double>& w_raw) const {
    if (sortedVerticesByCost.empty()) return 0;

    // Create a reference direction (e.g., [1/sqrt(D), ..., 1/sqrt(D)])
    std::vector<double> ref_direction(D);
    double inv_sqrt_d = 1.0 / std::sqrt(static_cast<double>(D));
    for (size_t d = 0; d < D; ++d) {
        ref_direction[d] = inv_sqrt_d;
    }

    // Calculate dot product of w_raw with reference direction
    double referenceValue = 0.0;
    for (size_t i = 0; i < D; ++i) {
        referenceValue += w_raw[i] * ref_direction[i];
    }

    auto it = std::lower_bound(sortedVerticesByCost.begin(), sortedVerticesByCost.end(),
                               std::make_pair(referenceValue, 0u));

    if (it == sortedVerticesByCost.end()) {
        return sortedVerticesByCost.back().second;
    }
    // Consider also checking 'it-1' if 'it' points to something much larger
    if (it != sortedVerticesByCost.begin()) {
        if (std::abs(it->first - referenceValue) > std::abs((it - 1)->first - referenceValue)) {
            return (it - 1)->second;
        }
    }

    return it->second;
}

// =========================================================================
// cleanCache (No changes needed)
// =========================================================================
inline void PL::cleanCache() const {
    // 使用更小的缓存大小限制
    const size_t REDUCED_MAX_CACHE_SIZE = 128;

    if (queryCache.size() <= REDUCED_MAX_CACHE_SIZE) {
        return;
    }

    // 更激进的清理策略：先清理过期条目
    auto now = std::chrono::steady_clock::now();
    for (auto it = queryCache.begin(); it != queryCache.end();) {
        if (std::chrono::duration_cast<std::chrono::minutes>(now - it->second.timestamp).count() > 5) {
            it = queryCache.erase(it);
        } else {
            ++it;
        }
    }

    // 如果仍然过大，删除最旧的2/3
    if (queryCache.size() > REDUCED_MAX_CACHE_SIZE) {
        size_t toRemove = queryCache.size() - REDUCED_MAX_CACHE_SIZE / 3;

        std::vector<std::unordered_map<size_t, QueryCacheEntry>::iterator> iterators;
        iterators.reserve(queryCache.size());

        for (auto it = queryCache.begin(); it != queryCache.end(); ++it) {
            iterators.push_back(it);
        }

        std::partial_sort(iterators.begin(),
                          iterators.begin() + toRemove,
                          iterators.end(),
                          [](const auto& a, const auto& b) {
                              return a->second.timestamp < b->second.timestamp;
                          });

        for (size_t i = 0; i < toRemove; ++i) {
            queryCache.erase(iterators[i]);
        }
    }
}



size_t PL::getIndexSize() const {
    size_t total_size = 0;

    // 顶点点集内存
    for (const auto& point : vertexPoints) {
        total_size += point.capacity() * sizeof(double); // 使用capacity而不是size
    }
    total_size += vertexPoints.capacity() * sizeof(std::vector<double>);

    // 预计算结果内存 - 现在数量大幅减少
    for (const auto& result : precomputedResults) {
        total_size += result.direction.capacity() * sizeof(double);
        total_size += result.optimalPoint.capacity() * sizeof(double);
        total_size += sizeof(double) + sizeof(unsigned) + sizeof(double);
        total_size += sizeof(std::vector<double>) * 2;
    }
    total_size += precomputedResults.capacity() * sizeof(PrecomputedResult);

    // 邻居关系内存 - 现在更紧凑
    for (const auto& pair : vertexNeighbors) {
        total_size += sizeof(unsigned);
        total_size += pair.second.capacity() * sizeof(unsigned);
        total_size += sizeof(std::vector<unsigned>);
    }

    // 查询缓存内存 - 现在更小
    for (const auto& cache_pair : queryCache) {
        total_size += cache_pair.second.result_point.capacity() * sizeof(double);
        total_size += sizeof(double) + sizeof(std::chrono::steady_clock::time_point);
        total_size += sizeof(std::vector<double>);
    }

    // 其他数据结构
    total_size += boxid2vertexid.capacity() * sizeof(unsigned);
    total_size += segment_data.capacity() * sizeof(std::pair<double, std::vector<unsigned>>);

    for (const auto& seg : segment_data) {
        total_size += seg.second.capacity() * sizeof(unsigned);
    }

    // directionToVertices - 现在更小
    for (const auto& pair : directionToVertices) {
        total_size += pair.second.capacity() * sizeof(unsigned);
    }

    // representativeDirections - 现在大幅减少
    for (const auto& dir : representativeDirections) {
        total_size += dir.capacity() * sizeof(double);
    }

    // 临时缓冲区
    total_size += temp_buffer.capacity() * sizeof(double);
    total_size += temp_reduced_weight.capacity() * sizeof(double);
    total_size += temp_mins.capacity() * sizeof(double);
    total_size += temp_maxs.capacity() * sizeof(double);

    // R*-Tree估算（保持不变）
    if (storage != nullptr) {
        total_size += boxid2vertexid.capacity() * 64; // 减少估算开销
    }

    return total_size;
}



#endif // POP_POINTLOCATION_HPP
