cmake_minimum_required(VERSION 3.20)
project(POP)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 禁用 CGAL 的 Debug 模式性能警告
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

# 查找必需的核心库
find_package(OpenMP REQUIRED)

# 查找 GMP
find_library(GMP_LIBRARIES NAMES gmp)
find_path(GMP_INCLUDE_DIRS NAMES gmp.h)
if(NOT GMP_LIBRARIES OR NOT GMP_INCLUDE_DIRS)
    message(FATAL_ERROR "GMP not found! Install: sudo apt-get install libgmp-dev")
endif()

# 查找 MPFR
find_library(MPFR_LIBRARIES NAMES mpfr)
find_path(MPFR_INCLUDE_DIRS NAMES mpfr.h)
if(NOT MPFR_LIBRARIES OR NOT MPFR_INCLUDE_DIRS)
    message(FATAL_ERROR "MPFR not found! Install: sudo apt-get install libmpfr-dev")
endif()

# 查找 CGAL
find_package(CGAL REQUIRED)
if(NOT CGAL_FOUND)
    message(FATAL_ERROR "CGAL not found! Install: sudo apt-get install libcgal-dev")
endif()

# 查找 Boost
find_package(Boost REQUIRED COMPONENTS filesystem system)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost not found! Install: sudo apt-get install libboost-all-dev")
endif()

# 查找 Eigen3 - 修复目标名称问题
find_package(Eigen3 REQUIRED)
if(NOT Eigen3_FOUND)
    message(FATAL_ERROR "Eigen3 not found! Install: sudo apt-get install libeigen3-dev")
endif()

# 检查 Eigen3 目标是否存在，如果不存在则创建
if(NOT TARGET Eigen3::Eigen)
    message(STATUS "Eigen3::Eigen target not found, creating interface target")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIRS}"
    )
endif()

# 查找 Spatialindex
find_library(Spatialindex_LIBRARIES NAMES spatialindex spatialindex_c)
find_path(Spatialindex_INCLUDE_DIRS NAMES spatialindex/SpatialIndex.h)
if(NOT Spatialindex_LIBRARIES OR NOT Spatialindex_INCLUDE_DIRS)
    message(WARNING "Spatialindex not found - some features may be disabled")
endif()

# 查找 Qhull
find_library(QHULLCPP_LIBRARIES NAMES qhullcpp)
find_library(QHULL_R_LIBRARIES NAMES qhull_r)
find_path(QHULL_INCLUDE_DIRS NAMES libqhullcpp/Qhull.h)
if(NOT QHULLCPP_LIBRARIES OR NOT QHULL_R_LIBRARIES OR NOT QHULL_INCLUDE_DIRS)
    message(FATAL_ERROR "Qhull not found! Install: sudo apt-get install libqhull-dev")
endif()

# 添加可执行文件
add_executable(POP src/main.cpp src/CH.hpp src/common.hpp src/Graph.hpp src/TD.hpp src/PL.hpp)

# 包含目录
target_include_directories(POP PRIVATE
        ${GMP_INCLUDE_DIRS}
        ${MPFR_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${Spatialindex_INCLUDE_DIRS}
        ${QHULL_INCLUDE_DIRS}
)

# 链接库
target_link_libraries(POP PRIVATE
        CGAL::CGAL
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES}
        Boost::filesystem
        Boost::system
        Eigen3::Eigen
        ${Spatialindex_LIBRARIES}
        ${QHULLCPP_LIBRARIES}
        ${QHULL_R_LIBRARIES}
        OpenMP::OpenMP_CXX
)

# 优化设置
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(POP PRIVATE -O2 -march=native)
endif()

message(STATUS "=== Configuration complete ===")