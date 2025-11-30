# POP: Preference-Optimal Path Querying over Multi-Cost Graphs

[![PVLDB 2025](https://img.shields.io/badge/Publication-PVLDB_2025-blue)](https://www.pvldb.org/)
[![C++](https://img.shields.io/badge/C++-17-blue)](https://isocpp.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

This repository contains the official implementation of the paper:

"One-for-All: Efficient Preference-Optimal Path Querying over Multi-Cost Graphs"
Xinyu Xie, Yajun Yang, Jeffrey Xu Yu, Hong Gao
PVLDB, 14(1), 2025

=====================================================================

Abstract
--------

Preference-optimal path query aims to find a path with the minimum score under a user-specific linear function in multi-cost graphs. Existing methods on single-cost graphs cannot handle arbitrary linear preference functions. We propose a novel one-for-all index (POP-index) that enables efficient processing of preference-optimal path queries under any linear function. Our approach combines graph partitioning, tree-based contraction, and geometric dimensionality-reduction to achieve orders of magnitude improvement over state-of-the-art methods.

=====================================================================

Quick Start
-----------

Prerequisites:
- OS: Linux (Ubuntu 20.04+ recommended)
- Compiler: GCC 9.0+ or Clang 10.0+ (supporting C++17)
- CMake: 3.20+

Install dependencies (Ubuntu/Debian):

    sudo apt-get update
    sudo apt-get install -y \
        build-essential \
        cmake \
        libboost-all-dev \
        libcgal-dev \
        libgmp-dev \
        libmpfr-dev \
        libeigen3-dev \
        libqhull-dev \
        libspatialindex-dev

=====================================================================

Building POP
------------

    mkdir build && cd build
    cmake ..
    cmake --build .

After building, an executable named "POP" will appear inside build/.

=====================================================================

Program Usage
-------------

The POP program requires multiple input lines describing:
1. Dataset name (short code)
2. Partition node upper bound
3. Root directory of the project
4. Full dataset identifier
5. Fixed dimension count (always 2)
6. Number of additional (random) dimensions
7. Weight vector (normalized, length = 2 + random_dim)

Example input:

    BAY
    1000
    /home/xxy/POP
    BAY-0.32M-0.80M
    2
    1
    0.1
    0.1
    0.1

Explanation of each line:

1. Dataset code:
   Available options: BAY, COL, FLA, NW, NY

2. Partition node upper bound:
   Examples: 1000, 2000, 4000, 8000
   Controls the node limit per partition block.

3. Root directory:
   Path to the project directory containing graph data.

4. Dataset full name:
   Must be one of:
     - BAY-0.32M-0.80M
     - COL-0.44M-1.06M
     - FLA-1.07M-2.71M
     - NW-1.21M-2.84M
     - NY-0.26M-0.73M

5. Fixed dimension count:
   Always "2". All datasets come with exactly two real cost dimensions.

6. Additional dimension count:
   Set as (desired_total_dimension - 2).
   Example:
     - If total dimension = 3 → input "1"
     - If total dimension = 5 → input "3"

7. Weight vector:
   Number of weights = fixed_dim (2) + random_dim
   Requirements:
     - each weight >= 0
     - sum of weights = 1

Example:
If fixed_dim=2 and random_dim=1 → total_dim=3 → require 3 weights:

    0.1
    0.1
    0.8

=====================================================================

IMPORTANT NOTE (Very Important!)
--------------------------------

Before running the program, **you must ensure that the following constant matches the total dimension**:

File: `src/common.hpp`  
Line: **58**  
Code:

    const int D = 3;

This `D` must equal:

    total_dimension = fixed_dimension (2) + random_dimension

Examples:
- If fixed_dim = 2 and random_dim = 1 → total_dim = 3 → set `D = 3`
- If fixed_dim = 2 and random_dim = 3 → total_dim = 5 → set `D = 5`

If this value does not match, the algorithm may produce incorrect results.

=====================================================================

Dataset Format
--------------

Each dataset is a multi-cost road network.
Each edge contains exactly two real costs (fixed-cost dimensions).
POP can extend the cost vector by generating additional random cost dimensions.

=====================================================================

Citation
--------

If you use this repository, please cite:

Xinyu Xie, Yajun Yang, Jeffrey Xu Yu, Hong Gao.
"One-for-All: Efficient Preference-Optimal Path Querying over Multi-Cost Graphs."
PVLDB 14(1), 2025.

=====================================================================

License
-------

Released under the MIT License.

=====================================================================
