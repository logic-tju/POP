# POP: Preference-Optimal Path Querying over Multi-Cost Graphs

[![PVLDB 2025](https://img.shields.io/badge/Publication-PVLDB_2025-blue)](https://www.pvldb.org/)
[![C++](https://img.shields.io/badge/C++-17-blue)](https://isocpp.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

This repository contains the official implementation of the paper:

**"One-for-All: Efficient Preference-Optimal Path Querying over Multi-Cost Graphs"**  
*Xinyu Xie, Yajun Yang, Jeffrey Xu Yu, Hong Gao*  
PVLDB, 14(1), 2025

## 📖 Abstract

Preference-optimal path query aims to find a path with the minimum score under a user-specific linear function in multi-cost graphs. Existing methods on single-cost graphs cannot handle arbitrary linear preference functions. We propose a novel "one-for-all" index (POP-index) that enables efficient processing of preference-optimal path queries under any linear function. Our approach combines graph partitioning, tree-based contraction, and geometric dimensionality-reduction to achieve orders of magnitude improvement over state-of-the-art methods.

## 🚀 Quick Start

### Prerequisites

- **OS**: Linux (Ubuntu 20.04+ recommended)
- **Compiler**: GCC 9.0+ or Clang 10.0+ with C++17 support
- **CMake**: 3.20+

### Dependencies Installation

```bash
# Install system dependencies (Ubuntu/Debian)
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

