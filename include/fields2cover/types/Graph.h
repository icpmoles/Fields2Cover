//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#pragma once
#ifndef FIELDS2COVER_TYPES_GRAPH_H_
#define FIELDS2COVER_TYPES_GRAPH_H_

#include <cstdint>
#include <vector>
#include <set>
#include <functional>
#include <utility>
#include <unordered_map>

namespace f2c::types {

/**
 * Map:
 *
 * map2map2int[src_idx] = map2int from "src"
 *
 * map2map2int[src_idx][dst_idx] = cost of traversal from "src" to "dest"
 *
 */
typedef std::unordered_map<size_t, std::unordered_map<size_t, int64_t>>
  map_to_map_to_int;

  /**
   * Pair:
   *
   * .first = list of edges' indexes of the path
   *
   * .second = cost of the path
   */
typedef std::pair<std::vector<size_t>, int64_t> pair_vec_size__int;

class Graph {
 public:
  Graph& addDirectedEdge(size_t from, size_t to, int64_t cost);
  Graph& addEdge(size_t i, size_t j, int64_t cost);
  Graph& removeDirectedEdge(size_t from, size_t to);
  Graph& removeEdge(size_t i, size_t j);

  size_t numNodes() const;
  size_t numEdges() const;
  map_to_map_to_int getEdges() const;

  std::vector<size_t> getEdgesFrom(size_t s) const;

  int64_t getCostFromEdge(size_t from, size_t to,
      int64_t INF = 1e15) const;

  std::vector<std::vector<size_t>> allPathsBetween(
      size_t from, size_t to) const;

  std::vector<std::vector<pair_vec_size__int>>
    shortestPathsAndCosts(int64_t INF = 1e15);

  std::vector<size_t> shortestPath(size_t from, size_t to,
        int64_t INF = 1e15);

  int64_t shortestPathCost(size_t from, size_t to,
        int64_t INF = 1e15);

 void setNewGen(bool newGen) {
   this->new_gen_ = newGen;
 }
 protected:
  void DFS(size_t from, size_t to,
    std::vector<std::vector<size_t>>& routes,
    std::vector<bool>& visited,
    int& route_index) const;

 protected:
  map_to_map_to_int edges_;
  /**
   * sp_[src_idx][dest_idx] saves the shortest path from src to dest
   */
  std::vector<std::vector<pair_vec_size__int>> shortest_paths_;
  bool new_gen_=false;


  // nodes that represent the middle of the swaths
  std::set<size_t> swath_extremities_nodes_;
};

}  // namespace f2c::types

#endif  // FIELDS2COVER_TYPES_GRAPH_H_
