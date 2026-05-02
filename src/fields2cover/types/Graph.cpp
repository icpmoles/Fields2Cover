//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <numeric>
#include <algorithm>
#include "fields2cover/types/Graph.h"

#include <iostream>

namespace f2c::types {

Graph& Graph::addDirectedEdge(size_t from, size_t to, int64_t cost) {
  this->edges_[from][to] =  cost;
  if (this->shortest_paths_.size() > 0) {
    this->shortest_paths_.clear();
  }
  return *this;
}

Graph& Graph::addEdge(size_t i, size_t j, int64_t cost) {
  return this->addDirectedEdge(i, j, cost).addDirectedEdge(j, i, cost);
}


Graph& Graph::removeDirectedEdge(size_t from, size_t to) {
  if (this->edges_.count(from) == 0 || this->edges_.at(from).count(to) == 0) {
    return *this;
  }
  this->edges_.at(from).erase(to);
  return *this;
}
Graph& Graph::removeEdge(size_t i, size_t j) {
  return this->removeDirectedEdge(i, j).removeDirectedEdge(j, i);
}

size_t Graph::numNodes() const {
  return this->edges_.size();
}

size_t Graph::numEdges() const {
  return std::accumulate(this->edges_.begin(), this->edges_.end(), 0,
      [](int sum, const auto& e) {return sum + e.second.size();});
}

map_to_map_to_int Graph::getEdges() const {
  return this->edges_;
}

std::vector<size_t> Graph::getEdgesFrom(size_t s) const {
  if (this->edges_.count(s) == 0) {
    return {};
  }
  std::vector<size_t> connections;
  for (auto&& i : this->edges_.at(s)) {
    connections.emplace_back(i.first);
  }
  return connections;
}


// this is the only place where we can get 1E15 as output
// this is only used in the vrp callback
int64_t Graph::getCostFromEdge(size_t from, size_t to, int64_t INF) const {
  if (this->edges_.count(from) == 0 || this->edges_.at(from).count(to) == 0) {
    return INF;
  }
  return this->edges_.at(from).at(to);
}


/**
 * returns list of paths (defined as list of node's idx) between from and src (uses DFS)
 *
 * Nb: never used
*/
std::vector<std::vector<size_t>> Graph::allPathsBetween(
    size_t from, size_t to) const {
  // list of visited nodes
  std::vector<bool> visited(this->numNodes(), false);
  // list of routes
  std::vector<std::vector<size_t>> routes(1);
  int route_index = 0;
  this->DFS(from, to, routes, visited, route_index);

  // erases wrong paths:
  // paths with zero-length
  // paths that didn't reach destination
  routes.erase(std::remove_if(routes.begin(), routes.end(),
      [&to] (const std::vector<size_t> x) {
          return (x.size() < 1 || x.back() != to);
      }), routes.end());
  return routes;
}


// depth first search
// nb: never used
void Graph::DFS(
    size_t from, size_t to,
    std::vector<std::vector<size_t>>& routes,
    std::vector<bool>& visited,
    int& route_index) const {

  // copies the active path in a new path, this increases the size by 1
  routes.emplace_back(routes[route_index]);
  // in the new path: last element is the previous node in the exploration tree
  routes.back().emplace_back(from);
  // gets the increased size of the routes vector
  int i_route = routes.size() - 1;

  // if destination not reached
  if (from != to) {
    visited[from] = true;
    // iterate over all the edges
    for (auto&& i : this->getEdgesFrom(from)) {
      // if not visited already
      if (!visited[i]) {
        // recursive call
        this->DFS(i, to, routes, visited, i_route);
      }
    }
    // mark visited
    visited[from] = false;
  }
}

// computes the optimal distance matrices: INF: 1073741824 = 2^30
std::vector<std::vector<pair_vec_size__int>>
    Graph::shortestPathsAndCosts(int64_t INF) {
  const size_t N = this->numNodes();
  std::vector<std::vector<int64_t>> dist(N, std::vector<int64_t>(N, INF));
  std::vector<std::vector<int64_t>> next(N, std::vector<int64_t>(N, -1));


  std::cout << "--> --> shortestPathsAndCosts Init shortestPathsAndCosts, N = "<< N << std::endl;
  // Initialize distances and paths for direct connections
  // NB: dist matrix: symmetric if undirected

  if (new_gen_) {
    // for every src
    for (const auto& src : this->edges_) {
      // for every destination
      const size_t src_idx = src.first;
      for (const auto& dest : src.second) {
        // get idx of dest
        const size_t dest_idx = dest.first;
        // fills dist[][]
        dist[src_idx][dest_idx] = dest.second;
        // fills next[][]
        next[src_idx][dest_idx] = dest_idx;
      }
      // distance cost to itself = 0
      dist[src_idx][src_idx] = 0;
      next[src_idx][src_idx] = src_idx;
    }
  } else {
    // Initialize distances and paths for direct connections
    for (const auto& src : this->edges_) {
      for (const auto& dest : src.second) {
        dist[src.first][dest.first] = dest.second;
        next[src.first][dest.first] = dest.first;
      }
      dist[src.first][src.first] = 0;
    }
  }

  std::cout << "--> --> shortestPathsAndCosts Floyd-Warshall" << std::endl;
  // Floyd-Warshall
  for (size_t k = 0; k < N; ++k) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        if (dist[i][k] < INF && \
            dist[k][j] < INF && \
            dist[i][j] > dist[i][k] + dist[k][j]) {
          dist[i][j] = dist[i][k] + dist[k][j];
          next[i][j] = next[i][k];
        }
      }
    }
  }

  std::cout << "--> --> shortestPathsAndCosts Reconstruct paths" << std::endl;

  // Reconstruct paths
  std::vector<std::vector<pair_vec_size__int>>
      paths(N, std::vector<pair_vec_size__int>(N));

  if (new_gen_) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        if (i != j && next[i][j] != -1) {
          bool successfull_path = true;
          std::vector<size_t> path = {i};
          size_t current = i;
          while (current != j) {
            if (path.size() > N) {successfull_path=false; break;}
            current = next[current][j];
            path.push_back(current);
          }
          if (successfull_path) {
            paths[i][j] = std::make_pair(path, dist[i][j]);
          } else {
            paths[i][j].second = INF;
          }
        } else if (i != j && next[i][j] == -1) {
          paths[i][j].second = INF;
          // std::cout << i << " -> " << j << " not found" << std::endl;
        }
      }
    }
  } else {
     for (size_t i = 0; i < N; ++i) {
    //for (auto i: swath_extremities_nodes_) {
      // if (std::find(swath_extremities_nodes_.begin(), swath_extremities_nodes_.end(), i) == swath_extremities_nodes_.end()) { // value not found
      //   continue;
      // }
       for (size_t j = 0; j < N; ++j) {
   // for (auto j: swath_extremities_nodes_) {
        // if (std::find(swath_extremities_nodes_.begin(), swath_extremities_nodes_.end(), j) == swath_extremities_nodes_.end()) { // value not found
        //   continue;
        // }
        if (i != j && next[i][j] != -1) {
          std::vector<size_t> path = {i};
          size_t current = i;
          while (current != j) {
            current = next[current][j];
            path.push_back(current);
          }
          paths[i][j] = std::make_pair(path, dist[i][j]);
        } else if (i != j && next[i][j] == -1) {
          paths[i][j].second = INF;
        }
      }
    }
  }
  std::cout << "--> --> shortestPathsAndCosts move" << std::endl;

  this->shortest_paths_ = std::move(paths);
  return this->shortest_paths_;
}

std::vector<size_t> Graph::shortestPath(size_t from, size_t to, int64_t INF) {
  if (this->numNodes() > 0 && this->shortest_paths_.size() == 0) {
    std::cout << "Calculating shortest_paths_" << std::endl;
    this->shortestPathsAndCosts(INF);
  }
  return this->shortest_paths_[from][to].first;
}

int64_t Graph::shortestPathCost(size_t from, size_t to, int64_t INF) {
  if (this->numNodes() > 0 && this->shortest_paths_.size() == 0) {
    std::cout << "Calculating shortest_paths_" << std::endl;
    // initialize if not done already
    this->shortestPathsAndCosts(INF);
  }
  return this->shortest_paths_[from][to].second;
}

}  // namespace f2c::types


