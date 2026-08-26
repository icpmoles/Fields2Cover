//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#pragma once
#ifndef FIELDS2COVER_ROUTE_PLANNING_ROUTE_PLANNING_BASE_H_
#define FIELDS2COVER_ROUTE_PLANNING_ROUTE_PLANNING_BASE_H_

#include <vector>
#include <map>
#include <optional>
#include <limits>
#include <utility>
#include "fields2cover/types.h"
#include "fields2cover/objectives/rp_obj/rp_objective.h"
#include "fields2cover/objectives/rp_obj/direct_dist_path_obj.h"
#include "fields2cover/route_planning/single_cell_swaths_order_base.h"
#include <ortools/constraint_solver/routing_enums.pb.h>


namespace f2c::rp {

// taken from https://github.com/google/or-tools/blob/v9.9/ortools/constraint_solver/routing_enums.proto#L25-L110
inline std::vector<std::string> FirstSolutionStrategy = {
  "UNSET",
  "GLOBAL_CHEAPEST_ARC",
  "LOCAL_CHEAPEST_ARC",
  "PATH_CHEAPEST_ARC",
  "PATH_MOST_CONSTRAINED_ARC",
  "EVALUATOR_STRATEGY", // 5
  "ALL_UNPERFORMED",
  "BEST_INSERTION",
  "PARALLEL_CHEAPEST_INSERTION",
  "LOCAL_CHEAPEST_INSERTION",
  "SAVINGS", // 10
  "SWEEP",
  "FIRST_UNBOUND_MIN_VALUE",
  "CHRISTOFIDES",
  "SEQUENTIAL_CHEAPEST_INSERTION",
  "AUTOMATIC", // 15
  "LOCAL_CHEAPEST_COST_INSERTION",
};

// Shared read-only problem data across threads
struct CVrpData {
  F2CGraph2D Cov_Graph;
  std::vector<std::vector<bool>> constraints;
  long int time_limit_seconds;
  bool constrained = false;
  bool use_guided_local_search = true;
  bool show_log = true;

};

// Search configuration for a single thread worker
struct SearchConfig {
  operations_research::FirstSolutionStrategy::Value strategy;
  // int32_t random_seed;
};

// Result structure returned from worker threads
struct ThreadResult {
  int64_t cost = -1;
  bool success = false;
  operations_research::FirstSolutionStrategy::Value strategy;
  // int32_t random_seed;
  std::vector<long long int> route;
  F2CRoute p_route;
  int64_t path_length = -1;
  int collisions = -1;
};

class RoutePlannerBase {
 public:
  /// Generate route to cover the swaths on a field.
  ///   If two consecutive swaths are far away,
  ///   the route connects both through the headland.
  ///
  /// @param cells Headland swath rings used to travel through the headlands
  /// @param swaths_by_cells Swaths to be covered.
  /// @param show_log Show log from the optimizer
  /// @param d_tol Tolerance distance to consider if two points are the same.
  /// @param redirect_swaths Whether to allow redirecting swaths
  /// @param time_limit_seconds Maximum time to spend on optimization
  /// @param search_for_optimum If true, uses guided local search which may take longer
  ///        but can find more optimal solutions. If false, uses automatic search which is faster
  ///        but may find less optimal solutions.
  /// @param dist_exponent
  /// @param use_visibility
  /// @param visibility_factor
  /// @param visibility_use_crossing
  /// @param prefer_inter_rings_crossing
  /// @param free_space_planner
  /// @param constrained
  /// @param parallel
  /// @return Route that covers all the swaths
  virtual F2CRoute genRoute(const F2CCells& cells,
       const F2CSwathsByCells& swaths_by_cells,
       std::vector<ThreadResult> &output,
      ThreadResult &best,
       bool show_log = false,
       double d_tol = 1e-4,
       bool redirect_swaths = true,
       long int time_limit_seconds = 1,
       bool search_for_optimum = false,
       float dist_exponent = 2,
       bool use_visibility = false,
       float visibility_factor = 2,
       bool visibility_use_crossing = false,
       bool prefer_inter_rings_crossing = false,
       bool free_space_planner = true,
       bool constrained = false,
       bool parallel = false);

  /// Set the start and the end of the route.
  void setStartAndEndPoint(const F2CPoint& p);

  /// Create graph to compute the shortest path between two points
  ///   in the headlands.
  ///
  /// @param cells Headland swath rings used to travel through the headlands
  /// @param swaths_by_cells Swaths to be covered.
  /// @param d_tol Tolerance distance to consider if two points are the same.
  /// @param free_space_planner
  virtual F2CGraph2D createShortestGraph(const F2CCells& cells,
      const F2CSwathsByCells& swaths_by_cells,
      double d_tol,
      bool free_space_planner) const;

  /// Create graph to compute the cost of covering the swaths in a given order.
  ///
  /// @param cells Headland swath rings used to travel through the headlands
  /// @param swaths_by_cells Swaths to be covered.
  /// @param shortest_graph Graph to compute the shortest path
  ///          between two nodes.
  /// @param d_tol Tolerance distance to consider if two points are the same.
  /// @param dist_exponent
  /// @param use_visibility
  /// @param visibility_factor
  /// @param visibility_use_crossing
  /// @param prefer_crossings
  /// @param free_space_planner
  virtual F2CGraph2D createCoverageGraph(const F2CCells& cells,
      const F2CSwathsByCells& swaths_by_cells,
      F2CGraph2D& shortest_graph,
      double d_tol,
      bool redirect_swaths = true,
      float dist_exponent = 2,
      bool use_visibility = false,
      double visibility_factor = 2,
      bool visibility_use_crossing = false,
      bool prefer_crossings = false,
      bool free_space_planner = false,
      bool constrained = true) const;


  virtual ~RoutePlannerBase() = default;

  /// Use the optimizer to generate the index of the points of the best
  ///   coverage route.
  /// @param cov_graph Graph representing the coverage problem
  /// @param show_log Whether to show optimization logs
  /// @param time_limit_seconds Maximum time to spend on optimization
  /// @param use_guided_local_search If true, uses guided local search which may take longer
  ///        but can find more optimal solutions. If false, uses automatic search which is faster
  ///        but may find less optimal solutions.
  virtual std::vector<long long int> computeBestRoute(
      const F2CGraph2D& cov_graph, bool show_log,
      long int time_limit_seconds,
      bool use_guided_local_search = true,
      bool constrained = true) const;

  /// Tranform index of points to an actual Route.
  virtual F2CRoute transformSolutionToRoute(
      const std::vector<long long int>& route_ids,
      const F2CSwathsByCells& swaths_by_cells,
      const F2CGraph2D& coverage_graph,
      F2CGraph2D& shortest_graph) const;

  virtual std::vector<long long int> computeBestRouteParallel(F2CGraph2D& cov_graph,
      bool show_log,
      long int time_limit_seconds,
      std::vector<ThreadResult> &output,
      ThreadResult &best,
      bool use_guided_local_search = true,
      bool constrained = true,
      uint cores = 4) const;

  ThreadResult RunSingleInitialization(CVrpData& data, SearchConfig config) const;

 protected:
  std::optional<F2CPoint> r_start_end;

 private:
  double pessimistic_traversal_ = 0.0;

  //
  std::vector<operations_research::FirstSolutionStrategy::Value> possibleStrategies = {
    // operations_research::FirstSolutionStrategy::UNSET,
    operations_research::FirstSolutionStrategy::GLOBAL_CHEAPEST_ARC,
    operations_research::FirstSolutionStrategy::LOCAL_CHEAPEST_ARC,
    operations_research::FirstSolutionStrategy::PATH_CHEAPEST_ARC,
    operations_research::FirstSolutionStrategy::PATH_MOST_CONSTRAINED_ARC,
    operations_research::FirstSolutionStrategy::EVALUATOR_STRATEGY, // 5
    operations_research::FirstSolutionStrategy::ALL_UNPERFORMED,
    operations_research::FirstSolutionStrategy::BEST_INSERTION,
    operations_research::FirstSolutionStrategy::PARALLEL_CHEAPEST_INSERTION,
    operations_research::FirstSolutionStrategy::LOCAL_CHEAPEST_INSERTION,
    operations_research::FirstSolutionStrategy::SAVINGS, // 10
    operations_research::FirstSolutionStrategy::SWEEP,
    operations_research::FirstSolutionStrategy::FIRST_UNBOUND_MIN_VALUE,
    operations_research::FirstSolutionStrategy::CHRISTOFIDES,
    operations_research::FirstSolutionStrategy::SEQUENTIAL_CHEAPEST_INSERTION,
    operations_research::FirstSolutionStrategy::AUTOMATIC, // 15
    operations_research::FirstSolutionStrategy::LOCAL_CHEAPEST_COST_INSERTION,
  };
};



}  // namespace f2c::rp

#endif  // FIELDS2COVER_ROUTE_PLANNING_ROUTE_PLANNING_BASE_H_

