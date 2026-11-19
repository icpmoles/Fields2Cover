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

namespace f2c::rp {

/// First solution built from the swaths instead of by an or-tools heuristic.
///   The swaths are grouped by the line they were cut from: starting from the
///   line at the border of the field, every segment of a line is covered in
///   order, jumping over the obstacles that split it, then the route moves to
///   the next parallel line and runs it in the opposite direction (snake).
///   The route is given to the solver as the first solution, then improved
///   by the local search. The value is past the last or-tools strategy
///   (LOCAL_CHEAPEST_COST_INSERTION = 16), so it is never passed to the solver.
constexpr int kSwathSnakeStrategy = 17;

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
  /// @return Route that covers all the swaths
  virtual F2CRoute genRoute(
      const F2CCells& cells, const F2CSwathsByCells& swaths_by_cells,
      bool show_log = false, double d_tol = 1e-4, bool redirect_swaths = true,
      long int time_limit_seconds = 1, bool search_for_optimum = false);

  /// Set the start and the end of the route.
  void setStartAndEndPoint(const F2CPoint& p);

  /// Set the first solution strategy of the solver: either the value of an
  ///   or-tools operations_research::FirstSolutionStrategy::Value
  ///   (default: AUTOMATIC = 15) or kSwathSnakeStrategy.
  void setFirstSolutionStrategy(int strategy);
  int getFirstSolutionStrategy() const;

  /// Build the first solution of kSwathSnakeStrategy: the nodes of the
  ///   coverage graph in visiting order, depot excluded (same format as the
  ///   route returned by computeBestRoute).
  ///   Throws if the swaths and the coverage graph don't match.
  /// @param swaths_by_cells Swaths used to build `cov_graph`
  /// @param cov_graph Graph from createCoverageGraph (with redirect_swaths)
  std::vector<long long int> genSwathSnakeRoute(
      const F2CSwathsByCells& swaths_by_cells,
      const F2CGraph2D& cov_graph) const;

  /// Set the first solution used by computeBestRoute when the strategy is
  ///   kSwathSnakeStrategy. genRoute sets it by itself.
  void setInitialRoute(const std::vector<long long int>& route);
  const std::vector<long long int>& getInitialRoute() const;

  /// Create graph to compute the shortest path between two points
  ///   in the headlands.
  ///
  /// @param cells Headland swath rings used to travel through the headlands
  /// @param swaths_by_cells Swaths to be covered.
  /// @param d_tol Tolerance distance to consider if two points are the same.
  virtual F2CGraph2D createShortestGraph(
      const F2CCells& cells, const F2CSwathsByCells& swaths_by_cells,
      double d_tol) const;

  /// Create graph to compute the cost of covering the swaths in a given order.
  ///
  /// @param cells Headland swath rings used to travel through the headlands
  /// @param swaths_by_cells Swaths to be covered.
  /// @param shortest_graph Graph to compute the shortest path
  ///          between two nodes.
  /// @param d_tol Tolerance distance to consider if two points are the same.
  virtual F2CGraph2D createCoverageGraph(
      const F2CCells& cells, const F2CSwathsByCells& swaths_by_cells,
      F2CGraph2D& shortest_graph,
      double d_tol, bool redirect_swaths = true) const;


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
      const F2CGraph2D& cov_graph, bool show_log, long int time_limit_seconds,
      bool use_guided_local_search = true) const;

  /// Tranform index of points to an actual Route.
  virtual F2CRoute transformSolutionToRoute(
      const std::vector<long long int>& route_ids,
      const F2CSwathsByCells& swaths_by_cells,
      const F2CGraph2D& coverage_graph,
      F2CGraph2D& shortest_graph) const;

 protected:
  std::optional<F2CPoint> r_start_end;
  int first_solution_strategy_ {15};  // or-tools AUTOMATIC
  std::vector<long long int> initial_route_;
};



}  // namespace f2c::rp

#endif  // FIELDS2COVER_ROUTE_PLANNING_ROUTE_PLANNING_BASE_H_

