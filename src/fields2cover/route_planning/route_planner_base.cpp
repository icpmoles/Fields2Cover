//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <ortools/constraint_solver/routing.h>
#include <ortools/constraint_solver/routing_enums.pb.h>
#include <ortools/constraint_solver/routing_index_manager.h>
#include <ortools/constraint_solver/routing_parameters.h>
#include <math.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <limits>
#include "fields2cover/route_planning/route_planner_base.h"

namespace f2c::rp {

namespace ortools = operations_research;

F2CRoute RoutePlannerBase::genRoute(
    const F2CCells& cells, const F2CSwathsByCells& swaths,
    bool show_log, double d_tol, bool redirect_swaths,
    long int time_limit_seconds, bool search_for_optimum) {
  F2CGraph2D shortest_graph = createShortestGraph(cells, swaths, d_tol);

  F2CGraph2D cov_graph = createCoverageGraph(
      cells, swaths, shortest_graph, d_tol, redirect_swaths);

  if (first_solution_strategy_ == kSwathSnakeStrategy) {
    initial_route_ = genSwathSnakeRoute(swaths, cov_graph);
  }

  std::vector<long long int> v_route = computeBestRoute(
      cov_graph, show_log, time_limit_seconds, search_for_optimum);
  return transformSolutionToRoute(
      v_route, swaths, cov_graph, shortest_graph);
}

void RoutePlannerBase::setStartAndEndPoint(const F2CPoint& p) {
  this->r_start_end = p;
}

F2CGraph2D RoutePlannerBase::createShortestGraph(
    const F2CCells& cells, const F2CSwathsByCells& swaths_by_cells,
    double d_tol) const {
  F2CGraph2D g;
  // Connect each swath end to the border of the cell it belongs to. An end
  // outside every cell is left out: a line to the nearest border would cross
  // whatever lies between.
  auto addBorderEdge = [&g, &cells, d_tol] (const F2CPoint& p) {
    F2CCell cell = cells.getCellWherePoint(p, d_tol);
    if (!cell.isEmpty()) {
      g.addEdge(p, cell.closestPointOnBorderTo(p));
    }
  };
  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      addBorderEdge(s.startPoint());
      addBorderEdge(s.endPoint());
    }
  }

  // Add points in the border
  for (auto&& cell : cells) {
    for (auto&& ring : cell) {
      for (size_t i = 0; i < ring.size()-1; ++i) {
        g.addEdge(ring.getGeometry(i), ring.getGeometry(i+1));
      }
    }
  }

  // Add start and end point if they exists. Unlike a swath end, this point
  // is allowed to sit outside the field, so it always gets its border edge.
  if (this->r_start_end) {
    g.addEdge(*r_start_end, cells.closestPointOnBorderTo(*r_start_end));
  }

  std::vector<F2CPoint> nodes = g.getNodes();

  // Connect nodes that are near other edges.
  for (int i = 0; i < 2; ++i) {
    auto edges = g.getEdges();
    for (auto&& edge : edges) {
      size_t from = edge.first;
      for (auto&& e : edge.second) {
        size_t to = e.first;

        F2CLineString border {g.indexToNode(from), g.indexToNode(to)};
        for (auto&& n : nodes) {
          if (n != border[0] && n != border[1] && n.distance(border) < d_tol) {
            g.addEdge(border[0], n);
            g.addEdge(n, border[1]);
            g.removeEdge(border[0], border[1]);
          }
        }
      }
    }
  }
  return g;
}


F2CGraph2D RoutePlannerBase::createCoverageGraph(
    const F2CCells& cells, const F2CSwathsByCells& swaths_by_cells,
    F2CGraph2D& shortest_graph,
    double d_tol, bool redirect_swaths) const {
  F2CGraph2D g;
  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      F2CPoint mid_p {(s.startPoint() + s.endPoint()) * 0.5};
      if (redirect_swaths) {
        g.addEdge(s.startPoint(), mid_p, 0);
        g.addEdge(s.endPoint(), mid_p, 0);
      } else {
        g.addDirectedEdge(s.startPoint(), mid_p, 0);
        g.addDirectedEdge(mid_p, s.endPoint(), 0);
      }
    }
  }

  for (const auto& swaths1 : swaths_by_cells) {
    for (const auto& s1 : swaths1) {
      auto s1_s = s1.startPoint();
      auto s1_e = s1.endPoint();
      for (const auto& swaths2 : swaths_by_cells) {
        for (const auto& s2 : swaths2) {
          auto s2_s = s2.startPoint();
          auto s2_e = s2.endPoint();
          if (redirect_swaths) {
            g.addEdge(s1_s, s2_s, shortest_graph);
            g.addEdge(s1_e, s2_e, shortest_graph);
            g.addEdge(s1_s, s2_e, shortest_graph);
            g.addEdge(s1_e, s2_s, shortest_graph);
          } else {
            g.addDirectedEdge(s1_e, s2_s, shortest_graph);
          }
        }
      }
    }
  }

  F2CPoint deposit(-1e8, -1e8);  // Arbitrary point
  if (this->r_start_end) {
    deposit = *this->r_start_end;
  }

  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      if (this->r_start_end) {
        g.addEdge(s.startPoint(), deposit, shortest_graph);
        g.addEdge(s.endPoint(), deposit, shortest_graph);
      } else {
        g.addEdge(s.startPoint(), deposit, 0);
        g.addEdge(s.endPoint(), deposit, 0);
      }
    }
  }
  return g;
}

void RoutePlannerBase::setFirstSolutionStrategy(int strategy) {
  if (strategy != kSwathSnakeStrategy &&
      !ortools::FirstSolutionStrategy::Value_IsValid(strategy)) {
    throw std::invalid_argument(
        "Unknown first solution strategy: " + std::to_string(strategy));
  }
  this->first_solution_strategy_ = strategy;
}

int RoutePlannerBase::getFirstSolutionStrategy() const {
  return this->first_solution_strategy_;
}

void RoutePlannerBase::setInitialRoute(
    const std::vector<long long int>& route) {
  this->initial_route_ = route;
}

const std::vector<long long int>& RoutePlannerBase::getInitialRoute() const {
  return this->initial_route_;
}

std::vector<long long int> RoutePlannerBase::genSwathSnakeRoute(
    const F2CSwathsByCells& swaths_by_cells,
    const F2CGraph2D& cov_graph) const {
  const size_t n_swaths = swaths_by_cells.sizeTotal();
  if (n_swaths == 0) {
    return {};
  }

  // All the swaths are parallel: take the direction of the longest one,
  // pointing to +x (or +y if vertical), so the result does not depend on how
  // each swath is oriented
  size_t longest = 0;
  for (size_t j = 1; j < n_swaths; ++j) {
    if (swaths_by_cells.getSwath(j).length() >
        swaths_by_cells.getSwath(longest).length()) {
      longest = j;
    }
  }
  const F2CSwath& ref = swaths_by_cells.getSwath(longest);
  double dx = ref.endPoint().getX() - ref.startPoint().getX();
  double dy = ref.endPoint().getY() - ref.startPoint().getY();
  const double norm = std::hypot(dx, dy);
  if (norm <= 0.0) {
    throw std::runtime_error("SWATH_SNAKE: the swaths have no length");
  }
  dx /= norm;
  dy /= norm;
  if (dx < -1e-9 || (std::abs(dx) <= 1e-9 && dy < 0.0)) {
    dx = -dx;
    dy = -dy;
  }
  // along: position on the swath line; across: which parallel line
  auto along = [&](const F2CPoint& p) {
    return p.getX() * dx + p.getY() * dy;};
  auto across = [&](const F2CPoint& p) {
    return -p.getX() * dy + p.getY() * dx;};
  // nodeToIndex throws std::out_of_range on unknown points
  auto node_id = [&cov_graph](const F2CPoint& p) -> long long int {
    try {
      return static_cast<long long int>(cov_graph.nodeToIndex(p));
    } catch (const std::out_of_range&) {
      throw std::runtime_error(
          "SWATH_SNAKE: the swaths don't match the coverage graph");
    }
  };

  // Only plain values: sorting F2CPoints by swapping them is not safe
  struct Segment {
    double across;
    double along;        // lowest along coordinate of the two ends
    long long int low;   // node of the end with the lowest along coordinate
    long long int high;  // node of the other end
    long long int mid;   // node in the middle of the swath
  };
  std::vector<Segment> segments;
  segments.reserve(n_swaths);
  double min_width = std::numeric_limits<double>::max();
  for (size_t j = 0; j < n_swaths; ++j) {
    const F2CSwath& s = swaths_by_cells.getSwath(j);
    // same mid point as in createCoverageGraph
    const F2CPoint mid_p {(s.startPoint() + s.endPoint()) * 0.5};
    const bool start_is_low = along(s.startPoint()) <= along(s.endPoint());
    const long long int start_id = node_id(s.startPoint());
    const long long int end_id = node_id(s.endPoint());
    segments.push_back(Segment{across(mid_p),
        std::min(along(s.startPoint()), along(s.endPoint())),
        start_is_low ? start_id : end_id, start_is_low ? end_id : start_id,
        node_id(mid_p)});
    if (s.getWidth() > 0.0) {
      min_width = std::min(min_width, s.getWidth());
    }
  }
  // Segments cut from the same line have the same across coordinate,
  // parallel lines are one swath width apart
  const double tol = min_width < std::numeric_limits<double>::max() ?
      0.25 * min_width : 1e-6;

  std::sort(segments.begin(), segments.end(),
      [](const Segment& a, const Segment& b) {return a.across < b.across;});
  std::vector<std::vector<Segment>> lines;
  for (auto&& seg : segments) {
    if (lines.empty() || seg.across - lines.back().front().across > tol) {
      lines.emplace_back();
    }
    lines.back().push_back(seg);
  }

  // Start from the line on the border of the field, cover all its segments
  // (jumping over the obstacles between them), then go back along the next one
  std::vector<long long int> route;
  route.reserve(3 * n_swaths);
  for (size_t k = 0; k < lines.size(); ++k) {
    auto& line = lines[k];
    const bool forward = (k % 2 == 0);
    std::sort(line.begin(), line.end(),
        [forward](const Segment& a, const Segment& b) {
          return forward ? a.along < b.along : a.along > b.along;});
    for (auto&& seg : line) {
      route.push_back(forward ? seg.low : seg.high);
      route.push_back(seg.mid);
      route.push_back(forward ? seg.high : seg.low);
    }
  }

  // Every node but the depot (the last one) must be visited exactly once
  const size_t n_nodes = cov_graph.numNodes();
  std::vector<bool> seen(n_nodes, false);
  for (auto&& id : route) {
    if (id < 0 || static_cast<size_t>(id) + 1 >= n_nodes || seen[id]) {
      throw std::runtime_error(
          "SWATH_SNAKE: the swaths don't match the coverage graph");
    }
    seen[id] = true;
  }
  if (route.size() + 1 != n_nodes) {
    throw std::runtime_error("SWATH_SNAKE: the route does not visit every node");
  }
  return route;
}

std::vector<long long int> RoutePlannerBase::computeBestRoute(
    const F2CGraph2D& cov_graph, bool show_log, long int time_limit_seconds,
    bool use_guided_local_search) const {
  int depot_id = static_cast<int>(cov_graph.numNodes()-1);
  const ortools::RoutingIndexManager::NodeIndex depot{depot_id};
  ortools::RoutingIndexManager manager(cov_graph.numNodes(), 1, depot);
  ortools::RoutingModel routing(manager);

  const int transit_callback_index = routing.RegisterTransitCallback(
      [&cov_graph, &manager] (long long int from, long long int to) -> long long int {
        auto from_node = manager.IndexToNode(from).value();
        auto to_node = manager.IndexToNode(to).value();
        return cov_graph.getCostFromEdge(from_node, to_node);
      });
  routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);
  ortools::RoutingSearchParameters searchParameters =
    ortools::DefaultRoutingSearchParameters();
  searchParameters.set_use_full_propagation(false);
  const bool swath_snake = (first_solution_strategy_ == kSwathSnakeStrategy);
  // kSwathSnakeStrategy is not an or-tools strategy: the first solution is
  // given to the solver, the strategy in the parameters is not used
  searchParameters.set_first_solution_strategy(swath_snake ?
    ortools::FirstSolutionStrategy::AUTOMATIC :
    static_cast<ortools::FirstSolutionStrategy::Value>(
        first_solution_strategy_));
  if (use_guided_local_search) {
    searchParameters.set_local_search_metaheuristic(
      ortools::LocalSearchMetaheuristic::GUIDED_LOCAL_SEARCH);
  } else {
    searchParameters.set_local_search_metaheuristic(
      ortools::LocalSearchMetaheuristic::AUTOMATIC);
  }
  searchParameters.mutable_time_limit()->set_seconds(time_limit_seconds);
  searchParameters.set_log_search(show_log);
  const ortools::Assignment* solution = nullptr;
  if (swath_snake) {
    if (initial_route_.empty()) {
      throw std::runtime_error(
          "SWATH_SNAKE: no initial route, set it with setInitialRoute");
    }
    routing.CloseModelWithParameters(searchParameters);
    std::vector<std::vector<int64_t>> routes(1);
    routes[0].reserve(initial_route_.size());
    for (auto&& node : initial_route_) {
      routes[0].push_back(manager.NodeToIndex(
          ortools::RoutingIndexManager::NodeIndex(static_cast<int>(node))));
    }
    const ortools::Assignment* initial =
      routing.ReadAssignmentFromRoutes(routes, true);
    if (initial == nullptr) {
      throw std::runtime_error(
          "SWATH_SNAKE: the initial route is not feasible");
    }
    solution = routing.SolveFromAssignmentWithParameters(
        initial, searchParameters);
  } else {
    solution = routing.SolveWithParameters(searchParameters);
  }
  if (solution == nullptr) {
    throw std::runtime_error("No route found by the solver");
  }

  long long int index = routing.Start(0);
  std::vector<long long int> v_id;

  index = solution->Value(routing.NextVar(index));

  while (!routing.IsEnd(index)) {
    v_id.emplace_back(manager.IndexToNode(index).value());
    index = solution->Value(routing.NextVar(index));
  }
  return v_id;
}

F2CRoute RoutePlannerBase::transformSolutionToRoute(
    const std::vector<long long int>& route_ids,
    const F2CSwathsByCells& swaths_by_cells,
    const F2CGraph2D& coverage_graph,
    F2CGraph2D& shortest_graph) const {
  F2CRoute route;
  const size_t NS = swaths_by_cells.sizeTotal();
  for (int i = 0; i < route_ids.size()-2; ++i) {
    F2CPoint p_s = coverage_graph.indexToNode(route_ids[i]);
    F2CPoint p_e = coverage_graph.indexToNode(route_ids[i+2]);
    for (int j = 0; j < NS; ++j) {
      F2CSwath swath = swaths_by_cells.getSwath(j).clone();
      if (p_s == swath.startPoint() && p_e == swath.endPoint()) {
        if (route.isEmpty() && r_start_end) {
          route.addConnection(shortest_graph.shortestPath(
                *r_start_end, swath.startPoint()));
        }
        route.addSwath(swath, shortest_graph);
        break;
      } else if (p_e == swath.startPoint() && p_s == swath.endPoint()) {
        swath.reverse();
        if (route.isEmpty() && r_start_end) {
          route.addConnection(shortest_graph.shortestPath(
                *r_start_end, swath.startPoint()));
        }
        route.addSwath(swath, shortest_graph);
        break;
      }
    }
  }
  if (r_start_end) {
    route.addConnection(shortest_graph.shortestPath(
          route.endPoint(), *r_start_end));
  }
  return route;
}

}  // namespace f2c::rp


