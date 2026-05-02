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

  const bool new_gen=false;
  std::cout << "createShortestGraph" << std::endl;
  F2CGraph2D shortest_graph = createShortestGraph(cells, swaths, d_tol, new_gen);

  std::cout << "createCoverageGraph" << std::endl;
  F2CGraph2D cov_graph = createCoverageGraph(
      cells, swaths, shortest_graph, d_tol, redirect_swaths, new_gen);

  std::cout << "computeBestRoute" << std::endl;
  std::vector<long long int> v_route = computeBestRoute(
      cov_graph, show_log, time_limit_seconds, search_for_optimum, new_gen);
  return transformSolutionToRoute(
      v_route, swaths, cov_graph, shortest_graph);
}

void RoutePlannerBase::setStartAndEndPoint(const F2CPoint& p) {
  this->r_start_end = p;
}

F2CGraph2D RoutePlannerBase::createShortestGraph(const F2CCells &cells,
                                      const F2CSwathsByCells &swaths_by_cells,
                                      double d_tol, bool new_gen) const {
  F2CGraph2D g;
  g.setNewGen(new_gen);
  // Add points from swaths that touches border
  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      g.addEdge(s.startPoint(), cells.closestPointOnBorderTo(s.startPoint()));
      g.addEdge(s.endPoint(),   cells.closestPointOnBorderTo(s.endPoint()));
      if (false) {
        const int64_t swath_cost = g.getScale()*pow((s.length()*20.0+1000.0),1.5);
        g.addEdge(s.startPoint(), s.endPoint(),swath_cost, true);
      }
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

  // connect internal obstacles with the external ring if visible or with the closest
  // visible internal ring otherwise
  // for (auto&& cell : cells) {
  //   // first ring is always the external
  //   for (size_t int_ring_idx = 1; int_ring_idx < cell.size(); ++int_ring_idx) {
  //     const auto ring = cell.getInteriorRing(int_ring_idx);
  //     for (size_t i = 0; i < ring.size()-1; ++i) {
  //       g.addEdge(ring.getGeometry(i), ring.getGeometry(i+1));
  //     }
  //   }
  // }
  //

  // for (int cell_idx=0; cell_idx < cells.size(); ++cell_idx) {
  //   auto cell = cells[cell_idx];
  //   auto swaths = swaths_by_cells[cell_idx];
  //   for (auto&& s : swaths) {
  //     }
  // }

  // for (auto&& swaths : swaths_by_cells) {
  //   for (auto&& s : swaths) {
  //     g.addEdge(s.startPoint(), s.endPoint(),g.getScale()*pow((s.length()*20.0+1000),1.5));
  //
  //   }
  // }

  // Add start and end point if they exists
  if (this->r_start_end) {
    g.addEdge(*r_start_end, cells.closestPointOnBorderTo(*r_start_end));
  }

  std::vector<F2CPoint> nodes = g.getNodes();

  // Connect nodes that are near other edges.
  for (int i = 0; i < 2; ++i) {
    auto edges = g.getEdges();
    // For every "from"
    for (auto&& edge : edges) {
      size_t from = edge.first;
      //To every "to"
      for (auto&& e : edge.second) {
        size_t to = e.first;
        // create a tentative border
        F2CLineString border {g.indexToNode(from), g.indexToNode(to)};
        for (auto&& n : nodes) {
          // if there is a node in the graph close enough to the tentative
          // border which is not already assigned:
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
    const F2CCells &cells, const F2CSwathsByCells &swaths_by_cells,
    F2CGraph2D &shortest_graph, double d_tol, bool redirect_swaths,
    bool new_gen) const {
  F2CGraph2D g;
  g.setNewGen(new_gen);
  // connect mid swaths
  std::cout << "-->createCoverageGraph: connect mid swaths" << std::endl;
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

  // nb: here is where the swaths are interconnected though the headlands
  // if the swaths extremities are not connected directly their weights remain "INFINITE"
  std::cout << "-->createCoverageGraph: swaths extremities" << std::endl;
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

  std::cout << "-->createCoverageGraph: add deposit" << std::endl;
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

std::vector<long long int> RoutePlannerBase::computeBestRoute(
    const F2CGraph2D &cov_graph, bool show_log, long int time_limit_seconds,
    bool use_guided_local_search, bool new_gen) const {
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
  searchParameters.set_first_solution_strategy(
    ortools::FirstSolutionStrategy::AUTOMATIC);
  if (use_guided_local_search) {
    searchParameters.set_local_search_metaheuristic(
      ortools::LocalSearchMetaheuristic::GUIDED_LOCAL_SEARCH);
  } else {
    searchParameters.set_local_search_metaheuristic(
      ortools::LocalSearchMetaheuristic::AUTOMATIC);
  }
  searchParameters.mutable_time_limit()->set_seconds(time_limit_seconds);
  searchParameters.set_log_search(show_log);
  const ortools::Assignment* solution =
    routing.SolveWithParameters(searchParameters);

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
  F2CSwath swath_gen;
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


