//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <ortools/constraint_solver/routing.h>
#include <ortools/constraint_solver/routing_enums.pb.h>
#include <ortools/constraint_solver/routing_index_manager.h>
#include <ortools/constraint_solver/routing_parameters.h>
#include <cmath>
#include <utility>
#include <vector>
#include <limits>
#include "fields2cover/route_planning/route_planner_base.h"
#include <ctime>
#include <chrono>
using namespace std::chrono;

namespace f2c::rp {

namespace ortools = operations_research;

F2CRoute RoutePlannerBase::genRoute(const F2CCells& cells,
    const F2CSwathsByCells& swaths,
    std::string& timings_array_dest,
    bool show_log,
    double d_tol,
    bool redirect_swaths,
    long int time_limit_seconds,
    bool search_for_optimum,
    float dist_exponent,
    bool use_visibility,
    float visibility_factor,
    bool visibility_use_crossing,
    bool prefer_inter_rings_crossing,
    bool free_space_planner)
{
  for (auto cell : cells) {
    // F2CGraph2D dummy;
    pessimistic_traversal_ += cell.getExteriorRing().getMinSafeLength()*100;
  }
  std::cout<< "pessimistic: " << pessimistic_traversal_ << std::endl;

  auto start_sg = high_resolution_clock::now();
  F2CGraph2D shortest_graph = createShortestGraph(cells, swaths, d_tol, free_space_planner);
  auto end_sg_start_cg = high_resolution_clock::now();

  auto duration_1 = duration_cast<milliseconds>(end_sg_start_cg - start_sg);


  F2CGraph2D cov_graph = createCoverageGraph(cells, swaths, shortest_graph, d_tol, redirect_swaths,
      dist_exponent, use_visibility, visibility_factor, visibility_use_crossing,
      prefer_inter_rings_crossing, free_space_planner);
  auto end_cg_start_vpr = high_resolution_clock::now();


  auto duration_2 = duration_cast<milliseconds>(end_cg_start_vpr - end_sg_start_cg);

  std::vector<long long int> v_route = computeBestRoute(
      cov_graph, show_log, time_limit_seconds, search_for_optimum);

  auto end_vpr_start_trans = high_resolution_clock::now();

  auto duration_3 = duration_cast<milliseconds>(end_vpr_start_trans - end_cg_start_vpr);


  auto ret =  transformSolutionToRoute(
      v_route, swaths, cov_graph, shortest_graph);

  auto end_trans = high_resolution_clock::now();

  auto duration_4 = duration_cast<milliseconds>(end_trans - end_vpr_start_trans);

  timings_array_dest = "[ " + std::to_string(duration_1.count()) + ", " +
    std::to_string(duration_2.count()) + ", " +
    std::to_string(duration_3.count()) + ", " +
    std::to_string(duration_4.count()) + " ]";
  return ret;
}

void RoutePlannerBase::setStartAndEndPoint(const F2CPoint& p) {
  this->r_start_end = p;
}

F2CGraph2D RoutePlannerBase::createShortestGraph(const F2CCells& cells,
    const F2CSwathsByCells& swaths_by_cells,
    double d_tol,
    bool free_space_planner) const {
  F2CGraph2D g;
  // Add points from swaths that touches border
  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      g.addEdge(s.startPoint(), cells.closestPointOnBorderTo(s.startPoint()));
      g.addEdge(s.endPoint(),cells.closestPointOnBorderTo(s.endPoint()));
      // if we don't expect a free space planner, then we can route through swaths
      if (!free_space_planner) {
        const int64_t swath_cost = g.getScalingFactor()*(pow(s.length(),2) + pessimistic_traversal_);
        g.addEdge(s.startPoint(), s.endPoint(),swath_cost);
      }
    }
  }


  // double crossing_penalty = 2.0*pessimistic_traversal_;
  // if (crossing_factor > 0.0) {
  //   crossing_penalty = (2.0-crossing_factor)*pessimistic_traversal_;
  // }
  // for (auto&& s : swaths) {
  //   g.addEdge(s.startPoint(), cells.closestPointOnBorderTo(s.startPoint()));
  //   g.addEdge(s.endPoint(),   cells.closestPointOnBorderTo(s.endPoint()));
  //   if (crossing_factor > 0.0) {
  //     // std::cout  <<  "x-ing" << std::endl;
  //     double swath_cost = pow(s.length() + crossing_penalty,2) *g.getScalingFactor();
  //     g.addEdge(s.startPoint(), s.endPoint(), swath_cost);
  //   }

  // Add points in the border
  for (auto&& cell : cells) {
    for (auto&& ring : cell) {
      for (size_t i = 0; i < ring.size()-1; ++i) {
        g.addEdge(ring.getGeometry(i), ring.getGeometry(i+1));
      }
    }
  }

  // Add start and end point if they exists
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


F2CGraph2D RoutePlannerBase::createCoverageGraph(const F2CCells& cells,
    const F2CSwathsByCells& swaths_by_cells,
    F2CGraph2D& shortest_graph,
    double d_tol,
    bool redirect_swaths,
    float dist_exponent,
    bool use_visibility,
    double visibility_factor,
    bool visibility_use_crossing,
    bool prefer_crossings,
    bool free_space_planner) const {
  if (prefer_crossings == true && use_visibility == false) {
    throw std::invalid_argument("You can only prefer crossings if you allow for visbility checks");
  }

  F2CGraph2D g;
  // int64_t INF = 1<<29;
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
  // std::clock_t start_clock = 0;
  int64_t counter = 0;
  int64_t total_nodes = swaths_by_cells.sizeTotal()*swaths_by_cells.sizeTotal();
  for (const auto& swaths1 : swaths_by_cells) {
    for (const auto& s1 : swaths1) {
      auto s1_s = s1.startPoint();
      auto s1_e = s1.endPoint();
      for (const auto& swaths2 : swaths_by_cells) {
        for (const auto& s2 : swaths2) {
          if (counter%100==0) {
            std::cout << counter << "/" << total_nodes << " = " << (counter*1.0/total_nodes)*100.0 << "%"  << std::endl;
          }
          counter++;
          // if (counter==2) {
          //   start_clock = std::clock();
          // }
          auto s2_s = s2.startPoint();
          auto s2_e = s2.endPoint();
          if (redirect_swaths) {
            for (auto a: {s1_s, s1_e}) {
              for (auto b: {s2_s, s2_e}) {
                int64_t cost_function = shortest_graph.shortestPathCost(a,b);
                const int64_t l2_d = a.distance(b);

                // otherwise if we explicitly want to allow crossings we use a custom cost
                // function that will create crossings between rings
                if (prefer_crossings) {
                  // std::cout << "prefer crossings" << std::endl;
                  const int64_t collisions = cells.countCollisions(a, b, visibility_use_crossing);
                  if (l2_d < cost_function && collisions == 0) {
                    cost_function = l2_d*shortest_graph.getScalingFactor();
                  } else {
                    cost_function +=  pessimistic_traversal_ * ( collisions * visibility_factor) *shortest_graph.getScalingFactor();
                    cost_function += pow(l2_d*shortest_graph.getScalingFactor(), dist_exponent);
                  }
                }
                // we always check for visibility when connecting swaths extremities
                else if ( use_visibility && free_space_planner ) //cost_function > (pessimistic_traversal_ * shortest_graph.getScalingFactor()))
                  {
                  // std::cout << "disconn rings detected" << a << b << std::endl;

                  const int64_t collisions = cells.countCollisions(a, b, visibility_use_crossing);

                  cost_function +=pessimistic_traversal_ * ( collisions * visibility_factor) * shortest_graph.getScalingFactor();
                  cost_function += pow(l2_d, dist_exponent)*shortest_graph.getScalingFactor();
                }

                g.addEdge(a, b, cost_function);
              }
            }
          } else {
            g.addDirectedEdge(s1_e, s2_s, shortest_graph);
          }
        }
      }
    }
  }

  // const auto end_time = std::clock();
  // std::cout << "TIME" << float( end_time - start_clock ) /  CLOCKS_PER_SEC << std::endl;

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


