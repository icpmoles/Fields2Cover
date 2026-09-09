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
#include <thread>
#include <future>
#include <algorithm>
#include <csignal>
#include <mutex>
#include <queue>


using namespace std::chrono;

namespace f2c::rp {

namespace ortools = operations_research;

F2CRoute RoutePlannerBase::genRoute(const F2CCells& cells,
    const F2CSwathsByCells& swaths,
    std::vector<ThreadResult> &result_output,
    ThreadResult &best,
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
    bool free_space_planner,
    bool constrained,
    bool parallel)
{
  // bool parallel = true;
  F2CGraph2D shortest_graph = createShortestGraph(cells, swaths, d_tol, free_space_planner);

  F2CGraph2D cov_graph = createCoverageGraph(cells, swaths, shortest_graph, d_tol, redirect_swaths,
      dist_exponent, use_visibility, visibility_factor, visibility_use_crossing,
      prefer_inter_rings_crossing, free_space_planner, constrained);

  // std::vector<ThreadResult> result_output;

  std::vector<long long int> v_route;
  if (parallel) {
    v_route = computeBestRouteParallel(
        cov_graph, show_log, time_limit_seconds, result_output, best,search_for_optimum,
        constrained, 1);

    for (auto result : result_output) {
      if (result.success) {
        result.p_route = transformSolutionToRoute(
          result.route, swaths, cov_graph, shortest_graph);
      }
    }
  }
  else {
     v_route = computeBestRoute(
      cov_graph, show_log, time_limit_seconds, search_for_optimum, constrained);
  }

  return transformSolutionToRoute(
      v_route, swaths, cov_graph, shortest_graph);
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
    bool free_space_planner,
    bool constrained) const {
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
    const F2CGraph2D& cov_graph,
    bool show_log,
    long int time_limit_seconds,
    bool use_guided_local_search,
    bool constrained) const {
  const size_t n_nodes = cov_graph.numNodes();
  int depot_id = static_cast<int>(n_nodes-1);
  const ortools::RoutingIndexManager::NodeIndex depot{depot_id};
  ortools::RoutingIndexManager manager(n_nodes, 1, depot);
  ortools::RoutingModel routing(manager);
  const int transit_callback_index = routing.RegisterTransitCallback(
      [&cov_graph, &manager /*, &constrained*/] (long long int from, long long int to) -> long long int {
        const auto from_node = manager.IndexToNode(from).value();
        const auto to_node = manager.IndexToNode(to).value();
        return cov_graph.getCostFromEdge(from_node, to_node);
        // if (constrained && cost>1<<29) {
        //     std::cout << from_node<< " to "<< to_node  << " = "  << cost << std::endl;
        // }+
      });

  // get access to the underlying solver
  // operations_research::Solver* const solver = routing.solver();
  if (constrained) {
    long long int counter = 0;
    std::cout << "Inserting ROUTE constraints" << std::endl;
    for (int from_node = 0; from_node<n_nodes-1; ++from_node) {
      for (int to_node = from_node; to_node<n_nodes-1; ++to_node) {
        if (to_node == from_node) {continue;}
        const auto cost = cov_graph.getCostFromEdge(from_node, to_node);
        if (cost > 1<<29) {
          counter++;
            const auto from_idx = manager.NodeToIndex(ortools::RoutingIndexManager::NodeIndex(from_node));
            const auto to_idx = manager.NodeToIndex(ortools::RoutingIndexManager::NodeIndex(to_node));

            // solver->AddConstraint(solver->MakeNonEquality(routing.NextVar(from_idx), to_idx));
            // solver->AddConstraint(solver->MakeNonEquality(routing.NextVar(to_idx), from_idx));
          // MORE IDOMATIC WAY to FORBID connections
          // https://github.com/google/or-tools/blob/v9.9/examples/cpp/random_tsp.cc#L143
            routing.NextVar(from_idx)->RemoveValue(to_idx);
            routing.NextVar(to_idx)->RemoveValue(from_idx);
        }
      }
    }
    std::cout << "TOTAL: "<< counter << std::endl;
  }

  routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);
  ortools::RoutingSearchParameters searchParameters =
    ortools::DefaultRoutingSearchParameters();
  searchParameters.set_use_full_propagation(true);
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

  if (!solution) {
    return {};
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
std::vector<long long int> RoutePlannerBase::computeBestRouteParallel(F2CGraph2D& cov_graph,
    bool show_log,
    long int time_limit_seconds,
    std::vector<ThreadResult> &output, // breaking!!
    ThreadResult &best,
    bool use_guided_local_search,
    bool constrained,
    uint max_cores
    ) const {

    if (max_cores < 1) {
        throw std::logic_error("max_cores should be greater than 1");
    }

    CVrpData data;
    data.Cov_Graph = cov_graph;
    data.constrained = constrained;
    data.show_log = show_log;
    data.time_limit_seconds = time_limit_seconds;
    data.use_guided_local_search = use_guided_local_search;

    const size_t n_nodes = cov_graph.numNodes();

    if (constrained) {
        std::cout << "(PARALLEL) Inserting ROUTE constraints" << std::endl;

        for (size_t from_node = 0; from_node < n_nodes; ++from_node) {
            std::vector<bool> row_constraints;
            row_constraints.reserve(n_nodes);

            for (size_t to_node = 0; to_node < n_nodes; ++to_node) {
                const auto cost =
                    cov_graph.getCostFromEdge(from_node, to_node);

                const bool is_disallowed = cost > (1 << 29);

                row_constraints.push_back(is_disallowed);
            }

            data.constraints.push_back(std::move(row_constraints));
        }
    }

    // Create one configuration/job for each strategy.
    std::vector<SearchConfig> configs;
    configs.reserve(this->possibleStrategies.size());

    for (const auto& str : this->possibleStrategies) {
        configs.push_back({str});
    }

    if (configs.empty()) {
        return {};
    }

    // Limit the number of simultaneously running jobs.
    const unsigned int hw_threads =
        std::thread::hardware_concurrency();

    const size_t hardware_threads =
        hw_threads > 0
            ? hw_threads
            : static_cast<size_t>(max_cores);

    const size_t worker_count = std::min(
        {
            static_cast<size_t>(max_cores),
            hardware_threads,
            configs.size()
        });

    std::cout << "Launching "
              << configs.size()
              << " initializations across "
              << worker_count
              << " workers...\n";

    // ------------------------------------------------------------
    // Work queue
    // ------------------------------------------------------------

    std::queue<size_t> work_queue;

    for (size_t i = 0; i < configs.size(); ++i) {
        work_queue.push(i);
    }

    std::mutex queue_mutex;

    // Each worker writes to its own result slot.
    // Therefore, no result mutex is required.
    std::vector<ThreadResult> results(configs.size());

    // ------------------------------------------------------------
    // Worker function
    // ------------------------------------------------------------

    auto worker = [&]() {
        while (true) {
            size_t job_index;

            {
                std::lock_guard<std::mutex> lock(queue_mutex);

                if (work_queue.empty()) {
                    return;
                }

                job_index = work_queue.front();
                work_queue.pop();
            }

            // Do not hold queue_mutex while running the expensive job.
            results[job_index] =
                this->RunSingleInitialization(
                    data,
                    configs[job_index]);
        }
    };

    // ------------------------------------------------------------
    // Start workers
    // ------------------------------------------------------------

    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t i = 0; i < worker_count; ++i) {
        workers.emplace_back(worker);
    }

    // ------------------------------------------------------------
    // Wait for all workers
    // ------------------------------------------------------------

    for (auto& worker_thread : workers) {
        worker_thread.join();
    }

    // ------------------------------------------------------------
    // Find global best solution
    // ------------------------------------------------------------

    ThreadResult best_result{
        std::numeric_limits<int64_t>::max(),
      false,
        ortools::FirstSolutionStrategy::UNSET,
        {},
      {},
      -1,
      0
    };

    // save results
    output = results;


    for (const auto& result : results) {
        if (result.cost < best_result.cost) {
            best_result = result;
        }
    }
    best.strategy = best_result.strategy;

    std::cout << "best strategy: " << FirstSolutionStrategy.at(best_result.strategy) << std::endl;
    return best_result.route;
}

//
//   if (max_cores < 2 ) {
//     throw std::logic_error("max_cores should be greater than 2");
//   }
//
//   CVrpData data;
//   data.Cov_Graph = cov_graph;
//   data.constrained = constrained;
//   data.show_log = show_log;
//   data.time_limit_seconds = time_limit_seconds;
//   data.use_guided_local_search = use_guided_local_search;
//
//   // std::vector<std::vector<bool>> constraints;
//   const size_t n_nodes = cov_graph.numNodes();
//
//   if (constrained) {
//     long long int counter = 0;
//     std::cout << "(PARALLEL) Inserting ROUTE constraints" << std::endl;
//     for (int from_node = 0; from_node<n_nodes; ++from_node) {
//       std::vector<bool> row_constraints;
//       for (int to_node = 0; to_node<n_nodes; ++to_node) {
//         const auto cost = cov_graph.getCostFromEdge(from_node, to_node);
//         const bool is_disallowed = cost > 1<<29;
//         row_constraints.push_back(is_disallowed);
//       }
//       data.constraints.push_back(row_constraints);
//     }
//     std::cout << "TOTAL: "<< counter << std::endl;
//   }
//
//
//   std::vector<SearchConfig> configs;
//
//   for (auto str: this->possibleStrategies) {
//     configs.push_back({str});
//   }
//
//   const unsigned int hw_threads = std::thread::hardware_concurrency();
//   const size_t max_concurrency = hw_threads > 0 ? hw_threads : max_cores;
//
//   std::queue<size_t> work_queue;
//   std::mutex queue_mutex;
//   std::condition_variable cv;
//   // bool done = false;
//
//   std::vector<ThreadResult> results(configs.size());
//   const size_t worker_count = std::min(max_concurrency, configs.size());
//
//   std::vector<std::future<ThreadResult>> futures;
//   futures.reserve(configs.size());
//
//   std::cout << "Launching " << configs.size() << " initializations across "
//               << max_concurrency << " hardware cores...\n";
//
//   for (const auto& config : configs) {
//     futures.push_back(
//           std::async(std::launch::async, [this, &data, config]() {
//             return this->RunSingleInitialization(data, config);
//           })
//       );
//   }
//   // Collect results and determine the global best solution
//   ThreadResult best_result{std::numeric_limits<int64_t>::max(),
//     ortools::FirstSolutionStrategy::UNSET, {} };
//
//   for (auto& fut : futures) {
//     ThreadResult res = fut.get();
//     if (res.cost < best_result.cost) {
//       best_result = res;
//     }
//   }
//
//   return best_result.route;
// }

// Worker function: Executed concurrently across threads
ThreadResult RoutePlannerBase::RunSingleInitialization(CVrpData& data, SearchConfig config) const {

  const size_t n_nodes = data.Cov_Graph.numNodes();
  int depot_id = static_cast<int>(n_nodes-1);
  const ortools::RoutingIndexManager::NodeIndex depot{depot_id};

  ortools::RoutingIndexManager manager(
        n_nodes,
        1,
        depot
    );

  ortools::RoutingModel routing(manager);

  const int transit_callback_index = routing.RegisterTransitCallback(
    [&data, &manager /*, &constrained*/] (long long int from, long long int to) -> long long int {
      const auto from_node = manager.IndexToNode(from).value();
      const auto to_node = manager.IndexToNode(to).value();
      const auto cost = data.Cov_Graph.getCostFromEdge(from_node, to_node);
      // if (constrained && cost>1<<29) {
      //     std::cout << from_node<< " to "<< to_node  << " = "  << cost << std::endl;
      // }
      return cost;
    });

  // setup constraints
  if (data.constrained) {
    long long int counter = 0;
    std::cout << "(SINGLE) Inserting ROUTE constraints" << std::endl;
    for (int from_node = 0; from_node<n_nodes-1; ++from_node) {
      for (int to_node = from_node; to_node<n_nodes-1; ++to_node) {
        if (to_node == from_node) {continue;}
        // const auto is_constrained = data.constraints.at(from_node).at(to_node);
        if (data.constraints.at(from_node).at(to_node)) {
          counter++;
          const auto from_idx = manager.NodeToIndex(ortools::RoutingIndexManager::NodeIndex(from_node));
          const auto to_idx = manager.NodeToIndex(ortools::RoutingIndexManager::NodeIndex(to_node));

          // solver->AddConstraint(solver->MakeNonEquality(routing.NextVar(from_idx), to_idx));
          // solver->AddConstraint(solver->MakeNonEquality(routing.NextVar(to_idx), from_idx));
          routing.NextVar(from_idx)->RemoveValue(to_idx);
          routing.NextVar(to_idx)->RemoveValue(from_idx);
        }
      }
    }
    std::cout << "TOTAL: "<< counter << std::endl;
  }

  routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);

  ortools::RoutingSearchParameters searchParameters =
   ortools::DefaultRoutingSearchParameters();
  searchParameters.set_use_full_propagation(true);
  searchParameters.set_first_solution_strategy(config.strategy);
  searchParameters.set_report_intermediate_cp_sat_solutions(true);

    // ortools::FirstSolutionStrategy::FIRST_UNBOUND_MIN_VALUE);
  if (data.use_guided_local_search) {
    searchParameters.set_local_search_metaheuristic(
      ortools::LocalSearchMetaheuristic::GUIDED_LOCAL_SEARCH);
  } else {
    searchParameters.set_local_search_metaheuristic(
      ortools::LocalSearchMetaheuristic::AUTOMATIC);
  }
  searchParameters.mutable_time_limit()->set_seconds(data.time_limit_seconds);
  searchParameters.set_log_search(data.show_log);
  const ortools::Assignment* solution =
    routing.SolveWithParameters(searchParameters);
  // std::cout  << "Routing excuted" <<std::endl;

  if (solution != nullptr) {
    ThreadResult result;
    result.success = true;
    result.cost = solution->ObjectiveValue();
    result.strategy = config.strategy;
    // result.random_seed = config.random_seed;
    // std::vector<int64_t> route;
    std::vector<long long int> v_id;
    long long int index = routing.Start(0);

    index = solution->Value(routing.NextVar(index));

    // std::cout  << "getting value, new index = " << index <<std::endl;
    while (!routing.IsEnd(index)) {
      result.route.emplace_back(manager.IndexToNode(index).value());
      index = solution->Value(routing.NextVar(index));
    }
    // std::cout  << "solution exported" <<std::endl;
    std::cout << FirstSolutionStrategy.at(config.strategy) << " Solved successfully = " << result.cost <<std::endl;

    return result;
  } else {
    std::cout << FirstSolutionStrategy.at(config.strategy)  << " NOT Solved successfully" <<std::endl;
  }


  return {std::numeric_limits<int64_t>::max(), false, config.strategy,  {},{}};;
}

}  // namespace f2c::rp


