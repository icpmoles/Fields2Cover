//
// Created by icpmoles on 03/10/25.
//

#include "fields2cover.h"
#include <iostream>


struct sol_val
{
  float width;
  float cov_width;
  float multiplier;
} ;

int main() {
  // Import field
  F2CFields fields;
  //f2c::Parser::importJson(std::string(DATA_PATH) + "campo_full.geojson", fields);
  f2c::Parser::importJson(std::string(DATA_PATH) + "giurati_poles.geojson", fields);
  F2CField field = fields[0];
  // Transform into UTM to work in meters
  // f2c::Transform::transformToUTM(field);
  std::string epsg_system;
  field.setEPSGCoordSystem(7791);

  std::cout << "EPSG: " << field.getEPSGCoordSystem() << std::endl;

  sol_val success {.width = 0.2, .cov_width = 1.0, .multiplier = 1.5};
  sol_val wrong {.width = 1.1, .cov_width = 1.2, .multiplier=  0.6};
  sol_val corre {.width = 1.05, .cov_width = 1.4, .multiplier=  1.0};

  sol_val selected=wrong;


  F2CRobot robot (selected.width, selected.cov_width);
  f2c::hg::ConstHL const_hl;

  F2CCells field_cell = field.getField();
  F2CCells no_hl = const_hl.generateHeadlands(field_cell, 2.0*selected.multiplier * robot.getWidth());
  F2CCells mid_hl_c = const_hl.generateHeadlands(field_cell, selected.multiplier * robot.getWidth());
  f2c::sg::BruteForce bf;
  // F2CSwaths swaths = bf.generateSwaths(0.01, robot.getCovWidth(), no_hl.getGeometry(0));
  F2CSwathsByCells swathsbc = bf.generateSwaths(0.01, robot.getCovWidth(), no_hl);

  f2c::rp::BoustrophedonOrder sorter;

  std::cout << "Generated " << swathsbc[0].size() << " swaths" << std::endl;

  // f2c::Visualizer::figure();
  // f2c::Visualizer::plot(no_hl);
  // f2c::Visualizer::plot(swathsbc);
  //
  // f2c::Visualizer::show();

  // F2CSwaths boustrophedon_swaths = swaths; //sorter.genSortedSwaths(swaths);
  f2c::pp::PathPlanning path_planner;
  robot.setMinTurningRadius(0.02);  // m
  f2c::pp::DubinsCurves dubins;

  f2c::rp::RoutePlannerBase route_planner;
  F2CRoute route = route_planner.genRoute(mid_hl_c, swathsbc, true, 0.001,
    true, 20, true);
  F2CPath path = path_planner.planPath(robot, route, dubins);

  std::cout << "completed " << std::endl;

  f2c::Visualizer::figure();
  f2c::Visualizer::figure_size(1500, 1500);
  f2c::Visualizer::plot(field);

  // f2c::Visualizer::plot(swathsbc);
  // f2c::Visualizer::plot(no_hl);
  f2c::Visualizer::plot(route);
  f2c::Visualizer::save("campo_elaborato.png");

  f2c::Visualizer::figure();
  f2c::Visualizer::figure_size(1500, 1500);
  f2c::Visualizer::plot(field);

  // f2c::Visualizer::plot(swathsbc);
  f2c::Visualizer::plot(no_hl);
  f2c::Visualizer::plot(mid_hl_c);
  // f2c::Visualizer::plot(route);
  f2c::Visualizer::plot(path);
  f2c::Visualizer::save("campo_elaborato_path.png");

  // f2c::Transform::transformToUTM(field);
  f2c::Visualizer::figure();
  f2c::Visualizer::figure_size(1500, 1500);
  f2c::Visualizer::plot(field);
  f2c::Visualizer::save("campo_source_" + std::to_string(field.getEPSGCoordSystem()) + ".png");


  f2c::Transform::transformToUTM(field);
  f2c::Visualizer::figure();
  f2c::Visualizer::figure_size(1500, 1500);
  f2c::Visualizer::plot(field);
  f2c::Visualizer::save("campo_source_utm.png");

  auto path_gps = f2c::Transform::transformToPrevCRS(path, field);
  path_gps.saveToFile("path.csv",10);
  // std::string ss =route.exportToJson();


  // // Transform the generated path back to the previousa CRS.
  // F2CPath path_gps = f2c::Transform::transformToPrevCRS(path, field);
  // f2c::Transform::transformToPrevCRS(field);
  //
  // f2c::Visualizer::figure();
  // f2c::Visualizer::plot(field.getCellsAbsPosition());
  // f2c::Visualizer::plot(path_gps);
  // f2c::Visualizer::save("Tutorial_8_1_GPS.png");

  return 0;
}