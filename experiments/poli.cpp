//
// Created by icpmoles on 03/10/25.
//

#include "fields2cover.h"
#include <iostream>

int main() {
  // Import field
  F2CFields fields;
  f2c::Parser::importJson(std::string(DATA_PATH) + "giurati_poles.geojson", fields);
  F2CField field = fields[0];
  // Transform into UTM to work in meters
  // f2c::Transform::transformToUTM(field);
  std::string epsg_system;
  field.setEPSGCoordSystem(7791);

  std::cout << "EPSG: " << field.getEPSGCoordSystem() << std::endl;


  F2CRobot robot (0.20, 1.5);
  f2c::hg::ConstHL const_hl;
  F2CCells no_hl = const_hl.generateHeadlands(field.getField(), 3.0 * robot.getWidth());
  f2c::sg::BruteForce bf;
  F2CSwaths swaths = bf.generateSwaths(M_PI, robot.getCovWidth(), no_hl.getGeometry(0));
  f2c::rp::BoustrophedonOrder sorter;
  swaths = sorter.genSortedSwaths(swaths);
  f2c::pp::PathPlanning path_planner;
  robot.setMinTurningRadius(0.1);  // m
  f2c::pp::DubinsCurves dubins;
  F2CPath path = path_planner.planPath(robot, swaths, dubins);

  f2c::Visualizer::figure();
  f2c::Visualizer::figure_size(1500, 1500);
  f2c::Visualizer::plot(field);

  f2c::Visualizer::plot(swaths);
  f2c::Visualizer::plot(no_hl);
  f2c::Visualizer::plot(path);
  f2c::Visualizer::save("campo_elaborato.png");

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