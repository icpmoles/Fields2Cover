//
// Created by icpmoles on 03/10/25.
//

#include "fields2cover.h"
#include <iostream>

int main() {
  // Import field
  // F2CFields fields;
  F2CField field = f2c::Parser::importFieldGml(std::string(DATA_PATH) + "Envelope.xml");
  // F2CField field = fields[1];
  // Transform into UTM to work in meters
  // f2c::Transform::transformToUTM(field);
  std::string epsg_system;
  field.setEPSGCoordSystem(3857);

  std::cout << "EPSG: " << field.getEPSGCoordSystem() << std::endl;


  f2c::Visualizer::figure();
  f2c::Visualizer::plot(field);
  f2c::Visualizer::save("poli.png");

  f2c::Transform::transformToUTM(field);
  f2c::Visualizer::figure();
  f2c::Visualizer::plot(field);
  f2c::Visualizer::save("poli_utm.png");



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