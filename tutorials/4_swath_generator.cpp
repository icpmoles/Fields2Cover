//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================


#define ALLOW_PARALLELIZATION
#include "fields2cover.h"
#include <iostream>


int main() {
  f2c::Random rand(420);
  F2CRobot robot (2.0, 6.0);
  f2c::hg::ConstHL const_hl; // constant headland
  F2CCells cells = rand.generateRandField(10000, 5).getField();
  F2CCells no_hl = const_hl.generateHeadlands(cells, 2.0 * robot.getWidth());

  std::cout << "####### Tutorial 4.1 Brute force swath generator ######" << std::endl;

  f2c::sg::BruteForce bf_sw_gen;
  f2c::obj::NSwath n_swath_obj;
  // bf_sw_gen.setStepAngle(0.001);
  bf_sw_gen.setOffsetDivisions(50);
  // F2CSwaths swaths_bf_nswath = bf_sw_gen.generateBestSwaths(n_swath_obj, robot.getCovWidth(), no_hl.getGeometry(0));
  //
  //
  // f2c::Visualizer::figure();
  // f2c::Visualizer::plot(cells);
  // f2c::Visualizer::plot(no_hl);
  // f2c::Visualizer::plot(swaths_bf_nswath);
  // f2c::Visualizer::save("Tutorial_4_1_Brute_force_NSwath.png");




  F2CSwaths swaths_bf_angle = bf_sw_gen.generateSwaths(M_PI/2, -0.0,robot.getCovWidth(), no_hl.getGeometry(0));
  f2c::Visualizer::figure();
  f2c::Visualizer::plot(cells);
  f2c::Visualizer::plot(no_hl);
  f2c::Visualizer::plot(swaths_bf_angle);
  f2c::Visualizer::save("Tutorial_4_1_Brute_force_Angle.png");

  return 0;
}

