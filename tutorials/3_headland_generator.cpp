//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================


#include "fields2cover.h"
#include <iostream>

int main() {
  f2c::Random rand(45);
  F2CField field = rand.generateRandField(100*100, 8);
  F2CCells cells = field.getField();
  F2CRobot robot (2.0, 6.0);

  std::cout << "####### Tutorial 3.1 Constant width headland generator ######" << std::endl;
  f2c::hg::ConstHL const_hl; // headland generator
  F2CCells no_hl = const_hl.generateHeadlands(cells, 3.0 * robot.getWidth());
  std::cout << "The complete area is " << cells.area() <<
    ", and the area without headlands is " << no_hl.area() << std::endl;


  f2c::obj::RemArea rem_area;

  std::cout << "The efficency of the Headland Generator is: " <<
    rem_area.computeCost(cells, no_hl)   << std::endl;

  f2c::Visualizer::figure();
  f2c::Visualizer::plot(field);
  f2c::Visualizer::plot(no_hl);
  f2c::Visualizer::save("Tutorial_3_1_Const_width.png");

  return 0;
}

