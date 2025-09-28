//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include "fields2cover/swath_generator/swath_generator_base.h"

#include "fields2cover/utils/visualizer.h"

namespace f2c::sg {

bool SwathGeneratorBase::getAllowOverlap() const {
  return this->allow_overlap;
}

void SwathGeneratorBase::setAllowOverlap(bool value) {
  this->allow_overlap = value;
}

F2CSwaths SwathGeneratorBase::generateBestSwaths(
    f2c::obj::SGObjective& obj, double op_width, const F2CCell& poly) {
  double best_angle = computeBestAngle(obj, op_width, poly);
  std::cout << "best_angle: " << best_angle*180.0/3.14 << std::endl;
  return generateSwaths(best_angle, op_width, poly);
}

F2CSwathsByCells SwathGeneratorBase::generateBestSwaths(
    f2c::obj::SGObjective& obj, double op_width, const F2CCells& polys) {
  F2CSwathsByCells swaths;
  for (auto&& p : polys) {
    swaths.emplace_back(generateBestSwaths(obj, op_width, p));
  }
  return swaths;
}

F2CSwathsByCells SwathGeneratorBase::generateSwaths(double angle,
    double op_width, const F2CCells& polys) {
  F2CSwathsByCells swaths;
  for (auto&& p : polys) {
    swaths.emplace_back(generateSwaths(angle, op_width, p));
  }
  return swaths;
}


F2CSwaths SwathGeneratorBase::generateSwaths(double angle,
    double op_width, const F2CCell& poly) {
  auto rot_poly {F2CPoint(0.0, 0.0).rotateFromPoint(-angle, poly)};

  // double field_height {rot_poly.getHeight()};
  F2CPoint min_point(rot_poly.getDimMinX(), rot_poly.getDimMinY());
  F2CPoint centroid = rot_poly.GeometricCentre();
  double y_up = rot_poly.getDimMaxY() - centroid.getY();
  double y_down = centroid.getY() - rot_poly.getDimMinY();

  const F2CLineString seed_curve = rot_poly.createStraightLongLine(centroid, 0.0);
  double curve_y {-0.5 * op_width};
  F2CMultiLineString paths;
  while (std::max(y_up,-y_down) > curve_y + (allow_overlap ? 0.0 : 0.5 * op_width)) {
    curve_y += op_width;
    auto g_up = seed_curve + F2CPoint(0.0, curve_y);
    auto g_down = seed_curve + F2CPoint(0.0, -curve_y);
    paths.addGeometry(F2CPoint(0.0, 0.0).rotateFromPoint(angle,g_up));
    paths.addGeometry(F2CPoint(0.0, 0.0).rotateFromPoint(angle,g_down));

  }

  F2CSwaths swaths;
  swaths.append(paths, poly, op_width, f2c::types::SwathType::MAINLAND);
  return swaths;
}


double SwathGeneratorBase::computeCostOfAngle(f2c::obj::SGObjective& obj,
    double ang, double op_width, const F2CCell& poly) {
  if (obj.isFastCompAvailable()) {
    return obj.computeCostWithMinimizingSign(ang, op_width, poly);
  }
  return obj.computeCostWithMinimizingSign(poly,
      generateSwaths(ang, op_width, poly));
}

double SwathGeneratorBase::computeBestAngle(f2c::obj::SGObjective& obj,
    double width, const F2CCell& poly) {
  return computeBestAngle(obj, width, poly);
}

}  // namespace f2c::sg

