//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include "fields2cover/types/MultiPoint.h"
#include "fields2cover/types/Robot.h"

namespace f2c::types {

MultiPoint::MultiPoint() {
  data_ = std::shared_ptr<OGRMultiPoint>(
    static_cast<OGRMultiPoint*>(
      OGRGeometryFactory::createGeometry(wkbMultiPoint)),
    [](OGRMultiPoint* f) {OGRGeometryFactory::destroyGeometry(f);});
}

MultiPoint::MultiPoint(const std::vector<Point>& ps) {
  for (auto&& p : ps) {
    this->addGeometry(p);
  }
}

MultiPoint::MultiPoint(const std::initializer_list<Point>& ps) {
  for (auto&& p : ps) {
    this->addGeometry(p);
  }
}

size_t MultiPoint::size() const {
  return isEmpty() ? 0 : this->data_->getNumGeometries();
}

void MultiPoint::getGeometry(size_t i, Point& point) {
  if (i >= this->size()) {
    throw std::out_of_range(
        "MultiPoint does not contain point " + std::to_string(i));
  }
  point = Point(data_->getGeometryRef(i), EmptyDestructor());
}

void MultiPoint::getGeometry(size_t i, Point& point) const {
  if (i >= this->size()) {
    throw std::out_of_range(
        "MultiPoint does not contain point " + std::to_string(i));
  }
  point = Point(data_->getGeometryRef(i), EmptyDestructor());
}

Point MultiPoint::getGeometry(size_t i) {
  if (i >= this->size()) {
    throw std::out_of_range(
        "MultiPoint does not contain point " + std::to_string(i));
  }
  return Point(data_->getGeometryRef(i));
}

const Point MultiPoint::getGeometry(size_t i) const {
  if (i >= this->size()) {
    throw std::out_of_range(
        "MultiPoint does not contain point " + std::to_string(i));
  }
  return Point(data_->getGeometryRef(i));
}

const Point MultiPoint::getFirstPoint() const {
  return getGeometry(0);
}

const Point MultiPoint::getLastPoint() const {
  return getGeometry(size() - 1);
}


void MultiPoint::setGeometry(size_t i, const Point& p) {
    OGRPoint* point = downCast<OGRPoint*>(data_->getGeometryRef(i));
    point->setX(p.getX());
    point->setY(p.getY());
    point->setZ(p.getZ());
}


void MultiPoint::addGeometry(const Point& p) {
  this->addPoint(p);
}

void MultiPoint::addPoint(const Point& p) {
  this->addPoint(p.getX(), p.getY(), p.getZ());
}

void MultiPoint::addPoints(const MultiPoint& ps) {
  for (auto&& p : ps) {
    this->addPoint(p);
  }
}

void MultiPoint::addPoint(double x, double y, double z) {
  OGRPoint p(x, y, z);
  this->data_->addGeometry(&p);
}


void MultiPoint::operator*=(double b) {
  for (auto&& p : *this) {
    p *= b;
  }
}

double MultiPoint::getInAngle(size_t i) const {
  // First point does not have In Angle
  if (0 == i) {
    throw std::invalid_argument(
        "MultiPoint::getInAngle not defined for first point");
  }
  return (getGeometry(i)-getGeometry(i-1)).getAngleFromPoint();
}

double MultiPoint::getOutAngle(size_t i) const {
  // Last point does not have Out Angle
  if (i >= size() - 1) {
    throw std::invalid_argument(
        "MultiPoint::getOutAngle not defined for last point");
  }
  return (getGeometry(i+1)-getGeometry(i)).getAngleFromPoint();
}

double MultiPoint::getPointAngle(size_t i) const {
  if (2 > size()) {
    throw std::invalid_argument(
        "MultiPoint::getPointAngle not defined when size() < 2");
  }
  if (i == 0) {
    return getOutAngle(i);
  } else if (i == size() - 1) {
    return getInAngle(i);
  }
  return getAngleAvg(getInAngle(i), getOutAngle(i));
}


MultiPoint MultiPoint::densify(double length) const{
  MultiPoint mp;
  // Point prev_p = this->getGeometry(0);
  // std::cout << " first_point "<< prev_p << std::endl;
  // mp.addPoint(prev_p);

  double cumulated_length = 0.0;
  bool jump = false;
  for (int i = 0; i < this->size()-1; i++) {
    Point prev_p = this->getGeometry(i);
    Point p = this->getGeometry(i+1);
    double angle = this->getOutAngle(i);
    if (mp.size()==0) {
      mp.addPoint(prev_p);
    }
    const double segment_length = prev_p.distance(p);
    double remaining_length = segment_length;

    while (remaining_length > length) {
      if (!jump) {
        Point temp = mp.back().getPointFromAngle(angle, length);
        mp.addPoint(temp);
        remaining_length -= length;
      } else {
        Point temp = prev_p.getPointFromAngle(angle, length-cumulated_length);
        mp.addPoint(temp);
        jump = false;
        remaining_length -= (length-cumulated_length);
        cumulated_length = 0;
      }
    }

    cumulated_length = remaining_length;
    // end of last iteration, we can just assign even if it's shorter
    if (i == this->size() - 2 && mp.back().distance(p) > 1e-7) {
      mp.addPoint(p);
      continue;
    }

    // prev_p = p;
    jump = true;
  }

  return mp;
}

// taken from: https://en.wikipedia.org/wiki/Low-pass_filter#Discrete-time_realization
MultiPoint MultiPoint::low_pass_filter(double alpha) const {
  MultiPoint mp;
  mp.addPoint(this->getFirstPoint());
  for (int i = 1; i < this->size(); i++) {
    Point new_p(
      mp.back().getX()+ alpha * (this->getGeometry(i).getX() - mp.back().getX()),
      mp.back().getY()+ alpha * (this->getGeometry(i).getY() - mp.back().getY())
    );
    mp.addPoint(new_p);
  }
  return mp;
}
MultiPoint MultiPoint::smoothify(double length, double alpha) const{
  std::cout << "MultiPoint::smoothify-> starting from:" << this->size();
  auto ret =  this->densify(length).low_pass_filter(alpha);

  std::cout << " to "<< ret.size() << std::endl;
  return ret;
}

}  // namespace f2c::types

