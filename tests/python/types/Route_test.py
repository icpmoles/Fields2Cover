#==============================================================================
#     Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
#                      Author: Gonzalo Mier
#                         BSD-3 License
#==============================================================================

import pytest
import fields2cover as f2c

def near(a, b):
  assert pytest.approx(a) == pytest.approx(b)

def test_fields2cover_types_route_init():
  assert (f2c.Route().asLineString().isEmpty());
  route = f2c.Route();
  route.addConnection(f2c.MultiPoint(
      f2c.VectorPoint([f2c.Point(0, 0)])));
  swaths1 = f2c.Swaths();
  swaths1.push_back(f2c.Swath(f2c.LineString(f2c.VectorPoint(
      [f2c.Point(0, 0), f2c.Point(1, 0)])), 4));
  route.addSwaths(swaths1);
  route.addConnection(f2c.MultiPoint(f2c.VectorPoint(
      [f2c.Point(1, 0), f2c.Point(1, 1)])));
  swaths2 = f2c.Swaths();
  swaths2.push_back(f2c.Swath(f2c.LineString(f2c.VectorPoint(
      [f2c.Point(1, 1), f2c.Point(0, 1)])), 4));
  route.addSwaths(swaths2);
  route.addConnection(f2c.MultiPoint());
  assert not (route.asLineString().isEmpty());
  near(route.asLineString().length(), 3);

def test_fields2cover_types_route_to_path():
    route = f2c.Route();
    robot = f2c.Robot(1.0, 2.0);
    swath1 = f2c.Swath(f2c.LineString(f2c.VectorPoint([
        f2c.Point(1, 1), f2c.Point(1, 2), f2c.Point(1, 4)
    ])))
    conn12 = f2c.MultiPoint()
    [conn12.addGeometry(p) for p in
        [f2c.Point(1.1, 4), f2c.Point(1.5, 4), f2c.Point(2, 4)]];
    swath2 = f2c.Swath(f2c.LineString(f2c.VectorPoint([
        f2c.Point(2,4), f2c.Point(2, 2), f2c.Point(2,1)
    ])))

    route.addSwath(swath1);
    route.addConnection(conn12);
    route.addSwath(swath2);
    path = route.asPath(robot);

    assert (path.size()==7)
