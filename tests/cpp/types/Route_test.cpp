//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <gtest/gtest.h>
#include "fields2cover/types.h"

TEST(fields2cover_types_route, init) {
  EXPECT_TRUE(F2CRoute().asLineString().isEmpty());
  F2CRoute route;

  EXPECT_EQ(route.startPoint(), F2CPoint());
  route.addConnection(F2CMultiPoint({F2CPoint(-2, -1), F2CPoint(0, -1), F2CPoint(0, 0)}));
  EXPECT_EQ(route.getLastConnection().size(), 3);
  EXPECT_EQ(route.startPoint(), F2CPoint(-2, -1));
  F2CSwaths swaths1;
  swaths1.emplace_back(F2CLineString({F2CPoint(0, 0), F2CPoint(1, 0)}), 4);
  swaths1.emplace_back(F2CLineString({F2CPoint(1, 1), F2CPoint(0, 1)}), 4);
  route.addSwaths(swaths1);
  EXPECT_EQ(route.getLastSwaths().size(), 2);
  route.addConnection(F2CMultiPoint({F2CPoint(0, 1), F2CPoint(0, 2)}));
  EXPECT_EQ(route.getLastConnection().size(), 2);
  F2CSwaths swaths2;
  swaths2.emplace_back(F2CLineString({F2CPoint(0, 2), F2CPoint(1, 2)}), 6);
  swaths2.emplace_back(F2CLineString({F2CPoint(1, 3), F2CPoint(0, 3)}), 6);
  swaths2.emplace_back(F2CLineString({F2CPoint(0, 4), F2CPoint(1, 4)}), 6);
  swaths2.emplace_back(F2CLineString({F2CPoint(1, 5), F2CPoint(0, 5)}), 6);
  route.addSwaths(swaths2);

  const F2CRoute c_route {route};
  EXPECT_EQ(c_route.getLastSwaths().size(), 4);
  EXPECT_EQ(c_route.getLastConnection().size(), 2);

  route.addConnection();
  EXPECT_FALSE(route.isEmpty());
  EXPECT_FALSE(route.asLineString().isEmpty());
  EXPECT_EQ(route.sizeVectorSwaths(), 2);
  EXPECT_EQ(route.sizeConnections(), 3);
  EXPECT_EQ(route.asLineString().length(), 14);
  EXPECT_EQ(route.length(), 14);
  EXPECT_EQ(route.clone().length(), 14);
  EXPECT_EQ(route.getSwaths(0).size(), 2);
  EXPECT_EQ(route.getSwaths(1).size(), 4);

  route.addConnection(F2CMultiPoint({F2CPoint(0, 5), F2CPoint(0, 6)}));
  route.addConnection(F2CMultiPoint({F2CPoint(0, 8), F2CPoint(5, 8)}));
  route.addConnection(F2CMultiPoint({F2CPoint(5, 10), F2CPoint(10, 10)}));
  EXPECT_EQ(route.sizeVectorSwaths(), 2);
  EXPECT_EQ(route.sizeConnections(), 3);
  F2CSwaths swaths3;
  swaths3.emplace_back(F2CLineString({F2CPoint(20, 20), F2CPoint(21, 20)}), 6);

  route.addConnectedSwaths(F2CMultiPoint({F2CPoint(10, 10), F2CPoint(20, 10)}), swaths3);
  EXPECT_EQ(route.sizeVectorSwaths(), 3);
  EXPECT_EQ(route.sizeConnections(), 3);
  EXPECT_EQ(route.length(), 50);

  EXPECT_EQ(route.getConnection(0).size(), 3);
  EXPECT_EQ(route.getSwaths(0).size(),     2);
  EXPECT_EQ(route.getConnection(1).size(), 2);
  EXPECT_EQ(route.getSwaths(1).size(),     4);
  EXPECT_EQ(route.getConnection(2).size(), 8);
  EXPECT_EQ(route.getSwaths(2).size(),     1);

  EXPECT_EQ(route.getLastConnection().size(), 8);
  EXPECT_EQ(route.getLastSwaths().size(),     1);
  EXPECT_EQ(route.getLastConnection().back(), F2CPoint(20, 10));
  EXPECT_EQ(route.getLastSwaths().back().getPath().endPoint(), F2CPoint(21, 20));
  F2CPoint p_rand1 (33, 44), p_rand2 (-22, 11), p_rand3(5, 6);
  route.getLastSwaths().back().setPath(F2CLineString(p_rand1, p_rand2));
  size_t n = route.getLastConnection().size();
  route.getLastConnection().setGeometry(n-1, p_rand3);
  EXPECT_EQ(route.getLastSwaths().back().getPath().endPoint(), p_rand2);
  EXPECT_EQ(route.getLastConnection().back(), p_rand3);

  EXPECT_EQ(route.startPoint(), F2CPoint(-2, -1));


  F2CRoute route2;
  route2.addSwaths(c_route.getLastSwaths());
  route2.addConnection(route.getLastConnection());
  EXPECT_EQ(route2.startPoint(), F2CPoint(0, 2));
  EXPECT_EQ(route2.endPoint(), p_rand3);


}


TEST(fields2cover_types_route, convert_to_path) {
  F2CRoute route;
  F2CRobot robot(1., 2.);
  F2CLineString line1({F2CPoint(1.0, 1.0), F2CPoint(1.0, 2.0), F2CPoint(1.0, 4.0)});
  F2CSwath swath1(line1);
  F2CMultiPoint conn12({F2CPoint(1.1, 4.0), F2CPoint(1.5, 4.0), F2CPoint(2.0, 4.0)});

  F2CLineString line2({F2CPoint(2.0, 4.0), F2CPoint(2.0, 2.0), F2CPoint(2.0, 1.0)});
  F2CSwath swath2(line2);
  F2CMultiPoint conn23({F2CPoint(2.0, 1.0), F2CPoint(2.25, 0.75),
      F2CPoint(2.75, 0.75), F2CPoint(3.0, 1.0)});

  F2CLineString line3({F2CPoint(3.0, 1.0), F2CPoint(3.0, 2.0)});
  F2CSwath swath3(line3);
  F2CMultiPoint conn34({F2CPoint(3.0, 2.0), F2CPoint(4.0, 2.0)});

  F2CLineString line4({F2CPoint(4.0, 2.0), F2CPoint(4.0, 1.0)});
  F2CSwath swath4(line4);
  F2CMultiPoint conn4end({F2CPoint(4.0, 1.0), F2CPoint(5.0, 1.0)});


  route.addSwath(swath1);
  route.addConnection(conn12);
  route.addSwath(swath2);
  route.addConnection(conn23);
  route.addSwath(swath3);
  route.addConnection(conn34);
  route.addSwath(swath4);
  route.addConnection(conn4end);

  F2CPath path = route.asPath(robot, false, 0.1, 0.5);
  size_t n = path.size();
  EXPECT_EQ(n, 17);

  EXPECT_EQ(path[0].point.getX(), line1.getGeometry(0).getX());
  EXPECT_EQ(path[0].point.getY(), line1.getGeometry(0).getY());
  EXPECT_EQ(path[0].angle, .5*M_PI);
  EXPECT_EQ(path[0].type, f2c::types::PathSectionType::SWATH);

  EXPECT_EQ(path[1].point.getX(), line1.getGeometry(1).getX());
  EXPECT_EQ(path[1].point.getY(), line1.getGeometry(1).getY());
  EXPECT_EQ(path[1].angle, .5*M_PI);
  EXPECT_EQ(path[1].type, f2c::types::PathSectionType::SWATH);

  EXPECT_EQ(path[2].point.getX(), conn12.getGeometry(0).getX());
  EXPECT_EQ(path[2].point.getY(), conn12.getGeometry(0).getY());
  EXPECT_EQ(path[2].angle, 0);
  EXPECT_EQ(path[2].type, f2c::types::PathSectionType::TURN);

  EXPECT_EQ(path[n-3].point.getX(), line4.getGeometry(0).getX());
  EXPECT_EQ(path[n-3].point.getY(), line4.getGeometry(0).getY());
  EXPECT_EQ(path[n-3].angle,1.5*M_PI);
  EXPECT_EQ(path[n-3].type, f2c::types::PathSectionType::SWATH);

  EXPECT_EQ(path[n-2].point.getX(), conn4end.getGeometry(0).getX());
  EXPECT_EQ(path[n-2].point.getY(), conn4end.getGeometry(0).getY());
  EXPECT_EQ(path[n-2].angle,0);
  EXPECT_EQ(path[n-2].type, f2c::types::PathSectionType::TURN);

  EXPECT_EQ(path[n-1].point.getX(), conn4end.getGeometry(1).getX());
  EXPECT_EQ(path[n-1].point.getY(), conn4end.getGeometry(1).getY());
  EXPECT_EQ(path[n-1].angle,0);
  EXPECT_EQ(path[n-1].type, f2c::types::PathSectionType::TURN);


  F2CRoute route_with_deposit;

  F2CMultiPoint conndep1({F2CPoint(0.0, 0.0), F2CPoint(1.0, 1.0)});
  F2CLineString wd_line1({F2CPoint(1.0, 1.0), F2CPoint(1.0, 2.0), F2CPoint(1.0, 4.0)});
  F2CSwath wd_swath1(wd_line1);
  F2CMultiPoint wd_conn12({F2CPoint(1.1, 4.0), F2CPoint(1.5, 4.0), F2CPoint(2.0, 4.0)});

  F2CLineString wd_line2({F2CPoint(2.0, 4.0), F2CPoint(2.0, 2.0), F2CPoint(2.0, 1.0)});
  F2CSwath wd_swath2(wd_line2);
  F2CMultiPoint wd_conn2dep({F2CPoint(2.0, 1.0), F2CPoint(2.0, 0.0),
      F2CPoint(1.0, 0.0), F2CPoint(0.0, 0.0)});

  route_with_deposit.addConnection(conndep1);
  route_with_deposit.addSwath(wd_swath1);
  route_with_deposit.addConnection(wd_conn12);
  route_with_deposit.addSwath(wd_swath2);
  route_with_deposit.addConnection(wd_conn2dep);

  F2CPath path_wd = route_with_deposit.asPath(robot, false, 0.1, 0.5);
  size_t wd_n =path_wd.size();
  EXPECT_EQ(wd_n,2+2+3+2+4);
  EXPECT_EQ(path_wd[0].point.getX(), path_wd[wd_n-1].point.getX());
  EXPECT_EQ(path_wd[0].point.getY(), path_wd[wd_n-1].point.getY());

}