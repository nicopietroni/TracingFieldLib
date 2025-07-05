#ifndef POLYGON3_H
#define POLYGON3_H

#include <assert.h>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

#include "best_rotation.h"
#include "line3.h"
#include "plane3.h"
#include "point3.h"
#include "segment3.h"
#include "polyline3.h"

namespace Geo {
  template <class ScalarType>
ScalarType BestAlignedPoly(const std::vector<Geo::Point3<ScalarType>> &FacePos,
                           std::vector<Geo::Point3<ScalarType>> &bestPoly) {
  // find the perimeter
  Geo::Polyline3<ScalarType> Poly;
  Poly.PolyPos = FacePos;
  ScalarType Len = Poly.Lenght();

  ScalarType sideLength = Len / FacePos.size(); // Compute side length
  ScalarType radius =
      sideLength / (2 * sin(M_PI / FacePos.size())); // Compute circumradius

  bestPoly.clear();
  for (int i = 0; i < FacePos.size(); ++i) {
    ScalarType angle = (2.0 * M_PI * i) / FacePos.size(); // Angle in radians
    ScalarType x = radius * cos(angle);
    ScalarType y = radius * sin(angle);
    bestPoly.push_back(Geo::Point3<ScalarType>(x, y, 0));
  }

  Geo::PointMatch<ScalarType>::TrasformRigidMatch(FacePos, bestPoly);

  ScalarType Err = 0;
  for (size_t i = 0; i < FacePos.size(); i++)
    Err += pow((FacePos[i] - bestPoly[i]).Norm(), 2);
  return Err;
}


} // namespace Geo

#endif
