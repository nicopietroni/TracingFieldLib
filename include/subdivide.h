#ifndef SUBDIVIDE_H
#define SUBDIVIDE_H

#include <array>
#include <assert.h>
#include <vector>
#include "point3.h"

namespace Geo
{

    template <class ScalarType>
    class Subdivide
    {
    public:

       //return the subdivision position between Points[1] and Points[2]
       static Geo::Point3<ScalarType> PolylineEdge(std::array<Geo::Point3<ScalarType>,4> &Points,
                                                   ScalarType omega= ScalarType(1.0) / ScalarType (16.0))
       {
            Point3<ScalarType> new_point = (Points[0]*(-omega)) +
                                  (Points[1]*(ScalarType(0.5) + omega)) +
                                  (Points[2]*(ScalarType(0.5) + omega)) +
                                  (Points[3]*(-omega));
            return new_point;
       }

    };

}

#endif
