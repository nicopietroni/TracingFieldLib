#ifndef SEGMENT_2_H
#define SEGMENT_2_H

#include <assert.h>
#include "point2.h"

namespace Geo{

template <class ScalarType>
class Segment2{
public:
    typedef Geo::Point2<ScalarType> CoordType;

private:
    CoordType Pos[2];

public:

    const CoordType &P(const int &IndexV)const
    {
        assert(IndexV>=0);
        assert(IndexV<2);
        return Pos[IndexV];
    }


    Segment2(const CoordType &P0,const CoordType &P1)
    {
        Pos[0]=P0;
        Pos[1]=P1;
    }

};


template <class ScalarType>
Geo::Point2<ScalarType> ClosestPoint( const Segment2<ScalarType> &s,
                                      const Geo::Point2<ScalarType> &p)
{

   Geo::Point2<ScalarType> v1 = s.P(1) - s.P(0);
   Geo::Point2<ScalarType> v2 = p-s.P(0);

   ScalarType dot = Dot(v1, v2);
   ScalarType t = dot / (v1.X * v1.X + v1.Y * v1.Y);

       if (t < 0.0) {
           return s.P(0); // Closest point is the start point of the segment
       } else if (t > 1.0) {
           return s.P(1); // Closest point is the end point of the segment
       } else {
           // Closest point is on the segment
           Geo::Point2<ScalarType> closestPoint = s.P(0)+v1*t;
           return closestPoint;
       }
}

template <class ScalarType>
ScalarType Distance( const Segment2<ScalarType> &s,
                     const Geo::Point2<ScalarType> &p)
{
    Geo::Point2<ScalarType> ClosP=ClosestPoint<ScalarType>(s,p);
    return (ClosP-p).Norm();
}

} //namespace


#endif
