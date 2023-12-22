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

template <class ScalarType>
bool Intersect(const Segment2<ScalarType> &seg0,
               const Segment2<ScalarType> &seg1,
               ScalarType &lambda0)
{
    const ScalarType Eps= ScalarType(1e-8);
    ScalarType lambda1;
    const Geo::Point2<ScalarType> & p0 = seg0.P(0);
    const Geo::Point2<ScalarType> & p1 = seg0.P(1);
    const Geo::Point2<ScalarType> & p2 = seg1.P(0);
    const Geo::Point2<ScalarType> & p3 = seg1.P(1);

    ScalarType a = (p1-p0).X;
    ScalarType b = (p2-p3).X;
    ScalarType c = (p1-p0).Y;
    ScalarType d = (p2-p3).Y;

    ScalarType e = (p2-p0).X;
    ScalarType f = (p2-p0).Y;

    ScalarType det = a*d-b*c;

    lambda0 = (d*e-b*f)/det;
    lambda1 = (-c*e+a*f)/det;

    if (fabs(det)<Eps)
        return false;// parallell

    if (!(lambda0 >= 0.0 && lambda0 <= 1.0 && lambda1 >= 0.0 && lambda1 <= 1.0))
        return false;

    lambda0=1-lambda0;
    //p_inters = p0*(1-lambda0)+p1*lambda0;
    return true;
}

//const SCALAR_TYPE Eps= SCALAR_TYPE(1e-8);
//SCALAR_TYPE lambda0,lambda1;
//const Point2<SCALAR_TYPE> & p0 = seg0.P0();
//const Point2<SCALAR_TYPE> & p1 = seg0.P1();
//const Point2<SCALAR_TYPE> & p2 = seg1.P0();
//const Point2<SCALAR_TYPE> & p3 = seg1.P1();

//SCALAR_TYPE a = (p1-p0)[0];
//SCALAR_TYPE b = (p2-p3)[0];
//SCALAR_TYPE c = (p1-p0)[1];
//SCALAR_TYPE d = (p2-p3)[1];

//SCALAR_TYPE e = (p2-p0)[0];
//SCALAR_TYPE f = (p2-p0)[1];

//SCALAR_TYPE det = a*d-b*c;

//lambda0 = (d*e-b*f)/det;
//lambda1 = (-c*e+a*f)/det;
//if (fabs(det)<Eps)
//    return false;// they are parallell

//if (!(lambda0 >= 0.0 && lambda0 <= 1.0 && lambda1 >= 0.0 && lambda1 <= 1.0))
//    return false;
//p_inters = p0*(1-lambda0)+p1*lambda0;
//return true;
} //namespace


#endif
