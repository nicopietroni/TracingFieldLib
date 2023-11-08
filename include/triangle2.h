#ifndef TRIANGLE2_H
#define TRIANGLE2_H

#include <assert.h>
#include "point2.h"

namespace Geo{

template <class ScalarType>
class Traingle2{
public:
    typedef Geo::Point2<ScalarType> CoordType;

private:
    CoordType Pos[3];

public:

    const CoordType &P(const int &IndexV)const
    {
        assert(IndexV>=0);
        assert(IndexV<3);
        return Pos[IndexV];
    }

    Traingle2(const CoordType &P0,const CoordType &P1,const CoordType &P2)
    {
        Pos[0]=P0;
        Pos[1]=P1;
        Pos[2]=P2;
    }

};

template <class ScalarType>
void BarycentricCoords(const Traingle2<ScalarType> &T,
                       const Geo::Point2<ScalarType> &Test_pos,
                       ScalarType &alpha,ScalarType &beta,ScalarType &gamma)
{
    ScalarType x1=T.P(0).X;
    ScalarType x2=T.P(1).X;
    ScalarType x3=T.P(2).X;

    ScalarType y1=T.P(0).Y;
    ScalarType y2=T.P(1).Y;
    ScalarType y3=T.P(2).Y;

    ScalarType x=Test_pos.X;
    ScalarType y=Test_pos.Y;

    alpha=((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
    beta=((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/((y3-y1)*(x2-x3)+(x1-x3)*(y2-y3));
    gamma=1-alpha-beta;
}

} //namespace


#endif
