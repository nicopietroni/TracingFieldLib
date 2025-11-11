#ifndef GRID_BASICS_2D
#define GRID_BASICS_2D

#include "../segment2.h"
#include "../point2.h"
#include "../box2.h"

namespace Geo{

template<class ScalarType>
void BestDim(const int n_elems,
             const Geo::Point2<ScalarType> & size2D,
             Geo::Point2<int> & dim )
{
    const int mincells   = 1;   // NUMBER OF CELLS
    const ScalarType GFactor = 1;	// GridEntry = NumElem*GFactor
    ScalarType diag = size2D.Norm();// Diagonal OF THE BOUNDING BOX
    ScalarType eps  = diag*1e-4;	// TOLERANCE

    assert(n_elems>0);
    assert(size2D.X>=0);
    assert(size2D.Y>=0);


    int ncell = (int)(n_elems*GFactor);	// NUMBER OF VOXELS
    if(ncell<mincells)
        ncell = mincells;

    dim.X = 1;
    dim.Y = 1;

    if(size2D.X>eps)
    {
        if(size2D.Y>eps)
        {
            ScalarType k = pow((ScalarType)(ncell/(size2D.X*size2D.Y)),ScalarType(1.0/3.f));
             dim.X = int(size2D.X * k);
             dim.Y = int(size2D.Y * k);
        }
    }

    dim.X = std::max(dim.X,1);
    dim.Y = std::max(dim.Y,1);
};

template <class ScalarType>
Geo::Box2<ScalarType> TriangleBBox(const  typename std::vector< Geo::Point2<ScalarType> > &vert,
                                   const  std::vector<std::vector<int> > &face,
                                    const int IndexF)
{
    typedef Geo::Point2<ScalarType> MyCoordType;
    Geo::Box2<ScalarType> Ret;
    for (int j=0;j<3;j++)
    {
        int IndexV=face[IndexF][j];
        MyCoordType pos=vert[IndexV];
        Ret.Add(pos);
    }
    return Ret;
}

template <class ScalarType>
Geo::Box2<ScalarType> SegmentBBox(const Geo::Segment2<ScalarType> &seg)
{
    Geo::Box2<ScalarType> Ret;
    Ret.Add(seg.P(0));
    Ret.Add(seg.P(1));
    return Ret;
}

}

#endif
