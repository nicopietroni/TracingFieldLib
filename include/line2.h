#ifndef LINE2
#define LINE2

#include "point2.h"

namespace Geo{

template <class ScalarType>
class Line2
{
public:
  typedef Point2<ScalarType> CoordType;

private:

    //Origin
    CoordType lineOrigin;
    //Direction (not necessarily normalized, unless so specified by NORM)
    CoordType lineDirection;

public:

  inline const CoordType &Origin() const { return lineOrigin; }
  inline CoordType &Origin() { return lineOrigin; }

  inline const CoordType &Direction() const { return lineDirection; }
  inline  CoordType &Direction() { return lineDirection; }

  inline CoordType Pos( const ScalarType t ) const
  { return Origin() + Direction() * t; }

    Line2(const CoordType &_lineOrigin, const CoordType &_lineDirection)
    {
        lineOrigin=_lineOrigin;
        lineDirection=_lineDirection;
        lineDirection.Normalize();
    }

    Line2()
    {
        lineOrigin=CoordType(0,0,0);
        lineDirection=CoordType(0,0,0);
    }

    /// Projects a point on the line
    inline ScalarType Projection( const  CoordType &p ) const
    {
        return ScalarType(Dot(p-Origin(),Direction()));
        //else      return ScalarType((p-_ori).dot(_dir)/_dir.SquaredNorm());
    }

    inline CoordType P( const ScalarType t ) const
    { return Origin() + Direction() * t; }


    CoordType ClosestPoint(const CoordType & p) const
    {
        return P(Projection(p));
    }

    ScalarType Distance(const CoordType & p) const
    {
        CoordType ClosP=ClosestPoint(p);
        return (p-ClosP).Norm();
    }
};

} // end namespace
#endif
