#ifndef TRIANGLE3_H
#define TRIANGLE3_H

#include <assert.h>
#include "point3.h"
#include "line3.h"
#include "segment3.h"
#include "plane3.h"
#include "triangle2.h"
#include <vector>

namespace Geo{

template <class ScalarType>
class Traingle3{
public:
    typedef Geo::Point3<ScalarType> CoordType;

private:
    CoordType Pos[3];

public:

    CoordType &P(const int &IndexV)
    {
        assert(IndexV>=0);
        assert(IndexV<3);
        return Pos[IndexV];
    }

    const CoordType &P(const int &IndexV)const
    {
        assert(IndexV>=0);
        assert(IndexV<3);
        return Pos[IndexV];
    }

    CoordType N()const
    {
        CoordType nDir=Cross((P(1)-P(0)),(P(2)-P(0)));
        nDir.Normalize();
        return(nDir);
    }

    Traingle3(const CoordType &P0,const CoordType &P1,const CoordType &P2)
    {
        Pos[0]=P0;
        Pos[1]=P1;
        Pos[2]=P2;
    }

};

template <class ScalarType>
bool IsInside(const Traingle3<ScalarType> &T,
              const Geo::Point3<ScalarType> &p)
{
    typedef typename Geo::Point3<ScalarType> CoordType;

    CoordType n=Cross((T.P(1)-T.P(0)),(T.P(2)-T.P(0)));
    CoordType n0=Cross((T.P(0)-p),(T.P(1)-p));
    CoordType n1=Cross((T.P(1)-p),(T.P(2)-p));
    CoordType n2=Cross((T.P(2)-p),(T.P(0)-p));
    return ((Dot(n,n0)>=0)&&(Dot(n,n1)>=0)&&(Dot(n,n2)>=0));
}

template <class ScalarType>
Geo::Point3<ScalarType> ClosestPoint(const Traingle3<ScalarType> &T,
                                     const Geo::Point3<ScalarType> &Test_pos)
{
    //std::cout<<"closest check"<<std::endl;

    typedef typename Geo::Point3<ScalarType> CoordType;

    ScalarType distv[3];
    ScalarType distproj;

    //find distance on the plane
    Geo::Plane3<ScalarType> TrPlane(T.P(0),T.P(1),T.P(2));
    Geo::Point3<ScalarType> clos_proj=TrPlane.Project(Test_pos);

    if (IsInside<ScalarType>(T,clos_proj))
        return (clos_proj);

    //distance from the edges
    Segment3<ScalarType> e0=Segment3<ScalarType>(T.P(0),T.P(1));
    Segment3<ScalarType> e1=Segment3<ScalarType>(T.P(1),T.P(2));
    Segment3<ScalarType> e2=Segment3<ScalarType>(T.P(2),T.P(0));

    CoordType clos_e[3];
    clos_e[0]=ClosestPoint(e0,Test_pos);
    clos_e[1]=ClosestPoint(e1,Test_pos);
    clos_e[2]=ClosestPoint(e2,Test_pos);

    ScalarType dist_e[3];

    dist_e[0]=(clos_e[0]-Test_pos).Norm();
    dist_e[1]=(clos_e[1]-Test_pos).Norm();
    dist_e[2]=(clos_e[2]-Test_pos).Norm();

    if ((dist_e[0]<=dist_e[1])&&(dist_e[0]<=dist_e[2]))return clos_e[0];
    if ((dist_e[1]<=dist_e[0])&&(dist_e[1]<=dist_e[2]))return clos_e[1];
    /*if ((dist_e[2]<dist_e[0])&&(dist_e[2]<dist_e[1]))*/
    return clos_e[2];
}

template <class ScalarType>
Geo::Point3<ScalarType> ClosestVert(const Traingle3<ScalarType> &T,
                                    const Geo::Point3<ScalarType> &Test_pos)
{
    ScalarType D0=(T.P(0)-Test_pos).Norm();
    ScalarType D1=(T.P(1)-Test_pos).Norm();
    ScalarType D2=(T.P(2)-Test_pos).Norm();
    if ((D0<D1)&&(D0<D2))return T.P(0);
    if ((D1<D0)&&(D1<D2))return T.P(1);
    return T.P(2);
}

template <class ScalarType>
inline Geo::Point2<ScalarType> RemoveComponent( const Geo::Point3<ScalarType> &P,
                                                const int &compIndex)
{
    assert(compIndex>=0);
    assert(compIndex<3);
    if (compIndex==0)return Geo::Point2<ScalarType>(P.Y,P.Z);
    if (compIndex==1)return Geo::Point2<ScalarType>(P.X,P.Z);
    return Geo::Point2<ScalarType>(P.X,P.Y);
}

template <class ScalarType>
inline Geo::Traingle2<ScalarType> RemoveComponent( const Geo::Traingle3<ScalarType> &T,
                                                   const int &compIndex)
{
    Geo::Traingle2<ScalarType> T2(RemoveComponent(T.P(0),compIndex),
                                  RemoveComponent(T.P(1),compIndex),
                                  RemoveComponent(T.P(2),compIndex));
    return T2;
}

template <class ScalarType>
void BarycentricCoords(const Traingle3<ScalarType> &T,
                       const Geo::Point3<ScalarType> &Test_pos,
                       ScalarType &alpha,ScalarType &beta,ScalarType &gamma)
{
    Geo::Point3<ScalarType> N=T.N();
    ScalarType MaxV;
    size_t Axis;
    N.GetMaxAbsCoord(MaxV,Axis);
    Geo::Point2<ScalarType> P0,P1,P2;

    Geo::Traingle2<ScalarType> T2=RemoveComponent(T,Axis);
    Geo::Point2<ScalarType> Test_pos2D=RemoveComponent(Test_pos,Axis);
    Geo::BarycentricCoords(T2,Test_pos2D,alpha,beta,gamma);
}


} //namespace


#endif
