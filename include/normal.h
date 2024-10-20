#ifndef COMPUTE_NORMAL
#define COMPUTE_NORMAL

#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdio.h>
#include <iostream>
#include "point3.h"
#include "triangle3.h"
#include <assert.h>

namespace  Geo{

template <class ScalarType>
Point3<ScalarType> Normal(const Point3<ScalarType> &P0,
                          const Point3<ScalarType> &P1,
                          const Point3<ScalarType> &P2)
{
    Point3<ScalarType> Dir0=P1-P0;
    Point3<ScalarType> Dir1=P2-P0;
    Dir0.Normalize();
    Dir1.Normalize();
    Point3<ScalarType> Norm=Cross(Dir0,Dir1);
    Norm.Normalize();
    return(Norm);
}


template <class ScalarType>
Point3<ScalarType> Normal(const std::vector<Point3<ScalarType> > &Vert,
                          const std::vector<std::vector<int> > &Faces,
                          int IndexF)
{
    assert(IndexF>=0);
    assert(IndexF<(int)Faces.size());
    assert(Faces[IndexF].size()>2);
    int sizeF=(int)Faces[IndexF].size();
    if (sizeF==3)
    {
        int IndexV0=Faces[IndexF][0];
        int IndexV1=Faces[IndexF][1];
        int IndexV2=Faces[IndexF][2];
        assert(IndexV0>=0);
        assert(IndexV0<Vert.size());
        assert(IndexV1>=0);
        assert(IndexV1<Vert.size());
        assert(IndexV2>=0);
        assert(IndexV2<Vert.size());
        Point3<ScalarType> P0=Vert[IndexV0];
        Point3<ScalarType> P1=Vert[IndexV1];
        Point3<ScalarType> P2=Vert[IndexV2];
        Point3<ScalarType> Dir0=P1-P0;
        Point3<ScalarType> Dir1=P2-P0;
        Dir0.Normalize();
        Dir1.Normalize();
        Point3<ScalarType> NormDir=Cross(Dir0,Dir1);
        NormDir.Normalize();
        return(NormDir);
    }

    Point3<ScalarType> Norm(0,0,0);
    for (int i=0;i<sizeF;i++)
    {
        int IndexV0=Faces[IndexF][(i+sizeF-1)%sizeF];
        int IndexV1=Faces[IndexF][i];
        int IndexV2=Faces[IndexF][(i+1)%sizeF];
        Point3<ScalarType> P0=Vert[IndexV0];
        Point3<ScalarType> P1=Vert[IndexV1];
        Point3<ScalarType> P2=Vert[IndexV2];
        Point3<ScalarType> Dir0=P1-P0;
        Point3<ScalarType> Dir1=P2-P0;
        Dir0.Normalize();
        Dir1.Normalize();
        Norm+=Cross(Dir0,Dir1);
    }
    Norm.Normalize();
    return(Norm);
}

template <class ScalarType>
Point3<ScalarType> NormalVert(const std::vector<Point3<ScalarType> > &Vert,
                              const std::vector<std::vector<int> > &Faces,
                              const std::vector<std::vector<int> > &VF,
                              int IndexV)
{
    Point3<ScalarType> Norm(0,0,0);
    assert(IndexV>=0);
    assert(IndexV<VF.size());
    for (size_t i =0;i<VF[IndexV].size();i++)
    {
        int IndexF=VF[IndexV][i];
        Norm+=Normal(Vert,Faces,IndexF);
    }
    Norm.Normalize();
    return Norm;
}

template <class ScalarType>
void ComputeNormals(const std::vector<Point3<ScalarType> > &Vert,
                    const std::vector<std::vector<int> > &Faces,
                    std::vector<Point3<ScalarType> > &FaceNormals,
                    std::vector<Point3<ScalarType> > &VertNormals)
{
    FaceNormals.clear();
    FaceNormals.resize(Faces.size());


    for (size_t i=0;i<Faces.size();i++)
        FaceNormals[i]=Normal(Vert,Faces,i);

    VertNormals.clear();
    VertNormals.resize(Vert.size(),Point3<ScalarType>(0,0,0));
    for (size_t i=0;i<Faces.size();i++)
    {
        Point3<ScalarType> FaceN=FaceNormals[i];
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int IndexV=Faces[i][j];
            assert(IndexV>=0);
            assert(IndexV<(int)VertNormals.size());
            VertNormals[IndexV]+=FaceN;
        }
    }

    for (size_t i=0;i<VertNormals.size();i++)
    {
        if (VertNormals[i]==Point3<ScalarType>(0,0,0))continue;
        VertNormals[i].Normalize();
    }
}


}

#endif
