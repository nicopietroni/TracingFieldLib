#ifndef SMOOTH_MESH
#define SMOOTH_MESH

#include <cmath>
#include <vector>
#include "point3.h"

namespace Geo
{

template <class ScalarType>
void SmoothLaplacian(std::vector<Geo::Point3<ScalarType> > &VertPos,
                     const std::vector<std::vector<int> > &Faces,
                     int NumIte,ScalarType Damp=0.5)
{

    for (int s=0;s<NumIte;s++)
    {
        std::vector<Geo::Point3<ScalarType> > NewVertPos(VertPos.size(),Geo::Point3<ScalarType>(0,0,0));
        std::vector<int> SumNum(VertPos.size(),0);

        for (size_t i=0;i<Faces.size();i++)
            for (size_t j=0;j<Faces[i].size();j++)
            {
                int IndexV0=Faces[i][j];
                int IndexV1=Faces[i][(j+1)%Faces[i].size()];
                Geo::Point3<ScalarType> P0=VertPos[IndexV0];
                Geo::Point3<ScalarType> P1=VertPos[IndexV1];
                NewVertPos[IndexV0]+=P1;
                NewVertPos[IndexV1]+=P0;
                SumNum[IndexV0]++;
                SumNum[IndexV1]++;
            }
        for (size_t i=0;i<NewVertPos.size();i++)
        {
            NewVertPos[i]/=(ScalarType)SumNum[i];
            VertPos[i]=NewVertPos[i]*(1-Damp)+VertPos[i]*Damp;
        }
    }
}

}
#endif
