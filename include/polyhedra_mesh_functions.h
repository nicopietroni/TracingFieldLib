/*
 * //Non-Commercial License Agreement
 * //Copyright (c) 2025, Nico Pietroni
 * //
 * //1. Permission is hereby granted to use, copy, modify, and distribute
 * //   this software for academic research, teaching, and personal
 * //   non-commercial purposes, provided that this copyright notice is
 * //   included in all copies or substantial portions of the software.
 * //
 * //2. Commercial use of this software, including but not limited to use in
 * //   commercial products, services, tools, or consulting, is strictly
 * //   prohibited without prior written permission from the copyright holder.
 * //
 * //3. This software is provided "as is", without warranty of any kind,
 * //   express or implied, including but not limited to the warranties of
 * //   merchantability, fitness for a particular purpose, and
 * //   noninfringement. In no event shall the authors or copyright holders
 * //   be liable for any claim, damages, or other liability, arising from,
 * //   out of, or in connection with the software or the use or other
 * //   dealings in the software.
 * //
 * //4. For commercial licensing, please contact:
 * //   Nico Pietroni
 * //   Email: nico.pietroni@gmail.com
 */

#ifndef POLYHEDRA_MESH_FUNCTIONS
#define POLYHEDRA_MESH_FUNCTIONS

#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "point3.h"
#include "tetra3.h"
#include "polyhedra_connectivity.h"
// #include "box3.h"
// #include "connectivity.h"
// #include "normal.h"
// #include "polyhedra_connectivity.h"
// #include "tetra3.h"

namespace Tetra
{


template <class ScalarType>
void FindCenter(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                const std::vector<std::vector<int> > &Connectivity,
                std::vector<Geo::Point3<ScalarType> > &Centers)
{
    Centers.clear();
    for (size_t i=0;i<Connectivity.size();i++)
    {
        Geo::Point3<ScalarType> Center(0,0,0);
        for (size_t j=0;j<Connectivity[i].size();j++)
            Center+=VertPos[Connectivity[i][j]];

        Center/=Connectivity[i].size();
        Centers.push_back(Center);
    }
}


// template <class ScalarType>
// class Tetra3
// {
//     Geo::Point3<ScalarType> P[4];

// public:

//     Tetra3(const Geo::Point3<ScalarType> &_P0,
//            const Geo::Point3<ScalarType> &_P1,
//            const Geo::Point3<ScalarType> &_P2,
//            const Geo::Point3<ScalarType> &_P3)
//     {
//         P[0]=_P0;
//         P[1]=_P1;
//         P[2]=_P2;
//         P[3]=_P3;
//     }

//     Geo::Point3<ScalarType> Barycenter()
//     {
//         return ((P[0] + P[1] + P[2] + P[3]) / (ScalarType)4.0);
//     }

//     //    ScalarType Volume()
//     //    {
//     //        Point3<ScalarType> Vect0=Cross((P[2] - P[0]) , (P[1] - P[0]));
//     //        Point3<ScalarType> Vect1=(P[3] - P[0]);
//     //        ScalarType Vol=Dot(Vect0,Vect1) / 6.0;
//     //        return Vol;
//     //    }

//     ScalarType Volume()const
//     {
//         Geo::Point3<ScalarType> L0 = P[1] - P[0];
//         Geo::Point3<ScalarType> L2 = P[0] - P[2];
//         Geo::Point3<ScalarType> L3 = P[3] - P[0];

//         return (Dot(Cross(L2,L0),L3) / 6.0);
//     }

//     Geo::Box3<ScalarType> GetBox()const
//     {
//         std::vector<Geo::Point3<ScalarType> > VertPos;
//         VertPos.push_back(P[0]);
//         VertPos.push_back(P[1]);
//         VertPos.push_back(P[2]);
//         VertPos.push_back(P[3]);
//         Geo::Box3<ScalarType> BB;
//         BB.Init(VertPos);
//         return BB;
//     }

//     bool IsInside(const Geo::Point3<ScalarType> &TestP)const
//     {
//         Geo::Box3<ScalarType> BB=GetBox();
//         if (!BB.IsInside(TestP))return false;

//         for (size_t i=0;i<4;i++)
//         {
//             //swap because normals has to look inside here
//             int LocV0=VofF(i, 0);
//             int LocV1=VofF(i, 2);
//             int LocV2=VofF(i, 1);

//             assert(LocV0>=0);
//             assert(LocV0<4);
//             assert(LocV1>=0);
//             assert(LocV1<4);
//             assert(LocV2>=0);
//             assert(LocV2<4);

//             assert(LocV0!=LocV1);
//             assert(LocV1!=LocV2);

//             Tetra3<ScalarType> CurrT(P[LocV0],P[LocV1],P[LocV2],TestP);
//             ScalarType CurrV=CurrT.Volume();
//             if (CurrV<0)return false;
//         }
//         return true;
//     }

//     Geo::Point3<ScalarType> NormalF(const int IndexF)
//     {
//         assert(IndexF>=0);
//         assert(IndexF<4);

//         //swap because normals has to look inside here
//         int LocV0=VofF(IndexF, 0);
//         int LocV1=VofF(IndexF, 1);
//         int LocV2=VofF(IndexF, 2);

//         assert(LocV0>=0);
//         assert(LocV0<4);
//         assert(LocV1>=0);
//         assert(LocV1<4);
//         assert(LocV2>=0);
//         assert(LocV2<4);

//         assert(LocV0!=LocV1);
//         assert(LocV1!=LocV2);

//         Geo::Point3<ScalarType> PF0=P[LocV0];
//         Geo::Point3<ScalarType> PF1=P[LocV1];
//         Geo::Point3<ScalarType> PF2=P[LocV2];

//         return Geo::Normal(PF0,PF1,PF2);
//     }

//     ScalarType DiehedralAngle(const int IndexE)
//     {
//         int F0idx=FofE(IndexE,0);
//         int F1idx=FofE(IndexE,1);
//         Geo::Point3<ScalarType> N0=NormalF(F0idx);
//         Geo::Point3<ScalarType> N1=NormalF(F1idx);
//         return(Geo::AngleRad(N0,-N1));
//     }
// };

template <class ScalarType>
ScalarType DiehedralAngle(const std::vector<std::vector<int> > &Tetra,
                          const std::vector<Geo::Point3<ScalarType> > &TetraPos,
                          size_t IndexT,size_t IndexE)
{
    assert(IndexT<Tetra.size());
    assert(IndexE<6);

    Tetra3<ScalarType> T(TetraPos[Tetra[IndexT][0]],
            TetraPos[Tetra[IndexT][1]],
            TetraPos[Tetra[IndexT][2]],
            TetraPos[Tetra[IndexT][3]]);
    return (T.DiehedralAngle(IndexE));
}

template <class ScalarType>
Tetra::Tetra3<ScalarType> GetTetra(const std::vector<std::vector<int> > &Tetra,
                                   const std::vector<Geo::Point3<ScalarType> > &TetraPos,
                                   size_t IndexT)
{
    assert(IndexT<Tetra.size());
    assert(Tetra[IndexT].size()==4);

    int IndexV0=Tetra[IndexT][0];
    int IndexV1=Tetra[IndexT][1];
    int IndexV2=Tetra[IndexT][2];
    int IndexV3=Tetra[IndexT][3];

    assert(IndexV0>=0);
    assert(IndexV0<TetraPos.size());
    assert(IndexV1>=0);
    assert(IndexV1<TetraPos.size());
    assert(IndexV2>=0);
    assert(IndexV2<TetraPos.size());
    assert(IndexV3>=0);
    assert(IndexV3<TetraPos.size());

    Geo::Point3<ScalarType> P0=TetraPos[IndexV0];
    Geo::Point3<ScalarType> P1=TetraPos[IndexV1];
    Geo::Point3<ScalarType> P2=TetraPos[IndexV2];
    Geo::Point3<ScalarType> P3=TetraPos[IndexV3];

    return (Tetra::Tetra3<ScalarType>(P0,P1,P2,P3));
}

template <class ScalarType>
ScalarType VolumeTetMesh(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                         const std::vector<std::vector<int> > &Connectivity)
{
    ScalarType CurrVol=0;

    for (size_t i=0;i<Connectivity.size();i++)
    {
        //check is a tetrahedron
        assert(Connectivity[i].size()==4);

        int IndexV0=Connectivity[i][0];
        int IndexV1=Connectivity[i][1];
        int IndexV2=Connectivity[i][2];
        int IndexV3=Connectivity[i][3];

        assert(IndexV0>=0);
        assert(IndexV1>=0);
        assert(IndexV2>=0);
        assert(IndexV3>=0);

        assert(IndexV0<VertPos.size());
        assert(IndexV1<VertPos.size());
        assert(IndexV2<VertPos.size());
        assert(IndexV3<VertPos.size());
        Tetra3<ScalarType> T(VertPos[IndexV0],VertPos[IndexV1],VertPos[IndexV2],VertPos[IndexV3]);

        CurrVol+=T.Volume();
    }
    return CurrVol;
}

//template <class ScalarType>
//ScalarType AverageEdgeSize(const std::vector<Point3<ScalarType> > &VertPos,
//                           const std::vector<std::vector<int> > &Connectivity)
//{
//    ScalarType AdgeEdge=0;
//    size_t NumE=0;
//    for (size_t i=0;i<Connectivity.size();i++)
//    {
//        size_t numV=Connectivity[i].size();
//        assert(numV==4);
//        for (size_t j=0;j<numV;j++)
//            for (size_t k=(j+1);j<numV;j++)
//            {

//                int IndexV0=Connectivity[i][j];
//                int IndexV1=Connectivity[i][(j+1)%numV];

//                assert(IndexV0>=0);
//                assert(IndexV1>=0);

//                assert(IndexV0<VertPos.size());
//                assert(IndexV1<VertPos.size());

//                Point3<ScalarType> P0=VertPos[IndexV0];
//                Point3<ScalarType> P1=VertPos[IndexV1];

//                AdgeEdge+=(P0-P1).Norm();
//                NumE++;
//            }
//    }
//    return (AdgeEdge/NumE);
//}

template <class ScalarType>
ScalarType AverageEdgeSize(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                           const std::vector<std::vector<int> > &Connectivity)
{
    //get average tetra volume
    ScalarType Vol=VolumeTetMesh(VertPos,Connectivity);
    //std::cout<<"Volume"<<Vol<<std::endl;
    Vol/=Connectivity.size();
    ScalarType edge_length = std::cbrt(12/1.4142135*Vol);
    return edge_length;

    //    for (size_t i=0;i<Connectivity.size();i++)
    //    {
    //        size_t numV=Connectivity[i].size();
    //        assert(numV==4);
    //        for (size_t j=0;j<numV;j++)
    //            for (size_t k=(j+1);j<numV;j++)
    //            {

    //                int IndexV0=Connectivity[i][j];
    //                int IndexV1=Connectivity[i][(j+1)%numV];

    //                assert(IndexV0>=0);
    //                assert(IndexV1>=0);

    //                assert(IndexV0<VertPos.size());
    //                assert(IndexV1<VertPos.size());

    //                Point3<ScalarType> P0=VertPos[IndexV0];
    //                Point3<ScalarType> P1=VertPos[IndexV1];

    //                AdgeEdge+=(P0-P1).Norm();
    //                NumE++;
    //            }
    //    }
    //    return (AdgeEdge/NumE);
}


};

#endif
