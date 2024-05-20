#ifndef MESH_FUNCTIONS
#define MESH_FUNCTIONS

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "point3.h"
#include "point2.h"
#include "segment2.h"
#include "plane3.h"
#include "connectivity.h"
#include <pos3.h>
#include "Space_Query/point_grid_3D.h"
#include "smooth_mesh.h"

// #include <igl/is_vertex_manifold.h>
// #include <igl/is_edge_manifold.h>
// //#include "eigen_interface.h"
// #include "Space_Query/triangle_grid_2D.h"

template <class ScalarType>
bool HasUnfererencedVert(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces)
{
    std::vector<bool> Referred(VertPos.size(), false);
    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
            Referred[Faces[i][j]] = true;

    for (size_t i = 0; i < Referred.size(); i++)
        if (!Referred[i])
            return true;

    return false;
}

template <class ScalarType>
bool HasUnfererencedVert(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces)
{
    std::vector<bool> Referred(VertPos.size(), false);
    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
            Referred[Faces[i][j]] = true;

    for (size_t i = 0; i < Referred.size(); i++)
        if (!Referred[i])
            return true;

    return false;
}

template <class ScalarType>
void RemoveUnfererencedVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                            std::vector<std::vector<int>> &Faces)
{
    std::vector<bool> Referred(VertPos.size(), false);
    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
            Referred[Faces[i][j]] = true;

    std::vector<int> NewIdx(VertPos.size(), -1);
    std::vector<Geo::Point3<ScalarType>> NewVertPos;

    int currIdx = 0;
    for (size_t i = 0; i < Referred.size(); i++)
        if (Referred[i])
        {
            NewIdx[i] = currIdx;
            currIdx++;
            NewVertPos.push_back(VertPos[i]);
        }

    VertPos = NewVertPos;
    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
        {
            int currIdx = Faces[i][j];
            assert(currIdx < NewIdx.size());
            int newIdx = NewIdx[currIdx];
            assert(newIdx >= 0);
            assert(newIdx < VertPos.size());
            Faces[i][j] = newIdx;
        }
}

template <class ScalarType>
void RemoveDuplicatedVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                          std::vector<std::vector<int>> &Faces)
{
    std::map<Geo::Point3<ScalarType>, int> PosIdxMap;
    std::vector<Geo::Point3<ScalarType>> NewVertPos;
    for (size_t i = 0; i < VertPos.size(); i++)
    {
        if (PosIdxMap.count(VertPos[i]) > 0)
            continue;
        NewVertPos.push_back(VertPos[i]);
        PosIdxMap[VertPos[i]] = NewVertPos.size() - 1;
    }

    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < 3; j++)
        {
            int IndexV = Faces[i][j];
            Geo::Point3<ScalarType> Pos = VertPos[IndexV];
            assert(PosIdxMap.count(Pos) > 0);
            int NewIdx = PosIdxMap[Pos];
            Faces[i][j] = NewIdx;
        }
    VertPos = NewVertPos;
}

template <class ScalarType>
void RemoveDuplicatedVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                          std::vector<Geo::Point2<ScalarType>> &VertUVs,
                          std::vector<std::vector<int>> &Faces)
{
    typedef std::pair<Geo::Point2<ScalarType>, Geo::Point3<ScalarType>> UV_3DKey;
    std::map<UV_3DKey, int> PosUVIdxMap;
    std::vector<Geo::Point3<ScalarType>> NewVertPos;
    std::vector<Geo::Point2<ScalarType>> NewVertUV;
    // std::cout<<"Size0:"<<VertPos.size()<<std::endl;
    assert(VertPos.size() == VertUVs.size());
    for (size_t i = 0; i < VertPos.size(); i++)
    {
        Geo::Point3<ScalarType> Pos3D = VertPos[i];
        Geo::Point2<ScalarType> PosUV = VertUVs[i];
        UV_3DKey Key(PosUV, Pos3D);

        if (PosUVIdxMap.count(Key) > 0)
            continue;

        NewVertPos.push_back(VertPos[i]);
        NewVertUV.push_back(VertUVs[i]);
        PosUVIdxMap[Key] = NewVertPos.size() - 1;
    }

    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < 3; j++)
        {
            int IndexV = Faces[i][j];
            Geo::Point3<ScalarType> Pos3D = VertPos[IndexV];
            Geo::Point2<ScalarType> PosUV = VertUVs[IndexV];
            UV_3DKey Key(PosUV, Pos3D);
            assert(PosUVIdxMap.count(Key) > 0);
            int NewIdx = PosUVIdxMap[Key];
            Faces[i][j] = NewIdx;
        }
    VertUVs = NewVertUV;
    VertPos = NewVertPos;
    // std::cout<<"Size2:"<<VertPos.size()<<std::endl;
}

template <class ScalarType>
void ClusterCloseVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                      std::vector<std::vector<int>> &Faces,
                      const ScalarType &delta)
{
    Geo::PointSpatialIndex<ScalarType> PGrid;
    PGrid.Init(VertPos);
    for (size_t i = 0; i < VertPos.size(); i++)
    {
        Geo::Point3<ScalarType> TestP = VertPos[i];
        std::vector<int> PosIn;
        PGrid.PointsInSphere(VertPos, TestP, delta, PosIn);

        if (PosIn.size() <= 1)
            continue;

        Geo::Point3<ScalarType> MinP = VertPos[PosIn[0]];
        for (size_t j = 0; j < PosIn.size(); j++)
            MinP = std::min(MinP, VertPos[PosIn[j]]);

        VertPos[i] = MinP;
    }
}

template <class ScalarType>
void MergeCloseVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    std::vector<std::vector<int>> &Faces,
                    const ScalarType &delta)
{

    // for (size_t i = 0; i < VertPos.size() - 1; i++)
    //     for (size_t j = (i + 1); j < VertPos.size(); j++)
    //     {
    //         Geo::Point3<ScalarType> Pos0 = VertPos[i];
    //         Geo::Point3<ScalarType> Pos1 = VertPos[j];
    //         if ((Pos0 - Pos1).Norm() > delta)
    //             continue;
    //         VertPos[i] = VertPos[j];
    //     }
    ClusterCloseVert(VertPos, Faces, delta);
    RemoveDuplicatedVert(VertPos, Faces);
}

template <class ScalarType>
static void GetInternalBorders(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                               const std::vector<std::vector<int>> &Faces,
                               std::vector<int> &IndexF, std::vector<int> &IndexE,
                               bool UseCloseDelta = false)
{
    IndexF.clear();
    IndexE.clear();

    std::vector<Geo::Point3<ScalarType>> VertMerged = VertPos;
    std::vector<std::vector<int>> FacesMerged = Faces;

    ScalarType delta = std::numeric_limits<ScalarType>::epsilon() * 2;
    if (UseCloseDelta)
        MergeCloseVert(VertMerged, FacesMerged, delta);
    else
        RemoveDuplicatedVert(VertMerged, FacesMerged);

    std::vector<int> BorderF_m;
    std::vector<int> BorderE_m;
    Geo::getBorderEdges(FacesMerged, BorderF_m, BorderE_m);
    assert(BorderF_m.size() == BorderE_m.size());

    typedef std::pair<int, int> FaceEdgePair;
    std::set<FaceEdgePair> BorderFaceEdgeSet_m;
    for (size_t i = 0; i < BorderE_m.size(); i++)
    {
        int IndexF = BorderF_m[i];
        int IndexE = BorderE_m[i];
        BorderFaceEdgeSet_m.insert(FaceEdgePair(IndexF, IndexE));
    }

    std::vector<int> BorderF0;
    std::vector<int> BorderE0;
    Geo::getBorderEdges(Faces, BorderF0, BorderE0);
    assert(BorderF0.size() == BorderE0.size());
    for (size_t i = 0; i < BorderF0.size(); i++)
    {
        int CurrIndexF = BorderF0[i];
        int CurrIndexE = BorderE0[i];
        FaceEdgePair Key(CurrIndexF, CurrIndexE);
        if (BorderFaceEdgeSet_m.count(Key) > 0)
            continue;
        IndexF.push_back(CurrIndexF);
        IndexE.push_back(CurrIndexE);
    }
}

void NumFaceValences(const std::vector<std::vector<int>> &Faces,
                     std::map<int, int> &ValMap)
{
    ValMap.clear();
    for (size_t i = 0; i < Faces.size(); i++)
    {
        int Valence = Faces[i].size();
        if (ValMap.find(Valence) == ValMap.end())
            ValMap[Valence] = 0;

        ValMap[Valence]++;
    }
}

bool IsAllQuad(const std::vector<std::vector<int>> &Faces)
{
    std::map<int, int> ValMap;
    NumFaceValences(Faces, ValMap);
    return (ValMap[4] == (int)Faces.size());
}

bool IsAllTris(const std::vector<std::vector<int>> &Faces)
{
    std::map<int, int> ValMap;
    NumFaceValences(Faces, ValMap);
    return (ValMap[3] == (int)Faces.size());
}

void VertValences(const std::vector<std::vector<int>> &Faces,
                  const int NumV,
                  std::vector<int> &VertValence)
{
    VertValence.clear();
    VertValence.resize(NumV, 0);
    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
        {
            int IndexV = Faces[i][j];
            assert(IndexV >= 0);
            assert(IndexV < (int)VertValence.size());
            VertValence[IndexV]++;
        }
}

void SingularVertices(const std::vector<std::vector<int>> &Faces,
                      const std::vector<bool> &BorderV,
                      std::vector<int> &SingularVert)
{
    std::vector<int> VertValence;
    VertValences(Faces, BorderV.size(), VertValence);
    for (size_t i = 0; i < VertValence.size(); i++)
    {
        if (BorderV[i] && (VertValence[i] > 2))
            SingularVert.push_back(i);
        if (!BorderV[i] && (VertValence[i] != 4))
            SingularVert.push_back(i);
    }
}

enum MType
{
    MT_Triangle,
    MT_Tetra,
    MT_Other
};

MType MeshType(const std::vector<std::vector<int>> &Elems)
{
    bool AllTris = true;
    bool AllTetra = true;
    for (size_t i = 0; i < Elems.size(); i++)
    {
        AllTris &= (Elems[i].size() == 3);
        AllTetra &= (Elems[i].size() == 4);
        if ((!AllTris) && (!AllTetra))
            break;
    }
    assert(!(AllTris && AllTetra));
    if (AllTris)
        return MT_Triangle;
    if (AllTetra)
        return MT_Tetra;
    return MT_Other;
}

template <class ScalarType>
ScalarType Area(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                const std::vector<std::vector<int>> &Faces,
                const int &IndexF)
{
    assert(Faces[IndexF].size() == 3);
    Geo::Point3<ScalarType> P0 = VertPos[Faces[IndexF][0]];
    Geo::Point3<ScalarType> P1 = VertPos[Faces[IndexF][1]];
    Geo::Point3<ScalarType> P2 = VertPos[Faces[IndexF][2]];
    ScalarType FaceArea = (Geo::Cross((P1 - P0), (P2 - P0))).Norm() / 2;
    return FaceArea;
}

template <class ScalarType>
ScalarType Area(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                const std::vector<std::vector<int>> &Faces,
                const int &IndexF, bool unsignedArea = true)
{
    assert(Faces[IndexF].size() == 3);
    Geo::Point2<ScalarType> P0 = VertPos[Faces[IndexF][0]];
    Geo::Point2<ScalarType> P1 = VertPos[Faces[IndexF][1]];
    Geo::Point2<ScalarType> P2 = VertPos[Faces[IndexF][2]];
    ScalarType FaceArea = (Geo::Cross((P1 - P0), (P2 - P0)) / 2);

    if (unsignedArea)
        return fabs(FaceArea);
    else
        return (FaceArea);
}

template <class ScalarType>
ScalarType Area(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                const std::vector<std::vector<int>> &Faces)
{
    ScalarType AreaVal = 0;
    for (size_t i = 0; i < Faces.size(); i++)
        AreaVal += Area(VertPos, Faces, i);
    return AreaVal;
}

template <class ScalarType>
ScalarType Area(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                const std::vector<std::vector<int>> &Faces,
                bool unsignedArea = true)
{
    ScalarType AreaVal = 0;
    for (size_t i = 0; i < Faces.size(); i++)
        AreaVal += Area(VertPos, Faces, i, unsignedArea);
    return AreaVal;
}

// template <class ScalarType>
// int ClosestVert(const std::vector<Geo::Point3<ScalarType> > &VertPos,
//                 const Geo::Point3<ScalarType> &TargetPos)
//{
//     ScalarType DistMin=std::numeric_limits<ScalarType>::max();
//     int IndexV=-1;
//     for (size_t i=0;i<VertPos.size();i++)
//     {
//         ScalarType TestDist=(VertPos[i]-TargetPos).Norm();
//         if (TestDist>DistMin)continue;
//         DistMin=TestDist;
//         IndexV=i;
//     }
//     return IndexV;
// }

template <class ScalarType>
static int ClosestVertFromVect(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                               const Geo::Point3<ScalarType> &TargetPos)
{
    ScalarType DistMin = std::numeric_limits<ScalarType>::max();
    int IndexV = -1;
    for (size_t i = 0; i < VertPos.size(); i++)
    {
        ScalarType TestDist = (VertPos[i] - TargetPos).Norm();
        if (TestDist > DistMin)
            continue;
        DistMin = TestDist;
        IndexV = i;
    }
    return IndexV;
}

template <class ScalarType>
ScalarType EdgeLenght(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces,
                      const int &IndexF, const int &IndexE)
{
    assert((IndexF >= 0) && (IndexF < (int)Faces.size()));
    assert((IndexE >= 0) && (IndexE < (int)Faces[IndexF].size()));

    int IndexV0 = Faces[IndexF][IndexE];
    int IndexV1 = Faces[IndexF][(IndexE + 1) % Faces[IndexF].size()];
    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
    return (P1 - P0).Norm();
}

template <class ScalarType>
ScalarType EdgeLenght(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces,
                      const int &IndexF, const int &IndexE)
{
    assert((IndexF >= 0) && (IndexF < (int)Faces.size()));
    assert((IndexE >= 0) && (IndexE < (int)Faces[IndexF].size()));

    int IndexV0 = Faces[IndexF][IndexE];
    int IndexV1 = Faces[IndexF][(IndexE + 1) % Faces[IndexF].size()];
    Geo::Point2<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point2<ScalarType> P1 = VertPos[IndexV1];
    return (P1 - P0).Norm();
}

template <class ScalarType>
ScalarType BorderLenght(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces)
{
    ScalarType BorderL=0;
    std::vector<std::pair<int,int> > BorderE;
    Geo::getBorderEdges(Faces,BorderE);
    for (size_t i=0;i<BorderE.size();i++)
    {
        int IndexV0=BorderE[i].first;
        int IndexV1=BorderE[i].second;
        Geo::Point3<ScalarType> P0=VertPos[IndexV0];
        Geo::Point3<ScalarType> P1=VertPos[IndexV1];
        BorderL+=(P0-P1).Norm();
    }
    return BorderL;
}

template <class ScalarType>
void LongestEdge(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 int &IndexF, int &IndexE)
{
    ScalarType LargestELenght = 0;
    IndexF = -1;
    IndexE = -1;
    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
        {
            //            int IndexV0=Faces[i][j];
            //            int IndexV1=Faces[i][(j+1)%Faces[i].size()];
            //            Geo::Point3<ScalarType> P0=VertPos[IndexV0];
            //            Geo::Point3<ScalarType> P1=VertPos[IndexV1];
            //            ScalarType CurrELenght=(P0-P1).Norm();

            ScalarType CurrELenght = EdgeLenght(VertPos, Faces, i, j);

            if (LargestELenght > CurrELenght)
                continue;

            LargestELenght = CurrELenght;
            IndexF = i;
            IndexE = j;
        }
}

template <class ScalarType>
ScalarType AvgEdgeLen(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces)
{
    ScalarType EdgeSum = 0;
    int Num = 0;
    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
        {
            ScalarType CurrELenght = EdgeLenght(VertPos, Faces, i, j);
            // should consider also border but become slow
            EdgeSum += CurrELenght;
            Num++;
        }
    return EdgeSum / Num;
}

template <class ScalarType>
ScalarType AvgEdgeLen(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces)
{
    ScalarType EdgeSum = 0;
    int Num = 0;
    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
        {
            ScalarType CurrELenght = EdgeLenght(VertPos, Faces, i, j);
            // should consider also border but become slow
            EdgeSum += CurrELenght;
            Num++;
        }
    return EdgeSum / Num;
}

template <class ScalarType>
ScalarType FaceVertAngleDeg(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const int &IndexF, const int &IndexV)
{
    assert((IndexF >= 0) && (IndexF < (int)Faces.size()));
    assert((IndexV >= 0) && (IndexV < (int)Faces[IndexF].size()));

    size_t NumV = Faces[IndexF].size();

    int IndexV0 = Faces[IndexF][(IndexV + NumV - 1) % NumV];
    int IndexV1 = Faces[IndexF][IndexV];
    int IndexV2 = Faces[IndexF][(IndexV + 1) % NumV];

    Geo::Point2<ScalarType> p0 = VertPos[IndexV0];
    Geo::Point2<ScalarType> p1 = VertPos[IndexV1];
    Geo::Point2<ScalarType> p2 = VertPos[IndexV2];

    Geo::Point2<ScalarType> dir0 = p0 - p1;
    Geo::Point2<ScalarType> dir1 = p2 - p1;
    dir0.Normalize();
    dir1.Normalize();
    ScalarType angle = Geo::AngleDeg(dir0, dir1);
    return angle;
}

template <class ScalarType>
ScalarType FaceVertAngleDeg(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const int &IndexF, const int &IndexV)
{
    assert((IndexF >= 0) && (IndexF < (int)Faces.size()));
    assert((IndexV >= 0) && (IndexV < (int)Faces[IndexF].size()));

    size_t NumV = Faces[IndexF].size();

    int IndexV0 = Faces[IndexF][(IndexV + NumV - 1) % NumV];
    int IndexV1 = Faces[IndexF][IndexV];
    int IndexV2 = Faces[IndexF][(IndexV + 1) % NumV];

    Geo::Point3<ScalarType> p0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> p1 = VertPos[IndexV1];
    Geo::Point3<ScalarType> p2 = VertPos[IndexV2];

    Geo::Point3<ScalarType> dir0 = p0 - p1;
    Geo::Point3<ScalarType> dir1 = p2 - p1;
    dir0.Normalize();
    dir1.Normalize();
    ScalarType angle = Geo::AngleDeg(dir0, dir1);
    return angle;
}

template <class ScalarType>
void InterpolatePerFace(const std::vector<Geo::Point3<ScalarType>> &VertVal,
                        const std::vector<std::vector<int>> &Faces,
                        std::vector<Geo::Point3<ScalarType>> &FaceVal)
{
    FaceVal.clear();
    FaceVal.resize(Faces.size(), Geo::Point3<ScalarType>(0, 0, 0));
    for (size_t i = 0; i < Faces.size(); i++)
    {
        for (size_t j = 0; j < Faces[i].size(); j++)
        {
            int IndexV = Faces[i][j];
            assert(IndexV < (int)VertVal.size());
            FaceVal[i] += VertVal[IndexV];
        }
        FaceVal[i] /= Faces[i].size();
    }
}

template <class ScalarType>
void SplitQuadIntoTris(const std::vector<Geo::Point3<ScalarType>> &Pos,
                       const std::vector<std::vector<int>> &Quads,
                       std::vector<std::vector<int>> &Tris)
{
    Tris.clear();
    // 2 triangles for each quad
    Tris.resize((int)Quads.size() * 2, std::vector<int>(3, -1));
    for (size_t i = 0; i < Quads.size(); i++)
    {
        assert(Quads[i].size() == 4);
        int IndexV0 = Quads[i][0];
        int IndexV1 = Quads[i][1];
        int IndexV2 = Quads[i][2];
        int IndexV3 = Quads[i][3];

        Geo::Point3<ScalarType> PosV0 = Pos[IndexV0];
        Geo::Point3<ScalarType> PosV1 = Pos[IndexV1];
        Geo::Point3<ScalarType> PosV2 = Pos[IndexV2];
        Geo::Point3<ScalarType> PosV3 = Pos[IndexV3];

        // split along shortest distance
        ScalarType D0 = (PosV0 - PosV2).Norm();
        ScalarType D1 = (PosV1 - PosV3).Norm();

        int IndexTriF0 = i * 2;
        int IndexTriF1 = (i * 2) + 1;
        if (D0 < D1)
        {
            Tris[IndexTriF0][0] = IndexV0;
            Tris[IndexTriF0][1] = IndexV1;
            Tris[IndexTriF0][2] = IndexV2;
            Tris[IndexTriF1][0] = IndexV0;
            Tris[IndexTriF1][1] = IndexV2;
            Tris[IndexTriF1][2] = IndexV3;
        }
        else
        {
            Tris[IndexTriF0][0] = IndexV0;
            Tris[IndexTriF0][1] = IndexV1;
            Tris[IndexTriF0][2] = IndexV3;
            Tris[IndexTriF1][0] = IndexV1;
            Tris[IndexTriF1][1] = IndexV2;
            Tris[IndexTriF1][2] = IndexV3;
        }
    }
}

void FFConnectedComponenent(const std::vector<std::vector<int>> &Elements,
                            const std::vector<std::vector<int>> &NextP,
                            const std::vector<std::vector<int>> &NextE,
                            std::vector<std::vector<int>> &ComponentsF)
{
    ComponentsF.clear();

    // find partitions
    std::vector<bool> Has_Explored(Elements.size(), false);

    for (size_t i = 0; i < Elements.size(); i++)
    {
        if (Has_Explored[i])
            continue;

        ComponentsF.resize(ComponentsF.size() + 1);

        std::vector<size_t> StackP;
        StackP.push_back(i);
        do
        {
            int CurrP = StackP.back();
            StackP.pop_back();

            if (Has_Explored[CurrP])
                continue;

            Has_Explored[CurrP] = true;
            ComponentsF.back().push_back(CurrP);

            for (size_t n = 0; n < NextP[CurrP].size(); n++)
            {
                // get next elem
                int AdjP = NextP[CurrP][n];
                assert(AdjP != CurrP);

                // get edge as seen from the other face
                int AdjE = NextE[CurrP][n];

                if (AdjP == -1)
                    continue;
                assert(AdjE != -1);

                bool has_to_diffuse = (!Has_Explored[AdjP]);
                if (has_to_diffuse)
                {
                    StackP.push_back(AdjP);
                }
            }
        } while (StackP.size() > 0);
    }
}

void FFConnectedComponenent(const std::vector<std::vector<int>> &Elements,
                            const std::vector<std::vector<int>> &NextP,
                            const std::vector<std::vector<int>> &NextE,
                            std::vector<std::vector<int>> &ComponentsF,
                            std::vector<std::vector<int>> &ComponentsV)
{
    FFConnectedComponenent(Elements, NextP, NextE, ComponentsF);
    ComponentsV.resize(ComponentsF.size());

    for (size_t c = 0; c < ComponentsF.size(); c++)
    {
        for (size_t i = 0; i < ComponentsF[c].size(); i++)
        {
            int IndexF = ComponentsF[c][i];
            assert(IndexF >= 0);
            assert(IndexF < (int)Elements.size());
            for (size_t j = 0; j < Elements[IndexF].size(); j++)
                ComponentsV[c].push_back(Elements[IndexF][j]);
        }

        std::sort(ComponentsV[c].begin(), ComponentsV[c].end());
        std::vector<int>::iterator last = std::unique(ComponentsV[c].begin(), ComponentsV[c].end());
        ComponentsV[c].erase(last, ComponentsV[c].end());
    }
}

void FFConnectedComponenent(const std::vector<std::vector<int>> &Elements,
                            std::vector<std::vector<int>> &ComponentsF,
                            std::vector<std::vector<int>> &ComponentsV)
{
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Elements, NextF, NextE);
    FFConnectedComponenent(Elements, NextF, NextE, ComponentsF, ComponentsV);
}

template <class ScalarType>
void AppendMesh(std::vector<Geo::Point3<ScalarType>> &Vert,
                std::vector<std::vector<int>> &Elements,
                const std::vector<Geo::Point3<ScalarType>> &To_append_Vert,
                const std::vector<std::vector<int>> &To_append_Elements)
{
    int NumV0 = Vert.size();
    Vert.insert(Vert.end(), To_append_Vert.begin(), To_append_Vert.end());

    int NumF0 = Elements.size();
    Elements.insert(Elements.end(), To_append_Elements.begin(), To_append_Elements.end());

    for (size_t i = NumF0; i < Elements.size(); i++)
        for (size_t j = 0; j < Elements[i].size(); j++)
            Elements[i][j] += NumV0;
}

template <class ScalarType>
void CopySubMesh(const std::vector<Geo::Point3<ScalarType>> &Vert,
                 const std::vector<std::vector<int>> &Elements,
                 const std::vector<int> &ToCopyEl,
                 std::vector<Geo::Point3<ScalarType>> &CopiedMeshesVert,
                 std::vector<std::vector<int>> &CopiedMeshesElem,
                 std::vector<int> &OriginalVert,
                 std::vector<int> &OriginalElem)
{
    OriginalVert.clear();
    OriginalElem.clear();
    CopiedMeshesVert.clear();
    CopiedMeshesElem.clear();

    // copy faces
    for (size_t i = 0; i < ToCopyEl.size(); i++)
    {
        int IndexEl = ToCopyEl[i];
        assert(IndexEl >= 0);
        assert(IndexEl < Elements.size());

        CopiedMeshesElem.push_back(Elements[IndexEl]);
        OriginalElem.push_back(IndexEl);
    }

    // then copy/reassign vertices
    std::map<int, int> VertToVert;
    for (size_t i = 0; i < CopiedMeshesElem.size(); i++)
        for (size_t j = 0; j < CopiedMeshesElem[i].size(); j++)
        {
            int IndexV = CopiedMeshesElem[i][j];
            assert(IndexV >= 0);
            assert(IndexV < Vert.size());
            if (VertToVert.count(IndexV) == 0)
            {
                CopiedMeshesVert.push_back(Vert[IndexV]);
                OriginalVert.push_back(IndexV);
                VertToVert[IndexV] = CopiedMeshesVert.size() - 1;
            }
            int NewIndexV = VertToVert[IndexV];
            CopiedMeshesElem[i][j] = NewIndexV;
        }
}

template <class ScalarType>
void SplitByFFConnectedComponenent(const std::vector<Geo::Point3<ScalarType>> &Vert,
                                   const std::vector<std::vector<int>> &Elements,
                                   const std::vector<std::vector<int>> &NextP,
                                   const std::vector<std::vector<int>> &NextE,
                                   std::vector<std::vector<Geo::Point3<ScalarType>>> &SubMeshesV,
                                   std::vector<std::vector<std::vector<int>>> &SubMeshesF,
                                   std::vector<std::vector<int>> &OriginalV,
                                   std::vector<std::vector<int>> &OriginalF)
{
    std::vector<std::vector<int>> Components;
    FFConnectedComponenent(Elements, NextP, NextE, Components);

    SubMeshesV.clear();
    SubMeshesV.resize(Components.size());
    SubMeshesF.clear();
    SubMeshesF.resize(Components.size());

    OriginalV.clear();
    OriginalV.resize(Components.size());
    OriginalF.clear();
    OriginalF.resize(Components.size());
    for (size_t i = 0; i < Components.size(); i++)
        CopySubMesh(Vert, Elements, Components[i],
                    SubMeshesV[i], SubMeshesF[i],
                    OriginalV[i], OriginalF[i]);
}

template <class ScalarType>
void SplitByFFConnectedComponenent(const std::vector<Geo::Point3<ScalarType>> &Vert,
                                   const std::vector<std::vector<int>> &Elements,
                                   std::vector<std::vector<Geo::Point3<ScalarType>>> &SubMeshesV,
                                   std::vector<std::vector<std::vector<int>>> &SubMeshesF,
                                   std::vector<std::vector<int>> &OriginalV,
                                   std::vector<std::vector<int>> &OriginalF)
{
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Elements, NextF, NextE);
    SplitByFFConnectedComponenent(Vert, Elements, NextF, NextE, SubMeshesV,
                                  SubMeshesF, OriginalV, OriginalF);
}

template <class ScalarType>
void SplitByFFConnectedComponenent(const std::vector<Geo::Point2<ScalarType>> &Vert,
                                   const std::vector<std::vector<int>> &Elements,
                                   std::vector<std::vector<Geo::Point2<ScalarType>>> &SubMeshesV,
                                   std::vector<std::vector<std::vector<int>>> &SubMeshesF,
                                   std::vector<std::vector<int>> &OriginalV,
                                   std::vector<std::vector<int>> &OriginalF)
{
    std::vector<Geo::Point3<ScalarType>> Vert3D;
    for (size_t i = 0; i < Vert.size(); i++)
        Vert3D.push_back(Geo::Point3<ScalarType>(Vert[i].X, Vert[i].Y, 0));

    std::vector<std::vector<Geo::Point3<ScalarType>>> SubMeshesV3D;
    SplitByFFConnectedComponenent(Vert3D, Elements, SubMeshesV3D, SubMeshesF, OriginalV, OriginalF);

    SubMeshesV.clear();
    SubMeshesV.resize(SubMeshesV3D.size());
    for (size_t i = 0; i < SubMeshesV3D.size(); i++)
        for (size_t j = 0; j < SubMeshesV3D[i].size(); j++)
        {
            Geo::Point3<ScalarType> Pos = SubMeshesV3D[i][j];
            SubMeshesV[i].push_back(Geo::Point2<ScalarType>(Pos.X, Pos.Y));
        }
}

int PolyMeshGenus(int NVert, int nEdges, int NFaces, int NHoles, int NComponents)
{
    return -((NVert + NFaces - nEdges + NHoles - 2 * NComponents) / 2);
}

int NumEdges(const std::vector<std::vector<int>> &Connectivity)
{
    std::set<std::pair<int, int>> EdgeSet;
    for (size_t i = 0; i < Connectivity.size(); i++)
    {
        int sizeF = Connectivity[i].size();
        for (int j = 0; j < sizeF; j++)
        {
            int IndexV0 = Connectivity[i][j];
            int IndexV1 = Connectivity[i][(j + 1) % sizeF];
            std::pair<int, int> entry(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
            EdgeSet.insert(entry);
        }
    }
    return EdgeSet.size();
}

template <class ScalarType>
void ComputePerVertAngleDeg(const std::vector<Geo::Point3<ScalarType>> &Pos,
                            const std::vector<std::vector<int>> &Faces,
                            std::vector<ScalarType> &VertAngle)
{
    VertAngle.clear();
    VertAngle.resize(Pos.size(), 0);
    for (size_t i = 0; i < Faces.size(); i++)
    {
        size_t FaceS = Faces[i].size();
        for (size_t j = 0; j < Faces[i].size(); j++)
        {
            int IndexLoc0 = (j + FaceS - 1) % FaceS;
            int IndexLoc1 = j;
            int IndexLoc2 = (j + 1) % FaceS;

            int IndexV0 = Faces[i][IndexLoc0];
            int IndexV1 = Faces[i][IndexLoc1];
            int IndexV2 = Faces[i][IndexLoc2];

            Geo::Point3<ScalarType> P0 = Pos[IndexV0];
            Geo::Point3<ScalarType> P1 = Pos[IndexV1];
            Geo::Point3<ScalarType> P2 = Pos[IndexV2];

            Geo::Point3<ScalarType> Dir0 = P0 - P1;
            Geo::Point3<ScalarType> Dir1 = P2 - P1;
            Dir0.Normalize();
            Dir1.Normalize();
            VertAngle[IndexV1] += Geo::AngleDeg(Dir0, Dir1);
        }
    }
}

template <class ScalarType>
void ComputeFaceBarycenters(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            std::vector<Geo::Point3<ScalarType>> &BaryF)
{
    BaryF.clear();
    // go over all faces
    for (size_t i = 0; i < Faces.size(); i++)
    {
        size_t sizeF = Faces[i].size();

        // get vertex barycenter
        Geo::Point3<ScalarType> BaryVal(0, 0, 0);
        for (size_t j = 0; j < sizeF; j++)
        {
            int IndexV = Faces[i][j];
            BaryVal += VertPos[IndexV];
        }
        BaryVal /= sizeF;
        BaryF.push_back(BaryVal);
    }
}

template <class ScalarType>
Geo::Point3<ScalarType> InterpolatePosTriMesh(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                                              const std::vector<std::vector<int>> &Faces,
                                              const int &FaceIdx,
                                              const Geo::Point3<ScalarType> &bary)
{
    assert(FaceIdx >= 0);
    assert(FaceIdx < Faces.size());
    int IndexV0 = Faces[FaceIdx][0];
    int IndexV1 = Faces[FaceIdx][1];
    int IndexV2 = Faces[FaceIdx][2];

    Geo::Point3<ScalarType> FaceP0 = Pos3D[IndexV0];
    Geo::Point3<ScalarType> FaceP1 = Pos3D[IndexV1];
    Geo::Point3<ScalarType> FaceP2 = Pos3D[IndexV2];
    Geo::Point3<ScalarType> InterpP = FaceP0 * bary.X +
                                      FaceP1 * bary.Y +
                                      FaceP2 * bary.Z;
    return InterpP;
}

template <class ScalarType>
Geo::Point2<ScalarType> InterpolateUVTriMesh(const std::vector<Geo::Point2<ScalarType>> &PosUV,
                                             const std::vector<std::vector<int>> &Faces,
                                             const int &FaceIdx,
                                             const Geo::Point3<ScalarType> &bary)
{
    assert(FaceIdx >= 0);
    assert(FaceIdx < Faces.size());
    int IndexV0 = Faces[FaceIdx][0];
    int IndexV1 = Faces[FaceIdx][1];
    int IndexV2 = Faces[FaceIdx][2];

    Geo::Point2<ScalarType> FaceP0 = PosUV[IndexV0];
    Geo::Point2<ScalarType> FaceP1 = PosUV[IndexV1];
    Geo::Point2<ScalarType> FaceP2 = PosUV[IndexV2];
    Geo::Point2<ScalarType> InterpP = FaceP0 * bary.X +
                                      FaceP1 * bary.Y +
                                      FaceP2 * bary.Z;
    return InterpP;
}

template <class ScalarType>
void SplitMeshFromEdges(std::vector<Geo::Point3<ScalarType>> &Pos,
                        std::vector<std::vector<int>> &Faces,
                        const std::vector<std::vector<int>> &NextF,
                        const std::vector<std::vector<int>> &NextE,
                        const std::vector<std::pair<int, int>> &SplitEdges)
{
    std::set<std::pair<int, int>> SplitEdgeSet(SplitEdges.begin(), SplitEdges.end());
    //    std::vector<std::vector<int> > NextF,NextE;
    //    Geo::getFFAdjacency(Faces, NextF,NextE);

    std::vector<Geo::Point3<ScalarType>> NewPos = Pos;
    std::vector<std::vector<int>> NewFaces = Faces;

    std::vector<bool> CheckedV(Pos.size(), false);

    for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++)
        {
            // if border do nothing
            if (NextF[i][j] == -1)
                continue;

            int IndexV0 = Faces[i][j];
            // already checked

            int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
            std::pair<int, int> Entry0(IndexV0, IndexV1);
            std::pair<int, int> Entry1(IndexV1, IndexV0);
            bool HasSplit = ((SplitEdgeSet.count(Entry0) > 0) ||
                             (SplitEdgeSet.count(Entry1) > 0));

            // no split
            if (!HasSplit)
                continue;

            if (CheckedV[IndexV0])
                continue;
            CheckedV[IndexV0] = true;

            Pos3 StartP(i, j, IndexV0);
            std::vector<Pos3> PosFan;
            StartP.FlipE(Faces);
            bool HasBorder;

            GetSortedFacesFaceFun(StartP, Faces, NextF, NextE, PosFan, HasBorder);

            std::vector<std::vector<Pos3>> PosPartitions;
            PosPartitions.resize(1);
            for (size_t k = 0; k < PosFan.size(); k++)
            {
                PosPartitions.back().push_back(PosFan[k]);

                // do not check the last one because moight be back on first marked edge
                if ((!HasBorder) && (k == PosFan.size() - 1))
                    continue;

                // check if has a border on other side
                Pos3 TestPos = PosFan[k];
                // TestPos.FlipE(Faces);
                int IndexVTest0 = TestPos.V;
                TestPos.FlipV(Faces);
                int IndexVTest1 = TestPos.V;
                std::pair<int, int> Entry0(IndexVTest0, IndexVTest1);
                std::pair<int, int> Entry1(IndexVTest1, IndexVTest0);
                bool HasSplit = ((SplitEdgeSet.count(Entry0) > 0) ||
                                 (SplitEdgeSet.count(Entry1) > 0));

                if (HasSplit)
                    PosPartitions.resize(PosPartitions.size() + 1);
            }

            // no need to split the vertex
            if (PosPartitions.size() == 1)
                continue;

            Geo::Point3<ScalarType> VertPos = Pos[StartP.V];
            // first partition remain the same then from the secon to other change vertex
            for (size_t k = 1; k < PosPartitions.size(); k++)
            {
                NewPos.push_back(VertPos);
                int newVIndex = NewPos.size() - 1;
                for (size_t w = 0; w < PosPartitions[k].size(); w++)
                {
                    int IndexF = PosPartitions[k][w].F;
                    int IndexV = PosPartitions[k][w].V;
                    int LocV = Geo::WhichVert(Faces[IndexF], IndexV);
                    assert(LocV >= 0);
                    assert(LocV <= (int)Faces[IndexF].size());
                    NewFaces[IndexF][LocV] = newVIndex;
                }
            }
        }
    Pos = NewPos;
    Faces = NewFaces;
}

template <class ScalarType>
static void RemoveFacesOnNegativeSide(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                                      std::vector<std::vector<int>> &Faces,
                                      const Geo::Plane3<ScalarType> &SplitPl)
{
    std::vector<std::vector<int>> NewFaces;
    for (size_t i = 0; i < Faces.size(); i++)
    {
        Geo::Point3<ScalarType> Bary(0, 0, 0);
        for (size_t j = 0; j < Faces[i].size(); j++)
            Bary += Pos3D[Faces[i][j]];
        Bary /= Faces[i].size();
        if (SplitPl.SignedDistance(Bary) < 0)
            continue;
        NewFaces.push_back(Faces[i]);
    }
    Faces = NewFaces;
    RemoveUnfererencedVert(Pos3D, Faces);
}

template <class ScalarType>
static void MirrorMesh(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                       std::vector<std::vector<int>> &Faces,
                       const Geo::Plane3<ScalarType> &SplitPl,
                       const std::set<int> &MiddleV)
//,bool SnapClosePlane=true)
{
    std::vector<Geo::Point3<ScalarType>> MirroredPos3D = Pos3D;
    std::vector<std::vector<int>> MirroredFaces = Faces;
    std::vector<bool> IsB;
    Geo::getBorderVert(Faces, IsB);
    // ScalarType delta = 0.0000001;
    for (size_t i = 0; i < MirroredPos3D.size(); i++)
    {
        // no need to mirror on the place
        if (MiddleV.count(i) > 0)
            continue;
        ScalarType Dist = fabs(SplitPl.SignedDistance(MirroredPos3D[i]));
        // if ((IsB[i])&&(Dist<delta)&&(SnapClosePlane))continue;
        // ScalarType delta = std::numeric_limits<ScalarType>::epsilon() * 2;
        MirroredPos3D[i] = SplitPl.Mirror(MirroredPos3D[i]);
    }

    // invert the orientation of faces
    for (size_t i = 0; i < MirroredFaces.size(); i++)
        std::reverse(MirroredFaces[i].begin(), MirroredFaces[i].end());

    AppendMesh(Pos3D, Faces, MirroredPos3D, MirroredFaces);
}

template <class ScalarType>
static void RemoveDeleted(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                          std::vector<std::vector<int>> &Connectivity,
                          const std::vector<bool> &DeletedV,
                          const std::vector<bool> &DeletedEl)
{
    assert(Pos3D.size() == DeletedV.size());
    assert(Connectivity.size() == DeletedEl.size());

    std::vector<Geo::Point3<ScalarType>> NewPos3D;
    std::vector<int> NewIdx(Pos3D.size(), -1);
    for (size_t i = 0; i < Pos3D.size(); i++)
    {
        if (DeletedV[i])
            continue;
        NewPos3D.push_back(Pos3D[i]);
        NewIdx[i] = NewPos3D.size() - 1;
    }

    std::vector<std::vector<int>> NewConnectivity;
    for (size_t i = 0; i < Connectivity.size(); i++)
    {
        if (DeletedEl[i])
            continue;
        NewConnectivity.resize(NewConnectivity.size() + 1);
        for (size_t j = 0; j < Connectivity[i].size(); j++)
        {
            int IndexVOr = Connectivity[i][j];
            assert(IndexVOr >= 0);
            //            std::cout<<"IndexVOr "<<IndexVOr<<std::endl;
            //            std::cout<<"NewIdx.size() "<<NewIdx.size()<<std::endl;
            assert(IndexVOr < NewIdx.size());
            int IndexVNew = NewIdx[IndexVOr];
            assert(IndexVNew >= 0);
            assert(IndexVNew < NewPos3D.size());
            NewConnectivity.back().push_back(IndexVNew);
        }
    }
}

static void InvertOrientation(std::vector<std::vector<int>> &Faces)
{
    for (size_t i = 0; i < Faces.size(); i++)
        std::reverse(Faces[i].begin(), Faces[i].end());
}

template <class ScalarType>
static bool HasNaNVert(const std::vector<Geo::Point3<ScalarType>> &VertPos)
{
    for (size_t i = 0; i < VertPos.size(); i++)
    {
        if (std::isnan(VertPos[i].X))
            return true;
        if (std::isnan(VertPos[i].Y))
            return true;
        if (std::isnan(VertPos[i].Z))
            return true;
    }
    return false;
}

static void GetDegenerateConnectivityTriangles(const std::vector<std::vector<int>> &Connectivity,
                                               std::vector<int> &DegenerateF)
{
    DegenerateF.clear();
    for (size_t i = 0; i < Connectivity.size(); i++)
    {
        assert(Connectivity[i].size() == 3);
        std::set<int> FaceVSet(Connectivity[i].begin(), Connectivity[i].end());
        if (FaceVSet.size() == Connectivity[i].size())
            continue;
        DegenerateF.push_back(i);
    }
}

template <class ScalarType>
static void RemoveDegenerateConnectivityTriangles(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                                  std::vector<std::vector<int>> &Connectivity,
                                                  std::vector<Geo::Point3<ScalarType>> &FaceNormals)
{
    std::vector<int> DegenerateF;
    GetDegenerateConnectivityTriangles(Connectivity, DegenerateF);
    if (DegenerateF.size() == 0)
        return;

    std::set<int> DegenerateSet(DegenerateF.begin(), DegenerateF.end());
    std::vector<std::vector<int>> NewF;
    std::vector<Geo::Point3<ScalarType>> NewNorm;
    for (size_t i = 0; i < Connectivity.size(); i++)
    {
        if (DegenerateSet.count(i) > 0)
            continue;

        NewF.push_back(Connectivity[i]);
        NewNorm.push_back(FaceNormals[i]);
    }
    FaceNormals = NewNorm;
    Connectivity = NewF;
    RemoveUnfererencedVert(VertPos, Connectivity);
}

template <class ScalarType>
static void GetDegenerateTriangles(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                   const std::vector<std::vector<int>> &Connectivity,
                                   const ScalarType &angleDelta,
                                   std::vector<int> &DegenerateF)
{
    for (size_t i = 0; i < Connectivity.size(); i++)
    {
        assert(Connectivity[i].size() == 3);
        for (size_t j = 0; j < 3; j++)
        {
            int IndexV = Connectivity[i][j];
            int IndexV0 = Connectivity[i][(j + 2) % 3];
            int IndexV1 = Connectivity[i][(j + 1) % 3];
            Geo::Point3<ScalarType> P = VertPos[IndexV];
            Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
            Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
            Geo::Point3<ScalarType> Dir0 = P0 - P;
            Geo::Point3<ScalarType> Dir1 = P1 - P;
            Dir0.Normalize();
            Dir1.Normalize();
            ScalarType Angle = Geo::AngleDeg(Dir0, Dir1);
            if ((std::isnan(Angle)) ||
                (Angle < angleDelta) ||
                (Angle > (180 - angleDelta)))
            {
                DegenerateF.push_back(i);
                break;
            }
        }
    }
}

template <class ScalarType>
static void SmoothDegerateTriangles(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                    const std::vector<std::vector<int>> &Connectivity,
                                    const ScalarType &angleDelta)
{
    std::vector<int> DegenerateF;
    GetDegenerateTriangles(VertPos, Connectivity, angleDelta, DegenerateF);
    std::vector<bool> IsFix = std::vector<bool>(VertPos.size(), true);
    for (size_t i = 0; i < DegenerateF.size(); i++)
    {
        int IndexF = DegenerateF[i];
        for (size_t j = 0; j < 3; j++)
        {
            int IndexV = Connectivity[IndexF][j];
            IsFix[IndexV] = false;
        }
    }
    std::vector<int> FixedVert;
    int To_Smooth = 0;
    for (size_t i = 0; i < IsFix.size(); i++)
    {
        if (!IsFix[i])
        {
            To_Smooth++;
            continue;
        }
        FixedVert.push_back(i);
    }
    if (To_Smooth == 0)
        return;
    std::cout << "**have to smooth:" << To_Smooth << std::endl;
    const std::vector<std::pair<int, int>> Features;
    Geo::SmoothLaplacian<ScalarType>(VertPos, Connectivity, 1,
                                     FixedVert, Features, false);
}
// if (cleanF)
// {
//     std::vector<bool> IsFix = std::vector<bool>(SubMeshesV[i].size(), true);
//     ScalarType delta = AvgEdgeLen(SubMeshesV[i], SubMeshesF[i]) * 0.0001;
//     for (size_t j = 0; j < SubMeshesF[i].size(); j++)
//     {
//         for (size_t k = 0; k < 3; k++)
//         {
//             int IndexV0 = SubMeshesF[i][j][k];
//             int IndexV1 = SubMeshesF[i][j][(k + 1) % 3];
//             Geo::Point3<ScalarType> P0 = SubMeshesV[i][IndexV0];
//             Geo::Point3<ScalarType> P1 = SubMeshesV[i][IndexV1];
//             if ((P0 - P1).Norm() > delta)
//                 continue;
//             IsFix[IndexV0] = false;
//             IsFix[IndexV1] = false;
//         }
//     }
//     std::vector<int> FixedVert;
//     for (size_t j = 0; j < IsFix.size(); j++)
//     {
//         if (!IsFix[i])continue;
//         FixedVert.push_back(j);
//     }
//     SmoothLaplacian<ScalarType>(SubMeshesV[i],SubMeshesF[i],10,FixedVert);

// template <class ScalarType>
// static bool EdgeManifold(std::vector<Geo::Point3<ScalarType> > &Pos3D,
//                         std::vector<std::vector<int> > &Faces)
// {
//     Eigen::MatrixXi F;
//     typename EigenInterface<ScalarType>::MatrixXm Vf;

//     EigenInterface<ScalarType>::GetTriMeshData(Pos3D, Faces, F, Vf);
//     Eigen::MatrixXd V = Vf.template cast<double>();

//      bool is_manifold = igl::is_edge_manifold(F);
//      return is_manifold;
// }

// template <class ScalarType>
// static bool VertManifold(std::vector<Geo::Point3<ScalarType> > &Pos3D,
//                         std::vector<std::vector<int> > &Faces)
// {
//     Eigen::MatrixXi F;
//     typename EigenInterface<ScalarType>::MatrixXm Vf;

//     EigenInterface<ScalarType>::GetTriMeshData(Pos3D, Faces, F, Vf);
//     Eigen::MatrixXd V = Vf.template cast<double>();

//     bool is_manifold = igl::is_vertex_manifold(F);
//     return is_manifold;
// }

template <class ScalarType>
void PosInstances(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                  std::vector<int> &Istances)
{
    Istances = std::vector<int>(Pos3D.size(),0);
    std::map<Geo::Point3<ScalarType>, int> PosCount;
    for (size_t i = 0; i < Pos3D.size(); i++)
        PosCount[Pos3D[i]]++;
    for (size_t i = 0; i < Pos3D.size(); i++)
    {
        assert(PosCount.count(Pos3D[i]) > 0);
        Istances[i]=PosCount[Pos3D[i]];
    }
}

// template <class ScalarType>
// void DesambiguateSamePosVertices(std::vector<Geo::Point3<ScalarType>> &Pos3D,
//                                  const std::vector<std::vector<int>> &Faces,
//                                  ScalarType delta)
// {
//     // first find vertices with double positions
//     std::map<Geo::Point3<ScalarType>, std::vector<int>> PosToIdx;
//     for (size_t i = 0; i < Pos3D.size(); i++)
//         PosToIdx[Pos3D[i]].push_back(i);

//     // then compute laplacian smooth
//     std::vector<Geo::Point3<ScalarType>> OldPos = Pos3D;
//     std::vector<std::pair<int,int> > Features;
//     SmoothLaplacian<ScalarType>(Pos3D, Faces, 5,Features,false);
//     for (size_t i = 0; i < OldPos.size(); i++)
//     {
//         assert(PosToIdx.count(OldPos[i]) > 0);
//         if (PosToIdx[OldPos[i]].size() == 2)
//         {
//             Geo::Point3<ScalarType> DirMov = Pos3D[i] - OldPos[i];
//             DirMov.Normalize();
//             DirMov = DirMov * delta;
//             Pos3D[i] = OldPos[i] + DirMov;
//         }
//         else
//             Pos3D[i] = OldPos[i];
//     }
// }

template <class ScalarType>
int DesambiguateSamePosVertices(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                                 const std::vector<std::vector<int>> &Faces,
                                 ScalarType delta)
{
    // first find vertices with double positions
    std::map<Geo::Point3<ScalarType>, std::vector<int>> PosToIdx;
    for (size_t i = 0; i < Pos3D.size(); i++)
        PosToIdx[Pos3D[i]].push_back(i);

    std::vector<int> FixedV;
    for (size_t i = 0; i < Pos3D.size(); i++)
    {
        if (PosToIdx[Pos3D[i]].size() == 1)
            FixedV.push_back(i);
    }
    // then compute laplacian smooth
    std::vector<Geo::Point3<ScalarType>> OldPos = Pos3D;
    std::vector<std::pair<int,int> > Features;

    SmoothLaplacian<ScalarType>(Pos3D, Faces, 5,FixedV,Features,false);
    int Moved=0;
    for (size_t i = 0; i < OldPos.size(); i++)
    {
        assert(PosToIdx.count(OldPos[i]) > 0);
        if (PosToIdx[OldPos[i]].size() == 2)
        {
            Geo::Point3<ScalarType> DirMov = Pos3D[i] - OldPos[i];
            DirMov.Normalize();
            DirMov = DirMov * delta;
            Pos3D[i] = OldPos[i] + DirMov;
            Moved++;
        }
        else
            Pos3D[i] = OldPos[i];
    }
    return Moved;
}

#endif
