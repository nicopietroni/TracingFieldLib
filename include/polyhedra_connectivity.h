#ifndef POLYHEDRA_CONNECTIVITY
#define POLYHEDRA_CONNECTIVITY

#include "box3.h"
#include "point3.h"
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
// #include "connectivity.h"
// #include "normal.h"

namespace Tetra {

// Vertex of a Edge
static int VofE(const int &indexE, const int &indexV) {
  assert((indexE < 6) && (indexV < 2));
  static int edgevert[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  return (edgevert[indexE][indexV]);
}

// Vertex of a Face
static int VofF(const int &indexF, const int &indexV) {
  assert((indexF < 4) && (indexV < 3));
  static int facevert[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
  return (facevert[indexF][indexV]);
}

// static int EofV(const int &indexV, const int &indexE)
//{
//     assert((indexE < 3) && (indexV < 4));
//     static int vertedge[4][3] = {{0, 1, 2},
//                                  {0, 3, 4},
//                                  {5, 1, 3},
//                                  {4, 5, 2}};
//     return vertedge[indexV][indexE];
// }

// edge of a Face
static int EofF(const int &indexF, const int &indexE) {
  assert((indexF < 4) && (indexE < 3));
  static int faceedge[4][3] = {{0, 3, 1}, {2, 4, 0}, {1, 5, 2}, {4, 5, 3}};
  return faceedge[indexF][indexE];
}

// static int FofV(const int &indexV, const int &indexF)
//{
//     assert((indexV < 4) && (indexF < 3));
//     static int vertface[4][3] = {{0, 1, 2},
//                                  {0, 3, 1},
//                                  {0, 2, 3},
//                                  {1, 3, 2}};
//     return vertface[indexV][indexF];
// }

static int FofE(const int &indexE, const int &indexSide) {
  assert((indexE < 6) && (indexSide < 2));
  static int edgeface[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 3}, {2, 3}};
  return edgeface[indexE][indexSide];
}

// return the vertex between two edges
static int VofEE(const int &indexE0, const int &indexE1) {
  assert((indexE0 < 6) && (indexE0 >= 0));
  assert((indexE1 < 6) && (indexE1 >= 0));
  static int edgesvert[6][6] = {{-1, 0, 0, 1, 1, -1}, {0, -1, 0, 2, -1, 2},
                                {0, 0, -1, -1, 3, 3}, {1, 2, -1, -1, 1, 2},
                                {1, -1, 3, 1, -1, 3}, {-1, 2, 3, 2, 3, -1}};
  return (edgesvert[indexE0][indexE1]);
}

// static int VofFFF(const int &indexF0, const int &indexF1, const int &indexF2)
//{
//     assert((indexF0 < 4) && (indexF0 >= 0));
//     assert((indexF1 < 4) && (indexF1 >= 0));
//     assert((indexF2 < 4) && (indexF2 >= 0));
//     static int facesvert[4][4][4] = {
//         {//0
//          {-1, -1, -1, -1},
//          {-1, -1, 0, 1},
//          {-1, 0, -1, 2},
//          {-1, 1, 2, -1}},
//         {//1
//          {-1, -1, 0, 1},
//          {-1, -1, -1, -1},
//          {0, -1, -1, 3},
//          {1, -1, 3, -1}},
//         {//2
//          {-1, 0, -1, 2},
//          {0, -1, -1, 3},
//          {-1, -1, -1, -1},
//          {2, 3, -1, -1}},
//         {//3
//          {-1, 1, 2, -1},
//          {1, -1, 3, -1},
//          {2, 3, -1, -1},
//          {-1, -1, -1, -1}}};
//     return facesvert[indexF0][indexF1][indexF2];
// }

// static int EofFF(const int &indexF0, const int &indexF1)
//{
//     assert((indexF0 < 4) && (indexF0 >= 0));
//     assert((indexF1 < 4) && (indexF1 >= 0));
//     static int facesedge[4][4] = {{-1, 0, 1, 3},
//                                   {0, -1, 2, 4},
//                                   {1, 2, -1, 5},
//                                   {3, 4, 5, -1}};
//     return (facesedge[indexF0][indexF1]);
// }

// static int EofVV(const int &indexV0, const int &indexV1)
//{
//     assert((indexV0 < 4) && (indexV0 >= 0));
//     assert((indexV1 < 4) && (indexV1 >= 0));
//     static int verticesedge[4][4] = {{-1, 0, 1, 2},
//                                      {0, -1, 3, 4},
//                                      {1, 3, -1, 5},
//                                      {2, 4, 5, -1}};

//    return verticesedge[indexV0][indexV1];
//}

// static int FofVVV(const int &indexV0, const int &indexV1, const int &indexV2)
//{

//    assert((indexV0 < 4) && (indexV0 >= 0));
//    assert((indexV1 < 4) && (indexV1 >= 0));
//    assert((indexV2 < 4) && (indexV2 >= 0));

//    static int verticesface[4][4][4] = {
//        {//0
//         {-1, -1, -1, -1},
//         {-1, -1, 0, 1},
//         {-1, 0, -1, 2},
//         {-1, 1, 2, -1}},
//        {//1
//         {-1, -1, 0, 1},
//         {-1, -1, -1, -1},
//         {0, -1, -1, 3},
//         {1, -1, 3, -1}},
//        {//2
//         {-1, 0, -1, 2},
//         {0, -1, -1, 3},
//         {-1, -1, -1, -1},
//         {2, 3, -1, -1}},
//        {//3
//         {-1, 1, 2, -1},
//         {1, -1, 3, -1},
//         {2, 3, -1, -1},
//         {-1, -1, -1, -1}}};
//    return verticesface[indexV0][indexV1][indexV2];
//}

// static int FofEE(const int &indexE0, const int &indexE1)
//{
//     assert((indexE0 < 6) && (indexE0 >= 0));
//     assert((indexE1 < 6) && (indexE1 >= 0));
//     static int edgesface[6][6] = {{-1, 0, 1, 0, 1, -1},
//                                   {0, -1, 2, 0, -1, 2},
//                                   {1, 2, -1, -1, 1, 2},
//                                   {0, 0, -1, -1, 3, 3},
//                                   {1, -1, 1, 3, -1, 3},
//                                   {-1, 2, 2, 3, 3, -1}};

//    return edgesface[indexE0][indexE1];
//}

// static int FoppositeV (const int & indexV)
//{
//     assert(indexV < 4 && indexV >= 0);
//     static int oppFaces[4] = { 3, 2, 1, 0 };

//    return oppFaces[indexV];
//}

// static int VoppositeF (const int & indexF)
//{
//     assert(indexF < 4 && indexF >= 0);
//     static int oppVerts[4] = { 3, 2, 1, 0 };

//    return oppVerts[indexF];
//}

static int EoppositeE(const int &indexE) {
  assert(indexE < 6 && indexE >= 0);
  return 5 - indexE;
}

// void Get
// template <class ScalarType>
// bool GetBoundingBox(const std::vector<std::vector<ScalarType> > &VertPos,
//                     std::vector<ScalarType> &MinP,
//                     std::vector<ScalarType> &MaxP)
//{

//}

// template <class ScalarType>
// bool GetBoundingBox(const std::vector<std::vector<ScalarType> > &VertPos,
//                     std::vector<ScalarType> &MinP,
//                     std::vector<ScalarType> &MaxP)
//{

//}

static int WhichEdge(const std::vector<int> &TetraV, const int &indexV0,
                     const int &indexV1) {
  assert(TetraV.size() == 4);

  std::pair<int, int> EdgeTest(std::min(indexV0, indexV1),
                               std::max(indexV0, indexV1));
  for (size_t e = 0; e < 6; e++) {
    int IndexEV0 = TetraV[VofE(e, 0)];
    int IndexEV1 = TetraV[VofE(e, 1)];
    std::pair<int, int> EdgeV(std::min(IndexEV0, IndexEV1),
                              std::max(IndexEV0, IndexEV1));

    if (EdgeV == EdgeTest)
      return e;
  }
  return -1;
}

void FaceVert(const std::vector<std::vector<int>> &Connectivity,
              const size_t &IndexT, const size_t &IndexF,
              std::vector<int> &FaceV) {
  int LocV0 = VofF(IndexF, 0);
  int LocV1 = VofF(IndexF, 1);
  int LocV2 = VofF(IndexF, 2);

  int GlobV0 = Connectivity[IndexT][LocV0];
  assert(GlobV0 >= 0);
  int GlobV1 = Connectivity[IndexT][LocV1];
  assert(GlobV1 >= 0);
  int GlobV2 = Connectivity[IndexT][LocV2];
  assert(GlobV2 >= 0);

  FaceV.push_back(GlobV0);
  FaceV.push_back(GlobV1);
  FaceV.push_back(GlobV2);
}

void GetPPAdjacency(const std::vector<std::vector<int>> &Connectivity,
                    std::vector<std::vector<int>> &NextP,
                    std::vector<std::vector<int>> &NextF) {
  NextP.clear();
  NextP.resize(Connectivity.size(), std::vector<int>(4, -1));

  NextF = NextP;

  std::map<std::vector<int>, std::vector<std::pair<int, int>>> TetsFace;
  for (size_t i = 0; i < Connectivity.size(); i++) {
    assert(Connectivity[i].size() == 4);
    for (size_t f = 0; f < 4; f++) {
      std::vector<int> Key;
      FaceVert(Connectivity, i, f, Key);

      std::sort(Key.begin(), Key.end());

      // std::cout<<Key[0]<<","<<Key[1]<<","<<Key[2]<<","<<std::endl;
      TetsFace[Key].push_back(std::pair<int, int>(i, f));
    }
  }

  std::map<std::vector<int>, std::vector<std::pair<int, int>>>::iterator
      IteTetsFace;
  for (IteTetsFace = TetsFace.begin(); IteTetsFace != TetsFace.end();
       IteTetsFace++) {
    size_t sizeT = (*IteTetsFace).second.size();
    assert((sizeT == 1) || (sizeT == 2));
    if (sizeT == 1)
      continue; // no adjacency
    std::pair<int, int> TF0 = (*IteTetsFace).second[0];
    std::pair<int, int> TF1 = (*IteTetsFace).second[1];

    NextP[TF0.first][TF0.second] = TF1.first;
    NextF[TF0.first][TF0.second] = TF1.second;

    NextP[TF1.first][TF1.second] = TF0.first;
    NextF[TF1.first][TF1.second] = TF0.second;

    //        Adjacency[TF0.first][TF0.second].first=TF1.first;
    //        Adjacency[TF0.first][TF0.second].second=TF1.second;
    //        Adjacency[TF1.first][TF1.second].first=TF0.first;
    //        Adjacency[TF1.first][TF1.second].second=TF0.second;
  }
}

void GetVPAdjacency(
    const size_t numV, const std::vector<std::vector<int>> &Connectivity,
    std::vector<std::vector<std::pair<int, int>>> &VPAdjacency) {
  VPAdjacency.clear();
  VPAdjacency.resize(numV);

  for (size_t i = 0; i < Connectivity.size(); i++)
    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      int currV = Connectivity[i][j];
      assert(currV >= 0);
      assert(currV < (int)VPAdjacency.size());
      VPAdjacency[currV].push_back(std::pair<int, int>(i, j));
    }
}

void GetPPKernelAlongVertices(const size_t numV,
                              const std::vector<std::vector<int>> &Connectivity,
                              std::vector<std::vector<int>> &Adjacency,
                              int diffuse_step = 0) {

  std::vector<std::vector<std::pair<int, int>>> VPAdjacency;
  Tetra::GetVPAdjacency(numV, Connectivity, VPAdjacency);

  Adjacency.clear();
  Adjacency.resize(Connectivity.size());

  for (size_t i = 0; i < VPAdjacency.size(); i++) {
    std::vector<int> Polys;
    for (size_t j = 0; j < VPAdjacency[i].size(); j++)
      Polys.push_back(VPAdjacency[i][j].first);

    for (size_t k = 0; k < Polys.size(); k++)
      for (size_t w = 0; w < Polys.size(); w++) {
        if (k == w)
          continue;
        int IndexP0 = Polys[k];
        int IndexP1 = Polys[w];
        assert(IndexP0 != IndexP1);
        assert(IndexP0 >= 0);
        assert(IndexP1 >= 0);
        assert(IndexP0 < (int)Adjacency.size());
        assert(IndexP1 < (int)Adjacency.size());
        Adjacency[IndexP0].push_back(IndexP1);
        Adjacency[IndexP1].push_back(IndexP0);
      }
  }

  // remove duplicates
  for (size_t i = 0; i < Adjacency.size(); i++) {
    std::sort(Adjacency[i].begin(), Adjacency[i].end());
    std::vector<int>::iterator last =
        std::unique(Adjacency[i].begin(), Adjacency[i].end());
    Adjacency[i].erase(last, Adjacency[i].end());
  }

  if (diffuse_step > 0)
    Geo::ExtendStar(Adjacency, diffuse_step);
}

void GetPPKernelAlongVertices(const std::vector<std::vector<int>> &Connectivity,
                              std::vector<std::vector<int>> &Adjacency,
                              int diffuse_step = 0) {
  int NumIdx = Geo::FindMaxVertIdx(Connectivity);
  NumIdx++;
  GetPPKernelAlongVertices(NumIdx, Connectivity, Adjacency, diffuse_step);
}

void FFExternalMesh(const std::vector<std::vector<int>> &Connectivity,
                    const std::vector<std::vector<int>> &NextP,
                    const std::vector<std::vector<int>> &NextF,
                    std::vector<std::vector<int>> &Faces,
                    std::vector<int> &PolyIndex, std::vector<int> &PolyFace) {
  Faces.clear();
  PolyIndex.clear();
  PolyFace.clear();

  for (size_t i = 0; i < NextP.size(); i++)
    for (size_t f = 0; f < 4; f++) {
      assert(NextP[i].size() == 4);
      if (NextP[i][f] != -1)
        continue;
      assert(NextF[i][f] == -1);
      Faces.push_back(std::vector<int>());
      FaceVert(Connectivity, i, f, Faces.back());
      PolyIndex.push_back(i);
      PolyFace.push_back(f);
    }
}

void FFExternalPolys(const std::vector<std::vector<int>> &NextP,
                     const std::vector<std::vector<int>> &NextF,
                     std::vector<int> &ExtPoly) {
  ExtPoly.clear();
  for (size_t i = 0; i < NextP.size(); i++)
    for (size_t f = 0; f < 4; f++) {
      assert(NextP[i].size() == 4);
      if (NextP[i][f] != -1)
        continue;
      assert(NextF[i][f] == -1);
      ExtPoly.push_back(i);
    }
}

template <class ScalarType>
void TetraFaceMesh(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Connectivity,
                   const std::vector<int> &TetraIdx,
                   const std::vector<int> &TetraFaces,
                   std::vector<Geo::Point3<ScalarType>> &SubMeshVertPos,
                   std::vector<std::vector<int>> &SubMeshConnectivity) {
  assert(TetraIdx.size() == TetraFaces.size());

  SubMeshVertPos.clear();
  SubMeshConnectivity.clear();

  std::map<int, int> VertToVertMap;
  for (size_t i = 0; i < TetraIdx.size(); i++) {
    int CurrTet = TetraIdx[i];
    assert(CurrTet >= 0);
    assert(CurrTet < Connectivity.size());
    assert(Connectivity[CurrTet].size() == 4);

    int CurrF = TetraFaces[i];
    assert(CurrF >= 0);
    assert(CurrF < 4);

    SubMeshConnectivity.push_back(std::vector<int>(3, -1));

    for (size_t j = 0; j < 3; j++) {
      int IdxVj = Tetra::VofF(CurrF, j);
      IdxVj = Connectivity[CurrTet][IdxVj];

      // then map to the new vertex idx
      int NewVertIdx = -1;
      if (VertToVertMap.count(IdxVj) == 0) {
        SubMeshVertPos.push_back(VertPos[IdxVj]);
        NewVertIdx = SubMeshVertPos.size() - 1;
        VertToVertMap[IdxVj] = NewVertIdx;
      } else {
        NewVertIdx = VertToVertMap[IdxVj];
      }

      SubMeshConnectivity.back()[j] = NewVertIdx;
    }
  }
}

template <class ScalarType>
void BorderVertices(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &ExtFace,
                    std::vector<bool> &BorderV) {
  BorderV.clear();
  BorderV.resize(VertPos.size(), false);
  for (size_t i = 0; i < ExtFace.size(); i++)
    for (size_t j = 0; j < ExtFace[i].size(); j++)
      BorderV[ExtFace[i][j]] = true;
}

template <class ScalarType>
void TetraInternalFaces(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &TetraConnectivity,
                        std::vector<std::vector<int>> &TetraInternalFaces,
                        std::vector<std::vector<int>> &FT,
                        std::vector<std::vector<int>> &FTIndex) {
  TetraInternalFaces.clear();
  FT.clear();
  FTIndex.clear();
  std::vector<std::vector<int>> NextT, NextF;
  GetPPAdjacency(TetraConnectivity, NextT, NextF);
  //std::cout<<"Size NextP 0 "<<NextT.size()<<std::endl;
  for (size_t i = 0; i < NextT.size(); i++)
  {
    assert(NextT[i].size()==4);
    for (size_t j = 0; j < NextT[i].size(); j++) {
      int T0 = i;
      int IndexF0 = j;
      int T1 = NextT[i][j];
      int IndexF1 = NextF[i][j];
      if ((T1 != -1) && (T1 < T0))
        continue;
      assert(T1!=T0);
      std::vector<int> FaceV;
      FaceVert(TetraConnectivity, T0, IndexF0, FaceV);
      TetraInternalFaces.push_back(FaceV);
      FT.resize(FT.size() + 1);
      FT.back().push_back(T0);
      FTIndex.resize(FTIndex.size() + 1);
      FTIndex.back().push_back(IndexF0);
      if (T1 != -1) {
        assert(IndexF1 != -1);
        FT.back().push_back(T1);
        FTIndex.back().push_back(IndexF1);
      }
    }
  }
}

}; // namespace Tetra

#endif
