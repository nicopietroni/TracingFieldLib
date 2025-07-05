#ifndef MESH_FEATURES
#define MESH_FEATURES

#include "point3.h"
#include "pos3.h"
#include "mesh_features_base.h"
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <normal.h>
#include <set>
#include <string>
#include <sys/_types/_size_t.h>
#include <utility>
#include <vector>
// #include <tuple>

namespace Geo {

template <class ScalarType>
bool IsFeature(const std::vector<Geo::Point3<ScalarType>> &VertPos,
               const std::vector<std::vector<int>> &Faces, const int &IndexF0,
               const int &IndexF1, const ScalarType feature_angle) {
  assert(IndexF0 >= 0);
  assert(IndexF1 >= 0);
  assert(IndexF0 < Faces.size());
  assert(IndexF1 < Faces.size());
  Geo::Point3<ScalarType> Norm0 = Geo::Normal(VertPos, Faces, IndexF0);
  Geo::Point3<ScalarType> Norm1 = Geo::Normal(VertPos, Faces, IndexF1);
  return (Geo::AngleDeg(Norm0, Norm1) >= feature_angle);
}

template <class ScalarType>
bool IsCorner(const std::vector<Geo::Point3<ScalarType>> &VertPos,
              const int IndexV0, const int IndexVCentral, const int IndexV1,
              const ScalarType feature_angle) {
  assert(IndexV0 >= 0);
  assert(IndexVCentral >= 0);
  assert(IndexV1 >= 0);
  assert(IndexV0 != IndexVCentral);
  assert(IndexV1 != IndexVCentral);
  assert(IndexV0 != IndexV1);
  Geo::Point3<ScalarType> Dir0 = VertPos[IndexV0] - VertPos[IndexVCentral];
  Geo::Point3<ScalarType> Dir1 = VertPos[IndexV1] - VertPos[IndexVCentral];
  Dir0.Normalize();
  Dir1.Normalize();
  ScalarType Angle = Geo::AngleDeg(Dir0, Dir1);
  if (Angle < (180 - feature_angle))
    return true;
  if (Angle > (180 + feature_angle))
    return true;
  return (false);
}

template <class ScalarType>
void CornersByCount(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::pair<int, int>> &VertPair,
                    std::vector<int> &CornersIdx) {
  CornersIdx.clear();
  std::vector<int> VertCount(VertPos.size(),0);
  for (size_t i = 0; i < VertPair.size(); i++) {
    int IndexV0 = VertPair[i].first;
    int IndexV1 = VertPair[i].second;
    assert(IndexV0 >= 0);
    assert(IndexV1 >= 0);
    assert(IndexV0 != IndexV1);
    VertCount[IndexV0]++;
    VertCount[IndexV1]++;
  }

  for (size_t i = 0; i < VertCount.size(); i++)
  {
    if ((VertCount[i]==2)||(VertCount[i]==0))continue;
    CornersIdx.push_back(i);
  }
}

template <class ScalarType>
void Corners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
             const std::vector<std::pair<int, int>> &VertPair,
             const ScalarType feature_angle, std::vector<int> &CornersIdx) {
  CornersIdx.clear();
  std::map<int, std::vector<int>> VertCount;
  for (size_t i = 0; i < VertPair.size(); i++) {
    int IndexV0 = VertPair[i].first;
    int IndexV1 = VertPair[i].second;
    assert(IndexV0 >= 0);
    assert(IndexV1 >= 0);
    assert(IndexV0 != IndexV1);
    VertCount[IndexV0].push_back(IndexV1);
    VertCount[IndexV1].push_back(IndexV0);
  }

  std::map<int, std::vector<int>>::iterator IteMap;
  for (IteMap = VertCount.begin(); IteMap != VertCount.end(); IteMap++) {
    int IndexV = (*IteMap).first;
    int CountNum = (*IteMap).second.size();
    if (CountNum != 2) {
      CornersIdx.push_back(IndexV);
      continue;
    }
    {
      int IndexV0 = (*IteMap).second[0];
      int IndexVCentr = (*IteMap).first;
      int IndexV1 = (*IteMap).second[1];
      bool IsC =
          IsCorner(VertPos, IndexV0, IndexVCentr, IndexV1, feature_angle);
      if (IsC)
        CornersIdx.push_back(IndexVCentr);
    }
  }
}

template <class ScalarType>
void FeatureEdges(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces,
                  const std::vector<std::vector<int>> &NextF,
                  const std::vector<int> &FaceIDX,
                  const ScalarType feature_angle,
                  std::vector<std::pair<int, int>> &VertPairFeatures,
                  bool add_border = true) {
  assert(NextF.size() == Faces.size());
  VertPairFeatures.clear();

  for (size_t i = 0; i < FaceIDX.size(); i++) {
    int IndexF0 = FaceIDX[i];
    for (size_t j = 0; j < Faces[IndexF0].size(); j++) {
      int IndexV0 = Faces[IndexF0][j];
      int IndexV1 = Faces[IndexF0][(j + 1) % Faces[i].size()];
      std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));

      int IndexF1 = NextF[IndexF0][j];
      if (IndexF1 == -1) {
        if (add_border)
          VertPairFeatures.push_back(KeyE);
        continue;
      }
      if (IndexF0 > IndexF1)
        continue;
      if (!IsFeature(VertPos, Faces, IndexF0, IndexF1, feature_angle))
        continue;
      VertPairFeatures.push_back(KeyE);
    }
  }
}

template <class ScalarType>
void FeatureEdges(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces,
                  const std::vector<std::vector<int>> &NextF,
                  const ScalarType &feature_angle,
                  std::vector<std::pair<int, int>> &VertPairFeatures,
                  bool add_border = true) {
  std::vector<int> FaceIDX;
  for (size_t i = 0; i < Faces.size(); i++)
    FaceIDX.push_back(i);

  FeatureEdges<ScalarType>(VertPos, Faces, NextF, FaceIDX, feature_angle,
                           VertPairFeatures, add_border);
}

template <class ScalarType>
void FeatureCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    const std::vector<std::vector<int>> &NextF,
                    const std::vector<int> &FaceIDX,
                    const ScalarType feature_angle,
                    std::vector<std::pair<int, int>> &VertPairFeatures,
                    std::vector<int> &CornersIdx, bool add_border = true) {
  FeatureEdges<ScalarType>(VertPos, Faces, NextF, FaceIDX, feature_angle,
                           VertPairFeatures, add_border);
  Corners<ScalarType>(VertPos, VertPairFeatures, feature_angle, CornersIdx);
}

template <class ScalarType>
void FeatureCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    const std::vector<std::vector<int>> &NextF,
                    const ScalarType &feature_angle,
                    std::vector<std::pair<int, int>> &VertPairFeatures,
                    std::vector<int> &CornersIdx, bool add_border = true) {
  std::vector<int> FaceIDX;
  for (size_t i = 0; i < Faces.size(); i++)
    FaceIDX.push_back(i);

  FeatureCorners(VertPos, Faces, NextF, FaceIDX, feature_angle,
                 VertPairFeatures, CornersIdx, add_border);
}

template <class ScalarType>
void FeatureCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    const ScalarType &feature_angle,
                    std::vector<std::pair<int, int>> &VertPairFeatures,
                    std::vector<int> &CornersIdx, bool add_border = true) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  FeatureCorners(VertPos, Faces, NextF, feature_angle, VertPairFeatures,
                 CornersIdx, add_border);
}

template <class ScalarType>
void FeatureCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    const ScalarType &feature_angle,
                    std::vector<std::vector<bool>> &IsFaceEFeature,
                    std::vector<bool> &IsVertCorner, bool add_border = true) {

  std::vector<std::pair<int, int>> VertPairFeatures;
  std::vector<int> CornersIdx;
  FeatureCorners(VertPos, Faces, feature_angle, VertPairFeatures, CornersIdx,
                 add_border);

  std::set<std::pair<int, int>> VertPairFeaturesMap;
  for (size_t i = 0; i < VertPairFeatures.size(); i++) {
    int IndexV0 = VertPairFeatures[i].first;
    int IndexV1 = VertPairFeatures[i].second;
    std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                            std::max(IndexV0, IndexV1));
    VertPairFeaturesMap.insert(Key);
  }

  IsFaceEFeature.clear();
  IsVertCorner.clear();

  IsFaceEFeature.resize(Faces.size());
  for (size_t i = 0; i < Faces.size(); i++) {
    int VertNum = Faces[i].size();
    IsFaceEFeature[i] = std::vector<bool>(VertNum, false);
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % VertNum];
      std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      if (VertPairFeaturesMap.count(Key) == 0)
        continue;
      IsFaceEFeature[i][j] = true;
    }
  }

  IsVertCorner.resize(VertPos.size(), false);
  for (size_t i = 0; i < CornersIdx.size(); i++)
    IsVertCorner[CornersIdx[i]] = true;
}

template <class ScalarType>
void GetVertPairFromBoolFeatures(
    const std::vector<std::vector<int>> &Faces,
    const std::vector<std::vector<bool>> &IsFaceEFeature,
    std::vector<std::pair<int, int>> &VertPairFeatures) {

  assert(IsFaceEFeature.size() == Faces.size());
  VertPairFeatures.clear();
  for (size_t i = 0; i < Faces.size(); i++) {
    assert(IsFaceEFeature[i].size() == Faces[i].size());
    for (size_t j = 0; j < Faces[i].size(); j++) {
      if (!IsFaceEFeature[i][j])
        continue;
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      std::pair<int, int> Entry(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
      VertPairFeatures.push_back(Entry);
    }
  }
  std::sort(VertPairFeatures.begin(), VertPairFeatures.end());
  std::vector<std::pair<int, int>>::iterator last =
      std::unique(VertPairFeatures.begin(), VertPairFeatures.end());
  VertPairFeatures.erase(last, VertPairFeatures.end());
}

template <class ScalarType>
void GetVertIndexFromBoolCorners(const std::vector<bool> &IsVertCorner,
                                 std::vector<int> &CornersIdx) {
  CornersIdx.clear();
  for (size_t i = 0; i < IsVertCorner.size(); i++) {
    if (!IsVertCorner[i])
      continue;
    CornersIdx.push_back(i);
  }
}

// // this function get triple of consecutive feature vertices if they form
// // an angle which is Anle + - tolerance, also consider original borders if
// // they are not included in the features
// template <class ScalarType>
// void GetAdjFeatureAngle(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
//                         const std::vector<std::vector<int>> &Faces,
//                         const std::vector<std::pair<int, int>> &FeatureE,
//                         const std::vector<ScalarType> Angles,
//                         ScalarType Tolerance, std::vector<std::array<int, 3>
//                         > &TripleVert) {
//   TripleVert.clear();
//   if (FeatureE.size() == 0)
//     return;

//   // create a mep to retrieve original indexes
//   std::map<Geo::Point3<ScalarType>, int> PosToVert;
//   for (size_t i = 0; i < Pos3D.size(); i++)
//     PosToVert[Pos3D[i]] = i;

//   // make a copy of the mesh
//   std::vector<Geo::Point3<ScalarType>> Pos3DSplit = Pos3D;
//   std::vector<std::vector<int>> FacesSplit = Faces;

//   // open around features
//   SplitMeshFromEdges(Pos3DSplit, FacesSplit, FeatureE);

//   // then get the sequence of borders
//   std::vector<std::vector<int>> BorderSeq;
//   getBorderSequncesVert(FacesSplit, BorderSeq);

//   // then check the triple of consecutive ones whose angle
//   for (int i = 0; i < (int)BorderSeq.size(); i++) {
//     // those are always loops
//     int SizeL = BorderSeq[i].size();
//     for (int j = 0; j < (int)BorderSeq[i].size(); j++) {

//       int IndexV0 = ((j + SizeL - 1) % SizeL);
//       int IndexV1 = j;
//       int IndexV2 = ((j + 1) % SizeL);

//       IndexV0 = BorderSeq[i][IndexV0];
//       IndexV1 = BorderSeq[i][IndexV1];
//       IndexV2 = BorderSeq[i][IndexV2];

//       Geo::Point3<ScalarType> P0 = Pos3DSplit[IndexV0];
//       Geo::Point3<ScalarType> P1 = Pos3DSplit[IndexV1];
//       Geo::Point3<ScalarType> P2 = Pos3DSplit[IndexV2];

//       // then compute the angle
//       Geo::Point3<ScalarType> Dir0 = P0 - P1;
//       Geo::Point3<ScalarType> Dir1 = P2 - P1;
//       Dir0.Normalize();
//       Dir1.Normalize();
//       ScalarType CurrAngle = AngleDeg(Dir0, Dir1);

//       // then in this case add the triple
//       bool to_select = false;
//       for (size_t k =0;k < Angles.size(); k++)
//         to_select |= (fabs(Angles[k] - CurrAngle) < Tolerance);

//       if (!to_select)
//         continue;

//       assert(PosToVert.count(P0) > 0);
//       assert(PosToVert.count(P1) > 0);
//       assert(PosToVert.count(P2) > 0);
//       int IndexV0Or = PosToVert[P0];
//       int IndexV1Or = PosToVert[P1];
//       int IndexV2Or = PosToVert[P2];

//       // same order to avoid duplicates
//       if (IndexV0Or > IndexV2Or)
//         std::swap(IndexV0Or, IndexV2Or);

//       TripleVert.push_back({IndexV0Or, IndexV1Or,IndexV2Or});
//     }
//   }
//   // then sort and unique
//   std::sort(TripleVert.begin(), TripleVert.end());
//   std::vector<std::array<int, 3> >::iterator last =
//   std::unique(TripleVert.begin(), TripleVert.end()); TripleVert.erase(last,
//   TripleVert.end());
// }

// template <class ScalarType>
// void FeatureCorners(const std::vector<Geo::Point3<ScalarType> > &VertPos,
//                   const std::vector<std::vector<int> > &Faces,
//                   const std::vector<std::vector<int> > &NextF,
//                   const std::vector<int> &FaceIDX,
//                   const ScalarType feature_angle,
//                   std::vector<std::pair<int,int> > &VertPairFeatures,
//                   std::vector<int> &CornersIdx,
//                   bool add_border=true)
// {
//   for ()
// }

template <class ScalarType>
ScalarType
SignedDihedralAngleRad(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::vector<int>> &NextF,
                       const std::vector<std::vector<int>> &NextE, 
                       int IndexF,
                       int IndexE) {
  int IndexF0 = IndexF;
  assert(IndexF0 >= 0);
  assert(IndexE >= 0);
  assert(IndexF0 < Faces.size());
  assert(IndexF0 < NextF.size());
  assert(IndexF0 < NextE.size());
  assert(IndexE < Faces[IndexF0].size());
  assert(IndexE < NextF[IndexF0].size());
  assert(IndexE < NextE[IndexF0].size());

  int IndexV0=Faces[IndexF][IndexE];
  int IndexV1=Faces[IndexF][(IndexE+1)%Faces[IndexF].size()];
  Geo::Point3<ScalarType> DirE=Pos3D[IndexV1]-Pos3D[IndexV0];
  DirE.Normalize();

  int IndexF1 = NextF[IndexF0][IndexE];

  // in this case it is border
  if (IndexF1 == -1)
    return 0;
  Geo::Point3<ScalarType> Norm0 = Normal(Pos3D, Faces, IndexF0);
  Geo::Point3<ScalarType> Norm1 = Normal(Pos3D, Faces, IndexF1);
  Geo::Point3<ScalarType> c = Geo::Cross(Norm0, Norm1);
  ScalarType d = Geo::Dot(c, DirE);
  ScalarType AngleRad = Geo::AngleRad(Norm0,Norm1);
  if (d>0)
    return AngleRad;
  else
    return -AngleRad;
}

template <class ScalarType>
void 
SignedDihedralAngleRad(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                       const std::vector<std::vector<int>> &Faces,
                       std::vector<std::vector<ScalarType>> &Dihandle) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  Dihandle.clear();
  Dihandle.resize(Faces.size());
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexF = i;
      int IndexE = j;
      ScalarType DAngle =
          SignedDihedralAngleRad(Pos3D, Faces, NextF, NextE, IndexF, IndexE);
      Dihandle[i].push_back(DAngle);
    }
}

template <class ScalarType>
void getFeatureTypes(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const std::vector<std::pair<int, int>> &Edges,
                     std::vector<FeatureType> &FType) {
  FType.clear();

  std::vector<std::vector<ScalarType>> Dihandle;
  SignedDihedralAngleRad(VertPos, Faces, Dihandle);

  std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeMap;
  FindEdgeFacesCorrespondences(VertPos, Faces, EdgeMap);
  for (size_t i = 0; i < Edges.size(); i++) {
    std::pair<int, int> KeyE = Edges[i];
    assert(EdgeMap.count(KeyE) > 0);
    assert(EdgeMap[KeyE].size() > 0);

    int IndexF = EdgeMap[KeyE][0].first;
    int IndexE = EdgeMap[KeyE][0].second;

    assert(IndexF >= 0);
    assert(IndexE >= 0);
    assert(IndexF < Dihandle.size());
    assert(IndexE < Dihandle[IndexF].size());
    ScalarType CurrAngle = Dihandle[IndexF][IndexE];
    if (CurrAngle < 0)
      FType.push_back(FTConcave);
    else
      FType.push_back(FTConvex);
  }
}

} // namespace Geo

#endif
