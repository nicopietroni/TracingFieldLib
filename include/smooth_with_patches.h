#include "IO/io_obj.h"
#include "connectivity.h"
#include "geodesic.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "point3.h"
#include "pos3.h"
#include "quad3.h"
#include <cassert>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#ifndef SMOOTH_WITH_PATCHES
#define SMOOTH_WITH_PATCHES

namespace Geo {

template <class ScalarType> class QuadPatchFunctions {

public:
  static bool
  HasConcaveAngle(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                  const std::vector<std::vector<int>> &PatchFaces) {
    std::vector<int> BorderV;
    getBorderVert(PatchFaces, BorderV);
    std::vector<std::vector<int>> PatchVF, PatchVFIndex;
    getVFAdjacency(PatchFaces, PatchVF, PatchVFIndex);
    for (size_t i = 0; i < BorderV.size(); i++) {
      int IndexV = BorderV[i];
      assert(IndexV >= 0);
      assert(IndexV < PatchVF.size());
      if (PatchVF[IndexV].size() > 2)
        return true;
    }
    return false;
  }

  static void FindCorners(const std::vector<std::vector<int>> &PatchFaces,
                          std::vector<int> &CornersIdx) {
    std::vector<int> BorderIdx;
    getBorderVert(PatchFaces, BorderIdx);

    std::vector<std::vector<int>> VF, VFIndex;
    getVFAdjacency(PatchFaces, VF, VFIndex);

    CornersIdx.clear();
    for (size_t i = 0; i < BorderIdx.size(); i++) {
      int TestV = BorderIdx[i];
      if (VF[TestV].size() == 1)
        CornersIdx.push_back(TestV);
    }
  }

  static void FindCorners(const std::vector<Pos3> &Sides,
                          std::vector<int> &CornersIdx) {
    CornersIdx.clear();
    for (size_t i = 0; i < Sides.size(); i++) {
      int CurrI = i;
      int NextI = (i + 1) % Sides.size();
      if (Sides[CurrI].F != Sides[NextI].F)
        continue;
      CornersIdx.push_back(CurrI);
    }
  }

  static int NumCorners(const std::vector<std::vector<int>> &PatchFaces) {
    std::vector<int> CornersIdx;
    FindCorners(PatchFaces, CornersIdx);
    return CornersIdx.size();
  }

  static bool
  IsDiskPatch(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
              const std::vector<std::vector<int>> &PatchFaces) {

    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(PatchFaces, NextF, NextE);

    std::vector<std::vector<Pos3>> BorderPos;
    getBorderSequnces(PatchFaces, NextF, NextE, BorderPos);
    if (BorderPos.size() != 1)
      return false;

    int Genus = PolyMeshGenus(PatchVertPos, PatchFaces);
    if (Genus != 0)
      return false;

    return true;
  }

  static bool RemoveDoublets(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                             std::vector<std::vector<int>> &QuadFaces,
                             std::vector<int> &To_Faces) {
    bool SolvedDoublet = false;

    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(QuadFaces, NextF, NextE);

    std::set<int> ModifiedV;
    std::vector<int> ToDeleteF;
    std::vector<int> ToDeleteV;

    std::map<int, std::vector<int>> MappedF;

    int InitialFsize = QuadFaces.size();

    int num_doublets = 0;
    for (size_t i = 0; i < NextF.size(); i++) {

      for (size_t j = 0; j < NextF[i].size(); j++) {
        int sizeF = NextF[i].size();
        assert(sizeF == 4);
        int IndexPrev = (j + sizeF - 1) % sizeF;
        if (NextF[i][j] != NextF[i][IndexPrev])
          continue;

        // this is the one to be deleted
        int CurrV = QuadFaces[i][j];
        int IndexV[4] = {-1, -1, -1, -1};
        IndexV[0] = QuadFaces[i][(j + 1) % sizeF];
        IndexV[1] = QuadFaces[i][(j + 2) % sizeF];
        IndexV[2] = QuadFaces[i][(j + 3) % sizeF];
        assert(CurrV != IndexV[0]);
        assert(CurrV != IndexV[1]);
        assert(CurrV != IndexV[2]);
        int OppF = NextF[i][j];
        int OppE = NextE[i][j];
        IndexV[3] = QuadFaces[OppF][(OppE + 3) % sizeF];
        assert(CurrV != IndexV[3]);
        assert(IndexV[3] != IndexV[0]);
        assert(IndexV[3] != IndexV[1]);
        assert(IndexV[3] != IndexV[2]);

        bool HasModified = false;
        for (size_t i = 0; i < 4; i++)
          HasModified |= (ModifiedV.count(IndexV[i]) > 0);

        if (HasModified)
          continue;

        for (size_t i = 0; i < 4; i++)
          ModifiedV.insert(IndexV[i]);

        int IndexF0 = i;
        int IndexF1 = OppF;

        ToDeleteV.push_back(CurrV);
        ToDeleteF.push_back(IndexF0);
        ToDeleteF.push_back(IndexF1);

        // add new face
        // QuadFaces.push_back(VertFan);
        QuadFaces.resize(QuadFaces.size() + 1);
        for (size_t i = 0; i < 4; i++)
          QuadFaces.back().push_back(IndexV[i]);

        int IndexNewF = QuadFaces.size() - 1;
        MappedF[IndexNewF].push_back(IndexF0);
        MappedF[IndexNewF].push_back(IndexF1);
        SolvedDoublet = true;
        num_doublets++;
      }
    }
    std::cout << "Solved " << num_doublets << " Doublets" << std::endl;

    if (SolvedDoublet) {
      To_Faces = std::vector<int>(InitialFsize, -1);

      std::vector<int> OriginalFace, OriginalVert;
      RemoveDeleted(QuadVertPos, QuadFaces, ToDeleteV, ToDeleteF, OriginalVert,
                    OriginalFace);
      // then remap
      assert(OriginalFace.size() == QuadFaces.size());
      for (size_t i = 0; i < OriginalFace.size(); i++) {
        int IndexFOr = OriginalFace[i];
        // if (IndexFOr>=To_Faces.size())continue;//new added face

        if (MappedF.count(IndexFOr) == 0) {
          assert(IndexFOr >= 0);
          assert(IndexFOr < To_Faces.size());
          assert(To_Faces[IndexFOr] == -1);
          To_Faces[IndexFOr] = i;
        } else {
          std::vector<int> OrFacesIdx = MappedF[IndexFOr];
          assert(OrFacesIdx.size() == 2);
          int IndexFOr0 = OrFacesIdx[0];
          int IndexFOr1 = OrFacesIdx[1];
          assert(IndexFOr0 >= 0);
          assert(IndexFOr0 < To_Faces.size());
          assert(IndexFOr1 >= 0);
          assert(IndexFOr1 < To_Faces.size());
          assert(To_Faces[IndexFOr0] == -1);
          assert(To_Faces[IndexFOr1] == -1);
          To_Faces[IndexFOr0] = i;
          To_Faces[IndexFOr1] = i;
        }
      }
      // sanity check
      int WrongAss = 0;
      for (size_t i = 0; i < To_Faces.size(); i++) {
        if (To_Faces[i] < 0)
          WrongAss++;
      }
      assert(WrongAss == 0);
    }
    return SolvedDoublet;
  }

  static bool RemoveDoublets(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                             std::vector<std::vector<int>> &QuadFaces,
                             std::vector<std::vector<int>> &FromFace) {
    std::vector<int> To_Faces;
    int NumF0 = QuadFaces.size();
    bool had_doublet = RemoveDoublets(QuadVertPos, QuadFaces, To_Faces);
    FromFace.clear();
    if (!had_doublet)
      return false;

    assert(To_Faces.size() == NumF0);

    // get the from Faces
    FromFace.resize(QuadFaces.size());
    for (size_t i = 0; i < To_Faces.size(); i++) {
      int IndexNewF = To_Faces[i];
      assert(IndexNewF >= 0);
      assert(IndexNewF < FromFace.size());
      FromFace[IndexNewF].push_back(i);
    }

    // consistency check
    for (size_t i = 0; i < FromFace.size(); i++) {
      assert(FromFace[i].size() > 0);
    }

    return true;
  }

  static bool RemoveDoublets(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                             std::vector<std::vector<int>> &QuadFaces,
                             std::vector<std::pair<int, int>> &Features,
                             std::vector<int> &Corners,
                             std::vector<int> &PatchIdx) {
    assert(PatchIdx.size() == QuadFaces.size());

    // get the from Faces
    std::vector<std::vector<int>> FromFace;
    std::vector<Geo::Point3<ScalarType>> OldPos = QuadVertPos;
    bool had_doublets = RemoveDoublets(QuadVertPos, QuadFaces, FromFace);
    if (!had_doublets)
      return false;

    // map patches (get the first one , shoul be ok)
    std::vector<int> NewPatchIDX;
    for (size_t i = 0; i < QuadFaces.size(); i++) {
      int OrFace = FromFace[i][0];
      assert(OrFace >= 0);
      assert(OrFace < PatchIdx.size());
      NewPatchIDX.push_back(PatchIdx[OrFace]);
    }

    // map corners and Features
    std::vector<int> VertMap;
    int NumMap = FindVertCorrespondences(OldPos, QuadVertPos, VertMap);

    std::vector<std::pair<int, int>> NewFeatures;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 < VertMap.size());
      assert(IndexV1 < VertMap.size());
      int NewIndexV0 = VertMap[IndexV0];
      int NewIndexV1 = VertMap[IndexV1];
      if ((NewIndexV0 == -1) && (NewIndexV0 == -1))
        continue;

      std::pair<int, int> NewFeat(std::min(NewIndexV0, NewIndexV1),
                                  std::max(NewIndexV0, NewIndexV1));
      NewFeatures.push_back(NewFeat);
    }

    std::vector<int> NewCorners;
    for (size_t i = 0; i < Corners.size(); i++) {
      int IndexV = Corners[i];
      assert(IndexV >= 0);
      assert(IndexV < VertMap.size());
      int NewIndexV = VertMap[IndexV];
      if (NewIndexV == -1)
        continue;

      NewCorners.push_back(NewIndexV);
    }

    PatchIdx = NewPatchIDX;
    Corners = NewCorners;
    Features = NewFeatures;
  }
};

template <class ScalarType> class ConvexPatchFinder {

  const std::vector<Geo::Point3<ScalarType>> &QuadVertPos;
  const std::vector<std::vector<int>> &QuadFaces;
  // const std::vector<std::pair<int, int>> &QuadFeatures;

  std::vector<std::vector<int>> VF, VFIndex;
  std::vector<std::vector<int>> FF, FFEdge;
  std::vector<bool> IsB;
  // std::set<std::pair<int, int>> QuadFeaturesSet;

  // using clock_t = std::chrono::steady_clock;
  // using ns = std::chrono::nanoseconds;

  int ExpandingSideScore(const std::vector<Geo::Point3<ScalarType>> &PatchPos,
                         const std::vector<std::vector<int>> &PatchFaces,
                         const std::vector<std::vector<Pos3>> &Side,
                         const int CentralV,
                         std::vector<std::pair<ScalarType, int>> &SideDist) {

    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(PatchFaces, PatchPos.size(), VVAdj);
    Geo::TMark TempMark;
    typename Geodesic<ScalarType>::GeodesicParam GParam;
    GParam.Seeds.push_back(CentralV);
    GParam.unitary_dist = true;
    Geodesic<ScalarType>::NodeDijkstra(PatchPos, VVAdj, TempMark, GParam);

    SideDist.clear();
    for (int i = 0; i < Side.size(); i++) {
      ScalarType MinSideD = std::numeric_limits<ScalarType>::max();
      for (int j = 0; j < Side[i].size(); j++) {
        int IndexV0 = Side[i][j].V;
        int IndexV1 = Side[i][j].VFlip(PatchFaces);
        assert(IndexV0 >= 0);
        assert(IndexV1 >= 0);
        assert(IndexV0 < GParam.Dist.size());
        assert(IndexV1 < GParam.Dist.size());
        MinSideD = std::min(MinSideD, GParam.Dist[IndexV0]);
        MinSideD = std::min(MinSideD, GParam.Dist[IndexV1]);
      }
      assert(MinSideD != std::numeric_limits<ScalarType>::max());
      SideDist.push_back(std::pair<ScalarType, int>(MinSideD, i));
    }
    std::sort(SideDist.begin(), SideDist.end());
  }

  void Expand(std::vector<int> &CurrFaces, const std::vector<Pos3> &ExpandPos) {

    int ExpF0 = ExpandPos[0].FFlip(FF, FFEdge);
    int ExpF1 = ExpandPos.back().FFlip(FF, FFEdge);
    CurrFaces.push_back(ExpF0);
    CurrFaces.push_back(ExpF1);

    for (size_t i = 1; i < ExpandPos.size(); i++) {
      int IndexV = ExpandPos[i].VFlip(QuadFaces);
      CurrFaces.insert(CurrFaces.end(), VF[IndexV].begin(), VF[IndexV].end());
    }
    std::sort(CurrFaces.begin(), CurrFaces.end());
    CurrFaces.erase(std::unique(CurrFaces.begin(), CurrFaces.end()),
                    CurrFaces.end());
  }

  bool CanExpand(const std::vector<int> &CurrFaces,
                 const std::vector<Pos3> &ExpandPos) {

    // check boundary edges
    for (size_t i = 0; i < ExpandPos.size(); i++) {
      if (ExpandPos[i].IsBorder(FF))
        return false;

      // int IndexV0 = ExpandPos[i].V;
      // int IndexV1 = ExpandPos[i].VFlip(QuadFaces);
      // assert(IndexV0 >= 0);
      // assert(IndexV1 >= 0);
      // assert(IndexV0 < QuadVertPos.size());
      // assert(IndexV1 < QuadVertPos.size());
      // std::pair<int, int> Key(std::min(IndexV0, IndexV1),
      //                         std::max(IndexV0, IndexV1));
      // if (QuadFeaturesSet.count(Key) > 0)
      //   return false;
    }

    // check internal features (only vertex in between)
    for (size_t i = 1; i < ExpandPos.size() - 1; i++) {
      int IndexV0 = ExpandPos[i].V;
      int IndexV1 = ExpandPos[i].VFlip(QuadFaces);
      // if (CornersV.count(IndexV0) > 0)
      //   return false;
      // if (CornersV.count(IndexV1) > 0)
      //   return false;
      if (IsB[IndexV0])
        return false;
      if (IsB[IndexV1])
        return false;
    }

    // try to expand
    std::vector<int> NextFaces = CurrFaces;
    Expand(NextFaces, ExpandPos);

    // get the mesh
    std::vector<int> OriginalFaces, OriginalVert;
    std::vector<Geo::Point3<ScalarType>> PatchVertPos;
    std::vector<std::vector<int>> PatchFaces;
    CopySubMesh(QuadVertPos, QuadFaces, NextFaces, PatchVertPos, PatchFaces,
                OriginalVert, OriginalFaces);

    std::vector<std::vector<Pos3>> BorderPos;
    std::vector<std::vector<int>> PatchNextF, PatchNextE;
    getFFAdjacency(PatchFaces, PatchNextF, PatchNextE);

    getBorderSequnces(PatchFaces, PatchNextF, PatchNextE, BorderPos);
    if (BorderPos.size() != 1)
      return false;

    int Genus = PolyMeshGenus(PatchVertPos, PatchFaces);
    if (Genus != 0)
      return false;

    std::vector<std::vector<Pos3>> Sides;
    FindSides(PatchVertPos, PatchFaces, Sides);
    if (Sides.size() < 3)
      return false;

    return true;
  }

  void FindSides(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                 const std::vector<std::vector<int>> &QuadFaces,
                 std::vector<std::vector<Pos3>> &Sides) {
    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(QuadFaces, NextF, NextE);

    std::vector<std::vector<Pos3>> BorderPos;
    getBorderSequnces(QuadFaces, NextF, NextE, BorderPos);
    // if (BorderPos.size() != 1) {
    //   WriteOBJ("testMesh.obj", QuadVertPos, QuadFaces);
    //   assert(0);
    // }
    assert(BorderPos.size() == 1);

    // then find corners
    std::vector<int> CornersIdx;
    std::vector<Pos3> Border = BorderPos[0];
    QuadPatchFunctions<ScalarType>::FindCorners(Border, CornersIdx);
    // assert(CornersIdx.size() >= 2);
    // shift all by one
    for (size_t i = 0; i < CornersIdx.size(); i++)
      CornersIdx[i] = (CornersIdx[i] + 1) % Border.size();
    // first rotate to have first corner first
    std::set<int> CornerSet(CornersIdx.begin(), CornersIdx.end());
    int StartIdx = CornersIdx[0];
    Sides.clear();
    Sides.resize(1);
    Sides.back().push_back(Border[StartIdx]);
    for (int i = 1; i < Border.size(); i++) {
      int CurrIdx = (StartIdx + i) % Border.size();

      if (CornerSet.count(CurrIdx) > 0)
        Sides.resize(Sides.size() + 1);

      Sides.back().push_back(Border[CurrIdx]);
    }
  }

  int FindCentralVert(const std::vector<int> &SeedFaceIdx) {

    std::vector<Geo::Point3<ScalarType>> SeedVertPos;
    std::vector<std::vector<int>> SeedFaces;
    std::vector<int> OriginalFaces, OriginalVert;

    CopySubMesh(QuadVertPos, QuadFaces, SeedFaceIdx, SeedVertPos, SeedFaces,
                OriginalVert, OriginalFaces);

    std::vector<int> CornersIdx;
    QuadPatchFunctions<ScalarType>::FindCorners(SeedFaces, CornersIdx);

    typename Geodesic<ScalarType>::GeodesicParam GParam;
    GParam.Seeds = CornersIdx;
    GParam.unitary_dist = true;

    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(SeedFaces, SeedVertPos.size(), VVAdj);
    Geo::TMark TempMark;

    Geodesic<ScalarType>::NodeDijkstra(SeedVertPos, VVAdj, TempMark, GParam);

    ScalarType MaxD = -1;
    int IndexSeed = -1;
    for (int i = 0; i < SeedVertPos.size(); i++) {
      if (GParam.Dist[i] < MaxD)
        continue;
      IndexSeed = i;
      MaxD = GParam.Dist[i];
    }
    assert(IndexSeed >= 0);
    return OriginalVert[IndexSeed];
  }

public:
  void FindConvexPatchFrom(const std::vector<int> &SeedF,
                           std::vector<int> &FaceIdx,
                           const int MaxExpandSteps = -1,
                           bool debug_msg = false) {

    assert(SeedF.size() >= 0);

    // initialize with faces sharing vertex
    std::vector<int> CurrFIdx = SeedF;
    int SeedIndexV = FindCentralVert(SeedF);
    assert(SeedIndexV >= 0);
    assert(SeedIndexV < QuadVertPos.size());
    bool has_expanded = false;
    int curr_step = 0;
    do {

      if (debug_msg)
        std::cout << "Expanding Step " << curr_step << std::endl;
      // copy the mesh from current set
      std::vector<Geo::Point3<ScalarType>> PatchVertPos;
      std::vector<std::vector<int>> PatchFaces;
      std::vector<int> OriginalFaces, OriginalVert;

      CopySubMesh(QuadVertPos, QuadFaces, CurrFIdx, PatchVertPos, PatchFaces,
                  OriginalVert, OriginalFaces);

      if (!QuadPatchFunctions<ScalarType>::HasConcaveAngle(PatchVertPos,
                                                           PatchFaces)) {
        FaceIdx = CurrFIdx;
        // BestPatchOriginalVert = OriginalVert;
      }

      // retrieve the central vertex
      int CentralV = -1;
      for (size_t i = 0; i < OriginalVert.size(); i++) {
        if (OriginalVert[i] == SeedIndexV) {
          CentralV = i;
          break;
        }
      }
      assert(CentralV != -1);

      // then update sides
      std::vector<std::vector<Pos3>> Sides;
      FindSides(PatchVertPos, PatchFaces, Sides);

      // map sides to original mesh
      std::vector<std::vector<Pos3>> OriginalSides = Sides;
      for (size_t i = 0; i < OriginalSides.size(); i++)
        for (size_t j = 0; j < OriginalSides[i].size(); j++) {
          int CurrF = OriginalSides[i][j].F;
          int CurrV = OriginalSides[i][j].V;
          OriginalSides[i][j].F = OriginalFaces[CurrF];
          OriginalSides[i][j].V = OriginalVert[CurrV];
        }

      // get the score of expansion
      std::vector<std::pair<ScalarType, int>> SideDist;
      ExpandingSideScore(PatchVertPos, PatchFaces, Sides, CentralV, SideDist);

      // then expand the one that can
      has_expanded = false;
      for (size_t i = 0; i < SideDist.size(); i++) {
        int IndexSide = SideDist[i].second;
        int CurrDist = SideDist[i].first;
        if (!CanExpand(CurrFIdx, OriginalSides[IndexSide]))
          continue;

        // already expanded to the limit
        if ((MaxExpandSteps != -1) && (CurrDist >= MaxExpandSteps))
          return;

        Expand(CurrFIdx, OriginalSides[IndexSide]);
        has_expanded = true;
        curr_step++;
        break;
      }
    } while (has_expanded);
  }

  void Init() {
    Geo::getVFAdjacency(QuadFaces, VF, VFIndex);
    Geo::getFFAdjacency(QuadFaces, FF, FFEdge);
    Geo::getBorderVert(QuadFaces, IsB);

    // then initialize the domains
    assert(IsAllQuad(QuadFaces));
  }

  ConvexPatchFinder(const std::vector<Geo::Point3<ScalarType>> &_QuadVertPos,
                    const std::vector<std::vector<int>> &_QuadFaces)
      : QuadVertPos(_QuadVertPos), QuadFaces(_QuadFaces) {}
};

#define MIN_KERN 1
#define MAX_KERN 4

template <class ScalarType> class QuadPatchSmooth {

  std::vector<Geo::Point3<ScalarType>> &QuadVertPos;
  const std::vector<std::vector<int>> &QuadFaces;
  std::vector<int> QuadPatch;
  std::vector<std::vector<int>> PatchQuad;
  const std::vector<std::pair<int, int>> &QuadFeatures;
  const std::vector<int> &Constraints;

  std::vector<std::vector<int>> FF, FFEdge;
  std::vector<std::vector<int>> SmoothDomainFaces;

  // ConvexPatchFinder<ScalarType> CPatch;
  std::vector<std::vector<int>> VF, VFIndex;
  //    std::vector<std::vector<int>> FF, FFEdge;
  std::vector<bool> IsB;

  using clock_t = std::chrono::steady_clock;
  using ns = std::chrono::nanoseconds;

  void getBiggestComponent(std::vector<int> &FacesIdx) {

    std::vector<Geo::Point3<ScalarType>> PatchVertPos;
    std::vector<std::vector<int>> PatchFaces;
    std::vector<int> OriginalFaces, OriginalVert;

    CopySubMesh(QuadVertPos, QuadFaces, FacesIdx, PatchVertPos, PatchFaces,
                OriginalVert, OriginalFaces);

    std::vector<std::vector<int>> ComponentsF;
    FFConnectedComponenent(PatchFaces, ComponentsF);
    if (ComponentsF.size() == 1)
      return;

    // get biggest component
    int bestC = 0;
    int bestSize = ComponentsF[0].size();
    for (size_t i = 1; i < ComponentsF.size(); i++) {
      if (ComponentsF[i].size() <= bestSize)
        continue;
      bestSize = ComponentsF[i].size();
      bestC = i;
    }
    FacesIdx = ComponentsF[bestC];

    // then map back to original one
    for (size_t i = 0; i < FacesIdx.size(); i++)
      FacesIdx[i] = OriginalFaces[FacesIdx[i]];
  }

  void GetDomainMesh(const int &IndexP,
                     std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                     std::vector<std::vector<int>> &PatchFaces,
                     std::vector<int> &OriginalFaces,
                     std::vector<int> &OriginalVert) const {
    assert(IndexP >= 0);
    assert(IndexP < SmoothDomainFaces.size());
    CopySubMesh(QuadVertPos, QuadFaces, SmoothDomainFaces[IndexP], PatchVertPos,
                PatchFaces, OriginalVert, OriginalFaces);
  }

  enum PatchType { PTOk, PTConcave, PTNonDiskLike, PTNonCorners };

  PatchType
  GetPatchType(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
               const std::vector<std::vector<int>> &PatchFaces) {
    if (QuadPatchFunctions<ScalarType>::HasConcaveAngle(PatchVertPos,
                                                        PatchFaces))
      return PTConcave;

    if (!QuadPatchFunctions<ScalarType>::IsDiskPatch(PatchVertPos, PatchFaces))
      return PTNonDiskLike;

    if (QuadPatchFunctions<ScalarType>::NumCorners(PatchFaces) < 3)
      return PTNonCorners;

    return PTOk;
  }

  void CheckPatches() {
    int ConcaveP = 0;
    int NonDiskP = 0;
    int NonRightCorn = 0;
    for (size_t i = 0; i < SmoothDomainFaces.size(); i++) {
      if (SmoothDomainFaces[i].size() == 0)
        continue;

      std::vector<Geo::Point3<ScalarType>> PatchVertPos;
      std::vector<std::vector<int>> PatchFaces;
      std::vector<int> OriginalFaces, OriginalVert;

      GetDomainMesh(i, PatchVertPos, PatchFaces, OriginalFaces, OriginalVert);

      PatchType PT = GetPatchType(PatchVertPos, PatchFaces);
      switch (PT) {
      case PTConcave:
        ConcaveP++;
        break;
      case PTNonDiskLike:
        NonDiskP++;
        break;
      case PTNonCorners:
        NonRightCorn++;
        break;
      default:
        assert(PT == PTOk);
      }
    }
    std::cout << "Total Patches Num:" << SmoothDomainFaces.size() << std::endl;
    std::cout << "Patch With Concave Vert:" << ConcaveP << std::endl;
    std::cout << "Patch Non Disk Like:" << NonDiskP << std::endl;
    std::cout << "Patch Right Num Corners:" << NonRightCorn << std::endl;
  }

  void AddBasicPatches() {
    for (size_t i = 0; i < PatchQuad.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> PatchVertPos;
      std::vector<std::vector<int>> PatchFaces;
      std::vector<int> OriginalFaces, OriginalVert;

      CopySubMesh(QuadVertPos, QuadFaces, PatchQuad[i], PatchVertPos,
                  PatchFaces, OriginalVert, OriginalFaces);

      PatchType PT = GetPatchType(PatchVertPos, PatchFaces);

      if (PT != PTOk)
        continue;

      SmoothDomainFaces.push_back(PatchQuad[i]);
    }
  }

  // void AddVertexPatches()
  // {
  //   for (size_t i=0;i<)
  // }

  void SolvePatchIfNeeded(std::vector<int> &FaceIdx,
                          const std::vector<int> &StartingF) {
    std::vector<Geo::Point3<ScalarType>> PatchVertPos;
    std::vector<std::vector<int>> PatchFaces;
    std::vector<int> OriginalFaces, OriginalVert;

    CopySubMesh(QuadVertPos, QuadFaces, FaceIdx, PatchVertPos, PatchFaces,
                OriginalVert, OriginalFaces);
    // remap starting vert
    std::vector<int> StartingFPatch;
    for (size_t i = 0; i < StartingF.size(); i++) {
      int TargetF = StartingF[i];
      int NewIdxTargetF = -1;
      for (size_t j = 0; j < OriginalFaces.size(); j++) {
        if (OriginalFaces[j] != TargetF)
          continue;
        NewIdxTargetF = j;
        break;
      }
      assert(NewIdxTargetF != -1);
      StartingFPatch.push_back(NewIdxTargetF);
    }

    PatchType PT = GetPatchType(PatchVertPos, PatchFaces);
    if (PT == PTOk)
      return;

    ConvexPatchFinder<ScalarType> CPatch(PatchVertPos, PatchFaces);
    CPatch.Init();
    std::vector<int> ConvexFaceIdx;
    CPatch.FindConvexPatchFrom(StartingFPatch, ConvexFaceIdx);

    /// then maps back
    FaceIdx.clear();
    for (size_t i = 0; i < ConvexFaceIdx.size(); i++)
      FaceIdx.push_back(OriginalFaces[ConvexFaceIdx[i]]);
  }

  void AddEdgePatches() {
    // first find patch patch adjacency
    std::vector<std::pair<int, int>> PatchAdj;
    std::map<std::pair<int, int>, std::vector<int>> BorderF;

    for (size_t i = 0; i < FF.size(); i++)
      for (size_t j = 0; j < FF[i].size(); j++) {
        int IndexF0 = i;
        int IndexF1 = FF[i][j];
        if (IndexF1 == -1)
          continue;

        int IndexP0 = QuadPatch[IndexF0];
        int IndexP1 = QuadPatch[IndexF1];
        if (IndexP0 == IndexP1)
          continue;

        std::pair<int, int> KeyAdj(std::min(IndexP0, IndexP1),
                                   std::max(IndexP0, IndexP1));
        PatchAdj.push_back(KeyAdj);

        // int IndexV0 = QuadFaces[i][j];
        // int IndexV1 = QuadFaces[i][(j + 1) % QuadFaces[i].size()];
        BorderF[KeyAdj].push_back(IndexF0);
        BorderF[KeyAdj].push_back(IndexF1);
      }
    // unique all adjacence
    std::sort(PatchAdj.begin(), PatchAdj.end());
    PatchAdj.erase(std::unique(PatchAdj.begin(), PatchAdj.end()),
                   PatchAdj.end());

    // unique all in between vertices
    std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
    for (IteMap = BorderF.begin(); IteMap != BorderF.end(); IteMap++) {
      std::vector<int> FaceIdx = (*IteMap).second;
      std::sort(FaceIdx.begin(), FaceIdx.end());
      FaceIdx.erase(std::unique(FaceIdx.begin(), FaceIdx.end()), FaceIdx.end());
      (*IteMap).second = FaceIdx;
    }

    // add all adjacent patches
    int Size0 = SmoothDomainFaces.size();
    for (size_t i = 0; i < PatchAdj.size(); i++) {
      int IndexP0 = PatchAdj[i].first;
      int IndexP1 = PatchAdj[i].second;
      assert(IndexP0 >= 0);
      assert(IndexP1 >= 0);
      assert(IndexP0 < PatchQuad.size());
      assert(IndexP1 < PatchQuad.size());
      assert(IndexP0 != IndexP1);
      assert(IndexP1 >= 0);

      // add the two patches
      std::vector<int> PatchQ0 = PatchQuad[IndexP0];
      std::vector<int> PatchQ1 = PatchQuad[IndexP1];

      std::vector<int> EdgePatchQ = PatchQ0;
      EdgePatchQ.insert(EdgePatchQ.end(), PatchQ1.begin(), PatchQ1.end());

      std::pair<int, int> KeyAdj(std::min(IndexP0, IndexP1),
                                 std::max(IndexP0, IndexP1));
      assert(BorderF.count(KeyAdj) > 0);
      // border faces should be a single connected component
      std::vector<int> Single_Comp_BorderF = BorderF[KeyAdj];
      getBiggestComponent(Single_Comp_BorderF);
      SolvePatchIfNeeded(EdgePatchQ, Single_Comp_BorderF);
      SmoothDomainFaces.push_back(EdgePatchQ);
    }
    int Size1 = SmoothDomainFaces.size();
    std::cout << "ADDED " << Size1 - Size0 << " edge domain patches"
              << std::endl;
  }

  void AddVertPatches() {

    // first find vert patch adjacency
    std::vector<std::vector<int>> VFPatch;
    std::vector<std::vector<int>> SeedF;

    for (size_t i = 0; i < VF.size(); i++) {
      std::set<int> VertPatchSet;
      for (size_t j = 0; j < VF[i].size(); j++) {
        int IndexF = VF[i][j];
        VertPatchSet.insert(QuadPatch[IndexF]);
      }
      bool IsPatchVert = false;
      if (VertPatchSet.size() > 2)
        IsPatchVert = true;
      if ((VertPatchSet.size() == 2) && (IsB[i]))
        IsPatchVert = true;
      if (IsPatchVert) {
        std::vector<int> VertPVect(VertPatchSet.begin(), VertPatchSet.end());
        VFPatch.push_back(VertPVect);
        SeedF.push_back(VF[i]);
      }
    }

    int Size0 = SmoothDomainFaces.size();

    // // add all adjacent patches
    for (size_t i = 0; i < VFPatch.size(); i++) {
      std::vector<int> FaceDomainIDX;
      for (size_t j = 0; j < VFPatch[i].size(); j++) {
        int IndexP = VFPatch[i][j];
        assert(IndexP >= 0);
        assert(IndexP < PatchQuad.size());
        FaceDomainIDX.insert(FaceDomainIDX.end(), PatchQuad[IndexP].begin(),
                             PatchQuad[IndexP].end());
      }
      SolvePatchIfNeeded(FaceDomainIDX, SeedF[i]);
      SmoothDomainFaces.push_back(FaceDomainIDX);
    }
    int Size1 = SmoothDomainFaces.size();
    std::cout << "ADDED " << Size1 - Size0 << " vertex domain patches"
              << std::endl;
  }

  void GetDomainKernel(const int IndexDomain,
                       std::vector<std::pair<int, int>> &VertKernel) const {
    std::vector<Geo::Point3<ScalarType>> PatchVertPos;
    std::vector<std::vector<int>> PatchFaces;
    std::vector<int> OriginalFaces, OriginalVert;

    GetDomainMesh(IndexDomain, PatchVertPos, PatchFaces, OriginalFaces,
                  OriginalVert);

    // compute geodesic from border
    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(PatchFaces, PatchVertPos.size(), VVAdj);
    Geo::TMark TempMark;
    typename Geodesic<ScalarType>::GeodesicParam GParam;
    getBorderVert(PatchFaces, GParam.Seeds);
    GParam.unitary_dist = true;
    Geodesic<ScalarType>::NodeDijkstra(PatchVertPos, VVAdj, TempMark, GParam);

    // return result
    VertKernel.clear();
    for (size_t i = 0; i < PatchVertPos.size(); i++) {
      int IndexOrV = OriginalVert[i];
      int Dist = (int)GParam.Dist[i];
      VertKernel.push_back(std::pair<int, int>(IndexOrV, Dist));
    }
  }

  void GetBestVertexKernel(std::vector<int> &VertKernel) const {
    VertKernel = std::vector<int>(QuadVertPos.size(), 0);
    for (size_t i = 0; i < SmoothDomainFaces.size(); i++) {
      std::vector<std::pair<int, int>> CurrKernel;
      GetDomainKernel(i, CurrKernel);
      for (size_t j = 0; j < CurrKernel.size(); j++) {
        int IndexV = CurrKernel[j].first;
        int KernV = CurrKernel[j].second;
        VertKernel[IndexV] = std::max(VertKernel[IndexV], KernV);
      }
    }
  }

  void LowKernelVert(std::vector<int> &LowLern, int MinKern = MIN_KERN) const {
    std::vector<int> VertKernel;
    GetBestVertexKernel(VertKernel);
    LowLern.clear();
    for (size_t i = 0; i < VertKernel.size(); i++) {
      if (VertKernel[i] < MinKern)
        LowLern.push_back(i);
    }
  }

  void AddKernelVert(std::vector<int> &VertIdx, int MaxKern = MAX_KERN) {
    ConvexPatchFinder<ScalarType> CPatch(QuadVertPos, QuadFaces);
    CPatch.Init();

    for (size_t i = 0; i < VertIdx.size(); i++) {
      int CurrV = VertIdx[i];
      std::vector<int> ConvexFaceIdx;
      CPatch.FindConvexPatchFrom(VF[CurrV], ConvexFaceIdx, MaxKern, false);
      SmoothDomainFaces.push_back(ConvexFaceIdx);
    }
  }

  void FindSmoothDomains_UsingQuadPatches() {

    auto t0 = clock_t::now();

    SmoothDomainFaces.clear();
    std::cout << "*** ADDING BASIC PATCHES ***" << std::endl;
    AddBasicPatches();
    std::cout << "*** ADDING EDGE PATCHES ***" << std::endl;
    AddEdgePatches();
    std::cout << "*** ADDING VERT PATCHES ***" << std::endl;
    AddVertPatches();
    std::cout << "*** Checking Remaining Vert ***" << std::endl;
    std::vector<int> LowLernV;
    LowKernelVert(LowLernV);
    std::cout << "There are " << LowLernV.size() << " Low Kernel V"
              << std::endl;
    if (LowLernV.size() > 0) {
      AddKernelVert(LowLernV);
      LowKernelVert(LowLernV);
      std::cout << "There are " << LowLernV.size() << " Low Kernel V"
                << std::endl;
    }
    std::cout << "*** Done ***" << std::endl;

    auto t1 = clock_t::now();
    ns t_delta{0};
    t_delta = std::chrono::duration_cast<ns>(t1 - t0);
    auto to_ms = [](ns d) {
      return std::chrono::duration<double, std::milli>(d).count();
    };

    std::cout << "Time   : " << to_ms(t_delta) << " ms\n";
  }

  void FindSmoothDomains_PerVertex() {

    auto t0 = clock_t::now();

    SmoothDomainFaces.clear();
    std::cout << "*** ADDING VERTEX CONVEX PATCHES ***" << std::endl;
    std::vector<int> KernVert;
    for (size_t i = 0; i < QuadVertPos.size(); i++)
      KernVert.push_back(i);

    AddKernelVert(KernVert);
    LowKernelVert(KernVert);
    std::cout << "There are " << KernVert.size() << " Low Kernel V"
              << std::endl;

    std::cout << "*** Done ***" << std::endl;

    auto t1 = clock_t::now();
    ns t_delta{0};
    t_delta = std::chrono::duration_cast<ns>(t1 - t0);
    auto to_ms = [](ns d) {
      return std::chrono::duration<double, std::milli>(d).count();
    };

    std::cout << "Time   : " << to_ms(t_delta) << " ms\n";
  }

  void FindVertOrderByDistortion(std::vector<int> &VertSeq) {
    VertSeq.clear();
    std::vector<ScalarType> FaceAR;
    AspectRatioScaledJacobian<ScalarType>(QuadVertPos, QuadFaces, FaceAR);

    std::vector<ScalarType> VertAr(QuadVertPos.size(), 0);
    std::vector<int> Num(QuadVertPos.size(), 0);
    for (size_t i = 0; i < QuadFaces.size(); i++)
      for (size_t j = 0; j < QuadFaces[i].size(); j++) {
        int IndexV = QuadFaces[i][j];
        VertAr[IndexV] += FaceAR[i];
        Num[IndexV]++;
      }
    for (size_t i = 0; i < VertAr.size(); i++)
      VertAr[i] /= Num[i];

    std::vector<std::pair<ScalarType, int>> RatioF;
    for (size_t i = 0; i < VertAr.size(); i++)
      RatioF.push_back(std::pair<ScalarType, int>(VertAr[i], i));

    std::sort(RatioF.begin(), RatioF.end());

    for (size_t i = 0; i < RatioF.size(); i++)
      VertSeq.push_back(RatioF[i].second);
  }

  void FindSmoothDomains_PerVertexIncremental(bool Sort_Distortion = true,
                                              int MaxKern = MAX_KERN) {

    auto t0 = clock_t::now();

    std::vector<bool> Solved(QuadVertPos.size(), false);
    SmoothDomainFaces.clear();

    ConvexPatchFinder<ScalarType> CPatch(QuadVertPos, QuadFaces);
    CPatch.Init();

    std::vector<int> VertSeq;
    if (Sort_Distortion) {
      FindVertOrderByDistortion(VertSeq);
      // std::vector<ScalarType> FaceAR;
      // AspectRatioScaledJacobian<ScalarType>(QuadVertPos, QuadFaces, FaceAR);

      // std::vector<ScalarType> VertAr(QuadVertPos.size(),0);
      // std::vector<int> Num(QuadVertPos.size(),0);
      // for (size_t i=0;i<QuadFaces.size();i++)
      // for (size_t j=0;j<QuadFaces[i].size();j++)
      // {
      //   int IndexV=QuadFaces[i][j];
      //   VertAr[IndexV]+=FaceAR[i];
      //   Num[IndexV]++;
      // }
      // for (size_t i=0;i<VertAr.size();i++)
      //   VertAr[i]/=Num[i];

      // std::vector<std::pair<ScalarType, int>> RatioF;
      // for (size_t i = 0; i < VertAr.size(); i++)
      //   RatioF.push_back(std::pair<ScalarType, int>(VertAr[i], i));

      // std::sort(RatioF.begin(), RatioF.end());

      // for (size_t i = 0; i < RatioF.size(); i++)
      //   VertSeq.push_back(RatioF[i].second);
    } else {
      for (size_t i = 0; i < QuadVertPos.size(); i++)
        VertSeq.push_back(i);
    }

    std::cout << "*** ADDING VERTEX CONVEX PATCHES ***" << std::endl;
    for (size_t i = 0; i < VertSeq.size(); i++) {
      int CurrV = VertSeq[i];
      assert(CurrV >= 0);
      assert(CurrV < QuadVertPos.size());
      if (Solved[CurrV])
        continue;

      std::vector<int> ConvexFaceIdx;
      CPatch.FindConvexPatchFrom(VF[CurrV], ConvexFaceIdx, MaxKern, false);
      SmoothDomainFaces.push_back(ConvexFaceIdx);

      int LastSmoothD = SmoothDomainFaces.size() - 1;
      std::vector<std::pair<int, int>> CurrKernel;
      GetDomainKernel(LastSmoothD, CurrKernel);
      for (size_t j = 0; j < CurrKernel.size(); j++) {
        int IndexV = CurrKernel[j].first;
        int KernV = CurrKernel[j].second;
        if (KernV >= MIN_KERN)
          Solved[IndexV] = true;
      }
    }
    std::vector<int> KernVert;
    LowKernelVert(KernVert);
    std::cout << "There are " << KernVert.size() << " Low Kernel V"
              << std::endl;

    std::cout << "*** Done ***" << std::endl;

    auto t1 = clock_t::now();
    ns t_delta{0};
    t_delta = std::chrono::duration_cast<ns>(t1 - t0);
    auto to_ms = [](ns d) {
      return std::chrono::duration<double, std::milli>(d).count();
    };

    std::cout << "Time   : " << to_ms(t_delta) << " ms\n";
  }

  void InitPatchQuad() {
    PatchQuad.clear();
    int MaxPIdx = -1;
    for (size_t i = 0; i < QuadPatch.size(); i++)
      MaxPIdx = std::max(MaxPIdx, QuadPatch[i]);
    assert(MaxPIdx >= 0);

    PatchQuad.resize(MaxPIdx + 1);

    for (size_t i = 0; i < QuadPatch.size(); i++) {
      int IndexP = QuadPatch[i];
      assert(IndexP < PatchQuad.size());
      int IndexF = i;
      PatchQuad[IndexP].push_back(IndexF);
    }
  }

  void SmoothWithPatches()
  {
    for (size_t i=0;i<SmoothDomainFaces.size();i++)
    {
      std::vector<Geo::Point3<ScalarType>> PatchVertPos;
      std::vector<std::vector<int>> PatchFaces;
      std::vector<int> OriginalFaces;
      std::vector<int> OriginalVert;

      std::cout<<"Step 0"<<std::endl;
      GetDomainMesh(i,PatchVertPos,PatchFaces,OriginalFaces,OriginalVert);
      std::vector<int> CornersIdx;
      QuadPatchFunctions<ScalarType>::FindCorners(PatchFaces, CornersIdx);

      std::cout<<"Step 1"<<std::endl;
      std::vector<std::pair<int, int>> Features;
      std::vector<int> Fixed;
      UVPatchSmoother<ScalarType>::SmoothPatch(PatchVertPos,PatchFaces,
                                               CornersIdx,Features,Fixed);
      
      std::cout<<"Step 2"<<std::endl;
      //then update positions on original mesh
      for (size_t j=0;j<OriginalVert.size();j++)
      {
        int IndexV=OriginalVert[j];
        QuadVertPos[IndexV]=PatchVertPos[j];
      }
      std::cout<<"Step 3"<<std::endl;
    }
  }

public:

  QuadPatchSmooth(std::vector<Geo::Point3<ScalarType>> &_QuadVertPos,
                  const std::vector<std::vector<int>> &_QuadFaces,
                  const std::vector<std::pair<int, int>> &_QuadFeatures,
                  const std::vector<int> &_Constraints)
      : QuadVertPos(_QuadVertPos), QuadFaces(_QuadFaces),
        QuadFeatures(_QuadFeatures), Constraints(_Constraints) {}

  void InitViaPatch(std::vector<int> &_QuadPatch) {
    QuadPatch = _QuadPatch;
    assert(QuadPatch.size() == QuadFaces.size());

    Geo::getVFAdjacency(QuadFaces, VF, VFIndex);
    Geo::getFFAdjacency(QuadFaces, FF, FFEdge);
    Geo::getBorderVert(QuadFaces, IsB);

    // then initialize the domains
    assert(IsAllQuad(QuadFaces));

    InitPatchQuad();

    FindSmoothDomains_UsingQuadPatches();

    CheckPatches();
  }

  void InitPerVertex() {

    Geo::getVFAdjacency(QuadFaces, VF, VFIndex);
    Geo::getFFAdjacency(QuadFaces, FF, FFEdge);
    Geo::getBorderVert(QuadFaces, IsB);

    // then initialize the domains
    assert(IsAllQuad(QuadFaces));

    // FindSmoothDomains_PerVertex();
    FindSmoothDomains_PerVertexIncremental();

    CheckPatches();

    SmoothWithPatches();
  }
};

} // namespace Geo

#endif
