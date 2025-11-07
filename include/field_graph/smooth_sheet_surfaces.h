#ifndef SMOOTH_SHEET_SURFACES
#define SMOOTH_SHEET_SURFACES

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "smooth_mesh.h"
#include <vector>

template <class ScalarType> class FoilsSmooth {

  static void SmoothReprojectSubMeshes(
      std::vector<Geo::Point3<ScalarType>> &NonManifV,
      std::vector<std::vector<int>> &NonManifFaces,
      const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
      const std::vector<std::vector<int>> &ReprojectFaces,
      const std::vector<std::pair<int, int>> &FeatureE,
      const std::vector<int> &Corners, int SmoothSteps) {

    // get corner and features position
    std::set<Geo::Point3<ScalarType>> FeatureVPos, CornerVPos;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    for (size_t i = 0; i < FeatureE.size(); i++) {
      int IndexV0 = FeatureE[i].first;
      int IndexV1 = FeatureE[i].second;
      Geo::Point3<ScalarType> Pos0 = ReprojectV[IndexV0];
      Geo::Point3<ScalarType> Pos1 = ReprojectV[IndexV1];
      FeatureVPos.insert(Pos0);
      FeatureVPos.insert(Pos1);
      Geo::Segment3<ScalarType> Seg(Pos0, Pos1);
      FeatureSegs.push_back(Seg);
    }

    for (size_t i = 0; i < Corners.size(); i++) {
      int IndexV = Corners[i];
      CornerVPos.insert(ReprojectV[IndexV]);
    }

    // then get the indexes on tne non manif
    std::set<int> FeatureIdxSet, CornerIdxSet;
    for (size_t i = 0; i < NonManifV.size(); i++) {
      Geo::Point3<ScalarType> TestP = NonManifV[i];
      if (FeatureVPos.count(TestP) > 0)
        FeatureIdxSet.insert(i);
      if (CornerVPos.count(TestP) > 0)
        CornerIdxSet.insert(i);
    }

    std::vector<int> FeatureIdxVert(FeatureIdxSet.begin(), FeatureIdxSet.end());
    std::vector<int> CornerIdxVert(CornerIdxSet.begin(), CornerIdxSet.end());

    std::vector<Geo::Point3<ScalarType>> PrevPos = NonManifV;

    // create the reprojection data structure
    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    TriGrid.Init(ReprojectV, ReprojectFaces);

    // get the ones to be reprojected
    std::vector<size_t> BorderV;
    Geo::getBorderVertNonManifMesh(NonManifFaces, BorderV);

    for (size_t s = 0; s < SmoothSteps; s++) {

      Geo::SmoothLaplacianNonManif<ScalarType>(NonManifV, NonManifFaces, 1);

      // reproject regular
      for (size_t i = 0; i < BorderV.size(); i++) {
        int IndexV = BorderV[i];

        if (FeatureIdxSet.count(IndexV))
          continue;
        if (CornerIdxSet.count(IndexV))
          continue;

        assert(IndexV >= 0);
        assert(IndexV < NonManifV.size());
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        int ClosestF;
        Geo::Point3<ScalarType> ClosestPos;
        TriGrid.GridClosest(ReprojectV, ReprojectFaces, Pos,
                            TriGrid.bbox.Diag(), ClosestF, ClosestPos);
        NonManifV[IndexV] = ClosestPos;
      }

      // reproject features
      for (size_t i = 0; i < FeatureIdxVert.size(); i++) {
        int IndexV = FeatureIdxVert[i];
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        ScalarType MinD = std::numeric_limits<ScalarType>::max();
        int BestSeg = -1;
        Geo::Point3<ScalarType> BestClos = Pos;
        for (size_t j = 0; j < FeatureSegs.size(); j++) {
          Geo::Segment3<ScalarType> Seg = FeatureSegs[j];
          Geo::Point3<ScalarType> TestClos = Geo::ClosestPoint(Seg, Pos);
          ScalarType TestD = (TestClos - Pos).Norm();
          if (TestD > MinD)
            continue;
          MinD = TestD;
          BestSeg = j;
          BestClos = TestClos;
        }
        if (BestSeg == -1)
          continue;
        NonManifV[IndexV] = BestClos;
      }

      // put back corners on original positions
      for (size_t i = 0; i < CornerIdxVert.size(); i++) {
        int IndexV = CornerIdxVert[i];
        NonManifV[IndexV] = PrevPos[IndexV];
      }
    }
  }

public:
  static void
  SmoothFoils(std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
              std::vector<std::vector<std::vector<int>>> &FoilFaces,
              const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
              const std::vector<std::vector<int>> &ReprojectFaces,
              const std::vector<std::pair<int, int>> &FeatureE,
              const std::vector<int> &Corners, int SmoothSteps) {

    // concatenate to a single surface
    std::vector<Geo::Point3<ScalarType>> NonManifV;
    std::vector<std::vector<int>> NonManifFaces;
    AppendMeshes<ScalarType>(NonManifV, NonManifFaces, FoilVertices, FoilFaces);
    RemoveDuplicatedVert(NonManifV, NonManifFaces);

    std::map<Geo::Point3<ScalarType>, int> PosIdxMap;
    for (size_t i = 0; i < NonManifV.size(); i++)
      PosIdxMap[NonManifV[i]] = i;

    SmoothReprojectSubMeshes(NonManifV, NonManifFaces, ReprojectV,
                             ReprojectFaces, FeatureE, Corners, SmoothSteps);

    // copy back positions
    for (size_t i = 0; i < FoilVertices.size(); i++)
      for (size_t j = 0; j < FoilVertices[i].size(); j++) {
        Geo::Point3<ScalarType> VertP = FoilVertices[i][j];
        assert(PosIdxMap.count(VertP) > 0);
        int IndexNonManif = PosIdxMap[VertP];
        assert(IndexNonManif >= 0);
        assert(IndexNonManif < NonManifV.size());
        FoilVertices[i][j] = NonManifV[IndexNonManif];
      }
  }
};

#endif