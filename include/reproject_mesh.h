#ifndef REPROJECT_MESH
#define REPROJECT_MESH

#include "Space_Query/segment_grid_3D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "box3.h"
#include "point3.h"
#include "segment3.h"
#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
void AttractProject(const std::vector<Geo::Point3<ScalarType>> &TargetVert,
                    std::vector<Geo::Point3<ScalarType>> &VertProj,
                    const std::vector<std::vector<int>> &FacesProj) {

  std::vector<Geo::Point3<ScalarType>> AvgDiff(
      VertProj.size(), Geo::Point3<ScalarType>(0, 0, 0));

  std::vector<int> AvgNum(VertProj.size(), 0);

  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(VertProj, FacesProj);

  for (size_t i = 0; i < TargetVert.size(); i++) {
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(VertProj, FacesProj, TargetVert[i], TriGrid.bbox.Diag(),
                        closestF, closestPt);

    Geo::Point3<ScalarType> Diff = TargetVert[i] - closestPt;

    for (size_t j = 0; j < FacesProj[closestF].size(); j++) {
      int IndexV = FacesProj[closestF][j];
      // Geo::Point3<ScalarType> TargetP=VertProj[IndexV]+Diff;
      AvgDiff[IndexV] += Diff;
      AvgNum[IndexV]++;
    }
  }

  // average it
  for (size_t i = 0; i < AvgNum.size(); i++) {
    if (AvgNum[i] == 0)
      continue;
    AvgDiff[i] /= AvgNum[i];
  }

  // smooth the difference
  int stepDiff = 3;
  for (size_t i = 0; i < stepDiff; i++) {
    std::vector<Geo::Point3<ScalarType>> SmDiff(
        AvgDiff.size(), Geo::Point3<ScalarType>(0, 0, 0));

    std::vector<int> SmNum(AvgDiff.size(), 0);

    // cumulate
    for (size_t i = 0; i < FacesProj.size(); i++)
      for (size_t j = 0; j < FacesProj[i].size(); j++) {
        int IndexV0 = FacesProj[i][j];
        int IndexV1 = FacesProj[i][(j + 1) % FacesProj[i].size()];
        SmDiff[IndexV0] += AvgDiff[IndexV1];
        SmNum[IndexV0]++;
      }

    // then average
    for (size_t i = 0; i < SmDiff.size(); i++) {
      if (SmNum[i] == 0)
        continue;
      AvgDiff[i] = SmDiff[i] / SmNum[i];
    }
  }

  for (size_t i = 0; i < AvgDiff.size(); i++)
    VertProj[i] += AvgDiff[i];
}

template <class ScalarType>
void InitProjStructures(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        const std::vector<std::pair<int, int>> &Features,
                        Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                        std::vector<Geo::Segment3<ScalarType>> &Segs,
                        Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {
  TriGrid.Init(VertPos, Faces);
  Geo::InitFromMeshFeatures<ScalarType>(VertPos, Features, Segs, SegGrid);
}

template <class ScalarType>
void ProjectOn(std::vector<Geo::Point3<ScalarType>> &VertPos,
               const std::vector<std::vector<int>> &Faces,
               const std::vector<std::pair<int,int>> &Features,
               const std::vector<Geo::Point3<ScalarType>> &TargetPos,
               const std::vector<std::vector<int>> &TargetFaces,
               Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
               const std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
               Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {

  std::vector<bool> FeatureVert(VertPos.size(), false);

  for (size_t i = 0; i < Features.size(); i++) {
    FeatureVert[Features[i].first] = true;
    FeatureVert[Features[i].second] = true;
  }

  for (size_t i = 0; i < VertPos.size(); i++) {
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    if (!FeatureVert[i]) {
      TriGrid.GridClosest(TargetPos, TargetFaces,
                          VertPos[i], TriGrid.bbox.Diag(), 
                          closestF, closestPt);
    } else {
      int closestS=-1;
      SegGrid.GridClosest(FeatureSegs, VertPos[i], 
                    SegGrid.bbox.Diag(), closestS,
                          closestPt);
      // ScalarType minD=std::numeric_limits<ScalarType>::max();
      // for (size_t j=0;j<FeatureSegs.size();j++)
      // {
      //   Geo::Point3<ScalarType> closestTest=ClosestPoint(FeatureSegs[j], VertPos[i]);
      //   if ((closestTest-VertPos[i]).Norm()>minD)continue;
      //   minD=(closestTest-VertPos[i]).Norm();
      //   closestS=j;
      //   closestPt=closestTest;
      // }
    }
    VertPos[i] = closestPt;
  }
}


template <class ScalarType>
void ReprojectBasis(std::vector<Geo::Point3<ScalarType>> &VertPos,
               const std::vector<Geo::Point3<ScalarType>> &TargetPos,
               const std::vector<std::vector<int>> &TargetFaces,
               std::vector<int> &IndexF,
               std::vector<Geo::Point3<ScalarType>> &BaryCoord) 
{
  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(TargetPos,TargetFaces);

  IndexF.clear();
  BaryCoord.clear();

  for (size_t i = 0; i < VertPos.size(); i++) {
    int closestF=-1;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(TargetPos, TargetFaces,
                          VertPos[i], TriGrid.bbox.Diag(), 
                          closestF, closestPt);
    assert(closestF!=-1);
    
    Geo::Point3<ScalarType> CurrBary;
    BarycentricCoords(TargetPos,TargetFaces,closestF,
                      closestPt,CurrBary.X,
                      CurrBary.Y,CurrBary.Z);
    
    IndexF.push_back(closestF);
    BaryCoord.push_back(CurrBary);
  }
}

} // namespace Geo
#endif
