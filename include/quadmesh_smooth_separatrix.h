#include "IO/io_obj.h"
#include "connectivity.h"
#include "geodesic.h"
#include "laplacian_neighbors.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "point3.h"
#include "polyline3.h"
#include "pos3.h"
#include "quad3.h"
#include "smooth_mesh.h"
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>
#ifndef SEPARATRIX_QUAD_MESH_SMOOTH
#define SEPARATRIX_QUAD_MESH_SMOOTH

namespace Geo {

template <class ScalarType> class QuadMeshSeparatrixSmooth {

  static void
  SmoothPosSeparatrix(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                      const std::vector<int> &SeparatrixV,
                      const std::set<int> &FixedV,
                      std::vector<Geo::Point3<ScalarType>> &NewVertPos,
                      int smooth_steps = 3, ScalarType delta = 0.25) {
    Polyline3<ScalarType> poly_sep;
    for (size_t i = 0; i < SeparatrixV.size(); i++) {
      int IndexV = SeparatrixV[i];
      poly_sep.PolyPos.push_back(QuadVertPos[IndexV]);
    }

    poly_sep.Smooth(smooth_steps, delta);
    NewVertPos = poly_sep.PolyPos;
  }

  static void
  SmoothPosSeparatrixStep(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                          const std::vector<std::vector<int>> &QuadFaces,
                          const std::vector<std::vector<int>> &Separatrices,
                          const std::vector<std::pair<int, int>> &Features,
                          const std::set<int> &FixedV, int smooth_steps = 3,
                          ScalarType delta = 0.25) {
    std::vector<Geo::Point3<ScalarType>> AvgPos(
        QuadVertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));
    std::vector<ScalarType> SumN(QuadVertPos.size(), 0);

    std::vector<std::vector<int>> Neighbors;
    GetLaplacianNeighbours(QuadVertPos, QuadFaces, Neighbors, Features);

    // add separatrices term
    for (size_t i = 0; i < Separatrices.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> NewVertPos;
      SmoothPosSeparatrix(QuadVertPos, Separatrices[i], FixedV, NewVertPos,
                          smooth_steps, 0.25);

      bool isLoop = (Separatrices[i][0] == Separatrices[i].back());
      int startIdx = 0;
      int endIdx = Separatrices[i].size();
      if (!isLoop) {
        startIdx++;
        endIdx--;
      }
      for (size_t j = startIdx; j < endIdx; j++) {
        int IndexV = Separatrices[i][j];
        AvgPos[IndexV] += NewVertPos[j];
        SumN[IndexV]+=1;
      }
    }

    // add laplacian term if needed
    for (size_t i = 0; i < AvgPos.size(); i++) {
      // add separatrices term
      if (SumN[i] == 0) {
        for (size_t j = 0; j < Neighbors[i].size(); j++) {
          AvgPos[i] += QuadVertPos[Neighbors[i][j]];
          SumN[i]+=1;
        }
      }
    }

    // add the anti-shear term along diagonals
    std::vector<ScalarType> FaceAR;
    AspectRatioWorstAngle(QuadVertPos,QuadFaces,FaceAR);
    // for (size_t i=0;i<FaceAR.size();i++)
    //   FaceAR[i]=1-FaceAR[i];

    ScalarType MinWeightSh=0.05;
    for (size_t i = 0; i < QuadFaces.size(); i++) 
    for (size_t j = 0; j < QuadFaces[i].size(); j++) 
    {
      int IndexV0=QuadFaces[i][j];
      int IndexV1=QuadFaces[i][(j+1)%4];
      ScalarType CurrW=MinWeightSh+FaceAR[i];
      AvgPos[IndexV0] += QuadVertPos[IndexV1]*CurrW;
      SumN[IndexV0]+=CurrW;
    }

    for (size_t i = 0; i < AvgPos.size(); i++) {
      assert(SumN[i] > 0);
      if (FixedV.count(i) > 0)
        continue;
      Geo::Point3<ScalarType> TargetPos = AvgPos[i] / (ScalarType)SumN[i];
      QuadVertPos[i] = QuadVertPos[i] * (1 - delta) + TargetPos * delta;
    }
  }

public:
  struct SeparatrixSmoothParam {
    std::vector<int> Fixed;
    std::vector<std::pair<int, int>> FeatureLines;
    ScalarType Damp;
    size_t Steps;
    bool fixBorder;
    std::vector<Geo::Point3<ScalarType>> *ReprojVert;
    std::vector<std::vector<int>> *ReprojFaces;
    std::vector<std::pair<int, int>> ReprojFeatureLines;
    //int ReprojStep;
    bool post_project;

    void Clear() {
      Fixed.clear();
      FeatureLines.clear();
      ReprojVert = NULL;
      ReprojFaces = NULL;
      ReprojFeatureLines.clear();
    }

    SeparatrixSmoothParam() {
      Damp = 0.2;
      Steps = 5;
      fixBorder = false;
      //ReprojStep = 5;
      ReprojVert = NULL;
      ReprojFaces = NULL;
      post_project = true;
    }
  };

private:

  static void
  ReprojectIfNeeded(const SeparatrixSmoothParam &SParam,
                    std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                    const std::vector<int> &ReprojectVertIdx,
                    std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
                    Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {

    bool do_reproject =((SParam.ReprojVert != NULL) &&
                        (SParam.ReprojFaces != NULL));

    bool post_project = (do_reproject && SParam.post_project);
    
    if (post_project)
      Geo::AttractProject<ScalarType>(*SParam.ReprojVert, VertPos, Faces);

    if (do_reproject)
      Geo::ProjectOn<ScalarType>(VertPos, Faces, ReprojectVertIdx,
                                 *SParam.ReprojVert, *SParam.ReprojFaces,
                                 TriGrid, FeatureSegs, SegGrid);
  }

  static void
  InitGridIfNeeded(const SeparatrixSmoothParam &SParam,
                   Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                   std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
                   Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {

    bool do_reproject = ((SParam.ReprojVert != NULL) &&
                        (SParam.ReprojFaces != NULL));

    bool post_project = (do_reproject && SParam.post_project);

    if (!do_reproject)
      return;

    Geo::InitProjStructures<ScalarType>(*SParam.ReprojVert, *SParam.ReprojFaces,
                                        SParam.ReprojFeatureLines, TriGrid,
                                        FeatureSegs, SegGrid);
  }

public:
  static void Smooth(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                     const std::vector<std::vector<int>> &QuadFaces,
                     const SeparatrixSmoothParam &SParam) {
    std::vector<std::vector<int>> Separatrices;

    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    Geo::SegmentSpatialIndex<ScalarType> SegGrid;
    InitGridIfNeeded(SParam, TriGrid, FeatureSegs, SegGrid);
    
    std::vector<int> ReprojectVertIdx;
    Geo::GetToReprojectVert<ScalarType>(QuadVertPos, ReprojectVertIdx, FeatureSegs, SegGrid);

    QuadPatchFunctions<ScalarType>::FindSeparatrices(
        QuadVertPos, QuadFaces, SParam.FeatureLines, SParam.Fixed,
        Separatrices);

    std::set<int> FixedSet(SParam.Fixed.begin(), SParam.Fixed.end());

    if (SParam.fixBorder) {
      std::vector<int> BorderV;
      getBorderVert(QuadFaces, BorderV);
      FixedSet.insert(BorderV.begin(), BorderV.end());
    }

    for (size_t s = 0; s < SParam.Steps; s++) {
      SmoothPosSeparatrixStep(QuadVertPos, QuadFaces, Separatrices,
                              SParam.FeatureLines, FixedSet, 5, SParam.Damp);

      ReprojectIfNeeded(SParam, QuadVertPos, QuadFaces, TriGrid,ReprojectVertIdx, FeatureSegs,
                          SegGrid);
    }
    // ReprojectIfNeeded(SParam, QuadVertPos, QuadFaces, TriGrid, ReprojectVertIdx,FeatureSegs,
    //                   SegGrid);
  }
};

} // namespace Geo

#endif
