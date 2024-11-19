#ifndef PATCH_OPTIMIZE
#define PATCH_OPTIMIZE

#include "../tangent_space_smooth.h"
#include "path_functions.h"
#include "point3.h"
#include "segment3.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <limits>
#include <sys/_types/_size_t.h>
#include <vector>

namespace Geo {
template <class ScalarType> class PatchOptimize {

  // static Geo::Point3<ScalarType>
  // Reproject(const Geo::Point3<ScalarType> &TestP,
  //           const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //           const std::vector<std::pair<int, int>> &IndexE) {

  //   ScalarType minD = std::numeric_limits<ScalarType>::max();
  //   Geo::Point3<ScalarType> MinP = TestP;
  //   for (size_t i = 0; i < IndexE.size(); i++) {
  //     int IndexV0 = IndexE[i].first;
  //     int IndexV1 = IndexE[i].second;
  //     assert(IndexV0 != IndexV1);

  //     Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
  //     Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
  //     Geo::Segment3<ScalarType> Seg(P0, P1);
  //     Geo::Point3<ScalarType> Clos = Geo::ClosestPoint(Seg, TestP);
  //     ScalarType TestD = (Clos - TestP).Norm();
  //     if (TestD > minD)
  //       continue;
  //     minD = TestD;
  //     MinP = Clos;
  //   }
  //   return MinP;
  // }

  // static void SmoothSingPaths(
  //     std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const std::vector<std::vector<int>> &Connectivity,
  //     const std::vector<std::vector<int>> &NextF,
  //     const std::vector<std::vector<int>> &NextE,
  //     const std::vector<std::vector<int>> &BorderVertSeq,
  //     //   const std::vector<Geo::Segment3<ScalarType> > &BorderSegs,
  //     //   Geo::SegmentSpatialIndex<ScalarType> &BorderSegIndex,
  //     const std::vector<Geo::Point3<ScalarType>> &FaceNorm,
  //     const std::vector<Geo::Point3<ScalarType>> &VertNorm,
  //     const std::vector<std::vector<int>> &PathsV,
  //     std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>
  //         &ReprojMap,
  //     ScalarType Damp = 0.5, size_t smooth_steps = 3) {

  //   typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  //   // UVP.TJunctions = TJunctions;
  //   std::set<int> SmoothSet;

  //   for (size_t i = 0; i < PathsV.size(); i++) {
  //     std::vector<Geo::Point3<ScalarType>> OldVPos = VertPos;
  //     if (PathsV[i].size() == 0)
  //       continue;
  //     UVP.SmoothOnly.clear();
  //     UVP.FeatureLines.clear();
  //     for (int j = 0; j < (int)PathsV[i].size() - 1; j++) {
  //       int IndexV0 = PathsV[i][j];
  //       UVP.SmoothOnly.push_back(IndexV0);
  //       int IndexV1 = PathsV[i][j + 1];
  //       std::pair<int, int> key(std::min(IndexV0, IndexV1),
  //                               std::max(IndexV0, IndexV1));
  //       UVP.FeatureLines.push_back(key);
  //     }
  //     UVP.SmoothOnly.push_back((int)PathsV[i].back());
  //     UVP.Steps = smooth_steps;
  //     // std::cout<<"Smoothing Path: "<<i<<std::endl;

  //     //   Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity,
  //     //   NextF,NextE,
  //     //                                               BorderVertSeq,FaceNorm,
  //     //                                               VertNorm,
  //     //                                               BorderSegs,BorderSegIndex,UVP);

  //     Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity, NextF,
  //                                                 NextE, BorderVertSeq,
  //                                                 FaceNorm, VertNorm, UVP);

  //     // then reproject
  //     for (int j = 0; j < (int)PathsV[i].size(); j++) {
  //       int IndexV = PathsV[i][j];
  //       std::pair<int, int> VertPath(IndexV, i);
  //       if (ReprojMap.count(VertPath) == 0)
  //         continue;

  //       Geo::Point3<ScalarType> ToReproj = VertPos[IndexV];
  //       VertPos[IndexV] = Reproject(ToReproj, OldVPos, ReprojMap[VertPath]);
  //     }
  //   }
  // }

  // static void
  // SmoothInternal(std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                const std::vector<std::vector<int>> &Connectivity,
  //                const std::vector<std::vector<int>> &NextF,
  //                const std::vector<std::vector<int>> &NextE,
  //                const std::vector<std::vector<int>> &BorderVertSeq,
  //                //  const std::vector<Geo::Segment3<ScalarType> > &BorderSegs,
  //                //  Geo::SegmentSpatialIndex<ScalarType> &BorderSegIndex,
  //                const std::vector<int> &Corners,
  //                const std::vector<Geo::Point3<ScalarType>> &FaceNorm,
  //                const std::vector<Geo::Point3<ScalarType>> &VertNorm,
  //                const std::vector<std::vector<int>> &PathsV,
  //                ScalarType Damp = 0.5, size_t smooth_steps = 3) {

  //   typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  //   for (size_t i = 0; i < PathsV.size(); i++) {
  //     if (PathsV[i].size() == 0)
  //       continue;
  //     for (int j = 0; j < (int)PathsV[i].size() - 1; j++) {
  //       int IndexV0 = PathsV[i][j];
  //       int IndexV1 = PathsV[i][j + 1];
  //       std::pair<int, int> key(std::min(IndexV0, IndexV1),
  //                               std::max(IndexV0, IndexV1));
  //       UVP.FeatureLines.push_back(key);
  //     }
  //   }
  //   // for (size_t i = 0; i < PathsV.size(); i++)
  //   //   for (size_t j = 0; j < PathsV[i].size(); j++)
  //   //     UVP.Fixed.push_back(PathsV[i][j]);
  //   for (size_t i = 0; i < Corners.size(); i++)
  //     UVP.Fixed.push_back(Corners[i]);

  //   UVP.Steps = smooth_steps;
  //   // Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity, NextF,
  //   //                                             NextE,BorderVertSeq,FaceNorm,
  //   //                                             VertNorm,BorderSegs,BorderSegIndex,
  //   //                                             UVP);
  //   Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity, NextF,
  //                                               NextE, BorderVertSeq, FaceNorm,
  //                                               VertNorm, UVP);
  // }

  // static void FindReprojectionBases(
  //     const std::vector<std::vector<int>> &PathsV,
  //     std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>
  //         &ReprojMap) {

  //   std::map<int, std::vector<std::pair<int, int>>> VertEdges;
  //   std::map<std::pair<int, int>, int> EdgesPaths;

  //   for (int i = 0; i < (int)PathsV.size(); i++) {
  //     for (int j = 0; j < (int)PathsV[i].size() - 1; j++) {
  //       int IndexV0 = PathsV[i][j];
  //       int IndexV1 = PathsV[i][j + 1];
  //       std::pair<int, int> EdgeKey(std::min(IndexV0, IndexV1),
  //                                   std::max(IndexV0, IndexV1));
  //       assert(EdgesPaths.count(EdgeKey) == 0);
  //       EdgesPaths[EdgeKey] = i;
  //       VertEdges[IndexV0].push_back(EdgeKey);
  //       VertEdges[IndexV1].push_back(EdgeKey);
  //     }
  //   }

  //   ReprojMap.clear();
  //   for (size_t i = 0; i < PathsV.size(); i++)
  //     for (size_t j = 0; j < PathsV[i].size(); j++) {
  //       int IndexV0 = PathsV[i][j];
  //       assert(VertEdges.count(IndexV0) > 0);
  //       std::vector<std::pair<int, int>> Edges = VertEdges[IndexV0];
  //       if (Edges.size() < 4)
  //         continue;
  //       for (size_t e = 0; e < Edges.size(); e++) {
  //         std::pair<int, int> KeyE = Edges[e];
  //         assert(EdgesPaths.count(KeyE) > 0);
  //         int IndexP = EdgesPaths[KeyE];

  //         // same Path then not use
  //         if (IndexP == i)
  //           continue;

  //         std::pair<int, int> VertPath(IndexV0, i);
  //         ReprojMap[VertPath].push_back(KeyE);
  //       }
  //     }
  // }

public:

  static void SmoothPaths(std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Connectivity,
                          const std::vector<std::vector<int>> &PathsV,
                          const std::vector<std::pair<int,int> > &Features,
                          ScalarType Damp = 0.5, size_t Steps = 5) {

    typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
    std::vector<Geo::TJunction> TJunctions;
    Geo::PathFunctions<ScalarType>::FindTJunctions(PathsV, TJunctions);

    for (size_t i = 0; i < PathsV.size(); i++) {
      if (PathsV[i].size() == 0)
        continue;
      for (int j = 0; j < (int)PathsV[i].size() - 1; j++) {
        int IndexV0 = PathsV[i][j];
        int IndexV1 = PathsV[i][j + 1];
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        UVP.FeatureLines.push_back(key);
      }
    }
    UVP.Steps = Steps;
    UVP.TJunctions = TJunctions;
    UVP.Damp = Damp;
    UVP.fixBorder = true;
    //fix the features
    for (size_t i=0;i<Features.size();i++)
    {
      UVP.Fixed.push_back(Features[i].first);
      UVP.Fixed.push_back(Features[i].second);
    }
    Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity, UVP);
  }

  // static void TwoStepSmooth(std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                           const std::vector<std::vector<int>> &Connectivity,
  //                           const std::vector<std::vector<int>> &PathsV,
  //                           ScalarType Damp = 0.5, size_t Steps = 5) {
  //   // do a first smooth step
  //   //  auto t0 = std::chrono::high_resolution_clock::now();
  //   //  SmoothPaths(VertPos,Connectivity,PathsV,Damp,Steps);
  //   //  auto t1 = std::chrono::high_resolution_clock::now();
  //   //  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t1
  //   //  - t0); std::cout << "@@@@ Time : " << duration.count() << "
  //   //  milliseconds" << std::endl;

  //   std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> ReprojMap;
  //   FindReprojectionBases(PathsV, ReprojMap);

  //   std::vector<int> Corners;
  //   Geo::PathFunctions<ScalarType>::FindCorners(VertPos, PathsV, Corners);
  //   // std::cout<<"There are "<<Corners.size()<<" corners"<<std::endl;

  //   std::vector<std::vector<int>> NextF, NextE;
  //   Geo::getFFAdjacency(Connectivity, NextF, NextE);
  //   std::vector<std::vector<int>> BorderVertSeq;
  //   getBorderSequncesVert(Connectivity, NextF, NextE, BorderVertSeq);
  //   std::vector<Geo::Point3<ScalarType>> FaceNorm, VertNorm;
  //   ComputeNormals<ScalarType>(VertPos, Connectivity, FaceNorm, VertNorm);

  //   // std::vector<Geo::Segment3<ScalarType> > BorderSegs;
  //   // Geo::SegmentSpatialIndex<ScalarType> BorderSegIndex;
  //   // InitFromMeshBorder<ScalarType>(VertPos,Connectivity,NextF,
  //   //                                NextE,BorderSegs,BorderSegIndex);

  //   for (size_t s = 0; s < Steps; s++) {

  //     auto t0 = std::chrono::high_resolution_clock::now();

  //     //   SmoothSingPaths(VertPos, Connectivity, NextF,NextE,
  //     //                  BorderVertSeq,BorderSegs,BorderSegIndex,
  //     //                  FaceNorm, VertNorm,PathsV, ReprojMap, Damp);
  //     SmoothSingPaths(VertPos, Connectivity, NextF, NextE, BorderVertSeq,
  //                     FaceNorm, VertNorm, PathsV, ReprojMap, Damp);

  //     auto t1 = std::chrono::high_resolution_clock::now();

  //     //   SmoothInternal(VertPos, Connectivity, NextF,NextE,
  //     //                  BorderVertSeq,BorderSegs,BorderSegIndex,
  //     //                  Corners, FaceNorm, VertNorm,PathsV, Damp);
  //     SmoothInternal(VertPos, Connectivity, NextF, NextE, BorderVertSeq,
  //                    Corners, FaceNorm, VertNorm, PathsV, Damp);

  //     auto t2 = std::chrono::high_resolution_clock::now();
  //     auto duration0 =
  //         std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
  //     auto duration1 =
  //         std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

  //     std::cout << "Time 0: " << duration0.count() << " milliseconds"
  //               << std::endl;
  //     std::cout << "Time 1: " << duration1.count() << " milliseconds"
  //               << std::endl;
  //   }
  // }
};

} // namespace Geo
#endif
