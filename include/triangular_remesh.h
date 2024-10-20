#ifndef TRI_REMESH
#define TRI_REMESH

#include "./local_operations/density_map.h"
#include "./local_operations/iterative_optimizer.h"
#include "./local_operations/tri_remesh_parameter.h"
#include "box3.h"
#include "math_utils.h"
#include "mesh_manifold.h"
#include "polyhedra_mesh.h"
#include "quadric.h"
#include <algorithm>
#include <cassert>
#include <vector>

namespace Geo {

template <class ScalarType>
void CorrectDefaultValues(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          TriRemParam<ScalarType> &RemP) {
  if (RemP.steps <= 0)
    RemP.steps = 5;
  if (RemP.feature_angle <= 0)
    RemP.feature_angle = 180;

  if (RemP.max_error <= 0) {
    Geo::Box3<ScalarType> BB;
    BB.Init(VertPos);
    RemP.max_error = BB.Diag() * 0.01;
  }

  if (RemP.TargetL <= 0)
    RemP.TargetL = AvgEdgeLen(VertPos, Faces) * 1.5;
}

template <class ScalarType>
bool MatchTriRemeshPreconditions(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                 std::vector<std::vector<int>> &Faces) {
  bool ManifE = EdgeManifold<ScalarType>(VertPos, Faces);
  bool ManifV = VertManifold<ScalarType>(VertPos, Faces);
  return (ManifE && ManifV);
}

template <class ScalarType>
int TriRemesh(std::vector<Geo::Point3<ScalarType>> &VertPos,
              std::vector<std::vector<int>> &Faces,
              TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) {

  CorrectDefaultValues(VertPos, Faces, TRemParam);

  std::vector<std::pair<Geo::Quadric<ScalarType>, ScalarType>> VertQ;
  InitQuadrics(VertPos, Faces, VertQ);

  auto durationSmooth = 0;
  auto durationCollapseSmall = 0;
  auto durationRefineBig = 0;
  auto durationOptValence = 0;
  auto durationOptCurvFlip = 0;
  auto durationComputeAdaptLen = 0;

  // auto durationOprARatio = 0;
  // auto durationReproject = 0;

  std::vector<Geo::Point3<ScalarType>> VertPos0 = VertPos;
  std::vector<std::vector<int>> Faces0 = Faces;
  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(VertPos0, Faces0);

  int iter = 0;
  auto t0 = std::chrono::high_resolution_clock::now();
  auto t1 = std::chrono::high_resolution_clock::now();
  int all_steps = 0;

  std::vector<ScalarType> VertDens;
  if (TRemParam.adaptive) {
    t0 = std::chrono::high_resolution_clock::now();
    ScalarType RadSize = 2 * AvgEdgeLen(VertPos, Faces);
    FindDensityMap(VertPos, Faces, RadSize, VertDens);
    ScalarType MinLVal = TRemParam.TargetL * TRemParam.minEdgeRatio;
    ScalarType MaxLVal = TRemParam.TargetL * TRemParam.maxEdgeRatio;
    LinearMap(VertDens, MinLVal, MaxLVal);
    t1 = std::chrono::high_resolution_clock::now();
    durationComputeAdaptLen +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  }

  for (int s = 0; s < TRemParam.steps; s++) {

    int stepsCollapseSmall = 0;
    int stepsRefineBig = 0;
    int stepsOptValence = 0;
    int stepsOprARatio = 0;

    // COLLAPSE SMALL EDGES
    if (TRemParam.collapse_small) {
      t0 = std::chrono::high_resolution_clock::now();
      stepsCollapseSmall =
          CollapseSmallEdges(VertPos, Faces, VertQ, TRemParam, &VertDens);
      t1 = std::chrono::high_resolution_clock::now();
      durationCollapseSmall +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
    }
    
    if (TRemParam.adaptive) {
      assert(VertDens.size() == VertPos.size());
    }

    // REFINE BIG ONES
    if (TRemParam.refine_big) {
      t0 = std::chrono::high_resolution_clock::now();
      stepsRefineBig =
          RefineBigEdges(VertPos, Faces, VertQ, TRemParam, &VertDens);
      t1 = std::chrono::high_resolution_clock::now();
      durationRefineBig +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
    }

    if (TRemParam.adaptive) {
      assert(VertDens.size() == VertPos.size());
    }

    // OPTIMIZE VALENCE
    if (TRemParam.flip_valence) {
      t0 = std::chrono::high_resolution_clock::now();
      stepsOptValence = OptimizeValence(VertPos, Faces, VertQ, TRemParam);
      t1 = std::chrono::high_resolution_clock::now();
      durationOptValence +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
    }
    if (TRemParam.adaptive) {
      assert(VertDens.size() == VertPos.size());
    }

    // SMOOTH STEP
    if (TRemParam.smooth_mesh) {
      t0 = std::chrono::high_resolution_clock::now();
      // int min_smooth_steps = 1;
      // int max_smooth_steps = 3;
      int smooth_steps = 3;
      // if (!TRemParam.iterative) {
      //   smooth_steps = max_smooth_steps - s;
      //   smooth_steps = std::max(smooth_steps, min_smooth_steps);
      // }

      SmoothReproject3DStep(VertPos, Faces, smooth_steps,
                            TRemParam.feature_angle,
                            TRemParam.IsFaceEFeature,
                            TRemParam.IsVertCorner,
                            VertPos0, Faces0, TriGrid);

      t1 = std::chrono::high_resolution_clock::now();
      durationSmooth +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
    }

    int total_steps =
        stepsCollapseSmall + stepsRefineBig + stepsOprARatio + stepsOptValence;
    all_steps += total_steps;
    iter = s + 1;

    if (TRemParam.writeMSG) {
      std::cout << "** Performed Collapse Small: " << stepsCollapseSmall
                << std::endl;
      std::cout << "** Performed Refine Big: " << stepsRefineBig << std::endl;
      std::cout << "** Performed Optimize Ratio: " << stepsOprARatio
                << std::endl;
      std::cout << "** Performed Optimize Valence: " << stepsOptValence
                << std::endl;
    }

    if (total_steps == 0)
      break;
  }

  if (TRemParam.flip_quadrics) {
    t0 = std::chrono::high_resolution_clock::now();
    all_steps += OptimizeFlipQuadric(VertPos, Faces, TRemParam);
    t1 = std::chrono::high_resolution_clock::now();
    durationOptCurvFlip +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  }

  if (TRemParam.smooth_mesh) {
    t0 = std::chrono::high_resolution_clock::now();
    SmoothReproject3DStep(VertPos, Faces, 3, 
                        TRemParam.feature_angle,
                        TRemParam.IsFaceEFeature,
                          TRemParam.IsVertCorner, VertPos0,
                          Faces0, TriGrid);
    t1 = std::chrono::high_resolution_clock::now();
    durationSmooth +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  }
  // Output the duration
  if (TRemParam.writeMSG) {
    std::cout << "** Performed Iterations: " << iter << std::endl;

    std::cout << "** Time Smooth: " << (durationSmooth) / (ScalarType)1000000
              << " seconds" << std::endl;

    std::cout << "** Time Collapse Small: "
              << (durationCollapseSmall) / (ScalarType)1000000 << " seconds"
              << std::endl;

    std::cout << "** Time Refine Big: "
              << (durationRefineBig) / (ScalarType)1000000 << " seconds"
              << std::endl;

    // std::cout << "** Time Opt A Ratio: "
    //           << (durationOprARatio) / (ScalarType)1000000 << " seconds"
    //           << std::endl;

    std::cout << "** Time Opt Valence: "
              << (durationOptValence) / (ScalarType)1000000 << " seconds"
              << std::endl;

    std::cout << "** Time Flip Quadric: "
              << (durationOptCurvFlip) / (ScalarType)1000000 << " seconds"
              << std::endl;

    if (TRemParam.adaptive) {
      std::cout << "** Time Compute Adaptive: "
                << (durationComputeAdaptLen) / (ScalarType)1000000 << " seconds"
                << std::endl;
    }
  }
  return (all_steps);
}

} // namespace Geo

#endif