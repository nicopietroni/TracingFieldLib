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

#ifndef TRI_REMESH
#define TRI_REMESH

#include "./local_operations/density_map.h"
#include "./local_operations/iterative_optimizer.h"
#include "./local_operations/tri_remesh_parameter.h"
#include "box3.h"
#include "eigen_interface.h"
#include "math_utils.h"
#include "mesh_manifold.h"
// #include "polyhedra_mesh.h"
#include "point3.h"
#include "quadric.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace Geo {

template <class ScalarType>
void AddAllBoundaryAsCorner(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            TriRemParam<ScalarType> &RemP) {
  std::vector<int> BoundaryV;
  Geo::getBorderVert(Faces, BoundaryV);
  for (size_t i = 0; i < BoundaryV.size(); i++)
    RemP.IsVertCorner[BoundaryV[i]] = true;
}

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

  if (RemP.IsFaceEFeature.size() == 0)
    RemP.IsFaceEFeature.resize(Faces.size(), std::vector<bool>(3, false));

  if (RemP.IsVertCorner.size() == 0)
    RemP.IsVertCorner.resize(VertPos.size(), false);
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

  if (TRemParam.fix_borders) {
    AddAllBoundaryAsCorner(VertPos, Faces, TRemParam);
  }

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

  // std::vector<ScalarType> VertDens;
  if (TRemParam.adaptive) {
    assert(TRemParam.VertSize.size() == VertPos.size());
    // t0 = std::chrono::high_resolution_clock::now();
    // ScalarType RadSize = 2 * AvgEdgeLen(VertPos, Faces);
    // FindDensityMap(VertPos, Faces, RadSize, VertDens);
    // ScalarType MinLVal = TRemParam.TargetL * TRemParam.minEdgeRatio;
    // ScalarType MaxLVal = TRemParam.TargetL * TRemParam.maxEdgeRatio;
    // LinearMap(VertDens, MinLVal, MaxLVal);
    // t1 = std::chrono::high_resolution_clock::now();
    // durationComputeAdaptLen +=
    //     std::chrono::duration_cast<std::chrono::microseconds>(t1 -
    //     t0).count();
  }
  for (int s = 0; s < TRemParam.steps; s++) {

    int stepsCollapseSmall = 0;
    int stepsRefineBig = 0;
    int stepsOptValence = 0;
    int stepsOprARatio = 0;

    // COLLAPSE SMALL EDGES
    if (TRemParam.collapse_small) {
      t0 = std::chrono::high_resolution_clock::now();
      // std::cout << "COLLAPSE SMALL EDGES" << std::endl;
      stepsCollapseSmall = CollapseSmallEdges(VertPos, Faces, VertQ, TRemParam,
                                              &TRemParam.VertSize);

      t1 = std::chrono::high_resolution_clock::now();
      durationCollapseSmall +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
    }

    if (TRemParam.adaptive) {
      assert(TRemParam.VertSize.size() == VertPos.size());
    }

    // REFINE BIG ONES
    if (TRemParam.refine_big) {
      t0 = std::chrono::high_resolution_clock::now();
      // std::cout << "REFINE BIG EDGES" << std::endl;
      stepsRefineBig =
          RefineBigEdges(VertPos, Faces, VertQ, TRemParam, &TRemParam.VertSize);
      t1 = std::chrono::high_resolution_clock::now();
      durationRefineBig +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
    }

    if (TRemParam.adaptive) {
      assert(TRemParam.VertSize.size() == VertPos.size());
    }

    // OPTIMIZE VALENCE
    if (TRemParam.flip_valence) {
      t0 = std::chrono::high_resolution_clock::now();
      // std::cout << "OPTIMIZE VALENCE" << std::endl;
      stepsOptValence = OptimizeValence(VertPos, Faces, VertQ, TRemParam);
      t1 = std::chrono::high_resolution_clock::now();
      durationOptValence +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
    }
    if (TRemParam.adaptive) {
      assert(TRemParam.VertSize.size() == VertPos.size());
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
                            TRemParam.feature_angle, TRemParam.IsFaceEFeature,
                            TRemParam.IsVertCorner, VertPos0, Faces0, TriGrid);

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
    SmoothReproject3DStep(VertPos, Faces, 3, TRemParam.feature_angle,
                          TRemParam.IsFaceEFeature, TRemParam.IsVertCorner,
                          VertPos0, Faces0, TriGrid);
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

template <class ScalarType>
int TriRemesh(Eigen::MatrixXd &vert,
              Eigen::MatrixXi &faces,
              TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) 
  {
    std::vector<Geo::Point3<ScalarType>> VertPos;
    std::vector<std::vector<int>> FacesIndex;
    // get the trimesh data from eigen matrices
    EigenInterface<ScalarType>::GetTriMeshData(faces, vert, VertPos,
                                               FacesIndex);
    // then remesh
    int steps = TriRemesh<ScalarType>(VertPos, FacesIndex, TRemParam);
    // transfer back to eigen
    EigenInterface<ScalarType>::GetTriMeshData(VertPos, FacesIndex, faces, vert);
    return steps;
  }

  template <class ScalarType>
int TriRemesh(Eigen::MatrixXf &vert,
              Eigen::MatrixXi &faces,
              TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) 
  {
    std::vector<Geo::Point3<ScalarType>> VertPos;
    std::vector<std::vector<int>> FacesIndex;
    // get the trimesh data from eigen matrices
    EigenInterface<ScalarType>::GetTriMeshData(faces, vert, VertPos,
                                               FacesIndex);
    // then remesh
    int steps = TriRemesh<ScalarType>(VertPos, FacesIndex, TRemParam);
    // transfer back to eigen
    EigenInterface<ScalarType>::GetTriMeshData(VertPos, FacesIndex, faces, vert);
    return steps;
  }

  template <class ScalarType>
  int TriRemesh(std::vector<std::vector<ScalarType>> & VertPos,
                std::vector<std::vector<int>> & Faces,
                TriRemParam<ScalarType> &TRemParam =
                    TriRemParam<ScalarType>()) 
    {
      std::vector<Geo::Point3<ScalarType>> VertPosVect;
      for (size_t i = 0; i < VertPos.size(); i++) {
        Geo::Point3<ScalarType> P;
        assert(VertPos[i].size() == 3);
        P.X = VertPos[i][0];
        P.Y = VertPos[i][1];
        P.Z = VertPos[i][2];
        VertPosVect.push_back(P);
      }

      // then remesh
      int steps = TriRemesh<ScalarType>(VertPosVect, Faces, TRemParam);

      // transfer back to vector of vectors
      VertPos.clear();
      for (size_t i = 0; i < VertPosVect.size(); i++) {
        Geo::Point3<ScalarType> P;
        P = VertPosVect[i];
        std::vector<ScalarType> PVect(3, 0);
        PVect[0] = P.X;
        PVect[1] = P.Y;
        PVect[2] = P.Z;
        VertPos.push_back(PVect);
      }
      return steps;
    }
  } // namespace Geo

#endif