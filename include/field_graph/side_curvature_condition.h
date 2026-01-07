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

#ifndef SIDE_CURVATURE_CONDITION
#define SIDE_CURVATURE_CONDITION

#include "../smooth_mesh.h"
#include "patch_managing.h"
#include "../Space_Query/triangle_grid_3D.h"
#include "../polyline3.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace Geo {

template <class ScalarType>
class SideCurvatureCondition : public PatchCondition<ScalarType> {
  std::vector<Geo::Point3<ScalarType>> TargetPos;
  std::vector<std::vector<int>> TargetFaces;
  Geo::TriangleSpatialIndex<ScalarType> *TargetTriGrid;

public:
  ScalarType MaxSideAngle;
  bool smooth_paths;

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;
    std::vector<std::vector<int>> SideV;
    PatchM.GetAllSideVert(IndexPatch, SideV);
    for (size_t i = 0; i < SideV.size(); i++) {

      Geo::Polyline3<ScalarType> SidePoly;
      for (size_t j = 0; j < SideV[i].size(); j++) {
        int IndexV = SideV[i][j];
        SidePoly.AddPoint(PatchM.PData.SubPatchVertPos[IndexPatch][IndexV]);
      }

// smooth and reproject
      if (smooth_paths)
      Geo::SmoothReproject<ScalarType>(SidePoly, 10, TargetPos, TargetFaces, *TargetTriGrid, 0.5);
      
      ScalarType TotA = SidePoly.TotAngle();
      if (fabs(TotA) > MaxSideAngle)
        return false;
    }
    return true;
  }

  std::string Name() const override { return std::string("side_curvature"); }

  void Init(const std::vector<Geo::Point3<ScalarType>> &_TargetPos,
            const std::vector<std::vector<int>> &_TargetFaces) {
    TargetPos = _TargetPos;
    TargetFaces = _TargetFaces;
    TargetTriGrid = new Geo::TriangleSpatialIndex<ScalarType>();
    TargetTriGrid->Init(TargetPos, TargetFaces);
  }

  //destructor
  ~SideCurvatureCondition() {
    if (TargetTriGrid)
      delete TargetTriGrid;
  }

  SideCurvatureCondition(ScalarType _MaxSideAngle) {
    MaxSideAngle = _MaxSideAngle;
    smooth_paths = false;
  }
};



template <class ScalarType>
class SideCurvatureUniformCondition : public PatchCondition<ScalarType> {
  std::vector<Geo::Point3<ScalarType>> TargetPos;
  std::vector<std::vector<int>> TargetFaces;
  Geo::TriangleSpatialIndex<ScalarType> *TargetTriGrid;

public:
  ScalarType MaxSideAngle;
  bool smooth_paths;

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;
    std::vector<std::vector<int>> SideV;
    PatchM.GetAllSideVert(IndexPatch, SideV);
    for (size_t i = 0; i < SideV.size(); i++) {

      Geo::Polyline3<ScalarType> SidePoly;
      for (size_t j = 0; j < SideV[i].size(); j++) {
        int IndexV = SideV[i][j];
        SidePoly.AddPoint(PatchM.PData.SubPatchVertPos[IndexPatch][IndexV]);
      }

// smooth and reproject
      if (smooth_paths)
      Geo::SmoothReproject<ScalarType>(SidePoly, 10, TargetPos, TargetFaces, *TargetTriGrid, 0.5);
      
      ScalarType TotA = SidePoly.TotAngle();
      if (fabs(TotA) > MaxSideAngle)
        return false;
    }
    return true;
  }

  std::string Name() const override { return std::string("side_curvature"); }

  void Init(const std::vector<Geo::Point3<ScalarType>> &_TargetPos,
            const std::vector<std::vector<int>> &_TargetFaces) {
    TargetPos = _TargetPos;
    TargetFaces = _TargetFaces;
    TargetTriGrid = new Geo::TriangleSpatialIndex<ScalarType>();
    TargetTriGrid->Init(TargetPos, TargetFaces);
  }

  //destructor
  ~SideCurvatureUniformCondition() {
    if (TargetTriGrid)
      delete TargetTriGrid;
  }

  SideCurvatureUniformCondition(ScalarType _MaxSideAngle) {
    MaxSideAngle = _MaxSideAngle;
    smooth_paths = false;
  }
};

} // namespace Geo
#endif
