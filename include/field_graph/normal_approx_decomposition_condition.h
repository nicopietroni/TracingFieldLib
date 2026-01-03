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

#ifndef NORMAL_APPROX_DECOMPOSITION_CONDITION
#define NORMAL_APPROX_DECOMPOSITION_CONDITION

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
class NormalApproxCondition : public PatchCondition<ScalarType> {

public:
  ScalarType MaxErr;
  ScalarType MaxErrPercent;
  ScalarType MaxBorderVariance;

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    // should have at least one border
    if (PatchM.NumBorders(IndexPatch) == 0)
      return false;

    std::vector<Geo::Point3<ScalarType>> SubPatchVertPos =
        PatchM.PData.SubPatchVertPos[IndexPatch];

    std::vector<std::vector<int>> SubPatchFaces =
        PatchM.PData.SubPatchFaces[IndexPatch];

    std::vector<Geo::Point3<ScalarType>> FaceNormalsTarget, VertNormalsTarget;

    // compute target normals
    ComputeNormals<ScalarType>(SubPatchVertPos, SubPatchFaces,
                               FaceNormalsTarget, VertNormalsTarget);
    std::vector<int> BorderV;
    getBorderVert(SubPatchFaces, BorderV);

    // check maximum angle between every possible pairs of border normals
    ScalarType MaxBorderAngle = 0;
    for (size_t i = 0; i < BorderV.size() - 1; i++) {
      for (size_t j = i + 1; j < BorderV.size(); j++) {
        ScalarType Ang = Geo::AngleDeg(VertNormalsTarget[BorderV[i]],
                                       VertNormalsTarget[BorderV[j]]);
        if (Ang > MaxBorderAngle)
          MaxBorderAngle = Ang;
      }
    }
    if (MaxBorderAngle > MaxBorderVariance)
      return false;

    // interpolate normals from boundary
    std::vector<Geo::Point3<ScalarType>> VertNormalsTest = VertNormalsTarget;
    Geo::SmoothLaplacianImplicit(VertNormalsTest, SubPatchFaces, BorderV);

    for (size_t i = 0; i < VertNormalsTest.size(); i++)
      VertNormalsTest[i].Normalize();

    // then find the number of faces with wrong normal
    int WrongVert = 0;
    for (size_t i = 0; i < VertNormalsTest.size(); i++) {
      if (Geo::AngleDeg(VertNormalsTest[i], VertNormalsTarget[i]) > MaxErr)
        WrongVert++;
    }

    if (MaxErrPercent == 0)
      return (WrongVert == 0);

    ScalarType PercWrong =
        (ScalarType)WrongVert / (ScalarType)SubPatchVertPos.size();
    return (PercWrong <= MaxErrPercent);
  }

  std::string Name() const override { return std::string("normal_approx"); }

  NormalApproxCondition(ScalarType _MaxErr, ScalarType _MaxErrPercent) {
    MaxErr = _MaxErr;
    MaxErrPercent = _MaxErrPercent;
    MaxBorderVariance = 80;
  }
};

} // namespace Geo
#endif
