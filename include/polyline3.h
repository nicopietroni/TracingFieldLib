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

#ifndef POLYLINE_3_H
#define POLYLINE_3_H

#include "Space_Query/triangle_grid_3D.h"
#include "best_rotation.h"
#include "point3.h"
#include "segment3.h"
#include <assert.h>
#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class Polyline3 {

public:
  std::vector<Geo::Point3<ScalarType>> PolyPos;

  void AddPoint(const Geo::Point3<ScalarType> &P) { PolyPos.push_back(P); }

  ScalarType Lenght() {
    ScalarType l = 0;
    for (size_t i = 0; i < PolyPos.size() - 1; i++) {
      Geo::Point3<ScalarType> P0 = PolyPos[i];
      Geo::Point3<ScalarType> P1 = PolyPos[i + 1];
      l += (P0 - P1).Norm();
    }
    return l;
  }

  void SmoothLaplacian(int iters = 3, ScalarType lambda = ScalarType(0.25)) {
    bool isLoop = (PolyPos.front() == PolyPos.back());
    if (isLoop && PolyPos.size() < 4)
      return;
    if ((!isLoop) && PolyPos.size() < 3)
      return;

    for (int it = 0; it < iters; ++it) {
      std::vector<Geo::Point3<ScalarType>> P = PolyPos;
      if (isLoop)
        P.pop_back(); // work on unique vertices

      const int n = (int)P.size();
      const int startIdx = isLoop ? 0 : 1;
      const int endIdx = isLoop ? n : n - 1;

      std::vector<Geo::Point3<ScalarType>> newP = P;

      for (int i = startIdx; i < endIdx; ++i) {
        int im = (i - 1 + n) % n;
        int ip = (i + 1) % n;

        Geo::Point3<ScalarType> lap = (P[im] + P[ip]) * ScalarType(0.5) - P[i];
        newP[i] = P[i] + lap * lambda;
      }

      // write back
      for (int i = 0; i < n; ++i)
        P[i] = newP[i];

      if (isLoop) {
        PolyPos.resize(n + 1);
        for (int i = 0; i < n; ++i)
          PolyPos[i] = P[i];
        PolyPos.back() = PolyPos.front();
      } else {
        // endpoints unchanged by construction
        for (int i = 0; i < n; ++i)
          PolyPos[i] = P[i];
      }
    }
  }

  void Smooth(int iters = 3, ScalarType damp = 0.25) {
    bool isLoop = (PolyPos[0] == PolyPos.back());
    if (isLoop && PolyPos.size() < 4)
      return;
    if ((!isLoop) && PolyPos.size() < 3)
      return;

    for (int it = 0; it < iters; ++it) {

      std::vector<Geo::Point3<ScalarType>> SmoothPos = PolyPos;

      // remove the last one if loops
      if (isLoop)
        SmoothPos.pop_back();

      std::vector<Geo::Point3<ScalarType>> Displ(
          SmoothPos.size(), Geo::Point3<ScalarType>(0, 0, 0));

      // Compute tangential displacements
      int n = SmoothPos.size();

      int startIdx = 1;
      int endIdx = SmoothPos.size() - 1;

      if (isLoop) {
        startIdx = 0;
        endIdx = SmoothPos.size();
      }

      for (int i = startIdx; i < endIdx; ++i) {
        int im = (i - 1 + n) % n; // previous index (circular)
        int ip = (i + 1) % n;     // next index (circular)

        // Classic 1D Laplacian target: midpoint of neighbors
        Geo::Point3<ScalarType> target =
            (SmoothPos[im] + SmoothPos[ip]) * ScalarType(0.5);
        Geo::Point3<ScalarType> d = target - SmoothPos[i];

        // Project displacement onto tangent direction (curve direction)
        Geo::Point3<ScalarType> t = (SmoothPos[ip] - SmoothPos[im]);
        t.Normalize();
        Geo::Point3<ScalarType> dpar = t * Dot(d, t);

        Displ[i] = dpar * damp;
      }

      // Compute mean displacement to remove global drift (zero-mean)
      Geo::Point3<ScalarType> meanD(0, 0, 0);
      for (int i = 0; i < n; ++i)
        meanD = meanD + Displ[i];

      meanD = meanD * (ScalarType(1) / ScalarType(n));

      // Apply displacement minus the mean
      for (int i = startIdx; i < endIdx; ++i) {
        PolyPos[i] = PolyPos[i] + (Displ[i] - meanD);
      }

      // Close the loop again (last == first)
      if (isLoop)
        PolyPos.back() = PolyPos[0];
    }
  }

  std::vector<Geo::Point3<ScalarType>> SampleUniformly(int n) {
    std::vector<Geo::Point3<ScalarType>> sampledPoints;
    if (PolyPos.empty() || n <= 0) {
      return sampledPoints;
    }

    // Compute segment lengths and total length
    std::vector<ScalarType> segmentLengths;
    ScalarType totalLength = 0;

    for (size_t i = 1; i < PolyPos.size(); ++i) {
      ScalarType len = (PolyPos[i] - PolyPos[i - 1]).Norm();
      segmentLengths.push_back(len);
      totalLength += len;
    }

    if (totalLength == 0) {
      // If all points are the same, return one of them.
      sampledPoints.push_back(PolyPos.front());
      return sampledPoints;
    }

    // Step size for uniform sampling
    ScalarType stepSize = totalLength / (n - 1);
    ScalarType accumulatedLength = 0;
    size_t segmentIndex = 0;
    ScalarType segmentStartLength = 0;

    // sampledPoints.push_back(PolyPos.front()); // First point

    for (int i = 1; i < n - 1; ++i) {
      ScalarType targetLength = i * stepSize;

      // Find the segment where the targetLength falls
      while (segmentIndex < segmentLengths.size() &&
             segmentStartLength + segmentLengths[segmentIndex] < targetLength) {
        segmentStartLength += segmentLengths[segmentIndex];
        ++segmentIndex;
      }

      if (segmentIndex >= segmentLengths.size()) {
        break;
      }

      // Interpolate between polyline.PolyPos[segmentIndex] and
      // polyline.PolyPos[segmentIndex + 1]
      ScalarType segmentFraction =
          (targetLength - segmentStartLength) / segmentLengths[segmentIndex];
      Geo::Point3<ScalarType> interpolatedPoint =
          PolyPos[segmentIndex] * (1 - segmentFraction) +
          PolyPos[segmentIndex + 1] * segmentFraction;

      sampledPoints.push_back(interpolatedPoint);
    }

    // sampledPoints.push_back(polyline.PolyPos.back()); // Last point
    return sampledPoints;
  }

  // compute the sum of angles between consecutive segments
  // considering the polyline as a sequence of segments or
  // it can be a loop when first and last point are the same
  ScalarType TotAngle() {
    ScalarType TotalAngle = 0;
    if (PolyPos.size() < 3)
      return TotalAngle;

    bool isLoop = (PolyPos[0] == PolyPos.back());
    size_t nSeg = PolyPos.size() - 2;
    if (isLoop)
      nSeg = PolyPos.size();

    for (size_t i = 0; i < nSeg; i++) {
      Geo::Point3<ScalarType> P0 = PolyPos[i];
      Geo::Point3<ScalarType> P1 = PolyPos[(i + 1) % PolyPos.size()];
      Geo::Point3<ScalarType> P2 = PolyPos[(i + 2) % PolyPos.size()];

      Geo::Point3<ScalarType> V1 = P1 - P0;
      Geo::Point3<ScalarType> V2 = P2 - P1;

      V1.Normalize();
      V2.Normalize();

      ScalarType Angle = fabs(Geo::AngleDeg(V1, V2));
      // std::cout << "Angle " << i << " : " << Angle << std::endl;
      TotalAngle += Angle;
    }
    return TotalAngle;
  }

  void DenoiseHighFrequencyClosed_TaubinBilateral(
      int iterations = 10,
      ScalarType lambda = ScalarType(0.33), // smoothing step
      ScalarType mu = ScalarType(
          -0.34), // inverse step (negative) -> zero shrink / zero drift
      ScalarType sigma_geom =
          ScalarType(0.05) // geometric bilateral weight (detail preservation)
  ) {
    const int n = (int)PolyPos.size();
    if (n < 4)
      return; // closed curve requires at least 3 distinct points + duplicated
              // last

    // Work only on the distinct vertices: [0 .. m-1]
    const int m = n - 1;

    auto sqr = [](ScalarType x) { return x * x; };

    std::vector<Geo::Point3<ScalarType>> temp(m);

    // Performs one bilateral Laplacian step with the given coefficient
    auto OneBilateralStep = [&](ScalarType step) {
      for (int i = 0; i < m; ++i) {
        // Circular indexing for closed curve
        int im1 = (i - 1 + m) % m;
        int ip1 = (i + 1) % m;

        // Geometric differences
        ScalarType d1 = (PolyPos[im1] - PolyPos[i]).Norm();
        ScalarType d2 = (PolyPos[ip1] - PolyPos[i]).Norm();

        // Bilateral weights: preserve sharp features and corners
        ScalarType w1 = std::exp(-sqr(d1) / sqr(sigma_geom));
        ScalarType w2 = std::exp(-sqr(d2) / sqr(sigma_geom));

        ScalarType wsum = w1 + w2;
        if (wsum <= ScalarType(0)) {
          temp[i] = PolyPos[i];
          continue;
        }

        // Bilateral Laplacian displacement
        Geo::Point3<ScalarType> lap = (((PolyPos[im1] - PolyPos[i]) * w1) +
                                       ((PolyPos[ip1] - PolyPos[i]) * w2)) /
                                      wsum;

        // Apply step
        temp[i] = PolyPos[i] + lap * step;
      }

      // Commit the update
      for (int i = 0; i < m; ++i)
        PolyPos[i] = temp[i];
    };

    // Taubin scheme: smooth (lambda) then unsmooth (mu) -> no shrinkage
    for (int it = 0; it < iterations; ++it) {
      OneBilateralStep(lambda);
      OneBilateralStep(mu);
    }

    // Restore closure explicitly
    PolyPos[n - 1] = PolyPos[0];
  }

  // void ParametrizeLenght(std::vector<ScalarType> &UVal)
  // {
  //     UVal.clear();
  //     ScalarType tot_l = Lenght();
  //     UVal.push_back(0);
  //     ScalarType curr_l = 0;
  //     for (size_t i = 0; i < PolyPos.size() - 1; i++)
  //     {
  //         Geo::Point2<ScalarType> P0 = PolyPos[i];
  //         Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
  //         curr_l += (P0 - P1).Norm();
  //         UVal.push_back(curr_l / tot_l);
  //     }
  //     UVal.back()=1;
  // }

  // bool Interpolate(const ScalarType &AlphaVal,
  //                  int &IndexSeg,
  //                  ScalarType &SubAlpha)
  // {
  //     assert(AlphaVal >= 0);
  //     assert(AlphaVal <= 1);
  //     assert(Pos.size() >= 2);
  //     if (AlphaVal == 0)
  //     {
  //         IndexSeg = 0;
  //         SubAlpha = 0;
  //         return;
  //     }
  //     if (AlphaVal == 1)
  //     {
  //         IndexSeg = PolyPos.size() - 1;
  //         SubAlpha = 1;
  //         return;
  //     }
  //     std::vector<ScalarType> UVal;
  //     ParametrizeLenght(UVal);

  //     //then search the right one
  //     for (size_t i = 0; i < UVal.size() - 1; i++)
  //     {
  //         ScalarType U0=UVal[i];
  //         ScalarType U1=UVal[i+1];
  //         if ((AlphaVal>=U0)&&(AlphaVal<=U1))
  //         {
  //             IndexSeg=i;
  //             ScalarType Int0=AlphaVal-U0;
  //             ScalarType Int1=U1-AlphaVal;
  //             assert(Int0>=0);
  //             assert(Int1>=0);
  //             SubAlpha=1-(Int0/(Int0+Int1));
  //             return true;
  //         }
  //     }
  //     return false;
  // }
};

template <class ScalarType>
void SmoothReproject(Polyline3<ScalarType> &Poly, int SmoothSteps,
                     const std::vector<Geo::Point3<ScalarType>> &TargetPos,
                     const std::vector<std::vector<int>> &TargetFaces,
                     Geo::TriangleSpatialIndex<ScalarType> &TargetTriGrid,
                     ScalarType damp = 0.5) {

  for (int it = 0; it < SmoothSteps; it++) {
    Poly.SmoothLaplacian(1, damp);
    for (size_t i = 0; i < Poly.PolyPos.size(); i++) {
      int closestF;
      Geo::Point3<ScalarType> closestPt = Poly.PolyPos[i];
      bool found = TargetTriGrid.GridClosest(
          TargetPos, TargetFaces, Poly.PolyPos[i], TargetTriGrid.bbox.Diag(),
          closestF, closestPt);
      if (found)
        Poly.PolyPos[i] = closestPt;
    }
  }
}
// template <class ScalarType>
//   Geo::Point3<ScalarType>
//   ClosestPoint(const Geo::Point3<ScalarType> &TestP) const {
//     ScalarType MinD = std::numeric_limits<ScalarType>::max();
//     Geo::Point3<ScalarType> ClosP(0, 0, 0);
//     for (size_t i = 0; i < PolyPos.size() - 1; i++) {
//       Geo::Segment3<ScalarType> CurrSeg(PolyPos[i], PolyPos[(i + 1)]);
//       Geo::Point3<ScalarType> ClosTest = ClosestPoint(CurrSeg, TestP);
//       ScalarType DistTest = (ClosTest - TestP).Norm();
//       if (DistTest > MinD)
//         continue;
//       MinD = DistTest;
//       ClosP = ClosTest;
//     }
//     return ClosP;
//   }
template <class ScalarType> class RotationHelpers {
public:
  static inline Eigen::Matrix<ScalarType, 3, 1>
  ToEigen(const Geo::Point3<ScalarType> &p) {
    return Eigen::Matrix<ScalarType, 3, 1>(p.X, p.Y, p.Z);
  }

  // Helper to compute best rotation matrix mapping a->b
  static inline Eigen::Matrix<ScalarType, 3, 3>
  RotationAroundAxis(const Geo::Point3<ScalarType> &axis_unit,
                     ScalarType angle_rad) {
    Eigen::Matrix<ScalarType, 3, 1> a = ToEigen(axis_unit);
    return Eigen::AngleAxis<ScalarType>(angle_rad, a).toRotationMatrix();
  }

  inline static ScalarType Clamp(const ScalarType &x, const ScalarType &a,
                                 const ScalarType &b) {
    return (x < a) ? a : (x > b) ? b : x;
  }

  // Helper to compute signed twist angle around tangent
  static ScalarType
  SignedTwistAngleAroundTangent(const Geo::Point3<ScalarType> &Nend_unit,
                                const Geo::Point3<ScalarType> &N1_unit,
                                const Geo::Point3<ScalarType> &Tend_unit) {
    // All inputs must be unit and Nend/N1 perpendicular to Tend
    // ScalarType c = Clamp<ScalarType>(Nend_unit.Dot(N1_unit), ScalarType(-1),
    //                                  ScalarType(1));
    // ScalarType s = Tend_unit.Dot(Nend_unit.Cross(N1_unit)); // signed sine

    ScalarType c =
        Clamp(Dot(Nend_unit, N1_unit), ScalarType(-1), ScalarType(1));
    ScalarType s = Dot(Tend_unit, Cross(Nend_unit, N1_unit)); // signed sine
    return std::atan2(s, c);                                  // [-pi, pi]
  }

  static bool OrthoUnitNormalToTangent(const Geo::Point3<ScalarType> &T_unit,
                                       Geo::Point3<ScalarType> &N_io,
                                       ScalarType eps = ScalarType(1e-12)) {
    // Remove tangent component
    N_io = N_io - T_unit * (Dot(T_unit, N_io));
    ScalarType lN = N_io.Norm();

    if (lN > eps) {
      // N_io /= lN;
      N_io.Normalize();
      return true;
    }

    // Fallback: build any vector orthogonal to T
    Geo::Point3<ScalarType> a = (std::abs(T_unit.X) < ScalarType(0.9))
                                    ? Geo::Point3<ScalarType>(1, 0, 0)
                                    : Geo::Point3<ScalarType>(0, 1, 0);

    N_io = Cross(T_unit, Cross(T_unit, a));
    lN = N_io.Norm();
    if (lN <= eps)
      return false;

    // N_io /= lN;
    N_io.Normalize();
    return true;
  }
};

template <class ScalarType>
inline bool InterpolateNormals_PT_WithTwistRedistribution(
    const std::vector<Geo::Point3<ScalarType>> &PolyPos,
    const Geo::Point3<ScalarType> &N0_in, const Geo::Point3<ScalarType> &N1_in,
    std::vector<Geo::Point3<ScalarType>> &OutNormals) {
  // Must be an open polyline
  assert(PolyPos.size() >= 2);
  if (PolyPos.front() == PolyPos.back())
    return false;

  const int n = (int)PolyPos.size();
  if (n < 2)
    return false;

  OutNormals.assign(n, Geo::Point3<ScalarType>(0, 0, 0));

  // ---------------------------------------------------------------------
  // Build unit tangents per segment
  // ---------------------------------------------------------------------
  std::vector<Geo::Point3<ScalarType>> Tseg(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    Geo::Point3<ScalarType> t = PolyPos[i + 1] - PolyPos[i];
    if (t.Norm() <= ScalarType(1e-12))
      return false;

    t.Normalize();
    Tseg[i] = t;
  }

  // ---------------------------------------------------------------------
  // Initialize normal at the start
  // ---------------------------------------------------------------------
  Geo::Point3<ScalarType> N0 = N0_in;
  if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(Tseg[0], N0))
    return false;

  OutNormals[0] = N0;

  // ---------------------------------------------------------------------
  // Parallel transport along the polyline
  // ---------------------------------------------------------------------
  for (int i = 0; i < n - 1; ++i) {
    Geo::Point3<ScalarType> t0 = Tseg[i];
    Geo::Point3<ScalarType> t1 = (i + 1 < n - 1) ? Tseg[i + 1] : Tseg[i];

    Eigen::Matrix3<ScalarType> R =
        Geo::PointMatch<ScalarType>::BestRotationMatrix(t0, t1);

    Geo::Point3<ScalarType> Ni = R * OutNormals[i];

    if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(t1, Ni))
      return false;

    OutNormals[i + 1] = Ni;
  }

  // ---------------------------------------------------------------------
  // Enforce end normal constraint
  // ---------------------------------------------------------------------
  Geo::Point3<ScalarType> Tend = Tseg[n - 2];

  Geo::Point3<ScalarType> N1 = N1_in;
  if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(Tend, N1))
    return false;

  Geo::Point3<ScalarType> Nend = OutNormals[n - 1];
  if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(Tend, Nend))
    return false;

  // Compute signed twist needed aSround Tend to align Nend to N1
  ScalarType twist = RotationHelpers<ScalarType>::SignedTwistAngleAroundTangent(
      Nend, N1, Tend);

  // ---------------------------------------------------------------------
  // Build arc-length parameterization
  // ---------------------------------------------------------------------
  std::vector<ScalarType> s(n, ScalarType(0));
  ScalarType total = ScalarType(0);

  for (int i = 1; i < n; ++i) {
    ScalarType ds = (PolyPos[i] - PolyPos[i - 1]).Norm();
    total += ds;
    s[i] = total;
  }

  if (total <= ScalarType(1e-12))
    return false;

  // ---------------------------------------------------------------------
  // Redistribute twist smoothly along the polyline
  // ---------------------------------------------------------------------
  for (int i = 0; i < n; ++i) {
    ScalarType u = s[i] / total; // from 0 at start to 1 at end
    ScalarType a = twist * u;

    Geo::Point3<ScalarType> ti = (i == n - 1) ? Tend : Tseg[std::min(i, n - 2)];

    Eigen::Matrix3<ScalarType> Rtw =
        RotationHelpers<ScalarType>::RotationAroundAxis(ti, a);

    // OutNormals[i] = RotationHelpers<ScalarType>::ApplyRotation(Rtw,
    // OutNormals[i]);
    OutNormals[i] = Rtw * OutNormals[i];

    if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(ti,
                                                               OutNormals[i]))
      return false;
  }

  return true;
}

template <class ScalarType> class CompositeCurve3 {

  std::vector<Geo::Polyline3<ScalarType>> Curves;

  bool IsLoop() const {
    if (Curves.size() < 1)
      return false;
    Geo::Point3<ScalarType> FirstP = Curves.front().PolyPos.front();
    Geo::Point3<ScalarType> LastP = Curves.back().PolyPos.back();
    return (FirstP == LastP);
  }

  // check last point of each curve matches first point of next curve
  // and each curve has at least 2 points
  bool IsProper() const {

    int numCurves = (int)Curves.size();

    for (int i = 0; i < numCurves; i++) {
      if (Curves[i].PolyPos.size() < 2)
        return false;
    }
    for (int i = 0; i < numCurves - 1; i++) {
      Geo::Point3<ScalarType> LastP = Curves[i].PolyPos.back();
      Geo::Point3<ScalarType> NextFirstP = Curves[i + 1].PolyPos.front();
      if (!(LastP == NextFirstP))
        return false;
    }
    return true;
  }

  void GetNextIndex(int &currCurveIdx, int &currPointIdx, int &nextCurveIdx,
                    int &nextPointIdx) const {
    assert(currCurveIdx >= 0);
    assert(currCurveIdx < (int)Curves.size());
    assert(currPointIdx >= 0);
    int sizeCurrCurve = (int)Curves[currCurveIdx].PolyPos.size();
    assert(currPointIdx < sizeCurrCurve);

    int numCurves = (int)Curves.size();
    if (currPointIdx == (sizeCurrCurve - 1)) {
      // move to next curve
      nextCurveIdx = (currCurveIdx + 1) % numCurves;
      // should not be the first as last matches first
      nextPointIdx = 1;
    } else {
      nextCurveIdx = currCurveIdx;
      nextPointIdx = currPointIdx + 1;
    }
  }

  void GetPrevIndex(int &currCurveIdx, int &currPointIdx, 
                    int &prevCurveIdx,int &prevPointIdx) const 
  {
    assert(currCurveIdx >= 0);
    assert(currCurveIdx < (int)Curves.size());
    assert(currPointIdx >= 0);
    int sizeCurrCurve = (int)Curves[currCurveIdx].PolyPos.size();
    assert(currPointIdx < sizeCurrCurve);

    int numCurves = (int)Curves.size();
    if (currPointIdx == 0) {
      // move to next curve
      prevCurveIdx = (currCurveIdx + numCurves - 1) % numCurves;
      // should not be the last as last matches first
      prevPointIdx = Curves[prevCurveIdx].PolyPos.size() - 2;
    } else {
      prevCurveIdx = currCurveIdx;
      prevPointIdx = currPointIdx - 1;
    }
  }

  Geo::Point3<ScalarType> GetPointAtIndex(int curveIdx, int pointIdx) const {
    assert(curveIdx >= 0);
    assert(curveIdx < (int)Curves.size());
    assert(pointIdx >= 0);
    int sizeCurrCurve = (int)Curves[curveIdx].PolyPos.size();
    assert(pointIdx < sizeCurrCurve);
    return Curves[curveIdx].PolyPos[pointIdx];
  }

  Geo::Point3<ScalarType> EsteemNormalAtIndex(int curveIdx,
                                              int pointIdx) const {
    // get all indexes
    int nextCurveIdx, nextPointIdx;
    GetNextIndex(curveIdx, pointIdx, nextCurveIdx, nextPointIdx);
    int prevCurveIdx, prevPointIdx;
    GetPrevIndex(curveIdx, pointIdx, prevCurveIdx, prevPointIdx);

    // then get points
    Geo::Point3<ScalarType> Pcurr = GetPointAtIndex(curveIdx, pointIdx);
    Geo::Point3<ScalarType> Pnext = GetPointAtIndex(nextCurveIdx, nextPointIdx);
    Geo::Point3<ScalarType> Pprev = GetPointAtIndex(prevCurveIdx, prevPointIdx);

    Geo::Point3<ScalarType> V1 = Pnext - Pcurr;
    V1.Normalize();

    Geo::Point3<ScalarType> V2 = Pcurr - Pprev;
    V2.Normalize();

    Geo::Point3<ScalarType> N = Cross(V1, V2);
    N.Normalize();
    return N;
  }
  // // compute normals using cross product of vectors between last of previous
  // // segment between pairs of vectors around each Vertex index
  // bool GetExtremesNormal(std::vector<Geo::Point3<ScalarType>> &VertexNormals)
  // {
  //   if (!IsLoop())
  //     return false;

  //   VertexNormals.clear();
  //   MakeCompleteVertexIndexes();

  //   for (size_t i = 0; i < VertexIndexes.size(); i++) {
  //     int vidx = VertexIndexes[i];

  //     // get the previous and the next index in the entire point list
  //     int prevIdx = vidx - 1;
  //     if (prevIdx < 0)
  //       prevIdx = Points.size() - 2; // loop case

  //     int nextIdx = vidx + 1;
  //     if (nextIdx >= (int)Points.size())
  //       nextIdx = 1; // loop case

  //     Geo::Point3<ScalarType> Pcurr = Points[vidx];
  //     Geo::Point3<ScalarType> Pprev = Points[prevIdx];
  //     Geo::Point3<ScalarType> Pnext = Points[nextIdx];

  //     Geo::Point3<ScalarType> V1 = Pcurr - Pprev;
  //     Geo::Point3<ScalarType> V2 = Pnext - Pcurr;

  //     V1.Normalize();
  //     V2.Normalize();

  //     Geo::Point3<ScalarType> N = Cross(V1, V2);
  //     N.Normalize();
  //     VertexNormals.push_back(N);
  //   }
  // }

  // compute normals using cross product of vectors between last of previous
  // segment between pairs of vectors around each Vertex index
  bool GetExtremesNormal(std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> &ExtremesNormals) const 
  {

    if (!IsLoop())
      return false;

    if (!IsProper())
      return false;

    ExtremesNormals.clear();

    for (size_t i = 0; i < Curves.size(); i++) {
      Geo::Point3<ScalarType> N0 = EsteemNormalAtIndex(i, 0);
      Geo::Point3<ScalarType> N1 = EsteemNormalAtIndex(i, (int)Curves[i].PolyPos.size() - 1);
      ExtremesNormals.push_back(std::make_pair(N0, N1));
    }
    return true;
  }

public:
  void AddPolyline(Geo::Polyline3<ScalarType> &PolyL) {
    Curves.push_back(PolyL);
  }

  bool ComputeVertexNormals(std::vector<std::vector<Geo::Point3<ScalarType>>> &OutNormals) const 
  {
    OutNormals.clear();
    std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        ExtremesNormals;
    bool res = GetExtremesNormal(ExtremesNormals);
    if (!res)
      return false;

    // then interpolate normals for each curve
    int numCurves = (int)Curves.size();
    OutNormals.resize(numCurves);
    for (int i = 0; i < numCurves; i++) {
      Geo::Point3<ScalarType> N0 = ExtremesNormals[i].first;
      Geo::Point3<ScalarType> N1 = ExtremesNormals[i].second;
      // OutNormals[i]=std::vector<Geo::Point3<ScalarType>>(Curves[i].PolyPos.size(),
      //                                                    Geo::Point3<ScalarType>(0,0,0));
      // OutNormals[i][0] = N0;
      // OutNormals[i].back() = N1;
      bool res2 = InterpolateNormals_PT_WithTwistRedistribution<ScalarType>(
          Curves[i].PolyPos, N0, N1, OutNormals[i]);
      if (!res2)
        return false;
    }
    return true;
  }
};

} // namespace Geo

#endif
