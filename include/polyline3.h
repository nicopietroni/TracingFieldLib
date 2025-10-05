#ifndef POLYLINE_3_H
#define POLYLINE_3_H

#include "point3.h"
#include "segment3.h"
#include <assert.h>
#include <cstddef>
#include <limits>
#include <vector>

namespace Geo {

template <class ScalarType> class Polyline3 {

public:
  std::vector<Geo::Point3<ScalarType>> PolyPos;

  ScalarType Lenght() {
    ScalarType l = 0;
    for (size_t i = 0; i < PolyPos.size() - 1; i++) {
      Geo::Point3<ScalarType> P0 = PolyPos[i];
      Geo::Point3<ScalarType> P1 = PolyPos[i + 1];
      l += (P0 - P1).Norm();
    }
    return l;
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

      std::vector<Geo::Point3<ScalarType>> Displ(SmoothPos.size(), Geo::Point3<ScalarType>(0, 0, 0));

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

} // namespace Geo

#endif
