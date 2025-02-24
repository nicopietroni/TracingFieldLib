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

  std::vector<Geo::Point3<ScalarType> > PolyPos;
  // bool IsLoop;

  //void Add(const Geo::Point3<ScalarType> &P) { PolyPos.push_back(P); }

  // void SmoothBilaplacian() {
  //   // for (size_t i=0;i<BorderVIdx.size();i++)
  //   //     for (size_t j=0;j<BorderVIdx[i].size();j++)
  //   //     {
  //   //         int SizeB=BorderVIdx[i].size();
  //   //         int Ind0=(j+(SizeB-2))%SizeB;
  //   //         int Ind1=(j+(SizeB-1))%SizeB;
  //   //         int Ind2=j;
  //   //         int Ind3=(j+1)%SizeB;
  //   //         int Ind4=(j+2)%SizeB;

  //   //         int VInd0=BorderVIdx[i][Ind0];
  //   //         int VInd1=BorderVIdx[i][Ind1];
  //   //         int VInd2=BorderVIdx[i][Ind2];
  //   //         int VInd3=BorderVIdx[i][Ind3];
  //   //         int VInd4=BorderVIdx[i][Ind4];

  //   //         //Geo::Point3<ScalarType> Pos0=TargetVertPos[VInd1];
  //   //         Geo::Point3<ScalarType> Pos1=TargetVertPos[VInd2];
  //   //         //Geo::Point3<ScalarType> Pos2=TargetVertPos[VInd3];
  //   //         //if (FixedPosSet.count(Pos0)>0)continue;
  //   //         if (FixedPosSet.count(Pos1)>0)continue;
  //   //         //if (FixedPosSet.count(Pos2)>0)continue;

  //   //         for (size_t c=0;c<3;c++)
  //   //         {
  //   //             Geo::LinearConstraint<ScalarType> Constr;
  //   //             Constr.Weight=WeightBorderSmooth;

  //   //             ScalarType TotW=0;
  //   //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd0,c,0.1));
  //   //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd1,c,0.4));
  //   //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd2,c,-1));
  //   //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd3,c,0.4));
  //   //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd4,c,0.1));
  //   //             Constr.RHSide=0;
  //   //             LinConstr.push_back(Constr);
  //   //             numAdd++;
  //   //         }
  //   //     }
  // }

  // // Function to subdivide a polyline using Cubic B-Spline Subdivision
  // void Subdivide(bool IsLoop,
  //                 std::vector<bool> &ToSplit) {
  //     // Subdivide each edge
  //     std::vector<Geo::Point3<ScalarType> > new_points;
  //     for (size_t i = 0; i < PolyPos.size(); ++i) {

  //         // Original point P_i
  //         Geo::Point3<ScalarType> p0 = PolyPos[i];
  //         if (!ToSplit[i])continue;

  //         // Add the original point to new_points
  //         new_points.push_back(p0);

  //         // Create new points only if we are not at the last point
  //         int LastP=PolyPos.size()-1;
  //         if (IsLoop)LastP=PolyPos.size();
  //         if (i < LastP) {
  //             Geo::Point3<ScalarType> p1 = PolyPos[(i + 1)%PolyPos.size()];
  //             Geo::Point3<ScalarType> prev = (i > 0) ? PolyPos[i - 1] : p0;
  //             Geo::Point3<ScalarType> next = (i < result.size() - 2) ?
  //             result[i + 2] : p1;  // Handle boundary

  //             // New point: weighted average using B-Spline mask
  //             Point3<T> new_point = (prev * T(1) + p0 * T(6) + p1 * T(1)) /
  //             T(8);

  //             // Add the new point to the polyline
  //             new_points.push_back(new_point);
  //         }
  //     }

  //     // Update the polyline with new points
  //     result = new_points;
  // }

  // return result;
  // }

  ScalarType Lenght()
  {
      ScalarType l = 0;
      for (size_t i = 0; i < PolyPos.size() - 1; i++)
      {
          Geo::Point3<ScalarType> P0 = PolyPos[i];
          Geo::Point3<ScalarType> P1 = PolyPos[i + 1];
          l += (P0 - P1).Norm();
      }
      return l;
  }

  std::vector<Geo::Point3<ScalarType>> SampleUniformly( int n) {
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

    //sampledPoints.push_back(PolyPos.front()); // First point

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

        // Interpolate between polyline.PolyPos[segmentIndex] and polyline.PolyPos[segmentIndex + 1]
        ScalarType segmentFraction = (targetLength - segmentStartLength) / segmentLengths[segmentIndex];
        Geo::Point3<ScalarType> interpolatedPoint = 
            PolyPos[segmentIndex] * (1 - segmentFraction) + 
            PolyPos[segmentIndex + 1] * segmentFraction;

        sampledPoints.push_back(interpolatedPoint);
    }

    //sampledPoints.push_back(polyline.PolyPos.back()); // Last point
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
