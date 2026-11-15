#ifndef FIXED_POLYLINE_SIDES_H
#define FIXED_POLYLINE_SIDES_H

// Decides which variables (sub-sides) of a SideDiscretization are fixed, by
// matching them with given polylines of samples (e.g. the border vertices of
// a quad region to be matched exactly; a closed polyline repeats its first
// point at the end). SideDiscretization only offers the mechanism
// (SetFixedVariable); it does not know about polylines, regions or sketches.
// A variable is fixed if all its points lie on one polyline (within
// TolRatio * average segment); with AllBorderSides a border variable (one
// patch side only) is also taken if within half a segment (this only tells
// WHICH polyline it lies on: the samples are still exactly the polyline
// vertices; inner borders opened along features are farther and stay free).
// Fixed value = number of polyline segments between its two extremes,
// samples = the polyline points in between. The extremes (patch corners) are
// snapped to the closest polyline sample, never two on the same one:
// CornerSnap (corner position -> sample position) must be applied to the
// resulting quad mesh (MeshQuadrangulator::QuadVertSnap).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include "patch_side_discretization.h"

namespace Geo {

template <class ScalarType> class FixedPolylineSides {
  typedef Geo::Point3<ScalarType> CoordType;

public:
  // returns the number of fixed variables
  static int
  Apply(Geo::SideDiscretization<ScalarType> &SideD,
        const std::vector<std::vector<CoordType>> &FixedPolylines,
        const ScalarType TolRatio, const bool AllBorderSides,
        std::map<CoordType, CoordType> &CornerSnap) {
    SideD.ClearFixedVariables();
    CornerSnap.clear();
    const int NumV = SideD.NumVars();
    if (FixedPolylines.empty())
      return 0;

    // arc length of the polylines, closed if the first point is repeated
    size_t NumPoly = FixedPolylines.size();
    std::vector<std::vector<ScalarType>> ArcL(NumPoly);
    std::vector<bool> Closed(NumPoly, false);
    ScalarType SumSeg = 0;
    int NumSeg = 0;
    for (size_t j = 0; j < NumPoly; j++) {
      const std::vector<CoordType> &Poly = FixedPolylines[j];
      ArcL[j].assign(Poly.size(), 0);
      for (size_t k = 1; k < Poly.size(); k++) {
        ScalarType L = (Poly[k] - Poly[k - 1]).Norm();
        ArcL[j][k] = ArcL[j][k - 1] + L;
        SumSeg += L;
        NumSeg++;
      }
      Closed[j] = (Poly.size() > 3) && (Poly.front() == Poly.back());
    }
    if (NumSeg == 0)
      return 0;
    ScalarType Tol = TolRatio * SumSeg / NumSeg;
    // to tell the region border from the inner borders along the features
    ScalarType BorderTol = std::max(Tol, (ScalarType)0.5 * SumSeg / NumSeg);

    // 1. the variables lying on a polyline
    std::vector<int> VarPoly(NumV, -1);
    int NumBorderVar = 0, NumBorderFixed = 0;
    for (int i = 0; i < NumV; i++) {
      std::vector<CoordType> SubS = SideD.PointsOfVariable(i);
      if (SubS.size() < 2)
        continue;
      // (a variable with coincident extremes cannot be measured)
      if (SubS.front() == SubS.back())
        continue;
      bool IsBorderVar = SideD.IsBorderVariable(i);
      if (IsBorderVar)
        NumBorderVar++;
      // the polyline with the smallest maximum distance
      int BestJ = -1;
      ScalarType BestMaxD = std::numeric_limits<ScalarType>::max();
      for (size_t j = 0; j < NumPoly; j++) {
        ScalarType MaxD = 0;
        for (size_t k = 0; k < SubS.size(); k++) {
          ScalarType D, Par;
          ProjectOnPolyline(FixedPolylines[j], ArcL[j], SubS[k], D, Par);
          MaxD = std::max(MaxD, D);
          if (MaxD >= BestMaxD)
            break;
        }
        if (MaxD < BestMaxD) {
          BestMaxD = MaxD;
          BestJ = (int)j;
        }
      }
      bool Take = (BestJ >= 0) &&
                  ((AllBorderSides && IsBorderVar && (BestMaxD <= BorderTol)) ||
                   (BestMaxD <= Tol));
      if (AllBorderSides && IsBorderVar && !Take)
        std::cout << "Fixed sides: border side " << i << " left free, "
                  << SubS.size() << " points, max distance from the region "
                     "border "
                  << BestMaxD << " (half segment " << BorderTol << ")"
                  << std::endl;
      if (Take) {
        VarPoly[i] = BestJ;
        if (IsBorderVar)
          NumBorderFixed++;
      }
    }
    if (AllBorderSides && (NumBorderFixed < NumBorderVar))
      std::cout << "Fixed sides: " << NumBorderVar - NumBorderFixed << " / "
                << NumBorderVar
                << " border sides not on the region border (inner borders, "
                   "e.g. along features), left free"
                << std::endl;

    // 2. the corners (extremes of these variables) snapped to a sample
    std::vector<std::map<CoordType, ScalarType>> PolyCorners(NumPoly);
    for (int i = 0; i < NumV; i++) {
      if (VarPoly[i] < 0)
        continue;
      int j = VarPoly[i];
      std::vector<CoordType> SubS = SideD.PointsOfVariable(i);
      CoordType Ext[2] = {SubS.front(), SubS.back()};
      for (int e = 0; e < 2; e++) {
        ScalarType D, Par;
        ProjectOnPolyline(FixedPolylines[j], ArcL[j], Ext[e], D, Par);
        PolyCorners[j][Ext[e]] = Par;
      }
    }
    std::vector<std::map<CoordType, int>> CornerSample(NumPoly);
    std::vector<bool> PolyOk(NumPoly, true);
    for (size_t j = 0; j < NumPoly; j++) {
      if (PolyCorners[j].empty())
        continue;
      const std::vector<CoordType> &Poly = FixedPolylines[j];
      // number of distinct samples
      int M = Closed[j] ? (int)Poly.size() - 1 : (int)Poly.size();
      // corners along the polyline
      std::vector<std::pair<ScalarType, CoordType>> Sorted;
      typename std::map<CoordType, ScalarType>::iterator ItC;
      for (ItC = PolyCorners[j].begin(); ItC != PolyCorners[j].end(); ItC++)
        Sorted.push_back(std::make_pair((*ItC).second, (*ItC).first));
      std::sort(Sorted.begin(), Sorted.end());
      // arc length of the M samples, as an open line. A closed loop is cut
      // halfway before the sample closest to its first corner: no corner can
      // fall across the cut, so a loop is handled exactly as an open line
      // (its samples renumbered from Shift)
      std::vector<ScalarType> SA(ArcL[j].begin(), ArcL[j].begin() + M);
      int Shift = 0;
      if (Closed[j]) {
        ScalarType L = ArcL[j].back();
        ScalarType BestD = std::numeric_limits<ScalarType>::max();
        for (int k = 0; k < M; k++) {
          ScalarType D = std::fabs(ArcL[j][k] - Sorted.front().first);
          D = std::min(D, L - D);
          if (D < BestD) {
            BestD = D;
            Shift = k;
          }
        }
        ScalarType PrevLen = (Shift == 0) ? (L - ArcL[j][M - 1])
                                          : (ArcL[j][Shift] - ArcL[j][Shift - 1]);
        ScalarType Origin = ArcL[j][Shift] - (ScalarType)0.5 * PrevLen;
        for (int k = 0; k < M; k++)
          SA[k] = std::fmod(ArcL[j][(Shift + k) % M] - Origin + 2 * L, L);
        for (size_t c = 0; c < Sorted.size(); c++)
          Sorted[c].first = std::fmod(Sorted[c].first - Origin + 2 * L, L);
        std::sort(Sorted.begin(), Sorted.end());
      }
      // closest sample (by arc length)
      std::vector<int> Idx(Sorted.size(), 0);
      for (size_t c = 0; c < Sorted.size(); c++) {
        ScalarType BestD = std::numeric_limits<ScalarType>::max();
        for (int k = 0; k < M; k++) {
          ScalarType D = std::fabs(SA[k] - Sorted[c].first);
          if (D < BestD) {
            BestD = D;
            Idx[c] = k;
          }
        }
      }
      // never two corners on the same sample: on a collision the one that
      // moves less goes to the free sample next to it (the previous one back,
      // or this one forward)
      auto ArcAt = [&](int K) -> ScalarType {
        return SA[std::max(0, std::min(K, M - 1))];
      };
      for (size_t c = 1; c < Idx.size(); c++) {
        if (Idx[c] > Idx[c - 1])
          continue;
        int Fwd = Idx[c - 1] + 1;
        ScalarType CostFwd = std::fabs(ArcAt(Fwd) - Sorted[c].first) -
                             std::fabs(ArcAt(Idx[c]) - Sorted[c].first);
        int Back = Idx[c] - 1;
        bool CanBack = (Idx[c] == Idx[c - 1]) &&
                       ((c >= 2) ? (Back > Idx[c - 2]) : (Back >= 0));
        if (CanBack) {
          ScalarType CostBack =
              std::fabs(ArcAt(Back) - Sorted[c - 1].first) -
              std::fabs(ArcAt(Idx[c - 1]) - Sorted[c - 1].first);
          if (CostBack < CostFwd) {
            Idx[c - 1] = Back;
            continue;
          }
        }
        Idx[c] = Fwd;
      }
      // pushed past the last sample: the last corners moved back
      if (Idx.back() > M - 1) {
        Idx.back() = M - 1;
        for (int c = (int)Idx.size() - 2; c >= 0; c--)
          if (Idx[c] >= Idx[c + 1])
            Idx[c] = Idx[c + 1] - 1;
        if (Idx.front() < 0)
          PolyOk[j] = false;
      }
      // back to the numbering of the polyline
      for (size_t c = 0; c < Idx.size(); c++)
        Idx[c] = (Idx[c] + Shift) % M;
      if (!PolyOk[j]) {
        std::cout << "WARNING: fixed polyline " << j << " has more corners ("
                  << Idx.size() << ") than samples, left free" << std::endl;
        continue;
      }
      for (size_t c = 0; c < Sorted.size(); c++) {
        CornerSample[j][Sorted[c].second] = Idx[c];
        CornerSnap[Sorted[c].second] = Poly[Idx[c]];
      }
    }

    // 3. value and samples of each variable
    for (int i = 0; i < NumV; i++) {
      if (VarPoly[i] < 0)
        continue;
      int j = VarPoly[i];
      if (!PolyOk[j])
        continue;
      const std::vector<CoordType> &Poly = FixedPolylines[j];
      int M = Closed[j] ? (int)Poly.size() - 1 : (int)Poly.size();
      std::vector<CoordType> SubS = SideD.PointsOfVariable(i);
      int K0 = CornerSample[j][SubS.front()];
      int K1 = CornerSample[j][SubS.back()];
      bool Forward = (K1 > K0);
      int Count = std::abs(K1 - K0);
      if (Closed[j]) {
        // direction: the middle of the side comes before its end going
        // forward from its start
        ScalarType L = ArcL[j].back();
        ScalarType D, S0, S1, SM;
        ProjectOnPolyline(Poly, ArcL[j], SubS.front(), D, S0);
        ProjectOnPolyline(Poly, ArcL[j], SubS.back(), D, S1);
        CoordType Mid = (SubS.size() > 2)
                            ? SubS[SubS.size() / 2]
                            : (SubS.front() + SubS.back()) * (ScalarType)0.5;
        ProjectOnPolyline(Poly, ArcL[j], Mid, D, SM);
        ScalarType DMid = std::fmod(SM - S0 + L, L);
        ScalarType DEnd = std::fmod(S1 - S0 + L, L);
        Forward = (DMid < DEnd);
        Count = Forward ? ((K1 - K0 + M) % M) : ((K0 - K1 + M) % M);
      }
      if (Count < 1) {
        std::cout << "WARNING: fixed variable " << i
                  << " with no segment, left free" << std::endl;
        continue;
      }
      std::vector<CoordType> Samples;
      Samples.push_back(SubS.front());
      for (int step = 1; step < Count; step++) {
        int K = Forward ? (K0 + step) : (K0 - step);
        if (Closed[j])
          K = ((K % M) + M) % M;
        Samples.push_back(Poly[K]);
      }
      Samples.push_back(SubS.back());
      SideD.SetFixedVariable(i, Count, Samples);
    }
    std::cout << "Fixed variables: " << SideD.NumFixedVariables() << " / "
              << NumV << ", snapped corners: " << CornerSnap.size()
              << std::endl;
    return SideD.NumFixedVariables();
  }

private:
  // distance of P from the polyline and arc length parameter of its closest
  // point (ArcL: cumulative length of the points)
  static void ProjectOnPolyline(const std::vector<CoordType> &Poly,
                                const std::vector<ScalarType> &ArcL,
                                const CoordType &P, ScalarType &Dist,
                                ScalarType &Param) {
    Dist = std::numeric_limits<ScalarType>::max();
    Param = 0;
    if (Poly.size() == 1) {
      Dist = (P - Poly[0]).Norm();
      return;
    }
    for (size_t k = 0; k + 1 < Poly.size(); k++) {
      CoordType A = Poly[k];
      CoordType Dir = Poly[k + 1] - A;
      ScalarType L2 = Geo::Dot(Dir, Dir);
      ScalarType t = 0;
      if (L2 > 0)
        t = Geo::Dot(P - A, Dir) / L2;
      t = std::max((ScalarType)0, std::min((ScalarType)1, t));
      ScalarType D = (P - (A + Dir * t)).Norm();
      if (D < Dist) {
        Dist = D;
        Param = ArcL[k] + t * (ArcL[k + 1] - ArcL[k]);
      }
    }
  }
};

} // namespace Geo

#endif
