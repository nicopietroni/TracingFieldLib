#ifndef PATCH_SIDE_DISCRETIZATION
#define PATCH_SIDE_DISCRETIZATION

#include "./solve/iterative_integer_lin_constraint_solver.h"
#include "mesh_functions.h"
#include "point3.h"
#include "polyline3.h"
#include "pos3.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class SideDiscretization {

public:
  enum typeMatching { NoQuad, OddSubd, NoMatch, OkMatch };

  enum energyMode { EM_ForceRegular, EM_Blend };

  enum solveMode { SM_Base, SM_Base_Regular, SM_Iterative, SM_GlobalSolve };

private:
  energyMode EM;
  ScalarType RegVsIso;

  struct PatchSide {
    std::vector<std::vector<Geo::Point3<ScalarType>>> SubSides;
    std::vector<std::array<int, 3>> PPAdj;
    std::vector<int> Variables;
  };

  struct PatchInfo {
    std::vector<PatchSide> Sides;
  };

  std::vector<PatchInfo> PatchSides;

  std::vector<std::vector<std::array<int, 3>>> VarSide;
  std::vector<int> VarTargetValue;
  std::vector<int> VarValue;
  bool force_parity;
  bool check_parity;
  bool current_doubleVal;
  // std::vector<bool> VarFeature;

  std::vector<std::vector<int>> VarSequences;
  int NumVariables;

  bool UseGeoAVG;

  std::vector<std::pair<std::vector<int>, std::vector<int>>> OppEquivalences;

  std::vector<std::vector<int>> PatchFacesIdx;
  std::vector<std::vector<int>> PatchCornersIdx;

public:
  // std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
  //     GeoFeatures;
  std::vector<std::vector<Geo::Point3<ScalarType>>> SampleSidePoints;
  std::vector<typeMatching> VarMatching;

  int NumPatches() const { return PatchSides.size(); }

  int NumSides(const int &PatchIdx) const {
    assert(PatchIdx >= 0);
    assert(PatchIdx < NumPatches());
    return PatchSides[PatchIdx].Sides.size();
  }

  int NumSubSides(const int &PatchIdx, const int &SideIdx) const {
    assert(PatchIdx >= 0);
    assert(PatchIdx < NumPatches());
    assert(SideIdx >= 0);
    assert(SideIdx < NumSides(PatchIdx));
    return (PatchSides[PatchIdx].Sides[SideIdx].SubSides.size());
  }

  int &Variable(int IndexP, int IndexS, int IndexSubS) {
    assert(IndexP >= 0);
    assert(IndexS >= 0);
    assert(IndexSubS >= 0);
    assert(IndexP < NumPatches());
    assert(IndexS < NumSides(IndexP));
    assert(IndexSubS < NumSubSides(IndexP, IndexS));
    return (PatchSides[IndexP].Sides[IndexS].Variables[IndexSubS]);
  }

  const int &Variable(int IndexP, int IndexS, int IndexSubS) const {
    assert(IndexP >= 0);
    assert(IndexS >= 0);
    assert(IndexSubS >= 0);
    assert(IndexP < NumPatches());
    assert(IndexS < NumSides(IndexP));
    assert(IndexSubS < NumSubSides(IndexP, IndexS));
    return (PatchSides[IndexP].Sides[IndexS].Variables[IndexSubS]);
  }

  std::vector<int> Variables(int IndexP, int IndexS) const {
    assert(IndexP >= 0);
    assert(IndexS >= 0);
    assert(IndexP < NumPatches());
    assert(IndexS < NumSides(IndexP));
    return (PatchSides[IndexP].Sides[IndexS].Variables);
  }

  std::vector<std::pair<int, int>> PatchesSideOfVariable(int IndexVar) const {
    assert(IndexVar >= 0);
    assert(IndexVar < VarSide.size());
    std::vector<std::pair<int, int>> PatchesSide;
    for (size_t i = 0; i < VarSide[IndexVar].size(); i++) {
      int IndexPatch = VarSide[IndexVar][i][0];
      int IndexSide = VarSide[IndexVar][i][1];
      std::pair<int, int> Entry(IndexPatch, IndexSide);
      PatchesSide.push_back(Entry);
    }
    return (PatchesSide);
  }

  void SideOfVariable(const int &VarIndex,
                      std::vector<std::array<int, 3>> &PatchSides) {
    PatchSides.clear();
    assert(VarIndex >= 0);
    assert(VarIndex < VarSide.size());
    PatchSides = VarSide[VarIndex];
  }

  void GetSidePoints(
      const int &PatchIdx, const int &SideIdx,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &SubSides) const {
    assert(PatchIdx >= 0);
    assert(PatchIdx < NumPatches());
    assert(SideIdx >= 0);
    assert(SideIdx < NumSides(PatchIdx));
    SubSides = PatchSides[PatchIdx].Sides[SideIdx].SubSides;
  }

  void GetSubSidePoints(const int &PatchIdx, const int &SideIdx,
                        const int &SubSideIdx,
                        std::vector<Geo::Point3<ScalarType>> &SubSide) const {
    assert(PatchIdx >= 0);
    assert(PatchIdx < NumPatches());
    assert(SideIdx >= 0);
    assert(SideIdx < NumSides(PatchIdx));
    SubSide = PatchSides[PatchIdx].Sides[SideIdx].SubSides[SubSideIdx];
  }

  void GetSubSidePoints(const std::array<int, 3> &PatchSide,
                        std::vector<Geo::Point3<ScalarType>> &SubSide) const {
    GetSubSidePoints(PatchSide[0], PatchSide[1], PatchSide[2], SubSide);
  }

  void GetVariablesMidP(std::vector<Geo::Point3<ScalarType>> &VarMidP) const {
    VarMidP.clear();
    for (size_t i = 0; i < VarSide.size(); i++) {
      assert(VarSide[i].size() > 0);
      std::array<int, 3> SideVar = VarSide[i][0];
      int IndexPatch = SideVar[0];
      int IndexSide = SideVar[1];
      int IndexSubSide = SideVar[2];
      std::vector<Geo::Point3<ScalarType>> SubSide;
      GetSubSidePoints(IndexPatch, IndexSide, IndexSubSide, SubSide);

      int IndexM0 = floor(SubSide.size() / 2) - 1;
      IndexM0 = std::max(IndexM0, 0);

      int IndexM1 = IndexM0 + 1;
      IndexM1 = std::min(IndexM1, (int)SubSide.size() - 1);

      Geo::Point3<ScalarType> P0 = SubSide[IndexM0];
      Geo::Point3<ScalarType> P1 = SubSide[IndexM1];
      VarMidP.push_back((P0 + P1) / 2);
    }
  }

  void GetVariablesValues(std::vector<int> &_VarValue) const {
    _VarValue = VarValue;
  }

private:
  void SortFromSeparator(std::vector<Geo::Point3<ScalarType>> &SideSeq,
                         const std::set<Geo::Point3<ScalarType>> &Separator) {
    // sort the sequence from the first corner
    int FistCorner = -1;
    for (size_t j = 0; j < SideSeq.size(); j++) {
      Geo::Point3<ScalarType> VertPos = SideSeq[j];
      if (Separator.count(VertPos) == 0)
        continue;
      FistCorner = j;
    }
    assert(FistCorner != -1);
    std::rotate(SideSeq.begin(), SideSeq.begin() + FistCorner, SideSeq.end());
  }

  void SplitIntoSubSequences(
      const std::vector<Geo::Point3<ScalarType>> &SideSeq,
      const std::set<Geo::Point3<ScalarType>> &Separator,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &SplitSeq,
      bool circular = true) {

    assert(Separator.count(SideSeq[0]) > 0);

    if (SideSeq.empty())
      return;

    std::vector<Geo::Point3<ScalarType>> TestSeq = SideSeq;
    // repeat the first to the enc
    TestSeq.push_back(TestSeq[0]);

    // add one sequence
    SplitSeq.resize(1);
    // set the first by defualt
    SplitSeq.back().push_back(SideSeq[0]);
    // start from second
    for (size_t i = 1; i < SideSeq.size(); ++i) {
      Geo::Point3<ScalarType> point = SideSeq[i];

      bool AddNewSide = (Separator.count(point) > 0);

      if (!circular)
        AddNewSide &= (i != (SideSeq.size() - 1));

      if (AddNewSide) {
        // repeat the separator at the end
        SplitSeq.back().push_back(point);
        // add new sequence
        SplitSeq.resize(SplitSeq.size() + 1);
      }

      // add the point
      SplitSeq.back().push_back(point);
    }
    // re add the first one at the end
    if (circular)
      SplitSeq.back().push_back(SideSeq[0]);
    // // If the first element is a separator, make it circular
    // if (!SplitSeq.empty() && Separator.find(SideSeq.front()) !=
    // Separator.end()) {
    //     SplitSeq.front().insert(SplitSeq.front().begin(), SideSeq.front());
    // }
  }

  void InitPatchAdj() {
    VarSide.clear();
    VarSide.resize(NumVariables);

    // resize with max value
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++) {
        int CurrSubSides = NumSubSides(i, j);
        PatchSides[i].Sides[j].PPAdj.resize(CurrSubSides,
                                            std::array<int, 3>{-1, -1, -1});

        for (size_t k = 0; k < CurrSubSides; k++) {
          int CurrVar = Variable(i, j, k);
          std::array<int, 3> Entry;
          Entry[0] = i;
          Entry[1] = j;
          Entry[2] = k;
          VarSide[CurrVar].push_back(Entry);
        }
      }

    for (size_t i = 0; i < VarSide.size(); i++) {
      std::vector<std::array<int, 3>> CurrSide = VarSide[i];

      // manifold check
      assert((CurrSide.size() == 1) || (CurrSide.size() == 2));

      // in this case is border
      if (CurrSide.size() == 1)
        continue;

      std::array<int, 3> Pside0 = CurrSide[0];
      int IndexPatch0 = Pside0[0];
      int IndexSide0 = Pside0[1];
      int IndexSubS0 = Pside0[2];

      std::array<int, 3> Pside1 = CurrSide[1];
      int IndexPatch1 = Pside1[0];
      int IndexSide1 = Pside1[1];
      int IndexSubS1 = Pside1[2];

      PatchSides[IndexPatch0].Sides[IndexSide0].PPAdj[IndexSubS0] = Pside1;
      PatchSides[IndexPatch1].Sides[IndexSide1].PPAdj[IndexSubS1] = Pside0;
    }

    // write a bit of stats
    int BorderSide = 0;
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++)
        for (size_t k = 0; k < NumSubSides(i, j); k++) {
          std::array<int, 3> Pside = PatchSides[i].Sides[j].PPAdj[k];
          int IndexPatch = Pside[0];
          int IndexSide = Pside[1];
          int IndexSubS = Pside[2];
          if (Pside == std::array<int, 3>{-1, -1, -1})
            BorderSide++;
        }
    std::cout << "Num Border Patch:" << BorderSide << std::endl;
  }

  void InitTargetValue(const ScalarType &TargetL) {
    VarValue.clear();
    for (size_t i = 0; i < NumVariables; i++) {
      std::vector<std::array<int, 3>> PatchSides;
      SideOfVariable(i, PatchSides);
      assert(PatchSides.size() > 0);
      std::array<int, 3> CurrPSide = PatchSides[0];
      std::vector<Geo::Point3<ScalarType>> SubSidePoints;
      GetSubSidePoints(CurrPSide, SubSidePoints);

      Geo::Polyline3<ScalarType> Poly3;
      Poly3.PolyPos = SubSidePoints;
      ScalarType L = Poly3.Lenght();
      int SubdNum = floor(0.5 + L / TargetL);
      if (force_parity)
        SubdNum = floor(0.5 + L / (TargetL * 2));
      SubdNum = std::max(SubdNum, 1);
      if (force_parity)
        VarValue.push_back(SubdNum * 2);
      else
        VarValue.push_back(SubdNum);
    }
    VarTargetValue = VarValue;
  }

  void InitVariables() {
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++) {
        int CurrSubNum = NumSubSides(i, j);
        PatchSides[i].Sides[j].Variables = std::vector<int>(CurrSubNum, -1);
      }

    std::map<std::vector<Geo::Point3<ScalarType>>, int> SideVarMap;
    int CurrVar = 0;
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++) {
        std::vector<std::vector<Geo::Point3<ScalarType>>> SubSides;
        GetSidePoints(i, j, SubSides);
        for (size_t k = 0; k < SubSides.size(); k++) {
          std::vector<Geo::Point3<ScalarType>> KeySide = SubSides[k];
          std::sort(KeySide.begin(), KeySide.end());
          if (SideVarMap.count(KeySide) == 0) {
            SideVarMap[KeySide] = CurrVar;
            Variable(i, j, k) = CurrVar;
            CurrVar++;
          } else {
            Variable(i, j, k) = SideVarMap[KeySide];
          }
        }
      }
    NumVariables = CurrVar;
    std::cout << "There are " << NumVariables << " Variables" << std::endl;
  }

  void FindDiscretizationPoint(const int IndexVariable) {
    std::vector<std::array<int, 3>> PatchSides;
    SideOfVariable(IndexVariable, PatchSides);
    assert(PatchSides.size() > 0);
    std::array<int, 3> CurrPSide = PatchSides[0];

    std::vector<Geo::Point3<ScalarType>> SubSide;
    GetSubSidePoints(CurrPSide[0], CurrPSide[1], CurrPSide[2], SubSide);

    Geo::Polyline3<ScalarType> PolyP;

    // //sample in the same way so I can recover the order
    // if (SubSide.back()<SubSide[0])
    //   std::reverse(SubSide.begin(),SubSide.begin());

    PolyP.PolyPos = SubSide;
    std::vector<Geo::Point3<ScalarType>> SidePos =
        PolyP.SampleUniformly(VarValue[IndexVariable] + 1);
    // PolyP.SampleUniformly(VarValue[IndexVariable] - 1);

    SampleSidePoints[IndexVariable].push_back(SubSide[0]);
    SampleSidePoints[IndexVariable].insert(
        SampleSidePoints[IndexVariable].end(), SidePos.begin(), SidePos.end());
    SampleSidePoints[IndexVariable].push_back(SubSide.back());
    // SideD.
  }

  void FindDiscretizationPoint() {
    SampleSidePoints.clear();
    SampleSidePoints.resize(NumVariables);
    for (size_t i = 0; i < NumVariables; i++)
      FindDiscretizationPoint(i);
  }

  bool IsQuadPatch(const int IndexP) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    return (PatchSides[IndexP].Sides.size() == 4);
  }

  int SideSize(const int IndexP, const int IndexSide) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    assert(IndexSide >= 0);
    assert(IndexSide < NumSides(IndexP));
    int CurrVal = 0;
    for (size_t i = 0; i < PatchSides[IndexP].Sides[IndexSide].Variables.size();
         i++) {
      int IndexVar = PatchSides[IndexP].Sides[IndexSide].Variables[i];
      assert(IndexVar >= 0);
      assert(IndexVar < VarValue.size());
      CurrVal += VarValue[IndexVar];
    }
    return CurrVal;
  }

  int PatchSubsides(const int IndexP) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    int SubVal = 0;
    for (size_t i = 0; i < NumSides(IndexP); i++)
      SubVal += SideSize(IndexP, i);
    return SubVal;
  }

  void PatchSideVariables(const int IndexP, int IndexSide,
                          std::vector<int> &VarList) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    assert(IndexSide >= 0);
    assert(IndexSide < NumSides(IndexP));
    for (size_t j = 0; j < NumSubSides(IndexP, IndexSide); j++)
      VarList.push_back(Variable(IndexP, IndexSide, j));
  }

  void PatchVariables(const int IndexP, std::vector<int> &VarList) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    int SubVal = 0;
    VarList.clear();
    for (size_t i = 0; i < NumSides(IndexP); i++)
      for (size_t j = 0; j < NumSubSides(IndexP, i); j++)
        VarList.push_back(Variable(IndexP, i, j));
  }

  void InitMatchingType() {
    VarMatching.clear();
    VarMatching.resize(NumVariables, OkMatch);

    for (size_t i = 0; i < NumPatches(); i++) {
      if (NumSides(i) != 4)
        continue;

      int Size0 = SideSize(i, 0);
      int Size1 = SideSize(i, 1);
      int Size2 = SideSize(i, 2);
      int Size3 = SideSize(i, 3);

      if (Size0 != Size2) {
        std::vector<int> VarList0;
        PatchSideVariables(i, 0, VarList0);
        std::vector<int> VarList2;
        PatchSideVariables(i, 2, VarList2);
        for (size_t j = 0; j < VarList0.size(); j++)
          VarMatching[VarList0[j]] = NoMatch;
        for (size_t j = 0; j < VarList2.size(); j++)
          VarMatching[VarList2[j]] = NoMatch;
      }

      if (Size1 != Size3) {
        std::vector<int> VarList1;
        PatchSideVariables(i, 1, VarList1);
        std::vector<int> VarList3;
        PatchSideVariables(i, 3, VarList3);
        for (size_t j = 0; j < VarList1.size(); j++)
          VarMatching[VarList1[j]] = NoMatch;
        for (size_t j = 0; j < VarList3.size(); j++)
          VarMatching[VarList3[j]] = NoMatch;
      }
    }

    for (size_t i = 0; i < NumPatches(); i++) {
      if (PatchSubsides(i) % 2 == 1) {
        std::vector<int> VarList;
        PatchVariables(i, VarList);
        for (size_t j = 0; j < VarList.size(); j++) {
          VarMatching[VarList[j]] = OddSubd;
          continue;
        }
      }
    }

    for (size_t i = 0; i < NumPatches(); i++) {
      if (!IsQuadPatch(i)) {
        std::vector<int> VarList;
        PatchVariables(i, VarList);
        for (size_t j = 0; j < VarList.size(); j++)
          VarMatching[VarList[i]] = NoQuad;
        continue;
      }
    }
  }

  void InitOppositeEquivalences() {
    OppEquivalences.clear();
    std::vector<bool> Visited(PatchSides.size(), false);
    for (size_t i = 0; i < PatchSides.size(); i++) {
      if (Visited[i])
        continue;
      Visited[i] = true;
      int SideSize = NumSides(i);
      if (SideSize != 4)
        continue;

      std::vector<int> VarSide0 = Variables(i, 0);
      std::vector<int> VarSide1 = Variables(i, 1);
      std::vector<int> VarSide2 = Variables(i, 2);
      std::vector<int> VarSide3 = Variables(i, 3);

      std::pair<std::vector<int>, std::vector<int>> Entry0(VarSide0, VarSide2);
      std::pair<std::vector<int>, std::vector<int>> Entry1(VarSide1, VarSide3);

      OppEquivalences.push_back(Entry0);
      OppEquivalences.push_back(Entry1);
    }
  }

  void InitVarEqSets() {
    VarSequences.clear();
    std::vector<std::vector<int>> EqGraph(NumVariables);

    for (size_t i = 0; i < OppEquivalences.size(); i++) {

      int Size0 = OppEquivalences[i].first.size();
      int Size1 = OppEquivalences[i].second.size();

      if (Size0 != 1)
        continue;
      if (Size1 != 1)
        continue;

      int IndexVar0 = OppEquivalences[i].first[0];
      int IndexVar1 = OppEquivalences[i].second[0];

      assert(IndexVar0 != IndexVar1);
      assert(IndexVar0 >= 0);
      assert(IndexVar1 >= 0);
      assert(IndexVar0 < NumVariables);
      assert(IndexVar1 < NumVariables);
      EqGraph[IndexVar0].push_back(IndexVar1);
      EqGraph[IndexVar1].push_back(IndexVar0);
    }

    std::vector<bool> Visited(EqGraph.size(), false);
    for (size_t i = 0; i < EqGraph.size(); i++) {
      if (Visited[i])
        continue;
      VarSequences.resize(VarSequences.size() + 1);
      std::vector<int> VertStack;
      VertStack.push_back(i);
      do {
        int IndexVar = VertStack.back();
        VertStack.pop_back();
        if (Visited[IndexVar])
          continue;

        VarSequences.back().push_back(IndexVar);
        Visited[IndexVar] = true;

        for (size_t j = 0; j < EqGraph[IndexVar].size(); j++) {
          int IndexNVar = EqGraph[IndexVar][j];
          if (Visited[IndexNVar])
            continue;
          VertStack.push_back(IndexNVar);
        }
      } while (!VertStack.empty());

      // std::set<int>
      // TestSet(VarSequences.back().begin(),VarSequences.back().end());
      //  if(VarSequences.back().size()!=TestSet.size())
      //  {
      //    std::cout<<"DEDE"<<std::endl;
      //    exit(0);
      //  }
    }

    // check non visited
    // int NonVisited = 0;
    for (size_t i = 0; i < Visited.size(); i++) {
      if (!Visited[i]) {
        // NonVisited++;
        VarSequences.push_back(std::vector<int>(1, i));
      }
    }
    // std::cout << "Non Visited Num " << NonVisited << std::endl;
  }

  void UpdateValueEqSet() {
    for (size_t i = 0; i < VarSequences.size(); i++) {
      int AvgSum = 0;
      ScalarType GeoAvg = 1;
      for (size_t j = 0; j < VarSequences[i].size(); j++) {
        int IndexVar = VarSequences[i][j];
        AvgSum += VarTargetValue[IndexVar];
        GeoAvg *= (ScalarType)VarTargetValue[IndexVar];
      }
      int VarNum = VarSequences[i].size();
      AvgSum = floor(0.5 + ((ScalarType)AvgSum / (ScalarType)VarNum));
      AvgSum = std::max(AvgSum, 1);
      GeoAvg = pow(GeoAvg, 1 / (ScalarType)VarNum);
      int GeoAvgi = floor(GeoAvg + 0.5);
      GeoAvgi = std::max(GeoAvgi, 1);

      if ((force_parity) && ((GeoAvgi % 2) == 1))
        GeoAvgi++;
      if ((force_parity) && ((AvgSum % 2) == 1))
        AvgSum++;
      for (size_t j = 0; j < VarSequences[i].size(); j++) {
        int IndexVar = VarSequences[i][j];
        if (UseGeoAVG)
          VarValue[IndexVar] = GeoAvgi;
        else
          VarValue[IndexVar] = AvgSum;
      }
    }
  }

  struct ChangeSubd {
    // the variables that are changed
    std::vector<int> IndexVar;

    // the difference of value
    int delta;

    // difference of patches that solve the parity
    int ParitySolve;
    // difference of patches that solve the matches
    int MatchSolve;
    // total missmatch
    int MissMatchDiff;
    // difference wth respect to ideal value
    ScalarType TargetSizeDiff;
    // difference wth respect to blend enedgy
    ScalarType BlendEnergyDiff;
    // use or not blend energy for priority
    bool UseBlend;

    inline bool operator==(const ChangeSubd &CSub) {
      if (IndexVar != CSub.IndexVar)
        return false;
      if (delta != CSub.delta)
        return false;
      return true;
    }

    inline bool operator<(const ChangeSubd &CSub) const {

      if (UseBlend)
        return (BlendEnergyDiff < BlendEnergyDiff);

      if (ParitySolve != CSub.ParitySolve)
        return (ParitySolve < CSub.ParitySolve);

      if (MatchSolve != CSub.MatchSolve)
        return (MatchSolve < CSub.MatchSolve);

      if (MissMatchDiff != CSub.MissMatchDiff)
        return (MissMatchDiff > CSub.MissMatchDiff);

      return (TargetSizeDiff > CSub.TargetSizeDiff);
    }

    ChangeSubd() {
      delta = 0;
      ParitySolve = 0;
      MatchSolve = 0;
      TargetSizeDiff = -1;
      BlendEnergyDiff = 0;
      UseBlend = false;
    }
  };

  int MissMatchValue(const int IndexP) const {
    int Size0 = SideSize(IndexP, 0);
    int Size1 = SideSize(IndexP, 1);
    int Size2 = SideSize(IndexP, 2);
    int Size3 = SideSize(IndexP, 3);
    int MissMatchVal = 0;
    MissMatchVal += fabs(Size0 - Size2);
    MissMatchVal += fabs(Size1 - Size3);
    return MissMatchVal;
  }

  int NumMissMatch(const int IndexP) const {
    int Size0 = SideSize(IndexP, 0);
    int Size1 = SideSize(IndexP, 1);
    int Size2 = SideSize(IndexP, 2);
    int Size3 = SideSize(IndexP, 3);
    int MissMatchNum = 0;
    if (Size0 != Size2)
      MissMatchNum++;
    if (Size1 != Size3)
      MissMatchNum++;
    return MissMatchNum;
  }

  int SizeDiff(const int IndexP) const {
    std::vector<int> VarList;
    PatchVariables(IndexP, VarList);
    int DiffVar = 0;
    for (size_t i = 0; i < VarList.size(); i++) {
      int CurrVar = VarList[i];
      DiffVar += fabs(VarTargetValue[CurrVar] - VarValue[CurrVar]);
    }
    return DiffVar;
  }

  bool IsOdd(const int IndexP) const {
    int Size0 = SideSize(IndexP, 0);
    int Size1 = SideSize(IndexP, 1);
    int Size2 = SideSize(IndexP, 2);
    int Size3 = SideSize(IndexP, 3);
    int SubdNum = Size0 + Size1 + Size2 + Size3;
    return ((SubdNum % 2) == 1);
  }

  int NumOdd(const std::vector<int> &Patches) const {
    int NumOdd = 0;
    for (size_t i = 0; i < Patches.size(); i++)
      if (IsOdd(Patches[i]))
        NumOdd++;
    return NumOdd;
  }

  int NumOdd() const {
    int NumOdd = 0;
    for (size_t i = 0; i < NumPatches(); i++)
      if (IsOdd(i))
        NumOdd++;
    return NumOdd;
  }

  int NumMissMatch() const {
    int NumMMatch = 0;
    for (size_t i = 0; i < NumPatches(); i++)
      NumMMatch += NumMissMatch(i);
    return NumMMatch;
  }

  int MissMatchValue() const {
    int MMValue = 0;
    for (size_t i = 0; i < NumPatches(); i++)
      MMValue += MissMatchValue(i);
    return MMValue;
  }

  int SizeDiff() const {
    int SizeDVal = 0;
    for (size_t i = 0; i < NumPatches(); i++)
      SizeDVal += SizeDiff(i);
    return SizeDVal;
  }

  int SizeDiff(const std::vector<int> &Patches) const {
    int DiffVal = 0;
    for (size_t i = 0; i < Patches.size(); i++)
      DiffVal += SizeDiff(Patches[i]);
    return DiffVal;
  }

  int NumMissMatch(const std::vector<int> &Patches) const {
    int NumMMatch = 0;
    for (size_t i = 0; i < Patches.size(); i++)
      NumMMatch += NumMissMatch(Patches[i]);
    return NumMMatch;
  }

  int MissMatchValue(const std::vector<int> &Patches) const {
    int MMVal = 0;
    for (size_t i = 0; i < Patches.size(); i++)
      MMVal += MissMatchValue(Patches[i]);
    return MMVal;
  }

  ScalarType BlendEnergy(const std::vector<int> &Patches) const {
    ScalarType EVal = 0;
    for (size_t i = 0; i < Patches.size(); i++)
      EVal += pow(MissMatchValue(Patches[i]), 2) * RegVsIso;

    for (size_t i = 0; i < Patches.size(); i++)
      EVal += pow(SizeDiff(Patches[i]), 2) * (1 - RegVsIso);

    return EVal;
  }

  ScalarType BlendEnergy() const {
    std::vector<int> PatchesIdx;
    for (size_t i = 0; i < NumPatches(); i++)
      PatchesIdx.push_back(i);

    return BlendEnergy(PatchesIdx);
  }

  std::vector<int> InvolvedPatches(ChangeSubd &ChOp) {
    std::vector<int> AllPatches;
    for (size_t i = 0; i < ChOp.IndexVar.size(); i++) {
      int CurrVar = ChOp.IndexVar[i];
      std::vector<std::pair<int, int>> PSide = PatchesSideOfVariable(CurrVar);
      for (size_t j = 0; j < PSide.size(); j++)
        AllPatches.push_back(PSide[j].first);
    }
    std::sort(AllPatches.begin(), AllPatches.end());
    AllPatches.erase(std::unique(AllPatches.begin(), AllPatches.end()),
                     AllPatches.end());
    return AllPatches;
  }

  void UpdatePriorityRegular(ChangeSubd &ChOp) {
    // get previous values
    std::vector<int> PatchIdx = InvolvedPatches(ChOp);

    int Parity0 = 0;
    if (check_parity)
      Parity0 = NumOdd(PatchIdx);

    int Match0 = NumMissMatch(PatchIdx);
    int MMatchDiff0 = MissMatchValue(PatchIdx);
    int SizeDiff0 = SizeDiff(PatchIdx);
    // get involved variables
    std::vector<int> OldVal; //=VarValue;
    // std::set<int> TestVarSet;
    for (size_t i = 0; i < ChOp.IndexVar.size(); i++) {
      int CurrVar = ChOp.IndexVar[i];
      // if (TestVarSet.count(CurrVar)>0)
      // {
      //   std::cout<<"ERROR"<<std::endl;
      //   exit(0);
      // }
      // TestVarSet.insert(CurrVar);
      OldVal.push_back(VarValue[CurrVar]);

      VarValue[CurrVar] += ChOp.delta;
      VarValue[CurrVar] = std::max(1, VarValue[CurrVar]);
    }

    int Parity1 = 0;
    if (check_parity)
      Parity1 = NumOdd(PatchIdx);
    // int Parity1=NumOdd();
    int Match1 = NumMissMatch(PatchIdx);
    int MMatchDiff1 = MissMatchValue(PatchIdx);
    int SizeDiff1 = SizeDiff(PatchIdx);

    // VarValue=OldVal;
    for (size_t i = 0; i < ChOp.IndexVar.size(); i++) {
      int CurrVar = ChOp.IndexVar[i];
      VarValue[CurrVar] = OldVal[i];
    }

    ChOp.ParitySolve = Parity0 - Parity1;
    ChOp.MatchSolve = Match0 - Match1;
    ChOp.MissMatchDiff = MMatchDiff0 - MMatchDiff1;
    ChOp.TargetSizeDiff = SizeDiff0 - SizeDiff1;
  }

  void UpdatePriorityBlend(ChangeSubd &ChOp) {
    // get previous values
    std::vector<int> PatchIdx = InvolvedPatches(ChOp);

    // int Parity0 = 0;
    // if (!force_parity)
    //   Parity0 = NumOdd(PatchIdx);

    ScalarType Energy0 = BlendEnergy(PatchIdx);

    // get involved variables
    std::vector<int> OldVal; //=VarValue;
    // std::set<int> TestVarSet;
    for (size_t i = 0; i < ChOp.IndexVar.size(); i++) {
      int CurrVar = ChOp.IndexVar[i];

      OldVal.push_back(VarValue[CurrVar]);

      VarValue[CurrVar] += ChOp.delta;
      VarValue[CurrVar] = std::max(1, VarValue[CurrVar]);
    }

    // int Parity1 = 0;
    // if (!force_parity)
    //   Parity1 = NumOdd(PatchIdx);

    ScalarType Energy1 = BlendEnergy(PatchIdx);

    for (size_t i = 0; i < ChOp.IndexVar.size(); i++) {
      int CurrVar = ChOp.IndexVar[i];
      VarValue[CurrVar] = OldVal[i];
    }

    ChOp.BlendEnergyDiff = Energy0 - Energy1;
    ChOp.UseBlend = true;
  }

  bool IsProfitable(ChangeSubd &ChOp) const {
    if (ChOp.UseBlend)
      return (ChOp.BlendEnergyDiff > 0);

    if (ChOp.ParitySolve < 0)
      return false;
    if (ChOp.ParitySolve > 0)
      return true;

    if (ChOp.MatchSolve < 0)
      return false;
    if (ChOp.MatchSolve > 0)
      return true;

    if (ChOp.MissMatchDiff < 0)
      return false;
    if (ChOp.MissMatchDiff > 0)
      return true;

    if (ChOp.TargetSizeDiff < 0)
      return false;
    if (ChOp.TargetSizeDiff > 0)
      return true;

    return false;
  }

  std::vector<ChangeSubd> ChangeOp;

  void DoOperation(const ChangeSubd &ChOp) {
    // std::cout << "Doing Operation" << std::endl;
    // std::cout << "Energy Reduction" << ChOp.BlendEnergyDiff<<std::endl;

    // std::cout << "- parity solve:" << ChOp.ParitySolve << std::endl;
    // std::cout<<"Doing Operation"<<ChOp.MatchSolve<<std::endl;
    // std::cout<<"Doing Operation"<<std::endl;
    // std::cout << "Energy 0: " << BlendEnergy() << std::endl;

    for (size_t i = 0; i < ChOp.IndexVar.size(); i++) {
      int CurrVar = ChOp.IndexVar[i];
      VarValue[CurrVar] += ChOp.delta;
      VarValue[CurrVar] = std::max(1, VarValue[CurrVar]);
    }
    // std::cout << "Energy 1: " << BlendEnergy() << std::endl;
  }

  void InitOpStack() {
    ChangeOp.clear();
    std::vector<int> DeltaOpt;
    DeltaOpt.push_back(1);
    DeltaOpt.push_back(-1);
    DeltaOpt.push_back(2);
    DeltaOpt.push_back(-2);

    if (EM == EM_ForceRegular) {
      for (size_t d = 0; d < DeltaOpt.size(); d++)
        for (size_t i = 0; i < VarSequences.size(); i++) {
          ChangeSubd CurrOp;
          CurrOp.IndexVar = VarSequences[i];
          CurrOp.delta = DeltaOpt[d];

          UpdatePriorityRegular(CurrOp);

          if (IsProfitable(CurrOp))
            ChangeOp.push_back(CurrOp);
        }
    } else {
      for (size_t d = 0; d < DeltaOpt.size(); d++)
        for (size_t i = 0; i < VarValue.size(); i++) {
          ChangeSubd CurrOp;
          CurrOp.IndexVar.push_back(i);
          CurrOp.delta = DeltaOpt[d];

          UpdatePriorityBlend(CurrOp);

          if (IsProfitable(CurrOp))
            ChangeOp.push_back(CurrOp);
        }
    }
    std::sort(ChangeOp.begin(), ChangeOp.end());
  }

  // void SetBaseIntegerValues(const ScalarType &TargetL) {
  //   InitTargetValue(TargetL);
  //   InitMatchingType();
  //   FindDiscretizationPoint();
  // }

  void MatchIntegerValues() {

    InitOppositeEquivalences();
    InitVarEqSets();
    UpdateValueEqSet();

    // if (force_parity) {
    //   for (size_t i = 0; i < VarValue.size(); i++)
    //     VarValue[i] *= 2;
    // }
    // InitMatchingType();
    // FindDiscretizationPoint();
  }

  void SolveIntegerValues() {
    int Odd0 = NumOdd();
    int MMatch0 = NumMissMatch();
    int MMVal0 = MissMatchValue();
    int DiffVal0 = SizeDiff();

    ScalarType En0 = BlendEnergy();

    if (force_parity) {
      for (size_t i = 0; i < VarValue.size(); i++) {
        assert(VarValue[i] % 2 == 0);
        VarValue[i] /= 2;
        VarTargetValue[i] /= 2;
      }
    }

    do {
      InitOpStack();
      // UpdatePriority(ChangeOp.back());
      check_parity = (!force_parity);
      if (ChangeOp.size() > 0)
        DoOperation(ChangeOp.back());
    } while (ChangeOp.size() > 0);

    if (force_parity) {
      for (size_t i = 0; i < VarValue.size(); i++)
        VarValue[i] *= 2;

      // force_parity = false;
      for (size_t i = 0; i < VarTargetValue.size(); i++)
        VarTargetValue[i] *= 2;

      if (EM == EM_ForceRegular) {
        // do one last step of optimization with no parity force
        check_parity = true;
        do {
          InitOpStack();
          // UpdatePriority(ChangeOp.back());
          if (ChangeOp.size() > 0)
            DoOperation(ChangeOp.back());
        } while (ChangeOp.size() > 0);
      }
    }

    InitMatchingType();
    FindDiscretizationPoint();
    ScalarType En1 = BlendEnergy();
    int Odd1 = NumOdd();
    int MMatch1 = NumMissMatch();
    int MMVal1 = MissMatchValue();
    int DiffVal1 = SizeDiff();
    if (EM == EM_ForceRegular) {
      std::cout << "* SIZE STATISTICS *" << std::endl;
      std::cout << "- Num Odd " << Odd0 << " - > " << Odd1 << std::endl;
      std::cout << "- Num MMatch " << MMatch0 << " - > " << MMatch1
                << std::endl;
      std::cout << "- Val MMatch " << MMVal0 << " - > " << MMVal1 << std::endl;
      std::cout << "- Diff Value " << DiffVal0 << " - > " << DiffVal1
                << std::endl;
    } else {
      std::cout << "* ENERGY STATISTICS *" << std::endl;
      std::cout << "- Energy " << En0 << " - > " << En1 << std::endl;
    }
  }

  // void SolveSequencesIntegerValues() {
  //   EM = EM_ForceRegular;
  //   SolveIntegerValues();
  // }

  // void SolveBlendIntegerValues(ScalarType _RegVsIso = 1) {
  //   EM = EM_Blend;
  //   RegVsIso = _RegVsIso;
  //   SolveIntegerValues();
  // }

  void SolveGlobal() {
    std::vector<LinearConstraint<ScalarType>> LSQCostr;

     if (force_parity) {
      for (size_t i = 0; i < VarValue.size(); i++) {
        assert(VarValue[i] % 2 == 0);
        VarValue[i] /= 2;
        VarTargetValue[i] /= 2;
      }
    }
    // add equality constraints
    //int MaxVal=0;
    for (size_t i = 0; i < OppEquivalences.size(); i++) {
      std::vector<int> Var0 = OppEquivalences[i].first;
      std::vector<int> Var1 = OppEquivalences[i].second;

      LinearConstraint<ScalarType> EqConstr;
      EqConstr.RHSide = 0;
      EqConstr.Weight = RegVsIso;

      for (size_t j = 0; j < Var0.size(); j++)
      {
        EqConstr.LHSide.push_back(VarTerm<ScalarType>(Var0[j], 0, 1));
        //MaxVal=std::max(MaxVal,Var0[j]);
      }

      for (size_t j = 0; j < Var1.size(); j++)
      {
        EqConstr.LHSide.push_back(VarTerm<ScalarType>(Var1[j], 0, -1));
        //MaxVal=std::max(MaxVal,Var0[j]);
      }
      LSQCostr.push_back(EqConstr);
    }

    // add isotropy constraint
    ScalarType delta_val=0.00001;
    for (size_t i = 0; i < VarTargetValue.size(); i++) {
      LinearConstraint<ScalarType> IsoConstr;
      IsoConstr.RHSide = VarTargetValue[i];
      IsoConstr.Weight = delta_val+(1 - RegVsIso);
      IsoConstr.LHSide.push_back(VarTerm<ScalarType>(i, 0, 1));
      LSQCostr.push_back(IsoConstr);
    }
    
    Geo::IntegerLinearConstraintLeastSQSolver<ScalarType> IntLSQ;
    IntLSQ.Init(LSQCostr, 1);
    //std::cout<<"DE 1 "<<VarTargetValue.size()<<std::endl;
    IntLSQ.Solve(VarValue);
    if (force_parity) {
      for (size_t i = 0; i < VarValue.size(); i++)
        VarValue[i] *= 2;
    }
    //std::cout<<"DE 2 "<<VarValue.size()<<std::endl;
    // EqConstr.RHSide=0;
    // EqConstr.Weight=RegVsIso;
    // for (size_t j=0;j<Var0.size();j++)
    //   EqConstr.LHSide.push_back(VarTerm<ScalarType>(Var0[j],0,1));
    // for (size_t j=0;j<Var1.size();j++)
    //   EqConstr.LHSide.push_back(VarTerm<ScalarType>(Var1[j],0,-1));
  }

public:
  void InitFromLabels(const std::vector<Geo::Point3<ScalarType>> &Vert,
                      const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::vector<int>> &_PatchFacesIdx,
                      const std::vector<std::vector<int>> &_PatchCornersIdx) {

    PatchSides.clear();
    VarSide.clear();
    VarTargetValue.clear();
    VarValue.clear();
    VarSequences.clear();
    NumVariables = 0;

    OppEquivalences.clear();
    PatchFacesIdx.clear();
    PatchCornersIdx.clear();

    // filters the ones that are not quads
    assert(_PatchFacesIdx.size() == _PatchCornersIdx.size());
    for (size_t i = 0; i < _PatchFacesIdx.size(); i++) {
      if (_PatchCornersIdx[i].size() != 4)
        continue;
      PatchFacesIdx.push_back(_PatchFacesIdx[i]);
      PatchCornersIdx.push_back(_PatchCornersIdx[i]);
    }
    std::set<Geo::Point3<ScalarType>> AllCornerPos;
    std::vector<std::set<Geo::Point3<ScalarType>>> PatchCornerPos;
    PatchCornerPos.resize(PatchCornersIdx.size());

    for (size_t i = 0; i < PatchCornersIdx.size(); i++)
      for (size_t j = 0; j < PatchCornersIdx[i].size(); j++) {
        int IndexV = PatchCornersIdx[i][j];
        assert(IndexV >= 0);
        assert(IndexV < Vert.size());
        Geo::Point3<ScalarType> CornPos = Vert[IndexV];
        AllCornerPos.insert(CornPos);
        PatchCornerPos[i].insert(CornPos);
      }

    PatchSides.resize(PatchFacesIdx.size());

    // then get all boundaries
    int NuMTJunctions = 0;
    for (size_t i = 0; i < PatchFacesIdx.size(); i++) {
      std::vector<std::vector<int>> SubMeshFaces;
      std::vector<Geo::Point3<ScalarType>> SubMeshVert;
      std::vector<int> OriginalF;
      std::vector<int> OriginalV;

      // make a copy of the mesh
      CopySubMesh<ScalarType>(Vert, Faces, PatchFacesIdx[i], SubMeshVert,
                              SubMeshFaces, OriginalV, OriginalF);

      // then get the boundaries
      std::vector<std::vector<int>> BorderVIdxAll;
      getBorderSequncesVert(SubMeshFaces, BorderVIdxAll);
      if (BorderVIdxAll.size() != 1) {
        std::cout << "MULTIPLE BORDERS" << std::endl;
        continue;
      }
      std::vector<int> BorderVIdx = BorderVIdxAll[0];

      std::vector<Geo::Point3<ScalarType>> BorderVPos;
      for (size_t j = 0; j < BorderVIdx.size(); j++)
        BorderVPos.push_back(SubMeshVert[BorderVIdx[j]]);

      // sort from the first
      SortFromSeparator(BorderVPos, PatchCornerPos[i]);

      // split in sequences
      std::vector<std::vector<Geo::Point3<ScalarType>>> SideSeq;
      SplitIntoSubSequences(BorderVPos, PatchCornerPos[i], SideSeq);

      if (SideSeq.size() != 4) {
        std::cout << "WARNING Side Size: " << SideSeq.size() << std::endl;
        std::cout << "Corner Size: " << PatchCornerPos[i].size() << std::endl;
      }

      // then split subsequences
      PatchSides[i].Sides.resize(SideSeq.size());

      for (size_t j = 0; j < SideSeq.size(); j++) {
        std::vector<std::vector<Geo::Point3<ScalarType>>> SubSideSeq;
        SplitIntoSubSequences(SideSeq[j], AllCornerPos, SubSideSeq, false);
        PatchSides[i].Sides[j].SubSides = SubSideSeq;
        if (SubSideSeq.size() != 1)
          NuMTJunctions++;
      }
      //   std::cout<<"Seq Size: "<<SideSeq.size()<<std::endl;
    }
    std::cout << "Num TJunctions: " << NuMTJunctions << std::endl;
    InitVariables();

    // VarFeature.resize(NumVariables, false);

    InitPatchAdj();

    // InitVarFeatures(Vert,_Features);
  }

  void GetSubSideExtremes(const int IndexP, const int IndexS,
                          const int IndexSubS, Geo::Point3<ScalarType> &P0,
                          Geo::Point3<ScalarType> &P1) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    assert(IndexS >= 0);
    assert(IndexS < NumSides(IndexP));
    assert(IndexSubS >= 0);
    assert(IndexSubS < NumSubSides(IndexP, IndexS));

    std::vector<Geo::Point3<ScalarType>> SubSideP;
    GetSubSidePoints(IndexP, IndexS, IndexSubS, SubSideP);
    assert(SubSideP.size() > 1);
    P0 = SubSideP[0];
    P1 = SubSideP.back();
  }

  void
  GetSubSideExtremes(std::vector<Geo::Point3<ScalarType>> &SubSExtremes) const {
    SubSExtremes.clear();
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++)
        for (size_t k = 0; k < NumSubSides(i, j); k++) {
          Geo::Point3<ScalarType> P0, P1;
          GetSubSideExtremes(i, j, k, P0, P1);
          SubSExtremes.push_back(P0);
          SubSExtremes.push_back(P1);
        }
  }

  void
  GetSortedSideSamples(const int &IndexPatch, const int &IndexSide,
                       std::vector<Geo::Point3<ScalarType>> &SidePoints) const {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchSides.size());
    assert(IndexSide >= 0);
    assert(IndexSide < NumSides(IndexPatch));
    SidePoints.clear();
    // IsFeature.clear();
    for (size_t i = 0; i < NumSubSides(IndexPatch, IndexSide); i++) {

      int IndexVar = PatchSides[IndexPatch].Sides[IndexSide].Variables[i];
      std::vector<Geo::Point3<ScalarType>> SeqP =
          PatchSides[IndexPatch].Sides[IndexSide].SubSides[i];
      assert(SeqP.size() >= 2);

      Geo::Point3<ScalarType> Pos0 = SeqP[0];
      Geo::Point3<ScalarType> Pos1 = SeqP.back();
      std::vector<Geo::Point3<ScalarType>> side_samples =
          SampleSidePoints[IndexVar];

      if (Pos0 == side_samples.back()) {
        assert(Pos1 == side_samples[0]);
        std::reverse(side_samples.begin(), side_samples.end());
      } else {
        assert(Pos0 == side_samples[0]);
        assert(Pos1 == side_samples.back());
      }

      SidePoints.insert(SidePoints.end(), side_samples.begin(),
                        side_samples.end());
    }
    // erase duplicates along TJunctions
    SidePoints.erase(std::unique(SidePoints.begin(), SidePoints.end()),
                     SidePoints.end());
  }

  void GetSortedSideSamples(
      const int &IndexPatch,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &SidePoints) const {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchSides.size());
    SidePoints.clear();
    SidePoints.resize(NumSides(IndexPatch));
    // IsFeature.resize(NumSides(IndexPatch));
    for (size_t i = 0; i < NumSides(IndexPatch); i++)
      GetSortedSideSamples(IndexPatch, i, SidePoints[i]);
  }

  void GetSortedSideSamples(
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &SidePoints) const {
    SidePoints.clear();

    if (SampleSidePoints.size() != NumVariables)
      return;

    SidePoints.resize(PatchSides.size());
    for (size_t i = 0; i < PatchSides.size(); i++)
      GetSortedSideSamples(i, SidePoints[i]);
  }

  void GetFeaturesQuadrData(
      const std::vector<Geo::Point3<ScalarType>> &Vert,
      const std::vector<std::pair<int, int>> &Features,
      std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
          &GeoFeatures) {
    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> EdgePos;
    assert(SampleSidePoints.size() == NumVariables);
    GeoFeatures.clear();
    std::set<EdgePos> FeaturesPos;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      Geo::Point3<ScalarType> Pos0 = Vert[IndexV0];
      Geo::Point3<ScalarType> Pos1 = Vert[IndexV1];
      EdgePos Entry(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
      FeaturesPos.insert(Entry);
    }

    for (size_t i = 0; i < NumVariables; i++) {
      std::vector<std::array<int, 3>> PatchSides;
      SideOfVariable(i, PatchSides);
      assert(PatchSides.size() > 0);
      std::array<int, 3> CurrPSide = PatchSides[0];
      std::vector<Geo::Point3<ScalarType>> SubSidePoints;
      GetSubSidePoints(CurrPSide, SubSidePoints);
      bool IsFeature = false;
      for (size_t j = 0; j < SubSidePoints.size() - 1; j++) {
        Geo::Point3<ScalarType> Pos0 = SubSidePoints[j];
        Geo::Point3<ScalarType> Pos1 = SubSidePoints[j + 1];
        EdgePos Entry(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        IsFeature |= (FeaturesPos.count(Entry) > 0);
        if (IsFeature)
          break;
      }
      if (IsFeature) {
        assert(SampleSidePoints[i].size() >= 2);
        for (size_t j = 0; j < SampleSidePoints[i].size() - 1; j++) {
          Geo::Point3<ScalarType> Pos0 = SampleSidePoints[i][j];
          Geo::Point3<ScalarType> Pos1 = SampleSidePoints[i][j + 1];
          EdgePos Entry(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
          GeoFeatures.push_back(Entry);
        }
      }
    }

    std::sort(GeoFeatures.begin(), GeoFeatures.end());
    GeoFeatures.erase(std::unique(GeoFeatures.begin(), GeoFeatures.end()),
                      GeoFeatures.end());
  }

  void GetQuadrangulationData(
      std::vector<std::vector<int>> &PatchFaces,
      std::vector<std::vector<int>> &PatchCorners,
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &SideSamples) const {
    PatchFaces = PatchFacesIdx;
    PatchCorners = PatchCornersIdx;
    GetSortedSideSamples(SideSamples);
    //_GeoFeatures = GeoFeatures;
  }

  void Solve(ScalarType _TargetL, const solveMode &_SM,
             bool _force_parity = true, ScalarType _RegVsIso = 0.9) {

    RegVsIso = _RegVsIso;
    force_parity = _force_parity;

    switch (_SM) {
    case SM_Base:
      InitTargetValue(_TargetL);
      break;
    case SM_Base_Regular:
      InitTargetValue(_TargetL);
      MatchIntegerValues();
      break;
    case SM_Iterative:
      InitTargetValue(_TargetL);
      MatchIntegerValues();
      EM = EM_ForceRegular;
      if (RegVsIso != 1)
        EM = EM_Blend;
      SolveIntegerValues();
      break;
    // case SM_Ite_Blend:
    //   InitTargetValue(_TargetL);
    //   MatchIntegerValues();
    //   EM = EM_Blend;
    //   SolveIntegerValues();
    //   break;
    default:
      InitTargetValue(_TargetL);
      InitOppositeEquivalences();
      SolveGlobal();
      assert(_SM == SM_GlobalSolve);
      // code block
    }
    InitMatchingType();
    FindDiscretizationPoint();
  }

  SideDiscretization() {
    NumVariables = 0;
    UseGeoAVG = false;
    force_parity = true;
    current_doubleVal = false;
    check_parity = true;
    EM = EM_ForceRegular;
  }
};

} // namespace Geo
#endif
