#ifndef SMOOTH_MESH_SAMPLING_DIR
#define SMOOTH_MESH_SAMPLING_DIR

#include "connectivity.h"
#include "fitting_points.h"
#include "point3.h"
#include "quad3.h"
#include "temp_mark.h"
#include "triangle3.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sys/_types/_size_t.h>
#include <vector>

namespace Geo {

template <class ScalarType> class QuadSmoothSampling {

  std::vector<Geo::Point3<ScalarType>> &QuadVertPos;
  const std::vector<std::vector<int>> &QuadFaces;
  const std::vector<std::pair<int, int>> &Quad_Features;
  std::vector<std::vector<int>> VF, VFIndex;
  int SampleDensity;

  std::vector<int> VertTMark;
  int GlobalMark;

  struct OptimizationOp {
    int IndexV;
    int TMark;
    Geo::Point3<ScalarType> TargetP;
    ScalarType Priority;

    inline bool operator<(const OptimizationOp &oth) const {
      return (Priority > oth.Priority);
    }
  };

  std::vector<OptimizationOp> OpHeap;

  ScalarType AspectRatio(const int IndexF) {
    bool HasConcave = HasConcaveV<ScalarType>(QuadVertPos, QuadFaces, IndexF);
    ScalarType Q = ScaledJacobianQuad(QuadVertPos, QuadFaces, IndexF);
    if (HasConcave)
      Q -= 1;
    return Q;
  }

  ScalarType WorstAspect(const int IndexV) {
    assert(VF[IndexV].size() > 0);
    ScalarType WorstA = AspectRatio(VF[IndexV][0]);
    for (size_t i = 1; i < VF[IndexV].size(); i++)
      WorstA = std::min(WorstA, AspectRatio(VF[IndexV][i]));

    return WorstA;
  }

  ScalarType BestPosition(const int IndexV, Geo::Point3<ScalarType> &BestPos) {
    std::vector<Geo::Point3<ScalarType>> SamplePos;
    assert(IndexV>=0);
    assert(IndexV<VF.size());
    for (size_t i = 0; i < VF[IndexV].size(); i++) {
      int IndexF = VF[IndexV][i];
      std::vector<Geo::Point3<ScalarType>> CurrSamples;
      SampleBilinear(QuadVertPos, QuadFaces, IndexF, SampleDensity, CurrSamples);
                    
      SamplePos.insert(SamplePos.end(), CurrSamples.begin(), CurrSamples.end());
    }
    // then get current aspect ratio
    ScalarType CurrRatio = WorstAspect(IndexV);
    Geo::Point3<ScalarType> OldPos = QuadVertPos[IndexV];
    BestPos = OldPos;
    ScalarType BestRatio = CurrRatio;
    for (size_t i = 0; i < SamplePos.size(); i++) {
      QuadVertPos[IndexV] = SamplePos[i];
      ScalarType TestRatio = WorstAspect(IndexV);
      if (TestRatio <= BestRatio)
        continue;
      BestRatio = TestRatio;
      BestPos = SamplePos[i];
    }
    QuadVertPos[IndexV] = OldPos;
    return (BestRatio - CurrRatio);
  }

  std::vector<int> InvolvedVert(const OptimizationOp &Op)const {
    std::vector<int> VertKern;
    int IndexV = Op.IndexV;
    for (size_t i = 0; i < VF[IndexV].size(); i++) {
      int CurrF = VF[IndexV][i];
      for (size_t j = 0; j < QuadFaces[CurrF].size(); j++)
        VertKern.push_back(QuadFaces[CurrF][j]);
    }
    std::sort(VertKern.begin(), VertKern.end());
    VertKern.erase(std::unique(VertKern.begin(), VertKern.end()),
                   VertKern.end());
    return VertKern;
  }

  bool IsUpToDate(const OptimizationOp &Op) const {
    std::vector<int> InvolvV = InvolvedVert(Op);
    for (size_t i = 0; i < InvolvV.size(); i++) {
      int IndexV = InvolvV[i];
      assert(IndexV >= 0);
      assert(IndexV < VertTMark.size());
      if (VertTMark[IndexV] > Op.TMark)
        return false;
    }
    return true;
  }

  void Apply(const OptimizationOp &Op) {
    QuadVertPos[Op.IndexV] = Op.TargetP;
    std::vector<int> InvolvV = InvolvedVert(Op);
    GlobalMark++;
    for (size_t i = 0; i < InvolvV.size(); i++)
      VertTMark[InvolvV[i]] = GlobalMark;
  }

  void UpdatePriority(OptimizationOp &Op) {
    Geo::Point3<ScalarType> BestPos;
    ScalarType Prio = BestPosition(Op.IndexV, BestPos);
    Op.Priority = Prio;
    Op.TargetP = BestPos;
  }

  std::vector<OptimizationOp> GetNewOp(const OptimizationOp &Op) {
    std::vector<int> InvolvV = InvolvedVert(Op);
    std::vector<OptimizationOp> NewOps;
    for (size_t i = 0; i < InvolvV.size(); i++) {
      int CurrIndexV = InvolvV[i];

      OptimizationOp Op;
      Op.IndexV = CurrIndexV;
      Op.TMark = GlobalMark;
      UpdatePriority(Op);
      NewOps.push_back(Op);
    }
    return NewOps;
  }

  void InitHeap() {
    VertTMark = std::vector<int>(QuadVertPos.size(), 0); 
    GlobalMark = 0;
    OpHeap.clear();
    for (size_t i = 0; i < QuadVertPos.size(); i++) {
     
      OptimizationOp Op;
      Op.IndexV = i;
      UpdatePriority(Op);
      if (Op.Priority<=0)continue;
      OpHeap.push_back(Op);
    }
    std::sort(OpHeap.begin(), OpHeap.end());
    std::cout<<"Initialized "<<OpHeap.size()<<" Operations"<<std::endl;
  }

public:

  void Optimize() {
    getVFAdjacency(QuadFaces,VF,VFIndex);
    InitHeap();
    do {
      std::pop_heap(OpHeap.begin(), OpHeap.end());
      OptimizationOp Op = OpHeap.back();

      OpHeap.pop_back();

      bool isUptoDate = IsUpToDate(Op);
      if (!isUptoDate)
        continue;

      Apply(Op);

      std::vector<OptimizationOp> NewOp = GetNewOp(Op);
      
      for (size_t i=0;i<NewOp.size();i++)
      {
        if (NewOp[i].Priority<=0)continue;
        OpHeap.push_back(NewOp[i]);
        std::push_heap(OpHeap.begin(), OpHeap.end());
      }
    } while (OpHeap.size() > 0);
  }
  
  QuadSmoothSampling(std::vector<Geo::Point3<ScalarType>> &_QuadVertPos,
                     const std::vector<std::vector<int>> &_QuadFaces,
                     const std::vector<std::pair<int, int>> &_Quad_Features)
      : QuadVertPos(_QuadVertPos), QuadFaces(_QuadFaces),
        Quad_Features(_Quad_Features) {
    SampleDensity = 10;
  }
};

} // namespace Geo
#endif
