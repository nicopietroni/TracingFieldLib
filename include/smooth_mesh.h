#ifndef SMOOTH_MESH
#define SMOOTH_MESH

#include "Space_Query/triangle_grid_3D.h"
#include "box3.h"
#include <cassert>
#include <cmath>
#include <connectivity.h>
#include <cstddef>
#include <iostream>
#include <set>
#include <solve/least_squares_lin_constraint_solver.h>
#include <stack>
#include <utility>
#include <vector>
// #include "mesh_mapping.h"
#include "Space_Query/segment_grid_3D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "point3.h"

// #include "polyhedra_connectivity.h"
// #include "header"

namespace Geo {

template <class ScalarType>
void SmoothLaplacianVertNormals(
    std::vector<Geo::Point3<ScalarType>> &VertNorm,
    const std::vector<std::vector<int>> &Faces, int NumIte,
    const std::vector<int> &FixedVert = std::vector<int>(),
    ScalarType Damp = 0.5) {
  std::set<int> FixedVertSet(FixedVert.begin(), FixedVert.end());
  for (int s = 0; s < NumIte; s++) {
    std::vector<Geo::Point3<ScalarType>> NewNorm(
        VertNorm.size(), Geo::Point3<ScalarType>(0, 0, 0));
    std::vector<int> SumNum(VertNorm.size(), 0);

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        Geo::Point3<ScalarType> N0 = VertNorm[IndexV0];
        Geo::Point3<ScalarType> N1 = VertNorm[IndexV1];

        NewNorm[IndexV0] += N1;
        SumNum[IndexV0]++;

        NewNorm[IndexV1] += N0;
        SumNum[IndexV1]++;
      }
    for (size_t i = 0; i < NewNorm.size(); i++) {
      if (FixedVertSet.count(i) > 0)
        continue;
      NewNorm[i] /= (ScalarType)SumNum[i];
      NewNorm[i].Normalize();
      VertNorm[i] = NewNorm[i] * (1 - Damp) + VertNorm[i] * Damp;
      VertNorm[i].Normalize();
    }
  }
}

// Bilateral filter function for mesh signals
template <class ScalarType>
void SmoothLBilateralVertSignal(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    std::vector<ScalarType> &VertSignal, ScalarType sigma_s = 1.0,
    ScalarType sigma_r = 0.5, int iterations = 3) {
  int num_vertices = VertPos.size();
  // std::vector<double> filtered_signal = signal;

  std::vector<std::vector<int>> VVAdj;
  FindVertexStar(Faces, VertPos.size(), VVAdj);
  // Loop over the number of iterations
  for (int iter = 0; iter < iterations; ++iter) {
    std::vector<ScalarType> new_signal(num_vertices, 0);

    // Iterate over all vertices
    for (int i = 0; i < num_vertices; ++i) {
      Geo::Point3<ScalarType> v_i = VertPos[i];
      ScalarType f_i = VertSignal[i];

      ScalarType numerator = 0.0;
      ScalarType denominator = 0.0;

      // Iterate over the neighboring vertices
      for (size_t j = 0; j < VVAdj[i].size(); j++) {
        int indexN = VVAdj[i][j];
        Geo::Point3<ScalarType> v_j = VertPos[indexN];
        ScalarType f_j = VertSignal[indexN];

        // Compute spatial weight (distance between vertices)
        ScalarType spatial_weight =
            std::exp(std::sqrt((v_i - v_j).Norm()) / (2 * sigma_s * sigma_s));

        // Compute range weight (difference in signal values)
        ScalarType range_weight =
            std::exp(-(f_i - f_j) * (f_i - f_j) / (2 * sigma_r * sigma_r));

        // Total weight
        ScalarType weight = spatial_weight * range_weight;

        // Accumulate numerator and denominator for the weighted average
        numerator += weight * f_j;
        denominator += weight;
      }

      // Update the new signal value at vertex i
      if (denominator != 0.0) {
        new_signal[i] = numerator / denominator;
      } else {
        new_signal[i] =
            f_i; // In case of zero denominator, retain the original signal
      }
    }

    // Update the filtered signal with the new values
    VertSignal = new_signal;
  }
}

// Bilateral filter function for mesh signals
template <class ScalarType>
void SmoothLaplacianlSignal(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            std::vector<ScalarType> &VertSignal,
                            int iterations = 3) {
  for (int s = 0; s < iterations; s++) {
    std::vector<ScalarType> VertSignalNew = VertSignal;
    std::vector<int> VertSignalNum;
    VertSignalNum.resize(VertPos.size(), 1);
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        assert(IndexV0 != IndexV1);
        VertSignalNum[IndexV0]++;
        VertSignalNum[IndexV1]++;
        VertSignalNew[IndexV0] += VertSignal[IndexV1];
        VertSignalNew[IndexV1] += VertSignal[IndexV0];
      }
    for (size_t i = 0; i < VertSignalNew.size(); i++)
      VertSignal[i] = VertSignalNew[i] / (ScalarType)VertSignalNum[i];
  }
}

template <class ScalarType> struct SmoothParam {
  int NumIte;
  std::vector<int> FixedVert;
  std::vector<std::pair<int, int>> Features;
  bool line_border;
  ScalarType Damp;
  std::vector<Geo::Point3<ScalarType>> *ReprojVert;
  std::vector<std::vector<int>> *ReprojFaces;

  SmoothParam() {
    NumIte = 10;
    line_border = true;
    Damp = 0.5;
    ReprojVert = NULL;
    ReprojFaces = NULL;
  }
};

// template <class ScalarType>
// void SmoothLaplacian(std::vector<Geo::Point3<ScalarType>> &VertPos,
//                      const std::vector<std::vector<int>> &Faces, int NumIte,
//                      const std::vector<int> &FixedVert = std::vector<int>(),
//                      const std::vector<std::pair<int, int>> &Features =
//                          std::vector<std::pair<int, int>>(),
//                      bool line_border = true, ScalarType Damp = 0.5) {

template <class ScalarType>
void SmoothLaplacian(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const SmoothParam<ScalarType> &SParam = SmoothParam<ScalarType>()) {
  
  //std::cout<<"Features NUM:"<<SParam.Features.size()<<std::endl;
  //exit(0);
  // save the features
  bool has_features = (SParam.Features.size() > 0);
  Geo::SegmentSpatialIndex<ScalarType> SegGrid;
  std::vector<Geo::Segment3<ScalarType>> Segs;
  if (has_features)
    Geo::InitFromMeshFeatures<ScalarType>(VertPos, SParam.Features, Segs,
                                          SegGrid);

  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  bool reproject = false;
  if ((SParam.ReprojVert != NULL) && (SParam.ReprojFaces != NULL)) {
    TriGrid.Init((*SParam.ReprojVert), (*SParam.ReprojFaces));
    reproject = true;
  }

  std::set<int> FixedVertSet(SParam.FixedVert.begin(), SParam.FixedVert.end());
  std::vector<bool> FeatureVert(VertPos.size(), false);
  for (size_t i = 0; i < SParam.Features.size(); i++) {
    FeatureVert[SParam.Features[i].first] = true;
    FeatureVert[SParam.Features[i].second] = true;
  }

  std::vector<bool> IsB;
  Geo::getBorderVert(Faces, IsB);
  for (int s = 0; s < SParam.NumIte; s++) {
    std::vector<Geo::Point3<ScalarType>> NewVertPos(
        VertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));
    std::vector<int> SumNum(VertPos.size(), 0);

    // cumulate features
    for (size_t i = 0; i < SParam.Features.size(); i++) {
      int IndexV0 = SParam.Features[i].first;
      int IndexV1 = SParam.Features[i].second;
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      NewVertPos[IndexV1] += P0;
      SumNum[IndexV1]++;
      NewVertPos[IndexV0] += P1;
      SumNum[IndexV0]++;
    }

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];

        if (SParam.line_border) {
          if ((!IsB[IndexV0]) && (!FeatureVert[IndexV0])) {
            NewVertPos[IndexV0] += P1;
            SumNum[IndexV0]++;
          }
          if ((!IsB[IndexV1]) && (!FeatureVert[IndexV1])) {
            NewVertPos[IndexV1] += P0;
            SumNum[IndexV1]++;
          }
        } else {
          if (!FeatureVert[IndexV0]) {
            NewVertPos[IndexV0] += P1;
            SumNum[IndexV0]++;
          }
          if (!FeatureVert[IndexV1]) {
            NewVertPos[IndexV1] += P0;
            SumNum[IndexV1]++;
          }
        }
      }
    if (SParam.line_border) {
      std::vector<std::pair<int, int>> BorderE;
      Geo::getBorderEdges(Faces, BorderE);
      for (size_t i = 0; i < BorderE.size(); i++) {
        int IndexV0 = BorderE[i].first;
        int IndexV1 = BorderE[i].second;
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];

        if (!FeatureVert[IndexV0]) {
          NewVertPos[IndexV0] += P1;
          SumNum[IndexV0]++;
        }
        if (!FeatureVert[IndexV1]) {
          NewVertPos[IndexV1] += P0;
          SumNum[IndexV1]++;
        }
      }
    }
    for (size_t i = 0; i < NewVertPos.size(); i++) {
      if (FixedVertSet.count(i) > 0)
        continue;
      if (SumNum[i] < 2)
        continue;
      NewVertPos[i] /= (ScalarType)SumNum[i];
      VertPos[i] = NewVertPos[i] * (1 - SParam.Damp) + VertPos[i] * SParam.Damp;
    }
    if (reproject) {
      for (size_t i = 0; i < VertPos.size(); i++) {
        int closestF;
        Geo::Point3<ScalarType> closestPt;
        if (!FeatureVert[i]) {
          TriGrid.GridClosest((*SParam.ReprojVert), (*SParam.ReprojFaces),
                              VertPos[i], TriGrid.bbox.Diag(), closestF,
                              closestPt);
        } else {
          int closestS;
          SegGrid.GridClosest(Segs, VertPos[i], SegGrid.bbox.Diag(), closestS,
                              closestPt);
        }
        VertPos[i] = closestPt;
      }
    }
  }
}

template <class ScalarType>
void SmoothLaplacian(
    std::vector<Geo::Point2<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const SmoothParam<ScalarType> &SParam = SmoothParam<ScalarType>()) {
      std::vector<Geo::Point3<ScalarType>> VertPos3D;
  for (size_t i=0;i<VertPos.size();i++)
    VertPos3D.push_back(Geo::Point3<ScalarType>(VertPos[i].X,VertPos[i].Y,0));
  
  SmoothLaplacian(VertPos3D,Faces,SParam);

  for (size_t i=0;i<VertPos3D.size();i++)
  {
    VertPos[i].X=VertPos3D[i].X;
    VertPos[i].Y=VertPos3D[i].Y;
  }
}

template <class ScalarType>
void SmoothLaplacianStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::set<int> &FixedV,
                         const std::vector<std::pair<int, int>> &LineSmooth =
                             std::vector<std::pair<int, int>>(),
                         ScalarType Damp = 0.5) {
  std::vector<Geo::Point3<ScalarType>> AvgPos(VertPos.size(),
                                              Geo::Point3<ScalarType>(0, 0, 0));
  std::vector<int> AvgNum(VertPos.size(), 0);

  std::set<int> LineSmoothV;
  for (size_t i = 0; i < LineSmooth.size(); i++) {
    LineSmoothV.insert(LineSmooth[i].first);
    LineSmoothV.insert(LineSmooth[i].second);
  }

  // cumulate the non line
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
      bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
      bool IsLineSmoothV0 = (LineSmoothV.count(IndexV0) > 0);
      bool IsLineSmoothV1 = (LineSmoothV.count(IndexV1) > 0);
      if ((!IsFixedV0) && (!IsLineSmoothV0)) {
        AvgPos[IndexV0] += VertPos[IndexV1];
        AvgNum[IndexV0]++;
      }
      if ((!IsFixedV1) && (!IsLineSmoothV1)) {
        AvgPos[IndexV1] += VertPos[IndexV0];
        AvgNum[IndexV1]++;
      }
    }

  // cumulate the line
  for (size_t i = 0; i < LineSmooth.size(); i++) {
    int IndexV0 = LineSmooth[i].first;
    int IndexV1 = LineSmooth[i].second;
    bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
    bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
    if (!IsFixedV0) {
      AvgPos[IndexV0] += VertPos[IndexV1];
      AvgNum[IndexV0]++;
    }
    if (!IsFixedV1) {
      AvgPos[IndexV1] += VertPos[IndexV0];
      AvgNum[IndexV1]++;
    }
  }

  // then average
  for (size_t i = 0; i < AvgNum.size(); i++) {
    if (AvgNum[i] == 0) {
      AvgPos[i] = VertPos[i];
      continue;
    }
    AvgPos[i] /= AvgNum[i];
  }
  for (size_t i = 0; i < AvgPos.size(); i++)
    VertPos[i] = VertPos[i] * Damp + AvgPos[i] * (1 - Damp);
}

template <class ScalarType>
void SmoothLaplacianNonManif(std::vector<Geo::Point3<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             int NumIte, ScalarType Damp = 0.5) {

  std::map<std::pair<int, int>, int> EdgeNum;
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      if (EdgeNum.count(KeyE) == 0)
        EdgeNum[KeyE] = 1;
      else
        EdgeNum[KeyE]++;
    }

  std::vector<std::pair<int, int>> NotManifE;
  std::vector<std::pair<int, int>> BorderE;
  std::map<std::pair<int, int>, int>::iterator IteEdgeMap;
  std::vector<int> BorderV;
  for (IteEdgeMap = EdgeNum.begin(); IteEdgeMap != EdgeNum.end();
       IteEdgeMap++) {
    int NuME = (*IteEdgeMap).second;
    std::pair<int, int> EdgeK = (*IteEdgeMap).first;
    if (NuME == 1) {
      BorderE.push_back(EdgeK);
      BorderV.push_back(EdgeK.first);
      BorderV.push_back(EdgeK.second);
    }
    if (NuME > 2) {
      NotManifE.push_back(EdgeK);
    }
  }

  std::sort(BorderE.begin(), BorderE.end());
  BorderE.erase(std::unique(BorderE.begin(), BorderE.end()), BorderE.end());

  std::sort(NotManifE.begin(), NotManifE.end());
  NotManifE.erase(std::unique(NotManifE.begin(), NotManifE.end()),
                  NotManifE.end());

  std::sort(BorderV.begin(), BorderV.end());
  BorderV.erase(std::unique(BorderV.begin(), BorderV.end()), BorderV.end());

  // find the fixed vertx
  std::set<int> FixedV;

  // then get the vertices tobe fixed
  std::map<int, int> VertNonManifE;
  for (size_t i = 0; i < NotManifE.size(); i++) {
    int IndexV0 = NotManifE[i].first;
    int IndexV1 = NotManifE[i].second;
    VertNonManifE[IndexV0]++;
    VertNonManifE[IndexV1]++;
    // if (BorderV.count(IndexV0))
    //   FixedV.insert(IndexV0);
    // if (BorderV.count(IndexV1))
    //   FixedV.insert(IndexV1);
  }

  // std::map<int, int>::iterator IteMapV;
  // for (IteMapV = VertNonManifE.begin(); IteMapV != VertNonManifE.end();
  //      IteMapV++)
  //   if ((*IteMapV).second > 2)
  //     FixedV.insert((*IteMapV).first);

  // then smooth
  std::vector<std::pair<int, int>> LineSmooth = NotManifE;
  LineSmooth.insert(LineSmooth.end(), BorderE.begin(), BorderE.end());
  for (size_t s = 0; s < NumIte; s++) {
    SmoothLaplacianStep(VertPos, Faces, FixedV, LineSmooth, Damp);

    // //reproject
    // for (size_t i=0;i<BorderV.size();i++)
    // {
    //   int IndexV=BorderV[i];
    //   assert(IndexV>=0);
    //   assert(IndexV<VertPos.size());
    //   Geo::Point3<ScalarType> Pos=VertPos[IndexV];
    //   int ClosestF;
    //   Geo::Point3<ScalarType> ClosestPos;
    //   TriGrid.GridClosest(ReprojectV,ReprojectFaces,Pos,
    //   TriGrid.bbox.Diag(),ClosestF,ClosestPos); VertPos[IndexV]=ClosestPos;
    // }
  }
}

template <class ScalarType>
void SmoothLaplacianImplicit(std::vector<Geo::Point3<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &FixedVert) {
  ScalarType WeightFixed = 1000000;
  ScalarType WeightLap = 0.00001;

  std::vector<Geo::LinearConstraint<ScalarType>> LinConstr;
  std::set<int> FixedVertSet(FixedVert.begin(), FixedVert.end());
  // std::cout<<"there are "<<FixedVertSet.size()<<" fixed"<<std::endl;
  // add fixed velue
  for (size_t i = 0; i < VertPos.size(); i++) {
    if (FixedVertSet.count((int)i) > 0) {
      for (size_t c = 0; c < 3; c++) {
        Geo::LinearConstraint<ScalarType> Constr;
        Constr.RHSide = VertPos[i].V(c);
        Constr.Weight = WeightFixed;
        Geo::VarTerm<ScalarType> Var(i, c, 1);
        Constr.LHSide.push_back(Var);
        LinConstr.push_back(Constr);
      }
    }
  }

  // add laplacian value
  std::vector<std::vector<int>> KernelV;
  Geo::FindVertexStar(Faces, VertPos.size(), KernelV);

  //    std::vector<bool> IsB;
  //    getBorderVert(Faces,IsB);

  for (size_t i = 0; i < KernelV.size(); i++) {
    int IndexV0 = i;
    //        bool IsBV0=IsB[IndexV0];
    //        if (IsBV0)continue;
    for (size_t c = 0; c < 3; c++) {
      Geo::LinearConstraint<ScalarType> Constr;
      Constr.Weight = WeightLap;
      Constr.RHSide = 0;
      Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV0, c, -1));
      ScalarType OtherW = (ScalarType)1 / (ScalarType)KernelV[i].size();

      for (size_t j = 0; j < KernelV[i].size(); j++) {
        int IndexV1 = KernelV[i][j];
        assert(IndexV0 != IndexV1);

        //                bool IsBV1=IsB[IndexV1];
        //                if (IsBV0 && (!IsBV1))continue;
        //                if (IsBV0 && IsBV1)
        //                    Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1,c,0.5));
        //                else
        Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1, c, OtherW));
      }
      LinConstr.push_back(Constr);
    }
  }

  Geo::LinearConstraintLeastSQSolver<ScalarType> LinSqSolver;
  LinSqSolver.Init(LinConstr);
  std::vector<ScalarType> Solution;
  LinSqSolver.Solve(Solution);

  for (size_t i = 0; i < VertPos.size(); i++) {
    int IndexX = (i * 3);
    int IndexY = IndexX + 1;
    int IndexZ = IndexX + 2;
    ScalarType ValX = Solution[IndexX];
    ScalarType ValY = Solution[IndexY];
    ScalarType ValZ = Solution[IndexZ];
    VertPos[i] = Geo::Point3<ScalarType>(ValX, ValY, ValZ);
  }
}

template <class ScalarType>
void SmoothLaplacianImplicit(std::vector<Geo::Point2<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &FixedVert)
{
  std::vector<Geo::Point3<ScalarType>> VertPos3D;
  for (size_t i=0;i<VertPos.size();i++)
    VertPos3D.push_back(Geo::Point3<ScalarType>(VertPos[i].X,VertPos[i].Y,0));
  
  SmoothLaplacianImplicit(VertPos3D,Faces,FixedVert);

  for (size_t i=0;i<VertPos3D.size();i++)
  {
    VertPos[i].X=VertPos3D[i].X;
    VertPos[i].Y=VertPos3D[i].Y;
  }
}

template <class ScalarType>
static void
SmoothDegerateTriangles(std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Connectivity,
                        const ScalarType &angleDelta) {
  std::vector<int> DegenerateF;
  GetDegenerateTriangles(VertPos, Connectivity, angleDelta, DegenerateF);
  std::vector<bool> IsFix = std::vector<bool>(VertPos.size(), true);
  for (size_t i = 0; i < DegenerateF.size(); i++) {
    int IndexF = DegenerateF[i];
    for (size_t j = 0; j < 3; j++) {
      int IndexV = Connectivity[IndexF][j];
      IsFix[IndexV] = false;
    }
  }
  std::vector<int> FixedVert;
  int To_Smooth = 0;
  for (size_t i = 0; i < IsFix.size(); i++) {
    if (!IsFix[i]) {
      To_Smooth++;
      continue;
    }
    FixedVert.push_back(i);
  }
  if (To_Smooth == 0)
    return;
  std::cout << "**have to smooth:" << To_Smooth << std::endl;
  const std::vector<std::pair<int, int>> Features;
  Geo::SmoothLaplacian<ScalarType>(VertPos, Connectivity, 1, FixedVert,
                                   Features, false);
}

template <class ScalarType>
void PosInstances(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                  std::vector<int> &Istances) {
  Istances = std::vector<int>(Pos3D.size(), 0);
  std::map<Geo::Point3<ScalarType>, int> PosCount;
  for (size_t i = 0; i < Pos3D.size(); i++)
    PosCount[Pos3D[i]]++;
  for (size_t i = 0; i < Pos3D.size(); i++) {
    assert(PosCount.count(Pos3D[i]) > 0);
    Istances[i] = PosCount[Pos3D[i]];
  }
}

template <class ScalarType>
int DesambiguateSamePosVertices(
    std::vector<Geo::Point3<ScalarType>> &Pos3D,
    const std::vector<std::vector<int>> &Faces, ScalarType delta,
    std::vector<int> *ModVert = NULL,
    std::vector<Geo::Point3<ScalarType>> *OrigPos = NULL) {
  if (ModVert != NULL) {
    assert(OrigPos != NULL);
    (*ModVert).clear();
    (*OrigPos).clear();
  }

  // first find vertices with double positions
  std::map<Geo::Point3<ScalarType>, std::vector<int>> PosToIdx;
  for (size_t i = 0; i < Pos3D.size(); i++)
    PosToIdx[Pos3D[i]].push_back(i);

  std::vector<int> FixedV;
  for (size_t i = 0; i < Pos3D.size(); i++) {
    if (PosToIdx[Pos3D[i]].size() == 1)
      FixedV.push_back(i);
  }
  // then compute laplacian smooth
  std::vector<Geo::Point3<ScalarType>> OldPos = Pos3D;
  // std::vector<std::pair<int, int>> Features;
  Geo::SmoothParam<ScalarType> SParam;
  SParam.NumIte = 5;
  SParam.line_border = false;
  // Geo::SmoothLaplacian<ScalarType>(Pos3D, Faces, 5, FixedV, Features, false);
  Geo::SmoothLaplacian<ScalarType>(Pos3D, Faces, SParam);
  int Moved = 0;
  for (size_t i = 0; i < OldPos.size(); i++) {
    assert(PosToIdx.count(OldPos[i]) > 0);
    // if (PosToIdx[OldPos[i]].size() == 2) {
    if (PosToIdx[OldPos[i]].size() > 1) {
      if (ModVert != NULL) {
        (*ModVert).push_back(i);
        (*OrigPos).push_back(OldPos[i]);
        assert(OrigPos != NULL);
      }
      Geo::Point3<ScalarType> DirMov = Pos3D[i] - OldPos[i];
      DirMov.Normalize();
      DirMov = DirMov * delta;
      Pos3D[i] = OldPos[i] + DirMov;
      Moved++;
    } else
      Pos3D[i] = OldPos[i];
  }
  return Moved;
}

// Coherently orient the triangle mesh
void OrientCoherently(std::vector<std::vector<int>> &faces) {
  std::map<std::pair<int, int>, int> edgeToFace;
  std::vector<bool> visited(faces.size(), false);

  // Populate edge-to-face map
  for (int i = 0; i < faces.size(); ++i) {
    const auto &face = faces[i];
    for (int j = 0; j < faces[i].size(); ++j) {
      int v1 = face[j];
      int v2 = face[(j + 1) % faces[i].size()];
      edgeToFace[{v1, v2}] = i;
    }
  }

  // BFS to propagate orientation
  std::queue<int> toVisit;
  toVisit.push(0);
  visited[0] = true;

  while (!toVisit.empty()) {
    int currentFace = toVisit.front();
    toVisit.pop();

    const auto &face = faces[currentFace];
    for (int j = 0; j < faces[currentFace].size(); ++j) {
      int v1 = face[j];
      int v2 = face[(j + 1) % faces[currentFace].size()];

      auto it = edgeToFace.find({v2, v1});
      if (it != edgeToFace.end()) {
        int neighborFace = it->second;

        if (!visited[neighborFace]) {
          visited[neighborFace] = true;

          // Ensure consistent orientation
          auto &neighbor = faces[neighborFace];
          if (neighbor[0] == v2 && neighbor[1] == v1) {
            std::swap(neighbor[1], neighbor[2]);
          } else if (neighbor[1] == v2 && neighbor[2] == v1) {
            std::swap(neighbor[0], neighbor[2]);
          } else if (neighbor[2] == v2 && neighbor[0] == v1) {
            std::swap(neighbor[0], neighbor[1]);
          }

          toVisit.push(neighborFace);
        }
      }
    }
  }
}

} // namespace Geo
#endif
