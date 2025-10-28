#ifndef SMOOTH_MESH
#define SMOOTH_MESH

#include "IO/io_obj.h"
#include "Space_Query/triangle_grid_3D.h"
#include "box3.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <connectivity.h>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <set>
#include <solve/least_squares_lin_constraint_solver.h>
#include <stack>
#include <sys/_types/_size_t.h>
#include <utility>
#include <vector>
// #include "mesh_mapping.h"
#include "Space_Query/segment_grid_3D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "best_rotation.h"
#include "edge_mesh_functions.h"
#include "laplacian_neighbors.h"
#include "mesh_features.h"
#include "mesh_functions.h"
#include "normal.h"
#include "point3.h"
#include "polygon3.h"
#include "polyline3.h"
#include "quad3.h"
#include "reproject_mesh.h"
#include "triangle3.h"
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
                            int iterations = 3,
                            std::vector<int> FixedV = std::vector<int>()) {
  std::set<int> FixedVSet(FixedV.begin(), FixedV.end());
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
    for (size_t i = 0; i < VertSignalNew.size(); i++) {
      if (FixedVSet.count(i) > 0)
        continue;
      VertSignal[i] = VertSignalNew[i] / (ScalarType)VertSignalNum[i];
    }
  }
}

// template <class ScalarType>
// ScalarType BestAlignedPoly(const std::vector<Geo::Point3<ScalarType>>
// &FacePos,
//                            std::vector<Geo::Point3<ScalarType>> &bestPoly) {
//   // find the perimeter
//   Geo::Polyline3<ScalarType> Poly;
//   Poly.PolyPos = FacePos;
//   ScalarType Len = Poly.Lenght();

//   ScalarType sideLength = Len / FacePos.size(); // Compute side length
//   ScalarType radius =
//       sideLength / (2 * sin(M_PI / FacePos.size())); // Compute circumradius

//   bestPoly.clear();
//   for (int i = 0; i < FacePos.size(); ++i) {
//     ScalarType angle = (2.0 * M_PI * i) / FacePos.size(); // Angle in radians
//     ScalarType x = radius * cos(angle);
//     ScalarType y = radius * sin(angle);
//     bestPoly.push_back(Geo::Point3<ScalarType>(x, y, 0));
//   }

//   Geo::PointMatch<ScalarType>::TrasformRigidMatch(FacePos, bestPoly);

//   ScalarType Err = 0;
//   for (size_t i = 0; i < FacePos.size(); i++)
//     Err += pow((FacePos[i] - bestPoly[i]).Norm(), 2);
//   return Err;
// }

// template <class ScalarType>
// Geo::Point3<ScalarType>
// NormalPoly(const std::vector<Geo::Point3<ScalarType>> &PolyPos) {

//   // compute normal of polygonal face
//   Point3<ScalarType> Norm(0, 0, 0);
//   int sizeF = PolyPos.size();
//   for (int i = 0; i < sizeF; i++) {
//     int IndexV0 = (i + sizeF - 1) % sizeF;
//     int IndexV1 = i;
//     int IndexV2 = (i + 1) % sizeF;
//     Point3<ScalarType> P0 = PolyPos[IndexV0];
//     Point3<ScalarType> P1 = PolyPos[IndexV1];
//     Point3<ScalarType> P2 = PolyPos[IndexV2];
//     Point3<ScalarType> Dir0 = P1 - P0;
//     Point3<ScalarType> Dir1 = P2 - P0;
//     Dir0.Normalize();
//     Dir1.Normalize();
//     Norm += Cross(Dir0, Dir1);
//   }
//   Norm.Normalize();
//   return Norm;
// }

// template <class ScalarType>
// void SmoothLaplacian(std::vector<Geo::Point3<ScalarType>> &VertPos,
//                      const std::vector<std::vector<int>> &Faces, int NumIte,
//                      const std::vector<int> &FixedVert = std::vector<int>(),
//                      const std::vector<std::pair<int, int>> &Features =
//                          std::vector<std::pair<int, int>>(),
//                      bool line_border = true, ScalarType Damp = 0.5) {

template <class ScalarType>
void GetLaplacianValues(const std::vector<Geo::Point3<ScalarType>> &Values,
                        const std::vector<std::vector<int>> &Neighbors,
                        std::vector<Geo::Point3<ScalarType>> &LaplVal) {

  LaplVal = std::vector<Geo::Point3<ScalarType>>(
      Values.size(), Geo::Point3<ScalarType>(0, 0, 0));

  assert(Neighbors.size() == Values.size());
  for (size_t i = 0; i < Neighbors.size(); i++) {
    for (size_t j = 0; j < Neighbors[i].size(); j++) {
      int IndexN = Neighbors[i][j];
      assert(IndexN >= 0);
      assert(IndexN < Values.size());
      LaplVal[i] += Values[IndexN];
    }
    LaplVal[i] /= Neighbors[i].size();
  }
}

// template <class ScalarType>
// void SmoothLaplacianStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
//                          const std::vector<std::vector<int>> &Faces,
//                          const std::set<int> &FixedV,
//                          const std::vector<std::pair<int, int>> &LineSmooth =
//                              std::vector<std::pair<int, int>>(),
//                          ScalarType Damp = 0.5) {
//   std::vector<Geo::Point3<ScalarType>> AvgPos(VertPos.size(),
//                                               Geo::Point3<ScalarType>(0, 0,
//                                               0));
//   std::vector<int> AvgNum(VertPos.size(), 0);

//   std::set<int> LineSmoothV;
//   for (size_t i = 0; i < LineSmooth.size(); i++) {
//     LineSmoothV.insert(LineSmooth[i].first);
//     LineSmoothV.insert(LineSmooth[i].second);
//   }

//   // cumulate the non line
//   for (size_t i = 0; i < Faces.size(); i++)
//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int IndexV0 = Faces[i][j];
//       int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
//       bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
//       bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
//       bool IsLineSmoothV0 = (LineSmoothV.count(IndexV0) > 0);
//       bool IsLineSmoothV1 = (LineSmoothV.count(IndexV1) > 0);
//       if ((!IsFixedV0) && (!IsLineSmoothV0)) {
//         AvgPos[IndexV0] += VertPos[IndexV1];
//         AvgNum[IndexV0]++;
//       }
//       if ((!IsFixedV1) && (!IsLineSmoothV1)) {
//         AvgPos[IndexV1] += VertPos[IndexV0];
//         AvgNum[IndexV1]++;
//       }
//     }

//   // cumulate the line
//   for (size_t i = 0; i < LineSmooth.size(); i++) {
//     int IndexV0 = LineSmooth[i].first;
//     int IndexV1 = LineSmooth[i].second;
//     bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
//     bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
//     if (!IsFixedV0) {
//       AvgPos[IndexV0] += VertPos[IndexV1];
//       AvgNum[IndexV0]++;
//     }
//     if (!IsFixedV1) {
//       AvgPos[IndexV1] += VertPos[IndexV0];
//       AvgNum[IndexV1]++;
//     }
//   }

//   // then average
//   for (size_t i = 0; i < AvgNum.size(); i++) {
//     if (AvgNum[i] == 0) {
//       AvgPos[i] = VertPos[i];
//       continue;
//     }
//     AvgPos[i] /= AvgNum[i];
//   }
//   for (size_t i = 0; i < AvgPos.size(); i++)
//     VertPos[i] = VertPos[i] * Damp + AvgPos[i] * (1 - Damp);
// }

// template <class ScalarType>
// void GetLaplacianNeighbours(std::vector<Geo::Point3<ScalarType>> &VertPos,
//                             const std::vector<std::vector<int>> &Faces,
//                             std::vector<std::vector<int>> &Neighbors,
//                             const std::vector<std::pair<int, int>>
//                             &LineSmooth =
//                                 std::vector<std::pair<int, int>>()) {
//   Neighbors.clear();
//   Neighbors.resize(VertPos.size());

//   for (size_t i = 0; i < Faces.size(); i++) {
//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int IndexV = Faces[i][j];
//       Neighbors[IndexV].insert(Neighbors[IndexV].end(), Faces[i].begin(),
//                                Faces[i].end());
//     }
//   }

//   // for (size_t i = 0; i < Faces.size(); i++) {
//   //   for (size_t j = 0; j < Faces[i].size(); j++) {
//   //     int IndexV0 = Faces[i][j];
//   //     int IndexV1 = Faces[i][(j+1)%Faces[i].size()];
//   //     Neighbors[IndexV1].push_back(IndexV0);
//   //     Neighbors[IndexV0].push_back(IndexV1);
//   //   }
//   // }

//   // set all features
//   std::set<std::pair<int, int>> LineSmoothSet;
//   std::set<int> FeatureV;
//   for (size_t i = 0; i < LineSmooth.size(); i++) {
//     std::pair<int, int> KeyE = LineSmooth[i];
//     if (KeyE.first > KeyE.second)
//       std::swap(KeyE.first, KeyE.second);
//     LineSmoothSet.insert(KeyE);
//     FeatureV.insert(KeyE.first);
//     FeatureV.insert(KeyE.second);
//   }

//   // otherwise filter neighbours
//   for (size_t i = 0; i < Neighbors.size(); i++) {
//     // if not a feature all neighbours
//     if (FeatureV.count(i) == 0)
//       continue;

//     // filters with features
//     std::vector<int> CurrN = Neighbors[i];
//     Neighbors[i].clear();
//     for (size_t j = 0; j < CurrN.size(); j++) {
//       int IndexV0 = i;
//       int IndexV1 = CurrN[j];
//       std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
//                                std::max(IndexV0, IndexV1));
//       if (LineSmoothSet.count(KeyE) == 0)
//         continue;
//       Neighbors[i].push_back(IndexV1);
//     }
//   }
// }

template <class ScalarType>
void SmoothLaplacianStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::set<int> &FixedV,
                         const std::vector<std::pair<int, int>> &LineSmooth =
                             std::vector<std::pair<int, int>>(),
                         ScalarType Damp = 0.5,
                         bool line_neight_features = true) {
  std::vector<std::vector<int>> Neighbors;
  if (line_neight_features)
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmooth);
  else {
    std::vector<std::pair<int, int>> LineSmoothTest;
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmoothTest);
  }

  // then average the values
  std::vector<Geo::Point3<ScalarType>> LaplVal;
  GetLaplacianValues(VertPos, Neighbors, LaplVal);
  assert(LaplVal.size() == VertPos.size());

  // if (Bilaplacian) {
  //   // encode as difference
  //   for (size_t i = 0; i < LaplVal.size(); i++)
  //     LaplVal[i] = LaplVal[i] - VertPos[i];

  //   // smooth the laplacian to get bilaplacian
  //   std::vector<Geo::Point3<ScalarType>> BiLaplVal;
  //   GetLaplacianValues(LaplVal, Neighbors, BiLaplVal);

  //   // and readd to original vert
  //   for (size_t i = 0; i < LaplVal.size(); i++) {
  //     LaplVal[i] = VertPos[i] + BiLaplVal[i] * 0.001; // BiLaplVal[i];
  //   }
  // }
  for (size_t i = 0; i < LaplVal.size(); i++) {
    if (FixedV.count(i) > 0)
      continue;
    VertPos[i] = VertPos[i] * Damp + LaplVal[i] * (1 - Damp);
  }
}

template <class ScalarType>
void SmoothTaubinStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces,
                      const std::set<int> &FixedV,
                      const std::vector<std::pair<int, int>> &LineSmooth =
                          std::vector<std::pair<int, int>>(),
                      ScalarType Damp = 0.5, bool line_neight_features = true,
                      ScalarType Lambda = 0.6307, ScalarType Mu = -0.6732) {

  std::vector<std::vector<int>> Neighbors;
  if (line_neight_features)
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmooth);
  else {
    std::vector<std::pair<int, int>> LineSmoothTest;
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmoothTest);
  }
  // then average the values
  std::vector<Geo::Point3<ScalarType>> L1;
  GetLaplacianValues(VertPos, Neighbors, L1);
  assert(L1.size() == VertPos.size());

  // smoothing
  std::vector<Geo::Point3<ScalarType>> V1(VertPos.size(),
                                          Geo::Point3<ScalarType>(0, 0, 0));
  for (size_t i = 0; i < L1.size(); i++)
    V1[i] = VertPos[i] + (L1[i] - VertPos[i]) * Lambda;

  // inflation
  std::vector<Geo::Point3<ScalarType>> L2;
  GetLaplacianValues(V1, Neighbors, L2);
  for (size_t i = 0; i < L1.size(); i++)
    L2[i] = L2[i] - V1[i];

  for (size_t i = 0; i < VertPos.size(); i++) {
    if (FixedV.count(i) > 0)
      continue;
    Geo::Point3<ScalarType> TargetP = V1[i] + L2[i] * Mu;
    VertPos[i] = VertPos[i] * Damp + TargetP * (1 - Damp);
  }
}

template <class ScalarType>
void SmoothCentroidStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        const std::set<int> &FixedV,
                        const std::vector<std::pair<int, int>> &LineSmooth =
                            std::vector<std::pair<int, int>>(),
                        ScalarType Damp = 0.5,
                        bool line_neight_features = true) {
  std::vector<Geo::Point3<ScalarType>> AvgPos(VertPos.size(),
                                              Geo::Point3<ScalarType>(0, 0, 0));
  std::vector<int> AvgNum(VertPos.size(), 0);

  std::map<int, int> LineSmoothV;
  if (line_neight_features) {
    for (size_t i = 0; i < LineSmooth.size(); i++) {
      int IndexV0 = LineSmooth[i].first;
      int IndexV1 = LineSmooth[i].second;

      if (LineSmoothV.count(IndexV0) == 0)
        LineSmoothV[IndexV0] = 0;
      if (LineSmoothV.count(IndexV1) == 0)
        LineSmoothV[IndexV1] = 0;
      LineSmoothV[IndexV0]++;
      LineSmoothV[IndexV1]++;
    }
  }

  std::vector<Geo::Point3<ScalarType>> BaryF;
  ComputeFaceBarycenters<ScalarType>(VertPos, Faces, BaryF);

  // cumulate the non line
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      bool IsFixedV = (FixedV.count(IndexV) > 0);
      bool IsLineSmoothV = false;

      if ((LineSmoothV.count(IndexV) > 0) && (LineSmoothV[IndexV] > 1))
        IsLineSmoothV = true;

      if ((!IsFixedV) && (!IsLineSmoothV)) {
        AvgPos[IndexV] += BaryF[i];
        AvgNum[IndexV]++;
      }
    }

  // cumulate the line
  if (line_neight_features) {
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

// template <class ScalarType>
// void SmoothPushStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
//                          const std::vector<std::vector<int>> &Faces,
//                          const std::set<int> &FixedV,
//                          const std::vector<std::pair<int, int>> &LineSmooth =
//                              std::vector<std::pair<int, int>>(),
//                          ScalarType Damp = 0.5) {
//   std::vector<Geo::Point3<ScalarType>> AvgPos(VertPos.size(),
//                                               Geo::Point3<ScalarType>(0, 0,
//                                               0));
//   std::vector<int> AvgNum(VertPos.size(), 0);

//   std::set<int> LineSmoothV;
//   for (size_t i = 0; i < LineSmooth.size(); i++) {
//     LineSmoothV.insert(LineSmooth[i].first);
//     LineSmoothV.insert(LineSmooth[i].second);
//   }

//   std::vector<ScalarType> VertL(VertPos.size(),0);
//   std::vector<int> VertLNum(VertPos.size(),0);
//   std::vector<std::vector<int> > VertN(VertPos.size(),std::vector<int>());
//    for (size_t i = 0; i < Faces.size(); i++)
//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int IndexV0 = Faces[i][j];
//       int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
//       ScalarType Len=(VertPos[IndexV0]-VertPos[IndexV1]).Norm();
//       VertL[IndexV0]+=Len;
//       VertL[IndexV1]+=Len;
//       VertLNum[IndexV0]++;
//       VertLNum[IndexV1]++;
//       VertN[IndexV0].push_back(IndexV1);
//       VertN[IndexV1].push_back(IndexV0);
//     }
//   for (size_t i = 0; i < VertL.size(); i++)
//     VertL[i]/=VertLNum[i];

//   for (size_t i = 0; i < VertN.size(); i++)
//   {
//     std::sort(VertN[i].begin(),VertN[i].end());
//     VertN[i].erase(std::unique(VertN[i].begin(), VertN[i].end()),
//     VertN[i].end());
//   }

//   //compute the position
//   for (size_t i = 0; i < VertN.size(); i++)
//     for (size_t j = 0; j < VertN[i].size(); j++) {
//       int IndexV0 = i;
//       int IndexV1 = VertN[i][j];
//       assert(IndexV0!=IndexV1);
//       Geo::Point3<ScalarType> Dir=VertPos[IndexV1]-VertPos[IndexV0];
//       Dir.Normalize();
//       Dir=Dir*VertL[IndexV0];
//       Geo::Point3<ScalarType> TargetP=VertPos[IndexV0]+Dir;
//       bool IsFixedV = (FixedV.count(IndexV1) > 0);
//       bool IsLineSmoothV = (LineSmoothV.count(IndexV1) > 0);
//       if ((!IsFixedV) && (!IsLineSmoothV)) {
//         AvgPos[IndexV1] += TargetP;
//         AvgNum[IndexV1]++;
//       }
//     }

//   // cumulate the line
//   for (size_t i = 0; i < LineSmooth.size(); i++) {
//     int IndexV0 = LineSmooth[i].first;
//     int IndexV1 = LineSmooth[i].second;
//     bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
//     bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
//     if (!IsFixedV0) {
//       AvgPos[IndexV0] += VertPos[IndexV1];
//       AvgNum[IndexV0]++;
//     }
//     if (!IsFixedV1) {
//       AvgPos[IndexV1] += VertPos[IndexV0];
//       AvgNum[IndexV1]++;
//     }
//   }

//   // then average
//   for (size_t i = 0; i < AvgNum.size(); i++) {
//     if (AvgNum[i] == 0) {
//       AvgPos[i] = VertPos[i];
//       continue;
//     }
//     AvgPos[i] /= AvgNum[i];
//   }
//   for (size_t i = 0; i < AvgPos.size(); i++)
//     VertPos[i] = VertPos[i] * Damp + AvgPos[i] * (1 - Damp);
// }

// template <class ScalarType>
// void GetTargetPoly(const std::vector<Geo::Point3<ScalarType>> &VertPos,
//                    const std::vector<std::vector<int>> &Faces,
//                    const std::set<int> &LineSmoothV,
//                    std::vector<std::vector<Geo::Point3<ScalarType> > >
//                    &TargetPoly, std::vector<ScalarType> &Weights,
//                    std::vector<Geo::Point3<ScalarType>> &TargetN =
//                    std::vector<Geo::Point3<ScalarType>>())
// {
//   ScalarType WCorrect=100;
//   ScalarType WDefault=1;
//   Weights.clear();
//   Weights.resize(Faces.size(),WDefault);

//   bool InitializedPoly=(TargetPoly.size()==Faces.size());
//   bool HasTargetN=(TargetN.size()==Faces.size());

//   if (!InitializedPoly)TargetPoly.clear();
//   if (!HasTargetN)TargetPoly.clear();

//   int PreFolded=0;
//   // cumulate the non line
//   for (size_t i = 0; i < Faces.size(); i++)
//   {

//     std::vector<Geo::Point3<ScalarType>> FacePos;
//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int IndexV = Faces[i][j];
//       FacePos.push_back(VertPos[IndexV]);
//     }

//     std::vector<Geo::Point3<ScalarType>> CurrPoly = BestAlignedPoly(FacePos);

//     //std::cout<<"Test 0"<<std::endl;
//     //keep the old one if flip
//     if (InitializedPoly && HasTargetN)
//     {
//       Geo::Point3 < ScalarType > TestN = TargetN[i];
//       Geo::Point3<ScalarType> PolyN=NormalPoly(CurrPoly);

//       //std::cout<<"TEST NORMALs"<<std::endl;
//       //update the target Norm if needed
//       if (Dot(TestN,PolyN)>0)
//         TargetPoly[i]=CurrPoly;
//       else
//       {
//         std::cout<<"Folded"<<std::endl;
//         Geo::Point3<ScalarType> PolyN1=NormalPoly(TargetPoly[i]);
//         if (Dot(TestN,PolyN1)<0)
//         {
//           assert(0);
//           std::cout<<"Folded 1"<<std::endl;
//         }
//         Weights[i]=WCorrect;
//       }
//     }else
//     {
//       TargetPoly.push_back(CurrPoly);

//       if (HasTargetN)
//       {
//         Geo::Point3 < ScalarType > TestN = TargetN[i];
//         Geo::Point3<ScalarType> PolyN1=NormalPoly(CurrPoly);
//         if (Dot(TestN,PolyN1)<0)
//         PreFolded++;
//       }
//     }
//     //std::cout<<"Test 1"<<std::endl;

//     // translate if on a border
//     Geo::Point3<ScalarType> TranslateV(0, 0, 0);
//     int TranslateNum = 0;
//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int IndexV = Faces[i][j];
//       bool IsLineSmoothV = (LineSmoothV.count(IndexV) > 0);
//       if (IsLineSmoothV) {
//         TranslateV += VertPos[IndexV] - TargetPoly[i][j];
//         TranslateNum++;
//       }
//     }
//     //bool HasB = false;
//     if (TranslateNum > 0) {
//       //HasB = true;
//       TranslateV /= TranslateNum;
//       for (size_t j = 0; j < TargetPoly[i].size(); j++)
//       TargetPoly[i][j] += TranslateV;
//     }

//     assert(Faces[i].size() == TargetPoly[i].size());
//   }
//   std::cout<<"Prefolded:"<<PreFolded<<" out of
//   "<<TargetPoly.size()<<std::endl;
// }

template <class ScalarType>
void GetTargetPoly(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::set<int> &LineSmoothV,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &TargetPoly,
    std::vector<ScalarType> &Weights, bool weight_by_Aspect = false) {
  ScalarType WDefault = 1;
  Weights.clear();
  Weights.resize(Faces.size(), WDefault);

  TargetPoly.clear();

  // cumulate the non line
  for (size_t i = 0; i < Faces.size(); i++) {

    std::vector<Geo::Point3<ScalarType>> FacePos;
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      FacePos.push_back(VertPos[IndexV]);
    }

    std::vector<Geo::Point3<ScalarType>> CurrPoly;
    ScalarType W = BestAlignedPoly(FacePos, CurrPoly);
    TargetPoly.push_back(CurrPoly);

    // translate if on a border
    Geo::Point3<ScalarType> TranslateV(0, 0, 0);
    int TranslateNum = 0;
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      bool IsLineSmoothV = (LineSmoothV.count(IndexV) > 0);
      if (IsLineSmoothV) {
        TranslateV += VertPos[IndexV] - TargetPoly[i][j];
        TranslateNum++;
      }
    }
    // bool HasB = false;
    if (TranslateNum > 0) {
      // HasB = true;
      TranslateV /= TranslateNum;
      for (size_t j = 0; j < TargetPoly[i].size(); j++)
        TargetPoly[i][j] += TranslateV;
    }
    // set the weight as the distance
    if (weight_by_Aspect) {
      ScalarType Err = 0;
      for (size_t j = 0; j < FacePos.size(); j++)
        Err += pow((FacePos[j] - TargetPoly[i][j]).Norm(), 2);
      Weights[i] = Err;
    }
  }
}

template <class ScalarType>
void SmoothPolyAlignStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::set<int> &FixedV,
                         const std::vector<std::pair<int, int>> &LineSmooth =
                             std::vector<std::pair<int, int>>(),
                         ScalarType Damp = 0.5,
                         bool line_neight_features = true) {
  std::vector<Geo::Point3<ScalarType>> AvgPos(VertPos.size(),
                                              Geo::Point3<ScalarType>(0, 0, 0));
  std::vector<ScalarType> AvgW(VertPos.size(), 0);
  // std::vector<std::vector<bool>> FromBorder(VertPos.size());

  std::map<int, int> LineSmoothV;
  std::set<int> LineSmoothVSet;

  //if (line_neight_features) {
    for (size_t i = 0; i < LineSmooth.size(); i++) {
      int IndexV0 = LineSmooth[i].first;
      int IndexV1 = LineSmooth[i].second;

      if (LineSmoothV.count(IndexV0) == 0)
        LineSmoothV[IndexV0] = 0;
      if (LineSmoothV.count(IndexV1) == 0)
        LineSmoothV[IndexV1] = 0;

      LineSmoothV[IndexV0]++;
      LineSmoothV[IndexV1]++;

      LineSmoothVSet.insert(IndexV0);
      LineSmoothVSet.insert(IndexV1);
    }
  //}

  // std::vector<std::vector<Geo::Point3<ScalarType> > > TargetPoly;
  std::vector<ScalarType> Weights;
  std::vector<std::vector<Geo::Point3<ScalarType>>> TargetPoly;
  GetTargetPoly<ScalarType>(VertPos, Faces, LineSmoothVSet, TargetPoly,
                            Weights); //,TargetN );
  for (size_t i = 0; i < Faces.size(); i++) {
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      bool IsFixedV = (FixedV.count(IndexV) > 0);
      bool IsLineSmoothV = false;

      if ((line_neight_features)&&(LineSmoothV.count(IndexV) > 0) && (LineSmoothV[IndexV] > 1))
        IsLineSmoothV = true;

      if ((!IsFixedV) && (!IsLineSmoothV)) {
        AvgPos[IndexV] += (TargetPoly[i][j] * Weights[i]);
        AvgW[IndexV] += Weights[i];
        // FromBorder[IndexV].push_back(HasB);
      }
    }
  }

  // cumulate the line
  if (line_neight_features) {
    for (size_t i = 0; i < LineSmooth.size(); i++) {
      int IndexV0 = LineSmooth[i].first;
      int IndexV1 = LineSmooth[i].second;
      bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
      bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
      if (!IsFixedV0) {
        AvgPos[IndexV0] += (VertPos[IndexV1] * Weights[i]);
        // AvgNum[IndexV0]++;
        AvgW[IndexV0] += Weights[i];
      }
      if (!IsFixedV1) {
        AvgPos[IndexV1] += (VertPos[IndexV0] * Weights[i]);
        // AvgNum[IndexV1]++;
        AvgW[IndexV1] += Weights[i];
      }
    }
  }

  // add the laplacian term
  std::vector<std::vector<int>> Neighbors;
  if (line_neight_features)
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmooth);
  else {
    std::vector<std::pair<int, int>> LineSmoothTest;
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmoothTest);
  }
  // get laplaxian values
  std::vector<Geo::Point3<ScalarType>> LaplVal;
  GetLaplacianValues(VertPos, Neighbors, LaplVal);

  ScalarType LaplW = 0.1;
  // then average
  for (size_t i = 0; i < AvgW.size(); i++) {
    if (AvgW[i] == 0) {
      AvgPos[i] = VertPos[i];
      continue;
    }
    // AvgPos[i] /= AvgNum[i];
    AvgPos[i] /= AvgW[i];
    AvgPos[i] = AvgPos[i] * (1 - LaplW) + LaplVal[i] * LaplW;
  }
  for (size_t i = 0; i < AvgPos.size(); i++) {
    Geo::Point3<ScalarType> OldP = VertPos[i];
    VertPos[i] = VertPos[i] * Damp + AvgPos[i] * (1 - Damp);
  }
}

// template <class ScalarType>
// void SmoothPolyAlignStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
//                          const std::vector<std::vector<int>> &Faces,
//                          const std::set<int> &FixedV,
//                          const std::vector<std::pair<int, int>> &LineSmooth =
//                              std::vector<std::pair<int, int>>(),
//                          ScalarType Damp = 0.5) {

//   std::vector<std::vector<Geo::Point3<ScalarType>>> NeighPos(VertPos.size());

//   // std::vector<int> AvgNum(VertPos.size(), 0);
//   //std::vector<std::vector<bool>> FromBorder(VertPos.size());

//   std::set<int> LineSmoothV;
//   for (size_t i = 0; i < LineSmooth.size(); i++) {
//     LineSmoothV.insert(LineSmooth[i].first);
//     LineSmoothV.insert(LineSmooth[i].second);
//   }

//   // cumulate the non line
//   for (size_t i = 0; i < Faces.size(); i++) {

//     std::vector<Geo::Point3<ScalarType>> FacePos;
//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int IndexV = Faces[i][j];
//       FacePos.push_back(VertPos[IndexV]);
//     }

//     std::vector<Geo::Point3<ScalarType>> TargetF = BestAlignedPoly(FacePos);

//     // translate if on a border
//     Geo::Point3<ScalarType> TranslateV(0, 0, 0);
//     int TranslateNum = 0;
//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int IndexV = Faces[i][j];
//       bool IsLineSmoothV = (LineSmoothV.count(IndexV) > 0);
//       if (IsLineSmoothV) {
//         TranslateV += VertPos[IndexV] - TargetF[j];
//         TranslateNum++;
//       }
//     }

//     //bool HasB = false;
//     if (TranslateNum > 0) {
//       //HasB = true;
//       TranslateV /= TranslateNum;
//       for (size_t j = 0; j < TargetF.size(); j++)
//         TargetF[j] += TranslateV;
//     }

//     assert(Faces[i].size() == TargetF.size());

//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int IndexV = Faces[i][j];
//       bool IsFixedV = (FixedV.count(IndexV) > 0);
//       bool IsLineSmoothV = (LineSmoothV.count(IndexV) > 0);
//       if ((!IsFixedV) && (!IsLineSmoothV)) {
//         NeighPos[IndexV].push_back(TargetF[j]);
//         // AvgPos[IndexV] += TargetF[j];
//         // AvgNum[IndexV]++;
//         //FromBorder[IndexV].push_back(HasB);
//       }
//     }
//   }

//   // cumulate the line
//   for (size_t i = 0; i < LineSmooth.size(); i++) {
//     int IndexV0 = LineSmooth[i].first;
//     int IndexV1 = LineSmooth[i].second;
//     bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
//     bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
//     if (!IsFixedV0) {
//       // AvgPos[IndexV0] = VertPos[IndexV1];
//       // AvgNum[IndexV0]++;
//       NeighPos[IndexV0].push_back(VertPos[IndexV1]);
//     }
//     if (!IsFixedV1) {
//       // AvgPos[IndexV1] += VertPos[IndexV0];
//       // AvgNum[IndexV1]++;
//       NeighPos[IndexV1].push_back(VertPos[IndexV0]);
//     }
//   }

//   std::vector<Geo::Point3<ScalarType>> TargetPos = VertPos;
//   // then average
//   for (size_t i = 0; i < NeighPos.size(); i++) {
//     int IndexV = i;
//     bool IsLineSmoothV = (LineSmoothV.count(IndexV) > 0);

//     Geo::Point3<ScalarType> AvgTarget(0, 0, 0);

//     bool IsFixedV = (FixedV.count(IndexV) > 0);

//     if (!IsFixedV) {
//       if (IsLineSmoothV) {
//         for (size_t j = 0; j < NeighPos[IndexV].size(); j++)
//           AvgTarget += NeighPos[IndexV][j];

//         AvgTarget /= NeighPos[IndexV].size();
//         TargetPos[IndexV] = AvgTarget;
//       } else {
//         // bool HasB=false;
//         // for (size_t j = 0; j < FromBorder[IndexV].size(); j++)
//         //   HasB|=FromBorder[IndexV][j];

//         // if (!HasB)
//         // {
//         for (size_t j = 0; j < NeighPos[IndexV].size(); j++)
//           AvgTarget += NeighPos[IndexV][j];

//         AvgTarget /= NeighPos[IndexV].size();
//         // }
//         // else {
//         //   int Num=0;
//         //   for (size_t j = 0; j < NeighPos[IndexV].size(); j++)
//         //   {
//         //     if (!FromBorder[IndexV][j])continue;
//         //     AvgTarget += NeighPos[IndexV][j];
//         //     Num++;
//         //   }

//         //   AvgTarget /= Num;
//         // }

//         TargetPos[IndexV] = AvgTarget;
//       }
//     }
//     // if (AvgNum[i] == 0) {
//     //   AvgPos[i] = VertPos[i];
//     //   continue;
//     // }
//     // AvgPos[i] /= AvgNum[i];
//   }
//   for (size_t i = 0; i < TargetPos.size(); i++) {
//     VertPos[i] = VertPos[i] * Damp + TargetPos[i] * (1 - Damp);
//   }
// }

enum SmoothType {
  SMLaplacian,
  // SMBiLaplacian,
  SMTaubin,
  SMBestPoly,
  SMCentroid
};

template <class ScalarType> struct SmoothParam {
  int NumIte;
  std::vector<int> FixedVert;
  std::vector<std::pair<int, int>> Features;
  bool line_border;
  bool post_project;
  ScalarType Damp;
  std::vector<Geo::Point3<ScalarType>> *ReprojVert;
  std::vector<std::vector<int>> *ReprojFaces;
  std::vector<std::pair<int, int>> ReprojFeatures;
  int ReprojStep;
  SmoothType SType;
  bool check_AspectR;
  bool no_reproj_feat_endpoints;
  bool line_neight_features;

  void Clear() {
    FixedVert.clear();
    Features.clear();
  }

  SmoothParam() {
    NumIte = 10;
    line_border = true;

    Damp = 0.5;
    ReprojVert = NULL;
    ReprojFaces = NULL;
    post_project = false;
    ReprojStep = 1;
    SType = SMLaplacian;
    check_AspectR = false;
    no_reproj_feat_endpoints = true;
    line_neight_features = true;
  }
};

template <class ScalarType>
void PostProject(std::vector<Geo::Point3<ScalarType>> &VertPosTarget,
                 std::vector<Geo::Point3<ScalarType>> &VertProj,
                 const std::vector<std::vector<int>> &FacesProj) {
  std::vector<Geo::Point3<ScalarType>> AvgDiff(
      VertProj.size(), Geo::Point3<ScalarType>(0, 0, 0));
  std::vector<int> AvgNum(VertProj.size(), 0);

  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(VertProj, FacesProj);

  for (size_t i = 0; i < VertPosTarget.size(); i++) {
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(VertProj, FacesProj, VertPosTarget[i],
                        TriGrid.bbox.Diag(), closestF, closestPt);

    Geo::Point3<ScalarType> Diff = VertPosTarget[i] - closestPt;

    for (size_t j = 0; j < FacesProj[closestF].size(); j++) {
      int IndexV = FacesProj[closestF][j];
      // Geo::Point3<ScalarType> TargetP=VertProj[IndexV]+Diff;
      AvgDiff[IndexV] += Diff;
      AvgNum[IndexV]++;
    }
  }

  std::set<int> FixedV;
  int stepDiff = 3;
  for (size_t i = 0; i < stepDiff; i++)
    SmoothLaplacianStep(AvgDiff, FacesProj, FixedV);

  // then smooth the displacement
  for (size_t i = 0; i < AvgNum.size(); i++) {
    if (AvgNum[i] == 0) {
      // AvgDiff[i] = 0;
      continue;
    }
    AvgDiff[i] /= AvgNum[i];
  }
  for (size_t i = 0; i < AvgDiff.size(); i++)
    VertProj[i] += AvgDiff[i]; // VertProj[i] * Damp + AvgPos[i] * (1 - Damp);
}

template <class ScalarType>
void SampleTargetNormals(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                         const SmoothParam<ScalarType> &SParam,
                         std::vector<Geo::Point3<ScalarType>> &TargetN) {
  std::vector<Geo::Point3<ScalarType>> BaryF;
  ComputeFaceBarycenters(VertPos, Faces, BaryF);

  assert((*SParam.ReprojVert).size() > 0);
  assert((*SParam.ReprojFaces).size() > 0);

  std::vector<Geo::Point3<ScalarType>> ReprojectN;
  ComputeFaceNormals((*SParam.ReprojVert), (*SParam.ReprojFaces), ReprojectN);

  TargetN.clear();

  for (size_t i = 0; i < BaryF.size(); i++) {
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest((*SParam.ReprojVert), (*SParam.ReprojFaces), BaryF[i],
                        TriGrid.bbox.Diag(), closestF, closestPt);
    assert(closestF >= 0);
    assert(closestF < ReprojectN.size());
    TargetN.push_back(ReprojectN[closestF]);
  }
}

template <class ScalarType>
ScalarType
QualityQuadsOnVert(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const std::vector<std::vector<int>> &VFAdj,
                   const int IndexV) {
  assert(IndexV >= 0);
  assert(IndexV < VFAdj.size());
  assert(VFAdj[IndexV].size() > 0);
  ScalarType WorstQ = 1;
  // get the worst of neighboring faces
  for (size_t i = 0; i < VFAdj[IndexV].size(); i++)
    WorstQ =
        std::min(WorstQ, ScaledJacobianQuad(VertPos, Faces, VFAdj[IndexV][i]));

  // WorstQ = std::min(WorstQ,
  //                   AspectRatioAlignPoly(VertPos, Faces, VFAdj[IndexV][i]));
  return WorstQ;
}

template <class ScalarType>
int ConcaveQuadsOnVert(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::vector<int>> &VFAdj,
                       const int IndexV) {
  assert(IndexV >= 0);
  assert(IndexV < VFAdj.size());
  assert(VFAdj[IndexV].size() > 0);
  int FlippedQ = 0;
  // get the worst of neighboring faces
  for (size_t i = 0; i < VFAdj[IndexV].size(); i++)
    if (HasConcaveV(VertPos, Faces, VFAdj[IndexV][i]))
      FlippedQ++;

  return FlippedQ;
}

template <class ScalarType>
void SetPositionIfBetter(const std::vector<Geo::Point3<ScalarType>> &OldPos,
                         std::vector<Geo::Point3<ScalarType>> &NewPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::vector<std::vector<int>> &VFAdj,
                         const std::vector<Point3<ScalarType>> OldFaceNormals) {

  // compute the original quality per vert
  std::vector<std::pair<ScalarType, int>> VertAR;
  for (size_t i = 0; i < VFAdj.size(); i++) {
    ScalarType VertQ = QualityQuadsOnVert(OldPos, Faces, VFAdj, i);
    VertAR.push_back(std::pair<ScalarType, int>(VertQ, i));
  }

  // set the old position
  std::vector<Geo::Point3<ScalarType>> SolvedPos = OldPos;

  // then sort to start from the worse
  std::sort(VertAR.begin(), VertAR.end());
  for (size_t i = 0; i < VertAR.size(); i++) {
    int IndexV = VertAR[i].second;

    // compute the old quality
    ScalarType OldQ = QualityQuadsOnVert(SolvedPos, Faces, VFAdj, IndexV);
    int oldFlip = ConcaveQuadsOnVert(SolvedPos, Faces, VFAdj, IndexV);

    // then try to move
    SolvedPos[IndexV] = NewPos[IndexV];

    // see if flipped faces
    bool has_flip = false;
    for (size_t j = 0; j < VFAdj[IndexV].size(); j++) {
      int IndexF = VFAdj[IndexV][j];
      Quad3<ScalarType> Q = QuadOfFace(SolvedPos, Faces, IndexF);
      Geo::Point3<ScalarType> N = Q.N();
      has_flip |= (Dot(N, OldFaceNormals[IndexF]) <= 0);
    }
    if (has_flip) {
      // restore old pos
      // std::cout<<"RESTORE BECAUSE OF FLIP"<<std::endl;
      SolvedPos[IndexV] = OldPos[IndexV];
      continue;
    }

    ScalarType NewQ = QualityQuadsOnVert(SolvedPos, Faces, VFAdj, IndexV);
    int newFlip = ConcaveQuadsOnVert(SolvedPos, Faces, VFAdj, IndexV);

    OldQ -= oldFlip;
    NewQ -= newFlip;
    // if (newFlip > oldFlip) {
    //   // std::cout<<"RESTORE BECAUSE OF CONCAVE"<<std::endl;
    //   //  restore old pos
    //   SolvedPos[IndexV] = OldPos[IndexV];
    //   continue;
    // }
    // if ((NewQ < 0.25) && (NewQ < OldQ)) {
    if ((OldQ < 0.25) && (NewQ < OldQ)) {
      // restore old pos
      SolvedPos[IndexV] = OldPos[IndexV];
      continue;
    }
  }
  NewPos = SolvedPos;
}

// template <class ScalarType>
// void ReprojectIfNeeded2(const SmoothParam<ScalarType> &SParam,
//                        std::vector<Geo::Point3<ScalarType>> &VertPos,
//                        const std::vector<std::vector<int>> &Faces,
//                        const std::vector<int> &ReprojFeatureV,
//                        Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
//                        std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
//                        Geo::SegmentSpatialIndex<ScalarType> &SegGrid,
//                        const std::set<Geo::Point3<ScalarType>> &NoReprPoints)
//                        {

//   bool do_reproject = ((SParam.ReprojStep > 0) && (SParam.ReprojVert != NULL)
//   &&
//                        (SParam.ReprojFaces != NULL));

//   bool post_project = (do_reproject && SParam.post_project);

//   // save the fixed one
//   std::vector<Geo::Point3<ScalarType>> FixedPos;
//   for (size_t i = 0; i < SParam.FixedVert.size(); i++)
//     FixedPos.push_back(VertPos[SParam.FixedVert[i]]);

//   if (post_project)
//     Geo::AttractProject<ScalarType>(*SParam.ReprojVert, VertPos, Faces);

//   if (do_reproject) {

//     Geo::ProjectOn<ScalarType>(VertPos, Faces, ReprojFeatureV,
//                                *SParam.ReprojVert, *SParam.ReprojFaces,
//                                TriGrid, FeatureSegs, SegGrid,NoReprPoints);
//   }
//   // restore the fixed one
//   for (size_t i = 0; i < SParam.FixedVert.size(); i++)
//     VertPos[SParam.FixedVert[i]] = FixedPos[i];
// }

template <class ScalarType>
void ReprojectIfNeeded(const SmoothParam<ScalarType> &SParam,
                       std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                       Geo::FeatureReprojector<ScalarType> *FeatureReproj) {

  bool do_reproject = ((SParam.ReprojStep > 0) && (SParam.ReprojVert != NULL) &&
                       (SParam.ReprojFaces != NULL));

  bool post_project = (do_reproject && SParam.post_project);

  // save the fixed one
  std::vector<Geo::Point3<ScalarType>> FixedPos;
  for (size_t i = 0; i < SParam.FixedVert.size(); i++)
    FixedPos.push_back(VertPos[SParam.FixedVert[i]]);

  if (post_project)
    Geo::AttractProject<ScalarType>(*SParam.ReprojVert, VertPos, Faces);

  if (do_reproject) {

    // Geo::ProjectOn<ScalarType>(VertPos, Faces, ReprojFeatureV,
    //                            *SParam.ReprojVert, *SParam.ReprojFaces,
    //                            TriGrid, FeatureSegs, SegGrid,NoReprPoints);
    Geo::ProjectOn<ScalarType>(VertPos, Faces, *SParam.ReprojVert,
                               *SParam.ReprojFaces, TriGrid);
    if (FeatureReproj != NULL)
      (*FeatureReproj).Reproject();
  }
  // restore the fixed one
  for (size_t i = 0; i < SParam.FixedVert.size(); i++)
    VertPos[SParam.FixedVert[i]] = FixedPos[i];
}

// template <class ScalarType>
// void SmoothLaplacian2(
//     std::vector<Geo::Point3<ScalarType>> &VertPos,
//     const std::vector<std::vector<int>> &Faces,
//     const SmoothParam<ScalarType> &SParam = SmoothParam<ScalarType>()) {

//   //  save the features
//   std::vector<Point3<ScalarType>> OldFaceNormals;
//   ComputeFaceNormals(VertPos, Faces, OldFaceNormals);

//   bool has_reproj_features = ((SParam.Features.size() >
//   0)&&(SParam.ReprojFeatures.size() > 0));
//   Geo::SegmentSpatialIndex<ScalarType> SegGrid;
//   std::vector<Geo::Segment3<ScalarType>> FeatureSegs;

//   Geo::TriangleSpatialIndex<ScalarType> TriGrid;
//   bool reproject = false;

//   if ((SParam.ReprojVert != NULL) && (SParam.ReprojFaces != NULL)) {
//     TriGrid.Init(*SParam.ReprojVert, *SParam.ReprojFaces);
//     reproject = true;
//   }

//   std::vector<int> ReprojectVertIdx;

//   std::set<int> FixedVertSet(SParam.FixedVert.begin(),
//   SParam.FixedVert.end());

//   std::set<Geo::Point3<ScalarType> > NoReprPoints;
//   if ((SParam.no_reproj_feat_endpoints)&&(SParam.ReprojVert!=NULL))
//   {
//     std::vector<int> Valence1Idx;
//     Geo::MeshFeatures<ScalarType>::VerticesWithValence1(SParam.ReprojFeatures,
//     Valence1Idx);
//     for (size_t j=0;j<Valence1Idx.size();j++)
//     {
//       int IndexV=Valence1Idx[j];
//       assert(IndexV>=0);
//       assert(IndexV<(*SParam.ReprojVert).size());
//       NoReprPoints.insert((*SParam.ReprojVert)[IndexV]);
//     }
//     //std::sort(Valence1Idx.begin(), Valence1Idx.end());
//   }

//   if ((has_reproj_features) && (SParam.ReprojVert != NULL) &&
//       (SParam.ReprojFaces != NULL)) {

//     Geo::InitFromMeshFeatures<ScalarType>(
//         *SParam.ReprojVert, SParam.ReprojFeatures, FeatureSegs, SegGrid);

//     for (size_t i = 0; i < SParam.Features.size(); i++) {
//       int IndexV0=SParam.Features[i].first;
//       int IndexV1=SParam.Features[i].second;
//       assert(IndexV0>=0);
//       assert(IndexV1>=0);
//       ReprojectVertIdx.push_back(IndexV0);
//       ReprojectVertIdx.push_back(IndexV1);
//     }
//     std::sort(ReprojectVertIdx.begin(), ReprojectVertIdx.end());
//     ReprojectVertIdx.erase(std::unique(ReprojectVertIdx.begin(),
//     ReprojectVertIdx.end()), ReprojectVertIdx.end());

//     GetToReprojectVert(VertPos, ReprojectVertIdx, FeatureSegs, SegGrid);
//   }

//   std::vector<std::vector<int>> VF, VFIndex;
//   if (SParam.check_AspectR)
//     getVFAdjacency(Faces, VF, VFIndex);

//   for (int s = 0; s < SParam.NumIte; s++) {

//     std::vector<Geo::Point3<ScalarType>> OldPos;
//     if (SParam.check_AspectR)
//       OldPos = VertPos;

//     if (SParam.SType == SMLaplacian)
//       SmoothLaplacianStep(VertPos, Faces, FixedVertSet, SParam.Features,
//                           SParam.Damp);

//     // if (SParam.SType == SMBiLaplacian)
//     //   SmoothLaplacianStep(VertPos, Faces, FixedVertSet, SParam.Features,
//     //                       SParam.Damp, true);

//     if (SParam.SType == SMTaubin)
//       SmoothTaubinStep(VertPos, Faces, FixedVertSet, SParam.Features,
//                        SParam.Damp);

//     if (SParam.SType == SMBestPoly)
//       SmoothPolyAlignStep(VertPos, Faces, FixedVertSet, SParam.Features,
//                           SParam.Damp);

//     if (SParam.SType == SMCentroid)
//       SmoothCentroidStep(VertPos, Faces, FixedVertSet, SParam.Features,
//                          SParam.Damp);

//     ReprojectIfNeeded2(SParam, VertPos, Faces, ReprojectVertIdx, TriGrid,
//                       FeatureSegs, SegGrid, NoReprPoints);

//     if (SParam.check_AspectR)
//       SetPositionIfBetter(OldPos, VertPos, Faces, VF, OldFaceNormals);
//   }

// }

template <class ScalarType>
void SmoothLaplacian(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    SmoothParam<ScalarType> &SParam = SmoothParam<ScalarType>()) {

  // make unique
  Geo::EdgeMeshFunctions<ScalarType>::MakeUnique(SParam.Features);
  Geo::EdgeMeshFunctions<ScalarType>::MakeUnique(SParam.ReprojFeatures);

  //  save the features
  std::vector<Point3<ScalarType>> OldFaceNormals;
  if (SParam.check_AspectR)
    ComputeFaceNormals(VertPos, Faces, OldFaceNormals);

  bool has_reproj_features =
      ((SParam.Features.size() > 0) && (SParam.ReprojFeatures.size() > 0));
  // Geo::SegmentSpatialIndex<ScalarType> SegGrid;
  // std::vector<Geo::Segment3<ScalarType>> FeatureSegs;

  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  bool reproject = false;

  if ((SParam.ReprojVert != NULL) && (SParam.ReprojFaces != NULL)) {
    TriGrid.Init(*SParam.ReprojVert, *SParam.ReprojFaces);
    reproject = true;
  }

  // std::vector<int> ReprojectVertIdx;

  std::set<int> FixedVertSet(SParam.FixedVert.begin(), SParam.FixedVert.end());

  Geo::FeatureReprojector<ScalarType> *FeatureReproj = NULL;

  if ((has_reproj_features) && (SParam.ReprojVert != NULL) &&
      (SParam.ReprojFaces != NULL)) {
    FeatureReproj = new Geo::FeatureReprojector<ScalarType>(
        VertPos, SParam.Features, (*SParam.ReprojVert), SParam.ReprojFeatures);
    (*FeatureReproj).Init();
  }

  std::vector<std::vector<int>> VF, VFIndex;
  if (SParam.check_AspectR)
    getVFAdjacency(Faces, VF, VFIndex);

  for (int s = 0; s < SParam.NumIte; s++) {

    std::vector<Geo::Point3<ScalarType>> OldPos;
    if (SParam.check_AspectR)
      OldPos = VertPos;

    if (SParam.SType == SMLaplacian)
      SmoothLaplacianStep(VertPos, Faces, FixedVertSet, SParam.Features,
                          SParam.Damp, SParam.line_neight_features);

    // if (SParam.SType == SMBiLaplacian)
    //   SmoothLaplacianStep(VertPos, Faces, FixedVertSet, SParam.Features,
    //                       SParam.Damp, true);

    if (SParam.SType == SMTaubin)
      SmoothTaubinStep(VertPos, Faces, FixedVertSet, SParam.Features,
                       SParam.Damp, SParam.line_neight_features);

    if (SParam.SType == SMBestPoly)
      SmoothPolyAlignStep(VertPos, Faces, FixedVertSet, SParam.Features,
                          SParam.Damp, SParam.line_neight_features);

    if (SParam.SType == SMCentroid)
      SmoothCentroidStep(VertPos, Faces, FixedVertSet, SParam.Features,
                         SParam.Damp,
                         SParam.line_neight_features); //);

    ReprojectIfNeeded(SParam, VertPos, Faces, TriGrid, FeatureReproj);

    if (SParam.check_AspectR)
      SetPositionIfBetter(OldPos, VertPos, Faces, VF, OldFaceNormals);
  }

  if (FeatureReproj != NULL)
    delete (FeatureReproj);
}

template <class ScalarType>
void SmoothLaplacian(
    std::vector<Geo::Point2<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const SmoothParam<ScalarType> &SParam = SmoothParam<ScalarType>()) {
  std::vector<Geo::Point3<ScalarType>> VertPos3D;
  for (size_t i = 0; i < VertPos.size(); i++)
    VertPos3D.push_back(Geo::Point3<ScalarType>(VertPos[i].X, VertPos[i].Y, 0));

  SmoothParam<ScalarType> SParam1 = SParam;
  SParam1.ReprojVert = NULL;
  SParam1.ReprojFaces = NULL;
  SParam1.post_project = false;
  SmoothLaplacian(VertPos3D, Faces, SParam1);

  for (size_t i = 0; i < VertPos3D.size(); i++) {
    VertPos[i].X = VertPos3D[i].X;
    VertPos[i].Y = VertPos3D[i].Y;
  }
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
                             const std::vector<int> &FixedVert,
                             const std::vector<std::pair<int, int>> &Features =
                                 std::vector<std::pair<int, int>>()) {

  ScalarType WeightFixed = 1000000;
  ScalarType WeightLap = 1;

  std::vector<std::vector<int>> KernelV;
  GetLaplacianNeighbours(VertPos, Faces, KernelV, Features);

  std::vector<Geo::LinearConstraint<ScalarType>> LinConstr;

  std::set<int> FixedSet(FixedVert.begin(), FixedVert.end());
  for (size_t i = 0; i < KernelV.size(); i++) {
    int IndexV0 = i;
    if (FixedSet.count(IndexV0) > 0)
      continue;
    for (size_t c = 0; c < 3; c++) {
      Geo::LinearConstraint<ScalarType> Constr;
      Constr.Weight = WeightLap;
      Constr.RHSide = 0;
      Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV0, c, -1));
      ScalarType OtherW = (ScalarType)1 / (ScalarType)KernelV[i].size();

      for (size_t j = 0; j < KernelV[i].size(); j++) {
        int IndexV1 = KernelV[i][j];

        assert(IndexV0 != IndexV1);

        Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1, c, OtherW));
      }
      LinConstr.push_back(Constr);
    }
  }

  for (size_t i = 0; i < FixedVert.size(); i++) {
    int IndexV = FixedVert[i];
    for (size_t c = 0; c < 3; c++) {
      Geo::LinearConstraint<ScalarType> Constr;
      Constr.RHSide = VertPos[IndexV].V(c);
      Constr.Weight = WeightFixed;
      Geo::VarTerm<ScalarType> Var(IndexV, c, 1);
      Constr.LHSide.push_back(Var);
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

// template <class ScalarType>
// void SmoothLaplacianImplicit(std::vector<ScalarType> &Signal,
//                              const std::vector<std::vector<int>> &Faces,
//                              const std::vector<int> &FixedVert,
//                              const std::vector<std::pair<int, int>> &Features
//                              =
//                                  std::vector<std::pair<int, int>>()) {

//   ScalarType WeightFixed = 1000000;
//   ScalarType WeightLap = 1;

//   std::vector<std::vector<int>> KernelV;
//   GetLaplacianNeighbours(VertPos, Faces, KernelV, Features);

//   std::vector<Geo::LinearConstraint<ScalarType>> LinConstr;

//   std::set<int> FixedSet(FixedVert.begin(), FixedVert.end());
//   for (size_t i = 0; i < KernelV.size(); i++) {
//     int IndexV0 = i;
//     if (FixedSet.count(IndexV0) > 0)
//       continue;
//     for (size_t c = 0; c < 3; c++) {
//       Geo::LinearConstraint<ScalarType> Constr;
//       Constr.Weight = WeightLap;
//       Constr.RHSide = 0;
//       Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV0, c, -1));
//       ScalarType OtherW = (ScalarType)1 / (ScalarType)KernelV[i].size();

//       for (size_t j = 0; j < KernelV[i].size(); j++) {
//         int IndexV1 = KernelV[i][j];

//         assert(IndexV0 != IndexV1);

//         Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1, c,
//         OtherW));
//       }
//       LinConstr.push_back(Constr);
//     }
//   }

//   for (size_t i = 0; i < FixedVert.size(); i++) {
//     int IndexV = FixedVert[i];
//     for (size_t c = 0; c < 3; c++) {
//       Geo::LinearConstraint<ScalarType> Constr;
//       Constr.RHSide = VertPos[IndexV].V(c);
//       Constr.Weight = WeightFixed;
//       Geo::VarTerm<ScalarType> Var(IndexV, c, 1);
//       Constr.LHSide.push_back(Var);
//       LinConstr.push_back(Constr);
//     }
//   }

//   Geo::LinearConstraintLeastSQSolver<ScalarType> LinSqSolver;
//   LinSqSolver.Init(LinConstr);
//   std::vector<ScalarType> Solution;
//   LinSqSolver.Solve(Solution);

//   for (size_t i = 0; i < VertPos.size(); i++) {
//     int IndexX = (i * 3);
//     int IndexY = IndexX + 1;
//     int IndexZ = IndexX + 2;
//     ScalarType ValX = Solution[IndexX];
//     ScalarType ValY = Solution[IndexY];
//     ScalarType ValZ = Solution[IndexZ];
//     VertPos[i] = Geo::Point3<ScalarType>(ValX, ValY, ValZ);
//   }
// }

template <class ScalarType>
void SmoothLaplacianImplicit(std::vector<Geo::Point2<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &FixedVert) {
  std::vector<Geo::Point3<ScalarType>> VertPos3D;
  for (size_t i = 0; i < VertPos.size(); i++)
    VertPos3D.push_back(Geo::Point3<ScalarType>(VertPos[i].X, VertPos[i].Y, 0));

  SmoothLaplacianImplicit(VertPos3D, Faces, FixedVert);

  for (size_t i = 0; i < VertPos3D.size(); i++) {
    VertPos[i].X = VertPos3D[i].X;
    VertPos[i].Y = VertPos3D[i].Y;
  }
}

template <class ScalarType>
void SmoothLaplacianImplicit(std::vector<ScalarType> &Signal,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &FixedVert) {
  std::vector<Geo::Point3<ScalarType>> VertPos3D;
  for (size_t i = 0; i < Signal.size(); i++)
    VertPos3D.push_back(Geo::Point3<ScalarType>(Signal[i], 0, 0));

  SmoothLaplacianImplicit(VertPos3D, Faces, FixedVert);

  for (size_t i = 0; i < VertPos3D.size(); i++)
    Signal[i] = VertPos3D[i].X;
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
  // Geo::SmoothLaplacian<ScalarType>(Pos3D, Faces, 5, FixedV, Features,
  // false);
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
