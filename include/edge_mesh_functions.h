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

#ifndef EDGE_MESH_FUNCTIONS
#define EDGE_MESH_FUNCTIONS

#include "connectivity.h"
#include "fitting_points.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "point3.h"
#include "polyline3.h"
#include "segment3.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <plane3.h>
#include <pos3.h>
#include <set>
#include <string>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <unordered_set>
#include <utility>
#include <vector>
// #include "Space_Query/triangle_grid_2D.h"

namespace Geo {
template <class ScalarType> class EdgeMeshFunctions {

  static ScalarType Curvature(const std::vector<Point3<ScalarType>> &VertPos) {
    ScalarType Curv = 0;

    Geo::Plane3<ScalarType> FitPlane;
    FittingPoint<ScalarType>::FitPlaneToPointSet(VertPos, FitPlane);
    Geo::Point3<ScalarType> up = FitPlane.Direction();

    for (int j = 1; j < VertPos.size() - 1; j++) {
      int IndexV0 = j - 1;
      int IndexV1 = j;
      int IndexV2 = j + 1;

      Geo::Point3<ScalarType> Dir0 = (VertPos[IndexV1] - VertPos[IndexV0]);
      Geo::Point3<ScalarType> Dir1 = (VertPos[IndexV2] - VertPos[IndexV1]);
      Dir0.Normalize();
      Dir1.Normalize();

      Curv += AngleDeg(Dir0, Dir1);
    }
    return fabs(Curv);
  }

public:
  static std::vector<int>
  getBranchIndexes(const std::vector<std::pair<int, int>> &Edges) {
    std::vector<std::vector<int>> AdjList;
    getVVAdjacency(Edges, AdjList);

    // get places to split
    std::vector<int> SplitP;
    for (size_t i = 0; i < AdjList.size(); i++) {
      if (AdjList[i].size() <= 2)
        continue;
      SplitP.push_back(i);
    }
    return SplitP;
  }

  static bool IsProperEdgeM(const std::vector<std::pair<int, int>> &Edges) {
    for (size_t i = 0; i < Edges.size(); i++)
      if (Edges[i].first == Edges[i].second)
        return false;
    return true;
  }

  static void SplitAlong(const std::vector<int> &SplitVert,
                         std::vector<Geo::Point3<ScalarType>> &EdgeVertPos,
                         std::vector<std::pair<int, int>> &Edges) {
    std::set<int> SplitSet(SplitVert.begin(), SplitVert.end());
    std::vector<bool> Linked(EdgeVertPos.size(), false);
    // std::cout<<"Num V0 "<<EdgeVertPos.size()<<std::endl;
    for (size_t i = 0; i < Edges.size(); i++) {
      int IndexV0 = Edges[i].first;
      int IndexV1 = Edges[i].second;

      assert(IndexV0 >= 0);
      assert(IndexV0 < Linked.size());

      assert(IndexV1 >= 0);
      assert(IndexV1 < Linked.size());

      if (SplitSet.count(IndexV0) > 0) {
        if (!Linked[IndexV0])
          Linked[IndexV0] = true;
        else {
          // duplicate if needed
          EdgeVertPos.push_back(EdgeVertPos[IndexV0]);
          Edges[i].first = EdgeVertPos.size() - 1;
        }
      }
      if (SplitSet.count(IndexV1) > 0) {
        if (!Linked[IndexV1])
          Linked[IndexV1] = true;
        else {
          // duplicate if needed
          EdgeVertPos.push_back(EdgeVertPos[IndexV1]);
          Edges[i].second = EdgeVertPos.size() - 1;
        }
      }
    }

    // std::cout<<"Num V1 "<<EdgeVertPos.size()<<std::endl;
  }

  static void
  ExtractFromMeshBoundary(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          std::vector<Geo::Point3<ScalarType>> &EdgeVertPos,
                          std::vector<std::pair<int, int>> &Edges) {
    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
    EdgeVertPos.clear();
    Edges.clear();
    std::vector<std::vector<int>> BorderSeq;
    getBorderSequncesVert(Faces, BorderSeq);
    std::map<Geo::Point3<ScalarType>, int> PosIDX;
    std::set<PosPair> AddedPos;
    for (size_t i = 0; i < BorderSeq.size(); i++) {
      int sizeB = BorderSeq[i].size();
      for (size_t j = 0; j < BorderSeq[i].size(); j++) {
        int IndexV0 = BorderSeq[i][j];
        int IndexV1 = BorderSeq[i][(j + 1) % sizeB];
        Geo::Point3<ScalarType> Pos0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1 = VertPos[IndexV1];
        int NewIndexV0 = -1;
        int NewIndexV1 = -1;
        if (PosIDX.count(Pos0) > 0)
          NewIndexV0 = PosIDX[Pos0];
        else {
          EdgeVertPos.push_back(Pos0);
          NewIndexV0 = EdgeVertPos.size() - 1;
        }

        if (PosIDX.count(Pos1) > 0)
          NewIndexV1 = PosIDX[Pos1];
        else {
          EdgeVertPos.push_back(Pos1);
          NewIndexV1 = EdgeVertPos.size() - 1;
        }
        assert(NewIndexV0 >= 0);
        assert(NewIndexV1 >= 0);

        PosPair Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        if (AddedPos.count(Key) > 0)
          continue;
        Edges.push_back(std::pair<int, int>(NewIndexV0, NewIndexV1));
        AddedPos.insert(Key);
      }
    }
  }

  static void
  ExtractFromFeatures(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      // const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::pair<int, int>> &Features,
                      std::vector<Geo::Point3<ScalarType>> &EdgeVertPos,
                      std::vector<std::pair<int, int>> &Edges) {
    EdgeVertPos.clear();
    Edges.clear();

    // make the features unique
    std::vector<std::pair<int, int>> UniqFFeatures;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      std::pair<int, int> Entry(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
      UniqFFeatures.push_back(Entry);
    }
    std::sort(UniqFFeatures.begin(), UniqFFeatures.end());
    UniqFFeatures.erase(std::unique(UniqFFeatures.begin(), UniqFFeatures.end()),
                        UniqFFeatures.end());

    std::vector<std::vector<int>> EdgeVect;
    for (size_t i = 0; i < UniqFFeatures.size(); i++) {
      int IndexV0 = UniqFFeatures[i].first;
      int IndexV1 = UniqFFeatures[i].second;
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      EdgeVect.resize(EdgeVect.size() + 1);
      EdgeVertPos.push_back(P0);
      EdgeVect.back().push_back(EdgeVertPos.size() - 1);
      EdgeVertPos.push_back(P1);
      EdgeVect.back().push_back(EdgeVertPos.size() - 1);
    }
    // std::cout<<"size V0:"<<EdgeVertPos.size()<<std::endl;
    RemoveDuplicatedVert(EdgeVertPos, EdgeVect);
    // std::cout<<"size V1:"<<EdgeVertPos.size()<<std::endl;
    RemoveUnfererencedVert(EdgeVertPos, EdgeVect);
    // std::cout<<"size V2:"<<EdgeVertPos.size()<<std::endl;
    Edges.resize(EdgeVect.size(), std::pair<int, int>(-1, -1));
    for (size_t i = 0; i < EdgeVect.size(); i++) {
      Edges[i].first = EdgeVect[i][0];
      Edges[i].second = EdgeVect[i][1];
    }
    // exit(0);
  }

  // // Function to find all sorted sequences based on branch points
  // static void ExtractSortedNodesSequences(
  //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const std::vector<std::pair<int, int>> &Edges,
  //     std::vector<std::vector<int>> &subsequences) {

  //   assert(IsProperEdgeM(Edges));
  //   // Build adjacency list
  //   std::vector<std::vector<int>> AdjList;
  //   getVVAdjacency(Edges, AdjList);

  //   // get places to split
  //   std::vector<int> SplitP;
  //   for (size_t i = 0; i < AdjList.size(); i++) {
  //     if (AdjList[i].size() <= 2)
  //       continue;
  //     SplitP.push_back(i);
  //   }

  //   // split the mesh in multiple components
  //   std::vector<Geo::Point3<ScalarType>> VertPos1 = VertPos;
  //   std::vector<std::pair<int, int>> Edges1 = Edges;
  //   SplitAlong(SplitP, VertPos1, Edges1);
  //   assert(IsProperEdgeM(Edges));

  //   // get new VV adj
  //   getVVAdjacency(Edges1, AdjList);

  //   std::vector<std::vector<int>> subsequences1;
  //   std::unordered_set<int> visited;

  //   // IDENTIFY SEGMENTS
  //   for (size_t startV = 0; startV < AdjList.size(); startV++) {
  //     if (visited.count(startV) > 0)
  //       continue;
  //     // should be already split at this point
  //     assert(AdjList[startV].size() <= 2);
  //     // no disconnected vertices
  //     // assert(AdjList[startV].size() > 0);
  //     if (AdjList[startV].size() == 0)
  //       continue;

  //     //  Start a new sequence if valence == 1
  //     if (AdjList[startV].size() == 1) {
  //       std::vector<int> sequence;
  //       int current = startV;
  //       int prev = -1;

  //       // Traverse the path starting from `current`
  //       while (true) {
  //         sequence.push_back(current);
  //         visited.insert(current);

  //         // Find the next vertex in the sequence
  //         int next = -1;
  //         for (int j = 0; j < AdjList[current].size(); j++) {
  //           int neighbor = AdjList[current][j];
  //           if (neighbor != prev) { //} && !visited.count(neighbor)) {
  //             next = neighbor;
  //             break;
  //           }
  //         }

  //         if (next == -1 || AdjList[next].size() != 2) {
  //           if (next != -1) {
  //             sequence.push_back(next);
  //             visited.insert(next);
  //           }
  //           break; // End of the sequence
  //         }

  //         prev = current;
  //         current = next;
  //       }
  //       subsequences1.push_back(sequence);
  //     }
  //   }

  //   // IDENTIFY LOOPS
  //   for (size_t startV = 0; startV < AdjList.size(); startV++) {
  //     if (visited.count(startV) > 0)
  //       continue;

  //     if (AdjList[startV].size() == 0)
  //       continue;

  //     std::vector<int> loop;
  //     int current = startV;
  //     int prev = -1;

  //     // Traverse the loop
  //     do {
  //       loop.push_back(current);
  //       visited.insert(current);

  //       // Find the next vertex in the loop
  //       int next = -1;
  //       assert(AdjList[current].size() == 2);
  //       for (int j = 0; j < AdjList[current].size(); j++) {
  //         int neighbor = AdjList[current][j];
  //         if (neighbor != prev) {
  //           next = neighbor;
  //           break;
  //         }
  //       }
  //       prev = current;
  //       current = next;
  //     } while (current != startV);
  //     // then add back the first one
  //     loop.push_back(startV);
  //     // at least 3 elements and the first repeated equal to 4
  //     assert(loop.size() >= 4);
  //     subsequences1.push_back(loop);
  //   }

  //   // report to the original sequence
  //   std::map<Geo::Point3<ScalarType>, int> OrigMap;
  //   for (size_t i = 0; i < VertPos.size(); i++)
  //     OrigMap[VertPos[i]] = i;

  //   // then update the sequences
  //   subsequences.clear();
  //   for (size_t i = 0; i < subsequences1.size(); i++) {
  //     subsequences.resize(subsequences.size() + 1);
  //     for (size_t j = 0; j < subsequences1[i].size(); j++) {
  //       int IndexV1 = subsequences1[i][j];
  //       Geo::Point3<ScalarType> Pos1 = VertPos1[IndexV1];
  //       assert(OrigMap.count(Pos1) > 0);
  //       subsequences.back().push_back(OrigMap[Pos1]);
  //     }
  //   }
  // }

  static void MakeUnique(std::vector<std::pair<int, int>> &Edges) {
    for (size_t i = 0; i < Edges.size(); i++) {
      if (Edges[i].first < Edges[i].second)
        continue;
      std::swap(Edges[i].first, Edges[i].second);
    }
    std::sort(Edges.begin(), Edges.end());
    Edges.erase(std::unique(Edges.begin(), Edges.end()), Edges.end());
  }

  // Function to find all sorted sequences based on branch points
  // this function require NO unreferenced vertices
  static void ExtractSortedNodesSequences(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::pair<int, int>> &Edges,
      std::vector<std::vector<int>> &subsequences) {

    assert(IsProperEdgeM(Edges));
    // Build adjacency list
    std::vector<std::vector<int>> AdjList;
    getVVAdjacency(Edges, AdjList);

    // get places to split
    std::vector<int> SplitP;
    for (size_t i = 0; i < AdjList.size(); i++) {
      if (AdjList[i].size() <= 2)
        continue;
      SplitP.push_back(i);
    }

    // split the mesh in multiple components
    std::vector<Geo::Point3<ScalarType>> VertPos1 = VertPos;
    std::vector<std::pair<int, int>> Edges1 = Edges;
    SplitAlong(SplitP, VertPos1, Edges1);
    assert(IsProperEdgeM(Edges));

    // get new VV adj
    getVVAdjacency(Edges1, AdjList);

    std::vector<std::vector<int>> subsequences1;
    std::unordered_set<int> visited;

    // IDENTIFY SEGMENTS
    for (size_t startV = 0; startV < AdjList.size(); startV++) {
      if (visited.count(startV) > 0)
        continue;
      // should be already split at this point
      assert(AdjList[startV].size() <= 2);
      // no disconnected vertices
      // assert(AdjList[startV].size() > 0);
      if (AdjList[startV].size() == 0)
        continue;

      //  Start a new sequence if valence == 1
      if (AdjList[startV].size() == 1) {
        std::vector<int> sequence;
        int current = startV;
        int prev = -1;

        // Traverse the path starting from `current`
        while (true) {
          sequence.push_back(current);
          visited.insert(current);

          // Find the next vertex in the sequence
          int next = -1;
          for (int j = 0; j < AdjList[current].size(); j++) {
            int neighbor = AdjList[current][j];
            if (neighbor != prev) { //} && !visited.count(neighbor)) {
              next = neighbor;
              break;
            }
          }

          if (next == -1 || AdjList[next].size() != 2) {
            if (next != -1) {
              sequence.push_back(next);
              visited.insert(next);
            }
            break; // End of the sequence
          }

          prev = current;
          current = next;
        }
        subsequences1.push_back(sequence);
      }
    }

    // IDENTIFY LOOPS
    for (size_t startV = 0; startV < AdjList.size(); startV++) {
      if (visited.count(startV) > 0)
        continue;

      if (AdjList[startV].size() == 0)
        continue;

      std::vector<int> loop;
      int current = startV;
      int prev = -1;

      // Traverse the loop
      do {
        loop.push_back(current);
        visited.insert(current);

        // Find the next vertex in the loop
        int next = -1;
        assert(AdjList[current].size() == 2);
        for (int j = 0; j < AdjList[current].size(); j++) {
          int neighbor = AdjList[current][j];
          if (neighbor != prev) {
            next = neighbor;
            break;
          }
        }
        prev = current;
        current = next;
      } while (current != startV);
      // then add back the first one
      loop.push_back(startV);
      // at least 3 elements and the first repeated equal to 4
      assert(loop.size() >= 4);
      subsequences1.push_back(loop);
    }

    // report to the original sequence
    std::map<Geo::Point3<ScalarType>, int> OrigMap;
    for (size_t i = 0; i < VertPos.size(); i++)
      OrigMap[VertPos[i]] = i;

    // then update the sequences
    subsequences.clear();
    for (size_t i = 0; i < subsequences1.size(); i++) {
      subsequences.resize(subsequences.size() + 1);
      for (size_t j = 0; j < subsequences1[i].size(); j++) {
        int IndexV1 = subsequences1[i][j];
        Geo::Point3<ScalarType> Pos1 = VertPos1[IndexV1];
        assert(OrigMap.count(Pos1) > 0);
        subsequences.back().push_back(OrigMap[Pos1]);
      }
    }
  }

  // Function to find all sorted sequences based on branch points
  // this function works with unreferenced vertices
  static void ExtractSortedNodesSequencesMap(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::pair<int, int>> &Edges,
      std::vector<std::vector<int>> &subsequences) {

    // get the submesh
    std::vector<Geo::Point3<ScalarType>> SelVertPos;
    std::vector<std::pair<int, int>> SelEdges;
    Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
        VertPos, Edges, SelVertPos, SelEdges);

    // get the vert-to-vert map
    std::vector<int> VertMap;
    Geo::FindVertCorrespondences(SelVertPos, VertPos, VertMap);
    assert(VertMap.size() == SelVertPos.size());

    // get the sequences
    subsequences.clear();
    ExtractSortedNodesSequences(SelVertPos, SelEdges, subsequences);

    // then remap back
    for (size_t i = 0; i < subsequences.size(); i++)
      for (size_t j = 0; j < subsequences[i].size(); j++) {
        int IndexV = subsequences[i][j];
        assert(IndexV >= 0);
        assert(IndexV < SelVertPos.size());
        assert(IndexV < VertMap.size());
        assert(VertMap[IndexV] >= 0);
        assert(VertMap[IndexV] < VertPos.size());
        subsequences[i][j] = VertMap[subsequences[i][j]];
      }
  }

  // resample positions on the edge mesh using branch points.
  static void
  ResampleNodeSequences(std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::pair<int, int>> &Edges) {
    // first get the sequences
    std::vector<std::vector<int>> subsequences;
    ExtractSortedNodesSequencesMap(VertPos, Edges, subsequences);
    // resample each sequence
    for (size_t i = 0; i < subsequences.size(); i++) {
      // can't be a single point but at least an edge
      assert(subsequences[i].size() >= 2);

      // in this case no need to resample
      if (subsequences[i].size() == 2)
        continue;

      // get current polyline
      Geo::Polyline3<ScalarType> CurrPoly;
      for (size_t j = 0; j < subsequences[i].size(); j++) {
        Geo::Point3<ScalarType> CurrPos = VertPos[subsequences[i][j]];
        CurrPoly.PolyPos.push_back(CurrPos);
      }
      int num_points = subsequences[i].size();

      // resample
      std::vector<Geo::Point3<ScalarType>> ResampledPoly;
      int numSamplePoints = subsequences[i].size();
      assert(numSamplePoints > 1);
      ResampledPoly = CurrPoly.SampleUniformly(numSamplePoints);
      // std::cout<<"Resampling target "<<numSamplePoints
      //          <<" points to "<<ResampledPoly.size()<<" points."<<std::endl;

      // update positions (not considering the last one to avoid numerical
      // issues)
      assert(ResampledPoly.size() == subsequences[i].size() - 2);
      for (size_t j = 0; j < ResampledPoly.size(); j++)
        VertPos[subsequences[i][j + 1]] = ResampledPoly[j];
    }
  }

  static bool SampleByCurvature(const std::vector<Point3<ScalarType>> &Points,
                                const int n_internal,
                                std::vector<int> &Selected) {
    if (n_internal <= 0)
      return false;

    bool IsCircle = (Points[0] == Points.back());

    ScalarType CurvValue = Curvature(Points);
    ScalarType CurvInt = CurvValue / n_internal;
    if (!IsCircle) {
      CurvInt = CurvValue / (n_internal + 1);
      CurvInt *= 1.01;
    }

    ScalarType CurrA = 0;
    Selected.clear();

    Geo::Plane3<ScalarType> FitPlane;
    FittingPoint<ScalarType>::FitPlaneToPointSet(Points, FitPlane);
    Geo::Point3<ScalarType> up = FitPlane.Direction();

    for (size_t i = 1; i < (Points.size() - 1); i++) {

      int IndexV0 = i - 1;
      int IndexV1 = i;
      int IndexV2 = i + 1;

      Geo::Point3<ScalarType> Dir0 = (Points[IndexV1] - Points[IndexV0]);
      Geo::Point3<ScalarType> Dir1 = (Points[IndexV2] - Points[IndexV1]);
      Dir0.Normalize();
      Dir1.Normalize();

      // ScalarType dot = Dot(Dir0, Dir1);
      // Point3<ScalarType> cross = Cross(Dir0, Dir1);
      // ScalarType sin_theta = cross.Norm();
      // ScalarType angle = fabs(std::atan2(sin_theta, dot));
      // angle *= (ScalarType)180.0 / M_PI;
      // ScalarType sign = (Dot(cross, up) >= 0) ? 1 : -1;
      // CurrA += (sign * angle);
      CurrA += AngleDeg(Dir0, Dir1);
      if (fabs(CurrA) >= CurvInt) {
        CurrA = 0;
        Selected.push_back(IndexV1);
      }
    }

    // add the last one in case is a circle
    if (IsCircle && (Selected.size() < n_internal))
      Selected.push_back(0);
  }

  static void
  SubSampleRegularCurvature(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::pair<int, int>> &Edges,
                            const std::vector<int> &Corners,
                            std::vector<int> &Selected,
                            const ScalarType &CurvatureIntervalCircle = 80,
                            const ScalarType &CurvatureIntervalPoly = 40) {
    // save indexes
    std::map<Geo::Point3<ScalarType>, int> PosToIDX;
    for (size_t i = 0; i < VertPos.size(); i++)
      PosToIDX[VertPos[i]] = i;

    // split along existing ones
    std::vector<Geo::Point3<ScalarType>> VertPos1 = VertPos;
    std::vector<std::pair<int, int>> Edges1 = Edges;
    SplitAlong(Corners, VertPos1, Edges1);

    // then get the sequences
    std::vector<std::vector<int>> subsequences;
    ExtractSortedNodesSequences(VertPos1, Edges1, subsequences);

    // then compute the curvature for each piece
    std::vector<Geo::Point3<ScalarType>> SelPos;
    for (size_t i = 0; i < subsequences.size(); i++) {
      // get the current point sequences
      std::vector<Point3<ScalarType>> CurrPoly;
      for (size_t j = 0; j < subsequences[i].size(); j++)
        CurrPoly.push_back(VertPos1[subsequences[i][j]]);

      ScalarType currCurv = Curvature(CurrPoly);
      bool isCycle = (CurrPoly[0] == CurrPoly.back());
      // if (isCycle)
      //   std::cout << "Cycle Curv:" << currCurv << std::endl;

      int InterValNum = 0;
      if (isCycle)
        InterValNum = floor(0.5 + (currCurv / CurvatureIntervalCircle));
      else
        InterValNum = floor(0.5 + (currCurv / CurvatureIntervalPoly));

      if (InterValNum == 0)
        continue;

      std::vector<int> Selected;
      SampleByCurvature(CurrPoly, InterValNum, Selected);
      // if (isCycle)
      //   std::cout << "Sel Pos:" << Selected.size() << std::endl;

      for (size_t j = 0; j < Selected.size(); j++)
        SelPos.push_back(CurrPoly[Selected[j]]);
    }

    // retrieve pos
    Selected.clear();
    for (size_t i = 0; i < SelPos.size(); i++) {
      assert(PosToIDX.count(SelPos[i]) > 0);
      Selected.push_back(PosToIDX[SelPos[i]]);
    }
  }

  static void ExtractEdgeMeshFromVertPairs(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::pair<int, int>> &Edges,
      std::vector<Geo::Point3<ScalarType>> &EdgeMeshVertPos,
      std::vector<std::vector<int>> &ConnectivityEdges) {

    // get the ones touched
    std::vector<bool> Used(VertPos.size(), false);
    for (size_t i = 0; i < Edges.size(); i++) {
      assert(Edges[i].first >= 0);
      assert(Edges[i].second >= 0);
      assert(Edges[i].first < VertPos.size());
      assert(Edges[i].second < VertPos.size());

      Used[Edges[i].first] = true;
      Used[Edges[i].second] = true;
    }

    std::map<int, int> VertMap;
    EdgeMeshVertPos.clear();
    for (size_t i = 0; i < VertPos.size(); i++) {
      if (!Used[i])
        continue;
      EdgeMeshVertPos.push_back(VertPos[i]);
      VertMap[i] = EdgeMeshVertPos.size() - 1;
    }

    ConnectivityEdges.clear();
    for (size_t i = 0; i < Edges.size(); i++) {
      int IndexV0 = Edges[i].first;
      int IndexV1 = Edges[i].second;
      assert(VertMap.count(IndexV0) > 0);
      assert(VertMap.count(IndexV1) > 0);
      IndexV0 = VertMap[IndexV0];
      IndexV1 = VertMap[IndexV1];
      assert(IndexV0 >= 0);
      assert(IndexV0 < EdgeMeshVertPos.size());
      assert(IndexV1 >= 0);
      assert(IndexV1 < EdgeMeshVertPos.size());
      ConnectivityEdges.resize(ConnectivityEdges.size() + 1);
      ConnectivityEdges.back().push_back(IndexV0);
      ConnectivityEdges.back().push_back(IndexV1);
    }
  }

  static void ExtractEdgeMeshFromVertPairs(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<std::pair<int, int>>> &Edges,
      std::vector<Geo::Point3<ScalarType>> &EdgeMeshVertPos,
      std::vector<std::vector<int>> &ConnectivityEdges)
  {
    EdgeMeshVertPos.clear();
    ConnectivityEdges.clear();
    for (size_t i = 0; i < Edges.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> SubEdgeMeshVertPos;
      std::vector<std::vector<int>> SubConnectivityEdges;
      ExtractEdgeMeshFromVertPairs(VertPos, Faces, Edges[i],
                                   SubEdgeMeshVertPos, SubConnectivityEdges);
      AppendMesh<ScalarType>(EdgeMeshVertPos, ConnectivityEdges, SubEdgeMeshVertPos, SubConnectivityEdges);
    }
  }

  static void ExtractEdgeMeshFromPosPairs(
      const std::vector<std::pair<Geo::Point3<ScalarType>,
                                  Geo::Point3<ScalarType>>> &EdgesPos,
      std::vector<Geo::Point3<ScalarType>> &EdgeMeshVertPos,
      std::vector<std::vector<int>> &ConnectivityEdges) {

    std::vector<std::pair<int, int>> Edges;
    std::map<Geo::Point3<ScalarType>, int> VertIdx;
    EdgeMeshVertPos.clear();
    ConnectivityEdges.clear();
    for (size_t i = 0; i < EdgesPos.size(); i++) {
      Geo::Point3<ScalarType> P0 = EdgesPos[i].first;
      if (VertIdx.count(P0) == 0) {
        EdgeMeshVertPos.push_back(P0);
        VertIdx[P0] = EdgeMeshVertPos.size() - 1;
      }
      Geo::Point3<ScalarType> P1 = EdgesPos[i].second;
      if (VertIdx.count(P1) == 0) {
        EdgeMeshVertPos.push_back(P1);
        VertIdx[P1] = EdgeMeshVertPos.size() - 1;
      }
      assert(VertIdx.count(P0) > 0);
      assert(VertIdx.count(P1) > 0);
      int IndexV0 = VertIdx[P0];
      int IndexV1 = VertIdx[P1];
      ConnectivityEdges.resize(ConnectivityEdges.size() + 1);
      ConnectivityEdges.back().push_back(IndexV0);
      ConnectivityEdges.back().push_back(IndexV1);
    }
  }

  static void ExtractEdgeMeshFromIdxSequences(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Seuquences,
      std::vector<Geo::Point3<ScalarType>> &EdgeMeshVertPos,
      std::vector<std::vector<int>> &ConnectivityEdges) {
    std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        EdgesPos;
    for (size_t i = 0; i < Seuquences.size(); i++) {
      for (size_t j = 0; j < Seuquences[i].size() - 1; j++) {
        int IndexV0 = Seuquences[i][j];
        int IndexV1 = Seuquences[i][j + 1];
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        EdgesPos.push_back(
            std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>(P0,
                                                                        P1));
      }
    }
    ExtractEdgeMeshFromPosPairs(EdgesPos, EdgeMeshVertPos, ConnectivityEdges);
  }

  static void ExtractEdgeMeshFromPolylines(
      const std::vector<Polyline3<ScalarType>> &Polylines,
      std::vector<Geo::Point3<ScalarType>> &EdgeMeshVertPos,
      std::vector<std::vector<int>> &ConnectivityEdges) {

    std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> EdgesPos;
    for (size_t i = 0; i < Polylines.size(); i++) {
      for (size_t j = 0; j < Polylines[i].PolyPos.size() - 1; j++) {
        Geo::Point3<ScalarType> P0 = Polylines[i].PolyPos[j];
        Geo::Point3<ScalarType> P1 = Polylines[i].PolyPos[j + 1];
        EdgesPos.push_back(std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>(P0, P1));
      }
    }
    ExtractEdgeMeshFromPosPairs(EdgesPos, EdgeMeshVertPos, ConnectivityEdges);
  }

  static void ExtractEdgePairFromIdxSequences(
      const std::vector<std::vector<int>> &Sequences,
      std::vector<std::vector<std::pair<int, int>>> &ConnectivityEdges) {
    ConnectivityEdges.clear();
    ConnectivityEdges.resize(Sequences.size());

    for (size_t i = 0; i < Sequences.size(); i++) {
      for (size_t j = 0; j < Sequences[i].size() - 1; j++) {
        int IndexV0 = Sequences[i][j];
        int IndexV1 = Sequences[i][j + 1];
        ConnectivityEdges[i].push_back(std::pair<int, int>(IndexV0, IndexV1));
      }
    }
  }
};
}; // namespace Geo
#endif
