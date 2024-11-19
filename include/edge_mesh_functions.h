#ifndef EDGE_MESH_FUNCTIONS
#define EDGE_MESH_FUNCTIONS

#include "connectivity.h"
#include "mesh_functions.h"
#include "point3.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <plane3.h>
#include <pos3.h>
#include <set>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// #include "Space_Query/triangle_grid_2D.h"

namespace Geo {
template <class ScalarType> class EdgeMeshFunctions {

  //   // Function to traverse from a start node until a branch point or end
  //   static void
  //   traverseSequence(int start,
  //                    const std::vector<std::vector<int> > &adjList,
  //                    std::unordered_set<int> &visited,
  //                    std::vector<int> &sequence,
  //                    const std::unordered_set<int> &branchPoints) {
  //     int current = start;
  //     int prev = -1; // To keep track of previous node and avoid cycles
  //     sequence.clear();
  //     while (true) {
  //       sequence.push_back(current);
  //       visited.insert(current);

  //       // Check if current node is a branch point (end of sequence)
  //       if (branchPoints.find(current) != branchPoints.end() &&
  //           current != start) {
  //         break;
  //       }

  //       // Move to the next unvisited neighbor
  //       bool foundNext = false;
  //     //   for (int neighbor : adjList.at(current)) {
  //     //     if (neighbor != prev && visited.find(neighbor) == visited.end())
  //     {
  //      for (int i=0;i<adjList[current].size();i++) {
  //         int neighbor=adjList[current][i];
  //         if (neighbor != prev && visited.find(neighbor) == visited.end()) {
  //           prev = current;
  //           current = neighbor;
  //           foundNext = true;
  //           break;
  //         }
  //       }

  //       if (!foundNext) { // Reached an endpoint (no further unvisited
  //       neighbors)
  //         break;
  //       }
  //     }
  //   }

public:
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
                      const std::vector<std::vector<int>> &Faces,
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

  // Function to find all sorted sequences based on branch points
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
      assert(AdjList[startV].size() <=2); 
      // no disconnected vertices
      assert(AdjList[startV].size() > 0); 
      
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
            if (neighbor != prev){//} && !visited.count(neighbor)) {
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

      std::vector<int> loop;
      int current = startV;
      int prev = -1;

      // Traverse the loop
      do {
        loop.push_back(current);
        visited.insert(current);

        // Find the next vertex in the loop
        int next = -1;
        assert(AdjList[current].size()==2);
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
      //then add back the first one 
      loop.push_back(startV);
      //at least 3 elements and the first repeated equal to 4
      assert(loop.size()>=4);
      subsequences1.push_back(loop);
    }

    //report to the original sequence
    std::map<Geo::Point3<ScalarType>, int> OrigMap;
    for (size_t i=0;i<VertPos.size();i++)
        OrigMap[VertPos[i]]=i;
    
    //then update the sequences
    subsequences.clear();
    for (size_t i=0;i<subsequences1.size();i++)
    {
        subsequences.resize(subsequences.size()+1);
        for (size_t j=0;j<subsequences1[i].size();j++)
        {
            int IndexV1=subsequences1[i][j];
            Geo::Point3<ScalarType> Pos1=VertPos1[IndexV1];
            assert(OrigMap.count(Pos1)>0);
            subsequences.back().push_back(OrigMap[Pos1]);
        }
    }
  }
  
};

}; // namespace Geo
#endif
