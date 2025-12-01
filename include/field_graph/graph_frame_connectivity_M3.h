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

#ifndef FRAME_FIELD_GRAPH_CONNECTIVITY_M3
#define FRAME_FIELD_GRAPH_CONNECTIVITY_M3

#include "../connectivity.h"
#include "../geodesic.h"
#include "../mesh_functions.h"
#include "graph.h"
#include "graph_frame_field.h"
#include "graph_frame_field_sheet_solve_SAT_M3.h"
#include "label_graph_functions.h"
#include "smooth_mesh.h"
#include "smooth_sheet_surfaces.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>
#include <chrono>

template <class ScalarType> class LabelGraph_ConnectivityM3 {
public:
  typedef enum {
    MANIFOLD,
    NON_MANIFOLD,
    INCORRECT_EXIT,
    UNDETERMINED
  } ManifoldType;

  ScalarType TimeSolvingSat = 0;
  ScalarType TimeSettingConstraints = 0;
  ScalarType TimeCollectingResults = 0;
  ScalarType TimeDeactivatingNodes = 0;

private:
  struct SubGraph {
    std::vector<int> NodesIdx;
    int IndexV;
    std::vector<std::vector<int>> Neigh;
    std::vector<std::vector<ScalarType>> Weight;
    std::vector<std::vector<bool>> Active;
    std::vector<bool> IsExit;
    std::vector<bool> IsDoubleExit;

    int LocalIndex(const int IndexN) {
      int LocIdx = -1;
      // return the index of IndexB in NodesIdx
      for (size_t i = 0; i < NodesIdx.size(); i++) {
        if (NodesIdx[i] == IndexN) {
          LocIdx = i;
          break;
        }
      }
      return LocIdx;
    }

    void LocalToGlobalIndex(const std::vector<int> LocalIdx,
                            std::vector<int> &GlobalIdx) const {
      GlobalIdx.clear();
      for (size_t i = 0; i < LocalIdx.size(); i++) {
        int LIdx = LocalIdx[i];
        assert(LIdx >= 0);
        assert(LIdx < NodesIdx.size());
        GlobalIdx.push_back(NodesIdx[LIdx]);
      }
    }

    // // ------------------------------------------------------------
    // // Find a cycle passing through `start`
    // // Returns true if such a cycle exists, and fills `cycleNodes`
    // // ------------------------------------------------------------
    // bool FindCycleThroughNode(int start, std::vector<int> &cycleNodes) {
    //   int N = Neigh.size();

    //   // List of active neighbors of start
    //   std::vector<int> neighStart;
    //   for (int k = 0; k < (int)Neigh[start].size(); ++k)
    //     if (Active[start][k])
    //       neighStart.push_back(Neigh[start][k]);

    //   // We need at least 2 neighbors to form a cycle
    //   if (neighStart.size() < 2)
    //     return false;

    //   // Try all ordered pairs (a, b)
    //   for (int a : neighStart) {
    //     for (int b : neighStart) {
    //       if (a == b)
    //         continue;

    //       // BFS from a to b avoiding edge (start, a)
    //       std::vector<int> parent(N, -1);
    //       std::vector<bool> visited(N, false);
    //       std::queue<int> q;

    //       visited[a] = true;
    //       q.push(a);

    //       while (!q.empty()) {
    //         int u = q.front();
    //         q.pop();

    //         if (u == b)
    //           break;

    //         for (int k = 0; k < (int)Neigh[u].size(); ++k) {
    //           int v = Neigh[u][k];
    //           if (!Active[u][k])
    //             continue;

    //           // avoid using the edge (start,a)
    //           if ((u == start && v == a) || (u == a && v == start))
    //             continue;

    //           if (!visited[v]) {
    //             visited[v] = true;
    //             parent[v] = u;
    //             q.push(v);
    //           }
    //         }
    //       }

    //       // If we reached b, we found a cycle through start
    //       if (!visited[b])
    //         continue;

    //       // Reconstruct path a → ... → b
    //       std::vector<int> path;
    //       for (int cur = b; cur != -1; cur = parent[cur])
    //         path.push_back(cur);
    //       std::reverse(path.begin(), path.end());

    //       // Complete cycle: start → a → ... → b → start
    //       cycleNodes.clear();
    //       cycleNodes.push_back(start);
    //       cycleNodes.insert(cycleNodes.end(), path.begin(), path.end());
    //       cycleNodes.push_back(start);

    //       return true;
    //     }
    //   }

    //   return false;
    // }

    // bool FindCycleThroughGlobalNode(int start, std::vector<int> &cycleNodes)
    // {
    //   int localStart = LocalIndex(start);
    //   assert(localStart != -1);
    //   std::vector<int> localcycleNodes;
    //   bool found = FindCycleThroughNode(localStart, localcycleNodes);
    //   if (found) {
    //     LocalToGlobalIndex(localcycleNodes, cycleNodes);
    //   }
    //   return found;
    // }

    SubGraph(int _IndexV, const std::vector<int> &_NodesIdx) {
      IndexV = _IndexV;
      NodesIdx = _NodesIdx;
    }
  };

  // // A function to check if a node can reach two exit nodes
  // // and return true if it can as well as the list of nodes in the path
  // // consider also that the initial node can be exit node and should be
  // counted bool HasTwoDifferentExitPaths(const int &IndexN, const int &IndexV,
  //                               std::vector<int> &PathNodes) {
  //   PathNodes.clear();

  //   assert(IndexN >= 0);
  //   assert(IndexN < SubGraphIndex.size());
  //   // get the subgraph
  //   int IndexGraph = SubGraphIndex[IndexN][IndexV];
  //   assert(IndexGraph >= 0);
  //   assert(IndexGraph < SubGraphs.size());

  //   // then get the local index
  //   int LocalIdx = SubGraphs[IndexGraph].LocalIndex(IndexN);
  //   assert(LocalIdx >= 0);
  //   assert(LocalIdx < SubGraphs[IndexGraph].NodesIdx.size());

  //   // BFS to find exit nodes
  //   std::queue<int> Q;
  //   std::set<int> visited;
  //   Q.push(LocalIdx);
  //   visited.insert(LocalIdx);

  //   int exit_count = 0;

  //   while (!Q.empty()) {
  //     int curr = Q.front();
  //     Q.pop();

  //     // Check if current node is an exit
  //     if (SubGraphs[IndexGraph].IsExit[curr]) {
  //       exit_count++;
  //       if (exit_count >= 2) {
  //         // set the path nodes
  //         for (const int &v : visited)
  //           PathNodes.push_back(SubGraphs[IndexGraph].NodesIdx[v]);

  //         return true; // Found two different exits
  //       }
  //     }

  //     // Explore neighbors
  //     for (size_t n = 0; n < SubGraphs[IndexGraph].Neigh[curr].size(); n++) {
  //       if (!SubGraphs[IndexGraph].Active[curr][n])
  //         continue; // Skip inactive edges

  //       int neighbor = SubGraphs[IndexGraph].Neigh[curr][n];
  //       if (visited.count(neighbor) == 0) {
  //         visited.insert(neighbor);
  //         Q.push(neighbor);
  //       }
  //     }
  //   }
  //   return false; // Less than two exits found
  // }

  // // A function to check if there is a cycle around a node vertex
  // // considering that the connections are symmetric,
  // // reurns true if there is a cycle and the list of nodes in the cycle
  // bool HasCycleAround(const int &IndexN, const int &IndexV,
  //                     std::vector<int> &CycleNodes) {
  //   CycleNodes.clear();
  //   assert(IndexN >= 0);
  //   assert(IndexN < SubGraphIndex.size());
  //   // get the subgraph
  //   int IndexGraph = SubGraphIndex[IndexN][IndexV];
  //   assert(IndexGraph >= 0);
  //   assert(IndexGraph < SubGraphs.size());
  //   // then get the local index
  //   int LocalIdx = SubGraphs[IndexGraph].LocalIndex(IndexN);
  //   assert(LocalIdx >= 0);
  //   assert(LocalIdx < SubGraphs[IndexGraph].NodesIdx.size());
  //   // DFS to find cycles
  //   std::set<int> visited;
  //   std::set<int> recStack;
  //   std::function<bool(int)> dfs = [&](int node) {
  //     if (recStack.count(node))
  //       return true; // Cycle detected

  //     if (visited.count(node))
  //       return false; // Already visited

  //     visited.insert(node);
  //     recStack.insert(node);

  //     for (size_t n = 0; n < SubGraphs[IndexGraph].Neigh[node].size(); n++) {
  //       if (!SubGraphs[IndexGraph].Active[node][n])
  //         continue; // Skip inactive edges

  //       int neighbor = SubGraphs[IndexGraph].Neigh[node][n];
  //       if (dfs(neighbor))
  //         return true;
  //     }

  //     recStack.erase(node);
  //     return false;
  //   };
  //   if (dfs(LocalIdx)) {
  //     // If a cycle is detected, populate CycleNodes
  //     for (const int &n : recStack)
  //       CycleNodes.push_back(SubGraphs[IndexGraph].NodesIdx[n]);
  //     return true;
  //   }
  //   return false; // No cycle found
  // }

  // // return the manifold type around a node vertex also return all the nodes
  // //  and the corresponding vertices are satisfied by checking the paths
  // ManifoldType
  // CheckManifoldAround(const int &IndexN, const int &IndexV,
  //                     std::vector<std::pair<int, int>> &PathNodesVertices) {
  //   assert(IndexN >= 0);
  //   assert(IndexN < SubGraphIndex.size());
  //   assert(IndexV >= 0);
  //   assert(IndexV < SubGraphIndex[IndexN].size());

  //   // get the subgraph
  //   int IndexGraph = SubGraphIndex[IndexN][IndexV];
  //   assert(IndexGraph >= 0);
  //   assert(IndexGraph < SubGraphs.size());

  //   // get the corresponding face
  //   int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
  //   assert(IndexF >= 0);
  //   assert(IndexF < InternalFaces.size());

  //   // and the global vertex
  //   int GlobalV = InternalFaces[IndexF][IndexV];
  //   assert(GlobalV >= 0);
  //   assert(GlobalV < VF.size());

  //   // then get the local index
  //   int LocalIdx = SubGraphs[IndexGraph].LocalIndex(IndexN);
  //   assert(LocalIdx >= 0);
  //   assert(LocalIdx < SubGraphs[IndexGraph].NodesIdx.size());

  //   // check if is double exit then is manifold
  //   if (SubGraphs[IndexGraph].IsDoubleExit[LocalIdx])
  //     return MANIFOLD;

  //   // else check the number of exit edges can't be one
  //   int NumExitEdges = 0;
  //   for (size_t n = 0; n < SubGraphs[IndexGraph].Neigh[LocalIdx].size(); n++)
  //   {
  //     if (SubGraphs[IndexGraph].Active[LocalIdx][n])
  //       continue;
  //     if (SubGraphs[IndexGraph]
  //             .IsExit[SubGraphs[IndexGraph].Neigh[LocalIdx][n]])
  //       NumExitEdges++;
  //   }
  //   if (NumExitEdges == 1)
  //     return INCORRECT_EXIT;

  //   std::vector<int> PathNodes;
  //   ManifoldType TypeResult = UNDETERMINED;

  //   // then should have two path to exit nodes
  //   if (NumExitEdges >= 2) {
  //     // else check if can reach two different exit nodes
  //     bool HasTwoExits = HasTwoDifferentExitPaths(IndexN, IndexV, PathNodes);
  //     if (HasTwoExits)
  //       TypeResult = MANIFOLD;
  //     else
  //       TypeResult = NON_MANIFOLD;
  //   } else {
  //     bool HasCycle =
  //         SubGraphs[IndexGraph].FindCycleThroughGlobalNode(IndexN,
  //         PathNodes);

  //     // bool HasCycle = HasCycleAround(IndexN, IndexV, PathNodes);
  //     if (HasCycle)
  //       TypeResult = MANIFOLD;
  //     else
  //       TypeResult = NON_MANIFOLD;
  //   }

  //   // set the path nodes and vertices
  //   PathNodesVertices.clear();
  //   if (TypeResult == MANIFOLD) {
  //     for (const int &PN : PathNodes) {
  //       // get the corresponding face
  //       int IndexF_PN = Field::FrameFieldM3Graph::ElemOfNode(PN);
  //       assert(IndexF_PN >= 0);
  //       assert(IndexF_PN < InternalFaces.size());
  //       // and the index of vertex
  //       int IndexLocV_PN = Geo::WhichVert(InternalFaces[IndexF_PN], GlobalV);
  //       assert(IndexLocV_PN >= 0);
  //       PathNodesVertices.push_back(std::make_pair(PN, IndexLocV_PN));
  //     }
  //   }

  //   return TypeResult;
  // }

  // return the manifold type around a node vertex also return all the nodes
  //  and the corresponding vertices are satisfied by checking the paths
  ManifoldType
  CheckManifoldAround(const int &IndexN, const int &IndexV,
                      std::vector<std::pair<int, int>> &PathNodesVertices) {
    using clock = std::chrono::steady_clock;

    assert(IndexN >= 0);
    assert(IndexN < SubGraphIndex.size());
    assert(IndexV >= 0);
    assert(IndexV < SubGraphIndex[IndexN].size());

    // get the subgraph
    int IndexGraph = SubGraphIndex[IndexN][IndexV];
    assert(IndexGraph >= 0);
    assert(IndexGraph < SubGraphs.size());

    // get the corresponding face
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF >= 0);
    assert(IndexF < InternalFaces.size());

    // and the global vertex
    int GlobalV = InternalFaces[IndexF][IndexV];
    assert(GlobalV >= 0);
    assert(GlobalV < VF.size());

    // temporary disable all vertices that are not in the kernel
    auto t0 = clock::now();
    std::set<int> OnlyActive(SubGraphs[IndexGraph].NodesIdx.begin(),
                             SubGraphs[IndexGraph].NodesIdx.end());

    std::vector<int> DeactivatedNodes;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      if (OnlyActive.count(i) == 0 && M3Graph.IsEnabled(i)) {
        M3Graph.DisableNode(i);
        DeactivatedNodes.push_back(i);
      }
    }
    auto t1 = clock::now();

    SheetSolve_ConstraintsM3<ScalarType> ConstrM3(M3Graph, InternalFaces);
    SheetSolveSATM3<ScalarType> ShSATM3(ConstrM3, M3Graph, InternalFaces);
    std::vector<int> Seeds(1, IndexN);
    ConstrM3.InitConstraints(Seeds, true, 0, false);
    auto t2 = clock::now();

    typename SheetSolveSATM3<ScalarType>::SAT_ExtractParam ExtrP;
    ExtrP.HardNodes.push_back(IndexN);
    std::vector<bool> MatchedSeed;
    ShSATM3.print_msg = false;
    bool solved = ShSATM3.Solve(ExtrP, MatchedSeed);
    auto t3 = clock::now();

    PathNodesVertices.clear();
    if (solved) {
      // extract all the nodes that are active in the solution
      for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
        if (!M3Graph.IsEnabled(i))
          continue;
        if (ConstrM3.M3Solution[i]) {
          // get the corresponding face
          int IndexF_PN = Field::FrameFieldM3Graph::ElemOfNode(i);
          assert(IndexF_PN >= 0);
          assert(IndexF_PN < InternalFaces.size());
          // and the index of vertex
          int IndexLocV_PN = Geo::WhichVert(InternalFaces[IndexF_PN], GlobalV);
          assert(IndexLocV_PN >= 0);
          PathNodesVertices.push_back(std::make_pair(i, IndexLocV_PN));
        }
      }
    }
    auto t4 = clock::now();
    for (size_t i = 0; i < DeactivatedNodes.size(); i++)
      M3Graph.EnableNode(DeactivatedNodes[i]);
    auto t5 = clock::now();

    TimeDeactivatingNodes +=
        std::chrono::duration<ScalarType, std::milli>(t1 - t0).count();
    TimeSettingConstraints +=
        std::chrono::duration<ScalarType, std::milli>(t2 - t1).count();
    TimeSolvingSat +=
        std::chrono::duration<ScalarType, std::milli>(t3 - t2).count();
    TimeCollectingResults +=
        std::chrono::duration<ScalarType, std::milli>(t4 - t3).count();
    TimeDeactivatingNodes +=
        std::chrono::duration<ScalarType, std::milli>(t5 - t4).count(); 
    if (solved)
      return MANIFOLD;
    else
      return NON_MANIFOLD;
  }

  void UpdateManifoldNodes() {
    TimeSolvingSat = 0;
    TimeSettingConstraints = 0;
    TimeCollectingResults = 0;
    TimeDeactivatingNodes = 0;

    ManifoldNodes.clear();

    ManifoldNodevType = std::vector<std::array<ManifoldType, 3>>(
        M3Graph.NumNodes(), {UNDETERMINED, UNDETERMINED, UNDETERMINED});

    ManifoldNodesIndex =
        std::vector<std::array<int, 3>>(M3Graph.NumNodes(), {-1, -1, -1});

    std::cout << "Number of nodes: " << M3Graph.NumNodes() << std::endl;
    std::cout << "Number of nodes to check manifoldness: "
              << M3Graph.NumNodes() * 3 << std::endl;

    int curr_interval = 1000;

    int iterated_count = 0;
    int non_enabled = 0;
    int already_computed = 0;
    int has_computed = 0;
    int manifold_found = 0;
    int non_manifold_found = 0;

    std::set<std::pair<int, int>> checked_pairs;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {

      iterated_count += 3;
      if (!M3Graph.IsEnabled(i)) {
        non_enabled += 3;
        continue;
      }

      for (size_t j = 0; j < 3; j++) {

        assert(SubGraphIndex[i][j] != -1);

        if (ManifoldNodevType[i][j] != UNDETERMINED)
        // if (checked_pairs.count(std::make_pair(i,j))>0)
        {
          already_computed++;
          continue; // already computed
        }

        has_computed++;
        checked_pairs.insert(std::make_pair(i, j));
        std::vector<std::pair<int, int>> PathNodesVertices;
        ManifoldType type = CheckManifoldAround(i, j, PathNodesVertices);
        ManifoldNodevType[i][j] = type;
        // if namifold propagate to all the path nodes and vertices
        if (type == MANIFOLD) {

          // push back new manifold node the node
          ManifoldNodes.resize(ManifoldNodes.size() + 1);
          for (size_t k = 0; k < PathNodesVertices.size(); k++) {
            int PN = PathNodesVertices[k].first;
            ManifoldNodes.back().push_back(PN);
          }

          ManifoldNodesIndex[i][j] = ManifoldNodes.size() - 1;
          // std::cout << "Manifold node found: " << ManifoldNodes.size()
          //           << " with " << PathNodesVertices.size() << " nodes."
          //           << std::endl;

          for (const auto &PNV : PathNodesVertices) {
            int PN = PNV.first;
            int PV = PNV.second;
            assert(PN < ManifoldNodevType.size());
            assert(PN >= 0);
            assert(PV >= 0);
            assert(PV < 3);

            checked_pairs.insert(std::make_pair(PN, PV));
            if (ManifoldNodevType[PN][PV] == UNDETERMINED) {
              manifold_found++;
              ManifoldNodevType[PN][PV] = MANIFOLD;
              ManifoldNodesIndex[PN][PV] = ManifoldNodes.size() - 1;
            }
          }
        } else {
          non_manifold_found++;
        }
        // print progress after every 1000 nodes
        if (iterated_count >= curr_interval) {
          std::cout << "Progress Report:" << std::endl;

          std::cout << "  - Iterated over: " << iterated_count << std::endl;

          std::cout << "  - Non-enabled nodes skipped: " << non_enabled
                    << std::endl;

          std::cout << "  - Has Computed " << has_computed << " nodes."
                    << std::endl;

          std::cout << "  - Already computed: " << already_computed << " nodes."
                    << std::endl;
          std::cout << " * Time Setting Constraints: " << TimeSettingConstraints
                    << " ms." << std::endl;
          std::cout << " * Time Solving SAT: " << TimeSolvingSat << " ms."
                    << std::endl;
          std::cout << " * Time Collecting Results: " << TimeCollectingResults
                    << " ms." << std::endl;
          std::cout << " * Time Deactivating/Activating Nodes: "
                    << TimeDeactivatingNodes << " ms." << std::endl;  
          // std::cout << "  - Checked pairs " << checked_pairs.size()<<
          // std::endl;

          // std::cout << "    Manifold found: " << manifold_found
          //           << " nodes." << std::endl;
          // std::cout << "    Non-manifold found: " << non_manifold_found
          //           << " nodes." << std::endl;
          curr_interval += 1000;
        }
      }
    }
    // then count and print the nomber of each type
    int NumManifold = 0;
    int NumNonManifold = 0;
    int NumUndetermined = 0;
    int NumIncorrectExit = 0;
    for (size_t i = 0; i < ManifoldNodevType.size(); i++) {
      if (!M3Graph.IsEnabled(i))
        continue;
      for (size_t j = 0; j < 3; j++) {
        // rewrite is as switch case
        switch (ManifoldNodevType[i][j]) {
        case MANIFOLD:
          NumManifold++;
          break;
        case NON_MANIFOLD:
          NumNonManifold++;
          break;
        case UNDETERMINED:
          NumUndetermined++;
          break;
        case INCORRECT_EXIT:
          NumIncorrectExit++;
          break;
        }
      }
    }
    std::cout << "Num Manifold: " << NumManifold << std::endl;
    std::cout << "Num NonManifold: " << NumNonManifold << std::endl;
    std::cout << "Num IncorrectExit: " << NumIncorrectExit << std::endl;
    std::cout << "Num Undetermined: " << NumUndetermined << std::endl;
  }

  void InitSubGraphConnections() {

    for (size_t i = 0; i < SubGraphs.size(); i++) {
      // allocate
      int NumNodes = SubGraphs[i].NodesIdx.size();
      SubGraphs[i].Neigh.resize(NumNodes);
      SubGraphs[i].Weight.resize(NumNodes);
      SubGraphs[i].Active.resize(NumNodes);

      // initialize mapping
      std::map<int, int> global_to_localMap;
      for (size_t j = 0; j < SubGraphs[i].NodesIdx.size(); j++) {
        int IndexN = SubGraphs[i].NodesIdx[j];
        global_to_localMap[IndexN] = j;
      }

      // init connections
      for (size_t j = 0; j < SubGraphs[i].NodesIdx.size(); j++) {
        int IndexLoc0 = j;
        int IndexGlob0 = SubGraphs[i].NodesIdx[j];

        assert(global_to_localMap.count(IndexGlob0) > 0);
        assert(global_to_localMap[IndexGlob0] == j);

        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        std::vector<bool> ActiveLink;
        M3Graph.AllNeigh(IndexGlob0, IndexNeigh, WeightNeigh, ActiveLink);

        for (size_t n = 0; n < IndexNeigh.size(); n++) {
          int IndexGlob1 = IndexNeigh[n];
          if (global_to_localMap.count(IndexGlob1) == 0)
            continue;

          int IndexLoc1 = global_to_localMap[IndexGlob1];
          assert(IndexLoc1 != IndexLoc0);
          assert(IndexLoc1 >= 0);
          assert(IndexLoc1 < SubGraphs[i].NodesIdx.size());

          SubGraphs[i].Neigh[IndexLoc0].push_back(IndexLoc1);
          SubGraphs[i].Weight[IndexLoc0].push_back(WeightNeigh[n]);
          SubGraphs[i].Active[IndexLoc0].push_back(ActiveLink[n]);
        }
      }
    }
  }

  void InitSubGraphExitNodes() {
    for (size_t i = 0; i < SubGraphs.size(); i++) {
      int IndexV = SubGraphs[i].IndexV;
      int NumNodes = SubGraphs[i].NodesIdx.size();

      SubGraphs[i].IsExit = std::vector<bool>(NumNodes, false);
      SubGraphs[i].IsDoubleExit = std::vector<bool>(NumNodes, false);

      for (size_t j = 0; j < SubGraphs[i].NodesIdx.size(); j++) {
        int IndexN = SubGraphs[i].NodesIdx[j];
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
        int IndexLocV = Geo::WhichVert(InternalFaces[IndexF], IndexV);
        assert(IndexLocV >= 0);
        int IndexLab0 = (IndexLocV + 2) % 3;
        int IndexLab1 = IndexLocV;
        bool IsExit0 = M3Graph.IsExitLabel(IndexN, IndexLab0);
        bool IsExit1 = M3Graph.IsExitLabel(IndexN, IndexLab1);
        if (IsExit0 || IsExit1)
          SubGraphs[i].IsExit[j] = true;
        if (IsExit0 && IsExit1)
          SubGraphs[i].IsDoubleExit[j] = true;
      }
    }
  }

  void InitSubGraphNodes() {
    SubGraphIndex.clear();
    SubGraphIndex.resize(M3Graph.NumNodes(), std::array<int, 3>{-1, -1, -1});

    int all_combinations = 0;
    int computed = 0;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      if (!M3Graph.IsEnabled(i))
        continue;
      for (size_t j = 0; j < 3; j++) {
        int IndexF0 = Field::FrameFieldM3Graph::ElemOfNode(i);
        int IndexV0 = InternalFaces[IndexF0][j];
        assert(IndexV0 >= 0);
        assert(IndexV0 < VF.size());

        // std::cout << "Checking Node " << i << " vert " << j << std::endl;

        all_combinations++;

        if (SubGraphIndex[i][j] != -1)
          continue; // alreadt computed

        computed++;
        std::vector<int> KernN = GetSurroundingNodeVert(i, j);

        SubGraphs.push_back(SubGraph(IndexV0, KernN));
        int CurrGraphIDX = SubGraphs.size() - 1;

        for (size_t n = 0; n < KernN.size(); n++) {
          // get the node
          int IndexN = KernN[n];

          // get corresponding face
          int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);

          // and the index of vertex
          int LocalV = Geo::WhichVert(InternalFaces[IndexF], IndexV0);
          assert(LocalV >= 0);

          assert(SubGraphIndex[IndexN][LocalV] == -1);
          SubGraphIndex[IndexN][LocalV] = CurrGraphIDX;
        }
      }
    }
    std::cout << "Computed:" << computed << " out of " << all_combinations
              << std::endl;
  }

public:
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  std::vector<std::vector<int>> &InternalFaces;
  Geo::TMark TempMark;

  // the subgraphs
  std::vector<SubGraph> SubGraphs;
  std::vector<std::vector<int>> ManifoldNodes;

  // std::vector < std::vector<std::pair<int, int>> NodeVertIndex;

  std::vector<std::array<int, 3>> SubGraphIndex;
  std::vector<std::array<ManifoldType, 3>> ManifoldNodevType;
  std::vector<std::array<int, 3>> ManifoldNodesIndex;

  std::vector<std::vector<int>> VF, VFIndex;

  std::vector<bool> HasSolution;
  std::vector<bool> Tested;

  // bool CheckHasSolution(int IndexN, const std::vector<int> &SubGraphNodes)
  // {
  //   SheetSolve_ConstraintsM3<ScalarType> ConstrM3(M3Graph, InternalFaces);
  //   SheetSolveSATM3<ScalarType> ShSATM3(ConstrM3, M3Graph, InternalFaces);
  //   std::vector<int> Seeds(1, IndexN);
  //   ConstrM3.InitConstraints(Seeds, true, 0);

  //   typename SheetSolveSATM3<ScalarType>::SAT_ExtractParam ExtrP;
  //   ExtrP.HardNodes.push_back(IndexN);
  //   std::vector<int> MatchedSeed;
  //   bool solved = ShSATM3.Solve(ExtrP, MatchedSeed);
  // }

  // bool CheckHasSolution(int IndexN, const std::vector<int> &SubGraphNodes) {
  //   SheetSolve_ConstraintsM3<ScalarType> ConstrM3(M3Graph, InternalFaces);
  //   SheetSolveSATM3<ScalarType> ShSATM3(ConstrM3, M3Graph, InternalFaces);
  //   std::vector<int> Seeds(1, IndexN);
  //   ConstrM3.InitConstraints(Seeds, true, 0);

  //   typename SheetSolveSATM3<ScalarType>::SAT_ExtractParam ExtrP;
  //   ExtrP.HardNodes.push_back(IndexN);
  //   std::vector<int> MatchedSeed;
  //   bool solved = ShSATM3.Solve(ExtrP, MatchedSeed);
  // }

  void UpdateHasSolution() {
    Tested = std::vector<bool>(M3Graph.NumNodes(), false);
    HasSolution = std::vector<bool>(M3Graph.NumNodes(), false);
  }

  std::vector<int>
  GetConnectedComponentFrom(const std::vector<int> &SeedN,
                            const std::set<int> &PossibleNodes) {

    std::vector<int> component;
    component.reserve(PossibleNodes.size());

    if (SeedN.empty())
      return component;

    std::queue<int> Q;
    std::unordered_set<int> visited;
    visited.reserve(PossibleNodes.size());

    // Push valid seeds
    for (int s : SeedN) {
      if (PossibleNodes.count(s)) {
        Q.push(s);
        visited.insert(s);
      }
    }

    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();

      component.push_back(u);

      // ---- Replace this with your neighbor access logic ----
      std::vector<int> nbrs; // = M3Graph.Neighbors(u);
      std::vector<ScalarType> WeightNeigh;
      M3Graph.AllNeigh(u, nbrs, WeightNeigh);
      // -------------------------------------------------------

      for (int v : nbrs) {
        if (PossibleNodes.count(v) && !visited.count(v)) {
          visited.insert(v);
          Q.push(v);
        }
      }
    }

    return component;
  }

  std::vector<int> GetSurroundingNodeVert(const int IndexN, const int IndexV) {
    assert(IndexN >= 0);
    assert(IndexN < M3Graph.NumNodes());
    assert(IndexV >= 0);
    assert(IndexV < 3);

    // get the corresponding face
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF >= 0);
    assert(IndexF < InternalFaces.size());
    assert(InternalFaces[IndexF].size() == 3);

    std::vector<int> AllNeighF;

    int GlobalV = InternalFaces[IndexF][IndexV];
    assert(GlobalV >= 0);
    assert(GlobalV < VF.size());
    std::vector<int> NeighF = VF[GlobalV];
    AllNeighF.insert(AllNeighF.end(), NeighF.begin(), NeighF.end());

    // get all possible neighbouring nodes
    std::vector<int> AllNeighN;
    Field::FrameFieldM3Graph::NodesOfElems(AllNeighF, AllNeighN);
    std::set<int> AllNeighNSet(AllNeighN.begin(), AllNeighN.end());

    std::vector<int> SeedN(1, IndexN);
    std::vector<int> ReachN = GetConnectedComponentFrom(SeedN, AllNeighNSet);
    return ReachN;
  }

  void GetNonManfoldNodes(std::vector<std::pair<int, int>> &NonManifoldNodev) {
    NonManifoldNodev.clear();
    for (size_t i = 0; i < ManifoldNodevType.size(); i++) {
      if (!M3Graph.IsEnabled(i))
        continue;
      for (size_t j = 0; j < 3; j++) {
        if (ManifoldNodevType[i][j] == NON_MANIFOLD) {
          NonManifoldNodev.push_back(std::make_pair(i, j));
        }
      }
    }
  }

  std::vector<int> GetManifoldNodesIndex(int IndexN, int IndexV) {
    std::vector<int> Result;
    if (ManifoldNodesIndex[IndexN][IndexV] == -1)
      return Result;

    int ManifoldIdx = ManifoldNodesIndex[IndexN][IndexV];
    Result = ManifoldNodes[ManifoldIdx];
    return Result;
    // std::vector<int> Result;
    // int SubGraphIdx = SubGraphIndex[IndexN][IndexV];
    // if (SubGraphIdx == -1)
    //   return Result;
    // return SubGraphs[SubGraphIdx].NodesIdx;
  }

  void Init() {
    // initialize connectivity

    Geo::getVFAdjacency(InternalFaces, VF, VFIndex);

    InitSubGraphNodes();

    InitSubGraphConnections();

    InitSubGraphExitNodes();

    UpdateManifoldNodes();
  }

  LabelGraph_ConnectivityM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                            std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif