#ifndef CONSTRAINT_GRAPH
#define CONSTRAINT_GRAPH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>

template <class ScalarType> class OrientFLowGraph {

  std::vector<std::pair<int, int>> Node_to_VarPatch;
  std::map<std::pair<int, int>, int> VarPatch_to_Node;
  // std::vector<std::vector<int> > ExitPatchNode,EnterPatchNode;
  std::vector<std::pair<int, int>> VarPatch;
  std::vector<std::vector<std::pair<int, int>>> VarPatchSide;

  std::vector<std::vector<int>> GraphAdj;
  std::vector<std::vector<ScalarType>> GraphWeight;
  std::vector<ScalarType> GraphNodeWeight;

  std::vector<std::vector<bool>> Enabled;
  std::vector<bool> EnabledN;

  std::vector<int> BoundarySinks;
  std::vector<int> SourceIncr;
  std::vector<int> SourceDecr;
  std::vector<int> TargetIncr;
  std::vector<int> TargetDecr;
  std::vector<int> SourceParity;
  std::vector<int> TargetParity;
  const std::vector<std::vector<std::vector<int>>> &PatchSideVariables;
  std::vector<int> &VarValue;
  std::vector<int> &VarTargetValue;

  bool WriteDebug;

  std::vector<std::pair<int, int>>
  GetArcsBetween(const int CurrPatch, const std::vector<int> &Side0,
                 const std::vector<int> &Side1) {

    std::vector<std::pair<int, int>> Arcs;
    for (size_t j = 0; j < Side0.size(); j++) {
      int IndexVar0 = Side0[j];
      // get if is first or second istance
      for (size_t k = 0; k < Side1.size(); k++) {
        int IndexVar1 = Side1[k];

        // connect them two ways
        int P0IndexV0 = VarPatch[IndexVar0].first;
        int P1IndexV0 = VarPatch[IndexVar0].second;
        assert(P0IndexV0 != P1IndexV0);

        int P0IndexV1 = VarPatch[IndexVar1].first;
        int P1IndexV1 = VarPatch[IndexVar1].second;
        assert(P0IndexV1 != P1IndexV1);

        assert((P0IndexV0 == CurrPatch) || (P1IndexV0 == CurrPatch));
        assert((P0IndexV1 == CurrPatch) || (P1IndexV1 == CurrPatch));

        int OtherPatch0 = P0IndexV0;
        if (P0IndexV0 == CurrPatch)
          OtherPatch0 = P1IndexV0;

        int OtherPatch1 = P0IndexV1;
        if (P0IndexV1 == CurrPatch)
          OtherPatch1 = P1IndexV1;

        // std::cout<<"Other Patch 0:"<<OtherPatch0<<std::endl;
        // std::cout<<"Other Patch 1:"<<OtherPatch1<<std::endl;
        // assert(((OtherPatch0 == -1) && (OtherPatch1 == -1)) ||
        //        ((OtherPatch0 != OtherPatch1)));

        // then add the two connectsions
        std::pair<int, int> KeyV0Patch(IndexVar0, CurrPatch);
        std::pair<int, int> KeyV1Patch(IndexVar1, CurrPatch);
        assert(VarPatch_to_Node.count(KeyV0Patch) > 0);
        assert(VarPatch_to_Node.count(KeyV1Patch) > 0);

        int NodePatchV0 = VarPatch_to_Node[KeyV0Patch];
        int NodePatchV1 = VarPatch_to_Node[KeyV1Patch];

        std::pair<int, int> KeyV0OthPatch(IndexVar0, OtherPatch0);
        std::pair<int, int> KeyV1OthPatch(IndexVar1, OtherPatch1);

        assert(VarPatch_to_Node.count(KeyV0OthPatch) > 0);
        assert(VarPatch_to_Node.count(KeyV1OthPatch) > 0);

        int NodeOthPV0 = VarPatch_to_Node[KeyV0OthPatch];
        int NodeOthPV1 = VarPatch_to_Node[KeyV1OthPatch];

        assert(((NodeOthPV0 == -1) && (NodeOthPV1 == -1)) ||
               ((NodeOthPV0 != NodeOthPV1)));

        Arcs.push_back(std::pair<int, int>(NodePatchV0, NodeOthPV1));
        Arcs.push_back(std::pair<int, int>(NodePatchV1, NodeOthPV0));
      }
    }
    return Arcs;
  }

  void
  InitAdjGraphFromLinks(const std::vector<std::pair<int, int>> &Connections,
                        const std::vector<ScalarType> &ConnectionsW) {
    assert(Connections.size() == ConnectionsW.size());

    int MaxNodes = 0;
    for (size_t i = 0; i < Connections.size(); i++) {
      assert(Connections[i].first >= 0);
      assert(Connections[i].second >= 0);
      MaxNodes = std::max(MaxNodes, Connections[i].first);
      MaxNodes = std::max(MaxNodes, Connections[i].second);
    }

    GraphAdj.clear();
    GraphWeight.clear();
    GraphAdj.resize(MaxNodes + 1);
    GraphWeight.resize(MaxNodes + 1);

    for (size_t i = 0; i < Connections.size(); i++) {
      int IndexN0 = Connections[i].first;
      int IndexN1 = Connections[i].second;
      assert(IndexN0 >= 0);
      assert(IndexN1 >= 0);
      assert(IndexN0 < GraphAdj.size());
      assert(IndexN1 < GraphWeight.size());

      GraphAdj[IndexN0].push_back(IndexN1);
      GraphWeight[IndexN0].push_back(ConnectionsW[i]);
    }
  }

  int SideValue(int IndexP, int IndexSide) {
    assert(IndexP >= 0);
    assert(IndexP < PatchSideVariables.size());
    assert(IndexSide >= 0);
    assert(IndexSide < PatchSideVariables[IndexP].size());
    int SumVal = 0;
    for (size_t i = 0; i < PatchSideVariables[IndexP][IndexSide].size(); i++) {
      int IndexVar = PatchSideVariables[IndexP][IndexSide][i];
      SumVal += VarValue[IndexVar];
    }
    return SumVal;
  }

  std::vector<int> FindSinks() {
    std::vector<int> Sinks;
    for (size_t i = 0; i < GraphAdj.size(); i++)
      if (GraphAdj[i].size() == 0)
        Sinks.push_back(i);
    return Sinks;
  }

  int OtherPatch(const int IndexVar, const int IndexPatch) {
    assert(IndexVar >= 0);
    assert(IndexVar < VarPatch.size());
    assert((VarPatch[IndexVar].first == IndexPatch) ||
           (VarPatch[IndexVar].second == IndexPatch));

    if (VarPatch[IndexVar].first == IndexPatch)
      return VarPatch[IndexVar].second;
    else
      return VarPatch[IndexVar].first;
  }

  std::vector<int> OtherPatch(const std::vector<int> &IndexVar,
                              const int IndexPatch) {
    std::vector<int> resPatch;
    for (size_t i = 0; i < IndexVar.size(); i++) {
      int PatchIdx = OtherPatch(IndexVar[i], IndexPatch);
      resPatch.push_back(PatchIdx);
    }
    return resPatch;
  }

  std::vector<int> PatchNodes(const std::vector<int> &IndexVar,
                              const int IndexPatch) {
    std::vector<int> resNodes;
    for (size_t i = 0; i < IndexVar.size(); i++) {
      std::pair<int, int> VarPKey(IndexVar[i], IndexPatch);
      int IndexN = VarPatch_to_Node[VarPKey];
      resNodes.push_back(IndexN);
    }
    return resNodes;
  }

  std::vector<int> OtherPatchNodes(const std::vector<int> &IndexVar,
                                   const int IndexPatch) {
    std::vector<int> resNodes;
    std::vector<int> othPatches = OtherPatch(IndexVar, IndexPatch);
    assert(othPatches.size() == IndexVar.size());

    for (size_t i = 0; i < IndexVar.size(); i++) {
      std::pair<int, int> VarPKey(IndexVar[i], othPatches[i]);
      assert(VarPatch_to_Node.count(VarPKey) > 0);

      int IndexN = VarPatch_to_Node[VarPKey];
      assert(IndexN >= 0);
      assert(IndexN < Node_to_VarPatch.size());
      resNodes.push_back(IndexN);
    }
    return resNodes;
  }

  ScalarType Energy(int IndexVar, int TestValue) {
    assert(TestValue >=
           0); // in realigy should be at lest 1 need to reorganize code
    assert(IndexVar >= 0);
    assert(IndexVar < VarTargetValue.size());
    ScalarType TargetVal = VarTargetValue[IndexVar];
    ScalarType Diff = pow((TestValue - TargetVal), 2);
    return Diff;
  }

  ScalarType Energy(int IndexVar) {
    assert(IndexVar >= 0);
    assert(IndexVar < VarValue.size());
    int TestValue = VarValue[IndexVar];
    return Energy(IndexVar, TestValue);
  }
  // ScalarType PriorityChange(int IndexN, int Diff, bool testDiff = false) {
  //   int IndexVar = VarOfNode(IndexN);
  //   ScalarType CurrVal = VarValue[IndexVar];
  //   ScalarType NextVal = CurrVal + Diff;
  //   ScalarType TargetVal = VarTargetValue[IndexVar];

  //   ScalarType Diff0 = pow((CurrVal - TargetVal), 2);
  //   ScalarType Diff1 = pow((NextVal - TargetVal), 2);
  //   // ScalarType AvgDiff = (Diff0 + Diff1) / 2;

  //   // if (Diff1 < Diff0)
  //   //   return -AvgDiff;
  //   // return AvgDiff;
  //   return (Diff1-Diff0);
  // }

  ScalarType PriorityChange(int IndexN, int Diff, bool testDiff = false) {
    int IndexVar = VarOfNode(IndexN);
    ScalarType CurrVal = VarValue[IndexVar];
    ScalarType NextVal = CurrVal + Diff;

    ScalarType Energy0 = Energy(IndexVar, CurrVal);
    ScalarType Energy1 = Energy(IndexVar, NextVal);
    return (Energy1 - Energy0);
  }

  void SetupNodeWeightChange(int Diff, bool make_positive) {
    // compute the change
    GraphNodeWeight = std::vector<ScalarType>(GraphAdj.size(), 0);
    for (size_t i = 0; i < GraphNodeWeight.size(); i++) {
      GraphNodeWeight[i] = PriorityChange(i, Diff);
    }
    if (!make_positive)
      return;

    ScalarType minValue =
        *std::min_element(GraphNodeWeight.begin(), GraphNodeWeight.end());

    // put strictly positive
    for (size_t i = 0; i < GraphNodeWeight.size(); i++) {
      GraphNodeWeight[i] -= minValue;
      GraphNodeWeight[i] += 0.00001;
      assert(GraphNodeWeight[i] > 0);
    }
  }

  void SetupWeightChange(int Diff, bool make_positive) {

    SetupNodeWeightChange(Diff, make_positive);

    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int Node0 = i;
      ScalarType W0 = GraphNodeWeight[Node0];
      for (size_t j = 0; j < GraphAdj[i].size(); j++) {
        int Node1 = GraphAdj[i][j];
        ScalarType W1 = GraphNodeWeight[Node1];
        int IndexVar1 = VarOfNode(GraphAdj[i][j]);
        GraphWeight[i][j] = (W0 + W1) / 2;
      }
    }
  }

  void SetupUniformWeight(ScalarType W) {
    for (size_t i = 0; i < GraphAdj.size(); i++)
      for (size_t j = 0; j < GraphAdj[i].size(); j++) {
        GraphWeight[i][j] = W;
      }
  }

  bool SetupSolve(bool equality, bool parity) {
    assert(equality || parity);
    assert(!(equality && parity));

    BoundarySinks = FindSinks();
    SourceIncr.clear();
    SourceDecr.clear();
    TargetDecr.clear();
    TargetIncr.clear();
    SourceParity.clear();
    TargetParity.clear();

    bool has_to_solve = false;
    for (size_t i = 0; i < PatchSideVariables.size(); i++) {
      assert(PatchSideVariables[i].size() == 4);
      int Size0 = SideValue(i, 0);
      int Size1 = SideValue(i, 1);
      int Size2 = SideValue(i, 2);
      int Size3 = SideValue(i, 3);
      std::vector<int> DecrVar0, IncrVar0, DecrVar1, IncrVar1;

      if (equality) {
        if ((Size0 == Size2) && (Size1 == Size3))
          continue;
      }
      if (parity) {
        if (((Size0 + Size1 + Size2 + Size3) % 2) == 0)
          continue;
      }

      has_to_solve = true;

      if (Size0 < Size2) {
        DecrVar0 = PatchSideVariables[i][2];
        IncrVar0 = PatchSideVariables[i][0];
      }
      if (Size0 > Size2) {
        DecrVar0 = PatchSideVariables[i][0];
        IncrVar0 = PatchSideVariables[i][2];
      }

      if (Size1 < Size3) {
        DecrVar1 = PatchSideVariables[i][3];
        IncrVar1 = PatchSideVariables[i][1];
      }
      if (Size1 > Size3) {
        DecrVar1 = PatchSideVariables[i][1];
        IncrVar1 = PatchSideVariables[i][3];
      }

      std::vector<int> DecrVar, IncrVar;
      DecrVar = DecrVar0;
      DecrVar.insert(DecrVar.end(), DecrVar1.begin(), DecrVar1.end());

      IncrVar = IncrVar0;
      IncrVar.insert(IncrVar.end(), IncrVar1.begin(), IncrVar1.end());

      std::vector<int> CurrTargetIncr = PatchNodes(IncrVar, i);
      std::vector<int> CurrTargetDecr = PatchNodes(DecrVar, i);

      TargetIncr.insert(TargetIncr.end(), CurrTargetIncr.begin(),
                        CurrTargetIncr.end());
      TargetDecr.insert(TargetDecr.end(), CurrTargetDecr.begin(),
                        CurrTargetDecr.end());

      std::vector<int> CurrSourceIncr = OtherPatchNodes(IncrVar, i);
      std::vector<int> CurrSourceDecr = OtherPatchNodes(DecrVar, i);

      SourceIncr.insert(SourceIncr.end(), CurrSourceIncr.begin(),
                        CurrSourceIncr.end());
      SourceDecr.insert(SourceDecr.end(), CurrSourceDecr.begin(),
                        CurrSourceDecr.end());
    }

    // that is the same can increase or decrease
    SourceParity = SourceIncr;
    SourceParity.insert(SourceParity.end(), SourceDecr.begin(),
                        SourceDecr.end());
    TargetParity = TargetIncr;
    TargetParity.insert(TargetParity.end(), TargetIncr.begin(),
                        TargetIncr.end());

    for (size_t i = 0; i < SourceIncr.size(); i++) {
      assert(SourceIncr[i] >= 0);
      assert(SourceIncr[i] < Node_to_VarPatch.size());
    }

    for (size_t i = 0; i < TargetIncr.size(); i++) {
      assert(TargetIncr[i] >= 0);
      assert(TargetIncr[i] < Node_to_VarPatch.size());
    }

    for (size_t i = 0; i < SourceDecr.size(); i++) {
      assert(SourceDecr[i] >= 0);
      assert(SourceDecr[i] < Node_to_VarPatch.size());
    }

    for (size_t i = 0; i < TargetDecr.size(); i++) {
      assert(TargetDecr[i] >= 0);
      assert(TargetDecr[i] < Node_to_VarPatch.size());
    }

    // std::cout << "Source Increase Num:" << SourceIncr.size() << std::endl;
    // std::cout << "Source Decrease Num:" << SourceDecr.size() << std::endl;

    // std::cout << "Target Increase Num:" << TargetIncr.size() << std::endl;
    // std::cout << "Targed Decrease Num:" << TargetDecr.size() << std::endl;

    // enable all
    Enabled.clear();
    Enabled.resize(GraphAdj.size());
    for (size_t i = 0; i < GraphAdj.size(); i++)
      Enabled[i] = std::vector<bool>(GraphAdj[i].size(), true);

    // disable the sink
    for (size_t i = 0; i < TargetIncr.size(); i++) {
      int IndexN = TargetIncr[i];
      assert(IndexN >= 0);
      assert(IndexN < Enabled.size());
      Enabled[IndexN] = std::vector<bool>(GraphAdj[IndexN].size(), false);
    }

    for (size_t i = 0; i < TargetDecr.size(); i++) {
      int IndexN = TargetDecr[i];
      assert(IndexN >= 0);
      assert(IndexN < Enabled.size());
      Enabled[IndexN] = std::vector<bool>(GraphAdj[IndexN].size(), false);
    }

    return has_to_solve;
  }

  void InitGraph() {

    // count number of variables
    int MaxV = 0;
    for (size_t i = 0; i < PatchSideVariables.size(); i++)
      for (size_t j = 0; j < PatchSideVariables[i].size(); j++)
        for (size_t k = 0; k < PatchSideVariables[i][j].size(); k++)
          MaxV = std::max(MaxV, PatchSideVariables[i][j][k]);

    VarPatch.clear();
    VarPatch.resize(MaxV + 1, std::pair<int, int>(-1, -1));

    VarPatchSide.clear();
    VarPatchSide.resize(MaxV + 1);

    for (size_t i = 0; i < PatchSideVariables.size(); i++)
      for (size_t j = 0; j < PatchSideVariables[i].size(); j++)
        for (size_t k = 0; k < PatchSideVariables[i][j].size(); k++) {
          int IndexVar = PatchSideVariables[i][j][k];
          VarPatchSide[IndexVar].push_back(std::pair<int, int>(i, j));

          assert(IndexVar >= 0);
          assert(IndexVar < VarPatch.size());
          assert((VarPatch[IndexVar].first == -1) ||
                 (VarPatch[IndexVar].second == -1));
          if (VarPatch[IndexVar].first == -1)
            VarPatch[IndexVar].first = i;
          else
            VarPatch[IndexVar].second = i;
        }

    // there are two nodes for each variable (one for each patch)
    Node_to_VarPatch.clear();
    VarPatch_to_Node.clear();

    Node_to_VarPatch.resize((MaxV + 1) * 2);
    for (size_t i = 0; i < VarPatch.size(); i++) {
      int IndexP0 = VarPatch[i].first;
      int IndexP1 = VarPatch[i].second;
      assert(IndexP0 != IndexP1);

      int IndexN0 = i * 2;
      int IndexVar = i;
      Node_to_VarPatch[IndexN0].first = IndexVar;
      Node_to_VarPatch[IndexN0].second = IndexP0;
      std::pair<int, int> VarPatch0(IndexVar, IndexP0);
      assert(VarPatch_to_Node.count(VarPatch0) == 0);
      VarPatch_to_Node[VarPatch0] = IndexN0;

      int IndexN1 = (i * 2) + 1;
      Node_to_VarPatch[IndexN1].first = IndexVar;
      Node_to_VarPatch[IndexN1].second = IndexP1;
      std::pair<int, int> VarPatch1(IndexVar, IndexP1);
      assert(VarPatch_to_Node.count(VarPatch1) == 0);
      VarPatch_to_Node[VarPatch1] = IndexN1;
    }

    // now find the connections
    std::vector<std::pair<int, int>> Connections;
    std::vector<ScalarType> ConnectionsW;
    for (size_t i = 0; i < PatchSideVariables.size(); i++) {
      int CurrPatch = i;
      assert(PatchSideVariables[i].size() == 4);
      std::vector<int> Side0 = PatchSideVariables[CurrPatch][0];
      std::vector<int> Side1 = PatchSideVariables[CurrPatch][1];
      std::vector<int> Side2 = PatchSideVariables[CurrPatch][2];
      std::vector<int> Side3 = PatchSideVariables[CurrPatch][3];
      std::vector<std::pair<int, int>> Arcs0 = GetArcsBetween(i, Side0, Side2);
      std::vector<std::pair<int, int>> Arcs1 = GetArcsBetween(i, Side1, Side3);
      Connections.insert(Connections.end(), Arcs0.begin(), Arcs0.end());
      Connections.insert(Connections.end(), Arcs1.begin(), Arcs1.end());
    }
    ConnectionsW.resize(Connections.size(), 1);
    std::cout << "There are " << Connections.size() << " Connections "
              << std::endl;

    InitAdjGraphFromLinks(Connections, ConnectionsW);
  }

  struct ShortestPathResult {
    ScalarType length;     // sum of GraphWeight along the path
    std::vector<int> path; // nodes from source to target (inclusive)
    bool found;            // false if no path exists
  };

  int VarOfNode(const int &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < Node_to_VarPatch.size());
    return (Node_to_VarPatch[IndexN].first);
  }

  void CheckPath(const ShortestPathResult &SRes) {
    // for (size_t i=0;i<SRes.path.size();i++)
    // assert(EnabledN[SRes.path[i]]);

    int Var0 = VarOfNode(SRes.path[0]);
    int Var1 = VarOfNode(SRes.path.back());
    if (Var0 == Var1) {
      std::cout << "Cycle" << std::endl;
      assert(0);
    }
  }

  void CheckPath(const std::vector<int> &PathN) {
    int Var0 = VarOfNode(PathN[0]);
    int Var1 = VarOfNode(PathN.back());
    if (Var0 == Var1) {
      std::cout << "Cycle" << std::endl;
      assert(0);
    }
  }

  ShortestPathResult ShortestPath(const std::vector<int> &sourceN,
                                  const std::vector<int> &targetN) {
    const int n = (int)GraphAdj.size();
    const ScalarType INF =
        std::numeric_limits<ScalarType>::max() / ScalarType(4);

    // Basic consistency checks
    if ((int)GraphWeight.size() != n || (int)Enabled.size() != n ||
        (int)EnabledN.size() != n) {
      assert(false && "Inconsistent sizes");
    }
    for (int u = 0; u < n; ++u) {
      if (GraphAdj[u].size() != GraphWeight[u].size() ||
          GraphAdj[u].size() != Enabled[u].size()) {
        assert(false && "Adj/Weight/Enabled row size mismatch");
      }
    }

    // Mark targets (only enabled nodes are valid targets)
    std::vector<char> isTarget(n, 0);
    for (int t : targetN) {
      if (t >= 0 && t < n && EnabledN[t])
        isTarget[t] = 1;
    }

    // Dijkstra init (multi-source), only for enabled nodes
    std::vector<ScalarType> dist(n, std::numeric_limits<ScalarType>::max());
    std::vector<int> parent(n, -1);

    using QItem = std::pair<ScalarType, int>; // (distance, node)
    std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;

    for (int s : sourceN) {
      if (s >= 0 && s < n && EnabledN[s]) {
        dist[s] = ScalarType(0);
        parent[s] = -1;
        pq.emplace(ScalarType(0), s);
      }
    }

    int reachedTarget = -1;

    // Dijkstra loop
    while (!pq.empty()) {
      auto [du, u] = pq.top();
      pq.pop();
      // Use ">" to avoid equality issues with floating ScalarType
      if (du > dist[u])
        continue;
      if (!EnabledN[u])
        continue;

      if (isTarget[u]) { // early exit on first settled target
        reachedTarget = u;
        break;
      }

      const auto &nbrs = GraphAdj[u];
      const auto &wts = GraphWeight[u];
      const auto &en = Enabled[u];

      for (size_t k = 0; k < nbrs.size(); ++k) {
        if (!en[k])
          continue; // edge disabled
        const int v = nbrs[k];
        if (v < 0 || v >= n || !EnabledN[v])
          continue; // invalid/disabled node

        const ScalarType w = static_cast<ScalarType>(wts[k]);
        if (w < ScalarType(0))
          continue; // Dijkstra needs non-negative weights

        const ScalarType nd = du + w;
        if (nd < dist[v]) {
          dist[v] = nd;
          parent[v] = u;
          pq.emplace(nd, v);
        }
      }
    }

    if (reachedTarget == -1) {
      return {INF, {}, false};
    }

    // Reconstruct path
    std::vector<int> path;
    for (int v = reachedTarget; v != -1; v = parent[v])
      path.push_back(v);
    std::reverse(path.begin(), path.end());

    return {dist[reachedTarget], path, true};
  }

  void ApplyIncrease(const ShortestPathResult &PathNodes) {
    CheckPath(PathNodes);
    // std::cout << "* Increase PATH SIZE:" << PathNodes.size() << std::endl;
    for (size_t i = 0; i < PathNodes.path.size(); i++) {
      int IndexN = PathNodes.path[i];
      int VarNum = VarOfNode(IndexN);
      VarValue[VarNum]++;
    }
  }

  void ApplyDecrease(const ShortestPathResult &PathNodes) {
    CheckPath(PathNodes);
    // std::cout << "* Decrease PATH SIZE:" << PathNodes.size() << std::endl;
    for (size_t i = 0; i < PathNodes.path.size(); i++) {
      int IndexN = PathNodes.path[i];
      int VarNum = VarOfNode(IndexN);
      VarValue[VarNum]--;
    }
  }

  // bool CanDecrease(std::vector<int> &NodeSeq) {
  //   for (size_t j = 0; j < NodeSeq.size(); j++) {
  //     int IndexN = NodeSeq[j];
  //     assert(EnabledN[IndexN]);
  //     int IndexVar = VarOfNode(IndexN);
  //     if (!VarValue[IndexVar] > 1)
  //       return false;
  //   }
  //   return true;
  // }

  bool HasSameVarMultiple(std::vector<int> &NodeSeq) {
    for (size_t j = 0; j < NodeSeq.size(); j++) {
      int IndexN = NodeSeq[j];
      assert(EnabledN[IndexN]);
      int IndexVar = VarOfNode(IndexN);
      if (!VarValue[IndexVar] > 1)
        return false;
    }
    return true;
  }

  bool CanDecrease(std::vector<int> &PathNodes) {

    std::vector<int> NodeIdx = PathNodes;
    std::sort(NodeIdx.begin(), NodeIdx.end());
    NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()), NodeIdx.end());

    std::vector<int> TempVal = VarValue;

    for (size_t j = 0; j < NodeIdx.size(); j++) {
      int IndexN = NodeIdx[j];
      int IndexVar = VarOfNode(IndexN);
      TempVal[IndexVar]--;
      if (TempVal[IndexVar] < 1)
        return false;
    }

    return true;
  }

  void ApplyDecrease(const std::vector<int> &PathNodes) {
    // check if cycle
    std::vector<int> NodeIdx = PathNodes;
    std::sort(NodeIdx.begin(), NodeIdx.end());
    NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()), NodeIdx.end());

    assert(CanDecrease(NodeIdx));

    for (size_t i = 0; i < NodeIdx.size(); i++) {
      int IndexN = NodeIdx[i];
      int VarNum = VarOfNode(IndexN);
      assert(VarValue[VarNum] > 1);
      VarValue[VarNum]--;
    }
  }

  void ApplyIncrease(const std::vector<int> &PathNodes) {
    // int lastN=PathNodes.size();
    // if (PathNodes[0]==PathNodes.back())lastN--;
    std::vector<int> NodeIdx = PathNodes;
    std::sort(NodeIdx.begin(), NodeIdx.end());
    NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()), NodeIdx.end());

    for (size_t i = 0; i < NodeIdx.size(); i++) {
      int IndexN = NodeIdx[i];
      int VarNum = VarOfNode(IndexN);
      VarValue[VarNum]++;
    }
  }

  std::vector<int> AdmissibleTargets(const int SourceN,
                                     const std::vector<int> &TargetCand,
                                     bool check_common_patch = false) {

    int var_source = VarOfNode(SourceN);

    std::vector<int> curr_target;

    int SourceVar = VarOfNode(SourceN);
    std::vector<std::pair<int, int>> source_p_side = VarPatchSide[SourceVar];
    // PatchSideVariables

    // set the current target filtering the
    // same variable of source
    for (size_t i = 0; i < TargetCand.size(); i++) {
      int test_target = TargetCand[i];

      int var_test_target = VarOfNode(test_target);

      if (var_test_target == var_source)
        continue;

      std::vector<std::pair<int, int>> target_p_side =
          VarPatchSide[var_test_target];
      bool has_common_sideP = false;
      for (size_t j = 0; j < source_p_side.size(); j++)
        for (size_t k = 0; k < target_p_side.size(); k++)
          if (target_p_side[k] == source_p_side[j])
            has_common_sideP = true;

      if (has_common_sideP)
        continue;

      if (check_common_patch) {
        int SourceP0 = VarPatch[SourceVar].first;
        int SourceP1 = VarPatch[SourceVar].second;
        int TargetP0 = VarPatch[var_test_target].first;
        int TargetP1 = VarPatch[var_test_target].second;
        if ((SourceP0 != -1) && (SourceP0 == TargetP0))
          continue;
        if ((SourceP0 != -1) && (SourceP0 == TargetP1))
          continue;
        if ((SourceP1 != -1) && (SourceP1 == TargetP0))
          continue;
        if ((SourceP1 != -1) && (SourceP1 == TargetP1))
          continue;
      }
      curr_target.push_back(test_target);
    }
    return curr_target;
  }

  std::vector<int> GetUnitaryNodes() {
    std::vector<int> UnitNodes;
    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int VarIndex = VarOfNode(i);
      assert(VarIndex >= 0);
      assert(VarIndex < VarValue.size());
      if (VarValue[VarIndex] == 1)
        UnitNodes.push_back(i);
    }
    return UnitNodes;
  }

  void DisableUnitaryNodes() {
    EnabledN = std::vector<bool>(GraphAdj.size(), true);
    std::vector<int> UniNodes = GetUnitaryNodes();
    for (size_t i = 0; i < UniNodes.size(); i++)
      EnabledN[UniNodes[i]] = false;
  }

  std::vector<ShortestPathResult> FindPaths(const std::vector<int> &SourceN,
                                            const std::vector<int> &TargetN,
                                            bool DisableUnitary,
                                            bool check_common_patch = false) {

    EnabledN = std::vector<bool>(GraphAdj.size(), true);
    if (DisableUnitary)
      DisableUnitaryNodes();

    std::vector<ShortestPathResult> resPaths;

    for (size_t i = 0; i < SourceN.size(); i++) {
      int curr_sourceN = SourceN[i];
      // std::cout<<"Node Idx:"<<curr_sourceN<<std::endl;
      assert(curr_sourceN >= 0);
      assert(curr_sourceN < Node_to_VarPatch.size());

      std::vector<int> curt_target =
          AdmissibleTargets(curr_sourceN, TargetN, check_common_patch);
      std::vector<int> curt_source(1, curr_sourceN);

      // std::cout << "Curr Source:" << curr_sourceN << std::endl;
      // std::cout << "Curr Target Size:" << curt_target.size() << std::endl;

      ShortestPathResult SRes = ShortestPath(curt_source, curt_target);

      if (!SRes.found)
        continue;

      CheckPath(SRes);
      resPaths.push_back(SRes);
    }
    return resPaths;
  }

  bool SolveEqStep() {
    bool HasToSolve = true;

    HasToSolve = SetupSolve(true, false);

    if (!HasToSolve)
      return false;

    TargetDecr.insert(TargetDecr.end(), BoundarySinks.begin(),
                      BoundarySinks.end());
    TargetIncr.insert(TargetIncr.end(), BoundarySinks.begin(),
                      BoundarySinks.end());

    std::vector<ShortestPathResult> IncreasePaths, DecreasePaths;

    for (size_t i = 0; i < SourceIncr.size(); i++) {
      int curr_sourceN = SourceIncr[i];
      assert(curr_sourceN >= 0);
      assert(curr_sourceN < Node_to_VarPatch.size());
    }

    SetupWeightChange(1, true);

    // std::cout << "*** FINDING INCREASE PATHS ***" << std::endl;
    IncreasePaths = FindPaths(SourceIncr, TargetIncr, false);
    // std::cout << "*** DONE ***" << std::endl;

    // find the best path increase
    int BestIdxIncr = -1;
    ScalarType BestWIncr = std::numeric_limits<ScalarType>::max();
    for (size_t i = 0; i < IncreasePaths.size(); i++) {

      for (size_t j = 0; j < IncreasePaths[i].path.size(); j++) {
        assert(EnabledN[IncreasePaths[i].path[j]]);
      }

      ScalarType CurrErr = Priority(IncreasePaths[i].path, 1);
      // IncreasePaths[i].length; // ModifyError(IncreasePaths[i], 1);
      if (CurrErr > BestWIncr)
        continue;
      BestWIncr = CurrErr;
      BestIdxIncr = i;
    }

    SetupWeightChange(-1, true);

    // std::cout << "*** FINDING DECREASE PATHS ***" << std::endl;
    DecreasePaths = FindPaths(SourceDecr, TargetDecr, true);
    // std::cout << "*** DONE ***" << std::endl;

    if ((IncreasePaths.size() == 0) && (DecreasePaths.size() == 0))
      return false;

    int BestIdxDecr = -1;
    ScalarType BestWDecr = std::numeric_limits<ScalarType>::max();
    for (size_t i = 0; i < DecreasePaths.size(); i++) {

      for (size_t j = 0; j < DecreasePaths[i].path.size(); j++) {
        assert(EnabledN[DecreasePaths[i].path[j]]);
      }
      if (!CanDecrease(DecreasePaths[i].path))
        continue;

      ScalarType CurrErr = Priority(DecreasePaths[i].path, -1);
      // DecreasePaths[i].length; // ModifyError(DecreasePaths[i], -1);
      if (CurrErr > BestWDecr)
        continue;
      BestWDecr = CurrErr;
      BestIdxDecr = i;
    }

    if ((BestIdxIncr == -1) && (BestIdxDecr == -1))
      return false;

    if (BestWDecr < BestWIncr) {
      // std::cout << "*** APPLYING Decrease ***" << std::endl;
      ApplyDecrease(DecreasePaths[BestIdxDecr]);
    } else {
      // std::cout << "*** APPLYING Increase ***" << std::endl;
      ApplyIncrease(IncreasePaths[BestIdxIncr]);
    }

    // std::cout << "*** DONE ***" << std::endl;

    return true;
  }

  bool SolveParityStep() {
    bool HasToSolve = true;

    HasToSolve = SetupSolve(false, true);

    if (!HasToSolve)
      return false;

    // add sinks
    TargetParity.insert(TargetParity.end(), BoundarySinks.begin(),
                        BoundarySinks.end());

    std::vector<ShortestPathResult> IncreasePaths, DecreasePaths;

    // check
    for (size_t i = 0; i < SourceParity.size(); i++) {
      int curr_sourceN = SourceParity[i];
      assert(curr_sourceN >= 0);
      assert(curr_sourceN < Node_to_VarPatch.size());
    }

    // find increase candidate
    SetupNodeWeightChange(1, true);
    IncreasePaths = FindPaths(SourceParity, TargetParity, false, true);

    // find the best path increase
    int BestIdxIncr = -1;
    ScalarType BestWIncr = std::numeric_limits<ScalarType>::max();
    for (size_t i = 0; i < IncreasePaths.size(); i++) {

      for (size_t j = 0; j < IncreasePaths[i].path.size(); j++) {
        assert(EnabledN[IncreasePaths[i].path[j]]);
      }

      ScalarType CurrErr = Priority(IncreasePaths[i].path, 1);
      if (CurrErr > BestWIncr)
        continue;
      BestWIncr = CurrErr;
      BestIdxIncr = i;
    }

    // find the best path decrease
    SetupWeightChange(-1, true);
    DecreasePaths = FindPaths(SourceDecr, TargetDecr, true, true);
    if ((IncreasePaths.size() == 0) && (DecreasePaths.size() == 0))
      return false;

    int BestIdxDecr = -1;
    ScalarType BestWDecr = std::numeric_limits<ScalarType>::max();
    for (size_t i = 0; i < DecreasePaths.size(); i++) {

      for (size_t j = 0; j < DecreasePaths[i].path.size(); j++) {
        assert(EnabledN[DecreasePaths[i].path[j]]);
      }
      if (!CanDecrease(DecreasePaths[i].path))
        continue;

      ScalarType CurrErr = Priority(DecreasePaths[i].path, -1);
      if (CurrErr > BestWDecr)
        continue;
      BestWDecr = CurrErr;
      BestIdxDecr = i;
    }

    if ((BestIdxIncr == -1) && (BestIdxDecr == -1))
      return false;

    if (BestWDecr < BestWIncr) {
      // std::cout << "*** APPLYING Decrease ***" << std::endl;
      ApplyDecrease(DecreasePaths[BestIdxDecr]);
    } else {
      // std::cout << "*** APPLYING Increase ***" << std::endl;
      ApplyIncrease(IncreasePaths[BestIdxIncr]);
    }

    // std::cout << "*** DONE ***" << std::endl;

    return true;
  }


  struct NodeRoute {
    enum Kind { NONE, CYCLE, SRC_TO_SINK } kind = NONE;
    ScalarType cost = std::numeric_limits<ScalarType>::infinity();
    std::vector<int> nodes; // cycles: [u, ..., u]
  };

  struct PerNodeResult {
    std::vector<std::pair<int, NodeRoute>> results; // (nodeId, NodeRoute)
    void Add(const PerNodeResult &other) {
      results.insert(results.end(), other.results.begin(), other.results.end());
    }
  };

  // ---------- Helpers ----------
  // Scan the parent-chain and tell if 'varToCheck' already appears (strict: no
  // exceptions)
  static inline bool
  appears_in_chain_var_strict(const std::vector<int> &parent,
                              const std::vector<int> &VarNode, int curNode,
                              int varToCheck) {
    for (int x = curNode; x != -1; x = parent[x]) {
      if (VarNode[x] == varToCheck)
        return true;
    }
    return false;
  }

  static std::vector<std::vector<int>>
  buildReverse(const std::vector<std::vector<int>> &G,
               const std::vector<bool> &EnabledN) {
    int n = (int)G.size();
    std::vector<std::vector<int>> R(n);
    for (int u = 0; u < n; ++u)
      if (EnabledN[u]) {
        for (int v : G[u])
          if (v >= 0 && v < n && EnabledN[v])
            R[v].push_back(u);
      }
    return R;
  }

  // static void findSourcesSinks(const std::vector<std::vector<int>> &G,
  //                              const std::vector<bool> &EnabledN,
  //                              std::vector<int> &sources,
  //                              std::vector<int> &sinks) {
  //   int n = (int)G.size();
  //   std::vector<int> indeg(n, 0), outdeg(n, 0);
  //   for (int u = 0; u < n; ++u)
  //     if (EnabledN[u]) {
  //       for (int v : G[u])
  //         if (v >= 0 && v < n && EnabledN[v]) {
  //           ++indeg[v];
  //           ++outdeg[u];
  //         }
  //     }
  //   for (int u = 0; u < n; ++u)
  //     if (EnabledN[u]) {
  //       if (indeg[u] == 0)
  //         sources.push_back(u);
  //       if (outdeg[u] == 0)
  //         sinks.push_back(u);
  //     }
  // }

  static void findSourcesSinks(const std::vector<std::vector<int>> &G,
                               const std::vector<bool> &EnabledN,
                               std::vector<int> &sources,
                               std::vector<int> &sinks) {
    int n = (int)G.size();
    std::vector<int> indeg(n, 0), outdeg(n, 0);
    for (int u = 0; u < n; ++u)
      for (int v : G[u]) {
        assert(v >= 0);
        assert(v < indeg.size());
        assert(v < outdeg.size());
        ++indeg[v];
        ++outdeg[u];
      }
    for (int u = 0; u < n; ++u)
      if (EnabledN[u]) {
        if (indeg[u] == 0)
          sources.push_back(u);
        if (outdeg[u] == 0)
          sinks.push_back(u);
      }
  }

  static std::vector<int> rebuild_path_src_to(int t,
                                              const std::vector<int> &parent) {
    std::vector<int> path;
    for (int x = t; x != -1; x = parent[x])
      path.push_back(x);
    std::reverse(path.begin(), path.end());
    return path;
  }

  static std::vector<int>
  rebuild_cycle_u_v_u(int u, int v, const std::vector<int> &parentRev) {
    std::vector<int> tail; // v ... u following parentRev toward u
    for (int x = v; x != -1; x = parentRev[x]) {
      tail.push_back(x);
      if (x == u)
        break;
    }
    std::vector<int> cyc;
    cyc.reserve(tail.size() + 1);
    cyc.push_back(u);
    cyc.insert(cyc.end(), tail.begin(), tail.end());
    if (cyc.back() != u)
      cyc.push_back(u);
    return cyc;
  }

  // Uniqueness on a path: no variable repeats anywhere.
  static bool path_vars_unique(const std::vector<int> &seq,
                               const std::vector<int> &VarNode) {
    // O(n^2) is fine for short paths; for long paths use a bitset/hashset.
    for (int i = 0; i < (int)seq.size(); ++i) {
      int vi = VarNode[seq[i]];
      for (int j = i + 1; j < (int)seq.size(); ++j) {
        if (VarNode[seq[j]] == vi)
          return false;
      }
    }
    return true;
  }

  // Uniqueness on a cycle [u ... u]: allow only the closing u duplication.
  static bool cycle_vars_unique_allow_last_u(const std::vector<int> &cyc,
                                             const std::vector<int> &VarNode) {
    if (cyc.size() >= 2 && cyc.front() == cyc.back()) {
      // check [0 .. size-2]
      for (int i = 0; i + 1 < (int)cyc.size(); ++i) {
        int vi = VarNode[cyc[i]];
        for (int j = i + 1; j + 1 < (int)cyc.size(); ++j) {
          if (VarNode[cyc[j]] == vi)
            return false;
        }
      }
      return true;
    } else {
      return path_vars_unique(cyc, VarNode);
    }
  }

  // ---------- Constrained Dijkstra variants (STRICT uniqueness: no exceptions)
  // ---------- Forward, multi-source -> target (stop when target is settled).
  // Edge cost a->b = Wnode[b]; include Wnode[source] at init.
  static bool constrained_multi_source_to_target(
      const std::vector<std::vector<int>> &G,
      const std::vector<ScalarType> &Wnode, // >= 0
      const std::vector<bool> &EnabledN, const std::vector<int> &VarNode,
      const std::vector<int> &sources, int target,
      std::vector<ScalarType> &dist, // out
      std::vector<int> &parent       // out
  ) {
    const int n = (int)G.size();
    dist.assign(n, std::numeric_limits<ScalarType>::infinity());
    parent.assign(n, -1);

    using Item = std::pair<ScalarType, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    for (int s : sources)
      if (s >= 0 && s < n && EnabledN[s]) {
        dist[s] = Wnode[s];
        parent[s] = -1;
        pq.emplace(dist[s], s);
      }

    while (!pq.empty()) {
      auto [du, u] = pq.top();
      pq.pop();
      if (du > dist[u])
        continue;
      if (u == target)
        return true;

      for (int v : G[u])
        if (v >= 0 && v < n && EnabledN[v]) {
          int varV = VarNode[v];
          if (appears_in_chain_var_strict(parent, VarNode, u, varV))
            continue; // strict: no repeats allowed

          ScalarType nd = du + Wnode[v];
          if (nd < dist[v]) {
            dist[v] = nd;
            parent[v] = u;
            pq.emplace(nd, v);
          }
        }
    }
    return false;
  }

  // Forward, single-source u -> ANY sink (cheapest). Uniqueness strict.
  static int constrained_single_source_to_any_sink(
      const std::vector<std::vector<int>> &G,
      const std::vector<ScalarType> &Wnode, // >= 0
      const std::vector<bool> &EnabledN, const std::vector<int> &VarNode, int u,
      const std::vector<int> &sinks,
      std::vector<ScalarType> &dist, // out
      std::vector<int> &parent       // out
  ) {
    const int n = (int)G.size();
    dist.assign(n, std::numeric_limits<ScalarType>::infinity());
    parent.assign(n, -1);

    if (!(u >= 0 && u < n && EnabledN[u]))
      return -1;

    // (Optional) speed-up: premark sinks to avoid scanning the vector each
    // time.
    std::vector<char> isSink(n, 0);
    for (int t : sinks)
      if (t >= 0 && t < n && EnabledN[t])
        isSink[t] = 1;

    using Item = std::pair<ScalarType, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    dist[u] = Wnode[u];
    pq.emplace(dist[u], u);

    while (!pq.empty()) {
      auto [du, x] = pq.top();
      pq.pop();
      if (du > dist[x])
        continue;

      if (isSink[x])
        return x; // cheapest sink reached

      for (int v : G[x])
        if (EnabledN[v]) {
          int varV = VarNode[v];
          if (appears_in_chain_var_strict(parent, VarNode, x, varV))
            continue;

          ScalarType nd = du + Wnode[v];
          if (nd < dist[v]) {
            dist[v] = nd;
            parent[v] = x;
            pq.emplace(nd, v);
          }
        }
    }
    return -1;
  }

  // Reverse, single-source from u on reverse graph (for cycle through u).
  // Uniqueness strict (no exceptions).
  static void constrained_reverse_from_u(
      const std::vector<std::vector<int>> &R, // reverse graph
      const std::vector<ScalarType> &Wnode,   // >= 0
      const std::vector<bool> &EnabledN, const std::vector<int> &VarNode, int u,
      std::vector<ScalarType> &distRev, // out
      std::vector<int> &parentRev       // out
  ) {
    const int n = (int)R.size();
    distRev.assign(n, std::numeric_limits<ScalarType>::infinity());
    parentRev.assign(n, -1);

    if (!(u >= 0 && u < n && EnabledN[u]))
      return;

    using Item = std::pair<ScalarType, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    distRev[u] = Wnode[u]; // cost to "close" at u includes W[u]
    pq.emplace(distRev[u], u);

    while (!pq.empty()) {
      auto [du, x] = pq.top();
      pq.pop();
      if (du > distRev[x])
        continue;

      for (int p : R[x])
        if (EnabledN[p]) {
          int varP = VarNode[p];
          if (appears_in_chain_var_strict(parentRev, VarNode, x, varP))
            continue;

          ScalarType nd = du + Wnode[p];
          if (nd < distRev[p]) {
            distRev[p] = nd;
            parentRev[p] = x;
            pq.emplace(nd, p);
          }
        }
    }
  }

  // ---------- Main API ----------
  // For each u in routeNodes: take min between
  //  (1) shortest cycle through u    (strict var-uniqueness; only closing u is
  //  duplicated) (2) shortest src→…→u→…→sink     (strict var-uniqueness on both
  //  halves and merged path)
  // If either part missing or uniqueness fails, candidate is discarded.
  PerNodeResult computePerNodeShortestRouteVarUnique(
      const std::vector<int> &VarNode, const std::vector<int> &routeNodes,
      bool incremental_filter = true // if true, skip later u's that share a
                                     // variable with an accepted best.nodes
  ) {
    const int n = (int)GraphAdj.size();
    assert((int)GraphNodeWeight.size() == n);
    assert((int)EnabledN.size() == n);
    assert((int)VarNode.size() == n);

    // Precompute reverse graph, sources and sinks
    auto RevAdj = buildReverse(GraphAdj, EnabledN);
    std::vector<int> sources, sinks;
    findSourcesSinks(GraphAdj, EnabledN, sources, sinks);

    PerNodeResult out;
    out.results.reserve(routeNodes.size());

    // Reusable buffers
    std::vector<ScalarType> dist, distRev, distL, distR;
    std::vector<int> parent, parentRev, parentL, parentR;

    // Optional: incremental filter so we don't solve two u's sharing the same
    // variable
    std::unordered_set<int> exploredVar;
    exploredVar.reserve(routeNodes.size() * 2);

    for (int u : routeNodes) {
      NodeRoute best; // NONE
      if (u < 0 || u >= n || !EnabledN[u]) {
        out.results.emplace_back(u, best);
        continue;
      }

      int varU = VarNode[u];
      if (incremental_filter && exploredVar.count(varU)) {
        // skip this u because a previous accepted route already used its
        // variable
        out.results.emplace_back(u, best);
        continue;
      }

      // ----- Candidate 1: cycle through u (strict uniqueness) -----
      constrained_reverse_from_u(RevAdj, GraphNodeWeight, EnabledN, VarNode, u,
                                 distRev, parentRev);

      ScalarType bestCycleCost = std::numeric_limits<ScalarType>::infinity();
      int bestV = -1;
      for (int v : GraphAdj[u])
        if (v >= 0 && v < n && EnabledN[v]) {
          ScalarType cand =
              GraphNodeWeight[v] + distRev[v]; // enter v + (v..u inclusive)
          if (cand < bestCycleCost) {
            bestCycleCost = cand;
            bestV = v;
          }
        }
      if (bestV != -1 && std::isfinite(bestCycleCost)) {
        auto cyc = rebuild_cycle_u_v_u(u, bestV, parentRev); // [u ... u]
        if (cycle_vars_unique_allow_last_u(cyc, VarNode)) {
          best.kind = NodeRoute::CYCLE;
          best.cost = bestCycleCost;
          best.nodes = std::move(cyc);
        }
      }

      // ----- Candidate 2: src→…→u and u→…→sink (strict uniqueness on both
      // halves) -----
      bool hasLeft = constrained_multi_source_to_target(
          GraphAdj, GraphNodeWeight, EnabledN, VarNode, sources, u, distL,
          parentL);

      int sink = constrained_single_source_to_any_sink(
          GraphAdj, GraphNodeWeight, EnabledN, VarNode, u, sinks, distR,
          parentR);

      if (hasLeft && sink != -1) {
        auto left = rebuild_path_src_to(u, parentL);     // [src ... u]
        auto right = rebuild_path_src_to(sink, parentR); // [u ... sink]
        if (!right.empty() && !left.empty() && right.front() == u)
          right.erase(right.begin()); // avoid duplicating 'u'

        std::vector<int> merged = left;
        merged.insert(merged.end(), right.begin(), right.end());

        if (path_vars_unique(merged, VarNode)) {
          ScalarType cost = 0;
          for (int v : merged)
            cost += GraphNodeWeight[v];
          if (cost < best.cost) {
            best.kind = NodeRoute::SRC_TO_SINK;
            best.cost = cost;
            best.nodes = std::move(merged);
          }
        }
      }

      if (incremental_filter && !best.nodes.empty()) {
        for (int x : best.nodes)
          exploredVar.insert(VarNode[x]);
      }

      out.results.emplace_back(u, std::move(best));
    }

    return out;
  }

  ScalarType Priority(const std::vector<int> &NodeSeq, int Diff) {

    if (NodeSeq.size() == 0)
      return 0;

    std::vector<int> NodeIdx = NodeSeq;
    std::sort(NodeIdx.begin(), NodeIdx.end());
    NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()), NodeIdx.end());

    ScalarType Prio = 0;
    // checkPrio=false;
    for (size_t i = 0; i < NodeIdx.size(); i++) {
      assert(EnabledN[NodeIdx[i]]);
      ScalarType CurrPrio = PriorityChange(NodeIdx[i], Diff, true);
      Prio += PriorityChange(NodeIdx[i], Diff);
    }
    return Prio;
  }

  ScalarType Priority(const NodeRoute &NodeR, int Diff) {
    return (Priority(NodeR.nodes, Diff));
  }

  std::vector<int> GetStartingNode(bool check_unitary = false) {
    std::vector<int> StartN;
    std::set<int> CheckedVars;
    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int IndexVar = VarOfNode(i);
      if (VarValue[IndexVar] == VarTargetValue[IndexVar])
        continue;
      if (check_unitary && (VarValue[IndexVar] <= 1))
        continue;
      if (CheckedVars.count(IndexVar) > 0)
        continue;

      CheckedVars.insert(IndexVar);
      StartN.push_back(i);
    }
    return StartN;
  }

  bool UniqueVariable(const std::vector<int> &NodeIdx) {
    std::set<int> Variables;
    for (size_t i = 0; i < NodeIdx.size(); i++) {
      int IndexVar = VarOfNode(NodeIdx[i]);
      if (Variables.count(IndexVar) > 0)
        return false;
      Variables.insert(IndexVar);
    }
    return true;
  }

  void CleanCyclePaths(PerNodeResult &CycleRes) {
    // clean the cycle so that one
    // variable is inremented
    // only once
    for (size_t i = 0; i < CycleRes.results.size(); i++) {
      NodeRoute CurrR = CycleRes.results[i].second;
      assert(CurrR.kind != NodeRoute::SRC_TO_SINK);
      if (CurrR.kind == NodeRoute::CYCLE) {
        assert(CurrR.nodes[0] == CurrR.nodes.back());
        CurrR.nodes.pop_back();
      }
      CycleRes.results[i].second = CurrR;
      assert(UniqueVariable(CurrR.nodes));
    }
  }

  bool ImproveSizeStep() {

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    bool HasToSolve = true;

    // set all start nodes
    std::vector<int> StartNIncr = GetStartingNode();
    std::vector<int> StartNDecr = GetStartingNode(true);
    // StartNIncr=FilterIncrease(StartNIncr);
    // StartNDecr=FilteDecrease(StartNDecr);

    std::vector<int> VarNode;
    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int IndexVar = VarOfNode(i);
      VarNode.push_back(IndexVar);
    }
    // get the increase paths
    EnabledN = std::vector<bool>(GraphAdj.size(), true);
    SetupNodeWeightChange(1, true);

    std::cout << "*** Computing Increment Resizing Cycles.." << std::endl;
    PerNodeResult IncrRes = computePerNodeShortestRouteVarUnique(
        VarNode, StartNIncr); // computePerNodeShortestRoute(StartNIncr);
    std::cout << "Computed " << IncrRes.results.size() << " increments Paths"
              << std::endl;
    // get the decrease paths
    EnabledN = std::vector<bool>(GraphAdj.size(), true);
    DisableUnitaryNodes();
    SetupNodeWeightChange(-1, true);
    // PerNodeResult DecrRes = computePerNodeShortestRoute(StartNDecr);

    std::cout << "*** Computing Decrement Resizing Cycles.." << std::endl;
    PerNodeResult DecrRes =
        computePerNodeShortestRouteVarUnique(VarNode, StartNDecr);
    std::cout << "Computed " << DecrRes.results.size() << " decrement Paths"
              << std::endl;
    PerNodeResult AllRes = IncrRes;
    AllRes.Add(DecrRes);

    CleanCyclePaths(AllRes);

    bool modified = false;
    bool one_step = false;
    EnabledN = std::vector<bool>(GraphAdj.size(), true);
    int steps = 0;
    do {
      // find the best
      int BestModIdx = -1;
      ScalarType BestPrio = 0;
      bool Increase = true;
      modified = false;
      for (size_t i = 0; i < AllRes.results.size(); i++) {
        std::vector<int> NodeIdx = AllRes.results[i].second.nodes;
        // CheckUniqueVariable(NodeIdx);
        // std::sort(NodeIdx.begin(), NodeIdx.end());
        // NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()),
        //               NodeIdx.end());

        ScalarType CurrPrioAdd = Priority(NodeIdx, 1);
        ScalarType CurrPrioDec = 0;

        if (CanDecrease(NodeIdx)) {
          CurrPrioDec = Priority(NodeIdx, -1);
        }

        if (CurrPrioAdd < BestPrio) {
          BestPrio = CurrPrioAdd;
          Increase = true;
          BestModIdx = i;
        }

        if (CurrPrioDec < BestPrio) {
          BestPrio = CurrPrioDec;
          Increase = false;
          BestModIdx = i;
        }
      }

      if (BestModIdx != -1) {
        steps++;

        std::vector<int> NodeIdx = AllRes.results[BestModIdx].second.nodes;

        if (Increase) {
          ApplyIncrease(NodeIdx);
          for (size_t i = 0; i < VarValue.size(); i++) {
            assert(VarValue[i] >= 1);
          }
        } else {
          assert(CanDecrease(NodeIdx));
          ApplyDecrease(NodeIdx);

          for (size_t i = 0; i < VarValue.size(); i++) {
            assert(VarValue[i] >= 1);
          }
        }
        modified = true;
        one_step = true;
      }

    } while (modified);

    return one_step;
  }

  ScalarType CurrentSizeEnergy() {
    ScalarType CurrE = 0;
    for (size_t i = 0; i < VarValue.size(); i++)
      CurrE += Energy(i);
    return CurrE;
  }

public:
  void SolveEquality() {

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    InitGraph();
    std::cout << "*** Solving Equality ***" << std::endl;
    bool changed = false;
    int steps = 0;
    do {
      changed = SolveEqStep();
      steps++;
      // if (WriteDebug)
      //   std::cout << "Performed Equality size step" << std::endl;
    } while (changed);
    std::cout << "Performed " << steps << " Equality size step" << std::endl;

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }
  }

  void SolveParity() {

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    InitGraph();
    bool changed = false;

    std::cout << "*** Solving Parity ***" << std::endl;
    int steps = 0;
    do {
      changed = SolveParityStep();
      steps++;
      // if (WriteDebug)
      //   std::cout << "Performed Parity step" << std::endl;
    } while (changed);
    std::cout << "Performed " << steps << " Parity size step" << std::endl;

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }
  }

  void ImproveSize() {

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    InitGraph();

    std::vector<ScalarType> GraphNodeW(GraphAdj.size(), 1);
    bool HasImproved = false;
    if (WriteDebug)
      std::cout << "Energy Value:" << CurrentSizeEnergy() << std::endl;
    do {
      HasImproved = ImproveSizeStep();
      // HasImproved = ImproveSizeStep2();
      if (WriteDebug) {
        if (HasImproved) {
          std::cout << "Performed Improve size step" << std::endl;
          std::cout << "Energy Value:" << CurrentSizeEnergy() << std::endl;
        }
      }
    } while (HasImproved);

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }
  }

  OrientFLowGraph(
      std::vector<int> &_VarValue, std::vector<int> &_VarTargetValue,
      const std::vector<std::vector<std::vector<int>>> &_PatchSideVariables)
      : VarTargetValue(_VarTargetValue), VarValue(_VarValue),
        PatchSideVariables(_PatchSideVariables) {
    WriteDebug = true;
  }
};

template <class ScalarType> class ConstraintGraph {

  struct VarConstrL {
    int IndexC0;
    int IndexC1;
    VarConstrL() {
      IndexC0 = -1;
      IndexC1 = -1;
    }

    void AddLink(const int &IndexConstr) {
      if (IndexC0 == -1) {
        IndexC0 = IndexConstr;
        return;
      }
      if (IndexC1 == -1) {
        IndexC1 = IndexConstr;
        return;
      }
      assert(0);
    }
  };

  struct AdjConstr {
    int AdjC;
    int IndexVar;

    AdjConstr(int _AdjC, int _IndexVar) {
      AdjC = _AdjC;
      IndexVar = _IndexVar;
    }
  };

  // patch parity constraints and graph
  std::vector<std::vector<int>> ParityConstr;
  std::vector<std::vector<AdjConstr>> ParityConstraintsAdj;
  std::vector<VarConstrL> ParityConstrLinks;
  std::vector<int> NotParityNodes;

  // opposite Equality constraints and graph
  std::vector<std::pair<std::vector<int>, std::vector<int>>> OppEqConstr;
  std::vector<VarConstrL> OppEqConstrConstrLinks;
  std::vector<int> IncreaseNodes, DereaseNodes;
  std::vector<std::vector<AdjConstr>> EqConstraintsAdj;

  // std::vector<std::vector<AdjConstr>> EqConstraintsIncrAdj,
  //     EqConstraintsDecrAdj;
  // std::vector<VarConstrL> OppEqConstrConstrLinksIncr,
  //     OppEqConstrConstrLinksDecr;
  std::vector<int> NotOppEqNodes;
  // std::set<int> IncreaseVar, DecreaseVar;

  std::vector<int> &VarTargetValue;
  std::vector<int> &VarValue;

  void InitAdjFromLinks(std::vector<std::vector<AdjConstr>> &ConstraintsAdj,
                        const std::vector<VarConstrL> &ConstrLinks) {
    // then set adjacency between constraints
    ConstraintsAdj.clear();
    ConstraintsAdj.resize(ConstrLinks.size());
    for (size_t i = 0; i < ConstrLinks.size(); i++) {
      int IndexC0 = ConstrLinks[i].IndexC0;
      int IndexC1 = ConstrLinks[i].IndexC1;
      int IndexVar = i;
      if ((IndexC0 == -1) && (IndexC1 == -1))
        continue;

      if (IndexC1 == -1) {
        assert(IndexC0 != -1);
        assert(IndexC0 < ConstraintsAdj.size());
        ConstraintsAdj[IndexC0].push_back(AdjConstr(-1, IndexVar));
        continue;
      }
      if (IndexC0 == -1) {
        assert(IndexC1 != -1);
        assert(IndexC1 < ConstraintsAdj.size());
        ConstraintsAdj[IndexC1].push_back(AdjConstr(-1, IndexVar));
        continue;
      }
      assert(IndexC0 != -1);
      assert(IndexC0 < ConstraintsAdj.size());
      assert(IndexC1 != -1);
      assert(IndexC1 < ConstraintsAdj.size());
      ConstraintsAdj[IndexC0].push_back(AdjConstr(IndexC1, IndexVar));
      ConstraintsAdj[IndexC1].push_back(AdjConstr(IndexC0, IndexVar));
    }
  }

  void InitParityConstr(
      std::vector<std::vector<std::vector<int>>> &_PatchSideVariables) {
    ParityConstr.clear();
    for (size_t i = 0; i < _PatchSideVariables.size(); i++) {
      std::vector<int> PatchVars;
      for (size_t j = 0; j < _PatchSideVariables[i].size(); j++)
        for (size_t k = 0; k < _PatchSideVariables[i][j].size(); k++)
          PatchVars.push_back(_PatchSideVariables[i][j][k]);
      assert(PatchVars.size() >= 4);
      ParityConstr.push_back(PatchVars);
    }

    NotParityNodes.clear();
    for (size_t i = 0; i < ParityConstr.size(); i++) {
      size_t CurrSum = 0;
      for (size_t j = 0; j < ParityConstr[i].size(); j++) {
        int IndexVar = ParityConstr[i][j];
        assert(IndexVar >= 0);
        assert(IndexVar < VarValue.size());
        CurrSum += VarValue[IndexVar];
      }
      if (i == 82)
        std::cout << std::endl;
      if ((CurrSum % 2) == 1)
        NotParityNodes.push_back(i);
    }
    std::cout << "Num Not Parity Constraints " << NotParityNodes.size()
              << std::endl;
    // allocate all possible variables and links
    ParityConstrLinks.resize(VarValue.size());
    for (size_t i = 0; i < ParityConstr.size(); i++) {
      for (size_t j = 0; j < ParityConstr[i].size(); j++) {
        {
          int IndexVar = ParityConstr[i][j];
          int IndexConstr = i;
          assert(IndexVar >= 0);
          assert(IndexVar < ParityConstrLinks.size());
          ParityConstrLinks[IndexVar].AddLink(IndexConstr);
        }
      }
    }

    InitAdjFromLinks(ParityConstraintsAdj, ParityConstrLinks);
    std::cout << "Allocated Graph " << std::endl;
  }

  // BSF COMPUTATION

  struct BFSData {
    std::vector<int> dist;       // BFS distance in number of internal edges
    std::vector<int> parentNode; // previous constraint on the BFS tree
    std::vector<int>
        parentVar; // IndexVar used to reach this node from parentNode
  };

  // BFS over constraints graph using ONLY internal edges (AdjC != -1).
  BFSData bfs_from(const std::vector<std::vector<AdjConstr>> &ConstraintGraph,
                   int src) {
    const int n = (int)ConstraintGraph.size();
    BFSData out;
    out.dist.assign(n, INT_MAX);
    out.parentNode.assign(n, -1);
    out.parentVar.assign(n, -1);

    std::queue<int> q;
    out.dist[src] = 0;
    q.push(src);

    while (!q.empty()) {
      int u = q.front();
      q.pop();
      for (const auto &e : ConstraintGraph[u]) {
        if (e.AdjC == -1)
          continue; // skip leaves: they don't connect constraints
        int v = e.AdjC;
        if (out.dist[v] != INT_MAX)
          continue; // already visited
        out.dist[v] = out.dist[u] + 1;
        out.parentNode[v] = u;
        out.parentVar[v] = e.IndexVar; // variable that connects u -> v
        q.push(v);
      }
    }
    return out;
  }

  // Reconstruct IndexVar sequence along the shortest path src->dst using BFS
  // parents.
  std::vector<int> reconstruct_path_vars(const BFSData &bfs, int src, int dst) {
    std::vector<int> vars;
    if (bfs.dist[dst] == INT_MAX)
      return vars; // unreachable
    int v = dst;
    while (v != src) {
      int pv = bfs.parentVar[v];
      if (pv < 0) {
        vars.clear();
        return vars;
      } // safety
      vars.push_back(pv);
      v = bfs.parentNode[v];
    }
    reverse(vars.begin(), vars.end());
    return vars;
  }

  struct PCJoinResult {
    bool ok = true; // DP ran; true even if nothing could be fixed
    std::vector<int>
        varsToFlip; // IndexVar to flip (XOR: each appears at most once)
    std::vector<int>
        fixedConstraints; // constraints from NotEvenNodes that were fixed
    std::vector<int>
        unresolvedConstraints; // constraints from NotEvenNodes left odd
    int fixedCount = 0;        // how many constraints we fixed
    int matchedPairs = 0;      // how many pairs used
    int leafUsed = 0;          // how many leaves used
    int totalPathLength = 0;   // total number of internal edges flipped
    int totalCost = 0; // total flips (path edges + leaves), uniform weights
  };

  // Helper to compare (coveredNodes, cost): maximize coveredNodes, tie-break
  // by min cost
  static inline bool better_pair(int curCover, int curCost, int newCover,
                                 int newCost) {
    return (newCover > curCover) || (newCover == curCover && newCost < curCost);
  }

  // Parent action used to backtrack the DP
  struct ParentAct {
    enum Kind { NONE, LEAF, PAIR } kind = NONE;
    int i = -1, j = -1; // indices in T (for LEAF only i is used)
  };

  /**
   * Prize-collecting T-join (uniform weights) with partial coverage.
   * - Uses internal edges (AdjC>=0) as unit-cost connectors between
   * constraints.
   * - Uses leaves (AdjC==-1), if present on a constraint, as unit-cost single
   * fixes. Goal: maximize how many odd constraints (in NotEvenNodes) are
   * fixed; among those, minimize total number of flips (path length +
   * leaves).
   *
   * Inputs:
   *  - Constraints: graph on constraints with AdjConstr entries.
   *  - NotEvenNodes: list of odd constraints (need fixing).
   *
   * Output:
   *  - PCJoinResult: varsToFlip (IndexVar, XORed), fixed/unresolved
   * constraints, stats.
   *
   * Notes:
   *  - Complexity of the DP is O(2^|T| * |T|^2). Practical for |T| ≲ 22–24.
   *  - Distances are computed with BFS (unit weights).
   */
  PCJoinResult prize_collecting_tjoin_maxcover(
      const std::vector<std::vector<AdjConstr>> &Constraints,
      const std::vector<int> &ToPairNodes) {
    PCJoinResult res;
    const int k = (int)ToPairNodes.size();
    if (k == 0)
      return res; // nothing to do

    // --- Collect one leaf per constraint (uniform cost = 1) ---
    // leafVar[c] = IndexVar of some leaf on constraint c, or -1 if none.
    std::vector<int> leafVar(Constraints.size(), -1);
    for (int c = 0; c < (int)Constraints.size(); ++c) {
      for (const auto &e : Constraints[c]) {
        if (e.AdjC == -1) {
          leafVar[c] = e.IndexVar;
          break;
        }
      }
    }

    // --- BFS from each odd node to build metric on T ---
    std::vector<BFSData> bfsList;
    bfsList.reserve(k);
    for (int i = 0; i < k; ++i) {
      bfsList.push_back(bfs_from(Constraints, ToPairNodes[i]));
    }

    // --- Distance matrix among odd nodes (uniform weights) ---
    const int INF = INT_MAX / 4;
    std::vector<std::vector<int>> distT(k, std::vector<int>(k, INF));
    for (int i = 0; i < k; ++i) {
      distT[i][i] = 0;
      for (int j = i + 1; j < k; ++j) {
        int u = ToPairNodes[i], v = ToPairNodes[j];
        int d = bfsList[i].dist[v];
        if (d == INT_MAX)
          d = INF; // unreachable
        distT[i][j] = distT[j][i] = d;
      }
    }

    // --- DP over subsets (64-bit masks) ---
    // dp[mask]   = minimum cost to fix exactly nodes in 'mask'
    // cover[mask]= how many nodes are fixed in 'mask' (used as primary
    // objective) parent[mask]= how we reached this mask (PAIR i-j or LEAF i)
    using u64 = uint64_t;
    const uint32_t k32 = (uint32_t)k;

    // Safety guard: DP size = 2^k
    if (k32 >= 63) {
      // Too large for 64-bit masks / memory; fallback policy:
      // here we just return "nothing fixed". You can plug a Blossom fallback
      // if needed.
      res.ok = false;
      res.fixedCount = 0;
      res.totalCost = 0;
      res.unresolvedConstraints = ToPairNodes;
      return res;
    }

    const u64 FULL = (k32 == 64 ? ~0ULL : ((1ULL << k32) - 1ULL));
    const u64 DP_N = FULL + 1ULL;

    std::vector<int> dp(DP_N, INF);
    std::vector<int> cover(DP_N, -1);
    std::vector<ParentAct> parent(DP_N);

    dp[0] = 0;
    cover[0] = 0;

    for (u64 mask = 0; mask <= FULL; ++mask) {
      if (cover[mask] < 0)
        continue; // unreachable state

      // Find first free i (odd node not yet fixed in this mask)
      uint32_t i = 0;
      while (i < k32 && (mask & (1ULL << i)))
        ++i;
      if (i == k32)
        continue; // nothing left to add

      // Option A: fix i alone with a leaf (if available)
      {
        int ci = ToPairNodes[i];
        if (leafVar[ci] != -1) {
          u64 nmask = mask | (1ULL << i);
          int ncover = cover[mask] + 1;
          int ncost = dp[mask] + 1; // uniform leaf cost
          if (cover[nmask] < 0 ||
              better_pair(cover[nmask], dp[nmask], ncover, ncost)) {
            cover[nmask] = ncover;
            dp[nmask] = ncost;
            parent[nmask] = {ParentAct::LEAF, (int)i, -1};
          }
        }
      }

      // Option B: pair i with any j>i not yet fixed
      for (uint32_t j = i + 1; j < k32; ++j)
        if (!(mask & (1ULL << j))) {
          if (distT[i][j] >= INF)
            continue; // unreachable pair
          u64 nmask = mask | (1ULL << i) | (1ULL << j);
          int ncover = cover[mask] + 2;
          int ncost = dp[mask] + distT[i][j];
          if (cover[nmask] < 0 ||
              better_pair(cover[nmask], dp[nmask], ncover, ncost)) {
            cover[nmask] = ncover;
            dp[nmask] = ncost;
            parent[nmask] = {ParentAct::PAIR, (int)i, (int)j};
          }
        }
    }

    // --- Choose best mask: maximize cover, then minimize cost ---
    u64 bestMask = 0;
    int bestCover = -1, bestCost = INF;
    bool found = false;
    for (u64 mask = 0; mask <= FULL; ++mask) {
      if (cover[mask] < 0)
        continue;
      if (!found || better_pair(bestCover, bestCost, cover[mask], dp[mask])) {
        found = true;
        bestMask = mask;
        bestCover = cover[mask];
        bestCost = dp[mask];
      }
    }
    if (!found) {
      res.ok = false;
      return res;
    }

    // --- Backtrack actions to obtain pairs and leaves ---
    std::vector<std::pair<int, int>> pairs; // indices in T
    std::vector<int> leaves;                // indices in T
    u64 mask = bestMask;
    while (mask) {
      ParentAct act = parent[mask];
      if (act.kind == ParentAct::PAIR) {
        pairs.push_back({act.i, act.j});
        mask &= ~(1ULL << act.i);
        mask &= ~(1ULL << act.j);
      } else if (act.kind == ParentAct::LEAF) {
        leaves.push_back(act.i);
        mask &= ~(1ULL << act.i);
      } else {
        // Should not happen; break to avoid infinite loop
        break;
      }
    }

    // --- XOR selected path variables and leaf variables ---
    std::unordered_set<int> toggle;
    toggle.reserve(1024);
    auto toggle_var = [&](int varId) {
      if (varId < 0)
        return;
      auto it = toggle.find(varId);
      if (it == toggle.end())
        toggle.insert(varId);
      else
        toggle.erase(it); // XOR: twice cancels
    };

    int totalLen = 0, usedLeaves = 0;
    for (auto [ii, jj] : pairs) {
      int u = ToPairNodes[ii];
      int v = ToPairNodes[jj];
      const BFSData &bfsU = bfsList[ii]; // BFS rooted at u
      std::vector<int> pathVars = reconstruct_path_vars(bfsU, u, v);
      for (int varId : pathVars)
        toggle_var(varId);
      totalLen += (int)pathVars.size();
    }
    for (int ii : leaves) {
      int c = ToPairNodes[ii];
      int leaf = leafVar[c];
      if (leaf != -1) {
        toggle_var(leaf);
        usedLeaves++;
      }
    }

    // --- Build fixed / unresolved lists ---
    std::vector<char> fixed(k, 0);
    for (auto [ii, jj] : pairs) {
      fixed[ii] = fixed[jj] = 1;
    }
    for (int ii : leaves)
      fixed[ii] = 1;

    std::vector<int> fixedCons, unresolved;
    fixedCons.reserve(k);
    unresolved.reserve(k);
    for (int t = 0; t < k; ++t) {
      int c = ToPairNodes[t];
      if (fixed[t])
        fixedCons.push_back(c);
      else
        unresolved.push_back(c);
    }

    // --- Fill result ---
    res.ok = true;
    res.fixedCount = (int)fixedCons.size();
    res.matchedPairs = (int)pairs.size();
    res.leafUsed = usedLeaves;
    res.totalPathLength = totalLen;
    res.totalCost = totalLen + usedLeaves;
    res.fixedConstraints = std::move(fixedCons);
    res.unresolvedConstraints = std::move(unresolved);
    res.varsToFlip.assign(toggle.begin(), toggle.end());
    std::sort(res.varsToFlip.begin(), res.varsToFlip.end()); // deterministic
    return res;
  }

  struct HeuResult {
    std::vector<int> varsToFlip;            // IndexVar (unique, XORed)
    std::vector<int> fixedConstraints;      // odd constraints that became even
    std::vector<int> unresolvedConstraints; // odd constraints left odd
    int usedLeaves = 0;
    int usedPairs = 0;
    int totalCost = 0; // = usedLeaves + total internal-edge flips
  };

  /**
   * Greedy partial solver using both leaves and internal edges.
   * See explanation in the message.
   */
  HeuResult
  heuristic_maxcover(const std::vector<std::vector<AdjConstr>> &Constraints,
                     const std::vector<int> &NotEvenNodes) {
    HeuResult R;
    const int nC = (int)Constraints.size();

    // One leaf per constraint if any (uniform cost = 1)
    std::vector<int> leafVar(nC, -1);
    for (int c = 0; c < nC; ++c)
      for (auto &e : Constraints[c])
        if (e.AdjC == -1) {
          leafVar[c] = e.IndexVar;
          break;
        }

    // Connected components over internal edges
    std::vector<int> comp(nC, -1);
    int compCnt = 0;
    for (int s = 0; s < nC; ++s)
      if (comp[s] == -1) {
        std::queue<int> q;
        q.push(s);
        comp[s] = compCnt;
        while (!q.empty()) {
          int u = q.front();
          q.pop();
          for (auto &e : Constraints[u])
            if (e.AdjC != -1) {
              int v = e.AdjC;
              if (comp[v] == -1) {
                comp[v] = compCnt;
                q.push(v);
              }
            }
        }
        ++compCnt;
      }

    // Odds by component
    std::vector<std::vector<int>> oddsInComp(compCnt);
    std::vector<char> isOdd(nC, 0);
    for (int c : NotEvenNodes)
      if (0 <= c && c < nC) {
        isOdd[c] = 1;
        oddsInComp[comp[c]].push_back(c);
      }

    std::unordered_set<int> toggle;
    toggle.reserve(2048);
    auto toggle_var = [&](int varId) {
      auto it = toggle.find(varId);
      if (it == toggle.end())
        toggle.insert(varId);
      else
        toggle.erase(it); // XOR
    };

    std::vector<char> fixedFlag(nC, 0); // mark constraints we manage to fix

    for (int cid = 0; cid < compCnt; ++cid) {
      auto &odds = oddsInComp[cid];
      if (odds.empty())
        continue;
      std::unordered_set<int> oddSet(odds.begin(), odds.end());

      // Phase 1: use leaves whenever possible (fix as many singles as we can)
      bool changed = true;
      while (changed) {
        changed = false;
        std::vector<int> pick;
        for (int c : oddSet)
          if (leafVar[c] != -1)
            pick.push_back(c);
        if (pick.empty())
          break;
        for (int c : pick) {
          toggle_var(leafVar[c]);
          oddSet.erase(c);
          fixedFlag[c] = 1;
          R.usedLeaves++;
          R.totalCost += 1;
          changed = true;
        }
      }

      // Phase 2: pair nearest odds by shortest paths (greedy)
      while (oddSet.size() >= 2) {
        int u = *oddSet.begin();
        BFSData bfs = bfs_from(Constraints, u);
        int bestV = -1, bestD = INT_MAX;
        for (int v : oddSet)
          if (v != u && bfs.dist[v] < bestD) {
            bestD = bfs.dist[v];
            bestV = v;
          }
        if (bestV == -1 || bestD == INT_MAX)
          break; // no reachable partner via internal edges

        std::vector<int> pathVars = reconstruct_path_vars(bfs, u, bestV);
        for (int vid : pathVars)
          toggle_var(vid);
        fixedFlag[u] = fixedFlag[bestV] = 1;
        oddSet.erase(u);
        oddSet.erase(bestV);
        R.usedPairs++;
        R.totalCost += (int)pathVars.size();
      }

      // Whatever remains in oddSet is unresolved in this component
    }

    // Build output lists from original NotEvenNodes
    for (int c : NotEvenNodes) {
      if (fixedFlag[c])
        R.fixedConstraints.push_back(c);
      else
        R.unresolvedConstraints.push_back(c);
    }
    R.varsToFlip.assign(toggle.begin(), toggle.end());
    std::sort(R.varsToFlip.begin(), R.varsToFlip.end());
    return R;
  }

  void SolveConstraints(const std::vector<std::vector<AdjConstr>> &Constraints,
                        const std::vector<int> &ToPairNodes,
                        std::vector<int> &varsToFlip,
                        std::vector<int> &fixedConstraints,
                        std::vector<int> &unresolvedConstraints,
                        const int MaxPairNum = 24) {
    if (ToPairNodes.size() <= 24) {
      PCJoinResult res =
          prize_collecting_tjoin_maxcover(Constraints, ToPairNodes);
      if (!res.ok) {
        std::cout << "Failure" << std::endl;
        return;
      }
      varsToFlip = res.varsToFlip;
      fixedConstraints = res.fixedConstraints;
      unresolvedConstraints = res.unresolvedConstraints;
    } else {
      HeuResult res = heuristic_maxcover(Constraints, ToPairNodes);
      varsToFlip = res.varsToFlip;
      fixedConstraints = res.fixedConstraints;
      unresolvedConstraints = res.unresolvedConstraints;
    }
    std::cout << "Solved constraints Num: " << fixedConstraints.size()
              << std::endl;
    std::cout << "Unsolved constraints Num:" << unresolvedConstraints.size()
              << std::endl;
    std::cout << "Change variables Num:" << varsToFlip.size() << std::endl;
  }

public:
  void
  SolveParity(std::vector<std::vector<std::vector<int>>> &_PatchSideVariables) {

    InitParityConstr(_PatchSideVariables);

    if (NotParityNodes.size() > 0) {
      // if (NotEvenNodes.size())
      std::vector<int> varsToFlip, fixedConstraints, unresolvedConstraints;

      SolveConstraints(ParityConstraintsAdj, NotParityNodes, varsToFlip,
                       fixedConstraints, unresolvedConstraints);

      // flip variables
      for (size_t i = 0; i < varsToFlip.size(); i++) {
        int IndexVar = varsToFlip[i];
        if (VarValue[IndexVar] == 1) {
          VarValue[IndexVar]++;
          continue;
        }
        if (VarValue[IndexVar] > VarTargetValue[IndexVar]) {
          VarValue[IndexVar]--;
          continue;
        }
        VarValue[IndexVar]++;
      }
    }
  }

  ConstraintGraph(std::vector<int> &_VarValue,
                  std::vector<int> &_VarTargetValue)
      : VarTargetValue(_VarTargetValue), VarValue(_VarValue) {}
};

#endif
