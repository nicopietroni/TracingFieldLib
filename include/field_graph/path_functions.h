#ifndef PATH_FUNCTIONS
#define PATH_FUNCTIONS

#include "../math_utils.h"
#include "connectivity.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "graph_tracing.h"
#include "path_sampling.h"
#include "point3.h"
#include "t_junction.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <list>
#include <set>
#include <sys/_types/_size_t.h>
#include <unordered_set>
#include <vector>

namespace Geo {

template <class ScalarType> class PathFunctions {
public:
  static bool IsProperPath(const std::vector<int> &PathN) {
    for (size_t i = 0; i < PathN.size() - 1; i++)
      if (PathN[i] == PathN[(i + 1)])
        return false;
    return true;
  }

  static bool IsLoop(const std::vector<int> &PathN) {
    return (PathN[0] == PathN.back());
  }

  static bool IsConnected(const TracingGraph<ScalarType> &TGraph,
                          const std::vector<int> &PathN) {

    assert(PathN.size() >= 2);
    int startIDX = 0;
    int endIDX = PathN.size();
    if (!PathFunctions<ScalarType>::IsLoop(PathN)) {
      startIDX++;
      endIDX--;
    }
    // the first and the last can be
    // non enabled but still used to trace
    for (size_t i = startIDX; i < endIDX; i++) {
      if (!TGraph.IsEnabled(PathN[i])) {
        std::cout << "Non Enabled Node" << std::endl;
        return false;
      }
    }

    for (size_t i = 0; i < PathN.size() - 1; i++) {
      int IndexN0 = PathN[i];
      int IndexN1 = PathN[(i + 1)];

      if (TGraph.WhichConnection(IndexN0, IndexN1) == -1) {
        std::cout << "Non Symmetric Arc" << std::endl;
        return false;
      }
      if (!TGraph.ActiveArcNodes(IndexN0, IndexN1)) {
        std::cout << "Non Active Arc" << std::endl;
        return false;
      }
    }
    return true;
  }

  static bool ConnectDirect(TracingGraph<ScalarType> &TGraph,
                            const int &IndexN0, const int &IndexN1,
                            std::vector<int> &DirectPath) {
    typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    SParam.avoid_zero_path = true;
    SParam.StopAtSel[0] = true;

    TGraph.DeSelectAll(0);
    TGraph.Select(IndexN1, 0);
    SParam.onlyDirect = true;
    SParam.onlyVisible = false;
    // TO CHECK
    SParam.uniformW = true;
    SParam.MaxSteps = TGraph.diffuse_step + 1;
    std::vector<int> SourceNodes(1, IndexN0);
    ScalarType PathL;
    bool has_connected = Geo::GraphQuery<ScalarType>::ShortestPaths(
        TGraph, SourceNodes, DirectPath, PathL, SParam);
    //   for (size_t i=0;i<DirectPath.size();i++)
    //   {
    //     if (!TGraph.IsEnabled(DirectPath[i]))
    //     {
    //       std::cout<<"Error Enabled Path 1"<<std::endl;
    //       exit(0);
    //     }

    //     if (TGraph.NodeType(DirectPath[i])==Geo::TangentDomain)
    //     {
    //       std::cout<<"Error Tangent Path 1"<<std::endl;
    //       exit(0);
    //     }
    //  }
    return has_connected;
  }

  // simplify the connections to have one edge per triangle/tetra
  // only considering the direct one
  static bool SimplifyStep(TracingGraph<ScalarType> &TGraph,
                           std::vector<int> &PathN) {
    if (PathN.size() < 3)
      return false;

    std::vector<int> NewPath;
    bool isLoop = IsLoop(PathN);
    NewPath = PathN;

    // remove the last one which is equal to the first one
    // this is to make the right cycle then eventually put back at the end
    if (isLoop) {
      assert(NewPath[0] == NewPath.back());
      NewPath.pop_back();
    }

    // first and last threated separately if needed (only if loop)
    int StartIdx = 0;
    int EndIdx = NewPath.size();
    if (!isLoop) {
      StartIdx = 1;
      EndIdx = NewPath.size() - 1;
    }
    for (size_t i = StartIdx; i < EndIdx; i++) {
      int SizeP = NewPath.size();
      int IndexN0 = PathN[(i + SizeP - 1) % SizeP];
      int IndexN1 = PathN[i];
      int IndexN2 = PathN[(i + 1) % SizeP];

      // see if there is a direct connection
      // then cancel the one in the middle
      if (TGraph.IsDirectPairE(IndexN0, IndexN2)) {
        NewPath.erase(NewPath.begin() + i);
        // if loops reinsert the first one to the end
        if (isLoop)
          NewPath.push_back(NewPath[0]);
        // copy back to return
        PathN = NewPath;
        return true;
      }
    }
    return false;
  }

  // simplify the connections to have one edge per triangle/tetra
  // only considering the direct one
  static void SimplifyConnection(TracingGraph<ScalarType> &TGraph,
                                 std::vector<int> &PathN) {
    assert(IsProperPath(PathN));

    if (PathN.size() < 3)
      return;

    // std::cout<<"Start Simpplify"<<std::endl;
    while (SimplifyStep(TGraph, PathN)) {
    };
    // std::cout<<"Done"<<std::endl;
  }

  static bool MakeDirectConnection(TracingGraph<ScalarType> &TGraph,
                                   std::vector<int> &PathN,
                                   bool UpdateArc = true) {
    //TGraph.CurrDistIndex() = 1;
    assert(IsProperPath(PathN));
    std::vector<int> NewPath;
    bool IsLoop = PathFunctions<ScalarType>::IsLoop(PathN);
    for (size_t i = 0; i < PathN.size() - 1; i++) {
      int IndexN0 = PathN[i];
      int IndexN1 = PathN[(i + 1)];

      // add the first one
      NewPath.push_back(IndexN0);

      if (TGraph.IsDirectPairE(IndexN0, IndexN1))
        continue;

      std::vector<int> DirectPath;
      // first and last can be disabled if it is not a loop
      if (IsLoop) {
        assert(TGraph.IsEnabled(IndexN0));
        assert(TGraph.IsEnabled(IndexN1));
      } else {
        if ((i != 0) && (i != (PathN.size() - 1))) {
          assert(TGraph.IsEnabled(IndexN0));
          assert(TGraph.IsEnabled(IndexN1));
        }
      }

      assert(TGraph.WhichConnection(IndexN0, IndexN1) != -1);
      assert(TGraph.ActiveArcNodes(IndexN0, IndexN1));
      // assert(TGraph.WhichConnection(IndexN1,IndexN0)!=-1);
      // assert(TGraph.ActiveArcNodes(IndexN1, IndexN0));

      bool connected = ConnectDirect(TGraph, IndexN0, IndexN1, DirectPath);
      if (!connected) {
        if (UpdateArc) {
          TGraph.DectivateArc(IndexN0, IndexN1);
          int OppN1 = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN1);
          int OppN0 = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN0);
          TGraph.DectivateArc(OppN1, OppN0);
        }
        return false;
      }
      // the new path can be only 3 + diffusion steps (should be 2 but one
      // because of singularities)
      // if (DirectPath.size() > (TGraph.diffuse_step + 3))
      //   return false;
      // then from second to one before the last
      for (size_t j = 1; j < DirectPath.size() - 1; j++)
        NewPath.push_back(DirectPath[j]);

      assert(DirectPath[0] == IndexN0);
      assert(DirectPath.back() == IndexN1);
    }
    // then add the last
    NewPath.push_back(PathN.back());
    PathN = NewPath;
    return true;
  }

  static int OppositeNode(const Geo::TracingGraph<ScalarType> &TGraph,
                          const int &indexN) {
    // if the opposite is selected then is on the boundary
    Geo::GraphType GType = TGraph.TypeofGraph();
    int OppN = -1;
    if (GType == Geo::GT_CrossField)
      OppN = Field::CrossFieldGraph::OppositeNodeCrossField(indexN);
    if (GType == Geo::GT_FrameField)
      OppN = Field::FrameFieldGraph::OppositeNodeFrameField(indexN);
    assert(OppN != -1);
    return OppN;
  }

  static void OppositeNodes(const TracingGraph<ScalarType> &TGraph,
                            const std::vector<int> &PathN,
                            std::vector<int> &OppPathN) {
    OppPathN.clear();
    for (size_t i = 0; i < PathN.size(); i++) {
      int IndexN = PathN[i];
      assert(IndexN >= 0);
      assert(IndexN < TGraph.NumNodes());

      int OppN = OppositeNode(TGraph, IndexN);

      assert(OppN >= 0);
      assert(OppN < (int)TGraph.NumNodes());

      OppPathN.push_back(OppN);
    }
  }

  static void OppositeNodesCross(const std::vector<int> &PathN,
                                 std::vector<int> &OppPathN) {
    OppPathN.clear();
    for (size_t i = 0; i < PathN.size(); i++) {
      int IndexN = PathN[i];
      int OppN = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
      OppPathN.push_back(OppN);
    }
  }

  static bool SelfIntersectTangentially(const TracingGraph<ScalarType> &TGraph,
                                        const std::vector<int> &PathN) {
    if (PathN.size() < 2)
      return false;
    std::unordered_set<int> TestNSet;

    std::vector<int> OppPathN;
    OppositeNodes(TGraph, PathN, OppPathN);
    assert(OppPathN.size() == PathN.size());

    int LastTestN = PathN.size();
    if (IsLoop(PathN))
      LastTestN--;
    for (int i = 0; i < LastTestN; i++) {
      int IndexN = PathN[i];
      int OppN = OppPathN[i];
      assert(IndexN >= 0);
      assert(IndexN < (int)TGraph.NumNodes());
      assert(OppN >= 0);
      assert(OppN < (int)TGraph.NumNodes());

      if (TestNSet.count(IndexN) > 0)
        return true;
      if (TestNSet.count(OppN) > 0)
        return true;
      TestNSet.insert(IndexN);
      TestNSet.insert(OppN);
    }
    return false;
  }

  static bool IntersectTangentially(const TracingGraph<ScalarType> &TGraph,
                                    const std::vector<int> &PathN0,
                                    const std::vector<int> &PathN1) {

    std::vector<int> OppPathN0;
    OppositeNodes(TGraph, PathN0, OppPathN0);
    assert(OppPathN0.size() == PathN0.size());

    std::unordered_set<int> TestNSet(PathN0.begin(), PathN0.end());
    TestNSet.insert(OppPathN0.begin(), OppPathN0.end());

    std::vector<int> OppPathN1;
    OppositeNodes(TGraph, PathN1, OppPathN1);
    assert(OppPathN1.size() == PathN1.size());

    int LastTestN = PathN1.size();
    if (IsLoop(PathN1))
      LastTestN--;
    for (int i = 0; i < LastTestN; i++) {
      int IndexN = PathN1[i];
      int OppN = OppPathN1[i];
      assert(IndexN >= 0);
      assert(IndexN < TGraph.NumNodes());
      assert(OppN >= 0);
      assert(OppN < TGraph.NumNodes());

      if (TestNSet.count(IndexN) > 0)
        return true;
      if (TestNSet.count(OppN) > 0)
        return true;
    }
    return false;
  }

  static bool HasSelfEdgeTwice(const std::vector<int> &PathN) {
    std::unordered_set<std::pair<int, int>, Geo::hashPairIntFunction> TestNSet;
    for (int i = 0; i < PathN.size() - 1; i++) {
      int IndexN0 = PathN[i];
      int IndexN1 = PathN[i + 1];
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      if (TestNSet.count(key) > 0)
        return true;
      TestNSet.insert(key);
    }
    return false;
  }

  static void ReversePath(const TracingGraph<ScalarType> &TGraph,
                          std::vector<int> &PathN) {
    std::vector<int> OppPathN;
    OppositeNodes(TGraph, PathN, OppPathN);
    std::reverse(OppPathN.begin(), OppPathN.end());
    PathN = OppPathN;
  }

  static bool BestConnectionFrom(const TracingGraph<ScalarType> &TGraph,
                                 const int &Node0, const int &VertIdx1,
                                 int &NodeIdx1, ScalarType &bestW) {
    NodeIdx1 = -1;

    bool found = false;

    bestW = std::numeric_limits<ScalarType>::max();

    std::vector<int> Neigh;
    std::vector<ScalarType> NeighW;
    TGraph.GetNeigh(Node0, true, true, Neigh, NeighW);

    for (size_t j = 0; j < Neigh.size(); j++) {
      int NeighV = -1;

      if (TGraph.TypeofGraph() == GT_CrossField)
        NeighV = Field::CrossFieldGraph::ElemOfNode(Neigh[j]);
      else
        NeighV = Field::CrossFieldGraph::ElemOfNode(Neigh[j]);

      if (NeighV != VertIdx1)
        continue;

      found = true;

      if (NeighW[j] > bestW)
        continue;

      bestW = NeighW[j];
      NodeIdx1 = Neigh[j];
    }

    return found;
  }

  static bool BestConnectionNodePair(const TracingGraph<ScalarType> &TGraph,
                                     const int &VertIdx0, const int &VertIdx1,
                                     int &NodeIdx0, int &NodeIdx1) {
    NodeIdx0 = -1;
    NodeIdx1 = -1;

    assert(VertIdx0 != VertIdx1);
    std::vector<int> IndexN0;

    // if (TGraph.TypeofGraph()==GT_CrossField)
    Field::CrossFieldGraph::NodesOfElem(VertIdx0, IndexN0);
    //        else
    //            Field::FrameFieldGraph::NodesOfElem(VertIdx0,IndexN0);

    ScalarType bestW = std::numeric_limits<ScalarType>::max();
    bool found = false;
    for (size_t i = 0; i < IndexN0.size(); i++) {
      int TestNodeIdx1;
      ScalarType TestbestW;
      bool found_this = BestConnectionFrom(TGraph, IndexN0[i], VertIdx1,
                                           TestNodeIdx1, TestbestW);
      if (!found_this)
        continue;

      if (TestbestW > bestW)
        continue;

      found = true;
      bestW = TestbestW;
      NodeIdx0 = IndexN0[i];
      NodeIdx1 = TestNodeIdx1;
    }
    return found;
  }

  // should work generic but to be tested
  static bool ReconstructPathFromElem(const TracingGraph<ScalarType> &Tgraph,
                                      const std::vector<int> &Elems,
                                      std::vector<int> &Nodes) {
    Nodes.clear();
    assert(Elems.size() >= 2);

    int E0 = Elems[0];
    int E1 = Elems[1];
    int NodeIdx0, NodeIdx1;
    bool found = BestConnectionNodePair(Tgraph, E0, E1, NodeIdx0, NodeIdx1);
    if (!found) {
      return false;
    }

    Nodes.push_back(NodeIdx0);
    for (size_t i = 1; i < Elems.size(); i++) {
      int NodeIdx1;
      ScalarType bestW;
      bool found =
          BestConnectionFrom(Tgraph, Nodes.back(), Elems[i], NodeIdx1, bestW);
      if (!found)
        return true;
      Nodes.push_back(NodeIdx1);
    }
    return true;
  }

  // should work generic but to be tested
  static bool ReconstructPathsFromElem(const TracingGraph<ScalarType> &Tgraph,
                                       const std::vector<int> &Elems,
                                       std::vector<std::vector<int>> &Nodes) {
    if (Elems.size() < 2)
      return false;

    Nodes.clear();

    std::vector<std::pair<int, int>> NodesE;
    for (size_t i = 0; i < Elems.size() - 1; i++) {
      // get the first direction Node
      int E0 = Elems[i];
      int E1 = Elems[i + 1];

      int NodeIdx0, NodeIdx1;
      bool found = BestConnectionNodePair(Tgraph, E0, E1, NodeIdx0, NodeIdx1);
      if (!found) {
        return false;
      }
      NodesE.push_back(std::pair<int, int>(NodeIdx0, NodeIdx1));
    }

    Nodes.resize(1);
    for (size_t i = 0; i < NodesE.size() - 1; i++) {
      Nodes.back().push_back(NodesE[i].first);
      int N0 = NodesE[i].second;
      int N1 = NodesE[(i + 1)].first;
      if (N0 != N1) {
        Nodes.back().push_back(N0);
        Nodes.resize(Nodes.size() + 1);
      }
    }
    Nodes.back().push_back(NodesE.back().first);
    Nodes.back().push_back(NodesE.back().second);
    return true;
  }

  static void CrossLinkBetweenVert(const Geo::TracingGraph<ScalarType> &Graph,
                                   const int &V0, const int &V1,
                                   std::vector<int> &IndexN,
                                   std::vector<int> &IndexL,
                                   bool onlydirect = true) {
    IndexN.clear();
    IndexL.clear();
    std::vector<int> IndexN0, IndexN1;
    Field::CrossFieldGraph::NodesOfElem(V0, IndexN0);
    Field::CrossFieldGraph::NodesOfElem(V1, IndexN1);
    std::set<int> IndexN1Set(IndexN1.begin(), IndexN1.end());
    for (size_t i = 0; i < IndexN0.size(); i++) {
      int Node0 = IndexN0[i];
      for (size_t j = 0; j < Graph.Neigh(IndexN0[i]).size(); j++) {
        int Node1 = Graph.Neigh(IndexN0[i])[j].N;
        bool isDirect = Graph.Neigh(IndexN0[i])[j].Direct;
        if ((!isDirect) && (onlydirect))
          continue;
        if (IndexN1Set.count(Node1) == 0)
          continue;
        IndexN.push_back(Node0);
        IndexL.push_back(j);
      }
    }
  }

  static void CrossEnableArcsOfPath(const std::vector<int> &PathN,
                                    Geo::TracingGraph<ScalarType> &Graph,
                                    bool enable) {
    int lastIDX = PathN.size();
    if (IsLoop(PathN))
      lastIDX--;
    for (int i = 0; i < lastIDX - 1; i++) {
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(PathN[i]);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(PathN[i + 1]);
      std::vector<int> IndexN0, IndexL0, IndexN1, IndexL1;
      CrossLinkBetweenVert(Graph, IndexV0, IndexV1, IndexN0, IndexL0);
      CrossLinkBetweenVert(Graph, IndexV1, IndexV0, IndexN1, IndexL1);

      if (enable) {
        for (size_t j = 0; j < IndexN0.size(); j++)
          Graph.EnableArc(IndexN0[j], IndexL0[j]);
        for (size_t j = 0; j < IndexN0.size(); j++)
          Graph.EnableArc(IndexN1[j], IndexL1[j]);
      } else {
        for (size_t j = 0; j < IndexN0.size(); j++)
          Graph.DisableArc(IndexN0[j], IndexL0[j]);
        for (size_t j = 0; j < IndexN0.size(); j++)
          Graph.DisableArc(IndexN1[j], IndexL1[j]);
      }
    }
  }

  static void FindEdgeM(const std::vector<std::vector<int>> &Paths,
                        std::vector<int> &EdgePathIdx,
                        std::vector<std::vector<int>> &Edges) {
    Edges.clear();
    EdgePathIdx.clear();
    for (size_t i = 0; i < Paths.size(); i++)
      for (size_t j = 0; j < Paths[i].size() - 1; j++) {
        int V0 = Paths[i][j];
        int V1 = Paths[i][j + 1];
        assert(V0 != V1);
        Edges.resize(Edges.size() + 1);
        Edges.back().push_back(V0);
        Edges.back().push_back(V1);
        EdgePathIdx.push_back(i);
      }
  }

  static void FindTJunctions(const std::vector<std::vector<int>> &Paths,
                             std::vector<TJunction> &TJunctions) {
    std::vector<int> EdgePathIdx;
    std::vector<std::vector<int>> Edges;

    FindEdgeM(Paths, EdgePathIdx, Edges);

    std::vector<std::vector<int>> Kernel;
    FindVertexStar(Edges, Kernel);

    std::map<std::pair<int, int>, int> VertPairToPath;
    for (size_t i = 0; i < Edges.size(); i++) {
      assert(Edges[i].size() == 2);
      int V0 = Edges[i][0];
      int V1 = Edges[i][1];
      int IndexP = EdgePathIdx[i];
      std::pair<int, int> Key(std::min(V0, V1), std::max(V0, V1));
      assert(VertPairToPath.count(Key) == 0);
      VertPairToPath[Key] = IndexP;
    }

    TJunctions.clear();
    for (size_t i = 0; i < Kernel.size(); i++) {
      int V = i;
      // different than 3 cannot be a t-junction
      if (Kernel[i].size() != 3)
        continue;
      int Vk0 = Kernel[i][0];
      int Vk1 = Kernel[i][1];
      int Vk2 = Kernel[i][2];

      assert(Vk0 >= 0);
      assert(Vk0 < Kernel.size());
      assert(Vk1 >= 0);
      assert(Vk1 < Kernel.size());
      assert(Vk2 >= 0);
      assert(Vk2 < Kernel.size());

      std::pair<int, int> E0(std::min(V, Vk0), std::max(V, Vk0));
      std::pair<int, int> E1(std::min(V, Vk1), std::max(V, Vk1));
      std::pair<int, int> E2(std::min(V, Vk2), std::max(V, Vk2));
      assert(VertPairToPath.count(E0) > 0);
      assert(VertPairToPath.count(E1) > 0);
      assert(VertPairToPath.count(E2) > 0);
      int PIdx0 = VertPairToPath[E0];
      int PIdx1 = VertPairToPath[E1];
      int PIdx2 = VertPairToPath[E2];

      // three different path, no t-junction
      if ((PIdx0 != PIdx1) && (PIdx0 != PIdx2) && (PIdx1 != PIdx2))
        continue;
      // three belongs to the same path, MUST test geometrically?
      if ((PIdx0 == PIdx1) && (PIdx0 == PIdx2))
        continue;

      int V0 = -1;
      int V1 = -1;
      int VT = -1;
      // three belongs to the same path, then test geometrically
      if (PIdx0 == PIdx1) {
        assert(PIdx0 != PIdx2);
        V0 = Vk0;
        V1 = Vk1;
        VT = Vk2;
      }
      if (PIdx1 == PIdx2) {
        assert(PIdx0 != PIdx1);
        V0 = Vk1;
        V1 = Vk2;
        VT = Vk0;
      }
      if ((PIdx0 == PIdx2) && (PIdx0 != PIdx1)) {
        assert(PIdx0 != PIdx1);
        V0 = Vk0;
        V1 = Vk2;
        VT = Vk1;
      }
      assert(V0 != -1);
      assert(V1 != -1);
      assert(VT != -1);
      TJunctions.push_back(TJunction(V, V0, V1, VT));
    }
  }

  static void FindCorners(const std::vector<Geo::Point3<ScalarType>> &Pos,
                          const std::vector<std::vector<int>> &PathsVert,
                          std::vector<int> &Corners) {
    std::vector<int> VertN(Pos.size(), 0);
    for (size_t i = 0; i < PathsVert.size(); i++)
      for (size_t j = 0; j < PathsVert[i].size() - 1; j++) {

        int V0 = PathsVert[i][j];
        int V1 = PathsVert[i][(j + 1)];

        VertN[V0]++;
        VertN[V1]++;
      }

    for (size_t i = 0; i < VertN.size(); i++) {
      if (VertN[i] <= 2)
        continue;
      Corners.push_back(i);
      // std::cout<<"size e "<<VertN[i]<<std::endl;
    }
  }

  static bool JoinIfPossible(std::vector<int> &PathN0,
                             std::vector<int> &PathN1) {
    if (PathN0.size() == 0)
      return false;
    if (PathN1.size() == 0)
      return false;

    assert(PathN0.size() >= 2);
    assert(PathN1.size() >= 2);

    if (IsLoop(PathN0))
      return false;
    if (IsLoop(PathN1))
      return false;

    std::pair<int, int> Extr0;
    std::pair<int, int> ExtrDir0;
    Extr0.first = Field::CrossFieldGraph::ElemOfNode(PathN0[0]);
    ExtrDir0.first = Field::CrossFieldGraph::DirOfNode(PathN0[0]) % 2;
    Extr0.second = Field::CrossFieldGraph::ElemOfNode(PathN0.back());
    ExtrDir0.second = Field::CrossFieldGraph::DirOfNode(PathN0.back()) % 2;
    assert(Extr0.first != Extr0.second);

    std::pair<int, int> Extr1;
    std::pair<int, int> ExtrDir1;
    Extr1.first = Field::CrossFieldGraph::ElemOfNode(PathN1[0]);
    ExtrDir1.first = Field::CrossFieldGraph::DirOfNode(PathN1[0]) % 2;
    Extr1.second = Field::CrossFieldGraph::ElemOfNode(PathN1.back());
    ExtrDir1.second = Field::CrossFieldGraph::DirOfNode(PathN1.back()) % 2;
    assert(Extr1.first != Extr1.second);

    // this is a loop
    if ((Extr0.first == Extr1.second) && (Extr0.second == Extr1.first) &&
        (ExtrDir0.first == ExtrDir1.second) &&
        (ExtrDir0.second == ExtrDir1.first)) {
      // remove the last one on the first
      // that might be duplicated after merge
      PathN0.pop_back();
      PathN0.insert(PathN0.end(), PathN1.begin(), PathN1.end());
      PathN1.clear();
      return true;
    }

    // again a loop
    if ((Extr0.first == Extr1.first) && (Extr0.second == Extr1.second) &&
        (ExtrDir0.first == ExtrDir1.first) &&
        (ExtrDir0.second == ExtrDir1.second)) {
      std::vector<int> OppN1;
      OppositeNodesCross(PathN1, OppN1);
      std::reverse(OppN1.begin(), OppN1.end());
      PathN0.pop_back();
      PathN0.insert(PathN0.end(), OppN1.begin(), OppN1.end());
      OppN1.clear();
      return true;
    }

    // join contiguous paths Path1 after Path0
    if ((Extr0.second == Extr1.first) && (ExtrDir0.second == ExtrDir1.first)) {
      PathN0.pop_back();
      PathN0.insert(PathN0.end(), PathN1.begin(), PathN1.end());
      PathN1.clear();
      return true;
    }

    // join contiguous paths Path0 after Path1
    if ((Extr0.first == Extr1.second) && (ExtrDir0.first == ExtrDir1.second)) {
      PathN1.pop_back();
      PathN1.insert(PathN1.end(), PathN0.begin(), PathN0.end());
      PathN0.clear();
      return true;
    }

    // join contiguous paths Path0 after Path1
    if ((Extr0.first == Extr1.first) && (ExtrDir0.first == ExtrDir1.first)) {
      std::vector<int> OppN0;
      OppositeNodesCross(PathN0, OppN0);
      std::reverse(OppN0.begin(), OppN0.end());
      PathN0 = OppN0;
      PathN0.pop_back();
      PathN0.insert(PathN0.end(), PathN1.begin(), PathN1.end());
      PathN1.clear();
      return true;
    }

    // join contiguous paths Path0 after Path1
    if ((Extr0.second == Extr1.second) &&
        (ExtrDir0.second == ExtrDir1.second)) {
      std::vector<int> OppN1;
      OppositeNodesCross(PathN1, OppN1);
      std::reverse(OppN1.begin(), OppN1.end());
      PathN0.pop_back();
      PathN0.insert(PathN0.end(), PathN1.begin(), PathN1.end());
      PathN1.clear();
      return true;
    }
    return false;
  }

  static ScalarType
  CrossEuclideanLenght(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<int> &PathN) {
    ScalarType CurrL = 0;
    for (size_t i = 0; i < PathN.size() - 1; i++) {
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(PathN[i]);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(PathN[i + 1]);
      assert(IndexV0 >= 0);
      assert(IndexV0 < VertPos.size());
      assert(IndexV1 >= 0);
      assert(IndexV1 < VertPos.size());
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      CurrL += (P0 - P1).Norm();
    }
    return CurrL;
  }
};

} // namespace Geo
#endif
