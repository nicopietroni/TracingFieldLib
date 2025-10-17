#ifndef LABEL_GRAPH_FUNCTIONS
#define LABEL_GRAPH_FUNCTIONS

#include "../cross_field.h"
#include "../frame_field.h"
#include "../point3.h"
#include "best_rotation.h"
#include "connectivity.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "label_graph.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "patch_basic.h"
#include "polyhedra_mesh_functions.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <random>
#include <stdio.h>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class LabelGraphFuctions {
  

public:
  static bool IsDeadEnd(LabelNeighOrientedGraph<ScalarType> &Facegraph,
                        int IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < Facegraph.NumNodes());
    if (!Facegraph.IsEnabled(IndexN))
      return true;

    for (size_t e = 0; e < 3; e++) {

      if (Facegraph.IsExitLabel(IndexN, e))
        continue;

      // std::vector<int> IndexNeigh;
      // std::vector<ScalarType> WeightNeigh;
      // Facegraph.ActiveNeigh(IndexN, e, IndexNeigh, WeightNeigh);

      // if (IndexNeigh.size() == 0)
      //   return true;
      if (!Facegraph.HasActiveNeigh(IndexN,e))return true;
    }
    return false;
  }
  static void
  DeactivateTangentNodes(LabelNeighOrientedGraph<ScalarType> &Facegraph) {
    std::vector<int> IdxNode;
    Facegraph.GetNodeOfType(TangentDomain, IdxNode);
    for (size_t i = 0; i < IdxNode.size(); i++)
      Facegraph.DisableNode(IdxNode[i]);
    std::cout<<"Disabled "<<IdxNode.size()<<" out of "<<Facegraph.NumNodes()<<std::endl;

  }

  static void
  DeactivateInternalDeadEnds(LabelNeighOrientedGraph<ScalarType> &Facegraph) {
    int DeactivatedNum = 0;
    bool has_deactivated = false;
    int NumBefore = ActiveInternalNodes(Facegraph);
    do {
      has_deactivated = false;
      for (size_t i = 0; i < Facegraph.NumNodes(); i++) {

        int IndexN = i;

        if (!Facegraph.IsEnabled(IndexN))
          continue;

        if (!IsDeadEnd(Facegraph, IndexN))
          continue;

        Facegraph.DisableNode(IndexN);
        has_deactivated = true;
        DeactivatedNum++;
      }
    } while (has_deactivated);
    int NumAfter = ActiveInternalNodes(Facegraph);
    std::cout << "DISABLED "
              << 100 * (ScalarType)(NumBefore - NumAfter) /
                     (ScalarType)NumBefore
              << " perc of Internal nodes that were dead ends" << std::endl;
  }

  //this function is used to symmetrize the connections
  static void
  SymmetrizeM3(Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
               const std::vector<std::vector<int>> &TetraInternalFaces) {
    
    //structure to keep info on connections
    struct ConnInfo {
      ScalarType W;
      bool IsActive;

      ConnInfo(ScalarType _W, bool _IsActive) {
        W = _W;
        IsActive = _IsActive;
      }

      ConnInfo() {
        W = -1;
        IsActive = false;
      }
    };

    // for each pair of nodes saves the connections
    std::map<std::pair<int, int>, ConnInfo> NodeNodeInfo;

    //for each node
    for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
      int IndexN0 = i;

      //for each label
      for (size_t l = 0; l < 3; l++) {
        int IndexLabel = l;
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        std::vector<bool> ActiveLink;

        //get all neighbours info
        Facegraph.AllNeigh(IndexN0, IndexLabel, IndexNeigh, WeightNeigh,
                           ActiveLink);
        
        //then for each neighbours
        for (size_t j = 0; j < IndexNeigh.size(); j++) {
          int IndexN1 = IndexNeigh[j];
          //save the connection info
          ConnInfo CI(WeightNeigh[j], ActiveLink[j]);
          std::pair<int, int> VertPair(IndexN0, IndexN1);
          NodeNodeInfo[VertPair] = CI;
        }
      }
    }

    //THIS HAS TO BE COMPLETED
  }

  static void
  GetM3ConflictsConstr(const LabelNeighOrientedGraph<ScalarType> &Graph,
                       std::vector<std::vector<int>> &InternalConstraints,
                       std::vector<std::vector<int>> &ExitConstraints) {
    InternalConstraints.clear();
    ExitConstraints.clear();

    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!Graph.IsEnabled(IndexN))
        continue;

      for (size_t e = 0; e < 3; e++) {
        std::vector<int> IndexNeigh;
        Graph.ActiveNeigh(IndexN, e, IndexNeigh);
        // add the initial one
        IndexNeigh.push_back(IndexN);
        // sort to reduce the total amoung
        std::sort(IndexNeigh.begin(), IndexNeigh.end());

        if (Graph.IsExitLabel(IndexN, e))
          ExitConstraints.push_back(IndexNeigh);
        else
          InternalConstraints.push_back(IndexNeigh);
      }
    }
    std::cout << "Exit Constraint Before Reduction " << ExitConstraints.size()
              << std::endl;
    std::sort(ExitConstraints.begin(), ExitConstraints.end());
    ExitConstraints.resize(
        std::unique(ExitConstraints.begin(), ExitConstraints.end()) -
        ExitConstraints.begin());
    std::cout << "Exit Constraint After Reduction " << ExitConstraints.size()
              << std::endl;

    std::cout << "Internal Constraint Before Reduction "
              << InternalConstraints.size() << std::endl;
    std::sort(InternalConstraints.begin(), InternalConstraints.end());
    InternalConstraints.resize(
        std::unique(InternalConstraints.begin(), InternalConstraints.end()) -
        InternalConstraints.begin());
    std::cout << "Internal Constraint After Reduction "
              << InternalConstraints.size() << std::endl;
  }

  static void GetInternalNodesWithEmptyLabel(
      const LabelNeighOrientedGraph<ScalarType> &Graph,
      std::vector<int> &SelectedNodes) {
    
    SelectedNodes.clear();
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!Graph.IsEnabled(IndexN))
        continue;
      if (Graph.NodeType(IndexN) == TangentDomain)
        continue;
      bool has_empty = false;
      for (size_t e = 0; e < 3; e++) {
        // if border edge then does not matters
        if (Graph.NodeBorderE[IndexN][e])
          continue;

        // std::vector<int> IndexNeigh;
        // Graph.ActiveNeigh(IndexN, e, IndexNeigh);
        // has_empty |= (IndexNeigh.size() == 0);
        has_empty |=(!Graph.HasActiveNeigh(IndexN,e));
        if (has_empty)break;
      }
      if (has_empty)
        SelectedNodes.push_back(i);
    }
  }

  static int
  ActiveInternalNodes(const LabelNeighOrientedGraph<ScalarType> &Graph) {
    int NumNodes = 0;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (!Graph.IsEnabled(i))
        continue;
      if (Graph.NodeType(i) == TangentDomain)
        continue;
      NumNodes++;
    }
    return NumNodes;
  }

  static void DisableInternalNodesWithEmptyLabels(
      LabelNeighOrientedGraph<ScalarType> &Graph) {
    int NumBefore = ActiveInternalNodes(Graph);
    bool had_disabled = false;
    int NumDisabled = 0;
    int perf_steps=0;
    // int t0=clock();
    do {
      std::vector<int> SelectedNodes;
      GetInternalNodesWithEmptyLabel(Graph, SelectedNodes);
      had_disabled = (SelectedNodes.size() > 0);
      NumDisabled += SelectedNodes.size();
      perf_steps++;
      for (size_t i = 0; i < SelectedNodes.size(); i++)
        Graph.DisableNode(SelectedNodes[i]);
     
    } while (had_disabled);
    //  int t1=clock();
    // std::cout<<"DEDEED "<<t1-t0<<std::endl;
    std::cout << "DISABLED "
              << 100 * ((ScalarType)NumDisabled / (ScalarType)NumBefore)
              << " perc of Internal nodes with empty labels" << std::endl;
  }

  static void GetActiveFrom(const LabelNeighOrientedGraph<ScalarType> &Graph,
                            const std::vector<int> &SeedNodes,
                            std::set<int> &ActiveNodes) {
    ActiveNodes.clear();
    std::vector<int> GraphStack;

    for (size_t i = 0; i < SeedNodes.size(); i++)
      GraphStack.push_back(SeedNodes[i]);

    do {
      int CurrN = GraphStack.back();
      GraphStack.pop_back();
      if (!Graph.IsEnabled(CurrN))
        continue;

      ActiveNodes.insert(CurrN);
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      Graph.ActiveNeigh(CurrN, IndexNeigh, WeightNeigh);
      for (size_t i = 0; i < IndexNeigh.size(); i++) {
        int IndexNeighN = IndexNeigh[i];
        if (ActiveNodes.count(IndexNeighN) > 0)
          continue;
        GraphStack.push_back(IndexNeighN);
      }
    } while (!GraphStack.empty());
  }

  static void
  DisableNotUsefulFromSeed(LabelNeighOrientedGraph<ScalarType> &Graph,
                           const std::vector<int> &SeedNodes,
                           std::vector<int> &DisabledN) {
    std::set<int> ActiveNodes;
    DisabledN.clear();
    GetActiveFrom(Graph, SeedNodes, ActiveNodes);
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (ActiveNodes.count(i) > 0)
        continue;
      if (Graph.IsEnabled(i))
        continue;
      Graph.DisableNode(i);
      DisabledN.push_back(i);
    }
  }

  static void RemoveNodeLowQuality(LabelNeighOrientedGraph<ScalarType> &Graph,
                                   ScalarType MaxW = 0.66,
                                   int MinRedoundancy = -1) {

    // first get connectivity map, for each node the
    // node pointing to it and the label
    std::vector<std::vector<std::pair<int, int>>> ReferredFrom(
        Graph.NumNodes());

    if (MinRedoundancy != -1) {
      for (size_t i = 0; i < Graph.NumNodes(); i++) {
        int IndexN = i;

        if (!Graph.IsEnabled(IndexN))
          continue;
        if (Graph.NodeType(IndexN) == TangentDomain)
          continue;

        for (size_t l = 0; l < 3; l++) {
          std::vector<int> IndexNeigh;
          int IndexLabel = l;
          Graph.ActiveNeigh(IndexN, IndexLabel, IndexNeigh);
          for (size_t j = 0; j < IndexNeigh.size(); j++) {
            int NeighNode = IndexNeigh[j];
            assert(NeighNode >= 0);
            assert(NeighNode < ReferredFrom.size());
            ReferredFrom[NeighNode].push_back(
                std::pair<int, int>(IndexN, IndexLabel));
          }
        }
      }
    }

    // then sort nodes by considering the quality
    std::vector<std::pair<ScalarType, int>> VertWIdx;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;

      if (!Graph.IsEnabled(IndexN))
        continue;
      if (Graph.NodeType(IndexN) == TangentDomain)
        continue;
      if (Graph.HasExitLabel(IndexN))
        continue;
      ScalarType Q = Graph.NodeW[IndexN];
      assert(Q >= 0);

      VertWIdx.push_back(std::pair<ScalarType, int>(Q, IndexN));
    }

    // sort and reverse
    std::sort(VertWIdx.begin(), VertWIdx.end());
    std::reverse(VertWIdx.begin(), VertWIdx.end());

    int NumDisabled = 0;
    bool has_disabled = false;
    do {
      has_disabled = false;

      for (size_t i = 0; i < VertWIdx.size(); i++) {

        int NodeIdx = VertWIdx[i].second;
        ScalarType CurrW = VertWIdx[i].first;
        if (CurrW <= MaxW)
          continue;

        if (!Graph.IsEnabled(NodeIdx))
          continue;
        if (Graph.NodeType(NodeIdx) == TangentDomain)
          continue;
        if (Graph.HasExitLabel(NodeIdx))
          continue;

        // get the ones pointing to it
        if (MinRedoundancy != -1) {
          int MinN = std::numeric_limits<int>::max();
          for (size_t j = 0; j < ReferredFrom[NodeIdx].size(); j++) {
            int IndexNFrom = ReferredFrom[NodeIdx][j].first;
            int IndexLFrom = ReferredFrom[NodeIdx][j].second;

            if (!Graph.IsEnabled(IndexNFrom))
              continue;
            if (Graph.NodeType(IndexNFrom) == TangentDomain)
              continue;

            int NumNeigh=Graph.NumActiveNeigh(IndexNFrom,IndexLFrom);
            MinN = std::min(MinN,NumNeigh-1);
          }
          if (MinN < MinRedoundancy)
            continue;
        }
        Graph.DisableNode(NodeIdx);
        NumDisabled++;
        has_disabled = true;
      }
      std::cout << "** Disabled " << NumDisabled << " Nodes **" << std::endl;
    } while (has_disabled);
  }

  //this one resize the number of neighbour 
  // for each node for each label to a maximum number
  static void ResizeNeigh(LabelNeighOrientedGraph<ScalarType> &Graph,
                          int NumNeigh = 2) {
    //for each node
    for (size_t i = 0; i < Graph.NumNodes(); i++)
    //for each label
      for (size_t l = 0; l < 3; l++) {
        std::vector<std::pair<ScalarType, int>> WeightNeigh;
        std::vector<int> IndexN;
        std::vector<ScalarType> WeightN;
        //if the node has exit label try to keep
        //maximum connectivity
        if (Graph.HasExitLabel(i))
          continue;
        //if the node is tangent doesn't care
        if (Graph.NodeType(i) == TangentDomain)
          continue;
        //if the node is not enabled doesn't care
        if (!Graph.IsEnabled(i))
          continue;

        //get the active neighbour for that 
        //node for that label
        Graph.ActiveNeigh(i, l, IndexN, WeightN);

        //if the number of neighbours is less than
        //the thresold than does nothing
        if (IndexN.size() <= NumNeigh)
          continue;

        //then for each neighbours push weight and neighbour index
        for (size_t j = 0; j < IndexN.size(); j++) {
          WeightNeigh.push_back(
              std::pair<ScalarType, int>(WeightN[j], IndexN[j]));
        }

        //sort them
        std::sort(WeightNeigh.begin(), WeightNeigh.end());

        //resize to the maximum
        WeightNeigh.resize(NumNeigh);

        //clear the set of neighbours for that node
        // and that specific label (edge of face)
        Graph.ClearNeighbor(i, l);

        //then re-add the neighbours
        for (size_t j = 0; j < WeightNeigh.size(); j++)
          Graph.AddNeighbour(i, WeightNeigh[j].second, l, WeightNeigh[j].first);
      }
  }

  static void IncreaseWeight(LabelNeighOrientedGraph<ScalarType> &Graph,
                             ScalarType &QualityTHR, ScalarType &Multiplier) {
    int NumDisable = 0;
    int ToTal = 0;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;
      ScalarType WeightN = Graph.NodeW[IndexN];
      if (!Graph.IsEnabled(IndexN))
        continue;
      if (Graph.HasExitLabel(IndexN))
        continue;
      ToTal++;
      if (WeightN > QualityTHR) {
        Graph.NodeW[IndexN] *= Multiplier;
      }
    }
    // std::cout<<"** Disabled "<<NumDisable<<" out of"<<ToTal<<std::endl;
  }

  static void ReactivateNodesOnConcaveFeatures(LabelNeighOrientedGraph<ScalarType> &Graph)
  {
    for (size_t i=0;i<Graph.NodesOnConcaveE.size();i++)
      Graph.EnableNode(Graph.NodesOnConcaveE[i]);

     for (size_t i=0;i<Graph.NodesOnConcaveEByV.size();i++)
      Graph.EnableNode(Graph.NodesOnConcaveEByV[i]);
  }


  static void GetM3ExitDrawPosition(const LabelNeighOrientedGraph<ScalarType> &Graph,
                                    const std::vector<int> &IndexN,
                                    const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                    const std::vector<std::vector<int> > &InternalFaces,
                                    std::vector<Geo::Point3<ScalarType> > &MidPointExit) {
    MidPointExit.clear();
    for (size_t i = 0; i < IndexN.size(); i++) {
     int CurrN=IndexN[i];
     for (size_t l = 0; l < Graph.N[CurrN].size(); l++) {
        if (!Graph.IsExitLabel(CurrN, l)) continue;
        int IndexF= Field::FrameFieldM3Graph::ElemOfNode(CurrN);
        int IndexV0= InternalFaces[IndexF][l];
        int IndexV1= InternalFaces[IndexF][(l+1)%InternalFaces[IndexF].size()];
        Geo::Point3<ScalarType> Pos0=VertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1=VertPos[IndexV1];
        Geo::Point3<ScalarType> MidP=Pos0*0.5+Pos1*0.5;
        MidPointExit.push_back(MidP);
      }
    }
  }
};

} // namespace Geo
#endif
