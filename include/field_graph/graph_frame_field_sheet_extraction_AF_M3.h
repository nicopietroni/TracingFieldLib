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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_AF_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_AF_M3

#include "../connectivity.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

template <class ScalarType> class SheetExtractorAFM3 {
public:
  // different types of energy minimization
  enum ObjType { FaceOBJ, ConnOBJ, FaceConnOBJ };

  struct AF_ExtractParam {
    std::vector<int> Seeds;
    ObjType OType;

    AF_ExtractParam() { OType = FaceOBJ; }
  };

  std::vector<bool> IsSelectedN;
  std::vector<int> SelectedN;

private:
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  std::vector<std::vector<int>> &InternalFaces;
  std::vector<ScalarType> DistFunct;

  void InitSeedDistance(const AF_ExtractParam &ExtrP) {
    Geo::TMark TempMark;
    typename Geo::Geodesic<ScalarType>::GeodesicParam GeoP;
    GeoP.Seeds = ExtrP.Seeds;
    Geo::LabelGraphFuctions<ScalarType>::Dijkstra(M3Graph, TempMark, GeoP);
    DistFunct = GeoP.Dist;
  }

  std::set<std::pair<int, int>> DisabledNodeLabel;
  std::set<int> DisabledFaces;
  std::map<std::pair<int, int>, int> InternalECount;
  std::map<std::pair<int, int>, int> ExternalECount;

  void InitDataStructures() {
    InternalECount.clear();
    ExternalECount.clear();
    std::vector<std::pair<int, int>> InternalTangentConstrEdges;

    // for each node the list of faces it is tangent with TB seen
    for (size_t i = 0; i < M3Graph.InternalTangentConstrEdges.size(); i++)
      InternalECount[M3Graph.InternalTangentConstrEdges[i]] = 0;
    for (size_t i = 0; i < M3Graph.BorderTangentConstrEdges.size(); i++)
      ExternalECount[M3Graph.BorderTangentConstrEdges[i]] = 0;
  }

  bool HasDisabledFace(const int &IndexN) {
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF > 0);
    assert(IndexF < InternalFaces.size());
    if (DisabledFaces.count(IndexF) > 0)
      return true;
    return false;
  }

  bool HasDisabledLabel(const int &IndexN) {
    for (size_t l = 0; l < 3; l++) {
      if (DisabledNodeLabel.count(std::pair<int, int>(IndexN, l)) > 0)
        return true;
    }
    return false;
  }

  std::vector<std::pair<int, int>> GetNodeEdgeVPairs(const int &IndexN) {
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF > 0);
    assert(IndexF < InternalFaces.size());
    std::vector<std::pair<int, int>> result;
    for (size_t e = 0; e < 3; e++) {
      int IndexV0 = InternalFaces[IndexF][e];
      int IndexV1 = InternalFaces[IndexF][(e + 1) % 3];
      assert(IndexV0 != IndexV1);
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      result.push_back(key);
    }
    return result;
  }

  bool HasDisabledEdge(const int &IndexN) {
    std::vector<std::pair<int, int>> res = GetNodeEdgeVPairs(IndexN);
    for (size_t i = 0; i < res.size(); i++) {
      std::pair<int, int> key = res[i];

      bool IsInternal = (InternalECount.count(key) > 0);
      bool IsExternal = (ExternalECount.count(key) > 0);
      assert(!(IsInternal && IsExternal));
      assert(IsInternal || IsExternal);
      if (IsInternal) {
        assert(InternalECount[key] <= 4);
        return (InternalECount[key] != 4);
      }
      if (IsExternal) {
        assert(InternalECount[key] <= 2);
        return (InternalECount[key] != 2);
      }
    }
  }

  bool CanAddNode(const int &IndexN) {
    if (HasDisabledEdge(IndexN))
      return false;
    if (HasDisabledFace(IndexN))
      return false;
    if (HasDisabledLabel(IndexN))
      return false;
    return true;
  }

  void UpdateDisabledLabels(const int &AddedNode) {
    Field::FrameFieldM3Graph::ElemOfNode(AddedNode);
    assert(M3Graph.IsEnabled(AddedNode));

    // this one should have been already added
    assert(AddedNode >= 0);
    assert(AddedNode < IsSelectedN.size());
    assert(IsSelectedN[AddedNode]);

    // 3 labels for each node
    for (size_t l = 0; l < 3; l++) {
      // get all active neighbours
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      std::vector<bool> Reverse;
      M3Graph.ActiveNeigh(AddedNode, l, IndexNeigh, WeightNeigh, Reverse);

      // check if there is already another one active in M3
      bool to_disable = false;
      for (size_t i = 0; i < IndexNeigh.size(); i++) {
        int NextNode = IndexNeigh[i];
        assert(NextNode != AddedNode);
        if (IsSelectedN[AddedNode]) {
          to_disable = true;
          break;
        }
      }
      if (!to_disable)
        continue;

      // disable itself and neighbours
      std::pair<int, int> NodeL0(AddedNode, l);
      assert(DisabledNodeLabel.count(NodeL0) == 0);
      DisabledNodeLabel.insert(NodeL0);

      for (size_t i = 0; i < IndexNeigh.size(); i++) {
        int NextNode = IndexNeigh[i];
        assert(NextNode != AddedNode);
        int IndexL0, IndexL1;
        M3Graph.WhichLabel(AddedNode, NextNode, IndexL0, IndexL1);
        assert(IndexL0 == l);
        assert(IndexL1 != -1);
        std::pair<int, int> NodeL1(NextNode, IndexL1);
        DisabledNodeLabel.insert(NodeL1);
      }
    }
  }

  void UpdateDisabledFaces(const int &AddedNode) {
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(AddedNode);
    assert(IndexF > 0);
    assert(IndexF < InternalFaces.size());
    DisabledFaces.insert(IndexF);
  }

  void UpdateDisabledEdges(const int &AddedNode) {
    std::vector<std::pair<int, int>> res = GetNodeEdgeVPairs(AddedNode);
    for (size_t i = 0; i < res.size(); i++) {
      std::pair<int, int> key = res[i];
      bool IsInternal = (InternalECount.count(key) > 0);
      bool IsExternal = (ExternalECount.count(key) > 0);
      assert(!(IsInternal && IsExternal));
      assert(IsInternal || IsExternal);
      if (IsInternal)
        InternalECount[key]++;
      if (IsExternal)
        ExternalECount[key]++;
    }
  }

  void AddNodeToFoil(const int &AddedNode) {

    assert(AddedNode>=0);
    assert(AddedNode<IsSelectedN.size());
    assert(CanAddNode(AddedNode));

    IsSelectedN[AddedNode]=true;
    SelectedN.push_back(AddedNode);

    // update the disabled labels
    UpdateDisabledLabels(AddedNode);

    // then update the disabled faces
    UpdateDisabledFaces(AddedNode);

    // update the edge count, 3 labels
    UpdateDisabledEdges(AddedNode);
  }

  bool ExpandStep() {
    // for (size_t i=0;i<SelectedN.size();)
  }

public:

  void Init() { InitDataStructures(); }

  void Solve(const AF_ExtractParam &ExtrP, std::vector<int> &ChosenF,
             std::vector<bool> &MatchedSeed) {

    MatchedSeed = std::vector<bool>(ExtrP.Seeds.size(), false);
    ChosenF.clear();
    InitSeedDistance(ExtrP);

    IsSelectedN = std::vector<bool>(M3Graph.NumNodes(), false);
    SelectedN.clear();

    bool has_expanded = false;

    for (size_t i=0;i<ExtrP.Seeds.size();i++)
      AddNodeToFoil(ExtrP.Seeds[i]);

    do {
      has_expanded = ExpandStep();
    } while (has_expanded);
  }

  SheetExtractorAFM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                     std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif