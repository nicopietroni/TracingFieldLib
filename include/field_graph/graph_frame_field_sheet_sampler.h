#ifndef FRAME_FIELD_GRAPH_SHEET_SAMPLER
#define FRAME_FIELD_GRAPH_SHEET_SAMPLER

#include "../connectivity.h"
#include "../frame_field.h"
#include "../point3.h"
#include "graph.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>
// #include "graph.h"
// #include "graph_tracing.h"
// #include "./frame_field.h"
#include "../mesh_sampling.h"
#include "graph_frame_field.h"
// #include <connectivity.h>
// #define TRACE_DEBUG

template <class ScalarType> class SheetSampler {
  typedef typename Geo::UndirectedNodeWeightedGraph<ScalarType> GraphType;

public:
  std::vector<int> ChosenNodes;

  void SampleInitialPoints(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraph, const ScalarType &Radius) {
    std::vector<int> Selected;
    Geo::MeshSampling<ScalarType>::SampleRadius3D(VertPos, Radius, Selected);
    // then for each vertex get the internal faces
    assert(VertFrameField.size() == VertPos.size());
    assert(InternalFFrameField.size() == InternalFaces.size());
    ChosenNodes.clear();

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(InternalFaces, VF, VFIndex);
    // for each (selected) vertex
    for (size_t i = 0; i < Selected.size(); i++) {
      int IndexV = Selected[i];
      assert(IndexV >= 0);
      assert(IndexV < VF.size());
      assert(IndexV < VertPos.size());

      for (size_t j = 0; j < 3; j++) {
        int ChosenNode = -1;
        ScalarType BestW = 1;
        Geo::Point3<ScalarType> DirTest =
            VertFrameField[IndexV].GetDirection(j);
        for (size_t f = 0; f < VF[IndexV].size(); f++) {
          int IndexF = VF[IndexV][f];
          assert(IndexF >= 0);
          assert(IndexF < InternalFaces.size());
          assert(IndexF < InternalFFrameField.size());
          int IndexDir = InternalFFrameField[IndexF].BestM3DirectionI(DirTest);
          int IndexNode =
              Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
          if (!InternalFacesGraph.IsEnabled(IndexNode))
            continue;
          if (InternalFacesGraph.NodeType(IndexNode) == Geo::TangentDomain)
            continue;
          ScalarType CurrW = InternalFacesGraph.Weight(IndexNode);
          if (CurrW > BestW)
            continue;
          BestW = CurrW;
          ChosenNode = IndexNode;
        }
        if (ChosenNode == -1)
          continue;
        ChosenNodes.push_back(ChosenNode);
      }
    }
    std::cout << "** Selected " << ChosenNodes.size() << " initial Seeds"
              << std::endl;
  }
};
#endif
