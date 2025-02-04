#ifndef FRAME_FIELD_GRAPH_SHEET_SAMPLER
#define FRAME_FIELD_GRAPH_SHEET_SAMPLER

#include "../connectivity.h"
#include "../frame_field.h"
#include "../point3.h"
#include "graph.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>
// #include "graph.h"
// #include "graph_tracing.h"
// #include "./frame_field.h"
#include "../mesh_sampling.h"
#include "../normal.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
// #include <connectivity.h>
// #define TRACE_DEBUG

// template <class ScalarType> class SheetSampler {
//   typedef typename Geo::UndirectedNodeWeightedGraph<ScalarType> GraphType;

//   void GetBestNodes(
//       const std::vector<int> &InternalFacesKernel,
//       const Field::FrameField<ScalarType> &VertFrame,
//       const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
//       const GraphType &InternalFacesGraph, std::vector<int> &ChosenNodes) {
//     ChosenNodes.clear();
//     for (size_t j = 0; j < 3; j++) {
//       int ChosenNode = -1;
//       ScalarType BestW = 1;
//       Geo::Point3<ScalarType> DirTest = VertFrame.GetDirection(j);
//       for (size_t f = 0; f < InternalFacesKernel.size(); f++) {
//         int IndexF = InternalFacesKernel[f];
//         assert(IndexF >= 0);
//         assert(IndexF < InternalFFrameField.size());
//         int IndexDir = InternalFFrameField[IndexF].BestM3DirectionI(DirTest);
//         int IndexNode =
//             Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
//         if (!InternalFacesGraph.IsEnabled(IndexNode))
//           continue;
//         if (InternalFacesGraph.NodeType(IndexNode) == Geo::TangentDomain)
//           continue;
//         ScalarType CurrW = InternalFacesGraph.Weight(IndexNode);
//         if (CurrW > BestW)
//           continue;
//         BestW = CurrW;
//         ChosenNode = IndexNode;
//       }
//       if (ChosenNode == -1)
//         continue;
//       ChosenNodes.push_back(ChosenNode);
//     }
//   }

// public:
//   std::vector<int> ChosenNodes;

//   void SampleInitialPoints(
//       const std::vector<Geo::Point3<ScalarType>> &ExtVertPos,
//       const std::vector<Geo::Point3<ScalarType>> &VertPos,
//       const std::vector<std::vector<int>> &InternalFaces,
//       const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
//       const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
//       const GraphType &InternalFacesGraph, const ScalarType &Radius,
//       bool FilterComplete = true) {
//     std::vector<int> Selected;
//     // Geo::MeshSampling<ScalarType>::SampleRadius3D(VertPos, Radius, Selected);
//     Geo::MeshSampling<ScalarType>::SampleVoronoiRadius3D(VertPos, Radius,
//                                                          Selected, 30);
//     std::cout << "** Selected " << Selected.size() << " initial Pos"
//               << std::endl;
//     // then for each vertex get the internal faces
//     assert(VertFrameField.size() == VertPos.size());
//     assert(InternalFFrameField.size() == InternalFaces.size());
//     ChosenNodes.clear();

//     std::vector<std::vector<int>> VF, VFIndex;
//     Geo::getVFAdjacency(InternalFaces, VF, VFIndex);

//     if (FilterComplete) {
//       std::vector<Geo::Point3<ScalarType>> Filtered;
//       std::vector<std::vector<int>> BestNodes;

//       std::set<Geo::Point3<ScalarType>> ExtPosSet(ExtVertPos.begin(),
//                                                   ExtVertPos.end());
//       int InternalNum = 0;
//       for (size_t i = 0; i < VertPos.size(); i++) {
//         if (ExtPosSet.count(VertPos[i]) > 0)
//           continue;
//         InternalNum++;

//         std::vector<int> CurrBest;
//         GetBestNodes(VF[i], VertFrameField[i], InternalFFrameField,
//                      InternalFacesGraph, CurrBest);
//         if (CurrBest.size() < 3)
//           continue;
//         Filtered.push_back(VertPos[i]);
//         BestNodes.push_back(CurrBest);
//       }
//       std::cout << "Filtered " << Filtered.size() << " out of " << InternalNum
//                 << " total " << VertPos.size() << std::endl;
//       ;
//     }

//     // for each (selected) vertex
//     for (size_t i = 0; i < Selected.size(); i++) {
//       int IndexV = Selected[i];
//       std::vector<int> Curr_ChosenNodes;
//       GetBestNodes(VF[IndexV], VertFrameField[IndexV], InternalFFrameField,
//                    InternalFacesGraph, Curr_ChosenNodes);
//       ChosenNodes.insert(ChosenNodes.end(), Curr_ChosenNodes.begin(),
//                          Curr_ChosenNodes.end());
//       // assert(IndexV >= 0);
//       // assert(IndexV < VF.size());
//       // assert(IndexV < VertPos.size());

//       // for (size_t j = 0; j < 3; j++) {
//       //   int ChosenNode = -1;
//       //   ScalarType BestW = 1;
//       //   Geo::Point3<ScalarType> DirTest =
//       //       VertFrameField[IndexV].GetDirection(j);
//       //   for (size_t f = 0; f < VF[IndexV].size(); f++) {
//       //     int IndexF = VF[IndexV][f];
//       //     assert(IndexF >= 0);
//       //     assert(IndexF < InternalFaces.size());
//       //     assert(IndexF < InternalFFrameField.size());
//       //     int IndexDir =
//       //     InternalFFrameField[IndexF].BestM3DirectionI(DirTest); int
//       //     IndexNode =
//       //         Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
//       //     if (!InternalFacesGraph.IsEnabled(IndexNode))
//       //       continue;
//       //     if (InternalFacesGraph.NodeType(IndexNode) == Geo::TangentDomain)
//       //       continue;
//       //     ScalarType CurrW = InternalFacesGraph.Weight(IndexNode);
//       //     if (CurrW > BestW)
//       //       continue;
//       //     BestW = CurrW;
//       //     ChosenNode = IndexNode;
//       //   }
//       //   if (ChosenNode == -1)
//       //     continue;
//       //   ChosenNodes.push_back(ChosenNode);
//       // }
//     }
//     std::cout << "** Selected " << ChosenNodes.size() << " initial Seeds"
//               << std::endl;
//   }
// };

template <class ScalarType> class GraphLabelSheetSampler {
  typedef typename Geo::LabelNeighOrientedGraph<ScalarType> GraphType;

  // int SampleBestM3Dir(
  //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const std::vector<std::vector<int>> &InternalFaces,
  //     const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
  //     const int &M6IndexN) {

  //   std::vector<int> M3Nodes;
  //   Field::FrameFieldGraph::ToM3(M6IndexN,M3Nodes);
  //   assert(M3Nodes.size())
  //   int IndexF = Field::FrameFieldGraph::PolyOfNode(M6IndexN);
  //   int IndexDir = Field::FrameFieldGraph::FrameOfNode(M6IndexN);
  //   assert(IndexF >= 0);
  //   assert(IndexF < InternalFaces.size());
  //   assert(IndexF < InternalFFrameField.size());
  //   Geo::Point3<ScalarType> FaceNorm =
  //       Geo::Normal(VertPos, InternalFaces, IndexF);
  //   FaceNorm.Normalize();

  //   // // get the two other ortho direction
  //   // int IndexDir0 = (IndexDir + 1) % 3;
  //   // int IndexDir1 = (IndexDir + 2) % 3;

  //   // Geo::Point3<ScalarType> DirFr0 =
  //   //     InternalFFrameField[IndexF].GetDirection(IndexDir0);
  //   // Geo::Point3<ScalarType> DirFr1 =
  //   //     InternalFFrameField[IndexF].GetDirection(IndexDir1);

  //   // ScalarType Dot0 = fabs(Geo::Dot(DirFr0, FaceNorm));
  //   // ScalarType Dot1 = fabs(Geo::Dot(DirFr1, FaceNorm));
  //   // if (Dot0 > Dot1)
  //   //   return IndexDir0;
  //   // return IndexDir1;
  // }

  int GetBestM3Node(
      const std::vector<Geo::Point3<ScalarType> > &FaceNorms,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const int &M6IndexN) {

    //get the face and the normal
    int IndexF = Field::FrameFieldGraph::PolyOfNode(M6IndexN);
    assert(IndexF < FaceNorms.size());
    assert(IndexF < InternalFFrameField.size());
    Geo::Point3<ScalarType> FaceNorm = FaceNorms[IndexF];
    FaceNorm.Normalize();

    //get the M3 candidates
    std::vector<int> M3Nodes;
    Field::FrameFieldGraph::ToM3(M6IndexN,M3Nodes);
    assert(M3Nodes.size()==2);

    // get the two other ortho direction
    int IndexDir0 = Field::FrameFieldM3Graph::DirOfNode(M3Nodes[0]);
    int IndexDir1 = Field::FrameFieldM3Graph::DirOfNode(M3Nodes[1]);

    //get the two candidate directions
    Geo::Point3<ScalarType> DirFr0 =
        InternalFFrameField[IndexF].GetDirection(IndexDir0);
    Geo::Point3<ScalarType> DirFr1 =
        InternalFFrameField[IndexF].GetDirection(IndexDir1);

    ScalarType Dot0 = fabs(Geo::Dot(DirFr0, FaceNorm));
    ScalarType Dot1 = fabs(Geo::Dot(DirFr1, FaceNorm));
    if (Dot0 > Dot1)
      return M3Nodes[0];
    return M3Nodes[1];
  }

  // void GetBestNodes(
  //     const std::vector<int> &InternalFacesKernel,
  //     const Field::FrameField<ScalarType> &VertFrame,
  //     const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
  //     const GraphType &InternalFacesGraphM3,
  //     const GraphType &InternalFacesGraphM6, 
  //     std::vector<int> &ChosenNodes) {

  //   ChosenNodes.clear();

  //   std::set<int> AddedF;
  //   for (size_t j = 0; j < 6; j++) {
  //     int ChosenNode = -1;
  //     ScalarType BestW = 1;
  //     Geo::Point3<ScalarType> DirTest = VertFrame.GetDirection(j);
  //     for (size_t f = 0; f < InternalFacesKernel.size(); f++) {
  //       int IndexF = InternalFacesKernel[f];
  //       assert(IndexF >= 0);
  //       assert(IndexF < InternalFFrameField.size());
  //       int IndexDir = InternalFFrameField[IndexF].BestDirectionI(DirTest);
  //       int IndexNode =
  //           Field::FrameFieldGraph::NodeOfPolyFrame(IndexF, IndexDir);

  //       if (!InternalFacesGraph.IsEnabled(IndexNode))
  //         continue;
  //       if (InternalFacesGraph.NodeType(IndexNode) == Geo::TangentDomain)
  //         continue;

  //       // Geo::Point3<ScalarType>
  //       // DirF=InternalFFrameField[IndexF].getDirection(IndexDir);
  //       // Geo::Point3<ScalarType> NormF=Geo::Normal()
  //       // ScalarType CurrW = InternalFFrameField[]
  //       ScalarType CurrW = InternalFacesGraph.NodeAlign[IndexNode];
  //       if (CurrW > BestW)
  //         continue;

  //       if (AddedF.count(IndexF) > 0)
  //         continue;

  //       BestW = CurrW;
  //       ChosenNode = IndexNode;
  //     }
  //     if (ChosenNode == -1)
  //       continue;
  //     AddedF.insert(Field::FrameFieldGraph::PolyOfNode(ChosenNode));
  //     ChosenNodes.push_back(ChosenNode);
  //   }
  // }

  void GetBestNodesM6(
      const std::vector<int> &InternalFacesKernel,
      const Field::FrameField<ScalarType> &VertFrame,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM6, 
      std::vector<int> &ChosenNodes) {

    ChosenNodes.clear();

    std::set<int> AddedF;
    for (size_t j = 0; j < 6; j++) {
      int ChosenNode = -1;
      ScalarType BestW = 1;
      Geo::Point3<ScalarType> DirTest = VertFrame.GetDirection(j);
      for (size_t f = 0; f < InternalFacesKernel.size(); f++) {
        int IndexF = InternalFacesKernel[f];
        assert(IndexF >= 0);
        assert(IndexF < InternalFFrameField.size());
        int IndexDir = InternalFFrameField[IndexF].BestDirectionI(DirTest);
        int IndexNode =
            Field::FrameFieldGraph::NodeOfPolyFrame(IndexF, IndexDir);

        if (!InternalFacesGraphM6.IsEnabled(IndexNode))
          continue;
        if (InternalFacesGraphM6.NodeType(IndexNode) == Geo::TangentDomain)
          continue;
        
        ScalarType CurrW = InternalFacesGraphM6.NodeW[IndexNode];
        if (CurrW > BestW)
          continue;

        if (AddedF.count(IndexF) > 0)
          continue;

        BestW = CurrW;
        ChosenNode = IndexNode;
      }
      if (ChosenNode == -1)
        continue;
      AddedF.insert(Field::FrameFieldGraph::PolyOfNode(ChosenNode));
      ChosenNodes.push_back(ChosenNode);
    }
  }
  
  void GetBestNodesM3(
      const std::vector<int> &InternalFacesKernel,
      const Field::FrameField<ScalarType> &VertFrame,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM3, 
      std::vector<int> &ChosenNodes) {

    ChosenNodes.clear();

    std::set<int> AddedF;
    for (size_t j = 0; j < 3; j++) {
      int ChosenNode = -1;
      ScalarType BestW = std::numeric_limits<ScalarType>::max();
      Geo::Point3<ScalarType> DirTest = VertFrame.GetDirection(j);
      for (size_t f = 0; f < InternalFacesKernel.size(); f++) {
        int IndexF = InternalFacesKernel[f];
        assert(IndexF >= 0);
        assert(IndexF < InternalFFrameField.size());
        int IndexDir = InternalFFrameField[IndexF].BestM3DirectionI(DirTest);
        int IndexNode = Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);

        if (!InternalFacesGraphM3.IsEnabled(IndexNode))
          continue;
        if (InternalFacesGraphM3.NodeType(IndexNode) == Geo::TangentDomain)
          continue;
        
        ScalarType CurrW = InternalFacesGraphM3.NodeW[IndexNode];
        if (CurrW > BestW)
          continue;

        if (AddedF.count(IndexF) > 0)
          continue;

        BestW = CurrW;
        ChosenNode = IndexNode;
      }
      if (ChosenNode == -1)
        continue;
      AddedF.insert(Field::FrameFieldM3Graph::ElemOfNode(ChosenNode));
      ChosenNodes.push_back(ChosenNode);
    }
  }

  void GetCompleteSamplesM6M3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM3,
      const GraphType &InternalFacesGraphM6,
      std::vector<Geo::Point3<ScalarType>> &Filtered,
      std::vector<std::vector<std::pair<int, int>>> &BestNodes) {

    Filtered.clear();
    BestNodes.clear();

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(InternalFaces, VF, VFIndex);
    
    std::vector<Geo::Point3<ScalarType> > FaceNormals;
    Geo::ComputeFaceNormals(VertPos,InternalFaces,FaceNormals);

    // int InternalNum=0;
    for (size_t i = 0; i < VertPos.size(); i++) {
      // //in this case is border
      // if (VerteMap[i] != -1)
      //   continue;
      // InternalNum++;

      std::vector<int> CurrBestM6Nodes;
      GetBestNodesM6(VF[i], VertFrameField[i], 
                   InternalFFrameField,
                   InternalFacesGraphM6,
                   CurrBestM6Nodes);

      if (CurrBestM6Nodes.size() < 6)
        continue;

      std::vector<int> CurrBestM3Nodes;
      for (size_t j = 0; j < CurrBestM6Nodes.size();j++) {
        int BestM3Node=GetBestM3Node(FaceNormals,InternalFFrameField,CurrBestM6Nodes[j]);

        assert(InternalFacesGraphM3.IsEnabled(BestM3Node));
        assert(InternalFacesGraphM3.NodeType(BestM3Node) != Geo::TangentDomain);

        CurrBestM3Nodes.push_back(BestM3Node);
      }

      Filtered.push_back(VertPos[i]);
      assert(CurrBestM3Nodes.size()==CurrBestM6Nodes.size());
      BestNodes.resize(BestNodes.size()+1);
      for (size_t j=0;j<CurrBestM6Nodes.size();j++)
        BestNodes.back().push_back(std::pair<int,int>(CurrBestM6Nodes[j],CurrBestM3Nodes[j]));
    }

    std::cout << "Filtered " << Filtered.size() << " out of " << VertPos.size()
              << std::endl;
  }

  void GetCompleteSamplesM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM3,
      std::vector<Geo::Point3<ScalarType>> &Filtered,
      std::vector<std::vector<int>> &BestM3Nodes,
      bool FilterComplete) {

    Filtered.clear();
    BestM3Nodes.clear();

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(InternalFaces, VF, VFIndex);
    
    // std::vector<Geo::Point3<ScalarType> > FaceNormals;
    // Geo::ComputeFaceNormals(VertPos,InternalFaces,FaceNormals);

    // int InternalNum=0;
    for (size_t i = 0; i < VertPos.size(); i++) {
      // //in this case is border
      // if (VerteMap[i] != -1)
      //   continue;
      // InternalNum++;

      std::vector<int> CurrBestM3Nodes;
      GetBestNodesM3(VF[i], VertFrameField[i], 
                          InternalFFrameField,
                          InternalFacesGraphM3,
                          CurrBestM3Nodes);

      if ((FilterComplete)&&(CurrBestM3Nodes.size() < 3))continue;

      if ((!FilterComplete)&&(CurrBestM3Nodes.size() == 0))continue;

      //if (CurrBestM3Nodes.size() < 3)
      // if (CurrBestM3Nodes.size() ==0)
      //   continue;

      Filtered.push_back(VertPos[i]);
      BestM3Nodes.push_back(CurrBestM3Nodes);
    }

    std::cout << "Filtered " << Filtered.size() << " out of " << VertPos.size()
              << std::endl;
  }

public:

  std::vector<int> ChosenNodesM6;
  std::vector<int> ChosenNodesM3;

  void SampleInitialPointsM6M3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM3,
      const GraphType &InternalFacesGraphM6, 
      const ScalarType &Radius,
      bool FilterComplete = true) {

    std::vector<Geo::Point3<ScalarType>> Filtered;
    std::vector<std::vector<std::pair<int, int>>> BestNodes;
    GetCompleteSamplesM6M3(VertPos, InternalFaces, VertFrameField,
                       InternalFFrameField, InternalFacesGraphM3,
                       InternalFacesGraphM6, 
                       Filtered,
                       BestNodes);
    ChosenNodesM6.clear();
    ChosenNodesM3.clear();

    assert(Filtered.size() == BestNodes.size());
    std::vector<int> Selected;
    // Geo::MeshSampling<ScalarType>::SampleRadius3D(VertPos, Radius, Selected);
    Geo::MeshSampling<ScalarType>::SampleVoronoiRadius3D(VertPos, Radius,
                                                         Selected, 30);
    std::cout << "** Selected " << Selected.size() << " initial Pos"
              << std::endl;

    Geo::PointSpatialIndex<ScalarType> PIndex;
    PIndex.Init(Filtered);

    // for each (selected) vertex
    Geo::Box3<ScalarType> BBox;
    BBox.Init(VertPos);
    ScalarType MaxD = BBox.Diag();

    for (size_t i = 0; i < Selected.size(); i++) {
      Geo::Point3<ScalarType> SelPos = VertPos[Selected[i]];
      int closestP_idx;
      PIndex.GridClosest(Filtered, SelPos, MaxD, closestP_idx);
      assert(closestP_idx >= 0);
      assert(closestP_idx < BestNodes.size());

      for (size_t j=0;j<BestNodes[closestP_idx].size();j++)
      {
        ChosenNodesM6.push_back(BestNodes[closestP_idx][j].first);
        ChosenNodesM3.push_back(BestNodes[closestP_idx][j].second);
      }
      // ChosenNodes.insert(ChosenNodes.end(), BestNodes[closestP].begin(),
      //                    BestNodes[closestP].end());
    }
    assert(ChosenNodesM6.size()==ChosenNodesM3.size());
    std::cout << "** Selected " << ChosenNodesM6.size() << " initial Seeds"
              << std::endl;
  }


  void SampleInitialPointsM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM3,
      const ScalarType &Radius,
      bool FilterComplete = true) {

    std::vector<Geo::Point3<ScalarType>> Filtered;
    std::vector<std::vector<int>> BestNodes;
    GetCompleteSamplesM3(VertPos, InternalFaces, VertFrameField,
                       InternalFFrameField,InternalFacesGraphM3, 
                       Filtered, BestNodes,FilterComplete);
    ChosenNodesM6.clear();
    ChosenNodesM3.clear();

    assert(Filtered.size() == BestNodes.size());
    std::vector<int> Selected;
    // Geo::MeshSampling<ScalarType>::SampleRadius3D(VertPos, Radius, Selected);
    Geo::MeshSampling<ScalarType>::SampleVoronoiRadius3D(VertPos, Radius,
                                                         Selected, 30);
    std::cout << "** Selected " << Selected.size() << " initial Pos"
              << std::endl;

    Geo::PointSpatialIndex<ScalarType> PIndex;
    PIndex.Init(Filtered);

    // for each (selected) vertex
    Geo::Box3<ScalarType> BBox;
    BBox.Init(VertPos);
    ScalarType MaxD = BBox.Diag();

    for (size_t i = 0; i < Selected.size(); i++) {
      Geo::Point3<ScalarType> SelPos = VertPos[Selected[i]];
      int closestP_idx;
      PIndex.GridClosest(Filtered, SelPos, MaxD, closestP_idx);
      assert(closestP_idx >= 0);
      assert(closestP_idx < BestNodes.size());

      ChosenNodesM3.insert(ChosenNodesM3.end(), 
                           BestNodes[closestP_idx].begin(),
                           BestNodes[closestP_idx].end());
    }
    //check
    for (size_t i=0;i<ChosenNodesM3.size();i++)
    {
      if ((!InternalFacesGraphM3.IsEnabled(ChosenNodesM3[i]))||
          (InternalFacesGraphM3.NodeType(ChosenNodesM3[i]) == Geo::TangentDomain))
      {
        std::cout<<"WARNING SELECTED TANGED OR NOT ENABLED NODE"<<std::endl;
        exit(0);
      }
    }

    std::cout << "** Selected " << ChosenNodesM3.size() << " initial Seeds"
              << std::endl;
  }
};
#endif
