#ifndef INITIALIZE_M3_GRAPH
#define INITIALIZE_M3_GRAPH

#include "../frame_field.h"
#include "../point3.h"
#include "best_rotation.h"
#include "connectivity.h"
#include "graph.h"
#include "graph_frame_field.h"
#include "m3_frame_graph.h"
#include "label_graph.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "polyhedra_mesh_functions.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <set>
#include <stdio.h>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class UndirectedGraphFuctions {

  static bool CanConnect(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<int> &Face0,
                         const std::vector<int> &Face1) //,
  // const int E0,const int E1)
  {
    Geo::Point3<ScalarType> NF0, NF1;
    ComputeCoherentNormal(VertPos, Face0, Face1, NF0, NF1);
    // Point3<ScalarType> PF0_0=VertPos[Face0[0]];
    // Point3<ScalarType> PF0_1=VertPos[Face0[1]];
    // Point3<ScalarType> PF0_2=VertPos[Face0[2]];

    // Point3<ScalarType> PF1_0=VertPos[Face1[0]];
    // Point3<ScalarType> PF1_1=VertPos[Face1[1]];
    // Point3<ScalarType> PF1_2=VertPos[Face1[2]];

    // //orient coherently
    // assert(Face0.size()==3);
    // assert(Face1.size()==3);

    // int VE0_F0=Face0[E0];
    // int VE1_F0=Face0[(E0+1)%3];

    // int VE0_F1=Face1[E1];
    // int VE1_F1=Face1[(E1+1)%3];
    // if (VE0_F0==VE0_F1)
    // {
    //   assert(VE1_F0==VE1_F1);
    //   std::swap(PF1_0,PF1_1);
    // }
    // else {
    //   assert(VE0_F0==VE1_F1);
    //   assert(VE1_F0==VE0_F1);
    // }
    // Point3<ScalarType> NF0=Normal(PF0_0,PF0_1,PF0_2);
    // Point3<ScalarType> NF1=Normal(PF1_0,PF1_1,PF1_2);
    NF0.Normalize();
    NF1.Normalize();
    return (Geo::AngleDeg(NF0, NF1) < 90);
  }

  static bool IsDeadEnd(UndirectedNodeWeightedGraph<ScalarType> &Facegraph,
                        int IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < Facegraph.NumNodes());
    std::vector<std::vector<int>> LabelNodes;
    // one label for each edge
    LabelNodes.resize(3);

    for (size_t j = 0; j < Facegraph.NumNeigh(IndexN); j++) {

      if (!Facegraph.IsEnabledL(IndexN, j))
        continue;
      int IndexNeigh = Facegraph.NeighNode(IndexN, j);

      if (!Facegraph.IsEnabled(IndexNeigh))
        continue;

      // get label adjacency
      int IndexLab = Facegraph.NeighLabel(IndexN, j);
      assert(IndexLab >= 0);
      assert(IndexLab < 3);

      LabelNodes[IndexLab].push_back(IndexNeigh);
    }
    for (size_t j = 0; j < 3; j++) {
      if (Facegraph.IsExitLabel(IndexN, j))
        continue;
      if (LabelNodes[j].size() == 0)
        return true;
    }
    return false;
  }

public:
  static void
  DeactivateTangentNodes(UndirectedNodeWeightedGraph<ScalarType> &Facegraph) {
    std::vector<int> IdxNode;
    Facegraph.GetNodeOfType(TangentDomain, IdxNode);
    for (size_t i = 0; i < IdxNode.size(); i++)
      Facegraph.DisableNode(IdxNode[i]);
  }

  static void DeactivateInternalDeadEnds(
      UndirectedNodeWeightedGraph<ScalarType> &Facegraph) {
    int DeactivatedNum = 0;
    bool has_deactivated = false;
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
    std::cout << "Deactivated " << DeactivatedNum << " nodes" << std::endl;
  }

  static void InitFromInternalFacesFrameField(
      UndirectedNodeWeightedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<std::vector<int>> &FT,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      ScalarType MaxAngle = 60) {
    std::vector<Point3<ScalarType>> FaceNormals;
    Geo::ComputeFaceNormals(VertPos, TetraInternalFaces, FaceNormals);

    // create three instances for each face
    Facegraph.AllocateNodes(TetraInternalFaces.size() * 3);

    // set the faces on the surface as tangent
    std::vector<int> ExtFaces;
    for (size_t i = 0; i < FT.size(); i++) {
      assert(FT[i].size() > 0);
      if (FT[i].size() > 1)
        continue;

      // otherwise is an external face
      ExtFaces.push_back(i);
      std::vector<int> IndexN;
      Field::FrameFieldM3Graph::NodesOfElem(i, IndexN);
      for (size_t j = 0; j < IndexN.size(); j++)
        Facegraph.NodeType(IndexN[j]) = TangentDomain;
    }

    std::map<std::pair<int, int>, std::vector<Geo::Point3<ScalarType>>>
        ExtEdgesDir;
    for (size_t i = 0; i < ExtFaces.size(); i++) {
      int IndexF = ExtFaces[i];
      assert(IndexF >= 0);
      assert(IndexF < TetraInternalFaces.size());
      int sizeF = TetraInternalFaces[IndexF].size();
      Geo::Point3<ScalarType> ExitDir = FaceNormals[IndexF];

      for (size_t j = 0; j < sizeF; j++) {
        int IndexV0 = TetraInternalFaces[IndexF][j];
        int IndexV1 = TetraInternalFaces[IndexF][(j + 1) % sizeF];
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));
        ExtEdgesDir[KeyE].push_back(ExitDir);
      }
    }

    // then cycle over all faces
    for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
      int IndexF = i;
      int sizeF = TetraInternalFaces[IndexF].size();

      // for each edge
      for (size_t j = 0; j < sizeF; j++) {
        int IndexV0 = TetraInternalFaces[IndexF][j];
        int IndexV1 = TetraInternalFaces[IndexF][(j + 1) % sizeF];
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));

        // if not on the border
        if (ExtEdgesDir.count(KeyE) == 0)
          continue;

        /// get all exit directions
        for (size_t j = 0; j < ExtEdgesDir[KeyE].size(); j++) {
          // get the direction
          Geo::Point3<ScalarType> ExitDIr = ExtEdgesDir[KeyE][j];

          // then get best oriented node in M3
          assert(IndexF < FacesFrameField.size());
          int IndexDir = FacesFrameField[IndexF].BestM3DirectionI(ExitDIr);

          // get the index of node
          int IndexN =
              Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
          assert(IndexN >= 0);
          assert(IndexN < Facegraph.NumNodes());
          Facegraph.NodeType(IndexN) = TangentDomain;
        }
      }
    }

    // set the exit label
    for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
      int IndexF = i;
      int sizeF = TetraInternalFaces[IndexF].size();

      // for each edge
      for (size_t e = 0; e < sizeF; e++) {
        int IndexV0 = TetraInternalFaces[IndexF][e];
        int IndexV1 = TetraInternalFaces[IndexF][(e + 1) % sizeF];
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));

        // if not on the border
        if (ExtEdgesDir.count(KeyE) == 0)
          continue;

        // for each node
        std::vector<int> IndexN;
        Field::FrameFieldM3Graph::NodesOfElem(IndexF, IndexN);
        for (size_t i = 0; i < IndexN.size(); i++) {
          if (Facegraph.NodeType(IndexN[i]) == TangentDomain)
            continue;
          Facegraph.SetExitLabel(IndexN[i], e);
        }
      }
    }

    // initialize the weight of nodes and disable in case
    // thy exceed a thresold
    for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
      int DirIdx = Field::FrameFieldM3Graph::DirOfNode(i);
      int ElemIdx = Field::FrameFieldM3Graph::ElemOfNode(i);
      assert(ElemIdx >= 0);
      assert(ElemIdx < FacesFrameField.size());
      assert(ElemIdx < FaceNormals.size());
      assert(DirIdx >= 0);
      assert(DirIdx < 3);
      Geo::Point3<ScalarType> DirFrame =
          FacesFrameField[ElemIdx].GetDirection(DirIdx);
      Geo::Point3<ScalarType> NormFace = FaceNormals[ElemIdx];

      // invert if needed
      if (Geo::Dot(DirFrame, NormFace) < 0)
        DirFrame = -DirFrame;

      ScalarType CurrAngle = Geo::AngleDeg(DirFrame, NormFace);
      // std::cout<<"Angle "<<CurrAngle<<std::endl;
      assert(CurrAngle <= 90);
      Facegraph.Weight(i) = pow(CurrAngle / (ScalarType)90, 4);

      if (CurrAngle > MaxAngle)
        Facegraph.DisableNode(i);
    }

    // compute the FF
    std::vector<std::vector<std::vector<int>>> NextF, NextE;
    Geo::getFFAdjacencyNonManifold(TetraInternalFaces, NextF, NextE);

    // for each face
    for (size_t i = 0; i < NextF.size(); i++) {
      // for each M3 direction
      std::vector<int> Nodes;
      Field::FrameFieldM3Graph::NodesOfElem(i, Nodes);
      int IndexF0 = i;

      for (size_t n = 0; n < Nodes.size(); n++) {
        // for each surrounding edge
        int IndexN0 = Nodes[n];
        assert(IndexN0 >= 0);
        assert(IndexN0 < Facegraph.NumNodes());
        int IndexDir0 = Field::FrameFieldM3Graph::DirOfNode(IndexN0);

        Field::FrameField<ScalarType> Frame0 = FacesFrameField[i];

        for (size_t e = 0; e < NextF[i].size(); e++) {
          // for each adjacent face along that edge
          int Label0 = (int)e;
          for (size_t f = 0; f < NextF[i][e].size(); f++) {
            int IndexF1 = NextF[i][e][f];
            int Label1 = NextE[i][e][f];

            if (IndexF1 == -1)
              continue; // no neighbours
            Geo::Point3<ScalarType> TargetD = Frame0.GetDirection(IndexDir0);
            Field::FrameField<ScalarType> Frame1 = FacesFrameField[IndexF1];
            int IndexDir1 = Frame1.BestM3DirectionI(TargetD);
            int IndexN1 =
                Field::FrameFieldM3Graph::NodeOfElemDir(IndexF1, IndexDir1);

            // check if they are connected
            bool IsEnabled = CanConnect(VertPos, TetraInternalFaces[IndexF0],
                                        TetraInternalFaces[IndexF1]); //,
            // Label0,Label1);
            // add the connection
            Facegraph.AddSymmLink(IndexN0, IndexN1, Label0, Label1, IsEnabled);
          }
        }
      }
    }
    // then make it symmetric
    Facegraph.MakeUniqueConnection();

    // print statistics
    Facegraph.WriteStats();
  }

  // // static void GetExternalFaceNodes(Geo::LabelNeighOrientedGraph<ScalarType>
  // // &Facegraph,
  // //                                       const std::vector<int> &FaceMap)
  // // {
  // //   // set all nodes of external faces as tangent
  // //   for (size_t i = 0; i < FaceMap.size(); i++) {
  // //     if (FaceMap[i] == -1)
  // //       continue;
  // //     std::vector<int> IndexN;
  // //     Field::FrameFieldM3Graph::NodesOfElem(i, IndexN);

  // //     // set the one on surfaces as tangent
  // //     for (size_t j = 0; j < IndexN.size(); j++) {
  // //       assert(IndexN[j] >= 0);
  // //       assert(IndexN[j] < Facegraph.NumNodes());
  // //       Facegraph.NodeType(IndexN[j]) = TangentDomain;
  // //     }
  // //   }
  // // }

  // static void
  // GetExternalEdges(const std::vector<std::vector<int>> &TetraInternalFaces,
  //                  const std::vector<int> &FaceMap,
  //                  std::set<std::pair<int, int>> &ExtEdges) {
  //   for (size_t i = 0; i < FaceMap.size(); i++) {
  //     if (FaceMap[i] == -1)
  //       continue;
  //     // then set directions for edges
  //     int sizeF = TetraInternalFaces[i].size();
  //     assert(sizeF == 3);
  //     for (size_t j = 0; j < sizeF; j++) {
  //       int IndexV0 = TetraInternalFaces[i][j];
  //       int IndexV1 = TetraInternalFaces[i][(j + 1) % sizeF];
  //       std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
  //                                std::max(IndexV0, IndexV1));
  //       ExtEdges.insert(KeyE);
  //     }
  //   }
  // }

  // static void
  // GetTangentNodes(const std::vector<int> &VertMap,
  //                 const std::vector<Geo::Point3<ScalarType>> &VertNorm,
  //                 const std::vector<int> &FaceMap,
  //                 const std::vector<std::vector<int>> &TetraInternalFaces,
  //                 const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
  //                 std::vector<int> &TangentNodes) {
  //   TangentNodes.clear();

  //   // first get the nodes encoding external faces
  //   for (size_t i = 0; i < FaceMap.size(); i++) {
  //     if (FaceMap[i] == -1)
  //       continue;
  //     std::vector<int> IndexN;
  //     Field::FrameFieldM3Graph::NodesOfElem(i, IndexN);
  //     TangentNodes.insert(TangentNodes.end(), IndexN.begin(), IndexN.end());
  //   }

  //   // std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNorm;
  //   // ComputeNormals<ScalarType>(SurfVertPos, SurfFaces, FaceNormals,
  //   // VertNorm);

  //   // then get all nodes having one vertex on surface
  //   for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
  //     int IndexF = i;
  //     int sizeF = TetraInternalFaces[IndexF].size();
  //     assert(sizeF == 3);
  //     for (size_t j = 0; j < sizeF; j++) {
  //       int IndexV = TetraInternalFaces[IndexF][j];

  //       if (VertMap[IndexV] == -1)
  //         continue; // not an external vert

  //       // get the normal on the surface
  //       int IndexVSurf = VertMap[IndexV];
  //       assert(IndexVSurf >= 0);
  //       assert(IndexVSurf < VertNorm.size());
  //       Geo::Point3<ScalarType> ExitDIr = VertNorm[IndexVSurf];

  //       // get the most aligned direction
  //       int IndexDir = FacesFrameField[IndexF].BestDirectionI(ExitDIr);

  //       // map between 0 and 2 for M3
  //       if (IndexDir >= 3)
  //         IndexDir = Field::OppositeFrameDir(IndexDir);

  //       // then find the node
  //       int IndexNode =
  //           Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
  //       TangentNodes.push_back(IndexNode);
  //     }
  //   }
  // }

  // static void
  // GetExitNodes(Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
  //              const std::vector<int> &FaceMap,
  //              const std::vector<std::vector<int>> &TetraInternalFaces,
  //              const std::vector<std::vector<std::vector<int>>> &NextF,
  //              std::vector<std::pair<int, int>> &ExitNodeLabel) {

  //   ExitNodeLabel.clear();
  //   // get the external edges
  //   std::set<std::pair<int, int>> ExtEdges;
  //   GetExternalEdges(TetraInternalFaces, FaceMap, ExtEdges);

  //   // cycle over all faces
  //   for (size_t i = 0; i < TetraInternalFaces.size(); i++) {

  //     int IndexF = i;

  //     // if is an external surface then
  //     //  it is already set as tangent
  //     if (FaceMap[IndexF] != -1)
  //       continue;

  //     // it can have maximum 3 edges
  //     assert(NextF[IndexF].size() == 3);

  //     // for each edge
  //     int sizeF = NextF[IndexF].size();
  //     for (size_t e = 0; e < NextF[IndexF].size(); e++) {
  //       int IndexLabel = e;
  //       int IndexV0 = TetraInternalFaces[IndexF][e];
  //       int IndexV1 = TetraInternalFaces[IndexF][(e + 1) % sizeF];
  //       std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
  //                                std::max(IndexV0, IndexV1));

  //       // check if the edge is on the border
  //       if (ExtEdges.count(KeyE) == 0)
  //         continue;

  //       // if on the border simply set as exit if not already tangent
  //       for (size_t d = 0; d < 3; d++) {
  //         int NodeIdx = Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, d);
  //         if (Facegraph.NodeType(NodeIdx) == TangentDomain)
  //           continue;
  //         // Facegraph.SetExit(NodeIdx, IndexLabel);
  //         ExitNodeLabel.push_back(std::pair<int, int>(NodeIdx, IndexLabel));
  //       }
  //     }
  //   }
  // }

  // static void
  // GetNodeNodeConnection(const std::vector<int> &FaceMap,
  //                       const std::vector<std::vector<int>> &TetraInternalFaces,
  //                       const std::vector<std::vector<std::vector<int>>> &NextF,
  //                       const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
  //                       std::vector<std::pair<int, int> > &Connections,
  //                       std::vector<int> &ConnectionsLabel) {
  //   Connections.clear();
  //   ConnectionsLabel.clear();

  //   // then set the exit nodes/label
  //   for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
  //     int IndexF0 = i;
  //     assert(IndexF0 < FacesFrameField.size());

  //     // if is an external surface then
  //     //  it has no neigh
  //     if (FaceMap[IndexF0] != -1)
  //       continue;

  //     // it can have maximum 3 edges
  //     assert(NextF[IndexF0].size() == 3);

  //     // for each edge
  //     int sizeF = NextF[IndexF0].size();
  //     for (size_t e = 0; e < NextF[IndexF0].size(); e++) {
  //       int IndexLabel = e;
        
  //       // cycle over all the other faces
  //       for (size_t j = 0; j < NextF[IndexF0][e].size(); j++) {
  //         int IndexF1 = NextF[IndexF0][e][j];

  //         // can't be border we tested before
  //         assert(IndexF1 >= 0);
  //         assert(IndexF1 < TetraInternalFaces.size());
  //         assert(IndexF1 < FacesFrameField.size());

  //         // if is an external surface then
  //         // can't ne neigh
  //         if (FaceMap[IndexF1] != -1)
  //           continue;

  //         // then cycle all over possible directions
  //         for (size_t d = 0; d < 3; d++) {
            
  //           // get best direction to collow
  //           int NeighD = FacesFrameField[IndexF0].FollowDirection(
  //                         FacesFrameField[IndexF1], d);

  //           if (NeighD >= 3)
  //             NeighD = Field::OppositeFrameDir(NeighD);

  //           int Node0 = Field::FrameFieldGraph::NodeOfPolyFrame(IndexF0, d);
  //           int Node1 = Field::FrameFieldGraph::NodeOfPolyFrame(IndexF1, NeighD);

  //           // ScalarType Weight=1;
  //           // Facegraph.AddNeighbour(Node0, Node1, IndexLabel, Weight);
  //           std::pair<int,int> KeyConn(std::min(Node0,Node1),
  //                                      std::max(Node0,Node1));

  //           Connections.push_back(KeyConn);
  //         }
  //       }
  //     }
  //   }

  //   //then symmetrize and make unique the connections
  //   std::sort(Connections.begin(), Connections.end());
  //   Connections.resize(std::unique(Connections.begin(), Connections.end()) - Connections.begin());

  //   //then find the labels and add the symmetric connnections
  //   std::vector<std::pair<int, int> > SymmConn;
  //   for (size_t i=0;i<Connections.size();i++)
  //   {
  //     int Node0 = Connections[i].first;
  //     int Node1 = Connections[i].second;
  //     int IndexF0=Field::FrameFieldM3Graph::ElemOfNode(Node0);
  //     int IndexF1=Field::FrameFieldM3Graph::ElemOfNode(Node1);
  //     int IndexE0=-1;
  //     int IndexE1=-1;
  //     Geo::SharedE(TetraInternalFaces[IndexF0],
  //                  TetraInternalFaces[IndexF1],
  //                  IndexE0,IndexE1);
  //     assert(IndexE0!=-1);
  //     assert(IndexE1!=-1);

  //     std::pair<int,int> CurrConn=Connections[i];
  //     SymmConn.push_back(CurrConn);
  //     ConnectionsLabel.push_back(IndexE0);

  //     std::swap(CurrConn.first, CurrConn.second);
  //     SymmConn.push_back(CurrConn);
  //     ConnectionsLabel.push_back(IndexE1);     
  //   }
  // }

  // static void
  // InitM3(Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
  //        const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //        const std::vector<std::vector<int>> &TetraInternalFaces,
  //        const std::vector<Geo::Point3<ScalarType>> &SurfVertPos,
  //        const std::vector<std::vector<int>> &SurfFaces,
  //        // const std::vector<std::pair<int,int> > &SurfConvexE,
  //        // const std::vector<std::pair<int,int> > &SurfConcaveE,
  //        const std::vector<Field::FrameField<ScalarType>> &FacesFrameField) {
  //   // get the mapping from face to face
  //   std::vector<int> FaceMap;
  //   int MappedFNum = Geo::FindFacesCorrespondences(
  //       VertPos, TetraInternalFaces, SurfVertPos, SurfFaces, FaceMap);
  //   assert(MappedFNum == SurfFaces.size());

  //   // get the mapping from vert to vert
  //   std::vector<int> VertMap;
  //   int MappedVNum = FindVertCorrespondences(VertPos, SurfVertPos, VertMap);
  //   assert(MappedVNum == SurfVertPos.size());

  //   // get the tangent nodes
  //   std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNorm;
  //   ComputeNormals<ScalarType>(SurfVertPos, SurfFaces, FaceNormals, VertNorm);
  //   std::vector<int> TangentNodes;
  //   GetTangentNodes(VertMap, VertNorm, FaceMap, TetraInternalFaces,
  //                   FacesFrameField,TangentNodes);

  //   // set as tangent nodes
  //   for (size_t i = 0; i < TangentNodes.size(); i++) {
  //     int IndexN = TangentNodes[i];
  //     Facegraph.NodeType(IndexN) = TangentDomain;
  //   }

  //   std::vector<std::vector<std::vector<int>>> NextF, NextE;
  //   Geo::getFFAdjacencyNonManifold(TetraInternalFaces, NextF, NextE);

  //   // get the exit label nodes
  //   std::vector<std::pair<int, int>> ExitNodeLabel;
  //   GetExitNodes(FaceMap, TetraInternalFaces, NextF, ExitNodeLabel);

  //   // set the exit label
  //   for (size_t i = 0; i < ExitNodeLabel.size(); i++) {
  //     Facegraph.SetExit(ExitNodeLabel[i].first, ExitNodeLabel[i].second);
  //   }

  //   //then get the connections
  //   std::vector<std::pair<int, int>> Connections;
  //   std::vector<int> ConnectionsLabel;
  //   GetNodeNodeConnection(FaceMap,TetraInternalFaces,
  //                         NextF,FacesFrameField,
  //                         Connections,ConnectionsLabel);

  //   assert(Connections.size()==ConnectionsLabel.size());
  //   for (size_t i=0;i<Connections.size();i++)
  //   {
  //     int IndexN0=Connections[i].first;
  //     int IndexN1=Connections[i].second;
  //     int IndexLab=ConnectionsLabel[i];
  //     ScalarType W=1;
  //     AddNeighbour(IndexN0,IndexN1,IndexLab,W); 
  //   }
  // }
};

} // namespace Geo
#endif