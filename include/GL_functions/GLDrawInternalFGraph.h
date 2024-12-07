#ifndef GL_DRAW_INTERNAL_F_GRAPH
#define GL_DRAW_INTERNAL_F_GRAPH

// #include "../field_graph/path_functions.h"
#include "../field_graph/graph.h"
#include "../field_graph/graph_frame_field.h"
// #include "../frame_field.h"
#include "../point3.h"
// #include "GLDrawField.h"
#include "../color_functions.h"
#include "GLDrawBase.h"
#include "GLDrawField.h"
#include "GLDrawMesh.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

namespace GLGraphDraw {

template <class ScalarType> class GLFrameUndirectedInternalFGraph {

  static void DrawNodeFace(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &InternalFaces,
                           const int Node_To_Draw,
                           const Geo::Point3<ScalarType> PolyColorLine,
                           const ScalarType PolyLineWidth) {
    int ElemIdx = Field::FrameFieldM3Graph::ElemOfNode(Node_To_Draw);

    assert(ElemIdx >= 0);
    assert(ElemIdx < InternalFaces.size());

    if (PolyLineWidth > 0) {
      std::vector<Geo::Point3<ScalarType>> Poly;
      for (size_t i = 0; i < InternalFaces[ElemIdx].size(); i++) {
        int IndexV = InternalFaces[ElemIdx][i];
        Geo::Point3<ScalarType> Pos = VertPos[IndexV];
        Poly.push_back(Pos);
      }
      // add the first one to close the polyline
      Poly.push_back(Poly[0]);
      GLDraw::GLDrawPolyline(Poly, PolyColorLine, PolyLineWidth);
    }
  }

  static void DrawNode(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const int Node_To_Draw, const Geo::Point3<ScalarType> ColFace,
      const Geo::Point3<ScalarType> ColDir, const ScalarType ThickFace,
      const ScalarType ThickDir) {
    int DirIdx = Field::FrameFieldM3Graph::DirOfNode(Node_To_Draw);
    int ElemIdx = Field::FrameFieldM3Graph::ElemOfNode(Node_To_Draw);

    // draw the face
    DrawNodeFace(VertPos, InternalFaces, Node_To_Draw, ColFace, ThickFace);

    // draw the direction of the frame
    ScalarType scaleVal =
        AvgEdgeLen<ScalarType>(VertPos, InternalFaces, ElemIdx);
    GLDraw::DrawDirectionFrame<ScalarType>(InternalFacesFrameField[ElemIdx],
                                           DirIdx, scaleVal, 0, ColDir, ColDir,
                                           ThickDir, ThickDir);
  }

  static Geo::Point3<ScalarType>
  GetColorNode(int Node_To_Draw,
               const Geo::UndirectedNodeWeightedGraph<ScalarType> &FaceGraph) {

    Geo::Point3<ScalarType> Col0 = Geo::ColorRGBRamp<ScalarType>(
        (ScalarType)1, (ScalarType)0, FaceGraph.Weight(Node_To_Draw));

    if (!FaceGraph.IsEnabled(Node_To_Draw))
      Col0 = Geo::Point3<ScalarType>(0, 0, 0);

    if (FaceGraph.NodeType(Node_To_Draw) == Geo::TangentDomain)
      Col0 = Geo::Point3<ScalarType>(0.8, 0.8, 0.8);

    return Col0;
  }

  static void DrawNeighbor(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::UndirectedNodeWeightedGraph<ScalarType> &FaceGraph,
      int IndexNode, int LocalN, bool DrawDir) {
    bool IsActiveL = FaceGraph.IsEnabledL(IndexNode, LocalN);
    int IndexNeigh = FaceGraph.NeighNode(IndexNode, LocalN);
    Geo::Point3<ScalarType> Col = GetColorNode(IndexNeigh, FaceGraph);
    if (!IsActiveL)
      Col = Geo::Point3<ScalarType>(1, 0, 0);

    if (DrawDir)
      DrawNode(VertPos, InternalFaces, InternalFacesFrameField, IndexNeigh, Col,
               Col, 5, 5);
    else
      DrawNode(VertPos, InternalFaces, InternalFacesFrameField, IndexNeigh, Col,
               Col, 5, 0);
  }

public:
  struct InternalFDrawParam {
    ScalarType lineFaceSize;
    ScalarType lineDirSize;
    Geo::Point3<ScalarType> ColFaceLine;
    Geo::Point3<ScalarType> ColDirLine;
    bool DrawNeigh;
    int CurrNeigh;
    bool DrawNeghDir;
    // bool DrawWeight;

    InternalFDrawParam() {
      lineFaceSize = 10;
      lineDirSize = 10;
      ColFaceLine = Geo::Point3<ScalarType>(0, 1, 0);
      ColDirLine = Geo::Point3<ScalarType>(0, 0, 1);
      DrawNeigh = true;
      CurrNeigh = -1;
      DrawNeghDir = false;
      // DrawWeight = true;
    }
  };

  static void DrawNode(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::UndirectedNodeWeightedGraph<ScalarType> &FaceGraph,
      const int Node_To_Draw, InternalFDrawParam &DParam) {
    // std::cout<<"W:"<<FaceGraph.Weight(Node_To_Draw)<<std::endl;

    Geo::Point3<ScalarType> Col0 = GetColorNode(Node_To_Draw, FaceGraph);

    DrawNode(VertPos, InternalFaces, InternalFacesFrameField, Node_To_Draw,
             Col0, Col0, DParam.lineFaceSize, DParam.lineDirSize);
    //compute normals to show

    int FaceIdx0 = Field::FrameFieldM3Graph::ElemOfNode(Node_To_Draw);
    GLDraw::DrawNormal<ScalarType>(VertPos,InternalFaces,FaceIdx0);

    std::vector<Geo::Segment3<ScalarType>> SegNorm;

    for (size_t i = 0; i < FaceGraph.NumNeigh(Node_To_Draw); i++)
    {
      int IndexNNode=FaceGraph.NeighNode(Node_To_Draw,i);
      int FaceIdx1 = Field::FrameFieldM3Graph::ElemOfNode(IndexNNode);
      Geo::Point3<ScalarType> NF0,NF1;
      Geo::ComputeCoherentNormal(VertPos,InternalFaces[FaceIdx0],InternalFaces[FaceIdx1],NF0,NF1);
      //NF1=Geo::Normal(VertPos,InternalFaces,FaceIdx1);
      ScalarType AvgL1=AvgEdgeLen(VertPos,InternalFaces,FaceIdx1);
      Geo::Point3<ScalarType> Bary1=ComputeFaceBarycenter<ScalarType>(VertPos,InternalFaces,FaceIdx1);
    
      Geo::Segment3<ScalarType> SegNorm1(Bary1,Bary1+NF1*AvgL1*0.3);
      SegNorm.push_back(SegNorm1);
    }

    if (DParam.DrawNeigh) {
      if (DParam.CurrNeigh == -1) {
        for (size_t i = 0; i < FaceGraph.NumNeigh(Node_To_Draw); i++) {
          DrawNeighbor(VertPos, InternalFaces, InternalFacesFrameField,
                       FaceGraph, Node_To_Draw, i, false);
        }
      } else {
        int IndexNeigh=DParam.CurrNeigh;
        int NeighNode= FaceGraph.NeighNode(Node_To_Draw,IndexNeigh);
        //int FaceIdx1 = Field::FrameFieldM3Graph::ElemOfNode(NeighNode);
        //GLDraw::DrawNormal<ScalarType>(VertPos,InternalFaces,FaceIdx1);

        GLDraw::GLDrawSegment<ScalarType>(SegNorm[IndexNeigh].P(0),SegNorm[IndexNeigh].P(1),100);

        // Geo::Point3<ScalarType> Dir0=Geo::Normal(VertPos,InternalFaces,FaceIdx0);
        // Geo::Point3<ScalarType> Dir1=SegNorm[IndexNeigh].P(1)-SegNorm[IndexNeigh].P(0);
        // Dir1.Normalize();
        // std::cout<<"Angle:"<<Geo::AngleDeg(Dir0,Dir1)<<std::endl;

        DrawNeighbor(VertPos, InternalFaces, InternalFacesFrameField, FaceGraph,
                     Node_To_Draw, DParam.CurrNeigh, true);
      }
    }
  }
};

}; // namespace GLGraphDraw

#endif
