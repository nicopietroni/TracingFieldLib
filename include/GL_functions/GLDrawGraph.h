#ifndef GL_DRAW_GRAPH
#define GL_DRAW_GRAPH

#include <vector>
#include <GL/glew.h>
#include "../point3.h"
#include "../frame_field.h"
#include "GLDrawMesh.h"
#include "GLDrawField.h"
#include "../field_graph/graph.h"
#include "../field_graph/graph_frame_field.h"
#include "../field_graph/graph_cross_field.h"

namespace GLGraphDraw
{

    template <class ScalarType>
    class GLCrossFieldGraph
    {

    public:
        struct PathDrawParam
        {
            bool DrawPath;
            ScalarType pointSize;
            ScalarType lineSize;
            Geo::Point3<ScalarType> ColPoint;
            Geo::Point3<ScalarType> ColLine;
            ScalarType moveVal;

            PathDrawParam()
            {
                DrawPath = true;
                pointSize = 0;
                lineSize = 5;
                ColPoint = Geo::Point3<ScalarType>(0, 0, 1);
                ColLine = Geo::Point3<ScalarType>(0, 0, 0);
                moveVal = (ScalarType)0.1;
            }
        };

        static void DrawNodes(const Geo::TracingGraph<ScalarType> &TGraph,
                              const std::vector<Field::CrossF<ScalarType>> &Cross,
                              const std::vector<int> &Nodes_To_Draw,
                              const ScalarType scaleVal,
                              const PathDrawParam &PDParam = PathDrawParam(),
                              const Geo::Point3<ScalarType> &ColLine = Geo::Point3<ScalarType>(0, 1, 1))
        {
            glPushAttrib(GL_ALL_ATTRIB_BITS);

            glDepthRange(0.0, 0.9999);
            if (PDParam.DrawPath)
            {
                glDisable(GL_LIGHTING);
                glLineWidth(10);
                GLDraw::glColor(ColLine);
                glBegin(GL_LINE_STRIP);
                for (size_t i = 0; i < Nodes_To_Draw.size(); i++)
                {
                    int currNode = Nodes_To_Draw[i];
                    int currVert = Field::CrossFieldGraph::ElemOfNode(currNode);
                    GLDraw::glVertex(Cross[currVert].Pos);
                }
                glEnd();
            }
            for (size_t i = 0; i < Nodes_To_Draw.size(); i++)
            {
                int currNode = Nodes_To_Draw[i];
                assert(currNode >= 0);
                // assert(currNode<(int)TGraph.N.size());
                assert(currNode < (int)TGraph.NumNodes());
                int CrossD = Field::CrossFieldGraph::DirOfNode(currNode);
                int CrossI = Field::CrossFieldGraph::ElemOfNode(currNode);
                assert(CrossD < 4);
                assert(CrossI < (int)Cross.size());

                float Psize = (float)PDParam.pointSize;
                float Lsize = (float)PDParam.lineSize;
                // if (i == 0)
                // {
                //     Lsize *= 2;
                //     Psize = Lsize * 2;
                // }
                GLDraw::DrawDirectionCross<ScalarType>(Cross[CrossI], CrossD, scaleVal, PDParam.moveVal,
                                                       PDParam.ColPoint, PDParam.ColLine,
                                                       Psize, Lsize);
            }
            glDepthRange(0.0, 1.0);
            glPopAttrib();
        }

        static void DrawPaths(const Geo::TracingGraph<ScalarType> &TGraph,
                              const std::vector<Field::CrossF<ScalarType>> &Cross,
                              const std::vector<std::vector<int>> &Paths_To_Draw,
                              const ScalarType scaleVal,
                              const PathDrawParam &PDParam = PathDrawParam(),
                              const std::vector<Geo::Point3<ScalarType>> &PathColors =
                                  std::vector<Geo::Point3<ScalarType>>())
        {
            for (size_t i = 0; i < Paths_To_Draw.size(); i++)
            {
                if (PathColors.size() == Paths_To_Draw.size())
                    DrawNodes(TGraph, Cross, Paths_To_Draw[i], scaleVal, PDParam, PathColors[i]);
                else
                    DrawNodes(TGraph, Cross, Paths_To_Draw[i], scaleVal, PDParam);
            }
        }

        static void DrawSelectedNodes(const Geo::TracingGraph<ScalarType> &TGraph,
                                      const std::vector<Field::CrossF<ScalarType>> &Cross,
                                      const int IndexSel,
                                      const ScalarType scaleVal,
                                      const ScalarType lineSize = 5,
                                      Geo::Point3<ScalarType> ColLine = Geo::Point3<ScalarType>(0.5, 0.9, 0.5))
        {
            std::vector<int> Nodes_To_Draw;
            for (size_t i = 0; i < TGraph.NumNodes(); i++)
                if (TGraph.IsSelected(i, IndexSel))
                    Nodes_To_Draw.push_back(i);

            DrawNodes(TGraph, Cross, Nodes_To_Draw, scaleVal, lineSize * 2, lineSize, Geo::Point3<ScalarType>(0, 0, 1), ColLine);
        }

        static void DrawNodesOfType(const Geo::TracingGraph<ScalarType> &TGraph,
                                    const std::vector<Field::CrossF<ScalarType>> &Cross,
                                    const typename Geo::NType TNode,
                                    const ScalarType scaleVal, const PathDrawParam &PDParam = PathDrawParam(),
                                    const Geo::Point3<ScalarType> &ColLine = Geo::Point3<ScalarType>(0, 1, 1))
        {
            std::vector<int> Nodes_To_Draw;
            for (size_t i = 0; i < TGraph.NumNodes(); i++)
            {
                if (TGraph.NodeType(i) != TNode)
                    continue;
                Nodes_To_Draw.push_back(i);
            }
            PathDrawParam PDParam1 = PDParam;
            PDParam1.DrawPath = false;
            DrawNodes(TGraph, Cross, Nodes_To_Draw, scaleVal, PDParam1, ColLine);
        }
    };

    template <class ScalarType>
    class GLFrameFieldGraph
    {

    public:
        static void DrawNodeTetras(const std::vector<std::vector<int>> &TetConnectivity,
                                   const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                   const std::vector<int> &Nodes_To_Draw,
                                   Geo::Point3<ScalarType> ColTet = Geo::Point3<ScalarType>(0, 0, 1))
        {
            std::vector<int> TetraIdx;
            for (size_t i = 0; i < Nodes_To_Draw.size(); i++)
            {
                int IndexTet = Field::FrameFieldGraph::PolyOfNode(Nodes_To_Draw[i]);
                assert(IndexTet >= 0);
                assert(IndexTet < (int)TetConnectivity.size());
                TetraIdx.push_back(IndexTet);
            }

            GLDraw::GLDrawTetras(TetConnectivity, VertPos, TetraIdx, ColTet);
        }

        static void DrawNodes(const Geo::TracingGraph<ScalarType> &TGraph,
                              const std::vector<Field::FrameField<ScalarType>> &Frames,
                              const std::vector<int> &Nodes_To_Draw,
                              const ScalarType scaleVal,
                              const ScalarType pointSize = 10,
                              const ScalarType lineSize = 5,
                              Geo::Point3<ScalarType> ColPoint = Geo::Point3<ScalarType>(0, 0, 1),
                              Geo::Point3<ScalarType> ColLine = Geo::Point3<ScalarType>(0.5, 0.9, 0.5),
                              const ScalarType moveVal = (ScalarType)0.1)
        {
            glPushAttrib(GL_ALL_ATTRIB_BITS);

            for (size_t i = 0; i < Nodes_To_Draw.size(); i++)
            {
                int currNode = Nodes_To_Draw[i];
                assert(currNode >= 0);
                // assert(currNode<(int)TGraph.N.size());
                assert(currNode < (int)TGraph.NumNodes());
                int FrameD = Field::FrameFieldGraph::FrameOfNode(currNode);
                int FrameI = Field::FrameFieldGraph::PolyOfNode(currNode);
                assert(FrameD < 6);
                assert(FrameI < (int)Frames.size());

                ScalarType curr_moveVal = moveVal;

                GLDraw::DrawDirectionFrame<ScalarType>(Frames[FrameI], FrameD, scaleVal, curr_moveVal,
                                                       ColPoint, ColLine, (float)pointSize, (float)lineSize);
            }
            glPopAttrib();
        }

        static void DrawPathNodes(const Geo::TracingGraph<ScalarType> &TGraph,
                                  const std::vector<Field::FrameField<ScalarType>> &Frames,
                                  const std::vector<int> &Path,
                                  const ScalarType scaleVal,
                                  const ScalarType LineWidth = 3,
                                  Geo::Point3<ScalarType> ColPath = Geo::Point3<ScalarType>(0, 0, 1))
        {
            std::vector<Geo::Point3<ScalarType>> PolyL;
            for (size_t i = 0; i < Path.size(); i++)
            {
                int IndexP = Field::FrameFieldGraph::PolyOfNode(Path[i]);
                PolyL.push_back(Frames[IndexP].Pos);
            }

            GLDraw::GLDrawPolyline<ScalarType>(PolyL, ColPath, LineWidth);

            DrawNodes(TGraph, Frames, Path, scaleVal); //,(ScalarType)0,false,false,false);

            glPopAttrib();
        }

        static void DrawGraphFrames(const Geo::TracingGraph<ScalarType> &TGraph,
                                    const std::vector<Field::FrameField<ScalarType>> &Frames,
                                    const std::vector<int> &To_DrawFrames,
                                    const ScalarType scaleVal,
                                    const ScalarType moveVal = (ScalarType)0.1,
                                    Geo::Point3<ScalarType> Col0 = Geo::Point3<ScalarType>(0.5, 0.9, 0.5),
                                    Geo::Point3<ScalarType> Col1 = Geo::Point3<ScalarType>(0.5, 0.5, 0.9),
                                    bool RedSink = true)
        {
            std::vector<int> Nodes_To_Draw;

            for (size_t i = 0; i < To_DrawFrames.size(); i++)
                for (size_t j = 0; j < 6; j++)
                {
                    int NodeI = Field::FrameFieldGraph::NodeOfPolyFrame((int)To_DrawFrames[i], (int)j);
                    Nodes_To_Draw.push_back(NodeI);
                }

            DrawNodes(TGraph, Frames, Nodes_To_Draw, scaleVal, moveVal, Col0, Col1, RedSink);
        }
    };

};

#endif
