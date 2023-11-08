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

namespace GLGraphDraw
{

template <class ScalarType>
void DrawNodeTetras(const std::vector<std::vector<int> > &TetConnectivity,
                    const std::vector<Geo::Point3<ScalarType> > &VertPos,
                    const std::vector<int> &Nodes_To_Draw,
                    Geo::Point3<ScalarType> ColTet=Geo::Point3<ScalarType>(0,0,1))
{
    std::vector<int> TetraIdx;
    for (size_t i=0;i<Nodes_To_Draw.size();i++)
    {
        int IndexTet=PolyOfNode(Nodes_To_Draw[i]);
        assert(IndexTet>=0);
        assert(IndexTet<TetConnectivity.size());
        TetraIdx.push_back(IndexTet);
    }

    GLDraw::GLDrawTetras(TetConnectivity,VertPos,TetraIdx,ColTet);

}

//template <class ScalarType>
//void DrawNodes(const Geo::TracingGraph<ScalarType> &TGraph,
//               const std::vector<Field::FrameField<ScalarType> > &Frames,
//               const std::vector<int> &Nodes_To_Draw,
//               const ScalarType scaleVal,
//               const ScalarType moveVal=(ScalarType)0.1,
//               bool RedSink=true,
//               bool Sel0=true,
//               bool Sel1=true)
//{
//    glPushAttrib(GL_ALL_ATTRIB_BITS);

//    for (size_t i=0;i<Nodes_To_Draw.size();i++)
//    {
//        int currNode=Nodes_To_Draw[i];
//        assert(currNode>=0);
//        assert(currNode<TGraph.N.size());
//        int FrameD=FrameOfNode(currNode);
//        int FrameI=PolyOfNode(currNode);
//        assert(FrameD<6);
//        assert(FrameI<(int)Frames.size());

//        Geo::Point3<ScalarType> Col0=Geo::Point3<ScalarType>(0,0,1);
//        Geo::Point3<ScalarType> Col1=Geo::Point3<ScalarType>(0.5,0.9,0.5);
//        ScalarType curr_moveVal=moveVal;
//        ScalarType pointSize=10.f;
//        ScalarType lineSize=5.f;

//        if ((Sel0)&&(TGraph.IsSelected(currNode,0)))
//        {
//            Col0=Geo::Point3<ScalarType>(1,0,0);
//            Col1=Geo::Point3<ScalarType>(1,0,0);
//            pointSize*=2;
//            lineSize*=2;
//        }
//        if ((Sel1)&&(TGraph.IsSelected(currNode,1)))
//        {
//            Col0=Geo::Point3<ScalarType>(1,0,1);
//            Col1=Geo::Point3<ScalarType>(1,0,1);
//            pointSize*=2;
//            lineSize*=2;
//        }
//        if ((RedSink)&&(TGraph.IsSink((size_t)currNode)))
//        {
//            Col0=Geo::Point3<ScalarType>(1,0,0);
//            Col1=Geo::Point3<ScalarType>(1,0,0);
//            pointSize/=2;
//            lineSize=1;
//            curr_moveVal/=2;
//        }
//        GLDraw::DrawDirection<ScalarType>(Frames[FrameI],FrameD,scaleVal,curr_moveVal,
//                                          Col0,Col1,pointSize,lineSize);

//    }
//    glPopAttrib();
//}

template <class ScalarType>
void DrawNodes(const Geo::TracingGraph<ScalarType> &TGraph,
               const std::vector<Field::FrameField<ScalarType> > &Frames,
               const std::vector<int> &Nodes_To_Draw,
               const ScalarType scaleVal,
               const ScalarType pointSize=10,
               const ScalarType lineSize=5,
               Geo::Point3<ScalarType> ColPoint=Geo::Point3<ScalarType>(0,0,1),
               Geo::Point3<ScalarType> ColLine=Geo::Point3<ScalarType>(0.5,0.9,0.5),
               const ScalarType moveVal=(ScalarType)0.1)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    for (size_t i=0;i<Nodes_To_Draw.size();i++)
    {
        int currNode=Nodes_To_Draw[i];
        assert(currNode>=0);
        assert(currNode<TGraph.N.size());
        int FrameD=FrameOfNode(currNode);
        int FrameI=PolyOfNode(currNode);
        assert(FrameD<6);
        assert(FrameI<(int)Frames.size());

        ScalarType curr_moveVal=moveVal;

        GLDraw::DrawDirection<ScalarType>(Frames[FrameI],FrameD,scaleVal,curr_moveVal,
                                          ColPoint,ColLine,(float)pointSize,(float)lineSize);

    }
    glPopAttrib();
}

template <class ScalarType>
void DrawPathNodes(const Geo::TracingGraph<ScalarType> &TGraph,
                   const std::vector<Field::FrameField<ScalarType> > &Frames,
                   const std::vector<int> &Path,
                   const ScalarType scaleVal,
                   const ScalarType LineWidth=3,
                   Geo::Point3<ScalarType> ColPath=Geo::Point3<ScalarType>(0,0,1))
{
    std::vector<Geo::Point3<ScalarType> > PolyL;
    for (size_t i=0;i<Path.size();i++)
    {
        int IndexP=PolyOfNode(Path[i]);
        PolyL.push_back(Frames[IndexP].Pos);
    }

    GLDraw::GLDrawPolyline<ScalarType>(PolyL,ColPath,LineWidth);

    DrawNodes<ScalarType>(TGraph,Frames,Path,scaleVal);//,(ScalarType)0,false,false,false);

    glPopAttrib();
}

template <class ScalarType>
void DrawGraphFrames(const Geo::TracingGraph<ScalarType> &TGraph,
                     const std::vector<Field::FrameField<ScalarType> > &Frames,
                     const std::vector<int> &To_DrawFrames,
                     const ScalarType scaleVal,
                     const ScalarType moveVal=(ScalarType)0.1,
                     Geo::Point3<ScalarType> Col0=Geo::Point3<ScalarType>(0.5,0.9,0.5),
                     Geo::Point3<ScalarType> Col1=Geo::Point3<ScalarType>(0.5,0.5,0.9),
                     bool RedSink=true)
{
    std::vector<int> Nodes_To_Draw;

    for (size_t i=0;i<To_DrawFrames.size();i++)
        for (size_t j=0;j<6;j++)
        {
            int NodeI=NodeOfPolyFrame((int)To_DrawFrames[i],(int)j);
            Nodes_To_Draw.push_back(NodeI);
        }

    DrawNodes(TGraph,Frames,Nodes_To_Draw,scaleVal,moveVal,Col0,Col1,RedSink);
    //        int currF=To_Draw[i];
    //        assert(currF>=0);
    //        assert(currF<Frames.size());
    //        for (size_t j=0;j<6;j++)
    //        {
    //            Point3<ScalarType> Col0=Point3<ScalarType>(0,0,1);
    //            Point3<ScalarType> Col1=Point3<ScalarType>(0.5,0.9,0.5);
    //            int NodeI=NodeOfPolyFrame((int)currF,(int)j);
    //            ScalarType curr_moveVal=moveVal;
    //            ScalarType pointSize=10.f;
    //            ScalarType lineSize=5.f;
    //            if ((RedSink)&&(TGraph.IsSink((size_t)NodeI)))
    //            {
    //                Col0=Point3<ScalarType>(1,0,0);
    //                Col1=Point3<ScalarType>(1,0,0);
    //                pointSize/=2;
    //                lineSize=1;
    //                curr_moveVal/=2;
    //            }
    //            GLDraw::DrawDirection<ScalarType>(Frames[currF],j,scaleVal,curr_moveVal,Col0,Col1,pointSize,lineSize);
    //        }
    //    }
    //    glPopAttrib();
}
};

#endif
