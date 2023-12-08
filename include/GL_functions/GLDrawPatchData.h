#ifndef GL_DRAW_PATCH_DATA
#define GL_DRAW_PATCH_DATA

#include <vector>
#include <GL/glew.h>
#include "../field_graph/patch_managing.h"
#include "GLDrawBase.h"

//#include "../point3.h"
//#include "../frame_field.h"
//#include "GLDrawMesh.h"
//#include "GLDrawField.h"
//#include "../field_graph/graph.h"
//#include "../field_graph/graph_frame_field.h"
//#include "../field_graph/graph_cross_field.h"

namespace GLGraphDraw
{

template <class ScalarType>
class GLDrawPatchData
{

public:

    //        static void DrawCorners(const std::vector<Geo::Point3<ScalarType> > &VertPos,
    //                                Geo::PatchManaging<ScalarType> &PM,
    //                                const int &IndexPatch,const ScalarType PointSize=20,
    //                                const Geo::Point3<ScalarType> &ColPoints=Geo::Point3<ScalarType>(0,1,1))
    //        {
    //            glPushAttrib(GL_ALL_ATTRIB_BITS);

    //            glDepthRange(0.0,0.9999);

    //            glDisable(GL_LIGHTING);
    //            glPointSize((GLfloat)PointSize);

    //            GLDraw::glColor(ColPoints);
    //            glBegin(GL_POINTS);
    //            for (size_t i=0;i<PM.PatchCorners[IndexPatch].size();i++)
    //                for (size_t j=0;j<PM.PatchCorners[IndexPatch][i].size();j++)
    //                {
    //                    int IndexV=PM.PatchCorners[IndexPatch][i][j];
    //                    Geo::Point3<ScalarType> Pos=VertPos[IndexV];
    //                    GLDraw::glVertex(Pos);
    //                }
    //            glEnd();

    //            glDepthRange(0.0,1.0);
    //            glPopAttrib();
    //        }

    //        static void DrawCorners(const std::vector<Geo::Point3<ScalarType> > &VertPos,
    //                                Geo::PatchManaging<ScalarType> &PM,const ScalarType PointSize=20,
    //                                const Geo::Point3<ScalarType> &ColPoints=Geo::Point3<ScalarType>(0,1,1))
    //        {
    //            for (size_t i=0;i<PM.PatchCorners.size();i++)
    //                DrawCorners(VertPos,PM,i,PointSize,ColPoints);
    //        }

    static void DrawCorners(Geo::PatchManaging<ScalarType> &PM,
                            const int &IndexPatch,const ScalarType PointSize=20,
                            const Geo::Point3<ScalarType> &ColPoints=Geo::Point3<ScalarType>(0,1,1))
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);

        glDepthRange(0.0,0.9999);

        glDisable(GL_LIGHTING);
        glPointSize((GLfloat)PointSize);

        GLDraw::glColor(ColPoints);
        glBegin(GL_POINTS);
        for (size_t i=0;i<PM.VisualCorners[IndexPatch].size();i++)
        {
            Geo::Point3<ScalarType> Pos=PM.VisualCorners[IndexPatch][i];
            GLDraw::glVertex(Pos);
        }
        glEnd();

        glDepthRange(0.0,1.0);
        glPopAttrib();
    }

    static void DrawCorners(Geo::PatchManaging<ScalarType> &PM,const ScalarType PointSize=20,
                            const Geo::Point3<ScalarType> &ColPoints=Geo::Point3<ScalarType>(0,1,1))
    {
        for (size_t i=0;i<PM.VisualCorners.size();i++)
            DrawCorners(PM,i,PointSize,ColPoints);
    }

    static void DrawSides(Geo::PatchManaging<ScalarType> &PM,
                          const int &IndexPatch,
                          const ScalarType LineSize=10,
                          const Geo::Point3<ScalarType> &ColLine=Geo::Point3<ScalarType>(0,1,0))
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);

        glDepthRange(0.0,0.9999);

        glDisable(GL_LIGHTING);
        glLineWidth((GLfloat)LineSize);

        GLDraw::glColor(ColLine);

        for (size_t i=0;i<PM.VisualPaths[IndexPatch].size();i++)
        {
            glBegin(GL_LINE_STRIP);
            for (size_t j=0;j<PM.VisualPaths[IndexPatch][i].size();j++)
            {
                Geo::Point3<ScalarType> Pos=PM.VisualPaths[IndexPatch][i];
                GLDraw::glVertex(Pos);
            }
             glEnd();
        }

        glDepthRange(0.0,1.0);
        glPopAttrib();
    }

    static void DrawSides(Geo::PatchManaging<ScalarType> &PM,const ScalarType LineSize=10,
                          const Geo::Point3<ScalarType> &ColLine=Geo::Point3<ScalarType>(0,1,0))
    {
        for (size_t i=0;i<PM.VisualCorners.size();i++)
            DrawSides(PM,i,LineSize,ColLine);
    }
};

}

#endif
