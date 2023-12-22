#ifndef GL_DRAW_PATCH_DATA
#define GL_DRAW_PATCH_DATA

#include <vector>
#include <GL/glew.h>
#include "../field_graph/patch_managing.h"
#include "GLDrawBase.h"
#include "color_functions.h"


namespace GLGraphDraw
{

template <class ScalarType>
class GLDrawPatchData
{

public:

    static void DrawCorners(Geo::PatchManaging<ScalarType> &PM,
                            const int &IndexPatch,const ScalarType PointSize=20)
    {
        assert(IndexPatch<PM.VisualCornerType.size());
        assert(PM.VisualCorners[IndexPatch].size()==PM.VisualCornerType[IndexPatch].size());

        glPushAttrib(GL_ALL_ATTRIB_BITS);

        glDepthRange(0.0,0.9999);

        glDisable(GL_LIGHTING);
        glPointSize((GLfloat)PointSize);


        glBegin(GL_POINTS);
        for (size_t i=0;i<PM.VisualCorners[IndexPatch].size();i++)
        {

            Geo::Point3<ScalarType> ColPoints(0,0,0);
            typename Geo::PatchManaging<ScalarType>::CornerType CType=PM.VisualCornerType[IndexPatch][i];
            switch(CType) {
              case Geo::PatchManaging<ScalarType>::CTypeConvex:
                ColPoints=Geo::Point3<ScalarType>(0,1,1);
                break;
              case Geo::PatchManaging<ScalarType>::CTypeConcave:
                ColPoints=Geo::Point3<ScalarType>(1,0,1);
                break;
             case Geo::PatchManaging<ScalarType>::CTypeSpike:
                ColPoints=Geo::Point3<ScalarType>(1,1,0);
                break;
             default:
                ColPoints=Geo::Point3<ScalarType>(0,10,0);
            }

            GLDraw::glColor(ColPoints);
            Geo::Point3<ScalarType> Pos=PM.VisualCorners[IndexPatch][i];
            GLDraw::glVertex(Pos);
        }
        glEnd();

        glDepthRange(0.0,1.0);
        glPopAttrib();
    }

    static void DrawCorners(Geo::PatchManaging<ScalarType> &PM,const ScalarType PointSize=20)
    {
        for (size_t i=0;i<PM.VisualCorners.size();i++)
            DrawCorners(PM,i,PointSize);
    }

    static void DrawSides(Geo::PatchManaging<ScalarType> &PM,
                          const int &IndexPatch,
                          const ScalarType LineSize=10)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);

        glDepthRange(0.0,0.9999);

        glDisable(GL_LIGHTING);
        glLineWidth((GLfloat)LineSize);


        for (size_t i=0;i<PM.VisualPaths[IndexPatch].size();i++)
        {
            Geo::Point3<ScalarType> ColLine=Geo::ScatterColor<ScalarType>((int)PM.VisualPaths[IndexPatch].size(),i);
            GLDraw::glColor(ColLine);
            glBegin(GL_LINE_STRIP);
            for (size_t j=0;j<PM.VisualPaths[IndexPatch][i].size();j++)
            {
                Geo::Point3<ScalarType> Pos=PM.VisualPaths[IndexPatch][i][j];
                GLDraw::glVertex(Pos);
            }
             glEnd();
        }

        glDepthRange(0.0,1.0);
        glPopAttrib();
    }

    static void DrawSides(Geo::PatchManaging<ScalarType> &PM,const ScalarType LineSize=10)
    {
        for (size_t i=0;i<PM.VisualCorners.size();i++)
        {
            DrawSides(PM,i,LineSize);
        }
    }
};

}

#endif
