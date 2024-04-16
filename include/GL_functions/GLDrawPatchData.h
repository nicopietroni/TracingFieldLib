#ifndef GL_DRAW_PATCH_DATA
#define GL_DRAW_PATCH_DATA

#include <vector>
#include <GL/glew.h>
#include "../field_graph/patch_managing.h"
#include "../field_graph/path_sampling.h"
#include "GLDrawBase.h"
#include "GLDrawGraph.h"

#include "color_functions.h"


namespace GLGraphDraw
{

template <class ScalarType>
class GLDrawPatchData
{

public:

    static void DrawCorners(Geo::PatchManaging<ScalarType> &PM,
                            const int &IndexPatch,
                            const ScalarType PointSize=20)
    {
        assert(IndexPatch<PM.PData.VisualCornerType.size());
        assert(PM.PData.VisualCorners[IndexPatch].size()==PM.PData.VisualCornerType[IndexPatch].size());

        glPushAttrib(GL_ALL_ATTRIB_BITS);

        glDepthRange(0.0,0.9999);

        glDisable(GL_LIGHTING);
        glPointSize((GLfloat)PointSize);


        glBegin(GL_POINTS);
        for (size_t i=0;i<PM.PData.VisualCorners[IndexPatch].size();i++)
        {

            Geo::Point3<ScalarType> ColPoints(0,0,0);
            typename Geo::CornerType CType=PM.PData.VisualCornerType[IndexPatch][i];
            switch(CType) {
              case Geo::CTypeConvex:
                ColPoints=Geo::Point3<ScalarType>(0,1,1);
                break;
              case Geo::CTypeConcave:
                ColPoints=Geo::Point3<ScalarType>(1,0,1);
                break;
             case Geo::CTypeSpike:
                ColPoints=Geo::Point3<ScalarType>(1,1,0);
                break;
             default:
                ColPoints=Geo::Point3<ScalarType>(0,1,0);
            }

            GLDraw::glColor(ColPoints);
            Geo::Point3<ScalarType> Pos=PM.PData.VisualCorners[IndexPatch][i];
            GLDraw::glVertex(Pos);
        }
        glEnd();

        glDepthRange(0.0,1.0);
        glPopAttrib();
    }

    static void DrawCorners(Geo::PatchManaging<ScalarType> &PM,
                            const ScalarType PointSize=20)
    {
        for (size_t i=0;i<PM.PData.VisualCorners.size();i++)
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


        for (size_t i=0;i<PM.PData.VisualPaths[IndexPatch].size();i++)
        {
            Geo::Point3<ScalarType> ColLine=Geo::ScatterColor<ScalarType>((int)PM.PData.VisualPaths[IndexPatch].size(),i);
            GLDraw::glColor(ColLine);
            glBegin(GL_LINE_STRIP);
            for (size_t j=0;j<PM.PData.VisualPaths[IndexPatch][i].size();j++)
            {
                Geo::Point3<ScalarType> Pos=PM.PData.VisualPaths[IndexPatch][i][j];
                GLDraw::glVertex(Pos);
            }
             glEnd();
        }

        glDepthRange(0.0,1.0);
        glPopAttrib();
    }

    static void DrawSides(Geo::PatchManaging<ScalarType> &PM,const ScalarType LineSize=10)
    {
        for (size_t i=0;i<PM.PData.VisualCorners.size();i++)
        {
            DrawSides(PM,i,LineSize);
        }
    }

    static void DrawCandidates(const Geo::CrossPathSampling<ScalarType> &PathSampl,
                                ScalarType ScaleSize,
                                int curr_cant=-1)
    {
        for (size_t i=0;i<PathSampl.Candidates.size();i++)
        {
            if ((curr_cant!=-1)&&(curr_cant!=i))continue;
            int IndexC=PathSampl.Candidates[i].second;
            assert(IndexC<PathSampl.SampledPaths.size());
            assert(IndexC<PathSampl.SampledPaths_Status.size());
            Geo::Point3<ScalarType> PathColors=Geo::ColorRGBRamp((ScalarType)0,(ScalarType)PathSampl.Candidates.size(),(ScalarType)i);
            typename GLGraphDraw::GLCrossFieldGraph<ScalarType>::PathDrawParam PDparam;
            GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodes(PathSampl.Graph(),PathSampl.CrossF(),PathSampl.
                                                                  SampledPaths[IndexC],ScaleSize,
                                                                  PDparam,PathColors);
        }
    }
};

}

#endif
