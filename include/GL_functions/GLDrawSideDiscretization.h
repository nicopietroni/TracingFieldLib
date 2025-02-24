#ifndef GL_DRAW_SIDE_DISCRETIZATION
#define GL_DRAW_SIDE_DISCRETIZATION

#include "GLDrawBase.h"
#include "GLDrawGraph.h"
#include <GL/glew.h>
#include "GLDrawMesh.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include "../patch_side_discretization.h"
#include "color_functions.h"
#include "GLDrawLabel_imgui.h"

namespace GLGraphDraw {

template <class ScalarType> class GLDrawPatchSideDiscretization {

public:

  static void DrawMultipleSide(const Geo::SideDiscretization<ScalarType> &SideD)
  {
    for (size_t i=0;i<SideD.NumPatches();i++)
    for (size_t j=0;j<SideD.NumSides(i);j++)
    {
      std::vector<std::vector<Geo::Point3<ScalarType> > > SubSides;
      SideD.GetSidePoints(i,j,SubSides);
      if (SubSides.size()<=1)continue;
      for (size_t k=0;k<SubSides.size();k++)
        GLDraw::GLDrawPolyline(SubSides[k],
                               Geo::Point3<ScalarType>(1,0,0),
                               (ScalarType)30);
    }
  }
  
  static void DrawSideValue(const Geo::SideDiscretization<ScalarType> &SideD)
  {
    std::vector<Geo::Point3<ScalarType> > VarMidP;
    SideD.GetVariablesMidP(VarMidP);
    
    std::vector<int> VarValue;
    SideD.GetVariablesValues(VarValue);

    std::vector<Geo::Point3<ScalarType> > LabelCol;
    LabelCol.resize(VarMidP.size(),Geo::Point3<ScalarType>(0,0,1));

    std::vector<std::string> LabelText;
    for (size_t i=0;i<VarValue.size();i++)
      LabelText.push_back(std::to_string(VarValue[i]));

    GLDraw::GLDrawLabel<ScalarType>::IMGuiDrawLabels3D(VarMidP,LabelCol,LabelText);
  }

  static void
  DrawSubdivision(const Geo::SideDiscretization<ScalarType> &SideD,
                  bool drawExtremes=true) 
  {
   for (size_t i=0;i<SideD.SampleSidePoints.size();i++)
   {
    Geo::Point3<ScalarType> pointCol = Geo::Point3<ScalarType>(0, 1, 0);
    if (SideD.VarMatching[i]==Geo::SideDiscretization<ScalarType>::NoQuad)
      pointCol = Geo::Point3<ScalarType>(0, 0, 0);
    if (SideD.VarMatching[i]==Geo::SideDiscretization<ScalarType>::OddSubd)
      pointCol = Geo::Point3<ScalarType>(1, 0, 0);
    if (SideD.VarMatching[i]==Geo::SideDiscretization<ScalarType>::NoMatch)
      pointCol = Geo::Point3<ScalarType>(1, 1, 0);
    GLDraw::DrawPoints(SideD.SampleSidePoints[i],pointCol);
   }
  }

  static void
  GLDrawSortedSubdivision(const Geo::SideDiscretization<ScalarType> &SideD
                          ) 
  {
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType> > > > SidePoints;
    SideD.GetSortedSideSamples(SidePoints);
    if (SidePoints.size()==0)return;
    for (size_t i=0;i<1;i++)
    {
      std::vector<Geo::Point3<ScalarType> > PatchSamples;
      for (size_t j=0;j<SidePoints[i].size();j++)
        PatchSamples.insert(PatchSamples.end(),
                            SidePoints[i][j].begin(),
                            SidePoints[i][j].end());

      std::vector<Geo::Point3<ScalarType> > pointCol(PatchSamples.size(),Geo::Point3<ScalarType>(0,0,0));

      for (size_t j=0;j<pointCol.size();j++)
        pointCol[j]=Geo::ColorRGBRamp<ScalarType>((ScalarType)0,
                                                  (ScalarType)PatchSamples.size(),
                                                  (ScalarType)j);

      GLDraw::DrawPoints(PatchSamples,pointCol,50);
    }
  }

};

} // namespace GLGraphDraw

#endif
