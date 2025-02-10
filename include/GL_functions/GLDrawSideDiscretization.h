#ifndef GL_DRAW_SIDE_DISCRETIZATION
#define GL_DRAW_SIDE_DISCRETIZATION

#include "GLDrawBase.h"
#include "GLDrawGraph.h"
#include <GL/glew.h>
#include "GLDrawMesh.h"
#include <iostream>
#include <vector>
#include "../patch_side_discretization.h"
#include "color_functions.h"

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
  
};

} // namespace GLGraphDraw

#endif
