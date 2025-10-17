#ifndef TRI_REMESH_PARAMETER
#define TRI_REMESH_PARAMETER

#include "../box3.h"
#include "../point3.h"
#include "../pos3.h"
#include "basic_priority_function.h"
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> struct TriRemParam {
  ScalarType TargetL;
  typename Geo::ImproveAspectRatioPriority<ScalarType>::AspectRatioType AType;
  ScalarType feature_angle;
  ScalarType max_error;
  int steps;
  bool iterative;
  bool adaptive;
  bool writeMSG;
  ScalarType minEdgeRatio;
  ScalarType maxEdgeRatio;

  std::vector<std::vector<bool> > IsFaceEFeature;
  std::vector<bool> IsVertCorner;
  std::vector<ScalarType> VertSize;
  //std::vector<bool> Exclude;
  // std::vector<std::pair<int,int> > VertPairFeatures;
  // std::vector<int> CornersIdx;

  bool collapse_small;
  bool refine_big;
  bool flip_valence;
  bool flip_quadrics;
  bool smooth_mesh;

  TriRemParam() {
    TargetL = -1;
    AType = Geo::ImproveAspectRatioPriority<ScalarType>::RadiusRatio;
    feature_angle = 60;
    max_error = -1;
    steps = 5;
    iterative = false;
    adaptive = false;
    writeMSG = false;
    minEdgeRatio = 0.75;
    maxEdgeRatio = 3;

    collapse_small=true;
    refine_big=true;
    flip_valence=true;
    flip_quadrics=true;
    smooth_mesh=true;
  }
};

} // namespace Geo

#endif