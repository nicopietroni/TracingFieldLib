#ifndef MESH_HAUSDORFF
#define MESH_HAUSDORFF

#include "./Space_Query/triangle_grid_3D.h"
#include "Space_Query/point_grid_3D.h"
#include "point3.h"
#include "segment2.h"

namespace Geo {

template <class ScalarType>
void OneWayHausdorff(const std::vector<Geo::Point3<ScalarType> > &Source_VertPos,
                     const std::vector<std::vector<int> > &Source_Faces,
                     const std::vector<Geo::Point3<ScalarType> > &Target_VertPos,
                     const std::vector<std::vector<int> > &Target_Faces,
                     std::vector<ScalarType> &VertDist)
{
  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(Target_VertPos, Target_Faces);

  VertDist.clear();
  ScalarType maxDistQ = TriGrid.bbox.Diag();
  for (size_t i = 0; i < Source_VertPos.size(); i++) {
    Geo::Point3<ScalarType> Pos = Source_VertPos[i];
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(Target_VertPos, Target_Faces, Pos, maxDistQ, closestF,
                        closestPt);
    ScalarType CurrDist = 0;
    if (closestF == -1)
      CurrDist = maxDistQ;
    else
      CurrDist = (closestPt - Pos).Norm();

    VertDist.push_back(CurrDist);
  }
}

template <class ScalarType>
void OneWayHausdorff(const std::vector<Geo::Point3<ScalarType> > &Source_VertPos,
                     const std::vector<std::vector<int> > &Source_Faces,
                     const std::vector<Geo::Point3<ScalarType> > &Target_VertPos,
                     const std::vector<std::vector<int> > &Target_Faces,
                     ScalarType &AvgDist, ScalarType &MaxDist) {

  std::vector<ScalarType> VertDist;
  OneWayHausdorff(Source_VertPos,Source_Faces,Target_VertPos,Target_Faces,VertDist);
  AvgDist = 0;
  MaxDist = 0;
  for (size_t i = 0; i < VertDist.size(); i++) {
    ScalarType CurrDist=VertDist[i];
    AvgDist += CurrDist;
    MaxDist = std::max(MaxDist, CurrDist);
  }
  AvgDist /= Source_VertPos.size();
}

template <class ScalarType>
void Hausdorff(const std::vector<Geo::Point3<ScalarType> > &VertPos0,
               const std::vector<std::vector<int> > &Faces0,
               const std::vector<Geo::Point3<ScalarType> > &VertPos1,
               const std::vector<std::vector<int> > &Faces1, 
               ScalarType &AvgDist,
               ScalarType &MaxDist) {
  ScalarType AvgDist0,MaxDist0,AvgDist1,MaxDist1;
  OneWayHausdorff(VertPos0,Faces0,VertPos1,Faces1,AvgDist0,MaxDist0);
  OneWayHausdorff(VertPos1,Faces1,VertPos0,Faces0,AvgDist1,MaxDist1);
  AvgDist=(AvgDist0+AvgDist1)/2;
  MaxDist=std::max(MaxDist0,MaxDist1);
}

}; // namespace Geo

#endif
