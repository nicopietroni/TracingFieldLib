#ifndef SMOOTH_MESH_SAMPLING_DIR
#define SMOOTH_MESH_SAMPLING_DIR


#include "point3.h"
#include "fitting_points.h"
#include "triangle3.h"
#include <sys/_types/_size_t.h>
#include <vector>

namespace Geo {

// template <class ScalarType>
// class SmoothSampleDirections
// {

//   void FindTargetPositions(const int IndexV,
//                            const std::vector<int> &KernelF,
//                            const std::vector<Geo::Point3<ScalarType> > &VertPos,
//                            const std::vector<std::vector<int> > &Faces,
//                            const Geo::Point3<ScalarType> &StartDir,
//                            const ScalarType &AngleInt,
//                            const int SampleStep,
//                            Geo::Point3<ScalarType> &TargetP)
//   {
//     //first collect all triangles to test
//     std::vector<Geo::Traingle3<ScalarType> > NeighTris;
//     std::vector<Geo::Point3<ScalarType> > FitPoints;
//     for (size_t i=0;i<KernelF.size();i++)
//     {
//       Geo::Traingle3<ScalarType> T=TriangleOfFace(VertPos,Faces,KernelF[i]);
//       NeighTris.push_back(T);
//       FitPoints.push_back(Barycenter(T));
//     }
//     FitPoints.push_back(VertPos[IndexV]);
//     Geo::Plane3<ScalarType> FitPlane;

//     //find the fitting plane
//     FittingPoint<ScalarType>::FitPlaneToPointSet(FitPoints,FitPlane);

//     //get the normal of the plane
//     Geo::Point3<ScalarType> N=FitPlane.Direction();
    
//     //flatten the direction to the plane

//     //then rotate on one and other direction

//     // reproject them on triangles

//   }

//   public:

//   void Smooth(
//     std::vector<Geo::Point3<ScalarType> > &VertNorm,
//     const std::vector<std::vector<int> > &Faces, 
//     int NumIte,
//     int SamplingDir=32,
//     int Level=4,
//     const std::vector<int> &FixedVert = std::vector<int>(),
//     const std::vector<std::pair<int, int> > &Features =
//     std::vector<std::pair<int, int> >()) 
// {
  
// }

//};

} // namespace Geo
#endif
