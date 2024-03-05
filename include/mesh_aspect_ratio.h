#ifndef MESH_ASPECT_RATIO
#define MESH_ASPECT_RATIO

#include <mesh_functions.h>
#include <triangle3.h>

namespace Geo{

enum QType{QTRadii,QTMeanRatio};

template<class ScalarType>
ScalarType TriMeshQuality(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                          const std::vector<std::vector<int> > &Faces,
                          const int &IndexF,QType QT=QTRadii)
{
   assert(IndexF>=0);
   assert(IndexF<Faces.size());
   assert(Faces[IndexF].size()==3);
   Geo::Traingle3<ScalarType> Tri=Geo::TriangleOfFace(VertPos,Faces,IndexF);
   if (QT==QTRadii)
       return (Geo::QualityRadii(Tri));
   else
       return (Geo::QualityMeanRatio(Tri));
}

template<class ScalarType>
ScalarType TriMeshAvgQuality(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                             const std::vector<std::vector<int> > &Faces,
                             QType QT=QTRadii)
{
    ScalarType AvgQ=0;
    for (size_t i=0;i<Faces.size();i++)
        AvgQ+=TriMeshQuality(VertPos,Faces,i,QT);
    AvgQ/=Faces.size();
    return (AvgQ);
}

template<class ScalarType>
ScalarType TriMeshAvgQuality(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                             const std::vector<std::vector<int> > &Faces,
                             const std::vector<int> &FacesIdx,QType QT=QTRadii)
{
    ScalarType AvgQ=0;
    for (size_t i=0;i<FacesIdx.size();i++)
        AvgQ+=TriMeshQuality(VertPos,Faces,FacesIdx[i],QT);
    AvgQ/=FacesIdx.size();
    return (AvgQ);
}

//   assert(IsAllTris(VertPos,Faces));
//   for (size)

}


#endif
