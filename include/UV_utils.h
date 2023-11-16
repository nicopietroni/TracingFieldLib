#ifndef UV_UTILS
#define UV_UTILS

#include <point2.h>
#include <cross_field.h>
#include <Space_Query/triangle_grid_2D.h>

namespace Geo{

template <class ScalarType >
class UVUtils{

//        static void SampleUVIntegralLines(const std::vector<Geo::Point2<ScalarType> > &vert,
//                                         const std::vector<std::vector<int> > &face,
//                                         const Geo::Point2<ScalarType> &Ori,
//                                         const Geo::Point2<ScalarType> &Dir,
//                                         std::vector<int> &baryFIdx,
//                                         std::vector<Geo::Point2<ScalarType> > &bary)
//        {
//            baryFIdx.clear();
//            bary.clear();

//            //            Geo::Box2<ScalarType> BB_2D;
//            //            BB_2D.Init(vert);
//            //            if (edgeSize<=0)
//            //                edgeSize=BB_2D.Diag()*0.05;
//            //            BB_2D.Inflate(edgeSize);

//            //            Geo::TriangleSpatialIndex2D<ScalarType> Grid_2D;
//            //            Grid_2D.Init(vert,face);

//            //            int DivX=ceil(BB_2D.DimX()/edgeSize);
//            //            int DivY=ceil(BB_2D.DimY()/edgeSize);
//            //            DivX=std::max(DivX,1);
//            //            DivY=std::max(DivY,1);

//        }


//    static void MapCrossFieldToUV(const std::vector<Geo::Point3<ScalarType> > &VertNormal3D,
//                                  const std::vector<Geo::Point2<ScalarType> > &vert,
//                                  const std::vector<Field::CrossF<ScalarType> > &VertField)
//    {

//    }

public:

    static bool MapBetween(const std::vector<Geo::Point2<ScalarType> > &vert0,
                           const std::vector<Geo::Point2<ScalarType> > &vert1,
                           const std::vector<std::vector<int> > &face1,
                           std::vector<int> &faceIDx,std::vector<Geo::Point3<ScalarType>> &bary)
    {
        Geo::TriangleSpatialIndex2D<ScalarType> Grid_2D;
        Grid_2D.Init(vert1,face1);
        faceIDx.resize(vert0.size(),-1);
        bary.resize(vert0.size(),Geo::Point3<ScalarType>(0,0,0));
        for (size_t i=0;i<vert0.size();i++)
        {
            bool found=Grid_2D.FirstTrianglesOverPoint(vert1,face1,vert0[i],faceIDx[i],bary[i]);
            if (!found)return false;
        }
        return true;
    }


};

}

#endif
