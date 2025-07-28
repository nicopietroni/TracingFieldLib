#ifndef REMESH_SIZE
#define REMESH_SIZE

#include "IO/io_obj.h"
#include "geodesic.h"
#include "point3.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include "Space_Query/point_grid_3D.h"

namespace Geo {

template <class ScalarType> class RemeshSize {

public:
  static ScalarType
  RemeshSizeNormal(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const ScalarType &MinSize, const ScalarType &MaxSize,
                   ScalarType MaxNormVar = 135) {
    assert(MaxSize > MinSize);

    std::vector<Geo::Point3<ScalarType>> FaceNorm, VertNorm;
    Geo::ComputeNormals<ScalarType>(VertPos, Faces, FaceNorm, VertNorm);

    // typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
    // GParam.MaxD = MaxSize;

    PointSpatialIndex<ScalarType> PIndex;
    PIndex.Init(VertPos);

    std::vector<ScalarType> RemSize =
        std::vector<ScalarType>(VertPos.size(), MinSize);

    for (size_t i = 0; i < VertPos.size(); i++) {
      int CurrV=i;
      // GParam.Seeds.clear();
      // GParam.Seeds.push_back(i);
      // GParam.InInterval.clear();
      Geo::TMark TempMark;

      Geo::Point3<ScalarType> CurrN = VertNorm[CurrV];
      // Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, Faces, TempMark, GParam);

      // // sort wrt distance
      std::vector<int> InInterval;
      std::vector<std::pair<ScalarType, int>> DistKernel;
      // std::vector<int> InInterval(GParam.InInterval.begin(),
      //                             GParam.InInterval.end());
      PIndex.PointsInSphere(VertPos,VertPos[CurrV],MaxSize,InInterval);
      for (size_t j = 0; j < InInterval.size(); j++) {
        int IndexV = InInterval[j];
        // ScalarType DistV = GParam.Dist[IndexV];
        ScalarType DistV = (VertPos[IndexV] - VertPos[CurrV]).Norm();
        DistKernel.push_back(std::pair<ScalarType, int>(DistV, IndexV));
      }

      std::sort(DistKernel.begin(), DistKernel.end());
      ScalarType CurrSize = MaxSize;
      for (size_t j = 0; j < DistKernel.size(); j++) {
        int IndexV = DistKernel[j].second;
        Geo::Point3<ScalarType> TestN = VertNorm[IndexV];

        if (fabs(Geo::AngleDeg(TestN, CurrN)) < MaxNormVar)
          continue;

        CurrSize = DistKernel[j].first;
        break;
      }
      RemSize[CurrV] = std::max(RemSize[CurrV], CurrSize);
    }
    std::sort(RemSize.begin(), RemSize.end());
    int lessIdx=0;
    for (size_t i = 0; i < RemSize.size(); i++) {
      if (RemSize[i]>=MaxSize)continue;
      lessIdx++;
    }
    // std::cout<<"Less are "<<lessIdx<<std::endl;
    // std::cout<<"Total "<<RemSize.size()<<std::endl;
    int IndexPerc = floor(0.5 + RemSize.size() * 0.05);
    IndexPerc = std::min(IndexPerc, (int)RemSize.size() - 1);
    IndexPerc = std::max(IndexPerc, 0);
    
    return RemSize[IndexPerc];
     //return RemSize[0];
  }
};
} // namespace Geo

#endif
