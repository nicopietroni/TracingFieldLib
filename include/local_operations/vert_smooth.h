#ifndef VERT_SMOOTH
#define VERT_SMOOTH

#include "../connectivity.h"
#include "../point3.h"
#include "../pos3.h"
#include "basic_local_operations.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <mesh_functions.h>
#include "../triangle3.h"
#include <vector>
//#include "../tangent_space_smooth.h"

namespace Geo {

template <class ScalarType> class VertSmooth {
  typedef Geo::Point3<ScalarType> CoordType;

  static CoordType AveragePos(const std::vector<Geo::Point3<ScalarType> > &Pos,
                              const std::vector<int> &VertStar) {
    CoordType AVGPos(0, 0, 0);
    for (size_t i = 0; i < VertStar.size(); i++) {
      int IndexV = VertStar[i];
      assert(IndexV >= 0);
      assert(IndexV < Pos.size());
      AVGPos += Pos[IndexV];
    }
    AVGPos /= VertStar.size();
    return AVGPos;
  }

    static void GetTriangles(const int &IndexV,
                             const std::vector<Geo::Point3<ScalarType> > &Pos,
                             const std::vector<std::vector<int> > &Faces,
                             const std::vector<std::vector<int> > &VF,
                             std::vector<Geo::Traingle3<ScalarType> > &TrisFan)
    {
        TrisFan.clear();
        for (size_t i=0;i<VF[IndexV].size();i++)
        {
            int CurrF=VF[IndexV][i];
            int IndexV0=Faces[CurrF][0];
            int IndexV1=Faces[CurrF][1];
            int IndexV2=Faces[CurrF][2];
            assert(IndexV0 >= 0);
            assert(IndexV0 < Pos.size());
            assert(IndexV1 >= 0);
            assert(IndexV1 < Pos.size());
            assert(IndexV2 >= 0);
            assert(IndexV2 < Pos.size());
            Geo::Point3<ScalarType>  P0=Pos[IndexV0];
            Geo::Point3<ScalarType>  P1=Pos[IndexV1];
            Geo::Point3<ScalarType>  P2=Pos[IndexV2];
            Geo::Traingle3<ScalarType> T(P0,P1,P2);
            TrisFan.push_back(T);
        }
    }

    static CoordType Reproject(const Geo::Point3<ScalarType> &Pos,
                               const std::vector<Geo::Traingle3<ScalarType> > &TrisFan)
    {
        ScalarType minD=std::numeric_limits<ScalarType>::max();
        Geo::Point3<ScalarType> minP=Pos;
        for (size_t i=0;i<TrisFan.size();i++)
        {
            Geo::Point3<ScalarType> TestP=ClosestPoint<ScalarType>(TrisFan[i],Pos);
            ScalarType TestD=(TestP-Pos).Norm();
            if (TestD>minD)continue;
            minP=TestP;
            minD=TestD;
        }
        return minP;
    }

public:

  static void GetInvolvedFaces(const Pos3 &SmoothP,
                               const std::vector<Geo::Point3<ScalarType> > &Pos,
                               const std::vector<std::vector<int> > &Faces,
                               const std::vector<std::vector<int> > &VF,
                               std::vector<int> &InvolvedF) {
    InvolvedF.clear();
    int IndexV = SmoothP.V;
    assert(IndexV >= 0);
    assert(IndexV < VF.size());
    InvolvedF = VF[IndexV];
  }

  static void Do(const Pos3 &SmoothP,
                 std::vector<Geo::Point3<ScalarType> > &Pos,
                 const std::vector<std::vector<int> > &Faces,
                 const std::vector<std::vector<int> > &VF, 
                 const std::vector<std::vector<int> > &VFIndex,
                 bool do_reproject = true) 
 {
    int IndexV = SmoothP.V;
    assert(IndexV >= 0);
    assert(IndexV < VF.size());
    std::vector<int> VertStar;
    
    Geo::getVVStar(Faces, VF, VFIndex,IndexV,VertStar);
    CoordType NewPos = AveragePos(Pos,VertStar);
     std::vector<Geo::Traingle3<ScalarType> > TriangleFan;
    GetTriangles(IndexV,Pos,Faces,VF,TriangleFan);
    if (do_reproject)
        Pos[IndexV]=Reproject(NewPos,TriangleFan);
    else
        Pos[IndexV]=NewPos;
  }
};

template <class ScalarType>
class VertSmoothOperation : public Geo::LocalOperation<ScalarType> {
  typedef Geo::Point3<ScalarType> CoordType;

public:
  bool Feasible() const override { return (true); }

  virtual void Do(bool UpdateAdj = true) override {
    VertSmooth<ScalarType>::GetInvolvedFaces(this->pos,this->Vert,
                                                       this->Faces,this->VF,
                                                       this->ModifiedF);
    VertSmooth<ScalarType>::Do(this->pos,
                              this->Vert,
                              this->Faces,
                              this->VF,
                              this->VFIndex); 
  }

  void GetInvolvedOLDFaces(std::vector<int> &involved_f) const override {
     VertSmooth<ScalarType>::GetInvolvedFaces(this->pos,this->Vert,
                                                       this->Faces,this->VF,
                                                       involved_f);
  }

  //    void GetInvolvedNEWFaces(std::vector<int> &involved_f)const override
  //    {
  // //       GetInvolvedNEWFaces(std::vector<int> &involved_f)const
  // // TriangleEdgeFlip<ScalarType>::GetInvolvedFaces(this->pos,this->Vert,
  // // this->Faces,this->NextF,
  // //                                                       this->NextE,
  // //                                                       involved_f);
  //    }
public:

  VertSmoothOperation(Pos3 _pos, std::vector<std::vector<int> > &_Faces,
                      std::vector<CoordType> &_Vert,
                      std::vector<bool> &_DeletedV,
                      std::vector<bool> &_DeletedF,
                      std::vector<std::vector<int> > &_VF,
                      std::vector<std::vector<int> > &_VFIndex,
                      std::vector<std::vector<int> > &_NextF,
                      std::vector<std::vector<int> > &_NextE,
                      std::vector<CoordType> &_VertNorm,
                      std::vector<CoordType> &_FaceNorm)
      : Geo::LocalOperation<ScalarType>(_pos, _Faces, _Vert, _DeletedV,
                                        _DeletedF, _VF, _VFIndex, _NextF,
                                        _NextE, _VertNorm, _FaceNorm) {}
};
} // namespace Geo

#endif
