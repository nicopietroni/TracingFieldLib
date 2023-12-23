#ifndef UV_DECOMPOSITION_CONDITION
#define UV_DECOMPOSITION_CONDITION

#include <vector>
#include "patch_managing.h"
#include "basic_decomposition_conditions.h"
#include "UV_param.h"
#include "UV_distortion.h"
#include "math_utils.h"

namespace Geo{

template <class ScalarType>
class UVCondition:public Geo::PatchCondition<ScalarType>
{
    std::vector<std::vector<Point2<ScalarType> > > PerPatchUVs;
    std::vector<std::vector<ScalarType > > PerPatchUVDists;

public:
    enum UVType{UVArap,UVLSqConformal};
    ScalarType MaxDist;
    UVType UVParamType;
    ScalarType Percentile;
    typename Geo::UVDistortion<ScalarType>::DistType DType;

    void UpdatePatchData(const PatchManaging<ScalarType> &PatchM,
                         const int &IndexPatch) override
    {
        assert(IndexPatch<PatchM.NumPatches());

        //allocate if needed
        if (IndexPatch>=PerPatchUVs.size())
        {
            PerPatchUVs.resize(IndexPatch+1);
            PerPatchUVDists.resize(IndexPatch+1);
        }

        //can't be parametrized in that case
        if (!PatchM.isDiskLike(IndexPatch))return;

        assert(PerPatchUVDists.size()==PerPatchUVs.size());
        //at least one face
        assert(PatchM.SubPatchVertPos[IndexPatch].size()>=3);
        assert(PatchM.SubPatchFaces[IndexPatch].size()>=1);

        if (UVParamType==UVArap)
            Geo::UVParam<ScalarType>::Arap(PatchM.SubPatchVertPos[IndexPatch],
                                           PatchM.SubPatchFaces[IndexPatch],
                                           PerPatchUVs[IndexPatch],true, 5);

        if (UVParamType==UVLSqConformal)
            Geo::UVParam<ScalarType>::LeastSquaresConformal(PatchM.SubPatchVertPos[IndexPatch],
                                                            PatchM.SubPatchFaces[IndexPatch],
                                                            PerPatchUVs[IndexPatch]);

        std::vector<ScalarType> FaceDistorion;
        Geo::UVDistortion<ScalarType>::ComputeDistortion(PatchM.SubPatchVertPos[IndexPatch],
                                                         PerPatchUVs[IndexPatch],
                                                         PatchM.SubPatchFaces[IndexPatch],
                                                         DType,PerPatchUVDists[IndexPatch]);
    }

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {
        assert(IndexPatch<PatchM.NumPatches());
        //can't be parametrized in that case
        //std::cout<<"test UV"<<std::endl;
        if (!PatchM.isDiskLike(IndexPatch))return false;

        assert(IndexPatch<PerPatchUVs.size());
        assert(IndexPatch<PerPatchUVDists.size());
        assert(PerPatchUVs[IndexPatch].size()==PatchM.SubPatchVertPos[IndexPatch].size());



        //then compute the max (considering the percentile)
        ScalarType testD=MaxPercentile(PerPatchUVDists[IndexPatch],Percentile);
        std::cout<<"test distortion:"<<testD<<std::endl;
        if (testD>MaxDist)return false;

        return true;
    }

    std::string Name()const override
    {return std::string("uv_distortion");}

    virtual void GetPerFaceScalar(const PatchManaging<ScalarType> &PatchM,
                                  const int &IndexPatch,
                                  std::vector<ScalarType> &FaceV)const override
    {
        //        assert(IndexPatch<PatchM.NumPatches());
        //        assert(IndexPatch<PerPatchUVs.size());
        //        assert(IndexPatch<PerPatchUVDists.size());
        if (IndexPatch>=PerPatchUVDists.size())
        {
            std::cout<<"Index Patch:"<<IndexPatch<<std::endl;
            std::cout<<"Size All Patches:"<<PerPatchUVDists.size()<<std::endl;
            int NumFPatch=PatchM.SubPatchFaces[IndexPatch].size();
            FaceV=std::vector<ScalarType>(NumFPatch,MaxDist);
            return;
        }
        if (PerPatchUVDists[IndexPatch].size()==0)
        {
            std::cout<<"Size All Data:"<<PerPatchUVDists[IndexPatch].size()<<std::endl;
            int NumFPatch=PatchM.SubPatchFaces[IndexPatch].size();
            FaceV=std::vector<ScalarType>(NumFPatch,MaxDist);
            return;
        }

        FaceV=PerPatchUVDists[IndexPatch];

        //        int NumFPatch=PatchM.SubPatchFaces[IndexPatch].size();
        //        FaceV=std::vector<ScalarType>(NumFPatch,MaxDist);

    }

    UVCondition(const UVType &_UVParamType,
                const typename Geo::UVDistortion<ScalarType>::DistType &_DType,
                const ScalarType &_MaxDist)
    {
        DType=_DType;
        UVParamType=_UVParamType;
        MaxDist=_MaxDist;
        Percentile=0.02;
    }
};



}
#endif
