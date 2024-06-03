#ifndef BASIC_DECOMPOSITION_CONDITION
#define BASIC_DECOMPOSITION_CONDITION

#include <vector>
#include "patch_managing.h"
#include "catmull_clarkability.h"

namespace Geo{

template <class ScalarType>
class PatchCondition
{
public:

    virtual void UpdatePatchIndex(std::map<int,int> &PatchIdxRemap)
    {(void)PatchIdxRemap;}

//    virtual bool MatchPreconditions(const PatchManaging<ScalarType> &PatchM,
//                                    const int &IndexPatch)const
//    {(void)PatchM;(void)IndexPatch;return true;}

    virtual void UpdatePatchData(const PatchManaging<ScalarType> &PatchM,
                                 const int &IndexPatch)
    {(void)PatchM;(void)IndexPatch;}

    virtual bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                           const int &IndexPatch)const
    {(void)PatchM;(void)IndexPatch;return true;}

    virtual bool IsCorrectAfterGlobalData(const Geo::PatchManaging<ScalarType> &PatchM,
                                          const int &IndexPatch)const
    {(void)PatchM;(void)IndexPatch;return true;}

    virtual std::string Name()const
    {return std::string("basic");}

    virtual void GetPerFaceScalar(const PatchManaging<ScalarType> &PatchM,
                                  const int &IndexPatch,
                                  std::vector<ScalarType> &FaceV)const    
    {(void)PatchM;(void)IndexPatch;(void)FaceV;}

    virtual bool UpdateGlobalData(const Geo::PatchManaging<ScalarType> &PatchM)
    {(void)PatchM;return true;}

    virtual void SaveStatus()
    {}

    virtual void RestoreStatus()
    {}
};

template <class ScalarType>
class SingleBoundaryCondition:public PatchCondition<ScalarType>
{

public:

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   return (PatchM.NumBorders(IndexPatch)==1);}

    std::string Name()const override
    {return std::string("single_bound");}
};

template <class ScalarType>
class DiskLikeCondition:public PatchCondition<ScalarType>
{

public:

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   return PatchM.isDiskLike(IndexPatch);}

    std::string Name()const override
    {return std::string("disk_like");}
};

template <class ScalarType>
class SingleSingularityCondition:public PatchCondition<ScalarType>
{

public:

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   return (!PatchM.HasMultipleSingularities(IndexPatch));}

    std::string Name()const override
    {return std::string("single_sing");}
};

template <class ScalarType>
class SelfConnectingPatchCondition:public PatchCondition<ScalarType>
{

public:

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    { 
        std::vector<int> Istances;
        PosInstances<ScalarType>(PatchM.PData.SubPatchVertPos[IndexPatch],Istances);
        for (size_t i=0;i<Istances.size();i++)
        {
            assert(Istances[i]>0);
            if (Istances[i]>1)return false;
        }
        return true;
    }

    std::string Name()const override
    {return std::string("dual_meshing");}
};

template <class ScalarType>
class MathchSingularityCondition:public PatchCondition<ScalarType>
{

public:

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   return (PatchM.MatchSingularity(IndexPatch));}

    std::string Name()const override
    {return std::string("match_sing");}
};

template <class ScalarType>
class NumSidesCondition:public PatchCondition<ScalarType>
{
public:
    int minCornerNum;
    int maxCornerNum;

    bool IsValidSetup()const
    {
        //one should be valid
        if ((minCornerNum<=0)&&(maxCornerNum<=0))return false;

        if (maxCornerNum>0)
            if (maxCornerNum<5)return false;
        if (minCornerNum>0)
            if (minCornerNum>3)return false;

        if ((maxCornerNum>0)&&(minCornerNum>0))
            if (minCornerNum>=maxCornerNum)return false;
        return true;
    }

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {
        if (!IsValidSetup())
        {
            std::cout<<"WARNING CONDITIONS ON EDGE SIZE INCORRECT"<<std::endl;
            return true;
        }
        assert(IndexPatch<PatchM.NumPatches());
        if (PatchM.NumBorders(IndexPatch)>1)return false;
        if ((maxCornerNum>0)&&(PatchM.NumSides(IndexPatch)>maxCornerNum))return false;
        if ((minCornerNum>0)&&(PatchM.NumSides(IndexPatch)<minCornerNum))return false;
        return true;
    }

    std::string Name()const override
    {return std::string("num_sides");}

    NumSidesCondition()
    {
        minCornerNum=3;
        maxCornerNum=6;
    }
};

template <class ScalarType>
class CCabilityCondition:public PatchCondition<ScalarType>
{

public:
    ScalarType SideThr;

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {
        assert(IndexPatch<PatchM.PData.SideLenghts.size());

        //if there is more than one hole return false
        //if (PatchM.PData.SideLenghts[IndexPatch].size()>1)return false;
        if (PatchM.NumBorders(IndexPatch)!=1)return false;
        return (IsCatmullClarkable(PatchM.PData.SideLenghts[IndexPatch][0],SideThr,true));
    }

    std::string Name()const override
    {return std::string("cc_ability");}

    CCabilityCondition(const ScalarType &_Thr)
    {
        SideThr=_Thr;
    }
};

template <class ScalarType>
class EqualSideCondition:public PatchCondition<ScalarType>
{

public:
    ScalarType SideThr;

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {
        //if there is more than one hole return false
        if (PatchM.SideLenghts[IndexPatch].size()>1)return false;

        if (PatchM.SideLenghts[IndexPatch][0].size()==4)return true;

        ScalarType AvgL=0;
        for (size_t i=0;i<PatchM.SideLenghts[IndexPatch][0].size();i++)
            AvgL+=PatchM.SideLenghts[IndexPatch][0][i];

        AvgL/=PatchM.SideLenghts[IndexPatch][0].size();

        for (size_t i=0;i<PatchM.SideLenghts[IndexPatch][0].size();i++)
            if (fabs(PatchM.SideLenghts[IndexPatch][0][i]-AvgL)>SideThr)return false;

        return true;
    }

    std::string Name()const override
    {return std::string("equal_side");}

    EqualSideCondition(const ScalarType &_Thr)
    {
        SideThr=_Thr;
    }
};


}
#endif
