#ifndef BASIC_DECOMPOSITION_CONDITION
#define BASIC_DECOMPOSITION_CONDITION

#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "patch_managing.h"
#include "catmull_clarkability.h"

namespace Geo{

template <class ScalarType>
class PatchCondition
{
    
public:

    //ScalarType elapsed_time;

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

    virtual bool Mandatory()const
    {return false;}

    virtual int Score(const PatchManaging<ScalarType> &PatchM,
                      const int &IndexPatch)const
    {return 0;}
    // PatchCondition()
    // {elapsed_time=0;}
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

    bool Mandatory()const override
    {return true;}
};

template <class ScalarType>
class DiskLikeCondition:public PatchCondition<ScalarType>
{

public:

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   
        //auto t0 = std::chrono::high_resolution_clock::now();
        if (PatchM.NumBorders(IndexPatch)!=1)return false;
        bool IsDisc= PatchM.isDiskLike(IndexPatch);
        // auto t1 = std::chrono::high_resolution_clock::now();
        // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        // //this->elapsed_time+=elapsed.count() * 1e-6;
        // std::cout<<"ELAPSED Disk Like:"<<elapsed.count()<<std::endl;
        return IsDisc;
    }

    std::string Name()const override
    {return std::string("disk_like");}

    bool Mandatory()const override
    {return true;}
};

template <class ScalarType>
class SingleSingularityCondition:public PatchCondition<ScalarType>
{

public:

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   
        //auto t0 = std::chrono::high_resolution_clock::now();
        bool HasMultiple= PatchM.HasMultipleSingularities(IndexPatch);
        // auto t1 = std::chrono::high_resolution_clock::now();
        // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        // //this->elapsed_time+=elapsed.count() * 1e-6;
        // std::cout<<"ELAPSED Multiple Sing:"<<elapsed.count()<<std::endl;
        return (!HasMultiple);
    }

    bool Mandatory()const override
    {return true;}

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
        // std::vector<int> Istances;
        // PosInstances<ScalarType>(PatchM.PData.SubPatchVertPos[IndexPatch],Istances);
        // for (size_t i=0;i<Istances.size();i++)
        // {
        //     assert(Istances[i]>0);
        //     if (Istances[i]>1)return false;
        // }
        // return true;
        //auto t0 = std::chrono::high_resolution_clock::now();
        bool HasDupl=HasDuplicateVert(PatchM.PData.SubPatchVertPos[IndexPatch]);
        // auto t1 = std::chrono::high_resolution_clock::now();
        // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        // //this->elapsed_time+=std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        // std::cout<<"ELAPSED Self Conn:"<<elapsed.count()<<std::endl;
        return (!HasDupl);
    }

    bool Mandatory()const override
    {return true;}
    
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
    {   
        //auto t0 = std::chrono::high_resolution_clock::now();
        bool MatchSing=PatchM.MatchSingularity(IndexPatch);
        // auto t1 = std::chrono::high_resolution_clock::now();
        // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        //this->elapsed_time+=elapsed.count() * 1e-6;
        //std::cout<<"ELAPSED Match Sing:"<<elapsed.count()<<std::endl;
        return (MatchSing);
    }

    std::string Name()const override
    {return std::string("match_sing");}

    //if has singularity and match then keep
    bool Mandatory()const override
    {
        return false;
    }

    //avoid to go to valence 2 and 6
    virtual int Score(const PatchManaging<ScalarType> &PatchM,
                     const int &IndexPatch)const override
    {
        int SideN=PatchM.NumSides(IndexPatch);
        if ((SideN==2)||(SideN==6))return -1;
        return 0;
    }
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

    bool Mandatory()const override
    {
        return true;
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
        // std::cout<<"##DEDE 0:"<<PatchM.NumBorders(IndexPatch)<<std::endl;
        // std::cout<<"##DEDE 1:"<<PatchM.NumSides(IndexPatch)<<std::endl;
        // std::cout<<"##DEDE 2:"<<PatchM.NumCorners(IndexPatch)<<std::endl;
        // std::cout<<"##DEDE 3:"<<PatchM.PData.SideLenghts[IndexPatch].size()<<std::endl;
        // PatchM.WriteMesh(IndexPatch,std::string("testCC.obj"));

        //auto t0 = std::chrono::high_resolution_clock::now();
        bool IsCC=IsCatmullClarkable(PatchM.PData.SideLenghts[IndexPatch][0],SideThr,true);
        // auto t1 = std::chrono::high_resolution_clock::now();
        // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        // this->elapsed_time+=elapsed.count() * 1e-6;
       
        return (IsCC);
        
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
        if (PatchM.PData.SideLenghts[IndexPatch].size()>1)return false;
        
        //patches made of 4 sides are ok
        if (PatchM.PData.SideLenghts[IndexPatch][0].size()==4)return true;

        ScalarType AvgL=0;
        ScalarType MinL=std::numeric_limits<ScalarType>::max();
        for (size_t i=0;i<PatchM.PData.SideLenghts[IndexPatch][0].size();i++)
        {   
            ScalarType CurrL=PatchM.PData.SideLenghts[IndexPatch][0][i];
            MinL=std::min(MinL,CurrL);
            AvgL+=CurrL;
        }

        AvgL/=PatchM.PData.SideLenghts[IndexPatch][0].size();

        ScalarType AbsThr=MinL*SideThr;
        for (size_t i=0;i<PatchM.PData.SideLenghts[IndexPatch][0].size();i++)
            if (fabs(PatchM.PData.SideLenghts[IndexPatch][0][i]-AvgL)>AbsThr)return false;

        return true;
    }

    std::string Name()const override
    {return std::string("equal_side");}

    EqualSideCondition(const ScalarType &_Thr)
    {
        SideThr=_Thr;
    }

    EqualSideCondition()
    {
        SideThr=1;
    }

    
};


}
#endif
