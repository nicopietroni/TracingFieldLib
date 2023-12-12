#ifndef PATCH_DECOMPOSER
#define PATCH_DECOMPOSER

#include <vector>
//#include "../point3.h"
//#include "graph.h"
#include "patch_managing.h"
#include "path_sampling.h"

namespace Geo{

template <class ScalarType>
class PatchCondition
{
public:
    virtual bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const
    {(void)PatchM;return true;}

    virtual bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                           const int &IndexPatch)const
    {(void)PatchM;(void)IndexPatch;return true;}
};

template <class ScalarType>
class DiskLikeCondition:public PatchCondition<ScalarType>
{

public:

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   return PatchM.isDiskLike(IndexPatch);}

};

template <class ScalarType>
class SingleSingularityCondition:public PatchCondition<ScalarType>
{

public:

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   return (!PatchM.HasMultipleSingularities(IndexPatch));}

};

template <class ScalarType>
class MathchSingularityCondition:public PatchCondition<ScalarType>
{

public:

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {   return (PatchM.MatchSingularity(IndexPatch));}

};

template <class ScalarType>
class PatchDecomposer
{
    PatchManaging<ScalarType> &PatchMan;
    Geo::CrossPathSampling<ScalarType> &PathSampl;
    std::vector<Geo::PatchCondition<ScalarType>* > PConditions;
    bool only_needed;


    std::set<int> UnsolvedPatches;

    void RemoveEmptyPatches()
    {
        std::vector<int> UnsolvedPatchesVect(UnsolvedPatches.begin(),UnsolvedPatches.end());
        for (size_t i=0;i<UnsolvedPatchesVect.size();i++)
        {
            int IndexP=UnsolvedPatchesVect[i];
            if (PatchMan.isEmpty(IndexP))
                UnsolvedPatches.erase(IndexP);
        }
    }

    void MarkUnsolvedPatches()
    {
        PathSampl.UnMarkAll(0);
        std::set<int>::iterator Ite;
        std::vector<int> to_mark_vert;
        for (Ite=UnsolvedPatches.begin();Ite!=UnsolvedPatches.end();Ite++)
        {
            int IndexPatch=(*Ite);
            for (int i=0;i<PatchMan.SubPatchVertToOriginal[IndexPatch].size();i++)
            {
                int IndexV=PatchMan.SubPatchVertToOriginal[IndexPatch][i];
                to_mark_vert.push_back(IndexV);
            }
        }
        PathSampl.MarkVertices(to_mark_vert,0);
    }

public:

    void ClearConditions()
    {
        PConditions.clear();
    }


    void AddCondition(Geo::PatchCondition<ScalarType> *PCond)
    {
        PConditions.push_back(PCond);
    }

    bool SatisfyConditions()
    {
        bool has_satisfied=true;
        for (size_t i=0;i<PConditions.size();i++)
        {
            has_satisfied&=PConditions[i]->IsCorrect(PatchMan);
            if (!has_satisfied)return false;
        }
        return true;
    }

    bool SatisfyConditions(const int &IndexPatch)
    {
        assert(!PatchMan.isEmpty(IndexPatch));

        bool has_satisfied=true;
        for (size_t i=0;i<PConditions.size();i++)
        {
            has_satisfied&=PConditions[i]->IsCorrect(PatchMan,IndexPatch);
            if (!has_satisfied)return false;
        }
        return true;
    }

    void Init()
    {
        for (size_t i=0;i<PatchMan.NumPatches();i++)
            if (!SatisfyConditions(i))UnsolvedPatches.insert(i);
    }

    void AddNextIfNeeded(bool &all_solved,
                         bool &added_sample)
    {
        all_solved=false;
        added_sample=false;
        if (only_needed)
           PathSampl.check_pass_on[0]=true;

        if (UnsolvedPatches.size()==0)
        {
            all_solved=true;
            added_sample=false;
            return;
        }

        added_sample=PathSampl.AddNextPath();
        if (!added_sample){return;}

        //verify local conditions
        std::vector<int> UpdatedPatches;
        PatchMan.Update(UpdatedPatches);

        for (size_t i=0;i<UpdatedPatches.size();i++)
        {
            int TestPatch=UpdatedPatches[i];

            if (PatchMan.isEmpty(TestPatch))continue;

            if (SatisfyConditions(TestPatch))
                UnsolvedPatches.erase(TestPatch);
            else
                UnsolvedPatches.insert(TestPatch);
        }

        RemoveEmptyPatches();

        if (only_needed)
            MarkUnsolvedPatches();
        all_solved=(UnsolvedPatches.size()==0);
    }

    bool RemoveIfPossible(int &Path_To_remove)
    {
        assert(UnsolvedPatches.size()==0);
        assert(Path_To_remove>=0);
        assert(Path_To_remove<PathSampl.SelectedPaths.size());
        std::vector<int> TestPath=PathSampl.SelectedPaths[Path_To_remove];
        if (PathSampl.SelectedPaths[Path_To_remove].size()==0)return false;
        PathSampl.SelectedPaths[Path_To_remove].clear();

        std::vector<int> UpdatedPatches;
        PatchMan.Update(UpdatedPatches);

        bool is_removable=true;
        std::cout<<"To_update:"<<UpdatedPatches.size()<<std::endl;
        for (size_t i=0;i<UpdatedPatches.size();i++)
        {
            int TestPatch=UpdatedPatches[i];

            if (PatchMan.isEmpty(TestPatch))continue;

            if (SatisfyConditions(TestPatch))
                continue;
            else
            {
                is_removable=false;
                break;
            }
        }
        if (!is_removable)
        {
            PathSampl.SelectedPaths[Path_To_remove]=TestPath;
            std::vector<int> UpdatedPatches;
            PatchMan.Update(UpdatedPatches);
            return false;
        }

        return true;
    }

    bool RemoveNextIfPossible()
    {
        bool has_removed=false;
        int Path_To_remove=PathSampl.SelectedPaths.size()-1;
        do{
            has_removed|=RemoveIfPossible(Path_To_remove);
            if (has_removed)return true;
            Path_To_remove--;
        }while ((!has_removed)&&(Path_To_remove>=0));
    }

//    void SatisfiedConditions(std::vector<bool> &Conditions)
//    {
//        bool has_satisfied=true;
//        for (size_t i=0;i<PConditions.size();i++)
//        {
//            has_satisfied&=PConditions[i].PConditions(PatchM,IndexPatch);
//            if (!has_satisfied)return false;
//        }
//        return true;
//    }

    PatchDecomposer(PatchManaging<ScalarType> &_PatchMan,
                    Geo::CrossPathSampling<ScalarType> &_PathSampl):
                    PatchMan(_PatchMan),PathSampl(_PathSampl)
    {
        only_needed=true;
    }

    ~PatchDecomposer()
    {
        for (size_t i=0;i<PConditions.size();i++)
           delete(PConditions[i]);
    }
};

}
#endif
