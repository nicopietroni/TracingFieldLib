#ifndef PATCH_DECOMPOSER
#define PATCH_DECOMPOSER

#include <vector>
//#include "../point3.h"
//#include "graph.h"
#include "patch_managing.h"
#include "path_sampling.h"
#include "mesh_functions.h"

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

    void GetUnsolvedMesh(std::vector<Geo::Point3<ScalarType> > &CopiedMeshesVert,
                         std::vector<std::vector<int> > &CopiedMeshesElem,
                         std::vector<int> &OriginalVert,
                         std::vector<int> &OriginalElem)
    {
        std::vector<int> ToCopyPatch(UnsolvedPatches.begin(),UnsolvedPatches.end());
        PatchMan.getMeshOfPatches(ToCopyPatch,CopiedMeshesVert,CopiedMeshesElem,
                                  OriginalVert,OriginalElem);
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

    void UpdateAllUnsolvedPatches()
    {
        for (size_t i=0;i<PatchMan.NumPatches();i++)
        {
            if (PatchMan.isEmpty(i))continue;

            if (!SatisfyConditions(i))UnsolvedPatches.insert(i);

        }
    }

    void GetBorderUnsolved(std::vector<int> &BorderUnsolved)
    {
        BorderUnsolved.clear();

        //get the non solved subpatches
        std::vector<Geo::Point3<ScalarType> > CopiedMeshesVert;
        std::vector<std::vector<int> > CopiedMeshesElem;
        std::vector<int> OriginalVert,OriginalElem;

        GetUnsolvedMesh(CopiedMeshesVert,CopiedMeshesElem,
                        OriginalVert,OriginalElem);

        //get BORDER VERTICES
        std::vector<std::vector<int> > NextF;
        std::vector<std::vector<int> > NextE;
        getFFAdjacency(CopiedMeshesElem,NextF,NextE);
        std::vector<size_t> BorderV;
        getBorderVert(CopiedMeshesElem,NextF,BorderV);

        for (size_t i=0;i<BorderV.size();i++)
        {
            int IndexLoc=BorderV[i];
            assert(IndexLoc<OriginalVert.size());
            assert(IndexLoc>=0);
            BorderUnsolved.push_back(OriginalVert[IndexLoc]);
        }
    }

    void RetraceOnUnsolved()
    {
        MarkUnsolvedPatches();

        PathSampl.sample_only_on[0]=true;
        PathSampl.limit_path_to[1]=true;

        //get the non solved subpatches
        std::vector<int> BorderUnsolved;
        GetBorderUnsolved(BorderUnsolved);

        //mark the ones on border
        PathSampl.UnMarkAll(1);
        PathSampl.MarkVertices(BorderUnsolved,1);
        PathSampl.SampleCandidates(false);

        PathSampl.limit_path_to[0]=false;
        PathSampl.sample_only_on[1]=false;
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

    //    bool RemoveIfPossible(int &Path_To_remove)
    //    {
    //        assert(UnsolvedPatches.size()==0);
    //        assert(Path_To_remove>=0);
    //        assert(Path_To_remove<PathSampl.SelectedPaths.size());
    //        std::vector<int> TestPath=PathSampl.SelectedPaths[Path_To_remove];
    //        if (PathSampl.SelectedPaths[Path_To_remove].size()==0)return false;
    //        PathSampl.SelectedPaths[Path_To_remove].clear();

    //        std::vector<int> UpdatedPatches;
    //        PatchMan.Update(UpdatedPatches);

    //        bool is_removable=true;
    //        std::cout<<"To_update:"<<UpdatedPatches.size()<<std::endl;
    //        for (size_t i=0;i<UpdatedPatches.size();i++)
    //        {
    //            int TestPatch=UpdatedPatches[i];

    //            if (PatchMan.isEmpty(TestPatch))continue;

    //            if (SatisfyConditions(TestPatch))
    //                continue;
    //            else
    //            {
    //                is_removable=false;
    //                break;
    //            }
    //        }
    //        if (!is_removable)
    //        {
    //            PathSampl.SelectedPaths[Path_To_remove]=TestPath;
    //            std::vector<int> UpdatedPatches;
    //            PatchMan.Update(UpdatedPatches);
    //            return false;
    //        }

    //        return true;
    //    }



    void SampleCandidatesUnsolvedPatches()
    {
        PathSampl.sample_only_on[0]=true;
        PathSampl.SampleCandidates();
        //PathSampl.SampleCandidates();
        PathSampl.sample_only_on[0]=false;
    }


    bool RemoveIfPossible(int &Path_To_remove)
    {

        assert(Path_To_remove>=0);
        assert(Path_To_remove<PathSampl.SelectedPaths.size());

        if (UnsolvedPatches.size()>0)
        {
            //check if passing from unsolved patch then can't do it
            std::set<int> PatchIdx;
            PatchMan.AddTouchedPatchesToSet(PathSampl.SelectedPaths[Path_To_remove],PatchIdx);
            std::set<int> intersectionIdx;
            std::set_intersection(PatchIdx.begin(), PatchIdx.end(),
                                  UnsolvedPatches.begin(), UnsolvedPatches.end(),
                                  std::inserter(intersectionIdx, intersectionIdx.begin()));

            // Check if the intersection is empty
            if (!intersectionIdx.empty())return false;
        }

        std::vector<int> TestPath=PathSampl.SelectedPaths[Path_To_remove];

        if (PathSampl.SelectedPaths[Path_To_remove].size()==0)return false;
        PathSampl.SelectedPaths[Path_To_remove].clear();

        std::vector<int> UpdatedPatches;
        PatchMan.Update(UpdatedPatches);

        bool is_removable=true;
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

    void GetColorBySolvedPatches(std::vector<Geo::Point3<ScalarType> > &FaceCol)
    {
        FaceCol.clear();
        FaceCol.resize(PatchMan.NumFaces(),Geo::Point3<ScalarType>(0,1,0));
        std::set<int>::iterator IteUns;
        for (IteUns=UnsolvedPatches.begin();IteUns!=UnsolvedPatches.end();IteUns++)
        {
            int IndexP=(*IteUns);
            for (size_t i=0;i<PatchMan.SubPatchFacesToOriginal[IndexP].size();i++)
            {
                int IndexF=PatchMan.SubPatchFacesToOriginal[IndexP][i];
                FaceCol[IndexF]=Geo::Point3<ScalarType>(1,0,0);
            }
        }
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
        return false;
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
