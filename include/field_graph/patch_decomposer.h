#ifndef PATCH_DECOMPOSER
#define PATCH_DECOMPOSER

#include <vector>
//#include "../point3.h"
//#include "graph.h"
#include "patch_managing.h"
#include "path_sampling.h"
#include "basic_decomposition_conditions.h"
//#include "mesh_functions.h"

namespace Geo{

template <class ScalarType>
class PatchDecomposer
{
    PatchManaging<ScalarType> &PatchMan;
    Geo::CrossPathSampling<ScalarType> &PathSampl;
    std::vector<Geo::PatchCondition<ScalarType>* > PConditions;
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
            assert(IndexPatch<PatchMan.NumPatches());
            for (int i=0;i<PatchMan.SubPatchVertToOriginal[IndexPatch].size();i++)
            {
                int IndexV=PatchMan.SubPatchVertToOriginal[IndexPatch][i];
                assert(IndexV>=0);
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

    bool add_only_needed;
    bool allow_darts;
    bool allow_concave;

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

    void UpdateCondition(const int &IndexPatch)
    {
        for (size_t i=0;i<PConditions.size();i++)
        {
            PConditions[i]->UpdatePatchData(PatchMan,IndexPatch);
        }
    }

    bool SatisfyConditions(const int &IndexPatch)
    {
        assert(!PatchMan.isEmpty(IndexPatch));

        bool has_satisfied=true;

        //check if has concave or darts
        if (!allow_darts)
        {
            std::vector<int> SpikeCorn;
            PatchMan.CornersOfType(IndexPatch,PatchManaging<ScalarType>::CTypeSpike,SpikeCorn);
            if (SpikeCorn.size()>0)return false;
        }
        if (!allow_concave)
        {
            std::vector<int> ConcaveCorn;
            PatchMan.CornersOfType(IndexPatch,PatchManaging<ScalarType>::CTypeConcave,ConcaveCorn);
            if (ConcaveCorn.size()>0)return false;
        }

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
        //remove empty selected paths
        PathSampl.CleanEmptySelectedPaths();
        PatchMan.CompactEmptyPatches();

        //rremove the space around singularities
        PathSampl.InitSingularityNodes(false);

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

    bool GetConditionPerFaceScalar(std::vector<ScalarType> &FaceV,
                                   const std::string &name_cond)
    {
        for (size_t i=0;i<PConditions.size();i++)
        {
            if (PConditions[i]->Name()!=name_cond)continue;
            FaceV.clear();
            FaceV=std::vector<ScalarType>(PatchMan.NumFaces(),0);
            //std::cout<<"Test 0 "<<std::endl;
            for (size_t pi=0;pi<PatchMan.NumPatches();pi++)
            {
                if (PatchMan.isEmpty(pi))continue;
                std::vector<ScalarType> FaceVPatch;
                PConditions[i]->GetPerFaceScalar(PatchMan,pi,FaceVPatch);
                if (FaceVPatch.size()==0)continue;
//                std::cout<<"Size 0 "<<FaceVPatch.size()<<std::endl;
//                std::cout<<"Size 1 "<<PatchMan.SubPatchFacesToOriginal[i].size()<<std::endl;
//                std::cout<<"Size 2 "<<PatchMan.SubPatchFaces[i].size()<<std::endl;
                assert(FaceVPatch.size()==PatchMan.SubPatchFacesToOriginal[pi].size());
                assert(FaceVPatch.size()==PatchMan.SubPatchFaces[pi].size());
                for (size_t j=0;j<FaceVPatch.size();j++)
                {
                    int IndexGLobF=PatchMan.MapFaceToOriginal(pi,j);
                    FaceV[IndexGLobF]=FaceVPatch[j];
                }
            }
            return true;
        }
        return false;
    }
    //    void AddNextIfNeeded(bool &all_solved,
    //                         bool &added_sample)
    //    {
    //        all_solved=false;
    //        added_sample=false;
    //        if (add_only_needed)
    //            PathSampl.check_pass_on[0]=true;

    //        if (UnsolvedPatches.size()==0)
    //        {
    //            all_solved=true;
    //            added_sample=false;
    //            return;
    //        }

    //        added_sample=PathSampl.AddNextPath();
    //        if (!added_sample){return;}

    //        //verify local conditions
    //        std::vector<int> UpdatedPatches;
    //        PatchMan.Update(UpdatedPatches);

    //        for (size_t i=0;i<UpdatedPatches.size();i++)
    //        {
    //            int TestPatch=UpdatedPatches[i];

    //            if (PatchMan.isEmpty(TestPatch))continue;

    //            if (SatisfyConditions(TestPatch))
    //                UnsolvedPatches.erase(TestPatch);
    //            else
    //                UnsolvedPatches.insert(TestPatch);
    //        }

    //        RemoveEmptyPatches();

    //        if (add_only_needed)
    //            MarkUnsolvedPatches();
    //        all_solved=(UnsolvedPatches.size()==0);
    //    }

    void AddNext(bool &all_solved,
                 bool &added_sample)
    {
        all_solved=false;
        added_sample=false;
        if (add_only_needed)
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
            assert(TestPatch<PatchMan.NumPatches());
            if (PatchMan.isEmpty(TestPatch))continue;

            UpdateCondition(TestPatch);

            if (SatisfyConditions(TestPatch))
                UnsolvedPatches.erase(TestPatch);
            else
                UnsolvedPatches.insert(TestPatch);
        }

        RemoveEmptyPatches();

        if (add_only_needed)
            MarkUnsolvedPatches();
        all_solved=(UnsolvedPatches.size()==0);
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

            UpdateCondition(TestPatch);

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
            for (size_t i=0;i<UpdatedPatches.size();i++)
                UpdateCondition(UpdatedPatches[i]);

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

    PatchDecomposer(PatchManaging<ScalarType> &_PatchMan,
                    Geo::CrossPathSampling<ScalarType> &_PathSampl):
        PatchMan(_PatchMan),PathSampl(_PathSampl)
    {
        add_only_needed=true;
        allow_darts=false;
        allow_concave=false;
    }

    ~PatchDecomposer()
    {
        for (size_t i=0;i<PConditions.size();i++)
            delete(PConditions[i]);
    }
};

}
#endif
