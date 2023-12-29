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
    bool stop_when_solved;
    bool add_only_where_needed;
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

    void CompactData()
    {
        PathSampl.CleanEmptySelectedPaths();
        std::map<int,int> PatchIdxRemap;
        PatchMan.CompactEmptyPatches(PatchIdxRemap);

        //remap unsolved one
        std::set<int>::iterator IteUnsolved;
        std::set<int> UnsolvedPatchesNew;
        for(IteUnsolved=UnsolvedPatches.begin();IteUnsolved!=UnsolvedPatches.end();IteUnsolved++)
        {
            int IndexPatch=(*IteUnsolved);
            assert(PatchIdxRemap.count(IndexPatch)>0);
            int newIndex=PatchIdxRemap[(*IteUnsolved)];
            UnsolvedPatchesNew.insert(newIndex);
        }
        UnsolvedPatches=UnsolvedPatchesNew;
        for (size_t i=0;i<PConditions.size();i++)
            (*PConditions[i]).UpdatePatchIndex(PatchIdxRemap);
    }

    void RetraceOnUnsolved()
    {
        //remove empty selected paths
//        PathSampl.CleanEmptySelectedPaths();
//        std::map<int,int> PatchIdxRemap;
//        PatchMan.CompactEmptyPatches(PatchIdxRemap);

//        //remap unsolved one
//        std::set<int>::iterator IteUnsolved;
//        std::set<int> UnsolvedPatchesNew;
//        for(IteUnsolved=UnsolvedPatches.begin();IteUnsolved!=UnsolvedPatches.end();IteUnsolved++)
//        {
//            int IndexPatch=(*IteUnsolved);
//            assert(PatchIdxRemap.count(IndexPatch)>0);
//            int newIndex=PatchIdxRemap[(*IteUnsolved)];
//            UnsolvedPatchesNew.insert(newIndex);
//        }
//        UnsolvedPatches=UnsolvedPatchesNew;
//        for (size_t i=0;i<PConditions.size();i++)
//            (*PConditions[i]).UpdatePatchIndex(PatchIdxRemap);

        CompactData();

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
        if (add_only_where_needed)
            PathSampl.check_pass_on[0]=true;

        if ((UnsolvedPatches.size()==0)&&(stop_when_solved))
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

        if (add_only_where_needed)
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
        add_only_where_needed=true;
        stop_when_solved=true;
        allow_darts=false;
        allow_concave=false;
    }

    bool SavePatchesLayout(const std::string meshFile)
    {
        CompactData();
        std::ofstream fout(meshFile.c_str());
        if(!fout) return false;

        //first set the number of patches
        fout<<"#p"<<PatchMan.EdgePatchSides.size()<<"\n";

        //for each patch
        for (size_t pi=0;pi<PatchMan.EdgePatchSides.size();pi++)
        {
            //write patch index
            fout<<"p "<<pi<<"\n";

            //for each loops
            fout<<"#l "<<PatchMan.EdgePatchSides[pi].size()<<"\n";
            assert(PatchMan.EdgePatchSides[pi].size()==PatchMan.PatchCorners[pi].size());

            //for each loop
            for (size_t li=0;li<PatchMan.PatchCorners[pi].size();li++)
            {
                fout<<"l "<<li<<"\n";
                //number of sides equal nummber of corners
                assert(PatchMan.PatchCorners[pi][li].size()==PatchMan.EdgePatchSides[pi][li].size());

                //write corners
                fout<<"#c "<<PatchMan.PatchCorners[pi][li].size()<<"\n";

                for (size_t ci=0;ci<PatchMan.PatchCorners[pi][li].size();ci++)
                {
                    int IndexC=PatchMan.PatchCorners[pi][li][ci];
                    IndexC=PatchMan.MapToOriginal(pi,IndexC);
                    assert(IndexC>=0);
                    assert(IndexC<PatchMan.VertPos.size());
                    fout<<"c "<<IndexC<<"\n";
                }

                //write number of sides
                fout<<"#s "<<PatchMan.EdgePatchSides[pi][li].size()<<"\n";
                for (size_t si=0;si<PatchMan.EdgePatchSides[pi][li].size();si++)
                {
                    fout<<"s "<<si<<"\n";
                    int sizeL=PatchMan.EdgePatchSides[pi][li][si].size();
                    fout<<"#e "<<sizeL<<"\n";

                    //for each edge
                    for (size_t ei=0;ei<PatchMan.EdgePatchSides[pi][li][si].size();ei++)
                    {
                        int IndexV0=PatchMan.EdgePatchSides[pi][li][si][ei].first;
                        int IndexV1=PatchMan.EdgePatchSides[pi][li][si][ei].second;
                        IndexV0=PatchMan.MapToOriginal(pi,IndexV0);
                        IndexV1=PatchMan.MapToOriginal(pi,IndexV1);
                        assert(IndexV0>=0);
                        assert(IndexV0<PatchMan.VertPos.size());
                        assert(IndexV1>=0);
                        assert(IndexV1<PatchMan.VertPos.size());
                        fout<<"e "<<IndexV0<<" "<<IndexV1<<"\n";
                    }
                }
            }
        }
        fout.close();
        return true;
    }


    void SaveSubMeshes(const std::string PathSave)
    {
        CompactData();

        std::vector<Geo::Point3<ScalarType> > PosVect;
        std::vector<std::vector<int> >  FaceVect;
        assert(PatchMan.SubPatchVertPos.size()==PatchMan.SubPatchFaces.size());

        for (size_t i=0;i<PatchMan.SubPatchVertPos.size();i++)
        {
            int baseIdx=PosVect.size();
            PosVect.insert(PosVect.end(),
                           PatchMan.SubPatchVertPos[i].begin(),
                           PatchMan.SubPatchVertPos[i].end());

            for (size_t j=0;j<PatchMan.SubPatchFaces[i].size();j++)
            {
                FaceVect.push_back(std::vector<int>(3,-1));
                assert(PatchMan.SubPatchFaces[i][j].size()==3);
                for (size_t k=0;k<3;k++)
                {
                    int IndexV=PatchMan.SubPatchFaces[i][j][k];
                    FaceVect.back()[k]=IndexV+baseIdx;
                }
            }
        }
        WriteOBJ(PathSave.c_str(),PosVect,FaceVect);
    }

    ~PatchDecomposer()
    {
        for (size_t i=0;i<PConditions.size();i++)
            delete(PConditions[i]);
    }
};

}
#endif
