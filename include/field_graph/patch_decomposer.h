#ifndef PATCH_DECOMPOSER
#define PATCH_DECOMPOSER

#include "graph.h"
#include "graph_cross_field.h"
#include "patch_basic.h"
#include <cassert>
#include <iostream>
#include <ostream>
#include <set>
#include <sys/_types/_size_t.h>
#define SMOOTH_KERNEL 2
#include "../edge_mesh_functions.h"
#include "../surface_cluster/face_cluster.h"
#include "basic_decomposition_conditions.h"
#include "connectivity.h"
#include "patch_managing.h"
#include "path_sampling.h"
#include "segment2.h"
#include "tangent_space_smooth.h"
#include <vector>

namespace Geo {

template <class ScalarType> class PatchDecomposer {
  PatchManaging<ScalarType> &PatchMan;
  Geo::CrossPathSampling<ScalarType> &PathSampl;
  std::vector<Geo::PatchCondition<ScalarType> *> PConditions;
  std::set<int> UnsolvedPatches;
  std::set<int> VertSplitForFeatures;
  
  std::vector<std::vector<int>> SmoothKernel;
  std::vector<std::vector<int>> VFAdj;
  std::vector<Field::CrossF<ScalarType>> &VertCurv;
  std::vector<Field::CrossF<ScalarType>> &FaceCurv;

  void RemoveEmptyPatches() {
    std::vector<int> UnsolvedPatchesVect(UnsolvedPatches.begin(),
                                         UnsolvedPatches.end());
    for (size_t i = 0; i < UnsolvedPatchesVect.size(); i++) {
      int IndexP = UnsolvedPatchesVect[i];
      if (PatchMan.isEmpty(IndexP))
        UnsolvedPatches.erase(IndexP);
    }
  }

  // void AssUnsolvedPatchesToSampleOn(bool use_unsolved)//
  // {
  // //bool must_pass, bool sample_on) {
  //   // if (must_pass)
  //   //   std::fill(PathSampl.must_pass_on.begin(),
  //   PathSampl.must_pass_on.end(),
  //   //             true);
  //     // std::fill(PathSampl.must_pass_on.begin(),
  //     PathSampl.must_pass_on.end(),
  //     //           false);

  //   // if (sample_on)
  //   //   std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(),
  //   false);

  //   std::set<int>::iterator Ite;
  //   std::vector<int> to_mark_vert;
  //   std::set<int> PatchToConsider = UnsolvedPatches;

  //   // if (UnsolvedMandatoryP.size() > 0)
  //   //   PatchToConsider = UnsolvedMandatoryP;

  //   PathSampl.must_cross_e.clear();

  //   for (Ite = PatchToConsider.begin(); Ite != PatchToConsider.end(); Ite++)
  //   {
  //     int IndexPatch = (*Ite);
  //     assert(IndexPatch < PatchMan.NumPatches());
  //     for (int i = 0;i <
  //     PatchMan.PData.SubPatchVertToOriginal[IndexPatch].size(); i++) {
  //       int IndexV = PatchMan.PData.SubPatchVertToOriginal[IndexPatch][i];
  //       assert(IndexV >= 0);
  //       //assert(IndexV < PathSampl.must_pass_on.size());
  //       assert(IndexV < PathSampl.sample_on.size());
  //       // if (must_pass)
  //       //   PathSampl.must_pass_on[IndexV] = true;
  //       if (sample_on)
  //         PathSampl.sample_on[IndexV] = true;
  //     }

  //     for (int i = 0;
  //          i < PatchMan.PData.SubPatchFacesToOriginal[IndexPatch].size();
  //          i++) {
  //       int IndexVF = PatchMan.PData.SubPatchFacesToOriginal[IndexPatch][i];
  //       assert(IndexVF >= 0);
  //       assert(IndexVF < PatchMan.Faces.size());
  //       int IndexV0 = PatchMan.Faces[IndexVF][0];
  //       int IndexV1 = PatchMan.Faces[IndexVF][1];
  //       int IndexV2 = PatchMan.Faces[IndexVF][2];
  //       std::pair<int, int> KeyE0(std::min(IndexV0, IndexV1),
  //                                 std::max(IndexV0, IndexV1));
  //       std::pair<int, int> KeyE1(std::min(IndexV1, IndexV2),
  //                                 std::max(IndexV1, IndexV2));
  //       std::pair<int, int> KeyE2(std::min(IndexV2, IndexV0),
  //                                 std::max(IndexV2, IndexV0));
  //       PathSampl.must_cross_e.insert(KeyE0);
  //       PathSampl.must_cross_e.insert(KeyE1);
  //       PathSampl.must_cross_e.insert(KeyE2);
  //     }
  //   }

  //   // PathSampl.MarkVertices(to_mark_vert, 0);
  // }
public:
  void InitMustPassEdges() {
    PathSampl.must_cross_e.clear();

    for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
      // for each patch
      if (UnsolvedPatches.count(i) == 0)
        continue;

      for (int j = 0; j < PatchMan.PData.SubPatchFacesToOriginal[i].size();
           j++) {
        int IndexVF = PatchMan.PData.SubPatchFacesToOriginal[i][j];
        assert(IndexVF >= 0);
        assert(IndexVF < PatchMan.Faces.size());
        int IndexV0 = PatchMan.Faces[IndexVF][0];
        int IndexV1 = PatchMan.Faces[IndexVF][1];
        int IndexV2 = PatchMan.Faces[IndexVF][2];
        std::pair<int, int> KeyE0(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
        std::pair<int, int> KeyE1(std::min(IndexV1, IndexV2),
                                  std::max(IndexV1, IndexV2));
        std::pair<int, int> KeyE2(std::min(IndexV2, IndexV0),
                                  std::max(IndexV2, IndexV0));
        PathSampl.must_cross_e.insert(KeyE0);
        PathSampl.must_cross_e.insert(KeyE1);
        PathSampl.must_cross_e.insert(KeyE2);
      }
    }
  }

  void InitSampling(bool only_unsolved) //
  {
    std::set<int>::iterator Ite;
    std::vector<int> to_mark_vert;
    std::set<int> PatchToConsider = UnsolvedPatches;

    std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(), false);
    std::fill(PathSampl.force_sample.begin(), PathSampl.force_sample.end(),
              false);

    for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
      // for each patch check if is
      if (only_unsolved && (UnsolvedPatches.count(i) == 0))
        continue;

      // set the sampling space
      for (int j = 0; j < PatchMan.PData.SubPatchVertToOriginal[i].size();
           j++) {
        int IndexV = PatchMan.PData.SubPatchVertToOriginal[i][j];
        assert(IndexV >= 0);
        assert(IndexV < PathSampl.sample_on.size());

        if (start_internal)
          PathSampl.sample_on[IndexV] = true;
      }

      if (start_concave) {
        std::vector<int> currSel;
        PatchMan.GetPatchGlobalVertOfType(i, CTypeConcave, currSel);
        for (size_t j = 0; j < currSel.size(); j++)
          PathSampl.force_sample[currSel[j]] = true;
      }

      if (start_spike) {
        std::vector<int> currSel;
        PatchMan.GetPatchGlobalVertOfType(i, CTypeSpike, currSel);
        for (size_t j = 0; j < currSel.size(); j++)
          PathSampl.force_sample[currSel[j]] = true;
      }
    }
  }

private:
  void GetUnsolvedMesh(std::vector<Geo::Point3<ScalarType>> &CopiedMeshesVert,
                       std::vector<std::vector<int>> &CopiedMeshesElem,
                       std::vector<int> &OriginalVert,
                       std::vector<int> &OriginalElem)
  //,
  // bool OnlyMandatory = false)
  {
    std::vector<int> ToCopyPatch;
    // if (!OnlyMandatory)
    ToCopyPatch =
        std::vector<int>(UnsolvedPatches.begin(), UnsolvedPatches.end());
    // else
    //   ToCopyPatch = std::vector<int>(UnsolvedMandatoryP.begin(),
    //                                  UnsolvedMandatoryP.end());

    PatchMan.getMeshOfPatches(ToCopyPatch, CopiedMeshesVert, CopiedMeshesElem,
                              OriginalVert, OriginalElem);
  }

  bool to_init_distance;

public:
  bool stop_when_solved;
  bool add_only_where_needed;
  bool allow_darts;
  bool allow_concave;
  bool dynamic_smooth;
  bool dynamic_updates;
  bool split_removal;

  bool start_internal;
  bool start_orginal_border;
  bool start_concave;
  bool start_spike;

  // bool end_internal;
  bool end_original_border;
  bool end_patch_border;
  bool end_concave;
  bool end_spike;

  std::set<int> FixedPaths;

  void ClearConditions() { PConditions.clear(); }

  void AddCondition(Geo::PatchCondition<ScalarType> *PCond) {
    PConditions.push_back(PCond);
  }

  void InitForceSamplingSpace() // bool and_sample) {
  {
    // force sampling on all points only for the case
    // of tracing from concave or spike
    // bool NoForce = !(start_concave || start_spike);
    // if (NoForce)
    //   std::fill(PathSampl.force_sample.begin(),
    //   PathSampl.force_sample.end(),false);

    std::vector<int> to_Mark;
    if (start_concave) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeConcave, currSel);
      to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    }
    if (start_spike) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeSpike, currSel);
      to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    }

    // initialize force sampling to false
    std::fill(PathSampl.force_sample.begin(), PathSampl.force_sample.end(),
              false);

    // then set to force sampling on selected categories
    for (size_t i = 0; i < to_Mark.size(); i++) {
      int IndexV = to_Mark[i];
      assert(IndexV >= 0);
      assert(IndexV < PatchMan.VertPos.size());
      assert(IndexV < PathSampl.force_sample.size());
      // if (!and_sample)
      PathSampl.force_sample[IndexV] = true;
      // else {
      // assert(IndexV < PathSampl.sample_on.size());
      // PathSampl.force_sample[IndexV] = PathSampl.sample_on[IndexV];
      //}
    }
  }

  // init the space where sampling the paths
  // the flag use_previous is used to check if
  // the sample set must AND (intersect) with previous one
  void InitSamplingSpace(bool use_previous) {
    bool allStart =
        (start_internal && start_orginal_border && start_concave && start_spike);

    if (allStart) {
      if (use_previous)
        return;
      std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(), true);
    }

    std::vector<int> to_Mark;
    if (start_internal) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeNone, currSel);
      to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    }
    if (start_orginal_border) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeFlat, currSel);
      //remove the one that have been inserted
      for (size_t i=0;i<currSel.size();i++)
      {
        if (VertSplitForFeatures.count(currSel[i]))continue;
        to_Mark.push_back(currSel[i]);
      }
    }
    if (start_concave) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeConcave, currSel);
      to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    }
    if (start_spike) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeSpike, currSel);
      to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    }

    if (!use_previous) {
      std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(), false);
      for (size_t i = 0; i < to_Mark.size(); i++) {
        int IndexV = to_Mark[i];
        assert(IndexV >= 0);
        assert(IndexV < PatchMan.VertPos.size());
        assert(IndexV < PathSampl.sample_on.size());
        PathSampl.sample_on[IndexV] = true;
      }
    } else {
      std::vector<bool> new_sample_on =
          std::vector<bool>(PathSampl.sample_on.size(), false);
      for (size_t i = 0; i < to_Mark.size(); i++) {
        int IndexV = to_Mark[i];
        assert(IndexV >= 0);
        assert(IndexV < new_sample_on.size());
        new_sample_on[IndexV] = true;
      }
      for (size_t i = 0; i < PathSampl.sample_on.size(); i++)
        PathSampl.sample_on[i] = PathSampl.sample_on[i] && new_sample_on[i];
    }
  }

  // this one
  void InitArrivalPoint(bool and_allow_pass) {
    bool allArrival = (end_original_border && end_patch_border && end_concave && end_spike);

    if (allArrival) {
      // if (and_allow_pass)
      //   PathSampl.allow_arrive_to = PathSampl.allow_path_on;
      // else
      std::fill(PathSampl.allow_arrive_to.begin(),
                PathSampl.allow_arrive_to.end(), true);
      return;
    }

    std::vector<int> to_Mark;
    if (end_concave) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeConcave, currSel);
      to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    }
    if (end_spike) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeSpike, currSel);
      to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    }
    if (end_original_border) {
      std::vector<int> currSel;
      PatchMan.GetAllGlobalVertOfType(CTypeFlat, currSel);
      for (size_t i=0;i<currSel.size();i++)
      {
        if (VertSplitForFeatures.count(currSel[i]))continue;
        to_Mark.push_back(currSel[i]);
      }
      //to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    }

    std::fill(PathSampl.allow_arrive_to.begin(),
              PathSampl.allow_arrive_to.end(), false);

    for (size_t i = 0; i < to_Mark.size(); i++) {
      int IndexV = to_Mark[i];
      assert(IndexV >= 0);
      assert(IndexV < PatchMan.VertPos.size());
      assert(IndexV < PathSampl.allow_arrive_to.size());
      // if (!and_allow_pass)
      PathSampl.allow_arrive_to[IndexV] = true;
      // else {
      //   assert(IndexV < PathSampl.allow_arrive_to.size());
      //   PathSampl.allow_arrive_to[IndexV] = PathSampl.allow_path_on[IndexV];
      // }
    }
  }

  void UpdateCondition(const int &IndexPatch) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < PConditions.size(); i++) {
      PConditions[i]->UpdatePatchData(PatchMan, IndexPatch);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed0 =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    TimeUpdateConditions += elapsed0.count() * 1e-6;
  }

  bool SatisfyConditions(const int &IndexPatch) {
    assert(!PatchMan.isEmpty(IndexPatch));

    bool has_satisfied = true;

    // check if has concave or darts
    if (!allow_darts) {
      std::vector<int> SpikeCorn;
      PatchMan.CornersOfType(IndexPatch, Geo::CTypeSpike, SpikeCorn);
      if (SpikeCorn.size() > 0)
        return false;
    }
    if (!allow_concave) {
      std::vector<int> ConcaveCorn;
      PatchMan.CornersOfType(IndexPatch, Geo::CTypeConcave, ConcaveCorn);
      if (ConcaveCorn.size() > 0)
        return false;
    }

    for (size_t i = 0; i < PConditions.size(); i++) {
      has_satisfied &= PConditions[i]->IsCorrect(PatchMan, IndexPatch);
      if (!has_satisfied)
        return false;
    }
    return true;
  }

  bool SatisfyConditionsAfterGlobalUpdate(const int &IndexPatch) {
    assert(!PatchMan.isEmpty(IndexPatch));

    bool has_satisfied = true;
    for (size_t i = 0; i < PConditions.size(); i++) {
      has_satisfied &=
          PConditions[i]->IsCorrectAfterGlobalData(PatchMan, IndexPatch);
      if (!has_satisfied)
        return false;
    }
    return true;
  }

  void Init() {
    PathSampl.Init();
    to_init_distance = true;
  }

  void UpdateAllUnsolvedPatches() {
    UnsolvedPatches.clear();
    for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
      std::cout << "Checking Patch " << i << std::endl;

      if (PatchMan.isEmpty(i))
        continue;
      UpdateCondition(i);
      if (!SatisfyConditions(i))
        UnsolvedPatches.insert(i);
      if (!SatisfyConditionsAfterGlobalUpdate(i))
        UnsolvedPatches.insert(i);
    }
  }

  void GetBorderUnsolved(std::vector<int> &BorderUnsolved) {
    BorderUnsolved.clear();

    // get the non solved subpatches
    std::vector<Geo::Point3<ScalarType>> CopiedMeshesVert;
    std::vector<std::vector<int>> CopiedMeshesElem;
    std::vector<int> OriginalVert, OriginalElem;

    GetUnsolvedMesh(CopiedMeshesVert, CopiedMeshesElem, OriginalVert,
                    OriginalElem); //, OnlyMandatory);

    // get BORDER VERTICES
    std::vector<std::vector<int>> NextF;
    std::vector<std::vector<int>> NextE;
    getFFAdjacency(CopiedMeshesElem, NextF, NextE);
    std::vector<size_t> BorderV;
    getBorderVert(CopiedMeshesElem, NextF, BorderV);

    for (size_t i = 0; i < BorderV.size(); i++) {
      int IndexLoc = BorderV[i];
      assert(IndexLoc < OriginalVert.size());
      assert(IndexLoc >= 0);
      BorderUnsolved.push_back(OriginalVert[IndexLoc]);
    }
  }

  void CompactData() {
    // // remap unsolved one
    std::set<int>::iterator IteUnsolved;
    std::set<int> UnsolvedPatchesNew;

    for (IteUnsolved = UnsolvedPatches.begin();
         IteUnsolved != UnsolvedPatches.end(); IteUnsolved++) {
      int IndexPatch = (*IteUnsolved);
      if (PatchMan.PData.PatchOriginalFaces[IndexPatch].size() == 0)
        continue;
      UnsolvedPatchesNew.insert(IndexPatch);
    }
    UnsolvedPatches = UnsolvedPatchesNew;

    PathSampl.CleanEmptySelectedPaths();
    std::map<int, int> PatchIdxRemap;
    PatchMan.CompactEmptyPatches(PatchIdxRemap);

    UnsolvedPatchesNew.clear();
    // std::set<int> UnsolvedMandatoryPatchesNew;

    for (IteUnsolved = UnsolvedPatches.begin();
         IteUnsolved != UnsolvedPatches.end(); IteUnsolved++) {
      int IndexPatch = (*IteUnsolved);
      assert(PatchIdxRemap.count(IndexPatch) > 0);
      int newIndex = PatchIdxRemap[(*IteUnsolved)];
      UnsolvedPatchesNew.insert(newIndex);
    }

    UnsolvedPatches = UnsolvedPatchesNew;
    // UnsolvedMandatoryP = UnsolvedMandatoryPatchesNew;

    for (size_t i = 0; i < PConditions.size(); i++)
      (*PConditions[i]).UpdatePatchIndex(PatchIdxRemap);
  }

  bool IsSolvingFeatures() {
    // if only start from particular ones
    if ((start_concave || start_spike) && (!start_internal) && (!start_orginal_border))
      return true;
    return false;
  }

  // void SetLengthPrioType()
  // {
  //   PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioLen;
  // }

  // void SetDistPrioType()
  // {
  //   PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioDist;
  // }

  void RetraceCandidates() {

    CompactData();

    // InitForceSamplingSpace();

    // std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(), false);

    // remove the space around singularities
    PathSampl.InitSingularityNodes(false);

    // mark the patches that are not solved
    // the first one says that useful path must pass on it
    // the second says that it sample on them only
    // if (add_only_where_needed)
    //  MarkUnsolvedPatches(true, true);
    // else {
    //   std::fill(PathSampl.must_pass_on.begin(), PathSampl.must_pass_on.end(),
    //             true);
    //   std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(),
    //   true);
    // }

    // allows paths to pass on the one that must solve only
    // PathSampl.allow_path_on = PathSampl.must_pass_on;

    // AssUnsolvedPatchesToSampleOn();
    // InitSamplingSpace(true);
    //  if (force_sample==NULL)
    //    InitForceSamplingSpace(false);
    //  else
    //    PathSampl.force_sample=(*force_sample);

    // bool doing_mandatory = (UnsolvedMandatoryP.size() > 0);
    // if (!doing_mandatory) {
    InitArrivalPoint(true);

    /// this one stops to previously found borders, useful for
    /// recursive call
    if (end_patch_border && add_only_where_needed) {
      // get the non solved subpatches
      std::vector<int> BorderUnsolved;
      GetBorderUnsolved(BorderUnsolved); //, false);
      // set artificially to exit so it stops to borders
      PathSampl.SetBorderCanPassToExit(BorderUnsolved);
    }

    // if ((start_concave || start_spike)&&(!start_internal)&&(!start_border)
    if (IsSolvingFeatures())
      PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioLen;
    else
      PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioDist;

    auto t0 = std::chrono::high_resolution_clock::now();
    // sample the candidate
    // PathSampl.Graph().EnableAllArcs();
    // PathSampl.Graph().EnableAllNodes();
    // PathSampl.Graph().PrintData();

    // std::vector<int> ForcedSamplVert;
    // PathSampl.GetForcedSampleVert(ForcedSamplVert);
    // std::cout<<"FORCING "<<ForcedSamplVert.size()<<std::endl;
    // if (!start_concave)
    // {
    //   std::cout<<"TEST NAVACCHIO 0"<<std::endl;
    //   PathSampl.sample_on=std::vector<bool>(PathSampl.sample_on.size(),
    //   false);
    //   PathSampl.force_sample=std::vector<bool>(PathSampl.force_sample.size(),
    //   false);
    // }

    PathSampl.SampleCandidates(to_init_distance);
    // if (!start_concave)
    // {
    //   exit(0);
    //   // std::cout<<"TEST NAVACCHIO 1"<<std::endl;
    //   //  PathSampl.sample_on=std::vector<bool>(PathSampl.sample_on.size(),
    //   false);
    // }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed1 =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    TimeTracing += elapsed1.count() * 1e-6;

    // if not initialized at this point will be
    to_init_distance = false;

    // restore the node type for borders that have been changed to exit
    // for the recursive call
    if (end_patch_border && add_only_where_needed) // || doing_mandatory) {
      PathSampl.RestoreNodeType();
    // }

    // auto t3 = std::chrono::high_resolution_clock::now();
    // auto elapsed2 =
    //     std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    // TimeAddUpdatePatch += elapsed2.count() * 1e-6;
  }

  bool GetConditionPerFaceScalar(std::vector<ScalarType> &FaceV,
                                 const std::string &name_cond) {
    for (size_t i = 0; i < PConditions.size(); i++) {
      if (PConditions[i]->Name() != name_cond)
        continue;
      FaceV.clear();
      FaceV = std::vector<ScalarType>(PatchMan.NumFaces(), 0);
      // std::cout<<"Test 0 "<<std::endl;
      for (size_t pi = 0; pi < PatchMan.NumPatches(); pi++) {
        if (PatchMan.isEmpty(pi))
          continue;
        std::vector<ScalarType> FaceVPatch;
        PConditions[i]->GetPerFaceScalar(PatchMan, pi, FaceVPatch);
        if (FaceVPatch.size() == 0)
          continue;

        assert(FaceVPatch.size() ==
               PatchMan.PData.SubPatchFacesToOriginal[pi].size());
        assert(FaceVPatch.size() == PatchMan.PData.SubPatchFaces[pi].size());
        for (size_t j = 0; j < FaceVPatch.size(); j++) {
          int IndexGLobF = PatchMan.MapFaceToOriginal(pi, j);
          FaceV[IndexGLobF] = FaceVPatch[j];
        }
      }
      return true;
    }
    return false;
  }

  void SmoothAround(const std::vector<int> &PathsIDX,
                    const std::set<int> InvolvedPatches, int smooth_steps = 5) {
    // initialize the kernel if needed
    if (SmoothKernel.size() != PatchMan.VertPos.size())
      Geo::FindVertexStar(PatchMan.Faces, PatchMan.VertPos.size(), SmoothKernel,
                          SMOOTH_KERNEL);

    std::vector<std::vector<int>> VFIndex;
    if (VFAdj.size() != PatchMan.VertPos.size())
      getVFAdjacency(PatchMan.Faces, VFAdj, VFIndex);

    // then find vertex paths
    std::vector<std::vector<int>> ModVertPaths;
    PathSampl.GetVertexSelectedPaths(PathsIDX, ModVertPaths);

    std::vector<int> ToSmoothVert;
    for (size_t i = 0; i < ModVertPaths.size(); i++)
      for (size_t j = 0; j < ModVertPaths[i].size(); j++) {
        int IndexV = ModVertPaths[i][j];
        assert(IndexV >= 0);
        assert(IndexV < SmoothKernel.size());
        std::vector<int> CurrKern = SmoothKernel[IndexV];
        ToSmoothVert.insert(ToSmoothVert.end(), CurrKern.begin(),
                            CurrKern.end());
      }

    std::sort(ToSmoothVert.begin(), ToSmoothVert.end());
    std::vector<int>::iterator last =
        std::unique(ToSmoothVert.begin(), ToSmoothVert.end());
    ToSmoothVert.erase(last, ToSmoothVert.end());

    // then filter the one belonging to other patches
    std::vector<int> ToSmoothVertSwap;
    for (size_t i = 0; i < ToSmoothVert.size(); i++) {
      int IndexV = ToSmoothVert[i];
      assert(IndexV >= 0);
      assert(IndexV < VFAdj.size());
      bool has_other_patch = false;
      for (size_t j = 0; j < VFAdj[IndexV].size(); j++) {
        int IndexF = VFAdj[IndexV][j];

        assert(IndexF >= 0);
        assert(IndexF < PatchMan.Faces.size());
        assert(IndexF < PatchMan.PData.OriginalFaceToPatch.size());

        int IndexP = PatchMan.PData.OriginalFaceToPatch[IndexF];

        assert(IndexP >= 0);
        assert(IndexP < PatchMan.NumPatches());

        if (InvolvedPatches.count(IndexP) == 0) {
          has_other_patch = true;
          break;
        }
      }
      if (has_other_patch)
        continue;
      ToSmoothVertSwap.push_back(IndexV);
    }
    ToSmoothVert = ToSmoothVertSwap;

    //        for (size_t i=0;i<InvolvedPatches.size();i++)
    //        {
    //            int IndexP=InvolvedPatches[i];
    //            assert(IndexP>=0);
    //            assert(IndexP<PatchMan.NumPatches());
    //        }

    // then add the kernel to only smooth ones
    typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
    UVP.SmoothOnly = ToSmoothVert;
    UVP.fixBorder = true;
    // get all paths
    std::vector<std::vector<int>> VertPaths;
    PathSampl.GetVertexSelectedPaths(VertPaths);
    for (size_t i = 0; i < VertPaths.size(); i++) {
      if (VertPaths[i].size() == 0)
        continue;
      for (int j = 0; j < (int)VertPaths[i].size() - 1; j++) {
        int IndexV0 = VertPaths[i][j];
        int IndexV1 = VertPaths[i][j + 1];
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        UVP.FeatureLines.push_back(key);
      }
    }

    UVP.Steps = smooth_steps;
    Geo::Local_Param_Smooth<ScalarType>::Smooth(PatchMan.VertPos,
                                                PatchMan.Faces, UVP);

    AdaptFieldToMesh(PatchMan.VertPos, PatchMan.Faces, VertCurv, FaceCurv);

    PatchMan.UpdateSubPatchPos();
    PatchMan.InitVisualCornersEdges();
  }

  void AddNext(bool &all_solved, bool &added_sample) {
    all_solved = false;
    added_sample = false;

    // if (add_only_where_needed)
    //     PathSampl.check_pass_on[0] = true;

    if ((UnsolvedPatches.size() == 0) && (stop_when_solved)) {
      all_solved = true;
      added_sample = false;
      std::cout << "All Solved!" << std::endl;
      return;
    }

    double timeReTrace, TimeUpdDist, TimeConflict;

    // check in case start on spike or concave that
    // if the vertex is no more concave or spike remove the candidate
    // if (start_spike && end_spike)

    // GetAllGlobalVertOfType(const CornerType CType,
    //                           std::vector<int> &SelectedVert)

    // see if already solved that features
    // in such a case no need to keep adding
    bool checkMustStartEnd = IsSolvingFeatures() && add_only_where_needed;
    if (checkMustStartEnd) {
      PathSampl.MustStartEnd.clear();
      std::vector<int> SelectedVertConcave, SelectedSpike;

      if (start_concave || end_concave)
        PatchMan.GetAllGlobalVertOfType(CTypeConcave, SelectedVertConcave);

      if (start_spike || end_spike)
        PatchMan.GetAllGlobalVertOfType(CTypeSpike, SelectedSpike);

      PathSampl.MustStartEnd = SelectedVertConcave;
      PathSampl.MustStartEnd.insert(PathSampl.MustStartEnd.end(),
                                    SelectedSpike.begin(), SelectedSpike.end());
    }

    added_sample =
        PathSampl.AddNextPath(add_only_where_needed, checkMustStartEnd,
                              timeReTrace, TimeUpdDist, TimeConflict);

    // auto t1 = std::chrono::high_resolution_clock::now();
    // auto elapsed0 =
    //     std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    TimeTracing += timeReTrace;
    TimeUpdateDistances += TimeUpdDist;
    TimeCheckConflicts += TimeConflict;

    if (!added_sample) {
      return;
    }

    // verify local conditions
    std::vector<int> ToAdd, ToRemove;
    // std::set<int> ToUpdPatches;
    PatchMan.FindAddRemovePaths(ToAdd, ToRemove);

    if (dynamic_smooth) {
      std::vector<int> PathsIDX = ToAdd;
      // PathsIDX.insert(PathsIDX.end(), ToRemove.begin(), ToRemove.end());

      std::set<int> TouchedPatches;
      PatchMan.GetTouchedPatchesFromPathIndex(PathsIDX, TouchedPatches);
      SmoothAround(PathsIDX, TouchedPatches);
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<int> UpdatedPatches;

    PatchMan.Update(ToAdd, ToRemove, UpdatedPatches);

    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed0 =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    TimeUpdateUpdatePatches += elapsed0.count() * 1e-6;

    for (size_t i = 0; i < UpdatedPatches.size(); i++) {
      int TestPatch = UpdatedPatches[i];
      assert(TestPatch < PatchMan.NumPatches());
      if (PatchMan.isEmpty(TestPatch))
        continue;

      UpdateCondition(TestPatch);
    }

    UnsolvedPatches.clear();
    for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
      if (PatchMan.isEmpty(i))
        continue;

      if (!SatisfyConditions(i))
        UnsolvedPatches.insert(i);
    }

    for (size_t i = 0; i < PConditions.size(); i++) {
      PConditions[i]->UpdateGlobalData(PatchMan);
    }

    for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
      if (PatchMan.isEmpty(i))
        continue;

      auto t4 = std::chrono::high_resolution_clock::now();
      if (!SatisfyConditionsAfterGlobalUpdate(i))
        UnsolvedPatches.insert(i);
      auto t5 = std::chrono::high_resolution_clock::now();
      auto elapsed2 =
          std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4);
      TimeUpdateConditions += elapsed2.count() * 1e-6;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed1 =
        std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    TimeUpdateConditions += elapsed1.count() * 1e-6;

    RemoveEmptyPatches();

    // if (add_only_where_needed)
    //   MarkUnsolvedPatches(true, false);

    if (stop_when_solved)
      all_solved = (UnsolvedPatches.size() == 0);

    // auto t2 = std::chrono::high_resolution_clock::now();
    // auto elapsed1 =
    //     std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    // TimeAddUpdatePatch += elapsed1.count() * 1e-6;
  }

  // int time_savSt;
  // int time_rem0;
  // int time_rem1;
  // int time_rem2;
  // int time_rem3;
  // int time_rem4;

  bool RemoveIfPossible(int &Path_To_remove, bool debug_msg = false) {
    if (debug_msg)
      std::cout << "Checking Path:" << Path_To_remove << std::endl;

    assert(Path_To_remove >= 0);
    assert(Path_To_remove < PathSampl.SelectedPaths.size());

    // check if already removed
    if (PathSampl.SelectedPaths[Path_To_remove].size() == 0)
      return false;

    if (UnsolvedPatches.size() > 0) {
      // check if passing from unsolved patch then can't do it
      std::set<int> PatchIdx;
      PatchMan.AddTouchedPatchesToSet(PathSampl.SelectedPaths[Path_To_remove],
                                      PatchIdx);
      std::set<int> intersectionIdx;
      std::set_intersection(
          PatchIdx.begin(), PatchIdx.end(), UnsolvedPatches.begin(),
          UnsolvedPatches.end(),
          std::inserter(intersectionIdx, intersectionIdx.begin()));

      // Check if the intersection is empty
      if (!intersectionIdx.empty())
        return false;
    }

    std::vector<int> TestPath = PathSampl.SelectedPaths[Path_To_remove];

    if (PathSampl.SelectedPaths[Path_To_remove].size() == 0)
      return false;

// SAVING STATUS
#ifdef SAVE_STATUS
    auto timesave_0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < PConditions.size(); i++)
      PConditions[i]->SaveStatus();
    PatchMan.SaveStatus();

    auto timesave_1 = std::chrono::high_resolution_clock::now();
    time_savSt += (std::chrono::duration_cast<std::chrono::milliseconds>(
                       timesave_1 - timesave_0))
                      .count();
// END SAVING STATUS
#endif

    PathSampl.SelectedPaths[Path_To_remove].clear();

    // Get the current time point before the first instruction
    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<int> UpdatedPatches;
    PatchMan.Update(UpdatedPatches);

    // std::cout<<"Updating:"<<UpdatedPatches.size()<<std::endl<<std::endl;
    bool is_removable = true;

    // Get the current time point before the first instruction
    auto time1 = std::chrono::high_resolution_clock::now();

    // LOCAL UPDATE
    for (size_t i = 0; i < UpdatedPatches.size(); i++) {
      int TestPatch = UpdatedPatches[i];

      if (PatchMan.isEmpty(TestPatch))
        continue;

      if (debug_msg)
        PatchMan.WriteStats(TestPatch);

      UpdateCondition(TestPatch);
    }

    // check correctness locally
    bool IsOkLocally = true;
    for (size_t i = 0; i < UpdatedPatches.size(); i++) {
      if (PatchMan.isEmpty(UpdatedPatches[i]))
        continue;

      IsOkLocally &= SatisfyConditions(UpdatedPatches[i]);
      if (!IsOkLocally)
        break;
    }

    // Get the current time point before the first instruction
    auto time2 = std::chrono::high_resolution_clock::now();

    // GLOBAL UPDATE
    if (IsOkLocally) {
      for (size_t i = 0; i < PConditions.size(); i++)
        PConditions[i]->UpdateGlobalData(PatchMan);
    }

    // Get the current time point before the first instruction
    auto time3 = std::chrono::high_resolution_clock::now();

    // CHECK GLOBAL CONDITIONS
    if (!IsOkLocally)
      is_removable = false;
    else {
      for (size_t i = 0; i < UpdatedPatches.size(); i++) {
        int TestPatch = UpdatedPatches[i];
        if (SatisfyConditionsAfterGlobalUpdate(TestPatch))
          continue;
        else {
          is_removable = false;
          break;
        }
      }
    }
    // Get the current time point before the first instruction
    auto time4 = std::chrono::high_resolution_clock::now();

    if (!is_removable) {
      PathSampl.SelectedPaths[Path_To_remove] = TestPath;

// // RESTORING STATUS
#ifdef SAVE_STATUS
      auto timesave_0 = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < PConditions.size(); i++)
        PConditions[i]->RestoreStatus();
      PatchMan.RestoreStatus();
      auto timesave_1 = std::chrono::high_resolution_clock::now();
      time_savSt += (std::chrono::duration_cast<std::chrono::milliseconds>(
                         timesave_1 - timesave_0))
                        .count();
// END RESTORING STATUS
#else
      auto timesave_0 = std::chrono::high_resolution_clock::now();
      std::vector<int> UpdatedPatches;
      PatchMan.Update(UpdatedPatches);
      //  std::cout<<"Updating
      //  Rollback:"<<UpdatedPatches.size()<<std::endl<<std::endl;
      for (size_t i = 0; i < UpdatedPatches.size(); i++)
        UpdateCondition(UpdatedPatches[i]);

      auto timesave_1 = std::chrono::high_resolution_clock::now();
      // time_savSt += (std::chrono::duration_cast<std::chrono::milliseconds>(
      //                    timesave_1 - timesave_0))
      //                   .count();
#endif
    }

    // Get the current time point before the first instruction
    // auto time5 = std::chrono::high_resolution_clock::now();

    // time_rem0 +=
    //     (std::chrono::duration_cast<std::chrono::milliseconds>(time1 -
    //     time0))
    //         .count();
    // time_rem1 +=
    //     (std::chrono::duration_cast<std::chrono::milliseconds>(time2 -
    //     time1))
    //         .count();
    // time_rem2 +=
    //     (std::chrono::duration_cast<std::chrono::milliseconds>(time3 -
    //     time2))
    //         .count();
    // time_rem3 +=
    //     (std::chrono::duration_cast<std::chrono::milliseconds>(time4 -
    //     time3))
    //         .count();
    // time_rem4 +=
    //     (std::chrono::duration_cast<std::chrono::milliseconds>(time5 -
    //     time4))
    //         .count();

    // Print the duration in microseconds
    // if (debug_msg) {

    //   std::cout << "Time Save/Restore St: " << time_rem0
    //             << " milli microseconds" << std::endl;
    //   std::cout << "Time 0: " << time_rem0 << " milli microseconds"
    //             << std::endl;
    //   std::cout << "Time 1: " << time_rem1 << " milli microseconds"
    //             << std::endl;
    //   std::cout << "Time 2: " << time_rem2 << " milli microseconds"
    //             << std::endl;
    //   std::cout << "Time 3: " << time_rem3 << " milli microseconds"
    //             << std::endl;
    //   std::cout << "Time 4: " << time_rem4 << " milli microseconds"
    //             << std::endl;
    // }
    if (is_removable) {
      // enable the nodes
      PathSampl.EnableGraphPathNodes(TestPath);
      PathSampl.EnableGraphArches(TestPath);
      PathSampl.RemoveToConflictNodes(TestPath);
      PathSampl.RemoveToConflictEdges(TestPath);
    }
    return is_removable;
  }

  void GetColorBySolvedPatches(std::vector<Geo::Point3<ScalarType>> &FaceCol) {
    FaceCol.clear();
    FaceCol.resize(PatchMan.NumFaces(), Geo::Point3<ScalarType>(0, 1, 0));
    std::set<int>::iterator IteUns;
    for (IteUns = UnsolvedPatches.begin(); IteUns != UnsolvedPatches.end();
         IteUns++) {
      int IndexP = (*IteUns);
      for (size_t i = 0;
           i < PatchMan.PData.SubPatchFacesToOriginal[IndexP].size(); i++) {
        int IndexF = PatchMan.PData.SubPatchFacesToOriginal[IndexP][i];
        FaceCol[IndexF] = Geo::Point3<ScalarType>(1, 0, 0);
      }
    }
  }

  int current_rem_cand;

  bool RemoveNextIfPossible() {
    if (current_rem_cand > PathSampl.SelectedPaths.size())
      current_rem_cand = PathSampl.SelectedPaths.size() - 1;

    bool has_removed = false;
    do {
      if (FixedPaths.count(current_rem_cand) > 0) {
        current_rem_cand--;
        continue;
      }
      has_removed |= RemoveIfPossible(current_rem_cand);

      // TestedRemovedPath.insert(current_rem_cand);
      current_rem_cand--;

      if (has_removed) {
        return (current_rem_cand >= 0);
      }

    } while ((!has_removed) && (current_rem_cand >= 0));
    return false;
  }

  // void SampleInitialCandidates() {

  //   InitSamplingSpace(false);
  //   InitForceSamplingSpace(); // false);
  //   InitArrivalPoint(false);

  //   // //NEW cannot pass from concave or spike can only depart or arrive
  //   // std::vector<int> ConcaveV;
  //   // std::vector<int> SpikeV;
  //   // PatchMan.GetAllGlobalVertOfType(CTypeConcave,ConcaveV);
  //   // PatchMan.GetAllGlobalVertOfType(CTypeSpike,SpikeV);
  //   // for (size_t i=0;i<ConcaveV.size();i++)
  //   //   PathSampl.allow_path_on[ConcaveV[i]]=false;
  //   //  for (size_t i=0;i<SpikeV.size();i++)
  //   //   PathSampl.allow_path_on[SpikeV[i]]=false;

  //   // MarkStartSamples(0,false);
  //   //  PathSampl.sample_only_on[0] = true;

  //   if (IsSolvingFeatures()) {
  //     PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioLen;
  //     std::cout << "Using Lenght" << std::endl;
  //   } else {
  //     PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioDist;
  //     std::cout << "Using Dist" << std::endl;
  //   }
  //   PathSampl.SampleCandidates(true);
  //   to_init_distance = false;
  //   // PathSampl.sample_only_on[0] = false;
  // }

  PatchDecomposer(PatchManaging<ScalarType> &_PatchMan,
                  Geo::CrossPathSampling<ScalarType> &_PathSampl,
                  std::vector<Field::CrossF<ScalarType>> &_VertCurv,
                  std::vector<Field::CrossF<ScalarType>> &_FaceCurv)
      : PatchMan(_PatchMan), PathSampl(_PathSampl), VertCurv(_VertCurv),
        FaceCurv(_FaceCurv) {
    add_only_where_needed = true;
    stop_when_solved = true;

    allow_darts = false;
    allow_concave = false;
    current_rem_cand = -1;
    dynamic_smooth = false;
    dynamic_updates = true;
    split_removal = false;
    // time_savSt = 0;
    // time_rem0 = 0;
    // time_rem1 = 0;
    // time_rem2 = 0;
    // time_rem3 = 0;
    // time_rem4 = 0;

    // StartTraceType=Geo::CTNone;
    // EndTraceType=Geo::CTNone;
    start_internal = true;
    start_orginal_border = true;
    start_concave = true;
    start_spike = true;

    // end_internal = true;
    end_original_border =true;
    end_patch_border = false;
    end_concave = true;
    end_spike = true;

    to_init_distance = true;
    first_ite = true;
  }

  //    std::set<int> TestedRemovedPath;
  //    std::map<int,std::vector<int> > TestedRemovedPathNeigh;

  void SetSpikeSpikeConfiguration() {
    start_internal = false;
    start_orginal_border = false;
    start_concave = false;
    start_spike = true;

    end_original_border=false;
    end_patch_border = false;
    end_concave = false;
    end_spike = true;
  }

  void SetFeatureFeatureConfiguration() {
    start_internal = false;
    start_orginal_border = false;
    start_concave = true;
    start_spike = true;

    end_original_border=false;
    end_patch_border = false;
    end_concave = true;
    end_spike = true;
  }

  void SetFeatureAllConfiguration() {
    start_internal = false;
    start_orginal_border = false;
    start_concave = true;
    start_spike = true;

    end_original_border=true;
    end_patch_border = true;
    end_concave = true;
    end_spike = true;
  }

  void SetBaseContinuumConfiguration() {
    start_internal = true;
    start_orginal_border = true;
    start_concave = false;
    start_spike = false;

    end_original_border=true;
    end_patch_border = false;
    end_concave = false;
    end_spike = false;
  }

  void SetBaseSplitConfiguration() {
    start_internal = true;
    start_orginal_border = true;
    start_concave = false;
    start_spike = false;

    end_original_border=true;
    end_patch_border = true;
    end_concave = false;
    end_spike = false;
  }

  bool AddAllPaths() {
    // auto t0 = std::chrono::high_resolution_clock::now();

    bool all_solved = false;
    bool added_sample = false;
    bool added_once = false;
    std::cout << "ADDING PATHS" << std::endl;

    do {
      AddNext(all_solved, added_sample);
      // update dart if needed
      UpdateOrthoDartNodeTypeIfNeeded();
      added_once |= added_sample;
    } while ((!all_solved) && (added_sample));

    CompactData();

    // auto t1 = std::chrono::high_resolution_clock::now();

    // auto elapsed0 =
    //     std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    // TimeAdd += elapsed0.count() * 1e-6;
    return (added_once);
  }

  void RemovePaths() {
    current_rem_cand = PathSampl.SelectedPaths.size() - 1;
    std::cout << "REMOVING PATHS" << std::endl;
    int NumPath0 = PathSampl.SelectedPaths.size();
    while (RemoveNextIfPossible()) {
      // update dart if needed
      // UpdateOrthoDartNodeTypeIfNeeded();
    };

    CompactData();

    int NumPath1 = PathSampl.SelectedPaths.size();
    // distance should be reinitialized with remaining
    if (NumPath0 != NumPath1)
      PathSampl.InitializeDistance();
  }

  void RemoveSplitted(bool mergeCollinear = true) {
    PathSampl.SplitPaths(FixedPaths);

    // TestedRemovedPath.clear();

    PatchMan.PData.OLDPathNodes = PathSampl.SelectedPaths;
    std::cout << "After Split: " << PathSampl.SelectedPaths.size() << " paths"
              << std::endl;

    std::cout << "REMOVING SPlitted PATHS" << std::endl;

    bool has_removed = false;
    do {
      std::cout << "REMOVING cycle" << std::endl;
      current_rem_cand = PathSampl.SelectedPaths.size() - 1;
      has_removed = false;
      while (RemoveNextIfPossible()) {
        has_removed = true;
      }
      std::cout << "Merging Collinear" << std::endl;
      if (mergeCollinear && has_removed)
        PathSampl.MergeCollinearPath();
      std::cout << "Merged" << std::endl;
    } while (has_removed);

    CompactData();
    // PathSampl.MergeCollinearPath();
  }

  bool SolveAllFeatures() // bool UpdateForce = true) {
  {
    // std::cout<<"STEP 0"<<std::endl;
    // PathSampl.Graph().PrintData();

    bool has_added = false;

    // SetSpikeSpikeConfiguration();

    // // then trace in the one unsolved
    // RetraceOnUnsolved();

    // has_added |=AddAllPaths();

    // std::vector<bool> force_sample_swap;
    PathSampl.force_direct = true;
    // first do feature to feature
    SetFeatureFeatureConfiguration();
    InitSampling(add_only_where_needed);
    if (add_only_where_needed)
      InitMustPassEdges();

    // if (UpdateForce)
    //  std::cout<<"STEP 0.1"<<std::endl;
    //  PathSampl.Graph().PrintData();

    // then trace in the one unsolved, only the first time then use borders
    // if (first_ite)
    std::cout << "**TRACING first feature step **" << std::endl;
    RetraceCandidates();

    // std::cout<<"STEP 1"<<std::endl;
    // PathSampl.Graph().PrintData();

    has_added |= AddAllPaths();

    // std::cout<<"STEP 2"<<std::endl;
    // PathSampl.Graph().PrintData();

    // set features to all possible arrival points
    SetFeatureAllConfiguration();
    InitSampling(add_only_where_needed);
    if (add_only_where_needed)
      InitMustPassEdges();
    // if (!end_border)
    // {
    // std::cout<<"STEP 3"<<std::endl;
    // PathSampl.Graph().PrintData();

    std::cout << "**TRACING second feature step **" << std::endl;
    RetraceCandidates();

    // std::cout<<"STEP 4"<<std::endl;
    // PathSampl.Graph().PrintData();
    has_added |= AddAllPaths();

    // std::cout<<"STEP 5"<<std::endl;
    // PathSampl.Graph().PrintData();
    //}
    PathSampl.force_direct = false;
    return has_added;
  }

  bool SolveAllMesh() {
    // UnsolvedPatches.clear();
    PathSampl.force_direct = false;
    bool has_added = false;

    if (first_ite) {
      SetBaseContinuumConfiguration();
      InitSampling(add_only_where_needed);
      if (add_only_where_needed)
        InitMustPassEdges();

      RetraceCandidates();
      has_added |= AddAllPaths();
      first_ite=false;
    }

    SetBaseSplitConfiguration();
    InitSampling(add_only_where_needed);
    if (add_only_where_needed)
      InitMustPassEdges();

    RetraceCandidates();
    has_added |= AddAllPaths();
    return has_added;
  }

  void Decompose(bool Do_remove) {

    AddAllPaths();

    if (Do_remove) {
      RemovePaths();
    }

    if (split_removal) {
      RemoveSplitted();
    }
  }

  double TimeTracing, TimeUpdateConditions, TimeUpdateUpdatePatches,
      TimeUpdateDistances, TimeCheckConflicts, TimeRemovePaths;

  bool ForceSolveSpikeConcave(int IndexP) // bool Do_features,
  {
    std::vector<int> ToStartV;
    std::vector<int> ConcaveV;

    if (PathSampl.avoidSingularities)
      PathSampl.EnableSingularities();

    PatchMan.GetPatchGlobalVertOfType(IndexP, Geo::CTypeConcave, ConcaveV);

    std::vector<int> SpikeV;
    PatchMan.GetPatchGlobalVertOfType(IndexP, Geo::CTypeSpike, SpikeV);

    ToStartV = ConcaveV;
    ToStartV.insert(ToStartV.end(), SpikeV.begin(), SpikeV.end());
    // std::cout << "Found " << ToStartV.size() << " Starting Vertices"
    //           << std::endl;

    // sample on vertices that are within this set
    PathSampl.sample_on = std::vector<bool>(PatchMan.VertPos.size(), false);

    // constrain to add the sample from this set
    PathSampl.force_sample = std::vector<bool>(PatchMan.VertPos.size(), false);
    for (size_t i = 0; i < ToStartV.size(); i++)
      PathSampl.force_sample[ToStartV[i]] = true;

    // limit the path to pass only on those ones
    // PathSampl.allow_path_on = std::vector<bool>(PatchMan.VertPos.size(),
    // false);

    // int sizeV = PatchMan.PData.SubPatchVertToOriginal[IndexP].size();
    // for (size_t i = 0; i < sizeV; i++) {
    //   int CurrV = PatchMan.PData.SubPatchVertToOriginal[IndexP][i];
    //   PathSampl.allow_path_on[CurrV] = true;
    // }

    // // limit the path to arrive only on those ones
    // PathSampl.must_pass_on = PathSampl.allow_path_on;

    // get all boundary vertices and set as possible arrive
    std::vector<int> BorderV;
    getBorderVert(PatchMan.PData.SubPatchFaces[IndexP], BorderV);
    PatchMan.MapToOriginal(IndexP, BorderV);
    PathSampl.allow_arrive_to =
        std::vector<bool>(PatchMan.VertPos.size(), false);
    for (size_t i = 0; i < BorderV.size(); i++)
      PathSampl.allow_arrive_to[BorderV[i]] = true;

    // then set the arrival as exit (except the one that are forced to sample )
    std::vector<std::pair<int, NType>> SwapNodeType;

    // set as exit node
    std::vector<int> SetAsExit;
    // std::cout << "BorderV " << BorderV.size() << std::endl;
    for (size_t i = 0; i < BorderV.size(); i++) {
      if (PathSampl.force_sample[BorderV[i]])
        continue;
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(BorderV[i], IndexN);
      for (size_t j = 0; j < IndexN.size(); j++)
        SetAsExit.push_back(IndexN[j]);
    }
    // std::cout << "Setting as exit " << SetAsExit.size() << std::endl;
    PathSampl.SetNodesToExit(SetAsExit);

    PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioLen;
    PathSampl.only_visible = false;
    // sample the candidate
    PathSampl.SampleCandidates(false);
    if (PathSampl.avoidSingularities)
      PathSampl.DisableSingularities();
    PathSampl.only_visible = true;
    bool has_added = AddAllPaths();
    PathSampl.RestoreNodeType();
    return has_added;
  }

  void ForceSolveSpikeConcave() {
    std::cout << "*** DOING FORCE SPIKE/CONCAVE PASS ***" << std::endl;
    bool has_added_path = false;
    do {
      has_added_path = false;
      for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
        if (PatchMan.isEmpty(i))
          continue;
        std::vector<int> SpikeCorn, ConcaveCorn;
        PatchMan.CornersOfType(i, Geo::CTypeSpike, SpikeCorn);
        PatchMan.CornersOfType(i, Geo::CTypeConcave, ConcaveCorn);
        if ((SpikeCorn.size() == 0) && (ConcaveCorn.size() == 0))
          continue;

        has_added_path |= ForceSolveSpikeConcave(i);
      }
    } while (has_added_path);
  }

  bool ForceSolveTopology(int IndexP) {
    std::vector<Geo::Point3<ScalarType>> PatchVertPos =
        PatchMan.PData.SubPatchVertPos[IndexP];
    std::vector<std::vector<int>> PatchFaces =
        PatchMan.PData.SubPatchFaces[IndexP];

    // split into patch with Voronoi
    std::vector<std::vector<int>> ClusterFaces;
    Geo::FindVoronoiClusters<ScalarType>(PatchVertPos, PatchFaces,
                                         ClusterFaces);

    // then convert border into paths
    // then add the path
    return true;
  }

  void ForceSolveTopology() {
    std::cout << "*** DOING FORCE TOPOLOGY PASS ***" << std::endl;
    bool has_added_path = false;
    do {
      has_added_path = false;
      for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
        if (PatchMan.isEmpty(i))
          continue;
        std::vector<int> SpikeCorn, ConcaveCorn;
        int NumB = PatchMan.NumBorders(i);
        int GenusP = PatchMan.Genus(i);
        bool must_spl = (NumB != 1) || ((NumB == 1) && (GenusP != 0));
        if (must_spl)
          has_added_path |= ForceSolveTopology(i);
      }
    } while (has_added_path);
  }

private:
  bool first_ite;

  void DecomposeCycle(bool Do_remove) {
    PathSampl.InitWithSelected();
    bool has_added = false;
    first_ite = true;
    do {
      // this firs step solve all the features (concave, spikes)
      has_added = false;
      std::cout << "* solving all features *" << std::endl;
      bool has_added_f = SolveAllFeatures();

      if (has_added_f)
        std::cout << "Added Features" << std::endl;

      has_added |= has_added_f;

      auto t0 = std::chrono::high_resolution_clock::now();
      // remove unnecessary
      if (has_added_f && Do_remove)
        RemovePaths();

      auto t1 = std::chrono::high_resolution_clock::now();
      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeRemovePaths += elapsed.count() * 1e-6;

      // add the remaining to solve all the rest of the objetive function
      bool has_added_m = SolveAllMesh();
      // if (has_added_m)
      //   std::cout << "Added Mesh" << std::endl;
      has_added |= has_added_m;

      // then remove unnecessary
      t0 = std::chrono::high_resolution_clock::now();
      if (has_added_m && Do_remove)
        RemovePaths();
      t1 = std::chrono::high_resolution_clock::now();
      elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeRemovePaths += elapsed.count() * 1e-6;

      // first_ite = false;
    } while (has_added);
  }

  void UpdateGraphSingleLineSpike(
      Geo::TracingGraph<ScalarType> &TGraph,
      const std::vector<Field::CrossF<ScalarType>> &VertCurv,
      const std::vector<std::pair<int, int>> &Features) {
    // update the spike vertices with only one neigh
    std::map<int, std::vector<int>> AdjList;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      AdjList[IndexV0].push_back(IndexV1);
      AdjList[IndexV1].push_back(IndexV0);
    }

    std::vector<int> SelectedVert;
    PatchMan.GetAllGlobalVertOfType(Geo::CTypeSpike, SelectedVert);

    for (size_t i = 0; i < SelectedVert.size(); i++) {
      int IndexV = SelectedVert[i];
      // opened spike
      if (AdjList[IndexV].size() == 0)
        continue;
      // assert(AdjList[IndexV].size() > 0);
      // std::cout<<"size "<<AdjList[IndexV].size()<<std::endl;
      // std::cout<<"TestN "<<IndexV<<std::endl;
      // // std::cout<<"N0 "<<AdjList[IndexV][0]<<std::endl;
      // // std::cout<<"N1 "<<AdjList[IndexV][1]<<std::endl;
      // exit(0);
      if (AdjList[IndexV].size() != 1)
        continue;
      // std::cout<<"DEEDEDE"<<std::endl;
      int IndexV0 = AdjList[IndexV][0];
      Geo::Point3<ScalarType> DirSpike =
          PatchMan.VertPos[IndexV] - PatchMan.VertPos[IndexV0];
      DirSpike.Normalize();
      // get the direction
      assert(IndexV < VertCurv.size());
      int DirIdx = VertCurv[IndexV].BestDirectionI(DirSpike);
      int NodeIdx = Field::CrossFieldGraph::NodeOfElemDir(IndexV, DirIdx);
      int OppNodeIdx = Field::CrossFieldGraph::OppositeNodeCrossField(NodeIdx);
      std::vector<int> OrthoN;
      Field::CrossFieldGraph::OrthoNodesCrossField(NodeIdx, OrthoN);

      // std::cout<<"Node 0 "<<NodeIdx<<std::endl;
      // std::cout<<"Node 1 "<<OppNodeIdx<<std::endl;
      // std::cout<<"Node 2 "<<OrthoN[0]<<std::endl;
      // std::cout<<"Node 3 "<<OrthoN[1]<<std::endl;
      assert(NodeIdx >= 0);
      assert(NodeIdx < TGraph.NumNodes());
      TGraph.NodeType(NodeIdx) = Geo::Internal;
      assert(OppNodeIdx >= 0);
      assert(OppNodeIdx < TGraph.NumNodes());
      TGraph.NodeType(OppNodeIdx) = Geo::ExitDomain;
      assert(OrthoN[0] >= 0);
      assert(OrthoN[0] < TGraph.NumNodes());
      TGraph.NodeType(OrthoN[0]) = Geo::DartOrthoDomain;
      assert(OrthoN[1] >= 0);
      assert(OrthoN[1] < TGraph.NumNodes());
      TGraph.NodeType(OrthoN[1]) = Geo::DartOrthoDomain;
    }
  }

  void UpdateOrthoDartNodeTypeIfNeeded() {
    std::vector<int> DartIdxNode;
    PathSampl.Graph().GetNodeType(Geo::DartOrthoDomain, DartIdxNode);
    if (DartIdxNode.size() == 0)
      return;

    // then get all spike nodes
    std::vector<int> SpikeVert;
    PatchMan.GetAllGlobalVertOfType(CTypeSpike, SpikeVert);

    std::set<int> SpikeVertSet(SpikeVert.begin(), SpikeVert.end());

    for (size_t i = 0; i < DartIdxNode.size(); i++) {
      int IndexN = DartIdxNode[i];
      int IndexV = Field::CrossFieldGraph::ElemOfNode(IndexN);
      // then no more spike change the type of node
      if (SpikeVertSet.count(IndexV) == 0)
        PathSampl.Graph().NodeType(IndexN) = Geo::Internal;
    }
  }

public:
  void
  AddFeaturesAsPaths(Geo::TracingGraph<ScalarType> &TGraph,
                     const std::vector<Field::CrossF<ScalarType>> &VertCurv,
                     const std::vector<std::pair<int, int>> &Features,
                     const std::vector<int> &NewBorder) 
  {
    VertSplitForFeatures=std::set<int>(NewBorder.begin(),NewBorder.end());

    if (Features.size() == 0)
      return;
    std::vector<Geo::Point3<ScalarType>> EdgeVertPos;
    std::vector<std::pair<int, int>> Edges;
    Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
        PatchMan.VertPos, PatchMan.Faces, Features, EdgeVertPos, Edges);
    std::vector<std::vector<int>> Sequences;
    Geo::EdgeMeshFunctions<ScalarType>::ExtractSortedNodesSequences(
        EdgeVertPos, Edges, Sequences);

    // then transform vertices into the original vertices
    std::map<Geo::Point3<ScalarType>, int> PosMap;
    for (size_t i = 0; i < PatchMan.VertPos.size(); i++)
      PosMap[PatchMan.VertPos[i]] = i;

    for (size_t i = 0; i < Sequences.size(); i++)
      for (size_t j = 0; j < Sequences[i].size(); j++) {
        int IndexVE = Sequences[i][j];
        Geo::Point3<ScalarType> PosVE = EdgeVertPos[IndexVE];
        assert(PosMap.count(PosVE) > 0);
        Sequences[i][j] = PosMap[PosVE];
      }

    // transform into sequences
    int Index0 = PathSampl.SelectedPaths.size();
    PathSampl.AddFeaturePathsFromVert(Sequences);
    int Index1 = PathSampl.SelectedPaths.size();

    // set those as fixed
    for (size_t i = Index0; i < Index1; i++)
      FixedPaths.insert(i);

    // then update the patch Manager
    std::vector<int> UpdatedPatches;
    std::vector<int> ToAdd, ToRemove;
    for (size_t i = Index0; i < Index1; i++)
      ToAdd.push_back(i);

    PatchMan.Update(ToAdd, ToRemove, UpdatedPatches);

    // update the spike vertices with only one neigh
    std::map<int, std::vector<int>> AdjList;
    for (size_t i = 0; i < UpdatedPatches.size(); i++) {
      int TestPatch = UpdatedPatches[i];
      std::vector<int> CornerVert;
      std::vector<CornerType> CurrType;
      PatchMan.GetAllCorners(TestPatch, CornerVert, CurrType);
      for (size_t j = 0; j < CurrType.size(); j++) {
        if (CurrType[j] != Geo::CTypeSpike)
          continue;
        int IndexV = CornerVert[j];
        if (AdjList[IndexV].size() != 1)
          continue;
        int IndexV0 = AdjList[IndexV][0];
      }
    }

    // then update conditions
    for (size_t i = 0; i < UpdatedPatches.size(); i++) {
      int TestPatch = UpdatedPatches[i];
      assert(TestPatch < PatchMan.NumPatches());
      if (PatchMan.isEmpty(TestPatch))
        continue;
      UpdateCondition(TestPatch);
    }

    // update unsolved patches
    UnsolvedPatches.clear();
    for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
      if (PatchMan.isEmpty(i))
        continue;

      if (!SatisfyConditions(i))
        UnsolvedPatches.insert(i);
    }

    for (size_t i = 0; i < PConditions.size(); i++) {
      PConditions[i]->UpdateGlobalData(PatchMan);
    }

    for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
      if (PatchMan.isEmpty(i))
        continue;

      if (!SatisfyConditionsAfterGlobalUpdate(i))
        UnsolvedPatches.insert(i);
    }

    RemoveEmptyPatches();

    // if (add_only_where_needed)
    //   MarkUnsolvedPatches(true, false);

    UpdateGraphSingleLineSpike(TGraph, VertCurv, Features);
    // if (stop_when_solved)
    //   all_solved = (UnsolvedPatches.size() == 0);
  }

  void BatchDecompose(bool Do_remove) {
    // bool has_added = false;

    PathSampl.dynamicUpdates = dynamic_updates;

    PathSampl.Candidates.clear();
    PathSampl.SampledPaths_Status.clear();
    PathSampl.Init();

    // Init();
    // PathSampl.Init();

    TimeTracing = 0;
    TimeUpdateConditions = 0;
    TimeUpdateDistances = 0;
    TimeUpdateUpdatePatches = 0;
    TimeCheckConflicts = 0;
    TimeRemovePaths = 0;

    DecomposeCycle(Do_remove);

    // final step, split and remove if possible the unnecessary
    auto t0 = std::chrono::high_resolution_clock::now();
    if (split_removal) {
      RemoveSplitted();
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    TimeRemovePaths += elapsed.count() * 1e-6;

    // PatchMan.InitVisualCornersEdges();

    std::cout << "<<< Decompose TIMINGS (sec) >>> " << std::endl;
    std::cout << "* Tracing: " << TimeTracing << std::endl;
    std::cout << "* Updating Conditions: " << TimeUpdateConditions << std::endl;
    std::cout << "* Updating Patches: " << TimeUpdateUpdatePatches << std::endl;
    std::cout << "* Updating Distances: " << TimeUpdateDistances << std::endl;
    std::cout << "* Check Conflicts: " << TimeCheckConflicts << std::endl;
    std::cout << "* Removing: " << TimeRemovePaths << std::endl;
  }

  void BatchMaxDecompose(std::vector<int> &RankPaths,
                         bool add_incomplete = true) {
    // bool has_added = false;

    PathSampl.dynamicUpdates = dynamic_updates;

    PathSampl.Candidates.clear();
    PathSampl.SampledPaths_Status.clear();
    PathSampl.Init();

    TimeTracing = 0;
    TimeUpdateConditions = 0;
    TimeUpdateUpdatePatches = 0;
    TimeUpdateDistances = 0;
    TimeCheckConflicts = 0;
    TimeRemovePaths = 0;

    bool has_added = false;
    bool added_incomplete = false;
    RankPaths.clear();

    do {

      has_added = false;
      SetFeatureFeatureConfiguration();
      InitForceSamplingSpace();
      std::vector<bool> force_sample_swap = PathSampl.force_sample;

      // std::cout<<"*** STEP 1 ***"<<std::endl;
      bool has_added_feature = false;
      // add all features
      add_only_where_needed = false;
      stop_when_solved = false;
      do {
        PathSampl.force_sample = force_sample_swap;
        has_added_feature = SolveAllFeatures();
        has_added |= has_added_feature;
      } while (has_added_feature);
      RankPaths.resize(PathSampl.SelectedPaths.size(), 0);

      // previously inserted cannot be removed
      for (size_t i = 0; i < PathSampl.SelectedPaths.size(); i++)
        FixedPaths.insert(i);

      // then add the mandatory
      add_only_where_needed = true;
      stop_when_solved = true;
      // std::cout<<"*** STEP 2 ***"<<std::endl;
      has_added |= SolveAllMesh();
      // then remove unnecessary
      RemovePaths();
      PathSampl.CleanEmptySelectedPaths();

      RankPaths.resize(PathSampl.SelectedPaths.size(), 1);

      std::vector<int> BoundaryP;
      for (size_t i = 0; i < PathSampl.SelectedPaths.size(); i++)
        BoundaryP.push_back(i);

      std::cout << "*** STEP 3 ***" << std::endl;
      add_only_where_needed = false;
      stop_when_solved = false;
      SolveAllMesh();
      RankPaths.resize(PathSampl.SelectedPaths.size(), 2);

      if (add_incomplete && (!added_incomplete)) {
        PathSampl.SetOrthoCanPassToExit(BoundaryP);
        SolveAllMesh();
        PathSampl.RestoreNodeType();
        RankPaths.resize(PathSampl.SelectedPaths.size(), 2);
        added_incomplete = true;
      }

    } while (has_added);
    // then add the rest of candidates

    std::cout << "<<< Decompose TIMINGS (sec) >>> " << std::endl;
    std::cout << "* Tracing: " << TimeTracing << std::endl;
    std::cout << "* Updating Conditions: " << TimeUpdateConditions << std::endl;
    std::cout << "* Updating Patches: " << TimeUpdateUpdatePatches << std::endl;
    std::cout << "* Updating Distances: " << TimeUpdateDistances << std::endl;
    std::cout << "* Check Conflicts: " << TimeCheckConflicts << std::endl;
    std::cout << "* Removing: " << TimeRemovePaths << std::endl;

    add_only_where_needed = true;
    stop_when_solved = true;
  }

  bool SaveFPatchesLayout(const std::string fPatchesFile) {
    CompactData();

    PatchMan.UpdateFaceToPatch();
    std::ofstream fout(fPatchesFile.c_str());
    if (!fout)
      return false;

    // write number of faces
    fout << PatchMan.NumFaces() << "\n";

    for (size_t i = 0; i < PatchMan.PData.OriginalFaceToPatch.size(); i++)
      fout << PatchMan.PData.OriginalFaceToPatch[i] << "\n";

    fout.close();
    return true;
  }

  bool SaveCorners(const std::string cPatchesFile) {
    CompactData();

    std::ofstream fout(cPatchesFile.c_str());
    if (!fout)
      return false;

    // write number of patches
    fout << PatchMan.NumPatches() << "\n";

    for (size_t pi = 0; pi < PatchMan.NumPatches(); pi++) {
      assert(PatchMan.PData.PatchCorners[pi].size() == 1);
      int NumC = PatchMan.PData.PatchCorners[pi][0].size();
      fout << NumC << "\n";
      for (size_t ci = 0; ci < PatchMan.PData.PatchCorners[pi][0].size();
           ci++) {
        int IndexC = PatchMan.PData.PatchCorners[pi][0][ci];
        IndexC = PatchMan.MapToOriginal(pi, IndexC);
        assert(IndexC >= 0);
        assert(IndexC < PatchMan.VertPos.size());
        fout << IndexC << "\n";
      }
    }
    fout.close();
    return true;
  }

  // bool SaveMeshColorByLayout(const std::string meshFile)
  // {

  // }

  bool SavePatchesLayout(const std::string meshFile) {
    CompactData();
    std::ofstream fout(meshFile.c_str());
    if (!fout)
      return false;

    // first set the number of patches
    fout << "#p" << PatchMan.PData.EdgePatchSides.size() << "\n";

    // for each patch
    for (size_t pi = 0; pi < PatchMan.PData.EdgePatchSides.size(); pi++) {
      // write patch index
      fout << "p " << pi << "\n";

      // for each loops
      fout << "#l " << PatchMan.PData.EdgePatchSides[pi].size() << "\n";
      assert(PatchMan.PData.EdgePatchSides[pi].size() ==
             PatchMan.PData.PatchCorners[pi].size());

      // for each loop
      for (size_t li = 0; li < PatchMan.PData.PatchCorners[pi].size(); li++) {
        fout << "l " << li << "\n";
        // number of sides equal nummber of corners
        assert(PatchMan.PData.PatchCorners[pi][li].size() ==
               PatchMan.PData.EdgePatchSides[pi][li].size());

        // write corners
        fout << "#c " << PatchMan.PData.PatchCorners[pi][li].size() << "\n";

        for (size_t ci = 0; ci < PatchMan.PData.PatchCorners[pi][li].size();
             ci++) {
          int IndexC = PatchMan.PData.PatchCorners[pi][li][ci];
          IndexC = PatchMan.MapToOriginal(pi, IndexC);
          assert(IndexC >= 0);
          assert(IndexC < PatchMan.VertPos.size());
          fout << "c " << IndexC << "\n";
        }

        // write number of sides
        fout << "#s " << PatchMan.PData.EdgePatchSides[pi][li].size() << "\n";
        for (size_t si = 0; si < PatchMan.PData.EdgePatchSides[pi][li].size();
             si++) {
          fout << "s " << si << "\n";
          int sizeL = PatchMan.PData.EdgePatchSides[pi][li][si].size();
          fout << "#e " << sizeL << "\n";

          // for each edge
          for (size_t ei = 0;
               ei < PatchMan.PData.EdgePatchSides[pi][li][si].size(); ei++) {
            int IndexV0 = PatchMan.PData.EdgePatchSides[pi][li][si][ei].first;
            int IndexV1 = PatchMan.PData.EdgePatchSides[pi][li][si][ei].second;
            IndexV0 = PatchMan.MapToOriginal(pi, IndexV0);
            IndexV1 = PatchMan.MapToOriginal(pi, IndexV1);
            assert(IndexV0 >= 0);
            assert(IndexV0 < PatchMan.VertPos.size());
            assert(IndexV1 >= 0);
            assert(IndexV1 < PatchMan.VertPos.size());
            fout << "e " << IndexV0 << " " << IndexV1 << "\n";
          }
        }
      }
    }
    fout.close();
    return true;
  }

  void WriteStats() {
    for (size_t i = 0; i < PConditions.size(); i++) {
      int NonSatisfied = 0;
      for (size_t j = 0; j < PatchMan.NumPatches(); j++) {
        if (PatchMan.isEmpty(j))
          continue;
        bool has_satisfied = PConditions[i]->IsCorrect(PatchMan, j);
        if (!has_satisfied)
          NonSatisfied++;
      }
      std::cout << "Unsatistied Condition " << PConditions[i]->Name().c_str()
                << " " << NonSatisfied << std::endl;
    }

    // std::set<int> MandatoryP;
    // GetMandatoryPatch(MandatoryP);
    // if (MandatoryP.size() > 0)
    //   std::cout << "Unsatistied MANDATAORY Condition " << MandatoryP.size()
    //             << std::endl;
  }

  void SaveSubMeshes(const std::string PathSave) {
    CompactData();

    std::vector<Geo::Point3<ScalarType>> PosVect;
    std::vector<std::vector<int>> FaceVect;
    assert(PatchMan.PData.SubPatchVertPos.size() ==
           PatchMan.PData.SubPatchFaces.size());

    for (size_t i = 0; i < PatchMan.PData.SubPatchVertPos.size(); i++) {
      int baseIdx = PosVect.size();
      PosVect.insert(PosVect.end(), PatchMan.PData.SubPatchVertPos[i].begin(),
                     PatchMan.PData.SubPatchVertPos[i].end());

      for (size_t j = 0; j < PatchMan.PData.SubPatchFaces[i].size(); j++) {
        FaceVect.push_back(std::vector<int>(3, -1));
        assert(PatchMan.PData.SubPatchFaces[i][j].size() == 3);
        for (size_t k = 0; k < 3; k++) {
          int IndexV = PatchMan.PData.SubPatchFaces[i][j][k];
          FaceVect.back()[k] = IndexV + baseIdx;
        }
      }
    }
    WriteOBJ(PathSave.c_str(), PosVect, FaceVect);
  }

  //    ~PatchDecomposer()
  //    {
  //        for (size_t i=0;i<PConditions.size();i++)
  //            delete(PConditions[i]);
  //    }
};

} // namespace Geo
#endif