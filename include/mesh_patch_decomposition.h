/*
 * //Non-Commercial License Agreement
 * //Copyright (c) 2025, Nico Pietroni
 * //
 * //1. Permission is hereby granted to use, copy, modify, and distribute
 * //   this software for academic research, teaching, and personal
 * //   non-commercial purposes, provided that this copyright notice is
 * //   included in all copies or substantial portions of the software.
 * //
 * //2. Commercial use of this software, including but not limited to use in
 * //   commercial products, services, tools, or consulting, is strictly
 * //   prohibited without prior written permission from the copyright holder.
 * //
 * //3. This software is provided "as is", without warranty of any kind,
 * //   express or implied, including but not limited to the warranties of
 * //   merchantability, fitness for a particular purpose, and
 * //   noninfringement. In no event shall the authors or copyright holders
 * //   be liable for any claim, damages, or other liability, arising from,
 * //   out of, or in connection with the software or the use or other
 * //   dealings in the software.
 * //
 * //4. For commercial licensing, please contact:
 * //   Nico Pietroni
 * //   Email: nico.pietroni@gmail.com
 */

#ifndef MESH_PATCH_DECOMPOSITION
#define MESH_PATCH_DECOMPOSITION

#include "IO/io_obj.h"
#include "IO/io_path.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include <connectivity.h>
#include <cross_field.h>
#include <field_graph/mesh_preprocess.h>
#include <mesh_features.h>
#include <point3.h>

#include <field_graph/dual_patches.h>
#include <field_graph/graph.h>
#include <field_graph/graph_tracing.h>
#include <field_graph/initialize_graph.h>
#include <field_graph/mesh_preprocess.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_managing.h>
#include <field_graph/patch_optimize.h>
#include <field_graph/path_functions.h>
#include <field_graph/path_sampling.h>
#include <field_graph/tracing_interface.h>

#include <patch_side_discretization.h>

#include "edge_mesh_functions.h"
#include "mesh_mapping.h"
#include "smooth_mesh.h"
#include <quad_mesh_functions.h>
#include <quadmesh_smooth_separatrix.h>

template <class ScalarType> class MeshPatchDecomposition {
  // MESH DATA
  std::vector<Geo::Point3<ScalarType>> &TriVertPos;
  std::vector<std::vector<int>> &TriConnectivity;
  std::vector<Geo::Point3<ScalarType>> TriFaceNormals;
  std::vector<Geo::Point3<ScalarType>> TriVertNormals;
  std::vector<std::vector<int>> NextF, NextE;
  Geo::Box3<ScalarType> MeshBox;

  // FIELD DATA
  std::vector<Field::CrossF<ScalarType>> &FaceCurv;
  std::vector<Field::CrossF<ScalarType>> VertCurv;
  std::vector<int> SingIndex, SingValue;

  // SETUP DATA
  int min_sing_distance;

  // // DUAL DATA
  // std::vector<Geo::Point3<ScalarType>> DualVertPos;
  // std::vector<std::vector<int>> DualConnectivity;
  // std::vector<std::pair<int, int>> DualFeatures;
  // std::vector<int> DualCorners;

  // GRAPH SETUP DATA
  ScalarType graph_drift_penalty;
  ScalarType graph_maxAngle;
  int graph_diffuse_step;

  // OTHER DATA
  std::vector<std::pair<int, int>> &TriFeatures;
  ScalarType feature_corner_angle;
  std::vector<int> TriCornersIdx;
  std::vector<int> FeatureSingIndex, FeatureSingValue;
  std::vector<std::pair<int, int>> RemainingF;
  std::vector<int> NewBorder;

  // DATA FOR REMOVE COINCIDENT VERT
  std::vector<int> ModSameBVert;
  std::vector<Geo::Point3<ScalarType>> ModSameBOrigPos;

  // DIFFERENT CLASS TO PATCH DECOMPOSITIN
  Geo::TracingGraph<ScalarType> TGraph;
  Geo::CrossPathSampling<ScalarType> PathSampl;
  Geo::PatchManaging<ScalarType> PatchM;
  Geo::PatchDecomposer<ScalarType> PDeco;

  void CleanMeshArtifacts() {
    // ScalarType AVLen = AvgEdgeLen(TriVertPos, TriConnectivity);
    Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
        TriVertPos, TriConnectivity);
  }

  void InitMesh() {
    MeshBox.Init(TriVertPos);
    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);
    CleanMeshArtifacts();
  }

  void InitField() {
    SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormals, FaceCurv,
                         VertCurv);

    Geo::MeshFeatures<ScalarType>::Corners(TriVertPos, TriFeatures,
                                           feature_corner_angle, TriCornersIdx);
  }

  void InitPatches() {
    // initialize patches corners
    // and all remaining data

    PatchM.CheckConsistentData();
    PatchM.Init();

    PatchM.CheckConsistentData();

    // then initialize the stopping conditions
    // InitConditions();
    PatchM.CheckConsistentData();

    // PDeco.InitBorderVertType();

    // initialise the decomposer
    PDeco.Init();
    PatchM.CheckConsistentData();

    // then check all patches and saved the unsolved ones
    PatchM.CheckConsistentData();
    PDeco.UpdateAllUnsolvedPatches();

    PDeco.AddFeaturesAsPaths(TGraph, VertCurv, RemainingF, NewBorder);
  }

  void UpdateAfterRemesh() {
    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);
    // ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(TriConnectivity, NextF, NextE);

    InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

    SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormals, FaceCurv,
                         VertCurv);

    GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                         SingIndex, SingValue);

    PatchM.SetSingularities(SingIndex, SingValue);
    PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
  }

  void PreprocessMesh() {

    Geo::RefineForTracing<ScalarType>::RefineForConnection(
        TriVertPos, TriConnectivity, FaceCurv, TriFeatures);

    Geo::getFFAdjacency(TriConnectivity, NextF, NextE);

    GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                         SingIndex, SingValue);

    Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(
        TriVertPos, TriConnectivity, FaceCurv, TriFeatures, SingIndex,
        min_sing_distance);

    UpdateAfterRemesh();

    std::vector<Field::CrossF<ScalarType>> VertCurv0;

    std::cout << "*** SPLIT FEATURES *** " << std::endl;

    RemainingF.clear();
    Geo::TracingPreprocess<ScalarType>::SplitCornerFeatures(
        TriVertPos, TriConnectivity, VertCurv, FaceCurv, TriFeatures,
        RemainingF, SingIndex, SingValue, FeatureSingIndex, FeatureSingValue,
        NewBorder);

    Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
        TriVertPos, TriConnectivity, &ModSameBVert, &ModSameBOrigPos);

    UpdateAfterRemesh();
  }

  void InitGraph() {
    Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(
        TGraph, TriVertPos, TriConnectivity, VertCurv, graph_drift_penalty,
        graph_maxAngle, graph_diffuse_step);

    Geo::PathFunctions<ScalarType>::DisableNonDirectableConnections(TGraph);
  }

  void SmoothPaths() {
    typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
    std::vector<std::vector<int>> VertPaths;
    std::vector<Geo::TJunction> TJunctions;

    PathSampl.GetVertexSelectedPaths(VertPaths);
    Geo::PathFunctions<ScalarType>::FindTJunctions(VertPaths, TJunctions);

    Geo::PatchOptimize<ScalarType>::SmoothPaths(TriVertPos, TriConnectivity,
                                                VertPaths, TriFeatures);

    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);

    // update cross F poistion
    for (size_t i = 0; i < VertCurv.size(); i++)
      VertCurv[i].Pos = TriVertPos[i];

    AdaptFieldToMesh(TriVertPos, TriConnectivity, VertCurv, FaceCurv);

    PatchM.UpdateSubPatchPos();
    PatchM.InitVisualCornersEdges();
  }

  void RestoreOriginalPosOnSplitBorders() {
    assert(ModSameBVert.size() == ModSameBOrigPos.size());
    for (size_t i = 0; i < ModSameBVert.size(); i++) {
      int IndexV = ModSameBVert[i];
      assert(IndexV >= 0);
      assert(IndexV < TriVertPos.size());
      TriVertPos[IndexV] = ModSameBOrigPos[i];
    }
  }

  // void ExtractDualPatches() {
  //   RestoreOriginalPosOnSplitBorders();
  //   Geo::DualPatch<ScalarType>::DisambiguateForDual(PatchM);

  //   typename Geo::DualPatch<ScalarType>::PrimalPatchData PData;
  //   PData.PatchFaces = PatchM.PData.SubPatchFaces;
  //   PData.PatchVertPos = PatchM.PData.SubPatchVertPos;
  //   PData.PatchCorners = PatchM.PData.PatchCorners;
  //   PData.EdgePatchSides = PatchM.PData.EdgePatchSides;
  //   for (size_t i = 0; i < PatchM.NumPatches(); i++) {
  //     std::vector<int> PatchSing;
  //     PatchM.PatchSingularities(i, PatchSing);
  //     PData.PatchSingularities.push_back(PatchSing);
  //   }

  //   std::vector<std::vector<int>> PerDualPatchFaceIDX;
  //   std::vector<std::vector<int>> PerDualPatchCornerIDX;

  //   int num_problematic_patches =
  //   Geo::DualPatch<ScalarType>::ExtractDual(
  //       PData, DualParam, DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
  //       PerDualPatchCornerIDX);

  //   has_problematic_patches|=(num_problematic_patches!=0);

  //   std::vector<int> PrimalCorners;
  //   Geo::DualPatch<ScalarType>::RetrieveDualFeatures(
  //       TriVertPos, TriConnectivity, TriFeatures, PrimalCorners, DualVertPos,
  //       DualConnectivity, DualFeatures, DualCorners);

  //   std::vector<std::vector<int>> DualNextF, DualNextE;
  //   Geo::getFFAdjacency(DualConnectivity, DualNextF, DualNextE);

  //   // INITIALIZE SIDE SUBDIVISION
  //   SideD.InitFromLabels(DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
  //                        PerDualPatchCornerIDX, DualFeatures);
  // }

  // void InitDefaultDualParameters() {
  //   DualParam.OnlyNonQuad = true;
  //   DualParam.InitSingPos = true;
  //   DualParam.FixSingPos = false;
  //   DualParam.JointBorder = true;
  // }

  // void UpdateQuadFeatures() {

  //   std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
  //   FeatureData;

  //   SideD.GetFeaturesQuadrData(DualVertPos, DualFeatures, FeatureData);

  //   QuadFeatures.clear();

  //   bool has_mapped = Geo::FindEdgePosCorrespondences(QuadVert, QuadFaces,
  //                                                FeatureData, QuadFeatures);

  //   if (!has_problematic_patches)
  //   {
  //     assert(has_mapped);
  //   }
  //   // find the corners
  //   QuadCorners.clear();

  //   Geo::MeshFeatures<ScalarType>::CornersByCount(QuadVert, QuadFeatures,
  //   QuadCorners,true);

  //   //add the geometric ones
  //   std::set<Geo::Point3<ScalarType>> CornerPosSet;
  //   for (size_t i = 0; i < TriCornersIdx.size(); i++)
  //     CornerPosSet.insert(TriVertPos[TriCornersIdx[i]]);

  //   for (size_t i = 0; i < QuadVert.size(); i++) {
  //     if (CornerPosSet.count(QuadVert[i])== 0)continue;
  //     QuadCorners.push_back(i);
  //   }

  // }

  void InitializeData() {
    std::cout << "*** INITIALIZATION *** " << std::endl;
    InitMesh();
    InitField();

    std::cout << "*** REFINE AND SPLIT *** " << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    PreprocessMesh();
    auto end0 = std::chrono::high_resolution_clock::now();
    auto elapsed0 =
        std::chrono::duration_cast<std::chrono::microseconds>(end0 - start);

    std::cout << "*** INIT GRAPHS *** " << std::endl;
    InitGraph();

    InitPatches();
  }

  void DecomposePatches() {
    bool remove_batch = true;
    PDeco.add_only_where_needed = true;
    PDeco.stop_when_solved = true;

    PDeco.BatchDecompose(remove_batch);

    SmoothPaths();

    RestoreOriginalPosOnSplitBorders();

    SmoothPaths();

    PatchM.WriteStats();
    PDeco.WriteStats();

    // ExtractDualPatches();
  }

public:
  bool single_boundary_cond;
  bool disk_like_cond;
  bool self_connection_cond;
  bool single_sing_cond;
  bool match_sing_cond;
  int MinSides;
  int MaxSides;
  ScalarType CCAbility;
  bool prefer_feature_features;
  bool dynamicUpdates;
  int num_samples;

  void AddExtraCondition(Geo::PatchCondition<ScalarType> *PCond) {
    PDeco.AddCondition(PCond);
  }

  void InitConditions() {
    PDeco.ClearConditions();
    if (single_boundary_cond)
      PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
    if (disk_like_cond)
      PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
    if (self_connection_cond)
      PDeco.AddCondition(new Geo::SelfConnectingPatchCondition<ScalarType>());
    if (single_sing_cond)
      PDeco.AddCondition(new Geo::SingleSingularityCondition<ScalarType>());
    if (match_sing_cond)
      PDeco.AddCondition(new Geo::MathchSingularityCondition<ScalarType>());
    if ((MinSides > 0) || (MaxSides > 0))
      PDeco.AddCondition(
          new Geo::NumSidesCondition<ScalarType>(MinSides, MaxSides));

    if (CCAbility > 0)
      PDeco.AddCondition(new Geo::CCabilityCondition<ScalarType>(CCAbility));
  }

  MeshPatchDecomposition(std::vector<Geo::Point3<ScalarType>> &_TriVertPos,
                         std::vector<std::vector<int>> &_TriConnectivity,
                         std::vector<Field::CrossF<ScalarType>> &_FaceCurv,
                        std::vector<std::pair<int, int>> &_TriFeatures
                      )
      : TriVertPos(_TriVertPos), TriConnectivity(_TriConnectivity),
        FaceCurv(_FaceCurv), TGraph(Geo::GT_CrossField),
        PathSampl(TriVertPos, TriConnectivity, TGraph, VertCurv, SingIndex,
                  SingValue),
        PatchM(TriVertPos, TriConnectivity, VertCurv, TGraph,
               PathSampl.SelectedPaths),
        PDeco(PatchM, PathSampl, VertCurv, FaceCurv),
        TriFeatures(_TriFeatures)

  {

    graph_drift_penalty = 100;
    graph_maxAngle = 45;
    graph_diffuse_step = 1;
    min_sing_distance = 1;

    dynamicUpdates = false;
    num_samples = 1000;

    feature_corner_angle = 45;

    single_boundary_cond = true;
    disk_like_cond = true;
    self_connection_cond = true;
    single_sing_cond = true;
    match_sing_cond = true;
    prefer_feature_features = false;

    MinSides = 3;
    MaxSides = 6;
    CCAbility = -1;
  }

  Geo::PatchManaging<ScalarType> &PatchManager() { return PatchM; }

  void GetVertexPaths(std::vector<std::vector<int>> &VertPaths) {
    PathSampl.GetVertexSelectedPaths(VertPaths);
  }

  void ExtractPatches(std::vector<std::vector<int>> &VertPaths) {
    //TriFeatures = _TriFeatures;
    PDeco.prefer_feature_features = prefer_feature_features;

    std::sort(TriFeatures.begin(), TriFeatures.end());
    auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
    TriFeatures.erase(last, TriFeatures.end());

    PathSampl.dynamicUpdates = dynamicUpdates;
    PathSampl.NumSamples = num_samples;

    InitializeData();

    DecomposePatches();

    PatchM.GetVertexPaths(VertPaths);
  }

  // void ExtractPatches(std::vector<std::pair<int, int>> &_TriFeatures,
  //                     std::vector<std::pair<int, int>> &PathBoundaries,
  //                     std::vector<int> &PathCorners) {
  //   std::vector<std::vector<int>> VertPaths;
  //   ExtractPatches(_TriFeatures, VertPaths);
  //   PatchM.GetVertexPaths(VertPaths);
  //   PatchM.GetAllSideGlobalEdges(PathBoundaries);
  //   PatchM.GetAllCorners(PathCorners);
  // }
  void ExtractPatches(std::vector<std::pair<int, int>> &PathBoundaries,
                      std::vector<int> &PathCorners) {
    std::vector<std::vector<int>> VertPaths;
    ExtractPatches(VertPaths);
    PatchM.GetVertexPaths(VertPaths);
    PatchM.GetAllSideGlobalEdges(PathBoundaries);
    PatchM.GetAllCorners(PathCorners);
  }
};

#endif