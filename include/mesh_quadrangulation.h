#ifndef MESH_QUADRANGULATION
#define MESH_QUADRANGULATION

#include "IO/io_obj.h"
#include "IO/io_path.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <set>
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
#include "smooth_mesh.h"
#include <quad_mesh_functions.h>
#include <quadmesh_smooth_separatrix.h>
#include <quadrangulation.h>
#include <quadrangulationLUT.h>

template <class ScalarType> class MeshQuadrangulator {
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

  // QUAD DATA
  std::vector<Geo::Point3<ScalarType>> &QuadVert;
  std::vector<std::vector<int>> &QuadFaces;
  // std::vector<Geo::Point3<ScalarType>> QuadFaceNormals, QuadVertNormals;
  std::vector<std::pair<int, int>> QuadFeatures;
  std::vector<int> QuadCorners;

  // SETUP DATA
  int min_sing_distance;

  // DUAL DATA
  std::vector<Geo::Point3<ScalarType>> DualVertPos;
  std::vector<std::vector<int>> DualConnectivity;
  // std::vector<Geo::Point3<ScalarType>> DualFaceNormals;
  // std::vector<Geo::Point3<ScalarType>> DualVertNormals;
  // std::vector<std::vector<int>> PerDualPatchFaceIDX;
  // std::vector<std::vector<int>> PerDualPatchCornerIDX;
  std::vector<std::pair<int, int>> DualFeatures;
  std::vector<int> DualCorners;

  // GRAPH SETUP DATA
  ScalarType graph_drift_penalty;
  ScalarType graph_maxAngle;
  int graph_diffuse_step;

  bool force_parity;
  bool OptimizeSize;

  // OTHER DATA
  std::vector<std::pair<int, int>> TriFeatures;
  ScalarType feature_corner_angle;
  std::vector<int> TriCornersIdx;
  std::vector<int> FeatureSingIndex, FeatureSingValue;
  std::vector<std::pair<int, int>> RemainingF;
  std::vector<int> NewBorder;

  // DATA FOR REMOVE COINCIDENT VERT
  std::vector<int> ModSameBVert;
  std::vector<Geo::Point3<ScalarType>> ModSameBOrigPos;

  // DIFFERENT CLASS TO PATCH DECOMPOSITIN AND QUAD
  Geo::TracingGraph<ScalarType> TGraph;
  Geo::SideDiscretization<ScalarType> SideD;
  Geo::CrossPathSampling<ScalarType> PathSampl;
  Geo::PatchManaging<ScalarType> PatchM;
  Geo::PatchDecomposer<ScalarType> PDeco;
  typename Geo::DualPatch<ScalarType>::DualParameter DualParam;
  Geo::SmoothParam<ScalarType> SParamQuadSmooth;
  Geo::QuadrangulationLUT<ScalarType> Q_LUT;
  bool has_problematic_patches;
  std::vector<int> QuadPatch;

  void CleanMeshArtifacts() {
    ScalarType AVLen = AvgEdgeLen(TriVertPos, TriConnectivity);
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

    Geo::MeshFeatures<ScalarType>::Corners(TriVertPos, TriFeatures, feature_corner_angle,
                             TriCornersIdx);
  }

  void InitConditions() {
    PDeco.ClearConditions();
    // if (single_boundary_cond)
    PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
    // if (disk_like_cond)
    PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
    // if (self_connection_cond)
    PDeco.AddCondition(new Geo::SelfConnectingPatchCondition<ScalarType>());
    // if (single_sing_cond)
    //   PDeco.AddCondition(new Geo::SingleSingularityCondition<ScalarType>());
    // if (match_sing_cond)
    PDeco.AddCondition(new Geo::MathchSingularityCondition<ScalarType>());
    // if ((MinSides > 0) || (MaxSides > 0)) {
    PDeco.AddCondition(new Geo::NumSidesCondition<ScalarType>());

    if (CCAbility > 0)
      PDeco.AddCondition(new Geo::CCabilityCondition<ScalarType>(CCAbility));
  }

  void InitPatches() {
    // initialize patches corners
    // and all remaining data

    PatchM.CheckConsistentData();
    PatchM.Init();

    PatchM.CheckConsistentData();

    // then initialize the stopping conditions
    InitConditions();
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

    // std::cout<<"Features size:"<<TriFeatures.size()<<std::endl;
    // std::cout<<"VertPath size:"<<VertPaths.size()<<std::endl;
    // WriteFeatures("feature.feat",TriFeatures);
    // WritePath("path.path",VertPaths);
    Geo::PatchOptimize<ScalarType>::SmoothPaths(TriVertPos, TriConnectivity,
                                                VertPaths, TriFeatures);
    // WriteOBJ("Primal.obj",TriVertPos, TriConnectivity);
    // exit(0);

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

  void ExtractDualPatches() {
    RestoreOriginalPosOnSplitBorders();
    Geo::DualPatch<ScalarType>::DisambiguateForDual(PatchM);

    typename Geo::DualPatch<ScalarType>::PrimalPatchData PData;
    PData.PatchFaces = PatchM.PData.SubPatchFaces;
    PData.PatchVertPos = PatchM.PData.SubPatchVertPos;
    PData.PatchCorners = PatchM.PData.PatchCorners;
    PData.EdgePatchSides = PatchM.PData.EdgePatchSides;
    for (size_t i = 0; i < PatchM.NumPatches(); i++) {
      std::vector<int> PatchSing;
      PatchM.PatchSingularities(i, PatchSing);
      PData.PatchSingularities.push_back(PatchSing);
    }

    std::vector<std::vector<int>> PerDualPatchFaceIDX;
    std::vector<std::vector<int>> PerDualPatchCornerIDX;

    Geo::DualPatch<ScalarType>::ExtractDual(
        PData, DualParam, DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
        PerDualPatchCornerIDX);

    std::vector<int> PrimalCorners;
    Geo::DualPatch<ScalarType>::RetrieveDualFeatures(
        TriVertPos, TriConnectivity, TriFeatures, PrimalCorners, DualVertPos,
        DualConnectivity, DualFeatures, DualCorners);

    // UPDATE MESH DATA
    // ComputeNormals(DualVertPos, DualConnectivity, DualFaceNormals,
    //                DualVertNormals);

    std::vector<std::vector<int>> DualNextF, DualNextE;
    Geo::getFFAdjacency(DualConnectivity, DualNextF, DualNextE);

    // INITIALIZE SIDE SUBDIVISION
    SideD.InitFromLabels(DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
                         PerDualPatchCornerIDX, DualFeatures);
  }

  void InitDefaultDualParameters() {
    DualParam.OnlyNonQuad = true;
    DualParam.InitSingPos = true;
    DualParam.FixSingPos = false;
    DualParam.JointBorder = true;
    // DualParam.RType=Geo::DualPatch<ScalarType>::RT_CCability;
  }

  void UpdateQuadFeatures() {
    std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        FeatureData;
    SideD.GetFeaturesQuadrData(DualVertPos, DualFeatures, FeatureData);

    QuadFeatures.clear();
    bool has_mapped = FindEdgePosCorrespondences(QuadVert, QuadFaces,
                                                 FeatureData, QuadFeatures);
    assert(has_mapped);

    // find the corners
    QuadCorners.clear();
   
    Geo::MeshFeatures<ScalarType>::CornersByCount(QuadVert, QuadFeatures, QuadCorners,true);

    //add the geometric ones
    std::set<Geo::Point3<ScalarType>> CornerPosSet;
    for (size_t i = 0; i < TriCornersIdx.size(); i++)
      CornerPosSet.insert(TriVertPos[TriCornersIdx[i]]);

    for (size_t i = 0; i < QuadVert.size(); i++) {
      if (CornerPosSet.count(QuadVert[i])== 0)continue;  
      QuadCorners.push_back(i);
    }

    //  //remove valence 1 vertices
    // std::vector<int> Valence1Idx;
    // Geo::MeshFeatures<ScalarType>::VerticesWithValence1(QuadFeatures,Valence1Idx);
    // std::set<int> Valence1IdxSet(Valence1Idx.begin(),Valence1Idx.end());
    // std::vector<int> NewQuadCorners;
    // for (size_t i = 0; i < QuadCorners.size(); i++) {
    //   if (Valence1IdxSet.count(QuadCorners[i])>0)continue;  
    //   NewQuadCorners.push_back(i);
    // }
    // QuadCorners=NewQuadCorners;
  }

  void Init3DSmooth() {
    SParamQuadSmooth.Clear();

    SParamQuadSmooth.Features = QuadFeatures;
    SParamQuadSmooth.FixedVert = QuadCorners;

    std::vector<std::pair<int, int>> BorderE;
    Geo::getBorderEdges(QuadFaces, BorderE);
    for (size_t i = 0; i < BorderE.size(); i++) {
      SParamQuadSmooth.FixedVert.push_back(BorderE[i].first);
      SParamQuadSmooth.FixedVert.push_back(BorderE[i].second);
    }

    SParamQuadSmooth.ReprojVert = &TriVertPos;
    SParamQuadSmooth.ReprojFaces = &TriConnectivity;
    SParamQuadSmooth.Damp = 0.5;
    SParamQuadSmooth.post_project = true;
  }

  // void SmoothCentroidMesh3D() {

  //   SParamQuadSmooth.NumIte = QuadSmoothSteps;
  //   SParamQuadSmooth.SType = Geo::SMCentroid;
  //   SParamQuadSmooth.check_AspectR = true;
  //   Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

  //   // Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
  //   //                                 QuadVertNormals);
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

    ExtractDualPatches();
  }

  void Quadrangulate() {

    ScalarType EdgeS = MeshBox.Diag() * QuadSize;

    // typename Geo::SideDiscretization<ScalarType>::solveMode SM;

    // SM = Geo::SideDiscretization<ScalarType>::SM_GlobalSolve;

    // SideD.Solve(EdgeS, SM, force_parity, RegVsIso);

    std::cout << "Edge Size:" << EdgeS << std::endl;
    SideD.PrioritizeTJunction = true;

    SideD.FindDiscretization(EdgeS, RegVsIso, OptimizeSize);
    int OddNum=SideD.NumOdd();
    //then reoptimize making them 
    if (OddNum>0)
    {
      std::cout << "Recompute with Even Variables" << std::endl;
      SideD.PMode=Geo::SideDiscretization<ScalarType>::EvenVar;
      SideD.FindDiscretization(EdgeS, RegVsIso, OptimizeSize);
    }
    std::cout << "GET QUADRANGULATION DATA" << std::endl;

    std::vector<std::vector<int>> PatchFaces, PatchCorners;
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> SideSamples;
    SideD.GetQuadrangulationData(PatchFaces, PatchCorners, SideSamples);

    int num_problematic_patches =
        Geo::Quadrangulation<ScalarType>::QuadrangulateQuadPatchsWithSubs(
            DualVertPos, DualConnectivity, PatchFaces, PatchCorners,
            SideSamples, Q_LUT, QuadVert, QuadFaces, QuadPatch);
    
    std::cout << "INIT SMOOTHING" << std::endl;

    has_problematic_patches=(num_problematic_patches>0);

    RemoveDuplicatedVert(QuadVert, QuadFaces);

    RemoveUnfererencedVert(QuadVert, QuadFaces);

    // Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
    //                                 QuadVertNormals);
    UpdateQuadFeatures();
    Init3DSmooth();
  }

  void SolveDoublets() {
    Geo::QuadPatchFunctions<ScalarType>::RemoveDoublets(
        QuadVert, QuadFaces, QuadFeatures, QuadCorners, QuadPatch);
  }

  void SmoothQuadMeshSeparatrix() {

    typename Geo::QuadMeshSeparatrixSmooth<ScalarType>::SeparatrixSmoothParam
        Sep_SParam;
    Sep_SParam.Fixed = QuadCorners;
    Sep_SParam.FeatureLines = QuadFeatures;
    Sep_SParam.Damp = 0.5;
    Sep_SParam.Steps = QuadSmoothSteps / 2;
    Sep_SParam.fixBorder = true;
    Sep_SParam.post_project=!(has_problematic_patches);
    Sep_SParam.ReprojVert = &TriVertPos;
    Sep_SParam.ReprojFaces = &TriConnectivity;
    Sep_SParam.ReprojFeatureLines = TriFeatures;
    Geo::QuadMeshSeparatrixSmooth<ScalarType>::Smooth(QuadVert, QuadFaces,
                                                      Sep_SParam);
  }

  void SmoothLaplacian() {
    Geo::SmoothParam<ScalarType> SParamQuadSmooth;

    SParamQuadSmooth.Clear();

    SParamQuadSmooth.Features = QuadFeatures;
    SParamQuadSmooth.FixedVert = QuadCorners;

    std::vector<std::pair<int, int>> BorderE;
    Geo::getBorderEdges(QuadFaces, BorderE);
    for (size_t i = 0; i < BorderE.size(); i++) {
      SParamQuadSmooth.FixedVert.push_back(BorderE[i].first);
      SParamQuadSmooth.FixedVert.push_back(BorderE[i].second);
    }

    SParamQuadSmooth.ReprojVert = &TriVertPos;
    SParamQuadSmooth.ReprojFaces = &TriConnectivity;
    SParamQuadSmooth.Damp = 0.5;
    SParamQuadSmooth.post_project = !(has_problematic_patches);
    SParamQuadSmooth.ReprojFeatures = TriFeatures;
    SParamQuadSmooth.NumIte = QuadSmoothSteps / 2;
    //SParamQuadSmooth.SType = Geo::SMLaplacian;
    SParamQuadSmooth.SType = Geo::SMBestPoly;
    SParamQuadSmooth.check_AspectR = false;

    Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);
  }

  void SmoothQuadrangulation() {
    // keep more field
    //SmoothQuadMeshSeparatrix();

    // laplacian smooth to solve some artifacts
    SmoothLaplacian();
  }

public:
  ScalarType QuadSize;
  int QuadSmoothSteps;
  ScalarType RegVsIso;
  ScalarType CCAbility;
  bool dynamicUpdates;
  int num_samples;

  MeshQuadrangulator(std::vector<Geo::Point3<ScalarType>> &_TriVertPos,
                     std::vector<std::vector<int>> &_TriConnectivity,
                     std::vector<Field::CrossF<ScalarType>> &_FaceCurv,
                     std::vector<Geo::Point3<ScalarType>> &_QuadVert,
                     std::vector<std::vector<int>> &_QuadFaces)
      : TriVertPos(_TriVertPos), TriConnectivity(_TriConnectivity),
        FaceCurv(_FaceCurv), QuadVert(_QuadVert), QuadFaces(_QuadFaces),
        TGraph(Geo::GT_CrossField),
        PathSampl(TriVertPos, TriConnectivity, TGraph, VertCurv, SingIndex,
                  SingValue),
        PatchM(TriVertPos, TriConnectivity, VertCurv, TGraph,
               PathSampl.SelectedPaths),
        PDeco(PatchM, PathSampl, VertCurv, FaceCurv)

  {

    graph_drift_penalty = 100;
    graph_maxAngle = 45;
    graph_diffuse_step = 1;
    min_sing_distance = 1;

    CCAbility = -1;
    dynamicUpdates = false;
    num_samples = 1000;

    feature_corner_angle = 45;
    QuadSize = 0.02;

    QuadSmoothSteps = 10;
    force_parity = true;
    RegVsIso = 1.0;

    OptimizeSize = true;

    has_problematic_patches = false;
  }

  void Init() {
    Q_LUT.Initialize();
    InitDefaultDualParameters();
    has_problematic_patches = false;
  }

  void ExtractQuad(std::vector<std::pair<int, int>> &_TriFeatures) {
    TriFeatures = _TriFeatures;

    std::sort(TriFeatures.begin(), TriFeatures.end());
    auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
    TriFeatures.erase(last, TriFeatures.end());

    PathSampl.dynamicUpdates = dynamicUpdates;
    PathSampl.NumSamples = num_samples;

    InitializeData();

    DecomposePatches();

    Quadrangulate();

    SolveDoublets();

    SmoothQuadrangulation();
  }
};

#endif