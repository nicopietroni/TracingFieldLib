#include <IO/file_path_functions.h>
#include <cassert>
#include <iostream>
#include <sstream>

// #include <point2.h>
#include "IO/io_features.h"
#include "IO/io_field.h"
#include <IO/io_obj.h>
#include <connectivity.h>
#include <cross_field.h>
#include <field_graph/mesh_preprocess.h>
#include <mesh_features.h>
#include <point3.h>

#include <field_graph/graph.h>
#include <field_graph/graph_tracing.h>
#include <field_graph/initialize_graph.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_managing.h>
#include <field_graph/patch_optimize.h>
#include <field_graph/path_functions.h>
#include <field_graph/path_sampling.h>
#include <field_graph/tracing_interface.h>

// #include <mesh_functions.h>
// #include <cross_field_functions.h>
// #include <smooth_cross_field.h>

#include <field_graph/dual_patches.h>
#include <field_graph/mesh_preprocess.h>
#include <patch_side_discretization.h>

#include <quadrangulation.h>
#include <quadrangulationLUT.h>
#include <vector>

#include <mesh_quadrangulation.h>

// // GLFWwindow *window;
std::string PathMesh;
std::string PathField;
// std::string ProjectName;

typedef double ScalarType;

// MESH DATA
bool loadedMesh = false;
std::vector<Geo::Point3<ScalarType>> TriVertPos;
std::vector<std::vector<int>> TriConnectivity;
std::vector<Geo::Point3<ScalarType>> TriFaceNormals;
std::vector<Geo::Point3<ScalarType>> TriVertNormals;
std::vector<std::vector<int>> NextF, NextE;
Geo::Box3<ScalarType> MeshBox;

// // // FIELD DATA
std::vector<Field::CrossF<ScalarType>> VertCurv;
std::vector<Field::CrossF<ScalarType>> FaceCurv;
std::vector<int> SingIndex, SingValue;

bool has_loaded_mesh = false;
bool has_cross_field = false;
bool has_features = false;
int min_sing_distance = 1;

std::vector<Geo::TJunction> TJunctions;

std::vector<Geo::Point3<ScalarType>> DualVertPos;
std::vector<std::vector<int>> DualConnectivity;
std::vector<Geo::Point3<ScalarType>> DualFaceNormals;
std::vector<Geo::Point3<ScalarType>> DualVertNormals;
std::vector<std::vector<int>> DualNextF, DualNextE;
std::vector<std::vector<int>> PerDualPatchFaceIDX;
std::vector<std::vector<int>> PerDualPatchCornerIDX;
std::vector<std::pair<int, int>> DualFeatures;
std::vector<int> DualCorners;

// bool use_curvature = true;
// bool dualize = false;
// bool smooth_tangent = false;

// bool single_boundary_cond = true;
// bool disk_like_cond = true;
// bool single_sing_cond = false;
// bool match_sing_cond = true;
// bool face_patch_cond = true;
// int MinSides = 3;
// int MaxSides = 6;

// ScalarType minQCrossVert, maxQCrossVert;

// std::vector<int> SingIndex, SingValue;

// Geo::TracingGraph<ScalarType> TGraph(Geo::GT_CrossField);

// Geo::CrossPathSampling<ScalarType> PathSampl(VertPos, Connectivity, TGraph,
// VertCurv, SingIndex, SingValue);

Geo::TracingGraph<ScalarType> TGraph(Geo::GT_CrossField);

float drift_penalty = 100.f;
float maxAngle = 45.f;
int diffuse_step = 1;
int smooth_steps = 10;
// int CurrN = 0;
// int OldN = -1;
// std::vector<int> PathN;

Geo::SideDiscretization<ScalarType> SideD;

Geo::CrossPathSampling<ScalarType> PathSampl(TriVertPos, TriConnectivity,
                                             TGraph, VertCurv, SingIndex,
                                             SingValue);

// std::vector<Geo::Point3<ScalarType>> PathColors;

// std::vector<Geo::Point3<ScalarType>> FaceCol;
// std::vector<Geo::Point3<ScalarType>> FaceColDual;
Geo::PatchManaging<ScalarType> PatchM(TriVertPos, TriConnectivity, VertCurv,
                                      TGraph, PathSampl.SelectedPaths);
// Geo::PatchDecomposer<ScalarType> PDeco(PatchM, PathSampl, VertCurv,
// FaceCurv);
Geo::PatchDecomposer<ScalarType> PDeco(PatchM, PathSampl, VertCurv, FaceCurv);

// ScalarType drift_penalty = 100.f;
// ScalarType maxAngle = 45.f;
// int diffuse_step = 1;
// int smooth_steps = 10;

// bool remove_batch = true;
// bool remove_split = false;

// int Nring = 3;

// Geo::DualPatch<ScalarType>::DualParameter DualParam;
// std::vector<Geo::Point3<ScalarType>> DualVertPos;
// std::vector<std::vector<int>> DualConnectivity;

// SmoothParam<ScalarType> SParam;
// int min_sing_distance = 1;

std::vector<std::pair<int, int>> TriFeatures;
std::vector<int> FeatureSingIndex, FeatureSingValue;
std::vector<std::pair<int, int>> RemainingF;
std::vector<int> NewBorder;

ScalarType Tri_feature_angle = 45;
std::vector<int> TriCornersIdx;

std::vector<int> ModSameBVert;
std::vector<Geo::Point3<ScalarType>> ModSameBOrigPos;

Geo::DualPatch<ScalarType>::DualParameter DualParam;

ScalarType QuadSize = 0.02;

std::vector<Geo::Point3<ScalarType>> QuadVert;
std::vector<std::vector<int>> QuadFaces;
std::vector<int> QuadPatch;
std::vector<Geo::Point3<ScalarType>> QuadFaceNormals, QuadVertNormals;
// std::vector<Geo::Point3<ScalarType>> QuadFaceColors;
std::vector<std::pair<int, int>> QuadFeatures;
std::vector<int> QuadCorners;

Geo::SmoothParam<ScalarType> SParamQuadSmooth;
Geo::QuadrangulationLUT<ScalarType> Q_LUT;
int QuadSmoothSteps = 10;

bool LoadMesh(const std::string &path) {
  //    oss<<"Loading "<<path.c_str()<<std::endl;
  loadedMesh = LoadOBJ(path, TriVertPos, TriConnectivity, false);
  if (loadedMesh) {
    std::cout << "Loaded " << TriConnectivity.size() << " faces "
              << TriVertPos.size() << " vertices" << std::endl;
    MeshBox.Init(TriVertPos);
    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);
    // Geo::getFFAdjacency(Connectivity, NextF, NextE);

    // std::cout<<"FNorm 0 "<<FaceNormals.size()<<std::endl;
  } else {
    std::cout << "Mesh not loaded" << std::endl;
  }
  return loadedMesh;
  // UpdateColor();
}

bool LoadField(const std::string &path) {
  // std::cout<<"Loading Field"<<path.c_str()<<std::endl;
  bool has_loaded =
      ReadField(TriVertPos, TriConnectivity, path, VertCurv, FaceCurv);

  if (!has_loaded)
    return false;

  // if the field is only defined on faces
  if ((VertCurv.size() == 0) || (VertCurv.size() != TriVertPos.size())) {
    assert(FaceCurv.size() == TriConnectivity.size());
    std::cout << "Init Vert Curv From Faces" << std::endl;
    SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormals, FaceCurv,
                         VertCurv);
    // InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);
  }

  if ((FaceCurv.size() == 0) || (FaceCurv.size() != TriConnectivity.size())) {
    assert(VertCurv.size() == TriVertPos.size());
    SetFaceCrossVectorFromVert(TriConnectivity, TriFaceNormals, VertCurv,
                               FaceCurv);

    // InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);
  }
  //   if (has_cross_field) {
  //     std::cout << "Loaded Field" << std::endl;
  //     showCross = true;
  //     showSing = true;

  //     // Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert,
  //     //                                       maxQCrossVert);
  //     Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert,
  //                                           maxQCrossVert);
  //     Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace,
  //                                           maxQCrossFace);

  //     GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE,
  //                          SingIndex, SingValue);

  //     PatchM.SetSingularities(SingIndex, SingValue);
  //     PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
  //     return true;
  //   } else {
  //     std::cout << "Field not loaded" << std::endl;
  //     return false;
  //   }
  return true;
}

// void SaveSubMeshes()
// {
//     std::string subMeshFilename = folderProjectName;
//     subMeshFilename.append("/bordered_layout.obj");
//     std::cout << "Saving on:" << subMeshFilename.c_str() << std::endl;

//     // subMeshFilename.append("_submesh.obj");
//     PDeco.SaveSubMeshes(subMeshFilename);
// }

// void SavePatches()
// {
//     std::string PatchFilename = ProjectName;
//     PatchFilename.append("/patch_layout.patch");
//     PDeco.SavePatchesLayout(PatchFilename);
// }

// void SaveField()
// {
//     std::string FinalFieldFilename = ProjectName;
//     FinalFieldFilename.append("/mesh_layout.ffield");
//     WriteField(FinalFieldFilename.c_str(), VertCurv, FaceCurv);
// }

// void SaveMesh()
// {
//     std::string FinalMFilename = ProjectName;
//     FinalMFilename.append("/mesh_layout.obj");
//     WriteOBJ(FinalMFilename.c_str(), VertPos, Connectivity);
// }

// void SaveDual()
// {
//     if (!dualize)
//         return;
//     std::string FinalMFilename = ProjectName;
//     FinalMFilename.append("/mesh_dual_layout.obj");
//     WriteOBJ(FinalMFilename.c_str(), DualVertPos, DualConnectivity);
// }

// void SaveFPatch()
// {
//     if (!face_patch_cond)
//         return;
//     std::string FPatchFilename = ProjectName;
//     FPatchFilename.append("/patch_layout.fpatch");
//     PDeco.SaveFPatchesLayout(FPatchFilename);

//     std::string CPatchFilename = ProjectName;
//     CPatchFilename.append("/patch_layout.corners");
//     PDeco.SaveCorners(CPatchFilename);
// }

// void SaveAll()
// {
//     SaveSubMeshes();
//     SavePatches();
//     SaveField();
//     SaveMesh();
//     SaveDual();
//     SaveFPatch();
// }

// void CreateOutputFolder()
// {
//     folderProjectName = ProjectName;

//     try
//     {
//         // Check if the folder already exists
//         if (!std::filesystem::exists(folderProjectName))
//         {
//             // Create the folder
//             std::filesystem::create_directory(folderProjectName);
//             std::cout << "Folder created successfully: " << folderProjectName
//             << std::endl;
//         }
//         else
//         {
//             std::cout << "Folder already exists: " << folderProjectName <<
//             std::endl;
//         }
//     }
//     catch (const std::exception &e)
//     {
//         std::cerr << "Error creating folder: " << e.what() << std::endl;
//     }
// }

// void UpdateAfterRemesh()
// {
//     ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
//     Geo::getFFAdjacency(Connectivity, NextF, NextE);

//     InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

//     SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv,
//     VertCurv); InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

//     Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert,
//     maxQCrossVert); GetFaceSingularities(FaceCurv, VertPos, Connectivity,
//     NextF, NextE, SingIndex, SingValue); PatchM.SetSingularities(SingIndex,
//     SingValue);
//     // UpdateColor();
// }

// bool single_boundary_cond = true;
// bool self_connection_cond = true;
// bool disk_like_cond = true;
// bool single_sing_cond = false;
// bool match_sing_cond = true;
// bool cc_ability_cond = false;
// bool equal_side_cond = false;
// bool use_uv = false;
// bool isometric_uv = false;
// bool conformal_uv = false;
// ScalarType max_uvDist = 0.05;

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
  //}
  // if (CCability >= 1) {
  // if (cc_ability_cond) {
  //   // ScalarType AvVAl = AvgEdgeLen(VertPos, Connectivity) * CCability;
  //   PDeco.AddCondition(new Geo::CCabilityCondition<ScalarType>(CCability));
  // }

  // if (isometric_uv) {
  //   Geo::UVCondition<ScalarType>::UVType uv_type =
  //       Geo::UVCondition<ScalarType>::UVArap;
  //   Geo::UVDistortion<ScalarType>::DistType dist_type =
  //       Geo::UVDistortion<ScalarType>::EdgeDist;
  //   PDeco.AddCondition(
  //       new Geo::UVCondition<ScalarType>(uv_type, dist_type, max_uvDist));
  // }
  // if (conformal_uv) {
  //   Geo::UVCondition<ScalarType>::UVType uv_type =
  //       Geo::UVCondition<ScalarType>::UVLSqConformal;
  //   Geo::UVDistortion<ScalarType>::DistType dist_type =
  //       Geo::UVDistortion<ScalarType>::AngleDist;
  //   PDeco.AddCondition(
  //       new Geo::UVCondition<ScalarType>(uv_type, dist_type, max_uvDist));
  // }
  // if (equal_side_cond)
  //   PDeco.AddCondition(new Geo::EqualSideCondition<ScalarType>(EqSideThr));
}

// void InitConditions()
// {
//     PDeco.ClearConditions();
//     if (single_boundary_cond)
//         PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
//     if (disk_like_cond)
//         PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
//     if (single_sing_cond)
//         PDeco.AddCondition(new
//         Geo::SingleSingularityCondition<ScalarType>());
//     if (match_sing_cond)
//         PDeco.AddCondition(new
//         Geo::MathchSingularityCondition<ScalarType>());
//     if (dualize || face_patch_cond)
//         PDeco.AddCondition(new
//         Geo::SelfConnectingPatchCondition<ScalarType>());
//     if ((MinSides > 0) || (MaxSides > 0))
//     {
//         PDeco.AddCondition(new Geo::NumSidesCondition<ScalarType>());
//     }
// }

// void InitPatches()
// {
//     PatchM.Init();
//     InitConditions();
//     // PDeco.InitBorderVertType();
//     PDeco.Init();
//     PDeco.UpdateAllUnsolvedPatches();
// }

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

  // std::cout << "Size Remaining:" << RemainingF.size() << std::endl;
  PDeco.AddFeaturesAsPaths(TGraph, VertCurv, RemainingF, NewBorder);
  // PathColors.resize(PathSampl.SelectedPaths.size());
  // Geo::ScatterColor<ScalarType>(PathColors);
  // UpdateColor();
}

// void InitGraph()
// {
//     Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(TGraph,
//     VertPos, Connectivity, VertCurv,
//                                                                    drift_penalty,
//                                                                    maxAngle,
//                                                                    diffuse_step);
// }

// void SmoothPaths()
// {
//     Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
//     std::vector<std::vector<int>> VertPaths;
//     PathSampl.GetVertexSelectedPaths(VertPaths);

//     for (size_t i = 0; i < VertPaths.size(); i++)
//     {
//         if (VertPaths[i].size() == 0)
//             continue;
//         for (int j = 0; j < (int)VertPaths[i].size() - 1; j++)
//         {
//             int IndexV0 = VertPaths[i][j];
//             int IndexV1 = VertPaths[i][j + 1];
//             std::pair<int, int> key(std::min(IndexV0, IndexV1),
//                                     std::max(IndexV0, IndexV1));
//             UVP.FeatureLines.push_back(key);
//         }
//     }
//     UVP.Steps = smooth_steps;
//     Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity, UVP);

//     //    //update cross F poistion
//     //    for (size_t i=0;i<VertCurv.size();i++)
//     //        VertCurv[i].Pos=VertPos[i];

//     AdaptFieldToMesh(VertPos, Connectivity, VertCurv, FaceCurv);

//     PatchM.UpdateSubPatchPos();
//     PatchM.InitVisualCornersEdges();
// }

// void SmoothFieldGlobal(bool iterative = false)
// {
//     if (!has_cross_field)
//         return;

//     std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;
//     std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;

//     if (!iterative)
//         SmoothGlobalPolyvector(VertPos, Connectivity, FaceNormals, FaceCurv,
//         SParam);
//     else
//         SmoothCrossIterative(VertPos, Connectivity, FaceNormals, FaceCurv,
//         SParam);

//     GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE,
//     SingIndex, SingValue); SetVertCrossFromFace(VertPos, Connectivity,
//     VertNormals, FaceCurv, VertCurv);

//     RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
//     RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

//     InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
//     InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
// }

// void InitFieldByCurvature()
// {
//     SParam.FixedCross.clear();
//     ComputeCurvatureField<ScalarType>(VertPos, Connectivity, FaceNormals,
//     VertNormals, VertCurv, FaceCurv, Nring);
//     Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert,
//     maxQCrossVert); GetFaceSingularities(FaceCurv, VertPos, Connectivity,
//     NextF, NextE, SingIndex, SingValue);

//     has_cross_field = true;
// }

// void BatchProcessFieldCurv()
// {
//     SParam.FixedCross.clear();
//     InitFieldByCurvature();
//     SmoothFieldGlobal();
// }

// void BatchProcessFieldBord()
// {
//     InitFieldByCurvature();
//     SParam.align_borders = true;
//     SParam.TopIntevalFix = 0;
//     SParam.GlobalSmoothVal = -1;
//     SmoothFieldGlobal();
//     std::cout << "There are " << SingIndex.size() << " singularities" <<
//     std::endl;
// }

// void SmoothTessellation()
// {
//     Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity);
//     ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
//     // has_remeshed = true;
// }

// void SplitBorderEdges()
// {
//     RefineMesh<ScalarType>::SpliToHaveOneBorderPerTris(VertPos,
//     Connectivity); Geo::getFFAdjacency(Connectivity, NextF, NextE);
//     ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
// }

// void FindDual()
// {
//     std::vector<std::vector<int>> PerPatchVertIDX;
//     std::vector<std::vector<int>> PerPatchFaceIDX;
//     std::vector<std::vector<int>> PerPatchCornerIDX;
//     Geo::DualPatch<ScalarType>::DisambiguateForDual(PatchM);
//     Geo::DualPatch<ScalarType>::ExtractDual(VertPos, Connectivity, PatchM,
//                                             DualParam, DualVertPos,
//                                             DualConnectivity,
//                                             PerPatchVertIDX, PerPatchFaceIDX,
//                                             PerPatchCornerIDX);
// }

// void UpdateFieldData()
// {

//     PatchM.SetSingularities(SingIndex, SingValue);
//     PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);

// }

// void UpdateAfterRemesh() {
//   ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
//   Geo::getFFAdjacency(Connectivity, NextF, NextE);

//   InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

//   SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv,
//   VertCurv); InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

//   Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert,
//   maxQCrossVert);

//   Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace,
//   maxQCrossFace);
//   // GetFaceSingularities(FaceCurv, VertPos, Connectivity,
//   //                      NextF, NextE, SingIndex, SingValue);

//   PatchM.SetSingularities(SingIndex, SingValue);
//   PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
//   UpdateColor();
// }

void UpdateAfterRemesh() {
  ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);
  // ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(TriConnectivity, NextF, NextE);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

  SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormals, FaceCurv,
                       VertCurv);
  // InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

  // Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert,
  // maxQCrossVert);

  // Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace,
  // maxQCrossFace);
  // // GetFaceSingularities(FaceCurv, VertPos, Connectivity,
  //                      NextF, NextE, SingIndex, SingValue);

  GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                       SingIndex, SingValue);

  PatchM.SetSingularities(SingIndex, SingValue);
  PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
  // UpdateColor();
}

void PreprocessMesh() {
  Geo::getFFAdjacency(TriConnectivity, NextF, NextE);

  GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                       SingIndex, SingValue);

  Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(
      TriVertPos, TriConnectivity, FaceCurv, TriFeatures, SingIndex,
      min_sing_distance);

  // MeshBox.Init(VertPos);
  UpdateAfterRemesh();

  // WriteOBJ("after_remesh.obj",TriVertPos, TriConnectivity);

  std::vector<Field::CrossF<ScalarType>> VertCurv0;
  // WriteField("after_remesh.ffield",VertCurv0,FaceCurv);
  // WriteFeatures("after_remesh.features",TriFeatures);
  // exit(0);

  std::cout << "*** SPLIT FEATURES *** " << std::endl;

  RemainingF.clear();
  Geo::TracingPreprocess<ScalarType>::SplitCornerFeatures(
      TriVertPos, TriConnectivity, VertCurv, FaceCurv, TriFeatures, RemainingF,
      SingIndex, SingValue, FeatureSingIndex, FeatureSingValue, NewBorder);

  // Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert,
  // maxQCrossVert); Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv,
  // minQCrossFace, maxQCrossFace);

  Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
      TriVertPos, TriConnectivity, &ModSameBVert, &ModSameBOrigPos);

  UpdateAfterRemesh();

  WriteOBJ("testBeforeGraph.obj", TriVertPos, TriConnectivity);

  // exit(0);
}

void InitGraph() {
  Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(
      TGraph, TriVertPos, TriConnectivity, VertCurv, drift_penalty, maxAngle,
      diffuse_step);

  Geo::PathFunctions<ScalarType>::DisableNonDirectableConnections(TGraph);
}

void SmoothPaths() {
  Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  std::vector<std::vector<int>> VertPaths;
  PathSampl.GetVertexSelectedPaths(VertPaths);
  Geo::PathFunctions<ScalarType>::FindTJunctions(VertPaths, TJunctions);

  // Geo::PatchOptimize<ScalarType>::TwoStepSmooth(VertPos, Connectivity,
  //                                               VertPaths);

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

void ExtractDualPatches() {
  RestoreOriginalPosOnSplitBorders();
  Geo::DualPatch<ScalarType>::DisambiguateForDual(PatchM);

  // Geo::DualPatch<ScalarType>::ExtractDual(
  //     VertPos, Connectivity, PatchM, DualParam, DualVertPos,
  //     DualConnectivity, PerDualPatchFaceIDX, PerDualPatchCornerIDX);

  Geo::DualPatch<ScalarType>::PrimalPatchData PData;
  PData.PatchFaces = PatchM.PData.SubPatchFaces;
  PData.PatchVertPos = PatchM.PData.SubPatchVertPos;
  PData.PatchCorners = PatchM.PData.PatchCorners;
  PData.EdgePatchSides = PatchM.PData.EdgePatchSides;
  for (size_t i = 0; i < PatchM.NumPatches(); i++) {
    std::vector<int> PatchSing;
    PatchM.PatchSingularities(i, PatchSing);
    PData.PatchSingularities.push_back(PatchSing);
  }

  Geo::DualPatch<ScalarType>::ExtractDual(PData, DualParam, DualVertPos,
                                          DualConnectivity, PerDualPatchFaceIDX,
                                          PerDualPatchCornerIDX);

  // std::cout << "Size F:" << PerDualPatchFaceIDX.size() << std::endl;
  // std::cout << "Vert FIdx:" << PerDualPatchVertIDX.size() << std::endl;
  // std::cout << "Corner FIdx:" << PerDualPatchCornerIDX.size() << std::endl;
  // // exit(0);

  // // UPDATE COLOR
  // FaceColDual.clear();
  // FaceColDual.resize(DualConnectivity.size(),
  //                    Geo::Point3<ScalarType>(0.5, 1, 0.83));

  // for (size_t i = 0; i < PerDualPatchFaceIDX.size(); i++) {
  //   Geo::Point3<ScalarType> CurrCol =
  //       Geo::ScatterColor<ScalarType>(PerDualPatchFaceIDX.size(), i);

  //   if (PerDualPatchCornerIDX[i].size() != 4)
  //     CurrCol = Geo::Point3<ScalarType>(1, 0, 0);

  //   for (size_t j = 0; j < PerDualPatchFaceIDX[i].size(); j++) {
  //     int IndexF = PerDualPatchFaceIDX[i][j];
  //     FaceColDual[IndexF] = CurrCol;
  //   }
  // }
  // // UPDATE CORNERS
  // Geo::DualPatch<ScalarType>::ExtractCornerVisPos(
  //     DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
  //     PerDualPatchCornerIDX, VisDualCornerPos);

  std::vector<int> PrimalCorners;
  Geo::DualPatch<ScalarType>::RetrieveDualFeatures(
      TriVertPos, TriConnectivity, TriFeatures, PrimalCorners, DualVertPos,
      DualConnectivity, DualFeatures, DualCorners);

  // has_dual_mesh = true;

  // UPDATE MESH DATA
  ComputeNormals(DualVertPos, DualConnectivity, DualFaceNormals,
                 DualVertNormals);

  Geo::getFFAdjacency(DualConnectivity, DualNextF, DualNextE);

  // CONTROLS FOR VISUALIZATION
  // showMesh = false;
  // showDual = true;
  // DrawPathMode = 2;

  // INITIALIZE SIDE SUBDIVISION
  SideD.InitFromLabels(DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
                       PerDualPatchCornerIDX, DualFeatures);
}

void InitDefaultDualParameters() {
  DualParam.OnlyNonQuad = true;
  DualParam.InitSingPos = true;
  DualParam.FixSingPos = false;
  DualParam.JointBorder = true;
}

bool force_parity = true;
ScalarType RegVsIso = 1.0;

// void InitTangentSmooth() {
//   std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
//       FeatureData;
//   SideD.GetFeaturesQuadrData(DualVertPos, DualFeatures, FeatureData);

//   bool has_mapped = FindEdgePosCorrespondences(QuadVert, QuadFaces,
//   FeatureData,
//                                                UVP.FeatureLines);
//   assert(has_mapped);

//   std::set<Geo::Point3<ScalarType>> CornerPosSet;
//   for (size_t i = 0; i < CornersIdx.size(); i++)
//     CornerPosSet.insert(VertPos[CornersIdx[i]]);

//   for (size_t i = 0; i < QuadVert.size(); i++) {
//     if (CornerPosSet.count(QuadVert[i]) > 0)
//       UVP.Fixed.push_back(i);
//   }
//   UVP.Damp = 0.5;
//   UVP.Steps = 10;
//   UVP.fixBorder = true;
//   UVP.post_project = true;
//   UVP.ReprojVert = &VertPos;
//   UVP.ReprojFaces = &Connectivity;
//   UVP.ReprojFeatureLines = Features;
// }

void UpdateQuadFeatures() {
  std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
      FeatureData;
  SideD.GetFeaturesQuadrData(DualVertPos, DualFeatures, FeatureData);

  QuadFeatures.clear();
  bool has_mapped = FindEdgePosCorrespondences(QuadVert, QuadFaces, FeatureData,
                                               QuadFeatures);
  assert(has_mapped);

  // find the corners
  QuadCorners.clear();
  std::set<Geo::Point3<ScalarType>> CornerPosSet;
  for (size_t i = 0; i < TriCornersIdx.size(); i++)
    CornerPosSet.insert(TriVertPos[TriCornersIdx[i]]);

  CornersByCount(QuadVert, QuadFeatures, QuadCorners);

  for (size_t i = 0; i < QuadVert.size(); i++) {
    if (CornerPosSet.count(QuadVert[i]) > 0)
      QuadCorners.push_back(i);
  }
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

void SmoothCentroidMesh3D() {

  SParamQuadSmooth.NumIte = QuadSmoothSteps;
  SParamQuadSmooth.SType = Geo::SMCentroid;
  SParamQuadSmooth.check_AspectR = true;
  Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

  Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                  QuadVertNormals);
  
  WriteOBJ("quad_mesh.obj",QuadVert,QuadFaces);
}

void BatchProcessQuad() {
  std::cout << "*** REFINE AND SPLIT *** " << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  PreprocessMesh();
  auto end0 = std::chrono::high_resolution_clock::now();
  auto elapsed0 =
      std::chrono::duration_cast<std::chrono::microseconds>(end0 - start);

  std::cout << "*** INIT GRAPHS *** " << std::endl;
  InitGraph();
  // auto end1 = std::chrono::high_resolution_clock::now();
  // auto elapsed1 =
  //     std::chrono::duration_cast<std::chrono::microseconds>(end1 - end0);

  // std::cout << "*** INIT PATCHES *** " << std::endl;
  InitPatches();
  // auto end2 = std::chrono::high_resolution_clock::now();
  // auto elapsed2 =
  //     std::chrono::duration_cast<std::chrono::microseconds>(end2 - end1);

  // std::cout << "*** ADD ALL PATHS *** " << std::endl;
  // // PDeco.dynamic_smooth=true;

  bool remove_batch = true;
  PDeco.add_only_where_needed = true;
  PDeco.stop_when_solved = true;
  // // PDeco.dynamic_smooth = false;
  // //  PDeco.SolveAllFeatures();
  PDeco.BatchDecompose(remove_batch);
  // // SmoothPaths();
  // // PDeco.SolveAllFeatures();
  // // SmoothPaths();
  // // PDeco.SolveAllMesh();

  SmoothPaths();

  // auto end3 = std::chrono::high_resolution_clock::now();
  // auto elapsed3 =
  //     std::chrono::duration_cast<std::chrono::microseconds>(end3 - end2);

  RestoreOriginalPosOnSplitBorders();

  // std::cout << "*** SMOOTH PATCHES *** " << std::endl;
  SmoothPaths();
  // auto end4 = std::chrono::high_resolution_clock::now();
  // auto elapsed4 =
  //     std::chrono::duration_cast<std::chrono::microseconds>(end4 - end3);

  // double elapsed_seconds0 =
  //     elapsed0.count() * 1e-6; // Convert microseconds to seconds
  // double elapsed_seconds1 =
  //     elapsed1.count() * 1e-6; // Convert microseconds to seconds
  // double elapsed_seconds2 =
  //     elapsed2.count() * 1e-6; // Convert microseconds to seconds
  // double elapsed_seconds3 =
  //     elapsed3.count() * 1e-6; // Convert microseconds to seconds
  // double elapsed_seconds4 =
  //     elapsed4.count() * 1e-6; // Convert microseconds to seconds

  // showCross = false;

  // auto total_elapsed =
  //     std::chrono::duration_cast<std::chrono::microseconds>(end4 - start);

  // double total_elapsed_seconds5 =
  //     total_elapsed.count() * 1e-6; // Convert microseconds to seconds

  // // std::cout << std::endl << std::endl;
  // // std::cout << "*** FINAL TIMINGS (sec) *** " << std::endl;
  // // std::cout << "* Remeshing time: " << elapsed_seconds0 << std::endl;
  // // std::cout << "* Init Graph Elapsed time: " << elapsed_seconds1 <<
  // // std::endl; std::cout << "* Init Patches: " << elapsed_seconds2 <<
  // // std::endl; std::cout << "* Adding Paths: " << elapsed_seconds3 <<
  // // std::endl; std::cout << "* Smooth Patches: " << elapsed_seconds4 <<
  // // std::endl; std::cout << "* TOTAL TIME: " << total_elapsed_seconds5 <<
  // // std::endl;

  // std::cout << std::endl << std::endl;
  // std::cout << "*** STATISTICS *** " << std::endl;
  // std::cout << "*Primal* " << std::endl;
  PatchM.WriteStats();
  PDeco.WriteStats();

  ExtractDualPatches();

  ScalarType EdgeS = MeshBox.Diag() * QuadSize;

  Geo::SideDiscretization<ScalarType>::solveMode SM;

  SM = Geo::SideDiscretization<ScalarType>::SM_GlobalSolve;

  SideD.Solve(EdgeS, SM, force_parity, RegVsIso);

  std::vector<std::vector<int>> PatchFaces, PatchCorners;
  std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> SideSamples;

  // std::vector<std::pair<Geo::Point3<ScalarType>,
  // Geo::Point3<ScalarType>>>
  //     FeatureData;
  // SideD.GetFeaturesQuadrData(DualVertPos, DualFeatures, FeatureData);
  SideD.GetQuadrangulationData(PatchFaces, PatchCorners, SideSamples);

  Geo::Quadrangulation<ScalarType>::QuadrangulateQuadPatchsWithSubs(
      DualVertPos, DualConnectivity, PatchFaces, PatchCorners, SideSamples,
      Q_LUT, QuadVert, QuadFaces, QuadPatch);

  RemoveDuplicatedVert(QuadVert, QuadFaces);
  RemoveUnfererencedVert(QuadVert, QuadFaces);

  // SmoothQuadMesh();

  // ElemColorByCluster<ScalarType>(QuadPatch, QuadFaceColors);

  // QuadFaceColors = std::vector<Geo::Point3<ScalarType>>(
  //     QuadFaces.size(), Geo::Point3<ScalarType>(0.8, 0.8, 0.8));

  Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                  QuadVertNormals);

  // InitTangentSmooth();
  UpdateQuadFeatures();
  Init3DSmooth();

  SmoothCentroidMesh3D();
}

// void BatchProcess()
// {
//     auto t0 = std::chrono::high_resolution_clock::now();
//     if (!has_cross_field)
//     {
//         if (smooth_tangent)
//             SmoothTessellation();

//         std::cout << "*** COMPUTING FIELD *** " << std::endl;
//         SplitBorderEdges();

//         if (use_curvature)
//             std::cout << "usign curvature " << std::endl;
//         else
//             std::cout << "using border " << std::endl;
//         if (use_curvature)
//             BatchProcessFieldCurv();
//         else
//             BatchProcessFieldBord();
//     }
//     std::cout << "*** REFINE STEP *** " << std::endl;
//     Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(VertPos,
//     Connectivity,
//                                                              FaceCurv,
//                                                              SingIndex,
//                                                              min_sing_distance);
//     MeshBox.Init(VertPos);
//     UpdateAfterRemesh();
//     auto t1 = std::chrono::high_resolution_clock::now();
//     auto elapsed0 = std::chrono::duration_cast<std::chrono::microseconds>(t1
//     - t0);

//     std::cout << "*** INIT GRAPHS *** " << std::endl;
//     InitGraph();
//     auto t2 = std::chrono::high_resolution_clock::now();
//     auto elapsed1 = std::chrono::duration_cast<std::chrono::microseconds>(t2
//     - t1);

//     std::cout << "*** INIT PATCHES *** " << std::endl;
//     InitPatches();
//     auto t3 = std::chrono::high_resolution_clock::now();
//     auto elapsed2 = std::chrono::duration_cast<std::chrono::microseconds>(t3
//     - t2);

//     std::cout << "*** ADD ALL PATHS *** " << std::endl;
//     PDeco.BatchDecompose(remove_batch, remove_split);
//     auto t4 = std::chrono::high_resolution_clock::now();
//     auto elapsed3 = std::chrono::duration_cast<std::chrono::microseconds>(t4
//     - t3);

//     std::cout << "*** SMOOTH PATCHES *** " << std::endl;
//     SmoothPaths();
//     auto t5 = std::chrono::high_resolution_clock::now();
//     auto elapsed4 = std::chrono::duration_cast<std::chrono::microseconds>(t5
//     - t4);

//     if (dualize)
//     {
//         std::cout << "*** DUALIZE *** " << std::endl;
//         FindDual();
//     }
//     auto t6 = std::chrono::high_resolution_clock::now();
//     auto elapsed5 = std::chrono::duration_cast<std::chrono::microseconds>(t6
//     - t5);

//     // Convert microseconds to seconds
//     double elapsed_seconds0 = elapsed0.count() * 1e-6;
//     double elapsed_seconds1 = elapsed1.count() * 1e-6;
//     double elapsed_seconds2 = elapsed2.count() * 1e-6;
//     double elapsed_seconds3 = elapsed3.count() * 1e-6;
//     double elapsed_seconds4 = elapsed4.count() * 1e-6;
//     double elapsed_seconds5 = elapsed4.count() * 1e-6;

//     std::cout << std::endl<< std::endl;
//     std::cout << "*** TIMINGS (sec) *** " << std::endl;
//     std::cout << "* Field & Remd time: " << elapsed_seconds0 << std::endl;
//     std::cout << "* Init Graph Elapsed time: " << elapsed_seconds1 <<
//     std::endl; std::cout << "* Init Patches: " << elapsed_seconds2 <<
//     std::endl; std::cout << "* Adding Paths: " << elapsed_seconds3 <<
//     std::endl; std::cout << "* Smooth Patches: " << elapsed_seconds4 <<
//     std::endl; std::cout << "* Dualize Patches: " << elapsed_seconds5 <<
//     std::endl;

//     std::cout << std::endl<< std::endl;
//     std::cout << "*** STATISTICS *** " << std::endl;
//     std::cout << "*Primal* " << std::endl;
//     PatchM.WriteStats();
//     PDeco.WriteStats();
// }

void CleanMeshArtifacts() {
  ScalarType AVLen = AvgEdgeLen(TriVertPos, TriConnectivity);
  Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
      TriVertPos, TriConnectivity);
}

std::string folderProjectName;
std::string commandLine;

int main(int argc, char *argv[]) {
  std::stringstream commandLineStr;
  for (int i = 0; i < argc; ++i) {
    commandLineStr << argv[i] << " ";
  }
  commandLine = commandLineStr.str();

  // // oss<<"Initialized Log Window"<<std::endl;
  // // AScene.SetOutputStream(&oss);

  if (argc < 2) {
    std::cout << "You should pass at least a mesh" << std::endl;
    exit(0);
  }

  for (int i = 1; i < argc; i++) {
    if (GetFileExtension(std::string(argv[i])) == std::string("obj")) {
      PathMesh = std::string(argv[1]);
      has_loaded_mesh = LoadMesh(PathMesh);
      if (!has_loaded_mesh) {
        std::cout << "ERROR LOADING INPUT MESH" << std::endl;
        exit(0);
      }
      CleanMeshArtifacts();
      continue;
    }
    if ((GetFileExtension(std::string(argv[i])) == std::string("ffield")) ||
        (GetFileExtension(std::string(argv[i])) == std::string("field"))) {
      PathField = std::string(argv[i]);
      std::cout << "Loading Field" << PathField.c_str() << std::endl;
      has_cross_field = LoadField(PathField);
      if (!has_cross_field) {
        std::cout << "ERROR LOADING FIELD" << std::endl;
        exit(0);
      }
      continue;
    }
    if (GetFileExtension(std::string(argv[i])) == std::string("features")) {
      std::string PathFeat = std::string(argv[i]);
      std::cout << "Loading Features" << std::endl;

      bool load_feat = ReadFeatures(PathFeat, TriFeatures, false);
      if (!load_feat) {
        std::cout << "ERROR LOADING FEATURES" << std::endl;
        exit(0);
      }

      Geo::Corners<ScalarType>(TriVertPos, TriFeatures, Tri_feature_angle,
                               TriCornersIdx);
      has_features = true;
      continue;
    }

    // if (std::string(argv[i]) == std::string("fpatch"))
    // {
    //     face_patch_cond = true;
    // }
    // if (std::string(argv[i]) == std::string("dual"))
    // {
    //     dualize = true;
    // }
    // if (std::string(argv[i]) == std::string("even"))
    // {
    //     DualParam.RType = Geo::DualPatch<ScalarType>::RTypeEven;
    // }
    // if (std::string(argv[i]) == std::string("ortho"))
    // {
    //     DualParam.RType = Geo::DualPatch<ScalarType>::RTypeOrtho;
    // }
    // if (std::string(argv[i]) == std::string("no_split_q"))
    // {
    //     DualParam.OnlyNonQuad = true;
    // }
    // if (std::string(argv[i]) == std::string("border"))
    // {
    //     use_curvature = false;
    // }
    // if (std::string(argv[i]) == std::string("smoothtangent"))
    // {
    //     smooth_tangent = true;
    // }
    // if (std::string(argv[i]) == std::string("no_match_sing_val"))
    // {
    //     single_sing_cond = true;
    //     match_sing_cond = false;
    // }
  }
  // if (has_features && has_cross_field && has_loaded_mesh) {
  //   std::cout << "BATCH PROCESSING" << std::endl;
  //   Q_LUT.Initialize();
  //   InitDefaultDualParameters();
  //   BatchProcessQuad();
  // }
  //if (has_features && has_cross_field && has_loaded_mesh) {
  if (has_cross_field && has_loaded_mesh) {
    std::cout << "BATCH PROCESSING" << std::endl;
    MeshQuadrangulator<ScalarType> MQuad(TriVertPos,TriConnectivity,
                                         VertCurv,FaceCurv,
                                        QuadVert,QuadFaces);
    MQuad.Init();
    MQuad.ExtractQuad(TriFeatures);
    WriteOBJ("quad_mesh.obj",QuadVert,QuadFaces);

    // Q_LUT.Initialize();
    // InitDefaultDualParameters();
    // BatchProcessQuad();
  }
  // ProjectName = RemoveExtension(PathMesh);

  // CreateOutputFolder();

  // SaveAll();

  return 0;
}
