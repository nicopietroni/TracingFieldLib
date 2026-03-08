// //Non-Commercial License Agreement
// //Copyright (c) 2025, Nico Pietroni
// //
// //1. Permission is hereby granted to use, copy, modify, and distribute
// //   this software for academic research, teaching, and personal
// //   non-commercial purposes, provided that this copyright notice is
// //   included in all copies or substantial portions of the software.
// //
// //2. Commercial use of this software, including but not limited to use in
// //   commercial products, services, tools, or consulting, is strictly
// //   prohibited without prior written permission from the copyright holder.
// //
// //3. This software is provided "as is", without warranty of any kind,
// //   express or implied, including but not limited to the warranties of
// //   merchantability, fitness for a particular purpose, and
// //   noninfringement. In no event shall the authors or copyright holders
// //   be liable for any claim, damages, or other liability, arising from,
// //   out of, or in connection with the software or the use or other
// //   dealings in the software.
// //
// //4. For commercial licensing, please contact:
// //   Nico Pietroni
// //   Email: nico.pietroni@gmail.com

#include <IO/file_path_functions.h>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <utility>
#include <vector>

// #include <point2.h>
#include "IO/io_features.h"
#include "IO/io_field.h"
#include "mesh_functions.h"
#include "plane3.h"
#include <IO/io_obj.h>
#include <point3.h>

#include <field_graph/refine_for_tracing.h>
#include <init_cross_field.h>
#include <mesh_mapping.h>
#include <mesh_quadrangulation.h>
#include <remesh_size.h>
#include <triangular_remesh.h>
#include "mesh_symmetry.h"

std::string PathMesh;
std::string PathField;
std::string ProjectName;
std::string Outfolder;

typedef double ScalarType;

bool extended_print = false;
// MESH DATA
bool loadedMesh = false;
std::vector<Geo::Point3<ScalarType>> TriVertPos;
std::vector<std::vector<int>> TriConnectivity;
std::vector<std::pair<int, int>> TriFeatures;
std::vector<int> TriCorners;

std::vector<Geo::Point3<ScalarType>> TriVertNormal, TriFaceNormal;

std::vector<std::vector<int>> NextF, NextE;

// // // FIELD DATA
std::vector<Field::CrossF<ScalarType>> FaceCurv, VertCurv;
FieldSmoothParam<ScalarType> SParam;
int Nring = 3;
ScalarType minQCrossVert, maxQCrossVert;
std::vector<int> SingIndex, SingValue;

bool has_loaded_mesh = false;
bool has_cross_field = false;
bool has_features = false;
bool to_split_symmetry = false;

ScalarType QuadSize;
std::vector<Geo::Point3<ScalarType>> QuadVert;
std::vector<std::vector<int>> QuadFaces;

ScalarType SizeMult = 1.0;

Geo::Plane3<ScalarType> SymmetryPlane;

bool LoadMesh(const std::string &path) {
  //    oss<<"Loading "<<path.c_str()<<std::endl;
  loadedMesh = LoadOBJ(path, TriVertPos, TriConnectivity, false);
  if (loadedMesh) {
    std::cout << "Loaded " << TriConnectivity.size() << " faces "
              << TriVertPos.size() << " vertices" << std::endl;
  } else {
    std::cout << "Mesh not loaded" << std::endl;
  }
  return loadedMesh;
}

bool LoadField(const std::string &path) {
  std::vector<Field::CrossF<ScalarType>> VertCurv;

  bool has_loaded =
      ReadField(TriVertPos, TriConnectivity, path, VertCurv, FaceCurv);

  if (!has_loaded)
    return false;

  return true;
}

void InitFieldByCurvature() {
  SParam.FixedCross.clear();
  ComputeCurvatureField<ScalarType>(TriVertPos, TriConnectivity, TriFaceNormal,
                                    TriVertNormal, VertCurv, FaceCurv, Nring);

  Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert, maxQCrossVert);
  GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                       SingIndex, SingValue);

  has_cross_field = true;
}

// void SaveMesh() {
//   assert(!PathMesh.empty());
//   std::string PathRemMesh = RemoveExtension(PathMesh);
//   PathRemMesh += std::string("_remesh.obj");
//   std::cout << "Saving Field on file " << PathRemMesh.c_str() << std::endl;
//   bool Saved = WriteOBJ(PathRemMesh, TriVertPos, TriConnectivity);
//   if (!Saved) {
//     std::cout << "ERROR SAVING FIELD FILE " << std::endl;
//     exit(0);
//   }
// }

void SaveField() {
  assert(!PathMesh.empty());
  std::string PathField = RemoveExtension(PathMesh);
  PathField += std::string("_remesh.field");
  std::cout << "Saving Field on file " << PathField.c_str() << std::endl;
  bool Saved = WriteField(PathField, VertCurv, FaceCurv);
  if (!Saved) {
    std::cout << "ERROR SAVING FIELD FILE " << std::endl;
    exit(0);
  }
  // ReadField(VertPos, Connectivity, PathField, VertCurv, FaceCurv);
}

void SaveFeatures() {
  if (TriFeatures.size() == 0)
    return;
  assert(!PathMesh.empty());
  std::string PathFeatures = RemoveExtension(PathMesh);
  PathFeatures += std::string("_remesh.features");
  std::cout << "Saving Features on file " << PathFeatures.c_str() << std::endl;
  bool Saved = WriteFeatures(PathFeatures, TriFeatures);
  if (!Saved) {
    std::cout << "ERROR SAVING FIELD FILE " << std::endl;
    exit(0);
  }
}

void SmoothGlobal(bool iterative = false) {
  if (!has_cross_field)
    return;

  std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;
  std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;

  SParam.Features = TriFeatures;

  if (!iterative)
    SmoothGlobalPolyvector(TriVertPos, TriConnectivity, TriFaceNormal, FaceCurv,
                           SParam);
  // else
  //   SmoothCrossIterative(VertPos, Connectivity, FaceNormals, FaceCurv,
  //   SParam);

  GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                       SingIndex, SingValue);

  ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormal, TriVertNormal);

  SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormal, FaceCurv,
                       VertCurv);

  RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
  RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
}

void RefineForTracing() {
  Geo::RefineForTracing<ScalarType>::RefineForFieldComputation(
      TriVertPos, TriConnectivity, TriFeatures);
  Geo::getFFAdjacency(TriConnectivity, NextF, NextE);
  ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormal, TriVertNormal);
  // has_remeshed = true;
}

void RefineForTracingConnectivity() {
  Geo::RefineForTracing<ScalarType>::RefineForConnection(
      TriVertPos, TriConnectivity, FaceCurv, TriFeatures);
  Geo::getFFAdjacency(TriConnectivity, NextF, NextE);
  ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormal, TriVertNormal);
}

void BatchProcessField() {

  RefineForTracing();
  SParam.FixedCross.clear();
  InitFieldByCurvature();

  SmoothGlobal();

  RefineForTracingConnectivity();

  for (size_t i = 0; i < TriFeatures.size(); i++) {
    if (TriFeatures[i].first > TriFeatures[i].second)
      std::swap(TriFeatures[i].first, TriFeatures[i].second);
  }
  std::sort(TriFeatures.begin(), TriFeatures.end());
  auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
  TriFeatures.erase(last, TriFeatures.end());
}

ScalarType MinSize = 0.005;
ScalarType MaxSize = 0.05;

Geo::Box3<ScalarType> TriMeshBox;

ScalarType initial_feature_angle = 60;
int ero_steps = 3;
int dil_steps = 3;
bool compute_features = true;
bool only_border_corners = false;
bool add_features = false;
int smooth_steps = 15;
std::string folderProjectName;
std::string commandLine;

bool remesh_near_symmetry = true;

// to be consistent with edge size always
// compute as it use 30 degrees features
void AutoQuadEdgeSize() {

  std::vector<std::pair<int, int>> TriFeatures_30;
  std::vector<int> TriCorners_30;

  if (has_features)
    Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
        TriVertPos, TriConnectivity, 30, TriFeatures_30, TriCorners_30,
        ero_steps, dil_steps);

  // std::vector<ScalarType> VertSize;
  ScalarType MinAggrSize = TriMeshBox.Diag() * MinSize;
  ScalarType MaxAggrSize = TriMeshBox.Diag() * MaxSize;
  ScalarType CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeAggregated(
      TriVertPos, TriConnectivity, TriFeatures_30, MinAggrSize, MaxAggrSize,
      NULL);
  // ScalarType CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeFeatures(
  //     TriVertPos, TriConnectivity, TriFeatures_30, MinAggrSize, MaxAggrSize,
  //     VertSize);
  QuadSize = CurrSize / TriMeshBox.Diag();
  QuadSize /= 2;
  QuadSize *= SizeMult;

  if (extended_print)
    std::cout << "Quad Size: " << QuadSize << std::endl;
}

void CreateOutputFolder() {

  Outfolder = GetDirectory(PathMesh);
  Outfolder += std::string("/quadrangulated/");

  try {
    // Check if the folder already exists
    if (!std::filesystem::exists(Outfolder))
      std::filesystem::create_directory(Outfolder);
  } catch (const std::exception &e) {
    std::cerr << "Error creating folder: " << e.what() << std::endl;
  }
  // exit(0);
}

void SaveMesh() {
  std::string RemFile = Outfolder;
  RemFile += std::string(GetFilanameNoExtension(PathMesh));
  RemFile.append("_quad.obj");
  std::cout << std::endl;
  std::cout << "Saving on:" << RemFile.c_str() << std::endl;
  WriteOBJ(RemFile, QuadVert, QuadFaces);
}

Geo::Point3<ScalarType> DirSymm(1, 0, 0);

// void SplitSymmetry() {
//   Geo::Point3<ScalarType> Center(0, 0, 0); // MeshBox.Center();
//   SymmetryPlane.Init(Center, DirSymm);

//   std::vector<Geo::Point3<ScalarType>> OldVertPos = TriVertPos;
//   std::vector<std::vector<int>> OldConnectivity = TriConnectivity;
//   RefineMesh<ScalarType>::SplitMeshWithPlane(TriVertPos, TriConnectivity,
//                                              SymmetryPlane);
//   RemoveFacesOnNegativeSide<ScalarType>(TriVertPos, TriConnectivity,
//                                         SymmetryPlane);

//   if (has_features) {
//     std::vector<int> VerteMap;
//     FindVertCorrespondences<ScalarType>(OldVertPos, TriVertPos, VerteMap);
//     // remap the features
//     std::vector<std::pair<int, int>> NewFeatures;
//     for (size_t i = 0; i < TriFeatures.size(); i++) {
//       int IndexV0 = TriFeatures[i].first;
//       int IndexV1 = TriFeatures[i].second;
//       assert(IndexV0 >= 0);
//       assert(IndexV1 >= 0);
//       assert(IndexV0 < VerteMap.size());
//       assert(IndexV1 < VerteMap.size());
//       int NewIndexV0 = VerteMap[IndexV0];
//       int NewIndexV1 = VerteMap[IndexV1];
//       if (NewIndexV0 == -1 || NewIndexV1 == -1) {
//         continue; // this feature is lost as one of the vertices is on the
//                   // symmetry plane
//       }
//       NewFeatures.push_back(std::pair<int, int>(
//           std::min(NewIndexV0, NewIndexV1), std::max(NewIndexV0, NewIndexV1)));
//     }
//     TriFeatures = NewFeatures;
//     std::sort(TriFeatures.begin(), TriFeatures.end());
//     auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
//     TriFeatures.erase(last, TriFeatures.end());
//   }

//   Geo::TriRemParam<ScalarType> RemP;
//   RemP.TargetL = AvgEdgeLen(TriVertPos, TriConnectivity); // QuadSize;
//   RemP.feature_angle = initial_feature_angle;
//   RemP.collapse_small = true;
//   RemP.refine_big = true;
//   RemP.flip_valence = true;
//   RemP.flip_quadrics = true;
//   RemP.smooth_mesh = true;
//   std::cout << "Remeshing near symmetry plane, target size: " << RemP.TargetL
//             << std::endl;
//   if (has_features)
//     Geo::MeshFeatures<ScalarType>::BoolVariablesFromPairs(
//         TriVertPos, TriConnectivity, TriFeatures, TriCorners,
//         RemP.IsFaceEFeature, RemP.IsVertCorner);
  
//   WriteOBJ("split_symmetry_before.obj", TriVertPos, TriConnectivity);
//   if (remesh_near_symmetry)
//     Geo::TriRemeshCloseToPlane<ScalarType>(TriVertPos, TriConnectivity,
//                                            SymmetryPlane, 0.05, RemP);
//   WriteOBJ("split_symmetry.obj", TriVertPos, TriConnectivity);

//   if (has_features) {
//     Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(
//         RemP.IsVertCorner, TriCorners);
//     Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
//         TriConnectivity, RemP.IsFaceEFeature, TriFeatures);
//   }

//   if (has_features)
//     Geo::MeshFeatures<ScalarType>::BoolVariablesFromPairs(
//         TriVertPos, TriConnectivity, TriFeatures, TriCorners,
//         RemP.IsFaceEFeature, RemP.IsVertCorner);

//   Geo::TriRemeshCloseToPlane<ScalarType>(TriVertPos, TriConnectivity,
//                                          SymmetryPlane, 0.025, RemP);

//   if (has_features) {
//     Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(
//         RemP.IsVertCorner, TriCorners);
//     Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
//         TriConnectivity, RemP.IsFaceEFeature, TriFeatures);
//   }

//   Geo::getFFAdjacency(TriConnectivity, NextF, NextE);
//   ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormal, TriVertNormal);

//   std::cout << "Mesh splitted with symmetry plane, new num faces: "
//             << TriConnectivity.size() << std::endl;
// }


void SplitSymmetry() {
  Geo::Point3<ScalarType> Center(0, 0, 0); // MeshBox.Center();
  SymmetryPlane.Init(Center, DirSymm);

  ScalarType Tolerance = AvgEdgeLen(TriVertPos, TriConnectivity) * 0.1;
  Geo::MeshSymmetry<ScalarType>::CutMeshWithSymmetryPlane(TriVertPos, TriConnectivity,
                                              TriFeatures, TriCorners,
                                              SymmetryPlane, Tolerance); 

  Geo::getFFAdjacency(TriConnectivity, NextF, NextE);
  ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormal, TriVertNormal);
  // WriteOBJ("split_symmetry.obj", TriVertPos, TriConnectivity);
  // exit(0);
  // std::cout << "Mesh splitted with symmetry plane, new num faces: "
  //           << TriConnectivity.size() << std::endl;
}

int main(int argc, char *argv[]) {
  std::stringstream commandLineStr;

  MeshQuadrangulator<ScalarType> MQuad(TriVertPos, TriConnectivity, FaceCurv,
                                       QuadVert, QuadFaces);
  MQuad.Init();

  for (int i = 0; i < argc; ++i) {
    commandLineStr << argv[i] << " ";
  }
  commandLine = commandLineStr.str();

  if (argc < 2) {
    std::cout << "You should pass at least a mesh as first argument"
              << std::endl;
    exit(0);
  }

  for (int i = 1; i < argc; i++) {
    if (GetFileExtension(std::string(argv[i])) == std::string("obj")) {
      PathMesh = std::string(argv[1]);
      has_loaded_mesh = LoadMesh(PathMesh);
      if (!MQuad.MatchPreconditions())
      {
        exit(0);
      }
      // has_loaded_mesh = LoadMesh(PathMesh);
      if (!has_loaded_mesh) {
        std::cout << "ERROR LOADING INPUT MESH" << std::endl;
        exit(0);
      }
      continue;
    }
    if ((GetFileExtension(std::string(argv[i])) == std::string("ffield")) ||
        (GetFileExtension(std::string(argv[i])) == std::string("field"))) {
      PathField = std::string(argv[i]);
      std::cout << "Loading Field" << PathField.c_str() << std::endl;
      has_cross_field = LoadField(PathField);
      // has_cross_field = LoadField(PathField);
      if (!has_cross_field) {
        std::cout << "ERROR LOADING FIELD" << std::endl;
        exit(0);
      }
      continue;
    }
    if (std::string(argv[i]) == std::string("-nofeatures")) {
      // only boundary features will be computed, no internal features
      compute_features = false;
      only_border_corners = true;
      continue;
    }

    if (std::string(argv[i]) == std::string("-angle")) {
      initial_feature_angle = static_cast<ScalarType>(std::stod(argv[i + 1]));
      std::cout << "Initial Feature Angle: " << initial_feature_angle
                << std::endl;
      add_features = true;
      continue;
    }

    if (std::string(argv[i]) == std::string("-size")) {
      SizeMult = static_cast<ScalarType>(std::stod(argv[i + 1]));
      std::cout << "Size Multiplier:" << SizeMult << std::endl;
      continue;
    }

    if (std::string(argv[i]) == std::string("-smooth")) {
      smooth_steps = std::atoi(argv[i + 1]);
      std::cout << "Smooth Steps:" << smooth_steps << std::endl;
      continue;
    }

    if (std::string(argv[i]) == std::string("-nosmooth")) {
      smooth_steps = 0;
      std::cout << "No Smooth Steps" << std::endl;
      continue;
    }

    if (std::string(argv[i]) == std::string("-fastsmooth")) {
      MQuad.reliable_smooth = false;
      std::cout << "Using fast but less reliable smooth" << std::endl;
      continue;
    }

    if (std::string(argv[i]) == std::string("-symmetry")) {
      to_split_symmetry = true;
      std::cout << "Splitting mesh with symmetry plane" << std::endl;
      continue;
    }

    if (std::string(argv[i]) == std::string("-symmetryr")) {
      to_split_symmetry = true;
      DirSymm = Geo::Point3<ScalarType>(-1, 0, 0);
      std::cout << "Splitting mesh with symmetry plane" << std::endl;
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

      has_features = true;
      continue;
    }

    if (GetFileExtension(std::string(argv[i])) == std::string("sharp")) {
      std::string PathFeat = std::string(argv[i]);
      std::cout << "Loading Features" << std::endl;

      bool load_feat =
          ReadSharpFile(PathFeat, TriConnectivity, TriFeatures, false);
      if (!load_feat) {
        std::cout << "ERROR LOADING FEATURES" << std::endl;
        exit(0);
      }

      has_features = true;
      continue;
    }
  }

  if (has_loaded_mesh) {
    if (to_split_symmetry)
      SplitSymmetry();
    if ((!has_features) && (!has_cross_field)) {
      //std::cout << "Step 0 " << std::endl;

      if (extended_print)
        std::cout << "CLEANING MESH" << std::endl;
      Geo::CleanMesh<ScalarType>::CollapseBadShapeFaces(TriVertPos,
                                                        TriConnectivity, 0.05);

      if (extended_print)
        std::cout << "CLEANING FOLDS" << std::endl;
      int steps_clean =
          Geo::CleanMesh<ScalarType>::SolveFolds(TriVertPos, TriConnectivity);

      if (extended_print)
        std::cout << "Performed:" << steps_clean << " unfold steps"
                  << std::endl;

      // WriteOBJ("clean.obj",TriVertPos, TriConnectivity);
      // exit(0);
    }

    if (extended_print)
      std::cout << "BATCH PROCESSING" << std::endl;

    if ((has_features) && (add_features)) {
      std::cout << "** ADDING FEATURES **" << std::endl;

      std::vector<std::pair<int, int>> AddedFeatures;
      std::vector<int> AddedCorners;
      Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
          TriVertPos, TriConnectivity, initial_feature_angle, AddedFeatures,
          AddedCorners, ero_steps, dil_steps);
      TriFeatures.insert(TriFeatures.end(), AddedFeatures.begin(),
                         AddedFeatures.end());
      TriCorners.insert(TriCorners.end(), AddedCorners.begin(),
                        AddedCorners.end());
      std::sort(TriFeatures.begin(), TriFeatures.end());
      auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
      TriFeatures.erase(last, TriFeatures.end());
      std::sort(TriCorners.begin(), TriCorners.end());
      auto last_corners = std::unique(TriCorners.begin(), TriCorners.end());
      TriCorners.erase(last_corners, TriCorners.end());
      has_features = true;
      compute_features = false;
    }
    if ((!has_features) && (compute_features)) {
      if (extended_print)
        std::cout << "** COMPUTING FEATURES **" << std::endl;
      Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
          TriVertPos, TriConnectivity, initial_feature_angle, TriFeatures,
          TriCorners, ero_steps, dil_steps);
      has_features = true;
    }

    if (only_border_corners) {
      Geo::MeshFeatures<ScalarType>::BorderCorners(TriVertPos, TriConnectivity,
                                                   45, TriFeatures, TriCorners);
      has_features = true;
      std::cout << "Only Border Corners, found " << TriCorners.size()
                << " corners" << std::endl;
    }

    //  ComputeNormals(TriVertPos, TriConnectivity, T`riFaceNormal,
    //  TriVertNormal);

    //  Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(TriCorners,
    //  TriCorners);
    //  Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(TriConnectivity,
    //  TriFeatures, TriFeatures); std::cout << "Mapped " << TriFeatures.size()
    //  << " features after cleaning"
    //           << std::endl;

    //  if (extended_print)
    //    std::cout << "GETTING DUAL" << std::endl;

    //  // get the dual
    //  std::vector<int> PerDualPatchFaceIDX, PerDualPatchCornerIDX;
    //  GetDualMesh(TriVertPos, TriConnectivity, TriFeatures, TriCorners,
    if (!has_cross_field) {
      BatchProcessField();
    }

    // std::sort(TriFeatures.begin(), TriFeatures.end());
    // auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
    // TriFeatures.erase(last, TriFeatures.end());

    if (extended_print)
      std::cout << "ESTEEM REMESH SIZE" << std::endl;
    TriMeshBox.Init(TriVertPos);
    AutoQuadEdgeSize();

    if (extended_print)
      std::cout << "QUADRANGULATE" << std::endl;

    MQuad.QuadSize = QuadSize; /// TriMeshBox.Diag();

    MQuad.QuadSmoothSteps = smooth_steps;
    // MQuad.RegVsIso=0.1;
    MQuad.CCAbility = 0.5; // 25;
    MQuad.num_samples = 3000;

    MQuad.ExtractQuad(TriFeatures);

    CreateOutputFolder();
    if (to_split_symmetry) {
      std::set<int> MiddleV;
      MirrorMesh<ScalarType>(QuadVert, QuadFaces, SymmetryPlane, MiddleV);
      RemoveDuplicatedVert(QuadVert, QuadFaces);
    }
    SaveMesh();
    // ProjectName = GetFilanameNoExtension(PathMesh);
    // ProjectName += std::string("_quad.obj");
    // std::cout << "Writing:" << ProjectName << std::endl;

    // WriteOBJ(ProjectName.c_str(), QuadVert, QuadFaces);
  }

  return 0;
}
