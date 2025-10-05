#include <IO/file_path_functions.h>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/_types/_size_t.h>
#include <vector>

// #include <point2.h>
#include "IO/io_features.h"
#include "IO/io_field.h"
#include <IO/io_obj.h>
#include <point3.h>

#include <field_graph/refine_for_tracing.h>
#include <init_cross_field.h>
#include <mesh_quadrangulation.h>
#include <remesh_size.h>

std::string PathMesh;
std::string PathField;
std::string ProjectName;

typedef double ScalarType;

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

ScalarType QuadSize;
std::vector<Geo::Point3<ScalarType>> QuadVert;
std::vector<std::vector<int>> QuadFaces;

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

void SaveMesh() {
  assert(!PathMesh.empty());
  std::string PathRemMesh = RemoveExtension(PathMesh);
  PathRemMesh += std::string("_remesh.obj");
  std::cout << "Saving Field on file " << PathRemMesh.c_str() << std::endl;
  bool Saved = WriteOBJ(PathRemMesh, TriVertPos, TriConnectivity);
  if (!Saved) {
    std::cout << "ERROR SAVING FIELD FILE " << std::endl;
    exit(0);
  }
}

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

  SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormal, FaceCurv,
                       VertCurv);

  RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
  RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
}

void RefineForTracing() {
  Geo::RefineForTracing<ScalarType>::RefineRequiredEdges(
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
}

ScalarType MinSize = 0.005;
ScalarType MaxSize = 0.05;

Geo::Box3<ScalarType> TriMeshBox;

void AutoQuadEdgeSize() {

  std::vector<ScalarType> VertSize;
  ScalarType MinAggrSize = TriMeshBox.Diag() * MinSize;
  ScalarType MaxAggrSize = TriMeshBox.Diag() * MaxSize;
  ScalarType CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeFeatures(
      TriVertPos, TriConnectivity, TriFeatures, MinAggrSize, MaxAggrSize,
      VertSize);
  QuadSize = CurrSize / TriMeshBox.Diag();
  QuadSize /= 4;

  std::cout << "Quad Size: " << QuadSize << std::endl;
}

ScalarType initial_feature_angle = 30;
int ero_steps=3;
int dil_steps=3;
bool compute_features=true;

void UpdateFeatureCorners() {
  Geo::MeshFeatures<ScalarType>::FeatureCorners(TriVertPos, TriConnectivity, initial_feature_angle,
                      TriFeatures, TriCorners);

  Geo::MeshFeatures<ScalarType>::RemoveSingleFaceFeatures(TriConnectivity,TriFeatures);

  if ((ero_steps>0)||(dil_steps>0))
  {
    Geo::MeshFeatures<ScalarType>::ErodeDilateFeatures(ero_steps,dil_steps,TriFeatures) ;
    Geo::MeshFeatures<ScalarType>::Corners(TriVertPos, TriFeatures, initial_feature_angle, TriCorners);
  }
}

std::string folderProjectName;
std::string commandLine;


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
    std::cout << "You should pass at least a mesh as first argument" << std::endl;
    exit(0);
  }

  for (int i = 1; i < argc; i++) {
    if (GetFileExtension(std::string(argv[i])) == std::string("obj")) {
      PathMesh = std::string(argv[1]);
      has_loaded_mesh = LoadMesh(PathMesh);
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
    if (std::string(argv[i]) == std::string("-nofeatures"))
    {
      compute_features=false;
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
  }

  if (has_loaded_mesh) {

    std::cout << "BATCH PROCESSING" << std::endl;
    if ((!has_features)&&(compute_features)) {
      std::cout << "** COMPUTING FEATURES **" << std::endl;
      UpdateFeatureCorners();
    }

    if (!has_cross_field) {
      std::cout << "** COMPUTING FIELD **" << std::endl;
      BatchProcessField();
    }

    // std::sort(TriFeatures.begin(), TriFeatures.end());
    // auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
    // TriFeatures.erase(last, TriFeatures.end());

    std::cout << "ESTEEM REMESH SIZE" << std::endl;
    TriMeshBox.Init(TriVertPos);
    AutoQuadEdgeSize();

    std::cout << "QUADRANGULATE" << std::endl;

    MQuad.QuadSize = QuadSize; /// TriMeshBox.Diag();

    MQuad.QuadSmoothSteps = 100;
    // MQuad.RegVsIso=0.1;
    MQuad.CCAbility = 0.5;
    // MQuad.num_samples=10000;
    MQuad.ExtractQuad(TriFeatures);

    ProjectName = GetFilanameNoExtension(PathMesh);
    ProjectName += std::string("_quad.obj");
    std::cout << "Writing:" << ProjectName << std::endl;

    WriteOBJ(ProjectName.c_str(), QuadVert, QuadFaces);
  }

  return 0;
}
