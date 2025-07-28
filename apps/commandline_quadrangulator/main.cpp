#include <IO/file_path_functions.h>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
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

// // // FIELD DATA
std::vector<Field::CrossF<ScalarType>> FaceCurv; //,VertCurv;
FieldSmoothParam<ScalarType> SParam;

bool has_loaded_mesh = false;
bool has_cross_field = false;
bool has_features = false;

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


void BatchProcessField() {
  Geo::RefineForTracing<ScalarType>::RefineRequiredEdges(
      TriVertPos, TriConnectivity, TriFeatures);

   SParam.FixedCross.clear();
  //  InitFieldByCurvature();
  //  SmoothGlobal();
  InitCrossFieldParam<ScalarType> InitFParam;
  std::vector<Field::CrossF<ScalarType>> VertCurv;
  InitMeshCrossField<ScalarType>(TriVertPos, TriConnectivity, TriFeatures,
                                 VertCurv, FaceCurv, InitFParam);

  Geo::RefineForTracing<ScalarType>::RefineForConnection(
      TriVertPos, TriConnectivity, FaceCurv, TriFeatures);
}

std::string folderProjectName;
std::string commandLine;

ScalarType feature_angle = 45;

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

      has_features = true;
      continue;
    }
  }

  if ((has_loaded_mesh)&&(has_cross_field)) {
    std::cout << "BATCH PROCESSING" << std::endl;
    // if (!has_features) {
    //   std::vector<int> Corners;
    //   std::cout << "** COMPUTING FEATURES **" << std::endl;
    //   Geo::FeatureCorners(TriVertPos, TriConnectivity, feature_angle,
    //                       TriFeatures, Corners);

    //   // WriteFeatures("test.features", TriFeatures);
    // }
    // if (!has_cross_field) {
    //   std::cout << "** COMPUTING FIELD **" << std::endl;
    //   BatchProcessField();
    //   // bool Saved = WriteField("test.field", VertCurv, FaceCurv);
    // }
    std::cout << "ESTEEM REMESH SIZE" << std::endl;
    Geo::Box3<ScalarType> MeshBox;
    MeshBox.Init(TriVertPos);
    ScalarType MinSize=MeshBox.Diag() * 0.005;
    ScalarType MaxSize=MeshBox.Diag() * 0.02;
    ScalarType CurrSize=Geo::RemeshSize<ScalarType>::RemeshSizeNormal(TriVertPos,TriConnectivity, 
                                                                      MinSize,MaxSize);  
    MQuad.QuadSize = CurrSize/MeshBox.Diag();//0.02;                                                                 
    // std::cout << "Size Computed " << CurrSize<<std::endl;
    // std::cout << "Size Base " << MQuad.QuadSize*MeshBox.Diag()<<std::endl;
    // exit(0);                       
    MQuad.QuadSmoothSteps=100;
    // MQuad.RegVsIso=0.1;
    MQuad.CCAbility=1;
    //MQuad.num_samples=10000;
    MQuad.ExtractQuad(TriFeatures);

    ProjectName=GetFilanameNoExtension(PathMesh);
    ProjectName+=std::string("_quad.obj");
    WriteOBJ(ProjectName.c_str(), QuadVert, QuadFaces);
  }

  return 0;
}
