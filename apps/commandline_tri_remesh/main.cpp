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
#include <IO/io_obj.h>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <box3.h>
#include <color_functions.h>
#include <connectivity.h>
#include <hausdorff.h>
#include <mesh_functions.h>
#include <point3.h>
#include <remesh_size.h>
#include <string>
#include <triangular_remesh.h>
#include <vector>
// GLFWwindow *window;
std::string PathMesh;
// std::string ProjectName;

typedef double ScalarType;

/// MESH DATA
bool loadedMesh = false;

// positions in 3D and UV
std::vector<typename Geo::Point3<ScalarType>> VertPos, VertPos0;
std::vector<std::vector<int>> Connectivity, Connectivity0;

std::vector<typename Geo::Point3<ScalarType>> FaceNormals;
std::vector<typename Geo::Point3<ScalarType>> VertNormals;

// Bounding Box in 3D
Geo::Box3<ScalarType> MeshBox;

ScalarType BBDiag = 0;
bool saveStats = false;

Geo::TriRemParam<ScalarType> RemP;

std::string Outfolder;

bool LoadMesh(const std::string &path) {
  loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  if (loadedMesh) {
    std::cout << "Loaded " << Connectivity.size() << " faces " << VertPos.size()
              << " vertices" << std::endl;

    MeshBox.Init(VertPos);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    // AVGL = AvgEdgeLen(VertPos, Connectivity);
    BBDiag = MeshBox.Diag();

    VertPos0 = VertPos;
    Connectivity0 = Connectivity;
    return true;
    // IteOpt.NextCollapseTest(NextV0, NextV1);
  } else {
    std::cout << "Mesh not loaded" << std::endl;
    return false;
  }
}

// std::string folderProjectName;
std::string commandLine;

void SaveMesh() {
  std::string RemFile = Outfolder;
  RemFile += std::string(GetFilanameNoExtension(PathMesh));
  RemFile.append("_remesh.obj");
  std::cout << "Saving on:" << RemFile.c_str() << std::endl;
  WriteOBJ(RemFile, VertPos, Connectivity);
}

void appendToCSV(const std::string &filename,
                 const std::vector<std::string> &str_field) {
  // Open the file in append mode (creates file if it doesn't exist)
  std::ofstream file;
  file.open(filename, std::ios::app);

  // Check if the file is open
  if (file.is_open()) {
    for (size_t i = 0; i < str_field.size() - 1; i++)
      file << str_field[i] << ",";

    file << str_field[str_field.size() - 1] << "\n";
    file.close();
    // std::cout << "Row appended to " << filename << std::endl;
  }
  // } else {
  //     std::cerr << "Unable to open file: " << filename << std::endl;
  // }
}

void SaveStats() {
  // std::cout<<"SAVING STATS"<<std::endl;
  std::vector<std::string> ToSaveStr;

  ToSaveStr.push_back(GetFilaname(PathMesh));
  auto start = std::chrono::high_resolution_clock::now();
  Geo::TriRemesh<ScalarType>(VertPos, Connectivity, RemP);

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count();

  ScalarType ElapsedT = duration / (ScalarType)1000000;
  ToSaveStr.push_back(std::to_string(ElapsedT));

  ToSaveStr.push_back(std::to_string(RemP.TargetL));
  ToSaveStr.push_back(std::to_string(AvgEdgeLen(VertPos, Connectivity)));

  ScalarType minLen, maxLen;
  ExtremenEdgeLen(VertPos, Connectivity, minLen, maxLen);
  ToSaveStr.push_back(std::to_string(minLen));
  ToSaveStr.push_back(std::to_string(maxLen));

  ToSaveStr.push_back(std::to_string(VertPos.size()));
  ToSaveStr.push_back(std::to_string(Connectivity.size()));

  ScalarType AvgQ = Geo::TriMeshAvgQuality<ScalarType>(VertPos, Connectivity);
  ScalarType WorstQ =
      Geo::TriMeshWorstQuality<ScalarType>(VertPos, Connectivity);

  ToSaveStr.push_back(std::to_string(AvgQ));
  ToSaveStr.push_back(std::to_string(WorstQ));

  ScalarType AvgDist, MaxDist;
  Geo::Hausdorff(VertPos0, Connectivity0, VertPos, Connectivity, AvgDist,
                 MaxDist);

  ToSaveStr.push_back(std::to_string(100 * AvgDist / BBDiag));
  ToSaveStr.push_back(std::to_string(100 * MaxDist / BBDiag));

  std::string StatFile = Outfolder;
  StatFile += std::string("statistics.csv");
  std::cout << "Saving statistics on " << StatFile.c_str() << std::endl;
  appendToCSV(StatFile.c_str(), ToSaveStr);
}

void CreateOutputFolder() {

  Outfolder = GetDirectory(PathMesh);
  Outfolder += std::string("/remeshed/");
  // std::cout<<"Creating Folder"<<Outfolder<<std::endl;

  try {
    // Check if the folder already exists
    if (!std::filesystem::exists(Outfolder))
      std::filesystem::create_directory(Outfolder);
    //{
    // Create the folder
    //     std::cout << "Folder created successfully: " << folderProjectName <<
    //     std::endl;
    // }
    // else
    // {
    //     std::cout << "Folder already exists: " << folderProjectName <<
    //     std::endl;
    // }
  } catch (const std::exception &e) {
    std::cerr << "Error creating folder: " << e.what() << std::endl;
  }
  // exit(0);
}

ScalarType MinSize = 0.005;
ScalarType MaxSize = 0.02;
int ero_steps = 3;
int dil_steps = 3;
ScalarType ratio = 1;

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
      bool has_load = LoadMesh(PathMesh);
      if (!has_load)
        exit(0);
      continue;
    }

    if (std::string(argv[i]) == std::string("-res")) {
      if (argc > i) {
        ratio = (ScalarType)std::stof(argv[i + 1]);
        //% of bounding box diagonal
        // RemP.TargetL = (BBDiag * ratio) / (ScalarType)100;
      }
    }

    if (std::string(argv[i]) == std::string("-ite")) {
      RemP.iterative = true;
    }

    if (std::string(argv[i]) == std::string("-angle")) {
      if (argc > i) {
        ScalarType angle = (ScalarType)std::stof(argv[i + 1]);
        RemP.feature_angle = angle;
      }
    }

    if (std::string(argv[i]) == std::string("-steps")) {
      if (argc > i) {
        RemP.steps = std::stoi(argv[i + 1]);
      }
    }

    if (std::string(argv[i]) == std::string("-fixborders")) {
      if (argc > i) {
        RemP.fix_borders = true;
      }
    }

    if (std::string(argv[i]) == std::string("-adapt")) {
      // RemP.max_error=BBDiag*0.001;
      RemP.adaptive = true;
    }

    if (std::string(argv[i]) == std::string("-stats")) {
      saveStats = true;
    }
  }

  if (!MatchTriRemeshPreconditions(VertPos, Connectivity)) {
    std::cout << "Error: Input mesh does not match preconditions" << std::endl;
    exit(0);
  }

  // Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
  //       VertPos, Connectivity, RemP.feature_angle, TriFeatures_30,
  //       TriCorners_30, ero_steps, dil_steps);

  // INITIALIZE FEATURES
  if ((RemP.feature_angle > 0) && (RemP.feature_angle < 180)) {
    Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
        VertPos, Connectivity, RemP.feature_angle, RemP.IsFaceEFeature,
        RemP.IsVertCorner, ero_steps, dil_steps);
  }
  else {
     RemP.IsVertCorner=std::vector<bool>(VertPos.size(),false);
     RemP.IsFaceEFeature=std::vector<std::vector<bool>>(Connectivity.size(),std::vector<bool>(3,false));
  }

  // FIND REMESHING SIZE
  ScalarType MinAggrSize = MeshBox.Diag() * MinSize;
  ScalarType MaxAggrSize = MeshBox.Diag() * MaxSize;

  ScalarType CurrSize = -1;
  std::cout<<"Finding Size"<<std::endl;
  if (RemP.adaptive)
    CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeAggregated(
        VertPos, Connectivity, RemP.IsFaceEFeature, MinAggrSize, MaxAggrSize,
        &RemP.VertSize);
  else
    CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeAggregated(
        VertPos, Connectivity, RemP.IsFaceEFeature, MinAggrSize, MaxAggrSize,
        NULL);

  assert(CurrSize > 0);
  std::cout << "Found Remesh Size of " << CurrSize / MeshBox.Diag()
            << " BB Diagonal" << std::endl;
  CurrSize/=2;
  CurrSize=std::max(CurrSize,MinSize);
  RemP.TargetL = ratio * CurrSize;

  std::cout << "REMESHING" << std::endl;
  
  Geo::TriRemesh<ScalarType>(VertPos, Connectivity, RemP);
  
  std::cout << "DONE" << std::endl;
  CreateOutputFolder();

  SaveMesh();

  if (saveStats)
    SaveStats();

  return 0;
}
