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
#include <mesh_clean.h>
#include <edge_mesh_functions.h>
#include <string>
#include <vector>
// GLFWwindow *window;
std::string PathMesh;
// std::string ProjectName;

typedef double ScalarType;

/// MESH DATA
bool loadedMesh = false;

// positions in 3D and UV
std::vector<typename Geo::Point3<ScalarType>> VertPos;
std::vector<std::vector<int>> Connectivity;

std::vector<Geo::Point3<ScalarType>> EdgeMeshVertPos;
std::vector<std::vector<int>> ConnectivityEdges;

std::vector<typename Geo::Point3<ScalarType>> FaceNormals;
std::vector<typename Geo::Point3<ScalarType>> VertNormals;

// Bounding Box in 3D
Geo::Box3<ScalarType> MeshBox;

ScalarType BBDiag = 0;
bool saveStats = false;

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
  RemFile.append("_edge.obj");
  std::cout << "Saving on:" << RemFile.c_str() << std::endl;
  WriteOBJ(RemFile, EdgeMeshVertPos,  ConnectivityEdges);
}


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
      PathMesh = std::string(argv[i]);
      bool has_load = LoadMesh(PathMesh);
      if (!has_load)
        exit(0);
      continue;
    }
  }

  Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromSurfaceMesh(VertPos, Connectivity,
                                                              EdgeMeshVertPos,ConnectivityEdges);
 
  SaveMesh();

  return 0;
}
