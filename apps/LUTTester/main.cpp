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
#include <mesh_clean.h>
#include <mesh_functions.h>
#include <point3.h>
#include <color_functions.h>
#include <string>
#include <vector>
#include <quadrangulation.h>
typedef double ScalarType;

/// MESH DATA
bool loadedMesh = false;

std::string PathMesh = "../test_mesh.obj";
std::string OutMesh = "../test_mesh_quad.obj";
// positions in 3D 
std::vector<typename Geo::Point3<ScalarType>> VertPosTris;
std::vector<typename Geo::Point2<ScalarType>> VertPosUV;
std::vector<std::vector<int>> ConnectivityTris;

std::vector<typename Geo::Point3<ScalarType>> VertPosQuad;
std::vector<std::vector<int>> ConnectivityQuad;


bool LoadMesh(const std::string &path,
              std::vector<typename Geo::Point3<ScalarType>> &CurrVertPos,
              std::vector<std::vector<int>> &CurrConnectivity) {
  loadedMesh = LoadOBJ(path, CurrVertPos, CurrConnectivity, true);
  if (loadedMesh) {
    std::cout << "Loaded " << CurrConnectivity.size() << " faces "
              << CurrVertPos.size() << " vertices" << std::endl;
    return true;
    // IteOpt.NextCollapseTest(NextV0, NextV1);
  } else {
    std::cout << "Mesh not loaded" << std::endl;
    return false;
  }
}

// std::string folderProjectName;
std::string commandLine;

void SaveQuadMesh() {
  WriteOBJ<ScalarType>(OutMesh.c_str(), VertPosQuad, ConnectivityQuad);
}

int Div[4] = {4, 4, 4, 4};

int main(int argc, char *argv[]) {

  Geo::QuadrangulationLUT<ScalarType> Q_LUT;
  Q_LUT.Initialize();

  if (argc > 4) {
      Div[0] = std::stoi(argv[1]);
      Div[1] = std::stoi(argv[2]);
      Div[2] = std::stoi(argv[3]);
      Div[3] = std::stoi(argv[4]);
  }

  bool loaded=  LoadMesh(PathMesh.c_str(), VertPosTris, ConnectivityTris);
  if (!loaded) {
    std::cout << "Mesh not loaded" << std::endl;
    exit(0);
  }

  //then find the vertices of the quadrangulation
  
  //ScalarType Threshold = 0.01;
  
  // MeshBox.Init(VertPos0);

  // ScalarType Diag = MeshBox.Diag();
  // Threshold = Diag * Threshold;

  // TwoWayFaceHausdorff(VertPos0, Connectivity0, VertPos1, Connectivity1,
  //                     ErrorFace0, ErrorFace1);

  // // for (size_t i = 0; i < ErrorFace0.size(); i++)
  // //   ErrorFace0[i] = std::min(ErrorFace0[i], Threshold);
  // // for (size_t i = 0; i < ErrorFace1.size(); i++)
  // //   ErrorFace1[i] = std::min(ErrorFace1[i], Threshold);
  // // for (size_t i = 0; i < ErrorFace0.size(); i++)
  // //   ErrorFace0[i] = std::max(ErrorFace0[i], (ScalarType)0);
  // // for (size_t i = 0; i < ErrorFace1.size(); i++)
  // //   ErrorFace1[i] = std::max(ErrorFace1[i], (ScalarType)0);

  // // check if there is nan
  // for (size_t i = 0; i < ErrorFace0.size(); i++) {
  //   if (ErrorFace0[i] < 0) {
  //     std::cout << "ErrorFace0 is < 0 at index " << i << std::endl;
  //     exit(0);
  //   }

  //   if (std::isnan(ErrorFace0[i])) {
  //     std::cout << "ErrorFace0 has nan at index " << i << std::endl;
  //     exit(0);
  //   }
  // }
  // for (size_t i = 0; i < ErrorFace1.size(); i++) {
  //   if (ErrorFace1[i] < 0) {
  //     std::cout << "ErrorFace1 is < 0 at index " << i << std::endl;
  //     exit(0);
  //   }

  //   if (std::isnan(ErrorFace1[i])) {
  //     std::cout << "ErrorFace1 has nan at index " << i << std::endl;
  //     exit(0);
  //   }
  // }
  // // WriteScalarVector("error_face0.txt", ErrorFace0);
  // // WriteScalarVector("error_face1.txt", ErrorFace1);


  // Geo::GetColorByScalar<ScalarType>(ErrorFace0, Threshold, 0, ColorFace0);
  // Geo::GetColorByScalar<ScalarType>(ErrorFace1, Threshold, 0, ColorFace1);


  // // WriteColorFace("color_face0.txt", ColorFace0);
  // // WriteColorFace("color_face1.txt", ColorFace1);
  
  // // then compute average and max error
  // ScalarType AvgError0 = 0, MaxError0 = 0;
  // for (size_t i = 0; i < ErrorFace0.size(); i++) {
  //   AvgError0 += ErrorFace0[i];
  //   MaxError0 = std::max(MaxError0, ErrorFace0[i]);
  // }
  // AvgError0 /= ErrorFace0.size();
  // ScalarType AvgError1 = 0, MaxError1 = 0;
  // for (size_t i = 0; i < ErrorFace1.size(); i++) {
  //   AvgError1 += ErrorFace1[i];
  //   MaxError1 = std::max(MaxError1, ErrorFace1[i]);
  // }
  // AvgError1 /= ErrorFace1.size();
  // // scale wrt the threshold
  // AvgError0 /= Diag;
  // MaxError0 /= Diag;
  // AvgError1 /= Diag;
  // MaxError1 /= Diag;
  // // in percentage
  // AvgError0 *= 100;
  // MaxError0 *= 100;
  // AvgError1 *= 100;
  // MaxError1 *= 100;
  // std::cout << "Mesh 0: Avg Error: " << AvgError0
  //           << "% Max Error: " << MaxError0 << "%" << std::endl;
  // std::cout << "Mesh 1: Avg Error: " << AvgError1
  //           << "% Max Error: " << MaxError1 << "%" << std::endl;

  // SaveMeshes();

  return 0;
}
