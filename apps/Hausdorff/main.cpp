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
#include <mesh_clean.h>
#include <mesh_functions.h>
#include <point3.h>
#include <color_functions.h>
#include <string>
#include <vector>
// GLFWwindow *window;
std::string PathMesh0, PathMesh1;
// std::string ProjectName;

typedef double ScalarType;

/// MESH DATA
bool loadedMesh = false;

// positions in 3D and UV
std::vector<typename Geo::Point3<ScalarType>> VertPos0, VertPos1;
std::vector<std::vector<int>> Connectivity0, Connectivity1;

std::vector<ScalarType> ErrorFace0, ErrorFace1;
std::vector<Geo::Point3<ScalarType>> ColorFace0, ColorFace1;
std::vector<Geo::Point3<ScalarType>> VertColor0, VertColor1;

// std::vector<typename Geo::Point3<ScalarType>> FaceNormals;
// std::vector<typename Geo::Point3<ScalarType>> VertNormals;

// Bounding Box in 3D
Geo::Box3<ScalarType> MeshBox;

ScalarType BBDiag = 0;
bool saveStats = false;

std::string Outfolder;

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

void SaveMeshes() {
  std::string ErrorFile0 = Outfolder;
  ErrorFile0 += std::string(GetFilanameNoExtension(PathMesh0));
  std::string ErrorFile1 = Outfolder;
  ErrorFile1 += std::string(GetFilanameNoExtension(PathMesh1));

  ErrorFile0.append("_error.obj");
  ErrorFile1.append("_error.obj");
  
  //copy face color to vertex color
  std::cout << "Getting vertex color from face color..." << std::endl;
  Geo::VertColorFromElem(Connectivity0, ColorFace0, VertColor0);
  Geo::VertColorFromElem(Connectivity1, ColorFace1, VertColor1);

  std::cout << "Saving meshes with error color..." << std::endl;
  WriteOBJWithVertColor<ScalarType>(ErrorFile0, VertPos0, Connectivity0, VertColor0);
  WriteOBJWithVertColor<ScalarType>(ErrorFile1, VertPos1, Connectivity1, VertColor1);
}

template <class ScalarType>
bool WriteScalarVector(const std::string &filename,
                       const std::vector<ScalarType> &values)
{
    std::ofstream out(filename.c_str());
    if (!out) return false;

    out << std::setprecision(17);
    for (size_t i = 0; i < values.size(); i++)
        out << values[i] << "\n";

    out.close();
    return true;
}

template <class ScalarType>
bool WriteColorFace(const std::string &filename,
                    const std::vector<Geo::Point3<ScalarType>> &ColorFace)
{
    std::ofstream out(filename.c_str());
    if (!out) return false;

    out << std::setprecision(17);
    for (size_t i = 0; i < ColorFace.size(); i++)
        out << ColorFace[i].X << " " << ColorFace[i].Y << " " << ColorFace[i].Z << "\n";

    out.close();
    return true;
}

int main(int argc, char *argv[]) {

  std::stringstream commandLineStr;
  for (int i = 0; i < argc; ++i) {
    commandLineStr << argv[i] << " ";
  }
  commandLine = commandLineStr.str();

  // // oss<<"Initialized Log Window"<<std::endl;
  // // AScene.SetOutputStream(&oss);

  if (argc < 3) {
    std::cout << "You should pass at least two mesh" << std::endl;
    exit(0);
  }

  PathMesh0 = std::string(argv[1]);
  bool has_load = LoadMesh(PathMesh0, VertPos0, Connectivity0);
  if (!has_load) {
    std::cout << "Mesh 0 not loaded" << std::endl;
    exit(0);
  }

  PathMesh1 = std::string(argv[2]);
  bool has_load1 = LoadMesh(PathMesh1, VertPos1, Connectivity1);
  if (!has_load1) {
    std::cout << "Mesh 1 not loaded" << std::endl;
    exit(0);
  }

  ScalarType Threshold = 0.01;
  if (argc > 3) {
    Threshold = (ScalarType)std::stof(argv[3]);
  }

  MeshBox.Init(VertPos0);

  ScalarType Diag = MeshBox.Diag();
  Threshold = Diag * Threshold;

  TwoWayFaceHausdorff(VertPos0, Connectivity0, VertPos1, Connectivity1,
                      ErrorFace0, ErrorFace1);

  // for (size_t i = 0; i < ErrorFace0.size(); i++)
  //   ErrorFace0[i] = std::min(ErrorFace0[i], Threshold);
  // for (size_t i = 0; i < ErrorFace1.size(); i++)
  //   ErrorFace1[i] = std::min(ErrorFace1[i], Threshold);
  // for (size_t i = 0; i < ErrorFace0.size(); i++)
  //   ErrorFace0[i] = std::max(ErrorFace0[i], (ScalarType)0);
  // for (size_t i = 0; i < ErrorFace1.size(); i++)
  //   ErrorFace1[i] = std::max(ErrorFace1[i], (ScalarType)0);

  // check if there is nan
  for (size_t i = 0; i < ErrorFace0.size(); i++) {
    if (ErrorFace0[i] < 0) {
      std::cout << "ErrorFace0 is < 0 at index " << i << std::endl;
      exit(0);
    }

    if (std::isnan(ErrorFace0[i])) {
      std::cout << "ErrorFace0 has nan at index " << i << std::endl;
      exit(0);
    }
  }
  for (size_t i = 0; i < ErrorFace1.size(); i++) {
    if (ErrorFace1[i] < 0) {
      std::cout << "ErrorFace1 is < 0 at index " << i << std::endl;
      exit(0);
    }

    if (std::isnan(ErrorFace1[i])) {
      std::cout << "ErrorFace1 has nan at index " << i << std::endl;
      exit(0);
    }
  }
  // WriteScalarVector("error_face0.txt", ErrorFace0);
  // WriteScalarVector("error_face1.txt", ErrorFace1);


  Geo::GetColorByScalar<ScalarType>(ErrorFace0, Threshold, 0, ColorFace0);
  Geo::GetColorByScalar<ScalarType>(ErrorFace1, Threshold, 0, ColorFace1);


  // WriteColorFace("color_face0.txt", ColorFace0);
  // WriteColorFace("color_face1.txt", ColorFace1);
  
  // then compute average and max error
  ScalarType AvgError0 = 0, MaxError0 = 0;
  for (size_t i = 0; i < ErrorFace0.size(); i++) {
    AvgError0 += ErrorFace0[i];
    MaxError0 = std::max(MaxError0, ErrorFace0[i]);
  }
  AvgError0 /= ErrorFace0.size();
  ScalarType AvgError1 = 0, MaxError1 = 0;
  for (size_t i = 0; i < ErrorFace1.size(); i++) {
    AvgError1 += ErrorFace1[i];
    MaxError1 = std::max(MaxError1, ErrorFace1[i]);
  }
  AvgError1 /= ErrorFace1.size();
  // scale wrt the threshold
  AvgError0 /= Diag;
  MaxError0 /= Diag;
  AvgError1 /= Diag;
  MaxError1 /= Diag;
  // in percentage
  AvgError0 *= 100;
  MaxError0 *= 100;
  AvgError1 *= 100;
  MaxError1 *= 100;
  std::cout << "Mesh 0: Avg Error: " << AvgError0
            << "% Max Error: " << MaxError0 << "%" << std::endl;
  std::cout << "Mesh 1: Avg Error: " << AvgError1
            << "% Max Error: " << MaxError1 << "%" << std::endl;

  SaveMeshes();

  return 0;
}
