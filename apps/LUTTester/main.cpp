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
#include <quadrangulation.h>
#include <string>
#include <utility>
#include <vector>
typedef double ScalarType;

/// MESH DATA
bool loadedMesh = false;

std::string PathMesh = "../test_mesh.obj";
std::string OutMesh = "test_mesh_quad.obj";
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

  bool loaded = LoadMesh(PathMesh.c_str(), VertPosTris, ConnectivityTris);
  //Geo::CleanMesh<ScalarType>::RemoveDegenerateFaces(VertPosTris, ConnectivityTris);

  if (!loaded) {
    std::cout << "Mesh not loaded" << std::endl;
    exit(0);
  }

  // then find the vertices of the quadrangulation
  
  std::vector<ScalarType> VertAngle;
  ComputePerVertAngleDeg(VertPosTris, ConnectivityTris, VertAngle);
  std::vector<std::pair<int, ScalarType>> VertAnglePair;
  // get the 4 corners of the quadrangulation, which are the vertices with
  // smaller angle
  for (size_t i = 0; i < VertAngle.size(); i++)
    VertAnglePair.push_back({VertAngle[i], (int)i});

  std::sort(VertAnglePair.begin(), VertAnglePair.end());
  std::vector<int> CornerIdx(4, -1);
  CornerIdx[0] = VertAnglePair[0].second;
  CornerIdx[1] = VertAnglePair[1].second;
  CornerIdx[2] = VertAnglePair[2].second;
  CornerIdx[3] = VertAnglePair[3].second;

  // then sort the vertices in counterclockwise order
  std::vector<std::pair<ScalarType, int>> Angles(4, {0, -1});
  for (size_t i = 0; i < 4; i++) {
    Geo::Point3<ScalarType> P = VertPosTris[CornerIdx[i]];
    ScalarType Angle = atan2(P.Y, P.X);
    Angles[i] = {Angle, CornerIdx[i]};
  }
  std::sort(Angles.begin(), Angles.end());
  std::reverse(Angles.begin(), Angles.end());
  for (size_t i = 0; i < 4; i++)
    CornerIdx[i] = Angles[i].second;

  // CornerIdx[0] = 73;
  // CornerIdx[1] = 2;
  // CornerIdx[2] = 21;
  // CornerIdx[3] = 42;
  // then set the vertices position by interpolating the position of the 4
  // corners on the original mesh
  std::vector<std::vector<Geo::Point3<ScalarType>>> SideSubd;
  for (size_t i = 0; i < 4; i++) {
    Geo::Point3<ScalarType> P0 = VertPosTris[CornerIdx[i]];
    Geo::Point3<ScalarType> P1 = VertPosTris[CornerIdx[(i + 1) % 4]];
    int NumSubd = Div[i];
    // then push number of subdivision points on the side
    std::vector<Geo::Point3<ScalarType>> CurrSideSubd;
    CurrSideSubd.push_back(P0);
    for (size_t j = 1; j < NumSubd; j++) {
      ScalarType t = (ScalarType)(j) / (NumSubd);
      Geo::Point3<ScalarType> P = P0 * (1 - t) + P1 * t;
      CurrSideSubd.push_back(P);
    }
    CurrSideSubd.push_back(P1);
    SideSubd.push_back(CurrSideSubd);
  }
  
  // CornerIdx[0] = 73;
  // CornerIdx[1] = 2;
  // CornerIdx[2] = 21;
  // CornerIdx[3] = 42;

  bool done=Geo::Quadrangulation<ScalarType>::QuadrangulateQuadPatchWithSubsRecursive(
      VertPosTris, ConnectivityTris, CornerIdx, SideSubd, Q_LUT, VertPosQuad,
      ConnectivityQuad);
  if (!done) {
    std::cout << "Quadrangulation failed" << std::endl;
    return 0;
  }

  SaveQuadMesh();

  // // test split
  // int SubSplitEdgeIdx = 0;
  // int SubSplitNum0 = Div[SubSplitEdgeIdx] / 2;
  // int SubSplitNum1 = Div[SubSplitEdgeIdx + 2] / 2;
  // std::cout<<"SubSplitEdge "<<SubSplitEdgeIdx<<" SubSplitNum0 "<<SubSplitNum0<<" SubSplitNum1 "<<SubSplitNum1<<std::endl;
  // SubSplitNum0 = std::max(1, SubSplitNum0);
  // SubSplitNum1 = std::max(1, SubSplitNum1);

  // std::vector<std::vector<Geo::Point3<ScalarType>>> SubPatchVert;
  // std::vector<std::vector<std::vector<int>>> SubPatchFaces;
  // std::vector<std::vector<int>> SubPatchCorners;
  // std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> SubPatchSideSubd;

  // Geo::Quadrangulation<ScalarType>::SplitInfo SInfo;
  // SInfo.SplitEdge = SubSplitEdgeIdx;
  // SInfo.SplitNum0 = SubSplitNum0;
  // SInfo.SplitNum1 = SubSplitNum1;
  // SInfo.NewSplitNum = Div[(SubSplitEdgeIdx+1)%4];
  // std::cout<<"SplitInfo: SplitEdge "<<SInfo.SplitEdge<<" SplitNum0 "<<SInfo.SplitNum0<<" SplitNum1 "<<SInfo.SplitNum1<<" NewSplitNum "<<SInfo.NewSplitNum<<std::endl;
  
  // done=
  // Geo::Quadrangulation<ScalarType>::SplitSideIntoSubPatches(
  //     VertPosTris, ConnectivityTris, CornerIdx, 
  //     SideSubd, SInfo, SubPatchVert, SubPatchFaces, SubPatchCorners, SubPatchSideSubd);
  // if (!done)
  // {
  //   std::cout<<"Split failed"<<std::endl;
  //   return 0;
  // }

  // for (size_t i = 0; i < SubPatchVert.size(); i++)
  // {
  //   std::cout<<"SubPatch "<<i<<" Vert "<<SubPatchVert[i].size()<<" Faces "<<SubPatchFaces[i].size()<<" Corners "<<SubPatchCorners[i].size()<<std::endl;
  //   std::string out_file = std::string("sub_patch_") + std::to_string(i) + std::string(".obj");
  //   WriteOBJ(out_file.c_str(), SubPatchVert[i], SubPatchFaces[i]);

  //   //save the mesh of corners
  //   std::vector<Geo::Point3<ScalarType>> SubPatchCornerVert;
  //   for (size_t j = 0; j < SubPatchCorners[i].size(); j++)
  //     SubPatchCornerVert.push_back(SubPatchVert[i][SubPatchCorners[i][j]]);
  //   out_file = std::string("sub_patch_corner_") + std::to_string(i) + std::string(".obj");
  //   WriteOBJ(out_file.c_str(), SubPatchCornerVert, std::vector<std::vector<int>>());

  //   //save boundary samples
  //   std::vector<Geo::Point3<ScalarType>> SubPatchBoundaryVert;
  //   for (size_t j = 0; j < SubPatchSideSubd[i].size(); j++)
  //     for (size_t k = 0; k < SubPatchSideSubd[i][j].size(); k++)
  //       SubPatchBoundaryVert.push_back(SubPatchSideSubd[i][j][k]);
  //   out_file = std::string("sub_patch_boundary_") + std::to_string(i) + std::string(".obj");
  //   WriteOBJ(out_file.c_str(), SubPatchBoundaryVert, std::vector<std::vector<int>>());

  //   std::vector<Geo::Point3<ScalarType>> SubPatchVertQuad;
  //   std::vector<std::vector<int>> SubPatchFacesQuad;
  //   bool done=Geo::Quadrangulation<ScalarType>::QuadrangulateQuadPatchWithSubs(
  //     SubPatchVert[i], SubPatchFaces[i], SubPatchCorners[i], SubPatchSideSubd[i], 
  //     Q_LUT, SubPatchVertQuad,SubPatchFacesQuad);
  //   if (!done)
  //   {
  //     std::cout<<"Quadrangulation of subpatch "<<i<<" failed"<<std::endl;
  //     continue;
  //   }
  //    out_file = std::string("sub_patch_quad_") + std::to_string(i) + std::string(".obj");
  //    WriteOBJ(out_file.c_str(), SubPatchVertQuad, SubPatchFacesQuad);
  // }
  
  return 0;
} 
