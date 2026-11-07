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
#include <IO/io_stl.h> // STL: loader
#include <cassert>
#include <chrono>
#include <cctype> // STL: std::tolower
#include <cmath>  // MAX FACES: std::sqrt
#include <cstdlib>
#include <filesystem> // FIX: used by CreateOutputFolder, was only included transitively
#include <fstream>    // FIX: used by appendToCSV, was only included transitively
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
// CSV: folder of the executable, where the csv files are collected
std::string ExeFolder;

// STL: true for .stl / .STL files
bool IsSTLFile(const std::string &path) {
  std::string Ext = GetFileExtension(path);
  for (size_t i = 0; i < Ext.size(); i++)
    Ext[i] = (char)std::tolower((unsigned char)Ext[i]);
  return (Ext == std::string("stl"));
}

bool LoadMesh(const std::string &path) {
  // STL: choose the loader from the extension (STL vertices are unified)
  // loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  if (IsSTLFile(path))
    loadedMesh = LoadSTL(path, VertPos, Connectivity, true);
  else
    loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  // FIX: a file read without errors but with no faces (e.g. a free-form
  // NURBS OBJ: only control points and "surf" records) made the rest of the
  // pipeline crash on an empty mesh
  if (loadedMesh && Connectivity.empty()) {
    std::cout << "Error: the mesh has no faces (free-form or empty file)"
              << std::endl;
    loadedMesh = false;
    return false;
  }
  if (loadedMesh) {
    std::cout << "Loaded " << Connectivity.size() << " faces " << VertPos.size()
              << " vertices" << std::endl;

    //Geo::CleanMesh<ScalarType>::SolveNonManifoldByDeletingFaces(VertPos, Connectivity);
    if (!IsAllTris(Connectivity)) {
      std::cout << "Input mesh is not triangular : Triangularize" << std::endl;

      std::vector<typename Geo::Point3<ScalarType>> VertPosTris;
      std::vector<std::vector<int>> ConnectivityTris;

      SplitIntoTris(VertPos, Connectivity, VertPosTris, ConnectivityTris);
      VertPos = VertPosTris;
      Connectivity = ConnectivityTris;
    }
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

// FIX: SaveStats used to call TriRemesh a second time on the already
// remeshed mesh. This doubled the running time and the statistics described
// a mesh that was remeshed twice (not the one saved on disk), with timing of
// the second pass only. Now it only measures the final mesh; the elapsed
// time of the real pipeline is passed in from main.
// void SaveStats() {
void SaveStats(ScalarType ElapsedT) {
  // std::cout<<"SAVING STATS"<<std::endl;
  std::vector<std::string> ToSaveStr;

  ToSaveStr.push_back(GetFilaname(PathMesh));
  // auto start = std::chrono::high_resolution_clock::now();
  // Geo::TriRemesh<ScalarType>(VertPos, Connectivity, RemP);

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

  // auto end = std::chrono::high_resolution_clock::now();
  // auto duration =
  //     std::chrono::duration_cast<std::chrono::microseconds>(end - start)
  //         .count();
  //
  // ScalarType ElapsedT = duration / (ScalarType)1000000;
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

  // CSV: collected in the executable folder
  // std::string StatFile = Outfolder;
  std::string StatFile = ExeFolder;
  StatFile += std::string("statistics.csv");
  std::cout << "Saving statistics on " << StatFile.c_str() << std::endl;
  appendToCSV(StatFile.c_str(), ToSaveStr);
}

void CreateOutputFolder() {

  // OUTPUT: same folder as the input mesh, no "remeshed" subfolder
  Outfolder = GetDirectory(PathMesh);
  // mesh given without a folder: use the current one
  if (Outfolder.empty())
    Outfolder = std::string(".");
  // Outfolder += std::string("/remeshed/");
  Outfolder += std::string("/");
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

// ISSUES: thresholds, the same used by the final cleaning, so the reported
// faces are the ones the cleaning could not fix
ScalarType IssueMinQ = 0.05;      // as CollapseBadShapeFaces
ScalarType IssueFoldDot = -0.99;  // as SolveFolds

// ISSUES: if the final mesh still has very bad faces or folds, append a row
// to issues.csv in the executable folder (created with header if missing)
void SaveIssues() {
  std::vector<bool> DeletedF(Connectivity.size(), false);
  std::vector<int> BadFaces = Geo::CleanMesh<ScalarType>::GetBadShapeFaces(
      VertPos, Connectivity, DeletedF, IssueMinQ);
  std::vector<std::pair<int, int>> Folds;
  Geo::CleanMesh<ScalarType>::GetFolds(VertPos, Connectivity, Folds,
                                       IssueFoldDot);

  if ((BadFaces.size() == 0) && (Folds.size() == 0))
    return;

  ScalarType WorstQ =
      Geo::TriMeshWorstQuality<ScalarType>(VertPos, Connectivity);
  std::cout << "Warning: " << BadFaces.size() << " bad faces, " << Folds.size()
            << " folds" << std::endl;

  std::string IssueFile = ExeFolder; // CSV: executable folder
  IssueFile += std::string("issues.csv");
  if (!std::filesystem::exists(IssueFile)) {
    std::vector<std::string> Header = {"mesh", "faces", "bad_faces", "folds",
                                       "worst_quality", "command"};
    appendToCSV(IssueFile, Header);
  }

  std::vector<std::string> Row;
  Row.push_back(GetFilaname(PathMesh));
  Row.push_back(std::to_string(Connectivity.size()));
  Row.push_back(std::to_string(BadFaces.size()));
  Row.push_back(std::to_string(Folds.size()));
  Row.push_back(std::to_string(WorstQ));
  // command line in quotes: options with commas would break the columns
  Row.push_back(std::string("\"") + commandLine + std::string("\""));
  std::cout << "Saving issues on " << IssueFile.c_str() << std::endl;
  appendToCSV(IssueFile, Row);
}

// TIMING: print the time elapsed since t0 for a pipeline phase, then reset t0
void PrintPhaseTime(const std::string &Phase,
                    std::chrono::high_resolution_clock::time_point &t0) {
  auto t1 = std::chrono::high_resolution_clock::now();
  ScalarType Secs =
      std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() /
      (ScalarType)1000000;
  std::cout << "** TIME " << Phase << ": " << Secs << " s" << std::endl;
  t0 = t1;
}

// MAX FACES: bound on the expected output faces (-maxfaces, 0 disables)
int MaxFaces = 2000000;
// MAX FACES: hard bound on the faces during refinement, whatever the
// estimate says: when exceeded the mesh fails with an error (bounds memory)
#define TRI_REMESH_HARD_MAX_FACES 5000000

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

  // CSV: folder of the executable from argv[0]. If launched through PATH
  // (no folder in argv[0]) this falls back to the current folder.
  ExeFolder =
      std::filesystem::absolute(std::filesystem::path(argv[0]))
          .parent_path()
          .string() +
      std::string("/");

  // // oss<<"Initialized Log Window"<<std::endl;
  // // AScene.SetOutputStream(&oss);

  if (argc < 2) {
    std::cout << "You should pass at least a mesh" << std::endl;
    exit(0);
  }

  // FIX: options with a value used "if (argc > i)", which is always true,
  // so a trailing "-res" read argv[argc] (NULL) and crashed in std::stof.
  // Now the value is checked to exist and is consumed (i++).
  for (int i = 1; i < argc; i++) {
    // STL: accept also .stl meshes
    // if (GetFileExtension(std::string(argv[i])) == std::string("obj")) {
    if ((GetFileExtension(std::string(argv[i])) == std::string("obj")) ||
        IsSTLFile(std::string(argv[i]))) {
      // FIX: was argv[1], wrong if the mesh is not the first argument
      // PathMesh = std::string(argv[1]);
      PathMesh = std::string(argv[i]);
      bool has_load = LoadMesh(PathMesh);
      // if (!has_load)
      //   exit(0);
      if (!has_load) {
        std::cout << "Error: cannot load the mesh " << PathMesh << std::endl;
        exit(0);
      }
      continue;
    }

    if (std::string(argv[i]) == std::string("-res")) {
      // if (argc > i) {
      if (i < argc - 1) {
        ratio = (ScalarType)std::stof(argv[i + 1]);
        i++;
        //% of bounding box diagonal
        // RemP.TargetL = (BBDiag * ratio) / (ScalarType)100;
      }
      continue;
    }

    if (std::string(argv[i]) == std::string("-ite")) {
      RemP.iterative = true;
      continue;
    }

    if (std::string(argv[i]) == std::string("-angle")) {
      // if (argc > i) {
      if (i < argc - 1) {
        ScalarType angle = (ScalarType)std::stof(argv[i + 1]);
        RemP.feature_angle = angle;
        i++;
      }
      continue;
    }

    if (std::string(argv[i]) == std::string("-steps")) {
      // if (argc > i) {
      if (i < argc - 1) {
        RemP.steps = std::stoi(argv[i + 1]);
        i++;
      }
      continue;
    }

    if (std::string(argv[i]) == std::string("-fixborders")) {
      // if (argc > i) {
      RemP.fix_borders = true;
      // }
      continue;
    }

    if (std::string(argv[i]) == std::string("-adapt")) {
      // RemP.max_error=BBDiag*0.001;
      RemP.adaptive = true;
      continue;
    }

    if (std::string(argv[i]) == std::string("-stats")) {
      saveStats = true;
      continue;
    }

    // TIMING: print the per-step times of the remesher
    if (std::string(argv[i]) == std::string("-verbose")) {
      RemP.writeMSG = true;
      continue;
    }

    // MAX FACES: "-maxfaces N" bounds the expected output faces (0 = no bound)
    if (std::string(argv[i]) == std::string("-maxfaces")) {
      if (i < argc - 1) {
        MaxFaces = std::stoi(argv[i + 1]);
        i++;
      }
      continue;
    }

    // EARLY STOP: "-earlystop" uses 0.001 (0.1% of the faces),
    // "-earlystop <ratio>" sets it explicitly
    if (std::string(argv[i]) == std::string("-earlystop")) {
      RemP.early_stop_ratio = 0.001;
      if ((i < argc - 1) && (argv[i + 1][0] != '-') &&
          (GetFileExtension(std::string(argv[i + 1])) != std::string("obj"))) {
        RemP.early_stop_ratio = (ScalarType)std::stof(argv[i + 1]);
        i++;
      }
      continue;
    }

    // FIX: report unknown options instead of silently ignoring them
    std::cout << "Warning: unknown option " << argv[i] << std::endl;
  }

  // FIX: no mesh among the arguments
  if (!loadedMesh) {
    std::cout << "Error: no .obj mesh loaded" << std::endl;
    exit(0);
  }

  // NON MANIFOLD: remove duplicated faces, split non-manifold edges and
  // vertices, pull apart coincident vertices (1% of the average edge
  // length). A clean manifold mesh is not modified.
  Geo::CleanMesh<ScalarType>::SplitNonManifold(VertPos, Connectivity);

  if (!MatchTriRemeshPreconditions(VertPos, Connectivity)) {
    // DIAGNOSTIC: report which precondition fails (it ends up in fail.txt)
    // std::cout << "Error: Input mesh does not match preconditions" << std::endl;
    bool ManifE = EdgeManifold<ScalarType>(VertPos, Connectivity);
    bool ManifV = VertManifold<ScalarType>(VertPos, Connectivity);
    std::cout << "Error: Input mesh does not match preconditions (edge manifold "
              << ManifE << ", vertex manifold " << ManifV << ")" << std::endl;
    exit(0);
  }

  // FIX: timing now covers the whole pipeline actually used to produce the
  // output mesh (features, size estimation, remesh, cleaning)
  auto start = std::chrono::high_resolution_clock::now();
  auto tPhase = start; // TIMING

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
  PrintPhaseTime("features", tPhase); // TIMING

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

  // FIX: assert disappears in release builds, check explicitly
  // assert(CurrSize > 0);
  if (!(CurrSize > 0)) {
    std::cout << "Error: invalid remesh size " << CurrSize << std::endl;
    exit(0);
  }
  std::cout << "Found Remesh Size of " << CurrSize / MeshBox.Diag()
            << " BB Diagonal" << std::endl;

  // ADAPTIVE: add the curvature/feature term to the size field, bounded
  // below by the global size (same units, before halving). Not used without
  // -adapt, so the non-adaptive path is unchanged.
  if (RemP.adaptive)
    Geo::RemeshSize<ScalarType>::AdaptiveSizeField(
        VertPos, Connectivity, CurrSize, MinAggrSize, MaxAggrSize,
        RemP.VertSize);
  PrintPhaseTime("size estimation", tPhase); // TIMING

  CurrSize/=2;
  // FIX (behavior change, please verify): MinSize is a fraction of the
  // bounding box diagonal, while CurrSize is absolute. The old clamp mixed
  // units: harmless for meshes with diagonal ~1, but on small meshes it
  // forced a huge target length (losing detail/features) and on large ones
  // it did nothing.
  // CurrSize=std::max(CurrSize,MinSize);
  // FIX: clamping to MinAggrSize doubled the minimum size on meshes with
  // diagonal > 1 (triangles too big near fine details). No clamp needed:
  // RemeshSizeAggregated already bounds every sample to >= MinAggrSize, so
  // after halving the size is >= MinAggrSize/2 in any unit.
  // CurrSize = std::max(CurrSize, MinAggrSize);
  RemP.TargetL = ratio * CurrSize;

  // FIX: the adaptive per-vertex size went straight to the remesher, skipping
  // the halving, the clamp and -res. It now gets the same transformation as
  // the global size, so -adapt and non-adaptive runs are on the same scale.
  if (RemP.adaptive) {
    assert(RemP.VertSize.size() == VertPos.size());
    for (size_t i = 0; i < RemP.VertSize.size(); i++) {
      ScalarType CurrVSize = RemP.VertSize[i] / 2;
      // CurrVSize = std::max(CurrVSize, MinAggrSize); // FIX: see above
      RemP.VertSize[i] = ratio * CurrVSize;
    }
  }

  // MAX FACES: TargetL is a fraction of the bbox diagonal, but the output
  // size depends on the area: meshes with a lot of area for their diagonal
  // (lattices, thin walls, many parts) produced tens of millions of faces
  // and tens of GB. The expected output (area / equilateral triangle area)
  // is bounded by MaxFaces by enlarging TargetL (and the adaptive field by
  // the same factor). MaxFaces <= 0 disables the bound.
  if (MaxFaces > 0) {
    ScalarType MeshArea = Area(VertPos, Connectivity);
    ScalarType TriArea = (std::sqrt((ScalarType)3) / 4) * RemP.TargetL * RemP.TargetL;
    ScalarType ExpFaces = MeshArea / TriArea;
    if (ExpFaces > (ScalarType)MaxFaces) {
      ScalarType Scale = std::sqrt(ExpFaces / (ScalarType)MaxFaces);
      std::cout << "Warning: expected " << (long long)ExpFaces
                << " faces, target length enlarged by " << Scale << " to stay below "
                << MaxFaces << std::endl;
      RemP.TargetL *= Scale;
      if (RemP.adaptive)
        for (size_t i = 0; i < RemP.VertSize.size(); i++)
          RemP.VertSize[i] *= Scale;
    }
  }

  RemP.max_faces = TRI_REMESH_HARD_MAX_FACES; // MAX FACES
  RemP.check_manifold = true; // SAFETY: stop if a step breaks manifoldness

  std::cout << "REMESHING" << std::endl;

  // Geo::TriRemesh<ScalarType>(VertPos, Connectivity, RemP);
  int RemeshSteps = Geo::TriRemesh<ScalarType>(VertPos, Connectivity, RemP);
  // MAX FACES: the refinement hit the hard bound (message already printed)
  if (RemeshSteps < 0)
    exit(0);
  PrintPhaseTime("remesh", tPhase); // TIMING

  // FIX: the cleaning steps ran without feature information, collapsing
  // corners/feature edges and smoothing feature vertices. After TriRemesh the
  // feature arrays in RemP match the remeshed mesh, so they are passed on.
  bool FeaturesOK = (RemP.IsFaceEFeature.size() == Connectivity.size()) &&
                    (RemP.IsVertCorner.size() == VertPos.size());
  if (!FeaturesOK)
    std::cout << "Warning: feature data out of sync, cleaning without features"
              << std::endl;

  std::cout<<"CLEANING MESH"<<std::endl;
  // Geo::CleanMesh<ScalarType>::CollapseBadShapeFaces(VertPos, Connectivity,0.05);
  Geo::CleanMesh<ScalarType>::CollapseBadShapeFaces(
      VertPos, Connectivity, 0.05, FeaturesOK ? &RemP.IsFaceEFeature : NULL,
      FeaturesOK ? &RemP.IsVertCorner : NULL);
  PrintPhaseTime("clean collapse", tPhase); // TIMING

  FeaturesOK = FeaturesOK &&
               (RemP.IsFaceEFeature.size() == Connectivity.size()) &&
               (RemP.IsVertCorner.size() == VertPos.size());

  // fixed vertices for unfolding: corners and feature-line vertices
  std::vector<bool> IsFixedV(VertPos.size(), false);
  if (FeaturesOK) {
    std::vector<std::pair<int, int>> FeatPairs;
    Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
        Connectivity, RemP.IsFaceEFeature, FeatPairs);
    for (size_t i = 0; i < FeatPairs.size(); i++) {
      IsFixedV[FeatPairs[i].first] = true;
      IsFixedV[FeatPairs[i].second] = true;
    }
    for (size_t i = 0; i < VertPos.size(); i++)
      if (RemP.IsVertCorner[i])
        IsFixedV[i] = true;
  }

  std::cout<<"CLEANING FOLDS"<<std::endl;
  // int steps_clean=Geo::CleanMesh<ScalarType>::SolveFolds(VertPos, Connectivity);
  int steps_clean = Geo::CleanMesh<ScalarType>::SolveFolds(
      VertPos, Connectivity, -0.99, FeaturesOK ? &IsFixedV : NULL);
  std::cout<<"Performed:"<<steps_clean<<" unfold steps"<<std::endl;
  PrintPhaseTime("clean folds", tPhase); // TIMING

  auto end = std::chrono::high_resolution_clock::now();
  ScalarType ElapsedT =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count() /
      (ScalarType)1000000;
  std::cout << "Elapsed: " << ElapsedT << " s" << std::endl;

  std::cout << "DONE" << std::endl;
  CreateOutputFolder();

  SaveMesh();

  SaveIssues(); // ISSUES

  if (saveStats)
    // SaveStats();
    SaveStats(ElapsedT);

  return 0;
}
