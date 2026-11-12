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


// Interactive triangle remesher: same steps as the command line tool
// (commandline_tri_remesh): non-manifold split, features, target size,
// TriRemesh, cleaning. Every remesh starts from the loaded mesh, so the
// parameters can be tuned and compared; the view switches between the
// input and the result.

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <GL/glew.h>
#include <UI/orbitball.h>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <iostream>
#include <stdio.h>
#include <string>
#include <utility>
#include <set>
#include <vector>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "IO/io_obj.h"
#include "IO/io_stl.h"
#include <GLFW/glfw3.h>
#include <GL_functions/GLDrawMesh.h>
#include <IO/file_path_functions.h>
#include <box3.h>
#include <chrono>
#include <color_functions.h>
#include <connectivity.h>
#include <hausdorff.h>
#include <local_operations/edge_flip.h>
#include <local_operations/iterative_optimizer.h>
#include <normal.h>
#include <point3.h>
#include <mesh_clean.h>
#include <remesh_size.h>
#include <triangular_remesh.h>
#include <GL_functions/GLPolylinePick.h>
#include <mesh_embed_polyline.h>

typedef double ScalarType;

GLFWwindow *window;
std::string PathMesh;
UI::OrbitBall trackball; //,trackballUV;

// std::ostringstream oss;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
Geo::Point3<ScalarType> mesh_color(0.79, 0.9, 0.87);

float line_thick = 1.f;
int selected_rend_item = 5;
bool showMesh = true;
bool showFeatures = true;
bool showCorners = true;

// [Win32] Our example includes a copy of glfw3.lib pre-compiled with VS2010 to
// maximize ease of testing and compatibility with old VS compilers. To link
// with VS2010-era libraries, VS2015+ requires linking with
// legacy_stdio_definitions.lib, which we do using this pragma. Your own project
// should not be affected, as you are likely to link with a newer binary of GLFW
// that is adequate for your version of Visual Studio.
#if defined(_MSC_VER) && (_MSC_VER >= 1900) &&                                 \
    !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

// MESH DATA
bool loadedMesh = false;
std::string LoadMsg; // why the mesh cannot be remeshed (empty if it can)

// DISPLAYED mesh (input or result)
std::vector<typename Geo::Point3<ScalarType>> VertPos;
std::vector<std::vector<int>> Connectivity;
// INPUT mesh (after the non-manifold split) and its features
std::vector<typename Geo::Point3<ScalarType>> VertPos0;
std::vector<std::vector<int>> Connectivity0;
std::vector<std::vector<bool>> IsFaceEFeature0;
std::vector<bool> IsVertCorner0;
// RESULT of the last remesh
bool hasResult = false;
std::vector<typename Geo::Point3<ScalarType>> VertPosR;
std::vector<std::vector<int>> ConnectivityR;
std::vector<std::vector<bool>> IsFaceEFeatureR;
std::vector<bool> IsVertCornerR;
int viewResult = 0; // 0 input, 1 result

std::vector<Geo::Point3<ScalarType>> FaceNormals;
std::vector<Geo::Point3<ScalarType>> VertNormals;

// Bounding Box in 3D (of the input: the view does not jump)
Geo::Box3<ScalarType> MeshBox;

ScalarType AVGL = 0;
ScalarType BBDiag = 0;
ScalarType MeshArea = 0;
std::vector<Geo::Point3<ScalarType>> VertCol;
std::vector<Geo::Point3<ScalarType>> FaceCol;

void InitConstantColor() {
  VertCol.clear();
  VertCol.resize(VertPos.size(), Geo::Point3<ScalarType>(0.9, 0.9, 0.9));
  FaceCol.clear();
}

// USER LINES: drawn holding SPACE, embedded in the input as features.
// The input as loaded is the target of the reprojection; the lines are kept
// as position pairs so they survive a new feature detection.
GLDraw::GLPolylinePick<ScalarType> LinePick;
// drawing: SPACE and the left button held together
bool SpaceDown = false;
bool LeftDown = false;
Geo::EmbedPolylineParam<ScalarType> EmbedP;
std::vector<typename Geo::Point3<ScalarType>> VertPosOrig;
std::vector<std::vector<int>> ConnectivityOrig;
typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosEdge;
std::set<PosEdge> UserEdgesPos;
std::set<Geo::Point3<ScalarType>> UserCornersPos;
std::vector<std::pair<int, int>> RenderUserEdges; // on the input view
std::string LineMsg; // result of the last embedding

PosEdge MakePosEdge(const Geo::Point3<ScalarType> &P0,
                    const Geo::Point3<ScalarType> &P1) {
  return PosEdge(std::min(P0, P1), std::max(P0, P1));
}

// PARAMETERS
Geo::TriRemParam<ScalarType> RemP;

// FEATURES (the angle is RemP.feature_angle, as in the command line tool)
bool use_features = true;
bool clean_features = true; // erode/dilate the feature lines
int ero_steps = 3;
int dil_steps = 3;

// SIZE
enum SizeMode { SizeAuto, SizeDiagRatio, SizeAvgEdge };
int size_mode = SizeAuto;
ScalarType ratio = 1;            // SizeAuto: scale of the automatic size (-res)
ScalarType diag_ratio = 0.01;    // SizeDiagRatio: edge / bbox diagonal
ScalarType avg_edge_ratio = 1;   // SizeAvgEdge: edge / input average edge
ScalarType MinSize = 0.005;      // SizeAuto: bounds, fractions of the diagonal
ScalarType MaxSize = 0.02;
// MAX FACES: bound on the expected output faces (0 disables)
int MaxFaces = 2000000;
// MAX FACES: hard bound on the faces during refinement
#define TRI_REMESH_HARD_MAX_FACES 5000000

float max_error_perc = 1.f;   // max geometric error, % of the diagonal
bool use_early_stop = false;
float early_stop_ratio = 0.001f;

std::vector<std::pair<int, int>> RenderVertFeatures;
std::vector<int> RenderVertCorners;

// STATS of the last remesh (shown in the panel)
std::vector<std::string> StatLines;
std::string SizeInfo;

void UpdateRenderCorners(const std::vector<std::vector<bool>> &IsFaceEFeature,
                         const std::vector<bool> &IsVertCorner) {
  RenderVertFeatures.clear();
  RenderVertCorners.clear();
  if ((IsFaceEFeature.size() != Connectivity.size()) ||
      (IsVertCorner.size() != VertPos.size()))
    return;
  Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
      Connectivity, IsFaceEFeature, RenderVertFeatures);
  Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(IsVertCorner,
                                                            RenderVertCorners);
}

// copy the input or the result into the displayed mesh
void ShowView(int view) {
  viewResult = (view == 1 && hasResult) ? 1 : 0;
  if (viewResult == 1) {
    VertPos = VertPosR;
    Connectivity = ConnectivityR;
  } else {
    VertPos = VertPos0;
    Connectivity = Connectivity0;
  }
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  InitConstantColor();
  if (viewResult == 1)
    UpdateRenderCorners(IsFaceEFeatureR, IsVertCornerR);
  else
    UpdateRenderCorners(IsFaceEFeature0, IsVertCorner0);
}

// features of the input with the current parameters
void UpdateFeatureCorners() {
  if (!loadedMesh)
    return;
  if (use_features && (RemP.feature_angle > 0) && (RemP.feature_angle < 180)) {
    if (clean_features)
      Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
          VertPos0, Connectivity0, RemP.feature_angle, IsFaceEFeature0,
          IsVertCorner0, ero_steps, dil_steps);
    else
      Geo::MeshFeatures<ScalarType>::FeatureCorners(
          VertPos0, Connectivity0, RemP.feature_angle, IsFaceEFeature0,
          IsVertCorner0);
  } else {
    IsVertCorner0 = std::vector<bool>(VertPos0.size(), false);
    IsFaceEFeature0 = std::vector<std::vector<bool>>(
        Connectivity0.size(), std::vector<bool>(3, false));
  }
  // USER LINES: always features (and their extremes corners)
  RenderUserEdges.clear();
  for (size_t i = 0; i < Connectivity0.size(); i++)
    for (size_t j = 0; j < 3; j++) {
      int IndexV0 = Connectivity0[i][j];
      int IndexV1 = Connectivity0[i][(j + 1) % 3];
      if (UserEdgesPos.count(MakePosEdge(VertPos0[IndexV0], VertPos0[IndexV1])) ==
          0)
        continue;
      IsFaceEFeature0[i][j] = true;
      if (IndexV0 < IndexV1)
        RenderUserEdges.push_back(std::make_pair(IndexV0, IndexV1));
    }
  for (size_t i = 0; i < VertPos0.size(); i++)
    if (UserCornersPos.count(VertPos0[i]) > 0)
      IsVertCorner0[i] = true;
  SizeInfo.clear(); // the automatic size depends on the features
  ShowView(0);
}

// USER LINES: embed a drawn line in the input (reprojected on the input as
// loaded); its edges become features for the next remesh
void EmbedUserLine(const std::vector<Geo::Point3<ScalarType>> &Line) {
  if (!loadedMesh || (Line.size() < 2))
    return;
  size_t NumV0 = VertPos0.size();
  std::vector<bool> WasCorner = IsVertCorner0;
  Geo::MeshEmbedPolyline<ScalarType> Embedder(VertPos0, Connectivity0,
                                              IsFaceEFeature0, IsVertCorner0,
                                              VertPosOrig, ConnectivityOrig);
  std::vector<std::pair<int, int>> NewEdges;
  bool done = Embedder.Embed(Line, EmbedP, NewEdges);
  if (!done) {
    LineMsg = "line not embedded (too short or off the mesh)";
    return;
  }
  for (size_t i = 0; i < NewEdges.size(); i++)
    UserEdgesPos.insert(
        MakePosEdge(VertPos0[NewEdges[i].first], VertPos0[NewEdges[i].second]));
  for (size_t i = 0; i < IsVertCorner0.size(); i++)
    if (IsVertCorner0[i] && ((i >= WasCorner.size()) || !WasCorner[i]))
      UserCornersPos.insert(VertPos0[i]);
  char buf[128];
  snprintf(buf, sizeof(buf), "line embedded: %zu edges, %zu new vertices",
           NewEdges.size(), VertPos0.size() - NumV0);
  LineMsg = buf;
  AVGL = AvgEdgeLen(VertPos0, Connectivity0);
  MeshArea = Area(VertPos0, Connectivity0);
  // features again: the detection on the refined input plus the user lines
  UpdateFeatureCorners();
}

bool IsSTLFile(const std::string &path) {
  if (path.size() < 4)
    return false;
  std::string ext = path.substr(path.size() - 4);
  for (size_t i = 0; i < ext.size(); i++)
    ext[i] = (char)std::tolower((unsigned char)ext[i]);
  return (ext == std::string(".stl"));
}

void LoadMesh(const std::string &path) {
  // loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  loadedMesh = IsSTLFile(path) ? LoadSTL(path, VertPos0, Connectivity0, true)
                               : LoadOBJ(path, VertPos0, Connectivity0, true);
  if (!loadedMesh || Connectivity0.empty()) {
    loadedMesh = false;
    std::cout << "Mesh not loaded" << std::endl;
    return;
  }
  if (!IsAllTris(Connectivity0)) {
    std::vector<typename Geo::Point3<ScalarType>> VertPosTris;
    std::vector<std::vector<int>> ConnectivityTris;
    SplitIntoTris(VertPos0, Connectivity0, VertPosTris, ConnectivityTris);
    VertPos0 = VertPosTris;
    Connectivity0 = ConnectivityTris;
  }
  std::cout << "Loaded " << Connectivity0.size() << " faces " << VertPos0.size()
            << " vertices" << std::endl;
  // the remesher needs a manifold mesh
  LoadMsg.clear();
  Geo::CleanMesh<ScalarType>::SplitNonManifold(VertPos0, Connectivity0);
  if (!MatchTriRemeshPreconditions(VertPos0, Connectivity0)) {
    bool ManifE = EdgeManifold<ScalarType>(VertPos0, Connectivity0);
    bool ManifV = VertManifold<ScalarType>(VertPos0, Connectivity0);
    LoadMsg = std::string("not manifold (edge manifold ") +
              std::to_string((int)ManifE) + ", vertex manifold " +
              std::to_string((int)ManifV) + ")";
    std::cout << "Warning: " << LoadMsg << std::endl;
  }

  // USER LINES: the input as loaded, target of the reprojection
  VertPosOrig = VertPos0;
  ConnectivityOrig = Connectivity0;
  UserEdgesPos.clear();
  UserCornersPos.clear();

  MeshBox.Init(VertPos0);
  AVGL = AvgEdgeLen(VertPos0, Connectivity0);
  BBDiag = MeshBox.Diag();
  MeshArea = Area(VertPos0, Connectivity0);
  hasResult = false;
  StatLines.clear();
  SizeInfo.clear();
  UpdateFeatureCorners();
}

// expected faces for a target length (area / equilateral triangle area)
ScalarType ExpectedFaces(ScalarType TargetL) {
  if (!(TargetL > 0))
    return 0;
  return MeshArea / ((std::sqrt((ScalarType)3) / 4) * TargetL * TargetL);
}

// TargetL (and the adaptive field) for the input, as the command line tool;
// the features must be in RP. False (reason in Msg) if no valid size.
bool SetupSize(Geo::TriRemParam<ScalarType> &RP, std::string &Msg) {
  Msg.clear();
  RP.VertSize.clear();
  if (size_mode == SizeAuto) {
    // FIND REMESHING SIZE
    ScalarType MinAggrSize = BBDiag * MinSize;
    ScalarType MaxAggrSize = BBDiag * MaxSize;
    ScalarType CurrSize = -1;
    if (RP.adaptive)
      CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeAggregated(
          VertPos0, Connectivity0, RP.IsFaceEFeature, MinAggrSize, MaxAggrSize,
          &RP.VertSize);
    else
      CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeAggregated(
          VertPos0, Connectivity0, RP.IsFaceEFeature, MinAggrSize, MaxAggrSize,
          NULL);
    if (!(CurrSize > 0)) {
      Msg = "invalid remesh size";
      return false;
    }
    if (RP.adaptive)
      Geo::RemeshSize<ScalarType>::AdaptiveSizeField(
          VertPos0, Connectivity0, CurrSize, MinAggrSize, MaxAggrSize,
          RP.VertSize);
    CurrSize /= 2;
    RP.TargetL = ratio * CurrSize;
    if (RP.adaptive)
      for (size_t i = 0; i < RP.VertSize.size(); i++)
        RP.VertSize[i] = ratio * RP.VertSize[i] / 2;
  } else {
    // the adaptive field needs the automatic size
    RP.adaptive = false;
    if (size_mode == SizeDiagRatio)
      RP.TargetL = diag_ratio * BBDiag;
    else
      RP.TargetL = avg_edge_ratio * AVGL;
  }

  // MAX FACES: enlarge TargetL if more faces are expected
  if (MaxFaces > 0) {
    ScalarType ExpFaces = ExpectedFaces(RP.TargetL);
    if (ExpFaces > (ScalarType)MaxFaces) {
      ScalarType Scale = std::sqrt(ExpFaces / (ScalarType)MaxFaces);
      RP.TargetL *= Scale;
      for (size_t i = 0; i < RP.VertSize.size(); i++)
        RP.VertSize[i] *= Scale;
      Msg = std::string("expected ") + std::to_string((long long)ExpFaces) +
            " faces: edge enlarged x" + std::to_string(Scale);
    }
  }
  return true;
}

void SaveResult() {
  if (!hasResult)
    return;
  std::string Folder = GetDirectory(PathMesh);
  if (Folder.empty())
    Folder = std::string(".");
  std::string OutFile = Folder + std::string("/") +
                        GetFilanameNoExtension(PathMesh) +
                        std::string("_remesh.obj");
  WriteOBJ(OutFile, VertPosR, ConnectivityR);
  StatLines.push_back(std::string("saved ") + OutFile);
  std::cout << "Saved " << OutFile << std::endl;
}

void DoRemesh() {
  if (!loadedMesh)
    return;
  StatLines.clear();
  if (!LoadMsg.empty()) {
    StatLines.push_back(std::string("cannot remesh: ") + LoadMsg);
    return;
  }
  auto start = std::chrono::high_resolution_clock::now();

  // always from the input; RemP keeps only the panel values (SetupSize may
  // change some of them, e.g. adaptive off for the fixed sizes)
  Geo::TriRemParam<ScalarType> RP = RemP;
  VertPosR = VertPos0;
  ConnectivityR = Connectivity0;
  RP.IsFaceEFeature = IsFaceEFeature0;
  RP.IsVertCorner = IsVertCorner0;
  // without features the remesher must not detect new ones either
  if (!use_features)
    RP.feature_angle = 180;
  RP.max_error = BBDiag * max_error_perc / 100;
  RP.early_stop_ratio = use_early_stop ? early_stop_ratio : 0;
  RP.writeMSG = true;

  std::string Msg;
  if (!SetupSize(RP, Msg)) {
    StatLines.push_back(std::string("error: ") + Msg);
    return;
  }
  if (!Msg.empty())
    StatLines.push_back(Msg);

  RP.max_faces = TRI_REMESH_HARD_MAX_FACES; // MAX FACES
  RP.check_manifold = true; // SAFETY: stop if a step breaks manifoldness
  int res = Geo::TriRemesh<ScalarType>(VertPosR, ConnectivityR, RP);

  // CLEANING with the features, always (as the command line tool)
  if (res >= 0) {
    bool FeaturesOK = (RP.IsFaceEFeature.size() == ConnectivityR.size()) &&
                      (RP.IsVertCorner.size() == VertPosR.size());
    Geo::CleanMesh<ScalarType>::CollapseBadShapeFaces(
        VertPosR, ConnectivityR, 0.05, FeaturesOK ? &RP.IsFaceEFeature : NULL,
        FeaturesOK ? &RP.IsVertCorner : NULL);
    FeaturesOK = FeaturesOK &&
                 (RP.IsFaceEFeature.size() == ConnectivityR.size()) &&
                 (RP.IsVertCorner.size() == VertPosR.size());
    // fixed vertices for unfolding: corners and feature-line vertices
    std::vector<bool> IsFixedV(VertPosR.size(), false);
    if (FeaturesOK) {
      std::vector<std::pair<int, int>> FeatPairs;
      Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
          ConnectivityR, RP.IsFaceEFeature, FeatPairs);
      for (size_t i = 0; i < FeatPairs.size(); i++) {
        IsFixedV[FeatPairs[i].first] = true;
        IsFixedV[FeatPairs[i].second] = true;
      }
      for (size_t i = 0; i < VertPosR.size(); i++)
        if (RP.IsVertCorner[i])
          IsFixedV[i] = true;
    }
    Geo::CleanMesh<ScalarType>::SolveFolds(VertPosR, ConnectivityR, -0.99,
                                           FeaturesOK ? &IsFixedV : NULL);
  }
  auto end = std::chrono::high_resolution_clock::now();
  ScalarType secs =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start)
          .count() /
      (ScalarType)1000000;
  if (res < 0) {
    StatLines.push_back("error: remesh stopped (see the console)");
    hasResult = false;
    ShowView(0);
    return;
  }
  IsFaceEFeatureR = RP.IsFaceEFeature;
  IsVertCornerR = RP.IsVertCorner;
  hasResult = true;
  ShowView(1);

  // STATS
  char buf[256];
  snprintf(buf, sizeof(buf), "time %.2f s, %zu faces, %zu vertices", secs,
           ConnectivityR.size(), VertPosR.size());
  StatLines.push_back(buf);
  ScalarType minLen, maxLen;
  ExtremenEdgeLen(VertPosR, ConnectivityR, minLen, maxLen);
  ScalarType optLen = AvgEdgeLen(VertPosR, ConnectivityR);
  snprintf(buf, sizeof(buf),
           "edge: target %.3g%% diag, avg %.3g%% min %.3g%% max %.3g%% of target",
           100 * RP.TargetL / BBDiag, 100 * optLen / RP.TargetL,
           100 * minLen / RP.TargetL, 100 * maxLen / RP.TargetL);
  StatLines.push_back(buf);
  ScalarType AvgDist, MaxDist;
  Geo::Hausdorff(VertPos0, Connectivity0, VertPosR, ConnectivityR, AvgDist,
                 MaxDist);
  snprintf(buf, sizeof(buf), "Hausdorff: avg %.3g%% max %.3g%% of diag",
           100 * AvgDist / BBDiag, 100 * MaxDist / BBDiag);
  StatLines.push_back(buf);
  ScalarType WorstQ =
      Geo::TriMeshWorstQuality<ScalarType>(VertPosR, ConnectivityR);
  ScalarType AvgQ = Geo::TriMeshAvgQuality<ScalarType>(VertPosR, ConnectivityR);
  snprintf(buf, sizeof(buf), "quality: avg %.3f worst %.3f", AvgQ, WorstQ);
  StatLines.push_back(buf);
  for (size_t i = 0; i < StatLines.size(); i++)
    std::cout << StatLines[i] << std::endl;
}

// expected size with the current parameters (the automatic size needs the
// size estimation, so it is computed on request)
void EstimateSize() {
  if (!loadedMesh)
    return;
  Geo::TriRemParam<ScalarType> TestP = RemP;
  TestP.IsFaceEFeature = IsFaceEFeature0;
  TestP.IsVertCorner = IsVertCorner0;
  std::string Msg;
  if (!SetupSize(TestP, Msg)) {
    SizeInfo = std::string("error: ") + Msg;
    return;
  }
  char buf[256];
  snprintf(buf, sizeof(buf), "edge %.3g%% of diag, ~%lld faces%s",
           100 * TestP.TargetL / BBDiag,
           (long long)ExpectedFaces(TestP.TargetL),
           Msg.empty() ? "" : " (capped)");
  SizeInfo = buf;
}

static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

// cursor in framebuffer pixels, y from below (as the GL viewport)
void CursorPixels(GLFWwindow *window, int &x, int &y) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  double xpos, ypos;
  glfwGetCursorPos(window, &xpos, &ypos);
  x = (int)(xpos * xscale);
  y = (int)(display_h - ypos * yscale);
}

static void cursor_position_callback(GLFWwindow *window, double xpos,
                                     double ypos) {
  // USER LINES: while drawing the mouse draws, the view does not move
  if (LinePick.IsActive()) {
    int x, y;
    CursorPixels(window, x, y);
    LinePick.AddScreenPoint(x, y);
    return;
  }
  if (ImGui::GetIO().WantCaptureMouse)
    return;

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  xpos *= xscale;
  ypos *= yscale;
  ypos = (display_h)-ypos;

  trackball.MouseMove((int)xpos, (int)ypos);
}

void mouse_button_callback(GLFWwindow *window, int button, int action,
                           int mods) {
  // USER LINES: the left button ends the line wherever it is released
  if ((button == GLFW_MOUSE_BUTTON_LEFT) && (action == GLFW_RELEASE)) {
    LeftDown = false;
    if (LinePick.IsActive()) {
      LinePick.End();
      return;
    }
  }
  if (ImGui::GetIO().WantCaptureMouse)
    return;
  // USER LINES: SPACE held + left button: start drawing (no rotation)
  if ((button == GLFW_MOUSE_BUTTON_LEFT) && (action == GLFW_PRESS)) {
    LeftDown = true;
    if (SpaceDown && loadedMesh) {
      LinePick.Begin();
      int x, y;
      CursorPixels(window, x, y);
      LinePick.AddScreenPoint(x, y);
      return;
    }
  }
  if (LinePick.IsActive())
    return;

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  double xpos, ypos;
  glfwGetCursorPos(window, &xpos, &ypos);
  xpos *= xscale;
  ypos *= yscale;
  ypos = (display_h)-ypos;

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {

    trackball.MouseDown((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
  }
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {

    trackball.MouseUp((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
  }
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
  if (ImGui::GetIO().WantCaptureMouse)
    return;
  trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods) {

  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
    trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);

  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
    trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);

  // USER LINES: SPACE + left button held to draw on the mesh (releasing
  // either one ends the line)
  if ((key == GLFW_KEY_SPACE) && (action == GLFW_PRESS) &&
      (!ImGui::GetIO().WantCaptureKeyboard))
    SpaceDown = true;
  if ((key == GLFW_KEY_SPACE) && (action == GLFW_RELEASE)) {
    SpaceDown = false;
    LinePick.End();
  }
}

void InitGLFW_Window() {
  if (!glfwInit()) {
    std::cout << "Error Initializing GLFW" << std::endl;
    exit(0);
  }
  // return 1;

  // window = glfwCreateWindow(1280, 720, "UV Mapping", NULL, NULL);
  window = glfwCreateWindow(1280, 720, "Triangle Remesher", NULL, NULL);

  if (window == NULL) {
    std::cout << "Error Initializing GLFW Window" << std::endl;
    exit(0);
  }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  // Setup window
  glfwSetErrorCallback(glfw_error_callback);
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetKeyCallback(window, key_callback);
}

void InitIMGui() {
  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  (void)io;
  // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable
  // Keyboard Controls

  // Setup Dear ImGui style
  ImGui::StyleColorsClassic();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL2_Init();
}

// a hover tooltip on the last item
void Help(const char *text) {
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", text);
}

// RENDER window: what is shown and how
void SetRenderWindow() {
  // first time on the right side, not over the remesher window
  ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, 10),
                          ImGuiCond_FirstUseEver);
  ImGui::Begin("Render", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  const float W = 140.0f;
  int view = viewResult;
  if (ImGui::RadioButton("input", &view, 0))
    ShowView(0);
  if (hasResult) {
    ImGui::SameLine();
    if (ImGui::RadioButton("result", &view, 1))
      ShowView(1);
  }
  ImGui::Checkbox("mesh", &showMesh);
  ImGui::SameLine();
  ImGui::Checkbox("features", &showFeatures);
  ImGui::SameLine();
  ImGui::Checkbox("corners", &showCorners);
  ImGui::SetNextItemWidth(W);
  ImGui::SliderFloat("line width", &line_thick, 0.5f, 5.f, "%.1f");
  ImGui::End();
}

// DRAW LINES window: lines drawn on the mesh, embedded as features
void SetDrawLinesWindow() {
  // first time below the Render window
  ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, 200),
                          ImGuiCond_FirstUseEver);
  ImGui::Begin("Draw lines", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  const float W = 140.0f;
  ImGui::TextWrapped("hold SPACE and the left button and move the mouse on "
                     "the mesh: the line is embedded in the input and kept "
                     "as a feature (corners at its extremes and crossings)");
  ImGui::SetNextItemWidth(W);
  ImGui::SliderInt("line smoothing", &EmbedP.smooth_steps, 0, 50);
  Help("smoothing steps of the drawn line, each reprojected on the input");
  float damp = (float)EmbedP.smooth_damp;
  ImGui::SetNextItemWidth(W);
  ImGui::SliderFloat("line smoothing strength", &damp, 0.05f, 1.f, "%.2f");
  Help("how much each smoothing step moves the line (Laplacian damping)");
  EmbedP.smooth_damp = damp;
  ImGui::SetNextItemWidth(W);
  ImGui::SliderInt("smooth line on the mesh", &EmbedP.relax_steps, 0, 50);
  Help("after the embedding: the line vertices smoothed along the line and "
       "reprojected, the mesh around following (removes the zigzag)");
  ImGui::SetNextItemWidth(W);
  ImGui::SliderInt("smooth around (rings)", &EmbedP.smooth_rings, 0, 10);
  Help("after the embedding, the mesh around the line is smoothed and "
       "reprojected, the line fixed: rings of vertices (0 = off)");
  ImGui::SetNextItemWidth(W);
  ImGui::SliderInt("smooth around steps", &EmbedP.smooth_around_steps, 0,
                   20);
  float sr = (float)(100 * EmbedP.sample_ratio);
  ImGui::SetNextItemWidth(W);
  ImGui::SliderFloat("line sampling (% of diag)", &sr, 0.1f, 5.f, "%.2f",
                     ImGuiSliderFlags_Logarithmic);
  Help("distance between the line vertices (at least 10 per line); the "
       "faces under the line are refined down to about twice this");
  EmbedP.sample_ratio = sr / 100;
  bool refine = (EmbedP.refine_steps > 0);
  ImGui::Checkbox("refine under the line", &refine);
  EmbedP.refine_steps = refine ? 20 : 0;
  if (ImGui::Button("clear lines")) {
    // the vertices stay, they are no longer features
    UserEdgesPos.clear();
    UserCornersPos.clear();
    UpdateFeatureCorners();
  }
  ImGui::SameLine();
  ImGui::Text("%zu line edges", UserEdgesPos.size());
  if (!LineMsg.empty())
    ImGui::TextWrapped("%s", LineMsg.c_str());
  ImGui::End();
}

void SetRenderBar() {
  // ImGui::Begin("Mesh Optimizer", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  ImGui::Begin("Triangle Remesher", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  const float W = 140.0f;

  // ---- MESH
  if (loadedMesh) {
    ImGui::Text("%s", GetFilaname(PathMesh).c_str());
    ImGui::Text("%zu faces, %zu vertices", Connectivity0.size(),
                VertPos0.size());
    ImGui::Text("avg edge %.3g%% of diag", 100 * AVGL / BBDiag);
    if (!LoadMsg.empty())
      ImGui::TextColored(ImVec4(0.8f, 0, 0, 1), "cannot remesh: %s",
                         LoadMsg.c_str());
  } else
    ImGui::Text("no mesh (pass a .obj or .stl on the command line)");

  // ---- FEATURES (shown live on the input)
  if (ImGui::CollapsingHeader("Features", ImGuiTreeNodeFlags_DefaultOpen)) {
    bool changed = ImGui::Checkbox("preserve sharp features", &use_features);
    if (use_features) {
      int angle = (int)floor(0.5 + RemP.feature_angle);
      ImGui::SetNextItemWidth(W);
      changed |= ImGui::SliderInt("sharp edge angle (deg)", &angle, 1, 179);
      Help("edges where the surface bends more than this are kept sharp");
      RemP.feature_angle = angle;
      changed |= ImGui::Checkbox("remove tiny feature fragments",
                                 &clean_features);
      Help("erode/dilate the feature lines (as the command line tool)");
    }
    if (changed)
      UpdateFeatureCorners();
    ImGui::Text("%zu feature edges, %zu corners (on the view)",
                RenderVertFeatures.size(), RenderVertCorners.size());
  }

  // ---- SIZE
  if (ImGui::CollapsingHeader("Target size", ImGuiTreeNodeFlags_DefaultOpen)) {
    int mode = size_mode;
    bool changed = false;
    changed |= ImGui::RadioButton("automatic", &mode, (int)SizeAuto);
    Help("from features and curvature, as the command line tool");
    ImGui::SameLine();
    changed |= ImGui::RadioButton("% of diagonal", &mode, (int)SizeDiagRatio);
    ImGui::SameLine();
    changed |= ImGui::RadioButton("x input edge", &mode, (int)SizeAvgEdge);
    size_mode = mode;

    if (size_mode == SizeAuto) {
      float s = (float)ratio;
      ImGui::SetNextItemWidth(W);
      changed |= ImGui::SliderFloat("scale", &s, 0.25f, 4.f, "%.2f",
                                    ImGuiSliderFlags_Logarithmic);
      Help("< 1 finer, > 1 coarser (the -res option)");
      ratio = s;
      changed |= ImGui::Checkbox("adaptive (finer on features and curvature)",
                                 &RemP.adaptive);
    } else if (size_mode == SizeDiagRatio) {
      float p = (float)(100 * diag_ratio);
      ImGui::SetNextItemWidth(W);
      changed |= ImGui::SliderFloat("edge (% of diag)", &p, 0.05f, 10.f, "%.2f",
                                    ImGuiSliderFlags_Logarithmic);
      diag_ratio = p / 100;
    } else {
      float r = (float)avg_edge_ratio;
      ImGui::SetNextItemWidth(W);
      changed |= ImGui::SliderFloat("edge (x input avg)", &r, 0.1f, 20.f, "%.2f",
                                    ImGuiSliderFlags_Logarithmic);
      avg_edge_ratio = r;
    }
    ImGui::SetNextItemWidth(W);
    changed |= ImGui::InputInt("max faces", &MaxFaces, 100000);
    Help("the edge is enlarged if more faces are expected (0 = no bound)");
    MaxFaces = std::max(MaxFaces, 0);

    // the fixed sizes are estimated live, the automatic one on request
    if (changed)
      SizeInfo.clear();
    if ((size_mode != SizeAuto) && loadedMesh && SizeInfo.empty())
      EstimateSize();
    if (size_mode == SizeAuto) {
      if (ImGui::Button("estimate"))
        EstimateSize();
      ImGui::SameLine();
    }
    ImGui::Text("%s", SizeInfo.c_str());
  }

  // ---- REMESHING
  if (ImGui::CollapsingHeader("Remeshing", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::SetNextItemWidth(W);
    ImGui::SliderInt("iterations", &RemP.steps, 1, 20);
    ImGui::Checkbox("linear sweep (faster)", &RemP.iterative);
    Help("operations in sequence instead of by priority (the -ite option)");
    ImGui::Checkbox("keep borders", &RemP.fix_borders);
    ImGui::Checkbox("early stop", &use_early_stop);
    Help("stop when an iteration changes less than ratio x faces");
    if (use_early_stop) {
      ImGui::SameLine();
      ImGui::SetNextItemWidth(80.0f);
      ImGui::InputFloat("ratio", &early_stop_ratio, 0, 0, "%.4f");
      early_stop_ratio = std::max(early_stop_ratio, 0.f);
    }
    ImGui::Checkbox("final smooth & projection",
                    &RemP.post_project);
    Help("final steps of smoothing, reprojection on the input and post "
         "projection (the input attracts the remeshed surface); feature lines "
         "reprojected on the input features, corners fixed");
    if (RemP.post_project) {
      ImGui::SameLine();
      ImGui::SetNextItemWidth(80.0f);
      ImGui::InputInt("steps", &RemP.post_project_steps);
      RemP.post_project_steps =
          std::max(RemP.post_project_steps, 1);
    }
    ImGui::SetNextItemWidth(W);
    ImGui::SliderFloat("max error (% of diag)", &max_error_perc, 0.05f, 5.f,
                       "%.2f", ImGuiSliderFlags_Logarithmic);
    Help("collapses and flips moving the surface more than this are refused");

    if (ImGui::Button("REMESH"))
      DoRemesh();
    ImGui::SameLine();
    if (ImGui::Button("SAVE"))
      SaveResult();
    Help("<input>_remesh.obj next to the input");
    for (size_t i = 0; i < StatLines.size(); i++)
      ImGui::TextWrapped("%s", StatLines[i].c_str());
  }

  ImGui::End();
}
void GLDrawMesh() {

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w,
               clear_color.z * clear_color.w, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(40, (GLdouble)display_w / (GLdouble)display_h, 0.1, 100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, 3.5f, 0, 0, 0, 0, 1, 0);

  trackball.GetView();
  trackball.Apply();
  trackball.DrawPostApply();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());

  bool DrawSurface = true;
  if (showMesh) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glLineWidth(line_thick);

    GLDraw::glColor(Geo::Point3<ScalarType>(0.8, 0.8, 0.8));

    if (FaceCol.size() == Connectivity.size())
      GLDraw::DrawSurfaceMesh<ScalarType>(
          VertPos, Connectivity, FaceNormals, VertNormals, DrawSurface, true,
          false, false, GLDraw::TRTriangle, &FaceCol, NULL); //,&FaceColor);
    else {
      if (VertCol.size() == VertPos.size())
        GLDraw::DrawSurfaceMesh<ScalarType>(
            VertPos, Connectivity, FaceNormals, VertNormals, DrawSurface, true,
            false, false, GLDraw::TRTriangle, NULL, &VertCol); //,&FaceColor);
    }

    glPopAttrib();
  }
  // USER LINES: unproject on the depth of the mesh just drawn
  LinePick.ProcessPending();
  LinePick.Draw();
  if (viewResult == 0) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    GLDraw::GLDrawEdges<ScalarType>(VertPos, RenderUserEdges, 12,
                                    Geo::Point3<ScalarType>(1, 0.5, 0));
    glPopAttrib();
  }
  if (showFeatures) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    GLDraw::GLDrawEdges<ScalarType>(VertPos, RenderVertFeatures, 10,
                                    Geo::Point3<ScalarType>(1, 0, 1));

    glPopAttrib();
  }
  if (showCorners) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    GLDraw::DrawVertices<ScalarType>(VertPos, RenderVertCorners,
                                     Geo::Point3<ScalarType>(0, 0, 0), 30);

    glPopAttrib();
  }
}

int main(int argc, char *argv[]) {

  // POST PROJECTION: on by default here (off in the remesher)
  RemP.post_project = true;

  if (argc == 2) {
    PathMesh = std::string(argv[1]);
    LoadMesh(PathMesh);
  }

  InitGLFW_Window();

  InitIMGui();

  trackball.center = Geo::Point3<float>(0, 0, 0);
  trackball.radius = 1;

  glewInit();
  // Main loop
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    {
      SetRenderBar();
      SetRenderWindow();
      SetDrawLinesWindow();
    }

    // Rendering
    ImGui::Render();

    GLDrawMesh();

    // USER LINES: a completed line is embedded in the input
    std::vector<Geo::Point3<ScalarType>> Line;
    if (LinePick.TakeFinished(Line))
      EmbedUserLine(Line);

    // If you are using this code with non-legacy OpenGL header/contexts (which
    // you should not, prefer using imgui_impl_opengl3.cpp!!), you may need to
    // backup/reset/restore other state, e.g. for current shader using the
    // commented lines below.
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

    glfwMakeContextCurrent(window);
    glfwSwapBuffers(window);
  }

  // Cleanup
  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}
