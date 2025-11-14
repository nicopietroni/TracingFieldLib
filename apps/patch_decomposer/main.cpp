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

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <GL/glew.h>
#include <UI/orbitball.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <sys/_types/_size_t.h>
#include <utility>
#include <vector>
#include <IO/io_path.h>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include "IO/file_path_functions.h"
#include "IO/io_obj.h"
#include <GLFW/glfw3.h>
#include <GL_functions/GLDrawField.h>
#include <GL_functions/GLDrawGraph.h>
#include <GL_functions/GLDrawLabel_imgui.h>
#include <GL_functions/GLDrawMesh.h>
#include <GL_functions/GLDrawPatchData.h>
#include <GL_functions/GLDrawSideDiscretization.h>
#include <box3.h>
#include <point3.h>

#include "IO/io_features.h"
#include "IO/io_field.h"
#include "tangent_space_smooth.h"
#include <GL_functions/GLDrawOutline.h>
#include <chrono>
#include <connectivity.h>
#include <cross_field_functions.h>
#include <field_graph/UV_decomposition_conditions.h>
#include <field_graph/dual_patches.h>
#include <field_graph/graph.h>
#include <field_graph/graph_tracing.h>
#include <field_graph/initialize_graph.h>
#include <field_graph/mesh_preprocess.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_managing.h>
#include <field_graph/patch_optimize.h>
#include <field_graph/path_functions.h>
#include <field_graph/path_sampling.h>

#include <field_graph/t_junction.h>
#include <field_graph/tracing_interface.h>
#include <math_utils.h>
#include <mesh_features.h>
#include <mesh_functions.h>
#include <mesh_mapping.h>
#include <patch_side_discretization.h>
#include <quadrangulation.h>
#include <quadrangulationLUT.h>
#include <rectify_features.h>
#include <remesh_size.h>
#include <smooth_mesh_sampling_dir.h>
#include <smooth_with_patches.h>
#include <quadmesh_smooth_separatrix.h>
// #include <motorcycle_graph.h>

GLFWwindow *window;
std::string PathMesh;
std::string PathField;
std::string ProjectName;
UI::OrbitBall trackball;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
// ImVec4 mesh_color = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
float line_thick = 1.f;
int selected_rend_item = 5;

bool DrawInitalPos = false;
bool DrawCandidates = false;
int curr_cand = -1;

bool DrawSamplesPaths = false;

bool showMesh = true;
bool showFeatures = true;
// bool showPaths = true;
int DrawPathMode = 0;
int CurrSelP = -1;
bool showDual = true;
bool showUnsolved = false;
bool showCross = true;
bool showSing = true;
bool has_graph = false;
// bool has_init_patches=false;

// TangentDomain,
bool DrawTangentNode = false;
// ExitDomain,
bool DrawExitNode = false;
// Internal,
bool DrawInternalNode = false;
// TangentExitDomain,
bool DrawTangentExit = false;
// ConvexExitDomain
bool DrawConvexExit = false;
// DartOrthoDomain
bool DrawDartOrthoDomain = false;
// JumpNodes
bool DrawJumpNodes = false;

bool DrawCorners = false;
bool DrawDualCorners = false;
bool DrawTJunction = false;

bool DrawLoops = true;
bool DrawPolyL = true;

bool DrawSides = false;

bool Draw_must_pass_on = false;
bool Draw_sample_on = false;
bool Draw_force_sample = false;
bool Draw_allow_path_on = false;
bool Draw_allow_arrive_to = false;

// bool DrawAssignTJuctions = false;

bool DrawSubdivision = false;
bool DrawAssignValues = false;
bool DrawQuadMesh = false;

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

typedef double ScalarType;

// DUAL MESH DATA
std::vector<Geo::Point3<ScalarType>> DualVertPos;
std::vector<std::vector<int>> DualConnectivity;
std::vector<Geo::Point3<ScalarType>> DualFaceNormals;
std::vector<Geo::Point3<ScalarType>> DualVertNormals;
std::vector<std::vector<int>> DualNextF, DualNextE;

// MESH DATA
bool loadedMesh = false;
std::vector<Geo::Point3<ScalarType>> VertPos;
std::vector<std::vector<int>> Connectivity;
std::vector<Geo::Point3<ScalarType>> FaceNormals;
std::vector<Geo::Point3<ScalarType>> VertNormals;
std::vector<std::vector<int>> NextF, NextE;
Geo::Box3<ScalarType> MeshBox;

// FIELD DATA
std::vector<Field::CrossF<ScalarType>> VertCurv;
std::vector<Field::CrossF<ScalarType>> FaceCurv;
bool has_cross_field = false;
bool has_dual_mesh = false;
ScalarType minQCrossVert, maxQCrossVert;
ScalarType minQCrossFace, maxQCrossFace;

std::vector<int> SingIndex, SingValue;
std::vector<int> FeatureSingIndex, FeatureSingValue;

Geo::TracingGraph<ScalarType> TGraph(Geo::GT_CrossField);

float drift_penalty = 100.f;
float maxAngle = 45.f;
int diffuse_step = 1;
int smooth_steps = 10;
int CurrN = 0;
int OldN = -1;
// std::vector<int> PathN;

bool OptimizeSize = true;

Geo::CrossPathSampling<ScalarType> PathSampl(VertPos, Connectivity, TGraph,
                                             VertCurv, SingIndex, SingValue);

std::vector<Geo::Point3<ScalarType>> PathColors;

std::vector<Geo::Point3<ScalarType>> FaceCol;
std::vector<Geo::Point3<ScalarType>> FaceColDual;
Geo::PatchManaging<ScalarType> PatchM(VertPos, Connectivity, VertCurv, TGraph,
                                      PathSampl.SelectedPaths);

int CurrCol = 0;
int OldCol = CurrCol;

Geo::Point3<ScalarType> mesh_color(0.79, 0.9, 0.87);

Geo::PatchDecomposer<ScalarType> PDeco(PatchM, PathSampl, VertCurv, FaceCurv);

std::vector<Geo::Point3<ScalarType>> SidesMidP;
// for debugging the marked parts
// std::vector<int> MarkedV;

int failed_trace = 0;
int tested_trace = 0;

bool single_boundary_cond = true;
bool self_connection_cond = true;
bool disk_like_cond = true;
bool single_sing_cond = false;
bool match_sing_cond = true;
bool cc_ability_cond = true;
bool equal_side_cond = false;
bool use_uv = false;
bool isometric_uv = false;
bool conformal_uv = false;
ScalarType max_uvDist = 0.05;

int MinSides = 3;
int MaxSides = 6;
ScalarType CCability = 0.5;
ScalarType EqSideThr = 1;

bool old_single_boundary_cond = single_boundary_cond;
bool old_self_connection_cond = self_connection_cond;
bool old_disk_like_cond = disk_like_cond;
bool old_single_sing_cond = single_sing_cond;
bool old_match_sing_cond = match_sing_cond;
bool old_cc_ability_cond = cc_ability_cond;
bool old_isometric_uv = isometric_uv;
bool old_conformal_uv = conformal_uv;
bool old_equal_side_cond = equal_side_cond;
ScalarType old_max_uvDist = max_uvDist;

int old_MinSides = MinSides;
int old_MaxSides = MaxSides;
ScalarType old_CCability = CCability;
ScalarType old_EqSideThr = EqSideThr;

int min_sing_distance = 1;

int DrawMeshMode = 1;
int DrawFieldMode = 0;

ScalarType AVLen = 0;

std::vector<std::pair<int, int>> Features;
std::vector<int> CornersIdx;

std::vector<int> ModSameBVert;
std::vector<Geo::Point3<ScalarType>> ModSameBOrigPos;

std::vector<Geo::TJunction> TJunctions;

ScalarType ThickOutline = 8;
ScalarType ThickBorders = 8;
ScalarType ThickPrimaryP = 8;
ScalarType ThickSecondaryP = 5;
ScalarType ThickTertiaryP = 2;

std::vector<int> MaxPathRank;

std::vector<std::pair<int, int>> RemainingF;
std::vector<int> NewBorder;

std::vector<Geo::Point3<ScalarType>> VisDualCornerPos;
// std::vector<std::vector<int>> PerDualPatchVertIDX;
std::vector<std::vector<int>> PerDualPatchFaceIDX;
std::vector<std::vector<int>> PerDualPatchCornerIDX;

std::vector<std::pair<int, int>> DualFeatures;
std::vector<int> DualCorners;

Geo::SideDiscretization<ScalarType> SideD;
// Geo::SideDiscretization<ScalarType>::parityMode
// PMode=Geo::SideDiscretization<ScalarType>::EvenVar;
//int ParityMode = 1;
ScalarType QuadSize = 0.02;

std::vector<Geo::Point3<ScalarType>> QuadVert;
std::vector<std::vector<int>> QuadFaces;
std::vector<int> QuadPatch;
std::vector<Geo::Point3<ScalarType>> QuadFaceNormals, QuadVertNormals;
std::vector<Geo::Point3<ScalarType>> QuadFaceColors;
std::vector<std::pair<int, int>> QuadFeatures;
std::vector<int> QuadCorners;

Geo::QuadrangulationLUT<ScalarType> Q_LUT;

int QuadSmoothSteps = 10;

ScalarType RegVsIso = 1.0;
ScalarType MinSize = 0.005;
ScalarType MaxSize = 0.05;

std::vector<std::vector<std::pair<int, int>>> Separatrices;

void AutoQuadEdgeSize() {
  //std::vector<ScalarType> VertSize;
  ScalarType MinAggrSize = MeshBox.Diag() * MinSize;
  ScalarType MaxAggrSize = MeshBox.Diag() * MaxSize;
  ScalarType CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeFeatures(VertPos, Connectivity, Features,
                                                                         MinAggrSize, MaxAggrSize);// VertSize);
  QuadSize = CurrSize / MeshBox.Diag();
  QuadSize/=2;
  std::cout<<"QuadSize: "<<QuadSize<<std::endl;
}

void RestoreOriginalPosOnSplitBorders() {
  assert(ModSameBVert.size() == ModSameBOrigPos.size());
  for (size_t i = 0; i < ModSameBVert.size(); i++) {
    int IndexV = ModSameBVert[i];
    assert(IndexV >= 0);
    assert(IndexV < VertPos.size());
    VertPos[IndexV] = ModSameBOrigPos[i];
  }
}


void UpdateColor() {
  bool has_uv = false;
  std::vector<ScalarType> FaceVal;

  switch (CurrCol) {
  case 0:
    FaceCol =
        std::vector<Geo::Point3<ScalarType>>(Connectivity.size(), mesh_color);
    break;
  case 1:
    PatchM.GetColorFaceByPatch(FaceCol);
    break;
  case 2:
    PatchM.GetColorFaceIfDiskLike(FaceCol);
    break;
  case 3:
    PatchM.GetColorByValence(FaceCol);
    break;
  case 4:
    PatchM.GetColorByMatchSingularity(FaceCol);
    break;
  case 5:
    PDeco.GetColorBySolvedPatches(FaceCol);
    break;
  case 6:
    has_uv =
        PDeco.GetConditionPerFaceScalar(FaceVal, std::string("uv_distortion"));
    if (has_uv) {
      Geo::GetColorByScalar<ScalarType>(FaceVal, max_uvDist, 0, FaceCol);
    }
    break;
  default:
    FaceCol =
        std::vector<Geo::Point3<ScalarType>>(Connectivity.size(), mesh_color);
    // code block
  }
}

void InitGraph() {
  Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(
      TGraph, VertPos, Connectivity, VertCurv, drift_penalty, maxAngle,
      diffuse_step);

  Geo::PathFunctions<ScalarType>::DisableNonDirectableConnections(TGraph);

  // TGraph.DeSelectAll(0);
  //  TGraph.SelectNoteType(Geo::ExitDomain, 0);
  //  TGraph.SelectNoteType(Geo::TangentDomain, 0);
}

bool UpdateConditions() {
  bool has_changed = false;
  has_changed |= old_single_boundary_cond != single_boundary_cond;
  has_changed |= old_self_connection_cond != self_connection_cond;
  has_changed |= old_disk_like_cond != disk_like_cond;
  has_changed |= old_single_sing_cond != single_sing_cond;
  has_changed |= old_match_sing_cond != match_sing_cond;
  has_changed |= old_cc_ability_cond != cc_ability_cond;
  has_changed |= old_equal_side_cond != equal_side_cond;

  has_changed |= old_MinSides != MinSides;
  has_changed |= old_MaxSides != MaxSides;
  has_changed |= old_CCability != CCability;
  has_changed |= old_EqSideThr != EqSideThr;
  has_changed |= isometric_uv != old_isometric_uv;
  has_changed |= conformal_uv != old_conformal_uv;
  has_changed |= max_uvDist != old_max_uvDist;

  old_single_boundary_cond = single_boundary_cond;
  old_self_connection_cond = self_connection_cond;
  old_disk_like_cond = disk_like_cond;
  old_single_sing_cond = single_sing_cond;
  old_match_sing_cond = match_sing_cond;
  old_cc_ability_cond = cc_ability_cond;

  old_MinSides = MinSides;
  old_MaxSides = MaxSides;
  old_CCability = CCability;
  old_EqSideThr = EqSideThr;
  old_isometric_uv = isometric_uv;
  old_conformal_uv = conformal_uv;
  old_max_uvDist = max_uvDist;
  old_equal_side_cond = equal_side_cond;

  return has_changed;
}

void InitConditions() {
  PDeco.ClearConditions();
  if (single_boundary_cond)
    PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
  if (disk_like_cond)
    PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
  if (self_connection_cond)
    PDeco.AddCondition(new Geo::SelfConnectingPatchCondition<ScalarType>());
  if (single_sing_cond)
    PDeco.AddCondition(new Geo::SingleSingularityCondition<ScalarType>());
  if (match_sing_cond)
    PDeco.AddCondition(new Geo::MathchSingularityCondition<ScalarType>());
  if ((MinSides > 0) || (MaxSides > 0)) {
    PDeco.AddCondition(new Geo::NumSidesCondition<ScalarType>());
  }
  // if (CCability >= 1) {
  if (cc_ability_cond) {
    // ScalarType AvVAl = AvgEdgeLen(VertPos, Connectivity) * CCability;
    PDeco.AddCondition(new Geo::CCabilityCondition<ScalarType>(CCability));
  }

  if (isometric_uv) {
    Geo::UVCondition<ScalarType>::UVType uv_type =
        Geo::UVCondition<ScalarType>::UVArap;
    Geo::UVDistortion<ScalarType>::DistType dist_type =
        Geo::UVDistortion<ScalarType>::EdgeDist;
    PDeco.AddCondition(
        new Geo::UVCondition<ScalarType>(uv_type, dist_type, max_uvDist));
  }
  if (conformal_uv) {
    Geo::UVCondition<ScalarType>::UVType uv_type =
        Geo::UVCondition<ScalarType>::UVLSqConformal;
    Geo::UVDistortion<ScalarType>::DistType dist_type =
        Geo::UVDistortion<ScalarType>::AngleDist;
    PDeco.AddCondition(
        new Geo::UVCondition<ScalarType>(uv_type, dist_type, max_uvDist));
  }
  if (equal_side_cond)
    PDeco.AddCondition(new Geo::EqualSideCondition<ScalarType>(EqSideThr));
}


void SampleTracing() {
  // PathN.clear();
  // PathSampl.Init();
  // if ((PDeco.start_spike) || (PDeco.start_concave))
  //   PathSampl.force_direct = true;

  // PDeco.SampleInitialCandidates();
  PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioDist;
  PDeco.InitSampling(false); // PDeco.add_only_where_needed);
  PathSampl.force_direct = true;
  // if (PDeco.add_only_where_needed)
  //   PDeco.InitMustPassEdges();

  PDeco.RetraceCandidates();
}

void SmoothPaths() {
  Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  std::vector<std::vector<int>> VertPaths;
  PathSampl.GetVertexSelectedPaths(VertPaths);
  Geo::PathFunctions<ScalarType>::FindTJunctions(VertPaths, TJunctions);

  // Geo::PatchOptimize<ScalarType>::TwoStepSmooth(VertPos, Connectivity,
  //                                               VertPaths);
  // std::cout<<"Features size:"<<Features.size()<<std::endl;
  // std::cout<<"VertPath size:"<<VertPaths.size()<<std::endl;
  // WriteFeatures("feature1.feat",Features);
  // WritePath("path1.path",VertPaths);

  Geo::PatchOptimize<ScalarType>::SmoothPaths(VertPos, Connectivity, VertPaths,
                                              Features);
  //   WriteOBJ("Primal1.obj",VertPos, Connectivity);

  // exit(0);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

  // update cross F poistion
  for (size_t i = 0; i < VertCurv.size(); i++)
    VertCurv[i].Pos = VertPos[i];

  AdaptFieldToMesh(VertPos, Connectivity, VertCurv, FaceCurv);

  PatchM.UpdateSubPatchPos();
  PatchM.InitVisualCornersEdges();
}

void InitPatches() {
  // initialize patches corners
  // and all remaining data

  PatchM.CheckConsistentData();
  PatchM.Init();

  PatchM.CheckConsistentData();

  // then initialize the stopping conditions
  InitConditions();
  PatchM.CheckConsistentData();

  // PDeco.InitBorderVertType();

  // initialise the decomposer
  PDeco.Init();
  PatchM.CheckConsistentData();

  // then check all patches and saved the unsolved ones
  PatchM.CheckConsistentData();
  PDeco.UpdateAllUnsolvedPatches();

  std::cout << "Size Remaining:" << RemainingF.size() << std::endl;
  PDeco.AddFeaturesAsPaths(TGraph, VertCurv, RemainingF, NewBorder);
  PathColors.resize(PathSampl.SelectedPaths.size());
  Geo::ScatterColor<ScalarType>(PathColors);
  UpdateColor();
}

bool LoadMesh(const std::string &path) {
  //    oss<<"Loading "<<path.c_str()<<std::endl;
  loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  if (loadedMesh) {
    std::cout << "Loaded " << Connectivity.size() << " faces " << VertPos.size()
              << " vertices" << std::endl;
    MeshBox.Init(VertPos);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);
    UpdateColor();
    return true;
  } else {
    std::cout << "Mesh not loaded" << std::endl;
    return false;
  }
}

bool LoadField(const std::string &path) {
  // std::cout<<"Loading Field"<<path.c_str()<<std::endl;
  has_cross_field = ReadField(VertPos, Connectivity, path, VertCurv, FaceCurv);
  if ((VertCurv.size() == 0) || (VertCurv.size() != VertPos.size())) {
    std::cout << "Init Vert Curv" << std::endl;
    SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv,
                         VertCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);
  }
  if (FaceCurv.size() == 0) {
    SetFaceCrossVectorFromVert(Connectivity, FaceNormals, VertCurv, FaceCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);
  }
  if (has_cross_field) {
    std::cout << "Loaded Field" << std::endl;
    showCross = true;
    showSing = true;

    // Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert,
    //                                       maxQCrossVert);
    Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert,
                                          maxQCrossVert);
    Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace,
                                          maxQCrossFace);

    GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE,
                         SingIndex, SingValue);

    PatchM.SetSingularities(SingIndex, SingValue);
    PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
    return true;
  } else {
    std::cout << "Field not loaded" << std::endl;
    return false;
  }
}

static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

static void cursor_position_callback(GLFWwindow *window, double xpos,
                                     double ypos) {
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
  (void)mods;
  if (ImGui::GetIO().WantCaptureMouse)
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

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    trackball.MouseDown((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    trackball.MouseUp((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
  (void)window;
  (void)xoffset;
  trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods) {
  (void)window;
  (void)scancode;
  (void)mods;

  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
    trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);
  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
    trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);
}

void InitGLFW_Window() {
  if (!glfwInit()) {
    // std::cout<<"Error Initializing GLFW"<<std::endl;
    exit(0);
  }
  // return 1;

  window = glfwCreateWindow(1280, 720, "Dear ImGui GLFW+OpenGL2 example", NULL,
                            NULL);
  // if (window == NULL)
  //     return 1;
  if (window == NULL) {
    // std::cout<<"Error Initializing GLFW Window"<<std::endl;
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

void SetRenderBar() {
  // GLGraphDraw::GLDrawPatchSideDiscretization<ScalarType>::DrawSideValue(SideD);

  ImGui::Begin(
      "Render", NULL,
      ImGuiWindowFlags_AlwaysAutoResize); // Create a window called "Hello,
                                          // world!" and append into it.

  // ImGui::Separator();
  // ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW MESH");
  if (ImGui::CollapsingHeader("DRAW MESH")) {
    ImGui::Checkbox("Draw Mesh", &showMesh);
    ImGui::Checkbox("Draw Features", &showFeatures);
    static const char *itemsDrawMode[] = {"Smooth", "Flat", "Smooth Wire",
                                          "Flat Wire"};
    ImGui::Checkbox("Draw Dual", &showDual);
    ImGui::SetNextItemWidth(100);

    ImGui::Checkbox("Draw Unsolved", &showUnsolved);
    ImGui::SetNextItemWidth(100);

    ImGui::Combo("Mesh Mode", &DrawMeshMode, itemsDrawMode,
                 IM_ARRAYSIZE(itemsDrawMode));

    ImGui::Checkbox("Draw Cross", &showCross);
    static const char *itemsFDrawMode[] = {"Face Field", "Vert Field"};
    ImGui::SetNextItemWidth(100);
    ImGui::Combo("Field Mode", &DrawFieldMode, itemsFDrawMode,
                 IM_ARRAYSIZE(itemsFDrawMode));

    ImGui::Checkbox("Draw Singularities", &showSing);
  }
  // ImGui::Checkbox("Draw Paths", &showPaths);
  // ImGui::Separator();
  // ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW PATHS");
  if (ImGui::CollapsingHeader("DRAW PATHS")) {
    if (ImGui::RadioButton("Draw Colored", DrawPathMode == 0)) {
      DrawPathMode = 0;
    }
    if (ImGui::RadioButton("Draw Importance", DrawPathMode == 1)) {
      DrawPathMode = 1;
    }
    if (ImGui::RadioButton("No Draw", DrawPathMode == 2)) {
      DrawPathMode = 2;
    }
    ImGui::InputInt("Curr Path", &CurrSelP);
    CurrSelP = std::min(CurrSelP, (int)PathSampl.SelectedPaths.size() - 1);
    ImGui::Checkbox("Draw Loops", &DrawLoops);
    ImGui::Checkbox("Draw PolyL", &DrawPolyL);
  }
  // ImGui::Separator();
  // ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW GRAPH");
  if (ImGui::CollapsingHeader("DRAW GRAPH")) {
    //  ImGui::Checkbox("Draw Exit Node", &DrawExitNode);
    ImGui::Checkbox("Draw Exit Node", &DrawExitNode);
    ImGui::Checkbox("Draw Tangent Node", &DrawTangentNode);
    ImGui::Checkbox("Draw Internal Node", &DrawInternalNode);
    ImGui::Checkbox("Draw Tangent Exit", &DrawTangentExit);
    ImGui::Checkbox("Draw Convex Exit", &DrawConvexExit);
    ImGui::Checkbox("Draw Dart Ortho", &DrawDartOrthoDomain);
    ImGui::Checkbox("Draw Jump Nodes", &DrawJumpNodes);
  }
  // ImGui::Separator();
  if (ImGui::CollapsingHeader("DRAW SAMPLING")) {
    // ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW SAMPLING");
    ImGui::Checkbox("Draw Must Pass", &Draw_must_pass_on);
    ImGui::Checkbox("Draw Sample On", &Draw_sample_on);
    ImGui::Checkbox("Draw Force Sample", &Draw_force_sample);
    ImGui::Checkbox("Draw Allow Path", &Draw_allow_path_on);
    ImGui::Checkbox("Draw Allow Arrive", &Draw_allow_arrive_to);

    ImGui::Checkbox("Draw Initial Pos", &DrawInitalPos);
    ImGui::InputInt("Cand Idx", &curr_cand);
    curr_cand = std::min(curr_cand, PathSampl.NumCandidates() - 1);
    ImGui::Checkbox("Draw Candidates", &DrawCandidates);
  }
  // ImGui::Separator();
  // ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW PATCH");
  if (ImGui::CollapsingHeader("DRAW PATCH")) {
    ImGui::Checkbox("Draw Corners", &DrawCorners);
    ImGui::Checkbox("Draw DualCorners", &DrawDualCorners);
    ImGui::Checkbox("Draw TJunctions", &DrawTJunction);

    ImGui::Checkbox("Draw Sides", &DrawSides);
    // ImGui::TextColored(ImVec4(1, 0, 1, 1), "COLOR PATCHES");
    static const char *itemsShape[] = {"Uniform",      "By Patch",   "if Disk",
                                       "Valence",      "Match Sing", "Unsolved",
                                       "Uv Distortion"};
    ImGui::SetNextItemWidth(100);
    ImGui::Combo("Color Patches", &CurrCol, itemsShape,
                 IM_ARRAYSIZE(itemsShape));
    if (CurrCol != OldCol) {
      UpdateColor();
      OldCol = CurrCol;
    }
  }
  if (ImGui::CollapsingHeader("DRAW SUBDIVISION")) {
    // ImGui::Checkbox("Draw Assing TJunct", &DrawAssignTJuctions);
    // ImGui::Checkbox("Draw Side Samples Sample", &DrawSubdivision);
    ImGui::Checkbox("Draw Side Samples Sample", &DrawSubdivision);
    ImGui::Checkbox("Draw Assign Values", &DrawAssignValues);
  }
  ImGui::End();
}

void UpdateAfterRemesh() {
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert, maxQCrossVert);

  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace, maxQCrossFace);
  // GetFaceSingularities(FaceCurv, VertPos, Connectivity,
  //                      NextF, NextE, SingIndex, SingValue);

  PatchM.SetSingularities(SingIndex, SingValue);
  PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
  UpdateColor();
}

void UpdateFeatureSing() {
  FeatureSingIndex.clear();
  FeatureSingValue.clear();
  std::vector<int> NewSingValue;
  std::vector<int> NewSingIndex;

  std::set<int> FeatureVSet;
  for (size_t i = 0; i < Features.size(); i++) {
    FeatureVSet.insert(Features[i].first);
    FeatureVSet.insert(Features[i].second);
  }

  for (size_t i = 0; i < SingIndex.size(); i++) {
    if (FeatureVSet.count(SingIndex[i]) == 0) {
      NewSingIndex.push_back(SingIndex[i]);
      NewSingValue.push_back(SingValue[i]);
    } else {
      FeatureSingIndex.push_back(SingIndex[i]);
      FeatureSingValue.push_back(SingValue[i]);
    }
  }

  SingIndex = NewSingIndex;
  SingValue = NewSingValue;
}

bool remove_batch = true;

void PreprocessMesh() {
  Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(
      VertPos, Connectivity, FaceCurv, Features, SingIndex, min_sing_distance);

  MeshBox.Init(VertPos);
  UpdateAfterRemesh();

  // WriteOBJ("after_remesh.obj",VertPos, Connectivity);
  // std::vector<Field::CrossF<ScalarType>> VertCurv0;
  // WriteField("after_remesh.ffield",VertCurv0,FaceCurv);
  // WriteFeatures("after_remesh.features",Features);
  // exit(0);

  std::cout << "*** SPLIT FEATURES *** " << std::endl;

  RemainingF.clear();

  std::cout << "* split corners " << std::endl;
  Geo::TracingPreprocess<ScalarType>::SplitCornerFeatures(
      VertPos, Connectivity, VertCurv, FaceCurv, Features, RemainingF,
      SingIndex, SingValue, FeatureSingIndex, FeatureSingValue, NewBorder);

  // WriteOBJ("Test0.obj",VertPos, Connectivity);
  // exit(0);

  std::cout << "* get min max " << std::endl;
  Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert, maxQCrossVert);
  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace, maxQCrossFace);

  std::cout << "* desambiguate " << std::endl;
  Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
      VertPos, Connectivity, &ModSameBVert, &ModSameBOrigPos);

  std::cout << "* update after remesh " << std::endl;
  UpdateAfterRemesh();
}

void BatchProcessMinPatch() {

  std::cout << "*** REFINE AND SPLIT *** " << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  PreprocessMesh();
  auto end0 = std::chrono::high_resolution_clock::now();
  auto elapsed0 =
      std::chrono::duration_cast<std::chrono::microseconds>(end0 - start);

  // WriteOBJ("testBeforeGraph.obj",VertPos, Connectivity);

  std::cout << "*** INIT GRAPHS *** " << std::endl;
  InitGraph();
  auto end1 = std::chrono::high_resolution_clock::now();
  auto elapsed1 =
      std::chrono::duration_cast<std::chrono::microseconds>(end1 - end0);

  std::cout << "*** INIT PATCHES *** " << std::endl;
  InitPatches();
  auto end2 = std::chrono::high_resolution_clock::now();
  auto elapsed2 =
      std::chrono::duration_cast<std::chrono::microseconds>(end2 - end1);

  std::cout << "*** ADD ALL PATHS *** " << std::endl;
  // PDeco.dynamic_smooth=true;

  remove_batch = true;
  PDeco.add_only_where_needed = true;
  PDeco.stop_when_solved = true;
  // PDeco.dynamic_smooth = false;
  //  PDeco.SolveAllFeatures();
  PDeco.BatchDecompose(remove_batch);

  SmoothPaths();

  auto end3 = std::chrono::high_resolution_clock::now();
  auto elapsed3 =
      std::chrono::duration_cast<std::chrono::microseconds>(end3 - end2);

  RestoreOriginalPosOnSplitBorders();

  std::cout << "*** SMOOTH PATCHES *** " << std::endl;
  SmoothPaths();

  auto end4 = std::chrono::high_resolution_clock::now();
  auto elapsed4 =
      std::chrono::duration_cast<std::chrono::microseconds>(end4 - end3);

  double elapsed_seconds0 =
      elapsed0.count() * 1e-6; // Convert microseconds to seconds
  double elapsed_seconds1 =
      elapsed1.count() * 1e-6; // Convert microseconds to seconds
  double elapsed_seconds2 =
      elapsed2.count() * 1e-6; // Convert microseconds to seconds
  double elapsed_seconds3 =
      elapsed3.count() * 1e-6; // Convert microseconds to seconds
  double elapsed_seconds4 =
      elapsed4.count() * 1e-6; // Convert microseconds to seconds

  showCross = false;

  auto total_elapsed =
      std::chrono::duration_cast<std::chrono::microseconds>(end4 - start);

  double total_elapsed_seconds5 =
      total_elapsed.count() * 1e-6; // Convert microseconds to seconds

  // std::cout << std::endl << std::endl;
  // std::cout << "*** FINAL TIMINGS (sec) *** " << std::endl;
  // std::cout << "* Remeshing time: " << elapsed_seconds0 << std::endl;
  // std::cout << "* Init Graph Elapsed time: " << elapsed_seconds1 <<
  // std::endl; std::cout << "* Init Patches: " << elapsed_seconds2 <<
  // std::endl; std::cout << "* Adding Paths: " << elapsed_seconds3 <<
  // std::endl; std::cout << "* Smooth Patches: " << elapsed_seconds4 <<
  // std::endl; std::cout << "* TOTAL TIME: " << total_elapsed_seconds5 <<
  // std::endl;

  std::cout << std::endl << std::endl;
  std::cout << "*** STATISTICS *** " << std::endl;
  std::cout << "*Primal* " << std::endl;
  PatchM.WriteStats();
  PDeco.WriteStats();
}

void BatchProcessMaxPaths() {

  std::cout << "*** REFINE AND SPLIT *** " << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  PreprocessMesh();
  auto end0 = std::chrono::high_resolution_clock::now();
  auto elapsed0 =
      std::chrono::duration_cast<std::chrono::microseconds>(end0 - start);

  std::cout << "*** INIT GRAPHS *** " << std::endl;
  InitGraph();
  auto end1 = std::chrono::high_resolution_clock::now();
  auto elapsed1 =
      std::chrono::duration_cast<std::chrono::microseconds>(end1 - end0);

  std::cout << "*** INIT PATCHES *** " << std::endl;
  InitPatches();
  auto end2 = std::chrono::high_resolution_clock::now();
  auto elapsed2 =
      std::chrono::duration_cast<std::chrono::microseconds>(end2 - end1);

  PDeco.BatchMaxDecompose(MaxPathRank);

  auto end3 = std::chrono::high_resolution_clock::now();
  auto elapsed3 =
      std::chrono::duration_cast<std::chrono::microseconds>(end3 - end2);

  RestoreOriginalPosOnSplitBorders();

  std::cout << "*** SMOOTH PATCHES *** " << std::endl;
  SmoothPaths();
  auto end4 = std::chrono::high_resolution_clock::now();
  auto elapsed4 =
      std::chrono::duration_cast<std::chrono::microseconds>(end4 - end3);

  double elapsed_seconds0 =
      elapsed0.count() * 1e-6; // Convert microseconds to seconds
  double elapsed_seconds1 =
      elapsed1.count() * 1e-6; // Convert microseconds to seconds
  double elapsed_seconds2 =
      elapsed2.count() * 1e-6; // Convert microseconds to seconds
  double elapsed_seconds3 =
      elapsed3.count() * 1e-6; // Convert microseconds to seconds
  double elapsed_seconds4 =
      elapsed4.count() * 1e-6; // Convert microseconds to seconds

  showCross = false;

  auto total_elapsed =
      std::chrono::duration_cast<std::chrono::microseconds>(end4 - start);

  double total_elapsed_seconds5 =
      total_elapsed.count() * 1e-6; // Convert microseconds to seconds

  std::cout << "*** TIMINGS (sec) *** " << std::endl;
  std::cout << "* Remeshing time: " << elapsed_seconds0 << std::endl;
  std::cout << "* Init Graph Elapsed time: " << elapsed_seconds1 << std::endl;
  std::cout << "* Init Patches: " << elapsed_seconds2 << std::endl;
  std::cout << "* Adding Paths: " << elapsed_seconds3 << std::endl;
  std::cout << "* Smooth Patches: " << elapsed_seconds4 << std::endl;
  std::cout << "* TOTAL TIME: " << total_elapsed_seconds5 << std::endl;

  // std::cout << std::endl << std::endl;
  // std::cout << "*** STATISTICS *** " << std::endl;
  // std::cout << "*Primal* " << std::endl;
  // PatchM.WriteStats();
  // PDeco.WriteStats();
}

int PrioIdx = 0;

void SetDecomposeBar() {

  ImGui::Begin(
      "Decompose", NULL,
      ImGuiWindowFlags_AlwaysAutoResize); // Create a window called "Hello,
                                          // world!" and append into it.

  ImGui::TextColored(ImVec4(1, 0, 1, 1), "BATCH");

  ImGui::SetNextItemWidth(100);
  ImGui::InputInt("Density", &PathSampl.NumSamples);
  ImGui::Checkbox("Dynamic Update", &PDeco.dynamic_updates);
  // ImGui::SetNextItemWidth(120);
  ImGui::Checkbox("Force Features", &PDeco.force_all_features);
  ImGui::Checkbox("Split Remove", &PDeco.split_removal);
  ImGui::Checkbox("Dynamic Smooth", &PDeco.dynamic_smooth);
  // split_removal
  // ImGui::SetNextItemWidth(120);
  // ImGui::Checkbox("Dynamic Smooth", &PDeco.dynamic_smooth);

  if (ImGui::Button("MINIMAL PATCH")) {
    BatchProcessMinPatch();
    PathColors.resize(PathSampl.SelectedPaths.size());
    Geo::ScatterColor<ScalarType>(PathColors);
    UpdateColor();
  }

  if (ImGui::Button("MAX PATHS")) {
    BatchProcessMaxPaths();
    PathColors.resize(PathSampl.SelectedPaths.size());
    Geo::ScatterColor<ScalarType>(PathColors);
    UpdateColor();
  }

  // if (ImGui::Button("Solve Remaining Spike/Concave")) {
  //   PDeco.ForceSolveSpikeConcave();
  //   PathColors.resize(PathSampl.SelectedPaths.size());
  //   Geo::ScatterColor<ScalarType>(PathColors);
  //   UpdateColor();
  //   PatchM.WriteStats();
  //   PDeco.WriteStats();
  // }

  ImGui::Separator();

  if (ImGui::CollapsingHeader("PREPROCESS")) {
    // ImGui::TextColored(ImVec4(1, 0, 1, 1), "PREPROCESS");

    if (ImGui::Button("Refine Sing")) {
      Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(
          VertPos, Connectivity, FaceCurv, Features, SingIndex,
          min_sing_distance);
      MeshBox.Init(VertPos);
      UpdateAfterRemesh();
      // WriteOBJ("refined_test.obj",VertPos,Connectivity);
    }

    if (ImGui::Button("Split Corners")) {
      RemainingF.clear();
      Geo::TracingPreprocess<ScalarType>::SplitCornerFeatures(
          VertPos, Connectivity, VertCurv, FaceCurv, Features, RemainingF,
          SingIndex, SingValue, FeatureSingIndex, FeatureSingValue, NewBorder);

      Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
          VertPos, Connectivity, &ModSameBVert, &ModSameBOrigPos);

      UpdateAfterRemesh();
    }

    ImGui::SameLine();
    ImGui::SetNextItemWidth(80);
    ImGui::InputInt("Sing Dist", &min_sing_distance);
    // if (ImGui::Button("Test"))
    // {
    //     TestSampling();
    // }
  }

  if (ImGui::CollapsingHeader("INIT GRAPH")) {

    // ImGui::TextColored(ImVec4(1, 0, 1, 1), "INIT");
    ImGui::SetNextItemWidth(150);
    ImGui::SliderFloat("Drift", &drift_penalty, 0, 1000);
    ImGui::SetNextItemWidth(150);
    ImGui::SliderFloat("Max Angle", &maxAngle, 0, 90);
    ImGui::SetNextItemWidth(150);
    ImGui::InputInt("Diffuse", &diffuse_step);
    ImGui::SetNextItemWidth(150);
    if (ImGui::Button("Init Graph")) {

      std::cout << "INIT GRAPH" << std::endl;
      InitGraph();
      std::cout << "INIT PATCHES" << std::endl;
      InitPatches();

      has_graph = true;
      showCross = false;
    }
  }
  ImGui::Separator();
  if (ImGui::CollapsingHeader("CONDITIONS")) {

    // ImGui::TextColored(ImVec4(1, 0, 1, 1), "CONDITIONS");
    ImGui::TextColored(ImVec4(0, 0, 1, 1), "structure");
    // ImGui::Checkbox("Darts", &PDeco.allow_darts);
    // ImGui::SameLine();
    // ImGui::Checkbox("Concave", &PDeco.allow_concave);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "topology");
    ImGui::Checkbox("1 Boundary", &single_boundary_cond);
    ImGui::SameLine();
    ImGui::Checkbox("Disk Like", &disk_like_cond);
    ImGui::SameLine();
    ImGui::Checkbox("Not Self Connection ", &self_connection_cond);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "singularity");
    ImGui::Checkbox("1 Sing Patch", &single_sing_cond);
    ImGui::SameLine();
    ImGui::Checkbox("Sing Match", &match_sing_cond);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "sides");
    ImGui::SetNextItemWidth(50);
    ImGui::InputInt("Min Sides", &MinSides, 0);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(50);
    ImGui::InputInt("Max Sides", &MaxSides, 0);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "shape");
    ImGui::SetNextItemWidth(100);

    ImGui::Checkbox("CCAbility", &cc_ability_cond);

    float CCAbilityf = CCability;
    ImGui::InputFloat("CC Tolerance", &CCAbilityf);
    CCAbilityf = std::min(CCAbilityf, 1.f);
    CCAbilityf = std::max(CCAbilityf, 0.f);
    CCability = (ScalarType)CCAbilityf;

    ImGui::Checkbox("Equal Side", &equal_side_cond);

    float EqSideThrf = EqSideThr;
    ImGui::InputFloat("EQSide Tolerance", &EqSideThrf);
    EqSideThrf = std::min(EqSideThrf, 1.f);
    EqSideThrf = std::max(EqSideThrf, 0.f);
    EqSideThr = (ScalarType)EqSideThrf;

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "UV");

    // Use ImGui::RadioButton for each option
    ImGui::Checkbox("UV", &use_uv);
    if (!use_uv) {
      isometric_uv = false;
      conformal_uv = false;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Isometric", isometric_uv)) {
      isometric_uv = true;
      conformal_uv = false;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Conformal", conformal_uv)) {
      isometric_uv = false;
      conformal_uv = true;
    }
    float max_uvDistf = max_uvDist;
    ImGui::SetNextItemWidth(100);
    ImGui::InputFloat("Max UV Distortion", &max_uvDistf);
    max_uvDist = max_uvDistf;

    bool needUpdate = UpdateConditions();
    if (needUpdate) {
      InitConditions();
      PDeco.UpdateAllUnsolvedPatches();
    }
  }
  if (!has_graph) {
    ImGui::End();
    return;
  }

  ImGui::Separator();
  ImGui::TextColored(ImVec4(1, 0, 1, 1), "PATH SAMPLING");

  ImGui::SetNextItemWidth(80);
  ImGui::InputInt("Samples", &PathSampl.NumSamples);
  ImGui::SameLine();
  ImGui::SetNextItemWidth(80);
  ImGui::InputInt("Sing Dist", &PathSampl.minSingResample);

  if (ImGui::Button("FEATURES")) {
    PathSampl.dynamicUpdates = PDeco.dynamic_updates;
    // PDeco.ForceAllFeatures();
    PDeco.SolveAllFeatures();
    //  PDeco.SolveSpikeSpike();
    PathColors.resize(PathSampl.SelectedPaths.size());
    Geo::ScatterColor<ScalarType>(PathColors);
    UpdateColor();
    // ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    // SmoothPaths();
  }
  ImGui::SameLine();
  if (ImGui::Button("ALL")) {
    PDeco.SolveAllMesh();
    PathColors.resize(PathSampl.SelectedPaths.size());
    Geo::ScatterColor<ScalarType>(PathColors);
    UpdateColor();
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    // SmoothPaths();
    //  SampleTracing();
  }
  ImGui::SameLine();
  if (ImGui::Button("REMOVE")) {
    PDeco.RemovePaths();
    PathColors.resize(PathSampl.SelectedPaths.size());
    Geo::ScatterColor<ScalarType>(PathColors);
    UpdateColor();
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  }
  ImGui::SameLine();
  if (ImGui::Button("SPLIT REMOVE")) {
    PDeco.RemoveSplitted();
    PathColors.resize(PathSampl.SelectedPaths.size());
    Geo::ScatterColor<ScalarType>(PathColors);
    UpdateColor();
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  }
  if (ImGui::CollapsingHeader("Sample Step")) {
    ImGui::BeginTable("Involved Nodes", 5);
    ImGui::TableSetupColumn("");
    ImGui::TableSetupColumn("Dart");
    ImGui::TableSetupColumn("Concave");
    ImGui::TableSetupColumn("Border");
    ImGui::TableSetupColumn("Internal");
    ImGui::TableHeadersRow();

    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::Text("Start");

    ImGui::TableSetColumnIndex(1);
    ImGui::Checkbox("##foo1", &PDeco.start_spike);
    ImGui::TableSetColumnIndex(2);
    ImGui::Checkbox("##foo2", &PDeco.start_concave);
    ImGui::TableSetColumnIndex(3);
    ImGui::Checkbox("##foo3", &PDeco.start_orginal_border);
    ImGui::TableSetColumnIndex(4);
    ImGui::Checkbox("##foo4", &PDeco.start_internal);

    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::Text("End");

    ImGui::TableSetColumnIndex(1);
    ImGui::Checkbox("##foo5", &PDeco.end_spike);
    ImGui::TableSetColumnIndex(2);
    ImGui::Checkbox("##foo6", &PDeco.end_concave);
    ImGui::TableSetColumnIndex(3);
    ImGui::Checkbox("##foo7", &PDeco.end_patch_border);
    ImGui::TableSetColumnIndex(4);
    // ImGui::Checkbox("##foo8", &PDeco.end_internal);

    ImGui::EndTable();

    // ImGui::Begin("Priority Type");
    // // Radio buttons with labels
    // if (ImGui::RadioButton("Distribution", PrioIdx == 0))
    // {
    //   PrioIdx=0;
    //   PDeco.SetDistPrioType();
    // }
    // if (ImGui::RadioButton("Lenght", PrioIdx == 1))
    // {
    //   PDeco.SetLengthPrioType();
    //   PrioIdx=1;
    // }
    // ImGui::End();

    // ImGui::SameLine();
    if (ImGui::Button("Sample Candidates")) {
      SampleTracing();
    }
  }

  ImGui::Separator();
  ImGui::TextColored(ImVec4(1, 0, 1, 1), "PATH DISTRIBUTION");

 
  if (ImGui::CollapsingHeader("Distributing Step")) {
    if (ImGui::Button("Add Next")) {
      bool terminated = false;
      bool solved = false;

      // PathSampl.InitializeDistance();
      // PathSampl.InitCandidates();
      // PathSampl.UpdateCandidatePriority();

      PDeco.AddNext(solved, terminated);

      PathColors.resize(PathSampl.SelectedPaths.size());
      Geo::ScatterColor<ScalarType>(PathColors);
      UpdateColor();
    }

    ImGui::SameLine();
    if (ImGui::Button("Resample")) {
      // SampleTracing();
    }
    ImGui::SameLine();
    if (ImGui::Button("Remove")) {
      bool has_removed = PDeco.RemoveNextIfPossible();
      if (has_removed)
        PathSampl.MergeCollinearPath();

      PathColors.resize(PathSampl.SelectedPaths.size());
      Geo::ScatterColor<ScalarType>(PathColors);
      UpdateColor();
    }
    ImGui::SameLine();
    if (ImGui::Button("Split")) {
      PathSampl.SplitPaths();
      PathColors.resize(PathSampl.SelectedPaths.size());
      Geo::ScatterColor<ScalarType>(PathColors);
      UpdateColor();
    }

    // ImGui::SetNextItemWidth(120);
    // ImGui::Checkbox("Dynamic Upd", &PathSampl.dynamicUpdates);
    ImGui::SetNextItemWidth(120);
    ImGui::Checkbox("Only Needed", &PDeco.add_only_where_needed);
    ImGui::SetNextItemWidth(120);
    ImGui::Checkbox("Stop Solved", &PDeco.stop_when_solved);
    ImGui::SameLine();
    // ImGui::SetNextItemWidth(120);
    // ImGui::Checkbox("Dynamic Smth", &PDeco.dynamic_smooth);
    ImGui::SetNextItemWidth(120);
    ImGui::Checkbox("Only Visible", &PathSampl.only_visible);
  }

  ImGui::End();
}

void SaveSubMeshes() {
  std::string subMeshFilename = ProjectName;
  subMeshFilename.append("_submesh.obj");
  PDeco.SaveSubMeshes(subMeshFilename);
}

void SavePatches() {
  std::string PatchFilename = ProjectName;
  PatchFilename.append("_layout.patch");
  PDeco.SavePatchesLayout(PatchFilename);
}

void SaveField() {
  std::string FinalFieldFilename = ProjectName;
  FinalFieldFilename.append("_layout.ffield");
  WriteField(FinalFieldFilename.c_str(), VertCurv, FaceCurv);
}

void SaveMesh() {
  std::string FinalMFilename = ProjectName;
  FinalMFilename.append("_layout.obj");
  WriteOBJ(FinalMFilename.c_str(), VertPos, Connectivity);
}

void SaveAll() {
  SaveSubMeshes();
  SavePatches();
  SaveField();
  SaveMesh();
}

void SaveDualPatches() {
  std::vector<std::vector<Geo::Point3<ScalarType>>> SubMeshesV;
  std::vector<std::vector<std::vector<int>>> SubMeshesF;
  std::vector<std::vector<int>> OriginalV, OriginalF;
  SplitByFFConnectedComponenent(DualVertPos, DualConnectivity, SubMeshesV,
                                SubMeshesF, OriginalV, OriginalF);
  for (size_t i = 0; i < SubMeshesV.size(); i++) {
    std::string FinalMFilename = ProjectName;
    std::string IDX = std::to_string(i);
    FinalMFilename.append(IDX);
    FinalMFilename.append("_dual.obj");
    std::cout << "Saving to:" << FinalMFilename.c_str() << std::endl;
    WriteOBJ(FinalMFilename.c_str(), SubMeshesV[i], SubMeshesF[i]);
  }
  // ConnectedCO
}

Geo::SmoothParam<ScalarType> SParamQuadSmooth;
Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;

void UpdateQuadFeatures() {
  std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
      FeatureData;
  SideD.GetFeaturesQuadrData(DualVertPos, DualFeatures, FeatureData);

  QuadFeatures.clear();
  bool has_mapped = FindEdgePosCorrespondences(QuadVert, QuadFaces, FeatureData,
                                               QuadFeatures);
  assert(has_mapped);

  // find the corners
  QuadCorners.clear();
  std::set<Geo::Point3<ScalarType>> CornerPosSet;
  for (size_t i = 0; i < CornersIdx.size(); i++)
    CornerPosSet.insert(VertPos[CornersIdx[i]]);

  Geo::MeshFeatures<ScalarType>::CornersByCount(QuadVert, QuadFeatures, QuadCorners,true);

  for (size_t i = 0; i < QuadVert.size(); i++) {
    if (CornerPosSet.count(QuadVert[i]) > 0)
      QuadCorners.push_back(i);
  }
}

void Init3DSmooth() {
  SParamQuadSmooth.Clear();

  SParamQuadSmooth.Features = QuadFeatures;
  SParamQuadSmooth.FixedVert = QuadCorners;

  std::vector<std::pair<int, int>> BorderE;
  Geo::getBorderEdges(QuadFaces, BorderE);
  for (size_t i = 0; i < BorderE.size(); i++) {
    SParamQuadSmooth.FixedVert.push_back(BorderE[i].first);
    SParamQuadSmooth.FixedVert.push_back(BorderE[i].second);
  }

  SParamQuadSmooth.ReprojVert = &VertPos;
  SParamQuadSmooth.ReprojFaces = &Connectivity;
  SParamQuadSmooth.Damp = 0.5;
  SParamQuadSmooth.post_project = true;
}

void SmoothQuadMesh3D() {

  SParamQuadSmooth.NumIte = QuadSmoothSteps;
  SParamQuadSmooth.SType = Geo::SMLaplacian;
  SParamQuadSmooth.check_AspectR = false;

  Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

  Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                  QuadVertNormals);
}

void SmoothQuadMeshSeparatrix() {

  Geo::QuadMeshSeparatrixSmooth<ScalarType>::SeparatrixSmoothParam Sep_SParam;
    Sep_SParam.Fixed = QuadCorners;
    Sep_SParam.FeatureLines = QuadFeatures;
    Sep_SParam.Damp = 0.5;
    Sep_SParam.Steps = QuadSmoothSteps;
    Sep_SParam.fixBorder = true;
    Sep_SParam.ReprojVert = &VertPos;
    Sep_SParam.ReprojFaces = &Connectivity;
    Sep_SParam.ReprojFeatureLines = Features;
    Geo::QuadMeshSeparatrixSmooth<ScalarType>::Smooth(QuadVert, QuadFaces,
                                                      Sep_SParam);
  Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                  QuadVertNormals);
}

void SmoothQuadMeshGradient() {

  Geo::SmoothSampling<ScalarType>::Improve(QuadVert, QuadFaces,QuadFeatures,QuadCorners,QuadSmoothSteps);
  Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                  QuadVertNormals);
}

void SolveDoublets() {
  Geo::QuadPatchFunctions<ScalarType>::RemoveDoublets(
      QuadVert, QuadFaces, QuadFeatures, QuadCorners, QuadPatch);
  Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                  QuadVertNormals);
  //ElemColorByCluster<ScalarType>(QuadPatch, QuadFaceColors);
}

void SmoothBestPolyMesh3D() {

  SParamQuadSmooth.NumIte = QuadSmoothSteps;
  SParamQuadSmooth.SType = Geo::SMBestPoly;
  SParamQuadSmooth.check_AspectR = false;
  Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

  Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                  QuadVertNormals);
}


Geo::DualPatch<ScalarType>::DualParameter DualParam;

int EvenMode = 0;
int OnlyDualP = -1;

void SetExportBar() {

  ImGui::Begin(
      "Finalize", NULL,
      ImGuiWindowFlags_AlwaysAutoResize); // Create a window called "Hello,
                                          // world!" and append into it.

  ImGui::Separator();
  // ImGui::TextColored(ImVec4(1, 0, 1, 1), "SMOOTH");
  if (ImGui::CollapsingHeader("PREPROCESS")) {

    ImGui::SetNextItemWidth(80);
    ImGui::InputInt("Smooth Steps", &smooth_steps);
    ImGui::SameLine();
    if (ImGui::Button("Smooth"))
      SmoothPaths();
  }

  ImGui::Separator();
  // ImGui::TextColored(ImVec4(1, 0, 1, 1), "DUALIZE");

  if (ImGui::CollapsingHeader("DUAL")) {

    ImGui::Checkbox("Only Non-Quad", &DualParam.OnlyNonQuad);
    ImGui::Checkbox("Init Sing Pos", &DualParam.InitSingPos);
    ImGui::Checkbox("Fix Sing Pos", &DualParam.FixSingPos);
    ImGui::Checkbox("Joint Border", &DualParam.JointBorder);

    if (DualParam.RType == Geo::DualPatch<ScalarType>::RT_UVDist_Tutte)
      EvenMode = 0;
    if (DualParam.RType == Geo::DualPatch<ScalarType>::RT_UVDist_LSC)
      EvenMode = 1;
    if (DualParam.RType == Geo::DualPatch<ScalarType>::RT_UVDist_ARAP)
      EvenMode = 2;
    if (DualParam.RType == Geo::DualPatch<ScalarType>::RT_CCability)
      EvenMode = 3;

    static const char *itemsEvenMode[] = {"UV_Tutte", "UV_LSC", "UV_Arap",
                                          "UV_CCability"};

    ImGui::Combo("Dual Mode", &EvenMode, itemsEvenMode,
                 IM_ARRAYSIZE(itemsEvenMode));

    if (EvenMode == 0)
      DualParam.RType = Geo::DualPatch<ScalarType>::RT_UVDist_Tutte;
    if (EvenMode == 1)
      DualParam.RType = Geo::DualPatch<ScalarType>::RT_UVDist_LSC;
    if (EvenMode == 2)
      DualParam.RType = Geo::DualPatch<ScalarType>::RT_UVDist_ARAP;
    if (EvenMode == 3)
      DualParam.RType = Geo::DualPatch<ScalarType>::RT_CCability;

    ImGui::InputInt("Only Patch", &OnlyDualP);
    DualParam.OnlyPatch.clear();
    if (OnlyDualP < -1)
      OnlyDualP = -1;
    OnlyDualP = std::min(OnlyDualP, PatchM.NumPatches() - 1);

    if (OnlyDualP >= 0)
      DualParam.OnlyPatch.push_back(OnlyDualP);

    if (ImGui::Button("Smooth Dual Patches")) {
      Geo::PatchOptimize<ScalarType>::DualMeshSmooth(
          VertPos, Connectivity, Features, DualVertPos, DualConnectivity);

      Geo::DualPatch<ScalarType>::ExtractCornerVisPos(
          DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
          PerDualPatchCornerIDX, VisDualCornerPos);
    }
  }

  if (ImGui::Button("Extract Dual")) {

    RestoreOriginalPosOnSplitBorders();
    Geo::DualPatch<ScalarType>::DisambiguateForDual(PatchM);

    // Geo::DualPatch<ScalarType>::ExtractDual(
    //     VertPos, Connectivity, PatchM, DualParam, DualVertPos,
    //     DualConnectivity, PerDualPatchFaceIDX, PerDualPatchCornerIDX);

    Geo::DualPatch<ScalarType>::PrimalPatchData PData;
    PData.PatchFaces = PatchM.PData.SubPatchFaces;
    PData.PatchVertPos = PatchM.PData.SubPatchVertPos;
    PData.PatchCorners = PatchM.PData.PatchCorners;
    PData.EdgePatchSides = PatchM.PData.EdgePatchSides;
    for (size_t i = 0; i < PatchM.NumPatches(); i++) {
      std::vector<int> PatchSing;
      PatchM.PatchSingularities(i, PatchSing);
      PData.PatchSingularities.push_back(PatchSing);
    }

    Geo::DualPatch<ScalarType>::ExtractDual(
        PData, DualParam, DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
        PerDualPatchCornerIDX);

    // std::cout << "Size F:" << PerDualPatchFaceIDX.size() << std::endl;
    // std::cout << "Vert FIdx:" << PerDualPatchVertIDX.size() << std::endl;
    // std::cout << "Corner FIdx:" << PerDualPatchCornerIDX.size() << std::endl;
    // // exit(0);

    // UPDATE COLOR
    FaceColDual.clear();
    FaceColDual.resize(DualConnectivity.size(),
                       Geo::Point3<ScalarType>(0.5, 1, 0.83));

    for (size_t i = 0; i < PerDualPatchFaceIDX.size(); i++) {
      Geo::Point3<ScalarType> CurrCol =
          Geo::ScatterColor<ScalarType>(PerDualPatchFaceIDX.size(), i);

      if (PerDualPatchCornerIDX[i].size() != 4)
        CurrCol = Geo::Point3<ScalarType>(1, 0, 0);

      for (size_t j = 0; j < PerDualPatchFaceIDX[i].size(); j++) {
        int IndexF = PerDualPatchFaceIDX[i][j];
        FaceColDual[IndexF] = CurrCol;
      }
    }
    // UPDATE CORNERS
    Geo::DualPatch<ScalarType>::ExtractCornerVisPos(
        DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
        PerDualPatchCornerIDX, VisDualCornerPos);

    std::vector<int> PrimalCorners;
    Geo::DualPatch<ScalarType>::RetrieveDualFeatures(
        VertPos, Connectivity, Features, PrimalCorners, DualVertPos,
        DualConnectivity, DualFeatures, DualCorners);

    has_dual_mesh = true;

    // UPDATE MESH DATA
    ComputeNormals(DualVertPos, DualConnectivity, DualFaceNormals,
                   DualVertNormals);

    Geo::getFFAdjacency(DualConnectivity, DualNextF, DualNextE);

    // CONTROLS FOR VISUALIZATION
    showMesh = false;
    showDual = true;
    DrawPathMode = 2;

    // INITIALIZE SIDE SUBDIVISION
    SideD.InitFromLabels(DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
                         PerDualPatchCornerIDX, DualFeatures);
    
    WriteOBJ("Dual1.obj",DualVertPos, DualConnectivity);
  }

  // PatchM.ExtractDual();
  if (has_dual_mesh) {
    ImGui::Separator();
    // ImGui::TextColored(ImVec4(1, 0, 1, 1), "QUADRANGULATE");

    ScalarType minValue = 0.001;
    ScalarType maxValue = 0.2;

    //     ScalarType MinSize=0.005;
    // ScalarType MaxSize=0.1;
    float MinSizef = MinSize;
    float MaxSizef = MaxSize;
    ImGui::InputFloat("Min Interval", &MinSizef);
    ImGui::InputFloat("Max Interval", &MaxSizef);
    MinSize = MinSizef;
    MaxSize = MaxSizef;
    MaxSize = std::min(MaxSize, maxValue);
    MinSize = std::max(MinSize, minValue);
    MaxSize = std::max(MaxSize, MinSize);

    if (ImGui::Button("Auto Detect Size")) {
      AutoQuadEdgeSize();
    }

    ImGui::SliderScalar("Quad Size", ImGuiDataType_Double, &QuadSize, &minValue,
                        &maxValue, "%.3f");

    if (ImGui::CollapsingHeader("SUBDIVISION")) {

      ScalarType minValueRvsIso = 0;
      ScalarType maxValueRvsIso = 1;
      ImGui::SliderScalar("Reg Vs Iso", ImGuiDataType_Double, &RegVsIso,
                          &minValueRvsIso, &maxValueRvsIso, "%.2f");

      ScalarType minRatio = 2;
      ScalarType maxRatio = 8;
      ImGui::SetNextItemWidth(80);
      ImGui::SliderScalar("Min Rt", ImGuiDataType_Double, &SideD.MinRatio,
                          &minRatio, &maxRatio, "%.1f");
      ImGui::SameLine();
      ImGui::SetNextItemWidth(80);
      ImGui::SliderScalar("Max Rt", ImGuiDataType_Double, &SideD.MaxRatio,
                          &minRatio, &maxRatio, "%.1f");

      ImGui::Checkbox("Weight FaceRatio", &SideD.WeightRegFaceRatio);
      ImGui::Checkbox("Weight Features", &SideD.WeightFeatures);

      bool useClosestInt =
          (SideD.FixType == Geo::IntegerLinearConstraintLeastSQSolver<
                                ScalarType>::ST_ClosestInt);
      bool useLowError =
          (SideD.FixType ==
           Geo::IntegerLinearConstraintLeastSQSolver<ScalarType>::ST_LowError);
      bool useStQRPivoting =
          (SideD.FixType == Geo::IntegerLinearConstraintLeastSQSolver<
                                ScalarType>::ST_StQRPivoting);

      if (ImGui::RadioButton("Closest", useClosestInt)) {
        SideD.FixType = Geo::IntegerLinearConstraintLeastSQSolver<
            ScalarType>::ST_ClosestInt;
      }
      ImGui::SameLine();
      if (ImGui::RadioButton("Low Err", useLowError)) {
        SideD.FixType =
            Geo::IntegerLinearConstraintLeastSQSolver<ScalarType>::ST_LowError;
      }
      ImGui::SameLine();
      if (ImGui::RadioButton("SQR Piv", useStQRPivoting)) {
        SideD.FixType = Geo::IntegerLinearConstraintLeastSQSolver<
            ScalarType>::ST_StQRPivoting;
      }
    }

    if (SideD.PMode == Geo::SideDiscretization<ScalarType>::NoParity)
      EvenMode = 0;
    if (SideD.PMode == Geo::SideDiscretization<ScalarType>::EvenVar)
      EvenMode = 1;
    if (SideD.PMode == Geo::SideDiscretization<ScalarType>::SolvePost)
      EvenMode = 2;

    static const char *itemsParity[] = {"No Parity", "Even Variable",
                                        "Fix Parity"};

    ImGui::Combo("Parity Mode", &EvenMode, itemsParity,
                 IM_ARRAYSIZE(itemsParity));

    Geo::SideDiscretization<ScalarType>::parityMode PMode;

    if (EvenMode == 0)
      PMode = Geo::SideDiscretization<ScalarType>::NoParity;
    if (EvenMode == 1)
      PMode = Geo::SideDiscretization<ScalarType>::EvenVar;
    if (EvenMode == 2)
      PMode = Geo::SideDiscretization<ScalarType>::SolvePost;

    SideD.PMode = PMode;

    ImGui::Checkbox("Optimize Size", &OptimizeSize);

    if (ImGui::Button("Find Subdivision")) {
      DrawSubdivision = true;
      DrawQuadMesh = false;
      showMesh = true;
      showDual = true;
      ScalarType EdgeS = MeshBox.Diag() * QuadSize;

      SideD.PrioritizeTJunction = true;

      std::cout<<"Edge Size:"<<EdgeS<<std::endl;
      SideD.FindDiscretization(EdgeS, RegVsIso, OptimizeSize);
    }

    if (ImGui::CollapsingHeader("Solve Steps")) {

      if (ImGui::Button("Solve System")) {
        DrawSubdivision = true;
        DrawQuadMesh = false;
        showMesh = true;
        showDual = true;
        ScalarType EdgeS = MeshBox.Diag() * QuadSize;

        SideD.PrioritizeTJunction = true;

        SideD.SolveSystem(EdgeS, RegVsIso);
      }

      if (ImGui::Button("Set To Target")) {
        DrawSubdivision = true;
        DrawQuadMesh = false;
        showMesh = true;
        showDual = true;
        ScalarType EdgeS = MeshBox.Diag() * QuadSize;

        SideD.SetToTarget(EdgeS);
      }

      if (ImGui::Button("Solve Parity")) {
        DrawSubdivision = true;
        DrawQuadMesh = false;
        showMesh = true;
        showDual = true;
        SideD.SolveParityViaGraph();
      }
      if (ImGui::Button("Solve Equality")) {
        DrawSubdivision = true;
        DrawQuadMesh = false;
        showMesh = true;
        showDual = true;
        SideD.SolveEqualityViaGraph();
      }

      if (ImGui::Button("Solve Size with Cycle")) {
        DrawSubdivision = true;
        DrawQuadMesh = false;
        showMesh = true;
        showDual = true;
        SideD.SolveSizeViaGraphCycles();
      }
    }

    if (ImGui::Button("Quadrangulate")) {
      DrawSubdivision = false;
      DrawQuadMesh = true;
      showMesh = false;
      showDual = false;
      std::vector<std::vector<int>> PatchFaces, PatchCorners;
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          SideSamples;

      SideD.GetQuadrangulationData(PatchFaces, PatchCorners, SideSamples);

      Geo::Quadrangulation<ScalarType>::QuadrangulateQuadPatchsWithSubs(
          DualVertPos, DualConnectivity, PatchFaces, PatchCorners, SideSamples,
          Q_LUT, QuadVert, QuadFaces, QuadPatch);

      RemoveDuplicatedVert(QuadVert, QuadFaces);
      RemoveUnfererencedVert(QuadVert, QuadFaces);

      // SmoothQuadMesh();

      //ElemColorByCluster<ScalarType>(QuadPatch, QuadFaceColors);

      QuadFaceColors = std::vector<Geo::Point3<ScalarType>>(
          QuadFaces.size(), Geo::Point3<ScalarType>(0.8, 0.8, 0.8));

      Geo::ComputeNormals<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                      QuadVertNormals);

      UpdateQuadFeatures();
      SolveDoublets();
      Init3DSmooth();
      // InitTangentSmooth();
      //  WriteOBJ("testQuad.obj",QuadVert, QuadFaces);
    }
  }

  ImGui::InputInt("Quad Smooth Steps", &QuadSmoothSteps, 0, 200);
  if (ImGui::Button("Laplacian Smooth 3D")) {
    SmoothQuadMesh3D();
  }

  if (ImGui::Button("Best Poly Smooth 3D")) {
    SmoothBestPolyMesh3D();
  }

  if (ImGui::Button("Smooth Separatrices")) {
   SmoothQuadMeshSeparatrix();
  }

  if (ImGui::Button("Smooth Gradient")) {
   SmoothQuadMeshGradient();
  }

  ImGui::Separator();
  // ImGui::TextColored(ImVec4(1, 0, 1, 1), "SAVE");
  if (ImGui::CollapsingHeader("SAVE")) {

    if (ImGui::Button("Mesh"))
      SaveMesh();
    ImGui::SameLine();
    if (ImGui::Button("Field"))
      SaveField();

    if (ImGui::Button("SubMeshes"))
      SaveSubMeshes();

    ImGui::SameLine();
    if (ImGui::Button("Patches"))
      SavePatches();

    if (ImGui::Button("Dual Patches")) {
      SaveDualPatches();
    }
    ImGui::SameLine();
    if (ImGui::Button("Dual Mesh")) {
      std::string FinalMFilename = ProjectName;
      FinalMFilename.append("_dual.obj");
      WriteOBJ(FinalMFilename.c_str(), DualVertPos, DualConnectivity);
    }

    if (ImGui::Button("Save ALL"))
      SaveAll();
    if (ImGui::Button("Save Loops")) {
      std::string LoopsFilename = ProjectName;
      LoopsFilename.append("_loops.txt");
      WriteCrossGraphPath<ScalarType>(LoopsFilename,VertPos,PathSampl.SelectedPaths);
    }
    if (ImGui::Button("Save Cycles")) {
      std::string CyclesFilename = ProjectName;
      CyclesFilename.append("_cycles.txt");
      PDeco.SaveCycleData(CyclesFilename);
      PDeco.SaveCurveCycleData(ProjectName);
      //PDeco.SaveSeparateSubMeshes(std::string("test"), std::string("off"));
    }
  }

  ImGui::End();
}

int BorderSize = 15;


void GLDrawMesh() {
  assert((glGetError()) == GL_NO_ERROR);

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
  // if (showTrackBall)
  //     trackball.DrawPostApply();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());

  assert((glGetError()) == GL_NO_ERROR);

  if (showDual || showMesh) {
    // // glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w,
    // // mesh_color.z * mesh_color.w, mesh_color.w);
    // glPushAttrib(GL_ALL_ATTRIB_BITS);

    // glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w,
    // mesh_color.z * mesh_color.w, mesh_color.w);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    bool wire = true;
    bool smoothshade = true;
    switch (DrawMeshMode) {
    case 0:
      wire = false;
      smoothshade = true;
      break;
    case 1:
      wire = false;
      smoothshade = false;
      break;
    case 2:
      wire = true;
      smoothshade = true;
      break;
    default:
      wire = true;
      smoothshade = false;
      break;
    }
    if (showDual) {
      GLDraw::glColor(Geo::Point3<ScalarType>(0.8, 0.8, 0.8));
      GLDraw::GLDrawBorders<ScalarType>(DualVertPos, DualConnectivity,
                                        DualNextF, 5);
      glLineWidth(line_thick);
      GLDraw::DrawSurfaceMesh<ScalarType>(
          DualVertPos, DualConnectivity, DualFaceNormals, DualVertNormals, true,
          wire, smoothshade, false, GLDraw::TRTriangle, &FaceColDual);
    }
    if (showMesh) {

      GLDraw::GLDrawBorders<ScalarType>(VertPos, Connectivity, NextF,
                                        BorderSize);
      glLineWidth(line_thick);
      GLDraw::DrawSurfaceMesh<ScalarType>(VertPos, Connectivity, FaceNormals,
                                          VertNormals, true, wire, smoothshade,
                                          false, GLDraw::TRTriangle, &FaceCol);
      // glPopAttrib();
    }
    if (showFeatures) {
      GLDraw::DrawVertices(VertPos, CornersIdx,
                           Geo::Point3<ScalarType>(0, 0, 0), 30);

      if (showMesh)
        GLDraw::GLDrawEdges<ScalarType>(VertPos, Features, 15,
                                        Geo::Point3<ScalarType>(0, 1, 1));
      if (showDual)
        GLDraw::GLDrawEdges<ScalarType>(DualVertPos, DualFeatures, 15,
                                        Geo::Point3<ScalarType>(1, 0, 1));
    }
    glPopAttrib();
    assert((glGetError()) == GL_NO_ERROR);
  }

  assert((glGetError()) == GL_NO_ERROR);

  if (has_cross_field) {
    ScalarType scaleVal = AVLen * (ScalarType)0.5;

    assert((glGetError()) == GL_NO_ERROR);

    if (showCross) {
      if (DrawFieldMode == 0) {
        GLDraw::DrawCrossFields<ScalarType>(FaceCurv, scaleVal, maxQCrossFace,
                                            minQCrossFace);
      } else {
        GLDraw::DrawCrossFields<ScalarType>(VertCurv, scaleVal, maxQCrossVert,
                                            minQCrossVert);
      }
    }
    assert((glGetError()) == GL_NO_ERROR);

    if (showSing)
      GLDraw::DrawSingularity<ScalarType>(VertPos, SingIndex, SingValue);

    assert((glGetError()) == GL_NO_ERROR);
  }
  assert((glGetError()) == GL_NO_ERROR);
  GLGraphDraw::GLCrossFieldGraph<ScalarType>::PathDrawParam PDParam;
  PDParam.DrawPath = true;
  if (DrawPathMode == 0) {
    if (DrawTJunction)
      GLDraw::GLDrawTJunctions<ScalarType>(VertPos, TJunctions);

    assert((glGetError()) == GL_NO_ERROR);

    PDParam.lineSize = 10;
    PDParam.pointSize = 0;
    // std::vector<std::vector<int>> Loops;

    if (CurrSelP == -1)
      GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawPaths(
          TGraph, VertCurv, PathSampl.SelectedPaths, MeshBox.Diag() * 0.005,
          PDParam, PathColors, DrawLoops, DrawPolyL);
    else
      GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodes(
          TGraph, VertCurv, PathSampl.SelectedPaths[CurrSelP],
          MeshBox.Diag() * 0.005, PDParam, PathColors[CurrSelP]);
  }
  assert((glGetError()) == GL_NO_ERROR);
  if (DrawPathMode == 1) {
    std::vector<ScalarType> RankThick;
    RankThick.push_back(ThickPrimaryP);
    RankThick.push_back(ThickSecondaryP);
    RankThick.push_back(ThickTertiaryP);
    std::vector<Geo::Point3<ScalarType>> RankCol;
    RankCol.push_back(Geo::Point3<ScalarType>(0, 0, 0));
    RankCol.push_back(Geo::Point3<ScalarType>(0.3, 0.3, 0.3));
    RankCol.push_back(Geo::Point3<ScalarType>(0.3, 0.3, 0.3));
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawPathRanks(
        TGraph, VertCurv, PathSampl.SelectedPaths, MaxPathRank, RankThick,
        RankCol);

    // std::cout<<"Num Solved: "<<PDeco.NumSolved<<"total
    // "<<PathSampl.SelectedPaths.size()<<std::endl;
    //  GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawPathsTwoRanks( TGraph,
    //  VertCurv, PathSampl.SelectedPaths,
    //                                                                PDeco.NumSolved,ThickSecondaryP,ThickPrimaryP,
    //                                                                Geo::Point3<ScalarType>(0.3,
    //                                                                0.3, 0.3),
    //                                                                Geo::Point3<ScalarType>(0,
    //                                                                0, 0));

    // GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawPathsRank(
    //     TGraph, VertCurv, PathSampl.SelectedPaths, 2, BorderSize);

    GLDraw::GLDrawBorders<ScalarType>(VertPos, Connectivity, NextF,
                                      ThickBorders);

    GLDraw::DrawOutline<ScalarType>(VertPos, Connectivity, ThickOutline);
  }

  PDParam.DrawPath = false;

  // GLGraphDraw::GLDrawPatchData<ScalarType>::DrawSpikesConcave(PatchM,40);
  assert((glGetError()) == GL_NO_ERROR);
  if (DrawCorners)
    GLGraphDraw::GLDrawPatchData<ScalarType>::DrawCorners(PatchM);
  if (DrawDualCorners)
    GLDraw::DrawPoints(VisDualCornerPos);

  if (DrawSides)
    GLGraphDraw::GLDrawPatchData<ScalarType>::DrawSides(PatchM);
  assert((glGetError()) == GL_NO_ERROR);
  if (DrawInitalPos) {
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::PathDrawParam PDParam1;
    PDParam1.ColPoint = Geo::Point3<ScalarType>(1, 0, 1);
    PDParam1.pointSize = 20;
    PDParam1.lineSize = 10;
    PDParam1.DrawPath = false;
    std::vector<int> NodesIdx;
    PathSampl.GetStartingNodes(NodesIdx);
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodes(
        TGraph, VertCurv, NodesIdx, MeshBox.Diag() * 0.005, PDParam1);

  }

  // if (Draw_must_pass_on)
  //   GLDraw::DrawVertices<ScalarType>(VertPos, PathSampl.must_pass_on);
  if (Draw_sample_on)
    GLDraw::DrawVertices<ScalarType>(VertPos, PathSampl.sample_on);
  if (Draw_force_sample)
    GLDraw::DrawVertices<ScalarType>(VertPos, PathSampl.force_sample);
  // if (Draw_allow_path_on)
  //   GLDraw::DrawVertices<ScalarType>(VertPos, PathSampl.allow_path_on);
  if (Draw_allow_arrive_to) {
    GLDraw::DrawVertices<ScalarType>(VertPos, PathSampl.allow_arrive_to);
  }

  assert((glGetError()) == GL_NO_ERROR);

  if (DrawCandidates) {
    GLGraphDraw::GLDrawPatchData<ScalarType>::DrawCandidates(
        PathSampl, MeshBox.Diag() * 0.005, curr_cand);
  }
  assert((glGetError()) == GL_NO_ERROR);

  PDParam.lineSize = 5;
  if (DrawExitNode) {
    PDParam.ColLine = Geo::Point3<ScalarType>(0, 1, 0);
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodesOfType(
        TGraph, VertCurv, Geo::ExitDomain, MeshBox.Diag() * 0.01, PDParam);
  }

  if (DrawConvexExit) {
    PDParam.ColLine = Geo::Point3<ScalarType>(1, 0, 0);
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodesOfType(
        TGraph, VertCurv, Geo::ConvexExitDomain, MeshBox.Diag() * 0.01,
        PDParam);
  }

  if (DrawTangentNode) {
    PDParam.ColLine = Geo::Point3<ScalarType>(1, 0, 0);
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodesOfType(
        TGraph, VertCurv, Geo::TangentDomain, MeshBox.Diag() * 0.01, PDParam);
  }

  if (DrawTangentExit) {
    PDParam.ColLine = Geo::Point3<ScalarType>(0, 1, 1);
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodesOfType(
        TGraph, VertCurv, Geo::TangentExitDomain, MeshBox.Diag() * 0.01,
        PDParam);
  }

  if (DrawDartOrthoDomain) {
    PDParam.ColLine = Geo::Point3<ScalarType>(0, 1, 1);
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodesOfType(
        TGraph, VertCurv, Geo::DartOrthoDomain, MeshBox.Diag() * 0.01, PDParam);
  }

  if (DrawJumpNodes) {
    PDParam.ColLine = Geo::Point3<ScalarType>(1, 1, 0);
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodesOfType(
        TGraph, VertCurv, Geo::JumpNodes, MeshBox.Diag() * 0.01, PDParam);
  }

  if (DrawInternalNode) {
    PDParam.ColLine = Geo::Point3<ScalarType>(1, 0, 1);
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodesOfType(
        TGraph, VertCurv, Geo::Internal, MeshBox.Diag() * 0.01, PDParam);
  }

  if (DrawAssignValues)
    GLGraphDraw::GLDrawPatchSideDiscretization<ScalarType>::DrawSideValue(
        SideD);

  if (DrawSubdivision)
    GLGraphDraw::GLDrawPatchSideDiscretization<ScalarType>::DrawSubdivision(
        SideD);

  if (DrawQuadMesh) {
    GLDraw::DrawSurfaceMesh<ScalarType>(QuadVert, QuadFaces, QuadFaceNormals,
                                        QuadVertNormals, true, true, false,
                                        false, GLDraw::TrQuad, &QuadFaceColors);

    if (showFeatures)
      GLDraw::GLDrawEdges<ScalarType>(QuadVert, QuadFeatures, 15,
                                      Geo::Point3<ScalarType>(1, 0, 1));

    if (Separatrices.size() > 0)
      GLDraw::GLDrawEdgeGroups<ScalarType>(QuadVert, Separatrices, 15);
  }

  assert((glGetError()) == GL_NO_ERROR);

}

void InitDefaultParameters() {
  DualParam.OnlyNonQuad = true;
  DualParam.InitSingPos = true;
  DualParam.FixSingPos = false;
  DualParam.JointBorder = true;
}

int main(int argc, char *argv[]) {

  InitDefaultParameters();
  // oss<<"Initialized Log Window"<<std::endl;
  // AScene.SetOutputStream(&oss);

  assert(argc >= 3);
  if (argc < 3) {
    std::cout << "ERROR need mesh and field file" << std::endl;
    exit(0);
  }

  loadedMesh = false;
  for (int i = 1; i < argc; i++) {

    if (GetFileExtension(std::string(argv[i])) == std::string("obj")) {
      PathMesh = std::string(argv[i]);
      loadedMesh = LoadMesh(PathMesh);
      if (!loadedMesh)
        exit(0);
      AVLen = AvgEdgeLen(VertPos, Connectivity);
      Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
          VertPos, Connectivity);

      ProjectName = RemoveExtension(PathMesh);

      Q_LUT.Initialize();

      continue;
    }

    if (GetFileExtension(std::string(argv[i])) == std::string("field")) {
      PathField = std::string(argv[i]);
      has_cross_field = LoadField(PathField);
      if (!has_cross_field)
        exit(0);
      continue;
    }

    if (GetFileExtension(std::string(argv[i])) == std::string("features")) {
      std::string PathFeat = std::string(argv[i]);
      bool load_feat = ReadFeatures(PathFeat, Features, true);
      if (!load_feat)
        exit(0);
      std::cout << "Loaded Features" << std::endl;

      Geo::MeshFeatures<ScalarType>::Corners(VertPos, Features, 45, CornersIdx);
      // has_load_features = true;
      continue;
    }
  }

  AutoQuadEdgeSize();

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

    SetDecomposeBar();
    int windowWidth, windowHeight;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);

    SetRenderBar();

    glfwGetWindowSize(window, &windowWidth, &windowHeight);

    SetExportBar();

    GLDrawMesh();

    ImGui::Render();

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
