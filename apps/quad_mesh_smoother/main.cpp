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
#include "mesh_functions.h"
#include <GL/glew.h>
#include <UI/orbitball.h>
#include <array>
#include <cassert>
#include <cstddef> // C++
#include <iostream>
#include <stdio.h>
#include <vector>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "IO/file_path_functions.h"
#include "IO/io_obj.h"
#include <GLFW/glfw3.h>
#include <GL_functions/GLDrawMesh.h>
#include <box3.h>
#include <chrono>
#include <connectivity.h>
#include <iostream>
#include <mesh_features.h>
#include <mesh_functions.h>
#include <point3.h>
#include <smooth_mesh.h>

GLFWwindow *window;
std::string PathMesh;
std::string PathField;
UI::OrbitBall trackball;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
ImVec4 mesh_color = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
float line_thick = 1.f;
int selected_rend_item = 5;

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
ScalarType AvEdge = 0;

// MESH DATA
bool loadedMesh = false;
std::vector<Geo::Point3<ScalarType>> VertPos;
std::vector<std::vector<int>> Connectivity;
std::vector<Geo::Point3<ScalarType>> FaceNormals;
std::vector<Geo::Point3<ScalarType>> VertNormals;
std::vector<std::vector<int>> NextF, NextE;
Geo::Box3<ScalarType> MeshBox;

std::vector<std::pair<int, int>> VertPairFeatures;
std::vector<int> CornersIdx;
int ero_steps = 2;
int dil_steps = 2;
bool include_border = true;
ScalarType feature_angle = 30;

ScalarType old_feature_angle = feature_angle;
int old_ero_steps = ero_steps;
int old_dil_steps = dil_steps;

bool ShowCorners = true;
bool ShowFeatures = true;
bool ReprojectMesh = true;

void InitFeatures() {
  Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
      VertPos, Connectivity, feature_angle, VertPairFeatures, CornersIdx,
      ero_steps, dil_steps, include_border);
}

std::vector<Geo::Point3<ScalarType>> TriVertPos;
std::vector<std::vector<int>> TriConnectivity;

void InitReprojectMesh() {
  SplitQuadIntoTris<ScalarType>(VertPos, Connectivity, TriConnectivity);
  TriVertPos = VertPos;
}

void SaveMesh() {
  // assert(!PathMesh.empty());
  // std::string PathRemMesh = RemoveExtension(PathMesh);
  // PathRemMesh += std::string("_remesh.obj");
  // std::cout << "Saving Field on file " << PathRemMesh.c_str() << std::endl;
  // bool Saved = WriteOBJ(PathRemMesh, VertPos, Connectivity);
  // if (!Saved) {
  //   std::cout << "ERROR SAVING FIELD FILE " << std::endl;
  //   exit(0);
  // }
}

Geo::SmoothParam<ScalarType> SParamQuadSmooth;

void SmoothMesh() {

  SParamQuadSmooth.Clear();

  SParamQuadSmooth.Features = VertPairFeatures;

  // std::vector<std::pair<int, int>> BorderEQuad;
  // Geo::getBorderEdges(QuadFaces, BorderEQuad);
  // for (size_t i = 0; i < BorderEQuad.size(); i++) {
  //   SParamQuadSmooth.Features.push_back(BorderEQuad[i]);
  // }

  // SParamQuadSmooth.ReprojFeatures = TriFeatures0;

  // std::vector<std::pair<int, int>> BorderETris;
  // Geo::getBorderEdges(TriConnectivity0, BorderETris);
  // for (size_t i = 0; i < BorderETris.size(); i++) {
  //   SParamQuadSmooth.ReprojFeatures.push_back(BorderETris[i]);
  // }

  // // add the subdivided border as features
  // SParamQuadSmooth.ReprojVert = &TriVertPos0;
  // SParamQuadSmooth.ReprojFaces = &TriConnectivity0;

  SParamQuadSmooth.Damp = 0.5;

  // if (reliable_smooth)
  // SParamQuadSmooth.post_project = !(has_problematic_patches);
  // else
  //   SParamQuadSmooth.post_project = false;

  // SParamQuadSmooth.NumIte = QuadSmoothSteps;
  // if (reliable_smooth)
  //   SParamQuadSmooth.SType = Geo::SMBestPoly;
  // else
  //   SParamQuadSmooth.SType = Geo::SMLaplacian;

  // Geo::SMCentroid;//Geo::SMBestPoly;//Geo::SMLaplacian;//Geo::SMTaubin;/

  // SParamQuadSmooth.line_border = true;
  // SParamQuadSmooth.check_AspectR = false;
  // SParamQuadSmooth.line_neight_features = false;
  // SParamQuadSmooth.reproject_voronoi_partition = reliable_smooth;
  // SParamQuadSmooth.check_multi_pass_fold = false;

  SParamQuadSmooth.FixedVert = CornersIdx;
  // std::sort(SParamQuadSmooth.FixedVert.begin(),
  //           SParamQuadSmooth.FixedVert.end());

  // std::vector<int> Valence1Idx;
  // Geo::MeshFeatures<ScalarType>::VerticesWithValence1(QuadFeatures,
  //                                                     Valence1Idx);
  // std::sort(Valence1Idx.begin(), Valence1Idx.end());

  // // std::cout<<"There are "<<Valence1Idx.size()<<" valence 1
  // // vert"<<std::endl;
  // std::vector<int> resultDiff;
  // std::set_difference(SParamQuadSmooth.FixedVert.begin(),
  //                     SParamQuadSmooth.FixedVert.end(), Valence1Idx.begin(),
  //                     Valence1Idx.end(), std::back_inserter(resultDiff));
  // SParamQuadSmooth.FixedVert = std::move(resultDiff);

  // SetFixQuadVertices(SParamQuadSmooth);
  Geo::SmoothLaplacian<ScalarType>(VertPos, Connectivity, SParamQuadSmooth);

  // if (!reliable_smooth)
  //   return;

  // SParamQuadSmooth.SType = Geo::SMCentroid;
  // SParamQuadSmooth.reproject_voronoi_partition = true;
  // SParamQuadSmooth.check_multi_pass_fold = false;
  // SParamQuadSmooth.check_AspectR = true;
  // Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
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
    AvEdge = AvgEdgeLen(VertPos, Connectivity);
    // DesambiguateSamePosVertices(VertPos,Connectivity,0.0001);
    return true;
    // std::cout<<"FNorm 0 "<<FaceNormals.size()<<std::endl;
  } else {
    std::cout << "Mesh not loaded" << std::endl;
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
  trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods) {

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

// bool lineF = (SParam.NDir == 2);
// bool crossF = (SParam.NDir == 4);

// std::vector<std::array<int, 3>> TripleVert;
// std::vector<int> DisconnVert;
// std::vector<Geo::Point3<ScalarType>> DisconnDir;
void SetRenderBar() {
  ImGui::Begin(
      "QuadMesh Smoother", NULL,
      ImGuiWindowFlags_AlwaysAutoResize); // Create a window called "Hello,
                                          // world!" and append into it.

  // if (ImGui::Button("SMOOTH MESH")) {
  //   SmoothMesh();
  // }
  // ImGui::SameLine();
  // if (ImGui::Button("SMOOTH MESHING")) {
  //   SmoothTessellation();
  // }
  ImGui::Text("Render");

  ImGui::Checkbox("Show Features", &ShowFeatures);
  ImGui::Checkbox("Show Corners", &ShowCorners);

  ImGui::Separator();

  ImGui::Text("Features");

  int valueA0 = floor(0.5 + feature_angle); // Variable to store slider value
  ImGui::SetNextItemWidth(100.0f);
  ImGui::SliderInt("Inital features", &valueA0, 0, 180);
  feature_angle = valueA0;

  ImGui::SliderInt("Ero Steps", &ero_steps, 0, 5);
  ImGui::SliderInt("Dil Steps", &dil_steps, 0, 5);

  if ((old_feature_angle != feature_angle) || (old_ero_steps != ero_steps) ||
      (old_dil_steps != dil_steps)) {
    InitFeatures();
    old_feature_angle = feature_angle;
    old_ero_steps = ero_steps;
    old_dil_steps = dil_steps;
  }

  ImGui::Separator();

  ImGui::Text("Smooth");

  const char *smoothTypeNames[] = {"Laplacian", "Angle Based", "Taubin", "Best Poly",
                                   "Centroid"};

  int currentSmoothType = static_cast<int>(SParamQuadSmooth.SType);

  if (ImGui::Combo("Smooth Type", &currentSmoothType, smoothTypeNames,
                   IM_ARRAYSIZE(smoothTypeNames))) {
    SParamQuadSmooth.SType = static_cast<Geo::SmoothType>(currentSmoothType);
  }

  ImGui::Checkbox("Post Project", &SParamQuadSmooth.post_project);
  ImGui::Checkbox("Reproject Mesh", &ReprojectMesh);
  ImGui::Checkbox("Reproject Tangent", &SParamQuadSmooth.reproj_tangent);

  if (ReprojectMesh || SParamQuadSmooth.post_project) {
    SParamQuadSmooth.ReprojVert = &TriVertPos;
    SParamQuadSmooth.ReprojFaces = &TriConnectivity;
    SParamQuadSmooth.ReprojFeatures = VertPairFeatures;
  } else {
    SParamQuadSmooth.ReprojVert = NULL;
    SParamQuadSmooth.ReprojFaces = NULL;
    SParamQuadSmooth.ReprojFeatures.clear();
  }
  ImGui::Checkbox("Check Aspect Ratio", &SParamQuadSmooth.check_AspectR);
  ImGui::Checkbox("Check Multi-Pass Fold",&SParamQuadSmooth.check_multi_pass_fold);
  // ImGui::Checkbox("Line Neighbor Features",
  // &SParamQuadSmooth.line_neight_features);
  ImGui::Checkbox("Reproject Voronoi Partition", &SParamQuadSmooth.reproject_voronoi_partition);
  ImGui::SliderInt("Smooth Steps", &SParamQuadSmooth.NumIte, 1, 30);

  // template <class ScalarType> struct SmoothParam {
  // int NumIte;
  // std::vector<int> FixedVert;
  // std::vector<std::pair<int, int>> Features;
  // bool line_border;
  // bool post_project;
  // ScalarType Damp;
  // std::vector<Geo::Point3<ScalarType>> *ReprojVert;
  // std::vector<std::vector<int>> *ReprojFaces;
  // std::vector<std::pair<int, int>> ReprojFeatures;
  // int ReprojStep;
  // SmoothType SType;
  // bool check_AspectR;
  // bool no_reproj_feat_endpoints;
  // bool line_neight_features;
  // bool reproject_voronoi_partition;
  // bool check_multi_pass_fold;
  // ScalarType reproj_dist_THR;
  // std::vector<Geo::Point3<ScalarType>> PrevVertPos;

  if (ImGui::Button("SMOOTH")) {
    auto start = std::chrono::steady_clock::now();

    SmoothMesh();

    auto end = std::chrono::steady_clock::now();

    double seconds = std::chrono::duration<double>(end - start).count();

    std::cout << "Elapsed: " << seconds << " s\n";
  }
  if (ImGui::Button("RESET")) {
    for (size_t i = 0; i < TriVertPos.size(); i++)
      VertPos[i] = TriVertPos[i];
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  }

  // if (ImGui::Button("SPLIT CORNERS")) {
  //   std::vector<ScalarType> Angles;
  //   // 180-45
  //   // 180+45
  //   Angles.push_back(135);
  //   Angles.push_back(225);
  //   //
  //   Geo::GetAdjFeatureAngle<ScalarType>(VertPos,Connectivity,Features,Angles,5,TripleVert);
  //   RefineMesh<ScalarType>::SubdivideAdjFeatureAngle(VertPos, Connectivity,
  //                                                    Features, Angles, 5);
  //   Geo::getFFAdjacency(Connectivity, NextF, NextE);
  //   ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

  //   std::cout << "found " << TripleVert.size() << " triplets" << std::endl;
  // }
  // ImGui::SameLine();
  // if (has_cross_field && loadedMesh) {
  //   if (ImGui::Button("SAVE ALL")) {
  //     SaveMesh();
  //     SaveField();
  //     SaveFeatures();
  //   }
  // }

  // ImGui::Separator();
  // ImGui::Text("Mesh");

  // if (!has_load_features) {
  //   int valueA0 =
  //       floor(0.5 + initial_feature_angle); // Variable to store slider value
  //   ImGui::SetNextItemWidth(100.0f);
  //   ImGui::SliderInt("Inital features", &valueA0, 0, 180);
  //   initial_feature_angle = valueA0;

  //   ImGui::SliderInt("Ero Steps", &ero_steps, 0, 5);
  //   ImGui::SliderInt("Dil Steps", &dil_steps, 0, 5);

  //   if ((old_initial_feature_angle != initial_feature_angle)||
  //       (old_ero_steps != ero_steps)||
  //       (old_dil_steps != dil_steps))
  //   {
  //     UpdateFeatureCorners();
  //     old_initial_feature_angle = initial_feature_angle;
  //     old_ero_steps = ero_steps;
  //     old_dil_steps = dil_steps;
  //   }
  //   // if (old_ero_dil_steps != ero_dil_steps) {
  //   //   UpdateFeatureCorners();
  //   //   old_ero_dil_steps = ero_dil_steps;
  //   // }
  // }

  // if (ImGui::Button("Refine For Tracing")) {
  //   // SplitBorderEdges();
  //   RefineForTracing();
  // }

  // if (has_cross_field)
  //   if (ImGui::Button("Refine for Connectivity")) {
  //     RefineForTracingConnectivity();
  //     // Geo::RefineForTracing<ScalarType>::RefineForConnection(
  //     //     VertPos, Connectivity, FaceCurv, Features);
  //     // Geo::getFFAdjacency(Connectivity, NextF, NextE);
  //     // ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  //     //
  //     Geo::RefineForTracing<ScalarType>::FindNonConnectedVert(VertPos,Connectivity,FaceCurv,Features,
  //     // DisconnVert,DisconnDir);
  //   }

  // // ImGui::SameLine();
  // // if (ImGui::Button("BATCH PROCESS BORD")) {
  // //   BatchProcessBord();
  // // }

  // ImGui::Separator();
  // ImGui::Text("Init");

  // ImGui::SliderInt("Ring", &Nring, 2, 5, "%d", ImGuiSliderFlags_AlwaysClamp);
  // if (ImGui::Button("Init By Curvature")) {
  //   if (loadedMesh)
  //     InitFieldByCurvature();
  // }

  // ImGui::Separator();
  // ImGui::Text("Smooth");
  // if (ImGui::RadioButton("Line F", lineF)) {
  //   lineF = true;
  //   crossF = false;
  //   SParam.NDir = 2;
  // }
  // ImGui::SameLine();
  // if (ImGui::RadioButton("Cross F", crossF)) {
  //   lineF = false;
  //   crossF = true;
  //   SParam.NDir = 4;
  // }

  // ImGui::Checkbox("Fix Border", &SParam.align_borders);
  // // ImGui::Checkbox("Show Sing", &showSing);

  // // Create a slider for the float value
  // ImGui::SliderFloat("Top Fix", &TopIntevalFix, 0.0f, 0.3f);
  // SParam.TopIntevalFix = (ScalarType)TopIntevalFix;

  // //        ImGui::SliderFloat("Top Fix",&TopIntevalFix, 0, 0.3,
  // //        "%2.f",ImGuiSliderFlags_AlwaysClamp); SParam.TopIntevalFix=
  // //        (ScalarType) TopIntevalFix;

  // ImGui::SliderFloat("Smoothness", &GlobalSmoothVal, 50, 500, "%1.f",
  //                    ImGuiSliderFlags_AlwaysClamp);
  // SParam.GlobalSmoothVal = GlobalSmoothVal;

  // if (has_cross_field && loadedMesh) {
  //   if (ImGui::Button("Smooth Field")) {
  //     if (loadedMesh)
  //       SmoothGlobal();
  //   }
  //   // ImGui::SameLine();
  //   ImGui::Separator();
  //   ImGui::Text("Save");

  //   if (ImGui::Button("Save Mesh")) {
  //     SaveMesh();
  //   }

  //   ImGui::SameLine();
  //   if (ImGui::Button("Save Field")) {
  //     SaveField();
  //   }

  //   ImGui::SameLine();
  //   if (ImGui::Button("Save Features")) {
  //     SaveFeatures();
  //   }
  // }

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
  // if (showTrackBall)
  //     trackball.DrawPostApply();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());

  // bool DrawSurface = true;

  // mesh.GLDrawBox();
  // if (showMesh) {
  glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w,
            mesh_color.z * mesh_color.w, mesh_color.w);
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glLineWidth(line_thick);

  GLDraw::DrawSurfaceMesh<ScalarType>(VertPos, Connectivity, FaceNormals,
                                      VertNormals, true, true, false, false,
                                      GLDraw::TrQuad);

  GLDraw::GLDrawBorders<ScalarType>(VertPos, Connectivity, NextF, 5);

  if (ShowCorners)
    GLDraw::DrawVertices(VertPos, CornersIdx, Geo::Point3<ScalarType>(1, 0, 0),
                         30);

  if (ShowFeatures)
    GLDraw::GLDrawEdges<ScalarType>(VertPos, VertPairFeatures, 15,
                                    Geo::Point3<ScalarType>(1, 0, 1));

  glPopAttrib();
}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cout << "You should pass at least a mesh" << std::endl;
    exit(0);
  }

  bool has_load_mesh = false;
  if (GetFileExtension(std::string(argv[1])) == std::string("obj")) {
    PathMesh = std::string(argv[1]);
    has_load_mesh = LoadMesh(PathMesh);
    InitFeatures();
    InitReprojectMesh();
  }
  if (!has_load_mesh) {
    std::cout << "You should pass at least a mesh" << std::endl;
    exit(0);
  }

  InitGLFW_Window();

  InitIMGui();

  trackball.center = Geo::Point3<float>(0, 0, 0);
  trackball.radius = 1;
  // trackball.SetSpinnable(true);

  glewInit();
  // Main loop
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    SetRenderBar();

    // Rendering
    ImGui::Render();

    GLDrawMesh();

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
