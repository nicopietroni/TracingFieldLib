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
#include <cassert>
#include <cstddef>
#include <iostream>
#include <stdio.h>
#include <utility>
#include <vector>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "IO/io_obj.h"
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
#include <mesh_subdivide.h>
#include <normal.h>
#include <point3.h>
#include <triangular_remesh.h>
#include <color_functions.h>
#include "tangent_space_smooth.h"

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

// positions in 3D and UV
std::vector<typename Geo::Point3<ScalarType>> VertPos, VertPos0;
std::vector<std::vector<int>> Connectivity, Connectivity0;

std::vector<Geo::Point3<ScalarType>> FaceNormals;
std::vector<Geo::Point3<ScalarType>> VertNormals;

// Bounding Box in 3D
Geo::Box3<ScalarType> MeshBox;

// Geo::ImproveEnergyPriority<ScalarType> IEPrio;
// Geo::IterativeOptimizer<ScalarType> IteOpt(VertPos, Connectivity,&IEPrio);

ScalarType AVGL = 0;
ScalarType BBDiag = 0;
// int NextV0 = 0, NextV1 = 0;
std::vector<Geo::Point3<ScalarType>> VertCol;
std::vector<Geo::Point3<ScalarType>> FaceCol;

void InitConstantColor() {
  VertCol.clear();
  VertCol.resize(VertPos.size(), Geo::Point3<ScalarType>(0.9, 0.9, 0.9));
}

void ColorByIndex(std::vector<int> &Index) {
  assert(Index.size()==Connectivity.size());
  FaceCol.clear();
  int MaxIdx=0;
  for (size_t i=0;i<Index.size();i++)
  MaxIdx=std::max(MaxIdx,Index[i]);

  for (size_t i=0;i<Index.size();i++)
    FaceCol.push_back(Geo::ScatterColor<ScalarType>(MaxIdx,Index[i]));
}

Geo::TriRemParam<ScalarType> RemP;

ScalarType initial_feature_angle = 45;
ScalarType old_initial_feature_angle = initial_feature_angle;

std::vector<std::pair<int, int>> RenderVertFeatures;
std::vector<int> RenderVertCorners;

int SelKernelSize = 2;

void UpdateRenderCorners() {
  assert(RemP.IsFaceEFeature.size() == Connectivity.size());
  assert(RemP.IsVertCorner.size() == VertPos.size());

  RenderVertFeatures.clear();
  RenderVertCorners.clear();
  Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(Connectivity, RemP.IsFaceEFeature, RenderVertFeatures);
  Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(RemP.IsVertCorner,RenderVertCorners);
}

// void UpdateSelectedArea() {

//   std::vector<int> KernelF;
//   for (size_t i = 0; i < RemP.IsFaceEFeature.size(); i++)
//     for (size_t j = 0; j < RemP.IsFaceEFeature[i].size(); j++)
//       if (RemP.IsFaceEFeature[i][j])
//         KernelF.push_back(i);

//   std::sort(KernelF.begin(), KernelF.end());
//   KernelF.erase(std::unique(KernelF.begin(), KernelF.end()), KernelF.end());

//   std::vector<std::vector<int>> StarF;
//   FindFaceStar(Connectivity, KernelF, StarF, SelKernelSize);

//   //add the faces to
//   RemP.Exclude= std::vector<bool>(Connectivity.size(), true);

//   //then update selection
//   for (size_t i=0;i<KernelF.size();i++)
//   {
//     int indexF=KernelF[i];
//     for (size_t j=0;j<StarF[indexF].size();j++)
//       RemP.Exclude[StarF[indexF][j]]=false;
//   }

//   //finally update color
//   InitConstantColor();
//   for (size_t i=0;i<Connectivity.size();i++)
//   {
//     if (RemP.Exclude[i])continue;

//   }
// }

void UpdateFeatureCorners() {
  Geo::MeshFeatures<ScalarType>::FeatureCorners(VertPos, Connectivity, initial_feature_angle,
                      RemP.IsFaceEFeature, RemP.IsVertCorner);

  UpdateRenderCorners();
}

void LoadMesh(const std::string &path) {
  loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  if (loadedMesh) {
    std::cout << "Loaded " << Connectivity.size() << " faces " << VertPos.size()
              << " vertices" << std::endl;

    MeshBox.Init(VertPos);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    UpdateFeatureCorners();

    AVGL = AvgEdgeLen(VertPos, Connectivity);
    BBDiag = MeshBox.Diag();

    VertPos0 = VertPos;
    Connectivity0 = Connectivity;
    // IteOpt.NextCollapseTest(NextV0, NextV1);
    InitConstantColor();
    // RemP.Exclude = std::vector<bool>(Connectivity.size(), false);
  } else
    std::cout << "Mesh not loaded" << std::endl;
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

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {

    trackball.MouseDown((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
  }
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {

    trackball.MouseUp((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
  }
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
    std::cout << "Error Initializing GLFW" << std::endl;
    exit(0);
  }
  // return 1;

  window = glfwCreateWindow(1280, 720, "UV Mapping", NULL, NULL);

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

float TargetRatio = 1.f;
// int selected_target = 0;
// Geo::OperationPriorityType PType=Geo::EqRigidDist;
// Geo::ImproveAspectRatioPriority<ScalarType>::AspectRatioType AType;

void SetRenderBar() {
  ImGui::Begin("Mesh Optimizer", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  int valueA0 =
      floor(0.5 + initial_feature_angle); // Variable to store slider value
  ImGui::SetNextItemWidth(100.0f);
  ImGui::SliderInt("Inital features", &valueA0, 0, 180);
  initial_feature_angle = valueA0;
  if (old_initial_feature_angle != initial_feature_angle) {
    UpdateFeatureCorners();
    old_initial_feature_angle = initial_feature_angle;
  }

  ImGui::SetNextItemWidth(100.0f);
  ImGui::InputFloat("target res (wrt average edge)", &TargetRatio);
  TargetRatio = std::max(TargetRatio, 0.2f);
  TargetRatio = std::min(TargetRatio, 20.f);
  ImGui::SetNextItemWidth(100.0f);
  ImGui::InputInt("optimization steps", &RemP.steps);
  RemP.steps = std::max(RemP.steps, 1);
  RemP.steps = std::min(RemP.steps, 20);
  ImGui::SetNextItemWidth(100.0f);
  ImGui::Checkbox("iterative", &RemP.iterative);

  int valueA1 =
      floor(0.5 + RemP.feature_angle); // Variable to store slider value
  ImGui::SetNextItemWidth(100.0f);
  ImGui::SliderInt("dynamic features angle", &valueA1, 0, 180);
  RemP.feature_angle = valueA1;

  ImGui::SetNextItemWidth(100.0f);
  ImGui::Checkbox("adaptive", &RemP.adaptive);

  ImGui::SetNextItemWidth(100.0f);

  if (ImGui::Button("DENSITY MAP")) {
    std::vector<ScalarType> VertESize;
    ScalarType Radius = AVGL * 2;
    FindDensityMap(VertPos, Connectivity, Radius, VertESize);
    Geo::GetColorByScalar<ScalarType>(VertESize, 0, 1, VertCol);
  }

  if (ImGui::Button("REMESH")) {
    RemP.TargetL = AVGL * TargetRatio;
    // std::cout << "Target E Len: " << target_l << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    ScalarType Err = MeshBox.Diag() * 0.01;
    RemP.max_error = Err;
    RemP.writeMSG = true;
    Geo::TriRemesh<ScalarType>(VertPos, Connectivity,
                               RemP); //,target_l,AType,60,Err);

    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    InitConstantColor();
    UpdateRenderCorners();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count();

    // Output the duration
    std::cout << std::endl << "*** FINAL STATS ***" << std::endl << std::endl;

    std::cout << "Elapsed time: " << duration / (ScalarType)1000000
              << " seconds" << std::endl
              << std::endl;

    std::cout << "Vert Num: " << VertPos.size() << std::endl;
    std::cout << "Face Num: " << Connectivity.size() << std::endl << std::endl;

    ScalarType minLen, maxLen;
    ExtremenEdgeLen(VertPos, Connectivity, minLen, maxLen);

    ScalarType optLen = AvgEdgeLen(VertPos, Connectivity);
    std::cout << "Err AVG E Len (% of Target L): "
              << 100 * fabs(RemP.TargetL - optLen) / RemP.TargetL << " %"
              << std::endl;
    std::cout << "Err Min E Len (% of Target L): "
              << 100 * fabs(RemP.TargetL - minLen) / RemP.TargetL << " %"
              << std::endl;
    std::cout << "Err Max E Len (% of Target L): "
              << 100 * fabs(RemP.TargetL - maxLen) / RemP.TargetL << " %"
              << std::endl
              << std::endl;

    ScalarType AvgDist, MaxDist;
    Geo::Hausdorff(VertPos0, Connectivity0, VertPos, Connectivity, AvgDist,
                   MaxDist);

    std::cout << "Avg Haus Distance (% of BB Diag): " << 100 * AvgDist / BBDiag
              << " %" << std::endl;
    std::cout << "Max Haus Distance (% of BB Diag): " << 100 * MaxDist / BBDiag
              << " %" << std::endl
              << std::endl;

    ScalarType WorstQ =
        Geo::TriMeshWorstQuality<ScalarType>(VertPos, Connectivity);
    ScalarType AvgQ = Geo::TriMeshAvgQuality<ScalarType>(VertPos, Connectivity);

    std::cout << "Avg Quality: " << AvgQ << std::endl;
    std::cout << "Worst Quality: " << WorstQ << std::endl << std::endl;

    WriteOBJ("test.obj", VertPos, Connectivity);
  }

  if (ImGui::Button("REPROJECT")) {
    std::vector<int> OldF;
    Geo::SubdivideMesh<ScalarType>::MidPointReproject(
        VertPos0, Connectivity0, VertPos, Connectivity, OldF, 2);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    InitConstantColor();
    ColorByIndex(OldF);
    // UpdateRenderCorners();
    RenderVertFeatures.clear();
    RenderVertCorners.clear();
  }
  if (ImGui::Button("SMOOTH")) {
    Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam SMParam;
    SMParam.ReprojVert=&VertPos0;
    SMParam.ReprojFaces=&Connectivity0;
    SMParam.post_project=true;
    SMParam.Steps=10;
    Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos,Connectivity,SMParam);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
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
    { SetRenderBar(); }

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
