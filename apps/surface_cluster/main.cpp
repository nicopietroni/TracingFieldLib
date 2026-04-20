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
// #include "mesh_functions.h"
#include <GL/glew.h>
#include <UI/orbitball.h>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <stdio.h>
#include <string>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <utility>
#include <vector>
// #include "ImGuiFileBrowser.h"
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "IO/file_path_functions.h"
#include "IO/io_obj.h"
#include "field_graph/dual_patches.h"
#include <GLFW/glfw3.h>
#include <GL_functions/GLDrawField.h>
#include <GL_functions/GLDrawMesh.h>
#include <IO/file_path_functions.h>
#include <box3.h>
#include <color_functions.h>
#include <point3.h>
#include <surface_cluster/face_cluster.h>

GLFWwindow *window;
std::string PathMesh;
UI::OrbitBall trackball;
// imgui_addons::ImGuiFileBrowser file_dialog;
// std::ostringstream oss;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
ImVec4 mesh_color = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
float line_thick = 1.f;
int selected_rend_item = 5;

bool showMesh = true;

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
// ScalarType AvEdge = 0;

// MESH DATA
bool loadedMesh = false;
std::vector<Geo::Point3<ScalarType>> VertPos0;
std::vector<std::vector<int>> Faces0;
std::vector<Geo::Point3<ScalarType>> VertPos;
std::vector<std::vector<int>> Faces;
std::vector<Geo::Point3<ScalarType>> FaceNormals;
std::vector<Geo::Point3<ScalarType>> FaceColor;
std::vector<std::pair<int, int>> BorderE;

std::vector<Geo::Point3<ScalarType>> VertNormals;
std::vector<std::vector<int>> NextF, NextE;
std::vector<std::vector<int>> VF, VFIndex;
Geo::Box3<ScalarType> MeshBox;

std::vector<std::vector<int>> ClusterFaces;
std::vector<int> FacesCluster;
std::vector<int> AllCorners;
std::string ProjectPath;

// void SaveMesh() {
//   // assert(!PathMesh.empty());
//   // std::string PathRemMesh = RemoveExtension(PathMesh);
//   // PathRemMesh += std::string("_remesh.obj");
//   // std::cout << "Saving Field on file " << PathRemMesh.c_str() <<
//   std::endl;
//   // bool Saved = WriteOBJ(PathRemMesh, VertPos, Connectivity);
//   // if (!Saved) {
//   //   std::cout << "ERROR SAVING FIELD FILE " << std::endl;
//   //   exit(0);
//   // }
// }

// void SavePatchLayou() {
//   // assert(!PathMesh.empty());
//   // std::string PathRemMesh = RemoveExtension(PathMesh);
//   // PathRemMesh += std::string("_remesh.obj");
//   // std::cout << "Saving Field on file " << PathRemMesh.c_str() <<
//   std::endl;
//   // bool Saved = WriteOBJ(PathRemMesh, VertPos, Connectivity);
//   // if (!Saved) {
//   //   std::cout << "ERROR SAVING FIELD FILE " << std::endl;
//   //   exit(0);
//   // }
// }

void UpdateCorners() {
  std::vector<std::vector<int>> ClusterCorners;
  Geo::ClusterOperations<ScalarType>::FindCorners(VertPos, Faces, ClusterFaces,
                                                  ClusterCorners);
  AllCorners.clear();
  for (size_t i = 0; i < ClusterCorners.size(); i++)
    AllCorners.insert(AllCorners.end(), ClusterCorners[i].begin(),
                      ClusterCorners[i].end());

  std::sort(AllCorners.begin(), AllCorners.end());
  std::vector<int>::iterator last =
      std::unique(AllCorners.begin(), AllCorners.end());
  AllCorners.erase(last, AllCorners.end());
}

void ColorByPartition() {
  std::cout << "Number of Partitions:" << ClusterFaces.size() << std::endl;
  for (size_t i = 0; i < ClusterFaces.size(); i++) {
    Geo::Point3<ScalarType> Color =
        Geo::ScatterColor<ScalarType>(ClusterFaces.size(), i);

    for (size_t j = 0; j < ClusterFaces[i].size(); j++) {
      int IndexF = ClusterFaces[i][j];
      assert(IndexF >= 0);
      assert(IndexF < FaceColor.size());
      FaceColor[IndexF] = Color;
    }
  }

  FacesCluster.resize(Faces.size(), -1);
  Geo::ClusterOperations<ScalarType>::UpdateFaceClusterFromClusterFaces(
      ClusterFaces, FacesCluster);

  Geo::ClusterOperations<ScalarType>::GetBorderCluster(FacesCluster, Faces,
                                                       NextF, BorderE);
}

int EnergyPerc = 20;
int TargetCl = -1;

bool LoadMesh(const std::string &path) {
  //    oss<<"Loading "<<path.c_str()<<std::endl;
  loadedMesh = LoadOBJ(path, VertPos, Faces, true);
  if (loadedMesh) {
    MeshBox.Init(VertPos);
    ComputeNormals(VertPos, Faces, FaceNormals, VertNormals);
    Geo::getFFAdjacency(Faces, NextF, NextE);
    Geo::getVFAdjacency(Faces, VF, VFIndex);

    FaceColor.clear();
    FaceColor.resize(Faces.size(), Geo::Point3<ScalarType>(0.9, 0.9, 1.0));
    VertPos0 = VertPos;
    Faces0 = Faces;
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
int RelaxIte = 20;

void SetRenderBar() {
  ImGui::Begin("Mesh Cluster", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  ImGui::Separator();
  ImGui::Text("Process");

  // int myInt = 0;
  ImGui::SliderInt("Num Cluster", &TargetCl, -1, 100);
  // ScalarType CurrL=
  ImGui::SliderInt("Error %", &EnergyPerc, -1, 100);
  ImGui::SliderInt("Relax Ite %", &RelaxIte, 5, 100);

  // ImGui::Checkbox("Show Cross", &showCross);
  // ImGui::Checkbox("Show Sing", &showSing);

  if (ImGui::Button("Voro Cluster")) {
    VertPos = VertPos0;
    Faces = Faces0;
    FindVoronoiClusters(VertPos, Faces, ClusterFaces, TargetCl,
                        (ScalarType)EnergyPerc, (ScalarType)-1, RelaxIte);
    ColorByPartition();
    UpdateCorners();
  }

  if (ImGui::Button("VSA Decomposition")) {
    VertPos = VertPos0;
    Faces = Faces0;
    FindVSAClusters(VertPos, Faces, ClusterFaces, TargetCl,
                    (ScalarType)EnergyPerc, (ScalarType)-1, RelaxIte);
    ColorByPartition();
    UpdateCorners();
  }

  // if (ImGui::Button("Dualize")) {
  //   std::vector<std::vector<Geo::Point3<ScalarType>>> PerClusterMeshVert;
  //   std::vector<std::vector<std::vector<int>>> PerClusterMeshFaces;
  //   std::vector<std::vector<int>> PerClusterCorners;
  //   Geo::ClusterOperations<ScalarType>::AllClusterMeshes(
  //       VertPos, Faces, ClusterFaces, AllCorners, 
  //       PerClusterMeshVert,PerClusterMeshFaces, 
  //       PerClusterCorners);

  //   std::vector<std::vector<int> > CentridsV;
  //   CentridsV.resize(PerClusterMeshFaces.size());
  //   for (size_t i=0;i<PerClusterMeshFaces.size();i++)
  //   {
  //     int CurrCentroidV=Geo::Centroids<ScalarType>::FindCentroidsVertUV(PerClusterMeshVert[i],PerClusterMeshFaces[i]);
  //     CentridsV[i].push_back(CurrCentroidV);
  //   }

  //   Geo::DualPatch<ScalarType>::PrimalPatchData PData;
  //   PData.PatchVertPos=PerClusterMeshVert;
  //   PData.PatchFaces=PerClusterMeshFaces;
  //   PData.PatchSingularities=CentridsV;
  //   assert(0);
  //   // TO BE COMPLETED

  //   // Geo::DualPatch<ScalarType>::ExtractDual(
  //   //     PData, DualParam, DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
  //   //     PerDualPatchCornerIDX);
  // }

  if (ImGui::Button("Save Patch Layout")) {
    std::string patch_path = ProjectPath + std::string(".patch");
    SaveCornerEdges(VertPos, Faces, ClusterFaces, patch_path);
    // SaveMesh();
    // SaveField();
    // SaveFeatures();
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
  // if (showTrackBall)
  //     trackball.DrawPostApply();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());

  bool DrawSurface = true;

  // mesh.GLDrawBox();
  if (showMesh) {
    glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w,
              mesh_color.z * mesh_color.w, mesh_color.w);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glLineWidth(line_thick);

    GLDraw::DrawSurfaceMesh<ScalarType>(VertPos, Faces, FaceNormals,
                                        VertNormals, DrawSurface, true, false,
                                        false, GLDraw::TRTriangle, &FaceColor);

    GLDraw::GLDrawBorders<ScalarType>(VertPos, Faces, NextF, 5);

    GLDraw::GLDrawEdges<ScalarType>(VertPos, BorderE, 15,
                                    Geo::Point3<ScalarType>(0, 0, 0));

    GLDraw::DrawVertices(VertPos, AllCorners);

    glPopAttrib();
  }
}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cout << "You should pass at least a mesh" << std::endl;
    exit(0);
  }

  bool has_load_mesh = false;
  PathMesh = std::string(argv[1]);
  has_load_mesh = LoadMesh(PathMesh);
  if (!has_load_mesh)
    exit(0);

  ProjectPath = RemoveExtension(PathMesh);
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
