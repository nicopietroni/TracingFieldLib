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
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "IO/io_obj.h"
#include <GLFW/glfw3.h>
#include <GL_functions/GLDrawMesh.h>
#include <box3.h>
#include <deformable_icp.h>
#include <icp.h>
#include <point3.h>
#include <ransac_align.h>

GLFWwindow *window;
std::string PathMesh0, PathMesh1;
UI::OrbitBall trackball;

ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
ImVec4 mesh_color0 = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
ImVec4 mesh_color1 = ImVec4(1.0f, 0.62f, 0.87f, 1.00f);

// float line_thick = 1.f;
// int selected_rend_item = 5;

bool showMesh0 = true;
bool showMesh1 = true;

typedef double ScalarType;

// MESH DATA
bool loadedMesh = false;
std::vector<Geo::Point3<ScalarType>> VertPos0, VertPos1;
std::vector<std::vector<int>> Connectivity0, Connectivity1;

std::vector<Geo::Point3<ScalarType>> FaceNormals0, FaceNormals1;
std::vector<Geo::Point3<ScalarType>> VertNormals0, VertNormals1;

Geo::Box3<ScalarType> MeshBox0;

Geo::ICP<ScalarType> ICPSolver(VertPos0, VertPos1, Connectivity0,
                               Connectivity1);
Geo::Deformable_ICP<ScalarType>
    Deformable_ICPSolver(VertPos0, VertPos1, Connectivity0, Connectivity1);

Geo::AlignRansac<ScalarType> RansacSolve(VertPos0, VertPos1, Connectivity0,
                                         Connectivity1);

void AlignRigid() {
  ICPSolver.Init();
  ICPSolver.IterateICP();
  ComputeNormals(VertPos1, Connectivity1, FaceNormals1, VertNormals1);
}

void AlignDeformable() {
  Deformable_ICPSolver.Init();
  Deformable_ICPSolver.DeformableAlignStep();
  ComputeNormals(VertPos1, Connectivity1, FaceNormals1, VertNormals1);
}

// void ColorByDistance() {
//   std::vector<ScalarType> VertDist0;
//   Geo::OneWayHausdorff<ScalarType>(VertPos0, Connectivity0, VertPos1,
//                                    Connectivity1, VertDist0);

//   std::vector<ScalarType> VertDist1;
//   Geo::OneWayHausdorff<ScalarType>(VertPos1, Connectivity1, 
//                                    VertPos0,Connectivity0, 
//                                    VertDist1);

//   ScalarType MaxDistance = MeshBox.Diag() * maxErrRatio;
//   std::vector<Geo::Point3<ScalarType>> VertColor0, VertColor1;
//   GetColorByScalar<ScalarType>(VertDist0, MaxDistance, 0, VertColor0);
//   GetColorByScalar<ScalarType>(VertDist1, MaxDistance, 0, VertColor1);

//   ElemColorFromVert<ScalarType>(Connectivity0, VertColor0, FaceCol0);
//   ElemColorFromVert<ScalarType>(Connectivity1, VertColor1, FaceCol1);
// }

void SaveMeshes() {
  WriteOBJ("./mesh0.obj",VertPos0, Connectivity0);
  WriteOBJ("./mesh1.obj",VertPos1, Connectivity1);
}

void AlignRansac() {
  RansacSolve.Init();
  RansacSolve.Align(3, 10000);
  ComputeNormals(VertPos1, Connectivity1, FaceNormals1, VertNormals1);
}

void LoadMeshes(const std::string &path0, const std::string &path1) {
  //    oss<<"Loading "<<path.c_str()<<std::endl;
  loadedMesh = LoadOBJ(path0, VertPos0, Connectivity0, true);
  loadedMesh &= LoadOBJ(path1, VertPos1, Connectivity1, true);
  if (loadedMesh) {
    std::cout << "Mesh 0 Loaded " << Connectivity0.size() << " faces "
              << VertPos0.size() << " vertices" << std::endl;
    std::cout << "Mesh 1 Loaded " << Connectivity1.size() << " faces "
              << VertPos1.size() << " vertices" << std::endl;

    MeshBox0.Init(VertPos0);
    ComputeNormals(VertPos0, Connectivity0, FaceNormals0, VertNormals0);
    ComputeNormals(VertPos1, Connectivity1, FaceNormals1, VertNormals1);
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
  //        if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_PRESS)
  //            trackball.ButtonDown(vcg::Trackball::KEY_SHIFT);
  //        if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_RELEASE)
  //            trackball.ButtonUp(vcg::Trackball::KEY_SHIFT);

  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
    trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);
  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
    trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);

  //        if (key == GLFW_KEY_LEFT_ALT && action == GLFW_PRESS)
  //            trackball.ButtonDown(vcg::Trackball::KEY_ALT);
  //        if (key == GLFW_KEY_LEFT_ALT && action == GLFW_RELEASE)
  //            trackball.ButtonUp(vcg::Trackball::KEY_ALT);
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
  ImGui::Begin(
      "ICP", NULL,
      ImGuiWindowFlags_AlwaysAutoResize); // Create a window called "Hello,
                                          // world!" and append into it.

  // if (ImGui::Button("Init Icp")) {
  //   InitICP();
  // }
  ImGui::Checkbox("Show Mesh0",
                  &showMesh0); // Create a checkbox to toggle `my_bool`
  ImGui::Checkbox("Show Mesh1", &showMesh1);

  if (ImGui::Button("Align Ransac Rigid")) {
    AlignRansac();
  }

  if (ImGui::Button("Align ICP Rigid")) {
    AlignRigid();
  }

  // if (ImGui::Button("Init Laplacian Icp")) {
  //   InitLaplacianICP();
  // }

  if (ImGui::Button("Align Deformable")) {
    AlignDeformable();
  }

  if (ImGui::Button("Rigid + Deformable")) {
    AlignRigid();
    AlignDeformable();
  }

  if (ImGui::Button("Save Mesh"))
    SaveMeshes();
  
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

  GLDraw::glScale(3 / MeshBox0.Diag());
  GLDraw::glTranslate(-MeshBox0.Center());

  bool DrawSurface = true;

  // mesh.GLDrawBox();
  if (showMesh0) {
    glColor4f(mesh_color0.x * mesh_color0.w, mesh_color0.y * mesh_color0.w,
              mesh_color0.z * mesh_color0.w, mesh_color0.w);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    // glLineWidth(line_thick);
    GLDraw::DrawSurfaceMesh<ScalarType>(VertPos0, Connectivity0, FaceNormals0,
                                        VertNormals0, DrawSurface, true, true,
                                        false);
    glPopAttrib();
  }

  if (showMesh1) {
    glColor4f(mesh_color1.x * mesh_color1.w, mesh_color1.y * mesh_color1.w,
              mesh_color1.z * mesh_color1.w, mesh_color1.w);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    // glLineWidth(line_thick);
    GLDraw::DrawSurfaceMesh<ScalarType>(VertPos1, Connectivity1, FaceNormals1,
                                        VertNormals1, DrawSurface, true, true,
                                        false);
    glPopAttrib();
  }
}

int main(int argc, char *argv[]) {

  // oss<<"Initialized Log Window"<<std::endl;
  // AScene.SetOutputStream(&oss);
  if (argc != 3) {
    std::cout << "Need 2 Meshes as arguments" << std::endl;
    exit(0);
  }
  PathMesh0 = std::string(argv[1]);
  PathMesh1 = std::string(argv[2]);
  LoadMeshes(PathMesh0, PathMesh1);

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
