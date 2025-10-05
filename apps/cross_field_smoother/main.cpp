#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include "mesh_functions.h"
#include <GL/glew.h>
#include <UI/orbitball.h>
#include <array>
#include <cassert>
#include <iostream>
#include <stdio.h>
#include <sys/_types/_size_t.h>
#include <vector>
// #include "ImGuiFileBrowser.h"
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "IO/io_obj.h"
#include <GLFW/glfw3.h>
#include <GL_functions/GLDrawField.h>
#include <GL_functions/GLDrawMesh.h>
#include <box3.h>
#include <point3.h>
// include "eigen_interface.h"
#include "IO/file_path_functions.h"
#include "IO/io_field.h"
#include <IO/io_features.h>
#include <connectivity.h>
#include <cross_field_functions.h>
#include <field_graph/refine_for_tracing.h>
#include <math_utils.h>
#include <mesh_features.h>
#include <mesh_refine.h>
#include <smooth_cross_field.h>
#include <subdivide.h>
#include <tangent_space_smooth.h>

GLFWwindow *window;
std::string PathMesh;
std::string PathField;
UI::OrbitBall trackball;
// imgui_addons::ImGuiFileBrowser file_dialog;
// std::ostringstream oss;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
ImVec4 mesh_color = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
float line_thick = 1.f;
int selected_rend_item = 5;

bool showMesh = true;
bool showCross = true;
bool showSing = true;
bool has_remeshed = false;
bool loaded_features = false;
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

// FIELD DATA
std::vector<Field::CrossF<ScalarType>> VertCurv;
std::vector<Field::CrossF<ScalarType>> FaceCurv;
bool has_cross_field = false;
bool has_load_features = false;
int Nring = 3;
ScalarType minQCrossVert, maxQCrossVert;
std::vector<int> SingIndex, SingValue;
FieldSmoothParam<ScalarType> SParam;

// Interface values
float TopIntevalFix = (float)SParam.TopIntevalFix;
float GlobalSmoothVal = (float)SParam.GlobalSmoothVal;

std::vector<std::pair<int, int>> Features;
std::vector<int> Corners;

int ero_steps=3;
int old_ero_steps=0;
int dil_steps=3;
int old_dil_steps=0;

ScalarType initial_feature_angle = 45;
ScalarType old_initial_feature_angle = 0;

bool processForTracing = true;

void InitFieldByCurvature() {
  SParam.FixedCross.clear();
  ComputeCurvatureField<ScalarType>(VertPos, Connectivity, FaceNormals,
                                    VertNormals, VertCurv, FaceCurv, Nring);
  Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert, maxQCrossVert);
  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);

  has_cross_field = true;
}

void SaveMesh() {
  assert(!PathMesh.empty());
  std::string PathRemMesh = RemoveExtension(PathMesh);
  PathRemMesh += std::string("_remesh.obj");
  std::cout << "Saving Field on file " << PathRemMesh.c_str() << std::endl;
  bool Saved = WriteOBJ(PathRemMesh, VertPos, Connectivity);
  if (!Saved) {
    std::cout << "ERROR SAVING FIELD FILE " << std::endl;
    exit(0);
  }
}

void SaveField() {
  assert(!PathMesh.empty());
  std::string PathField = RemoveExtension(PathMesh);
  PathField += std::string("_remesh.field");
  std::cout << "Saving Field on file " << PathField.c_str() << std::endl;
  bool Saved = WriteField(PathField, VertCurv, FaceCurv);
  if (!Saved) {
    std::cout << "ERROR SAVING FIELD FILE " << std::endl;
    exit(0);
  }
  // ReadField(VertPos, Connectivity, PathField, VertCurv, FaceCurv);
}

void SaveFeatures() {
  if (Features.size() == 0)
    return;
  assert(!PathMesh.empty());
  std::string PathFeatures = RemoveExtension(PathMesh);
  PathFeatures += std::string("_remesh.features");
  std::cout << "Saving Features on file " << PathFeatures.c_str() << std::endl;
  bool Saved = WriteFeatures(PathFeatures, Features);
  if (!Saved) {
    std::cout << "ERROR SAVING FIELD FILE " << std::endl;
    exit(0);
  }
}

void SmoothGlobal(bool iterative = false) {
  if (!has_cross_field)
    return;

  std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;
  std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;

  SParam.Features = Features;

  if (!iterative)
    SmoothGlobalPolyvector(VertPos, Connectivity, FaceNormals, FaceCurv,
                           SParam);
  else
    SmoothCrossIterative(VertPos, Connectivity, FaceNormals, FaceCurv, SParam);

  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);
  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);

  RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
  RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
}

void RefineForTracing() {
  Geo::RefineForTracing<ScalarType>::RefineRequiredEdges(VertPos, Connectivity,
                                                         Features);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  has_remeshed = true;
}

void RefineForTracingConnectivity() {
  Geo::RefineForTracing<ScalarType>::RefineForConnection(VertPos, Connectivity,
                                                         FaceCurv, Features);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
}

void BatchProcessCurv() {
  if (processForTracing)
    RefineForTracing();

  //WriteOBJ("testMesh0.obj",VertPos, Connectivity);
  
  SParam.FixedCross.clear();
  InitFieldByCurvature();
  SmoothGlobal();

  if (processForTracing) {
    RefineForTracingConnectivity();
  }
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
    //DesambiguateSamePosVertices(VertPos,Connectivity,0.0001);
    return true;
    // std::cout<<"FNorm 0 "<<FaceNormals.size()<<std::endl;
  } else {
    std::cout << "Mesh not loaded" << std::endl;
    return false;
  }
}

bool LoadField(const std::string &path) {
  // std::cout<<"Loading Field"<<path.c_str()<<std::endl;
  bool has_loaded = ReadField(VertPos, Connectivity, path, VertCurv, FaceCurv);
  if (has_loaded) {
    std::cout << "Loaded Field" << std::endl;
    showCross = true;
    showSing = true;

    Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert,
                                          maxQCrossVert);
    GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE,
                         SingIndex, SingValue);
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

void UpdateFeatureCorners() {
  Geo::FeatureCorners(VertPos, Connectivity, 
        initial_feature_angle, 
    Features,Corners);

  Geo::RemoveSingleFaceFeatures(Connectivity,Features);

  if ((ero_steps>0)||(dil_steps>0))
  {
    Geo::ErodeDilateFeatures(ero_steps,dil_steps,Features) ;
    Geo::Corners<ScalarType>(VertPos, Features, initial_feature_angle, Corners);
  }
  
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

bool lineF = (SParam.NDir == 2);
bool crossF = (SParam.NDir == 4);

// std::vector<std::array<int, 3>> TripleVert;
std::vector<int> DisconnVert;
std::vector<Geo::Point3<ScalarType>> DisconnDir;
void SetRenderBar() {
  ImGui::Begin(
      "Field Smoother", NULL,
      ImGuiWindowFlags_AlwaysAutoResize); // Create a window called "Hello,
                                          // world!" and append into it.

  // if (ImGui::Button("SMOOTH MESH")) {
  //   SmoothMesh();
  // }
  // ImGui::SameLine();
  // if (ImGui::Button("SMOOTH MESHING")) {
  //   SmoothTessellation();
  // }
  ImGui::Separator();
  ImGui::Text("Render");
  ImGui::Checkbox("Show Cross", &showCross);
  ImGui::Checkbox("Show Sing", &showSing);
  ImGui::Checkbox("Process for Tracing", &processForTracing);

  if (ImGui::Button("BATCH PROCESS")) {
    BatchProcessCurv();
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
  ImGui::SameLine();
  if (has_cross_field && loadedMesh) {
    if (ImGui::Button("SAVE ALL")) {
      SaveMesh();
      SaveField();
      SaveFeatures();
    }
  }

  ImGui::Separator();
  ImGui::Text("Mesh");

  if (!has_load_features) {
    int valueA0 =
        floor(0.5 + initial_feature_angle); // Variable to store slider value
    ImGui::SetNextItemWidth(100.0f);
    ImGui::SliderInt("Inital features", &valueA0, 0, 180);
    initial_feature_angle = valueA0;
    

    ImGui::SliderInt("Ero Steps", &ero_steps, 0, 5);
    ImGui::SliderInt("Dil Steps", &dil_steps, 0, 5);

    if ((old_initial_feature_angle != initial_feature_angle)||
        (old_ero_steps != ero_steps)||
        (old_dil_steps != dil_steps)) 
    {
      UpdateFeatureCorners();
      old_initial_feature_angle = initial_feature_angle;
      old_ero_steps = ero_steps;
      old_dil_steps = dil_steps;
    }
    // if (old_ero_dil_steps != ero_dil_steps) {
    //   UpdateFeatureCorners();
    //   old_ero_dil_steps = ero_dil_steps;
    // }
  }

  if (ImGui::Button("Refine For Tracing")) {
    // SplitBorderEdges();
    RefineForTracing();
  }

  if (has_cross_field)
    if (ImGui::Button("Refine for Connectivity")) {
      RefineForTracingConnectivity();
      // Geo::RefineForTracing<ScalarType>::RefineForConnection(
      //     VertPos, Connectivity, FaceCurv, Features);
      // Geo::getFFAdjacency(Connectivity, NextF, NextE);
      // ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
      // Geo::RefineForTracing<ScalarType>::FindNonConnectedVert(VertPos,Connectivity,FaceCurv,Features,
      //                                                         DisconnVert,DisconnDir);
    }

  // ImGui::SameLine();
  // if (ImGui::Button("BATCH PROCESS BORD")) {
  //   BatchProcessBord();
  // }

  ImGui::Separator();
  ImGui::Text("Init");

  ImGui::SliderInt("Ring", &Nring, 2, 5, "%d", ImGuiSliderFlags_AlwaysClamp);
  if (ImGui::Button("Init By Curvature")) {
    if (loadedMesh)
      InitFieldByCurvature();
  }

  ImGui::Separator();
  ImGui::Text("Smooth");
  if (ImGui::RadioButton("Line F", lineF)) {
    lineF = true;
    crossF = false;
    SParam.NDir = 2;
  }
  ImGui::SameLine();
  if (ImGui::RadioButton("Cross F", crossF)) {
    lineF = false;
    crossF = true;
    SParam.NDir = 4;
  }

  ImGui::Checkbox("Fix Border", &SParam.align_borders);
  // ImGui::Checkbox("Show Sing", &showSing);

  // Create a slider for the float value
  ImGui::SliderFloat("Top Fix", &TopIntevalFix, 0.0f, 0.3f);
  SParam.TopIntevalFix = (ScalarType)TopIntevalFix;

  //        ImGui::SliderFloat("Top Fix",&TopIntevalFix, 0, 0.3,
  //        "%2.f",ImGuiSliderFlags_AlwaysClamp); SParam.TopIntevalFix=
  //        (ScalarType) TopIntevalFix;

  ImGui::SliderFloat("Smoothness", &GlobalSmoothVal, 50, 500, "%1.f",
                     ImGuiSliderFlags_AlwaysClamp);
  SParam.GlobalSmoothVal = GlobalSmoothVal;

  if (has_cross_field && loadedMesh) {
    if (ImGui::Button("Smooth Field")) {
      if (loadedMesh)
        SmoothGlobal();
    }
    // ImGui::SameLine();
    ImGui::Separator();
    ImGui::Text("Save");

    if (ImGui::Button("Save Mesh")) {
      SaveMesh();
    }

    ImGui::SameLine();
    if (ImGui::Button("Save Field")) {
      SaveField();
    }

    ImGui::SameLine();
    if (ImGui::Button("Save Features")) {
      SaveFeatures();
    }
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
    GLDraw::DrawSurfaceMesh<ScalarType>(VertPos, Connectivity, FaceNormals,
                                        VertNormals, DrawSurface, true, false,
                                        false);
    GLDraw::GLDrawBorders<ScalarType>(VertPos, Connectivity, NextF, 5);

    GLDraw::GLDrawEdges<ScalarType>(VertPos, Features, 15,
                                    Geo::Point3<ScalarType>(1, 0, 1));

    // std::vector<int> DrawV;
    // for (size_t i = 0; i < TripleVert.size(); i++) {
    //   DrawV.push_back(TripleVert[i][0]);
    //   DrawV.push_back(TripleVert[i][1]);
    //   DrawV.push_back(TripleVert[i][2]);
    // }
    GLDraw::DrawVertices(VertPos, DisconnVert, Geo::Point3<ScalarType>(1, 0, 0),
                         30);
    glPopAttrib();
  }

  if (has_cross_field) {
    ScalarType scaleVal = AvEdge * (ScalarType)0.5;
    if (showCross)
      GLDraw::DrawCrossFields<ScalarType>(FaceCurv, scaleVal, maxQCrossVert,
                                          minQCrossVert);

    if (showSing)
      GLDraw::DrawSingularity<ScalarType>(VertPos, SingIndex, SingValue);
  }
}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cout << "You should pass at least a mesh" << std::endl;
    exit(0);
  }

  bool has_load_mesh = false;
  for (int i = 1; i < argc; i++) {
    if (GetFileExtension(std::string(argv[i])) == std::string("obj")) {
      PathMesh = std::string(argv[i]);
      has_load_mesh = LoadMesh(PathMesh);
      if (!has_load_mesh)
        exit(0);
      continue;
    }

    if (GetFileExtension(std::string(argv[i])) == std::string("field")) {
      PathField = std::string(argv[i]);
      has_cross_field = LoadField(PathField);
      if (!has_cross_field)
        exit(0);
      continue;
    }

    if (GetFileExtension(std::string(argv[i])) == std::string("txt")) {
      std::string PathClist = std::string(argv[i]);
      bool load_clist = LoadClist(PathClist, Features, true);
      if (!load_clist)
        exit(0);
      std::cout << "Loaded Features" << std::endl;
      has_load_features = true;
      continue;
    }

    if (GetFileExtension(std::string(argv[i])) == std::string("features")) {
      std::string PathFeat = std::string(argv[i]);
      bool load_feat = ReadFeatures(PathFeat, Features, true);
      if (!load_feat)
        exit(0);
      std::cout << "Loaded Features" << std::endl;
      has_load_features = true;
      continue;
    }
  }

  // glWrap.m=&mesh;

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
