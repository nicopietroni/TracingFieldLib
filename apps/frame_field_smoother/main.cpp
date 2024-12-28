#include "ImGuiFileBrowser.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <GL/glew.h>
#include <UI/orbitball.h>
#include <iostream>
#include <stdio.h>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>
#include <GL_functions/GLDrawField.h>
#include <GL_functions/GLDrawMesh.h>
#include <IO/file_path_functions.h>
#include <IO/io_field.h>
#include <IO/io_obj.h>
#include <IO/io_poly.h>
#include <box3.h>
#include <point3.h>

#include <color_functions.h>
#include <frame_field.h>
#include <frame_field_singularities.h>
#include <polyhedra_mesh_functions.h>
#include <sample_for_vis_frame_field.h>
#include <smooth_cross_field.h>
#include <smooth_frame_field.h>

GLFWwindow *window;
std::string PathMesh;
std::string PathField;
UI::OrbitBall trackball;
imgui_addons::ImGuiFileBrowser file_dialog;
// std::ostringstream oss;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);

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

typedef double UsedScalarType;

std::vector<Geo::Point3<UsedScalarType>> VertPos;
std::vector<Geo::Point3<UsedScalarType>> FaceNormals;
std::vector<Geo::Point3<UsedScalarType>> VertNormals;

Geo::Box3<UsedScalarType> MeshBox;
std::vector<std::vector<int>> PolyConnectivity;

std::vector<std::vector<int>> NextPoly;
std::vector<std::vector<int>> NextPolyF;

// std::vector<std::vector<int> > NextPoly;
// std::vector<std::vector<int> > NextPolyF;
////std::vector<std::vector<std::pair<int,int> > > PPAdjacency;
////std::vector<std::vector<std::pair<int,int> > > VPAdjacency;

std::vector<std::vector<int>> ExtFace;
std::vector<int> ExtPoly;
std::vector<int> ExtPolyFaceIndex;

std::vector<Field::FrameField<UsedScalarType>> FrameField;

// std::string TetraMeshPath;
// std::string FieldPath;

bool has_frame_field = false;
bool has_cross_field = false;

// TracingGraph<UsedScalarType> TGraph;

// UsedScalarType drift_penalty=100;
// UsedScalarType maxAngle=45;
// size_t diffuse_step=0;

// size_t TestN=0;
// std::vector<int> PathN;

std::vector<Field::CrossF<UsedScalarType>> VertCurv;
std::vector<Field::CrossF<UsedScalarType>> FaceCurv;

UsedScalarType maxQCrossVert = 0;
UsedScalarType minQCrossVert = 0;

// external surface maps
std::vector<std::vector<int>> SurfExtFace;
std::vector<Geo::Point3<UsedScalarType>> SurfVertPos;
std::vector<Geo::Point3<UsedScalarType>> SurfFaceNormals, SurfVertNormals;
std::vector<int> SingIndex, SingValue;

std::vector<std::vector<int>> NextF, NextE;
int Nring = 3;

int anglePercThr = 5;
SmoothParam<UsedScalarType> SParam;
int GlobSmoothI = (int)SParam.GlobalSmoothVal;

// SmoothFrameParam<UsedScalarType> SParamFrame;

////sampled pos to test
////std::vector<Point3<UsedScalarType> > SampledPos;

bool DrawSurface = true;
bool DrawWire = false;
bool TransparentSurface = true;

bool DrawCross = true;
bool DrawCrossSing = true;

bool DrawFrame = true;
bool DrawFrameSing = true;

std::vector<Pos4> FrameSingPos;
std::vector<int> FrameSingValue;
std::vector<std::vector<Geo::Point3<UsedScalarType>>> PolyLSing;
std::vector<Geo::Point3<UsedScalarType>> PolyColorSing;

std::vector<std::vector<Geo::Point3<UsedScalarType>>> PolyL;
std::vector<Geo::Point3<UsedScalarType>> PolyColor;

SmoothFrameParam<UsedScalarType> SParamFrame;

// void UpdateConstraints()
//{
//     InitLinearConstr();
//     UpdateArap();
// }

UsedScalarType initial_feature_angle = 45;
UsedScalarType old_initial_feature_angle = 0;


void LoadMesh(const std::string &path) {

  LoadTetraMesh<UsedScalarType>(path.c_str(), VertPos, PolyConnectivity);

  MeshBox.Init(VertPos);

  Tetra::GetPPAdjacency(PolyConnectivity, NextPoly, NextPolyF);

  // extract only external surface
  Tetra::FFExternalMesh(PolyConnectivity, NextPoly, NextPolyF, ExtFace, ExtPoly,
                        ExtPolyFaceIndex);
  ComputeNormals(VertPos, ExtFace, FaceNormals, VertNormals);

  SurfExtFace = ExtFace;
  SurfVertPos = VertPos;

  RemoveUnfererencedVert<UsedScalarType>(SurfVertPos, SurfExtFace);
  ComputeNormals(SurfVertPos, SurfExtFace, SurfFaceNormals, SurfVertNormals);

  // update adjacency
  Geo::getFFAdjacency(SurfExtFace, NextF, NextE);

  //    if (!FieldPath.empty())
  //    {

  //        has_frame_field=true;

  //        std::vector<Geo::Point3<UsedScalarType> > Centers;
  //        Tetra::FindCenter(VertPos,PolyConnectivity,Centers);

  //        Field::LoadFrameFields(FieldPath.c_str(),Centers,FrameField);

  //        UpdateFrameFieldVisualization();

  //        UpdateCrossFieldFromFrame();
  //    }
}

void InitSurfField() {
  ComputeCurvatureField(SurfVertPos, SurfExtFace, SurfFaceNormals,
                        SurfVertNormals, VertCurv, FaceCurv, Nring);
  Field::CrossF<UsedScalarType>::getMinMaxQ(VertCurv, minQCrossVert,
                                            maxQCrossVert);
  GetFaceSingularities(FaceCurv, VertPos, SurfExtFace, NextF, NextE, SingIndex,
                       SingValue);
  has_cross_field = true;
}

void PropagateCrossFromFix() {
  if (!has_cross_field)
    return;
  PropagateCrossFromFixed(SurfVertPos, SurfExtFace, SurfFaceNormals, FaceCurv,
                          SParam);
  GetFaceSingularities(FaceCurv, VertPos, SurfExtFace, NextF, NextE, SingIndex,
                       SingValue);
}

void UpdateCrossFieldFromFrame() {
  CrossFieldFromFrame(FrameField, ExtPoly, SurfExtFace, SurfFaceNormals,
                      SurfVertPos, FaceCurv);
  SetVertCrossFromFace<UsedScalarType>(SurfVertPos, SurfExtFace,
                                       SurfVertNormals, FaceCurv, VertCurv);
  GetFaceSingularities(FaceCurv, VertPos, SurfExtFace, NextF, NextE, SingIndex,
                       SingValue);
  has_cross_field = true;
}

void SmoothCrossIte() {
  if (!has_cross_field)
    return;
  SmoothCrossIterative(SurfVertPos, SurfExtFace, SurfFaceNormals, FaceCurv,
                       SParam); //,VertCurv);
  GetFaceSingularities(FaceCurv, VertPos, SurfExtFace, NextF, NextE, SingIndex,
                       SingValue);
}

void SmoothCrossGlobal() {
  if (!has_cross_field)
    return;
  SmoothGlobalPolyvector(SurfVertPos, SurfExtFace, SurfFaceNormals, FaceCurv,
                         SParam);
  GetFaceSingularities(FaceCurv, VertPos, SurfExtFace, NextF, NextE, SingIndex,
                       SingValue);
}

void SaveCrossField() {
  assert(!PathMesh.empty());
  std::string PathCrossField;
  if (PathField.empty()) {
    PathCrossField = RemoveExtension(PathMesh);
    PathCrossField += std::string("_cross.field");
  } else {
    PathCrossField = RemoveExtension(PathField);
    PathCrossField += std::string("_cross.field");
  }
  std::cout << "Saving Field on file " << PathCrossField.c_str() << std::endl;
  bool Saved = WriteField(PathCrossField, VertCurv, FaceCurv);
  if (!Saved) {
    std::cout << "ERROR SAVING FIELD FILE " << std::endl;
    exit(0);
  }
}

void SaveFrameField() {
  assert(!PathMesh.empty());
  std::string PathFrameField;
  if (PathField.empty()) {
    PathFrameField = RemoveExtension(PathMesh);
    PathFrameField += std::string("_frame.field");
  } else {
    PathFrameField = RemoveExtension(PathField);
    PathFrameField += std::string("_frame.field");
  }
  std::cout << "Saving Field on file " << PathFrameField.c_str() << std::endl;

  std::vector<Field::FrameField<UsedScalarType>> VertField;
  bool Saved = WriteField(PathFrameField, VertField, FrameField);
  if (!Saved) {
    std::cout << "ERROR SAVING FIELD FILE " << std::endl;
    exit(0);
  }
}

void UpdateFrameFieldVisualization() {
  if (!has_frame_field)
    return;

  GetFrameSingularities<UsedScalarType>(FrameField, PolyConnectivity, NextPoly,
                                        NextPolyF, FrameSingPos,
                                        FrameSingValue);
  GetPolylineFromSing(VertPos, PolyConnectivity, FrameSingPos, FrameSingValue,
                      PolyLSing, PolyColorSing);

  Field::SampleCrossPolylines<UsedScalarType>(VertPos, PolyConnectivity,
                                              FrameField, PolyL, 200);
  PolyColor.resize(PolyL.size());
  for (size_t i = 0; i < PolyColor.size(); i++)
    PolyColor[i] = Geo::ScatterColor<UsedScalarType>(PolyColor.size(), i);
  // PolyColor[i].ScatterColor(PolyColor.size(),i);
}

void PropagateFrameFromBorder() {
  if (!has_cross_field)
    return;
  if (!Field::IsProperCrossField(FaceCurv))
    return;

  // std::cout<<"Initializing Graph.."<<std::endl;
  // SmoothFrameParam<UsedScalarType> SParam;
  PropagateFrameFromBorderCross(SurfExtFace, FaceCurv, ExtPoly, VertPos,
                                PolyConnectivity, FrameField, SParamFrame);

  has_frame_field = true;
  UpdateFrameFieldVisualization();
}

void SmoothFrameIterative() {
  if (!has_frame_field)
    return;
  // int t0=clock();
  SmoothFrameIterative(VertPos, PolyConnectivity,
                       FrameField,SParamFrame); //,VertCurv);
  // int t1=clock();
  UpdateFrameFieldVisualization();

  UpdateCrossFieldFromFrame();
}

void BatchProcess() {
  InitSurfField();
  SmoothCrossGlobal();
  PropagateFrameFromBorder();
  SmoothFrameIterative();
}

void LoadField(const std::string &path) {
  has_frame_field = true;

  std::vector<Geo::Point3<UsedScalarType>> Centers;
  Tetra::FindCenter(VertPos, PolyConnectivity, Centers);

  std::string ExtField = GetFileExtension(path);
  bool LoadedField = false;
  if (ExtField == std::string("txt"))
    LoadedField = LoadFrameFieldTxt(path.c_str(), Centers, FrameField);

  if (ExtField == std::string("field")) {
    std::vector<Field::FrameField<UsedScalarType>> VertField;
    LoadedField =
        ReadField(VertPos, PolyConnectivity, path, VertField, FrameField);
  }

  if (!LoadedField) {
    std::cout << "ERROR LOADING FIELD FILE " << std::endl;
    exit(0);
  }

  UpdateFrameFieldVisualization();

  UpdateCrossFieldFromFrame();
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

  //    xPick=xpos;
  //    yPick=ypos;

  //    if (CMode!=CMNone)return;

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

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    trackball.MouseUp((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
    // CMode=CMNone;
    // PickedV=-1;
  }

  //    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS &&
  //    (CMode!=CMNone))
  //    {
  //        xPick=xpos;
  //        yPick=ypos;
  //        to_pick=true;
  //    }
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

// Now inside any function
void showMainMenu() {
  bool open = false, save = false;
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem("Open", NULL))
        open = true;
      if (ImGui::MenuItem("Quit", NULL))
        exit(0);

      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  // Remember the name to ImGui::OpenPopup() and showFileDialog() must be
  // same...
  if (open)
    ImGui::OpenPopup("Open File");
  // if(save)
  //     ImGui::OpenPopup("Save File");

  /* Optional third parameter. Support opening only compressed rar/zip files.
   * Opening any other file will show error, return false and won't close the
   * dialog.
   */
  if (file_dialog.showFileDialog(
          "Open File", imgui_addons::ImGuiFileBrowser::DialogMode::OPEN,
          ImVec2(700, 310), ".obj,.ply,.off")) {
    LoadMesh(file_dialog.selected_path);
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

void SetRenderBar() {
  ImGui::Begin(
      "Field Computation", NULL,
      ImGuiWindowFlags_AlwaysAutoResize); // Create a window called "Hello,
                                          // world!" and append into it.

  ImGui::Separator();
  ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW SURFACE");
  ImGui::Checkbox("Surface", &DrawSurface);
  ImGui::SameLine();
  ImGui::Checkbox("Wire", &DrawWire);
  ImGui::SameLine();
  ImGui::Checkbox("Trasparent", &TransparentSurface);

  ImGui::Separator();
  ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW CROSS");
  ImGui::Checkbox("Cross", &DrawCross);
  ImGui::SameLine();
  ImGui::Checkbox("Cross Sing", &DrawCrossSing);
  ImGui::Separator();

  ImGui::Separator();
  ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW FRAME");
  ImGui::Checkbox("Frame", &DrawFrame);
  ImGui::SameLine();
  ImGui::Checkbox("Frame Sing", &DrawFrameSing);
  ImGui::Separator();

  ImGui::TextColored(ImVec4(1, 0, 1, 1), "COMPUTE CROSS");

   int valueA0 =
        floor(0.5 + initial_feature_angle); // Variable to store slider value
    ImGui::SetNextItemWidth(100.0f);
    ImGui::SliderInt("Inital features", &valueA0, 0, 180);
    initial_feature_angle = valueA0;
    if (old_initial_feature_angle != initial_feature_angle) {
      //UpdateFeatureCorners();
      old_initial_feature_angle = initial_feature_angle;
    }

  ImGui::PushItemWidth(100);
  ImGui::SliderInt("NRing", &Nring, 0,
                   5); //, "%2d",ImGuiSliderFlags_AlwaysClamp);
  ImGui::PopItemWidth();

  ImGui::SameLine();
  if (ImGui::Button("Init Curvature")) {
    InitSurfField();
  }

  // ImGui::InputDouble("Hard Curvature", &SParam.TopIntevalFix, 0, 0.005);
  //&SParam.TopIntevalFix
  if (has_cross_field) {
    ImGui::PushItemWidth(100);
    ImGui::SliderInt("% Fixed", &anglePercThr, 0,
                     10); //, "%2d",ImGuiSliderFlags_AlwaysClamp);
    SParam.TopIntevalFix = (UsedScalarType)1 / (UsedScalarType)anglePercThr;
    ImGui::PopItemWidth();

    ImGui::SameLine();

    if (ImGui::Button("Propagate Fix")) {
      PropagateCrossFromFix();
    }
    // ImGui::SameLine();
    // ImGui::SameLine();
    //ImGui::Checkbox("Align Border", &SParam.align_borders);

    //    TwAddButton(bar,"Propagate
    //    Fix",PropagateCrossFromFix,0,"label='Propagate From Fix'");
    ImGui::PushItemWidth(100);
    ImGui::SliderInt("Smooth Ite", &SParam.iterations, 10,
                     100); //, "%2d",ImGuiSliderFlags_AlwaysClamp);
    ImGui::PopItemWidth();

    ImGui::SameLine();
    if (ImGui::Button("Smooth Cross Ite")) {
      SmoothCrossIte();
    }

    ImGui::PushItemWidth(100);
    ImGui::SliderInt("Smooth Val", &GlobSmoothI, 10, 1000);
    ImGui::PopItemWidth();
    SParam.GlobalSmoothVal = (double)GlobSmoothI;
    ImGui::SameLine();
    if (ImGui::Button("Smooth Global")) {
      SmoothCrossGlobal();
    }

    ImGui::Separator();
    ImGui::TextColored(ImVec4(1, 0, 1, 1), "COMPUTE FRAME");
    if (ImGui::Button("Propagate Frame")) {
      PropagateFrameFromBorder();
    }
    // ImGui::SameLine();

    if (has_frame_field) {
      ImGui::Checkbox("Weight Q", &SParamFrame.weightByQuality);
      ImGui::SameLine();
      ImGui::Checkbox("Fix Border", &SParamFrame.fix_borders);
      ImGui::SameLine();
      ImGui::Checkbox("Align Border", &SParamFrame.align_borders);

      if (ImGui::Button("Smooth Frame Ite")) {
        SmoothFrameIterative();
      }
    }

    ImGui::Separator();
    ImGui::TextColored(ImVec4(1, 0, 1, 1), "BATCH PROCESS");
    if (ImGui::Button("Batch Process")) {
      BatchProcess();
    }
  }
  ImGui::Separator();
  ImGui::TextColored(ImVec4(1, 0, 1, 1), "SAVE");

  if (ImGui::Button("Save Cross")) {
    SaveCrossField();
  }
  if (has_frame_field) {
    ImGui::SameLine();
    if (ImGui::Button("Save Frame")) {
      SaveFrameField();
    }
  }

  ImGui::End();
}

void GLDrawMesh() {
  // glPushAttrib(GL_ALL_ATTRIB_BITS);
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
  trackball.DrawPostApply();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());

  UsedScalarType scaleVal = MeshBox.Diag() * (UsedScalarType)0.02;

  if (DrawCross) {
    GLDraw::DrawCrossFields<UsedScalarType>(FaceCurv, scaleVal, maxQCrossVert,
                                            minQCrossVert);
  }

  if (DrawCrossSing) {
    GLDraw::DrawSingularity<UsedScalarType>(SurfVertPos, SingIndex, SingValue);
  }

  if (DrawFrameSing && has_frame_field)
    GLDraw::GLDrawPolylines(PolyLSing, PolyColorSing);

  if (DrawFrame && has_frame_field)
    GLDraw::GLDrawPolylines(PolyL, PolyColor);

  glColor4f(.5f, 1.f, .8f, .3f);
  GLDraw::DrawSurfaceMesh<UsedScalarType>(VertPos, ExtFace, FaceNormals,
                                          VertNormals, DrawSurface, DrawWire,
                                          true, TransparentSurface);
}

int main(int argc, char *argv[]) {

  // oss<<"Initialized Log Window"<<std::endl;
  // AScene.SetOutputStream(&oss);

  if (argc >= 2) {
    PathMesh = std::string(argv[1]);
    LoadMesh(PathMesh);
  }

  if (argc >= 3) {
    PathField = std::string(argv[2]);
    LoadField(PathField);
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
    {
      SetRenderBar();

      // SetLogWindow();
      showMainMenu();
    }
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
