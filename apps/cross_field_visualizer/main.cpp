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

#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <stdio.h>
#include <UI/orbitball.h>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include <GLFW/glfw3.h>
#include "IO/io_obj.h"
#include <point3.h>
#include <box3.h>
#include <GL_functions/GLDrawMesh.h>
#include <GL_functions/GLDrawField.h>
#include "IO/io_field.h"
#include <connectivity.h>
#include <Space_Query/triangle_grid_2D.h>

GLFWwindow* window;
std::string PathMesh;
std::string PathField;
UI::OrbitBall trackball;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
ImVec4 mesh_color = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
float line_thick=1.f;
int selected_rend_item = 5;

bool showMesh=true;
bool showCross=true;
bool showSing=true;


#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

typedef double ScalarType;

//MESH DATA
bool loadedMesh=false;
std::vector<Geo::Point3<ScalarType > > VertPos;
std::vector<std::vector<int> > Connectivity;
std::vector<Geo::Point3<ScalarType> > FaceNormals;
std::vector<Geo::Point3<ScalarType> > VertNormals;
std::vector<std::vector<int> > NextF,NextE;
Geo::Box3<ScalarType> MeshBox;
Geo::Box2<ScalarType> UvBox;

//FIELD DATA
std::vector<Field::CrossF<ScalarType> > VertCurv;
std::vector<Field::CrossF<ScalarType> > FaceCurv;
bool has_cross_field=false;
int Nring=3;
ScalarType minQCrossVert,maxQCrossVert;
std::vector<int> SingIndex,SingValue;

void TestSample()
{
    Geo::TriangleSpatialIndex2D<ScalarType> Grid_2D;

    //get from 3D to 2D position
    std::vector<Geo::Point2<ScalarType > > VertPos_2D;
    for (size_t i=0;i<VertPos.size();i++)
        VertPos_2D.push_back(Geo::Point2<ScalarType >(VertPos[i].X,VertPos[i].Y));

    //initialize the grid
    Grid_2D.Init(VertPos_2D,Connectivity);

    //then do some random search
    int NumTest=100;
    UvBox.Init(VertPos_2D);
    int interval_sample=10000;
    for (int i=0;i<NumTest;i++)
    {
      //get random x and y
      ScalarType XRand=(ScalarType)(rand()%interval_sample )/(ScalarType)interval_sample;
      ScalarType YRand=(ScalarType)(rand()%interval_sample )/(ScalarType)interval_sample;
      XRand*=UvBox.DimX();
      YRand*=UvBox.DimY();
      XRand+=UvBox.Min.X;
      YRand+=UvBox.Min.Y;

      //test point
      Geo::Point2<ScalarType > TestP(XRand,YRand);

      //returning data
      int IndexF;
      Geo::Point3<ScalarType> bary;
      //then search wich face
      Grid_2D.FirstTrianglesOverPoint(VertPos_2D,Connectivity,TestP,IndexF,bary);
    }
}

void LoadMesh(const std::string &path)
{
    //    oss<<"Loading "<<path.c_str()<<std::endl;
    loadedMesh=LoadOBJ(path,VertPos,Connectivity,true);
    if (loadedMesh)
    {
        std::cout<<"Loaded "<<Connectivity.size()<<" faces "<<VertPos.size()<<" vertices"<<std::endl;
        MeshBox.Init(VertPos);
        ComputeNormals(VertPos,Connectivity,FaceNormals,VertNormals);
        Geo::getFFAdjacency(Connectivity,NextF,NextE);

        //std::cout<<"FNorm 0 "<<FaceNormals.size()<<std::endl;
    }
    else
        std::cout<<"Mesh not loaded"<<std::endl;
}

void LoadField(const std::string &path)
{
    //std::cout<<"Loading Field"<<path.c_str()<<std::endl;
    has_cross_field=ReadField(VertPos,Connectivity,path,VertCurv,FaceCurv);
    if (VertCurv.size()==0)
    {
        SetVertCrossFromFace(VertPos,Connectivity,VertNormals,FaceCurv,VertCurv);
        InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);
    }
    if (FaceCurv.size()==0)
    {
        SetFaceCrossVectorFromVert(Connectivity,FaceNormals,VertCurv,FaceCurv);
        InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);
    }
    if (has_cross_field)
    {
        std::cout<<"Loaded Field"<<std::endl;
        showCross=true;
        showSing=true;

        Field::CrossF<ScalarType>::getMinMaxQ(VertCurv,minQCrossVert,maxQCrossVert);
        GetFaceSingularities(FaceCurv,VertPos,Connectivity,NextF,NextE,SingIndex,SingValue);
    }
    else
        std::cout<<"Field not loaded"<<std::endl;
}

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if(ImGui::GetIO().WantCaptureMouse) return;
    
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    float xscale, yscale;
    glfwGetWindowContentScale(window, &xscale, &yscale);
    xpos*=xscale;
    ypos*=yscale;
    ypos=(display_h)-ypos;
    trackball.MouseMove((int)xpos,(int)ypos);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    (void)mods;

    if(ImGui::GetIO().WantCaptureMouse) return;

    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    float xscale, yscale;
    glfwGetWindowContentScale(window, &xscale, &yscale);
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    xpos*=xscale;
    ypos*=yscale;
    ypos=(display_h)-ypos;

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
        trackball.MouseDown((int)xpos,(int)ypos,UI::OrbitBall::BUTTON_LEFT);

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
        trackball.MouseUp((int)xpos,(int)ypos,UI::OrbitBall::BUTTON_LEFT);

}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    (void)window;
    (void)xoffset;
    trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    (void)mods;
    (void)window;
    (void)scancode;

    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
        trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);
    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
        trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);
}

void InitGLFW_Window()
{
    if (!glfwInit())
        exit(0);

    window = glfwCreateWindow(1280, 720, "Field Visualizer", NULL, NULL);

    if (window == NULL)
        exit(0);


    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glfwSetErrorCallback(glfw_error_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetKeyCallback(window, key_callback);
}

void InitIMGui()
{
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;

    // Setup Dear ImGui style
    ImGui::StyleColorsClassic();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();
}

void SetRenderBar()
{
    ImGui::Begin("Field Visualizer",NULL,ImGuiWindowFlags_AlwaysAutoResize);       // Create a window called "Hello, world!" and append into it.

    ImGui::Checkbox("Draw Mesh", &showMesh);
    ImGui::Checkbox("Draw Cross", &showCross);
    ImGui::Checkbox("Draw Singularities", &showSing);

    if (ImGui::Button("Test Sample"))
    {
        TestSample();
    }

    ImGui::End();
}

void GLDrawMesh()
{
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(40, (GLdouble)display_w/(GLdouble)display_h, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);

    trackball.GetView();

    trackball.Apply();
    //if (showTrackBall)
    //    trackball.DrawPostApply();

    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    GLDraw::glScale(3/MeshBox.Diag());
    GLDraw::glTranslate(-MeshBox.Center());

    bool DrawSurface=true;

    //mesh.GLDrawBox();
    if (showMesh)
    {
        glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w, mesh_color.z * mesh_color.w, mesh_color.w);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(line_thick);
        GLDraw::DrawSurfaceMesh<ScalarType>(VertPos,Connectivity,FaceNormals,VertNormals,DrawSurface,true,true,false);
        glPopAttrib();

    }

    ScalarType scaleVal=MeshBox.Diag()*(ScalarType)0.02;
    if (showCross)
        GLDraw::DrawCrossFields<ScalarType>(FaceCurv,scaleVal,maxQCrossVert,minQCrossVert);

    if (showSing)
        GLDraw::DrawSingularity<ScalarType>(VertPos,SingIndex,SingValue);
}

int main(int argc, char *argv[])
{

    //oss<<"Initialized Log Window"<<std::endl;
    //AScene.SetOutputStream(&oss);

    if (argc>=2){
        PathMesh=std::string(argv[1]);
        LoadMesh(PathMesh);

    }

    if (argc>=3){
        PathField=std::string(argv[2]);
        LoadField(PathField);
    }

    //glWrap.m=&mesh;
    
    InitGLFW_Window();

    InitIMGui();

    trackball.center=Geo::Point3<float>(0, 0, 0);
    trackball.radius= 1;
    //trackball.SetSpinnable(true);

    glewInit();
    // Main loop
    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        SetRenderBar();

        // Rendering
        ImGui::Render();

        GLDrawMesh();

        // If you are using this code with non-legacy OpenGL header/contexts (which you should not, prefer using imgui_impl_opengl3.cpp!!),
        // you may need to backup/reset/restore other state, e.g. for current shader using the commented lines below.
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
