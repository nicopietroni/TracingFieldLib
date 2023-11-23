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
#include <GL_functions/GLDrawGraph.h>

#include "IO/io_field.h"
#include <connectivity.h>
#include <field_graph/graph.h>
#include <field_graph/initialize_graph.h>
#include <field_graph/graph_tracing.h>

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
bool has_graph=false;
// [Win32] Our example includes a copy of glfw3.lib pre-compiled with VS2010 to maximize ease of testing and compatibility with old VS compilers.
// To link with VS2010-era libraries, VS2015+ requires linking with legacy_stdio_definitions.lib, which we do using this pragma.
// Your own project should not be affected, as you are likely to link with a newer binary of GLFW that is adequate for your version of Visual Studio.
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

//FIELD DATA
std::vector<Field::CrossF<ScalarType> > VertCurv;
std::vector<Field::CrossF<ScalarType> > FaceCurv;
bool has_cross_field=false;
ScalarType minQCrossVert,maxQCrossVert;
std::vector<int> SingIndex,SingValue;

Geo::TracingGraph<ScalarType> TGraph(Geo::GT_CrossField);

float drift_penalty=100.f;
float maxAngle=45.f;
int diffuse_step=1;

int CurrN=0;
int OldN=-1;
std::vector<int> PathN;

void InitGraph()
{
    Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(TGraph,VertPos,Connectivity,VertCurv,
                                                                   drift_penalty,maxAngle,diffuse_step);
    TGraph.DeSelectAll(0);
    TGraph.SelectNoteType(Geo::ExitDomain,0);
    TGraph.SelectNoteType(Geo::TangentDomain,0);
}

void TestTrace()
{
    int currI=CurrN;
    do{
        if (TGraph.NodeType(currI)==Geo::Internal)
        {
           CurrN=currI;
           break;
        }
        currI=(currI+1)%TGraph.NumNodes();
    }while (true);

    //select the same for loop
    TGraph.DeSelectAll(1);
    TGraph.SelectNoteType(Geo::ExitDomain,1);
    TGraph.Select(CurrN,1);

    Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    SParam.avoid_zero_path=true;
    SParam.stop_at_S0=false;
    SParam.stop_at_S1=true;

    ScalarType PathL;
    std::vector<int> SourceNodes;
    SourceNodes.push_back(CurrN);
    bool traced=Geo::GraphQuery<ScalarType>::ShortestPaths(TGraph,SourceNodes,PathN,PathL,SParam);
//    if (traced)
//        std::cout<<"Traced Correctly"<<std::endl;
//    else
//        std::cout<<"Not Traced"<<std::endl;
    //Geo::GraphQuery<ScalarType>::TraceFrom(TGraph,CurrN,PathN,PathL,SParam);
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
    (void)window;
    (void)scancode;
    (void)mods;

    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
        trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);
    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
        trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);
}


void InitGLFW_Window()
{
    if (!glfwInit())
    {
        //std::cout<<"Error Initializing GLFW"<<std::endl;
        exit(0);
    }
    //return 1;

    window = glfwCreateWindow(1280, 720, "Dear ImGui GLFW+OpenGL2 example", NULL, NULL);
    // if (window == NULL)
    //     return 1;
    if (window == NULL)
    {
        //std::cout<<"Error Initializing GLFW Window"<<std::endl;
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

void InitIMGui()
{
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsClassic();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();
}

void SetRenderBar()
{
    ImGui::Begin("Patch Decomposer",NULL,ImGuiWindowFlags_AlwaysAutoResize);       // Create a window called "Hello, world!" and append into it.

    ImGui::Separator();
    ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW");
    ImGui::Checkbox("Draw Mesh", &showMesh);
    ImGui::Checkbox("Draw Cross", &showCross);
    ImGui::Checkbox("Draw Singularities", &showSing);
    ImGui::Separator();
    ImGui::TextColored(ImVec4(1, 0, 1, 1), "GRAPH");
    ImGui::SliderFloat("Drift",&drift_penalty, 0, 1000);
    ImGui::SliderFloat("Max Angle",&maxAngle, 0, 90);
    ImGui::InputInt("Diffuse Step",&diffuse_step);
    if (ImGui::Button("Init Graph"))
    {
        InitGraph();
        has_graph=true;
    }

    ImGui::InputInt("test node",&CurrN);
    if ((CurrN!=OldN)&&(has_graph))
    {
        TestTrace();
        OldN=CurrN;
    }
//    if (ImGui::Button("Test Trace"))
//    {
//        TestTrace();
//    }

    //    ImGui::SliderInt("Ring",&Nring, 2, 5, "%d",ImGuiSliderFlags_AlwaysClamp);
    //    if (ImGui::Button("Init By Curvature"))
    //    {
    //        if (loadedMesh)
    //            InitFieldByCurvature();
    //    }

    //    if (has_cross_field)
    //    {
    //        if (ImGui::Button("Smooth Field"))
    //        {
    //            if (loadedMesh)
    //                SmoothGlobal();
    //        }
    //        //ImGui::SameLine();

    //        ImGui::SliderFloat("Top Fix",&TopIntevalFix, 0, 0.3, "%2.f",ImGuiSliderFlags_AlwaysClamp);
    //        SParam.TopIntevalFix= (ScalarType) TopIntevalFix;

    //        ImGui::SliderFloat("Smoothness",&GlobalSmoothVal, 50, 500, "%1.f",ImGuiSliderFlags_AlwaysClamp);
    //        SParam.GlobalSmoothVal=GlobalSmoothVal;

    //        if (ImGui::Button("Save Field"))
    //        {
    //            if (loadedMesh)
    //                SaveField();
    //        }
    //    }
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

    if (showMesh)
    {
        glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w, mesh_color.z * mesh_color.w, mesh_color.w);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(line_thick);
        GLDraw::DrawSurfaceMesh<ScalarType>(VertPos,Connectivity,FaceNormals,VertNormals,true,true,true,false);
        glPopAttrib();
    }

    if (has_cross_field)
    {
        ScalarType scaleVal=MeshBox.Diag()*(ScalarType)0.02;
        if (showCross)
            //GLDraw::DrawCrossFields<ScalarType>(VertCurv,scaleVal,maxQCrossVert,minQCrossVert);
            GLDraw::DrawCrossFields<ScalarType>(FaceCurv,scaleVal,maxQCrossVert,minQCrossVert);

        if (showSing)
            GLDraw::DrawSingularity<ScalarType>(VertPos,SingIndex,SingValue);
    }


    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodes(TGraph,VertCurv,PathN,MeshBox.Diag()*0.005);
    //GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawSelectedNodes(TGraph,VertCurv,0,MeshBox.Diag()*0.005);

}

int main(int argc, char *argv[])
{

    //oss<<"Initialized Log Window"<<std::endl;
    //AScene.SetOutputStream(&oss);

    assert(argc>=3);
    if (argc<3)
    {
        std::cout<<"ERROR need mesh and field file"<<std::endl;
        exit(0);
    }

    PathMesh=std::string(argv[1]);
    LoadMesh(PathMesh);

    PathField=std::string(argv[2]);
    LoadField(PathField);
    
    InitGLFW_Window();

    InitIMGui();

    trackball.center=Geo::Point3<float>(0, 0, 0);
    trackball.radius= 1;

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
