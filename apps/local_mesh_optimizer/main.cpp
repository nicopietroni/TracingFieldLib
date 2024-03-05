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
#include <IO/file_path_functions.h>
#include <color_functions.h>
#include <local_operations/edge_flip.h>
#include <local_operations/iterative_optimizer.h>
#include <connectivity.h>

typedef double ScalarType;

GLFWwindow* window;
std::string PathMesh;
UI::OrbitBall trackball;//,trackballUV;

//std::ostringstream oss;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
Geo::Point3<ScalarType> mesh_color(0.79, 0.9, 0.87);

float line_thick=1.f;
int selected_rend_item = 5;
bool showMesh=true;

// [Win32] Our example includes a copy of glfw3.lib pre-compiled with VS2010 to maximize ease of testing and compatibility with old VS compilers.
// To link with VS2010-era libraries, VS2015+ requires linking with legacy_stdio_definitions.lib, which we do using this pragma.
// Your own project should not be affected, as you are likely to link with a newer binary of GLFW that is adequate for your version of Visual Studio.
#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif


//MESH DATA
bool loadedMesh=false;

//positions in 3D and UV
std::vector<Geo::Point3<ScalarType > > VertPos;
//std::vector<Geo::Point2<ScalarType> >  VerPosUV;

//Connectivity of the mesh and Normals
std::vector<std::vector<int> > Connectivity;
std::vector<Geo::Point3<ScalarType> > FaceNormals;
std::vector<Geo::Point3<ScalarType> > VertNormals;

//Bounding Box in 3D
Geo::Box3<ScalarType> MeshBox;


Geo::IterativeOptimizer<ScalarType> IteOpt(VertPos,Connectivity);

void UpdateColorMesh()
{
    //    // Do something with the selected item (e.g., print it)
    //    switch(currentColorItem)
    //    {

    //    case 0:
    //        FaceColor=std::vector<Geo::Point3<ScalarType> > (Connectivity.size(),mesh_color);
    //        break;

    //    case 1:
    //        Geo::UVDistortion<ScalarType>::ComputeDistortion(VertPos,VerPosUV,Connectivity,
    //                                                         Geo::UVDistortion<ScalarType>::AngleDist,
    //                                                         FaceDistortion);
    //        Geo::GetColorByScalar<ScalarType>(FaceDistortion,MaxAngleDist,0,FaceColor);

    //        break;

    //    case 2:
    //        Geo::UVDistortion<ScalarType>::ComputeDistortion(VertPos,VerPosUV,Connectivity,
    //                                                         Geo::UVDistortion<ScalarType>::AreaDist,
    //                                                         FaceDistortion);
    //        Geo::GetColorByScalar<ScalarType>(FaceDistortion,MaxAreaDist,0,FaceColor);
    //        break;

    //    case 3:
    //        Geo::UVDistortion<ScalarType>::ComputeDistortion(VertPos,VerPosUV,Connectivity,
    //                                                         Geo::UVDistortion<ScalarType>::EdgeDist,
    //                                                         FaceDistortion);
    //        Geo::GetColorByScalar<ScalarType>(FaceDistortion,MaxEdgeDist,0,FaceColor);
    //        break;
    //    case 4:
    //        Geo::UVDistortion<ScalarType>::ComputeDistortion(VertPos,VerPosUV,Connectivity,
    //                                                         Geo::UVDistortion<ScalarType>::ScaledJacobianDist,
    //                                                         FaceDistortion);
    //        Geo::GetColorByScalar<ScalarType>(FaceDistortion,MaxJacDist,0,FaceColor);
    //        break;
    //    default:
    //        assert(0);
    //    }

}

void LoadMesh(const std::string &path)
{
    loadedMesh=LoadOBJ(path,VertPos,Connectivity,true);
    if (loadedMesh)
    {
        std::cout<<"Loaded "<<Connectivity.size()<<" faces "<<VertPos.size()<<" vertices"<<std::endl;
        MeshBox.Init(VertPos);
        ComputeNormals(VertPos,Connectivity,FaceNormals,VertNormals);
    }
    else
        std::cout<<"Mesh not loaded"<<std::endl;
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
    {

        trackball.MouseDown((int)xpos,(int)ypos,UI::OrbitBall::BUTTON_LEFT);
    }
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {

        trackball.MouseUp((int)xpos,(int)ypos,UI::OrbitBall::BUTTON_LEFT);
    }

}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{

    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
        trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);

    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
        trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);

}


void InitGLFW_Window()
{
    if (!glfwInit())
    {
        std::cout<<"Error Initializing GLFW"<<std::endl;
        exit(0);
    }
    //return 1;

    window = glfwCreateWindow(1280, 720, "UV Mapping", NULL, NULL);

    if (window == NULL)
    {
        std::cout<<"Error Initializing GLFW Window"<<std::endl;
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

void TestRandomFlipChange()
{
    int NumTest=100000;
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Connectivity,NextF,NextE);
    std::vector<std::vector<int> > VF,VFIndex;
    Geo::getVFAdjacency(Connectivity,VF,VFIndex);
    for (int s=0;s<NumTest;s++)
    {
        int IndexF= rand()%Connectivity.size();
        int IndexE= rand()%3;
        Pos3 CurrP(IndexF,IndexE,Connectivity[IndexF][IndexE]);
        if (!Geo::TriangleEdgeFlip<ScalarType>::Feasible(CurrP,VertPos,Connectivity,
                                                          NextF,NextE,VF))continue;
        Geo::TriangleEdgeFlip<ScalarType>::Do(CurrP,Connectivity,NextF,NextE,VF,VFIndex);

        //Geo::TriangleEdgeFlip<ScalarType>::Do(CurrP,Connectivity,NextF,NextE,VF,VFIndex);
    }
    bool IsOk=Geo::CheckVFAdjacency(Connectivity,VF,VFIndex);
    assert(IsOk);
    IsOk=Geo::CheckFFAdjacencyManifold(Connectivity,NextF,NextE);
    assert(IsOk);
    ComputeNormals(VertPos,Connectivity,FaceNormals,VertNormals);

}

void TestOptimize()
{
    IteOpt.Init();
}

void TestTopologyChange()
{
    std::vector<std::vector<int> > VF,VFIndex;
    Geo::getVFAdjacency(Connectivity,VF,VFIndex);
    bool IsOk=Geo::CheckVFAdjacency(Connectivity,VF,VFIndex);
    assert(IsOk);

    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Connectivity,NextF,NextE);
    IsOk=Geo::CheckFFAdjacencyManifold(Connectivity,NextF,NextE);
    assert(IsOk);

    for (int i=0;i<(int)Connectivity.size();i++)
    {
        std::vector<int> ToUpd;
        ToUpd.push_back(i);
        Geo::UpdateVF(Connectivity,ToUpd,VF,VFIndex);
//        Geo::VFDetach(Connectivity,ToUpd,VF,VFIndex);
//        Geo::VFAttach(Connectivity,ToUpd,VF,VFIndex);
    }

    IsOk=Geo::CheckVFAdjacency(Connectivity,VF,VFIndex);


    for (int i=0;i<(int)Connectivity.size();i++)
    {
        //std::cout<<"Testing:"<<i<<std::endl;
        std::vector<int> ToUpd;
        ToUpd.push_back(i);
        for (int j=0;j<(int)NextF[i].size();j++)
        {
            int nextFIDX=NextF[i][j];
            if (nextFIDX==-1)continue;
            ToUpd.push_back(nextFIDX);
        }
        Geo::UpdateFF(Connectivity,ToUpd,NextF,NextE);

    }
    IsOk=Geo::CheckFFAdjacencyManifold(Connectivity,NextF,NextE);
    assert(IsOk);
}

void SetRenderBar()
{
    ImGui::Begin("Mesh Optimizer",NULL,ImGuiWindowFlags_AlwaysAutoResize);

    if (ImGui::Button("Test Connectivity Chance"))
    {
        TestTopologyChange();
    }
    if (ImGui::Button("Test Random Flip"))
    {
        TestRandomFlipChange();
    }
    if (ImGui::Button("Test Optimize"))
    {
        IteOpt.Init();
        IteOpt.Optimize();
        //IteOpt.Init();
        //TestTopologyChange();
    }

    ImGui::TextColored(ImVec4(1, 0, 1, 1), "OPTIMIZE");
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
    trackball.DrawPostApply();


    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    GLDraw::glScale(3/MeshBox.Diag());
    GLDraw::glTranslate(-MeshBox.Center());

    bool DrawSurface=true;
    if (showMesh)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(line_thick);


        GLDraw::glColor(Geo::Point3<ScalarType>(0.8,0.8,0.8));

        GLDraw::DrawSurfaceMesh<ScalarType>(VertPos,Connectivity,FaceNormals,
                                            VertNormals,DrawSurface,true,true,false,
                                            GLDraw::TRTriangle);//,&FaceColor);


        glPopAttrib();
    }

}

int main(int argc, char *argv[])
{

    if (argc==2){
        PathMesh=std::string(argv[1]);
        LoadMesh(PathMesh);
    }

    
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
        {
            SetRenderBar();
        }

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
