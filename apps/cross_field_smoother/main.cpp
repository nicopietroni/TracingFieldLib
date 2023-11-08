#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <stdio.h>
#include <UI/orbitball.h>
#include "ImGuiFileBrowser.h"
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>
#include "IO/io_obj.h"
#include <point3.h>
#include <box3.h>
#include <GL_functions/GLDrawMesh.h>
#include <GL_functions/GLDrawField.h>
#//include "eigen_interface.h"
#include "IO/io_field.h"
#include "IO/path_functions.h"
#include <connectivity.h>
#include <smooth_cross_field.h>

GLFWwindow* window;
std::string PathMesh;
std::string PathField;
UI::OrbitBall trackball;
imgui_addons::ImGuiFileBrowser file_dialog; 
//std::ostringstream oss;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
ImVec4 mesh_color = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
float line_thick=1.f;
int selected_rend_item = 5;

bool showMesh=true;
bool showCross=true;
bool showSing=true;


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
int Nring=3;
ScalarType minQCrossVert,maxQCrossVert;
std::vector<int> SingIndex,SingValue;
SmoothParam<ScalarType> SParam;

// Interface values
float TopIntevalFix=(float)SParam.TopIntevalFix;
float GlobalSmoothVal=(float)SParam.GlobalSmoothVal;


void InitFieldByCurvature()
{
    ComputeCurvatureField<ScalarType>(VertPos,Connectivity,FaceNormals,VertNormals,VertCurv,FaceCurv,Nring);
    Field::CrossF<ScalarType>::getMinMaxQ(VertCurv,minQCrossVert,maxQCrossVert);
    GetFaceSingularities(FaceCurv,VertPos,Connectivity,NextF,NextE,SingIndex,SingValue);
    has_cross_field=true;
}

void SaveField()
{
    assert(!PathMesh.empty());
    std::string PathField=RemoveExtension(PathMesh);
    PathField+=std::string(".field");
    std::cout<<"Saving Field on file "<<PathField.c_str()<<std::endl;
    bool Saved=WriteField(PathField,VertCurv,FaceCurv);
    if (!Saved)
    {
        std::cout<<"ERROR SAVING FIELD FILE "<<std::endl;
        exit(0);
    }
}

void SmoothGlobal()
{
    if (!has_cross_field)return;
    SmoothGlobalPolyvector(VertPos,Connectivity,FaceNormals,FaceCurv,SParam);
    GetFaceSingularities(FaceCurv,VertPos,Connectivity,NextF,NextE,SingIndex,SingValue);
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
    trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
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

// Now inside any function
void showMainMenu()
{
    bool open = false, save = false;
    if(ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("File"))
        {
            if (ImGui::MenuItem("Open", NULL))
                open = true;
            if (ImGui::MenuItem("Quit", NULL))
                exit(0);
            
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
    
    //Remember the name to ImGui::OpenPopup() and showFileDialog() must be same...
    if(open)
        ImGui::OpenPopup("Open File");
    // if(save)
    //     ImGui::OpenPopup("Save File");

    /* Optional third parameter. Support opening only compressed rar/zip files.
     * Opening any other file will show error, return false and won't close the dialog.
     */
    if(file_dialog.showFileDialog("Open File", imgui_addons::ImGuiFileBrowser::DialogMode::OPEN, ImVec2(700, 310), ".obj,.ply,.off"))
    {
        LoadMesh(file_dialog.selected_path);
    }
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
    ImGui::Begin("Field Smoother",NULL,ImGuiWindowFlags_AlwaysAutoResize);       // Create a window called "Hello, world!" and append into it.

    ImGui::SliderInt("Ring",&Nring, 2, 5, "%d",ImGuiSliderFlags_AlwaysClamp);
    if (ImGui::Button("Init By Curvature"))
    {
        if (loadedMesh)
            InitFieldByCurvature();
    }

    if (has_cross_field)
    {
        if (ImGui::Button("Smooth Field"))
        {
            if (loadedMesh)
                SmoothGlobal();
        }
        //ImGui::SameLine();

        ImGui::SliderFloat("Top Fix",&TopIntevalFix, 0, 0.3, "%2.f",ImGuiSliderFlags_AlwaysClamp);
        SParam.TopIntevalFix= (ScalarType) TopIntevalFix;

        ImGui::SliderFloat("Smoothness",&GlobalSmoothVal, 50, 500, "%1.f",ImGuiSliderFlags_AlwaysClamp);
        SParam.GlobalSmoothVal=GlobalSmoothVal;

        if (ImGui::Button("Save Field"))
        {
            if (loadedMesh)
                SaveField();
        }
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

    //    UsedScalarType scaleVal=MeshBox.Diag()*(UsedScalarType)0.02;
    //    vcg::glScale(2.0f/mesh.bbox.Diag());
    //    vcg::glTranslate(-mesh.bbox.Center());
    bool DrawSurface=true;
    //    bool DrawWire=true;
    //    std::cout<<"Conn "<<Connectivity.size()<<std::endl;
    //    std::cout<<"FNorm "<<FaceNormals.size()<<std::endl;


    //mesh.GLDrawBox();
    if (showMesh)
    {
        glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w, mesh_color.z * mesh_color.w, mesh_color.w);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(line_thick);
        //        glWrap.SetHintParamf(vcg::GLW::HNPPointSize,line_thick);
        //        glWrap.Draw(DMode,vcg::GLW::CMNone,vcg::GLW::TMNone);
        GLDraw::DrawSurfaceMesh<ScalarType>(VertPos,Connectivity,FaceNormals,VertNormals,DrawSurface,false,true,false);
        glPopAttrib();

    }


    if (has_cross_field)
    {
        ScalarType scaleVal=MeshBox.Diag()*(ScalarType)0.02;
        if (showCross)
            GLDraw::DrawCrossFields<ScalarType>(FaceCurv,scaleVal,maxQCrossVert,minQCrossVert);

        if (showSing)
            GLDraw::DrawSingularity<ScalarType>(VertPos,SingIndex,SingValue);
    }
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
        {
            SetRenderBar();

            //SetLogWindow();

            showMainMenu();
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
