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
#include <GL_functions/GLPick.h>
#include <arap.h>
#include <pos3.h>
#include <solve/linear_constraint.h>

/*
#include <GLDrawField.h>
#include <smooth_cross_field.h>*/

GLFWwindow* window;
std::string PathMesh;//="../test_data/grid.obj";
//MyTriMesh mesh;
//vcg::GlTrimesh<MyTriMesh> glWrap;
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

//bool showTrackBall=true;
//vcg::GLW::DrawMode DMode=vcg::GLW::DMSmooth;

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
std::vector<Geo::Point3<ScalarType> > VertRestPos;

std::vector<std::vector<int> > Connectivity;
std::vector<Geo::Point3<ScalarType> > FaceNormals;
std::vector<Geo::Point3<ScalarType> > VertNormals;

//std::vector<std::vector<int> > NextF,NextE;
Geo::Box3<ScalarType> MeshBox;

//the constraints of ARAP
std::map<int,Geo::Point3<ScalarType> > ConstraintVertPos;
std::vector<int> OldFixed;

Geo::ARAP<ScalarType> ArapSolver;

bool to_pick=false;

enum ConstrMode{CMNone,CMFix,CMZero,CMDrag,CMRemove};

ConstrMode CMode=CMNone;

double xPick, yPick,dPick;
int PickedV=-1;

//std::vector<std::vector<int> > KernelV0;
std::vector<std::vector<int> > KernelV1;

bool useCotangent=true;
bool useRotInvadiant=true;

void GetFixedArray(std::vector<int> &fixedVect)
{
    fixedVect.clear();
    std::map<int,Geo::Point3<ScalarType> >::iterator iteMap;
    for (iteMap=ConstraintVertPos.begin();iteMap!=ConstraintVertPos.end();iteMap++)
        fixedVect.push_back((*iteMap).first);
    std::sort(fixedVect.begin(),fixedVect.end());
}

void UpdateArap(bool forceUpdate=false)
{
    std::vector<int> NewFixed;
    GetFixedArray(NewFixed);
    bool HasToUpdate=OldFixed!=NewFixed;
    OldFixed=NewFixed;

    if (HasToUpdate || forceUpdate)
    {
        std::vector<Geo::LinearConstraint<ScalarType> > LinConstr;
        if (!useCotangent)
        {
            std::vector<std::vector<ScalarType> > Weight;
            Geo::SetUniformLaplacianWeightForMesh(Connectivity,VertRestPos.size(),Weight);
            ArapSolver.Init(VertRestPos,Connectivity,NewFixed,LinConstr,Weight);
        }
        else
            ArapSolver.Init(VertRestPos,Connectivity,NewFixed,LinConstr);
        //TestRotConstr();
    }

    std::vector<Geo::Point3<ScalarType> > FixPos;
    for (int i=0;i<NewFixed.size();i++)
    {
        int IndexV=NewFixed[i];
        assert(ConstraintVertPos.count(IndexV)>0);
        Geo::Point3<ScalarType> ContrPos=ConstraintVertPos[IndexV];
        FixPos.push_back(ContrPos);
    }

    ArapSolver.IterateUntilConvergence(VertPos,FixPos,MeshBox.Diag()*0.001,8,false,useRotInvadiant);
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
        //GetFaceFaceConnectivity(Connectivity,NextF,NextE);

        VertRestPos=VertPos;
        //Geo::FindVertexStar(Connectivity,KernelV0,1);
        Geo::FindVertexStar(Connectivity,KernelV1,3);
        //std::cout<<"FNorm 0 "<<FaceNormals.size()<<std::endl;
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

    xPick=xpos;
    yPick=ypos;

    if (CMode!=CMNone)return;

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
    {
        trackball.MouseUp((int)xpos,(int)ypos,UI::OrbitBall::BUTTON_LEFT);
        //CMode=CMNone;
        PickedV=-1;
    }

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && (CMode!=CMNone))
    {
        xPick=xpos;
        yPick=ypos;
        to_pick=true;
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_TAB && action == GLFW_PRESS)
        CMode=CMRemove;

    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
        CMode=CMFix;

    if (key == GLFW_KEY_0 && action == GLFW_PRESS)
        CMode=CMZero;

    if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
        CMode=CMNone;

    if (key == GLFW_KEY_0 && action == GLFW_RELEASE)
        CMode=CMNone;

    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_PRESS)
        CMode=CMDrag;

    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_RELEASE)
        CMode=CMNone;

    if (key == GLFW_KEY_TAB && action == GLFW_RELEASE)
        CMode=CMNone;

    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
        trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);

    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
        trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);

    if (key== GLFW_KEY_S)
        WriteOBJ<ScalarType>(std::string("deformed.obj"),VertPos,Connectivity);
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
    ImGui::Begin("Arap Deformer",NULL,ImGuiWindowFlags_AlwaysAutoResize);       // Create a window called "Hello, world!" and append into it.
    bool changedDevelop=ImGui::Checkbox("Use Cotangent W", &useCotangent);
    bool changedRotation=ImGui::Checkbox("Use Rotation Invariant", &useRotInvadiant);
    if (changedDevelop || changedRotation)
        UpdateArap(true);

}

void ConstraintVertKernel(const size_t &IndexV)
{
    assert(IndexV>=0);
    ConstraintVertPos[IndexV]=VertRestPos[IndexV];
}

void ConstraintVertKernelZero(const size_t &IndexV)
{
    assert(IndexV>=0);
    assert(IndexV<KernelV1.size());
    for (size_t i=0;i<KernelV1[IndexV].size();i++)
    {
        int CurrV=KernelV1[IndexV][i];
        Geo::Point3<ScalarType> CurrP=VertPos[CurrV];
        CurrP.Y=0;
        ConstraintVertPos[CurrV]=CurrP;
    }
}

void RemoveVertKernel(const size_t &IndexV)
{
    assert(IndexV>=0);
    assert(IndexV<KernelV1.size());
    for (size_t i=0;i<KernelV1[IndexV].size();i++)
    {
        int CurrV=KernelV1[IndexV][i];
        ConstraintVertPos.erase(CurrV);
    }
}


void AddDisplVertKernel(const size_t &IndexV,
                        const Geo::Point3<ScalarType> &TargetP)
{
    assert(IndexV>=0);
    Geo::Point3<ScalarType> Displ=TargetP-VertRestPos[IndexV];
    ConstraintVertPos[IndexV]=VertRestPos[IndexV]+Displ;
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
    trackball.DrawPostApply();

    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    GLDraw::glScale(3/MeshBox.Diag());
    GLDraw::glTranslate(-MeshBox.Center());

    bool DrawSurface=true;

    if (showMesh)
    {
        glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w, mesh_color.z * mesh_color.w, mesh_color.w);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(line_thick);
        GLDraw::DrawSurfaceMesh<ScalarType>(VertPos,Connectivity,FaceNormals,VertNormals,DrawSurface,true,true,false);

        if ((to_pick)&&((CMode==CMFix)||(CMode==CMRemove)||(CMode==CMZero)))
        {
            ScalarType dPickTest;
            PickedV=GLDraw::GLPick<ScalarType>::PickVert(xPick,yPick,VertPos,dPickTest);
            if (PickedV>=0)
            {
                if (CMode==CMFix)
                    ConstraintVertKernel(PickedV);
                if (CMode==CMZero)
                    ConstraintVertKernelZero(PickedV);
                if (CMode==CMRemove)
                    RemoveVertKernel(PickedV);

                UpdateArap();
            }
            to_pick=false;
            PickedV=-1;
            //dPick=-1;
        }

        if ((to_pick)&&(CMode==CMDrag))
        {
            PickedV=GLDraw::GLPick<ScalarType>::PickVert(xPick,yPick,VertPos,dPick);
            if (PickedV>=0)
            {
                AddDisplVertKernel(PickedV,VertPos[PickedV]);
                UpdateArap();
            }
            to_pick=false;
        }

        if ((!to_pick)&&(CMode==CMDrag)&&(PickedV>=0))
        {
            Geo::Point3<ScalarType> TargetP=GLDraw::GLPick<ScalarType>::Unproject(xPick,yPick,dPick);
            //ConstraintVertPos[PickedV]=TargetP;
            AddDisplVertKernel(PickedV,TargetP);
            UpdateArap();
            //UpdateLaplacianMatrixIfNeeded();
        }
        glPopAttrib();

    }

    std::map<int,Geo::Point3<ScalarType> >::iterator IteConstr;
    for (IteConstr=ConstraintVertPos.begin();IteConstr!=ConstraintVertPos.end();IteConstr++)
    {

        int IndexV=(*IteConstr).first;
        assert(IndexV>=0);
        assert(IndexV<(int)VertPos.size());
        Geo::Point3<ScalarType> P0=VertPos[IndexV];
        Geo::Point3<ScalarType> P1=(*IteConstr).second;
        GLDraw::GLDrawSegment<ScalarType>(P0,P1,5,20,Geo::Point3<ScalarType>(1,0,0),Geo::Point3<ScalarType>(1,0,0));
    }
}

int main(int argc, char *argv[])
{

    //oss<<"Initialized Log Window"<<std::endl;
    //AScene.SetOutputStream(&oss);

    if (argc==2){
        PathMesh=std::string(argv[1]);
        LoadMesh(PathMesh);
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
