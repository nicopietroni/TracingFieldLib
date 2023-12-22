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
#include <GL_functions/GLDrawPatchData.h>

#include "IO/io_field.h"
#include <connectivity.h>
#include <field_graph/graph.h>
#include <field_graph/initialize_graph.h>
#include <field_graph/graph_tracing.h>
#include <field_graph/path_functions.h>
#include <chrono>
#include <field_graph/tracing_interface.h>
#include <field_graph/path_sampling.h>
#include "tangent_space_smooth.h"
#include <field_graph/patch_managing.h>
#include <field_graph/patch_decomposer.h>
#include <mesh_functions.h>
#include <field_graph/dual_patches.h>

GLFWwindow* window;
std::string PathMesh;
std::string PathField;
UI::OrbitBall trackball;
ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
//ImVec4 mesh_color = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
float line_thick=1.f;
int selected_rend_item = 5;

bool DrawInitalPos=false;
bool DrawSamplesPaths=false;


bool showMesh=true;
bool showCross=true;
bool showSing=true;
bool has_graph=false;
bool has_init_patches=false;
bool DrawTangentNode=false;
bool DrawExitNode=false;
bool DrawCorners=false;
bool DrawSides=false;

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
int smooth_steps=10;
int CurrN=0;
int OldN=-1;
//std::vector<int> PathN;

Geo::CrossPathSampling<ScalarType> PathSampl(VertPos,Connectivity,TGraph,VertCurv,SingIndex,SingValue);

std::vector<Geo::Point3<ScalarType> > PathColors;

std::vector<Geo::Point3<ScalarType> > FaceCol;
Geo::PatchManaging<ScalarType> PatchM(VertPos,Connectivity,TGraph,PathSampl.SelectedPaths);

int CurrCol=0;
int OldCol=CurrCol;

Geo::Point3<ScalarType> mesh_color(0.79, 0.9, 0.87);

Geo::PatchDecomposer<ScalarType> PDeco(PatchM,PathSampl);

std::vector<Geo::Point3<ScalarType> > SidesMidP;
//for debugging the marked parts
//std::vector<int> MarkedV;

void UpdateColor()
{
    switch(CurrCol) {
    case 0:
        FaceCol=std::vector<Geo::Point3<ScalarType> >(Connectivity.size(),mesh_color);
        break;
    case 1:
        PatchM.GetColorFaceByPatch(FaceCol);
        break;
    case 2:
        PatchM.GetColorFaceIfDiskLike(FaceCol);
        break;
    case 3:
        PatchM.GetColorByValence(FaceCol);
        break;
    case 4:
        PatchM.GetColorByMatchSingularity(FaceCol);
        break;
    case 5:
        PDeco.GetColorBySolvedPatches(FaceCol);
        break;
    default:
        FaceCol=std::vector<Geo::Point3<ScalarType> >(Connectivity.size(),mesh_color);
        // code block
    }
}

void InitGraph()
{
    Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(TGraph,VertPos,Connectivity,VertCurv,
                                                                   drift_penalty,maxAngle,diffuse_step);
    TGraph.DeSelectAll(0);
    TGraph.SelectNoteType(Geo::ExitDomain,0);
    TGraph.SelectNoteType(Geo::TangentDomain,0);
}

int failed_trace=0;
int tested_trace=0;

bool single_boundary_cond=true;
bool disk_like_cond=true;
bool single_sing_cond=false;
bool match_sing_cond=true;

int MinSides=3;
int MaxSides=6;
ScalarType CCability=0;

bool old_single_boundary_cond=single_boundary_cond;
bool old_disk_like_cond=disk_like_cond;
bool old_single_sing_cond=single_sing_cond;
bool old_match_sing_cond=match_sing_cond;

int old_MinSides=MinSides;
int old_MaxSides=MaxSides;
ScalarType old_CCability=CCability;

bool UpdateConditions()
{
    bool has_changed=false;
    has_changed|=old_single_boundary_cond!=single_boundary_cond;
    has_changed|=old_disk_like_cond!=disk_like_cond;
    has_changed|=old_single_sing_cond!=single_sing_cond;
    has_changed|=old_match_sing_cond!=match_sing_cond;

    has_changed|=old_MinSides!=MinSides;
    has_changed|=old_MaxSides!=MaxSides;
    has_changed|=old_CCability!=CCability;

    old_single_boundary_cond=single_boundary_cond;
    old_disk_like_cond=disk_like_cond;
    old_single_sing_cond=single_sing_cond;
    old_match_sing_cond=match_sing_cond;
    old_MinSides=MinSides;
    old_MaxSides=MaxSides;
    old_CCability=CCability;
    return has_changed;
}

void InitConditions()
{
    PDeco.ClearConditions();
    if (single_boundary_cond)
        PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
    if (disk_like_cond)
        PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
    if (single_sing_cond)
        PDeco.AddCondition(new Geo::SingleSingularityCondition<ScalarType>());
    if (match_sing_cond)
        PDeco.AddCondition(new Geo::MathchSingularityCondition<ScalarType>());
    if ((MinSides>0)||(MaxSides>0))
    {
        PDeco.AddCondition(new Geo::NumSidesCondition<ScalarType>());
    }
    if (CCability>=1)
    {
        ScalarType AvVAl=AvgEdgeLen(VertPos,Connectivity);
        PDeco.AddCondition(new Geo::CCabilityCondition<ScalarType>(AvVAl));
    }
}

void InitTracing()
{
    //PathN.clear();
    PathSampl.Init();
}

void AddAllPaths(bool upd_patches=true)
{
    //auto start = std::chrono::high_resolution_clock::now();
    PathSampl.InsertAllPaths();
    //    auto end = std::chrono::high_resolution_clock::now();
    //    // Calculate the elapsed time
    //    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    //    double elapsed_seconds = duration.count() * 1e-6; // Convert microseconds to seconds
    //    std::cout << "Adding Paths Elapsed time: " << elapsed_seconds << " seconds " <<std::endl;

    PathColors.resize(PathSampl.SelectedPaths.size());
    Geo::ScatterColor<ScalarType>(PathColors);

    std::vector<int> UpdatedPatches;
    if (upd_patches)
    {
        PatchM.Update(UpdatedPatches);
        PDeco.UpdateAllUnsolvedPatches();
        CurrCol=1;
        UpdateColor();
    }
}

void SmoothPaths()
{
    Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
    std::vector<std::vector<int> > VertPaths;
    PathSampl.GetVertexSelectedPaths(VertPaths);

    for (size_t i=0;i<VertPaths.size();i++)
    {
        if (VertPaths[i].size()==0)continue;
        for (int j=0;j<(int)VertPaths[i].size()-1;j++)
        {
            int IndexV0=VertPaths[i][j];
            int IndexV1=VertPaths[i][j+1];
            std::pair<int,int> key(std::min(IndexV0,IndexV1),
                                   std::max(IndexV0,IndexV1));
            UVP.FeatureLines.push_back(key);
        }
    }
    UVP.Steps=smooth_steps;
    Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos,Connectivity,UVP);

    //update cross F poistion
    for (size_t i=0;i<VertCurv.size();i++)
        VertCurv[i].Pos=VertPos[i];

    PatchM.UpdateSubPatchPos();
    PatchM.InitVisualCornersEdges();
}

void InitPatches()
{
    PatchM.Init();
    InitConditions();
    PDeco.UpdateAllUnsolvedPatches();
}

void BatchProcess()
{
    std::cout << "*** INIT GRAPHS *** " <<std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    InitGraph();
    auto end0 = std::chrono::high_resolution_clock::now();
    auto elapsed0 = std::chrono::duration_cast<std::chrono::microseconds>(end0 - start);

    std::cout << "*** INIT PATCHES *** " <<std::endl;
    InitPatches();
    auto end1 = std::chrono::high_resolution_clock::now();
    auto elapsed1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - end0);

    std::cout << "*** INIT TRACING *** " <<std::endl;
    InitTracing();
    auto end2 = std::chrono::high_resolution_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - end1);

    std::cout << "*** ADD ALL PATHS *** " <<std::endl;
    AddAllPaths(false);
    auto end3 = std::chrono::high_resolution_clock::now();
    auto elapsed3 = std::chrono::duration_cast<std::chrono::microseconds>(end3 - end2);

    std::cout << "*** UPDATING PATCHES *** " <<std::endl;
    std::vector<int> UpdatedPatches;
    PatchM.Update(UpdatedPatches);
    CurrCol=1;
    UpdateColor();
    PatchM.CompactEmptyPatches();

    auto end4 = std::chrono::high_resolution_clock::now();
    auto elapsed4 = std::chrono::duration_cast<std::chrono::microseconds>(end4 - end3);

    std::cout << "*** SMOOTH PATCHES *** " <<std::endl;
    SmoothPaths();
    auto end5 = std::chrono::high_resolution_clock::now();
    auto elapsed5 = std::chrono::duration_cast<std::chrono::microseconds>(end5 - end4);

    //PatchM.InitVisualCornersEdges();

    double elapsed_seconds0 = elapsed0.count() * 1e-6; // Convert microseconds to seconds
    double elapsed_seconds1 = elapsed1.count() * 1e-6; // Convert microseconds to seconds
    double elapsed_seconds2 = elapsed2.count() * 1e-6; // Convert microseconds to seconds
    double elapsed_seconds3 = elapsed3.count() * 1e-6; // Convert microseconds to seconds
    double elapsed_seconds4 = elapsed4.count() * 1e-6; // Convert microseconds to seconds
    double elapsed_seconds5 = elapsed5.count() * 1e-6; // Convert microseconds to seconds

    showCross=false;

    std::cout << "*** TIMINGS (sec) *** " <<std::endl;
    std::cout << "* Init Graph Elapsed time: " << elapsed_seconds0 <<std::endl;
    std::cout << "* Init Patches: " << elapsed_seconds1 <<std::endl;
    std::cout << "* Init Tracing: " << elapsed_seconds2  <<std::endl;
    std::cout << "* Adding Paths: " << elapsed_seconds3  <<std::endl;
    std::cout << "* Updating Patches: " << elapsed_seconds4  <<std::endl;
    std::cout << "* Smooth Patches: " << elapsed_seconds5  <<std::endl;
}

//void TestTrace()
//{
//    int currI=CurrN;
//    do{
//        if (TGraph.NodeType(currI)==Geo::Internal)
//        {
//           CurrN=currI;
//           break;
//        }
//        currI=(currI+1)%TGraph.NumNodes();
//    }while (true);
//    tested_trace++;
//    bool has_traced=Geo::GraphTracing<ScalarType>::TraceLoopOrBorderBorderFrom(TGraph,CurrN,PathN);
//    if (!has_traced) failed_trace++;
//}

//void TestSpeed()
//{
//    // Record the start time
//    auto start = std::chrono::high_resolution_clock::now();
//    failed_trace=0;
//    tested_trace=0;
//    for (CurrN=0;CurrN<(int)TGraph.NumNodes();CurrN++)
//    {
//        if ((CurrN % 100)==0)
//            std::cout<<"Testing "<<CurrN<<std::endl;
//       TestTrace();
//    }
//    // Record the end time
//    auto end = std::chrono::high_resolution_clock::now();

//    // Calculate the elapsed time
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    double elapsed_seconds = duration.count() * 1e-6; // Convert microseconds to seconds

//    std::cout << "All Trace Elapsed time: " << elapsed_seconds << " seconds " << tested_trace <<" start nodes"<<std::endl;
//    std::cout << "Failed trace " << failed_trace <<std::endl;

//    CurrN=0;
//}

//void TestSelfIntersect()
//{
//    // Record the start time
//    auto start = std::chrono::high_resolution_clock::now();

//    for (CurrN=0;CurrN<(int)TGraph.NumNodes();CurrN++)
//    {
//        if ((CurrN % 100)==0)
//            std::cout<<"Testing "<<CurrN<<std::endl;
//       TestTrace();
//       if (Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,PathN))
//       {
//           std::cout<<"Found Self Intersection Path "<<CurrN<<std::endl;
////           for (size_t i=0;i<PathN.size();i++)
////           {}
//           break;
//       }
//    }
//    // Record the end time
//    auto end = std::chrono::high_resolution_clock::now();

//    // Calculate the elapsed time
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    double elapsed_seconds = duration.count() * 1e-6; // Convert microseconds to seconds

//    std::cout << "All Trace Elapsed time: " << elapsed_seconds << " seconds " << TGraph.NumNodes() <<" start nodes"<<std::endl;

//    CurrN=0;
//}


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

    UpdateColor();
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

        Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv,minQCrossVert,maxQCrossVert);
        GetFaceSingularities(FaceCurv,VertPos,Connectivity,NextF,NextE,SingIndex,SingValue);
        PatchM.AddSingularities(SingIndex,SingValue);
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
    ImGui::Begin("Render",NULL,ImGuiWindowFlags_AlwaysAutoResize);       // Create a window called "Hello, world!" and append into it.

    ImGui::Separator();
    ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW MESH");
    ImGui::Checkbox("Draw Mesh", &showMesh);
    ImGui::Checkbox("Draw Cross", &showCross);
    ImGui::Checkbox("Draw Singularities", &showSing);
    ImGui::Checkbox("Draw Corners", &DrawCorners);
    ImGui::Checkbox("Draw Sides", &DrawSides);

    ImGui::Separator();
    ImGui::TextColored(ImVec4(1, 0, 1, 1), "DRAW GRAPH");
    ImGui::Checkbox("Draw Exit Node", &DrawExitNode);
    ImGui::Checkbox("Draw Tangent Node", &DrawTangentNode);
    ImGui::TextColored(ImVec4(1, 0, 1, 1), "PATCHES");
    static const char* itemsShape[] = { "Color Uniform",
                                        "Color By Patch",
                                        "Color if Disk",
                                        "Color Valence",
                                        "Color Match Sing",
                                        "Color Unsolved Patches"};

    ImGui::End();
}

void SetDecomposeBar()
{
    ImGui::Begin("Decompose",NULL,ImGuiWindowFlags_AlwaysAutoResize);       // Create a window called "Hello, world!" and append into it.

    ImGui::TextColored(ImVec4(1, 0, 1, 1), "INIT");
    ImGui::SliderFloat("Drift",&drift_penalty, 0, 1000);
    ImGui::SliderFloat("Max Angle",&maxAngle, 0, 90);
    ImGui::InputInt("Diffuse Step",&diffuse_step);
    if (ImGui::Button("Init Graph"))
    {
        InitGraph();
        has_graph=true;
        showCross=false;
    }
    if (!has_graph){
        ImGui::End();
        return;
    }

    ImGui::InputInt("Initial Sample",&PathSampl.NumSamples);
    if (ImGui::Button("Init Tracing"))
    {
        InitTracing();
        InitPatches();
        has_init_patches=true;
    }
    if (!has_init_patches){
        ImGui::End();
        return;
    }
    //    if (ImGui::Button("Init Patches"))
    //    {
    //        InitPatches();
    //    }

    ImGui::Separator();
    ImGui::TextColored(ImVec4(1, 0, 1, 1), "PATH SAMPLING");

    ImGui::Checkbox("Dynamic Update", &PathSampl.dynamicUpdates);
    ImGui::Checkbox("Only Needed", &PDeco.add_only_needed);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "structure");
    ImGui::Checkbox("Darts", &PDeco.allow_darts);
    ImGui::SameLine();
    ImGui::Checkbox("Concave", &PDeco.allow_concave);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "topology");
    ImGui::Checkbox("1 Boundary", &single_boundary_cond);
    ImGui::SameLine();
    ImGui::Checkbox("Disk Like", &disk_like_cond);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "singularity");
    ImGui::Checkbox("1 Sing Patch", &single_sing_cond);
    ImGui::SameLine();
    ImGui::Checkbox("Sing Match", &match_sing_cond);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "sides");
    ImGui::SetNextItemWidth(50);
    ImGui::InputInt("Min Sides",&MinSides,0);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(50);
    ImGui::InputInt("Max Sides",&MaxSides,0);

    ImGui::TextColored(ImVec4(0, 0, 1, 1), "shape");
    ImGui::SetNextItemWidth(100);
    float CCAbilityf=CCability;
    ImGui::InputFloat("CC Ability", &CCAbilityf);
    CCability=(ScalarType)CCAbilityf;

    bool needUpdate=UpdateConditions();
    if (needUpdate)
        InitConditions();
    //    if (ImGui::Button("Add One"))
    //    {
    //        PathSampl.AddNextPath();
    //        PathColors.resize(PathSampl.SelectedPaths.size());
    //        Geo::ScatterColor<ScalarType>(PathColors);
    //        std::vector<int> UpdatedPatches;
    //        PatchM.Update(UpdatedPatches);
    //        UpdateColor();
    //        //PatchM.InitVisualCornersEdges();
    //    }

    //    if (ImGui::Button("Add All Paths"))
    //        AddAllPaths();

    if (ImGui::Button("Add Next"))
    {
        bool terminated=false;
        bool solved=false;
        PDeco.AddNext(solved,terminated);

        PathColors.resize(PathSampl.SelectedPaths.size());
        Geo::ScatterColor<ScalarType>(PathColors);
        UpdateColor();
    }
    ImGui::SameLine();
    if (ImGui::Button("Resample"))
    {
        PDeco.RetraceOnUnsolved();
    }
    ImGui::SameLine();
    if (ImGui::Button("Remove"))
    {
        PDeco.RemoveNextIfPossible();

        PathColors.resize(PathSampl.SelectedPaths.size());
        Geo::ScatterColor<ScalarType>(PathColors);
        UpdateColor();
    }
    //    if (ImGui::Button("Batch Process"))
    //        BatchProcess();
    //    ImGui::Separator();
    //    ImGui::TextColored(ImVec4(1, 0, 1, 1), "FINALIZE");

    //    ImGui::InputInt("Smooth Steps",&smooth_steps);
    //    if (ImGui::Button("Smooth"))
    //        SmoothPaths();

    //    if (ImGui::Button("Extract Dual"))
    //        Geo::DualPatch<ScalarType>::ExtractDual(VertPos,Connectivity,PatchM,SidesMidP);
    //        //PatchM.ExtractDual();

    //    if (ImGui::Button("Save Patches"))
    //        PatchM.WriteSubMeshes(std::string("test.obj"));

    ImGui::End();
}

void SetExportBar()
{
    ImGui::Begin("Export",NULL,ImGuiWindowFlags_AlwaysAutoResize);       // Create a window called "Hello, world!" and append into it.

    ImGui::InputInt("Smooth Steps",&smooth_steps);
    if (ImGui::Button("Smooth"))
        SmoothPaths();

    if (ImGui::Button("Extract Dual"))
        Geo::DualPatch<ScalarType>::ExtractDual(VertPos,Connectivity,PatchM,SidesMidP);
    //PatchM.ExtractDual();

    if (ImGui::Button("Save Patches"))
        PatchM.WriteSubMeshes(std::string("test.obj"));

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
        //glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w, mesh_color.z * mesh_color.w, mesh_color.w);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(line_thick);
        GLDraw::DrawSurfaceMesh<ScalarType>(VertPos,Connectivity,FaceNormals,VertNormals,true,true,true,false,GLDraw::TRTriangle,&FaceCol);
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


    //GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodes(TGraph,VertCurv,PathN,MeshBox.Diag()*0.005);


    GLGraphDraw::GLCrossFieldGraph<ScalarType>::PathDrawParam PDParam;
    PDParam.DrawPath=true;
    GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawPaths(TGraph,VertCurv,PathSampl.
                                                          SelectedPaths,MeshBox.Diag()*0.005,
                                                          PDParam,PathColors);

    PDParam.DrawPath=false;

    if (DrawCorners)
        GLGraphDraw::GLDrawPatchData<ScalarType>::DrawCorners(PatchM);

    if (DrawSides)
        GLGraphDraw::GLDrawPatchData<ScalarType>::DrawSides(PatchM);

    if (DrawInitalPos)
    {
        std::vector<int> StartIDX;
        PathSampl.GetStartingNodes(StartIDX);
        GLDraw::DrawVertices(VertPos,StartIDX);
    }
    GLDraw::DrawPoints(SidesMidP);
    //     GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodesOfType(TGraph,VertCurv,
    //                                                                 Geo::ExitDomain,MeshBox.Diag()*0.005);
    //bool DrawSamplesPaths=false;

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
        //ImGui::SetNextWindowPos(ImVec2(20, 20));
        SetDecomposeBar();
        int windowWidth, windowHeight;
        glfwGetWindowSize(window, &windowWidth, &windowHeight);
        //ImGui::SetNextWindowPos(ImVec2(windowWidth-200, 20));
        SetRenderBar();
        //ImGui::Render();

        glfwGetWindowSize(window, &windowWidth, &windowHeight);
        //ImGui::SetNextWindowPos(ImVec2(windowWidth-400, windowHeight-200));
        SetExportBar();

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
