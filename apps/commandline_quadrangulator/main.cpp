#include <sstream>
#include <iostream>
#include <point3.h>
#include <IO/io_obj.h>
#include "IO/io_field.h"
#include <IO/file_path_functions.h>
#include <field_graph/graph.h>
#include <field_graph/initialize_graph.h>
#include <field_graph/graph_tracing.h>
#include <field_graph/path_functions.h>
#include <field_graph/tracing_interface.h>
#include <field_graph/path_sampling.h>
#include <field_graph/patch_managing.h>
#include <field_graph/patch_decomposer.h>
#include <mesh_functions.h>
#include <cross_field_functions.h>
#include <smooth_cross_field.h>
#include <cross_field.h>
#include <field_graph/mesh_preprocess.h>
#include <field_graph/dual_patches.h>

// GLFWwindow *window;
std::string PathMesh;
std::string PathField;
std::string ProjectName;

typedef double ScalarType;

// // MESH DATA
bool loadedMesh = false;
std::vector<Geo::Point3<ScalarType>> VertPos;
std::vector<std::vector<int>> Connectivity;
std::vector<Geo::Point3<ScalarType>> FaceNormals;
std::vector<Geo::Point3<ScalarType>> VertNormals;
std::vector<std::vector<int>> NextF, NextE;
Geo::Box3<ScalarType> MeshBox;

// // FIELD DATA
std::vector<Field::CrossF<ScalarType>> VertCurv;
std::vector<Field::CrossF<ScalarType>> FaceCurv;
bool has_cross_field = false;
bool use_curvature = true;
bool dualize = false;
bool smooth_tangent = false;

ScalarType minQCrossVert, maxQCrossVert;

std::vector<int> SingIndex, SingValue;

Geo::TracingGraph<ScalarType> TGraph(Geo::GT_CrossField);

Geo::CrossPathSampling<ScalarType> PathSampl(VertPos, Connectivity, TGraph, VertCurv, SingIndex, SingValue);

Geo::PatchManaging<ScalarType> PatchM(VertPos, Connectivity, VertCurv, TGraph, PathSampl.SelectedPaths);

Geo::PatchDecomposer<ScalarType> PDeco(PatchM, PathSampl, VertCurv, FaceCurv);

ScalarType drift_penalty = 100.f;
ScalarType maxAngle = 45.f;
int diffuse_step = 1;
int smooth_steps = 10;

bool remove_batch = true;
bool remove_split = false;

int Nring = 3;

Geo::DualPatch<ScalarType>::DualParameter DualParam;
std::vector<Geo::Point3<ScalarType>> DualVertPos;
std::vector<std::vector<int>> DualConnectivity;

SmoothParam<ScalarType> SParam;
int min_sing_distance = 1;

bool LoadMesh(const std::string &path)
{
    //    oss<<"Loading "<<path.c_str()<<std::endl;
    loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
    if (loadedMesh)
    {
        std::cout << "Loaded " << Connectivity.size() << " faces " << VertPos.size() << " vertices" << std::endl;
        MeshBox.Init(VertPos);
        ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
        Geo::getFFAdjacency(Connectivity, NextF, NextE);

        // std::cout<<"FNorm 0 "<<FaceNormals.size()<<std::endl;
    }
    else
    {
        std::cout << "Mesh not loaded" << std::endl;
    }
    return loadedMesh;
    // UpdateColor();
}

bool LoadField(const std::string &path)
{
    // std::cout<<"Loading Field"<<path.c_str()<<std::endl;
    has_cross_field = ReadField(VertPos, Connectivity, path, VertCurv, FaceCurv);
    if (!has_cross_field)
    {
        std::cout << "Field not loaded" << std::endl;
        return false;
    }

    if (VertCurv.size() == 0)
    {
        SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);
        InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);
    }
    if (FaceCurv.size() == 0)
    {
        SetFaceCrossVectorFromVert(Connectivity, FaceNormals, VertCurv, FaceCurv);
        InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);
    }
    if (has_cross_field)
    {
        std::cout << "Loaded Field" << std::endl;

        Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert, maxQCrossVert);
        GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex, SingValue);
        PatchM.SetSingularities(SingIndex, SingValue);
    }
    return true;
}

std::string folderProjectName;
std::string commandLine;

void SaveSubMeshes()
{
    std::string subMeshFilename = folderProjectName;
    subMeshFilename.append("/bordered_layout.obj");
    std::cout << "Saving on:" << subMeshFilename.c_str() << std::endl;

    // subMeshFilename.append("_submesh.obj");
    PDeco.SaveSubMeshes(subMeshFilename);
}

void SavePatches()
{
    std::string PatchFilename = ProjectName;
    PatchFilename.append("/patch_layout.patch");
    PDeco.SavePatchesLayout(PatchFilename);
}

void SaveField()
{
    std::string FinalFieldFilename = ProjectName;
    FinalFieldFilename.append("/mesh_layout.ffield");
    WriteField(FinalFieldFilename.c_str(), VertCurv, FaceCurv);
}

void SaveMesh()
{
    std::string FinalMFilename = ProjectName;
    FinalMFilename.append("/mesh_layout.obj");
    WriteOBJ(FinalMFilename.c_str(), VertPos, Connectivity);
}

void SaveDual()
{
    if (!dualize)
        return;
    std::string FinalMFilename = ProjectName;
    FinalMFilename.append("/mesh_dual_layout.obj");
    WriteOBJ(FinalMFilename.c_str(), DualVertPos, DualConnectivity);
}

void SaveAll()
{
    SaveSubMeshes();
    SavePatches();
    SaveField();
    SaveMesh();
    SaveDual();
}

void CreateOutputFolder()
{
    folderProjectName = ProjectName;

    try
    {
        // Check if the folder already exists
        if (!std::filesystem::exists(folderProjectName))
        {
            // Create the folder
            std::filesystem::create_directory(folderProjectName);
            std::cout << "Folder created successfully: " << folderProjectName << std::endl;
        }
        else
        {
            std::cout << "Folder already exists: " << folderProjectName << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error creating folder: " << e.what() << std::endl;
    }
}

void UpdateAfterRemesh()
{
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);

    InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

    SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

    Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert, maxQCrossVert);
    GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex, SingValue);
    PatchM.SetSingularities(SingIndex, SingValue);
    // UpdateColor();
}

bool single_boundary_cond = true;
bool disk_like_cond = true;
bool single_sing_cond = false;
bool match_sing_cond = true;
int MinSides = 3;
int MaxSides = 6;

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
    if (dualize)
        PDeco.AddCondition(new Geo::DualMeshingCondition<ScalarType>());
    if ((MinSides > 0) || (MaxSides > 0))
    {
        PDeco.AddCondition(new Geo::NumSidesCondition<ScalarType>());
    }
}

void InitPatches()
{
    PatchM.Init();
    InitConditions();
    // PDeco.InitBorderVertType();
    PDeco.Init();
    PDeco.UpdateAllUnsolvedPatches();
}

void InitGraph()
{
    Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(TGraph, VertPos, Connectivity, VertCurv,
                                                                   drift_penalty, maxAngle, diffuse_step);
}

void SmoothPaths()
{
    Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
    std::vector<std::vector<int>> VertPaths;
    PathSampl.GetVertexSelectedPaths(VertPaths);

    for (size_t i = 0; i < VertPaths.size(); i++)
    {
        if (VertPaths[i].size() == 0)
            continue;
        for (int j = 0; j < (int)VertPaths[i].size() - 1; j++)
        {
            int IndexV0 = VertPaths[i][j];
            int IndexV1 = VertPaths[i][j + 1];
            std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                    std::max(IndexV0, IndexV1));
            UVP.FeatureLines.push_back(key);
        }
    }
    UVP.Steps = smooth_steps;
    Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity, UVP);

    //    //update cross F poistion
    //    for (size_t i=0;i<VertCurv.size();i++)
    //        VertCurv[i].Pos=VertPos[i];

    AdaptFieldToMesh(VertPos, Connectivity, VertCurv, FaceCurv);

    PatchM.UpdateSubPatchPos();
    PatchM.InitVisualCornersEdges();
}

void SmoothFieldGlobal(bool iterative = false)
{
    if (!has_cross_field)
        return;

    std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;
    std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;

    if (!iterative)
        SmoothGlobalPolyvector(VertPos, Connectivity, FaceNormals, FaceCurv, SParam);
    else
        SmoothCrossIterative(VertPos, Connectivity, FaceNormals, FaceCurv, SParam);

    GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex, SingValue);
    SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);

    RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
    RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

    InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
}

void InitFieldByCurvature()
{
    SParam.FixedCross.clear();
    ComputeCurvatureField<ScalarType>(VertPos, Connectivity, FaceNormals, VertNormals, VertCurv, FaceCurv, Nring);
    Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert, maxQCrossVert);
    GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex, SingValue);

    has_cross_field = true;
}

void BatchProcessFieldCurv()
{
    SParam.FixedCross.clear();
    InitFieldByCurvature();
    SmoothFieldGlobal();
}

void BatchProcessFieldBord()
{
    InitFieldByCurvature();
    SParam.align_borders = true;
    SParam.TopIntevalFix = 0;
    SParam.GlobalSmoothVal = -1;
    SmoothFieldGlobal();
    std::cout << "There are " << SingIndex.size() << " singularities" << std::endl;
}

void SmoothTessellation()
{
    Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    // has_remeshed = true;
}

void SplitBorderEdges()
{
    RefineMesh<ScalarType>::SpliToHaveOneBorderPerTris(VertPos, Connectivity);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
}

void FindDual()
{
    std::vector<std::vector<int>> PerPatchVertIDX;
    std::vector<std::vector<int>> PerPatchFaceIDX;
    std::vector<std::vector<int>> PerPatchCornerIDX;
    Geo::DualPatch<ScalarType>::DisambiguateForDual(PatchM);
    Geo::DualPatch<ScalarType>::ExtractDual(VertPos, Connectivity, PatchM,
                                            DualParam, DualVertPos, DualConnectivity,
                                            PerPatchVertIDX, PerPatchFaceIDX,
                                            PerPatchCornerIDX);
}

void BatchProcess()
{
    auto t0 = std::chrono::high_resolution_clock::now();
    if (!has_cross_field)
    {
        if (smooth_tangent)
        SmoothTessellation();

        std::cout << "*** COMPUTING FIELD *** " << std::endl;
        SplitBorderEdges();

        if (use_curvature)
            std::cout << "usign curvature " << std::endl;
        else
            std::cout << "using border " << std::endl;
        if (use_curvature)
            BatchProcessFieldCurv();
        else
            BatchProcessFieldBord();
    }
    std::cout << "*** REFINE STEP *** " << std::endl;
    Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(VertPos, Connectivity,
                                                             FaceCurv, SingIndex,
                                                             min_sing_distance);
    MeshBox.Init(VertPos);
    UpdateAfterRemesh();
    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed0 = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

    std::cout << "*** INIT GRAPHS *** " << std::endl;
    InitGraph();
    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed1 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

    std::cout << "*** INIT PATCHES *** " << std::endl;
    InitPatches();
    auto t3 = std::chrono::high_resolution_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);

    std::cout << "*** ADD ALL PATHS *** " << std::endl;
    PDeco.BatchDecompose(remove_batch, remove_split);
    auto t4 = std::chrono::high_resolution_clock::now();
    auto elapsed3 = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);

    std::cout << "*** SMOOTH PATCHES *** " << std::endl;
    SmoothPaths();
    auto t5 = std::chrono::high_resolution_clock::now();
    auto elapsed4 = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4);

    if (dualize)
    {
        std::cout << "*** DUALIZE *** " << std::endl;
        FindDual();
    }
    auto t6 = std::chrono::high_resolution_clock::now();
    auto elapsed5 = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5);

    // Convert microseconds to seconds
    double elapsed_seconds0 = elapsed0.count() * 1e-6;
    double elapsed_seconds1 = elapsed1.count() * 1e-6;
    double elapsed_seconds2 = elapsed2.count() * 1e-6;
    double elapsed_seconds3 = elapsed3.count() * 1e-6;
    double elapsed_seconds4 = elapsed4.count() * 1e-6;
    double elapsed_seconds5 = elapsed4.count() * 1e-6;

    std::cout << "*** TIMINGS (sec) *** " << std::endl;
    std::cout << "* Field & Remd time: " << elapsed_seconds0 << std::endl;
    std::cout << "* Init Graph Elapsed time: " << elapsed_seconds1 << std::endl;
    std::cout << "* Init Patches: " << elapsed_seconds2 << std::endl;
    std::cout << "* Adding Paths: " << elapsed_seconds3 << std::endl;
    std::cout << "* Smooth Patches: " << elapsed_seconds4 << std::endl;
    std::cout << "* Dualize Patches: " << elapsed_seconds5 << std::endl;

}

int main(int argc, char *argv[])
{

    std::stringstream commandLineStr;
    for (int i = 0; i < argc; ++i)
    {
        commandLineStr << argv[i] << " ";
    }
    commandLine = commandLineStr.str();

    // // oss<<"Initialized Log Window"<<std::endl;
    // // AScene.SetOutputStream(&oss);

    if (argc < 2)
    {
        std::cout << "You should pass at least a mesh" << std::endl;
        exit(0);
    }

    for (int i = 1; i < argc; i++)
    {
        if (GetFileExtension(std::string(argv[i])) == std::string("obj"))
        {
            PathMesh = std::string(argv[1]);
            bool has_load = LoadMesh(PathMesh);
            Geo::TracingPreprocess<ScalarType>::DesambiguateInitialDarts(VertPos, Connectivity);

            if (!has_load)
                exit(0);
            continue;
        }
        if ((GetFileExtension(std::string(argv[i])) == std::string("ffield")) ||
            (GetFileExtension(std::string(argv[i])) == std::string("field")))
        {
            PathField = std::string(argv[i]);
            std::cout << "Loading Field" << PathField.c_str() << std::endl;
            LoadField(PathField);
            std::cout << "DONE" << std::endl;
            // Used_curv = false;
            if (!has_cross_field)
            {
                std::cout << "ERROR LOADING FIELD" << std::endl;
                exit(0);
            }
            continue;
        }
        if (std::string(argv[i]) == std::string("dual"))
        {
            dualize = true;
        }
        if (std::string(argv[i]) == std::string("even"))
        {
            DualParam.RType = Geo::DualPatch<ScalarType>::RTypeEven;
        }
        if (std::string(argv[i]) == std::string("ortho"))
        {
            DualParam.RType = Geo::DualPatch<ScalarType>::RTypeOrtho;
        }
        if (std::string(argv[i]) == std::string("no_split_q"))
        {
            DualParam.OnlyNonQuad = true;
        }
        if (std::string(argv[i]) == std::string("border"))
        {
            use_curvature=false;
        }
        if (std::string(argv[i]) == std::string("smoothtangent"))
        {
            smooth_tangent=true;
        }
        if (std::string(argv[i]) == std::string("no_match_sing_val"))
        {
            single_sing_cond = true;
            match_sing_cond = false;
        }
    }
    BatchProcess();

    ProjectName = RemoveExtension(PathMesh);

    CreateOutputFolder();

    SaveAll();

    return 0;
}
