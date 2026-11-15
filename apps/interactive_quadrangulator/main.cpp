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

// Interactive quadrangulator: same steps as the command line tool
// (quad_commandline): features, cross field, quadrangulation, with an
// optional triangle remeshing before (same steps as the interactive triangle
// remesher). The features are chosen at the start and are the same for the
// triangle remeshing and the quadrangulation.
// SKETCHES (SPACE + left button): on the triangle views the line is embedded
// in the input as a feature (used by the next remesh/quadrangulation); on the
// quad view the patches touched by the line (plus some rings of neighbour
// patches) are remeshed with the line as a new feature, the border of the
// region kept as it is (QuadRegionRemesher). Hovering a sketch shows an X:
// clicking it removes the sketch (the feature, the mesh stays).

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <GL/glew.h>
#include <UI/orbitball.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <cctype>
#include <cstddef>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include "IO/io_obj.h"
#include "IO/io_stl.h"
#include <GLFW/glfw3.h>
#include <GL_functions/GLDrawMesh.h>
#include <GL_functions/GLPolylinePick.h>
#include <IO/file_path_functions.h>
#include <Space_Query/triangle_grid_3D.h>
#include <box3.h>
#include <chrono>
#include <color_functions.h>
#include <connectivity.h>
#include <mesh_clean.h>
#include <mesh_embed_polyline.h>
#include <normal.h>
#include <point3.h>
#include <quad_region_remesh.h>
#include <quad_region_select.h>
#include <quad_layout_refit.h>
#include <quad_region_history.h>
#include <polyline_fitting.h>
#include <mesh_sketch_remesh.h>
#include <sketch_log.h>
#include <sketch_symmetry.h>
#include <remesh_size.h>
#include <triangular_remesh.h>

typedef double ScalarType;
typedef Geo::Point3<ScalarType> CoordType;

GLFWwindow *window;
std::string PathMesh;
UI::OrbitBall trackball;

ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);

float line_thick = 1.f;
bool showMesh = true;
bool showFeatures = true;
bool showCorners = true;
bool showSketches = true;
// the border of the mesh is shown as a feature only if asked (the pipeline
// marks it as feature, but it is the border)
bool showBorderFeatures = false;
// DEBUG: on the quad view, the triangle features the quadrangulation started
// from (green), to tell where a quad feature comes from
bool showSourceFeatures = false;
// the quads of the last region remesh (shown in light blue on the quad view)
bool showLastRegion = true;
std::vector<bool> LastRegionFace;
std::vector<std::pair<Geo::Point3<double>, Geo::Point3<double>>>
    SourceFeatSegs;
// (the full quadrangulation part; the last region remesh is added to it)
std::vector<std::pair<Geo::Point3<double>, Geo::Point3<double>>>
    SourceFeatSegsFull;

#if defined(_MSC_VER) && (_MSC_VER >= 1900) &&                                 \
    !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

// VIEWS ---------------------------------------------------------------------
enum ViewType { ViewInput = 0, ViewTri = 1, ViewQuad = 2 };
int viewMode = ViewInput;

// INPUT mesh (after the non-manifold split, with the sketches embedded)
bool loadedMesh = false;
std::string LoadMsg; // why the mesh cannot be processed (empty if it can)
std::vector<CoordType> VertPos0;
std::vector<std::vector<int>> Connectivity0;
std::vector<std::vector<bool>> IsFaceEFeature0;
std::vector<bool> IsVertCorner0;
// the input as loaded, target of the sketch reprojection
std::vector<CoordType> VertPosOrig;
std::vector<std::vector<int>> ConnectivityOrig;

// TRIANGLE REMESH result
bool hasTri = false;
std::vector<CoordType> VertPosR;
std::vector<std::vector<int>> ConnectivityR;
std::vector<std::vector<bool>> IsFaceEFeatureR;
std::vector<bool> IsVertCornerR;

// QUAD result
bool hasQuad = false;
std::vector<CoordType> QuadVert;
std::vector<std::vector<int>> QuadFaces;
std::vector<std::pair<int, int>> QuadFeatures;
std::vector<int> QuadPatch; // layout patch of every quad (-1 unknown)
// Per quad, the mask of its vertices that are corners
// of its layout patch (bit j = vertex j), from the quadrangulator; kept along
// QuadPatch so the layout refit uses the real patches
std::vector<int> QuadPatchCorner;
// Which region remesh every quad comes from, kept up
// to date by the region remesher and the layout refit; a sketch on the quads
// remembers the id of its region, its removal remeshes that area
Geo::QuadRegionHistory<ScalarType> RegionHistory;

// DISPLAYED mesh
std::vector<CoordType> VertPos;
std::vector<std::vector<int>> Connectivity;
std::vector<CoordType> FaceNormals, VertNormals;
std::vector<CoordType> VertCol;
std::vector<std::pair<int, int>> RenderEdges; // quad view: polygon edges
std::vector<std::pair<int, int>> RenderVertFeatures;
std::vector<int> RenderVertCorners;
std::vector<std::pair<int, int>> RenderSketchEdges;

Geo::Box3<ScalarType> MeshBox; // of the input: the view does not jump
ScalarType AVGL = 0;
ScalarType BBDiag = 0;
ScalarType MeshArea = 0;

// FEATURES (the same for the triangle remesh and the quadrangulation) -------
bool use_features = true;
ScalarType feature_angle = 60;
bool clean_features = true; // erode/dilate the feature lines
int ero_steps = 3;
int dil_steps = 3;

// TRIANGLE REMESH (optional first step) --------------------------------------
bool use_tri_remesh = false;
Geo::TriRemParam<ScalarType> RemP;
ScalarType tri_ratio = 1; // scale of the automatic size
ScalarType TriMinSize = 0.005, TriMaxSize = 0.02;
float max_error_perc = 1.f;
#define TRI_REMESH_HARD_MAX_FACES 5000000

// QUADRANGULATION -------------------------------------------------------------
ScalarType QuadSizeMult = 1;           // the -size option
ScalarType QuadMinSize = 0.005, QuadMaxSize = 0.05;
int quad_smooth_steps = 50;            // the -smooth option
bool reliable_smooth = true;           // -fastsmooth = false
float cc_ability = 0.25f;
int num_samples = 6000;
// EXACT FEATURES: only the quad edges on the triangle feature lines are
// features. The quadrangulation marks as feature the whole patch side with a
// feature edge (an open line prolonged to the next corner): on, the side
// samples are marked exactly (exact_feature_samples) and the quad feature
// edges still off the lines (e.g. past the end of an open line) trimmed
bool exact_features = true;
// absolute quad edge of the last full quadrangulation: the region remeshes
// use it, so the size does not drift where a region is remeshed again
ScalarType AbsQuadEdge = -1;
// The sketch only selects the region, the line is not added
// (the remesh should give back about the same quads)
bool debug_region_no_line = false;

// Save region_tri.obj, region_remesh.obj, region_field_tri.obj,
// region_features.obj, region_quad.obj
bool debug_region_save = false;
// After a region remesh the side subdivisions of the
// whole layout solved again (the old patches can change too), the changed
// patches tessellated again; rings of patches allowed to change (-1 = all)
bool region_refit_layout = true;
// Refit keeping every changed side (no cut to the region)
bool debug_refit_no_cut = false;
int region_refit_rings = -1;
// SKETCH ON THE QUADS: true = the line becomes a feature of the input and
// everything is quadrangulated again (as a sketch on the input); false = only
// the region around it is remeshed (QuadRegionRemesher)
bool requad_all = false;
// smoothing of the seam after a region remesh (the region plus rings of
// quads around it, rings and steps the defaults of QuadRegionRemeshParam)
bool seam_smooth = true;

// REGION REMESH (sketch on the quads)
// The quads within a distance from the
// line along the surface, the border regularized by dilate + erode; it does
// not depend on where the line falls in the layout patches
// Whole layout patches, for the layout refit
enum RegionMode { RegionRings = 0, RegionDistance = 1, RegionPatches = 2 };
int region_mode = RegionPatches;
int region_patch_rings = 0;
float region_radius = 3.f; // in quad edges
int region_face_rings = 3;  // rings of quads (RegionRings)

std::vector<std::string> StatLines;

// SKETCHES --------------------------------------------------------------------
GLDraw::GLPolylinePick<ScalarType> LinePick;
bool SpaceDown = false;
bool LeftDown = false;
Geo::EmbedPolylineParam<ScalarType> EmbedP;
// The drawn line is cleaned, closed if almost a loop
// and smoothed once (SketchPolyline); the embedding does not smooth it again
Geo::SketchPolylineParam<ScalarType> FitP;
// The fitting window in quad edges
float fit_window_quads = 3.f;
// Fitting mode (0 B-spline, 1 local Bezier) and the
// distance of the control points in quad edges
int fit_mode = 0;
float fit_spacing_quads = 10.f;
// DEBUG (fitting preview): a drawn line is only processed and shown (the
// line as drawn in blue, the fitted one in orange), nothing is embedded
bool debug_fit_only = false;
std::vector<CoordType> PreviewDrawn, PreviewFitted;
// DEBUG (embedding check): the fitted line is embedded in the input (as a
// sketch on the triangles) but nothing is quadrangulated; shown: the fitted
// line (orange), the embedded one (the sketch edges) and in red the
// vertices where the embedded line turns more than the corner angle of the
// quadrangulator (each one would become a patch corner)
bool debug_embed_only = false;
// DEBUG, the line only embedded (the
// patch around it not remeshed)
bool debug_sketch_no_remesh = false;
std::vector<CoordType> DebugSharpPts;
// the embedding parameters
Geo::EmbedPolylineParam<ScalarType> EmbedNoSmooth() { return EmbedP; }
typedef std::pair<CoordType, CoordType> PosEdge;
PosEdge MakePosEdge(const CoordType &P0, const CoordType &P1) {
  return PosEdge(std::min(P0, P1), std::max(P0, P1));
}

struct Sketch {
  bool on_quad = false;           // false: feature of the input
  std::vector<PosEdge> Segs;      // its feature edges (positions)
  std::set<CoordType> CornersPos; // corners it added (input only)
  // quads only: the line as embedded in the region, its feature edges are
  // the quad features along it (Segs, updated after every region remesh)
  std::vector<PosEdge> LineSegs;
  // Quads only, the id of the region remeshed for it
  int RegionId = -1;
  // Quads only, the line as drawn (moved to the input
  // as it is when switching to "quadrangulate all again")
  std::vector<CoordType> Input;
  // screen segments (ImGui coordinates), updated at every frame
  std::vector<std::pair<ImVec2, ImVec2>> Screen;
  // the sketches of one drawn line (a line and its mirror) share it, so they
  // are removed together; -1 = alone
  int Group = -1;
};
std::vector<Sketch> Sketches;
// the group given to the sketches being added, and the next free one
int NewSketchGroup = -1;
int NumSketchGroups = 0;
std::string LineMsg;

// hover X
int HoverSketch = -1;
ImVec2 HoverAnchor;      // point of the sketch where the hover started
const float XRadius = 9; // pixels
const float HoverDist = 8;
ImVec2 XPos() { return ImVec2(HoverAnchor.x + 14, HoverAnchor.y - 14); }

// ---------------------------------------------------------------------------
// UTILITIES
// ---------------------------------------------------------------------------

ScalarType PointSegDist(const CoordType &P, const CoordType &A,
                        const CoordType &B) {
  CoordType D = B - A;
  ScalarType L2 = Geo::Dot(D, D);
  ScalarType t = 0;
  if (L2 > 0)
    t = Geo::Dot(P - A, D) / L2;
  t = std::max((ScalarType)0, std::min((ScalarType)1, t));
  return (P - (A + D * t)).Norm();
}

ScalarType PointSegsDist(const CoordType &P, const std::vector<PosEdge> &Segs) {
  ScalarType Best = std::numeric_limits<ScalarType>::max();
  for (size_t i = 0; i < Segs.size(); i++)
    Best = std::min(Best, PointSegDist(P, Segs[i].first, Segs[i].second));
  return Best;
}

float ScreenSegDist(const ImVec2 &P, const ImVec2 &A, const ImVec2 &B) {
  float dx = B.x - A.x, dy = B.y - A.y;
  float L2 = dx * dx + dy * dy;
  float t = 0;
  if (L2 > 0)
    t = ((P.x - A.x) * dx + (P.y - A.y) * dy) / L2;
  t = std::max(0.f, std::min(1.f, t));
  float qx = A.x + dx * t - P.x, qy = A.y + dy * t - P.y;
  return std::sqrt(qx * qx + qy * qy);
}

ImVec2 ScreenSegClosest(const ImVec2 &P, const ImVec2 &A, const ImVec2 &B) {
  float dx = B.x - A.x, dy = B.y - A.y;
  float L2 = dx * dx + dy * dy;
  float t = 0;
  if (L2 > 0)
    t = ((P.x - A.x) * dx + (P.y - A.y) * dy) / L2;
  t = std::max(0.f, std::min(1.f, t));
  return ImVec2(A.x + dx * t, A.y + dy * t);
}

// feature edges as vertex pairs from the per face flags
void FeaturePairs(const std::vector<std::vector<int>> &Faces,
                  const std::vector<std::vector<bool>> &IsFaceEFeature,
                  std::vector<std::pair<int, int>> &Pairs) {
  Pairs.clear();
  if (IsFaceEFeature.size() != Faces.size())
    return;
  Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
      Faces, IsFaceEFeature, Pairs);
}

// polygon edges (unique) of a mesh
void PolyEdges(const std::vector<std::vector<int>> &Faces,
               std::vector<std::pair<int, int>> &Edges) {
  std::set<std::pair<int, int>> ES;
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % Faces[i].size()];
      ES.insert(std::make_pair(std::min(V0, V1), std::max(V0, V1)));
    }
  Edges.assign(ES.begin(), ES.end());
}

// face normals (Newell) and vertex normals of a polygon mesh
void PolyNormals(const std::vector<CoordType> &Vert,
                 const std::vector<std::vector<int>> &Faces,
                 std::vector<CoordType> &FN, std::vector<CoordType> &VN) {
  FN.assign(Faces.size(), CoordType(0, 0, 0));
  VN.assign(Vert.size(), CoordType(0, 0, 0));
  for (size_t i = 0; i < Faces.size(); i++) {
    CoordType N(0, 0, 0);
    for (size_t j = 0; j < Faces[i].size(); j++) {
      const CoordType &A = Vert[Faces[i][j]];
      const CoordType &B = Vert[Faces[i][(j + 1) % Faces[i].size()]];
      N.X += (A.Y - B.Y) * (A.Z + B.Z);
      N.Y += (A.Z - B.Z) * (A.X + B.X);
      N.Z += (A.X - B.X) * (A.Y + B.Y);
    }
    for (size_t j = 0; j < Faces[i].size(); j++)
      VN[Faces[i][j]] += N;
    if (N.Norm() > 0)
      N.Normalize();
    FN[i] = N;
  }
  for (size_t i = 0; i < VN.size(); i++)
    if (VN[i].Norm() > 0)
      VN[i].Normalize();
}

// ---------------------------------------------------------------------------
// VIEW
// ---------------------------------------------------------------------------

void UpdateRenderSketches() {
  RenderSketchEdges.clear();
  std::set<PosEdge> SketchPos;
  for (size_t i = 0; i < Sketches.size(); i++) {
    if (Sketches[i].on_quad != (viewMode == ViewQuad))
      continue;
    if ((!Sketches[i].on_quad) && (viewMode != ViewInput))
      continue;
    SketchPos.insert(Sketches[i].Segs.begin(), Sketches[i].Segs.end());
  }
  if (SketchPos.empty())
    return;
  std::vector<std::pair<int, int>> Edges;
  PolyEdges(Connectivity, Edges);
  for (size_t i = 0; i < Edges.size(); i++)
    if (SketchPos.count(MakePosEdge(VertPos[Edges[i].first],
                                    VertPos[Edges[i].second])) > 0)
      RenderSketchEdges.push_back(Edges[i]);
}

void ShowView(int view) {
  if ((view == ViewTri) && !hasTri)
    view = ViewInput;
  if ((view == ViewQuad) && !hasQuad)
    view = ViewInput;
  viewMode = view;
  RenderEdges.clear();
  RenderVertFeatures.clear();
  RenderVertCorners.clear();
  if (viewMode == ViewQuad) {
    VertPos = QuadVert;
    Connectivity = QuadFaces;
    PolyNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    PolyEdges(Connectivity, RenderEdges);
    RenderVertFeatures = QuadFeatures;
  } else if (viewMode == ViewTri) {
    VertPos = VertPosR;
    Connectivity = ConnectivityR;
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    FeaturePairs(Connectivity, IsFaceEFeatureR, RenderVertFeatures);
    if (IsVertCornerR.size() == VertPos.size())
      Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(
          IsVertCornerR, RenderVertCorners);
  } else {
    VertPos = VertPos0;
    Connectivity = Connectivity0;
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    FeaturePairs(Connectivity, IsFaceEFeature0, RenderVertFeatures);
    if (IsVertCorner0.size() == VertPos.size())
      Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(
          IsVertCorner0, RenderVertCorners);
  }
  // the border edges are not shown as features
  if (!showBorderFeatures) {
    std::map<std::pair<int, int>, int> EdgeCount;
    for (size_t i = 0; i < Connectivity.size(); i++)
      for (size_t j = 0; j < Connectivity[i].size(); j++) {
        int V0 = Connectivity[i][j];
        int V1 = Connectivity[i][(j + 1) % Connectivity[i].size()];
        EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
      }
    std::vector<std::pair<int, int>> Kept;
    for (size_t i = 0; i < RenderVertFeatures.size(); i++) {
      int V0 = RenderVertFeatures[i].first, V1 = RenderVertFeatures[i].second;
      std::pair<int, int> Key(std::min(V0, V1), std::max(V0, V1));
      if ((EdgeCount.count(Key) > 0) && (EdgeCount[Key] == 1))
        continue;
      Kept.push_back(RenderVertFeatures[i]);
    }
    RenderVertFeatures = Kept;
  }
  VertCol.assign(VertPos.size(), CoordType(0.9, 0.9, 0.9));
  UpdateRenderSketches();
  HoverSketch = -1;
}

// ---------------------------------------------------------------------------
// FEATURES OF THE INPUT: detection + the sketches
// ---------------------------------------------------------------------------

// The detection of the panel on any mesh
// (as UpdateFeatureCorners did on the input)
void DetectFeatures(std::vector<CoordType> &Pos,
                    std::vector<std::vector<int>> &Faces,
                    std::vector<std::vector<bool>> &IsFaceEFeature,
                    std::vector<bool> &IsVertCorner) {
  if (use_features && (feature_angle > 0) && (feature_angle < 180)) {
    if (clean_features)
      Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
          Pos, Faces, feature_angle, IsFaceEFeature, IsVertCorner, ero_steps,
          dil_steps);
    else
      Geo::MeshFeatures<ScalarType>::FeatureCorners(
          Pos, Faces, feature_angle, IsFaceEFeature, IsVertCorner);
  } else {
    IsVertCorner.assign(Pos.size(), false);
    IsFaceEFeature.assign(Faces.size(), std::vector<bool>(3, false));
  }
}

void UpdateFeatureCorners() {
  if (!loadedMesh)
    return;
  DetectFeatures(VertPos0, Connectivity0, IsFaceEFeature0, IsVertCorner0);
  // SKETCHES on the input: always features (and their corners)
  std::set<PosEdge> SketchPos;
  std::set<CoordType> SketchCorners;
  for (size_t i = 0; i < Sketches.size(); i++) {
    if (Sketches[i].on_quad)
      continue;
    SketchPos.insert(Sketches[i].Segs.begin(), Sketches[i].Segs.end());
    SketchCorners.insert(Sketches[i].CornersPos.begin(),
                         Sketches[i].CornersPos.end());
  }
  for (size_t i = 0; i < Connectivity0.size(); i++)
    for (size_t j = 0; j < 3; j++) {
      int V0 = Connectivity0[i][j], V1 = Connectivity0[i][(j + 1) % 3];
      if (SketchPos.count(MakePosEdge(VertPos0[V0], VertPos0[V1])) > 0)
        IsFaceEFeature0[i][j] = true;
    }
  for (size_t i = 0; i < VertPos0.size(); i++)
    if (SketchCorners.count(VertPos0[i]) > 0)
      IsVertCorner0[i] = true;
  if (viewMode == ViewInput)
    ShowView(ViewInput);
}

// ---------------------------------------------------------------------------
// LOAD
// ---------------------------------------------------------------------------

bool IsSTLFile(const std::string &path) {
  if (path.size() < 4)
    return false;
  std::string ext = path.substr(path.size() - 4);
  for (size_t i = 0; i < ext.size(); i++)
    ext[i] = (char)std::tolower((unsigned char)ext[i]);
  return (ext == std::string(".stl"));
}

// SKETCH LOG: every action recorded (with its settings) before it runs, to
// replay a session from the loaded mesh (e.g. up to a crash)
Geo::SketchLog<ScalarType> SLog;
// symmetric sketches
Geo::SketchSymmetry<ScalarType> Sym;
// the quads are only the half on the positive side of HalfSym (the plane as
// it was at the quadrangulation), shown and saved with their mirror
bool QuadIsHalf = false;
Geo::SketchSymmetry<ScalarType> HalfSym;
std::string SketchLogPath;

void LoadMesh(const std::string &path) {
  loadedMesh = IsSTLFile(path) ? LoadSTL(path, VertPos0, Connectivity0, true)
                               : LoadOBJ(path, VertPos0, Connectivity0, true);
  if (!loadedMesh || Connectivity0.empty()) {
    loadedMesh = false;
    std::cout << "Mesh not loaded" << std::endl;
    return;
  }
  if (!IsAllTris(Connectivity0)) {
    std::vector<CoordType> VertPosTris;
    std::vector<std::vector<int>> ConnectivityTris;
    SplitIntoTris(VertPos0, Connectivity0, VertPosTris, ConnectivityTris);
    VertPos0 = VertPosTris;
    Connectivity0 = ConnectivityTris;
  }
  std::cout << "Loaded " << Connectivity0.size() << " faces " << VertPos0.size()
            << " vertices" << std::endl;
  LoadMsg.clear();
  Geo::CleanMesh<ScalarType>::SplitNonManifold(VertPos0, Connectivity0);
  if (!MatchTriRemeshPreconditions(VertPos0, Connectivity0)) {
    bool ManifE = EdgeManifold<ScalarType>(VertPos0, Connectivity0);
    bool ManifV = VertManifold<ScalarType>(VertPos0, Connectivity0);
    LoadMsg = std::string("not manifold (edge manifold ") +
              std::to_string((int)ManifE) + ", vertex manifold " +
              std::to_string((int)ManifV) + ")";
    std::cout << "Warning: " << LoadMsg << std::endl;
  }
  VertPosOrig = VertPos0;
  ConnectivityOrig = Connectivity0;
  Sketches.clear();

  MeshBox.Init(VertPos0);
  AVGL = AvgEdgeLen(VertPos0, Connectivity0);
  BBDiag = MeshBox.Diag();
  MeshArea = Area(VertPos0, Connectivity0);
  hasTri = false;
  hasQuad = false;
  StatLines.clear();
  viewMode = ViewInput;
  UpdateFeatureCorners();
  ShowView(ViewInput);
  // a new sketch log, next to the mesh
  SketchLogPath = path + ".sketch_log.txt";
  SLog.Start(SketchLogPath);
}

// ---------------------------------------------------------------------------
// TRIANGLE REMESH (as the interactive triangle remesher, automatic size)
// ---------------------------------------------------------------------------

// The remesh parameters of the panel
// (as RunTriRemesh sets them), shared by the patch remesh of the sketches
Geo::TriRemParam<ScalarType> PanelRemParam() {
  Geo::TriRemParam<ScalarType> RP = RemP;
  RP.feature_angle = use_features ? feature_angle : 180;
  RP.max_error = BBDiag * max_error_perc / 100;
  RP.writeMSG = false;
  return RP;
}

bool RunTriRemesh() {
  if (!loadedMesh)
    return false;
  if (!LoadMsg.empty()) {
    StatLines.push_back(std::string("cannot remesh: ") + LoadMsg);
    return false;
  }
  auto start = std::chrono::high_resolution_clock::now();
  Geo::TriRemParam<ScalarType> RP = PanelRemParam();
  VertPosR = VertPos0;
  ConnectivityR = Connectivity0;
  RP.IsFaceEFeature = IsFaceEFeature0;
  RP.IsVertCorner = IsVertCorner0;

  // automatic size
  RP.VertSize.clear();
  ScalarType MinAggrSize = BBDiag * TriMinSize;
  ScalarType MaxAggrSize = BBDiag * TriMaxSize;
  // One size for the whole mesh
  RP.adaptive = false;
  // The size is a property of the surface
  // as loaded (its own features), not of the embedding of the lines (the
  // refinement under them made it smaller)
  std::vector<CoordType> SizePos = VertPosOrig;
  std::vector<std::vector<int>> SizeFaces = ConnectivityOrig;
  std::vector<std::vector<bool>> SizeFE;
  std::vector<bool> SizeCorner;
  DetectFeatures(SizePos, SizeFaces, SizeFE, SizeCorner);
  ScalarType CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeAggregated(
      SizePos, SizeFaces, SizeFE, MinAggrSize, MaxAggrSize, NULL);
  if (!(CurrSize > 0)) {
    StatLines.push_back("tri remesh: invalid remesh size");
    return false;
  }
  CurrSize /= 2;
  RP.TargetL = tri_ratio * CurrSize;

  RP.max_faces = TRI_REMESH_HARD_MAX_FACES;
  RP.check_manifold = true;
  int res = Geo::TriRemesh<ScalarType>(VertPosR, ConnectivityR, RP);
  if (res >= 0) {
    bool FeaturesOK = (RP.IsFaceEFeature.size() == ConnectivityR.size()) &&
                      (RP.IsVertCorner.size() == VertPosR.size());
    Geo::CleanMesh<ScalarType>::CollapseBadShapeFaces(
        VertPosR, ConnectivityR, 0.05, FeaturesOK ? &RP.IsFaceEFeature : NULL,
        FeaturesOK ? &RP.IsVertCorner : NULL);
    FeaturesOK = FeaturesOK &&
                 (RP.IsFaceEFeature.size() == ConnectivityR.size()) &&
                 (RP.IsVertCorner.size() == VertPosR.size());
    std::vector<bool> IsFixedV(VertPosR.size(), false);
    if (FeaturesOK) {
      std::vector<std::pair<int, int>> FeatPairs;
      FeaturePairs(ConnectivityR, RP.IsFaceEFeature, FeatPairs);
      for (size_t i = 0; i < FeatPairs.size(); i++) {
        IsFixedV[FeatPairs[i].first] = true;
        IsFixedV[FeatPairs[i].second] = true;
      }
      for (size_t i = 0; i < VertPosR.size(); i++)
        if (RP.IsVertCorner[i])
          IsFixedV[i] = true;
    }
    Geo::CleanMesh<ScalarType>::SolveFolds(VertPosR, ConnectivityR, -0.99,
                                           FeaturesOK ? &IsFixedV : NULL);
  }
  if (res < 0) {
    StatLines.push_back("tri remesh: stopped (see the console)");
    hasTri = false;
    return false;
  }
  IsFaceEFeatureR = RP.IsFaceEFeature;
  IsVertCornerR = RP.IsVertCorner;
  if (IsFaceEFeatureR.size() != ConnectivityR.size())
    IsFaceEFeatureR.assign(ConnectivityR.size(), std::vector<bool>(3, false));
  if (IsVertCornerR.size() != VertPosR.size())
    IsVertCornerR.assign(VertPosR.size(), false);
  hasTri = true;
  auto end = std::chrono::high_resolution_clock::now();
  char buf[256];
  snprintf(buf, sizeof(buf), "tri remesh: %.2f s, %zu faces",
           std::chrono::duration<double>(end - start).count(),
           ConnectivityR.size());
  StatLines.push_back(buf);
  return true;
}

// ---------------------------------------------------------------------------
// QUADRANGULATION (as the command line tool)
// ---------------------------------------------------------------------------

ScalarType AvgQuadEdge(); // (below)

void DoQuadrangulate() {
  if (!loadedMesh)
    return;
  StatLines.clear();
  if (!LoadMsg.empty()) {
    StatLines.push_back(std::string("cannot quadrangulate: ") + LoadMsg);
    return;
  }
  auto start = std::chrono::high_resolution_clock::now();

  // 1. the triangle mesh: remeshed or the input, with its features
  std::vector<CoordType> TriV;
  std::vector<std::vector<int>> TriF;
  std::vector<std::pair<int, int>> TriFeatures;
  // (the sketches are remeshed in the input already: see SketchOnTriangles)
  if (use_tri_remesh) {
    if (!RunTriRemesh())
      return;
    TriV = VertPosR;
    TriF = ConnectivityR;
    FeaturePairs(TriF, IsFaceEFeatureR, TriFeatures);
  } else {
    TriV = VertPos0;
    TriF = Connectivity0;
    FeaturePairs(TriF, IsFaceEFeature0, TriFeatures);
  }
  // with the symmetry only the half on the positive side of the plane (the
  // plane fixed from now on for these quads)
  QuadIsHalf = Sym.Enabled;
  if (QuadIsHalf) {
    HalfSym = Sym;
    HalfSym.CutHalf(TriV, TriF, TriFeatures, MeshBox.Center(), BBDiag);
  }

  // DEBUG: the triangle features the quadrangulation starts from
  SourceFeatSegs.clear();
  for (size_t i = 0; i < TriFeatures.size(); i++)
    SourceFeatSegs.push_back(
        std::make_pair(TriV[TriFeatures[i].first], TriV[TriFeatures[i].second]));
  SourceFeatSegsFull = SourceFeatSegs;

  // the true surface for the region remeshes (before the pipeline moves TriV)

  // 2. cross field (BatchProcessField)
  std::vector<Field::CrossF<ScalarType>> FaceCurv;
  Geo::QuadRegionRemesher<ScalarType>::ComputeField(TriV, TriF, TriFeatures,
                                                    FaceCurv, 3);

  // 3. quad size (AutoQuadEdgeSize: always with 30 degrees features)
  Geo::Box3<ScalarType> TriBox;
  TriBox.Init(TriV);
  std::vector<std::pair<int, int>> TriFeatures_30;
  std::vector<int> TriCorners_30;
  if (use_features)
    Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
        TriV, TriF, 30, TriFeatures_30, TriCorners_30, ero_steps, dil_steps);
  ScalarType CurrSize = Geo::RemeshSize<ScalarType>::RemeshSizeAggregated(
      TriV, TriF, TriFeatures_30, TriBox.Diag() * QuadMinSize,
      TriBox.Diag() * QuadMaxSize, NULL);
  ScalarType QuadSize = (CurrSize / TriBox.Diag()) / 2 * QuadSizeMult;
  // (the pipeline uses MeshBox.Diag() * QuadSize, MeshBox of this mesh)
  AbsQuadEdge = TriBox.Diag() * QuadSize;

  // the triangle features as segments (before the pipeline moves TriV)
  std::vector<PosEdge> TriFeatSegs;
  for (size_t i = 0; i < TriFeatures.size(); i++)
    TriFeatSegs.push_back(
        MakePosEdge(TriV[TriFeatures[i].first], TriV[TriFeatures[i].second]));

  // 4. quadrangulation
  std::vector<CoordType> NewQuadVert;
  std::vector<std::vector<int>> NewQuadFaces;
  bool problematic = false;
  {
    MeshQuadrangulator<ScalarType> MQuad(TriV, TriF, FaceCurv, NewQuadVert,
                                         NewQuadFaces);
    MQuad.Init();
    if (!MQuad.MatchPreconditions()) {
      StatLines.push_back("quad: triangle mesh not manifold");
      return;
    }
    MQuad.QuadSize = QuadSize;
    MQuad.QuadSmoothSteps = quad_smooth_steps;
    MQuad.CCAbility = cc_ability;
    MQuad.num_samples = num_samples;
    MQuad.reliable_smooth = reliable_smooth;
    MQuad.exact_feature_samples = exact_features;
    MQuad.ExtractQuad(TriFeatures);
    MQuad.GetQuadFeatures(QuadFeatures);
    MQuad.GetQuadPatch(QuadPatch);
    MQuad.GetQuadPatchCorner(QuadPatchCorner);
    // The layouts of the full quadrangulation too, to
    // compare with the ones of a region
    if (debug_region_save) {
      std::vector<CoordType> PrimalV, DualV;
      std::vector<std::vector<int>> PrimalF, DualF;
      MQuad.GetPrimalDualMesh(PrimalV, PrimalF, DualV, DualF);
      WriteOBJ("full_primal_layout.obj", PrimalV, PrimalF);
      WriteOBJ("full_dual_layout.obj", DualV, DualF);
    }
    problematic = MQuad.HasProblematicPatches();
  }
  QuadVert = NewQuadVert;
  QuadFaces = NewQuadFaces;
  RegionHistory.Reset(QuadFaces.size());
  hasQuad = true;
  LastRegionFace.clear(); // no region remesh on this quad mesh yet
  // no prolongation of the open feature lines
  if (exact_features) // (cleanup: trim and exact features together)
    Geo::QuadRegionRemesher<ScalarType>::TrimFeaturesToLines(
        QuadVert, QuadFeatures, TriFeatSegs, AvgQuadEdge() * (ScalarType)0.25);

  // the old quad sketches belong to the old quad mesh
  std::vector<Sketch> Kept;
  for (size_t i = 0; i < Sketches.size(); i++)
    if (!Sketches[i].on_quad)
      Kept.push_back(Sketches[i]);
  Sketches = Kept;

  auto end = std::chrono::high_resolution_clock::now();
  char buf[256];
  snprintf(buf, sizeof(buf), "quad: %.2f s, %zu quads, %zu feature edges%s",
           std::chrono::duration<double>(end - start).count(),
           QuadFaces.size(), QuadFeatures.size(),
           problematic ? " (problematic patches)" : "");
  StatLines.push_back(buf);
  std::vector<ScalarType> FaceAR;
  Geo::AspectRatioScaledJacobian(QuadVert, QuadFaces, FaceAR);
  if (!FaceAR.empty()) {
    ScalarType ARMin = *std::min_element(FaceAR.begin(), FaceAR.end());
    ScalarType ARAvg = 0;
    for (size_t i = 0; i < FaceAR.size(); i++)
      ARAvg += FaceAR[i];
    ARAvg /= FaceAR.size();
    snprintf(buf, sizeof(buf), "scaled jacobian: avg %.3f min %.3f", ARAvg,
             ARMin);
    StatLines.push_back(buf);
  }
  for (size_t i = 0; i < StatLines.size(); i++)
    std::cout << StatLines[i] << std::endl;
  ShowView(ViewQuad);
}

// ---------------------------------------------------------------------------
// SKETCH ON THE TRIANGLES: embedded in the input as a feature
// ---------------------------------------------------------------------------

// EParam, if given, replaces the panel parameters of the
// embedding (e.g. no smoothing for a line already smoothed)
void SketchOnTriangles(const std::vector<CoordType> &Line,
                       const Geo::EmbedPolylineParam<ScalarType> *EParam =
                           NULL) {
  size_t NumV0 = VertPos0.size();
  std::set<CoordType> WasCornerPos;
  for (size_t i = 0; i < IsVertCorner0.size(); i++)
    if (IsVertCorner0[i])
      WasCornerPos.insert(VertPos0[i]);
  // Embedded, then the patch changed by
  // the embedding remeshed with its border fixed
  // the edges of the input sketches (pairs of points), updated by Apply
  // (split by the embedding, replaced inside the remeshed patch), the new
  // line appended
  std::vector<int> InputSketch;
  std::vector<std::vector<PosEdge>> SketchSegs;
  for (size_t i = 0; i < Sketches.size(); i++)
    if (!Sketches[i].on_quad) {
      InputSketch.push_back((int)i);
      SketchSegs.push_back(Sketches[i].Segs);
    }
  std::string RemMsg;
  if (!Geo::MeshSketchRemesh<ScalarType>::Apply(
          VertPos0, Connectivity0, IsFaceEFeature0, IsVertCorner0, VertPosOrig,
          ConnectivityOrig, Line, (EParam != NULL) ? *EParam : EmbedNoSmooth(),
          PanelRemParam(), !debug_sketch_no_remesh, SketchSegs, RemMsg)) {
    LineMsg = "line not embedded (too short or off the mesh)";
    return;
  }
  for (size_t i = 0; i < InputSketch.size(); i++)
    Sketches[InputSketch[i]].Segs = SketchSegs[i];
  Sketch S;
  S.on_quad = false;
  S.Input = Line;
  S.Group = NewSketchGroup;
  S.Segs = SketchSegs.back();
  for (size_t i = 0; i < IsVertCorner0.size(); i++)
    if (IsVertCorner0[i] && (WasCornerPos.count(VertPos0[i]) == 0))
      S.CornersPos.insert(VertPos0[i]);
  Sketches.push_back(S);
  char buf[128];
  snprintf(buf, sizeof(buf), "line embedded: %zu edges, %d vertices more%s",
           S.Segs.size(), (int)VertPos0.size() - (int)NumV0,
           debug_sketch_no_remesh ? " (DEBUG: not remeshed)"
                                  : ", patch remeshed");
  LineMsg = buf;
  if (!RemMsg.empty())
    LineMsg += std::string(" (") + RemMsg + ")";
  AVGL = AvgEdgeLen(VertPos0, Connectivity0);
  MeshArea = Area(VertPos0, Connectivity0);
  // the triangle remesh / quads no longer have this feature: redo them
  UpdateFeatureCorners();
  ShowView(ViewInput);
}

// ---------------------------------------------------------------------------
// SKETCH ON THE QUADS: region remesh with the line as a new feature
// ---------------------------------------------------------------------------

ScalarType AvgQuadEdge() {
  ScalarType L = 0;
  int Num = 0;
  for (size_t i = 0; i < QuadFaces.size(); i++)
    for (size_t j = 0; j < QuadFaces[i].size(); j++) {
      L += (QuadVert[QuadFaces[i][j]] -
            QuadVert[QuadFaces[i][(j + 1) % QuadFaces[i].size()]])
               .Norm();
      Num++;
    }
  return (Num > 0) ? L / Num : 0;
}

// the region parameters from the panel
typename Geo::QuadRegionSelector<ScalarType>::Param RegionParam() {
  typename Geo::QuadRegionSelector<ScalarType>::Param P;
  P.Mode = (typename Geo::QuadRegionSelector<ScalarType>::GrowMode)region_mode;
  P.Radius = region_radius;
  P.FaceRings = region_face_rings;
  P.PatchRings = region_patch_rings;
  return P;
}

// the quad feature edges of a sketch (the whole feature line along its
// embedded line)
void QuadSketchEdges(const Sketch &S, std::vector<int> &FeatIdx) {
  Geo::QuadFeatureLines<ScalarType>::EdgesOfLine(
      QuadVert, QuadFeatures, S.LineSegs, AvgQuadEdge() * (ScalarType)0.25,
      FeatIdx);
}

// Segs of all the quad sketches from the current quad features
void UpdateQuadSketchSegs() {
  for (size_t s = 0; s < Sketches.size(); s++) {
    if (!Sketches[s].on_quad)
      continue;
    std::vector<int> FeatIdx;
    QuadSketchEdges(Sketches[s], FeatIdx);
    Sketches[s].Segs.clear();
    for (size_t i = 0; i < FeatIdx.size(); i++)
      Sketches[s].Segs.push_back(
          MakePosEdge(QuadVert[QuadFeatures[FeatIdx[i]].first],
                      QuadVert[QuadFeatures[FeatIdx[i]].second]));
  }
}

// remesh of the region grown around the touched quads, with the new lines
// as features; if it fails the region is grown (up to two more times).
// EmbeddedLines: the lines as embedded (one per new line).
// Irregular vertices (valence != 4) inside a set of faces
// (all their faces in the set, not on a feature)
int InnerIrregular(const std::vector<std::vector<int>> &Faces,
                   const std::vector<int> &FaceSet,
                   const std::vector<std::pair<int, int>> &Features) {
  std::vector<bool> In(Faces.size(), false);
  for (size_t i = 0; i < FaceSet.size(); i++)
    if ((FaceSet[i] >= 0) && (FaceSet[i] < (int)Faces.size()))
      In[FaceSet[i]] = true;
  std::map<int, std::pair<int, int>> VC; // vertex -> (faces in, faces)
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      VC[Faces[i][j]].second++;
      if (In[i])
        VC[Faces[i][j]].first++;
    }
  std::set<int> FeatV;
  for (size_t i = 0; i < Features.size(); i++) {
    FeatV.insert(Features[i].first);
    FeatV.insert(Features[i].second);
  }
  // (a vertex on the mesh border has faces only on one side: its count is
  // not a valence; such vertices are left out by the border of the region)
  int Num = 0;
  std::map<int, std::pair<int, int>>::iterator It;
  for (It = VC.begin(); It != VC.end(); It++)
    if (((*It).second.first == (*It).second.second) &&
        ((*It).second.second != 4) && (FeatV.count((*It).first) == 0))
      Num++;
  return Num;
}

// The irregular vertices of InnerIrregular split in the
// ones on the border of the layout patches (corners of the layout: no
// solve of the subdivisions removes them) and the ones inside a patch (its
// tessellation: opposite sides different)
void IrregularKinds(const std::vector<std::vector<int>> &Faces,
                    const std::vector<int> &FaceSet,
                    const std::vector<std::pair<int, int>> &Features,
                    const std::vector<int> &Patch, int &OnLayout,
                    int &InPatch) {
  OnLayout = InPatch = 0;
  std::vector<bool> In(Faces.size(), false);
  for (size_t i = 0; i < FaceSet.size(); i++)
    if ((FaceSet[i] >= 0) && (FaceSet[i] < (int)Faces.size()))
      In[FaceSet[i]] = true;
  std::map<int, std::pair<int, int>> VC; // vertex -> (faces in, faces)
  std::map<int, std::set<int>> VP;       // vertex -> patches
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      VC[Faces[i][j]].second++;
      if (In[i])
        VC[Faces[i][j]].first++;
      VP[Faces[i][j]].insert((i < Patch.size()) ? Patch[i] : -1);
    }
  std::set<int> FeatV;
  for (size_t i = 0; i < Features.size(); i++) {
    FeatV.insert(Features[i].first);
    FeatV.insert(Features[i].second);
  }
  std::map<int, std::pair<int, int>>::iterator It;
  for (It = VC.begin(); It != VC.end(); It++)
    if (((*It).second.first == (*It).second.second) &&
        ((*It).second.second != 4) && (FeatV.count((*It).first) == 0)) {
      if (VP[(*It).first].size() > 1)
        OnLayout++;
      else
        InPatch++;
    }
}

// ForcedRegion, if given, is the region of the first
// attempt (the area of a sketch being removed); the other attempts grow
// from Touched as usual
bool RemeshQuadRegion(const std::vector<bool> &Touched,
                      const std::vector<std::vector<CoordType>> &Lines,
                      std::vector<std::vector<PosEdge>> &EmbeddedLines,
                      std::string &Msg,
                      const std::vector<int> *ForcedRegion = NULL) {
  Geo::QuadRegionRemeshParam<ScalarType> RParam;
  RParam.EmbedParam = EmbedNoSmooth();
  RParam.QuadSmoothSteps = quad_smooth_steps;
  RParam.CCAbility = cc_ability;
  RParam.NumSamples = num_samples;
  RParam.reliable_smooth = reliable_smooth;
  RParam.trim_features = exact_features;
  // The true surface (the triangle mesh of the full
  // quadrangulation)
  RParam.save_debug = debug_region_save;
  if (!seam_smooth)
    RParam.SeamSmoothSteps = 0;
  // the size chosen at the start, not measured again on the region
  RParam.TargetQuadEdge =
      (AbsQuadEdge > 0) ? AbsQuadEdge : -1; // (cleanup: no region scale)
  char buf[256];
  std::string Reasons;
  // The region as selected, no growing
  for (int attempt = 0; attempt < 1; attempt++) {
    std::vector<int> Region;
    ScalarType EdgeL = (AbsQuadEdge > 0) ? AbsQuadEdge : AvgQuadEdge();
    auto TSel = std::chrono::high_resolution_clock::now();
    // The forced region first
    if ((attempt == 0) && (ForcedRegion != NULL) && !ForcedRegion->empty())
      Region = *ForcedRegion;
    else
      Geo::QuadRegionSelector<ScalarType>::Grow(QuadVert, QuadFaces,
                                                Touched, RegionParam(), EdgeL,
                                                attempt + 1, Region,
                                                &QuadPatch);
    // The selection of the region
    double SelSec = std::chrono::duration<double>(
                        std::chrono::high_resolution_clock::now() - TSel)
                        .count();
    std::cout << "Region selection: " << SelSec << " s, " << Region.size()
              << " quads of " << QuadFaces.size() << " (attempt "
              << attempt + 1 << ")" << std::endl;
    if (Region.empty()) {
      Msg = "empty region";
      return false;
    }
    auto start = std::chrono::high_resolution_clock::now();
    std::set<int> RegionSet(Region.begin(), Region.end());
    std::vector<int> OldPatch = QuadPatch;
    std::vector<int> OldPatchCorner = QuadPatchCorner;
    std::vector<int> NewRegionFaces, NewRegionPatch, UsedRegion;
    std::vector<int> NewRegionPatchCorner;
    // The old mesh, to count its singularities
    std::vector<std::vector<int>> OldFaces = QuadFaces;
    std::vector<std::pair<int, int>> OldFeatures = QuadFeatures;
    bool done = Geo::QuadRegionRemesher<ScalarType>::Remesh(
        QuadVert, QuadFaces, Region, QuadFeatures, Lines, RParam,
        &NewRegionFaces, &NewRegionPatch, &EmbeddedLines, &UsedRegion,
        &RegionHistory,
        &NewRegionPatchCorner);
    // the region may have been completed to be manifold
    RegionSet = std::set<int>(UsedRegion.begin(), UsedRegion.end());
    auto end = std::chrono::high_resolution_clock::now();
    if (!done) {
      std::string Why = Geo::QuadRegionRemesher<ScalarType>::LastMessage();
      std::cout << "Region of " << Region.size()
                << " quads not remeshed, attempt " << attempt + 1
                << " of 3: " << Why << std::endl;
      Reasons += std::string("\n  attempt ") + std::to_string(attempt + 1) +
                 " (" + std::to_string(Region.size()) + " quads): " + Why;
      continue;
    }
    // layout patches: the kept faces keep theirs (same order), the new ones
    // get new ids
    std::vector<int> NewPatch;
    int MaxP = -1;
    for (size_t i = 0; i < OldPatch.size(); i++)
      MaxP = std::max(MaxP, OldPatch[i]);
    if (OldPatch.size() ==
        RegionSet.size() + (QuadFaces.size() - NewRegionFaces.size())) {
      for (size_t i = 0; i < OldPatch.size(); i++)
        if (RegionSet.count((int)i) == 0)
          NewPatch.push_back(OldPatch[i]);
    }
    for (size_t i = 0; i < NewRegionFaces.size(); i++) {
      int P = (i < NewRegionPatch.size()) ? NewRegionPatch[i] : -1;
      NewPatch.push_back((P >= 0) ? (MaxP + 1 + P) : -1);
    }
    QuadPatch = NewPatch;
    // The same for the corner masks (the old patches
    // keep their corners, the new ones those of the region quadrangulation)
    {
      std::vector<int> NewPatchCorner;
      if (OldPatchCorner.size() == OldPatch.size())
        for (size_t i = 0; i < OldPatchCorner.size(); i++)
          if (RegionSet.count((int)i) == 0)
            NewPatchCorner.push_back(OldPatchCorner[i]);
      for (size_t i = 0; i < NewRegionFaces.size(); i++)
        NewPatchCorner.push_back(
            (i < NewRegionPatchCorner.size()) ? NewRegionPatchCorner[i] : 0);
      QuadPatchCorner = NewPatchCorner;
      if (QuadPatchCorner.size() != QuadFaces.size())
        QuadPatchCorner.clear(); // (unknown: the refit guesses them)
    }
    // (updated by the remesher)
    // the quads of this region (the merged mesh: the new ones at the end)
    LastRegionFace.assign(QuadFaces.size(), false);
    for (size_t i = 0; i < NewRegionFaces.size(); i++)
      if (NewRegionFaces[i] < (int)LastRegionFace.size())
        LastRegionFace[NewRegionFaces[i]] = true;
    if (QuadPatch.size() != QuadFaces.size())
      QuadPatch.assign(QuadFaces.size(), -1);
    // The subdivisions of the layout solved again, the
    // region border no longer fixed; the shown region becomes the region plus
    // the patches tessellated again
    std::string RefitMsg;
    // With a region cut through the old
    // patches these are not 4 sided any more and stay fixed: no propagation
    if (region_refit_layout && (region_mode == RegionPatches)) {
      std::vector<int> AllF(QuadFaces.size());
      for (size_t i = 0; i < AllF.size(); i++)
        AllF[i] = (int)i;
      int SingAll0 = InnerIrregular(QuadFaces, AllF, QuadFeatures);
      // Of the region, on the layout / inside patches
      int RLay0, RIn0;
      IrregularKinds(QuadFaces, NewRegionFaces, QuadFeatures, QuadPatch, RLay0,
                     RIn0);
      Geo::QuadLayoutRefitParam<ScalarType> LP;
      LP.PatchRings = region_refit_rings;
      // the target edge of the full quadrangulation
      LP.TargetEdge = AbsQuadEdge;
      LP.keep_connected_only = !debug_refit_no_cut;
      // Its steps and mode
      LP.SmoothSteps = quad_smooth_steps;
      LP.reliable_smooth = reliable_smooth;
      std::vector<int> OldFaceOf;
      auto TR = std::chrono::high_resolution_clock::now();
      bool refit = Geo::QuadLayoutRefit<ScalarType>::Refit(
          QuadVert, QuadFaces, QuadPatch, QuadFeatures, NewRegionFaces, LP,
          &OldFaceOf, &RegionHistory,
          &QuadPatchCorner);
      double RSec = std::chrono::duration<double>(
                        std::chrono::high_resolution_clock::now() - TR)
                        .count();
      if (refit) {
        // (updated by the refit)
        std::vector<bool> NewLast(QuadFaces.size(), false);
        std::vector<int> NewRF;
        for (size_t i = 0; i < QuadFaces.size(); i++) {
          int O = (i < OldFaceOf.size()) ? OldFaceOf[i] : -1;
          bool InR = (O < 0) || ((O < (int)LastRegionFace.size()) &&
                                 LastRegionFace[O]);
          NewLast[i] = InR;
          if (InR)
            NewRF.push_back((int)i);
        }
        LastRegionFace = NewLast;
        NewRegionFaces = NewRF;
        AllF.resize(QuadFaces.size());
        for (size_t i = 0; i < AllF.size(); i++)
          AllF[i] = (int)i;
      }
      snprintf(buf, sizeof(buf),
               "\nlayout refit (%.3f s): %s; singularities of the mesh %d -> %d",
               RSec, Geo::QuadLayoutRefit<ScalarType>::LastMessage().c_str(),
               SingAll0, InnerIrregular(QuadFaces, AllF, QuadFeatures));
      RefitMsg = buf;
      int RLay1, RIn1;
      IrregularKinds(QuadFaces, NewRegionFaces, QuadFeatures, QuadPatch, RLay1,
                     RIn1);
      snprintf(buf, sizeof(buf),
               "\nrefit, region singularities on the layout %d -> %d, inside "
               "patches %d -> %d%s",
               RLay0, RLay1, RIn0, RIn1,
               debug_refit_no_cut ? " (DEBUG: no cut)" : "");
      RefitMsg += buf;
    }
    // Attempt, selected vs used (manifold completion)
    snprintf(buf, sizeof(buf),
             "region of %zu quads (selected %zu) remeshed in %.2f s (%zu new "
             "quads), attempt %d%s, border turns %d",
             RegionSet.size(), Region.size(),
             std::chrono::duration<double>(end - start).count(),
             NewRegionFaces.size(), attempt + 1,
             (attempt > 0) ? " (grown region)" : "",
             Geo::QuadRegionSelector<ScalarType>::LastTurns());
    Msg = buf;
    snprintf(buf, sizeof(buf), "\nregion selection %.3f s (%zu quads of %zu)",
             SelSec, Region.size(), OldFaces.size());
    Msg += buf;
    // Singularities inside the region, before and after
    snprintf(buf, sizeof(buf), "\ninner singularities: %d before -> %d after",
             InnerIrregular(OldFaces, UsedRegion, OldFeatures),
             InnerIrregular(QuadFaces, NewRegionFaces, QuadFeatures));
    Msg += buf;
    // After (and after the refit, if done): on the layout
    // or inside the patches
    {
      int Lay, Ins;
      IrregularKinds(QuadFaces, NewRegionFaces, QuadFeatures, QuadPatch, Lay,
                     Ins);
      snprintf(buf, sizeof(buf),
               " (now %d on the layout, %d inside patches)", Lay, Ins);
      Msg += buf;
    }
    // DEBUG: the triangle features of this region added to the green ones
    {
      const std::vector<PosEdge> &RegSegs =
          Geo::QuadRegionRemesher<ScalarType>::LastTriFeatures();
      // (only the last region, the older ones may no longer exist)
      SourceFeatSegs = SourceFeatSegsFull;
      SourceFeatSegs.insert(SourceFeatSegs.end(), RegSegs.begin(),
                            RegSegs.end());
    }
    Msg += RefitMsg;
    // sizes, samples and times of the remesh
    Msg += std::string("\n") + Geo::QuadRegionRemesher<ScalarType>::LastInfo();
    return true;
  }
  Msg = std::string("region not remeshed:") + Reasons;
  return false;
}

void SketchOnQuads(const std::vector<CoordType> &Line) {
  std::vector<bool> Touched;
  Geo::QuadRegionSelector<ScalarType>::TouchedByLine(
      QuadVert, QuadFaces, Line, AvgQuadEdge() * (ScalarType)0.25, Touched);
  std::vector<std::vector<PosEdge>> Embedded;
  std::string Msg;
  // The region only, without the line (no sketch kept)
  if (debug_region_no_line) {
    bool ok = RemeshQuadRegion(Touched, std::vector<std::vector<CoordType>>(),
                               Embedded, Msg);
    LineMsg = std::string("DEBUG region remeshed without the line") +
              (ok ? ": " : " FAILED: ") + Msg;
    StatLines.push_back(LineMsg);
    UpdateQuadSketchSegs();
    ShowView(ViewQuad);
    return;
  }
  if (!RemeshQuadRegion(Touched, std::vector<std::vector<CoordType>>(1, Line),
                        Embedded, Msg)) {
    LineMsg = std::string("sketch NOT added: ") + Msg;
    StatLines.push_back(LineMsg);
    return;
  }
  if (Embedded.empty() || Embedded[0].empty()) {
    LineMsg = std::string("line not embedded in the region: ") + Msg;
    StatLines.push_back(LineMsg);
    UpdateQuadSketchSegs();
    ShowView(ViewQuad);
    return;
  }
  Sketch S;
  S.on_quad = true;
  S.LineSegs = Embedded[0];
  S.RegionId = RegionHistory.LastId();
  S.Input = Line;
  S.Group = NewSketchGroup;
  Sketches.push_back(S);
  UpdateQuadSketchSegs();
  LineMsg = Msg + ", " + std::to_string(Sketches.back().Segs.size()) +
            " line edges";
  ShowView(ViewQuad);
}

// switching to "quadrangulate all again": the sketches made on the quads
// (region mode) become features of the input, then everything is
// quadrangulated again with them
void QuadSketchesToInputAndRequad() {
  std::vector<std::vector<CoordType>> Lines;
  for (size_t s = 0; s < Sketches.size(); s++) {
    if (!Sketches[s].on_quad)
      continue;
    // The line as embedded in the region (the shape seen
    // in the local mode), embedded on the input without smoothing below
    // The fitted line, the same embedded locally
    if (Sketches[s].Input.size() >= 2) {
      Lines.push_back(Sketches[s].Input);
      continue;
    }
    // One polyline per sketch (sketches without Input)
    std::vector<CoordType> SLine;
    Geo::QuadFeatureLines<ScalarType>::ChainSegsSingle(Sketches[s].LineSegs,
                                                       SLine);
    if (SLine.size() >= 2)
      Lines.push_back(SLine);
  }
  if (Lines.empty() || !hasQuad)
    return;
  int NumOk = 0;
  // The lines are already fitted
  Geo::EmbedPolylineParam<ScalarType> NoSmooth = EmbedNoSmooth();
  for (size_t i = 0; i < Lines.size(); i++) {
    size_t NumSketches = Sketches.size();
    SketchOnTriangles(Lines[i], &NoSmooth);
    if (Sketches.size() > NumSketches)
      NumOk++;
  }
  DoQuadrangulate(); // (drops the old quad sketches)
  LineMsg = std::to_string(NumOk) + " / " + std::to_string(Lines.size()) +
            " quad sketches moved to the input, quadrangulated again";
}

void EmbedUserLine(const std::vector<CoordType> &DrawnLine) {
  if (!loadedMesh || (DrawnLine.size() < 2))
    return;
  // The line made ready once, used as it is below
  // (samples at the spacing of the embedding)
  std::vector<CoordType> Line;
  {
    Geo::SketchPolylineParam<ScalarType> P = FitP;
    P.Step = BBDiag * EmbedP.sample_ratio;
    // In quad edges (the quad edge of the last
    // quadrangulation, or 2% of the diagonal before it)
    ScalarType QuadE = (AbsQuadEdge > 0) ? AbsQuadEdge : BBDiag * (ScalarType)0.02;
    P.WindowLength = fit_window_quads * QuadE;
    P.Mode = (fit_mode == 0)
                 ? Geo::SketchPolylineParam<ScalarType>::FitBSpline
                 : Geo::SketchPolylineParam<ScalarType>::FitLocalBezier;
    P.ControlSpacing = fit_spacing_quads * QuadE;
    bool Closed = false;
    Geo::SketchPolyline<ScalarType>::Process(DrawnLine, P, Line, &Closed);
    if (Line.size() < 2)
      return;
    // The fitted points lie on chords, slightly
    // off the surface, and the embedding reprojects only inside its own
    // smoothing (off now): each point on the closest point of the input
    if (!ConnectivityOrig.empty()) {
      Geo::TriangleSpatialIndex<ScalarType> Grid;
      Grid.Init(VertPosOrig, ConnectivityOrig);
      for (size_t i = 0; i < Line.size(); i++) {
        int IndexF = -1;
        CoordType Q;
        if (Grid.GridClosest(VertPosOrig, ConnectivityOrig, Line[i],
                             Grid.bbox.Diag(), IndexF, Q) &&
            (IndexF >= 0))
          Line[i] = Q;
      }
      if (Closed)
        Line.back() = Line.front();
    }
    if (Closed)
      std::cout << "Sketch: almost a loop, closed" << std::endl;
  }
  // DEBUG (embedding check): embedded in the input only, then measured
  if (debug_embed_only) {
    PreviewDrawn = DrawnLine;
    PreviewFitted = Line;
    DebugSharpPts.clear();
    size_t NumSketches = Sketches.size();
    SketchOnTriangles(Line);
    if (Sketches.size() == NumSketches) {
      LineMsg = "DEBUG: line not embedded";
      return;
    }
    const Sketch &S = Sketches.back();
    // the embedded line as one polyline
    std::vector<CoordType> E;
    Geo::QuadFeatureLines<ScalarType>::ChainSegsSingle(S.Segs, E);
    // sharp turns (as MeshFeatures::Corners of the quadrangulator: the angle
    // between two consecutive feature edges)
    const ScalarType CornerAngle = 45;
    int NumSharp = 0;
    ScalarType MaxTurn = 0;
    for (size_t i = 1; i + 1 < E.size(); i++) {
      CoordType A = E[i] - E[i - 1], B = E[i + 1] - E[i];
      if (!(A.Norm() > 0) || !(B.Norm() > 0))
        continue;
      A.Normalize();
      B.Normalize();
      ScalarType Turn =
          std::acos(std::max((ScalarType)-1,
                             std::min((ScalarType)1, Geo::Dot(A, B)))) *
          180 / (ScalarType)M_PI;
      MaxTurn = std::max(MaxTurn, Turn);
      if (Turn > CornerAngle) {
        NumSharp++;
        DebugSharpPts.push_back(E[i]);
      }
    }
    // distance of the embedded vertices from the fitted line
    std::vector<PosEdge> FitSegs;
    for (size_t i = 0; i + 1 < Line.size(); i++)
      FitSegs.push_back(MakePosEdge(Line[i], Line[i + 1]));
    ScalarType MaxDist = 0;
    for (size_t i = 0; i < E.size(); i++)
      MaxDist = std::max(MaxDist, Geo::QuadFeatureLines<ScalarType>::
                                      PointSegsDist(E[i], FitSegs));
    char buf[256];
    snprintf(buf, sizeof(buf),
             "DEBUG embedded: %zu edges, %zu corners added, max turn %.1f "
             "deg, %d vertices turning > %.0f deg (red), max distance from "
             "the fitted line %.3g (%.2f quad edges); not quadrangulated",
             S.Segs.size(), S.CornersPos.size(), (double)MaxTurn, NumSharp,
             (double)CornerAngle, (double)MaxDist,
             (AbsQuadEdge > 0) ? (double)(MaxDist / AbsQuadEdge) : 0.0);
    LineMsg = buf;
    StatLines.push_back(LineMsg);
    ShowView(ViewInput);
    return;
  }
  // DEBUG (fitting preview): only shown
  if (debug_fit_only) {
    PreviewDrawn = DrawnLine;
    PreviewFitted = Line;
    char buf[160];
    snprintf(buf, sizeof(buf),
             "PREVIEW: drawn %zu points -> fitted %zu points%s (nothing "
             "embedded)",
             DrawnLine.size(), Line.size(),
             (Line.size() > 2) && (Line.front() == Line.back()) ? ", closed"
                                                                 : "");
    LineMsg = buf;
    return;
  }
  // On the quads the line can also go into the
  // input, then the whole quadrangulation is done again
  if ((viewMode == ViewQuad) && !requad_all) {
    SketchOnQuads(Line);
    return;
  }
  bool from_quads = (viewMode == ViewQuad);
  size_t NumSketches = Sketches.size();
  SketchOnTriangles(Line); // also from the remesh view: into the input
  if (from_quads && (Sketches.size() > NumSketches)) {
    std::string Msg = LineMsg;
    DoQuadrangulate();
    LineMsg = Msg + ", quadrangulated again";
  }
}

// ---------------------------------------------------------------------------
// SKETCH REMOVAL
// input: the feature is removed (the next remesh/quadrangulation is without
// it); quads: the feature is removed and the region around it remeshed
// ---------------------------------------------------------------------------

// The input back as loaded and the
// lines still there embedded again (the refinement of a removed line goes
// away with it); quad sketches untouched
void RebuildInputSketches() {
  std::vector<std::vector<CoordType>> Lines;
  std::vector<Sketch> Kept;
  for (size_t i = 0; i < Sketches.size(); i++) {
    if (Sketches[i].on_quad)
      Kept.push_back(Sketches[i]);
    else if (Sketches[i].Input.size() >= 2)
      Lines.push_back(Sketches[i].Input);
  }
  Sketches = Kept;
  VertPos0 = VertPosOrig;
  Connectivity0 = ConnectivityOrig;
  AVGL = AvgEdgeLen(VertPos0, Connectivity0);
  MeshArea = Area(VertPos0, Connectivity0);
  UpdateFeatureCorners();
  for (size_t i = 0; i < Lines.size(); i++)
    SketchOnTriangles(Lines[i]);
}

void RemoveSketch(int Index, bool remesh = true) {
  if ((Index < 0) || (Index >= (int)Sketches.size()))
    return;
  Sketch S = Sketches[Index];
  HoverSketch = -1;
  if (!S.on_quad) {
    // (NEW: the view before, the embeddings again show the input)
    int ViewBefore = viewMode;
    Sketches.erase(Sketches.begin() + Index);
    RebuildInputSketches();
    LineMsg = "sketch removed from the input features";
    // The quads without it
    if (remesh && requad_all && hasQuad && (ViewBefore == ViewQuad)) {
      DoQuadrangulate();
      LineMsg = "sketch removed, quadrangulated again";
      return;
    }
    ShowView(ViewBefore);
    return;
  }
  // its feature edges and the quads around them
  std::vector<int> FeatIdx;
  QuadSketchEdges(S, FeatIdx);
  Sketches.erase(Sketches.begin() + Index);
  std::set<int> Removed(FeatIdx.begin(), FeatIdx.end());
  std::set<int> SketchV;
  std::vector<std::pair<int, int>> Kept;
  for (size_t i = 0; i < QuadFeatures.size(); i++) {
    if (Removed.count((int)i) > 0) {
      SketchV.insert(QuadFeatures[i].first);
      SketchV.insert(QuadFeatures[i].second);
      continue;
    }
    Kept.push_back(QuadFeatures[i]);
  }
  std::vector<std::pair<int, int>> FeaturesBefore = QuadFeatures;
  QuadFeatures = Kept;
  if (remesh && !SketchV.empty()) {
    std::vector<bool> Touched(QuadFaces.size(), false);
    for (size_t i = 0; i < QuadFaces.size(); i++)
      for (size_t j = 0; j < QuadFaces[i].size(); j++)
        if (SketchV.count(QuadFaces[i][j]) > 0)
          Touched[i] = true;
    std::vector<std::vector<PosEdge>> Embedded;
    std::string Msg;
    // The quads still coming from its region remesh
    std::vector<int> SketchRegion = RegionHistory.Faces(S.RegionId);
    // (patches region: the patches it made, whole, plus a ring of patches:
    // the patches along a line alone can be a thin strip)
    if (region_mode == RegionPatches)
      Geo::QuadRegionSelector<ScalarType>::WholePatches(
          QuadVert, QuadFaces, QuadPatch, SketchRegion, 1);
    if (RemeshQuadRegion(Touched, std::vector<std::vector<CoordType>>(),
                         Embedded, Msg,
                         SketchRegion.empty() ? NULL : &SketchRegion))
      LineMsg = std::string("sketch removed, ") + Msg;
    else {
      // not remeshed: the sketch stays as it was (its features too)
      QuadFeatures = FeaturesBefore;
      Sketches.insert(Sketches.begin() + Index, S);
      LineMsg = std::string("sketch NOT removed: ") + Msg;
    }
  } else
    LineMsg = "sketch removed from the quad features";
  UpdateQuadSketchSegs();
  ShowView(viewMode);
}

// ---------------------------------------------------------------------------
// SAVE
// ---------------------------------------------------------------------------

std::string OutName(const std::string &Suffix) {
  std::string Folder = GetDirectory(PathMesh);
  if (Folder.empty())
    Folder = std::string(".");
  return Folder + std::string("/") + GetFilanameNoExtension(PathMesh) + Suffix;
}

void SaveQuad() {
  if (!hasQuad)
    return;
  std::string OutFile = OutName("_quad.obj");
  // the half with its mirror
  if (QuadIsHalf) {
    std::vector<CoordType> V = QuadVert;
    std::vector<std::vector<int>> F = QuadFaces;
    HalfSym.MirrorQuad(V, F, MeshBox.Center(), BBDiag);
    WriteOBJ(OutFile, V, F);
  } else
    WriteOBJ(OutFile, QuadVert, QuadFaces);
  StatLines.push_back(std::string("saved ") + OutFile);
}

void SaveTri() {
  if (!hasTri)
    return;
  std::string OutFile = OutName("_remesh.obj");
  WriteOBJ(OutFile, VertPosR, ConnectivityR);
  StatLines.push_back(std::string("saved ") + OutFile);
}

// ---------------------------------------------------------------------------
// INPUT CALLBACKS
// ---------------------------------------------------------------------------

static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

// cursor in framebuffer pixels, y from below (as the GL viewport)
void CursorPixels(GLFWwindow *window, int &x, int &y) {
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  double xpos, ypos;
  glfwGetCursorPos(window, &xpos, &ypos);
  x = (int)(xpos * xscale);
  y = (int)(display_h - ypos * yscale);
}

// true if the cursor is on the X of the hovered sketch
bool CursorOnX(GLFWwindow *window) {
  if (HoverSketch < 0)
    return false;
  double xpos, ypos;
  glfwGetCursorPos(window, &xpos, &ypos);
  ImVec2 C = XPos();
  float dx = (float)xpos - C.x, dy = (float)ypos - C.y;
  return (std::sqrt(dx * dx + dy * dy) <= XRadius + 2);
}

static void cursor_position_callback(GLFWwindow *window, double xpos,
                                     double ypos) {
  if (LinePick.IsActive()) {
    int x, y;
    CursorPixels(window, x, y);
    LinePick.AddScreenPoint(x, y);
    return;
  }
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

Geo::SketchAction<ScalarType> CurrentAction(char Type) {
  Geo::SketchAction<ScalarType> A;
  A.Type = Type;
  A.View = viewMode;
  A.RequadAll = requad_all;
  A.RegionMode = region_mode;
  A.Refit = region_refit_layout;
  A.PatchRings = region_patch_rings;
  return A;
}

// the action recorded, then done with its settings
void DoAction(const Geo::SketchAction<ScalarType> &A, bool record = true) {
  if (record)
    SLog.Add(A);
  if (A.View != viewMode)
    ShowView(A.View);
  requad_all = A.RequadAll;
  region_mode = A.RegionMode;
  region_refit_layout = A.Refit;
  region_patch_rings = A.PatchRings;
  if (A.Type == 'Q')
    DoQuadrangulate();
  else if (A.Type == 'L')
    EmbedUserLine(A.Line);
  else if (A.Type == 'R')
    RemoveSketch(A.Index);
  else if (A.Type == 'A') {
    // one quadrangulation at the end
    bool requad = requad_all && hasQuad && (viewMode == ViewQuad);
    for (int i = (int)Sketches.size() - 1; i >= 0; i--)
      RemoveSketch(i, !requad);
    if (requad)
      DoQuadrangulate();
  }
}

// the sketch actions of the log queued: one is done at every frame
// (ReplayStep), so the mesh is shown after each of them; the replayed actions
// are not recorded again, so the log replayed is kept (also after another
// crash)
std::vector<Geo::SketchAction<ScalarType>> ReplayQueue;
size_t ReplayNext = 0;

void ReplaySketchLog() {
  std::vector<Geo::SketchAction<ScalarType>> Actions;
  if (!SLog.LoadLast(Actions)) {
    std::cout << "no sketch log " << SketchLogPath << std::endl;
    return;
  }
  // only the sketches, on the mesh as it is now (quadrangulate first, with
  // the same settings): the quadrangulation is not recorded (a 'Q' of an old
  // log is skipped)
  ReplayQueue.clear();
  for (size_t i = 0; i < Actions.size(); i++)
    if (Actions[i].Type != 'Q')
      ReplayQueue.push_back(Actions[i]);
  ReplayNext = 0;
}

void ReplayStep() {
  if (ReplayNext >= ReplayQueue.size())
    return;
  std::cout << "replay action " << ReplayNext + 1 << " / "
            << ReplayQueue.size() << ": " << ReplayQueue[ReplayNext].Type
            << std::endl;
  // (not recorded: the log being replayed stays as it is)
  DoAction(ReplayQueue[ReplayNext++], false);
  ShowView(viewMode);
}

void mouse_button_callback(GLFWwindow *window, int button, int action,
                           int mods) {
  if ((button == GLFW_MOUSE_BUTTON_LEFT) && (action == GLFW_RELEASE)) {
    LeftDown = false;
    if (LinePick.IsActive()) {
      LinePick.End();
      return;
    }
  }
  if (ImGui::GetIO().WantCaptureMouse)
    return;
  if ((button == GLFW_MOUSE_BUTTON_LEFT) && (action == GLFW_PRESS)) {
    // SKETCH REMOVAL: click on the X
    if ((!SpaceDown) && CursorOnX(window)) {
      // with its mirror (the sketches of the same drawn line), from the
      // last one so that the indices of the others do not change
      std::vector<int> ToRemove(1, HoverSketch);
      int G = ((HoverSketch >= 0) && (HoverSketch < (int)Sketches.size()))
                  ? Sketches[HoverSketch].Group
                  : -1;
      if (G >= 0)
        for (int i = 0; i < (int)Sketches.size(); i++)
          if ((i != HoverSketch) && (Sketches[i].Group == G))
            ToRemove.push_back(i);
      std::sort(ToRemove.rbegin(), ToRemove.rend());
      for (size_t i = 0; i < ToRemove.size(); i++) {
        Geo::SketchAction<ScalarType> A = CurrentAction('R');
        A.Index = ToRemove[i];
        DoAction(A);
      }
      return;
    }
    LeftDown = true;
    if (SpaceDown && loadedMesh) {
      LinePick.Begin();
      int x, y;
      CursorPixels(window, x, y);
      LinePick.AddScreenPoint(x, y);
      return;
    }
  }
  if (LinePick.IsActive())
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
  if (ImGui::GetIO().WantCaptureMouse)
    return;
  trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods) {
  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
    trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);
  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
    trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);
  if ((key == GLFW_KEY_SPACE) && (action == GLFW_PRESS) &&
      (!ImGui::GetIO().WantCaptureKeyboard))
    SpaceDown = true;
  if ((key == GLFW_KEY_SPACE) && (action == GLFW_RELEASE)) {
    SpaceDown = false;
    LinePick.End();
  }
}

void InitGLFW_Window() {
  if (!glfwInit()) {
    std::cout << "Error Initializing GLFW" << std::endl;
    exit(0);
  }
  window = glfwCreateWindow(1280, 720, "Quadrangulator", NULL, NULL);
  if (window == NULL) {
    std::cout << "Error Initializing GLFW Window" << std::endl;
    exit(0);
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);
  glfwSetErrorCallback(glfw_error_callback);
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetKeyCallback(window, key_callback);
}

void InitIMGui() {
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui::StyleColorsClassic();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL2_Init();
}

void Help(const char *text) {
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("%s", text);
}

// ---------------------------------------------------------------------------
// PANELS
// ---------------------------------------------------------------------------

void SetRenderWindow() {
  ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, 10),
                          ImGuiCond_FirstUseEver);
  ImGui::Begin("Render", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  const float W = 140.0f;
  int view = viewMode;
  bool changed = ImGui::RadioButton("input", &view, ViewInput);
  if (hasTri) {
    ImGui::SameLine();
    changed |= ImGui::RadioButton("tri remesh", &view, ViewTri);
  }
  if (hasQuad) {
    ImGui::SameLine();
    changed |= ImGui::RadioButton("quads", &view, ViewQuad);
  }
  if (changed)
    ShowView(view);
  ImGui::Checkbox("mesh", &showMesh);
  ImGui::SameLine();
  ImGui::Checkbox("features", &showFeatures);
  ImGui::SameLine();
  ImGui::Checkbox("corners", &showCorners);
  ImGui::Checkbox("sketches", &showSketches);
  ImGui::SameLine();
  if (ImGui::Checkbox("border as feature", &showBorderFeatures))
    ShowView(viewMode);
  ImGui::Checkbox("source features (green, quad view)", &showSourceFeatures);
  ImGui::Checkbox("last region (light blue, quad view)", &showLastRegion);
  Help("the quads made by the last region remesh (sketch or removal)");
  Help("DEBUG: the triangle features the last quadrangulation started from, "
       "drawn over the quads: a quad feature without green under it was "
       "made by the pipeline");
  ImGui::SetNextItemWidth(W);
  ImGui::SliderFloat("line width", &line_thick, 0.5f, 5.f, "%.1f");
  ImGui::End();
}

void SetDrawLinesWindow() {
  ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, 200),
                          ImGuiCond_FirstUseEver);
  ImGui::Begin("Sketch", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  const float W = 140.0f;
  ImGui::TextWrapped(
      "hold SPACE and the left button and move the mouse on the mesh.\n"
      "triangle views: the line becomes a feature of the input.\n"
      "quad view: the patches around the line are remeshed with the line as "
      "a feature, their border kept.\n"
      "hover a sketch and click its X to remove it.");
  ImGui::SetNextItemWidth(W);
  // The drawn line fitted once
  // The fitting of the drawn line (once, when it is
  // drawn; the embedding does not smooth it again)
  ImGui::Text("line fitting");
  ImGui::SameLine();
  ImGui::RadioButton("B-spline", &fit_mode, 0);
  Help("one least squares cubic B-spline through all the points of the line, "
       "a control point every given distance (the extremes kept, a closed "
       "line periodic)");
  ImGui::SameLine();
  ImGui::RadioButton("local Bezier", &fit_mode, 1);
  Help("local Bezier fits on windows of the line, averaged and repeated");
  if (fit_mode == 0) {
    ImGui::SetNextItemWidth(W);
    ImGui::SliderFloat("control points every (quad edges)", &fit_spacing_quads,
                       0.5f, 20.f, "%.1f");
    Help("longer = smoother (fewer control points)");
  } else {
    ImGui::SetNextItemWidth(W);
    ImGui::SliderInt("line fitting steps", &FitP.SmoothIters, 0, 10);
    ImGui::SetNextItemWidth(W);
    ImGui::SliderFloat("line fitting window (quad edges)", &fit_window_quads,
                       0.5f, 10.f, "%.1f");
  }
  ImGui::Checkbox("close almost loops", &FitP.close_loops);
  Help("a line whose start and end come close (less than 10% of its length, "
       "also passing each other) going on in the same direction is cut "
       "there and closed");

  ImGui::Separator();
  ImGui::Text("sketch on the quads");
  int mode = requad_all ? 0 : 1;
  bool was_requad_all = requad_all;
  ImGui::RadioButton("quadrangulate all again", &mode, 0);
  Help("the line becomes a feature of the input and the whole mesh is "
       "quadrangulated again (best quality, as without sketches)");
  ImGui::SameLine();
  ImGui::RadioButton("remesh a region", &mode, 1);
  Help("only the patches around the line are remeshed, their border kept");
  requad_all = (mode == 0);
  // switching to "quadrangulate all again": the region sketches become input
  // features and everything is quadrangulated again (to compare)
  if (requad_all && !was_requad_all)
    QuadSketchesToInputAndRequad();
  ImGui::Text("region (remesh a region)");
  ImGui::RadioButton("distance", &region_mode, RegionDistance);
  Help("the quads within a distance from the line measured along the surface "
       "(in quad edge lengths), then dilate + erode and the regular border");
  ImGui::SameLine();
  ImGui::RadioButton("quad rings", &region_mode, RegionRings);
  Help("the quads touched by the line plus rings of quads (counted through "
       "the vertices, whatever their size), then the regular border");
  ImGui::SameLine();
  ImGui::RadioButton("patches", &region_mode, RegionPatches);
  Help("the whole layout patches touched by the line (plus rings of "
       "patches): the old patches left stay whole, so the layout refit "
       "solves the flow on old and new patches together");
  ImGui::SetNextItemWidth(W);
  if (region_mode == RegionDistance) {
    ImGui::SliderFloat("radius (quad edges)", &region_radius, 1.f, 15.f,
                       "%.1f");
  } else if (region_mode == RegionPatches) {
    ImGui::SliderInt("patch rings", &region_patch_rings, 0, 3);
    // (moved here)
    ImGui::Checkbox("refit the layout after", &region_refit_layout);
    Help("after the region remesh the side subdivisions of the whole layout "
         "are solved again with the current ones as targets (the region border "
         "is no longer fixed): the old patches along the strips can change and "
         "absorb a different subdivision at the two ends of a line; only the "
         "patches with a changed side are tessellated again");
    if (region_refit_layout) {
      ImGui::SetNextItemWidth(W);
      ImGui::SliderInt("refit patch rings (-1 = all)", &region_refit_rings, -1,
                       10);
    }
  } else {
    ImGui::SliderInt("quad rings", &region_face_rings, 1, 10);
  }

  ImGui::Checkbox("smooth the seam", &seam_smooth);
  Help("after a region remesh the region and some rings of quads around it "
       "are smoothed together (features and mesh border fixed)");

  ImGui::Separator();
  int NumIn = 0, NumQ = 0;
  for (size_t i = 0; i < Sketches.size(); i++)
    (Sketches[i].on_quad ? NumQ : NumIn)++;
  ImGui::Text("%d sketches on the input, %d on the quads", NumIn, NumQ);
  if (ImGui::Button("remove all sketches"))
    DoAction(CurrentAction('A'));
  ImGui::SameLine();
  if (ImGui::Button("replay sketch log"))
    ReplaySketchLog();
  Help("the sketches of the log done again, one per frame, with their "
       "settings, on the mesh as it is now (quadrangulate first, with the same "
       "settings): the ones of this session, or if none yet the ones of the "
       "previous one (e.g. after a crash); the log is written next to the mesh "
       "(.sketch_log.txt, the previous one .sketch_log.txt.prev), a replay "
       "does not change it");
  // All the debug options together, closed
  if (ImGui::CollapsingHeader("Debug")) {
    // DEBUG (fitting preview)
    ImGui::Checkbox("DEBUG: only show the fitted line", &debug_fit_only);
    Help("a drawn line is only processed (closure, fitting, projection) and "
         "shown: blue as drawn, orange fitted with its samples; nothing is "
         "embedded or remeshed");
    ImGui::Checkbox("DEBUG: embed in the input only", &debug_embed_only);
    Help("the fitted line is embedded in the input triangles (a sketch on the "
         "input) and not quadrangulated: orange the fitted line, the sketch "
         "edges the embedded one, red the vertices where it turns more than "
         "45 deg (each would be a patch corner); the numbers in the status");
    ImGui::Checkbox("DEBUG: sketch on triangles not remeshed",
                    &debug_sketch_no_remesh);
    Help("the line only embedded in the input, the patch around it not "
         "remeshed");
    ImGui::Checkbox("DEBUG: refit without cutting", &debug_refit_no_cut);
    Help("patches region with the refit: every side changed by the new solve "
         "kept, also far from the region (as the full quadrangulation would "
         "do)");
    ImGui::Checkbox("DEBUG: region without the line", &debug_region_no_line);
    Help("the sketch only selects the region, which is remeshed without the "
         "line: it should give back about the same quads");
    ImGui::Checkbox("DEBUG: save region meshes", &debug_region_save);
    Help("writes region_tri.obj (quads triangulated), region_remesh.obj "
         "(remeshed, lines embedded), region_field_tri.obj (what the "
         "quadrangulator works on), region_features.obj (its features as "
         "lines), region_quad.obj, region_primal_layout.obj, "
         "region_dual_layout.obj (and full_primal_layout.obj, "
         "full_dual_layout.obj at a full quadrangulation), in the working "
         "folder");
  }
  if (!LineMsg.empty()) {
    std::cout << LineMsg << std::endl;
    LineMsg.clear();
  }
  ImGui::End();
}

void SetMainBar() {
  ImGui::Begin("Quadrangulator", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  const float W = 140.0f;

  if (loadedMesh) {
    ImGui::Text("%s", GetFilaname(PathMesh).c_str());
    ImGui::Text("%zu faces, %zu vertices", Connectivity0.size(),
                VertPos0.size());
    if (!LoadMsg.empty())
      ImGui::TextColored(ImVec4(0.8f, 0, 0, 1), "cannot process: %s",
                         LoadMsg.c_str());
  } else
    ImGui::Text("no mesh (pass a .obj or .stl on the command line)");

  // ---- FEATURES (shared by the triangle remesh and the quadrangulation)
  if (ImGui::CollapsingHeader("Features", ImGuiTreeNodeFlags_DefaultOpen)) {
    bool changed = ImGui::Checkbox("sharp features", &use_features);
    if (use_features) {
      int angle = (int)floor(0.5 + feature_angle);
      ImGui::SetNextItemWidth(W);
      changed |= ImGui::SliderInt("sharp edge angle (deg)", &angle, 1, 179);
      Help("edges where the surface bends more than this are features (the "
           "-angle option)");
      feature_angle = angle;
      changed |= ImGui::Checkbox("remove tiny feature fragments",
                                 &clean_features);
    }
    if (changed)
      UpdateFeatureCorners();
    ImGui::TextWrapped("the same features (and the sketches on the input) "
                       "are used by the triangle remesh and the quads");
  }

  // ---- TRIANGLE REMESH (optional)
  if (ImGui::CollapsingHeader("Triangle remesh (optional)")) {
    ImGui::Checkbox("remesh before quadrangulating", &use_tri_remesh);
    float s = (float)tri_ratio;
    ImGui::SetNextItemWidth(W);
    ImGui::SliderFloat("tri size scale", &s, 0.25f, 4.f, "%.2f",
                       ImGuiSliderFlags_Logarithmic);
    Help("< 1 finer, > 1 coarser (automatic size)");
    tri_ratio = s;
    ImGui::SetNextItemWidth(W);
    ImGui::SliderInt("iterations", &RemP.steps, 1, 20);
    ImGui::Checkbox("keep borders", &RemP.fix_borders);
    ImGui::Checkbox("final smooth & projection", &RemP.post_project);
    ImGui::SetNextItemWidth(W);
    ImGui::SliderFloat("max error (% of diag)", &max_error_perc, 0.05f, 5.f,
                       "%.2f", ImGuiSliderFlags_Logarithmic);
    if (ImGui::Button("REMESH TRIANGLES")) {
      StatLines.clear();
      if (RunTriRemesh())
        ShowView(ViewTri);
    }
    ImGui::SameLine();
    if (ImGui::Button("SAVE TRIS"))
      SaveTri();
  }

  // ---- QUADRANGULATION
  if (ImGui::CollapsingHeader("Quadrangulation",
                              ImGuiTreeNodeFlags_DefaultOpen)) {
    float q = (float)QuadSizeMult;
    ImGui::SetNextItemWidth(W);
    ImGui::SliderFloat("quad size scale", &q, 0.25f, 4.f, "%.2f",
                       ImGuiSliderFlags_Logarithmic);
    Help("< 1 smaller quads, > 1 bigger (the -size option)");
    QuadSizeMult = q;
    ImGui::SetNextItemWidth(W);
    ImGui::SliderInt("smooth steps", &quad_smooth_steps, 0, 100);
    ImGui::Checkbox("reliable smooth", &reliable_smooth);
    Help("off = faster but less reliable smoothing (-fastsmooth)");
    ImGui::SetNextItemWidth(W);
    ImGui::SliderFloat("CC-ability", &cc_ability, 0.f, 1.f, "%.2f");
    Help("Catmull-Clark-ability of the patches (0 disables the condition)");
    ImGui::Checkbox("exact features", &exact_features);
    Help("only the quad edges on the triangle feature lines are features (no "
         "edge past the end of an open line); off, the whole patch side with "
         "a feature edge is feature (open lines prolonged to the next patch "
         "corner, as the command line)");

    ImGui::Separator();
    ImGui::Checkbox("symmetry", &Sym.Enabled);
    Help("as -symmetry: only the half on the positive side of the plane is "
         "quadrangulated (and remeshed by the sketches), shown and saved with "
         "its mirror; on the input every line is drawn also mirrored, a line "
         "crossing the plane becomes one symmetric line");
    if (Sym.Enabled) {
      ImGui::SameLine();
      ImGui::Checkbox("show plane", &Sym.ShowPlane);
      ImGui::Text("plane normal");
      ImGui::SameLine();
      ImGui::RadioButton("X", &Sym.Axis, 0);
      ImGui::SameLine();
      ImGui::RadioButton("Y", &Sym.Axis, 1);
      ImGui::SameLine();
      ImGui::RadioButton("Z", &Sym.Axis, 2);
      ImGui::SetNextItemWidth(W);
      ImGui::SliderFloat("plane offset", &Sym.Offset, -0.5f, 0.5f, "%.3f");
      Help("from the center of the mesh box along the normal, in box "
           "diagonals");
    }
    if (QuadIsHalf)
      ImGui::TextWrapped("quads: the half on the positive side, shown and "
                         "saved with its mirror (plane fixed at the "
                         "quadrangulation, quadrangulate again to change it)");

    if (ImGui::Button("QUADRANGULATE"))
      DoQuadrangulate();
    ImGui::SameLine();
    if (ImGui::Button("SAVE QUADS"))
      SaveQuad();
    Help("<input>_quad.obj next to the input");
  }
  for (size_t i = 0; i < StatLines.size(); i++)
    ImGui::TextWrapped("%s", StatLines[i].c_str());
  ImGui::End();
}

// ---------------------------------------------------------------------------
// DRAW
// ---------------------------------------------------------------------------

void DrawQuadMesh() {
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1, 1);
  glColor3d(0.8, 0.8, 0.8);
  bool ColorRegion = showLastRegion &&
                     (LastRegionFace.size() == Connectivity.size());
  for (size_t i = 0; i < Connectivity.size(); i++) {
    // the last remeshed region in light blue
    if (ColorRegion)
      glColor3d(LastRegionFace[i] ? 0.6 : 0.8, LastRegionFace[i] ? 0.8 : 0.8,
                LastRegionFace[i] ? 1.0 : 0.8);
    glBegin(GL_POLYGON);
    glNormal3d(FaceNormals[i].X, FaceNormals[i].Y, FaceNormals[i].Z);
    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      const CoordType &P = VertPos[Connectivity[i][j]];
      glVertex3d(P.X, P.Y, P.Z);
    }
    glEnd();
  }
  glDisable(GL_LIGHTING);
  glLineWidth(line_thick);
  glColor3d(0.2, 0.2, 0.2);
  glBegin(GL_LINES);
  for (size_t i = 0; i < RenderEdges.size(); i++) {
    const CoordType &P0 = VertPos[RenderEdges[i].first];
    const CoordType &P1 = VertPos[RenderEdges[i].second];
    glVertex3d(P0.X, P0.Y, P0.Z);
    glVertex3d(P1.X, P1.Y, P1.Z);
  }
  glEnd();
  glPopAttrib();
}

// screen projection of the visible sketches (ImGui coordinates), with the
// matrices of the mesh just drawn
void ProjectSketches() {
  GLdouble MV[16], PR[16];
  GLint VP[4];
  glGetDoublev(GL_MODELVIEW_MATRIX, MV);
  glGetDoublev(GL_PROJECTION_MATRIX, PR);
  glGetIntegerv(GL_VIEWPORT, VP);
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  for (size_t s = 0; s < Sketches.size(); s++) {
    Sketches[s].Screen.clear();
    // (the input sketches also on the quads when everything is
    // quadrangulated again: they are the features of those quads)
    bool visible = showSketches &&
                   (Sketches[s].on_quad
                        ? (viewMode == ViewQuad)
                        : ((viewMode == ViewInput) ||
                           ((viewMode == ViewQuad) && requad_all)));
    if (!visible)
      continue;
    for (size_t i = 0; i < Sketches[s].Segs.size(); i++) {
      GLdouble x0, y0, z0, x1, y1, z1;
      const CoordType &P0 = Sketches[s].Segs[i].first;
      const CoordType &P1 = Sketches[s].Segs[i].second;
      gluProject(P0.X, P0.Y, P0.Z, MV, PR, VP, &x0, &y0, &z0);
      gluProject(P1.X, P1.Y, P1.Z, MV, PR, VP, &x1, &y1, &z1);
      if ((z0 > 1) || (z1 > 1) || (z0 < 0) || (z1 < 0))
        continue;
      ImVec2 A((float)(x0 / xscale), (float)((display_h - y0) / yscale));
      ImVec2 B((float)(x1 / xscale), (float)((display_h - y1) / yscale));
      Sketches[s].Screen.push_back(std::make_pair(A, B));
    }
  }
}

// hovered sketch and its X (drawn over everything)
void UpdateHover() {
  if (LinePick.IsActive() || SpaceDown) {
    HoverSketch = -1;
    return;
  }
  ImVec2 M = ImGui::GetIO().MousePos;
  // keep the current one while the mouse is around the anchor or the X
  if ((HoverSketch >= 0) && (HoverSketch < (int)Sketches.size())) {
    ImVec2 C = XPos();
    float dA = std::sqrt((M.x - HoverAnchor.x) * (M.x - HoverAnchor.x) +
                         (M.y - HoverAnchor.y) * (M.y - HoverAnchor.y));
    float dX = std::sqrt((M.x - C.x) * (M.x - C.x) + (M.y - C.y) * (M.y - C.y));
    if ((dA > 3 * XRadius) && (dX > XRadius + 4))
      HoverSketch = -1;
  } else
    HoverSketch = -1;
  if (HoverSketch < 0) {
    float Best = HoverDist;
    for (size_t s = 0; s < Sketches.size(); s++)
      for (size_t i = 0; i < Sketches[s].Screen.size(); i++) {
        const ImVec2 &A = Sketches[s].Screen[i].first;
        const ImVec2 &B = Sketches[s].Screen[i].second;
        float D = ScreenSegDist(M, A, B);
        if (D < Best) {
          Best = D;
          HoverSketch = (int)s;
          HoverAnchor = ScreenSegClosest(M, A, B);
        }
      }
  }
  if (HoverSketch < 0)
    return;
  ImDrawList *DL = ImGui::GetForegroundDrawList();
  // the hovered sketch highlighted
  for (size_t i = 0; i < Sketches[HoverSketch].Screen.size(); i++)
    DL->AddLine(Sketches[HoverSketch].Screen[i].first,
                Sketches[HoverSketch].Screen[i].second,
                IM_COL32(255, 60, 60, 255), 3);
  ImVec2 C = XPos();
  DL->AddCircleFilled(C, XRadius, IM_COL32(220, 40, 40, 255));
  DL->AddCircle(C, XRadius, IM_COL32(255, 255, 255, 255), 0, 1.5f);
  float r = XRadius * 0.45f;
  DL->AddLine(ImVec2(C.x - r, C.y - r), ImVec2(C.x + r, C.y + r),
              IM_COL32(255, 255, 255, 255), 2);
  DL->AddLine(ImVec2(C.x - r, C.y + r), ImVec2(C.x + r, C.y - r),
              IM_COL32(255, 255, 255, 255), 2);
}

// the mirror of the half quads (only drawn: the quads are the half), with
// its features and sketches
void DrawMirrorHalf() {
  CoordType O = HalfSym.Origin(MeshBox.Center(), BBDiag);
  CoordType N = HalfSym.Normal();
  glPushMatrix();
  glTranslated(O.X, O.Y, O.Z);
  glScaled(1 - 2 * N.X, 1 - 2 * N.Y, 1 - 2 * N.Z);
  glTranslated(-O.X, -O.Y, -O.Z);
  // the mirror reverses the order of the vertices on the screen: without
  // this the faces are seen as back faces and lit with the normal flipped
  glFrontFace(GL_CW);
  if (showMesh)
    DrawQuadMesh();
  if (showFeatures) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    GLDraw::GLDrawEdges<ScalarType>(VertPos, RenderVertFeatures, 10,
                                    CoordType(1, 0, 1));
    glPopAttrib();
  }
  if (showSketches) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    GLDraw::GLDrawEdges<ScalarType>(VertPos, RenderSketchEdges, 12,
                                    CoordType(1, 0.5, 0));
    glPopAttrib();
  }
  glFrontFace(GL_CCW);
  glPopMatrix();
}

// the symmetry plane, transparent, with its outline
void DrawSymmetryPlane() {
  CoordType C[4];
  Sym.PlaneQuad(MeshBox.Center(), BBDiag, C);
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);
  glColor4d(0.2, 0.5, 1, 0.2);
  glBegin(GL_QUADS);
  for (int i = 0; i < 4; i++)
    glVertex3d(C[i].X, C[i].Y, C[i].Z);
  glEnd();
  glDepthMask(GL_TRUE);
  glLineWidth(line_thick * 2);
  glColor4d(0.2, 0.5, 1, 0.8);
  glBegin(GL_LINE_LOOP);
  for (int i = 0; i < 4; i++)
    glVertex3d(C[i].X, C[i].Y, C[i].Z);
  glEnd();
  glPopAttrib();
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
  trackball.DrawPostApply();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  if (!loadedMesh)
    return;
  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());
  // (before the half: the line picking reads the matrices of the half)
  if (QuadIsHalf && (viewMode == ViewQuad))
    DrawMirrorHalf();

  if (showMesh) {
    if (viewMode == ViewQuad)
      DrawQuadMesh();
    else {
      glPushAttrib(GL_ALL_ATTRIB_BITS);
      glLineWidth(line_thick);
      GLDraw::glColor(CoordType(0.8, 0.8, 0.8));
      if (VertCol.size() == VertPos.size())
        GLDraw::DrawSurfaceMesh<ScalarType>(
            VertPos, Connectivity, FaceNormals, VertNormals, true, true, false,
            false, GLDraw::TRTriangle, NULL, &VertCol);
      glPopAttrib();
    }
  }
  // SKETCH: unproject on the depth of the mesh just drawn
  LinePick.ProcessPending();
  LinePick.Draw();
  if (showFeatures) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    GLDraw::GLDrawEdges<ScalarType>(VertPos, RenderVertFeatures, 10,
                                    CoordType(1, 0, 1));
    glPopAttrib();
  }
  // DEBUG: source triangle features over the quads
  if (showSourceFeatures && (viewMode == ViewQuad) && !SourceFeatSegs.empty()) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);
    glLineWidth(line_thick * 3);
    glColor3d(0, 0.7, 0);
    glBegin(GL_LINES);
    for (size_t i = 0; i < SourceFeatSegs.size(); i++) {
      const CoordType &P0 = SourceFeatSegs[i].first;
      const CoordType &P1 = SourceFeatSegs[i].second;
      glVertex3d(P0.X, P0.Y, P0.Z);
      glVertex3d(P1.X, P1.Y, P1.Z);
    }
    glEnd();
    glPopAttrib();
  }
  // DEBUG (fitting preview): the drawn line (blue) and the fitted one
  // (orange, its samples as points)
  if ((debug_fit_only || debug_embed_only) && (PreviewDrawn.size() >= 2)) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);
    glLineWidth(line_thick * 2);
    glColor3d(0.1, 0.3, 1);
    glBegin(GL_LINE_STRIP);
    for (size_t i = 0; i < PreviewDrawn.size(); i++)
      glVertex3d(PreviewDrawn[i].X, PreviewDrawn[i].Y, PreviewDrawn[i].Z);
    glEnd();
    glLineWidth(line_thick * 3);
    glColor3d(1, 0.5, 0);
    glBegin(GL_LINE_STRIP);
    for (size_t i = 0; i < PreviewFitted.size(); i++)
      glVertex3d(PreviewFitted[i].X, PreviewFitted[i].Y, PreviewFitted[i].Z);
    glEnd();
    glPointSize(line_thick * 5);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < PreviewFitted.size(); i++)
      glVertex3d(PreviewFitted[i].X, PreviewFitted[i].Y, PreviewFitted[i].Z);
    glEnd();
    // DEBUG (embedding check): the sharp turns of the embedded line
    if (debug_embed_only && !DebugSharpPts.empty()) {
      glPointSize(line_thick * 12);
      glColor3d(1, 0, 0);
      glBegin(GL_POINTS);
      for (size_t i = 0; i < DebugSharpPts.size(); i++)
        glVertex3d(DebugSharpPts[i].X, DebugSharpPts[i].Y, DebugSharpPts[i].Z);
      glEnd();
    }
    glPopAttrib();
  }
  if (showSketches) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    GLDraw::GLDrawEdges<ScalarType>(VertPos, RenderSketchEdges, 12,
                                    CoordType(1, 0.5, 0));
    glPopAttrib();
  }
  if (showCorners) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    GLDraw::DrawVertices<ScalarType>(VertPos, RenderVertCorners,
                                     CoordType(0, 0, 0), 30);
    glPopAttrib();
  }
  if (Sym.Enabled && Sym.ShowPlane)
    DrawSymmetryPlane();
  ProjectSketches();
}

int main(int argc, char *argv[]) {
  RemP.post_project = true;
  RemP.steps = 10;

  if (argc == 2) {
    PathMesh = std::string(argv[1]);
    LoadMesh(PathMesh);
  }

  InitGLFW_Window();
  InitIMGui();
  trackball.center = Geo::Point3<float>(0, 0, 0);
  trackball.radius = 1;
  glewInit();

  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    {
      SetMainBar();
      SetRenderWindow();
      SetDrawLinesWindow();
      UpdateHover(); // screen positions of the previous frame
    }
    ImGui::Render();

    GLDrawMesh();

    // SKETCH: a completed line is embedded / remeshed
    std::vector<CoordType> Line;
    ReplayStep(); // (a replay in progress: its next action)
    if (LinePick.TakeFinished(Line)) {
      // with the symmetry the line and its mirror (or one symmetric line),
      // each one an action of its own (so the log replays them as drawn)
      std::vector<std::vector<CoordType>> Lines(1, Line);
      // a region of the half quads: the line brought on the half
      if (QuadIsHalf && (viewMode == ViewQuad) && !requad_all)
        HalfSym.HalfLine(Line, MeshBox.Center(), BBDiag, Lines[0]);
      else if (Sym.Enabled)
        Sym.Lines(Line, MeshBox.Center(), BBDiag, Lines);
      NewSketchGroup = (Lines.size() > 1) ? NumSketchGroups++ : -1;
      for (size_t i = 0; i < Lines.size(); i++) {
        Geo::SketchAction<ScalarType> A = CurrentAction('L');
        A.Line = Lines[i];
        DoAction(A);
      }
      NewSketchGroup = -1;
    }

    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
    glfwMakeContextCurrent(window);
    glfwSwapBuffers(window);
  }

  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
  glfwDestroyWindow(window);
  glfwTerminate();
  return 0;
}
