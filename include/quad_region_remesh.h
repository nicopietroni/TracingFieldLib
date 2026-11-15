/*
 * //Non-Commercial License Agreement
 * //Copyright (c) 2025, Nico Pietroni
 * //
 * //1. Permission is hereby granted to use, copy, modify, and distribute
 * //   this software for academic research, teaching, and personal
 * //   non-commercial purposes, provided that this copyright notice is
 * //   included in all copies or substantial portions of the software.
 * //
 * //2. Commercial use of this software, including but not limited to use in
 * //   commercial products, services, tools, or consulting, is strictly
 * //   prohibited without prior written permission from the copyright holder.
 * //
 * //3. This software is provided "as is", without warranty of any kind,
 * //   express or implied, including but not limited to the warranties of
 * //   merchantability, fitness for a particular purpose, and
 * //   noninfringement. In no event shall the authors or copyright holders
 * //   be liable for any claim, damages, or other liability, arising from,
 * //   out of, or in connection with the software or the use or other
 * //   dealings in the software.
 * //
 * //4. For commercial licensing, please contact:
 * //   Nico Pietroni
 * //   Email: nico.pietroni@gmail.com
 */

#ifndef QUAD_REGION_REMESH
#define QUAD_REGION_REMESH

// Remeshing of a region of a quad mesh keeping its border as it is:
//  1. the region faces are copied and triangulated (shorter diagonal), the
//     border loops of the region become fixed polylines of samples;
//  2. the triangulation is remeshed (TriRemesh) with the border fixed and the
//     old feature lines kept as feature edges (or re-embedded as polylines,
//     reembed_old_features);
//  3. the new feature lines (3D polylines near the surface) are embedded
//     with MeshEmbedPolyline, the triangulated region as target surface;
//  4. the cross field is computed as in the batch pipeline;
//  5. the quadrangulation pipeline (MeshQuadrangulator) runs with the border
//     loops as fixed sides: every patch side on the border gets exactly the
//     old border vertices, the patch corners on the border are snapped to
//     them and the border is not smoothed;
//  6. the new quads replace the region (shared border vertices merged by
//     position) and the feature edges of the whole mesh are updated.
// If the new border does not match the old one the input mesh
// is left unchanged and false is returned.

#include "IO/io_obj.h"
#include "IO/io_field.h"
#include "fixed_polyline_sides.h"
#include "quad_region_history.h"
#include "quad_region_select.h"
#include "mesh_embed_polyline.h"
#include "mesh_functions.h"
#include "point3.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <field_graph/refine_for_tracing.h>
#include <init_cross_field.h>
#include <iostream>
#include <fstream>
#include <functional>
#include <map>
#include <mesh_quadrangulation.h>
#include <set>
#include <triangular_remesh.h>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> struct QuadRegionRemeshParam {
  // the refinements for the field and the tracing (RefineForFieldComputation,
  // RefineForConnection) on the region; they never split the region border
  // edges: every new border vertex would be a place where a separatrix can
  // end in the middle of an old border edge (a layout corner that is not an
  // old vertex), so the border vertices stay only the old ones
  bool field_refinements = true;
  // target quad edge, times the average edge of the region quads (1 = same
  // density as before; the border is fixed anyway)
  ScalarType QuadSizeMult = 1;
  // Target quad edge in absolute units; when > 0 it is
  // used instead of the average edge of the region (times QuadSizeMult), so
  // that repeated region remeshes keep the size chosen at the start
  ScalarType TargetQuadEdge = -1;
  // target edge of the triangle remeshing, times the target quad edge
  ScalarType TriEdgeRatio = 0.5;
  int TriRemeshSteps = 10;
  ScalarType FeatureCornerAngle = 45;
  // true: the old features are not kept by the triangle remesher but
  // embedded again (as polylines) with MeshEmbedPolyline, as the new ones
  bool reembed_old_features = false;
  EmbedPolylineParam<ScalarType> EmbedParam;
  // cross field (as in BatchProcessField of the command line)
  int FieldNRing = 3;
  // quadrangulation
  int QuadSmoothSteps = 10;
  ScalarType CCAbility = 0.25;
  // Weight of the equality of the opposite patch sides
  // against the target size in the side subdivision (MeshQuadrangulator
  // RegVsIso, 0.5 as in the full quadrangulation); with the border fixed the
  // sizes are already given by it, a higher value prefers equal opposite
  // sides (no irregular tessellation of the patches) to the target size
  ScalarType RegVsIso = 0.5;
  int NumSamples = 6000;
  // NumSamples is the number of tracing samples for
  // the WHOLE mesh (as in the full quadrangulation); the region gets the
  // part of its area, so the density of the traced paths (and with it the
  // patch layout and the time) is the same as in the full quadrangulation.
  // false = NumSamples on the region alone (the old behaviour)
  bool scale_samples = true;
  int MinSamples = 300;
  bool reliable_smooth = true;
  // After the merge the new region and SeamRings rings
  // of quads around it are smoothed together (the region border was fixed in
  // the quadrangulation, the corners snapped on it), reprojected on the old
  // surface; the feature vertices and the mesh border stay. 0 steps = off
  int SeamRings = 3;
  int SeamSmoothSteps = 10;
  ScalarType FixedSideTolRatio = 0.1;
  // debug: save the intermediate meshes (region_tri.obj, region_remesh.obj,
  // region_quad.obj)
  bool save_debug = false;
  // The quadrangulation marks as feature the whole patch
  // side containing a feature edge, so an open line is prolonged up to the
  // next corner; true keeps only the quad feature edges on the triangle
  // features (old ones and new lines)
  bool trim_features = false;
  // The region is completed where it would not be a
  // manifold (two parts touching at a vertex: all the faces around it are
  // added) and the faces almost enclosed (all sides but one on the region)
  // are added; the region actually used is returned in UsedRegion
  bool fix_region = true;
};

template <class ScalarType> class QuadRegionRemesher {
  typedef Geo::Point3<ScalarType> CoordType;
  typedef std::pair<CoordType, CoordType> EdgePos;

  static EdgePos MakeEdgePos(const CoordType &P0, const CoordType &P1) {
    return EdgePos(std::min(P0, P1), std::max(P0, P1));
  }

  static ScalarType AvgPolyEdge(const std::vector<CoordType> &Vert,
                                const std::vector<std::vector<int>> &Faces) {
    ScalarType L = 0;
    int Num = 0;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        L += (Vert[Faces[i][j]] - Vert[Faces[i][(j + 1) % Faces[i].size()]])
                 .Norm();
        Num++;
      }
    return (Num > 0) ? (L / Num) : 0;
  }

  // copy of the faces FaceIdx, vertices re-indexed (OrigV: original index)
  static void CopyFaces(const std::vector<CoordType> &Vert,
                        const std::vector<std::vector<int>> &Faces,
                        const std::vector<int> &FaceIdx,
                        std::vector<CoordType> &SubVert,
                        std::vector<std::vector<int>> &SubFaces,
                        std::vector<int> &OrigV) {
    SubVert.clear();
    SubFaces.clear();
    OrigV.clear();
    std::map<int, int> NewIdx;
    for (size_t i = 0; i < FaceIdx.size(); i++) {
      const std::vector<int> &F = Faces[FaceIdx[i]];
      std::vector<int> NewF;
      for (size_t j = 0; j < F.size(); j++) {
        if (NewIdx.count(F[j]) == 0) {
          NewIdx[F[j]] = (int)SubVert.size();
          SubVert.push_back(Vert[F[j]]);
          OrigV.push_back(F[j]);
        }
        NewF.push_back(NewIdx[F[j]]);
      }
      SubFaces.push_back(NewF);
    }
  }

  // quads split along the shorter diagonal, other polygons as a fan
  // The faces with a vertex in the box of Pts enlarged
  // by Margin, as a mesh of its own (a big reference surface is not searched
  // for a small region: the spatial grids were built on the whole mesh at
  // every region remesh)
  static void CropToBox(const std::vector<CoordType> &Vert,
                        const std::vector<std::vector<int>> &Faces,
                        const std::vector<CoordType> &Pts, ScalarType Margin,
                        std::vector<CoordType> &OutVert,
                        std::vector<std::vector<int>> &OutFaces) {
    OutVert.clear();
    OutFaces.clear();
    if (Pts.empty())
      return;
    CoordType Lo = Pts[0], Hi = Pts[0];
    for (size_t i = 1; i < Pts.size(); i++) {
      Lo.X = std::min(Lo.X, Pts[i].X);
      Lo.Y = std::min(Lo.Y, Pts[i].Y);
      Lo.Z = std::min(Lo.Z, Pts[i].Z);
      Hi.X = std::max(Hi.X, Pts[i].X);
      Hi.Y = std::max(Hi.Y, Pts[i].Y);
      Hi.Z = std::max(Hi.Z, Pts[i].Z);
    }
    auto Inside = [&](const CoordType &P) {
      return (P.X >= Lo.X - Margin) && (P.X <= Hi.X + Margin) &&
             (P.Y >= Lo.Y - Margin) && (P.Y <= Hi.Y + Margin) &&
             (P.Z >= Lo.Z - Margin) && (P.Z <= Hi.Z + Margin);
    };
    std::vector<int> NewIdx(Vert.size(), -1);
    for (size_t i = 0; i < Faces.size(); i++) {
      bool Take = false;
      for (size_t j = 0; (j < Faces[i].size()) && !Take; j++)
        Take = Inside(Vert[Faces[i][j]]);
      if (!Take)
        continue;
      std::vector<int> F;
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int V = Faces[i][j];
        if (NewIdx[V] < 0) {
          NewIdx[V] = (int)OutVert.size();
          OutVert.push_back(Vert[V]);
        }
        F.push_back(NewIdx[V]);
      }
      OutFaces.push_back(F);
    }
  }

  static void Triangulate(const std::vector<CoordType> &Vert,
                          const std::vector<std::vector<int>> &PolyFaces,
                          std::vector<std::vector<int>> &TriFaces) {
    TriFaces.clear();
    for (size_t i = 0; i < PolyFaces.size(); i++) {
      const std::vector<int> &F = PolyFaces[i];
      if (F.size() < 3)
        continue;
      if (F.size() == 4) {
        ScalarType D02 = (Vert[F[0]] - Vert[F[2]]).Norm();
        ScalarType D13 = (Vert[F[1]] - Vert[F[3]]).Norm();
        if (D02 <= D13) {
          TriFaces.push_back({F[0], F[1], F[2]});
          TriFaces.push_back({F[0], F[2], F[3]});
        } else {
          TriFaces.push_back({F[0], F[1], F[3]});
          TriFaces.push_back({F[1], F[2], F[3]});
        }
        continue;
      }
      for (size_t j = 1; j + 1 < F.size(); j++)
        TriFaces.push_back({F[0], F[(int)j], F[(int)j + 1]});
    }
  }

  // feature edges (vertex pairs) chained into polylines of positions, cut at
  // the vertices that do not have two feature edges (closed lines repeat the
  // first point)
  static void ChainFeatures(const std::vector<CoordType> &Vert,
                            const std::vector<std::pair<int, int>> &FEdges,
                            std::vector<std::vector<CoordType>> &Lines) {
    Lines.clear();
    std::map<int, std::vector<int>> Adj;
    std::set<std::pair<int, int>> Used;
    for (size_t i = 0; i < FEdges.size(); i++) {
      int V0 = FEdges[i].first, V1 = FEdges[i].second;
      if (V0 == V1)
        continue;
      std::pair<int, int> Key(std::min(V0, V1), std::max(V0, V1));
      if (Used.count(Key) > 0)
        continue;
      Used.insert(Key);
      Adj[V0].push_back(V1);
      Adj[V1].push_back(V0);
    }
    std::set<std::pair<int, int>> Visited;
    // first from the extremes and junctions, then the closed loops
    for (int pass = 0; pass < 2; pass++) {
      typename std::map<int, std::vector<int>>::iterator It;
      for (It = Adj.begin(); It != Adj.end(); It++) {
        int Start = (*It).first;
        if ((pass == 0) && ((*It).second.size() == 2))
          continue;
        for (size_t k = 0; k < (*It).second.size(); k++) {
          int Next = (*It).second[k];
          std::pair<int, int> Key(std::min(Start, Next), std::max(Start, Next));
          if (Visited.count(Key) > 0)
            continue;
          std::vector<int> Seq(1, Start);
          int Prev = Start, Curr = Next;
          Visited.insert(Key);
          Seq.push_back(Curr);
          while ((Adj[Curr].size() == 2) && (Curr != Start)) {
            int Cand = (Adj[Curr][0] == Prev) ? Adj[Curr][1] : Adj[Curr][0];
            std::pair<int, int> KeyN(std::min(Curr, Cand), std::max(Curr, Cand));
            if (Visited.count(KeyN) > 0)
              break;
            Visited.insert(KeyN);
            Prev = Curr;
            Curr = Cand;
            Seq.push_back(Curr);
          }
          std::vector<CoordType> Line;
          for (size_t s = 0; s < Seq.size(); s++)
            Line.push_back(Vert[Seq[s]]);
          Lines.push_back(Line);
        }
      }
    }
  }

  // corners of the feature lines: not two feature edges, or a turn sharper
  // than AngleDeg
  static void FeatureCornersOf(const std::vector<CoordType> &Vert,
                               const std::vector<std::pair<int, int>> &FEdges,
                               ScalarType AngleDeg,
                               std::vector<bool> &IsCorner) {
    IsCorner.assign(Vert.size(), false);
    std::map<int, std::vector<int>> Adj;
    for (size_t i = 0; i < FEdges.size(); i++) {
      Adj[FEdges[i].first].push_back(FEdges[i].second);
      Adj[FEdges[i].second].push_back(FEdges[i].first);
    }
    ScalarType CosThr = std::cos(AngleDeg * M_PI / 180.0);
    typename std::map<int, std::vector<int>>::iterator It;
    for (It = Adj.begin(); It != Adj.end(); It++) {
      int V = (*It).first;
      if ((*It).second.size() != 2) {
        IsCorner[V] = true;
        continue;
      }
      CoordType D0 = Vert[V] - Vert[(*It).second[0]];
      CoordType D1 = Vert[(*It).second[1]] - Vert[V];
      ScalarType L0 = D0.Norm(), L1 = D1.Norm();
      if ((L0 <= 0) || (L1 <= 0))
        continue;
      // turn angle between the two directions
      if (Geo::Dot(D0, D1) < CosThr * L0 * L1)
        IsCorner[V] = true;
    }
  }

public:
  // Why the last Remesh returned false (empty if ok)
  // How many border sides of the region quadrangulation were
  // fixed by the border polylines (set by ConstrainSides)
  static std::string &LastFixInfo() {
    static std::string S;
    return S;
  }

  static std::string &LastMessage() {
    static std::string Msg;
    return Msg;
  }
  // The triangle features the last Remesh quadrangulated
  // (old features of the region and new lines, as segments)
  static std::vector<EdgePos> &LastTriFeatures() {
    static std::vector<EdgePos> Segs;
    return Segs;
  }

  // Sizes, samples and times of the last successful Remesh
  // The dual patches of the last remesh by sides
  static std::string &LastLayout() {
    static std::string L;
    return L;
  }

  static std::string &LastInfo() {
    static std::string Msg;
    return Msg;
  }

  // completes the region to a manifold (see fix_region; the code is in
  // QuadRegionSelector, used also by the patches region)
  static void MakeRegionManifold(const std::vector<CoordType> &Vert,
                                 const std::vector<std::vector<int>> &Faces,
                                 std::vector<int> &Region) {
    QuadRegionSelector<ScalarType>::MakeRegionManifold(Vert, Faces, Region);
  }

  // Area of a set of polygons (fan triangulation)
  static ScalarType PolyArea(const std::vector<CoordType> &Vert,
                             const std::vector<std::vector<int>> &Faces) {
    ScalarType A = 0;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 1; j + 1 < Faces[i].size(); j++)
        A += Geo::Cross(Vert[Faces[i][j]] - Vert[Faces[i][0]],
                        Vert[Faces[i][j + 1]] - Vert[Faces[i][0]])
                 .Norm() *
             (ScalarType)0.5;
    return A;
  }

  // Keeps only the feature edges QF (of the mesh QV)
  // whose extremes and middle are within Tol of the segments Lines (a
  // uniform grid on the segments, so it is fast also on a whole mesh)
  static void TrimFeaturesToLines(const std::vector<CoordType> &QV,
                                  std::vector<std::pair<int, int>> &QF,
                                  const std::vector<EdgePos> &Lines,
                                  const ScalarType Tol) {
    if (Lines.empty()) {
      QF.clear();
      return;
    }
    // cell size: at least Tol, about the average segment
    ScalarType AvgL = 0;
    CoordType MinP = Lines[0].first;
    for (size_t i = 0; i < Lines.size(); i++) {
      AvgL += (Lines[i].first - Lines[i].second).Norm();
      const CoordType *Pts[2] = {&Lines[i].first, &Lines[i].second};
      for (int e = 0; e < 2; e++) {
        MinP.X = std::min(MinP.X, Pts[e]->X);
        MinP.Y = std::min(MinP.Y, Pts[e]->Y);
        MinP.Z = std::min(MinP.Z, Pts[e]->Z);
      }
    }
    AvgL /= Lines.size();
    ScalarType Cell = std::max(Tol, AvgL);
    if (!(Cell > 0))
      Cell = 1;
    MinP = MinP - CoordType(Tol, Tol, Tol) * 2;
    auto CellOf = [&](const CoordType &P) {
      return std::make_tuple((long long)std::floor((P.X - MinP.X) / Cell),
                             (long long)std::floor((P.Y - MinP.Y) / Cell),
                             (long long)std::floor((P.Z - MinP.Z) / Cell));
    };
    std::map<std::tuple<long long, long long, long long>, std::vector<int>>
        Grid;
    for (size_t i = 0; i < Lines.size(); i++) {
      CoordType A = Lines[i].first, B = Lines[i].second;
      CoordType Lo(std::min(A.X, B.X) - Tol, std::min(A.Y, B.Y) - Tol,
                   std::min(A.Z, B.Z) - Tol);
      CoordType Hi(std::max(A.X, B.X) + Tol, std::max(A.Y, B.Y) + Tol,
                   std::max(A.Z, B.Z) + Tol);
      std::tuple<long long, long long, long long> C0 = CellOf(Lo), C1 = CellOf(Hi);
      for (long long x = std::get<0>(C0); x <= std::get<0>(C1); x++)
        for (long long y = std::get<1>(C0); y <= std::get<1>(C1); y++)
          for (long long z = std::get<2>(C0); z <= std::get<2>(C1); z++)
            Grid[std::make_tuple(x, y, z)].push_back((int)i);
    }
    auto Near = [&](const CoordType &P) {
      typename std::map<std::tuple<long long, long long, long long>,
                        std::vector<int>>::const_iterator It =
          Grid.find(CellOf(P));
      if (It == Grid.end())
        return false;
      for (size_t k = 0; k < (*It).second.size(); k++) {
        const EdgePos &E = Lines[(*It).second[k]];
        CoordType D = E.second - E.first;
        ScalarType L2 = Geo::Dot(D, D);
        ScalarType t = 0;
        if (L2 > 0)
          t = Geo::Dot(P - E.first, D) / L2;
        t = std::max((ScalarType)0, std::min((ScalarType)1, t));
        if ((P - (E.first + D * t)).Norm() <= Tol)
          return true;
      }
      return false;
    };
    std::vector<std::pair<int, int>> Kept;
    for (size_t i = 0; i < QF.size(); i++) {
      const CoordType &P0 = QV[QF[i].first];
      const CoordType &P1 = QV[QF[i].second];
      if (Near(P0) && Near(P1) && Near((P0 + P1) * (ScalarType)0.5))
        Kept.push_back(QF[i]);
    }
    QF = Kept;
  }

  // Public, used also by the interactive quadrangulator
  // cross field on the triangle mesh, as BatchProcessField of the command
  // line (refine, curvature field, global smoothing, refine for connection)
  static void ComputeField(std::vector<CoordType> &TriVert,
                           std::vector<std::vector<int>> &TriFaces,
                           std::vector<std::pair<int, int>> &TriFeatures,
                           std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                           int NRing, bool Refine = true,
                           bool KeepBorderEdges = false) {
    std::vector<CoordType> FaceNormal, VertNormal;
    std::vector<std::vector<int>> NextF, NextE;
    std::vector<Field::CrossF<ScalarType>> VertCurv;
    std::vector<int> SingIndex, SingValue;

    if (Refine)
      Geo::RefineForTracing<ScalarType>::RefineForFieldComputation(
          TriVert, TriFaces, TriFeatures, true, false, KeepBorderEdges);
    Geo::getFFAdjacency(TriFaces, NextF, NextE);
    ComputeNormals(TriVert, TriFaces, FaceNormal, VertNormal);

    FieldSmoothParam<ScalarType> SParam;
    SParam.FixedCross.clear();
    ComputeCurvatureField<ScalarType>(TriVert, TriFaces, FaceNormal,
                                      VertNormal, VertCurv, FaceCurv, NRing);

    std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;
    std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;
    SParam.Features = TriFeatures;
    SmoothGlobalPolyvector(TriVert, TriFaces, FaceNormal, FaceCurv, SParam);
    GetFaceSingularities(FaceCurv, TriVert, TriFaces, NextF, NextE, SingIndex,
                         SingValue);
    ComputeNormals(TriVert, TriFaces, FaceNormal, VertNormal);
    SetVertCrossFromFace(TriVert, TriFaces, VertNormal, FaceCurv, VertCurv);
    RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
    RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);

    if (Refine)
      Geo::RefineForTracing<ScalarType>::RefineForConnection(
          TriVert, TriFaces, FaceCurv, TriFeatures, 45, false,
          KeepBorderEdges);
    for (size_t i = 0; i < TriFeatures.size(); i++)
      if (TriFeatures[i].first > TriFeatures[i].second)
        std::swap(TriFeatures[i].first, TriFeatures[i].second);
    std::sort(TriFeatures.begin(), TriFeatures.end());
    TriFeatures.erase(std::unique(TriFeatures.begin(), TriFeatures.end()),
                      TriFeatures.end());
  }

public:
  // QuadVert, QuadFaces: the quad mesh, changed in place.
  // RegionFaces: the faces to remesh.
  // QuadFeatures: feature edges of the quad mesh (vertex pairs), updated.
  // NewFeatureLines: feature lines to add (3D polylines near the surface).
  // NewRegionFaces (optional): the indices of the new faces.
  static bool Remesh(std::vector<CoordType> &QuadVert,
                     std::vector<std::vector<int>> &QuadFaces,
                     const std::vector<int> &InRegionFaces,
                     std::vector<std::pair<int, int>> &QuadFeatures,
                     const std::vector<std::vector<CoordType>> &NewFeatureLines,
                     const QuadRegionRemeshParam<ScalarType> &Param =
                         QuadRegionRemeshParam<ScalarType>(),
                     std::vector<int> *NewRegionFaces = NULL,
                     std::vector<int> *NewRegionPatch = NULL,
                     std::vector<std::vector<EdgePos>> *EmbeddedLines = NULL,
                     std::vector<int> *UsedRegion = NULL,
                     QuadRegionHistory<ScalarType> *History = NULL,
                     std::vector<int> *NewRegionPatchCorner = NULL) {
    // NewRegionPatchCorner, if given, the corner mask
    // of every new quad (as MeshQuadrangulator::GetQuadPatchCorner), in the
    // order of NewRegionFaces
    // The region actually remeshed
    std::vector<int> RegionFaces = InRegionFaces;
    std::sort(RegionFaces.begin(), RegionFaces.end());
    RegionFaces.erase(std::unique(RegionFaces.begin(), RegionFaces.end()),
                      RegionFaces.end());
    if (Param.fix_region)
      MakeRegionManifold(QuadVert, QuadFaces, RegionFaces);
    if (UsedRegion != NULL)
      *UsedRegion = RegionFaces;
    LastMessage().clear();
    LastInfo().clear();
    auto TStart = std::chrono::high_resolution_clock::now();
    // Time and size of every stage, printed and in
    // LastInfo
    std::string StageTimes;
    LastLayout().clear();
    auto TLap = TStart;
    auto Lap = [&](const char *Name, size_t NV, size_t NF) {
      auto Now = std::chrono::high_resolution_clock::now();
      double Sec = std::chrono::duration<double>(Now - TLap).count();
      TLap = Now;
      char b[256];
      snprintf(b, sizeof(b), "\n  %-28s %8.3f s  (%zu vert, %zu faces)", Name,
               Sec, NV, NF);
      StageTimes += b;
      std::cout << "Region remesh stage:" << b << std::endl;
    };
    if (RegionFaces.empty()) {
      LastMessage() = "empty region";
      return false;
    }

    // 1. REGION ------------------------------------------------------------
    std::vector<CoordType> RegVert;
    std::vector<std::vector<int>> RegPolyFaces, RegTriFaces;
    std::vector<int> OrigV;
    CopyFaces(QuadVert, QuadFaces, RegionFaces, RegVert, RegPolyFaces, OrigV);
    Triangulate(RegVert, RegPolyFaces, RegTriFaces);
    // Which step breaks the manifoldness
    if (!(EdgeManifold<ScalarType>(RegVert, RegTriFaces) &&
          VertManifold<ScalarType>(RegVert, RegTriFaces))) {
      LastMessage() = "the region itself is not manifold";
      return false;
    }
    ScalarType AvgQuadE = AvgPolyEdge(RegVert, RegPolyFaces);
    if (!(AvgQuadE > 0)) {
      LastMessage() = "degenerate region";
      return false;
    }
    if (Param.save_debug) {
      WriteOBJ("region_tri.obj", RegVert, RegTriFaces);
      // The old quads of the region and the rest of the mesh
      // around it, to compare with region_quad.obj
      WriteOBJ("region_old_quad.obj", RegVert, RegPolyFaces);
      std::set<int> InRegion(RegionFaces.begin(), RegionFaces.end());
      std::vector<int> RestFaces;
      for (size_t i = 0; i < QuadFaces.size(); i++)
        if (InRegion.count((int)i) == 0)
          RestFaces.push_back((int)i);
      std::vector<CoordType> RestVert;
      std::vector<std::vector<int>> RestPolyFaces;
      std::vector<int> RestOrigV;
      CopyFaces(QuadVert, QuadFaces, RestFaces, RestVert, RestPolyFaces,
                RestOrigV);
      WriteOBJ("region_rest_quad.obj", RestVert, RestPolyFaces);
    }
    Lap("region extraction", RegVert.size(), RegTriFaces.size());

    // border loops of the region = fixed sides (closed, first point repeated)
    std::vector<std::vector<int>> BorderLoops;
    getBorderSequncesVert(RegTriFaces, BorderLoops);
    std::vector<std::vector<CoordType>> FixedPolylines;
    std::set<CoordType> FixedPos;
    for (size_t i = 0; i < BorderLoops.size(); i++) {
      if (BorderLoops[i].size() < 3)
        continue;
      std::vector<CoordType> Loop;
      for (size_t j = 0; j < BorderLoops[i].size(); j++) {
        Loop.push_back(RegVert[BorderLoops[i][j]]);
        FixedPos.insert(RegVert[BorderLoops[i][j]]);
      }
      Loop.push_back(Loop.front());
      FixedPolylines.push_back(Loop);
    }
    // The region is a whole closed
    // component, quadrangulated from scratch without fixed sides (the border
    // checks below are trivially satisfied: 0 edges)
    if (FixedPolylines.empty())
      std::cout << "Region remesh: the region has no border, quadrangulated "
                   "from scratch"
                << std::endl;
    int NumBorderE = 0;
    for (size_t i = 0; i < FixedPolylines.size(); i++)
      NumBorderE += (int)FixedPolylines[i].size() - 1;
    if ((NumBorderE % 2) != 0)
      std::cout << "Region remesh: WARNING odd number of border edges "
                << NumBorderE << ", no quad mesh can match it" << std::endl;

    // old feature edges inside the region (or on its border)
    std::set<std::pair<int, int>> OldFeatSet;
    for (size_t i = 0; i < QuadFeatures.size(); i++)
      OldFeatSet.insert(std::make_pair(
          std::min(QuadFeatures[i].first, QuadFeatures[i].second),
          std::max(QuadFeatures[i].first, QuadFeatures[i].second)));
    std::vector<std::pair<int, int>> RegFeat;
    {
      std::set<std::pair<int, int>> Added;
      for (size_t i = 0; i < RegPolyFaces.size(); i++)
        for (size_t j = 0; j < RegPolyFaces[i].size(); j++) {
          int V0 = RegPolyFaces[i][j];
          int V1 = RegPolyFaces[i][(j + 1) % RegPolyFaces[i].size()];
          std::pair<int, int> KeyO(std::min(OrigV[V0], OrigV[V1]),
                                   std::max(OrigV[V0], OrigV[V1]));
          if (OldFeatSet.count(KeyO) == 0)
            continue;
          std::pair<int, int> Key(std::min(V0, V1), std::max(V0, V1));
          if (Added.count(Key) > 0)
            continue;
          Added.insert(Key);
          RegFeat.push_back(Key);
        }
    }

    // 2. TRIANGLE REMESH (border fixed) --------------------------------------
    // target surface for the reprojections of the embedding
    // The region itself
    const std::vector<CoordType> &TargetVert = RegVert;
    const std::vector<std::vector<int>> &TargetFaces = RegTriFaces;

    std::vector<CoordType> TriVert = RegVert;
    std::vector<std::vector<int>> TriFaces = RegTriFaces;
    ScalarType TargetQuadE = (Param.TargetQuadEdge > 0)
                                 ? Param.TargetQuadEdge
                                 : AvgQuadE * Param.QuadSizeMult;

    TriRemParam<ScalarType> TRemParam;
    TRemParam.fix_borders = true;
    TRemParam.TargetL = TargetQuadE * Param.TriEdgeRatio;
    TRemParam.steps = Param.TriRemeshSteps;
    TRemParam.IsFaceEFeature.assign(TriFaces.size(), std::vector<bool>(3, false));
    TRemParam.IsVertCorner.assign(TriVert.size(), false);

    std::vector<std::vector<CoordType>> ToEmbed;
    if (Param.reembed_old_features) {
      std::vector<std::vector<CoordType>> OldLines;
      ChainFeatures(RegVert, RegFeat, OldLines);
      ToEmbed.insert(ToEmbed.end(), OldLines.begin(), OldLines.end());
    } else {
      // old features kept by the remesher as feature edges, their extremes
      // and sharp turns as corners
      std::set<std::pair<int, int>> RegFeatSet(RegFeat.begin(), RegFeat.end());
      for (size_t i = 0; i < TriFaces.size(); i++)
        for (size_t j = 0; j < 3; j++) {
          int V0 = TriFaces[i][j], V1 = TriFaces[i][(j + 1) % 3];
          if (RegFeatSet.count(std::make_pair(std::min(V0, V1),
                                              std::max(V0, V1))) > 0)
            TRemParam.IsFaceEFeature[i][j] = true;
        }
      std::vector<bool> IsCorner;
      FeatureCornersOf(TriVert, RegFeat, Param.FeatureCornerAngle, IsCorner);
      for (size_t i = 0; i < IsCorner.size(); i++)
        if (IsCorner[i])
          TRemParam.IsVertCorner[i] = true;
    }
    // -------------------------------------
    // the steps after a triangle remesh, as functions (used twice)
    // the border edges preserved by the remesher
    auto PreserveBorderEdges = [&](TriRemParam<ScalarType> &RP) {
      RP.PreserveEdges.clear();
      std::map<std::pair<int, int>, int> EdgeCount;
      for (size_t i = 0; i < TriFaces.size(); i++)
        for (size_t j = 0; j < 3; j++) {
          int V0 = TriFaces[i][j], V1 = TriFaces[i][(j + 1) % 3];
          EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      typename std::map<std::pair<int, int>, int>::iterator ItE;
      for (ItE = EdgeCount.begin(); ItE != EdgeCount.end(); ItE++)
        if ((*ItE).second == 1)
          RP.PreserveEdges.insert(MakeEdgePos(TriVert[(*ItE).first.first],
                                              TriVert[(*ItE).first.second]));
    };
    // the feature data of the remesher sized as the mesh
    auto FixFeatureSizes = [&](TriRemParam<ScalarType> &RP) {
      if (RP.IsFaceEFeature.size() != TriFaces.size()) {
        std::cout << "Region remesh: WARNING features lost by the remesher"
                  << std::endl;
        RP.IsFaceEFeature.assign(TriFaces.size(), std::vector<bool>(3, false));
      }
      if (RP.IsVertCorner.size() != TriVert.size())
        RP.IsVertCorner.resize(TriVert.size(), false);
    };
    // (exact border) every border vertex exactly on the fixed polylines (the
    // old border vertices are already there; a vertex added on a border edge
    // goes back onto it)
    auto SnapBorder = [&]() {
      std::vector<std::vector<int>> Loops;
      getBorderSequncesVert(TriFaces, Loops);
      int NumMoved = 0;
      for (size_t l = 0; l < Loops.size(); l++)
        for (size_t k = 0; k < Loops[l].size(); k++) {
          CoordType &P = TriVert[Loops[l][k]];
          if (FixedPos.count(P) > 0)
            continue;
          ScalarType BestD = std::numeric_limits<ScalarType>::max();
          CoordType BestQ = P;
          for (size_t j = 0; j < FixedPolylines.size(); j++)
            for (size_t q = 0; q + 1 < FixedPolylines[j].size(); q++) {
              const CoordType &A = FixedPolylines[j][q];
              CoordType D = FixedPolylines[j][q + 1] - A;
              ScalarType L2 = Geo::Dot(D, D);
              ScalarType t = (L2 > 0) ? Geo::Dot(P - A, D) / L2 : 0;
              t = std::max((ScalarType)0, std::min((ScalarType)1, t));
              CoordType Q = A + D * t;
              ScalarType Dist = (P - Q).Norm();
              if (Dist < BestD) {
                BestD = Dist;
                BestQ = Q;
              }
            }
          if (!(BestQ == P)) {
            P = BestQ;
            NumMoved++;
          }
        }
      if (NumMoved > 0)
        std::cout << "Region remesh: " << NumMoved
                  << " border vertices put back on the border" << std::endl;
    };

    // Index in ToEmbed of the first new line
    size_t FirstNewLine = ToEmbed.size();
    if (EmbeddedLines != NULL)
      EmbeddedLines->assign(NewFeatureLines.size(), std::vector<EdgePos>());
    ToEmbed.insert(ToEmbed.end(), NewFeatureLines.begin(),
                   NewFeatureLines.end());

    // The lines are embedded in the
    // triangulated region, then the remesh keeps them as features and makes
    // the triangles around them as uniform as elsewhere
    // 3. EMBED THE FEATURE LINES -----------------------------------------------
    {
      MeshEmbedPolyline<ScalarType> Embedder(
          TriVert, TriFaces, TRemParam.IsFaceEFeature, TRemParam.IsVertCorner,
          TargetVert, TargetFaces);
      // Sample_ratio is a fraction of the diagonal
      // of the mesh given to the embedder, here the region: on a small region
      // the samples (and the refinement under the line) became much denser
      // than the triangles. The spacing is the triangle target edge instead.
      EmbedPolylineParam<ScalarType> EParam = Param.EmbedParam;
      {
        Geo::Box3<ScalarType> EBox;
        EBox.Init(TriVert);
        if (EBox.Diag() > 0)
          EParam.sample_ratio = TRemParam.TargetL / EBox.Diag();
      }
      for (size_t i = 0; i < ToEmbed.size(); i++) {
        if (ToEmbed[i].size() < 2)
          continue;
        std::vector<std::pair<int, int>> NewEdges;
        bool done = Embedder.Embed(ToEmbed[i], EParam, NewEdges);
        // The embedded line (positions on the region surface)
        if (done && (EmbeddedLines != NULL) && (i >= FirstNewLine))
          for (size_t k = 0; k < NewEdges.size(); k++)
            (*EmbeddedLines)[i - FirstNewLine].push_back(
                MakeEdgePos(TriVert[NewEdges[k].first],
                            TriVert[NewEdges[k].second]));
        if (!done)
          std::cout << "Region remesh: WARNING line " << i << " not embedded"
                    << std::endl;
        // A new line that cannot be embedded fails the
        // remesh (the caller can try a bigger region) instead of a remesh
        // without it
        if ((!done) && (i >= FirstNewLine)) {
          LastMessage() = "the line cannot be embedded in the region";
          return false;
        }
      }
    }
    Lap("line embedding", TriVert.size(), TriFaces.size());

    // 2b. TRIANGLE REMESH (after the embedding) -------------------------------
    // -1 if a step broke the manifoldness
    TRemParam.check_manifold = true;
    // The border edges are preserved (as the shared
    // edges of TriRemesh with OnlyFaces): fix_borders alone keeps only the
    // border vertices, a border edge could still be split and the new vertex
    // smoothed off the border
    PreserveBorderEdges(TRemParam);
    if (TriRemesh<ScalarType>(TriVert, TriFaces, TRemParam) < 0) {
      LastMessage() = "the triangle remesh broke the manifoldness";
      return false;
    }
    FixFeatureSizes(TRemParam);
    SnapBorder();

    if (Param.save_debug)
      WriteOBJ("region_remesh.obj", TriVert, TriFaces);
    Lap("tri remesh after embedding", TriVert.size(),
        TriFaces.size());

    std::vector<std::pair<int, int>> TriFeatures;
    Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
        TriFaces, TRemParam.IsFaceEFeature, TriFeatures);

    // The triangle features as segments
    std::vector<EdgePos> TriFeatSegs;
    for (size_t i = 0; i < TriFeatures.size(); i++)
      TriFeatSegs.push_back(MakeEdgePos(TriVert[TriFeatures[i].first],
                                        TriVert[TriFeatures[i].second]));
    LastTriFeatures() = TriFeatSegs;

    auto TTri = std::chrono::high_resolution_clock::now();

    // 4. FIELD -----------------------------------------------------------------
    std::vector<Field::CrossF<ScalarType>> FaceCurv;
    // The border edges of the region are always
    // features for the quadrangulator, as the mesh border in the command
    // line: the patch corners come from the features (MeshFeatures::Corners),
    // so the turns of the region border become patch corners. Added before
    // the field (the refinements for the field and the tracing see them too)
    {
      std::map<std::pair<int, int>, int> EdgeCount;
      for (size_t i = 0; i < TriFaces.size(); i++)
        for (size_t j = 0; j < 3; j++) {
          int V0 = TriFaces[i][j], V1 = TriFaces[i][(j + 1) % 3];
          EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      std::set<std::pair<int, int>> FeatSet;
      for (size_t i = 0; i < TriFeatures.size(); i++)
        FeatSet.insert(std::make_pair(
            std::min(TriFeatures[i].first, TriFeatures[i].second),
            std::max(TriFeatures[i].first, TriFeatures[i].second)));
      int NumAdded = 0;
      typename std::map<std::pair<int, int>, int>::iterator ItE;
      for (ItE = EdgeCount.begin(); ItE != EdgeCount.end(); ItE++) {
        if (((*ItE).second != 1) || (FeatSet.count((*ItE).first) > 0))
          continue;
        TriFeatures.push_back((*ItE).first);
        NumAdded++;
      }
      std::cout << "Region remesh: " << NumAdded
                << " border edges added as features" << std::endl;
    }
    // (the border edges of the region are never split)
    ComputeField(TriVert, TriFaces, TriFeatures, FaceCurv, Param.FieldNRing,
                 Param.field_refinements, true);
    // The triangle mesh the quadrangulator works on (after the
    // refinements of the field) and its features as lines
    if (Param.save_debug) {
      WriteOBJ("region_field_tri.obj", TriVert, TriFaces);
      std::ofstream FL("region_features.obj");
      for (size_t i = 0; i < TriVert.size(); i++)
        FL << "v " << TriVert[i].X << " " << TriVert[i].Y << " "
           << TriVert[i].Z << "\n";
      for (size_t i = 0; i < TriFeatures.size(); i++)
        FL << "l " << TriFeatures[i].first + 1 << " "
           << TriFeatures[i].second + 1 << "\n";
    }

    Lap("field (with refinements)", TriVert.size(),
        TriFaces.size());

    // 5. QUADRANGULATION WITH FIXED BORDER ---------------------------------------
    std::vector<CoordType> NewQuadVert;
    std::vector<std::vector<int>> NewQuadFaces;
    std::vector<std::pair<int, int>> NewQuadFeatures;
    std::vector<int> NewQuadPatch;
    std::vector<int> NewQuadPatchCorner;
    int NumSamplesUsed = 0;
    auto TField = std::chrono::high_resolution_clock::now();
    {
      MeshQuadrangulator<ScalarType> MQuad(TriVert, TriFaces, FaceCurv,
                                           NewQuadVert, NewQuadFaces);
      MQuad.Init();
      if (!MQuad.MatchPreconditions()) {
        std::cout << "Region remesh: remeshed region not manifold" << std::endl;
        LastMessage() = "remeshed region not manifold (after embedding the "
                        "lines and refining for the field)";
        return false;
      }
      Geo::Box3<ScalarType> RegBox;
      RegBox.Init(TriVert);
      MQuad.QuadSize = TargetQuadE / RegBox.Diag();
      MQuad.QuadSmoothSteps = Param.QuadSmoothSteps;
      MQuad.CCAbility = Param.CCAbility;
      MQuad.RegVsIso = Param.RegVsIso;
      int RegSamples = Param.NumSamples;
      if (Param.scale_samples) {
        ScalarType AllA = PolyArea(QuadVert, QuadFaces);
        ScalarType RegA = PolyArea(RegVert, RegPolyFaces);
        if (AllA > 0)
          RegSamples = (int)std::floor(0.5 + Param.NumSamples * RegA / AllA);
        RegSamples = std::max(RegSamples, Param.MinSamples);
        RegSamples = std::min(RegSamples, Param.NumSamples);
      }
      NumSamplesUsed = RegSamples;
      MQuad.num_samples = RegSamples;
      TField = std::chrono::high_resolution_clock::now();
      MQuad.reliable_smooth = Param.reliable_smooth;
      // The region border is the constraint;
      // the sides on it are fixed to the old border vertices (the polylines
      // are the whole border, AllBorderSides = true), the corners on it are
      // snapped to a border vertex. The new lines are NOT fixed: they are
      // only features, their sampling is free.
      const ScalarType TolR = Param.FixedSideTolRatio;
      LastFixInfo() = "no border constraint";
      MQuad.ConstrainSides =
          [FixedPolylines, TolR](Geo::SideDiscretization<ScalarType> &SD,
                                 std::map<Geo::Point3<ScalarType>,
                                          Geo::Point3<ScalarType>> &Snap) {
            Geo::FixedPolylineSides<ScalarType>::Apply(SD, FixedPolylines,
                                                       TolR, true, Snap);
            // Border sides fixed / all border sides
            int NB = 0, NBF = 0;
            for (int i = 0; i < SD.NumVars(); i++)
              if (SD.IsBorderVariable(i)) {
                NB++;
                if (SD.IsFixedVariable(i))
                  NBF++;
              }
            LastFixInfo() = "border sides fixed " + std::to_string(NBF) +
                            " / " + std::to_string(NB) + ", snapped corners " +
                            std::to_string(Snap.size());
            std::cout << "Region remesh: " << LastFixInfo() << std::endl;
          };
      MQuad.fix_border_vert = true;
      // The dual vertices on the region border go on old border vertices
      // (the corners there are then old vertices, not moved afterwards)
      MQuad.BorderAnchorPos = FixedPos;
      // Only the old features and the new lines,
      // no prolongation along the patch sides
      MQuad.exact_feature_samples = true;
      // the ends of the sketched lines stay where the line ends
      MQuad.fix_open_line_ends = true;
      MQuad.ExtractQuad(TriFeatures);
      MQuad.GetQuadFeatures(NewQuadFeatures);
      // No prolongation of the open lines
      if (Param.trim_features)
        TrimFeaturesToLines(NewQuadVert, NewQuadFeatures, TriFeatSegs,
                            TargetQuadE * (ScalarType)0.25);
      MQuad.GetQuadPatch(NewQuadPatch);
      MQuad.GetQuadPatchCorner(NewQuadPatchCorner);
      // The primal and dual patch layouts of the region
      // (as SaveLayouts of the command line) and the dual patches by number
      // of sides (a dual patch with other than 4 sides holds a singularity)
      {
        std::vector<CoordType> PrimalV, DualV;
        std::vector<std::vector<int>> PrimalF, DualF;
        MQuad.GetPrimalDualMesh(PrimalV, PrimalF, DualV, DualF);
        // The patches of the quads (NewQuadPatch) by
        // number of corners (a vertex with one quad of the patch)
        std::map<int, std::map<int, int>> PatchVCount;
        for (size_t i = 0; i < NewQuadFaces.size(); i++)
          if ((i < NewQuadPatch.size()) && (NewQuadPatch[i] >= 0))
            for (size_t j = 0; j < NewQuadFaces[i].size(); j++)
              PatchVCount[NewQuadPatch[i]][NewQuadFaces[i][j]]++;
        // (corners as QuadLayoutRefit: less than half of the quads around a
        // vertex inside, one quad or shared on the border)
        std::vector<int> VQ(NewQuadVert.size(), 0);
        std::map<std::pair<int, int>, int> EC;
        for (size_t i = 0; i < NewQuadFaces.size(); i++)
          for (size_t j = 0; j < NewQuadFaces[i].size(); j++) {
            VQ[NewQuadFaces[i][j]]++;
            int V0 = NewQuadFaces[i][j];
            int V1 = NewQuadFaces[i][(j + 1) % NewQuadFaces[i].size()];
            EC[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
          }
        std::vector<bool> BV(NewQuadVert.size(), false);
        for (std::map<std::pair<int, int>, int>::iterator ItE = EC.begin();
             ItE != EC.end(); ItE++)
          if ((*ItE).second == 1) {
            BV[(*ItE).first.first] = true;
            BV[(*ItE).first.second] = true;
          }
        std::map<int, int> BySides;
        std::map<int, std::map<int, int>>::iterator ItP;
        for (ItP = PatchVCount.begin(); ItP != PatchVCount.end(); ItP++) {
          int NumC = 0;
          std::map<int, int>::iterator ItV;
          for (ItV = (*ItP).second.begin(); ItV != (*ItP).second.end(); ItV++) {
            int KP = (*ItV).second, N = VQ[(*ItV).first];
            if (BV[(*ItV).first] ? ((KP == 1) || (KP < N)) : (2 * KP < N))
              NumC++;
          }
          BySides[NumC]++;
        }
        std::string Hist;
        std::map<int, int>::iterator ItS;
        for (ItS = BySides.begin(); ItS != BySides.end(); ItS++)
          Hist += " " + std::to_string((*ItS).second) + "x" +
                  std::to_string((*ItS).first) + "-cornered";
        std::cout << "Region remesh: " << PatchVCount.size()
                  << " layout patches:" << Hist << std::endl;
        LastLayout() =
            std::to_string(PatchVCount.size()) + " layout patches:" + Hist;
        if (Param.save_debug) {
          WriteOBJ("region_primal_layout.obj", PrimalV, PrimalF);
          WriteOBJ("region_dual_layout.obj", DualV, DualF);
        }
      }
      // The quads are kept, only a
      // warning; the border check below still refuses a region that does
      // not match the old border
      if (MQuad.HasProblematicPatches())
        std::cout << "Region remesh: WARNING problematic patches in the "
                     "quadrangulation"
                  << std::endl;
    }
    if (Param.save_debug)
      WriteOBJ("region_quad.obj", NewQuadVert, NewQuadFaces);

    // the new border must be exactly the old one: as many border edges, and
    // as many border vertices, all old ones (counted on the edges, not by
    // walking the border loops: a broken quadrangulation can have a
    // non-manifold border, where the walk never ends)
    {
      std::map<std::pair<int, int>, int> EdgeCount;
      for (size_t i = 0; i < NewQuadFaces.size(); i++)
        for (size_t j = 0; j < NewQuadFaces[i].size(); j++) {
          int V0 = NewQuadFaces[i][j];
          int V1 = NewQuadFaces[i][(j + 1) % NewQuadFaces[i].size()];
          EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      int NewBorderE = 0;
      std::set<int> NewBorderV;
      typename std::map<std::pair<int, int>, int>::iterator ItE;
      for (ItE = EdgeCount.begin(); ItE != EdgeCount.end(); ItE++) {
        if ((*ItE).second != 1)
          continue;
        NewBorderE++;
        NewBorderV.insert((*ItE).first.first);
        NewBorderV.insert((*ItE).first.second);
      }
      int Matched = 0;
      std::set<int>::iterator ItV;
      for (ItV = NewBorderV.begin(); ItV != NewBorderV.end(); ItV++)
        if (FixedPos.count(NewQuadVert[*ItV]) > 0)
          Matched++;
      if ((NewBorderE != NumBorderE) || (Matched != NumBorderE) ||
          ((int)NewBorderV.size() != NumBorderE)) {
        std::cout << "Region remesh: border not matched (edges " << NewBorderE
                  << " / " << NumBorderE << ", vertices " << Matched << " / "
                  << NumBorderE << "), mesh unchanged" << std::endl;
        LastMessage() = "new border not matching the old one (edges " +
                        std::to_string(NewBorderE) + " / " +
                        std::to_string(NumBorderE) + ", vertices " +
                        std::to_string(Matched) + " / " +
                        std::to_string(NumBorderE) + "; " + LastFixInfo() +
                        ")";
        return false;
      }
    }

    // 6. MERGE -----------------------------------------------------------------
    // feature edges by position: the old ones (those left outside the region
    // survive) and the new ones
    std::set<EdgePos> FeatPos;
    std::set<int> RegionSet(RegionFaces.begin(), RegionFaces.end());
    for (size_t i = 0; i < QuadFeatures.size(); i++)
      FeatPos.insert(MakeEdgePos(QuadVert[QuadFeatures[i].first],
                                 QuadVert[QuadFeatures[i].second]));
    // The new features minus the border edges
    {
      std::map<std::pair<int, int>, int> EdgeCount;
      for (size_t i = 0; i < NewQuadFaces.size(); i++)
        for (size_t j = 0; j < NewQuadFaces[i].size(); j++) {
          int V0 = NewQuadFaces[i][j];
          int V1 = NewQuadFaces[i][(j + 1) % NewQuadFaces[i].size()];
          EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      for (size_t i = 0; i < NewQuadFeatures.size(); i++) {
        int V0 = NewQuadFeatures[i].first, V1 = NewQuadFeatures[i].second;
        std::pair<int, int> Key(std::min(V0, V1), std::max(V0, V1));
        if ((EdgeCount.count(Key) > 0) && (EdgeCount[Key] == 1))
          continue;
        FeatPos.insert(MakeEdgePos(NewQuadVert[V0], NewQuadVert[V1]));
      }
    }

    std::vector<CoordType> MergedVert = QuadVert;
    std::vector<std::vector<int>> MergedFaces;
    for (size_t i = 0; i < QuadFaces.size(); i++)
      if (RegionSet.count((int)i) == 0)
        MergedFaces.push_back(QuadFaces[i]);
    size_t FirstNewF = MergedFaces.size();
    AppendMesh(MergedVert, MergedFaces, NewQuadVert, NewQuadFaces);
    RemoveDuplicatedVert(MergedVert, MergedFaces);
    RemoveUnfererencedVert(MergedVert, MergedFaces);

    // features of the merged mesh: the edges whose positions are features
    // (old features inside the region are gone with their edges)
    std::vector<std::pair<int, int>> MergedFeatures;
    {
      std::set<std::pair<int, int>> Added;
      for (size_t i = 0; i < MergedFaces.size(); i++)
        for (size_t j = 0; j < MergedFaces[i].size(); j++) {
          int V0 = MergedFaces[i][j];
          int V1 = MergedFaces[i][(j + 1) % MergedFaces[i].size()];
          if (FeatPos.count(MakeEdgePos(MergedVert[V0], MergedVert[V1])) == 0)
            continue;
          std::pair<int, int> Key(std::min(V0, V1), std::max(V0, V1));
          if (Added.count(Key) > 0)
            continue;
          Added.insert(Key);
          MergedFeatures.push_back(Key);
        }
    }

    Lap("quadrangulation + merge", NewQuadVert.size(),
        NewQuadFaces.size());

    if ((Param.SeamSmoothSteps > 0) && (FirstNewF < MergedFaces.size())) {
      // target surface: the old quad mesh triangulated
      // The true surface if given
      // Built below, around the band
      // movable: the new faces plus SeamRings rings of faces
      std::vector<std::vector<int>> VF(MergedVert.size());
      for (size_t i = 0; i < MergedFaces.size(); i++)
        for (size_t j = 0; j < MergedFaces[i].size(); j++)
          VF[MergedFaces[i][j]].push_back((int)i);
      std::vector<bool> InBand(MergedFaces.size(), false);
      std::vector<int> Front;
      for (size_t i = FirstNewF; i < MergedFaces.size(); i++) {
        InBand[i] = true;
        Front.push_back((int)i);
      }
      for (int r = 0; r < Param.SeamRings; r++) {
        std::vector<int> NewFront;
        for (size_t k = 0; k < Front.size(); k++)
          for (size_t j = 0; j < MergedFaces[Front[k]].size(); j++) {
            int V = MergedFaces[Front[k]][j];
            for (size_t l = 0; l < VF[V].size(); l++)
              if (!InBand[VF[V][l]]) {
                InBand[VF[V][l]] = true;
                NewFront.push_back(VF[V][l]);
              }
          }
        Front = NewFront;
      }
      // The surface around the band only
      std::vector<CoordType> TargetQV;
      std::vector<std::vector<int>> TargetQF;
      {
        std::vector<CoordType> BandPts;
        for (size_t i = 0; i < MergedFaces.size(); i++)
          if (InBand[i])
            for (size_t j = 0; j < MergedFaces[i].size(); j++)
              BandPts.push_back(MergedVert[MergedFaces[i][j]]);
        ScalarType Margin = 2 * std::max(AvgQuadE, Param.TargetQuadEdge);
        // (the old quads as the surface)
        std::vector<std::vector<int>> AllTri;
        Triangulate(QuadVert, QuadFaces, AllTri);
        CropToBox(QuadVert, AllTri, BandPts, Margin, TargetQV, TargetQF);
      }
      // a vertex moves if all its faces are in the band, it is not on a
      // feature and not on the mesh border
      std::vector<bool> Movable(MergedVert.size(), true);
      for (size_t v = 0; v < VF.size(); v++)
        for (size_t l = 0; l < VF[v].size(); l++)
          if (!InBand[VF[v][l]])
            Movable[v] = false;
      for (size_t i = 0; i < MergedFeatures.size(); i++) {
        Movable[MergedFeatures[i].first] = false;
        Movable[MergedFeatures[i].second] = false;
      }
      std::map<std::pair<int, int>, int> EdgeCount;
      for (size_t i = 0; i < MergedFaces.size(); i++)
        for (size_t j = 0; j < MergedFaces[i].size(); j++) {
          int V0 = MergedFaces[i][j];
          int V1 = MergedFaces[i][(j + 1) % MergedFaces[i].size()];
          EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      typename std::map<std::pair<int, int>, int>::iterator ItE;
      for (ItE = EdgeCount.begin(); ItE != EdgeCount.end(); ItE++)
        if ((*ItE).second == 1) {
          Movable[(*ItE).first.first] = false;
          Movable[(*ItE).first.second] = false;
        }
      Geo::SmoothParam<ScalarType> SParam;
      SParam.Clear();
      for (size_t v = 0; v < Movable.size(); v++)
        if (!Movable[v])
          SParam.FixedVert.push_back((int)v);
      SParam.Features = MergedFeatures;
      SParam.ReprojVert = &TargetQV;
      SParam.ReprojFaces = &TargetQF;
      SParam.Damp = 0.5;
      // (no final projection: it made the region noisier)
      SParam.post_project = false;
      SParam.NumIte = Param.SeamSmoothSteps;
      SParam.SType =
          Param.reliable_smooth ? Geo::SMBestPoly : Geo::SMLaplacian;
      SParam.line_border = true;
      SParam.check_AspectR = false;
      SParam.line_neight_features = false;
      SParam.reproject_voronoi_partition = Param.reliable_smooth;
      SParam.check_multi_pass_fold = false;
      Geo::SmoothLaplacian<ScalarType>(MergedVert, MergedFaces, SParam);
    }

    Lap("seam smoothing", MergedVert.size(),
        MergedFaces.size());

    QuadVert = MergedVert;
    QuadFaces = MergedFaces;
    QuadFeatures = MergedFeatures;
    // The kept faces keep their region, the new ones
    // get a new id
    if (History != NULL)
      History->AfterRemesh(RegionSet, QuadFaces.size(),
                           QuadFaces.size() - FirstNewF);
    if (NewRegionFaces != NULL) {
      NewRegionFaces->clear();
      for (size_t i = FirstNewF; i < QuadFaces.size(); i++)
        NewRegionFaces->push_back((int)i);
    }
    // layout patch of the new faces (same order as the new quads)
    if (NewRegionPatchCorner != NULL) {
      *NewRegionPatchCorner = NewQuadPatchCorner;
      NewRegionPatchCorner->resize(NewQuadFaces.size(), 0);
    }
    if (NewRegionPatch != NULL) {
      *NewRegionPatch = NewQuadPatch;
      NewRegionPatch->resize(NewQuadFaces.size(), -1);
    }
    std::cout << "Region remesh: " << RegionFaces.size() << " faces replaced by "
              << NewQuadFaces.size() << std::endl;
    {
      auto TEnd = std::chrono::high_resolution_clock::now();
      auto Sec = [](const std::chrono::high_resolution_clock::time_point &A,
                    const std::chrono::high_resolution_clock::time_point &B) {
        return std::chrono::duration<double>(B - A).count();
      };
      char buf[512];
      snprintf(buf, sizeof(buf),
               "target edge %.4g, region avg edge before %.4g after %.4g, "
               "%d tracing samples; time tri %.2fs field %.2fs quad %.2fs",
               (double)TargetQuadE, (double)AvgQuadE,
               (double)AvgPolyEdge(NewQuadVert, NewQuadFaces), NumSamplesUsed,
               Sec(TStart, TTri), Sec(TTri, TField), Sec(TField, TEnd));
      LastInfo() = buf;
      LastInfo() += StageTimes;
      LastInfo() += "\nlayout: " + LastLayout();
      std::cout << "Region remesh: " << LastInfo() << std::endl;
    }
    return true;
  }
};

} // namespace Geo

#endif
