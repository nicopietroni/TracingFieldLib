#ifndef QUAD_LAYOUT_REFIT
#define QUAD_LAYOUT_REFIT

// The patch layout of a quad mesh rebuilt as a graph and its side
// subdivisions optimized again, globally, exactly as MeshQuadrangulator does
// (same target edge, same steps): where the layout did not change the same
// subdivisions come out; only the patches whose sides change are tessellated
// again, the rest of the quad mesh stays identical. After a region remesh the region border subdivisions were
// fixed: where a new feature line splits the region into two parts that need
// different subdivisions at its two ends (e.g. 2|3 at the top and 3|2 at the
// bottom) singularities were forced inside the region; here the old patches
// along the strips can change too and absorb the difference.
// The layout: the patches are the labels of the quads (QuadPatch, -1 =
// unknown), their corners the boundary vertices with one quad of the patch
// around them. The quads are triangulated (two triangles each, the patch
// sides are untouched) as the domain of SideDiscretization and of the
// tessellation. Patches with other than 4 corners and unknown patches are
// kept: their sides are fixed.

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <mesh_quadrangulation.h>
#include <patch_side_discretization.h>
#include <quadrangulation.h>
#include "quad_region_history.h"

namespace Geo {

template <class ScalarType> struct QuadLayoutRefitParam {
  // patches allowed to change: rings of patches around the patches of the
  // seed faces; -1 = the whole mesh
  int PatchRings = -1;
  // weight of the equality of the opposite sides against the target (the
  // current subdivisions), as MeshQuadrangulator::RegVsIso
  ScalarType RegVsIso = 0.5;
  // target quad edge (absolute), the one of the full quadrangulation: the
  // same optimization finds the same subdivisions where the layout did not
  // change; <= 0 = the average edge of the quad mesh
  ScalarType TargetEdge = -1;
  // Steps of smoothing of the patches tessellated again
  // (their vertices only: the features, the mesh border and the kept quads
  // stay), reprojected on the old quads; 0 = off
  int SmoothSteps = 10;
  bool reliable_smooth = true;
  // True = of the changed sides only the ones linked to the
  // seed patches are kept (the others back as before); false = every change
  // of the new solve kept, as the full quadrangulation would do
  bool keep_connected_only = true;
};

template <class ScalarType> class QuadLayoutRefit {
  typedef Geo::Point3<ScalarType> CoordType;
  typedef std::pair<CoordType, CoordType> EdgePos;

  static EdgePos MakeEdgePos(const CoordType &P0, const CoordType &P1) {
    return EdgePos(std::min(P0, P1), std::max(P0, P1));
  }

public:
  static std::string &LastMessage() {
    static std::string M;
    return M;
  }

  // QuadVert, QuadFaces, QuadPatch, QuadFeatures: changed in place.
  // SeedFaces: faces whose patches start the rings (e.g. the new faces of a
  // region remesh). OldFaceOf (optional): for every face of the result the
  // index of the same face before, -1 for a new face. Returns false (mesh
  // unchanged) if nothing changes or something fails (LastMessage says it).
  static bool Refit(std::vector<CoordType> &QuadVert,
                    std::vector<std::vector<int>> &QuadFaces,
                    std::vector<int> &QuadPatch,
                    std::vector<std::pair<int, int>> &QuadFeatures,
                    const std::vector<int> &SeedFaces,
                    const QuadLayoutRefitParam<ScalarType> &Param =
                        QuadLayoutRefitParam<ScalarType>(),
                    std::vector<int> *OldFaceOf = NULL,
                    QuadRegionHistory<ScalarType> *History = NULL,
                    std::vector<int> *QuadPatchCorner = NULL) {
    // QuadPatchCorner, if given (one mask per quad, bit
    // j = vertex j is a corner of the quad's patch), gives the corners of the
    // patches as the quadrangulator made them; it is updated with the
    // patches tessellated again. Without it the corners are guessed
    LastMessage().clear();
    size_t NF = QuadFaces.size();
    if (QuadPatch.size() != NF) {
      LastMessage() = "no patch layout";
      return false;
    }
    for (size_t i = 0; i < NF; i++)
      if (QuadFaces[i].size() != 4) {
        LastMessage() = "not a quad mesh";
        return false;
      }

    // 1. the domain: the quads triangulated, the patches -------------------
    std::vector<std::vector<int>> TriF;
    TriF.reserve(2 * NF);
    for (size_t i = 0; i < NF; i++) {
      const std::vector<int> &Q = QuadFaces[i];
      TriF.push_back({Q[0], Q[1], Q[2]});
      TriF.push_back({Q[0], Q[2], Q[3]});
    }
    std::map<int, std::vector<int>> LabelQuads;
    int NoLabel = 0;
    for (size_t i = 0; i < NF; i++)
      if (QuadPatch[i] >= 0)
        LabelQuads[QuadPatch[i]].push_back((int)i);
      else
        NoLabel++;
    if (NoLabel > 0)
      std::cout << "Layout refit: " << NoLabel
                << " quads without a layout patch, kept fixed" << std::endl;
    if ((QuadPatchCorner == NULL) || (QuadPatchCorner->size() != NF))
      std::cout << "Layout refit: corners of the patches not known, guessed"
                << std::endl;
    // A vertex is a corner of a patch if the
    // patch takes less than half of the quads around it (inside the mesh), or
    // on the mesh border if it has one quad (a corner of the mesh) or shares
    // the vertex with another patch: its angle there is less than a straight
    // one. One quad of the patch (the old rule) is the usual case;
    // the old rule missed the corners with two quads of the patch at an
    // irregular vertex (e.g. valence 3 on the border, valence 5 inside),
    // typical of the tessellations with different opposite sides, so exactly
    // the patches to improve were left out (and their neighbours fixed)
    std::vector<int> VertQuads(QuadVert.size(), 0);
    for (size_t i = 0; i < NF; i++)
      for (size_t j = 0; j < 4; j++)
        VertQuads[QuadFaces[i][j]]++;
    std::vector<bool> OnMeshBorderV(QuadVert.size(), false);
    {
      std::map<std::pair<int, int>, int> EC0;
      for (size_t i = 0; i < NF; i++)
        for (size_t j = 0; j < 4; j++) {
          int V0 = QuadFaces[i][j], V1 = QuadFaces[i][(j + 1) % 4];
          EC0[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      std::map<std::pair<int, int>, int>::iterator ItB;
      for (ItB = EC0.begin(); ItB != EC0.end(); ItB++)
        if ((*ItB).second == 1) {
          OnMeshBorderV[(*ItB).first.first] = true;
          OnMeshBorderV[(*ItB).first.second] = true;
        }
    }
    std::vector<int> Labels; // patch passed -> its label
    std::vector<std::vector<int>> PatchTris, PatchCorners;
    std::map<int, int> LabelIdx;
    // The patches left out (not 4 corners: e.g. an old
    // patch cut by the region border), with their quads
    int NumNot4 = 0, QuadsNot4 = 0;
    typename std::map<int, std::vector<int>>::iterator ItL;
    for (ItL = LabelQuads.begin(); ItL != LabelQuads.end(); ItL++) {
      std::map<int, int> VCount;
      for (size_t k = 0; k < (*ItL).second.size(); k++)
        for (size_t j = 0; j < 4; j++)
          VCount[QuadFaces[(*ItL).second[k]][j]]++;
      std::vector<int> Corners;
      std::map<int, int>::iterator ItV;
      // The known ones
      bool Known = (QuadPatchCorner != NULL) && (QuadPatchCorner->size() == NF);
      if (Known) {
        std::set<int> CS;
        for (size_t k = 0; k < (*ItL).second.size(); k++) {
          int Fq = (*ItL).second[k];
          for (size_t j = 0; j < 4; j++)
            if (((*QuadPatchCorner)[Fq] >> j) & 1)
              CS.insert(QuadFaces[Fq][j]);
        }
        Corners.assign(CS.begin(), CS.end());
      }
      for (ItV = VCount.begin(); (!Known) && (ItV != VCount.end()); ItV++) {
        int KP = (*ItV).second, N = VertQuads[(*ItV).first];
        // on the mesh border: a corner of the mesh (one quad) or a vertex
        // shared with another patch; inside: less than half of the quads
        bool corner = OnMeshBorderV[(*ItV).first] ? ((KP == 1) || (KP < N))
                                                  : (2 * KP < N);
        if (corner)
          Corners.push_back((*ItV).first);
      }
      if (Corners.size() != 4) {
        NumNot4++;
        QuadsNot4 += (int)(*ItL).second.size();
        // Which one and why
        std::cout << "Layout refit: patch " << (*ItL).first << " ("
                  << (*ItL).second.size() << " quads) has "
                  << Corners.size() << " corners ("
                  << (Known ? "known" : "guessed") << "), kept fixed"
                  << std::endl;
        continue; // kept as it is (its neighbours' sides fixed below)
      }
      std::vector<int> Tris;
      for (size_t k = 0; k < (*ItL).second.size(); k++) {
        Tris.push_back(2 * (*ItL).second[k]);
        Tris.push_back(2 * (*ItL).second[k] + 1);
      }
      LabelIdx[(*ItL).first] = (int)Labels.size();
      Labels.push_back((*ItL).first);
      PatchTris.push_back(Tris);
      PatchCorners.push_back(Corners);
    }
    if (Labels.empty()) {
      LastMessage() = "no 4 sided patch";
      return false;
    }

    // patches allowed to change
    std::vector<bool> Allowed(Labels.size(), true);
    if (Param.PatchRings >= 0) {
      Allowed.assign(Labels.size(), false);
      std::map<int, std::set<int>> Adj; // label -> adjacent labels
      {
        std::map<std::pair<int, int>, std::vector<int>> EF;
        for (size_t i = 0; i < NF; i++)
          for (size_t j = 0; j < 4; j++) {
            int V0 = QuadFaces[i][j], V1 = QuadFaces[i][(j + 1) % 4];
            EF[std::make_pair(std::min(V0, V1), std::max(V0, V1))].push_back(
                (int)i);
          }
        typename std::map<std::pair<int, int>, std::vector<int>>::iterator
            ItE;
        for (ItE = EF.begin(); ItE != EF.end(); ItE++) {
          if ((*ItE).second.size() != 2)
            continue;
          int L0 = QuadPatch[(*ItE).second[0]];
          int L1 = QuadPatch[(*ItE).second[1]];
          if (L0 != L1) {
            Adj[L0].insert(L1);
            Adj[L1].insert(L0);
          }
        }
      }
      std::set<int> In;
      for (size_t i = 0; i < SeedFaces.size(); i++)
        if ((SeedFaces[i] >= 0) && (SeedFaces[i] < (int)NF))
          In.insert(QuadPatch[SeedFaces[i]]);
      for (int r = 0; r < Param.PatchRings; r++) {
        std::set<int> Next = In;
        std::set<int>::iterator It;
        for (It = In.begin(); It != In.end(); It++)
          Next.insert(Adj[*It].begin(), Adj[*It].end());
        In = Next;
      }
      for (size_t p = 0; p < Labels.size(); p++)
        Allowed[p] = (In.count(Labels[p]) > 0);
    }

    // 2. the variables and their current subdivisions ------------------------
    Geo::SideDiscretization<ScalarType> SideD;
    SideD.InitFromLabels(QuadVert, TriF, PatchTris, PatchCorners,
                         QuadFeatures);
    int NumV = SideD.NumVars();
    if (NumV == 0) {
      LastMessage() = "no patch side";
      return false;
    }
    // InitFromLabels may discard some of the patches
    // passed (e.g. a patch whose corners do not split its border in 4
    // sides), so its patch indices are not the ones of Labels: the label of
    // each of its patches from its faces (two triangles per quad). Before,
    // a discarded patch shifted all the following ones: the wrong quads were
    // replaced and the merge failed ("edges not matched")
    std::vector<int> PsdLabel;
    std::vector<bool> PsdAllowed;
    {
      const std::vector<std::vector<int>> &PFI = SideD.GetPatchFacesIdx();
      for (size_t p = 0; p < PFI.size(); p++) {
        int L = PFI[p].empty() ? -1 : QuadPatch[PFI[p][0] / 2];
        PsdLabel.push_back(L);
        PsdAllowed.push_back((LabelIdx.count(L) > 0) ? Allowed[LabelIdx[L]]
                                                      : false);
      }
      if (PFI.size() != Labels.size()) {
        std::cout << "Layout refit: " << Labels.size() - PFI.size()
                  << " patches discarded by the side discretization (kept "
                     "as they are)"
                  << std::endl;
        NumNot4 += (int)(Labels.size() - PFI.size());
      }
    }
    // the edges on the mesh border
    std::set<std::pair<int, int>> MeshBorder;
    {
      std::map<std::pair<int, int>, int> EC;
      for (size_t i = 0; i < NF; i++)
        for (size_t j = 0; j < 4; j++) {
          int V0 = QuadFaces[i][j], V1 = QuadFaces[i][(j + 1) % 4];
          EC[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      std::map<std::pair<int, int>, int>::iterator ItE;
      for (ItE = EC.begin(); ItE != EC.end(); ItE++)
        if ((*ItE).second == 1)
          MeshBorder.insert((*ItE).first);
    }
    std::set<EdgePos> MeshBorderPos;
    {
      std::set<std::pair<int, int>>::iterator It;
      for (It = MeshBorder.begin(); It != MeshBorder.end(); It++)
        MeshBorderPos.insert(
            MakeEdgePos(QuadVert[(*It).first], QuadVert[(*It).second]));
    }
    // the current subdivisions (to know what changed) and samples
    std::vector<int> OldVal(NumV, 0);
    std::vector<std::vector<CoordType>> OldPts(NumV);
    int NumFixed = 0;
    // The variables on the mesh border (their new samples
    // are the new border edges)
    std::vector<bool> VarOnMeshBorder(NumV, false);
    for (int i = 0; i < NumV; i++) {
      OldPts[i] = SideD.PointsOfVariable(i);
      OldVal[i] = (int)OldPts[i].size() - 1;
      // fixed: a side of a kept patch (one patch only, not on the mesh
      // border), or of a patch not allowed to change
      bool Fix = false;
      if (SideD.IsBorderVariable(i)) {
        bool OnMeshBorder = (OldPts[i].size() >= 2);
        for (size_t k = 0; OnMeshBorder && (k + 1 < OldPts[i].size()); k++)
          OnMeshBorder = (MeshBorderPos.count(
                              MakeEdgePos(OldPts[i][k], OldPts[i][k + 1])) > 0);
        Fix = !OnMeshBorder;
        VarOnMeshBorder[i] = OnMeshBorder;
      }
      if (!Fix) {
        std::vector<std::array<int, 3>> PS;
        SideD.SideOfVariable(i, PS);
        for (size_t k = 0; k < PS.size(); k++)
          if ((PS[k][0] >= 0) && (PS[k][0] < (int)PsdAllowed.size()) &&
              !PsdAllowed[PS[k][0]])
            Fix = true;
      }
      if (Fix && (OldVal[i] >= 1)) {
        SideD.SetFixedVariable(i, OldVal[i], OldPts[i]);
        NumFixed++;
      }
    }
    std::cout << "Layout refit: " << Labels.size() << " patches, " << NumV
              << " sides (" << NumFixed << " fixed)" << std::endl;

    // 3. solve (as MeshQuadrangulator::FindSubdivisions) ---------------------
    ScalarType AvgE = 0;
    int NumE = 0;
    for (size_t i = 0; i < NF; i++)
      for (size_t j = 0; j < 4; j++) {
        AvgE += (QuadVert[QuadFaces[i][j]] - QuadVert[QuadFaces[i][(j + 1) % 4]])
                    .Norm();
        NumE++;
      }
    AvgE /= std::max(NumE, 1);
    SideD.ExactFeatureSamples = true;
    SideD.PrioritizeTJunction = true;
    // the same optimization as MeshQuadrangulator::FindSubdivisions
    ScalarType EdgeS = (Param.TargetEdge > 0) ? Param.TargetEdge : AvgE;
    SideD.FindDiscretization(EdgeS, Param.RegVsIso, true, true, true);
    if (SideD.NumOdd() > 0) {
      SideD.PMode = Geo::SideDiscretization<ScalarType>::EvenVar;
      SideD.FindDiscretization(EdgeS, Param.RegVsIso, true, true, true);
    }
    std::vector<int> Val;
    SideD.GetVariablesValues(Val);
    // The patches with a changed side,
    // linked when they share a changed side, form independent changes;
    // only the ones touching the patches of the seed are kept, the others
    // (e.g. a strip around the mesh changed by the rounding of the new
    // solve) go back exactly to their old subdivisions and vertices
    int NP = SideD.NumPatches();
    std::vector<int> Parent(NP);
    for (int p = 0; p < NP; p++)
      Parent[p] = p;
    std::function<int(int)> Root = [&](int p) -> int {
      while (Parent[p] != p) {
        Parent[p] = Parent[Parent[p]];
        p = Parent[p];
      }
      return p;
    };
    std::vector<std::vector<int>> VarPatches(NumV);
    for (int i = 0; i < NumV; i++) {
      if (Val[i] == OldVal[i])
        continue;
      std::vector<std::array<int, 3>> PS;
      SideD.SideOfVariable(i, PS);
      for (size_t k = 0; k < PS.size(); k++)
        if ((PS[k][0] >= 0) && (PS[k][0] < NP))
          VarPatches[i].push_back(PS[k][0]);
      for (size_t k = 1; k < VarPatches[i].size(); k++)
        Parent[Root(VarPatches[i][k])] = Root(VarPatches[i][0]);
    }
    std::set<int> KeepRoots;
    {
      std::set<int> SeedLabels;
      for (size_t i = 0; i < SeedFaces.size(); i++)
        if ((SeedFaces[i] >= 0) && (SeedFaces[i] < (int)NF))
          SeedLabels.insert(QuadPatch[SeedFaces[i]]);
      for (int p = 0; (p < NP) && (p < (int)PsdLabel.size()); p++)
        if (SeedLabels.count(PsdLabel[p]) > 0)
          KeepRoots.insert(Root(p));
    }
    std::set<int> Changed; // patches passed with a changed side (kept)
    int NumChangedV = 0, NumReverted = 0;
    for (int i = 0; i < NumV; i++) {
      // The cut can be switched off
      bool Revert = Param.keep_connected_only && (Val[i] != OldVal[i]) &&
                    !VarPatches[i].empty() &&
                    (KeepRoots.count(Root(VarPatches[i][0])) == 0);
      if ((Val[i] == OldVal[i]) || Revert) {
        // the unchanged (or reverted) sides keep exactly their old vertices
        if (!SideD.IsFixedVariable(i))
          SideD.SetFixedVariable(i, OldVal[i], OldPts[i]);
        if (Revert)
          NumReverted++;
        continue;
      }
      NumChangedV++;
      for (size_t k = 0; k < VarPatches[i].size(); k++)
        Changed.insert(VarPatches[i][k]);
    }
    if (NumReverted > 0)
      std::cout << "Layout refit: " << NumReverted
                << " changed sides not connected to the region, back as "
                   "before"
                << std::endl;
    if (Changed.empty()) {
      LastMessage() = "no side changed, " + std::to_string(NumNot4) +
                      " patches (" + std::to_string(QuadsNot4) +
                      " quads) not 4 sided, kept fixed";
      std::cout << "Layout refit: no side changed" << std::endl;
      return false;
    }
    // (the reverted sides: their values back to the fixed old ones)
    SideD.EnforceFixedValues(false);
    SideD.RecomputeSamples();
    // A mesh border side whose subdivision changed has new
    // edges, not in MeshBorderPos (the old ones): they are border edges too,
    // not cracks (before, they were counted as edges not matched)
    for (int i = 0; i < NumV; i++)
      if (VarOnMeshBorder[i] && (i < (int)SideD.SampleSidePoints.size()))
        for (size_t k = 0; k + 1 < SideD.SampleSidePoints[i].size(); k++)
          MeshBorderPos.insert(MakeEdgePos(SideD.SampleSidePoints[i][k],
                                           SideD.SampleSidePoints[i][k + 1]));

    // 4. the changed patches tessellated again ---------------------------------
    std::vector<std::vector<int>> PF, PC;
    std::vector<std::vector<std::vector<CoordType>>> SS;
    SideD.GetQuadrangulationData(PF, PC, SS);
    std::vector<std::vector<int>> PFc, PCc;
    std::vector<std::vector<std::vector<CoordType>>> SSc;
    std::vector<int> ChangedLabel;
    std::set<int>::iterator ItC;
    for (ItC = Changed.begin(); ItC != Changed.end(); ItC++) {
      if ((*ItC < 0) || (*ItC >= (int)PF.size()))
        continue;
      PFc.push_back(PF[*ItC]);
      PCc.push_back(PC[*ItC]);
      SSc.push_back(SS[*ItC]);
      ChangedLabel.push_back(PsdLabel[*ItC]);
    }
    std::vector<CoordType> NQV;
    std::vector<std::vector<int>> NQF;
    std::vector<int> NQP;
    int Problems =
        Geo::Quadrangulation<ScalarType>::QuadrangulateQuadPatchsWithSubs(
            QuadVert, TriF, PFc, PCc, SSc,
            MeshQuadrangulator<ScalarType>::GetSharedLUT(), NQV, NQF, NQP);
    if (Problems > 0) {
      LastMessage() = std::to_string(Problems) +
                      " patches could not be tessellated, mesh unchanged";
      std::cout << "Layout refit: " << LastMessage() << std::endl;
      return false;
    }
    // features of the new patches (sample segments on feature edges)
    std::vector<EdgePos> NewFeatData;
    SideD.GetFeaturesQuadrData(QuadVert, QuadFeatures, NewFeatData);

    // 5. merge: the kept quads, then the new ones (shared vertices by
    // position: the sides between a kept and a changed patch did not change)
    std::set<int> ChangedLabelSet(ChangedLabel.begin(), ChangedLabel.end());
    std::set<EdgePos> FeatPos;
    for (size_t i = 0; i < QuadFeatures.size(); i++)
      FeatPos.insert(MakeEdgePos(QuadVert[QuadFeatures[i].first],
                                 QuadVert[QuadFeatures[i].second]));
    for (size_t i = 0; i < NewFeatData.size(); i++)
      FeatPos.insert(MakeEdgePos(NewFeatData[i].first, NewFeatData[i].second));
    std::vector<CoordType> MV;
    std::vector<std::vector<int>> MF;
    std::vector<int> MP, MOld;
    std::vector<int> MC;
    std::map<CoordType, int> PosIdx;
    auto VertIdx = [&](const CoordType &P) -> int {
      typename std::map<CoordType, int>::iterator It = PosIdx.find(P);
      if (It != PosIdx.end())
        return (*It).second;
      int Idx = (int)MV.size();
      MV.push_back(P);
      PosIdx[P] = Idx;
      return Idx;
    };
    for (size_t i = 0; i < NF; i++) {
      if (ChangedLabelSet.count(QuadPatch[i]) > 0)
        continue;
      std::vector<int> F;
      for (size_t j = 0; j < 4; j++)
        F.push_back(VertIdx(QuadVert[QuadFaces[i][j]]));
      MF.push_back(F);
      MP.push_back(QuadPatch[i]);
      MOld.push_back((int)i);
      MC.push_back(((QuadPatchCorner != NULL) &&
                    (QuadPatchCorner->size() == NF))
                       ? (*QuadPatchCorner)[i]
                       : 0);
    }
    for (size_t i = 0; i < NQF.size(); i++) {
      std::vector<int> F;
      for (size_t j = 0; j < NQF[i].size(); j++)
        F.push_back(VertIdx(NQV[NQF[i][j]]));
      MF.push_back(F);
      int P = ((NQP[i] >= 0) && (NQP[i] < (int)ChangedLabel.size()))
                  ? ChangedLabel[NQP[i]]
                  : -1;
      MP.push_back(P);
      MOld.push_back(-1);
      // The ends of the sides of its patch
      int Mask = 0;
      if ((NQP[i] >= 0) && (NQP[i] < (int)SSc.size()))
        for (size_t j = 0; j < NQF[i].size(); j++)
          for (size_t k = 0; k < SSc[NQP[i]].size(); k++)
            if (!SSc[NQP[i]][k].empty() &&
                (SSc[NQP[i]][k].front() == NQV[NQF[i][j]]))
              Mask |= (1 << j);
      MC.push_back(Mask);
    }
    // the border of the changed part must match (each edge in two faces or
    // on the mesh border as before)
    {
      std::map<std::pair<int, int>, int> EC;
      for (size_t i = 0; i < MF.size(); i++)
        for (size_t j = 0; j < MF[i].size(); j++) {
          int V0 = MF[i][j], V1 = MF[i][(j + 1) % MF[i].size()];
          EC[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      int Open = 0;
      std::map<std::pair<int, int>, int>::iterator ItE;
      for (ItE = EC.begin(); ItE != EC.end(); ItE++)
        if (((*ItE).second == 1) &&
            (MeshBorderPos.count(MakeEdgePos(MV[(*ItE).first.first],
                                             MV[(*ItE).first.second])) == 0))
          Open++;
      if (Open > 0) {
        // Where the open edges are and why. An open edge
        // of a kept quad is a side of a kept patch whose new neighbour does
        // not have the same vertices there; an open edge of a new quad whose
        // ends are almost on kept vertices is a rounding difference; far
        // from them, a different subdivision (or the wrong patch replaced)
        {
          int OnKept = 0, OnNew = 0, NearKept = 0, NonManif = 0;
          ScalarType MaxNear = 0, AvgE2 = 0;
          std::set<int> KeptV;
          for (size_t i = 0; i < MF.size(); i++)
            if (MOld[i] >= 0)
              KeptV.insert(MF[i].begin(), MF[i].end());
          std::map<std::pair<int, int>, int> EFace; // edge -> a face
          for (size_t i = 0; i < MF.size(); i++)
            for (size_t j = 0; j < MF[i].size(); j++) {
              int V0 = MF[i][j], V1 = MF[i][(j + 1) % MF[i].size()];
              EFace[std::make_pair(std::min(V0, V1), std::max(V0, V1))] =
                  (int)i;
              AvgE2 += (MV[V0] - MV[V1]).Norm();
            }
          AvgE2 /= std::max((size_t)1, 4 * MF.size());
          for (ItE = EC.begin(); ItE != EC.end(); ItE++) {
            if ((*ItE).second > 2)
              NonManif++;
            if (((*ItE).second != 1) ||
                (MeshBorderPos.count(MakeEdgePos(MV[(*ItE).first.first],
                                                 MV[(*ItE).first.second])) > 0))
              continue;
            int Fc = EFace[(*ItE).first];
            if (MOld[Fc] >= 0) {
              OnKept++;
              continue;
            }
            OnNew++;
            int Ends[2] = {(*ItE).first.first, (*ItE).first.second};
            for (int e = 0; e < 2; e++) {
              if (KeptV.count(Ends[e]) > 0)
                continue;
              ScalarType Best = std::numeric_limits<ScalarType>::max();
              for (std::set<int>::iterator ItK = KeptV.begin();
                   ItK != KeptV.end(); ItK++)
                Best = std::min(Best, (MV[*ItK] - MV[Ends[e]]).Norm());
              if (Best < AvgE2 * (ScalarType)1e-3) {
                NearKept++;
                MaxNear = std::max(MaxNear, Best);
              }
            }
          }
          std::cout << "Layout refit: open edges " << OnKept
                    << " on kept quads, " << OnNew << " on new quads; "
                    << NearKept
                    << " ends of new quads almost on a kept vertex (max "
                       "distance "
                    << MaxNear << ", avg edge " << AvgE2 << "); " << NonManif
                    << " edges in more than 2 quads" << std::endl;
        }
        LastMessage() = std::to_string(Open) +
                        " edges not matched between kept and new patches, "
                        "mesh unchanged";
        std::cout << "Layout refit: " << LastMessage() << std::endl;
        return false;
      }
    }
    std::vector<std::pair<int, int>> MFeat;
    {
      std::set<std::pair<int, int>> Done;
      for (size_t i = 0; i < MF.size(); i++)
        for (size_t j = 0; j < MF[i].size(); j++) {
          int V0 = MF[i][j], V1 = MF[i][(j + 1) % MF[i].size()];
          std::pair<int, int> Key(std::min(V0, V1), std::max(V0, V1));
          if (Done.count(Key) > 0)
            continue;
          Done.insert(Key);
          if (FeatPos.count(MakeEdgePos(MV[V0], MV[V1])) > 0)
            MFeat.push_back(Key);
        }
    }
    std::cout << "Layout refit: " << NumChangedV << " sides changed, "
              << ChangedLabel.size() << " patches tessellated again ("
              << NF - (MF.size() - NQF.size()) << " quads -> " << NQF.size()
              << ")" << std::endl;
    LastMessage() = std::to_string(NumChangedV) + " sides changed, " +
                    std::to_string(ChangedLabel.size()) +
                    " patches tessellated again";
    LastMessage() += ", " + std::to_string(NumNot4) + " patches (" +
                     std::to_string(QuadsNot4) +
                     " quads) not 4 sided, kept fixed";
    // The
    // features and the mesh border are feature lines the vertices slide
    // along (reprojected on the old features and border), the corners of
    // the features fixed (not the open ends, as there); limited to the
    // quads tessellated again plus one ring of kept quads around them
    if ((Param.SmoothSteps > 0) && !NQF.empty()) {
      Geo::SmoothParam<ScalarType> SParam;
      SParam.Clear();
      SParam.Features = MFeat;
      std::vector<std::pair<int, int>> BorderEQuad;
      Geo::getBorderEdges(MF, BorderEQuad);
      SParam.Features.insert(SParam.Features.end(), BorderEQuad.begin(),
                             BorderEQuad.end());
      // the surface before: the old quads triangulated, their features and
      // border as the lines to reproject the feature vertices on
      SParam.ReprojVert = &QuadVert;
      SParam.ReprojFaces = &TriF;
      SParam.ReprojFeatures = QuadFeatures;
      std::vector<std::pair<int, int>> BorderETris;
      Geo::getBorderEdges(TriF, BorderETris);
      SParam.ReprojFeatures.insert(SParam.ReprojFeatures.end(),
                                   BorderETris.begin(), BorderETris.end());
      SParam.Damp = 0.5;
      SParam.post_project = true;
      SParam.NumIte = Param.SmoothSteps;
      SParam.SType = Param.reliable_smooth ? Geo::SMBestPoly : Geo::SMLaplacian;
      SParam.line_border = true;
      SParam.check_AspectR = false;
      SParam.line_neight_features = false;
      SParam.reproject_voronoi_partition = Param.reliable_smooth;
      SParam.check_multi_pass_fold = false;
      // fixed: the corners of the features, the open ends of the lines
      // included (a free end is pulled off its line, a kink), and the
      // vertices of the kept quads
      std::vector<int> Corners;
      Geo::MeshFeatures<ScalarType>::CornersByCount(MV, MFeat, Corners, true);
      std::vector<int> Valence1;
      Geo::MeshFeatures<ScalarType>::VerticesWithValence1(MFeat, Valence1);
      std::set<int> Fixed(Corners.begin(), Corners.end());
      Fixed.insert(Valence1.begin(), Valence1.end());
      // movable: the vertices of the new quads and of the kept quads
      // touching them (one ring, so the seam is smoothed too, as the seam
      // smoothing of the region remesh); the other kept vertices fixed
      std::vector<bool> Movable(MV.size(), false);
      for (size_t i = 0; i < MF.size(); i++)
        if (MOld[i] < 0)
          for (size_t j = 0; j < MF[i].size(); j++)
            Movable[MF[i][j]] = true;
      std::vector<bool> NewV = Movable;
      for (size_t i = 0; i < MF.size(); i++) {
        if (MOld[i] < 0)
          continue;
        bool touch = false;
        for (size_t j = 0; j < MF[i].size(); j++)
          touch = touch || NewV[MF[i][j]];
        if (touch)
          for (size_t j = 0; j < MF[i].size(); j++)
            Movable[MF[i][j]] = true;
      }
      for (size_t v = 0; v < MV.size(); v++)
        if (!Movable[v])
          Fixed.insert((int)v);
      SParam.FixedVert.assign(Fixed.begin(), Fixed.end());
      Geo::SmoothLaplacian<ScalarType>(MV, MF, SParam);
    }
    std::vector<int> OldQuadPatch = QuadPatch; // (for the region history)
    QuadVert = MV;
    QuadFaces = MF;
    QuadPatch = MP;
    QuadFeatures = MFeat;
    if (QuadPatchCorner != NULL)
      *QuadPatchCorner = MC;
    if (OldFaceOf != NULL)
      *OldFaceOf = MOld;
    // The kept faces keep their region, the faces
    // tessellated again take the one of their patch
    if (History != NULL)
      History->AfterRebuild(MOld, MP, &OldQuadPatch);
    return true;
  }
};

} // namespace Geo

#endif
