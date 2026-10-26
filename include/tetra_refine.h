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

// Refinement of tetrahedral meshes by edge bisection.
//
// - SplitMarkedEdges: splits a set of marked edges at their midpoint; all the
//   tetras around a split edge are split in two, so the mesh stays conforming
//   by construction (one edge at a time, longest first). Can be reused with
//   any marking.
// - MarkTetraEdges: marks all the edges of a set of tetras.
// - MarkEdgesBetweenSingularities: marks the regions where two different
//   singular lines of a frame field are closer than a given number of edges
//   (measured with a topological geodesic on the vertex graph), so that after
//   the split the sheets have enough room to pass between them.

#ifndef TETRA_REFINE
#define TETRA_REFINE

#include "geodesic.h"
#include "point3.h"
#include "temp_mark.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace TetraRefine {

typedef std::pair<int, int> EdgeKey;

inline EdgeKey Key(const int V0, const int V1) {
  return EdgeKey(std::min(V0, V1), std::max(V0, V1));
}

// local edges of a tetra
static const int TetEdges[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                   {1, 2}, {1, 3}, {2, 3}};

// marks all the edges of the given tetras
inline void MarkTetraEdges(const std::vector<std::vector<int>> &Tetra,
                           const std::vector<int> &TetIdx,
                           std::set<EdgeKey> &Marked) {
  for (size_t i = 0; i < TetIdx.size(); i++) {
    int t = TetIdx[i];
    assert(t >= 0);
    assert(t < (int)Tetra.size());
    assert(Tetra[t].size() == 4);
    for (int e = 0; e < 6; e++)
      Marked.insert(Key(Tetra[t][TetEdges[e][0]], Tetra[t][TetEdges[e][1]]));
  }
}

// Splits the marked edges (vertex pairs) at their midpoint. Each split edge
// splits all the tetras around it in two, so the result is always conforming.
// The edges are processed from the longest. TrackedEdges (optional, e.g.
// feature edges as vertex pairs) are updated: a split edge is replaced by its
// two halves. TetOrigin (optional): for each tetra the index of the tetra of
// the input mesh it comes from (to copy per-tetra data, e.g. the frames).
// SplitT (optional): position of the split along the edge, as parameter t in
// (0,1) from the first vertex of the key (the smaller index); the edges not
// in the map are split at the midpoint.
// Returns the number of split edges.
template <class ScalarType>
int SplitMarkedEdges(std::vector<Geo::Point3<ScalarType>> &VertPos,
                     std::vector<std::vector<int>> &Tetra,
                     const std::set<EdgeKey> &Marked,
                     std::vector<EdgeKey> *TrackedEdges = nullptr,
                     std::vector<int> *TetOrigin = nullptr,
                     const std::map<EdgeKey, ScalarType> *SplitT = nullptr) {
  if (TetOrigin != nullptr) {
    TetOrigin->resize(Tetra.size());
    for (size_t t = 0; t < Tetra.size(); t++)
      (*TetOrigin)[t] = (int)t;
  }
  // edge -> tetras
  std::map<EdgeKey, std::vector<int>> EdgeTet;
  for (size_t t = 0; t < Tetra.size(); t++) {
    assert(Tetra[t].size() == 4);
    for (int e = 0; e < 6; e++)
      EdgeTet[Key(Tetra[t][TetEdges[e][0]], Tetra[t][TetEdges[e][1]])]
          .push_back((int)t);
  }

  // longest first (more regular tetras)
  std::vector<std::pair<ScalarType, EdgeKey>> Order;
  for (std::set<EdgeKey>::const_iterator it = Marked.begin();
       it != Marked.end(); it++) {
    if (EdgeTet.count(*it) == 0)
      continue; // not an edge of the mesh
    ScalarType Len = (VertPos[it->first] - VertPos[it->second]).Norm();
    Order.push_back(std::make_pair(Len, *it));
  }
  std::sort(Order.begin(), Order.end());
  std::reverse(Order.begin(), Order.end());

  std::set<EdgeKey> Tracked;
  if (TrackedEdges != nullptr)
    for (size_t i = 0; i < TrackedEdges->size(); i++)
      Tracked.insert(Key((*TrackedEdges)[i].first, (*TrackedEdges)[i].second));

  auto RemoveTet = [&EdgeTet](const EdgeKey &E, const int t) {
    std::vector<int> &L = EdgeTet[E];
    L.erase(std::remove(L.begin(), L.end(), t), L.end());
  };

  int NumSplit = 0;
  for (size_t i = 0; i < Order.size(); i++) {
    EdgeKey E = Order[i].second;
    std::map<EdgeKey, std::vector<int>>::iterator ItE = EdgeTet.find(E);
    if (ItE == EdgeTet.end())
      continue;
    std::vector<int> Around = ItE->second;
    int a = E.first;
    int b = E.second;

    // new vertex at the midpoint (or at the given parameter)
    int m = (int)VertPos.size();
    // VertPos.push_back((VertPos[a] + VertPos[b]) / (ScalarType)2);
    ScalarType tE = (ScalarType)0.5;
    if (SplitT != nullptr) {
      typename std::map<EdgeKey, ScalarType>::const_iterator ItT =
          SplitT->find(E);
      if (ItT != SplitT->end())
        tE = ItT->second;
    }
    VertPos.push_back(VertPos[a] + (VertPos[b] - VertPos[a]) * tE);

    for (size_t k = 0; k < Around.size(); k++) {
      int t = Around[k];
      std::vector<int> T = Tetra[t];
      int c = -1, d = -1;
      for (int j = 0; j < 4; j++)
        if ((T[j] != a) && (T[j] != b)) {
          if (c == -1)
            c = T[j];
          else
            d = T[j];
        }
      assert((c != -1) && (d != -1));

      // T1 = T with b -> m (kept in place), T2 = T with a -> m (new);
      // replacing a vertex with a point of the segment keeps the orientation
      std::vector<int> T1 = T, T2 = T;
      for (int j = 0; j < 4; j++) {
        if (T1[j] == b)
          T1[j] = m;
        if (T2[j] == a)
          T2[j] = m;
      }
      int t2 = (int)Tetra.size();
      Tetra[t] = T1;
      Tetra.push_back(T2);
      if (TetOrigin != nullptr)
        TetOrigin->push_back((*TetOrigin)[t]);

      // update the edge map
      RemoveTet(Key(b, c), t);
      RemoveTet(Key(b, d), t);
      EdgeTet[Key(b, c)].push_back(t2);
      EdgeTet[Key(b, d)].push_back(t2);
      EdgeTet[Key(c, d)].push_back(t2);
      EdgeTet[Key(a, m)].push_back(t);
      EdgeTet[Key(m, c)].push_back(t);
      EdgeTet[Key(m, d)].push_back(t);
      EdgeTet[Key(m, b)].push_back(t2);
      EdgeTet[Key(m, c)].push_back(t2);
      EdgeTet[Key(m, d)].push_back(t2);
    }
    EdgeTet.erase(E);

    if (Tracked.count(E) > 0) {
      Tracked.erase(E);
      Tracked.insert(Key(a, m));
      Tracked.insert(Key(m, b));
    }
    NumSplit++;
  }

  if (TrackedEdges != nullptr)
    TrackedEdges->assign(Tracked.begin(), Tracked.end());
  return NumSplit;
}

// Marks the regions where two different singular lines are too close.
// SingEdges: the singular edges as vertex pairs (any grouping, e.g. the
// output of GetPolylineVertFromSing). A multi-source geodesic (number of
// edges) is run from all the singular vertices: on a mesh edge (u,v) whose
// two vertices are reached from different singular vertices s0, s1 the gap
// between them is Steps[u] + Steps[v] + 1. If it is below MinSteps and s0, s1
// are not close along the singular lines (so they are really two different
// lines, not two points of the same line or of two lines meeting at a node),
// the region is marked:
// - OnlyPathEdges = false: all the edges of the tetras around the vertices of
//   the two geodesic paths (more regular, but the mesh grows a lot when the
//   lines are close everywhere)
// - OnlyPathEdges = true: only the edges of the paths between the two lines
//   (each split adds a vertex across the corridor: much cheaper, may need
//   more rounds since the other short paths are found at the next round)
// Returns the number of narrow contacts found.
// Marks the edges of the tetras around the singular lines: the tetras with a
// vertex within Rings vertex rings of a singular vertex (Rings = 0: the
// tetras touching the singular lines). After the split the sheets have room
// to pass beside a line (around it, as in the dual of a hex mesh) instead of
// touching it. Returns the number of tetras whose edges are marked.
inline int MarkEdgesAroundSingularities(
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<std::vector<std::pair<int, int>>> &SingEdges, int Rings,
    std::set<EdgeKey> &Marked) {
  int NumV = 0;
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int v : Tetra[t])
      NumV = std::max(NumV, v + 1);
  std::vector<char> Near(NumV, 0);
  for (size_t i = 0; i < SingEdges.size(); i++)
    for (size_t j = 0; j < SingEdges[i].size(); j++) {
      Near[SingEdges[i][j].first] = 1;
      Near[SingEdges[i][j].second] = 1;
    }
  // grow by rings: the vertices of the tetras touching the current set
  for (int r = 0; r < Rings; r++) {
    std::vector<char> NewNear = Near;
    for (size_t t = 0; t < Tetra.size(); t++) {
      bool Touch = false;
      for (int v : Tetra[t])
        if (Near[v])
          Touch = true;
      if (Touch)
        for (int v : Tetra[t])
          NewNear[v] = 1;
    }
    Near = NewNear;
  }
  std::vector<int> TetIdx;
  for (size_t t = 0; t < Tetra.size(); t++) {
    bool Touch = false;
    for (int v : Tetra[t])
      if (Near[v])
        Touch = true;
    if (Touch)
      TetIdx.push_back((int)t);
  }
  MarkTetraEdges(Tetra, TetIdx, Marked);
  return (int)TetIdx.size();
}

template <class ScalarType>
int MarkEdgesBetweenSingularities(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<std::vector<std::pair<int, int>>> &SingEdges,
    const int MinSteps, std::set<EdgeKey> &Marked,
    const bool OnlyPathEdges = false) {

  // vertex graph of the tet mesh
  std::vector<std::set<int>> VVSet(VertPos.size());
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int e = 0; e < 6; e++) {
      int V0 = Tetra[t][TetEdges[e][0]];
      int V1 = Tetra[t][TetEdges[e][1]];
      VVSet[V0].insert(V1);
      VVSet[V1].insert(V0);
    }
  std::vector<std::vector<int>> VV(VertPos.size());
  for (size_t i = 0; i < VVSet.size(); i++)
    VV[i].assign(VVSet[i].begin(), VVSet[i].end());

  // singular graph
  std::map<int, std::vector<int>> SingAdj;
  for (size_t i = 0; i < SingEdges.size(); i++)
    for (size_t j = 0; j < SingEdges[i].size(); j++) {
      int V0 = SingEdges[i][j].first;
      int V1 = SingEdges[i][j].second;
      SingAdj[V0].push_back(V1);
      SingAdj[V1].push_back(V0);
    }
  if (SingAdj.empty())
    return 0;

  // multi-source geodesic (number of edges) from the singular vertices
  typename Geo::Geodesic<ScalarType>::GeodesicParam GeoP;
  for (std::map<int, std::vector<int>>::iterator it = SingAdj.begin();
       it != SingAdj.end(); it++)
    GeoP.Seeds.push_back(it->first);
  GeoP.unitary_dist = true;
  Geo::TMark TempMark;
  Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, VV, TempMark, GeoP);

  // true if S1 is within MaxDist edges from S0 along the singular lines
  const int MaxAlongDist = 3 * MinSteps;
  auto CloseAlongLines = [&SingAdj, MaxAlongDist](int S0, int S1) {
    std::map<int, int> Dist;
    std::queue<int> Q;
    Dist[S0] = 0;
    Q.push(S0);
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      if (u == S1)
        return true;
      if (Dist[u] >= MaxAlongDist)
        continue;
      const std::vector<int> &Adj = SingAdj[u];
      for (size_t i = 0; i < Adj.size(); i++)
        if (Dist.count(Adj[i]) == 0) {
          Dist[Adj[i]] = Dist[u] + 1;
          Q.push(Adj[i]);
        }
    }
    return false;
  };

  std::set<int> NarrowVert;
  std::set<std::pair<int, int>> CheckedPairs;
  int NumNarrow = 0;
  for (size_t u = 0; u < VV.size(); u++)
    for (size_t k = 0; k < VV[u].size(); k++) {
      int v = VV[u][k];
      if (v <= (int)u)
        continue;
      int S0 = GeoP.Source[u];
      int S1 = GeoP.Source[v];
      if ((S0 < 0) || (S1 < 0) || (S0 == S1))
        continue;
      int Gap = GeoP.Steps[u] + GeoP.Steps[v] + 1;
      if (Gap >= MinSteps)
        continue;
      std::pair<int, int> SKey(std::min(S0, S1), std::max(S0, S1));
      bool SameLine;
      if (CheckedPairs.count(SKey) > 0)
        continue;
      CheckedPairs.insert(SKey);
      SameLine = CloseAlongLines(S0, S1);
      if (SameLine)
        continue;
      NumNarrow++;
      // the two geodesic paths down to the singular vertices
      for (int x = (int)u; x != -1; x = GeoP.Father[x])
        NarrowVert.insert(x);
      for (int x = v; x != -1; x = GeoP.Father[x])
        NarrowVert.insert(x);
      if (OnlyPathEdges) {
        Marked.insert(Key((int)u, v));
        for (int x = (int)u; GeoP.Father[x] != -1; x = GeoP.Father[x])
          Marked.insert(Key(x, GeoP.Father[x]));
        for (int x = v; GeoP.Father[x] != -1; x = GeoP.Father[x])
          Marked.insert(Key(x, GeoP.Father[x]));
      }
    }
    if (OnlyPathEdges)
      return NumNarrow;

  // all the edges of the tetras around the narrow vertices
  std::vector<int> NarrowTet;
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int j = 0; j < 4; j++)
      if (NarrowVert.count(Tetra[t][j]) > 0) {
        NarrowTet.push_back((int)t);
        break;
      }
  MarkTetraEdges(Tetra, NarrowTet, Marked);
  return NumNarrow;
}

} // namespace TetraRefine

#endif
