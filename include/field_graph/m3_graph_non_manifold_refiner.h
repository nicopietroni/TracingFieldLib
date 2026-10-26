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

// Repair of the non manifold (node, vertex) pairs of the M3 graph by
// refinement instead of adding links.
// A pair (n, v) fails because around v there are no faces aligned enough
// with the sheet of n to close a local disk. The plane through v orthogonal
// to the direction of n is inserted in the star of v: every edge of the link
// of v (edges of the tetras around v opposite to v) crossed by the plane is
// split at the crossing point. In each tetra (v,a,b,c) the plane crosses two
// edges of abc sharing a vertex, and the conforming split always creates the
// face (v, p_ab, p_ac) lying in the plane: a node exactly aligned with the
// sheet. Only the tetras around the split edges change (no cascade).
// Several requests on the same edge: the ones with a similar direction (same
// sheet seen from close vertices) are averaged, the edge takes the largest
// group, the others wait for the next round (the edge is split by then).
// Splits too close to a vertex are skipped (that vertex is already about on
// the plane, and slivers are avoided).
// The graph must be rebuilt after each refinement (see main: rounds of
// refine, graph, check).

#ifndef M3_GRAPH_NON_MANIFOLD_REFINER
#define M3_GRAPH_NON_MANIFOLD_REFINER

#include "../frame_field.h"
#include "../frame_field_refine.h"
#include "../point3.h"
#include "../tetra_refine.h"
#include "graph_frame_connectivity_M3.h"
#include "graph_frame_field.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

template <class ScalarType> class M3GraphNonManifoldRefiner {
public:
  struct Param {
    // splits closer than this (parameter along the edge) to a vertex are
    // skipped
    ScalarType MinT = (ScalarType)0.2;
    // two requests on an edge are the same sheet if the angle between their
    // directions is below this (degrees, unsigned directions)
    ScalarType SameDirAngle = (ScalarType)30;
    // also the pairs with an incorrect exit
    bool IncludeIncorrectExit = true;
    // a tetra around v is not cut if one of its faces through v is already
    // about on the plane (angle between its normal and the direction below
    // this, degrees): the sheet can already use that face, only the gaps of
    // the disk are cut
    ScalarType GoodFaceAngle = (ScalarType)25;
  };
  Param P;

  struct Stats {
    int NumPairs = 0;      // failing pairs used
    int NumRequests = 0;   // edge crossings requested
    int NumEdges = 0;      // edges to split
    int NumConflicts = 0;  // edges with requests of different sheets
    int NumTooClose = 0;   // edges skipped (split too close to a vertex)
    int NumTetCut = 0;     // tetras around the failing vertices to cut
    int NumTetGood = 0;    // tetras skipped (a face already on the plane)
  };
  Stats LastStats;

private:
  struct Request {
    ScalarType T;                // from the first vertex of the key
    Geo::Point3<ScalarType> Dir; // normal of the plane
  };

public:
  // split positions for the failing pairs of the connectivity check
  // (Conn must be checked on the current graph)
  void ComputeSplits(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Tetra,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const LabelGraph_ConnectivityM3<ScalarType> &Conn,
      std::map<TetraRefine::EdgeKey, ScalarType> &SplitT) {
    SplitT.clear();
    LastStats = Stats();

    // vertex -> tetras
    std::vector<std::vector<int>> VT(VertPos.size());
    for (size_t t = 0; t < Tetra.size(); t++)
      for (int v : Tetra[t])
        VT[v].push_back((int)t);

    const ScalarType CosGood = std::cos(P.GoodFaceAngle * (ScalarType)M_PI /
                                        (ScalarType)180);
    // requests per edge
    std::map<TetraRefine::EdgeKey, std::vector<Request>> Req;
    for (size_t n = 0; n < Conn.ManifoldNodevType.size(); n++) {
      if (!Conn.IsCheckable((int)n))
        continue;
      for (int j = 0; j < 3; j++) {
        bool Fail =
            (Conn.ManifoldNodevType[n][j] ==
             LabelGraph_ConnectivityM3<ScalarType>::NON_MANIFOLD) ||
            (P.IncludeIncorrectExit &&
             (Conn.ManifoldNodevType[n][j] ==
              LabelGraph_ConnectivityM3<ScalarType>::INCORRECT_EXIT));
        if (!Fail)
          continue;
        LastStats.NumPairs++;
        int f = Field::FrameFieldM3Graph::ElemOfNode((int)n);
        int d = Field::FrameFieldM3Graph::DirOfNode((int)n);
        int v = Faces[f][j];
        Geo::Point3<ScalarType> Dir = FacesFrameField[f].GetDirection(d);
        ScalarType L = Dir.Norm();
        if (L <= 0)
          continue;
        Dir = Dir / L;
        std::set<TetraRefine::EdgeKey> Done;
        // the plane through v orthogonal to Dir, crossing the link of v:
        // only the tetras crossed by the plane without a face through v
        // already about on it (the gaps of the disk)
        for (int t : VT[v]) {
          std::vector<int> Opp;
          for (int w : Tetra[t])
            if (w != v)
              Opp.push_back(w);
          assert(Opp.size() == 3);
          std::vector<ScalarType> Side(3);
          bool Pos = false, Neg = false;
          for (int k = 0; k < 3; k++) {
            Side[k] = Geo::Dot(VertPos[Opp[k]] - VertPos[v], Dir);
            if (Side[k] > 0)
              Pos = true;
            if (Side[k] < 0)
              Neg = true;
          }
          if (!(Pos && Neg))
            continue;
          bool HasGood = false;
          for (int k = 0; k < 3; k++) {
            Geo::Point3<ScalarType> E0 = VertPos[Opp[k]] - VertPos[v];
            Geo::Point3<ScalarType> E1 = VertPos[Opp[(k + 1) % 3]] - VertPos[v];
            Geo::Point3<ScalarType> N = Geo::Cross(E0, E1);
            ScalarType LN = N.Norm();
            if (LN <= 0)
              continue;
            if (std::fabs(Geo::Dot(N / LN, Dir)) >= CosGood)
              HasGood = true;
          }
          if (HasGood) {
            LastStats.NumTetGood++;
            continue;
          }
          LastStats.NumTetCut++;
          for (int k = 0; k < 3; k++) {
            int a = Opp[k], b = Opp[(k + 1) % 3];
            ScalarType fa = Side[k], fb = Side[(k + 1) % 3];
            // strictly on opposite sides
            if (!(((fa < 0) && (fb > 0)) || ((fa > 0) && (fb < 0))))
              continue;
            TetraRefine::EdgeKey K = TetraRefine::Key(a, b);
            // the same edge from another tetra of the same pair
            if (Done.count(K) > 0)
              continue;
            Done.insert(K);
            if (K.first != a)
              std::swap(fa, fb);
            Request R;
            R.T = fa / (fa - fb);
            R.Dir = Dir;
            Req[K].push_back(R);
            LastStats.NumRequests++;
          }
        }
      }
    }

    // one position per edge: largest group of similar directions, averaged
    const ScalarType CosA = std::cos(P.SameDirAngle * (ScalarType)M_PI /
                                     (ScalarType)180);
    for (typename std::map<TetraRefine::EdgeKey,
                           std::vector<Request>>::iterator it = Req.begin();
         it != Req.end(); it++) {
      const std::vector<Request> &R = it->second;
      int Best = -1, BestSize = 0;
      for (size_t i = 0; i < R.size(); i++) {
        int Size = 0;
        for (size_t k = 0; k < R.size(); k++)
          if (std::fabs(Geo::Dot(R[i].Dir, R[k].Dir)) >= CosA)
            Size++;
        if (Size > BestSize) {
          BestSize = Size;
          Best = (int)i;
        }
      }
      if (BestSize < (int)R.size())
        LastStats.NumConflicts++;
      ScalarType SumT = 0;
      for (size_t k = 0; k < R.size(); k++)
        if (std::fabs(Geo::Dot(R[Best].Dir, R[k].Dir)) >= CosA)
          SumT += R[k].T;
      ScalarType T = SumT / (ScalarType)BestSize;
      if ((T < P.MinT) || (T > (1 - P.MinT))) {
        LastStats.NumTooClose++;
        continue;
      }
      SplitT[it->first] = T;
    }
    LastStats.NumEdges = (int)SplitT.size();
    std::cout << "** Non manifold refine: " << LastStats.NumPairs
              << " failing pairs, " << LastStats.NumRequests
              << " edge crossings, " << LastStats.NumEdges
              << " edges to split (" << LastStats.NumConflicts
              << " with different sheets, " << LastStats.NumTooClose
              << " skipped too close to a vertex), " << LastStats.NumTetCut
              << " tetras cut, " << LastStats.NumTetGood
              << " skipped (a face already on the plane) **" << std::endl;
  }

  // splits the edges at the given positions, the field is inherited.
  // Returns the number of split edges.
  static int Refine(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    std::vector<std::vector<int>> &Tetra,
                    std::vector<std::pair<int, int>> &Features,
                    std::vector<Field::FrameField<ScalarType>> &TetraFrames,
                    const std::map<TetraRefine::EdgeKey, ScalarType> &SplitT) {
    std::set<TetraRefine::EdgeKey> Marked;
    for (typename std::map<TetraRefine::EdgeKey, ScalarType>::const_iterator
             it = SplitT.begin();
         it != SplitT.end(); it++)
      Marked.insert(it->first);
    size_t NumTetBefore = Tetra.size();
    std::vector<int> TetOrigin;
    int NumSplit = TetraRefine::SplitMarkedEdges(VertPos, Tetra, Marked,
                                                 &Features, &TetOrigin, &SplitT);
    FrameFieldRefiner<ScalarType>::InheritFrames(VertPos, Tetra, TetOrigin,
                                                 TetraFrames);
    std::cout << "** Non manifold refine: " << NumSplit << " edges split, "
              << NumTetBefore << " -> " << Tetra.size()
              << " tetras (field inherited) **" << std::endl;
    return NumSplit;
  }
};

#endif
