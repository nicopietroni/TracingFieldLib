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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_SAT_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_SAT_M3

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "../solve/kissat_interface.h"
#include "graph.h"
#include "graph_frame_constraints_M3.h"
#include "graph_frame_field_sheet_extraction_M3.h"
#include "label_graph_functions.h"
#include "smooth_mesh.h"
#include "smooth_sheet_surfaces.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>


template <class ScalarType> class SheetSolveSATM3 {
public:

  struct SAT_ExtractParam {
    std::vector<int> HardNodes;
    typename SheetSolve_ConstraintsM3<ScalarType>::ContinuityType CType;
    SAT_ExtractParam() 
    {
      CType=SheetSolve_ConstraintsM3<ScalarType>::FaceBased;
    }
  };

private:
  //the set of constraints
  const SheetSolve_ConstraintsM3<ScalarType> &M3Constraints;
  //the graph
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  //internal face mesh
  const std::vector<std::vector<int>> &InternalFaces;


  bool SolveSAT(const SAT_ExtractParam &ExtrP, bool perLabel = true) {

    KissatExternal::Result res;

    std::vector<int> SatSeedTrue, SatSeedFalse;
    // map hard nodes to variables
    for (size_t i = 0; i < ExtrP.HardNodes.size(); i++) {
      int IndexN = ExtrP.HardNodes[i];
      if (!M3Constraints.IsVariableM3(IndexN))
        continue;
      int IndexVar = M3Constraints.M3NodeToVar[IndexN];
      SatSeedTrue.push_back(IndexVar);
    }

    // first add the at-most-one constraints
    std::vector<std::vector<int>> AtMostOne = M3Constraints.SingleF_ConstrM3;

    // then add all the pairs of avoided connections as at-most-one
    if (ExtrP.CType==SheetSolve_ConstraintsM3<ScalarType>::EdgeBased) {
      for (size_t i = 0; i < M3Constraints.Avoid_Connections.size(); i++) {
        int VarA = M3Constraints.Avoid_Connections[i].first;
        int VarB = M3Constraints.Avoid_Connections[i].second;
        AtMostOne.push_back({VarA, VarB});
      }
    }

    std::vector<std::vector<int>> AtMostTwo;
    if (ExtrP.CType==SheetSolve_ConstraintsM3<ScalarType>::FaceBased) 
      AtMostTwo= M3Constraints.Internal_TangentNodes_Constr;

    std::vector<std::pair<int, std::vector<int>>> ImplyOne;
    if (ExtrP.CType==SheetSolve_ConstraintsM3<ScalarType>::FaceBased) 
      ImplyOne=M3Constraints.EdgeLabel_Continuation_Constr;
    
    std::vector<std::vector<int>> ZeroOrTwo;
    if (ExtrP.CType==SheetSolve_ConstraintsM3<ScalarType>::EdgeBased)
      ZeroOrTwo = M3Constraints.Internal_TangentNodes_Constr;
    std::cout << "Zero or Two constraints: " << ZeroOrTwo.size() << std::endl;

    res = KissatExternal::solve(SatSeedTrue, SatSeedFalse, AtMostOne, AtMostTwo,
                                ZeroOrTwo, ImplyOne);
    bool solved = (res.sat);
    if (solved)
      std::cout << " ** GF2 M3 solved " << std::endl;
    else
      std::cout << " ** GF2 M3 not solved " << std::endl << std::endl;
    // map back the solution
    if (solved) {
      M3Solution = std::vector<bool>(M3Constraints.M3VarToNode.size(), false);
      for (size_t i = 0; i < M3Constraints.M3VarToNode.size(); i++) {
        int VarValue = res.assignment[i]; // + 1];
        if (VarValue == 1) {
          M3Solution[i] = true;
        } else {
          M3Solution[i] = false;
        }
      }
    } else {
      M3Solution = std::vector<bool>(M3Constraints.M3VarToNode.size(), false);
    }
    return solved;
  }

public:
  std::vector<bool> M3Solution;

  void ExtractSurfacesFormSolution(
      const std::vector<Geo::Point3<ScalarType>> &OriginalVertPos,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      std::vector<std::vector<std::vector<int>>> &FoilFaces,
      std::vector<std::vector<int>> &OriginalFaceIdx,
      std::vector<std::vector<int>> &OriginalVertIdx,
      const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
      const std::vector<std::vector<int>> &ReprojectFaces,
      const std::vector<std::pair<int, int>> &FeatureE,
      const std::vector<int> &Corners, int SmoothSteps = 20) {

    SheetExtractorM3<ScalarType>::ExtractSubMeshes(
        M3Graph, M3Solution, InternalFaces, OriginalVertPos, FoilVertices,
        FoilFaces, OriginalFaceIdx, OriginalVertIdx);

    if (SmoothSteps == 0)
      return;

    FoilsSmooth<ScalarType>::SmoothFoils(FoilVertices, FoilFaces, ReprojectV,
                                         ReprojectFaces, FeatureE, Corners,
                                         SmoothSteps);
  }

  int NumVariables() { return M3Constraints.M3VarToNode.size(); }

  bool Solve(const SAT_ExtractParam &ExtrP, std::vector<bool> &MatchedSeed) {

    auto time0 = std::chrono::high_resolution_clock::now();

    bool solved = SolveSAT(ExtrP);
    if (solved)
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), true);
    else
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);

    std::vector<bool> CompleteM3Solution(M3Graph.NumNodes(), false);
    for (size_t i = 0; i < M3Solution.size(); i++) {
      if (M3Solution[i]) {
        int IndexN = M3Constraints.M3VarToNode[i];
        CompleteM3Solution[IndexN] = true;
      }
    }

    M3Solution = CompleteM3Solution;
    Geo::LabelGraphFuctions<ScalarType>::FilterSelectedFromSeed(
        M3Graph, ExtrP.HardNodes, M3Solution);

    return solved;
  }

  SheetSolveSATM3(const SheetSolve_ConstraintsM3<ScalarType> &_M3Constraints,
                  Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                  std::vector<std::vector<int>> &_InternalFaces)
      : M3Constraints(_M3Constraints),M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif