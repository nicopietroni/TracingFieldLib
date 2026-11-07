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

// Orchestration of the sheet solvers:
// - the last SAT solution (and its seeds) is stored as start of the ILP and
//   of the cut refinement
// - StartFromSAT: sets as current solution the stored SAT solution of the
//   current seeds, or runs the SAT
// - SolveILP: ILP from the SAT solution, restricted to a band around it,
//   optionally repeated moving the band on the new solution until it does
//   not change

#ifndef M3_SHEET_TRACING_DRIVER
#define M3_SHEET_TRACING_DRIVER

#include "graph_frame_constraints_M3.h"
#include "graph_frame_field_sheet_solve_ILP_M3.h"
#include "graph_frame_field_sheet_solve_SAT_M3.h"
#include "label_graph.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>

template <class ScalarType> class M3SheetTracingDriver {
public:
  struct Param {
    // the ILP (and the cut) start from the SAT solution of the seeds
    bool StartFromSAT = true;
    // band around the starting solution, steps of the graph (-1 = none)
    int BandSteps = 3;
    // repeat the ILP moving the band until the solution does not change
    bool IterateBand = true;
    int BandMaxIter = 10;
  };
  Param P;

private:
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  SheetSolve_ConstraintsM3<ScalarType> &Constr;
  SheetSolveSATM3<ScalarType> &SAT;
  typename SheetSolveSATM3<ScalarType>::SAT_ExtractParam &SATParam;
  SheetSolveILPM3<ScalarType> &ILP;
  typename SheetSolveILPM3<ScalarType>::ExtractParam &ILPParam;

  // last SAT solution and its seeds (sorted)
  std::vector<bool> SATStartSolution;
  std::vector<int> SATStartSeeds;

  std::vector<int> AllILPSeeds() const {
    std::vector<int> Seeds = ILPParam.HardNodes;
    Seeds.insert(Seeds.end(), ILPParam.SoftNodes.begin(),
                 ILPParam.SoftNodes.end());
    return Seeds;
  }

public:
  // to be called when the graph changes (the stored solution may use nodes
  // not enabled any more)
  void ClearSATStart() {
    SATStartSolution.clear();
    SATStartSeeds.clear();
  }

  // to be called after a SAT solve (stores its solution if solved)
  void StoreSAT(bool Solved) {
    if (Solved && Constr.BeenSolved) {
      SATStartSolution = Constr.M3Solution;
      SATStartSeeds = SATParam.HardNodes;
      std::sort(SATStartSeeds.begin(), SATStartSeeds.end());
    } else
      ClearSATStart();
  }

  // the SAT solution of the current seeds becomes the current solution:
  // the stored one if computed for these seeds, otherwise the SAT is run
  void StartFromSAT() {
    std::vector<int> SortedSeeds = SATParam.HardNodes;
    std::sort(SortedSeeds.begin(), SortedSeeds.end());
    if ((SATStartSolution.size() == M3Graph.NumNodes()) &&
        (SATStartSeeds == SortedSeeds)) {
      Constr.M3Solution = SATStartSolution;
      Constr.BeenSolved = true;
      std::cout << "** Start: stored SAT solution for these seeds **"
                << std::endl;
      return;
    }
    std::cout << "** Start: running the SAT **" << std::endl;
    auto t0 = std::chrono::steady_clock::now();
    Constr.InitConstraints(SATParam.HardNodes, true);
    std::vector<bool> SATMatched;
    bool ok = SAT.Solve(SATParam, SATMatched);
    StoreSAT(ok);
    std::cout << "** Start: SAT "
              << (ok ? "solved" : "not solved (start from the seeds)")
              << " in "
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - t0)
                     .count()
              << " s **" << std::endl;
  }

  // ILP: start (SAT), constraints, single solve or band iterations
  bool SolveILP(std::vector<bool> &MatchedSeed) {
    ILPParam.FeasiblityOnly = false;
    if (P.StartFromSAT)
      StartFromSAT();
    std::vector<int> Seeds = AllILPSeeds();
    // the band (and its iterations) only around the SAT start: without it
    // the whole graph
    const int Band = P.StartFromSAT ? P.BandSteps : -1;
    // if (!(P.IterateBand && (P.BandSteps >= 0))) {
    if (!(P.StartFromSAT && P.IterateBand && (Band >= 0))) {
      Constr.InitConstraints(Seeds, false, Band, true, &ILPParam.HardNodes);
      return ILP.Solve(ILPParam, MatchedSeed);
    }
    bool Any = false;
    for (int it = 0; it < P.BandMaxIter; it++) {
      std::vector<bool> Prev = Constr.M3Solution;
      bool PrevValid = Constr.BeenSolved;
      // only the hard seeds must be in the starting solution
      Constr.InitConstraints(Seeds, false, Band, true, &ILPParam.HardNodes);
      bool FromSolution = Constr.HasStartSolution;
      bool solved = ILP.Solve(ILPParam, MatchedSeed);
      if (!solved) {
        std::cout << "** Band iteration " << it << ": ILP not solved, stop **"
                  << std::endl;
        // keep the last valid solution
        if (PrevValid && (Prev.size() == M3Graph.NumNodes())) {
          Constr.M3Solution = Prev;
          Constr.BeenSolved = true;
        }
        break;
      }
      Any = true;
      int NumNodes = 0;
      for (size_t i = 0; i < Constr.M3Solution.size(); i++)
        if (Constr.M3Solution[i])
          NumNodes++;
      std::cout << "** Band iteration " << it << ": " << NumNodes
                << " faces, energy " << Constr.SolutionW(Constr.M3Solution)
                << " **" << std::endl;
      // without a starting solution the band is not used: one solve
      if (!FromSolution)
        break;
      if (Constr.M3Solution == Prev) {
        std::cout << "** Band iteration converged (solution unchanged) **"
                  << std::endl;
        break;
      }
    }
    return Any;
  }

  // ILP feasibility check (same setup of the solve)
  bool CheckFeasibilityILP(std::vector<bool> &MatchedSeed) {
    Constr.InitConstraints(AllILPSeeds(), false,
                           P.StartFromSAT ? P.BandSteps : -1, true,
                           &ILPParam.HardNodes);
    ILPParam.FeasiblityOnly = true;
    bool solved = ILP.Solve(ILPParam, MatchedSeed);
    ILPParam.FeasiblityOnly = false;
    return solved;
  }

  M3SheetTracingDriver(
      Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
      SheetSolve_ConstraintsM3<ScalarType> &_Constr,
      SheetSolveSATM3<ScalarType> &_SAT,
      typename SheetSolveSATM3<ScalarType>::SAT_ExtractParam &_SATParam,
      SheetSolveILPM3<ScalarType> &_ILP,
      typename SheetSolveILPM3<ScalarType>::ExtractParam &_ILPParam)
      : M3Graph(_M3Graph), Constr(_Constr), SAT(_SAT), SATParam(_SATParam),
        ILP(_ILP), ILPParam(_ILPParam) {}
};

#endif
