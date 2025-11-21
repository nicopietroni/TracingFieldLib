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

#ifndef KISSAT_INTERFACE_H
#define KISSAT_INTERFACE_H

// #pragma once
// #include <vector>
// #include <cassert>
// #include <iostream>

// extern "C" {
// #include "kissat.h"
// }

// class KissatSatWrapperEmbedded
// {
// public:
//     struct Result {
//         bool sat;
//         std::vector<int> assignment;
//     };

//     static Result solve(
//         int numVars,
//         const std::vector<int> &seedTrue,
//         const std::vector<int> &seedFalse,
//         const std::vector<std::vector<int>> &atMostOne,
//         const std::vector<std::vector<int>> &zeroOrTwo)
//     {
//         kissat *solver = kissat_new();

//         // Add clause by adding literals, ending with 0
//         auto addClause = [&](const std::vector<int> &cl) {
//             for (int lit : cl)
//                 kissat_add(solver, lit);   // API corretta
//             kissat_add(solver, 0);          // end of clause
//         };

//         // 1. Seed TRUE
//         for (int v : seedTrue) {
//             addClause({ v });
//         }

//         // 2. Seed FALSE
//         for (int v : seedFalse) {
//             addClause({ -v });
//         }

//         // 3. At-most-one (pairwise)
//         for (const auto &g : atMostOne) {
//             for (size_t i=0;i<g.size();i++)
//                 for (size_t j=i+1;j<g.size();j++)
//                     addClause({ -g[i], -g[j] });
//         }

//         // 4. Zero-or-two (only pairs for now)
//         for (const auto &g : zeroOrTwo) {
//             if (g.size() == 2) {
//                 int a = g[0], b = g[1];
//                 // (¬a ∨  b)
//                 addClause({ -a, b });
//                 // ( a ∨ ¬b)
//                 addClause({  a, -b });
//             } else {
//                 std::cerr << "[WARNING] zeroOrTwo group > 2 not
//                 implemented.\n";
//             }
//         }

//         // Solve
//         int res = kissat_solve(solver);

//         Result R;
//         R.assignment.assign(numVars+1, 0);

//         if (res == 20) { // UNSAT
//             R.sat = false;
//             kissat_release(solver);
//             return R;
//         }

//         // SAT
//         R.sat = true;
//         for (int v=1; v<=numVars; v++) {
//             int val = kissat_value(solver, v);
//             R.assignment[v] = (val > 0 ? 1 : 0);
//         }

//         kissat_release(solver);
//         return R;
//     }
// };

#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/wait.h> // for WEXITSTATUS
#include <unordered_map>
#include <unordered_set>
#include <vector>

class KissatExternal {
public:
  struct Result {
    bool sat;
    std::vector<int> assignment; // 1-based, 0/1 (real IDs if remapped)
  };

  // Creazione nuova variabile SAT sicura
  static int newVar() {
    static int counter = 200000000;
    return counter++;
  }

  // AMO usando ladder encoding (molto veloce)
  static void addAtMostOne_Ladder(const std::vector<int> &g,
                                  std::vector<std::vector<int>> &C) {
    int n = g.size();
    if (n <= 1)
      return;

    // crea variabili s1 ... s(n-1)
    std::vector<int> s(n - 1);
    for (int i = 0; i < n - 1; i++)
      s[i] = newVar();

    // x1 → s1     = (¬x1 ∨ s1)
    C.push_back({-g[0], s[0]});

    // regole per i = 1..n-2
    for (int i = 1; i < n - 1; i++) {
      // xi+1 → si     = (¬g[i] ∨ s[i-1])
      C.push_back({-g[i], s[i]});

      // si → si+1
      C.push_back({-s[i - 1], s[i]});

      // si → ¬xi+1
      C.push_back({-s[i - 1], -g[i]});
    }

    // Ultimo vincolo: s(n-1) → ¬xn
    C.push_back({-s[n - 2], -g[n - 1]});
  }

  // Sequential counter AtMost2: sum(g) <= 2
  static void addAtMostTwo_Sequential(const std::vector<int> &g,
                                      std::vector<std::vector<int>> &C) {
    int n = g.size();
    if (n <= 2)
      return; // massimo 2 veri sempre

    // aux[i][1] e aux[i][2] per i = 0..n-1
    std::vector<int> s1(n), s2(n);

    for (int i = 0; i < n; i++) {
      s1[i] = newVar();
      s2[i] = newVar();
    }

    // i = 0
    // s1[0] = x0
    C.push_back({-g[0], s1[0]});
    C.push_back({-s1[0], g[0]});

    // s2[0] = false
    C.push_back({-s2[0]});

    // i >= 1
    for (int i = 1; i < n; i++) {
      // s1[i] = s1[i-1] OR g[i]
      C.push_back({-g[i], s1[i]});
      C.push_back({-s1[i - 1], s1[i]});

      // s2[i] = (s2[i-1]) OR (s1[i-1] AND g[i])
      // Parte AND
      int h = newVar();
      C.push_back({-s1[i - 1], -g[i], h});
      C.push_back({s1[i - 1], -h});
      C.push_back({g[i], -h});

      // s2[i] = h OR s2[i-1]
      C.push_back({-h, s2[i]});
      C.push_back({-s2[i - 1], s2[i]});

      // AtMost2: NON permettere terzo true:
      // se s2[i-1] già indica 2 true → g[i] deve essere false
      C.push_back({-s2[i - 1], -g[i]});
    }
  }

  static void addForbidExactlyOne(const std::vector<int> &g,
                                  std::vector<std::vector<int>> &C) {
    int n = g.size();
    if (n <= 1)
      return;

    for (int i = 0; i < n; i++) {
      std::vector<int> clause;
      clause.reserve(n);
      clause.push_back(-g[i]);
      for (int j = 0; j < n; j++)
        if (j != i)
          clause.push_back(g[j]);
      C.push_back(std::move(clause));
    }
  }

  static void addZeroOrTwo_Sequential(const std::vector<int> &g,
                                      std::vector<std::vector<int>> &C) {
    int n = g.size();
    if (n <= 1)
      return;

    // SUM != 1
    addForbidExactlyOne(g, C);

    // SUM <= 2
    addAtMostTwo_Sequential(g, C);
  }

  static void addAtMostOne_Simple(const std::vector<int> &g,
                                  std::vector<std::vector<int>> &C) {
    const int n = (int)g.size();
    if (n <= 1)
      return;

    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        C.push_back({-g[i], -g[j]});
  }

  static void addZeroOrTwo_Simple(const std::vector<int> &g,
                                  std::vector<std::vector<int>> &C) {
    int n = g.size();
    if (n <= 1)
      return;

    // 1) Forbid exactly one selected
    for (int i = 0; i < n; i++) {
      std::vector<int> clause;
      clause.reserve(n);
      clause.push_back(-g[i]); // ¬xi
      for (int j = 0; j < n; j++)
        if (j != i)
          clause.push_back(g[j]);
      C.push_back(std::move(clause));
    }

    // 2) Forbid >= 3 selected (AtMost2)
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++)
        for (int h = j + 1; h < n; h++)
          C.push_back({-g[i], -g[j], -g[h]});
  }

  // Per ogni pair (x0, group):
  //   se x0 = 1  → uno e uno solo tra group è TRUE
  //   se x0 = 0  → al massimo uno TRUE in group
  //
  // CNF:
  // (1) x0 -> (OR group)
  // (2) AtMostOne(group)

  static void addConditionalOneOf(
      const std::vector<std::pair<int, std::vector<int>>> &groups,
      std::vector<std::vector<int>> &C) {
    for (const auto &p : groups) {
      int x0 = p.first;
      const auto &g = p.second;

      int n = g.size();
      if (n == 0)
        continue;

      // ----------------------------------------------------
      // (1) x0 -> ( x1 OR x2 OR ... OR xn )
      //      (¬x0 ∨ x1 ∨ x2 ∨ ... ∨ xn)
      // ----------------------------------------------------
      {
        std::vector<int> clause;
        clause.reserve(n + 1);
        clause.push_back(-x0);
        for (int v : g)
          clause.push_back(v);
        C.push_back(std::move(clause));
      }

      // ----------------------------------------------------
      // (2) AtMostOne(group)
      //      ∀ i<j: (¬xi ∨ ¬xj)
      // ----------------------------------------------------
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
          C.push_back({-g[i], -g[j]});
    }
  }

  static void addAtMostTwo(const std::vector<int> &g,
                           std::vector<std::vector<int>> &C) {
    const int n = (int)g.size();
    if (n <= 2)
      return; // se <=2 sono già ammissibili

    // forbid any triple (xi AND xj AND xh)
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        for (int h = j + 1; h < n; ++h)
          C.push_back({-g[i], -g[j], -g[h]});
  }

  static Result solve(const std::vector<int> &seedTrue,
                      const std::vector<int> &seedFalse,
                      const std::vector<std::vector<int>> &atMostOne,
                      const std::vector<std::vector<int>> &atMostTwo,
                      const std::vector<std::vector<int>> &zeroOrTwo,
                      const std::vector<std::pair<int, std::vector<int>>> &ImplyOne,
                      bool useRemap = true,
                      const std::string &kissatPath =
                          "/Users/134750/Desktop/00_GeoProjects/"
                          "TracingFieldLib/external_lib/kissat/build/kissat") {
    //----------------------------------------------------------------------
    // TEMP FILES
    //----------------------------------------------------------------------
    const std::string cnfFile = "/tmp/tmp_kissat.cnf";
    const std::string outFile = "/tmp/tmp_kissat.out";

    //----------------------------------------------------------------------
    // (1) BUILD CLAUSES IN REAL-LITERAL SPACE
    //----------------------------------------------------------------------
    std::vector<std::vector<int>> clausesReal;

    // Seeds TRUE
    for (int v : seedTrue)
      clausesReal.push_back({v});

    // Seeds FALSE
    for (int v : seedFalse)
      clausesReal.push_back({-v});

    // At-Most-One constraints
    for (auto &g : atMostOne)
      addAtMostOne_Simple(g, clausesReal);

    // At-Most-Two constraints
     for (auto &g : atMostTwo)
      addAtMostTwo(g, clausesReal);

    // Zero-or-Two constraints
    for (auto &g : zeroOrTwo)
        addZeroOrTwo_Simple(g, clausesReal);                        
      //addZeroOrTwo_Sequential(g, clausesReal);

    // Imply-One constraints
    addConditionalOneOf(ImplyOne, clausesReal); 

    //----------------------------------------------------------------------
    // (2) REMAP VARIABLE IDS — FULLY DETERMINISTIC
    //----------------------------------------------------------------------
    std::vector<int> back;   // satVar → realVar
    std::vector<int> allIds; // sorted list of real IDs

    int satVars = 0;

    if (useRemap) {
      // Collect all real variable IDs
      for (auto &cl : clausesReal)
        for (int lit : cl)
          allIds.push_back(std::abs(lit));

      // Sort + unique = fully deterministic
      std::sort(allIds.begin(), allIds.end());
      allIds.erase(std::unique(allIds.begin(), allIds.end()), allIds.end());

      // Allocate remap
      satVars = (int)allIds.size();
      back.assign(satVars + 1, 0);

      // Deterministic real→sat map
      // No unordered_map needed — use binary search on sorted array
      auto mapRealToSat = [&](int real) {
        auto it = std::lower_bound(allIds.begin(), allIds.end(), real);
        int idx = int(it - allIds.begin()) + 1; // 1-based SAT variable
        return idx;
      };

      //------------------------------------------------------------------
      // (3) BUILD SAT CLAUSES WITH THE NEW MAPPING
      //------------------------------------------------------------------
      std::vector<std::vector<int>> clausesSat;
      clausesSat.reserve(clausesReal.size());

      for (auto &cl : clausesReal) {
        std::vector<int> satC;
        satC.reserve(cl.size());

        for (int lit : cl) {
          int real = std::abs(lit);
          int satv = mapRealToSat(real);
          int satLit = (lit > 0) ? satv : -satv;
          satC.push_back(satLit);

          back[satv] = real;
        }

        clausesSat.push_back(satC);
      }

      clausesReal.clear(); // free memory, unused now

      //------------------------------------------------------------------
      // (4) WRITE CNF — deterministically ordered
      //------------------------------------------------------------------
      std::ofstream cnf(cnfFile);
      if (!cnf) {
        std::cerr << "Error: cannot write CNF file\n";
        return {false, {}};
      }

      cnf << "p cnf " << satVars << " " << clausesSat.size() << "\n";

      for (auto &cl : clausesSat) {
        for (int lit : cl)
          cnf << lit << " ";
        cnf << "0\n";
      }
      cnf.close();
    } else {
      //------------------------------------------------------------------
      // (NO REMAP): find max var
      //------------------------------------------------------------------
      for (auto &cl : clausesReal)
        for (int lit : cl)
          satVars = std::max(satVars, std::abs(lit));

      std::ofstream cnf(cnfFile);
      if (!cnf) {
        std::cerr << "Error: cannot write CNF file\n";
        return {false, {}};
      }

      cnf << "p cnf " << satVars << " " << clausesReal.size() << "\n";
      for (auto &cl : clausesReal) {
        for (int lit : cl)
          cnf << lit << " ";
        cnf << "0\n";
      }
      cnf.close();
    }

    //----------------------------------------------------------------------
    // (5) RUN KISSAT (DETERMINISTIC IF CNF IS DETERMINISTIC!)
    //----------------------------------------------------------------------
    std::string kissatDir = kissatPath.substr(0, kissatPath.find_last_of('/'));
    std::string kissatExe = kissatPath.substr(kissatPath.find_last_of('/') + 1);

    // Fixed seed = deterministic model
    std::string cmd = "cd " + kissatDir + " && ./" + kissatExe +
                      " --seed=1234 --relaxed " + cnfFile;

    std::cout << "Running kissat: " << cmd << std::endl;
    FILE *pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
      std::cerr << "Error: cannot run kissat\n";
      return {false, {}};
    }
    std::cout << "Done" << std::endl;
    char buffer[4096];
    std::string solverOutput;
    while (fgets(buffer, sizeof(buffer), pipe))
      solverOutput += buffer;

    int status = pclose(pipe);
    int exitCode = WEXITSTATUS(status);

    if (exitCode != 10 && exitCode != 20) {
      std::cerr << "Error: kissat failed (exit " << exitCode << ")\n";
      std::cerr << "OUTPUT:\n" << solverOutput << "\n";
      return {false, {}};
    }

    //----------------------------------------------------------------------
    // (6) PARSE SOLUTION
    //----------------------------------------------------------------------
    Result R;
    R.sat = (exitCode == 10);
    R.assignment.assign(satVars + 1, 0);

    std::stringstream ss(solverOutput);
    std::string line;

    while (std::getline(ss, line)) {
      if (line.size() == 0)
        continue;
      if (line[0] != 'v')
        continue;

      std::stringstream ls(line.substr(2));
      int lit;
      while (ls >> lit) {
        if (lit == 0)
          continue;

        int sv = std::abs(lit);
        int val = (lit > 0);

        if (sv > satVars)
          continue;

        if (useRemap) {
          int real = back[sv];
          if (real >= 0 && real < (int)R.assignment.size())
            R.assignment[real] = val;
        } else {
          if (sv < (int)R.assignment.size())
            R.assignment[sv] = val;
        }
      }
    }

    return R;
  }
};

#endif
