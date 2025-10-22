#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_ILP_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_ILP_M3

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <include/gurobi_c++.h>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>
#include "../connectivity.h"
#include "../mesh_functions.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "smooth_mesh.h"

#define SOFT_SEED_CONST_M3 100000 // 10000 // 10000
#define POW_EXP 4

template <class ScalarType> class SheetExtractorILPM3 {
public:
  // different types of energy minimization
  enum ObjType { FaceOBJ, ConnOBJ, FaceConnOBJ };

private:
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;

  std::vector<std::vector<int>> &InternalFaces;

  // for each node the variable on the Gurobi solve
  std::vector<int> M3NodeToVar;

  // for each variable on the Gurobi solve to the node
  std::vector<int> M3VarToNode;

  // for each connection Var the list of nodes
  std::vector<std::pair<int, int>> M3LinkVarToNodes;

  // its weigth
  std::vector<ScalarType> M3LinkVarWeight;

  // for each pair of nodes the variable
  std::map<std::pair<int, int>, int> M3NodesToLinkVar;

  // only one single istance of the face can be chosen (one face might
  // correspont to multiple nodes)
  std::vector<std::vector<int>> SingleF_ConstrM3;

  // the edge can be transversed only twice maximum (SET CORRECT VALENCE with
  // Borders/Convex/Concave)
  std::vector<std::vector<int>> TwiceEdge_ConstrM3;

  // exit M3 should be only one for edge
  std::vector<std::vector<int>> M3_Exit_Constr;

  // Continuing M3 Label
  std::vector<std::vector<int>> M3_Label_Cont_Constr;

  //  // Prohibited Neighbors
  // std::vector<std::pair<int,int>> M3_Not_Legal_Connections;

  bool IsVariableM3(int IndexN) const {
    // return false;
    if (!M3Graph.IsEnabled(IndexN))
      return false;
    if (M3Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    return true;
  }

  void ConstraintStats(const std::vector<std::vector<int>> &Constraint,
                       int &MaxSize, ScalarType &AVGSize) {
    AVGSize = 0;
    MaxSize = 0;
    for (size_t i = 0; i < Constraint.size(); i++) {
      MaxSize = std::max(MaxSize, (int)Constraint[i].size());
      AVGSize += (ScalarType)Constraint[i].size();
    }
    AVGSize /= (ScalarType)Constraint.size();
  }

  void ConstraintStats(
      const std::vector<std::pair<int, std::vector<int>>> &Constraint,
      int &MaxSize, ScalarType &AVGSize) {
    AVGSize = 0;
    MaxSize = 0;
    for (size_t i = 0; i < Constraint.size(); i++) {
      MaxSize = std::max(MaxSize, (int)Constraint[i].second.size());
      AVGSize += (ScalarType)Constraint[i].second.size();
    }
    AVGSize /= (ScalarType)Constraint.size();
  }

  void InitVariableMap() {

    M3NodeToVar = std::vector<int>(M3Graph.NumNodes(), -1);
    M3VarToNode.clear();

    int CurrV = 0;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!IsVariableM3(IndexN))
        continue;

      M3NodeToVar[IndexN] = CurrV;
      M3VarToNode.push_back(IndexN);
      CurrV++;
    }

    M3LinkVarToNodes.clear();
    M3NodesToLinkVar.clear();
    M3LinkVarWeight.clear();

    // for each node
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!IsVariableM3(IndexN))
        continue;
      // get all active neighbours
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      M3Graph.ActiveNeigh(IndexN, IndexNeigh, WeightNeigh);
      for (size_t j = 0; j < IndexNeigh.size(); j++) {
        {
          int CurrNeigh = IndexNeigh[j];
          assert(CurrNeigh != IndexN);
          assert(CurrNeigh >= 0);
          assert(CurrNeigh < M3Graph.NumNodes());
          if (!IsVariableM3(CurrNeigh))
            continue;
          // then get the variables
          int IndexVar0 = M3NodeToVar[IndexN];
          int IndexVar1 = M3NodeToVar[CurrNeigh];
          assert(IndexVar0 >= 0);
          assert(IndexVar0 < M3VarToNode.size());
          assert(IndexVar1 >= 0);
          assert(IndexVar1 < M3VarToNode.size());
          std::pair<int, int> KeyLink(std::min(IndexVar0, IndexVar1),
                                      std::max(IndexVar0, IndexVar1));
          if (M3NodesToLinkVar.count(KeyLink) > 0)
            continue;
          // add new variable
          M3NodesToLinkVar[KeyLink] = M3LinkVarToNodes.size();
          M3LinkVarToNodes.push_back(KeyLink);
          M3LinkVarWeight.push_back(WeightNeigh[j]);
        }
      }
    }
  }

  // void InitNotLegalConnections_ConstrM3() {
  //   M3_Not_Legal_Connections.clear();
  //   for (size_t i=0;i<M3Graph.NumNodes();i++)
  //   {
  //     int IndexN0=i;
  //     if (!IsVariableM3(IndexN0))continue;
  //     //get active neighbours
  //     std::vector<int> ActiveNeigh;
  //     std::vector<ScalarType> WeightNeigh;
  //     M3Graph.ActiveNeigh(IndexN0, ActiveNeigh, WeightNeigh);
  //     std::set<int> ActiveNeighSet(ActiveNeigh.begin(),ActiveNeigh.end());
  //     assert(IndexN0<=M3Graph.InternalTangentConstr.size());
  //     std::vector<int> AllNeigh=M3Graph.InternalTangentConstr[IndexN0];
  //     // std::set<int> AllNeighSet(AllNeigh.begin(),AllNeigh.end());
      
  //     //  for (size_t j = 0; j < AllNeigh.size(); j++) {
  //     //    int IndexN1=AllNeigh[j];
  //     //    if (!IsVariableM3(IndexN1))continue;
  //     //    if (ActiveNeighSet.count(IndexN1))continue;
  //     //    std::pair<int,int> AvoidConn(std::min(IndexN0,IndexN1),std::max(IndexN0,IndexN1));
  //     //    M3_Not_Legal_Connections.push_back(AvoidConn);
  //     //  }
  //   }

  //   std::sort(M3_Not_Legal_Connections.begin(), M3_Not_Legal_Connections.end());
  //   M3_Not_Legal_Connections.erase(std::unique(M3_Not_Legal_Connections.begin(), 
  //                                             M3_Not_Legal_Connections.end()), 
  //                                             M3_Not_Legal_Connections.end());

  // }

  void InitSingleF_ConstrM3() {
    SingleF_ConstrM3.clear();

    for (size_t i = 0; i < InternalFaces.size(); i++) {
      std::vector<int> NodesIdx;
      Field::FrameFieldM3Graph::NodesOfElem(i, NodesIdx);
      std::vector<int> VarNodesForF;
      for (size_t j = 0; j < NodesIdx.size(); j++) {
        int testN = NodesIdx[j];
        if (!IsVariableM3(testN))
          continue;
        VarNodesForF.push_back(testN);
      }
      // no need if there is zero or none active
      // variable for such a face
      if (VarNodesForF.size() <= 1)
        continue;

      // otherwise add a new constraint
      SingleF_ConstrM3.resize(SingleF_ConstrM3.size() + 1);
      for (size_t j = 0; j < VarNodesForF.size(); j++) {
        int IndexN = VarNodesForF[j];

        assert(IsVariableM3(IndexN));
        assert(IndexN >= 0);
        assert(IndexN < M3NodeToVar.size());

        int IndexVar = M3NodeToVar[IndexN];
        assert(IndexVar >= 0);
        SingleF_ConstrM3.back().push_back(IndexVar);
      }
    }
  }

  void InitTwiceEdge_Constr() {

    TwiceEdge_ConstrM3.clear();

    std::map<std::pair<int, int>, std::vector<int>> EdgeNodesM3;

    for (size_t i = 0; i < InternalFaces.size(); i++)
      for (size_t j = 0; j < InternalFaces[i].size(); j++) {
        int IndexV0 = InternalFaces[i][j];
        int IndexV1 = InternalFaces[i][(j + 1) % InternalFaces[i].size()];

        std::vector<int> NodesIdxM3;
        Field::FrameFieldM3Graph::NodesOfElem(i, NodesIdxM3);
        assert(NodesIdxM3.size()==3);
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));

        for (size_t k = 0; k < NodesIdxM3.size(); k++) {
          int IndexN = NodesIdxM3[k];
          if (!IsVariableM3(IndexN))
            continue;

          int IndexVar = M3NodeToVar[IndexN];
          assert(IndexVar >= 0);
          EdgeNodesM3[Key].push_back(IndexVar);
        }
      }

    std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
    for (IteMap = EdgeNodesM3.begin(); IteMap != EdgeNodesM3.end(); IteMap++) {
      std::vector<int> NodesOfEdge = (*IteMap).second;
      std::sort(NodesOfEdge.begin(), NodesOfEdge.end());
      NodesOfEdge.erase(std::unique(NodesOfEdge.begin(), NodesOfEdge.end()),
                        NodesOfEdge.end());
      if (NodesOfEdge.size() < 3)
        continue;
      TwiceEdge_ConstrM3.push_back(NodesOfEdge);
    }
  }

  // only one times for each direction for each node
  std::vector<std::pair<int, std::vector<int>>> EdgeLabel_Constr;

  // // for the border edges a different kind of constraint
  // std::vector<std::vector<int>> EdgeLabel_ExitConstr;
  std::vector<std::vector<int>> Internal_TangentNodes_Constr,
      Border_TangentNodes_Constr;


  void InitInternalTangentNodes_Constr() {
    Internal_TangentNodes_Constr.clear();
    for (size_t i = 0; i < M3Graph.InternalTangentConstr.size(); i++) {
      std::vector<int> CurrTangVar;
      for (size_t j = 0; j < M3Graph.InternalTangentConstr[i].size(); j++) {
        int IndexN = M3Graph.InternalTangentConstr[i][j];
        if (!IsVariableM3(IndexN))
          continue;
        int IndexVar = M3NodeToVar[IndexN];
        assert(IndexVar >= 0);
        assert(IndexVar < M3VarToNode.size());
        CurrTangVar.push_back(IndexVar);
      }

      if (CurrTangVar.size() <= 1)
        continue;
      Internal_TangentNodes_Constr.push_back(CurrTangVar);
    }
  }

  void InitBorderTangentNodes_Constr() {
    Border_TangentNodes_Constr.clear();
    for (size_t i = 0; i < M3Graph.BorderTangentConstr.size(); i++) {
      std::vector<int> CurrTangVar;
      for (size_t j = 0; j < M3Graph.BorderTangentConstr[i].size(); j++) {
        int IndexN = M3Graph.BorderTangentConstr[i][j];
        if (!IsVariableM3(IndexN))
          continue;
        int IndexVar = M3NodeToVar[IndexN];
        assert(IndexVar >= 0);
        assert(IndexVar < M3VarToNode.size());
        CurrTangVar.push_back(IndexVar);
      }

      if (CurrTangVar.size() <= 1)
        continue;
      Border_TangentNodes_Constr.push_back(CurrTangVar);
    }
  }

  void InitLingNode_Constr() {}

  void InitSingleLabelF_Constr() {

    EdgeLabel_Constr.clear();

    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {

      int IndexN0 = i;
      if (!IsVariableM3(IndexN0))
        continue;

      std::vector<std::vector<int>> LabelVars;
      int IndexVar0 = M3NodeToVar[IndexN0];
      assert(IndexVar0 >= 0);

      // one label for each edge
      LabelVars.resize(3);

      for (size_t l = 0; l < 3; l++) {

        int IndexLab = l;
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        M3Graph.ActiveNeigh(IndexN0, IndexLab, IndexNeigh, WeightNeigh);

        for (size_t j = 0; j < IndexNeigh.size(); j++) {
          int IndexN1 = IndexNeigh[j];
          assert(IndexLab >= 0);
          assert(IndexLab < 3);
          // if (!IsVariableM3(IndexN1))
          //   continue;
          assert(IsVariableM3(IndexN1));

          int IndexVar1 = M3NodeToVar[IndexN1];
          assert(IndexVar1 >= 0);
          LabelVars[IndexLab].push_back(IndexVar1);
        }
      }

      for (size_t j = 0; j < LabelVars.size(); j++) {
        // if (LabelVars[j].size() < 2)
        //   continue;

        std::set<int> LabelVarsSet(LabelVars[j].begin(), LabelVars[j].end());
        assert(LabelVarsSet.size() == LabelVars[j].size());

        bool IsExitNodeLab = M3Graph.IsExitLabel(IndexN0, j);
        // if no exit separate them so different constraint
        if (!IsExitNodeLab) {
          EdgeLabel_Constr.resize(EdgeLabel_Constr.size() + 1);
          EdgeLabel_Constr.back().first = IndexVar0;
          EdgeLabel_Constr.back().second = LabelVars[j];
        }
      }
    }
  }

  bool PrintSolution() {
    std::cout << "** M3 **" << std::endl;
    for (size_t i = 0; i < M3Solution.size(); i++) {
      if (!M3Solution[i])
        continue;
      std::cout << "Chosen M3 Variable " << i << std::endl;
      std::cout << "Corresponding to M3 Node " << M3VarToNode[i] << std::endl;
    }
  }

  bool SolveGurobi(const std::vector<int> &HardConstrNodes,
                   const std::vector<int> &SoftConstrNodes,
                   const ScalarType &Gap, const ObjType &OType,
                   std::vector<int> &ChosenF, std::vector<int> &M3Nodes) {
    ChosenF.clear();
    M3Nodes.clear();

    ScalarType DeltaMin = 0.00000001;
    try {
      // Create Gurobi environment
      GRBEnv env = GRBEnv(true);
      // env.set("OutputFlag", "0"); // Set to 1 for detailed Gurobi output
      env.start();

      // Create Gurobi model
      GRBModel model = GRBModel(env);
      model.reset();

      std::set<int> HardConstrNodesSet(HardConstrNodes.begin(),
                                       HardConstrNodes.end());
      std::set<int> SoftConstrNodesSet(SoftConstrNodes.begin(),
                                       SoftConstrNodes.end());

      // *** VARIABLE DECLARATION ***
      // add M3 variables
      std::vector<GRBVar> M3NodeVars(M3VarToNode.size());
      for (int i = 0; i < M3VarToNode.size(); ++i) {
        int IndexN = M3VarToNode[i];
        ScalarType ObjW = pow(M3Graph.NodeW[IndexN], POW_EXP);

        // in this case only weight the links
        if (OType == ConnOBJ)
          ObjW = 0;

        if (SoftConstrNodesSet.count(IndexN) > 0)
          ObjW = -(SOFT_SEED_CONST_M3);

        // if (M3Graph.HasExitLabel(IndexN))
        //   ObjW = -100;
        // ObjW = -(SOFT_SEED_CONST_M3*M3Graph.NodeW[IndexN]);

        if (HardConstrNodesSet.count(IndexN) > 0)
        {
          //std::cout<<"Added Hard"<<std::endl;
          M3NodeVars[i] = model.addVar(1.0, 1.0, ObjW, GRB_BINARY,
                                       "M3_N" + std::to_string(i));
        }
        else
          M3NodeVars[i] = model.addVar(0.0, 1.0, ObjW, GRB_BINARY,
                                       "M3_N" + std::to_string(i));
      }

      // add the link nodes if needed
      std::vector<GRBVar> M3LinkVars;
      if ((OType == ConnOBJ) || (OType == FaceConnOBJ)) {
        M3LinkVars.resize(M3LinkVarToNodes.size());
        for (int i = 0; i < M3LinkVars.size(); ++i) {
          ScalarType LinkW = pow(M3LinkVarWeight[i], POW_EXP);
          M3LinkVars[i] = model.addVar(0.0, 1.0, LinkW, GRB_BINARY,
                                       "Link_W" + std::to_string(i));
        }
      }

      // // add M3 variables
      // std::vector<GRBVar> M3NodeVars(M3VarToNode.size());

      // // *** OBJECTIVE FUNCTION ***
      // // Set the objective: minimize the arc value
      // GRBLinExpr objective = 0;
      // for (int i = 0; i < M3NodeVars.size(); ++i) {
      //   int CurrN = M3VarToNode[i];
      //   assert(CurrN >= 0);
      //   assert(CurrN < M3Graph.NodeW.size());
      //   assert(M3Graph.IsEnabled(CurrN));
      //   ScalarType CurrW = M3Graph.NodeW[CurrN];
      //   assert(CurrW > 0);
      //   objective += (M3NodeVars[i] * (DeltaMin + CurrW));
      // }

      // // ADD SOFT INPUT CONSTRAINTS
      // // std::set<int> SoftSeed;
      // for (int i = 0; i < SoftConstrNodes.size(); ++i) {

      //   // add the ortho direction
      //   int IndexNM3 = SoftConstrNodes[i];
      //   if (IsVariableM3(IndexNM3)) {
      //     int IndexVarM3 = M3NodeToVar[IndexNM3];
      //     // SoftSeed.insert(IndexVarM3);
      //     assert(IndexVarM3 >= 0);
      //     assert(IndexVarM3 < M3NodeVars.size());
      //     objective += (-SOFT_SEED_CONST * M3NodeVars[IndexVarM3]);
      //   }
      // }

      // model.setObjective(objective, GRB_MINIMIZE);

      // // *** HARD INPUT CONSTRAINTS ***
      // // ADD HARD INPUT CONSTRAINTS
      // for (int i = 0; i < HardConstrNodes.size(); ++i) {

      //   // add the ortho direction
      //   int IndexNM3 = HardConstrNodes[i];
      //   if (IsVariableM3(IndexNM3)) {
      //     int IndexVarM3 = M3NodeToVar[IndexNM3];
      //     assert(IndexVarM3 >= 0);
      //     assert(IndexVarM3 < M3NodeVars.size());
      //     model.addConstr(M3NodeVars[IndexVarM3] == 1,
      //                     "Hard Input" + std::to_string(i));
      //   }
      // }

      // *** CONNECTION CONSTRAINTS ***
      if ((OType == ConnOBJ) || (OType == FaceConnOBJ)) {
        for (size_t i = 0; i < M3LinkVarToNodes.size(); i++) {
          int IndexVarM3Node0 = M3LinkVarToNodes[i].first;
          int IndexVarM3Node1 = M3LinkVarToNodes[i].second;
          int IndexVarLink = i;
          assert(IndexVarLink <= M3LinkVars.size());
          assert(IndexVarM3Node0 <= M3NodeVars.size());
          assert(IndexVarM3Node1 <= M3NodeVars.size());
          GRBLinExpr lin_constr;
          model.addConstr(M3LinkVars[IndexVarLink] <=
                              M3NodeVars[IndexVarM3Node0],
                          "And 0 " + std::to_string(i));
          model.addConstr(M3LinkVars[IndexVarLink] <=
                              M3NodeVars[IndexVarM3Node1],
                          "And 1 " + std::to_string(i));
          model.addConstr(M3LinkVars[IndexVarLink] >=
                              M3NodeVars[IndexVarM3Node0] +
                                  M3NodeVars[IndexVarM3Node1] - 1,
                          "And 0 " + std::to_string(i));
        }
      }
      // *** MESHING CONSTRAINTS ***

      // ADD SINGLE NODE FOR M3 FACE CONSTRAINT
      for (size_t i = 0; i < SingleF_ConstrM3.size(); i++) {
        GRBLinExpr lin_constr;
        for (size_t j = 0; j < SingleF_ConstrM3[i].size(); j++) {
          int IndexV = SingleF_ConstrM3[i][j];
          lin_constr += M3NodeVars[IndexV];
        }
        model.addConstr(lin_constr <= 1, "Single F M3" + std::to_string(i));
      }

      // ADD TWICE EDGE CONSTR M3
      for (size_t i = 0; i < TwiceEdge_ConstrM3.size(); i++) {
        GRBLinExpr twice_constr;
        for (size_t j = 0; j < TwiceEdge_ConstrM3[i].size(); j++) {
          int IndexVar = TwiceEdge_ConstrM3[i][j];
          assert(IndexVar >= 0);
          assert(IndexVar < M3NodeVars.size());
          twice_constr += M3NodeVars[IndexVar];
        }
        model.addConstr(twice_constr <= 4, "Twice M3 " + std::to_string(i));
      }

      // ADD EDGE LABEL CONSTRAINT
      for (size_t i = 0; i < EdgeLabel_Constr.size(); i++) {
        GRBLinExpr upper_bound;
        int IndexVar0 = EdgeLabel_Constr[i].first;
        assert(IndexVar0 < M3NodeVars.size());
        assert(IndexVar0 >= 0);
        //assert(EdgeLabel_Constr[i].second.size() > 0);

        for (size_t j = 0; j < EdgeLabel_Constr[i].second.size(); j++) {
          int IndexVar1 = EdgeLabel_Constr[i].second[j];
          assert(IndexVar0 != IndexVar1);
          assert(IndexVar1 < M3NodeVars.size());
          assert(IndexVar1 >= 0);
          upper_bound += M3NodeVars[IndexVar1];
        }

        //  model.addGenConstrIndicator(M3NodeVars[IndexVar0], 1, upper_bound ==
        //  1,
        //                             "UpperBound" + std::to_string(i));

        model.addConstr(M3NodeVars[IndexVar0] <= upper_bound,
                        "LowerBound" + std::to_string(i));
        // model.addConstr(upper_bound<=1,
        //                 "LowerBound" + std::to_string(i));
      }

      // ADD internal Label CONSTR
      for (size_t i = 0; i < Internal_TangentNodes_Constr.size(); i++) {

        // if (EdgeLabel_ExitConstr[i].size() <= 1)
        if (Internal_TangentNodes_Constr[i].size() <= 1)
          continue;
        GRBLinExpr upper_bound;

        for (size_t j = 0; j < Internal_TangentNodes_Constr[i].size(); j++) {
          int IndexVar = Internal_TangentNodes_Constr[i][j];

          assert(IndexVar < M3NodeVars.size());
          assert(IndexVar >= 0);
          upper_bound += M3NodeVars[IndexVar];
        }
        model.addConstr(upper_bound <= 2,
                        "Internal Tangent C " + std::to_string(i));
      }

      // ADD Exit Label CONSTR
      for (size_t i = 0; i < Border_TangentNodes_Constr.size(); i++) {

        // if (EdgeLabel_ExitConstr[i].size() <= 1)
        if (Border_TangentNodes_Constr[i].size() <= 1)
          continue;
        GRBLinExpr upper_bound;

        for (size_t j = 0; j < Border_TangentNodes_Constr[i].size(); j++) {
          int IndexVar = Border_TangentNodes_Constr[i][j];

          assert(IndexVar < M3NodeVars.size());
          assert(IndexVar >= 0);
          upper_bound += M3NodeVars[IndexVar];
        }
        model.addConstr(upper_bound <= 1,
                        "Border Tangent C " + std::to_string(i));
      }

      // exit(0);

      // Record the first timestamp
      // model.presolve();
      auto start = std::chrono::high_resolution_clock::now();

      // model.set(GRB_DoubleParam_Heuristics, 0.5);
      // model.set(GRB_DoubleParam_MIPGap, 0.05); // 1% gap
      model.set(GRB_DoubleParam_MIPGap, Gap); // 1% gap
      model.set(GRB_IntParam_NumericFocus, 1);
      model.set(GRB_IntParam_Method, 1);
      model.set(GRB_IntParam_MIPFocus, 1);
      // model.set(GRB_IntParam_NormAdjust, 2);

      // model.set(GRB_IntParam_Presolve, 1);  // Conservative presolve
      // model.tune();
      //    model.set(GRB_DoubleParam_MIPGapAbs, 10); // Absolute gap
      //     // Set time limit to 300 seconds (5 minutes)
      //  //model.set(GRB_DoubleParam_TimeLimit, 3.0);
      //    model.set(GRB_IntParam_MIPFocus, 1);
      //    model.set(GRB_IntParam_Threads, 8); // Use 4 threads
      //    model.set(GRB_IntParam_Cuts, 2); // Moderate cut aggressiveness
      //    model.set(GRB_IntParam_NodeMethod, 2);
      model.optimize();
      auto end = std::chrono::high_resolution_clock::now();
      // Calculate the duration
      auto duration =
          std::chrono::duration_cast<std::chrono::seconds>(end - start);
      // Output the duration
      std::cout << "*** TOTAL TIME ELAPSED : " << duration.count()
                << " *** seconds\n";

      // COLLECT THE RESULT
      if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        std::cout << "Optimal solution found:\n";

        // Read back the solution for the variables
        // std::vector<int> solution(n);
        M3Solution.clear();
        for (int i = 0; i < M3NodeVars.size(); ++i) {
          double ValueVar = M3NodeVars[i].get(GRB_DoubleAttr_X);

          if (ValueVar > 0.5) {
            M3Solution.push_back(true);
            int IndexN = M3VarToNode[i];
            M3Nodes.push_back(IndexN);
            int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
            ChosenF.push_back(IndexF);
          } else {
            M3Solution.push_back(false);
          }
        }

        std::cout << "Chosen " << ChosenF.size() << " faces" << std::endl;

        // PrintSolution();
        // CheckSolution();
        return true;
      } else {
        std::cout << "No optimal solution found.\n";
        return false;
      }
    }

    catch (GRBException &e) {
      std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (const std::exception &e) {
      std::cerr << "Standard exception: " << e.what() << "\n";
    }
    return false;
  }

  void ExtractSubMeshes(
      const std::vector<Geo::Point3<ScalarType>> &OriginalVertPos,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      std::vector<std::vector<std::vector<int>>> &FoilFaces,
      std::vector<std::vector<int>> &OriginalFaceIdx,
      std::vector<std::vector<int>> &OriginalVertIdx) {
    FoilFaces.clear();
    FoilVertices.clear();
    OriginalFaceIdx.clear();
    OriginalVertIdx.clear();

    // get the selected nodes
    std::set<int> SelectedNodes;
    for (size_t i = 0; i < M3Solution.size(); i++) {
      if (!M3Solution[i])
        continue;
      int IndexN = M3VarToNode[i];
      SelectedNodes.insert(IndexN);
    }

    // split into parts
    std::set<int> VisitedNodes;
    std::vector<std::vector<int>> Groups;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      // if not selected
      if (SelectedNodes.count(i) == 0)
        continue;

      // if already visited
      if (VisitedNodes.count(i) > 0)
        continue;

      // then create a new group
      Groups.resize(Groups.size() + 1);

      // then get connected components
      std::vector<int> ExplStack;
      ExplStack.push_back(i);
      do {
        int CurrN = ExplStack.back();
        ExplStack.pop_back();
        int CurrF = Field::FrameFieldM3Graph::ElemOfNode(CurrN);
        // push current node
        Groups.back().push_back(CurrF);
        VisitedNodes.insert(CurrN);

        // explore neighbours
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        // Field::FrameFieldGraph::OrthoNodesFrameField(int IndexN, int &OrthoN)
        M3Graph.ActiveNeigh(CurrN, IndexNeigh, WeightNeigh);
        for (size_t j = 0; j < IndexNeigh.size(); j++) {

          // // if not enabled link
          // if (!Graph.IsEnabledL(CurrN, j))
          //   continue;
          int NeighNode = IndexNeigh[j]; // Graph.NeighNode(CurrN, j);

          // if not active
          if (SelectedNodes.count(NeighNode) == 0)
            continue;

          // if already explored
          if (VisitedNodes.count(NeighNode) > 0)
            continue;

          // otherwise push in the stack
          ExplStack.push_back(NeighNode);
        }

      } while (!ExplStack.empty());
    }

    // transform groups into faces
    for (size_t i = 0; i < Groups.size(); i++) {
      FoilFaces.resize(FoilFaces.size() + 1);
      FoilVertices.resize(FoilVertices.size() + 1);
      OriginalFaceIdx.resize(OriginalFaceIdx.size() + 1);
      OriginalVertIdx.resize(OriginalVertIdx.size() + 1);

      std::sort(Groups[i].begin(), Groups[i].end());
      Groups[i].erase(std::unique(Groups[i].begin(), Groups[i].end()),
                      Groups[i].end());

      CopySubMesh<ScalarType>(OriginalVertPos, InternalFaces, Groups[i],
                              FoilVertices.back(), FoilFaces.back(),
                              OriginalVertIdx.back(), OriginalFaceIdx.back());
    }
  }

  void SmoothReprojectSubMeshes(
      std::vector<Geo::Point3<ScalarType>> &NonManifV,
      std::vector<std::vector<int>> &NonManifFaces,
      const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
      const std::vector<std::vector<int>> &ReprojectFaces,
      const std::vector<std::pair<int, int>> &FeatureE,
      const std::vector<int> &Corners, int SmoothSteps) {

    // get corner and features position
    std::set<Geo::Point3<ScalarType>> FeatureVPos, CornerVPos;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    for (size_t i = 0; i < FeatureE.size(); i++) {
      int IndexV0 = FeatureE[i].first;
      int IndexV1 = FeatureE[i].second;
      Geo::Point3<ScalarType> Pos0 = ReprojectV[IndexV0];
      Geo::Point3<ScalarType> Pos1 = ReprojectV[IndexV1];
      FeatureVPos.insert(Pos0);
      FeatureVPos.insert(Pos1);
      Geo::Segment3<ScalarType> Seg(Pos0, Pos1);
      FeatureSegs.push_back(Seg);
    }

    for (size_t i = 0; i < Corners.size(); i++) {
      int IndexV = Corners[i];
      CornerVPos.insert(ReprojectV[IndexV]);
    }

    // then get the indexes on tne non manif
    std::set<int> FeatureIdxSet, CornerIdxSet;
    for (size_t i = 0; i < NonManifV.size(); i++) {
      Geo::Point3<ScalarType> TestP = NonManifV[i];
      if (FeatureVPos.count(TestP) > 0)
        FeatureIdxSet.insert(i);
      if (CornerVPos.count(TestP) > 0)
        CornerIdxSet.insert(i);
    }

    std::vector<int> FeatureIdxVert(FeatureIdxSet.begin(), FeatureIdxSet.end());
    std::vector<int> CornerIdxVert(CornerIdxSet.begin(), CornerIdxSet.end());

    std::vector<Geo::Point3<ScalarType>> PrevPos = NonManifV;

    // create the reprojection data structure
    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    TriGrid.Init(ReprojectV, ReprojectFaces);

    // get the ones to be reprojected
    std::vector<size_t> BorderV;
    Geo::getBorderVertNonManifMesh(NonManifFaces, BorderV);

    for (size_t s = 0; s < SmoothSteps; s++) {

      Geo::SmoothLaplacianNonManif<ScalarType>(NonManifV, NonManifFaces, 1);

      // reproject regular
      for (size_t i = 0; i < BorderV.size(); i++) {
        int IndexV = BorderV[i];

        if (FeatureIdxSet.count(IndexV))
          continue;
        if (CornerIdxSet.count(IndexV))
          continue;

        assert(IndexV >= 0);
        assert(IndexV < NonManifV.size());
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        int ClosestF;
        Geo::Point3<ScalarType> ClosestPos;
        TriGrid.GridClosest(ReprojectV, ReprojectFaces, Pos,
                            TriGrid.bbox.Diag(), ClosestF, ClosestPos);
        NonManifV[IndexV] = ClosestPos;
      }

      // reproject features
      for (size_t i = 0; i < FeatureIdxVert.size(); i++) {
        int IndexV = FeatureIdxVert[i];
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        ScalarType MinD = std::numeric_limits<ScalarType>::max();
        int BestSeg = -1;
        Geo::Point3<ScalarType> BestClos = Pos;
        for (size_t j = 0; j < FeatureSegs.size(); j++) {
          Geo::Segment3<ScalarType> Seg = FeatureSegs[j];
          Geo::Point3<ScalarType> TestClos = Geo::ClosestPoint(Seg, Pos);
          ScalarType TestD = (TestClos - Pos).Norm();
          if (TestD > MinD)
            continue;
          MinD = TestD;
          BestSeg = j;
          BestClos = TestClos;
        }
        if (BestSeg == -1)
          continue;
        NonManifV[IndexV] = BestClos;
      }

      // put back corners on original positions
      for (size_t i = 0; i < CornerIdxVert.size(); i++) {
        int IndexV = CornerIdxVert[i];
        NonManifV[IndexV] = PrevPos[IndexV];
      }
    }
  }

public:
  struct ExtractParam {
    std::vector<int> HardNodes;
    std::vector<int> SoftNodes;
    bool SingleFConstrM3;
    bool TwiceEdgeConstrM3;
    bool InternalTangentM3;
    bool BorderTangentM3;
    bool LabelM3;
    //bool GlobalNeigh;
    ScalarType Gap;
    ObjType OType;
    // bool SingleLabelCostr;
    // bool IncompatibleNPairConstr;

    ExtractParam() {
      SingleFConstrM3 = true;
      TwiceEdgeConstrM3 = true;
      LabelM3 = true;
      InternalTangentM3 = true;
      BorderTangentM3 = true;
      //GlobalNeigh=true;
      Gap = 0.1;
      //OType = FaceOBJ;
      // OType = ConnOBJ;
      OType =FaceConnOBJ;
    }
  };

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

    ExtractSubMeshes(OriginalVertPos, FoilVertices, FoilFaces, OriginalFaceIdx,
                     OriginalVertIdx);

    if (SmoothSteps == 0)
      return;

    // concatenate to a single surface
    std::vector<Geo::Point3<ScalarType>> NonManifV;
    std::vector<std::vector<int>> NonManifFaces;
    AppendMeshes<ScalarType>(NonManifV, NonManifFaces, FoilVertices, FoilFaces);
    RemoveDuplicatedVert(NonManifV, NonManifFaces);

    std::map<Geo::Point3<ScalarType>, int> PosIdxMap;
    for (size_t i = 0; i < NonManifV.size(); i++)
      PosIdxMap[NonManifV[i]] = i;

    SmoothReprojectSubMeshes(NonManifV, NonManifFaces, ReprojectV,
                             ReprojectFaces, FeatureE, Corners, SmoothSteps);

    // copy back positions
    for (size_t i = 0; i < FoilVertices.size(); i++)
      for (size_t j = 0; j < FoilVertices[i].size(); j++) {
        Geo::Point3<ScalarType> VertP = FoilVertices[i][j];
        assert(PosIdxMap.count(VertP) > 0);
        int IndexNonManif = PosIdxMap[VertP];
        assert(IndexNonManif >= 0);
        assert(IndexNonManif < NonManifV.size());
        FoilVertices[i][j] = NonManifV[IndexNonManif];
      }
  }

  bool FindSurfaceSheet(const ExtractParam &ExtrP, 
                        std::vector<int> &ChosenF,
                        std::vector<bool> &MatchedSoftSeed,
                        bool cut_space = false) {
    std::vector<int> DisabledN;
    std::vector<int> Seeds = ExtrP.SoftNodes;
    Seeds.insert(Seeds.end(), ExtrP.HardNodes.begin(), ExtrP.HardNodes.end());
    MatchedSoftSeed.clear();
    MatchedSoftSeed=std::vector<bool>(ExtrP.SoftNodes.size(),false);
    // std::set<int> ActiveNodes;
    // Geo::LabelGraphFuctions<ScalarType>::GetActiveFrom(M3Graph,Seeds,ActiveNodes);
    if (cut_space) {
      Geo::LabelGraphFuctions<ScalarType>::DisableNotUsefulFromSeed(
          M3Graph, Seeds, DisabledN);
      std::cout << "Disabled From Seeds:" << DisabledN.size() << std::endl;
    }

    InitVariableMap();

    std::cout << "Num M3 Nodes Variables: " << M3VarToNode.size() << std::endl;

    SingleF_ConstrM3.clear();
    if (ExtrP.SingleFConstrM3) {
      InitSingleF_ConstrM3();
      std::cout << "Initialized " << SingleF_ConstrM3.size()
                << " single face constraint" << std::endl;
      int MaxSize;
      ScalarType AVGSize;
      ConstraintStats(SingleF_ConstrM3, MaxSize, AVGSize);
      std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize << std::endl;
    }

    // if (ExtrP.GlobalNeigh)
    // {
    //   InitNotLegalConnections_ConstrM3();
    //   std::cout << "Initialized " << M3_Not_Legal_Connections.size()
    //             << " global connection constraint" << std::endl;
    // }

    std::vector<std::vector<int>> SingleF_ConstrM6;
    TwiceEdge_ConstrM3.clear();
    if (ExtrP.TwiceEdgeConstrM3) {
      InitTwiceEdge_Constr();
      std::cout << "Initialized " << TwiceEdge_ConstrM3.size()
                << " M3 twice edge constraint" << std::endl;
      int MaxSize;
      ScalarType AVGSize;
      ConstraintStats(TwiceEdge_ConstrM3, MaxSize, AVGSize);
      std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize << std::endl;
    }

    EdgeLabel_Constr.clear();
    Internal_TangentNodes_Constr.clear();
    Border_TangentNodes_Constr.clear();
    // EdgeLabel_ExitConstr.clear();
    InitSingleLabelF_Constr();

    if (ExtrP.InternalTangentM3) {
      InitInternalTangentNodes_Constr();
      std::cout << "Initialized " << EdgeLabel_Constr.size()
                << " edge label constraint" << std::endl;
      int MaxSize;
      ScalarType AVGSize;
      ConstraintStats(EdgeLabel_Constr, MaxSize, AVGSize);
      std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize << std::endl;
    }

    if (ExtrP.BorderTangentM3) {
      InitBorderTangentNodes_Constr();

      std::cout << "Initialized " << Border_TangentNodes_Constr.size()
                << " border tangent constraint" << std::endl;
      int MaxSize;
      ScalarType AVGSize;
      ConstraintStats(Border_TangentNodes_Constr, MaxSize, AVGSize);
      std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize << std::endl;
    }

    // exit(0);
    std::vector<int> M3NodesSel;
    bool solved = SolveGurobi(ExtrP.HardNodes, ExtrP.SoftNodes, ExtrP.Gap,
                              ExtrP.OType, ChosenF, M3NodesSel);

    std::set<int> M3NodesSelSet(M3NodesSel.begin(),M3NodesSel.end());
    int NumMatched=0;
    for (size_t i=0;i<ExtrP.SoftNodes.size();i++)
    {
      int TestNode=ExtrP.SoftNodes[i];
      if (M3NodesSelSet.count(TestNode)==0)continue;
      MatchedSoftSeed[i]=true;
      NumMatched++;
    }
    std::cout<<"Matched "<<100*(ScalarType)NumMatched/(ScalarType)ExtrP.SoftNodes.size()<<" Perc of Soft Seeds"<<std::endl;
    if (cut_space) {
      for (size_t i = 0; i < DisabledN.size(); i++)
        M3Graph.EnableNode(DisabledN[i]);
    }

    return solved;
  }

  SheetExtractorILPM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                      std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif