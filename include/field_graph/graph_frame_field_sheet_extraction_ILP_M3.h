#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_ILP_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_ILP_M3

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "smooth_mesh.h"
#include "smooth_sheet_surfaces.h"
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

#define SOFT_SEED_CONST_M3 100000
#define SOFT_CONSTR_PENALTY 1000
#define POW_EXP 4

template <class ScalarType> class SheetExtractorILPM3 {
public:
  // different types of energy minimization
  enum ObjType { FaceOBJ, ConnOBJ, FaceConnOBJ };
  enum ContinuityType { FaceBased, EdgeBased };

  struct ExtractParam {
    std::vector<int> HardNodes;
    std::vector<int> SoftNodes;

    bool SoftSurfaceContrinuationM3;

    bool SingleFConstrM3;
    bool SoftSingleFConstrM3;

    bool TwiceEdgeConstrM3;
    bool SoftTwiceEdgeConstrM3;

    bool InternalTangentM3;
    bool SoftInternalTangentM3;

    bool BorderTangentM3;
    bool SoftBorderTangentM3;

    bool SoftAvoidConnectionsM3;

    ScalarType Gap;
    ObjType OType;
    bool FeasiblityOnly;
    bool OptimizeObj;

    ContinuityType CType;

    ExtractParam() {
      SingleFConstrM3 = true;
      TwiceEdgeConstrM3 = true;
      // LabelM3 = true;
      InternalTangentM3 = true;
      BorderTangentM3 = true;
      // GlobalNeigh=true;
      Gap = 0.1;
      OType = FaceOBJ;
      FeasiblityOnly = false;
      OptimizeObj = true;

      SoftSurfaceContrinuationM3 = false;
      SoftSingleFConstrM3 = false;
      SoftTwiceEdgeConstrM3 = false;
      SoftInternalTangentM3 = false;
      SoftBorderTangentM3 = false;
      SoftAvoidConnectionsM3 = false;

      CType = FaceBased;
    }
  };

private:
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;

  std::vector<std::vector<int>> &InternalFaces;

  // NODE VARIABLES
  //  for each node the variable on the Gurobi solve
  std::vector<int> M3NodeToVar;

  // for each variable on the Gurobi solve to the node
  std::vector<int> M3VarToNode;

  // LINK VARIABLES IF USED IN THE OBJECTIVE FUNCTION
  //  for each connection Var the list of nodes
  std::vector<std::pair<int, int>> M3LinkVarToNodes;

  // its weigth
  std::vector<ScalarType> M3LinkVarWeight;

  // for each pair of nodes the variable (this is when link variable used)
  std::map<std::pair<int, int>, int> M3NodesToLinkVar;

  // MANIFOLD CONSTRAINTS
  //  only one single istance of the face can be chosen (one face might
  //  correspont to multiple nodes)
  std::vector<std::vector<int>> SingleF_ConstrM3;

  // the edge can be transversed only twice maximum (SET CORRECT VALENCE with
  // Borders/Convex/Concave)
  std::vector<std::vector<int>> TwiceEdge_ConstrM3;

  // nodes which are tangent internally for each edge (valence 0 or 2)
  std::vector<std::vector<int>> Internal_TangentNodes_Constr;

  // nodes which are tangent on border for each edge (valence 0 or 1)
  std::vector<std::vector<int>> Border_TangentNodes_Constr;

  // APPROACH 1, EXPLICIT NEIGHBOURS AND CONSTRAINTS
  // list of pairs for each node the continuation, the first one
  // is repeated three times because the three labels
  std::vector<std::pair<int, std::vector<int>>> EdgeLabel_Continuation_Constr;

  // APPROACH 2, ON GROUP OF EDGES AND PROHIBITED NEIGH
  // the connections that are disabled on a group of faces sharing the edge
  std::vector<std::pair<int, int>> Avoid_Connections;

  // exit M3 should be only one for exit edge
  // std::vector<std::vector<int>> M3_Exit_Constr;

  // // Continuing M3 Label such that one face should trigger the next
  // std::vector<std::vector<int>> M3_Label_Cont_Constr;

  bool IsVariableM3(int IndexN) const {
    // return false;
    if (!M3Graph.IsEnabled(IndexN))
      return false;
    if (M3Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    return true;
  }

  // void ConstraintStats(const std::vector<std::vector<int>> &Constraint,
  //                      int &MaxSize, ScalarType &AVGSize) {
  //   AVGSize = 0;
  //   MaxSize = 0;
  //   for (size_t i = 0; i < Constraint.size(); i++) {
  //     MaxSize = std::max(MaxSize, (int)Constraint[i].size());
  //     AVGSize += (ScalarType)Constraint[i].size();
  //   }
  //   AVGSize /= (ScalarType)Constraint.size();
  // }

  // void ConstraintStats(
  //     const std::vector<std::pair<int, std::vector<int>>> &Constraint,
  //     int &MaxSize, ScalarType &AVGSize) {
  //   AVGSize = 0;
  //   MaxSize = 0;
  //   for (size_t i = 0; i < Constraint.size(); i++) {
  //     MaxSize = std::max(MaxSize, (int)Constraint[i].second.size());
  //     AVGSize += (ScalarType)Constraint[i].second.size();
  //   }
  //   AVGSize /= (ScalarType)Constraint.size();
  // }

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
  //     //    std::pair<int,int>
  //     AvoidConn(std::min(IndexN0,IndexN1),std::max(IndexN0,IndexN1));
  //     //    M3_Not_Legal_Connections.push_back(AvoidConn);
  //     //  }
  //   }

  //   std::sort(M3_Not_Legal_Connections.begin(),
  //   M3_Not_Legal_Connections.end());
  //   M3_Not_Legal_Connections.erase(std::unique(M3_Not_Legal_Connections.begin(),
  //                                             M3_Not_Legal_Connections.end()),
  //                                             M3_Not_Legal_Connections.end());

  // }

  //  for (size_t i = 0; i < M3Graph.InternalTangentConstr.size(); i++) {
  //     std::vector<int> CurrTangVar;
  //     for (size_t j = 0; j < M3Graph.InternalTangentConstr[i].size(); j++) {
  //       int IndexN = M3Graph.InternalTangentConstr[i][j];
  //       if (!IsVariableM3(IndexN))
  //         continue;

  void InitAvoidConnection_ConstrM3() {
    Avoid_Connections.clear();
    for (size_t i = 0; i < M3Graph.InternalTangentConstr.size(); i++) {
      // for each pair
      for (size_t j = 0; j < M3Graph.InternalTangentConstr[i].size() - 1; j++)
        for (size_t k = (j + 1); k < M3Graph.InternalTangentConstr[i].size();
             k++) {
          int IndexN0 = M3Graph.InternalTangentConstr[i][j];
          int IndexN1 = M3Graph.InternalTangentConstr[i][k];
          if (!IsVariableM3(IndexN0))
            continue;
          if (!IsVariableM3(IndexN1))
            continue;
          assert(IndexN0 != IndexN1);
          if (M3Graph.IsActiveNeighbour(IndexN0, IndexN1))
            continue;
          int IndexVar0 = M3NodeToVar[IndexN0];
          int IndexVar1 = M3NodeToVar[IndexN1];
          assert(IndexVar0 >= 0);
          assert(IndexVar1 >= 0);
          Avoid_Connections.push_back(std::pair<int, int>(
              std::min(IndexVar0, IndexVar1), std::max(IndexVar0, IndexVar1)));
        }
    }

    std::sort(Avoid_Connections.begin(), Avoid_Connections.end());
    Avoid_Connections.erase(
        std::unique(Avoid_Connections.begin(), Avoid_Connections.end()),
        Avoid_Connections.end());
  }

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
        assert(NodesIdxM3.size() == 3);
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

  // // for the border edges a different kind of constraint
  // std::vector<std::vector<int>> EdgeLabel_ExitConstr;

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

  void InitSingleLabelF_Constr() {

    EdgeLabel_Continuation_Constr.clear();

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
          EdgeLabel_Continuation_Constr.resize(
              EdgeLabel_Continuation_Constr.size() + 1);
          EdgeLabel_Continuation_Constr.back().first = IndexVar0;
          EdgeLabel_Continuation_Constr.back().second = LabelVars[j];
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

  void CheckBigM(GRBModel &model) {
    for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); i++) {
      GRBConstr c = model.getConstr(i);
      std::string name = c.get(GRB_StringAttr_ConstrName);

      // Ottieni la riga come espressione lineare
      GRBLinExpr expr = model.getRow(c);
      int nTerms = expr.size();

      for (int t = 0; t < nTerms; t++) {
        double coeff = expr.getCoeff(t);
        if (fabs(coeff) > 1e5) { // soglia da regolare
          std::cout << "Possibile Big-M in constraint " << name
                    << " coeff=" << coeff << std::endl;
        }
      }

      double rhs = c.get(GRB_DoubleAttr_RHS);
      if (fabs(rhs) > 1e5)
        std::cout << "Possibile Big-M in constraint " << name << " rhs=" << rhs
                  << std::endl;
    }
  }

  void CoeffDistribution(GRBModel &model, ScalarType &minC, ScalarType &maxC,
                         ScalarType &avgC) {
    minC = 0;
    maxC = 0;
    avgC = 0;

    // get if there are constraints
    int nConstrs = model.get(GRB_IntAttr_NumConstrs);
    if (nConstrs == 0) {
      std::cout << "⚠️  No Constraints." << std::endl;
      return;
    }
    std::vector<ScalarType> coeffs;

    GRBConstr *constrs = model.getConstrs(); // <-- get all constraints

    // for all constraints
    for (int i = 0; i < nConstrs; ++i) {
      GRBLinExpr expr = model.getRow(constrs[i]);
      int nTerms = expr.size();
      for (int t = 0; t < nTerms; ++t) {
        coeffs.push_back(fabs(expr.getCoeff(t)));
      }

      double rhs = constrs[i].get(GRB_DoubleAttr_RHS);
      if (fabs(rhs) > 0)
        coeffs.push_back(fabs(rhs)); // also consider RHS
    }

    delete[] constrs; // deallocate

    if (coeffs.empty()) {
      std::cout << "⚠️  No coefficient found." << std::endl;
      return;
    }

    // Statistiche
    minC = *std::min_element(coeffs.begin(), coeffs.end());
    maxC = *std::max_element(coeffs.begin(), coeffs.end());
    avgC = std::accumulate(coeffs.begin(), coeffs.end(), 0.0) / coeffs.size();
  }

  void ObjDistribution(GRBModel &model, ScalarType &minC, ScalarType &maxC,
                       ScalarType &avgC) {
    minC = 0;
    maxC = 0;
    avgC = 0;
    int nVars = model.get(GRB_IntAttr_NumVars);
    std::vector<double> objCoeffs;
    objCoeffs.reserve(nVars);

    for (int i = 0; i < nVars; ++i) {
      GRBVar v = model.getVar(i);
      double c = fabs(v.get(GRB_DoubleAttr_Obj));
      if (c > 0)
        objCoeffs.push_back(c);
    }

    if (objCoeffs.empty()) {
      std::cout << "⚠️ All coefficient of objective function are NULL!"
                << std::endl;
      return;
    }

    minC = *std::min_element(objCoeffs.begin(), objCoeffs.end());
    maxC = *std::max_element(objCoeffs.begin(), objCoeffs.end());
    avgC = std::accumulate(objCoeffs.begin(), objCoeffs.end(), 0.0) /
           objCoeffs.size();
  }

  void RescaleObjFunctions(GRBModel &model, const double maxConstraintCoeff,
                           const double maxObjCoeff,
                           bool breaks_symmetry = false) {

    // compute scaling factor
    double maxVal = std::max(maxConstraintCoeff, maxObjCoeff);
    if (maxVal < 1e-12) {
      std::cout << "⚠️  maxVal too small skip scaling." << std::endl;
      return;
    }

    double scale = 1000.0 / maxVal;
    std::cout << "🔧 Scaling model by factor " << scale << " (max coeff → 1000)"
              << std::endl;

    // get number of variables
    int nVars = model.get(GRB_IntAttr_NumVars);

    // rescale all objective values
    ScalarType SymmStep = (ScalarType)1 / (ScalarType)nVars;

    for (int i = 0; i < nVars; ++i) {
      GRBVar v = model.getVar(i);
      double c = v.get(GRB_DoubleAttr_Obj);
      if (!breaks_symmetry)
        v.set(GRB_DoubleAttr_Obj, 1 + c * scale);
      else
        v.set(GRB_DoubleAttr_Obj, SymmStep * i + c * scale);
    }

    // // 3️⃣ get all constraints and rescale RHS
    // int nConstrs = model.get(GRB_IntAttr_NumConstrs);
    // if (nConstrs > 0) {
    //   GRBConstr *constrs = model.getConstrs();
    //   for (int i = 0; i < nConstrs; ++i) {
    //     double rhs = constrs[i].get(GRB_DoubleAttr_RHS);
    //     constrs[i].set(GRB_DoubleAttr_RHS, rhs * scale);
    //   }
    //   delete[] constrs;
    // }

    // update model
    model.update();

    std::cout << "✅ Model rescaled successfully." << std::endl;
  }

  void InitNodeVariables(const ExtractParam &ExtrP, GRBModel &model,
                         std::vector<GRBVar> &M3NodeVars) {
    std::set<int> HardConstrNodesSet(ExtrP.HardNodes.begin(),
                                     ExtrP.HardNodes.end());
    std::set<int> SoftConstrNodesSet(ExtrP.SoftNodes.begin(),
                                     ExtrP.SoftNodes.end());

    // *** VARIABLE DECLARATION ***
    // add M3 variables
    M3NodeVars = std::vector<GRBVar>(M3VarToNode.size());
    for (int i = 0; i < M3VarToNode.size(); ++i) {
      int IndexN = M3VarToNode[i];
      ScalarType ObjW = pow(M3Graph.NodeW[IndexN], POW_EXP);

      // in this case only weight the links
      if (ExtrP.OType == ConnOBJ)
        ObjW = 0;

      if (SoftConstrNodesSet.count(IndexN) > 0)
        ObjW = -(SOFT_SEED_CONST_M3);

      // if (M3Graph.HasExitLabel(IndexN))
      //   ObjW = -100;
      // ObjW = -(SOFT_SEED_CONST_M3*M3Graph.NodeW[IndexN]);

      if (HardConstrNodesSet.count(IndexN) > 0) {
        // std::cout<<"Added Hard"<<std::endl;
        M3NodeVars[i] = model.addVar(1.0, 1.0, ObjW, GRB_BINARY,
                                     "M3_N" + std::to_string(i));
      } else
        M3NodeVars[i] = model.addVar(0.0, 1.0, ObjW, GRB_BINARY,
                                     "M3_N" + std::to_string(i));
    }
  }

  void InitLinkVariables(const ExtractParam &ExtrP, GRBModel &model,
                         const std::vector<GRBVar> &M3NodeVars,
                         std::vector<GRBVar> &M3LinkVars) {
    M3LinkVars.resize(M3LinkVarToNodes.size());
    for (int i = 0; i < M3LinkVars.size(); ++i) {
      ScalarType LinkW = pow(M3LinkVarWeight[i], POW_EXP);
      M3LinkVars[i] = model.addVar(0.0, 1.0, LinkW, GRB_BINARY,
                                   "Link_W" + std::to_string(i));
    }

    // *** link the variables to connections in case need for weighting ***
    for (size_t i = 0; i < M3LinkVarToNodes.size(); i++) {
      int IndexVarM3Node0 = M3LinkVarToNodes[i].first;
      int IndexVarM3Node1 = M3LinkVarToNodes[i].second;
      int IndexVarLink = i;
      assert(IndexVarLink <= M3LinkVars.size());
      assert(IndexVarM3Node0 <= M3NodeVars.size());
      assert(IndexVarM3Node1 <= M3NodeVars.size());
      GRBLinExpr lin_constr;
      model.addConstr(M3LinkVars[IndexVarLink] <= M3NodeVars[IndexVarM3Node0],
                      "And 0 " + std::to_string(i));
      model.addConstr(M3LinkVars[IndexVarLink] <= M3NodeVars[IndexVarM3Node1],
                      "And 1 " + std::to_string(i));
      model.addConstr(M3LinkVars[IndexVarLink] >=
                          M3NodeVars[IndexVarM3Node0] +
                              M3NodeVars[IndexVarM3Node1] - 1,
                      "And 0 " + std::to_string(i));
    }
  }

  void AddFoilFaceContinuityConstraint(GRBModel &model,
                                       const std::vector<GRBVar> &M3NodeVars,
                                       bool SoftConstr) {
    // ADD EDGE LABEL CONSTRAINT IF SELECTED A FACE SHOULD SELECT NEIGHBOURS
    for (size_t i = 0; i < EdgeLabel_Continuation_Constr.size(); i++) {
      GRBLinExpr upper_bound;
      int IndexVar0 = EdgeLabel_Continuation_Constr[i].first;
      assert(IndexVar0 < M3NodeVars.size());
      assert(IndexVar0 >= 0);

      for (size_t j = 0; j < EdgeLabel_Continuation_Constr[i].second.size();
           j++) {
        int IndexVar1 = EdgeLabel_Continuation_Constr[i].second[j];
        assert(IndexVar0 != IndexVar1);
        assert(IndexVar1 < M3NodeVars.size());
        assert(IndexVar1 >= 0);
        upper_bound += M3NodeVars[IndexVar1];
      }

      // SOFT_CONSTR_PENALTY
      //  value of a node less then neighbur for a label
      //  means if select it should select one neighbour
      if (!SoftConstr) {
        model.addConstr(M3NodeVars[IndexVar0] <= upper_bound,
                        "LowerBound" + std::to_string(i));
      } else {
        GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY, GRB_BINARY,
                                    "relax_" + std::to_string(i));
        model.addConstr(M3NodeVars[IndexVar0] <= upper_bound + relax);
      }
    }
  }

  void AddSingleFaceConstraint(GRBModel &model,
                               const std::vector<GRBVar> &M3NodeVars,
                               bool SoftConstr) {
    for (size_t i = 0; i < SingleF_ConstrM3.size(); i++) {
      GRBLinExpr lin_constr;
      for (size_t j = 0; j < SingleF_ConstrM3[i].size(); j++) {
        int IndexV = SingleF_ConstrM3[i][j];
        lin_constr += M3NodeVars[IndexV];
      }
      if (!SoftConstr) {
        model.addConstr(lin_constr <= 1, "Single F M3" + std::to_string(i));
      } else {
        // GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
        // GRB_BINARY,
        //                             "Single F M3 relax" + std::to_string(i));
        GRBVar relax =
            model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                         "Single F M3 relax" + std::to_string(i));

        model.addConstr(lin_constr <= 1 + relax,
                        "Single F M3 Soft" + std::to_string(i));
      }
    }
  }

  void AddTwiceEdgeConstraint(GRBModel &model,
                              const std::vector<GRBVar> &M3NodeVars,
                              bool SoftConstr) {
    for (size_t i = 0; i < TwiceEdge_ConstrM3.size(); i++) {
      GRBLinExpr twice_constr;
      for (size_t j = 0; j < TwiceEdge_ConstrM3[i].size(); j++) {
        int IndexVar = TwiceEdge_ConstrM3[i][j];
        assert(IndexVar >= 0);
        assert(IndexVar < M3NodeVars.size());
        twice_constr += M3NodeVars[IndexVar];
      }
      if (!SoftConstr) {
        model.addConstr(twice_constr <= 4, "Twice M3 " + std::to_string(i));
      } else {
        // GRBVar slack = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
        // GRB_BINARY,
        //                             "slack_twice_" + std::to_string(i));
        GRBVar relax =
            model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                         "Twice E M3 relax" + std::to_string(i));
        model.addConstr(twice_constr <= 4 + relax,
                        "Soft_TwiceM3_" + std::to_string(i));
      }
    }
  }

  void AddInternalTangentConstraint(GRBModel &model,
                                    const std::vector<GRBVar> &M3NodeVars,
                                    bool ImposeContinuity, bool SoftConstr) {

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
      // in this case also can be 1 ,continuity imposed elsewhere
      if (!ImposeContinuity) {
        if (!SoftConstr) {
          model.addConstr(upper_bound <= 2,
                          "Internal Tangent C " + std::to_string(i));
        } else {

          GRBVar relax =
              model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                           "Internal relax" + std::to_string(i));
          model.addConstr(upper_bound <= 2 + relax,
                          "Soft_InternalM3_" + std::to_string(i));
        }
      } else // in this case can be 0 or 2
      // because should ensure continuity (then need to avoid disabled
      // connections)
      {
        // add one more variable to ensure that the sum is even
        GRBVar evenVar =
            model.addVar(0.0, 1.0, 0, GRB_BINARY, "relax_" + std::to_string(i));
        if (!SoftConstr) {
          model.addConstr(upper_bound == 2 * evenVar,
                          "Even Continuation " + std::to_string(i));
        } else {

          GRBVar relax =
              model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                           "Even relax" + std::to_string(i));
          model.addConstr(upper_bound == 2 * evenVar + relax,
                          "Soft Even Continuation " + std::to_string(i));
        }
      }
    }
  }

  void AddBorderTangentConstraint(GRBModel &model,
                                  const std::vector<GRBVar> &M3NodeVars,
                                  bool SoftConstr) {
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
      if (!SoftConstr) {
        model.addConstr(upper_bound <= 1,
                        "Border Tangent C " + std::to_string(i));
      } else {
        // GRBVar slack = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
        // GRB_BINARY,
        //                             "slack_border_" + std::to_string(i));
        GRBVar relax =
            model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                         "Border relax" + std::to_string(i));
        model.addConstr(upper_bound <= 1 + relax,
                        "Soft_BorderM3_" + std::to_string(i));
      }
    }
  }

  void AddAvoidConnectionsConstraint(GRBModel &model,
                                     const std::vector<GRBVar> &M3NodeVars,
                                     bool SoftConstr) {

    for (size_t i = 0; i < Avoid_Connections.size(); i++) {
      int IndexVar0 = Avoid_Connections[i].first;
      int IndexVar1 = Avoid_Connections[i].second;

      assert(IndexVar0 >= 0);
      assert(IndexVar1 >= 0);
      assert(IndexVar0 < M3NodeVars.size());
      assert(IndexVar1 < M3NodeVars.size());

      GRBLinExpr no_twice_constr;
      no_twice_constr += M3NodeVars[IndexVar0];
      no_twice_constr += M3NodeVars[IndexVar1];

      if (!SoftConstr) {
        model.addConstr(no_twice_constr <= 1,
                        "Avoid Connection M3 " + std::to_string(i));
      } else {
        GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY, GRB_BINARY,
                                    "relax_" + std::to_string(i));
        model.addConstr(no_twice_constr <= 1 + relax,
                        "Soft Avoid Connection M3 " + std::to_string(i));
      }
    }
  }

  void OptimizeCoefficients(GRBModel &model) {

    model.update();

    CheckBigM(model);
    ScalarType minConstrC = 0;
    ScalarType maxConstrC = 0;
    ScalarType avgConstrC = 0;
    CoeffDistribution(model, minConstrC, maxConstrC, avgConstrC);

    ScalarType minObjC = 0;
    ScalarType maxObjC = 0;
    ScalarType avgObjC = 0;
    ObjDistribution(model, minObjC, maxObjC, avgObjC);

    RescaleObjFunctions(model, maxConstrC, maxObjC);

    CoeffDistribution(model, minConstrC, maxConstrC, avgConstrC);
    std::cout << "📊 Coeff Distribution: "
              << "min=" << minConstrC << "  max=" << maxConstrC
              << "  avg=" << avgConstrC
              << "  ratio=" << maxConstrC / std::max(minConstrC, 1e-9)
              << std::endl;

    ObjDistribution(model, minObjC, maxObjC, avgObjC);

    std::cout << "🎯 Obj coeff: min=" << minObjC << "  max=" << maxObjC
              << "  avg=" << avgObjC
              << "  ratio=" << maxObjC / std::max(minObjC, 1e-9) << std::endl;
  }

  bool SolveGurobi(const ExtractParam &ExtrP, std::vector<int> &ChosenF,
                   std::vector<int> &M3Nodes) {

    ChosenF.clear();
    M3Nodes.clear();

    ScalarType DeltaMin = 0.00000001;
    try {
      // auto time0 = std::chrono::high_resolution_clock::now();

      // Create Gurobi environment
      GRBEnv env = GRBEnv(true);
      // env.set("OutputFlag", "0"); // Set to 1 for detailed Gurobi output
      env.start();

      // Create Gurobi model
      GRBModel model = GRBModel(env);
      model.reset();

      std::vector<GRBVar> M3NodeVars;
      InitNodeVariables(ExtrP, model, M3NodeVars);

      // add the link nodes if needed
      std::vector<GRBVar> M3LinkVars;

      if ((ExtrP.OType == ConnOBJ) || (ExtrP.OType == FaceConnOBJ))
        InitLinkVariables(ExtrP, model, M3NodeVars, M3LinkVars);

      // ADD SINGLE NODE FOR M3 FACE CONSTRAINT
      if (ExtrP.SingleFConstrM3) {
        std::cout << " ** Adding Single F Constraint" << std::endl;
        AddSingleFaceConstraint(model, M3NodeVars, ExtrP.SoftSingleFConstrM3);
      }

      // ADD TWICE EDGE CONSTR M3
      if (ExtrP.TwiceEdgeConstrM3) {
        std::cout << " ** Adding Twice E Constraint" << std::endl;
        AddTwiceEdgeConstraint(model, M3NodeVars, ExtrP.SoftTwiceEdgeConstrM3);
      }

      // ADD internal Label CONSTR
      if (ExtrP.CType == FaceBased) {

        AddFoilFaceContinuityConstraint(model, M3NodeVars,
                                        ExtrP.SoftSurfaceContrinuationM3);

        if (ExtrP.InternalTangentM3) {
          std::cout << " ** Adding Internal Tangent Constraint" << std::endl;
          AddInternalTangentConstraint(model, M3NodeVars, false,
                                       ExtrP.SoftInternalTangentM3);
        }
      } else {
        std::cout << " ** Adding Internal Tangent Constraint with Continuity"
                  << std::endl;
        AddInternalTangentConstraint(model, M3NodeVars, true,
                                     ExtrP.SoftInternalTangentM3);

        std::cout << " ** Adding Avoid Connections" << std::endl;
        AddAvoidConnectionsConstraint(model, M3NodeVars,
                                      ExtrP.SoftAvoidConnectionsM3);
      }
      // ADD Exit Label CONSTR
      if (ExtrP.BorderTangentM3) {
        std::cout << " ** Adding Border Tangent Constraint" << std::endl;
        AddBorderTangentConstraint(model, M3NodeVars,
                                   ExtrP.SoftBorderTangentM3);
      }

      if (!ExtrP.FeasiblityOnly)
        model.set(GRB_DoubleParam_MIPGap, ExtrP.Gap); // gap
      else
        model.set(GRB_DoubleParam_MIPGap, 1); // gap

      model.set(GRB_IntParam_NumericFocus, 0); // 1
      model.set(GRB_IntParam_Method, 1);
      model.set(GRB_IntParam_MIPFocus, 1);
      // NEW
      model.set(GRB_DoubleParam_Heuristics, 0.05);
      model.set(GRB_IntParam_Threads, 8);
      model.set(GRB_IntParam_Seed, 42);
      model.set(GRB_IntParam_AggFill, 5); // riduce fill-in nel presolve

      if (ExtrP.OptimizeObj)
        OptimizeCoefficients(model);

      model.optimize();

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
        // Field::FrameFieldGraph::OrthoNodesFrameField(int IndexN, int
        // &OrthoN)
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

public:
  void InitConstraints() {
    InitVariableMap();
    M3Solution = std::vector<bool>(M3VarToNode.size(), false);
    std::cout << "Num M3 Nodes Variables: " << M3VarToNode.size() << std::endl;

    // int MaxSize;
    // ScalarType AVGSize;

    SingleF_ConstrM3.clear();
    InitSingleF_ConstrM3();
    std::cout << "Initialized " << SingleF_ConstrM3.size()
              << " single face constraint" << std::endl;
    // ConstraintStats(SingleF_ConstrM3, MaxSize, AVGSize);
    // std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize <<
    // std::endl;

    TwiceEdge_ConstrM3.clear();
    InitTwiceEdge_Constr();
    std::cout << "Initialized " << TwiceEdge_ConstrM3.size()
              << " M3 twice edge constraint" << std::endl;
    // ConstraintStats(TwiceEdge_ConstrM3, MaxSize, AVGSize);
    // std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize <<
    // std::endl;

    EdgeLabel_Continuation_Constr.clear();
    Internal_TangentNodes_Constr.clear();
    Border_TangentNodes_Constr.clear();

    InitSingleLabelF_Constr();

    InitInternalTangentNodes_Constr();
    std::cout << "Initialized " << EdgeLabel_Continuation_Constr.size()
              << " edge label constraint" << std::endl;

    InitAvoidConnection_ConstrM3();
    std::cout << "Initialized " << Avoid_Connections.size()
              << " connection avoid constraint" << std::endl;
    // ConstraintStats(EdgeLabel_Continuation_Constr, MaxSize, AVGSize);
    // std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize <<
    // std::endl;

    InitBorderTangentNodes_Constr();
    std::cout << "Initialized " << Border_TangentNodes_Constr.size()
              << " border tangent constraint" << std::endl;
    // ConstraintStats(Border_TangentNodes_Constr, MaxSize, AVGSize);
    // std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize <<
    // std::endl;
  }

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

    // // concatenate to a single surface
    // std::vector<Geo::Point3<ScalarType>> NonManifV;
    // std::vector<std::vector<int>> NonManifFaces;
    // AppendMeshes<ScalarType>(NonManifV, NonManifFaces, FoilVertices,
    // FoilFaces); RemoveDuplicatedVert(NonManifV, NonManifFaces);

    // std::map<Geo::Point3<ScalarType>, int> PosIdxMap;
    // for (size_t i = 0; i < NonManifV.size(); i++)
    //   PosIdxMap[NonManifV[i]] = i;

    // SmoothReprojectSubMeshes(NonManifV, NonManifFaces, ReprojectV,
    //                          ReprojectFaces, FeatureE, Corners,
    //                          SmoothSteps);

    // // copy back positions
    // for (size_t i = 0; i < FoilVertices.size(); i++)
    //   for (size_t j = 0; j < FoilVertices[i].size(); j++) {
    //     Geo::Point3<ScalarType> VertP = FoilVertices[i][j];
    //     assert(PosIdxMap.count(VertP) > 0);
    //     int IndexNonManif = PosIdxMap[VertP];
    //     assert(IndexNonManif >= 0);
    //     assert(IndexNonManif < NonManifV.size());
    //     FoilVertices[i][j] = NonManifV[IndexNonManif];
    //   }

    FoilsSmooth<ScalarType>::SmoothFoils(FoilVertices, FoilFaces, ReprojectV,
                                         ReprojectFaces, FeatureE, Corners,
                                         SmoothSteps);
  }

  //   // only one single istance of the face can be chosen (one face might
  // // correspont to multiple nodes)
  // std::vector<std::vector<int>> SingleF_ConstrM3;

  // // the edge can be transversed only twice maximum (SET CORRECT VALENCE with
  // // Borders/Convex/Concave)
  // std::vector<std::vector<int>> TwiceEdge_ConstrM3;

  // // exit M3 should be only one for exit edge
  // std::vector<std::vector<int>> M3_Exit_Constr;

  // // Continuing M3 Label such that one face should trigger the next
  // std::vector<std::vector<int>> M3_Label_Cont_Constr;

  // // only one times for each direction for each node
  // std::vector<std::pair<int, std::vector<int>>> EdgeLabel_Constr;

  int NumConstraints() {
    int size_constr = SingleF_ConstrM3.size();
    size_constr += TwiceEdge_ConstrM3.size();
    size_constr += EdgeLabel_Continuation_Constr.size();
    size_constr += Internal_TangentNodes_Constr.size();
    size_constr += Border_TangentNodes_Constr.size();
    size_constr += Avoid_Connections.size();
    return size_constr;
  }

  void ConstraintsStats(ScalarType &AvgVal, int &MaxVal) {
    MaxVal = 0;
    AvgVal = 0;
    ScalarType Num = NumConstraints();

    for (size_t i = 0; i < SingleF_ConstrM3.size(); i++) {
      MaxVal = std::max(MaxVal, (int)SingleF_ConstrM3[i].size());
      AvgVal += (ScalarType)SingleF_ConstrM3[i].size();
    }
    for (size_t i = 0; i < TwiceEdge_ConstrM3.size(); i++) {
      MaxVal = std::max(MaxVal, (int)TwiceEdge_ConstrM3[i].size());
      AvgVal += (ScalarType)TwiceEdge_ConstrM3[i].size();
    }
    for (size_t i = 0; i < Internal_TangentNodes_Constr.size(); i++) {
      MaxVal = std::max(MaxVal, (int)Internal_TangentNodes_Constr[i].size());
      AvgVal += (ScalarType)Internal_TangentNodes_Constr[i].size();
    }
    for (size_t i = 0; i < Border_TangentNodes_Constr.size(); i++) {
      MaxVal = std::max(MaxVal, (int)Border_TangentNodes_Constr[i].size());
      AvgVal += (ScalarType)Border_TangentNodes_Constr[i].size();
    }

    for (size_t i = 0; i < EdgeLabel_Continuation_Constr.size(); i++) {
      MaxVal =
          std::max(MaxVal, (int)EdgeLabel_Continuation_Constr[i].second.size());
      AvgVal += (ScalarType)EdgeLabel_Continuation_Constr[i].second.size();
    }

    for (size_t i = 0; i < Avoid_Connections.size(); i++) {
      MaxVal = std::max(MaxVal, 2);
      AvgVal += (ScalarType)2;
    }
    if (Num > 0)
      AvgVal /= Num;
  }

  int NumVariables() { return M3VarToNode.size(); }

  bool Solve(const ExtractParam &ExtrP, std::vector<int> &ChosenF,
             std::vector<bool> &MatchedSeed)
  // bool cut_space = false)
  {
    // InitConstraints();

    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<int> DisabledN;
    std::vector<int> Seeds = ExtrP.SoftNodes;
    Seeds.insert(Seeds.end(), ExtrP.HardNodes.begin(), ExtrP.HardNodes.end());
    MatchedSeed.clear();

    if (ExtrP.SoftNodes.size() > 0)
      MatchedSeed = std::vector<bool>(ExtrP.SoftNodes.size(), false);
    else
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);

    // exit(0);
    auto time1 = std::chrono::high_resolution_clock::now();

    std::vector<int> M3NodesSel;
    bool solved = SolveGurobi(ExtrP, ChosenF, M3NodesSel);

    auto time2 = std::chrono::high_resolution_clock::now();

    std::set<int> M3NodesSelSet(M3NodesSel.begin(), M3NodesSel.end());
    int NumMatched = 0;

    if (ExtrP.SoftNodes.size() > 0) {
      for (size_t i = 0; i < ExtrP.SoftNodes.size(); i++) {
        int TestNode = ExtrP.SoftNodes[i];
        if (M3NodesSelSet.count(TestNode) == 0)
          continue;
        MatchedSeed[i] = true;
        NumMatched++;
      }
    } else {
      for (size_t i = 0; i < ExtrP.HardNodes.size(); i++) {
        int TestNode = ExtrP.HardNodes[i];
        if (M3NodesSelSet.count(TestNode) == 0)
          continue;
        MatchedSeed[i] = true;
        NumMatched++;
      }
    }

    auto time3 = std::chrono::high_resolution_clock::now();

    std::cout << "Matched "
              << 100 * (ScalarType)NumMatched /
                     (ScalarType)ExtrP.SoftNodes.size()
              << " Perc of Soft Seeds" << std::endl;
    // if (cut_space) {
    //   for (size_t i = 0; i < DisabledN.size(); i++)
    //     M3Graph.EnableNode(DisabledN[i]);
    // }

    auto duration0 =
        std::chrono::duration_cast<std::chrono::seconds>(time1 - time0);
    auto duration1 =
        std::chrono::duration_cast<std::chrono::seconds>(time2 - time1);
    auto duration2 =
        std::chrono::duration_cast<std::chrono::seconds>(time3 - time2);

    // Output the duration
    // std::cout << "** Time Presolve : " << duration0.count()<< " **
    // seconds\n";
    std::cout << "** Time Setup : " << duration0.count() << " ** seconds\n";
    std::cout << "** Time Solve : " << duration1.count() << " ** seconds\n";
    std::cout << "** Time Collect : " << duration2.count() << " ** seconds\n";

    return solved;
  }

  SheetExtractorILPM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                      std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif