#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_ILP
#define FRAME_FIELD_GRAPH_SHEET_TRACING_ILP

#include "graph.h"
#include "graph_frame_field.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <include/gurobi_c++.h>
#include <iostream>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

template <class ScalarType> class SheetExtractorILP {

  Geo::UndirectedNodeWeightedGraph<ScalarType> &Graph;
  std::vector<std::vector<int>> &InternalFaces;

  std::vector<int> NodeToVar;
  std::vector<int> VarToNode;

  // only one single istance of the face can be chosen (one face might
  // correspont to multiple nodes)
  std::vector<std::vector<int>> SingleF_Constr;

  // only one times for each direction for each node
  std::vector<std::pair<int, std::vector<int>>> EdgeLabel_Constr;

  // for the border edges a different kind of constraint
  std::vector<std::vector<int>> EdgeLabel_ExitConstr;

  // the connections that cannot be done
  std::vector<std::pair<int, int>> IncompatibleNodePairs_Constr;

  // the edge can be transversed only twice maximum
  std::vector<std::vector<int>> TwiceEdge_Constr;

  bool IsVariable(int IndexN) const {
    if (!Graph.IsEnabled(IndexN))
      return false;
    if (Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    return true;
  }

  void InitVariableMap() {
    NodeToVar = std::vector<int>(Graph.NumNodes(), -1);
    VarToNode.clear();

    int CurrV = 0;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (!IsVariable(i))
        continue;

      NodeToVar[i] = CurrV;
      VarToNode.push_back(i);
      CurrV++;
    }
  }

  void InitSingleF_Constr() {
    SingleF_Constr.clear();
    for (size_t i = 0; i < InternalFaces.size(); i++) {
      std::vector<int> NodesIdx;
      Field::FrameFieldM3Graph::NodesOfElem(i, NodesIdx);
      std::vector<int> VarNodesForF;
      for (size_t j = 0; j < NodesIdx.size(); j++) {
        int testN = NodesIdx[j];
        if (!IsVariable(testN))
          continue;
        VarNodesForF.push_back(testN);
      }
      // no need if there is zero or none active
      // variable for such a face
      if (VarNodesForF.size() <= 1)
        continue;

      // otherwise add a new constraint
      SingleF_Constr.resize(SingleF_Constr.size() + 1);
      for (size_t j = 0; j < VarNodesForF.size(); j++) {
        int IndexN = VarNodesForF[j];
        assert(IsVariable(IndexN));
        assert(IndexN >= 0);
        assert(IndexN < NodeToVar.size());
        int IndexVar = NodeToVar[IndexN];
        assert(IndexVar >= 0);
        SingleF_Constr.back().push_back(IndexVar);
      }
    }
  }

  void InitIncompatibleNodePairs_Constr() {
    IncompatibleNodePairs_Constr.clear();
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (!IsVariable(i))
        continue;
      int IndexVar0 = NodeToVar[i];
      assert(IndexVar0 >= 0);
      for (size_t j = 0; j < Graph.NumNeigh(i); j++) {
        if (Graph.IsEnabledL(i, j))
          continue;
        int IndexN = Graph.NeighNode(i, j);
        if (!IsVariable(IndexN))
          continue;
        int IndexVar1 = NodeToVar[IndexN];
        assert(IndexVar1 >= 0);
        std::pair<int, int> ConstrPair(std::min(IndexVar0, IndexVar1),
                                       std::max(IndexVar0, IndexVar1));
        IncompatibleNodePairs_Constr.push_back(ConstrPair);
      }
    }
    std::sort(IncompatibleNodePairs_Constr.begin(),
              IncompatibleNodePairs_Constr.end());
    IncompatibleNodePairs_Constr.erase(
        std::unique(IncompatibleNodePairs_Constr.begin(),
                    IncompatibleNodePairs_Constr.end()),
        IncompatibleNodePairs_Constr.end());
  }

  void InitSingleLabelF_Constr() {

    // //get border edges
    // std::vector<int> BorderF,BorderE;
    // Geo::getBorderEdges(InternalFaces,BorderF,BorderE);
    // std::set<std::pair<int, int> > BorderESet;
    // for (size_t i=0;i<BorderF.size();i++)
    //     BorderESet.push_back(std::pair<int, int>(BorderF[i],BorderE[i]));

    EdgeLabel_Constr.clear();
    EdgeLabel_ExitConstr.clear();

    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (!IsVariable(i))
        continue;
      std::vector<std::vector<int>> LabelNodes;
      int IndexNode0 = i;
      int IndexVar0 = NodeToVar[i];
      assert(IndexVar0 >= 0);

      // one label for each edge
      LabelNodes.resize(3);
      for (size_t j = 0; j < Graph.NumNeigh(i); j++) {
        int IndexN = Graph.NeighNode(i, j);
        int IndexLab = Graph.NeighLabel(i, j);
        assert(IndexLab >= 0);
        assert(IndexLab < 3);
        if (!IsVariable(IndexN))
          continue;
        int IndexVar = NodeToVar[IndexN];
        assert(IndexVar >= 0);
        LabelNodes[IndexLab].push_back(IndexVar);
      }

      for (size_t j = 0; j < LabelNodes.size(); j++) {
        if (LabelNodes[j].size() == 0)
          continue;
        std::set<int> LabNodesSet(LabelNodes[j].begin(), LabelNodes[j].end());
        assert(LabNodesSet.size() == LabelNodes[j].size());

        bool IsExitNodeLab = Graph.IsExitLabel(IndexNode0, j);
        // if no exit separate them so different constraint
        if (!IsExitNodeLab) {
          EdgeLabel_Constr.resize(EdgeLabel_Constr.size() + 1);
          EdgeLabel_Constr.back().first = IndexVar0;
          EdgeLabel_Constr.back().second = LabelNodes[j];
        } else {
          // otherwise there should be only one
          EdgeLabel_ExitConstr.push_back(LabelNodes[j]);
          EdgeLabel_ExitConstr.back().push_back(IndexVar0);
        }
      }
    }
  }

  //   void InitExitLabel_Constr()
  //   {
  //     EdgeLabel_ExitConstr.clear();
  //     for (size_t i = 0; i < Graph.NumNodes(); i++) {
  //       if (!IsVariable(i))
  //         continue;
  //       std::vector<std::vector<int>> LabelNodes;
  //       int IndexVar0 = NodeToVar[i];
  //       assert(IndexVar0 >= 0);

  //       // one label for each edge
  //       LabelNodes.resize(3);
  //       for (size_t j = 0; j < Graph.NumNeigh(i); j++) {
  //         int IndexN = Graph.NeighNode(i, j);
  //         int IndexLab = Graph.NeighLabel(i, j);
  //         assert(IndexLab >= 0);
  //         assert(IndexLab < 3);
  //         if (!IsVariable(IndexN))
  //           continue;
  //         int IndexVar = NodeToVar[IndexN];
  //         assert(IndexVar >= 0);
  //         LabelNodes[IndexLab].push_back(IndexVar);
  //       }
  //   }

  void InitTwiceEdge_Constr() {
    std::map<std::pair<int, int>, std::vector<int>> EdgeNodes;
    TwiceEdge_Constr.clear();
    for (size_t i = 0; i < InternalFaces.size(); i++)
      for (size_t j = 0; j < InternalFaces[i].size(); j++) {
        int IndexV0 = InternalFaces[i][j];
        int IndexV1 = InternalFaces[i][(j + 1) % InternalFaces[i].size()];
        std::vector<int> NodesIdx;
        Field::FrameFieldM3Graph::NodesOfElem(i, NodesIdx);
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        for (size_t k = 0; k < NodesIdx.size(); k++) {
          int IndexN = NodesIdx[k];
          if (!IsVariable(IndexN))
            continue;
          int IndexVar = NodeToVar[IndexN];
          assert(IndexVar >= 0);
          EdgeNodes[Key].push_back(IndexVar);
        }
      }

    std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
    for (IteMap = EdgeNodes.begin(); IteMap != EdgeNodes.end(); IteMap++) {
      std::vector<int> NodesOfEdge = (*IteMap).second;
      std::sort(NodesOfEdge.begin(), NodesOfEdge.end());
      NodesOfEdge.erase(std::unique(NodesOfEdge.begin(), NodesOfEdge.end()),
                        NodesOfEdge.end());
      if (NodesOfEdge.size() < 3)
        continue;
      TwiceEdge_Constr.push_back(NodesOfEdge);
    }
  }

  // the initial set of contraints
  //   std::vector<GRBConstr> InitialConstraints;
  //   std::unique_ptr<GRBModel> model;
  //   std::vector<GRBVar> Vars;

  //   void InitGurobi() {
  //     ScalarType DeltaMin = 0.00000001;
  //     try {
  //       // Create Gurobi environment
  //       GRBEnv env = GRBEnv(true);
  //       env.set("OutputFlag", "0"); // Set to 1 for detailed Gurobi output
  //       env.start();
  //       // env.set("OutputFlag", "0"); // Set to 1 for detailed Gurobi output

  //       // Create Gurobi model
  //       // model = GRBModel(env);
  //       model = std::make_unique<GRBModel>(env);
  //       model->reset();
  //       // model->remove(model->getConstrs());

  //       // add the boolean variables
  //       Vars = std::vector<GRBVar>(VarToNode.size());
  //       for (int i = 0; i < VarToNode.size(); ++i) {
  //         Vars[i] =
  //             model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" +
  //             std::to_string(i));
  //       }

  //       // Set the objective: minimize the sum of the variables
  //       GRBLinExpr objective = 0;
  //       for (int i = 0; i < VarToNode.size(); ++i) {
  //         int IndexN = VarToNode[i];
  //         assert(IndexN >= 0);
  //         assert(IndexN < Graph.NumNodes());
  //         ScalarType CurrW = Graph.Weight(IndexN);
  //         objective += Vars[i] * (DeltaMin + CurrW);
  //       }
  //       model->setObjective(objective, GRB_MINIMIZE);

  //       // ADD SINGLE NODE FOR A FACE CONSTRAINT
  //       for (size_t i = 0; i < SingleF_Constr.size(); i++) {
  //         GRBLinExpr lin_constr;
  //         for (size_t j = 0; j < SingleF_Constr[i].size(); j++) {
  //           int IndexV = SingleF_Constr[i][j];
  //           lin_constr += Vars[IndexV];
  //         }
  //         model->addConstr(lin_constr <= 1, "Single F" + std::to_string(i));
  //       }

  //       // ADD EDGE LABEL CONSTRAINT
  //       for (size_t i = 0; i < EdgeLabel_Constr.size(); i++) {
  //         GRBLinExpr upper_bound;
  //         int IndexVar0 = EdgeLabel_Constr[i].first;
  //         assert(IndexVar0 < Vars.size());
  //         assert(IndexVar0 >= 0);
  //         assert(EdgeLabel_Constr[i].second.size() > 0);
  //         for (size_t j = 0; j < EdgeLabel_Constr[i].second.size(); j++) {
  //           int IndexVar1 = EdgeLabel_Constr[i].second[j];
  //           assert(IndexVar0 != IndexVar1);
  //           assert(IndexVar1 < Vars.size());
  //           assert(IndexVar1 >= 0);
  //           upper_bound += Vars[IndexVar1];
  //         }
  //         model->addGenConstrIndicator(Vars[IndexVar0], 1, upper_bound <= 1,
  //                                      "UpperBound" + std::to_string(i));
  //         // model.addConstr(upper_bound <= 1, "UpperBound"+
  //         std::to_string(i)); model->addConstr(Vars[IndexVar0] <=
  //         upper_bound,
  //                          "LowerBound" + std::to_string(i));
  //       }

  //       // ADD INCOMPATIBLE NODE CONSTRAINT
  //       for (size_t i = 0; i < IncompatibleNodePairs_Constr.size(); i++) {
  //         int IndexVar0 = IncompatibleNodePairs_Constr[i].first;
  //         int IndexVar1 = IncompatibleNodePairs_Constr[i].second;
  //         model->addConstr(Vars[IndexVar0] + Vars[IndexVar1] <= 1,
  //                          "Incompatible" + std::to_string(i));
  //       }

  //       // ADD TWICE EDGE CONSTR
  //       for (size_t i = 0; i < TwiceEdge_Constr.size(); i++) {
  //         GRBLinExpr twice_constr;
  //         for (size_t j = 0; j < TwiceEdge_Constr[i].size(); j++) {
  //           int IndexVar = TwiceEdge_Constr[i][j];
  //           assert(IndexVar >= 0);
  //           assert(IndexVar < Vars.size());
  //           twice_constr += Vars[IndexVar];
  //         }
  //         model->addConstr(twice_constr <= 2, "Twice" + std::to_string(i));
  //       }

  //       model->optimize();

  //     }

  //     catch (GRBException &e) {
  //       std::cerr << "Gurobi error: " << e.getMessage() << "\n";
  //     } catch (const std::exception &e) {
  //       std::cerr << "Standard exception: " << e.what() << "\n";
  //     }
  //   }

  //   void UpdateViolatedEdgeLabel()
  //   {
  //     if (Solution.size()<VarToNode.size())return;
  //     std::set<int> ChoosenSet(ChosenF.begin(),ChosenF.end());
  //     int violated_lower=0;
  //     int violated_upper=0;
  //     for (size_t i = 0; i < EdgeLabel_Constr.size(); i++) {

  //         int Lab0=0;
  //         int IndexVar0 = EdgeLabel_Constr[i].first;
  //         assert(IndexVar0>=0);
  //         assert(IndexVar0<Solution.size());

  //         if (Solution[IndexVar0])
  //             Lab0=1;
  //         else
  //             continue;

  //         int NumLab=0;
  //         for (size_t j = 0; j < EdgeLabel_Constr[i].second.size(); j++) {
  //           int IndexVar1 = EdgeLabel_Constr[i].second[j];
  //           assert(IndexVar1>=0);
  //           assert(IndexVar1<Solution.size());
  //           if (Solution[IndexVar1])NumLab++;
  //         }
  //         if (NumLab>1)
  //         std::cout<<"Label Num:"<<NumLab<<std::endl;
  //         if (NumLab>1)violated_upper++;
  //         if (Lab0>NumLab)violated_lower++;
  //       }
  //       std::cout<<"VIOLATED LABEL CONSTRAINT
  //       (Upperbound):"<<violated_upper<<std::endl; std::cout<<"VIOLATED LABEL
  //       CONSTRAINT (Lowerbound):"<<violated_lower<<std::endl;
  //   }

  bool SolveGurobi(const std::vector<int> &SeedNodes,
                   std::vector<int> &ChosenF) {
    ChosenF.clear();
    ScalarType DeltaMin = 0.00000001;
    try {
      // Create Gurobi environment
      GRBEnv env = GRBEnv(true);
      env.set("OutputFlag", "0"); // Set to 1 for detailed Gurobi output
      env.start();

      // Create Gurobi model
      GRBModel model = GRBModel(env);
      model.reset();

      // add the boolean variables
      std::vector<GRBVar> Vars(VarToNode.size());
      for (int i = 0; i < VarToNode.size(); ++i)
        Vars[i] =
            model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + std::to_string(i));

      // Set the objective: minimize the sum of the variables
      GRBLinExpr objective = 0;
      for (int i = 0; i < VarToNode.size(); ++i) {
        int IndexN = VarToNode[i];
        assert(IndexN >= 0);
        assert(IndexN < Graph.NumNodes());
        ScalarType CurrW = Graph.Weight(IndexN);
        objective += Vars[i] * (DeltaMin + CurrW);
      }
      model.setObjective(objective, GRB_MINIMIZE);

      // ADD SINGLE NODE FOR A FACE CONSTRAINT
      for (size_t i = 0; i < SingleF_Constr.size(); i++) {
        GRBLinExpr lin_constr;
        for (size_t j = 0; j < SingleF_Constr[i].size(); j++) {
          int IndexV = SingleF_Constr[i][j];
          lin_constr += Vars[IndexV];
        }
        model.addConstr(lin_constr <= 1, "Single F" + std::to_string(i));
      }

      // ADD EDGE LABEL CONSTRAINT
      for (size_t i = 0; i < EdgeLabel_Constr.size(); i++) {
        GRBLinExpr upper_bound;
        int IndexVar0 = EdgeLabel_Constr[i].first;
        assert(IndexVar0 >= 0);
        assert(IndexVar0 < VarToNode.size());
        int IndexN0 = VarToNode[IndexVar0];
        // int IndexF0=Field::FrameFieldM3Graph::ElemOfNode(IndexN0);

        assert(IndexVar0 < Vars.size());
        assert(IndexVar0 >= 0);
        assert(EdgeLabel_Constr[i].second.size() > 0);
        for (size_t j = 0; j < EdgeLabel_Constr[i].second.size(); j++) {
          int IndexVar1 = EdgeLabel_Constr[i].second[j];
          assert(IndexVar0 != IndexVar1);
          assert(IndexVar1 < Vars.size());
          assert(IndexVar1 >= 0);
          upper_bound += Vars[IndexVar1];
        }
        model.addGenConstrIndicator(Vars[IndexVar0], 1, upper_bound <= 1,
                                    "UpperBound" + std::to_string(i));
        // model.addConstr(upper_bound <= 1, "UpperBound"+ std::to_string(i));
        model.addConstr(Vars[IndexVar0] <= upper_bound,
                        "LowerBound" + std::to_string(i));
      }

      // ADD INCOMPATIBLE NODE CONSTRAINT
      for (size_t i = 0; i < IncompatibleNodePairs_Constr.size(); i++) {
        int IndexVar0 = IncompatibleNodePairs_Constr[i].first;
        int IndexVar1 = IncompatibleNodePairs_Constr[i].second;
        model.addConstr(Vars[IndexVar0] + Vars[IndexVar1] <= 1,
                        "Incompatible" + std::to_string(i));
      }

      // ADD TWICE EDGE CONSTR
      for (size_t i = 0; i < TwiceEdge_Constr.size(); i++) {
        GRBLinExpr twice_constr;
        for (size_t j = 0; j < TwiceEdge_Constr[i].size(); j++) {
          int IndexVar = TwiceEdge_Constr[i][j];
          assert(IndexVar >= 0);
          assert(IndexVar < Vars.size());
          twice_constr += Vars[IndexVar];
        }
        model.addConstr(twice_constr <= 4, "Twice" + std::to_string(i));
      }

      //   // ADD INITIAL SEED NODES
      //   for (size_t i = 0; i < SeedNodes.size(); i++) {
      //     if (!IsVariable(SeedNodes[i]))
      //       continue;
      //     int IndexVarTest = NodeToVar[SeedNodes[i]];
      //     assert(IndexVarTest >= 0);
      //     GRBConstr AddedC = model.addConstr(Vars[IndexVarTest] == 1,
      //                                        "Initial Seed" +
      //                                        std::to_string(i));
      //   }

      int IndexVarTest = NodeToVar[377];
      assert(IndexVarTest >= 0);
      GRBConstr AddedC0 = model.addConstr(Vars[IndexVarTest] == 1,
                                          "Initial Seed" + std::to_string(1));
      IndexVarTest = NodeToVar[50001];
      assert(IndexVarTest >= 0);
      GRBConstr AddedC1 = model.addConstr(Vars[IndexVarTest] == 1,
                                          "Initial Seed" + std::to_string(2));
     
     IndexVarTest = NodeToVar[99803];
       GRBConstr AddedC2 = model.addConstr(Vars[IndexVarTest] == 1,
                                          "Initial Seed" + std::to_string(3));

        
      // Record the first timestamp
      //model.presolve();
      auto start = std::chrono::high_resolution_clock::now();
      
    //   model.set(GRB_DoubleParam_Heuristics, 0.5);
    //   model.set(GRB_DoubleParam_MIPGap, 0.9); // 1% gap
    //   model.set(GRB_DoubleParam_MIPGapAbs, 10); // Absolute gap
    //    // Set time limit to 300 seconds (5 minutes)
    // //model.set(GRB_DoubleParam_TimeLimit, 3.0);
    //   model.set(GRB_IntParam_MIPFocus, 1);
    //   model.set(GRB_IntParam_Threads, 8); // Use 4 threads
    //   model.set(GRB_IntParam_Cuts, 2); // Moderate cut aggressiveness
    //   model.set(GRB_IntParam_NodeMethod, 2);
      model.optimize();
      auto end = std::chrono::high_resolution_clock::now();
      // Calculate the duration
      auto duration =
          std::chrono::duration_cast<std::chrono::seconds>(end - start);
      // Output the duration
      std::cout << "Time elapsed: " << duration.count() << " seconds\n";
      // COLLECT THE RESULT
      
      if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        std::cout << "Optimal solution found:\n";

        // Read back the solution for the variables
        // std::vector<int> solution(n);
        Solution.clear();
        for (int i = 0; i < Vars.size(); ++i) {
          double ValueVar = Vars[i].get(GRB_DoubleAttr_X);
          // int ValueVar = static_cast<int>(Vars[i].get(GRB_DoubleAttr_X));

          if (ValueVar > 0.5) {
            Solution.push_back(true);
            int IndexN = VarToNode[i];
            int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
            ChosenF.push_back(IndexF);
          } else {
            Solution.push_back(false);
          }
          // solution[i] = static_cast<int>(x[i].get(GRB_DoubleAttr_X));
          // std::cout << "x[" << i << "] = " << solution[i] << "\n";
        }
        std::cout << "Chosen " << ChosenF.size() << " faces" << std::endl;
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
  }

  //  void InitBorderLabelN()
  //  {
  //     BorderLabelN.clear();
  //     std::vector<std::vector<std::vector<int>>> NextF, NextE;
  //     Geo::getFFAdjacencyNonManifold(InternalFaces, NextF, NextE);
  //     for (size_t i=0;i<Graph.NumNodes();i++)
  //     {
  //         int IndexN=i;
  //         int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
  //         assert(IndexF>=0);
  //         assert(IndexF<NextF.size());
  //         assert(NextF[IndexF].size()==3);
  //         BorderLabelN.push_back(std::vector<bool>(3,false));
  //         for (size_t e=0;e<NextF[IndexF].size();e++)
  //             BorderLabelN.back()[e]=(NextF[IndexF][e]==-1);
  //     }
  //  }

public:
  // std::vector<std::vector<int>> ViolatedC;
  std::vector<bool> Solution;
  std::vector<std::vector<bool>> BorderLabelN;
  // std::vector<int> ChosenF;

  bool FindSurfaceSheet(std::vector<int> &SeedNodes, std::vector<int> &ChosenF,
                        bool SingleFConstr = true, bool ExitLabelConstr = true,
                        bool IncompatibleNPairConstr = true,
                        bool TwiceEdgeConstr = true) {

    SingleF_Constr.clear();
    EdgeLabel_Constr.clear();
    IncompatibleNodePairs_Constr.clear();
    TwiceEdge_Constr.clear();

    std::cout << "Initializing Variable Mapping" << std::endl;
    InitVariableMap();

    std::cout << "Num Variables: " << VarToNode.size() << std::endl;

    if (SingleFConstr) {
      InitSingleF_Constr();
      std::cout << "Initialized " << SingleF_Constr.size()
                << " single face constraint" << std::endl;
    }

    if (ExitLabelConstr) {
      InitSingleLabelF_Constr();
      std::cout << "Initialized " << EdgeLabel_Constr.size()
                << " face label constraint" << std::endl;
    }

    if (IncompatibleNPairConstr) {
      InitIncompatibleNodePairs_Constr();
      std::cout << "Initialized " << IncompatibleNodePairs_Constr.size()
                << " Incompatible faces constraint" << std::endl;
    }

    if (TwiceEdgeConstr) {
      InitTwiceEdge_Constr();
      std::cout << "Initialized " << TwiceEdge_Constr.size()
                << " twice edge constraint" << std::endl;
    }

    bool solved = SolveGurobi(SeedNodes, ChosenF);
    return solved;
  }

  //   bool SolveSurfaceFrom(std::vector<int> &SeedNodes,
  //                         std::vector<int> &ChosenF) {
  //     ChosenF.clear();
  //     if (SeedNodes.size() == 0)
  //       return false;

  //     // remove the old ones
  //     for (size_t i = 0; i < InitialConstraints.size(); i++)
  //       model->remove(InitialConstraints[i]);

  //     InitialConstraints.clear();

  //     // add the new ones
  //     for (size_t i = 0; i < SeedNodes.size(); i++) {
  //       if (!IsVariable(SeedNodes[i]))
  //         continue;
  //       int IndexVarTest = NodeToVar[SeedNodes[i]];
  //       assert(IndexVarTest >= 0);
  //       GRBConstr AddedC = model->addConstr(Vars[IndexVarTest] == 1,
  //                                           "Initial Seed" +
  //                                           std::to_string(i));
  //       InitialConstraints.push_back(AddedC);
  //     }
  //     if (SeedNodes.size() == 0)
  //       return true;
  //     // solve
  //     model->optimize();

  //     // collect the resul
  //     if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
  //       std::cout << "Optimal solution found:\n";

  //       // Read back the solution for the variables
  //       // std::vector<int> solution(n);
  //       Solution.clear();
  //       for (int i = 0; i < Vars.size(); ++i) {
  //         double ValueVar = Vars[i].get(GRB_DoubleAttr_X);
  //         // int ValueVar = static_cast<int>(Vars[i].get(GRB_DoubleAttr_X));

  //         if (ValueVar > 0.5) {
  //           Solution.push_back(true);
  //           int IndexN = VarToNode[i];
  //           int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
  //           ChosenF.push_back(IndexF);
  //         } else {
  //           Solution.push_back(false);
  //         }
  //         // solution[i] = static_cast<int>(x[i].get(GRB_DoubleAttr_X));
  //         // std::cout << "x[" << i << "] = " << solution[i] << "\n";
  //       }
  //       std::cout << "Chosen " << ChosenF.size() << " faces" << std::endl;
  //       return true;
  //     } else {
  //       std::cout << "No optimal solution found.\n";
  //       return false;
  //     }
  //   }

  //   void InitILP(bool SingleFConstr = true, bool ExitLabelConstr = true,
  //                bool IncompatibleNPairConstr = true,
  //                bool TwiceEdgeConstr = true) {
  //     std::cout << "Initializing Variable Mapping" << std::endl;
  //     InitVariableMap();

  //     std::cout << "Num Variables: " << VarToNode.size() << std::endl;

  //     if (SingleFConstr) {
  //       InitSingleF_Constr();
  //       std::cout << "Initialized " << SingleF_Constr.size()
  //                 << " single face constraint" << std::endl;
  //     }

  //     if (ExitLabelConstr) {
  //       InitSingleLabelF_Constr();
  //       std::cout << "Initialized " << EdgeLabel_Constr.size()
  //                 << " face label constraint" << std::endl;
  //     }

  //     if (IncompatibleNPairConstr) {
  //       InitIncompatibleNodePairs_Constr();
  //       std::cout << "Initialized " << IncompatibleNodePairs_Constr.size()
  //                 << " Incompatible faces constraint" << std::endl;
  //     }

  //     if (TwiceEdgeConstr) {
  //       InitTwiceEdge_Constr();
  //       std::cout << "Initialized " << TwiceEdge_Constr.size()
  //                 << " twice edge constraint" << std::endl;
  //     }

  //     InitGurobi();

  //     // UpdateViolatedEdgeLabel();
  //   }

  SheetExtractorILP(Geo::UndirectedNodeWeightedGraph<ScalarType> &_Graph,
                    std::vector<std::vector<int>> &_InternalFaces)
      : Graph(_Graph), InternalFaces(_InternalFaces) {}
};

#endif
