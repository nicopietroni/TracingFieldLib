#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_ILP
#define FRAME_FIELD_GRAPH_SHEET_TRACING_ILP

#include "../mesh_functions.h"
#include "../connectivity.h"
#include "graph.h"
#include "graph_frame_field.h"
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

#define SOFT_SEED_CONST 10000

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
        if (!Graph.IsEnabledL(i, j))continue;

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

  bool SolveGurobi(const std::vector<int> &HardConstrNodes,
                   const std::vector<int> &SoftConstrNodes,
                   std::vector<int> &ChosenF) {
    ChosenF.clear();
    ScalarType DeltaMin = 0.00000001;
    try {
      // Create Gurobi environment
      GRBEnv env = GRBEnv(true);
      // env.set("OutputFlag", "0"); // Set to 1 for detailed Gurobi output
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

      // ADD SOFT CONSTRAINTS
      for (int i = 0; i < SoftConstrNodes.size(); ++i) {
        int IndexN = SoftConstrNodes[i];
        if (!IsVariable(IndexN))
          continue;
        int IndexVar = NodeToVar[IndexN];
        assert(IndexVar >= 0);
        assert(IndexVar < Vars.size());
        objective += -SOFT_SEED_CONST * Vars[IndexVar];
        // GRBConstr AddedC = model.addConstr(Vars[IndexVarTest] == 1,
        //  assert(IndexN >= 0);
        //  assert(IndexN < Graph.NumNodes());
        //  ScalarType CurrW = Graph.Weight(IndexN);
        //  objective += Vars[i] * (DeltaMin + CurrW);
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

      // ADD EXIT LABEL CONSTR
      for (size_t i = 0; i < EdgeLabel_ExitConstr.size(); i++) {
        if (EdgeLabel_ExitConstr.size() <= 1)
          continue;
        GRBLinExpr upper_bound;

        for (size_t j = 0; j < EdgeLabel_ExitConstr[i].size(); j++) {
          int IndexVar = EdgeLabel_ExitConstr[i][j];
          assert(IndexVar < Vars.size());
          assert(IndexVar >= 0);
          upper_bound += Vars[IndexVar];
        }
        model.addConstr(upper_bound <= 1, "ExitB " + std::to_string(i));
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

      // ADD HARD CONSTRAINTS
      for (size_t i = 0; i < HardConstrNodes.size(); i++) {
        if (!IsVariable(HardConstrNodes[i]))
          continue;
        int IndexVarTest = NodeToVar[HardConstrNodes[i]];
        assert(IndexVarTest >= 0);
        GRBConstr AddedC = model.addConstr(Vars[IndexVarTest] == 1,
                                           "Hard Constr" + std::to_string(i));
      }

      // Record the first timestamp
      // model.presolve();
      auto start = std::chrono::high_resolution_clock::now();

      // model.set(GRB_DoubleParam_Heuristics, 0.5);
      model.set(GRB_DoubleParam_MIPGap, 0.05); // 1% gap
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
    for (size_t i = 0; i < Solution.size(); i++) {
      if (!Solution[i])
        continue;
      int IndexN = VarToNode[i];
      SelectedNodes.insert(IndexN);
    }

    // split into parts
    std::set<int> VisitedNodes;
    std::vector<std::vector<int>> Groups;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
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
        for (size_t j = 0; j < Graph.NumNeigh(CurrN); j++) {
          // if not enabled link
          if (!Graph.IsEnabledL(CurrN, j))
            continue;
          int NeighNode = Graph.NeighNode(CurrN, j);

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
      const std::vector<std::pair<int,int>> &FeatureE,
      const std::vector<int> &Corners,
      int SmoothSteps)
  {
   
    //get corner and features position
    std::set<Geo::Point3<ScalarType>> FeatureVPos,CornerVPos;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    for (size_t i=0;i<FeatureE.size();i++)
    {
      int IndexV0=FeatureE[i].first;
      int IndexV1=FeatureE[i].second;
      Geo::Point3<ScalarType> Pos0=ReprojectV[IndexV0];
      Geo::Point3<ScalarType> Pos1=ReprojectV[IndexV1];
      FeatureVPos.insert(Pos0);
      FeatureVPos.insert(Pos1);
      Geo::Segment3<ScalarType> Seg(Pos0,Pos1);
      FeatureSegs.push_back(Seg);
    }

    for (size_t i=0;i<Corners.size();i++)
    {
      int IndexV=Corners[i];
      CornerVPos.insert(ReprojectV[IndexV]);
    }

    //then get the indexes on tne non manif 
    std::set<int> FeatureIdxSet,CornerIdxSet;
    for (size_t i=0;i<NonManifV.size();i++)
    {
      Geo::Point3<ScalarType> TestP=NonManifV[i];
      if (FeatureVPos.count(TestP)>0)FeatureIdxSet.insert(i);
      if (CornerVPos.count(TestP)>0)CornerIdxSet.insert(i);
    }

    std::vector<int> FeatureIdxVert(FeatureIdxSet.begin(),FeatureIdxSet.end());
    std::vector<int> CornerIdxVert(CornerIdxSet.begin(),CornerIdxSet.end());

    std::vector<Geo::Point3<ScalarType>> PrevPos=NonManifV;
    
    //create the reprojection data structure
    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    TriGrid.Init(ReprojectV, ReprojectFaces);

    //get the ones to be reprojected
    std::vector<size_t> BorderV;
    Geo::getBorderVertNonManifMesh(NonManifFaces,BorderV);
   
    for (size_t s = 0; s < SmoothSteps; s++) {
      
      Geo::SmoothLaplacianNonManif<ScalarType>(NonManifV, NonManifFaces, 1);

      //reproject regular
      for (size_t i = 0; i < BorderV.size(); i++) {
        int IndexV = BorderV[i];

        if (FeatureIdxSet.count(IndexV))continue;
        if (CornerIdxSet.count(IndexV))continue;
        
        assert(IndexV >= 0);
        assert(IndexV < NonManifV.size());
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        int ClosestF;
        Geo::Point3<ScalarType> ClosestPos;
        TriGrid.GridClosest(ReprojectV, ReprojectFaces, Pos,
                            TriGrid.bbox.Diag(), ClosestF, ClosestPos);
        NonManifV[IndexV] = ClosestPos;
      }

      //reproject features
      for (size_t i=0;i<FeatureIdxVert.size();i++)
      {
        int IndexV=FeatureIdxVert[i];
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        ScalarType MinD=std::numeric_limits<ScalarType>::max();
        int BestSeg=-1;
        Geo::Point3<ScalarType> BestClos =Pos;
        for (size_t j=0;j<FeatureSegs.size();j++)
        {
          Geo::Segment3<ScalarType> Seg=FeatureSegs[j];
          Geo::Point3<ScalarType> TestClos=Geo::ClosestPoint(Seg,Pos);
          ScalarType TestD=(TestClos-Pos).Norm();
          if (TestD>MinD)continue;
          MinD=TestD;
          BestSeg=j;
          BestClos=TestClos;
        }
        if (BestSeg==-1)continue;
        NonManifV[IndexV]=BestClos;
      }

      //put back corners on original positions
      for (size_t i=0;i<CornerIdxVert.size();i++)
      {
        int IndexV=CornerIdxVert[i];
        NonManifV[IndexV]=PrevPos[IndexV];
      }
    }
  }

public:
  struct ExtractParam {
    std::vector<int> HardNodes;
    std::vector<int> SoftNodes;
    bool SingleFConstr;
    bool SingleLabelCostr;
    bool IncompatibleNPairConstr;
    bool TwiceEdgeConstr;

    ExtractParam() {
      SingleFConstr = true;
      SingleLabelCostr = true;
      IncompatibleNPairConstr = true;
      TwiceEdgeConstr = true;
    }
  };

  std::vector<bool> Solution;
  std::vector<std::vector<bool>> BorderLabelN;

  void ExtractSurfacesFormSolution(
      const std::vector<Geo::Point3<ScalarType>> &OriginalVertPos,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      std::vector<std::vector<std::vector<int>>> &FoilFaces,
      std::vector<std::vector<int>> &OriginalFaceIdx,
      std::vector<std::vector<int>> &OriginalVertIdx,
      const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
      const std::vector<std::vector<int>> &ReprojectFaces,
      const std::vector<std::pair<int,int>> &FeatureE,
      const std::vector<int> &Corners,
      int SmoothSteps = 20) {

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

    SmoothReprojectSubMeshes(NonManifV,NonManifFaces,
                                  ReprojectV,ReprojectFaces,
                                  FeatureE,Corners,SmoothSteps);

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

  bool FindSurfaceSheet(const ExtractParam &ExtrP, std::vector<int> &ChosenF) {
   
    SingleF_Constr.clear();
    EdgeLabel_Constr.clear();
    IncompatibleNodePairs_Constr.clear();
    TwiceEdge_Constr.clear();

    std::cout << "Initializing Variable Mapping" << std::endl;
    InitVariableMap();

    std::cout << "Num Variables: " << VarToNode.size() << std::endl;

    if (ExtrP.SingleFConstr) {
      InitSingleF_Constr();
      std::cout << "Initialized " << SingleF_Constr.size()
                << " single face constraint" << std::endl;
    }

    if (ExtrP.SingleLabelCostr) {
      InitSingleLabelF_Constr();
      std::cout << "Initialized " << EdgeLabel_Constr.size()
                << " face label constraint" << std::endl;
    }

    if (ExtrP.IncompatibleNPairConstr) {
      InitIncompatibleNodePairs_Constr();
      std::cout << "Initialized " << IncompatibleNodePairs_Constr.size()
                << " Incompatible faces constraint" << std::endl;
    }

    if (ExtrP.TwiceEdgeConstr) {
      InitTwiceEdge_Constr();
      std::cout << "Initialized " << TwiceEdge_Constr.size()
                << " twice edge constraint" << std::endl;
    }

    bool solved = SolveGurobi(ExtrP.HardNodes, ExtrP.SoftNodes, ChosenF);

    return solved;
  }


  SheetExtractorILP(Geo::UndirectedNodeWeightedGraph<ScalarType> &_Graph,
                    std::vector<std::vector<int>> &_InternalFaces)
      : Graph(_Graph), InternalFaces(_InternalFaces) {}
};

#endif
