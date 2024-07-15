#ifndef PATH_SAMPLING
#define PATH_SAMPLING

#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>
#define VERT_MARK_SIZE 2

#include "../mesh_sampling.h"
#include "./connectivity.h"
#include "geodesic.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "path_functions.h"
#include "tracing_interface.h"
#include <list>
#include <set>
#include <unordered_set>

namespace Geo {

template <class ScalarType> class CrossPathSampling {

  std::vector<Geo::Point3<ScalarType>> &VertPos;
  std::vector<std::vector<int>> &Connectivity;
  Geo::TracingGraph<ScalarType> &TGraph;
  std::vector<Field::CrossF<ScalarType>> &VertCurv;

  std::vector<int> &SingIndex, &SingValue;
  std::vector<int> SingNodes;

  //    std::vector<size_t> BorderV;
  //    std::vector<int> BorderNodes;

  std::vector<int> StartingNodes;
  // std::vector<ScalarType> NodeDistances;

  enum SP_Status { SP_Selected, SP_Discarded, SP_Evaluate };

  bool updatePatches;
  // bool orthoBorderDist;

  //    std::unordered_set<int> ConflictNodes;
  //    //these are used to detect when two path shared an edge and to solve it
  //    std::unordered_set<std::pair<int,int>,Geo::hashPairIntFunction >
  //    ConflictEdges;

  std::set<int> ConflictNodes;

  // these are used to detect when two path shared an edge and to solve it
  std::set<std::pair<int, int>> ConflictEdges;

  std::vector<std::vector<int>> VF;
  std::vector<std::vector<int>> VFIndex;
  std::vector<std::vector<int>> NextF;
  std::vector<std::vector<int>> NextE;
  std::vector<bool> IsBorder;

public:
  bool only_visible;
  bool avoidSingularities;

  std::vector<std::vector<int>> SampledPaths;
  std::vector<SP_Status> SampledPaths_Status;
  std::vector<std::pair<ScalarType, int>> Candidates;

  int NumSamples;

  bool dynamicUpdates;
  std::vector<std::vector<int>> SelectedPaths;

  // this determine the nodes that can emit
  std::vector<bool> CanEmit;
  // this determine the nodes that can receive
  std::vector<bool> CanReceive;

  // //this used to get some subdomains to sample into subdomain or
  // //to check if a path pass from some vertices/faces
  // std::vector<int> VertTMark[VERT_MARK_SIZE];
  // int GlobalTMark[VERT_MARK_SIZE];

  // //used to see if a path pass on that vertices or not
  // bool check_pass_on[VERT_MARK_SIZE];

  // //sample on nodes that are within this set
  // bool initial_sample_on[VERT_MARK_SIZE];

  // //limit the path to pass only on those marked one
  // bool allow_path_on[VERT_MARK_SIZE];

  //  //limit the path to pass only on those marked one
  // bool allow_path_on[VERT_MARK_SIZE];

  // this used to get some subdomains to sample into subdomain or

  // //to check if a path pass from some vertices/faces
  // std::vector<int> VertTMark[VERT_MARK_SIZE];
  // int GlobalTMark;

  // used to see if a path pass on one of that vertices or not
  std::vector<bool> must_pass_on;

  // sample on vertices that are within this set
  std::vector<bool> sample_on;

  // constrain to add the sample from this set
  std::vector<bool> force_sample;

  // limit the path to pass only on those ones
  std::vector<bool> allow_path_on;

  // limit the path to arrive only on those ones
  std::vector<bool> allow_arrive_to;

  // the minimum distance to add samples im between
  int minSingResample;

  Geo::TracingGraph<ScalarType> &Graph() const { return TGraph; }

  std::vector<Field::CrossF<ScalarType>> CrossF() const { return VertCurv; }

  // void InitMarkVector()
  // {
  //     for (size_t i=0;i<VERT_MARK_SIZE;i++)
  //     {
  //         VertTMark[i]=std::vector<int>(VertPos.size(),0);
  //         GlobalTMark[i]=0;
  //     }
  // }

  // void InitEmitReceive()
  // {
  //     CanEmit=std::vector<bool>(TGraph.NumNodes(),true);
  //     CanReceive=std::vector<bool>(TGraph.NumNodes(),true);
  // }

  void FindExitAllowPathAreaNodes(const std::vector<int> &TestVertex,
                                  std::vector<int> &ExitNodes) {
    ExitNodes.clear();
    for (size_t i = 0; i < TestVertex.size(); i++) {

      // get the nodes of the vertex
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(TestVertex[i], IndexN);

      // if already has exit node it is ok
      bool has_exit = false;
      for (size_t j = 0; j < IndexN.size(); j++) {
        has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::ExitDomain);
        has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::TangentExitDomain);
        has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::ConvexExitDomain);
      }
      if (has_exit)
        continue;
      for (size_t j = 0; j < IndexN.size(); j++) {
        int currN = IndexN[j];
        if (ConflictNodes.count(currN) > 0)
          continue;
        if (TGraph.NodeType(currN) != Geo::Internal)
          continue;
        if (!IsTraceable(currN))
          continue;

        int NumNeghTest = NumAllowPathNeigh(currN);
        if (NumNeghTest > 0)
          continue;
        ExitNodes.push_back(currN);
      }
    }
  }

  void SampleCloseSingTracingVert(std::vector<int> &SampledVert, int stepDist) {
    SampledVert.clear();

    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(Connectivity, VertPos.size(), VVAdj);
    typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
    GParam.MaxSteps = stepDist;
    Geo::TMark TempMark;
    for (size_t i = 0; i < SingIndex.size(); i++) {
      int IndexV0 = SingIndex[i];
      GParam.Seeds.clear();
      GParam.Seeds.push_back(IndexV0);
      Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, VVAdj, TempMark, GParam);
      // std::cout<<"In Interval "<<GParam.InInterval.size()<<std::endl;
      // check wich other singularities are in the area
      std::vector<int> OtherSing;

      for (size_t j = 0; j < SingIndex.size(); j++) {
        int IndexV1 = SingIndex[j];
        if (IndexV1 <= IndexV0)
          continue;
        if (GParam.InInterval.count(IndexV1) == 0)
          continue;
        OtherSing.push_back(IndexV1);
      }

      for (size_t j = 0; j < OtherSing.size(); j++) {
        std::vector<int> VertSeq;
        Geo::Geodesic<ScalarType>::RetrievePath(GParam, OtherSing[j], VertSeq);
        if (VertSeq.size() <= 2)
          continue;
        int num = VertSeq.size();
        int MidV = num / 2; //(num % 2 == 0) ? (num / 2) : ((num / 2) + 1);
        SampledVert.push_back(VertSeq[MidV]);
      }
    }
  }

  // void SampleTracingVert(int NumSources,std::vector<int> &SampledVert,
  //                        const int IndexMark=-1)
  // {
  //     //ScalarType AreaC=VertPos.size()/NumSources
  //     int
  //     RadiusStar=Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(VertPos.size(),NumSources);
  //     //int RadiusStar=VertPos.size()/NumSources;
  //     SampledVert.clear();
  //     if (IndexMark!=-1)
  //     {
  //         std::vector<int> MarkedV;
  //         GetMarkedVert(IndexMark,MarkedV);
  //         //           std::cout<<" There are "<< MarkedV.size()<<"places to
  //         sample "<<std::endl;
  //         //RadiusStar=MarkedV.size()/NumSources;
  //         RadiusStar=Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(MarkedV.size(),NumSources);
  //         Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,
  //                                                                RadiusStar,4,SampledVert,&MarkedV);
  //     }
  //     else
  //         Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,
  //                                                                RadiusStar,4,SampledVert);
  // }

  void SampleTracingVert(int NumSources, std::vector<int> &SampledVert,
                         const int IndexMark = -1) {
    // ScalarType AreaC=VertPos.size()/NumSources
    std::vector<int> SampleableV;
    GetSamplableVert(SampleableV);
    int RadiusStar = Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(
        SampleableV.size(), NumSources);

    // int RadiusStar=VertPos.size()/NumSources;

    SampledVert.clear();
    Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(
        VertPos.size(), Connectivity, RadiusStar, 4, SampledVert, &SampleableV);

    std::vector<int> ForcedV;
    GetForcedSampleVert(ForcedV);
    SampledVert.insert(SampledVert.end(), ForcedV.begin(), ForcedV.end());
    // then add the forced one
    //  if (IndexMark!=-1)
    //  {
    //      std::vector<int> MarkedV;
    //      GetMarkedVert(IndexMark,MarkedV);
    //      //           std::cout<<" There are "<< MarkedV.size()<<"places to
    //      sample "<<std::endl;
    //      //RadiusStar=MarkedV.size()/NumSources;
    //      RadiusStar=Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(MarkedV.size(),NumSources);
    //      Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,
    //                                                             RadiusStar,4,SampledVert,&MarkedV);
    //  }
    //  else
    //      Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,
    //                                                             RadiusStar,4,SampledVert);
  }

  int NumAllowPathNeigh(int IndexN) {
    std::vector<int> NeighN;
    std::vector<ScalarType> NeighW;
    TGraph.GetNeigh(IndexN, false, true, NeighN, NeighW);
    int numMarkN = 0;
    for (size_t i = 0; i < NeighN.size(); i++) {
      int IndexNeighNode = NeighN[i];
      int IndexV = Field::CrossFieldGraph::ElemOfNode(IndexNeighNode);
      if (allow_path_on[IndexV])
        numMarkN++;
    }
    return numMarkN;
  }

  // int NumMarkedNeigh(int IndexN,const int IndexMark)
  // {
  //     assert(IndexMark>=0);
  //     assert(IndexMark<VERT_MARK_SIZE);

  //     std::vector<int> NeighN;
  //     std::vector<ScalarType> NeighW;
  //     TGraph.GetNeigh(IndexN,false,true,NeighN,NeighW);
  //     int numMarkN=0;
  //     for (size_t i=0;i<NeighN.size();i++)
  //         if (IsMarkedNode(NeighN[i],IndexMark))numMarkN++;

  //     return numMarkN;
  // }

  // void FindSampledNodesFromVertex(const std::vector<int> &SampledVert,
  //                                 std::vector<int> &SelectedNodes,
  //                                 int IndexMark=-1)
  // {
  //     SelectedNodes.clear();
  //     std::set<int> SingSetId;

  //     if (avoidSingularities)
  //         SingSetId.insert(SingIndex.begin(),SingIndex.end());

  //     //get nodes from sampled vertices
  //     for (size_t i=0;i<SampledVert.size();i++)
  //     {
  //         //not start tracing from singularities
  //         if (avoidSingularities)
  //             if (SingSetId.count(SampledVert[i])>0)continue;

  //         //then get possible tracing nodes
  //         std::vector<int> IndexN;
  //         Field::CrossFieldGraph::NodesOfElem(SampledVert[i],IndexN);

  //         //get only one per direction (the one with more neigh)
  //         int StartN[2]={-1,-1};
  //         int NumNeigh[2]={-1,-1};
  //         for (size_t j=0;j<IndexN.size();j++)
  //         {
  //             //if (AddedDir2[j%2]==true)continue;//already added in this or
  //             opposite dir if (!IsTraceable(IndexN[j]))continue;//if not
  //             traceable then skip if
  //             (ConflictNodes.count(IndexN[j])>0)continue;//if conflict that
  //             node than skip

  //             int NumNeghTest=TGraph.NeighSize(IndexN[j]);
  //             if (IndexMark!=-1)
  //                 NumNeghTest=(NumMarkedNeigh(IndexN[j],IndexMark));

  //             //no neighbour continue
  //             if (NumNeghTest==0)continue;

  //             int IndexM2=(j%2);
  //             if (NumNeigh[IndexM2]>NumNeghTest)continue;

  //             NumNeigh[IndexM2]=NumNeghTest;
  //             StartN[IndexM2]=IndexN[j];
  //         }
  //         if (StartN[0]!=-1)SelectedNodes.push_back(StartN[0]);
  //         if (StartN[1]!=-1)SelectedNodes.push_back(StartN[1]);
  //     }
  // }

  void FindSampledNodesFromVertex(const std::vector<int> &SampledVert,
                                  std::vector<int> &SelectedNodes) {
    SelectedNodes.clear();
    std::set<int> SingSetId;

    if (avoidSingularities)
      SingSetId.insert(SingIndex.begin(), SingIndex.end());

    // get nodes from sampled vertices
    for (size_t i = 0; i < SampledVert.size(); i++) {
      // not start tracing from singularities
      if (avoidSingularities)
        if (SingSetId.count(SampledVert[i]) > 0)
          continue;

      // then get possible tracing nodes
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(SampledVert[i], IndexN);

      // get only one per direction (the one with more neigh)
      int StartN[2] = {-1, -1};
      int NumNeigh[2] = {-1, -1};
      for (size_t j = 0; j < IndexN.size(); j++) {
        // if (AddedDir2[j%2]==true)continue;//already added in this or opposite
        // dir
        if (!IsTraceable(IndexN[j]))
          continue; // if not traceable then skip
        if (ConflictNodes.count(IndexN[j]) > 0)
          continue; // if conflict that node than skip

        int NumNeghTest = TGraph.NeighSize(IndexN[j]);
        // if (IndexMark!=-1)
        //     NumNeghTest=(NumMarkedNeigh(IndexN[j],IndexMark));

        // no neighbour continue
        if (NumNeghTest == 0)
          continue;

        int IndexM2 = (j % 2);
        if (NumNeigh[IndexM2] > NumNeghTest)
          continue;

        NumNeigh[IndexM2] = NumNeghTest;
        StartN[IndexM2] = IndexN[j];
      }
      if (StartN[0] != -1)
        SelectedNodes.push_back(StartN[0]);
      if (StartN[1] != -1)
        SelectedNodes.push_back(StartN[1]);
    }
  }

private:
  bool IsTraceable(int &IndexN) {
    // std::cout<<"test A"<<std::endl;
    // std::cout<<"Test:"<<IndexN<<std::endl;
    assert(IndexN >= 0);
    // std::cout<<"test B"<<std::endl;
    // assert(IndexN < TGraph.NumNodes());

    if (TGraph.NodeType(IndexN) == Geo::ExitDomain)
      return false;
    if (TGraph.NodeType(IndexN) == Geo::TangentExitDomain)
      return false;
    if (TGraph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    if (TGraph.NodeType(IndexN) == Geo::ConvexExitDomain)
      return false;

    return true;
  }

  void DisableSingularities() {
    for (size_t i = 0; i < SingIndex.size(); i++) {
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(SingIndex[i], IndexN);
      for (size_t j = 0; j < IndexN.size(); j++)
        TGraph.DisableNode(IndexN[j]);
    }
  }

  // void DisableNonMarked(int IndexMark)
  // {
  //     assert(IndexMark>=0);
  //     assert(IndexMark<VERT_MARK_SIZE);
  //     for (size_t i=0;i<VertPos.size();i++)
  //     {
  //         if (IsMarkedVert(i,IndexMark))continue;
  //         std::vector<int> IndexN;
  //         Field::CrossFieldGraph::NodesOfElem(i,IndexN);
  //         for (size_t j=0;j<IndexN.size();j++)
  //             TGraph.DisableNode(IndexN[j]);
  //     }
  // }

  void DisableNonAllowPath() {
    for (size_t i = 0; i < VertPos.size(); i++) {
      if (allow_path_on[i])
        continue;
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(i, IndexN);
      for (size_t j = 0; j < IndexN.size(); j++)
        TGraph.DisableNode(IndexN[j]);
    }
  }

  void EnableSingularities() {
    for (size_t i = 0; i < SingIndex.size(); i++) {
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(SingIndex[i], IndexN);
      for (size_t j = 0; j < IndexN.size(); j++)
        TGraph.EnableNode(IndexN[j]);
    }
  }

  //    void InitializeDistance()
  //    {
  //        TGraph.CurrDistIndex()=0;

  //        std::vector<int> StartN;

  //        //first add singularities
  //        if (avoidSingularities)
  //            StartN.insert(StartN.end(),SingNodes.begin(),SingNodes.end());

  //        //add border nodes
  //        StartN.insert(StartN.end(),BorderNodes.begin(),BorderNodes.end());

  //        //then add all the previously added nodes
  //        for (size_t i=0;i<SelectedPaths.size();i++)
  //        {
  //            //add the path and opposite
  //            StartN.insert(StartN.end(),SelectedPaths[i].begin(),SelectedPaths[i].end());
  //            std::vector<int> OppNodes;
  //            Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths[i],OppNodes);
  //            StartN.insert(StartN.end(),OppNodes.begin(),OppNodes.end());
  //        }

  //        if (StartN.size()==0)
  //        {
  //            //add the ones with smaller X
  //            int IndexV=-1;
  //            ScalarType minX=std::numeric_limits<ScalarType>::max();
  //            for (size_t i=0;i<VertPos.size();i++)
  //            {
  //                if (VertPos[i].X>minX)continue;
  //                IndexV=i;
  //                minX=VertPos[i].X;
  //            }
  //            assert(IndexV>=0);
  //            StartN.push_back(IndexV);
  //        }

  //        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
  //        SParam.avoid_zero_path=true;
  //        for (size_t i=0;i<SELECTED_GRAPH_FLAG;i++)
  //            SParam.StopAtSel[i]=false;
  //        SParam.MaxL=std::numeric_limits<ScalarType>::max();
  //        SParam.onlyDirect=false;
  //        SParam.onlyVisible=false;
  //        SParam.stopAtPrevDist=false;
  //        SParam.stopAtPrevDist=false;

  //        int StopNode;
  //        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,StartN,SParam,StopNode);
  //    }

  void InitializeDistance() {
    TGraph.CurrDistIndex() = 0;
    // TGraph.SetMaxDistance();

    std::vector<int> StartN;

    // first add singularities
    if (avoidSingularities)
      StartN.insert(StartN.end(), SingNodes.begin(), SingNodes.end());

    //        GetNodeType(const NType &TypeNode,
    //                    std::vector<int> &IdxNode)

    //        //get also the orthogonal to border or not to compute distance
    //        if (!orthoBorderDist)
    //        {
    std::vector<int> BorderNodes;
    TGraph.GetNodeType(Geo::TangentDomain, BorderNodes);
    // add border nodes
    StartN.insert(StartN.end(), BorderNodes.begin(), BorderNodes.end());
    //        }else
    //        {
    //            std::vector<size_t> BorderV;
    //            std::vector<std::vector<int> > NextF,NextE;
    //            getFFAdjacency(Connectivity,NextF,NextE);
    //            getBorderVert(Connectivity,NextF,BorderV);
    //            for (size_t i=0;i<BorderV.size();i++)
    //            {
    //                std::vector<int> IndexN;
    //                Field::CrossFieldGraph::NodesOfElem(BorderV[i],IndexN);
    //                StartN.insert(StartN.end(),IndexN.begin(),IndexN.end());
    //            }
    //        }

    // then add all the previously added nodes
    for (size_t i = 0; i < SelectedPaths.size(); i++) {
      // add the path and opposite
      StartN.insert(StartN.end(), SelectedPaths[i].begin(),
                    SelectedPaths[i].end());
      std::vector<int> OppNodes;
      Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths[i],
                                                     OppNodes);
      StartN.insert(StartN.end(), OppNodes.begin(), OppNodes.end());
    }

    if (StartN.size() == 0) {
      // add the ones with smaller X
      int IndexV = -1;
      ScalarType minX = std::numeric_limits<ScalarType>::max();
      for (size_t i = 0; i < VertPos.size(); i++) {
        if (VertPos[i].X > minX)
          continue;
        IndexV = i;
        minX = VertPos[i].X;
      }
      assert(IndexV >= 0);
      StartN.push_back(IndexV);
    }

    typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    SParam.avoid_zero_path = true;
    for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
      SParam.StopAtSel[i] = false;
    SParam.MaxL = std::numeric_limits<ScalarType>::max();
    SParam.onlyDirect = false;
    SParam.onlyVisible = true;
    SParam.stopAtPrevDist = false;

    int StopNode;
    Geo::GraphQuery<ScalarType>::Dijkstra(TGraph, StartN, SParam, StopNode);

    //        //then check of some has not been reached
    //        for (size_t i=0;i<TGraph.NumNodes();i++)
    //        {
    //            if (TGraph.IsMarked(i))continue;
    //            std::vector<int> IndexN;
    //            Field::CrossFieldGraph::NodesOfElem(i,IndexN);
    //            ScalarType AvgDist=0;
    //            int Num=0;
    //            for (size_t j=0;j<IndexN.size();j++)
    //            {
    //                if (TGraph.IsMarked(IndexN[j]))continue;
    //                AvgDist+=TGraph.Distance(IndexN[j]);
    //                Num++;
    //            }
    //            if (Num>0)
    //            {
    //                AvgDist/=Num;
    //                TGraph.Distance(i)=AvgDist;
    //            }
    //        }
    //        //then get the maximum value and set to it the untouched
    ScalarType MaxVal = TGraph.MaxDistance();
    TGraph.SetUnMarkedToDistance(MaxVal);
  }

  ScalarType Distance(const std::vector<int> &TestPath) {
    ScalarType CurrD = 0;
    for (size_t i = 0; i < TestPath.size(); i++) {
      int N0 = TestPath[i];
      int N1 = Field::CrossFieldGraph::OppositeNodeCrossField(N0);
      assert(N0 >= 0);
      assert(N1 >= 0);
      assert(N0 < (int)TGraph.NumNodes());
      assert(N1 < (int)TGraph.NumNodes());
      CurrD += TGraph.Distance(N0);
      CurrD += TGraph.Distance(N1);
    }
    CurrD /= (ScalarType)TestPath.size() * 2;
    return CurrD;
  }

  void GetArrivalNodes(std::vector<int> &ArrivalNodes) {
    ArrivalNodes.clear();
    for (size_t i = 0; i < allow_arrive_to.size(); i++) {
      if (!allow_arrive_to[i])
        continue;
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(i, IndexN);
      for (size_t j = 0; j < IndexN.size(); j++) {
        int currN = IndexN[j];
        bool isExit = ((TGraph.NodeType(currN) == Geo::ExitDomain) ||
                       (TGraph.NodeType(currN) == Geo::TangentExitDomain));
        if (isExit)
          ArrivalNodes.push_back(currN);
      }
    }
  }

  void SamplePaths(int NumSamples) {
    StartingNodes.clear();
    SampledPaths.clear();

    std::cout << "Sampling starting Nodes" << std::endl;
    // int MarkSample=-1;
    // for (size_t i=0;i<VERT_MARK_SIZE;i++)
    //     if (sample_only_on[i])MarkSample=i;

    // if (MarkSample!=-1)
    //     std::cout<<"Curr Mark Sampling:"<<MarkSample<<std::endl;

    // InitTracingNodes(NumSamples,MarkSample);
    std::vector<int> SampledVert;
    SampleTracingVert(NumSamples, SampledVert);
    // SampleTracingVert(NumSamples,SampledVert,MarkSample);
    // std::cout<<"Sampled "<< SampledVert.size() << " starting Vert"
    // <<std::endl;

    // add in between singularity samples
    if (minSingResample > 0) {
      std::vector<int> SingSingSampledVert;
      SampleCloseSingTracingVert(SingSingSampledVert, minSingResample);
      for (size_t i = 0; i < SingSingSampledVert.size(); i++) {
        if (!sample_on[SingSingSampledVert[i]])
          continue;
        SampledVert.push_back(SingSingSampledVert[i]);
      }
      // SampledVert.insert(SampledVert.end(), SingSingSampledVert.begin(),
      // SingSingSampledVert.end());
    }

    // FindSampledNodesFromVertex(SampledVert,StartingNodes,MarkSample);
    FindSampledNodesFromVertex(SampledVert, StartingNodes);

    std::vector<int> NewStartingNodes;

    int non_traced = 0;
    int self_intersected = 0;

    // DisableNonAllowPath();

    std::vector<int> ArrivalNodes;
    GetArrivalNodes(ArrivalNodes);

    // FindArrivalNodesFromVertex(SampledVert,StartingNodes);

    for (size_t i = 0; i < StartingNodes.size(); i++) {
      assert(TGraph.EnabledN[StartingNodes[i]]);
      std::vector<int> CurrPath;
      assert(TGraph.NodeType(StartingNodes[i]) != ExitDomain);
      assert(TGraph.NodeType(StartingNodes[i]) != TangentDomain);
      assert(TGraph.IsEnabled(StartingNodes[i]));
      std::chrono::time_point<std::chrono::high_resolution_clock> start =
          std::chrono::high_resolution_clock::now();
      bool has_traced =
          Geo::GraphTracing<ScalarType>::TraceLoopOrBorderBorderFrom(
              TGraph, StartingNodes[i], CurrPath, true, &ArrivalNodes,
              only_visible);
      // std::cout<<"i:"<<i<<std::endl;
      // assert(TGraph.NodeType(CurrPath[0])!=ExitDomain);
      // assert(TGraph.NodeType(CurrPath.back())==ExitDomain);
      // assert(0);
      std::chrono::time_point<std::chrono::high_resolution_clock> end =
          std::chrono::high_resolution_clock::now();
      std::chrono::duration duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      // TimeTracing += duration.count() * 1e-6;

      if (!has_traced) {
        non_traced++;
        continue;
      }

      start = std::chrono::high_resolution_clock::now();

      bool self_intersect =
          Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,
                                                                    CurrPath);

      end = std::chrono::high_resolution_clock::now();
      duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start);

      // TimeConflict += duration.count() * 1e-6;

      if (self_intersect) {
        self_intersected++;
        continue;
      }
      SampledPaths.push_back(CurrPath);
      NewStartingNodes.push_back(StartingNodes[i]);
    }

    // //then restore
    // if (IndexLimit!=-1)
    // {
    //     for (size_t i=0;i<SwapNodeType.size();i++)
    //     {
    //         int IndexN=SwapNodeType[i].first;
    //         NType TypeNode=SwapNodeType[i].second;
    //         TGraph.NodeType(IndexN)=TypeNode;
    //     }
    // }

    StartingNodes = NewStartingNodes;
    SampledPaths_Status.clear();
    SampledPaths_Status.resize(SampledPaths.size(), SP_Evaluate);

    std::cout << "Sampled " << SampledPaths.size() << " Paths" << std::endl;
    std::cout << "self intersected " << self_intersected << std::endl;
    std::cout << "Non traced " << non_traced << std::endl;
  }

  void InitCandidates() {
    Candidates.clear();
    for (size_t i = 0; i < SampledPaths.size(); i++)
      Candidates.push_back(std::pair<ScalarType, int>(-1, i));
  }

  void UpdateCandidatePriority() {
    for (size_t i = 0; i < Candidates.size(); i++) {
      int IndexPath = Candidates[i].second;
      assert(IndexPath >= 0);
      assert(IndexPath < SampledPaths.size());
      ScalarType Dist = Distance(SampledPaths[IndexPath]);
      Candidates[i].first = Dist;
    }
    std::sort(Candidates.begin(), Candidates.end());
  }

  bool ConflictWithInserted(const std::vector<int> &TestPath) {

    for (size_t i = 0; i < TestPath.size(); i++) {
      int test_node = TestPath[i];
      if (ConflictNodes.count(test_node) > 0)
        return true;
    }
    return false;
  }

  bool CrossProperly(const std::vector<int> &TestPath) {
    // test if cross the the vertex might touch and stay on the
    // same side, if touch has to cross
    bool IsLoop = Geo::PathFunctions<ScalarType>::IsLoop(TestPath);
    int StartN = 1;
    if (IsLoop)
      StartN = 0;

    // if loops the last and the first are the same so do not consider
    if (TestPath.size() < 3)
      return true;
    for (size_t i = StartN; i < TestPath.size() - 1; i++) {
      int indexCurr = i;
      int indexPred = indexCurr - 1;
      int indexSucc = indexCurr + 1;

      // if was on zero then pred <0 then take the last
      if (indexPred < 0)
        indexPred = (TestPath.size() - 2);

      assert(indexCurr != indexPred);
      assert(indexCurr != indexSucc);
      assert(indexSucc != indexPred);

      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(TestPath[indexPred]);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(TestPath[indexCurr]);
      int IndexV2 = Field::CrossFieldGraph::ElemOfNode(TestPath[indexSucc]);

      assert(IndexV0 != IndexV1);
      assert(IndexV1 != IndexV2);
      assert(IndexV0 != IndexV2);

      assert(IndexV1 >= 0);
      assert(IndexV1 < IsBorder.size());
      // if border no problem
      if (IsBorder[IndexV1])
        continue;

      std::vector<int> NodesV1;
      Field::CrossFieldGraph::NodesOfElem(IndexV1, NodesV1);
      // check if intersect
      bool foundInt = false;
      for (size_t j = 0; j < NodesV1.size(); j++)
        if (ConflictNodes.count(NodesV1[j]) > 0) {
          foundInt = true;
          break;
        }
      if (!foundInt)
        continue;

      std::vector<int> VertFan;
      bool HasBorder;

      // then see if it cross it
      GetSortedVertStar(IndexV1, Connectivity, VF, VFIndex, NextF, NextE,
                        VertFan, HasBorder);
      assert(!HasBorder);

      std::vector<bool> NewPath;
      for (size_t j = 0; j < VertFan.size(); j++) {
        assert(VertFan[j] != IndexV1);
        int Ve0 = IndexV1;
        int Ve1 = VertFan[j];
        std::pair<int, int> key(std::min(Ve0, Ve1), std::max(Ve0, Ve1));
        bool IsInsertedE = (ConflictEdges.count(key) > 0);
        bool IsNewE = ((VertFan[j] == IndexV0) || (VertFan[j] == IndexV2));
        if (IsInsertedE)
          assert(!IsNewE);

        if (IsNewE)
          assert(!IsInsertedE);

        if (IsInsertedE)
          NewPath.push_back(false);

        if (IsNewE)
          NewPath.push_back(true);
      }
      //   if ((NewPath.size() != 4)) {
      //     std::cout << "WARINING Size:" << NewPath.size() << std::endl;
      //     for (size_t k = 0; k < NewPath.size(); k++)
      //       std::cout << NewPath[k] << " ";
      //     std::cout << std::endl;
      //   }

      assert(NewPath.size() == 4);
      for (size_t k = 0; k < NewPath.size(); k++)
        if (NewPath[k] == NewPath[(k + 1) % NewPath.size()])
          return false;
    }
    return true;
  }

  bool ConflictEdgeWithInsertedEdges(const std::vector<int> &TestPath) {
    // test if there is the same edge
    for (size_t i = 0; i < TestPath.size() - 1; i++) {
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(TestPath[i]);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(TestPath[i + 1]);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      if (ConflictEdges.count(key) > 0)
        return true;
    }

    return false;
  }

  //    void InitBorderNodes()
  //    {
  //        std::vector<std::vector<int> > NextF,NextE;
  //        Geo::getFFAdjacency(Connectivity,NextF,NextE);

  //        //then add the borders
  //        Geo::getBorderVert(Connectivity,NextF,BorderV);

  //        BorderNodes.clear();
  //        for (size_t j=0;j<BorderV.size();j++)
  //        {
  //            std::vector<int> IndexN;
  //            Field::CrossFieldGraph::NodesOfElem(BorderV[j],IndexN);
  //            for (size_t j=0;j<IndexN.size();j++)
  //            {
  //                if (!IsTraceable(IndexN[j]))continue;
  //                BorderNodes.push_back(IndexN[j]);
  //            }
  //        }
  //    }

  void AddToConflictNodes(const std::vector<int> &NodesPath) {
    std::vector<int> OppNodesPath;
    Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
    ConflictNodes.insert(NodesPath.begin(), NodesPath.end());
    ConflictNodes.insert(OppNodesPath.begin(), OppNodesPath.end());
  }

  void AddToConflictEdges(const std::vector<int> &NodesPath) {
    for (size_t i = 0; i < NodesPath.size() - 1; i++) {
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i]);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i + 1]);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      ConflictEdges.insert(key);
    }
  }

public:
  void EnableGraphNodes(const std::vector<int> &NodesPath) {
    std::vector<int> To_Enable = NodesPath;
    std::vector<int> OppNodesPath;
    Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
    To_Enable.insert(To_Enable.end(), OppNodesPath.begin(), OppNodesPath.end());
    TGraph.EnableNodes(To_Enable);
  }

  void DisableGraphNodes(const std::vector<int> &NodesPath) {
    std::vector<int> To_Disable = NodesPath;
    std::vector<int> OppNodesPath;
    Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
    To_Disable.insert(To_Disable.end(), OppNodesPath.begin(),
                      OppNodesPath.end());
    TGraph.DisableNodes(To_Disable);
  }

  void EnableGraphArches(const std::vector<int> &NodesPath) {
    Geo::PathFunctions<ScalarType>::CrossEnableArcsOfPath(NodesPath, TGraph,
                                                          true);
  }

  void DisableGraphArches(const std::vector<int> &NodesPath) {
    Geo::PathFunctions<ScalarType>::CrossEnableArcsOfPath(NodesPath, TGraph,
                                                          false);
  }

private:
  void UpdateDistancesWithLastInsertedPath() {
    std::vector<int> StartN;
    StartN.insert(StartN.end(), SelectedPaths.back().begin(),
                  SelectedPaths.back().end());
    std::vector<int> OppNodes;
    Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths.back(),
                                                   OppNodes);
    StartN.insert(StartN.end(), OppNodes.begin(), OppNodes.end());

    typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    SParam.avoid_zero_path = true;
    for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
      SParam.StopAtSel[i] = false;
    SParam.MaxL = std::numeric_limits<ScalarType>::max();
    SParam.onlyDirect = false;
    SParam.onlyVisible = false;
    SParam.stopAtPrevDist = true;

    int StopNode;
    Geo::GraphQuery<ScalarType>::Dijkstra(TGraph, StartN, SParam, StopNode);
  }

  void RetraceCandidatesConflictWithSelected() {
    assert(SampledPaths_Status.size() == SampledPaths.size());
    assert(StartingNodes.size() == SampledPaths.size());

    std::vector<int> ArrivalNodes;
    GetArrivalNodes(ArrivalNodes);

    for (size_t i = 0; i < SampledPaths.size(); i++) {
      if (SampledPaths_Status[i] != SP_Evaluate)
        continue;

      bool has_conflict_node = ConflictWithInserted(SampledPaths[i]);
      bool has_conflict_edge = ConflictEdgeWithInsertedEdges(SampledPaths[i]);
      if ((!has_conflict_node) && (!has_conflict_edge))
        continue;

      // check if the start nodes is avalable
      if (ConflictNodes.count(StartingNodes[i]) > 0) {
        assert(!TGraph.IsEnabled(StartingNodes[i]));
        SampledPaths_Status[i] = SP_Discarded;
        continue;
      }

      bool has_traced =
          Geo::GraphTracing<ScalarType>::TraceLoopOrBorderBorderFrom(
              TGraph, StartingNodes[i], SampledPaths[i], true, &ArrivalNodes,
              only_visible);
      if (!has_traced) {
        SampledPaths_Status[i] = SP_Discarded;
        continue;
      }
      bool self_intersect =
          Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(
              TGraph, SampledPaths[i]);
      self_intersect |=
          Geo::PathFunctions<ScalarType>::HasSelfEdgeTwice(SampledPaths[i]);
      if (self_intersect) {
        SampledPaths_Status[i] = SP_Discarded;
        continue;
      }
    }
  }

  void WriteTiming() {
    // std::cout<<"*** TIMINGS ***"<<std::endl;
    // std::cout<<"* Time Init:"<<TimeInitData<<std::endl;
    // std::cout<<"* Time tracing:"<<TimeTracing<<std::endl;
    // std::cout<<"* Time upd Dist:"<<TimeUpdDist<<std::endl;
    // std::cout<<"* Time Conflicts:"<<TimeConflict<<std::endl;
  }

  // bool IsMarkedVert(const int &IndexV,const int &IndexMark)
  // {
  //     assert(IndexMark>=0);
  //     assert(IndexMark<VERT_MARK_SIZE);
  //     assert(IndexV>=0);
  //     assert(IndexV<VertPos.size());
  //     return (VertTMark[IndexMark][IndexV]==GlobalTMark[IndexMark]);
  // }

  // bool IsMarkedNode(const int &IndexNode,int IndexMark)
  // {
  //     assert(IndexMark>=0);
  //     assert(IndexMark<VERT_MARK_SIZE);
  //     assert(IndexNode>=0);
  //     assert(IndexNode<TGraph.NumNodes());
  //     int IndexV=Field::CrossFieldGraph::ElemOfNode(IndexNode);
  //     return (IsMarkedVert(IndexV,IndexMark));
  // }

public:
  void CleanEmptySelectedPaths() {
    std::map<int, int> NewPathIdx;
    std::vector<std::vector<int>> SelectedPathsSwap;
    for (size_t i = 0; i < SelectedPaths.size(); i++) {
      if (SelectedPaths[i].size() == 0)
        continue;
      SelectedPathsSwap.push_back(SelectedPaths[i]);
      NewPathIdx[i] = SelectedPathsSwap.size() - 1;
    }
    SelectedPaths = SelectedPathsSwap;

    std::vector<std::pair<ScalarType, int>> CandidatesSwap;
    for (size_t i = 0; i < Candidates.size(); i++) {
      int IndexPath = Candidates[i].second;
      if (NewPathIdx.count(IndexPath) == 0)
        continue;

      ScalarType Priority = Candidates[i].first;
      IndexPath = NewPathIdx[IndexPath];

      CandidatesSwap.push_back(std::pair<ScalarType, int>(Priority, IndexPath));
    }
    Candidates = CandidatesSwap;
  }

  void InitSingularityNodes(bool use_star = true) {
    std::vector<int> VertToUse = SingIndex;
    // add sing star in case
    if (use_star) {
      std::vector<std::vector<int>> Kernel;
      Geo::FindVertexStar(Connectivity, VertPos.size(), Kernel);
      for (size_t j = 0; j < SingIndex.size(); j++) {
        int IndexV = SingIndex[j];
        assert(IndexV >= 0);
        assert(IndexV < Kernel.size());
        VertToUse.insert(VertToUse.end(), Kernel[IndexV].begin(),
                         Kernel[IndexV].end());
      }
    }
    // move one to
    SingNodes.clear();
    for (size_t j = 0; j < VertToUse.size(); j++) {

      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(VertToUse[j], IndexN);

      for (size_t k = 0; k < IndexN.size(); k++) {
        if (!IsTraceable(IndexN[k]))
          continue;

        SingNodes.push_back(IndexN[k]);

        // if (!IsTraceable(VertToUse[j]))
        //     continue;

        // SingNodes.push_back(VertToUse[j]);
      }
    }
  }

  void GetSamplableVert(std::vector<int> &SamplVert) {
    SamplVert.clear();
    assert(sample_on.size() == VertPos.size());
    for (size_t i = 0; i < sample_on.size(); i++) {
      if (!sample_on[i])
        continue;
      SamplVert.push_back(i);
    }
  }

  void GetForcedSampleVert(std::vector<int> &ForcedSamplVert) {
    ForcedSamplVert.clear();
    assert(force_sample.size() == VertPos.size());
    for (size_t i = 0; i < force_sample.size(); i++) {
      if (!force_sample[i])
        continue;
      ForcedSamplVert.push_back(i);
    }
  }

  // void GetMarkedVert(const int &IndexMark,std::vector<int> &MarkedV)
  // {
  //     assert(IndexMark>=0);
  //     assert(IndexMark<VERT_MARK_SIZE);
  //     MarkedV.clear();
  //     for (size_t i=0;i<VertPos.size();i++)
  //     {
  //         if (!IsMarkedVert(i,IndexMark))continue;
  //         MarkedV.push_back(i);
  //     }
  // }

  // void UnMarkAll(int IndexMark)
  // {
  //     assert(IndexMark>=0);
  //     assert(IndexMark<VERT_MARK_SIZE);
  //     GlobalTMark[IndexMark]++;
  // }

  // void MarkVertex(const int &IndexV,int IndexMark)
  // {

  //     assert(IndexV>=0);
  //     assert(IndexV<VertTMark[IndexMark].size());
  //     assert(IndexV<(int)VertPos.size());

  //     VertTMark[IndexMark][IndexV]=GlobalTMark[IndexMark];
  // }

  // void MarkVertices(const std::vector<int> &To_MarkVert,int IndexMark)
  // {
  //     assert(IndexMark>=0);
  //     assert(IndexMark<VERT_MARK_SIZE);
  //     for (size_t i=0;i<To_MarkVert.size();i++)
  //     {
  //         int IndexV=To_MarkVert[i];
  //         // assert(IndexV>=0);
  //         // assert(IndexV<VertTMark[IndexMark].size());
  //         // assert(IndexV<(int)VertPos.size());
  //         // VertTMark[IndexMark][IndexV]=GlobalTMark[IndexMark];
  //         MarkVertex(IndexV,IndexMark);
  //     }
  // }

  // bool HasMarkedVert(int &IndexCandidate,int IndexMark)
  // {
  //     assert(IndexCandidate<SampledPaths.size());
  //     for (size_t i=0;i<SampledPaths[IndexCandidate].size();i++)
  //         if (IsMarkedNode(SampledPaths[IndexCandidate][i],IndexMark))return
  //         true;
  //     return false;
  // }

  void ClearSampledPaths() {
    StartingNodes.clear();
    SampledPaths.clear();
    SampledPaths_Status.clear();
    Candidates.clear();
    // ConflictNodes.clear();
    // SelectedPaths.clear();
  }

  void WriteStats() {
    std::cout << "Added Total:" << SelectedPaths.size() << " paths"
              << std::endl;
    int NumNodes = 0;
    for (int i = 0; i < SelectedPaths.size(); i++)
      NumNodes += SelectedPaths[i].size();
    std::cout << "Nodes Total:" << NumNodes << std::endl;
  }

  // bool AddNextPath()
  // {
  //     if (Candidates.size()==0)return false;

  //     TGraph.CurrDistIndex()=0;
  //     //start from the last one in reverse order (add the farthest)
  //     for (int i=Candidates.size()-1;i>=0;i--)
  //     {
  //         int IndexPath=Candidates[i].second;

  //         if (SampledPaths_Status[IndexPath]!=SP_Evaluate)continue;

  //         //CHECK OF PASSING ON REQUIRED PATCHES
  //         bool pass_on_patches=true;
  //         for (size_t j=0;j<VERT_MARK_SIZE;j++)
  //             if (check_pass_on[j])pass_on_patches=false;

  //         if (!pass_on_patches)
  //         {
  //             //std::cout<<"checking useful path"<<std::endl;
  //             for (size_t j=0;j<VERT_MARK_SIZE;j++)
  //             {
  //                 if (check_pass_on[j] && (HasMarkedVert(IndexPath,j)))
  //                 {
  //                     pass_on_patches=true;
  //                     break;
  //                 }
  //             }
  //         }

  //         if (!pass_on_patches)
  //         {
  //             SampledPaths_Status[IndexPath]=SP_Discarded;
  //             continue;
  //         }

  //         //CHECK OF INTERSECTION
  //         bool has_conflict=false;

  //         auto start = std::chrono::high_resolution_clock::now();

  //         if (!dynamicUpdates)
  //             has_conflict=ConflictWithInserted(SampledPaths[IndexPath]);

  //         if (!has_conflict)
  //             has_conflict|=ConflictEdgeWithInsertedEdges(SampledPaths[IndexPath]);

  //         if (has_conflict)
  //         {
  //             SampledPaths_Status[IndexPath]=SP_Discarded;
  //             continue;
  //         }

  //         //add the new path
  //         SampledPaths_Status[IndexPath]=SP_Selected;
  //         SelectedPaths.push_back(SampledPaths[IndexPath]);

  //         AddToConflictNodes(SelectedPaths.back());
  //         AddToConflictEdges(SelectedPaths.back());

  //         auto end = std::chrono::high_resolution_clock::now();
  //         auto duration =
  //         std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  //         //TimeConflict += duration.count() * 1e-6;

  //         //UPDATE DISTANCES WITH LAST ADDED
  //         if (avoidSingularities)
  //             EnableSingularities();

  //         start = std::chrono::high_resolution_clock::now();
  //         UpdateDistancesWithLastInsertedPath();

  //         if (avoidSingularities)
  //             DisableSingularities();

  //         end = std::chrono::high_resolution_clock::now();
  //         duration =
  //         std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  //         //TimeUpdDist += duration.count() * 1e-6;

  //         //DisableGraphNodes(SelectedPaths.back());
  //         //UPDATE CANDIDATES IF NEEDED
  //         if (dynamicUpdates)
  //         {
  //             DisableGraphNodes(SelectedPaths.back());

  //             TGraph.CurrDistIndex()=1;
  //             start = std::chrono::high_resolution_clock::now();

  //             RetraceCandidatesConflictWithSelected();

  //             end = std::chrono::high_resolution_clock::now();
  //             duration =
  //             std::chrono::duration_cast<std::chrono::microseconds>(end -
  //             start);
  //             //TimeTracing += duration.count() * 1e-6;

  //             TGraph.CurrDistIndex()=0;
  //         }

  //         //UPDATE PRIORITY
  //         UpdateCandidatePriority();

  //         return true;
  //     }
  //     return false;
  // }

  bool PassOnMustVertex(int &IndexCandidate) {
    assert(IndexCandidate < SampledPaths.size());
    for (size_t i = 0; i < SampledPaths[IndexCandidate].size(); i++) {
      int IndexN = SampledPaths[IndexCandidate][i];
      int IndexV = Field::CrossFieldGraph::ElemOfNode(IndexN);
      assert(IndexV < must_pass_on.size());
      if (must_pass_on[IndexV])
        return true;
    }
    return false;
  }

  bool AddNextPath() {
    if (Candidates.size() == 0)
      return false;

    TGraph.CurrDistIndex() = 0;
    // start from the last one in reverse order (add the farthest)
    for (int i = Candidates.size() - 1; i >= 0; i--) {
      int IndexPath = Candidates[i].second;

      if (SampledPaths_Status[IndexPath] != SP_Evaluate)
        continue;

      // //CHECK OF PASSING ON REQUIRED PATCHES
      // bool pass_on_patches=true;
      // for (size_t j=0;j<VERT_MARK_SIZE;j++)
      //     if (check_pass_on[j])pass_on_patches=false;

      // if (!pass_on_patches)
      // {
      //     //std::cout<<"checking useful path"<<std::endl;
      //     for (size_t j=0;j<VERT_MARK_SIZE;j++)
      //     {
      //         if (check_pass_on[j] && (HasMarkedVert(IndexPath,j)))
      //         {
      //             pass_on_patches=true;
      //             break;
      //         }
      //     }
      // }

      // CHECK OF PASSING ON REQUIRED PATCHES
      bool pass_on_must_vertices = PassOnMustVertex(IndexPath);

      if (!pass_on_must_vertices) {
        SampledPaths_Status[IndexPath] = SP_Discarded;
        continue;
      }

      // CHECK OF INTERSECTION
      bool has_conflict = false;

      auto start = std::chrono::high_resolution_clock::now();

      if (!dynamicUpdates)
        has_conflict = ConflictWithInserted(SampledPaths[IndexPath]);

      if (!has_conflict)
        has_conflict |= ConflictEdgeWithInsertedEdges(SampledPaths[IndexPath]);

      if (has_conflict) {
        SampledPaths_Status[IndexPath] = SP_Discarded;
        continue;
      }

      bool cross_properly = CrossProperly(SampledPaths[IndexPath]);
      if (!cross_properly) {
        // std::cout << "Improper CROSS" << std::endl;
        SampledPaths_Status[IndexPath] = SP_Discarded;
        continue;
      }

      // add the new path
      SampledPaths_Status[IndexPath] = SP_Selected;
      SelectedPaths.push_back(SampledPaths[IndexPath]);

      AddToConflictNodes(SelectedPaths.back());
      AddToConflictEdges(SelectedPaths.back());

      auto end = std::chrono::high_resolution_clock::now();
      auto duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      // TimeConflict += duration.count() * 1e-6;

      // UPDATE DISTANCES WITH LAST ADDED
      if (avoidSingularities)
        EnableSingularities();

      start = std::chrono::high_resolution_clock::now();
      UpdateDistancesWithLastInsertedPath();

      if (avoidSingularities)
        DisableSingularities();

      end = std::chrono::high_resolution_clock::now();
      duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      // TimeUpdDist += duration.count() * 1e-6;

      // DisableGraphNodes(SelectedPaths.back());
      // UPDATE CANDIDATES IF NEEDED
      if (dynamicUpdates) {
        DisableGraphNodes(SelectedPaths.back());
        DisableGraphArches(SelectedPaths.back());
        TGraph.CurrDistIndex() = 1;
        start = std::chrono::high_resolution_clock::now();

        RetraceCandidatesConflictWithSelected();

        end = std::chrono::high_resolution_clock::now();
        duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        // TimeTracing += duration.count() * 1e-6;

        TGraph.CurrDistIndex() = 0;
      }

      // UPDATE PRIORITY
      UpdateCandidatePriority();

      return true;
    }
    return false;
  }

  // void GetMarkedVertices(std::vector<int> &MarkedV,
  //                        int IndexMark)
  // {
  //     assert(IndexMark>=0);
  //     assert(IndexMark<VERT_MARK_SIZE);
  //     MarkedV.clear();
  //     for (size_t i=0;i<VertPos.size();i++)
  //         if (IsMarkedVert(i,IndexMark))MarkedV.push_back(i);
  // }

  void InsertAllPaths() {
    while (AddNextPath()) {
    };
    WriteStats();
    WriteTiming();
  }

  //    void SampleCandidates(bool Init_Distances=true)
  //    {
  //        ClearSampledPaths();
  //
  //        if (avoidSingularities)
  //            DisableSingularities();
  //        SamplePaths(NumSamples);
  //
  //        if (avoidSingularities)
  //            EnableSingularities();
  //
  //        if (Init_Distances)
  //            InitializeDistance();
  //
  //        InitCandidates();
  //
  //        UpdateCandidatePriority();

  //    }

  std::vector<std::pair<int, NType>> SwapNodeType;

  void SetBorderCanPassToExit(const std::vector<int> &TestVertex) {
    SwapNodeType.clear();
    // set ones possible to exit
    std::vector<int> IndexN;
    FindExitAllowPathAreaNodes(TestVertex, IndexN);
    // set as exit node
    for (size_t i = 0; i < IndexN.size(); i++) {
      int currN = IndexN[i];
      std::pair<int, NType> Entry(currN, TGraph.NodeType(currN));
      SwapNodeType.push_back(Entry);
      TGraph.NodeType(currN) = Geo::ExitDomain;
    }
  }

  void RestoreNodeType() {
    for (size_t i = 0; i < SwapNodeType.size(); i++) {
      int IndexN = SwapNodeType[i].first;
      NType TypeNode = SwapNodeType[i].second;
      TGraph.NodeType(IndexN) = TypeNode;
    }
  }

  void SampleCandidates(bool Init_Distances = true) {

    ClearSampledPaths();

    if (avoidSingularities)
      DisableSingularities();

    TGraph.DisableNodeType(TangentDomain);

    // save the distance on another flag
    TGraph.CurrDistIndex() = 1;

    SamplePaths(NumSamples);

    TGraph.CurrDistIndex() = 0;

    if (avoidSingularities)
      EnableSingularities();

    if (Init_Distances) {
      TGraph.EnableAllNodes();
      InitializeDistance();
    }

    InitCandidates();
    UpdateCandidatePriority();

    WriteTiming();
  }

  void Init() {
    must_pass_on = std::vector<bool>(VertPos.size(), true);
    sample_on = std::vector<bool>(VertPos.size(), true);
    force_sample = std::vector<bool>(VertPos.size(), false);
    allow_path_on = std::vector<bool>(VertPos.size(), true);
    allow_arrive_to = std::vector<bool>(VertPos.size(), true);

    // InitMarkVector();

    // InitEmitReceive();

    TGraph.EnableAllNodes();

    // InitBorderNodes();
    InitSingularityNodes();

    Geo::getVFAdjacency(Connectivity, VF, VFIndex);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);
    Geo::getBorderVert(Connectivity, NextF, VertPos.size(), IsBorder);
  }
  //     SampleCandidates(true);

  //     WriteTiming();

  // }

  void GetStartingNodes(std::vector<int> &StartIDX) {
    StartIDX = StartingNodes;
    // StartIDX.clear();
    // for (size_t i=0;i<StartingNodes.size();i++)
    // {
    //     int IndexV=Field::CrossFieldGraph::ElemOfNode(StartingNodes[i]);
    //     StartIDX.push_back(IndexV);
    // }
  }

  void GetVertexSelectedPaths(std::vector<std::vector<int>> &VertPaths) {
    VertPaths.clear();
    VertPaths.resize(SelectedPaths.size());

    for (size_t i = 0; i < SelectedPaths.size(); i++)
      for (size_t j = 0; j < SelectedPaths[i].size(); j++) {
        int IndexV = Field::CrossFieldGraph::ElemOfNode(SelectedPaths[i][j]);
        VertPaths[i].push_back(IndexV);
      }
  }

  void GetVertexSelectedPaths(const std::vector<int> &IndexPaths,
                              std::vector<std::vector<int>> &VertPaths) {
    VertPaths.clear();
    VertPaths.resize(IndexPaths.size());

    for (size_t i = 0; i < IndexPaths.size(); i++) {
      int currPIdx = IndexPaths[i];
      assert(currPIdx >= 0);
      assert(currPIdx < SelectedPaths.size());
      for (size_t j = 0; j < SelectedPaths[currPIdx].size(); j++) {
        int IndexV =
            Field::CrossFieldGraph::ElemOfNode(SelectedPaths[currPIdx][j]);
        VertPaths[i].push_back(IndexV);
      }
    }
  }

  void SplitPaths(const std::set<int> &FixedPaths = std::set<int>()) {
    // for each vertex the path passing by
    std::map<int, std::set<int>> VertexPath;

    // find the splitting vertices
    for (int i = 0; i < (int)SelectedPaths.size(); i++)
      for (int j = 0; j < (int)SelectedPaths[i].size(); j++) {
        int PathN = SelectedPaths[i][j];
        int PathV = Field::CrossFieldGraph::ElemOfNode(PathN);
        VertexPath[PathV].insert(i);
      }

    std::set<int> SplitVert;
    std::map<int, std::set<int>>::iterator IteM;
    for (IteM = VertexPath.begin(); IteM != VertexPath.end(); IteM++) {
      if ((*IteM).second.size() == 1)
        continue;
      SplitVert.insert((*IteM).first);
    }

    std::vector<std::vector<int>> NewPaths;
    for (int i = 0; i < (int)SelectedPaths.size(); i++) {
      if (FixedPaths.count(i) > 0) {
        NewPaths.push_back(SelectedPaths[i]);
        continue;
      }
      bool is_loop = Geo::PathFunctions<ScalarType>::IsLoop(SelectedPaths[i]);
      int numInters = 0;
      int sizeP = (int)SelectedPaths[i].size() - 1;
      int first_inters = -1;
      for (int j = 0; j < sizeP; j++) {
        int PathN = SelectedPaths[i][j];
        int PathV = Field::CrossFieldGraph::ElemOfNode(PathN);

        // not loop not check first and last
        if (!is_loop) // no split begin and end
        {
          if ((j == 0) || (j == (sizeP - 1)))
            continue;
        }

        // if loop check only one among start/end
        if (is_loop) // no split begin and end
        {
          if (j == (sizeP - 1))
            continue;
        }

        if (SplitVert.count(PathV) == 0)
          continue;

        numInters++;
        first_inters = j;
      }

      // one split on the loop does nothing
      if ((is_loop) && (numInters < 2)) {
        //                NewPaths.push_back(NewPathsPrio(SelectedPaths[i],i));
        //                NewPaths.back().InitPriorityByY(VertPos);
        NewPaths.push_back(SelectedPaths[i]);
        continue;
      }

      // one split on the loop does nothing
      if ((!is_loop) && (numInters == 0)) {
        //                NewPaths.push_back(NewPathsPrio(SelectedPaths[i],i));
        //                NewPaths.back().InitPriorityByY(VertPos);
        NewPaths.push_back(SelectedPaths[i]);
        continue;
      }

      std::vector<int> to_splitPath = SelectedPaths[i];
      assert(first_inters >= 0);
      assert(first_inters < to_splitPath.size());
      if (is_loop) {
        // remove the last one
        assert(to_splitPath[0] == to_splitPath.back());
        to_splitPath.pop_back();

        if (first_inters > 0)
          std::rotate(to_splitPath.begin(), to_splitPath.begin() + first_inters,
                      to_splitPath.end());

        // replicate the last one
        to_splitPath.push_back(to_splitPath[0]);
      }

      std::vector<std::vector<int>> SplittedP;
      SplittedP.resize(1);
      for (int j = 0; j < to_splitPath.size(); j++) {
        int PathN = to_splitPath[j];
        int PathV = Field::CrossFieldGraph::ElemOfNode(PathN);
        SplittedP.back().push_back(PathN);
        if ((j == 0) || (j == (to_splitPath.size() - 1)))
          continue; // no split begin and end
        if (SplitVert.count(PathV) == 0)
          continue;
        // otherwise split
        SplittedP.resize(SplittedP.size() + 1);
        SplittedP.back().push_back(PathN);
      }
      for (size_t j = 0; j < SplittedP.size(); j++) {
        //                NewPathsPrio NPP(SplittedP[j],i);
        //                NPP.InitPriorityByY(VertPos);
        //                NewPaths.push_back(NPP);
        NewPaths.push_back(SplittedP[j]);
      }
      // NewPaths.insert(NewPaths.end(),SplittedP.begin(),SplittedP.end());
    }
    //        if (sort_path)
    //            std::sort(NewPaths.begin(),NewPaths.end());
    SelectedPaths = NewPaths;
    //        SelectedPaths.clear();
    //        for (size_t i=0;i<NewPaths.size();i++)
    //            SelectedPaths.push_back(NewPaths[i].Path);
  }

  void MergeCollinearPath() {
    // find the possible merging vertex

    // std::vector<std::vector<int> > NewPaths;
    // find the splitting vertices
    bool has_joined = false;
    do {

      has_joined = false;
      // std::map<int,std::set<int> > ExtrVertPath;
      std::map<int, std::vector<int>> ExtrVertPath;
      for (int i = 0; i < (int)SelectedPaths.size(); i++) {
        bool is_loop = Geo::PathFunctions<ScalarType>::IsLoop(SelectedPaths[i]);
        if (is_loop) {
          // NewPaths.push_back(SelectedPaths[i]);
          continue;
        }
        int Node0 = SelectedPaths[i][0];
        int Node1 = SelectedPaths[i].back();

        int PathV0 = Field::CrossFieldGraph::ElemOfNode(Node0);
        int PathV1 = Field::CrossFieldGraph::ElemOfNode(Node1);
        ExtrVertPath[PathV0].push_back(i);
        ExtrVertPath[PathV1].push_back(i);

        //                int
        //                Node0Opp=Field::CrossFieldGraph::OppositeNodeCrossField(Node0);
        //                int
        //                Node1Opp=Field::CrossFieldGraph::OppositeNodeCrossField(Node1);

        //                assert(Node0!=Node1);
        //                assert(Node0Opp!=Node1Opp);
        //                ExtrNodePath[Node0].push_back(i);
        //                ExtrNodePath[Node1].push_back(i);
        //                ExtrNodePath[Node0Opp].push_back(i);
        //                ExtrNodePath[Node1Opp].push_back(i);
      }

      std::map<int, std::vector<int>>::iterator IteM;
      for (IteM = ExtrVertPath.begin(); IteM != ExtrVertPath.end(); IteM++) {
        if ((*IteM).second.size() != 2)
          continue;

        int IndexP0 = (*IteM).second[0];
        int IndexP1 = (*IteM).second[1];

        if (SelectedPaths[IndexP0].size() == 0)
          continue;
        if (SelectedPaths[IndexP1].size() == 0)
          continue;

        assert(IndexP0 != IndexP1);

        std::vector<int> JoinedPath;
        std::vector<int> P0 = SelectedPaths[IndexP0];
        std::vector<int> P1 = SelectedPaths[IndexP1];
        has_joined = true;

        if (P0.back() == P1[0]) {

          JoinedPath = P0;
          JoinedPath.insert(JoinedPath.end(), P1.begin(), P1.end());
          SelectedPaths[IndexP0].clear();
          SelectedPaths[IndexP1].clear();
          SelectedPaths.push_back(JoinedPath);
          continue;
        }

        if (P0.back() ==
            Field::CrossFieldGraph::OppositeNodeCrossField(P1.back())) {
          JoinedPath = P0;
          Geo::PathFunctions<ScalarType>::ReversePath(TGraph, P1);
          JoinedPath.insert(JoinedPath.end(), P1.begin(), P1.end());
          SelectedPaths[IndexP0].clear();
          SelectedPaths[IndexP1].clear();
          SelectedPaths.push_back(JoinedPath);
          continue;
        }

        if (P1.back() == P0[0]) {
          JoinedPath = P1;
          JoinedPath.insert(JoinedPath.end(), P0.begin(), P0.end());
          SelectedPaths[IndexP0].clear();
          SelectedPaths[IndexP1].clear();
          SelectedPaths.push_back(JoinedPath);
          continue;
        }

        if (P1.back() ==
            Field::CrossFieldGraph::OppositeNodeCrossField(P0.back())) {
          JoinedPath = P1;
          Geo::PathFunctions<ScalarType>::ReversePath(TGraph, P0);
          JoinedPath.insert(JoinedPath.end(), P0.begin(), P0.end());
          SelectedPaths[IndexP0].clear();
          SelectedPaths[IndexP1].clear();
          SelectedPaths.push_back(JoinedPath);
          continue;
        }
        // assert(0);
      }
    } while (has_joined);
    // SelectedPaths=NewPaths;
    CleanEmptySelectedPaths();
  }

  CrossPathSampling(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                    std::vector<std::vector<int>> &_Connectivity,
                    Geo::TracingGraph<ScalarType> &_TGraph,
                    std::vector<Field::CrossF<ScalarType>> &_VertCurv,
                    std::vector<int> &_SingIndex, std::vector<int> &_SingValue)
      : VertPos(_VertPos), Connectivity(_Connectivity), TGraph(_TGraph),
        VertCurv(_VertCurv), SingIndex(_SingIndex), SingValue(_SingValue) {
    avoidSingularities = true;
    dynamicUpdates = true;
    updatePatches = true;
    // orthoBorderDist=true;
    NumSamples = 1000;
    minSingResample = 8;

    only_visible = true;
    // for (size_t i=0;i<VERT_MARK_SIZE;i++)
    // {
    //     check_pass_on[i]=false;
    //     limit_path_to[i]=false;
    //     sample_only_on[i]=false;
    // }
  }
};

} // namespace Geo
#endif
