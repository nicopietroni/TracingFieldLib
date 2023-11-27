#ifndef PATH_SAMPLING
#define PATH_SAMPLING

#include "graph.h"
#include <list>
#include <set>
#include <unordered_set>
#include "graph_cross_field.h"
#include "../mesh_sampling.h"
#include "tracing_interface.h"
#include "path_functions.h"
#include "./connectivity.h"

namespace Geo{

template <class ScalarType>
class CrossPathSampling
{

    std::vector<Geo::Point3<ScalarType > > &VertPos;
    std::vector<std::vector<int> > &Connectivity;
    Geo::TracingGraph<ScalarType> &TGraph;
    std::vector<Field::CrossF<ScalarType> > &VertCurv;

    std::vector<int> &SingIndex,&SingValue;
    std::vector<int> SingNodes;

    std::vector<size_t> BorderV;
    std::vector<int> BorderNodes;

    std::vector<int> StartingNodes;
    //std::vector<ScalarType> NodeDistances;
    int NumSamples;


    std::vector<std::vector<int> > SampledPaths;
    enum SP_Status{SP_Selected,SP_Discarded,SP_Evaluate};
    std::vector<SP_Status> SampledPaths_Status;

    bool avoidSingularities;
    bool dynamicUpdates;
    std::vector<std::pair<ScalarType,int> > Candidates;

    std::unordered_set<int> ConflictNodes;

public:

    std::vector<std::vector<int> > SelectedPaths;

private:

    bool IsTraceable(int &IndexN)
    {
        assert(IndexN>=0);
        assert(IndexN<TGraph.NumNodes());
        if (TGraph.NodeType(IndexN)==Geo::ExitDomain)return false;
        if (TGraph.NodeType(IndexN)==Geo::TangentDomain)return false;
        return true;
    }

    void InitTracingNodes(int NumSources)
    {
        int RadiusStar=VertPos.size()/NumSources;
        std::vector<int> SampledVert;
        Geo::MeshSampling<ScalarType>::SampleVert(Connectivity,VertPos.size(),
                                                  RadiusStar,SampledVert,true);

        std::cout<<"Sampled "<<SampledVert.size()<<" Starting vertices"<<std::endl;
        std::set<int> SingSetId(SingIndex.begin(),SingIndex.end());
        //get nodes from sampled vertices
        for (size_t i=0;i<SampledVert.size();i++)
        {
            //not start tracing from singularities
            if (avoidSingularities)
                if (SingSetId.count(SampledVert[i])>0)continue;

            std::vector<int> IndexN;
            Field::CrossFieldGraph::NodesOfElem(SampledVert[i],IndexN);
            bool AddedDir2[2]={false,false};
            //then put two orthogonal ones
            for (size_t j=0;j<IndexN.size();j++)
            {
                if (AddedDir2[j%2]==true)continue;//already added in this or opposite dir
                if (!IsTraceable(IndexN[j]))continue;
                AddedDir2[j%2]=true;
                StartingNodes.push_back(IndexN[j]);
            }
        }
    }

    void DisableSingularities()
    {
        for (size_t i=0;i<SingIndex.size();i++)
        {
            std::vector<int> IndexN;
            Field::CrossFieldGraph::NodesOfElem(SingIndex[i],IndexN);
            for (size_t j=0;j<IndexN.size();j++)
                TGraph.DisableNode(IndexN[j]);
        }
    }

    void EnableSingularities()
    {
        for (size_t i=0;i<SingIndex.size();i++)
        {
            std::vector<int> IndexN;
            Field::CrossFieldGraph::NodesOfElem(SingIndex[i],IndexN);
            for (size_t j=0;j<IndexN.size();j++)
                TGraph.EnableNode(IndexN[j]);
        }
    }

    void InitializeDistance()
    {
        TGraph.CurrDistIndex()=0;

        std::vector<int> StartN;

        //first add singularities
        if (avoidSingularities)
            StartN.insert(StartN.end(),SingNodes.begin(),SingNodes.end());

        //add border nodes
        StartN.insert(StartN.end(),BorderNodes.begin(),BorderNodes.end());

        //then add all the previously added nodes
        for (size_t i=0;i<SelectedPaths.size();i++)
        {
            //add the path and opposite
            StartN.insert(StartN.end(),SelectedPaths[i].begin(),SelectedPaths[i].end());
            std::vector<int> OppNodes;
            Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths[i],OppNodes);
            StartN.insert(StartN.end(),OppNodes.begin(),OppNodes.end());
        }

        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
        SParam.avoid_zero_path=true;
        for (size_t i=0;i<SELECTED_GRAPH_FLAG;i++)
            SParam.StopAtSel[i]=false;
        SParam.MaxL=std::numeric_limits<ScalarType>::max();
        SParam.onlyDirect=false;
        SParam.onlyVisible=false;
        SParam.stopAtPrevDist=false;
        SParam.stopAtPrevDist=false;

        int StopNode;
        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,StartN,SParam,StopNode);
    }

    ScalarType Distance(const std::vector<int> &TestPath)
    {
        ScalarType CurrD=0;
        for (size_t i=0;i<TestPath.size();i++)
        {
            int N0=TestPath[i];
            int N1=Field::CrossFieldGraph::OppositeNodeCrossField(N0);
            assert(N0>=0);
            assert(N1>=0);
            assert(N0<TGraph.NumNodes());
            assert(N1<TGraph.NumNodes());
            CurrD+=TGraph.Distance(N0);
            CurrD+=TGraph.Distance(N1);
        }
        CurrD/=(ScalarType)TestPath.size()*2;
        return CurrD;
    }


    void SamplePaths(int NumSamples)
    {
        std::cout<<"Sampling starting Nnodes"<<std::endl;
        InitTracingNodes(NumSamples);
        std::cout<<"Sampled "<<StartingNodes.size()<<" Starting Nodes"<<std::endl;

        std::vector<int> NewStartingNodes;
        for (size_t i=0;i<StartingNodes.size();i++)
        {
            std::vector<int> CurrPath;
            bool has_traced=Geo::GraphTracing<ScalarType>::TraceLoopOrBorderBorderFrom(TGraph,StartingNodes[i],CurrPath);
            if (!has_traced)continue;
            bool self_intersect=Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,CurrPath);
            if (self_intersect)continue;
            SampledPaths.push_back(CurrPath);
            NewStartingNodes.push_back(StartingNodes[i]);
        }
        StartingNodes=NewStartingNodes;
        SampledPaths_Status.clear();
        SampledPaths_Status.resize(SampledPaths.size(),SP_Evaluate);
        std::cout<<"Sampled "<<SampledPaths.size()<<" Paths"<<std::endl;
    }

    void InitCandidates()
    {
        Candidates.clear();
        for (size_t i=0;i<SampledPaths.size();i++)
            Candidates.push_back(std::pair<ScalarType,int>(-1,i));
    }

    void UpdateCandidatePriority()
    {
        for (size_t i=0;i<Candidates.size();i++)
        {
            int IndexPath=Candidates[i].second;
            assert(IndexPath>=0);
            assert(IndexPath<SampledPaths.size());
            ScalarType Dist=Distance(SampledPaths[IndexPath]);
            Candidates[i].first=Dist;
        }
        std::sort(Candidates.begin(),Candidates.end());
    }

    bool ConflictWithInserted(const std::vector<int> &TestPath)
    {

        for (size_t i=0;i<TestPath.size();i++)
        {
            int test_node=TestPath[i];
            if (ConflictNodes.count(test_node)>0)return true;
        }
        return false;
    }

    void InitBorderNodes()
    {
        std::vector<std::vector<int> > NextF,NextE;
        Geo::getFFAdjacency(Connectivity,NextF,NextE);

        //then add the borders
        Geo::getBorderVert(Connectivity,NextF,BorderV);

        BorderNodes.clear();
        for (size_t j=0;j<BorderV.size();j++)
        {
            std::vector<int> IndexN;
            Field::CrossFieldGraph::NodesOfElem(BorderV[j],IndexN);
            for (size_t j=0;j<IndexN.size();j++)
            {
                if (!IsTraceable(IndexN[j]))continue;
                BorderNodes.push_back(IndexN[j]);
            }
        }
    }

    void InitSingularityNodes(bool use_star=true)
    {
        std::vector<int> VertToUse=SingIndex;
        //add sing star in case
        if (use_star)
        {
            std::vector<std::vector<int> > Kernel;
            Geo::FindVertexStar(Connectivity,VertPos.size(),Kernel);
            for (size_t j=0;j<SingIndex.size();j++)
            {
                int IndexV=SingIndex[j];
                assert(IndexV>=0);
                assert(IndexV<Kernel.size());
                VertToUse.insert(VertToUse.end(),Kernel[IndexV].begin(),Kernel[IndexV].end());
            }
        }
        //move one to
        SingNodes.clear();
        for (size_t j=0;j<VertToUse.size();j++)
        {
            std::vector<int> IndexN;
            Field::CrossFieldGraph::NodesOfElem(VertToUse[j],IndexN);
            for (size_t j=0;j<IndexN.size();j++)
            {
                if (!IsTraceable(VertToUse[j]))continue;
                SingNodes.push_back(VertToUse[j]);
            }
        }
    }

    void AddToConflictNodes(const std::vector<int> &NodesPath)
    {
        std::vector<int> OppNodesPath;
        Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath,OppNodesPath);
        ConflictNodes.insert(NodesPath.begin(),NodesPath.end());
        ConflictNodes.insert(OppNodesPath.begin(),OppNodesPath.end());
    }

    void DisableGraphNodes(const std::vector<int> &NodesPath)
    {
        std::vector<int> To_Disable=NodesPath;
        std::vector<int> OppNodesPath;
        Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath,OppNodesPath);
        To_Disable.insert(To_Disable.end(),OppNodesPath.begin(),OppNodesPath.end());
        TGraph.DisableNodes(To_Disable);
    }

    void UpdateDistancesWithLastInsertedPath()
    {
        std::vector<int> StartN;
        StartN.insert(StartN.end(),SelectedPaths.back().begin(),SelectedPaths.back().end());
        std::vector<int> OppNodes;
        Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths.back(),OppNodes);
        StartN.insert(StartN.end(),OppNodes.begin(),OppNodes.end());

        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
        SParam.avoid_zero_path=true;
        for (size_t i=0;i<SELECTED_GRAPH_FLAG;i++)
            SParam.StopAtSel[i]=false;
        SParam.MaxL=std::numeric_limits<ScalarType>::max();
        SParam.onlyDirect=false;
        SParam.onlyVisible=false;
        SParam.stopAtPrevDist=true;

        int StopNode;
        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,StartN,SParam,StopNode);
    }

    void RetraceCandidatesConflictWithSelected()
    {
        assert(SampledPaths_Status.size()==SampledPaths.size());
        assert(StartingNodes.size()==SampledPaths.size());
        for (size_t i=0;i<SampledPaths.size();i++)
        {
            if (SampledPaths_Status[i]!=SP_Evaluate)continue;
            bool has_conflict=ConflictWithInserted(SampledPaths[i]);
            if (!has_conflict)continue;//fine no need to retrace

            //check if the start nodes is avalable
            if (ConflictNodes.count(StartingNodes[i])>0)
            {
                assert(!TGraph.IsEnabled(StartingNodes[i]));
                SampledPaths_Status[i]=SP_Discarded;
                continue;
            }

            bool has_traced=Geo::GraphTracing<ScalarType>::TraceLoopOrBorderBorderFrom(TGraph,StartingNodes[i],SampledPaths[i]);
            if (!has_traced)
            {
                SampledPaths_Status[i]=SP_Discarded;
                continue;
            }
            bool self_intersect=Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,SampledPaths[i]);
            if (self_intersect)
            {
                SampledPaths_Status[i]=SP_Discarded;
                continue;
            }
        }
    }

public:

    void WriteStats()
    {
        std::cout<<"Added Total:" <<SelectedPaths.size()<<" paths"<<std::endl;
        int NumNodes=0;
        for (int i=0;i<SelectedPaths.size();i++)
            NumNodes+=SelectedPaths[i].size();
        std::cout<<"Nodes Total:" <<NumNodes<<std::endl;
    }

    bool AddNextPath()
    {
        if (Candidates.size()==0)return false;

        TGraph.CurrDistIndex()=0;
        //start from the last one in reverse order (add the farthest)
        for (int i=Candidates.size()-1;i>=0;i--)
        {

            int IndexPath=Candidates[i].second;

            if (SampledPaths_Status[IndexPath]!=SP_Evaluate)continue;

            //in this case just check if intersect
            bool has_conflict=false;
            if (!dynamicUpdates)
                has_conflict=ConflictWithInserted(SampledPaths[IndexPath]);

            if (has_conflict)
            {
                SampledPaths_Status[IndexPath]=SP_Discarded;
                continue;
            }

            //add the new path
            SampledPaths_Status[IndexPath]=SP_Selected;
            SelectedPaths.push_back(SampledPaths[IndexPath]);

            AddToConflictNodes(SelectedPaths.back());

            //then update distance
            //InitializeDistance();

            if (avoidSingularities)
                EnableSingularities();

            UpdateDistancesWithLastInsertedPath();

            if (avoidSingularities)
                DisableSingularities();

            if (dynamicUpdates)
            {
                DisableGraphNodes(SelectedPaths.back());

                //              for (size_t i=0;i<SelectedPaths.size();i++)
                //                  DisableGraphNodes(SelectedPaths[i]);

                TGraph.CurrDistIndex()=1;
                RetraceCandidatesConflictWithSelected();
                TGraph.CurrDistIndex()=0;

                //                TGraph.EnableAllNodes();
                //                DisableSingularities();
            }

            //and update candiate priority
            UpdateCandidatePriority();

            return true;
        }
        return false;
    }

    void InsertAllPaths()
    {
        while(AddNextPath()){};
        WriteStats();
    }

    void Init()
    {
        InitBorderNodes();

        InitSingularityNodes();

        if (avoidSingularities)
            DisableSingularities();

        SamplePaths(NumSamples);

        if (avoidSingularities)
            EnableSingularities();

        InitializeDistance();

        InitCandidates();

        UpdateCandidatePriority();

    }

    CrossPathSampling(std::vector<Geo::Point3<ScalarType > > &_VertPos,
                      std::vector<std::vector<int> > &_Connectivity,
                      Geo::TracingGraph<ScalarType> &_TGraph,
                      std::vector<Field::CrossF<ScalarType> > &_VertCurv,
                      std::vector<int> &_SingIndex,
                      std::vector<int> &_SingValue)
        :VertPos(_VertPos),Connectivity(_Connectivity),
          TGraph(_TGraph),VertCurv(_VertCurv),
          SingIndex(_SingIndex),SingValue(_SingValue)
    {
        avoidSingularities=true;
        dynamicUpdates=true;
        NumSamples=1000;
    }
};

}
#endif
