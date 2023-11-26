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
        //get nodes from sampled vertices
        for (size_t i=0;i<SampledVert.size();i++)
        {

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

    std::vector<std::vector<int> > SampledPaths;
    bool avoidSingularities;
    bool dynamicUpdates;
    std::vector<std::pair<ScalarType,int> > Candidates;
    std::unordered_set<int> ConflictNodes;

    void SamplePaths(int NumSamples)
    {
        std::cout<<"Sampling starting Nnodes"<<std::endl;
        InitTracingNodes(NumSamples);
        std::cout<<"Sampled "<<StartingNodes.size()<<" Starting Nodes"<<std::endl;
        for (size_t i=0;i<StartingNodes.size();i++)
        {
            std::vector<int> CurrPath;
            bool has_traced=Geo::GraphTracing<ScalarType>::TraceLoopOrBorderBorderFrom(TGraph,StartingNodes[i],CurrPath);
            if (!has_traced)continue;
            bool self_intersect=Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,CurrPath);
            if (self_intersect)continue;
            SampledPaths.push_back(CurrPath);

        }
        std::cout<<"Sampled "<<SampledPaths.size()<<" Paths"<<std::endl;
    }

    void InitCandidates()
    {
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

    bool ConflictWithInserted(int IndexPath)
    {
        assert(IndexPath>=0);
        assert(IndexPath<(int)SampledPaths.size());
        for (size_t i=0;i<SampledPaths[IndexPath].size();i++)
        {
            int test_node=SampledPaths[IndexPath][i];
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

    void InitSingularityNodes()
    {
        SingNodes.clear();
        for (size_t j=0;j<SingIndex.size();j++)
        {
            std::vector<int> IndexN;
            Field::CrossFieldGraph::NodesOfElem(SingIndex[j],IndexN);
            for (size_t j=0;j<IndexN.size();j++)
            {
                if (!IsTraceable(IndexN[j]))continue;
                SingNodes.push_back(IndexN[j]);
            }
        }
    }

//    void InitInsertedPathNodes()
//    {
////        SingNodes.clear();
////        for (size_t j=0;j<SingIndex.size();j++)
////        {
////            std::vector<int> IndexN;
////            Field::CrossFieldGraph::NodesOfElem(SingIndex[j],IndexN);
////            for (size_t j=0;j<IndexN.size();j++)
////            {
////                if (!IsTraceable(IndexN[j]))continue;
////                SingNodes.push_back(IndexN[j]);
////            }
////        }
//    }

public:

    bool AddNextPath()
    {
        //start from the last one in reverse order (add the farthest)
        for (size_t i=Candidates.size()-1;i>=0;i--)
        {
            int IndexPath=Candidates[i].second;
            //in this case just check if intersect
            bool has_conflict=false;
            if (!dynamicUpdates)
                has_conflict=ConflictWithInserted(IndexPath);

            if (has_conflict)
                continue;

            //remove the one that have been marked and the inserted one
            int CurrSize=(i+1);
            Candidates.resize(CurrSize-1);

            //add the new path
            SelectedPaths.push_back(SampledPaths[IndexPath]);

            //then update distance
            InitializeDistance();

            //and update candiate priority
            UpdateCandidatePriority();

            return true;
        }
        return false;
    }

    void Init()
    {
        InitBorderNodes();

        InitSingularityNodes();

        InitializeDistance();
        if (avoidSingularities)
            DisableSingularities();

        SamplePaths(NumSamples);

        if (avoidSingularities)
            EnableSingularities();

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
        dynamicUpdates=false;
        NumSamples=1000;
    }
};

}
#endif
