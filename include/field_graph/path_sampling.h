#ifndef PATH_SAMPLING
#define PATH_SAMPLING

#define VERT_MARK_SIZE 2

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
    bool updatePatches;
    std::vector<std::pair<ScalarType,int> > Candidates;

    std::unordered_set<int> ConflictNodes;

    //these are used to detect when two path shared an edge and to solve it
    std::unordered_set<std::pair<int,int>,Geo::hashPairIntFunction > ConflictEdges;

public:

    bool dynamicUpdates;
    std::vector<std::vector<int> > SelectedPaths;

    //this used to get some subdomains to sample into subdomain or
    //to check if a path pass from some vertices/faces
    std::vector<int> VertTMark[VERT_MARK_SIZE];
    int GlobalTMark[VERT_MARK_SIZE];
    bool check_pass_on[VERT_MARK_SIZE];
    bool limit_path_to[VERT_MARK_SIZE];

    void InitMarkVector()
    {
        for (size_t i=0;i<VERT_MARK_SIZE;i++)
        {
            VertTMark[i]=std::vector<int>(VertPos.size(),0);
            GlobalTMark[i]=0;
        }
    }

private:

    ScalarType TimeTracing;
    ScalarType TimeUpdDist;
    ScalarType TimeConflict;
    ScalarType TimeInitData;
    ScalarType TimeUpdPatches;

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

        //        Geo::MeshSampling<ScalarType>::SampleVert(Connectivity,VertPos.size(),
        //                                                  RadiusStar,SampledVert,true);
        //        std::vector<std::vector<int> > NextF,NextE;
        //        Geo::getFFAdjacency(Connectivity,NextF,NextE);
        //        Geo::MeshSampling<ScalarType>::SampleBordersRegularly(Connectivity,NextF,NextE,RadiusStar,4,SampledVert);

        Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,RadiusStar,4,SampledVert);

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

        if (StartN.size()==0)
        {
            //add the ones with smaller X
            int IndexV=-1;
            ScalarType minX=std::numeric_limits<ScalarType>::max();
            for (size_t i=0;i<VertPos.size();i++)
            {
                if (VertPos[i].X>minX)continue;
                IndexV=i;
                minX=VertPos[i].X;
            }
            assert(IndexV>=0);
            StartN.push_back(IndexV);
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
        auto start = std::chrono::high_resolution_clock::now();

        std::cout<<"Sampling starting Nnodes"<<std::endl;
        InitTracingNodes(NumSamples);
        std::cout<<"Sampled "<<StartingNodes.size()<<" Starting Nodes"<<std::endl;

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        TimeInitData += duration.count() * 1e-6;

        std::vector<int> NewStartingNodes;
        for (size_t i=0;i<StartingNodes.size();i++)
        {
            std::vector<int> CurrPath;

            auto start = std::chrono::high_resolution_clock::now();

            bool has_traced=Geo::GraphTracing<ScalarType>::TraceLoopOrBorderBorderFrom(TGraph,StartingNodes[i],CurrPath);

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            TimeTracing += duration.count() * 1e-6;

            if (!has_traced)continue;

            start = std::chrono::high_resolution_clock::now();

            bool self_intersect=Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,CurrPath);

            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

            TimeConflict += duration.count() * 1e-6;

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

    bool ConflictEdgeWithInsertedEdges(const std::vector<int> &TestPath)
    {
        for (size_t i=0;i<TestPath.size()-1;i++)
        {
            int IndexV0=Field::CrossFieldGraph::ElemOfNode(TestPath[i]);
            int IndexV1=Field::CrossFieldGraph::ElemOfNode(TestPath[i+1]);
            std::pair<int,int> key(std::min(IndexV0,IndexV1),
                                   std::max(IndexV0,IndexV1));
            if(ConflictEdges.count(key)>0)return true;
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

    void AddToConflictEdges(const std::vector<int> &NodesPath)
    {
        for (size_t i=0;i<NodesPath.size()-1;i++)
        {
            int IndexV0=Field::CrossFieldGraph::ElemOfNode(NodesPath[i]);
            int IndexV1=Field::CrossFieldGraph::ElemOfNode(NodesPath[i+1]);
            std::pair<int,int> key(std::min(IndexV0,IndexV1),
                                   std::max(IndexV0,IndexV1));
            ConflictEdges.insert(key);
        }
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
            self_intersect|=Geo::PathFunctions<ScalarType>::HasSelfEdgeTwice(SampledPaths[i]);
            if (self_intersect)
            {
                SampledPaths_Status[i]=SP_Discarded;
                continue;
            }
        }
    }

    void WriteTiming()
    {
        std::cout<<"*** TIMINGS ***"<<std::endl;
        std::cout<<"* Time Init:"<<TimeInitData<<std::endl;
        std::cout<<"* Time tracing:"<<TimeTracing<<std::endl;
        std::cout<<"* Time upd Dist:"<<TimeUpdDist<<std::endl;
        std::cout<<"* Time Conflicts:"<<TimeConflict<<std::endl;
    }

    bool IsMarkeVert(const int &IndexV,const int &IndexMark)
    {
        assert(IndexMark>=0);
        assert(IndexMark<VERT_MARK_SIZE);
        assert(IndexV>=0);
        assert(IndexV<VertPos.size());
        return (VertTMark[IndexMark][IndexV]==GlobalTMark[IndexMark]);
    }

    bool IsMarkedNode(const int &IndexNode,int IndexMark)
    {
        assert(IndexMark>=0);
        assert(IndexMark<VERT_MARK_SIZE);
        assert(IndexNode>=0);
        assert(IndexNode<TGraph.NumNodes());
        int IndexV=Field::CrossFieldGraph::ElemOfNode(IndexNode);
        return (IsMarkeVert(IndexV,IndexMark));
    }


    //    void MarkNode(const int &IndexNode,int IndexMark)
    //    {
    //        assert(IndexMark>=0);
    //        assert(IndexMark<VERT_MARK_SIZE);
    //        assert(IndexNode>=0);
    //        assert(IndexNode<TGraph.NumNodes());
    //        int IndexV=Field::CrossFieldGraph::ElemOfNode(IndexNode);
    //        VertTMark[IndexMark][IndexV]=GlobalTMark[IndexV];
    //    }

public:

    void UnMarkAll(int IndexMark)
    {
        assert(IndexMark>=0);
        assert(IndexMark<VERT_MARK_SIZE);
        GlobalTMark[IndexMark]++;
    }

    void MarkVertices(const std::vector<int> &To_MarkVert,int IndexMark)
    {
        assert(IndexMark>=0);
        assert(IndexMark<VERT_MARK_SIZE);
        for (size_t i=0;i<To_MarkVert.size();i++)
        {
            int IndexV=To_MarkVert[i];
            assert(IndexV>=0);
            assert(IndexV<VertTMark[IndexMark].size());
            assert(IndexV<VertPos.size());
            VertTMark[IndexMark][IndexV]=GlobalTMark[IndexMark];
        }
    }

    bool HasMarkedVert(int &IndexCandidate,int IndexMark)
    {
        assert(IndexCandidate<SampledPaths.size());
        for (size_t i=0;i<SampledPaths[IndexCandidate].size();i++)
            if (IsMarkedNode(SampledPaths[IndexCandidate][i],IndexMark))return true;
        return false;
    }

    void ClearSamples()
    {
        StartingNodes.clear();
        SampledPaths.clear();
        SampledPaths_Status.clear();
        Candidates.clear();
        ConflictNodes.clear();
        SelectedPaths.clear();
    }

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

            //CHECK OF PASSING ON REQUIRED PATCHES
            bool pass_on_patches=true;
            for (size_t j=0;j<VERT_MARK_SIZE;j++)
                if (check_pass_on[j])pass_on_patches=false;

            if (!pass_on_patches)
            {
                //std::cout<<"checking useful path"<<std::endl;
                for (size_t j=0;j<VERT_MARK_SIZE;j++)
                {
                    if (check_pass_on[j] && (HasMarkedVert(IndexPath,j)))
                    {
                        pass_on_patches=true;
                        break;
                    }
                }
            }

            if (!pass_on_patches)
            {
                SampledPaths_Status[IndexPath]=SP_Discarded;
                continue;
            }

            //CHECK OF INTERSECTION
            bool has_conflict=false;

            auto start = std::chrono::high_resolution_clock::now();

            if (!dynamicUpdates)
                has_conflict=ConflictWithInserted(SampledPaths[IndexPath]);

            if (!has_conflict)
                has_conflict|=ConflictEdgeWithInsertedEdges(SampledPaths[IndexPath]);

            if (has_conflict)
            {
                SampledPaths_Status[IndexPath]=SP_Discarded;
                continue;
            }


            //add the new path
            SampledPaths_Status[IndexPath]=SP_Selected;
            SelectedPaths.push_back(SampledPaths[IndexPath]);

            AddToConflictNodes(SelectedPaths.back());
            AddToConflictEdges(SelectedPaths.back());

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            TimeConflict += duration.count() * 1e-6;

            //UPDATE DISTANCES WITH LAST ADDED
            if (avoidSingularities)
                EnableSingularities();

            start = std::chrono::high_resolution_clock::now();
            UpdateDistancesWithLastInsertedPath();

            if (avoidSingularities)
                DisableSingularities();

            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            TimeUpdDist += duration.count() * 1e-6;

            //UPDATE CANDIDATES IF NEEDED
            if (dynamicUpdates)
            {
                DisableGraphNodes(SelectedPaths.back());

                TGraph.CurrDistIndex()=1;
                start = std::chrono::high_resolution_clock::now();

                RetraceCandidatesConflictWithSelected();

                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                TimeTracing += duration.count() * 1e-6;

                TGraph.CurrDistIndex()=0;
            }

            //UPDATE PRIORITY
            UpdateCandidatePriority();

            return true;
        }
        return false;
    }

    void GetMarkedVertices(std::vector<int> &MarkedV,
                           int IndexMark)
    {
        assert(IndexMark>=0);
        assert(IndexMark<VERT_MARK_SIZE);
        MarkedV.clear();
        for (size_t i=0;i<VertPos.size();i++)
            if (IsMarkeVert(i,IndexMark))MarkedV.push_back(i);
    }

    void InsertAllPaths()
    {
        while(AddNextPath()){};
        WriteStats();
        WriteTiming();
    }

    void Init()
    {
        TimeTracing=0;
        TimeUpdDist=0;
        TimeConflict=0;
        TimeInitData=0;
        TimeUpdPatches=0;

        auto start = std::chrono::high_resolution_clock::now();

        TGraph.EnableAllNodes();

        ClearSamples();

        InitBorderNodes();

        InitSingularityNodes();

        if (avoidSingularities)
            DisableSingularities();

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        TimeInitData += duration.count() * 1e-6;

        SamplePaths(NumSamples);

        if (avoidSingularities)
            EnableSingularities();

        start = std::chrono::high_resolution_clock::now();

        InitializeDistance();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        TimeUpdDist += duration.count() * 1e-6;

        InitCandidates();

        UpdateCandidatePriority();

        start = std::chrono::high_resolution_clock::now();


        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        TimeUpdPatches += duration.count() * 1e-6;

        InitMarkVector();

        WriteTiming();

    }

    void GetVertexSelectedPaths(std::vector<std::vector<int> > &VertPaths)
    {
        VertPaths.clear();
        VertPaths.resize(SelectedPaths.size());

        for (size_t i=0;i<SelectedPaths.size();i++)
            for (size_t j=0;j<SelectedPaths[i].size();j++)
            {
                int IndexV=Field::CrossFieldGraph::ElemOfNode(SelectedPaths[i][j]);
                VertPaths[i].push_back(IndexV);
            }
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
        updatePatches=true;
        NumSamples=1000;
        TimeTracing=0;
        TimeUpdDist=0;
        TimeConflict=0;
        for (size_t i=0;i<VERT_MARK_SIZE;i++)
        {
            check_pass_on[i]=false;
            limit_path_to[i]=false;
        }

    }
};

}
#endif
