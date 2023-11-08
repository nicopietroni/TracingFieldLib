#ifndef GRAPH_TRACING
#define GRAPH_TRACING

#include "graph.h"
#include <list>
#include <set>
#include "graph_frame_field.h"

namespace Geo{

template <class ScalarType>
class GraphQuery
{
    struct HeapEntry
    {
        ScalarType Weight;
        size_t NodeI;
        size_t NodeSource;

        bool operator < (const HeapEntry &He) const
        {
            if( Weight != He.Weight)
                return (Weight > He.Weight);
            return (NodeI<He.NodeI);
        }

        HeapEntry(size_t _NodeI,
                  ScalarType _Weight,
                  size_t _NodeSource)
        {
            Weight=_Weight;
            NodeI=_NodeI;
            NodeSource=_NodeSource;
        }
    };

public:

    static bool RetrievePath(const TracingGraph<ScalarType> &TGraph,
                             const size_t endNode,
                             const std::set<int> startN,
                             std::vector<int> &path)
    {
        assert(endNode>=0);
        assert(endNode<TGraph.NumNodes());
        // get the nodes from
        size_t curr_node = endNode;
        std::list<size_t> nodeList;
        int step=0;
        nodeList.push_back(curr_node);
        while ((startN.count(curr_node)==0)&&(step<TGraph.NumNodes()))
        {
            //std::cout<<curr_node<<std::endl;
            curr_node = TGraph.Father(curr_node);
            assert(curr_node>=0);
            assert(curr_node<TGraph.NumNodes());
            nodeList.push_back(curr_node);
            step++;
        } //&&(curr_node != endNode)&&(step<1000));
        //if is a loop the put again the first
        //if (curr_node == endNode)
        //    nodeList.push_back(curr_node);

        //std::cout<<"step:"<<step<<std::endl;
        //std::cout<<curr_node<<std::endl;
        //std::cout<<"test end"<<std::endl;
        path.assign(nodeList.rbegin(), nodeList.rend());
        return (step<TGraph.NumNodes());
    }

    struct ShortestPathParam
    {
        bool avoid_zero_path;
        bool stop_at_S0;
        bool stop_at_S1;
        ScalarType MaxL;
        //bool AvoidIntersections;

        //bool AvoidIntersections
        ShortestPathParam()
        {
            avoid_zero_path=true;
            stop_at_S0=true;
            stop_at_S1=true;
            //AvoidIntersections=false;
            MaxL=std::numeric_limits<ScalarType>::max();
        }
    };

    static bool Dijkstra(TracingGraph<ScalarType> &TGraph,
                         const std::vector<int> &SourceNodes,
                         const ShortestPathParam &SParam,
                         int &StopNode)
    {
        StopNode=-1;
        //initialize the heap
        std::vector<HeapEntry> Heap;
        //std::set<int> StartN(SourceNodes.begin(),SourceNodes.end());
        TGraph.UnMarkAll();
        for (size_t i=0;i<SourceNodes.size();i++)
        {
            size_t IndexN0=SourceNodes[i];
            TGraph.Mark(IndexN0);
            TGraph.Distance(IndexN0)=0;
            TGraph.Father(IndexN0)=-1;
            TGraph.Source(IndexN0)=IndexN0;
            assert(TGraph.EnabledN[IndexN0]);
            Heap.push_back(HeapEntry(IndexN0,0,IndexN0));
        }

        //set initial
        std::make_heap(Heap.begin(),Heap.end());
        do
        {
            std::pop_heap(Heap.begin(),Heap.end());
            size_t CurrN=(Heap.back()).NodeI;
            ScalarType CurrWeight=TGraph.Distance(CurrN);
            size_t CurrSource=TGraph.Source(CurrN);
            //size_t CurrSource=(Heap.back()).NodeSource;

            Heap.pop_back();

            if (CurrWeight>SParam.MaxL)
                return false;

            assert(TGraph.IsMarked(CurrN));

            bool check_select=SParam.stop_at_S0|| SParam.stop_at_S1;
            if (check_select)
            {
                bool HasCompletedPath=false;

                if (SParam.stop_at_S0)
                    HasCompletedPath|=TGraph.IsSelected(CurrN,0);
                if (SParam.stop_at_S1)
                    HasCompletedPath|=TGraph.IsSelected(CurrN,1);

                //this is to also include loops
                if (SParam.avoid_zero_path)
                    HasCompletedPath&=(TGraph.Distance(CurrN)>0);

                if (HasCompletedPath)
                {
                    StopNode=CurrN;
                    //                    //std::cout<<"Retrieve"<<std::endl;
                    //                    RetrievePath(TGraph,CurrN,StartN,PathN);
                    //                    //std::cout<<"End Retrieve"<<std::endl;
                    //                    PathL=TGraph.Distance(CurrN);
                    return true;
                }
            }

            //go over neighbours
            for (size_t i=0;i<TGraph.N[CurrN].size();i++)
            {
                //check if the link is enabled
                if (!TGraph.EnabledA[CurrN][i])continue;
                size_t NextN=TGraph.N[CurrN][i];
                if (!TGraph.EnabledN[NextN])continue;

                //then get the next weight
                ScalarType NextWeight=CurrWeight+TGraph.W[CurrN][i];
                ScalarType OldWeight=std::numeric_limits<ScalarType>::max();

                //check if has been explored or not
                bool IsMarked=(TGraph.IsMarked(NextN));
                if (IsMarked)
                {
                    //check for loops
                    ScalarType testDist=TGraph.Distance(NextN);
                    if (testDist>0)
                        OldWeight=testDist;
                }


                if (NextWeight<OldWeight)
                {
                    TGraph.Mark(NextN);

//                    //mark also others in case would avoid self intersect
//                    if (SParam.AvoidIntersections)
//                    {
//                        Geo::GraphType GType=TGraph.TypeofGraph();
//                        if (GType==GT_CrossField)
//                        {
//                            std::cout<<"HAVE TO BE IMPLEMENTED"<<std::endl;
//                            exit(0);
//                        }
//                        if (GType==GT_FrameField)
//                        {
//                           std::vector<int> IndexNodes;
//                           NodesOfSamePoly(NextN,IndexNodes);
//                           for (size_t j=0;j<IndexNodes.size();j++)
//                               TGraph.Mark(IndexNodes[j]);
//                        }
//                    }
                    TGraph.Distance(NextN)=NextWeight;
                    TGraph.Father(NextN)=CurrN;
                    TGraph.Source(NextN)=CurrSource;

                    Heap.push_back(HeapEntry(NextN,NextWeight,CurrSource));
                    push_heap(Heap.begin(),Heap.end());
                }
            }

        }while (!Heap.empty());
        //std::cout<<"Heap Empty!"<<std::endl;
        return false;
    }


    static bool ShortestPaths(TracingGraph<ScalarType> &TGraph,
                              const std::vector<int> &SourceNodes,
                              std::vector<int> &PathN,
                              ScalarType &PathL,
                              const ShortestPathParam &SParam)
    {
        PathL=-1;
        std::set<int> StartN(SourceNodes.begin(),SourceNodes.end());

        int StopNode=-1;
        bool completed_trace=Dijkstra(TGraph,SourceNodes,SParam,StopNode);

        if (completed_trace)
        {
            assert(StopNode>=0);
            assert(StopNode<TGraph.NumNodes());
            //std::cout<<"Retrieve"<<std::endl;
            RetrievePath(TGraph,StopNode,StartN,PathN);
            //std::cout<<"End Retrieve"<<std::endl;
            PathL=TGraph.Distance(StopNode);
            return true;
        }

        return completed_trace;
    }

    static bool IsSelected(TracingGraph<ScalarType> &TGraph,
                           const int &TestNode,const ShortestPathParam &SParam)
    {
        if ((SParam.stop_at_S0) && (TGraph.IsSelected(TestNode,0)))
            return true;
        if ((SParam.stop_at_S1) && (TGraph.IsSelected(TestNode,1)))
            return true;
        return false;
    }

    static bool TraceFrom(TracingGraph<ScalarType> &TGraph,
                          const int &SourceNode,
                          std::vector<int> &PathN,
                          ScalarType &PathL,
                          const ShortestPathParam &SParam)
    {
        assert(TGraph.EnabledN[SourceNode]);
        assert(SourceNode>=0);
        assert(SourceNode<TGraph.NumNodes());

        PathN.clear();
        PathL=0;
        int CurrN=SourceNode;
        PathN.push_back(CurrN);
        //check the first one
        if (IsSelected(TGraph,CurrN,SParam))
            return true;
        while (true)
        {
            int NextI=-1;
            for (size_t i=0;i<TGraph.N[CurrN].size();i++)
            {
                //if (IsSelected(TGraph,CurrN,SParam))
                if (!TGraph.EnabledA[CurrN][i])continue;
                NextI=i;
                break;
            }
            if (NextI==-1)return false;
            assert(TGraph.EnabledN[NextI]);
            if (!TGraph.IsSink(CurrN))
            {
                PathL+=TGraph.W[CurrN][NextI];
                CurrN=TGraph.N[CurrN][NextI];

                if (PathL>SParam.MaxL)
                    return false;

                assert(CurrN>=0);
                assert(CurrN<TGraph.NumNodes());
                PathN.push_back(CurrN);

                if (IsSelected(TGraph,CurrN,SParam))
                    return true;
            }
            else
                return false;

        }
    }


};
}
#endif
