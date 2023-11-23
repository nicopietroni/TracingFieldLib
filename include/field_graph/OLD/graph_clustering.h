#ifndef GRAPH_CLUSTERING
#define GRAPH_CLUSTERING

#include "graph.h"
#include "graph_tracing.h"
#include "graph_frame_field.h"

template <class ScalarType>
class GraphClustering
{
    static void PossibleArrivalNodes(const size_t &startNode,
                                     std::vector<size_t> &ArrNodes)
    {
        ArrNodes.clear();
        ArrNodes.push_back(startNode);
        //OppositeNodeFrameField()
    }

public:
    static void SelectPossibleArrivalNodes(const size_t &startNode,
                                           TracingGraph<ScalarType> &TGraph,
                                           int IndexSel=0)
    {
        std::vector<int> OrthoN;
        OrthoNodesFrameField((int)startNode,OrthoN);

        TGraph.Select(startNode,IndexSel);
        for (size_t j=0;j<OrthoN.size();j++)
            TGraph.Select(OrthoN[j],IndexSel);
    }

    static void UnselectPossibleArrivalNodes(const size_t &startNode,
                                             TracingGraph<ScalarType> &TGraph,
                                             int IndexSel=0)
    {
        std::vector<int> OrthoN;
        OrthoNodesFrameField((int)startNode,OrthoN);

        TGraph.Unselect(startNode,IndexSel);
        for (size_t j=0;j<OrthoN.size();j++)
            TGraph.Unselect(OrthoN[j],IndexSel);
    }


    static void SetupTracing(TracingGraph<ScalarType> &TGraph)
    {
        //select 0 the end points
        TGraph.Unselect(ExitDomain,0);
        TGraph.Unselect(ExitDomain,1);
        TGraph.SelectNoteType(ExitDomain,1);

//        TGraph.EnableAllNodes();
//        TGraph.DisableNodeType(TangentDomain);
//        TGraph.DisableArconDisabledNodeType();
    }


    //    static bool FindPath(TracingGraph<ScalarType> &TGraph,
    //                        int IndexN,std::vector<int> &PathN)
    //    {
    //            assert(!TGraph.IsSink(IndexN));
    //            assert(IndexN<TGraph.NumNodes());

    //            PathN.clear();
    //            SelectPossibleArrivalNodes(IndexN,TGraph);
    //            std::vector<int> SourceNodes(1,IndexN);
    //            bool traced=GraphQuery<ScalarType>::ShortestPaths(TGraph,SourceNodes,PathN,true,true);
    //            UnselectPossibleArrivalNodes(IndexN,TGraph);
    //            if (traced)
    //                return true;
    //                //non_traced++;
    //            else
    //                return false;
    //                //has_traced++;
    //    }

//    static void FindCode(TracingGraph<ScalarType> &TGraph,int IndexPoly,
//                         std::vector<std::vector<int> > &Paths,
//                         std::vector<ScalarType > &Lengths,
//                         std::vector<int> &Code,
//                         std::vector<bool> &Traced,
//                         size_t printMsg=100)
//    {
//        if ((printMsg>0)&&((IndexPoly%printMsg)==0))
//            std::cout<<"Gatherig Code of Polygon:"<<IndexPoly<<std::endl;

//        assert(IndexPoly>=0);
//        std::vector<int> PolyNodes;
//        NodesOfPoly(IndexPoly,PolyNodes);
//        assert(PolyNodes.size()==6);

//        Paths.clear();
//        Paths.resize(6);
//        Lengths.clear();
//        Lengths.resize(6,-1);
//        Code.clear();
//        Code.resize(6,-2);
//        Traced.clear();
//        Traced.resize(6,false);
//        for (size_t i=0;i<6;i++)
//        {
//            int IndexN=PolyNodes[i];
//            //std::cout<<"Tracing from:"<<IndexN<<std::endl;
//            assert(IndexN<TGraph.NumNodes());
//            if (TGraph.IsSink(IndexN))
//                continue;

//            //            //trace possible loop
//            //            SelectPossibleArrivalNodes(IndexN,TGraph,1);
//            //            std::vector<int> SourceNodes(1,IndexN);
//            //            std::vector<int> Path0;
//            //            typename GraphQuery<ScalarType>::ShortestPathParam SParam0;
//            //            SParam0.stop_at_S0=false;
//            //            SParam0.stop_at_S1=true;
//            //            ScalarType Lenght0=0;
//            bool traced0=false;//GraphQuery<ScalarType>::ShortestPaths(TGraph,SourceNodes,Path0,Lenght0,SParam0);
//            //            UnselectPossibleArrivalNodes(IndexN,TGraph,1);
//            //std::cout<<"L0:"<<Lenght0<<std::endl;

//            std::vector<int> Path1;
//            typename GraphQuery<ScalarType>::ShortestPathParam SParam1;
//            SParam1.stop_at_S0=false;
//            SParam1.stop_at_S1=true;
//            //SParam1.MaxL=Lenght0;
//            ScalarType Lenght1=0;
//            bool traced1=GraphQuery<ScalarType>::TraceFrom(TGraph,IndexN,Path1,Lenght1,SParam1);
//            Paths.push_back(Path1);
//            Traced[i]=traced1;
//            //std::cout<<"L1:"<<Lenght1<<std::endl;

//            //            int AddPath=-1;
//            //            if ((traced0)&&(traced1))
//            //                AddPath = (Lenght0 < Lenght1) ? 0 : 1;
//            //            if ((traced0)&&(!traced1))
//            //                AddPath=0;
//            //            if ((!traced0)&&(traced1))
//            //                AddPath=1;

//            //            if (AddPath==0)
//            //            {
//            //                Paths.push_back(Path0);
//            //                Lengths.push_back(Lenght0);
//            //                Code[i]=0;
//            //            }
//            //            if (AddPath==1)
//            //            {
//            //                Paths.push_back(Path1);
//            //                Lengths.push_back(Lenght1);
//            //                Code[i]=1;
//            //            }
//            //            if (AddPath==-1)
//            //                Code[i]=-1;
//        }
//    }

    static void FindCodeLoop(TracingGraph<ScalarType> &TGraph,int IndexPoly,
                         std::vector<std::vector<int> > &Paths,
                         std::vector<ScalarType > &MaxLengths,
                         std::vector<int> &Code,
                         std::vector<bool> &Traced)
    {

        assert(MaxLengths.size()==6);
        assert(IndexPoly>=0);
        std::vector<int> PolyNodes;
        NodesOfPoly(IndexPoly,PolyNodes);
        assert(PolyNodes.size()==6);

        Paths.clear();
        Paths.resize(6);
        Code.clear();
        Code.resize(6,-2);
        Traced.clear();
        Traced.resize(6,false);

        for (size_t i=0;i<6;i++)
        {
            int IndexN=PolyNodes[i];
            if (!TGraph.EnabledN[IndexN])continue;

            assert(IndexN<TGraph.NumNodes());
            if (TGraph.IsSink(IndexN))
                continue;

            //trace possible loop
            SelectPossibleArrivalNodes(IndexN,TGraph,0);
            std::vector<int> SourceNodes(1,IndexN);
            std::vector<int> PathN;
            typename GraphQuery<ScalarType>::ShortestPathParam SParam;
            SParam.stop_at_S0=true;
            SParam.stop_at_S1=false;
            SParam.MaxL=MaxLengths[i];
            ScalarType LenghtTest=0;
            bool has_traced=GraphQuery<ScalarType>::ShortestPaths(TGraph,SourceNodes,PathN,LenghtTest,SParam);
            UnselectPossibleArrivalNodes(IndexN,TGraph,0);

            if (has_traced)
                assert(LenghtTest<SParam.MaxL);

            MaxLengths[i]=LenghtTest;
            Paths.push_back(PathN);
            Traced[i]=has_traced;
        }
    }

    static void FindCodeTrace(TracingGraph<ScalarType> &TGraph,
                              int IndexPoly,
                              std::vector<std::vector<int> > &Paths,
                              std::vector<ScalarType > &MaxLengths,
                              std::vector<int> &Code,
                              std::vector<bool> &Traced)
    {
        assert(MaxLengths.size()==6);
        assert(IndexPoly>=0);
        std::vector<int> PolyNodes;
        NodesOfPoly(IndexPoly,PolyNodes);
        assert(PolyNodes.size()==6);

        Paths.clear();
        Paths.resize(6);
        Code.clear();
        Code.resize(6,-2);
        Traced.clear();
        Traced.resize(6,false);
        for (size_t i=0;i<6;i++)
        {
            int IndexN=PolyNodes[i];
            if (!TGraph.EnabledN[IndexN])continue;

            assert(IndexN<TGraph.NumNodes());

            std::vector<int> PathN;
            typename GraphQuery<ScalarType>::ShortestPathParam SParam;

            SParam.stop_at_S0=false;
            SParam.stop_at_S1=true;
            SParam.MaxL=MaxLengths[i];
            ScalarType LenghtTest=0;
            bool has_traced=GraphQuery<ScalarType>::TraceFrom(TGraph,IndexN,PathN,LenghtTest,SParam);
            if (has_traced)
                assert(LenghtTest<SParam.MaxL);

            MaxLengths[i]=LenghtTest;
            Paths.push_back(PathN);
            Traced[i]=has_traced;
        }
    }

    static void FindCodes(TracingGraph<ScalarType> &TGraph,
                          int NumPoly)
    {
        size_t non_traced=0;
        size_t has_traced=0;

        //select 0 the end points
//        TGraph.Unselect(ExitDomain,0);
//        TGraph.Unselect(ExitDomain,1);
//        TGraph.SelectNoteType(ExitDomain,1);
        SetupTracing(TGraph);

        for (size_t i=0;i<NumPoly;i++)
        {
            if ((i%100)==0)
                std::cout<<"Gatherig Code of Polygon:"<<i<<std::endl;

            std::vector<std::vector<int> > Paths;
            std::vector<int> Code;
            std::vector<bool> Traced;
            std::vector<ScalarType> MaxLengths(6,std::numeric_limits<ScalarType>::max());
            //FindCodeLoop(TGraph,i,Paths,MaxLengths,Code,Traced);
            FindCodeTrace(TGraph,i,Paths,MaxLengths,Code,Traced);
            for (size_t j=0;j<Traced.size();j++)
            {
                if (Traced[j])
                    has_traced++;
                else
                    non_traced++;
            }
        }
        std::cout<<"Non Traced:"<<non_traced<<std::endl;
        std::cout<<"Traced:"<<has_traced<<std::endl;
    }

    //    static void retrievePath(TracingGraph<ScalarType> &TGraph,
    //                             const size_t endNode,
    //                             std::vector<size_t> &path)
    //    {
    //        assert(endNode>=0);
    //        assert(endNode<TGraph.NumNodes());
    //        // get the nodes from
    //        size_t curr_node = endNode;
    //        std::list<size_t> nodeList;
    //        do {
    //            nodeList.push_back(curr_node);
    //            curr_node = TGraph.Father(curr_node);
    //        } while ((curr_node != -1)&&(curr_node != endNode));

    //        path.assign(nodeList.rbegin(), nodeList.rend());
    //    }


    //    static bool ShortestPaths(TracingGraph<ScalarType> &TGraph,
    //                              const std::vector<int> &SourceNodes,
    //                              const std::vector<int> &DestNodes,
    //                              std::vector<int> &PathN,
    //                              bool avoid_zero_path=true)//,
    //                              //ScalarType &Len)
    //    {
    //        std::vector<HeapEntry> Heap;
    //        TGraph.UnMarkAll();
    //        for (size_t i=0;i<SourceNodes.size();i++)
    //        {
    //            size_t IndexN0=SourceNodes[i];
    //            TGraph.Mark(IndexN0);
    //            TGraph.Distance(IndexN0)=0;
    //            TGraph.Father(IndexN0)=-1;
    //            assert(TGraph.EnabledN[IndexN0]);
    //            Heap.push_back(HeapEntry(IndexN0,0));
    //        }

    //        for (size_t i=0;i<DestNodes.size();i++)
    //        {
    //            TGraph.MarkAsDest(DestNodes[i]);
    //            assert(TGraph.EnabledN[DestNodes[i]]);
    //        }

    //        //set initial
    //        std::make_heap(Heap.begin(),Heap.end());
    //        do
    //        {
    //            std::pop_heap(Heap.begin(),Heap.end());
    //            size_t CurrN=(Heap.back()).NodeI;
    //            ScalarType CurrWeight=TGraph.Distance(CurrN);
    //            Heap.pop_back();

    //            assert(TGraph.IsMarked(CurrN));


    //            bool HasCompletedPath=TGraph.IsDest(CurrN);
    //            //this is to also include loops
    //            if (avoid_zero_path)
    //                HasCompletedPath&=(TGraph.Distance(CurrN)>0);

    //            if (HasCompletedPath)
    //            {
    //                 RetrievePath(TGraph,CurrN,PathN);
    //            }
    //            for (size_t i=0;i<TGraph.N[CurrN].size();i++)
    //            {
    //                if (!TGraph.EnabledA[CurrN][i])continue;
    //                size_t NextN=TGraph.N[CurrN][i];
    //                if (!TGraph.EnabledN[NextN])continue;

    //                ScalarType NextWeight=CurrWeight+TGraph.W[CurrN][i];

    //                bool IsMarked=(TGraph.IsMarked(NextN));

    //                //check if has been explored or not
    //                ScalarType OldWeight=std::numeric_limits<ScalarType>::max();
    //                if (IsMarked)
    //                    OldWeight=TGraph.Distance(NextN);


    //                if (NextWeight<OldWeight)
    //                {
    //                    TGraph.Mark(NextN);
    //                    TGraph.Distance(NextN)=NextWeight;
    //                    TGraph.Father(NextN)=CurrN;

    //                    Heap.push_back(HeapEntry(NextN,NextWeight));
    //                    push_heap(Heap.begin(),Heap.end());
    //                }
    //            }

    //        }while (!Heap.empty());
    //        return false;
    //    }
};

#endif
