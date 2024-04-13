#ifndef GRAPH_TRACING_INTERFACE
#define GRAPH_TRACING_INTERFACE

#include "graph.h"
#include <list>
#include <set>
#include "path_functions.h"
//#include "../cross_field.h"
#include "graph_tracing.h"

namespace Geo{

template <class ScalarType>
class GraphTracing
{

    static bool CompleteTraceFrom(Geo::TracingGraph<ScalarType> &TGraph,
                                  const int &StartN,std::vector<int> &PathN,
                                  bool allow_loops,bool allow_exit_domains)
    {

        TGraph.DeSelectAll(0);

        assert(allow_loops || allow_exit_domains);

        //set destination borders if needed
        if (allow_exit_domains)
            TGraph.SelectNoteType(Geo::ExitDomain,0);

        //set destination loops if needed
        if (allow_loops)
            TGraph.Select(StartN,0);

        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
        SParam.avoid_zero_path=true;
        SParam.StopAtSel[0]=true;
        SParam.onlyDirect=false;
        SParam.stopAtPrevDist=false;

        ScalarType PathL;
        std::vector<int> SourceNodes;
        SourceNodes.push_back(StartN);
        bool traced=Geo::GraphQuery<ScalarType>::ShortestPaths(TGraph,SourceNodes,PathN,PathL,SParam);

        return (traced);

    }

    static bool IsOppositeExitDomain(const Geo::TracingGraph<ScalarType> &TGraph,
                                     const int &indexN)
    {
        int OppN=Geo::PathFunctions<ScalarType>::OppositeNode(TGraph,indexN);
        return(TGraph.NodeType(OppN)==(Geo::ExitDomain));
    }

    static bool IsOppositeTangent(const Geo::TracingGraph<ScalarType> &TGraph,
                                     const int &indexN)
    {
        int OppN=Geo::PathFunctions<ScalarType>::OppositeNode(TGraph,indexN);
        return(TGraph.NodeType(OppN)==(Geo::TangentDomain));
    }

public:


    static bool TraceLoopOrBorderBorderFrom(Geo::TracingGraph<ScalarType> &TGraph,
                                            const int &StartN,std::vector<int> &PathN,
                                            bool make_direct=true)
    {

        bool IsBorder=IsOppositeExitDomain(TGraph,StartN);
        IsBorder|=IsOppositeTangent(TGraph,StartN);
        bool has_traced=CompleteTraceFrom(TGraph,StartN,PathN,true,true);
        if (!has_traced)return false;
        bool IsLoop=Geo::PathFunctions<ScalarType>::IsLoop(PathN);

        bool has_to_complete=!(IsBorder || IsLoop);

        if (has_to_complete)
        {
            //std::cout<<"Completing"<<std::endl;
            int StartOpp=PathFunctions<ScalarType>::OppositeNode(TGraph,StartN);
            std::vector<int> PathNOpp;
            
            has_traced=CompleteTraceFrom(TGraph,StartOpp,PathNOpp,false,true);

            if (!has_traced)return false;

            //revert the path
            PathFunctions<ScalarType>::ReversePath(TGraph,PathNOpp);

            assert(PathNOpp.back()==PathN[0]);

            //remove the last one
            PathNOpp.pop_back();
            //then append
            PathNOpp.insert(PathNOpp.end(),PathN.begin(),PathN.end());
            PathN=PathNOpp;
        }

        if (make_direct)
        {
            bool has_direct=Geo::PathFunctions<ScalarType>::MakeDirectConnection(TGraph,PathN);
            if (!has_direct){return false;}
        }

        //FINAL CHECKS
        if (!IsLoop)
        {
            bool IsExit_0=TGraph.NodeType(PathN[0])==(Geo::ExitDomain);
            bool IsExitOpp_0=IsOppositeExitDomain(TGraph,PathN[0]);

            //in case joining with a dart or a concave
            bool IsTagent_0=TGraph.NodeType(PathN[0])==(Geo::TangentDomain);
            bool IsTagentOpp_0=IsOppositeTangent(TGraph,PathN[0]);
            assert(IsExit_0||IsExitOpp_0||IsTagent_0||IsTagentOpp_0);

            bool IsTagent_1=TGraph.NodeType(PathN[0])==(Geo::TangentDomain);
            bool IsExit_1=TGraph.NodeType(PathN.back())==(Geo::ExitDomain);

            bool IsExitOpp_1=IsOppositeExitDomain(TGraph,PathN.back());
            bool IsTagentOpp_1=IsOppositeTangent(TGraph,PathN.back());
            assert(IsExit_1||IsExitOpp_1||IsTagent_1||IsTagentOpp_1);
        }

        return true;
    }

};

}
#endif
