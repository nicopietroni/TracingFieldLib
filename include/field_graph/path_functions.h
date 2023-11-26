#ifndef PATH_FUNCTIONS
#define PATH_FUNCTIONS

#include "graph.h"
#include "graph_tracing.h"
#include <list>
#include <set>
#include <unordered_set>
#include "graph_cross_field.h"
#include "graph_frame_field.h"

namespace Geo{

template <class ScalarType>
class PathFunctions
{
public:

    static bool IsProperPath(const std::vector<int> &PathN)
    {
        for (size_t i=0;i<PathN.size()-1;i++)
            if (PathN[i]==PathN[(i+1)])return false;
        return true;
    }

    static bool IsLoop(const std::vector<int> &PathN)
    {return (PathN[0]==PathN.back());}

    static bool ConnectDirect(TracingGraph<ScalarType> &TGraph,
                              const int &IndexN0,
                              const int &IndexN1,
                              std::vector<int> &DirectPath)
    {
        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
        SParam.avoid_zero_path=true;
        SParam.StopAtSel[0]=true;

        TGraph.DeSelectAll(0);
        TGraph.Select(IndexN1,0);
        SParam.onlyDirect=true;
        SParam.onlyVisible=false;

        std::vector<int> SourceNodes(1,IndexN0);
        ScalarType PathL;
        bool has_connected=Geo::GraphQuery<ScalarType>::ShortestPaths(TGraph,SourceNodes,DirectPath,PathL,SParam);
        return has_connected;
    }

    static bool MakeDirectConnection(TracingGraph<ScalarType> &TGraph,
                                     std::vector<int> &PathN)
    {
        assert(IsProperPath(PathN));
        std::vector<int> NewPath;
        for (size_t i=0;i<PathN.size()-1;i++)
        {
            int IndexN0=PathN[i];
            int IndexN1=PathN[(i+1)];

            //add the first one
            NewPath.push_back(IndexN0);

            if (TGraph.IsDirectPairE(IndexN0,IndexN1))continue;
            std::vector<int> DirectPath;
            bool connected=ConnectDirect(TGraph,IndexN0,IndexN1,DirectPath);
            if (!connected)return false;

            //then from second to one before the last
            for (size_t j=1;j<DirectPath.size()-1;j++)
                NewPath.push_back(DirectPath[j]);

            assert(DirectPath[0]==IndexN0);
            assert(DirectPath.back()==IndexN1);
        }
        //then add the last
        NewPath.push_back(PathN.back());
        PathN=NewPath;
        return true;
    }

    static int OppositeNode(const Geo::TracingGraph<ScalarType> &TGraph,
                            const int &indexN)
    {
        //if the opposite is selected then is on the boundary
        Geo::GraphType GType=TGraph.TypeofGraph();
        int OppN=-1;
        if (GType==Geo::GT_CrossField)
            OppN=Field::CrossFieldGraph::OppositeNodeCrossField(indexN);
        if (GType==Geo::GT_FrameField)
            OppN=Field::FrameFieldGraph::OppositeNodeFrameField(indexN);
        assert(OppN!=-1);
        return OppN;
    }

    static void OppositeNodes(const TracingGraph<ScalarType> &TGraph,
                              const std::vector<int> &PathN,
                              std::vector<int> &OppPathN)
    {
        OppPathN.clear();
        for (size_t i=0;i<PathN.size();i++)
        {
            int IndexN=PathN[i];
            assert(IndexN>=0);
            assert(IndexN<TGraph.NumNodes());

            int OppN=OppositeNode(TGraph,IndexN);

            assert(OppN>=0);
            assert(OppN<(int)TGraph.NumNodes());

            OppPathN.push_back(OppN);
        }
    }

    static bool SelfIntersectTangentially(const TracingGraph<ScalarType> &TGraph,
                                          const std::vector<int> &PathN)
    {
        if (PathN.size()<2)return false;
        std::unordered_set<int> TestNSet;

        std::vector<int> OppPathN;
        OppositeNodes(TGraph,PathN,OppPathN);
        assert(OppPathN.size()==PathN.size());

        int LastTestN=PathN.size();
        if (IsLoop(PathN))LastTestN--;
        for (int i=0;i<LastTestN;i++)
        {
            int IndexN=PathN[i];
            int OppN=OppPathN[i];
            assert(IndexN>=0);
            assert(IndexN<(int)TGraph.NumNodes());
            assert(OppN>=0);
            assert(OppN<(int)TGraph.NumNodes());

            if (TestNSet.count(IndexN)>0)return true;
            if (TestNSet.count(OppN)>0)return true;
            TestNSet.insert(IndexN);
            TestNSet.insert(OppN);
        }
        return false;
    }

    static bool IntersectTangentially(const TracingGraph<ScalarType> &TGraph,
                                      const std::vector<int> &PathN0,
                                      const std::vector<int> &PathN1)
    {

        std::vector<int> OppPathN0;
        OppositeNodes(TGraph,PathN0,OppPathN0);
        assert(OppPathN0.size()==PathN0.size());

        std::unordered_set<int> TestNSet(PathN0.begin(),PathN0.end());
        TestNSet.insert(OppPathN0.begin(),OppPathN0.end());

        std::vector<int> OppPathN1;
        OppositeNodes(TGraph,PathN1,OppPathN1);
        assert(OppPathN1.size()==PathN1.size());

        int LastTestN=PathN1.size();
        if (IsLoop(PathN1))LastTestN--;
        for (int i=0;i<LastTestN;i++)
        {
            int IndexN=PathN1[i];
            int OppN=OppPathN1[i];
            assert(IndexN>=0);
            assert(IndexN<TGraph.NumNodes());
            assert(OppN>=0);
            assert(OppN<TGraph.NumNodes());

            if (TestNSet.count(IndexN)>0)return true;
            if (TestNSet.count(OppN)>0)return true;
        }
        return false;
    }


    static void ReversePath(const TracingGraph<ScalarType> &TGraph,
                            std::vector<int> &PathN)
    {
        std::vector<int> OppPathN;
        OppositeNodes(TGraph,PathN,OppPathN);
        std::reverse(OppPathN.begin(),OppPathN.end());
        PathN=OppPathN;
    }
};

}
#endif
