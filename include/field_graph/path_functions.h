#ifndef PATH_FUNCTIONS
#define PATH_FUNCTIONS

#include "graph.h"
#include "graph_tracing.h"
#include <list>
#include <set>
#include <unordered_set>
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "../math_utils.h"

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
        SParam.MaxSteps=TGraph.diffuse_step+1;
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
            assert(TGraph.IsEnabled(IndexN0));
            assert(TGraph.IsEnabled(IndexN1));
            bool connected=ConnectDirect(TGraph,IndexN0,IndexN1,DirectPath);
            if (!connected)return false;

            //the new path can be only 3 + diffusion steps (should be 2 but one because of singularities)
            if (DirectPath.size()>(TGraph.diffuse_step+3))
                return false;
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

    static bool HasSelfEdgeTwice(const std::vector<int> &PathN)
    {
        std::unordered_set<std::pair<int,int>,Geo::hashPairIntFunction > TestNSet;
        for (int i=0;i<PathN.size()-1;i++)
        {
            int IndexN0=PathN[i];
            int IndexN1=PathN[i+1];
            int IndexV0=Field::CrossFieldGraph::ElemOfNode(IndexN0);
            int IndexV1=Field::CrossFieldGraph::ElemOfNode(IndexN1);
            std::pair<int,int> key(std::min(IndexV0,IndexV1),
                                   std::max(IndexV0,IndexV1));
            if (TestNSet.count(key)>0)return true;
            TestNSet.insert(key);
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

    static bool BestConnectionFrom(const TracingGraph<ScalarType> &TGraph,
                                   const int &Node0,const int &VertIdx1,
                                   int &NodeIdx1,ScalarType &bestW)
    {
        NodeIdx1=-1;

        bool found=false;

        bestW=std::numeric_limits<ScalarType>::max();

        std::vector<int> Neigh;
        std::vector<ScalarType> NeighW;
        TGraph.GetNeigh(Node0,true,true,Neigh,NeighW);

        for (size_t j=0;j<Neigh.size();j++)
        {
            int NeighV=-1;

            if (TGraph.TypeofGraph()==GT_CrossField)
                NeighV=Field::CrossFieldGraph::ElemOfNode(Neigh[j]);
            else
                NeighV=Field::CrossFieldGraph::ElemOfNode(Neigh[j]);

            if (NeighV!=VertIdx1)continue;

            found=true;

            if (NeighW[j]>bestW)continue;

            bestW=NeighW[j];
            NodeIdx1=Neigh[j];
        }

        return found;
    }

    static bool BestConnectionNodePair(const TracingGraph<ScalarType> &TGraph,
                                       const int &VertIdx0,
                                       const int &VertIdx1,
                                       int &NodeIdx0,
                                       int &NodeIdx1)
    {
        NodeIdx0=-1;
        NodeIdx1=-1;

        assert(VertIdx0!=VertIdx1);
        std::vector<int> IndexN0;

        //if (TGraph.TypeofGraph()==GT_CrossField)
            Field::CrossFieldGraph::NodesOfElem(VertIdx0,IndexN0);
//        else
//            Field::FrameFieldGraph::NodesOfElem(VertIdx0,IndexN0);

        ScalarType bestW=std::numeric_limits<ScalarType>::max();
        bool found=false;
        for (size_t i=0;i<IndexN0.size();i++)
        {
            int TestNodeIdx1;
            ScalarType TestbestW;
            bool found_this=BestConnectionFrom(TGraph,IndexN0[i],VertIdx1,
                                               TestNodeIdx1,TestbestW);
            if (!found_this)continue;

            if (TestbestW>bestW)continue;

            found=true;
            bestW=TestbestW;
            NodeIdx0=IndexN0[i];
            NodeIdx1=TestNodeIdx1;
        }
        return found;
    }

    //should work generic but to be tested
    static bool ReconstructPathFromElem(const TracingGraph<ScalarType> &Tgraph,
                                             const std::vector<int> &Elems,
                                             std::vector<int> &Nodes)
    {
        Nodes.clear();
        assert(Elems.size()>=2);

        int E0=Elems[0];
        int E1=Elems[1];

        int NodeIdx0,NodeIdx1;
        bool found=BestConnectionNodePair(Tgraph,E0,E1,NodeIdx0,NodeIdx1);
        if (!found){return false;}

        Nodes.push_back(NodeIdx0);
        for (size_t i=1;i<Elems.size();i++)
        {
            int NodeIdx1;
            ScalarType bestW;
            bool found=BestConnectionFrom(Tgraph,Nodes.back(),
                                          Elems[i],NodeIdx1,
                                          bestW);
            if (!found)return true;
            Nodes.push_back(NodeIdx1);
        }
        return true;
    }
};

}
#endif
