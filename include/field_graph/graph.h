#ifndef FIELD_GRAPH
#define FIELD_GRAPH

#include <cmath>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <unordered_set>
#include <set>

//#include "../point3.h"

#define SELECTED_GRAPH_FLAG 4

namespace  Geo{

enum NType{
    ExitDomain, TangentDomain, Internal
};

enum GraphType
{
    GT_CrossField,GT_FrameField
};


template <class ScalarType>
class TracingGraph
{


    //TO BE USED IN DIJKSTRA
    //the temporary mark
    std::vector<int> TM;

    //array of fathers node
    std::vector<int> F;

    //the distance
    std::vector<ScalarType> D;

    //the global mark
    int GlobalMark;

    //selected marks used for Dijkstra
    std::vector<int> Sel[SELECTED_GRAPH_FLAG];
    int GlobalSelMark[SELECTED_GRAPH_FLAG];

    //source value for tracing short patsh
    std::vector<int> Src;
    //END OF DIJKSTRA

    GraphType GType;

    // Hash function
    struct hashPaitIntFunction
    {
      size_t operator()(const std::pair<int , int> &x) const
      {
        return x.first ^ x.second;
      }
    };

    std::unordered_set<std::pair<int,int>,hashPaitIntFunction > DirectE;
    //std::set<std::pair<int,int> > DirectE;

    bool AllocatedSelected()
    {
        for (size_t i=0;i<SELECTED_GRAPH_FLAG;i++)
            if(NT.size()!=Sel[i].size())return false;
        return true;
    }

    struct NeighInfo
    {
        //the next node
        int N;
        //the weight
        ScalarType W;
        //active or not
        bool Enabled;
        //direct or not
        bool Direct;
        //Visible or not
        bool Visible;

        NeighInfo()
        {
            N=-1;
            W=-1;
            Enabled=false;
            Direct=false;
            Visible=false;
        }

        inline bool operator ==(const NeighInfo &N1)const
        {
            return ((N==N1.N)&&(W==N1.W)&&
                    (Enabled==N1.Enabled)&&
                    (Direct==N1.Direct)&&
                    (Visible==N1.Visible));
        }

        inline bool operator <(const NeighInfo &N1)const
        {
            //assert(N!=N1.N);//can't have two times same neighbour
            if (Visible && !N1.Visible)return true;
            if (!Visible && N1.Visible)return false;
            if (W!=N1.W)return (W<N1.W);
            return (N<N1.N);
        }

        NeighInfo(const int _N,
                  const ScalarType _W,
                  const bool _Enabled,
                  const bool _Direct,
                  const bool _Visible)
        {
            N=_N;
            W=_W;
            Enabled=_Enabled;
            Direct=_Direct;
            Visible=_Visible;
        }
    };


    std::vector<std::vector<NeighInfo> > N;

public:

    void InitDirectEdge()
    {
        DirectE.clear();

        for (size_t i=0;i<NumNodes();i++)
            for (size_t j=0;j<N[i].size();j++)
            {
                int IndexN0=i;
                if (!DirectLink(i,j))continue;
                int IndexN1=NodeLink(i,j);
                std::pair<int,int> Key(std::min(IndexN0,IndexN1),
                                       std::max(IndexN1,IndexN0));
                DirectE.insert(Key);
            }
        std::cout<<"Direct Link "<<DirectE.size()<<std::endl;
    }

    //type of node
    std::vector<NType> NT;
    //set if a node is enabled or not
    std::vector<bool> EnabledN;

    bool IsDirectPairE(const int &IndexN0,
                       const int &IndexN1)
    {
        std::pair<int,int> Key(std::min(IndexN0,IndexN1),
                               std::max(IndexN0,IndexN1));
        return(DirectE.count(Key)>0);
    }

    bool ProperNeigh(const int &IndexN)
    {
        if (N[IndexN].size()<2)
            return true;

        for (size_t i=0;i<N[IndexN].size()-1;i++)
        {
            for (size_t j=(i+1);j<N[IndexN].size();j++)
                if (N[IndexN][i].N==N[IndexN][j].N)return false;
        }
        return true;
    }

    bool ProperNeigh()
    {
        for (size_t i=0;i<NumNodes();i++)
            if (!ProperNeigh(i))return false;
        return true;
    }

    void AddNeighbour(const int &IndexN0,
                      const int &IndexN1,
                      const ScalarType W,
                      const bool Enabled,
                      const bool Direct,
                      const bool Visible)
    {
        assert(IndexN0>=0);
        assert(IndexN0<(int)N.size());
        N[IndexN0].push_back(NeighInfo(IndexN1,W,Enabled,Direct,Visible));
    }

    void AllocateNodes(int NumNodes)
    {
        N.clear();
        N.resize(NumNodes);
        NT.clear();
        NT.resize(NumNodes,Internal);
        EnabledN.clear();
        EnabledN.resize(NumNodes,true);
    }

    int NeighSize(const int &IndexN0)
    {
        assert(IndexN0>=0);
        assert(IndexN0<(int)N.size());
        return ((int)N[IndexN0].size());
    }

    bool ActiveLink(const int &IndexN0,const int &LocalIdx)const
    {
        assert(IndexN0>=0);
        assert(IndexN0<(int)N.size());
        assert(LocalIdx>=0);
        assert(LocalIdx<(int)N[IndexN0].size());
        return (N[IndexN0][LocalIdx].Enabled);
    }

    bool DirectLink(const int &IndexN0,const int &LocalIdx)const
    {
        assert(IndexN0>=0);
        assert(IndexN0<(int)N.size());
        assert(LocalIdx>=0);
        assert(LocalIdx<(int)N[IndexN0].size());
        return (N[IndexN0][LocalIdx].Direct);
    }

    bool VisibleLink(const int &IndexN0,const int &LocalIdx)const
    {
        assert(IndexN0>=0);
        assert(IndexN0<(int)N.size());
        assert(LocalIdx>=0);
        assert(LocalIdx<(int)N[IndexN0].size());
        return (N[IndexN0][LocalIdx].Visible);
    }

    ScalarType WeightLink(const int &IndexN0,const int &LocalIdx)const
    {
        assert(IndexN0>=0);
        assert(IndexN0<(int)N.size());
        assert(LocalIdx>=0);
        assert(LocalIdx<(int)N[IndexN0].size());
        return (N[IndexN0][LocalIdx].W);
    }

    int NodeLink(const int &IndexN0,const int &LocalIdx)const
    {
        assert(IndexN0>=0);
        assert(IndexN0<(int)N.size());
        assert(LocalIdx>=0);
        assert(LocalIdx<(int)N[IndexN0].size());
        return (N[IndexN0][LocalIdx].N);
    }

    void SelectNoteType(NType TypeNode,int IndexS)
    {
        assert(AllocatedSelected());
        for (size_t i=0;i<NT.size();i++)
        {
            if (NT[i]==TypeNode)
                Select(i,IndexS);
        }
    }

    void EnableAllNodes()
    {
        EnabledN=std::vector<bool>(N.size(),true);
    }

    void EnableAllArcs()
    {
        for (size_t i=0;i<N.size();i++)
            for (size_t j=0;j<N[i].size();j++)
                N[i][j].Enabled=true;
    }

    void DisableNodeType(NType TypeNode)
    {
        assert(N.size()==NT.size());
        for (size_t i=0;i<NT.size();i++)
        {
            if (NT[i]==TypeNode)
                EnabledN[i]=false;
        }
    }

    void EnableNode(const int &IndexN)
    {
        assert(IndexN>=0);
        assert(IndexN<EnabledN.size());
        EnabledN[IndexN]=true;
    }

    void DisableNode(const int &IndexN)
    {
        assert(IndexN>=0);
        assert(IndexN<EnabledN.size());
        EnabledN[IndexN]=false;
    }

    void DisableArcsToDisabledNodes()
    {
        for (size_t i=0;i<N.size();i++)
        {
            for (size_t j=0;j<N[i].size();i++)
            {
                int IndexN=N[i][j];
                assert(IndexN>=0);
                assert(IndexN<EnabledN.size());
                if (EnabledN[IndexN])continue;
                N[i][j].Enabled=false;
            }
        }
    }

    void Clear()
    {
        EnabledN.clear();
        //EnabledA.clear();
        N.clear();

        //W.clear();
        TM.clear();
        F.clear();
        D.clear();
        NT.clear();
        Src.clear();
        GlobalMark=0;
        for (size_t i=0;i<SELECTED_GRAPH_FLAG;i++)
        {
            Sel[i].clear();
            GlobalSelMark[i]=0;
        }
    }

    bool IsSink(const size_t &IndexN)const
    {
        assert(IndexN<N.size());
        return(N[IndexN].size()==0);
    }

    bool AreConnected(const int &IndexN0,
                      const int &IndexN1)const
    {
        assert(EnabledN[IndexN0]);
        assert(EnabledN[IndexN1]);
        for (size_t i=0;i<N[IndexN0].size();i++)
        {
            if (!N[IndexN0][i].Enabled)continue;
            if (N[IndexN0][i].N==IndexN1)return true;
        }
        return false;
    }

    void RemoveNeighbours(const int &IndexN)
    {
        assert(IndexN>=0);
        assert(IndexN<(int)N.size());
        //assert(IndexN<(int)W.size());
        N[IndexN].clear();
        //W[IndexN].clear();
    }

    void UnMarkAll()
    {
        GlobalMark++;
    }

    void Mark(const size_t &IndexN)
    {
        assert(IndexN<TM.size());
        TM[IndexN]=GlobalMark;
    }

    bool IsMarked(const size_t &IndexN)const
    {
        assert(IndexN<TM.size());
        return(TM[IndexN]==GlobalMark);
    }

    void Select(const size_t &IndexN,const int IndexS)
    {
        assert((IndexS>=0)||(IndexS<=SELECTED_GRAPH_FLAG));
        assert(IndexN<Sel[IndexS].size());
        Sel[IndexS][IndexN]=GlobalSelMark[IndexS];
    }

    bool IsSelected(const size_t &IndexN,const int IndexS)const
    {
        assert((IndexS>=0)||(IndexS<=SELECTED_GRAPH_FLAG));
        assert(IndexN<Sel[IndexS].size());
        return (Sel[IndexS][IndexN]==GlobalSelMark[IndexS]);
    }

    void InitTracing()
    {
        TM.clear();
        TM.resize(N.size(),0);
        F.clear();
        F.resize(N.size(),-1);
        D.clear();
        D.resize(N.size(),0);

        for (size_t i=0;i<SELECTED_GRAPH_FLAG;i++)
        {
            GlobalSelMark[i]=0;
            Sel[i].resize(N.size(),-1);
        }

        Src.resize(N.size(),-1);
        GlobalMark=0;

        EnableAllNodes();

        InitDirectEdge();
    }

    int &Father(const size_t &IndexN)
    {
        assert(IndexN<F.size());
        return F[IndexN];
    }

    const int &Father(const size_t &IndexN)const
    {
        assert(IndexN<F.size());
        return F[IndexN];
    }

    int &Source(const size_t &IndexN)
    {
        assert(IndexN<Src.size());
        return Src[IndexN];
    }

    const int &Source(const size_t &IndexN)const
    {
        assert(IndexN<Src.size());
        return Src[IndexN];
    }

    const ScalarType &Distance(const size_t &IndexN)const
    {
        assert(IndexN<D.size());
        return D[IndexN];
    }

    ScalarType &Distance(const size_t &IndexN)
    {
        assert(IndexN<D.size());
        return D[IndexN];
    }

    const NType &NodeType(const size_t &IndexN)const
    {
        assert(IndexN<NT.size());
        return NT[IndexN];
    }

    NType &NodeType(const size_t &IndexN)
    {
        assert(IndexN<NT.size());
        return NT[IndexN];
    }

    size_t NumNodes()const
    {
        assert(TM.size()==F.size());
        assert(F.size()==D.size());
        return (N.size());
    }

    void SortNeigh(const size_t &IndexN)
    {
        assert(IndexN<N.size());
        std::sort(N[IndexN].begin(),N[IndexN].end());
    }


    void SortAllNeigh()
    {
        for (size_t i=0;i<NumNodes();i++)
            SortNeigh(i);
    }

    void DeSelectAll(const int IndexS)
    {
        assert((IndexS>=0)||(IndexS<=SELECTED_GRAPH_FLAG));
        GlobalSelMark[IndexS]++;
    }


    GraphType TypeofGraph()const
    {return GType;}

    int FirstActiveNeigh(const size_t &IndexN)
    {
        //then get the closest node
        for (size_t i=0;i<N[IndexN].size();i++)
        {
            //int NextN=N[IndexN][i];
            if (!ActiveLink(IndexN,i))continue;
            int NextN=NodeLink(IndexN,i);
            if (!EnabledN[NextN])continue;
            return NextN;
        }
        return -1;
    }

    TracingGraph(GraphType _GType)
    {GType=_GType;}
};

template <class ScalarType>
void FilterNodewOfType(const TracingGraph<ScalarType> &Tgraph,
                       const NType &TypeN,std::vector<int> &Nodes)
{
    std::vector<int> NodesTemp;
    for (size_t i=0;i<Nodes.size();i++)
    {
        int currN=Nodes[i];
        if (Tgraph.NodeType(currN)!=TypeN)continue;
        NodesTemp.push_back(currN);
    }
    Nodes=NodesTemp;
}

template <class ScalarType>
void FilterNonActiveNodes(const TracingGraph<ScalarType> &Tgraph,
                          std::vector<int> &Nodes)
{
    std::vector<int> NodesTemp;
    for (size_t i=0;i<Nodes.size();i++)
    {
        int currN=Nodes[i];
        if (!Tgraph.EnabledN[currN])continue;
        NodesTemp.push_back(currN);
    }
    Nodes=NodesTemp;
}

}
#endif
