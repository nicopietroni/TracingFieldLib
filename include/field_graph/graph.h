#ifndef FIELD_GRAPH
#define FIELD_GRAPH

#include <cmath>
#include <vector>
#include <stdio.h>
#include <iostream>
//#include "../point3.h"

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
    size_t GlobalMark;

    //std::vector<bool> Start;
    //selected flags
    std::vector<bool> S0,S1;
    //source value for tracing short patsh
    std::vector<int> Src;

    GraphType GType;
    //END OF DIJKSTRA

    //    //the visited mark to be used for other algotithms
    //    std::vector<int> VisitedM;
    //    //the visited global mark
    //    size_t VisitedGlobalMark;

public:

    //set if a node is enabled or not
    std::vector<bool> EnabledN;
    //set if an arc is enabled or not
    std::vector<std::vector<bool> > EnabledA;
    //neighbouring information
    std::vector<std::vector<int> > N;
    //weight of the arch
    std::vector<std::vector<ScalarType> > W;
    //type of node
    std::vector<NType> NT;

    void SelectNoteType(NType TypeNode,int IndexS)
    {
        assert(NT.size()==S0.size());
        assert(NT.size()==S1.size());
        assert((IndexS==0)||(IndexS==1));
        for (size_t i=0;i<NT.size();i++)
        {
            if (NT[i]==TypeNode)
                Select(i,IndexS);
        }
    }

    void EnableAllNodes()
    {
        EnabledN=std::vector<bool>(N.size(),true);
        EnabledA=std::vector<std::vector<bool> > (N.size(),std::vector<bool>());
        for (size_t i=0;i<N.size();i++)
            EnabledA[i]=std::vector<bool>(N[i].size(),true);
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

    void DisableArconDisabledNodeType()
    {
        assert(N.size()==EnabledA.size());
        for (size_t i=0;i<N.size();i++)
        {
            assert(EnabledA[i].size()==N[i].size());
            for (size_t j=0;j<N[i].size();i++)
            {
                int IndexN=N[i][j];
                assert(IndexN>=0);
                assert(IndexN<EnabledN.size());
                if (EnabledN[IndexN])continue;
                EnabledA[i][j]=false;
            }
        }
    }

    void Clear()
    {
        EnabledN.clear();
        EnabledA.clear();
        N.clear();
        W.clear();
        TM.clear();
        F.clear();
        D.clear();
        S0.clear();
        S1.clear();
        NT.clear();
        Src.clear();
        //VisitedM.clear();

        GlobalMark=0;
        //VisitedGlobalMark=0;
    }

    bool IsSink(const size_t &IndexN)const
    {
        assert(IndexN<N.size());
        assert(IndexN<W.size());
        assert(N[IndexN].size()==W[IndexN].size());
        return(N[IndexN].size()==0);
    }

    bool AreConnected(const int &IndexN0,
                      const int &IndexN1)const
    {
        assert(EnabledN[IndexN0]);
        assert(EnabledN[IndexN1]);
        for (size_t i=0;i<N[IndexN0].size();i++)
        {
            if (!EnabledA[IndexN0][i])continue;
            if (N[IndexN0][i]==IndexN1)return true;
        }
        return false;
    }

    void RemoveNeighbours(const int &IndexN)
    {
        assert(IndexN>=0);
        assert(IndexN<(int)N.size());
        assert(IndexN<(int)W.size());
        N[IndexN].clear();
        W[IndexN].clear();
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

    //    void SetAllUnVisited()
    //    {
    //        VisitedGlobalMark++;
    //    }

    //    void SetVisited(const size_t &IndexN)
    //    {
    //        assert(IndexN<VisitedM.size());
    //        VisitedM[IndexN]=VisitedGlobalMark;
    //    }

    //    bool IsVisited(const size_t &IndexN)
    //    {
    //        assert(IndexN<VisitedM.size());
    //        return(VisitedM[IndexN]==VisitedGlobalMark);
    //    }


    void Select(const size_t &IndexN,const int IndexS)
    {
        assert((IndexS==0)||(IndexS==1));
        assert(IndexN<S0.size());
        assert(IndexN<S1.size());
        if (IndexS==0)
            S0[IndexN]=true;
        else
            S1[IndexN]=true;
    }

    void Unselect(const size_t &IndexN,const int IndexS)
    {
        assert((IndexS==0)||(IndexS==1));
        assert(IndexN<S0.size());
        assert(IndexN<S1.size());
        if (IndexS==0)
            S0[IndexN]=false;
        else
            S1[IndexN]=false;
    }

    bool IsSelected(const size_t &IndexN,const int IndexS)const
    {
        assert(IndexN<S0.size());
        assert(IndexN<S1.size());
        if (IndexS==0)
            return(S0[IndexN]);
        else
            return(S1[IndexN]);
    }

    void InitTracing()
    {
        TM.clear();
        TM.resize(N.size(),0);
        F.clear();
        F.resize(N.size(),-1);
        D.clear();
        D.resize(N.size(),0);
        S0.clear();
        S0.resize(N.size(),false);
        S1.clear();
        S1.resize(N.size(),false);
        Src.clear();
        Src.resize(N.size(),-1);
        GlobalMark=0;
        //        VisitedM.clear();
        //        VisitedM.resize(N.size(),0);

        GlobalMark=0;
        //        VisitedGlobalMark=0;
        EnableAllNodes();
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
        assert(EnabledN.size()==EnabledA.size());
        assert(N.size()==W.size());
        assert(TM.size()==F.size());
        assert(F.size()==D.size());
        return (N.size());
    }

    void SortNeighbyW(const size_t &IndexN)
    {
        assert(N.size()==W.size());
        assert(IndexN<N.size());
        assert(IndexN<W.size());
        std::vector<std::pair<ScalarType,size_t> > WN;
        for (size_t i=0;i<N[IndexN].size();i++)
            WN.push_back(std::pair<ScalarType,size_t>(W[IndexN][i],N[IndexN][i]));

        std::sort(WN.begin(),WN.end());

        W[IndexN].clear();
        N[IndexN].clear();

        for (size_t i=0;i<WN.size();i++)
        {
            W[IndexN].push_back(WN[i].first);
            N[IndexN].push_back(WN[i].second);
        }

    }

    void SortAllNeighbyW()
    {
        for (size_t i=0;i<NumNodes();i++)
            SortNeighbyW(i);
    }

    void SelectAll(const int IndexS)
    {
        assert((IndexS==0)||(IndexS==1));
        if (IndexS==0)
            std::fill(S0.begin(), S0.end(), true);
        else
            std::fill(S1.begin(), S1.end(), true);
    }

    void DeSelectAll(const int IndexS)
    {
        assert((IndexS==0)||(IndexS==1));
        if (IndexS==0)
            std::fill(S0.begin(), S0.end(), false);
        else
            std::fill(S1.begin(), S1.end(), false);
    }

    GraphType TypeofGraph()const
    {return GType;}

    int FirstActiveNeigh(const size_t &IndexN)
    {
        //then get the closest node
        for (size_t i=0;i<N[IndexN].size();i++)
        {
            int NextN=N[IndexN][i];
            if (!EnabledN[NextN])continue;
            return NextN;
        }
        return -1;
    }

    //        void Check()
    //        {
    //            ScalarType L=0;
    //            for (size_t i=0;i<W.size();i++)
    //                for (size_t j=0;j<W[i].size();j++)
    //                {
    //                   L+=W[i][j];
    //                }
    //            std::cout<<"L"<<L<<std::endl;
    //        }

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
