#ifndef FRAME_FIELD_GRAPH
#define FRAME_FIELD_GRAPH

//#include "graph.h"
#include "../frame_field.h"

int FrameOfNode(int IndexN)
{
    return (IndexN%6);
}

int PolyOfNode(int IndexN)
{
    return (IndexN/6);
}

int NodeOfPolyFrame(int IndexP,int IndexD)
{
    assert(IndexD>=0);
    assert(IndexD<6);
    return (IndexP*6+IndexD);
    //return (IndexN%6);
}

void NodesOfPoly(int IndexP,std::vector<int> &IndexN)
{
    IndexN.clear();
    assert(IndexP>=0);
    for (size_t j=0;j<6;j++)
        IndexN.push_back(NodeOfPolyFrame(IndexP,j));
}

void NodesOfSamePoly(int IndexN,std::vector<int> &IndexNodes)
{
    IndexNodes.clear();
    assert(IndexN>=0);
    int PolyIDX=PolyOfNode(IndexN);
    NodesOfPoly(PolyIDX,IndexNodes);
}


int OppositeNodeFrameField(int IndexN)
{
    int FrameI=FrameOfNode(IndexN);
    assert(FrameI>=0);
    assert(FrameI<6);
    int OppDir=Field::OppositeFrameDir(FrameI);
    assert(OppDir>=0);
    assert(OppDir<6);
    int PolyI=PolyOfNode(IndexN);
    return NodeOfPolyFrame(PolyI,OppDir);
}

void OrthoNodesFrameField(int IndexN,std::vector<int> &OrthoN)
{
    int FrameI=FrameOfNode(IndexN);
    assert(FrameI>=0);
    assert(FrameI<6);
    std::vector<int> OrthoD;
    Field::OrthoFrameDirs(FrameI,OrthoD);
    OrthoN.clear();
    for (size_t i=0;i<OrthoD.size();i++)
    {
        assert(OrthoD[i]>=0);
        assert(OrthoD[i]<6);
        int PolyI=PolyOfNode(IndexN);
        int NodeI=NodeOfPolyFrame(PolyI,OrthoD[i]);
        OrthoN.push_back(NodeI);
    }
}

//template <class ScalarType>
//void SelectPossibleArrivalReturnNodes(const size_t &startNode,
//                                      Geo::TracingGraph<ScalarType> &TGraph,
//                                      int IndexSel=0)
//{
//    std::vector<int> OrthoN;
//    OrthoNodesFrameField((int)startNode,OrthoN);

//    TGraph.Select(startNode,IndexSel);
//    for (size_t j=0;j<OrthoN.size();j++)
//        TGraph.Select(OrthoN[j],IndexSel);
//}

//static void UnSelectPossibleArrivalReturnNodes(const size_t &startNode,
//                                         TracingGraph<ScalarType> &TGraph,
//                                         int IndexSel=0)
//{
//    std::vector<int> OrthoN;
//    OrthoNodesFrameField((int)startNode,OrthoN);

//    TGraph.Unselect(startNode,IndexSel);
//    for (size_t j=0;j<OrthoN.size();j++)
//        TGraph.Unselect(OrthoN[j],IndexSel);
//}

#endif
