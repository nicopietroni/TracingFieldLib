#ifndef CROSS_FIELD_GRAPH
#define CROSS_FIELD_GRAPH

#include "../cross_field.h"

namespace Field{

class CrossFieldGraph
{
public:

static int DirOfNode(int IndexN)
{
    return (IndexN%4);
}

static int ElemOfNode(int IndexN)
{
    return (IndexN/4);
}

static int NodeOfElemDir(int IndexE,int IndexD)
{
    assert(IndexD>=0);
    assert(IndexD<4);
    return (IndexE*4+IndexD);
}

static void NodesOfElem(int IndexE,std::vector<int> &IndexN)
{
    IndexN.clear();
    assert(IndexE>=0);
    for (size_t j=0;j<4;j++)
        IndexN.push_back(NodeOfElemDir(IndexE,j));
}

static void NodesOfSameElem(int IndexN,std::vector<int> &IndexNodes)
{
    IndexNodes.clear();
    assert(IndexN>=0);
    int ElemIDX=ElemOfNode(IndexN);
    NodesOfElem(ElemIDX,IndexNodes);
}

static int OppositeNodeCrossField(int IndexN)
{
    int DirI=DirOfNode(IndexN);
    assert(DirI>=0);
    assert(DirI<4);
    int OppDir=Field::OppositeCrossDir(DirI);
    assert(OppDir>=0);
    assert(OppDir<6);
    int ElemI=ElemOfNode(IndexN);
    return NodeOfElemDir(ElemI,OppDir);
}

void OrthoNodesCrossField(int IndexN,std::vector<int> &OrthoN)
{
    int DirI=DirOfNode(IndexN);
    assert(DirI>=0);
    assert(DirI<4);
    std::vector<int> OrthoD;
    Field::OrthoCrossDirs(DirI,OrthoD);
    OrthoN.clear();
    for (size_t i=0;i<OrthoD.size();i++)
    {
        assert(OrthoD[i]>=0);
        assert(OrthoD[i]<4);
        int ElemI=ElemOfNode(IndexN);
        int NodeI=NodeOfElemDir(ElemI,OrthoD[i]);
        OrthoN.push_back(NodeI);
    }
}

//template <class ScalarType>
//bool AreDirectionCoherent(const std::vector<Field::FrameField<ScalarType> > &FrameField,
//                          const int &IndexN0,const int &IndexN1)
//{
//    if (IndexN0==IndexN1)return true;
//    assert(IndexN0!=IndexN1);
//    assert(IndexN0>=0);
//    assert(IndexN1>=0);

//    int FrIndex0=FrameOfNode(IndexN0);
//    assert(FrIndex0<FrameField.size());
//    assert(FrIndex0>=0);
//    int DrIndex0=FrameOfNode(IndexN0);
//    assert(DrIndex0>=0);
//    assert(DrIndex0<6);

//    int FrIndex1=FrameOfNode(IndexN1);
//    assert(FrIndex1<FrameField.size());
//    assert(FrIndex1>=0);
//    int DrIndex1=FrameOfNode(IndexN1);
//    assert(DrIndex1>=0);
//    assert(DrIndex1<6);

//    if((FrIndex0==FrIndex1)&&(FrIndex0!=FrIndex1))return false;//same frame different direction

//    int DrTest1=FrameField[FrIndex1].FollowDirection(FrameField[FrIndex0],DrIndex0);

//    return (DrTest1==DrIndex1);
//}

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
};
}
#endif
