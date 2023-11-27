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
    assert(OppDir<4);
    int ElemI=ElemOfNode(IndexN);
    return NodeOfElemDir(ElemI,OppDir);
}

static void OppositeNodeCrossField(const std::vector<int> &IndexN,
                                   std::vector<int> &OppIndexN)
{
    OppIndexN.clear();
    for (size_t i=0;i<IndexN.size();i++)
        OppIndexN.push_back(OppositeNodeCrossField(IndexN[i]));
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

};
}
#endif
