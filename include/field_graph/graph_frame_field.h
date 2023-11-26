#ifndef FRAME_FIELD_GRAPH
#define FRAME_FIELD_GRAPH

#include "../frame_field.h"

namespace Field{

class FrameFieldGraph
{

public:

    static int FrameOfNode(int IndexN)
    {
        return (IndexN%6);
    }

    static int PolyOfNode(int IndexN)
    {
        return (IndexN/6);
    }

    static int NodeOfPolyFrame(int IndexP,int IndexD)
    {
        assert(IndexD>=0);
        assert(IndexD<6);
        return (IndexP*6+IndexD);
        //return (IndexN%6);
    }

    static void NodesOfPoly(int IndexP,std::vector<int> &IndexN)
    {
        IndexN.clear();
        assert(IndexP>=0);
        for (size_t j=0;j<6;j++)
            IndexN.push_back(NodeOfPolyFrame(IndexP,j));
    }

    static void NodesOfSamePoly(int IndexN,std::vector<int> &IndexNodes)
    {
        IndexNodes.clear();
        assert(IndexN>=0);
        int PolyIDX=PolyOfNode(IndexN);
        NodesOfPoly(PolyIDX,IndexNodes);
    }


    static int OppositeNodeFrameField(int IndexN)
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

    static void OrthoNodesFrameField(int IndexN,std::vector<int> &OrthoN)
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

    template <class ScalarType>
    static bool AreDirectionCoherent(const std::vector<Field::FrameField<ScalarType> > &FrameField,
                                     const int &IndexN0,const int &IndexN1)
    {
        if (IndexN0==IndexN1)return true;
        assert(IndexN0!=IndexN1);
        assert(IndexN0>=0);
        assert(IndexN1>=0);

        int FrIndex0=FrameOfNode(IndexN0);
        assert(FrIndex0<FrameField.size());
        assert(FrIndex0>=0);
        int DrIndex0=FrameOfNode(IndexN0);
        assert(DrIndex0>=0);
        assert(DrIndex0<6);

        int FrIndex1=FrameOfNode(IndexN1);
        assert(FrIndex1<(int)FrameField.size());
        assert(FrIndex1>=0);
        int DrIndex1=FrameOfNode(IndexN1);
        assert(DrIndex1>=0);
        assert(DrIndex1<6);

        if((FrIndex0==FrIndex1)&&(FrIndex0!=FrIndex1))return false;//same frame different direction

        int DrTest1=FrameField[FrIndex1].FollowDirection(FrameField[FrIndex0],DrIndex0);

        return (DrTest1==DrIndex1);
    }
};
}
#endif
