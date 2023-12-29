#ifndef POS3
#define POS3

#include <vector>
#include <map>
#include <iostream>
#include <connectivity.h>

class Pos3{

public:

    int F,E,V;

    inline bool operator == ( Pos3 const & p ) const {
        return (V==p.V && E==p.E && F==p.F);
    }

    inline bool operator != ( Pos3 const & p ) const {
        return (V!=p.V || E!=p.E || F!=p.F);
    }

    inline bool operator <= ( Pos3 const & p) const {
        return	(F!=p.F)?(F<p.F):
                         (E!=p.E)?(E<p.E):
                                  (V<=p.V);
    }

    Pos3()
    {
        V=-1;
        E=-1;
        F=-1;
    }

    Pos3(int _F,int _E, int _V)
    {
        V=_V;
        E=_E;
        F=_F;
        assert(V>=0);
        assert(E>=0);
        assert(F>=0);
    }

    void FlipV(const std::vector<std::vector<int> > &Faces)
    {
        assert(F<(int)Faces.size());
        assert(E<(int)Faces[F].size());
        int IndexE0=E;
        int IndexE1=(E+1)%Faces[F].size();
        int currV0=Faces[F][IndexE0];
        int currV1=Faces[F][IndexE1];
        assert((V==currV0)||(V==currV1));
        if (V==currV0)
            V=currV1;
        else
            V=currV0;
    }

    int VFlip(const std::vector<std::vector<int> > &Faces)
    {
        Pos3 TestP=(*this);
        TestP.FlipV(Faces);
        return TestP.V;
    }

    void FlipE(const std::vector<std::vector<int> > &Faces)
    {
        assert(F<(int)Faces.size());
        assert(E<(int)Faces[F].size());
        size_t sizeF=Faces[F].size();

        if(Faces[F][E]==V)
        {
            E=(E+sizeF-1)%sizeF;
        }else
        {
            E=(E+1)%sizeF;
        }
    }

    int CheckConsistency(const std::vector<std::vector<int> > &Faces,
                         const std::vector<std::vector<int> > &NextF,
                         const std::vector<std::vector<int> > &NextEd)const
    {
        if (F<0)return -1;
        if (E<0)return -2;
        if (F>=(int)Faces.size())return -3;
        if (E>=(int)Faces[F].size())return -4;

        size_t sizeF=Faces[F].size();
        if (!((Faces[F][E]==V)||(Faces[F][(E+1)%sizeF]==V)))return -5;

        int newF=NextF[F][E];
        if (newF!=-1)
        {
            int newE=NextEd[F][E];
            if (NextF[newF][newE]!=F)return -6;
        }
        return 0;
    }

    void FlipF(const std::vector<std::vector<int> > &NextF,
               const std::vector<std::vector<int> > &NextEd)
    {
        assert(F>=0);
        assert(E>=0);
        assert(F<(int)NextF.size());
        assert(E<(int)NextF[F].size());
        assert(F<(int)NextEd.size());
        assert(E<(int)NextEd[F].size());
        int newF=NextF[F][E];
        int newE=NextEd[F][E];
        F=newF;
        E=newE;
    }

    bool IsBorder(const std::vector<std::vector<int> > &NextF)
    {
        int newF=NextF[F][E];
        return (newF==-1);
    }

    void NextE(const std::vector<std::vector<int> > &Faces,
               const std::vector<std::vector<int> > &NextF,
               const std::vector<std::vector<int> > &NextEd)
    {
        FlipF(NextF,NextEd);
        FlipE(Faces);
        assert(F>=0);
        assert(E>=0);
        assert(F<(int)NextF.size());
        assert(E<(int)NextF[F].size());
    }

    void NextB(const std::vector<std::vector<int> > &Faces,
               const std::vector<std::vector<int> > &NextF,
               const std::vector<std::vector<int> > &NextEd)
    {
        assert(IsBorder(NextF));
        FlipE(Faces);
        while (!IsBorder(NextF))
        {
            NextE(Faces,NextF,NextEd);
        }
    }
};


void GetSortedFacesFaceFun(const Pos3 &StartPos,
                           const std::vector<std::vector<int> > &Faces,
                           const std::vector<std::vector<int> > &NextF,
                           const std::vector<std::vector<int> > &NextE,
                           std::vector<Pos3 > &PosFan,
                           bool &HasBorder)
{
    Pos3 currPos=StartPos;
    HasBorder=false;

    assert(currPos.CheckConsistency(Faces,NextF,NextE)==0);
    //    if (ret!=0)
    //    {
    //        std::cout<<"Error:"<<ret<<std::endl;
    //        assert(0);
    //    }

    do
    {
        PosFan.push_back(currPos);
        //if met the border the first time stop
        if (currPos.IsBorder(NextF))
            HasBorder=true;
        else
            currPos.NextE(Faces,NextF,NextE);
    }while ((currPos!=StartPos)&&(!HasBorder));
    assert(currPos.V==StartPos.V);

    if (!HasBorder)return;

    //otherwise search for next border
    PosFan.clear();
    currPos.FlipE(Faces);
    HasBorder=false;
    do
    {
        PosFan.push_back(currPos);

        //if met the border the first time stop
        if (currPos.IsBorder(NextF))
            HasBorder=true;
        else
        {
            currPos.NextE(Faces,NextF,NextE);
        }
    }while ((!HasBorder));
}

void GetSortedVertStar(const Pos3 &StartPos,
                       const std::vector<std::vector<int> > &Faces,
                       const std::vector<std::vector<int> > &NextF,
                       const std::vector<std::vector<int> > &NextE,
                       std::vector<int > &VertFan,
                       bool &HasBorder)
{
    VertFan.clear();

    std::vector<Pos3> PosFan;
    GetSortedFacesFaceFun(StartPos,Faces,NextF,NextE,PosFan,HasBorder);
    for (size_t i=0;i<PosFan.size();i++)
    {
        assert(PosFan[i].V==StartPos.V);
        PosFan[i].FlipV(Faces);
        VertFan.push_back(PosFan[i].V);
    }
}

void getBorderSequnces(const std::vector<std::vector<int> > &Faces,
                       const std::vector<std::vector<int> > &NextF,
                       const std::vector<std::vector<int> > &NextE,
                       std::vector<std::vector<Pos3 > > &BorderPos)
{
    BorderPos.clear();
    int MaxV=Geo::FindMaxVertIdx(Faces);
    assert(MaxV>=0);
    std::vector<bool> Tested(MaxV+1,false);

    for (size_t i=0;i<NextF.size();i++)
        for (size_t j=0;j<NextF[i].size();j++)
        {
            //not a border
            if (NextF[i][j]!=-1)continue;

            int FSize=Faces[i].size();
            int VertIdx0=Faces[i][j];
            int VertIdx1=Faces[i][(j+1)%FSize];
            if (Tested[VertIdx0]&&Tested[VertIdx1])continue;

            //otherwise start goind around the border
            BorderPos.resize(BorderPos.size()+1);
            Pos3 StartPos(i,j,VertIdx0);
            Pos3 CurrPos=StartPos;
            do{
                Tested[CurrPos.V]=true;
                BorderPos.back().push_back(CurrPos);
                //go next border
                CurrPos.NextB(Faces,NextF,NextE);
                //then on other side
                CurrPos.FlipV(Faces);
            }while (CurrPos!=StartPos);
        }
}

void getBorderSequncesVert(const std::vector<std::vector<int> > &Faces,
                           const std::vector<std::vector<int> > &NextF,
                           const std::vector<std::vector<int> > &NextE,
                           std::vector<std::vector<int > > &BorderVIdx)
{
    BorderVIdx.clear();
    std::vector<std::vector<Pos3 > > BorderPos;
    getBorderSequnces(Faces,NextF,NextE,BorderPos);
    BorderVIdx.resize(BorderPos.size());
    for (size_t i=0;i<BorderPos.size();i++)
    {
        for (size_t j=0;j<BorderPos[i].size();j++)
            BorderVIdx[i].push_back(BorderPos[i][j].V);

        assert(BorderVIdx[i].size()>1);
        assert(BorderVIdx[i][0]!=BorderVIdx[i].back());
    }
}

void getBorderSequncesVert(const std::vector<std::vector<int> > &Faces,
                           std::vector<std::vector<int > > &BorderVIdx)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Faces,NextF,NextE);
    getBorderSequncesVert(Faces,NextF,NextE,BorderVIdx);
}

#endif
