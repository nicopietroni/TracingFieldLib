#ifndef FRAME_FIELD_UTILS
#define FRAME_FIELD_UTILS

#include <cmath>
#include <vector>
#include <stdio.h>
#include "point3.h"
//#include "pos4.h"
#include "frame_field.h"
#include "cross_field.h"
#include <set>

//template <class ScalarType>
//void GetPolylineFromSing(const std::vector<Geo::Point3<ScalarType> > &VertPos,
//                         const std::vector<std::vector<int> > &Tetra,
//                         const std::vector<Pos4> &FrameSingPos,
//                         const std::vector<int> &FrameSingValue,
//                         std::vector<std::vector<Geo::Point3<ScalarType> > > &PolyLSing,
//                         std::vector<Geo::Point3<ScalarType> > &PolyColorSing)
//{
//    PolyLSing.clear();
//    PolyColorSing.clear();

//    PolyLSing.resize(FrameSingPos.size());
//    PolyColorSing.resize(FrameSingPos.size(),Geo::Point3<ScalarType>(0,0,0));

//    for (size_t i=0;i<FrameSingPos.size();i++)
//    {
//        Pos4 CurrPos=FrameSingPos[i];

//        int IndexV0=CurrPos.V;
//        CurrPos.FlipV(Tetra);
//        int IndexV1=CurrPos.V;
//        assert(IndexV1!=IndexV0);
//        Geo::Point3<ScalarType> P0=VertPos[IndexV0];
//        Geo::Point3<ScalarType> P1=VertPos[IndexV1];

//        PolyLSing[i].push_back(P0);
//        PolyLSing[i].push_back(P1);

//        Geo::Point3<ScalarType> Col(0,1,0);

//        if (FrameSingValue[i]>0)Col=Geo::Point3<ScalarType>(1,0,0);
//        if (FrameSingValue[i]<0)Col=Geo::Point3<ScalarType>(1,1,0);

//        PolyColorSing[i]=Col;
//    }
//}

//template <class ScalarType>
//bool FrameSingularValue(const std::vector<Field::FrameField<ScalarType> > &SortedFaceFrames)//,int &IndexDir)
//{
//    //SingVal=0;
//    //    IndexDir=-1;
//    //    bool IsSing=false
//    //std::cout<<"size fan:"<<SortedFaceFrames.size()<<std::endl;
//    for (size_t d=0;d<3;d++)
//    {
//        int CurrD=d;
//        int sizeCross=SortedFaceFrames.size();
//        for (size_t i=0;i<SortedFaceFrames.size();i++)
//        {
//            Field::FrameField<ScalarType> currF=SortedFaceFrames[i];
//            Field::FrameField<ScalarType> NextF=SortedFaceFrames[(i+1)%sizeCross];

//            Geo::Point3<ScalarType> TargetDir=currF.GetDirection(CurrD);
//            CurrD=NextF.BestDirectionI(TargetDir);
//        }

//        //get the singularity value and the one where it not variates
//        if (CurrD>=3)
//            CurrD-=3;

//        if (CurrD!=d)
//        {
//            //SingVal=(CurrD-d);
//            //IndexDir=d;
//            return true;
//        }
//    }
//    return false;
//}



template <class ScalarType>
void CrossFieldFromFrame(const std::vector<Field::FrameField<ScalarType> > &TetraFrames,
                         const std::vector<int> &TriTetraIndex,
                         const std::vector<std::vector<int> > &Tri,
                         const std::vector<Geo::Point3<ScalarType> > &TriNorm,
                         const std::vector<Geo::Point3<ScalarType> > &TriPos,
                         std::vector<Field::CrossF<ScalarType> > &FaceCross,
                         bool force_tangent_surface=true)
{
    assert(Field::IsProperTetraFrameField(TetraFrames));
    assert(TriTetraIndex.size()==Tri.size());
    const ScalarType eps=(ScalarType)0.00001;

    int warning_normals=0;
    FaceCross.resize(Tri.size());
    for (size_t i=0;i<TriTetraIndex.size();i++)
    {
        int IndexTetra=TriTetraIndex[i];
        assert(IndexTetra>=0);
        assert(IndexTetra<TetraFrames.size());

        Field::FrameField<ScalarType> CurrFrame=TetraFrames[IndexTetra];

        //get direction aligned with tri normal
        Geo::Point3<ScalarType> CurrTriNorm=TriNorm[i];
        int indexFrameNormal=CurrFrame.BestDirectionI(CurrTriNorm);
        Geo::Point3<ScalarType> NormFrame=CurrFrame.GetDirection(indexFrameNormal);

        //get the two other directions
        Geo::Point3<ScalarType> CrossD0=CurrFrame.GetDirection((indexFrameNormal+1)%6);
        Geo::Point3<ScalarType> CrossD1=CurrFrame.GetDirection((indexFrameNormal+2)%6);
        assert(Tri[i].size()==3);
        Geo::Point3<ScalarType> Pos0=TriPos[Tri[i][0]];
        Geo::Point3<ScalarType> Pos1=TriPos[Tri[i][1]];
        Geo::Point3<ScalarType> Pos2=TriPos[Tri[i][2]];
        Geo::Point3<ScalarType> Pos=(Pos0+Pos1+Pos2)/3;
        Field::CrossF<ScalarType> CopiedFrame(CrossD0,CrossD1,1,1,Pos);
        CopiedFrame.Q=CurrFrame.Q;
        CopiedFrame.MakeDirectionNormalCoherent(CurrTriNorm);

        //check
        if ((NormFrame-CurrTriNorm).Norm()>eps)
        {
            //std::cout<<"Diff:"<<(NormFrame-CurrTriNorm).Norm()<<std::endl;
            warning_normals++;
            //CopiedFrame.MakeTangentTo(CurrTriNorm);
        }
        if (force_tangent_surface)
            CopiedFrame.MakeTangentTo(CurrTriNorm);

        FaceCross[i]=CopiedFrame;
    }

    if (warning_normals>0)
        std::cout<<"* Warning Volume to Surface Mapping "<<warning_normals<<" faces"<<std::endl;
}

//template <class ScalarType>
//void GetFrameSingularities(const std::vector<Field::FrameField<ScalarType> > &TetraFrames,
//                           const std::vector<std::vector<int> > &Tetra,
//                           const std::vector<std::vector<int> > &NextT,
//                           const std::vector<std::vector<int> > &NextF,
//                           std::vector<Pos4 > &SingPos,
//                           std::vector<int> &SingValue)
//{
//    std::vector<std::pair<int,int> > BorderE;
//    GetBorderEdges(Tetra,NextT,NextF,BorderE);

//    std::set<std::pair<int,int> > BorderESet(BorderE.begin(),BorderE.end());

//    SingPos.clear();
//    SingValue.clear();

//    std::set<std::pair<int,int> > VisitedE;

//    //for each tetra
//    int currT,currF,currE,currV;
//    for (size_t t=0;t<Tetra.size();t++)
//    {
//        currT=t;
//        //for each face`
//        assert(Tetra[currT].size()==4);
//        for (size_t e=0;e<6;e++)
//        {
//            currE=e;

//            currF=Tetra::FofE(currE,0);

//            int currV0=Tetra::VofE(currE,0);
//            int currV1=Tetra::VofE(currE,1);

//            currV0=Tetra[currT][currV0];
//            currV1=Tetra[currT][currV1];

//            currV=currV0;
//            std::pair<int,int> EntryEV(std::min(currV0,currV1),std::max(currV0,currV1));

//            if (BorderESet.count(EntryEV)>0)continue;

//            //std::cout<<"EntryEV "<<EntryEV.first<<","<<EntryEV.second<<std::endl;

//            if (VisitedE.count(EntryEV)>0)continue;

//            VisitedE.insert(EntryEV);

//            Pos4 currPos(currT,currF,currE,currV);

//            std::vector<Pos4> PosFan;
//            bool HasBorder;
//            int err=currPos.CheckConsistency(Tetra,NextT,NextF);
//            assert(err==0);
//            GetSortedTetraFun(currPos,Tetra,NextT,NextF,PosFan,HasBorder);
//            assert(!HasBorder);

//            std::vector<Field::FrameField<ScalarType> > SortedFan;
//            for (size_t k=0;k<PosFan.size();k++)
//            {
//                int IndexT=PosFan[k].T;
//                SortedFan.push_back(TetraFrames[IndexT]);
//            }
//            //            int SingVal;
//            //            int IndexDir;
//            bool IsSing=FrameSingularValue(SortedFan);//,SingVal,IndexDir);

//            if (IsSing)
//            {
//                SingPos.push_back(currPos);
//                SingValue.push_back(1);
//            }
//            //            if (SingVal!=0)
//            //            {
//            //                SingPos.push_back(currPos);
//            //                SingValue.push_back(SingVal);
//            //            }
//        }
//    }
//}


template <class ScalarType>
void SetBorderFrameFromCross(const std::vector<std::vector<int> > &Tri,
                    const std::vector<Field::CrossF<ScalarType> > &FaceCross,
                    const std::vector<int> &TriTetraIndex,
                    const std::vector<Geo::Point3<ScalarType> > &TetraVertPos,
                    const std::vector<std::vector<int> > &Tetra,
                    std::vector<Field::FrameField<ScalarType> > &TetraFrames)
{
    assert(TriTetraIndex.size()>0);
    assert(TriTetraIndex.size()==Tri.size());
    assert(TriTetraIndex.size()==FaceCross.size());
//    TetraFrames.clear();
//    TetraFrames.resize(Tetra.size());
    for (size_t i=0;i<TriTetraIndex.size();i++)
    {
        int IndexT=TriTetraIndex[i];
        assert(IndexT>=0);
        assert(IndexT<TetraFrames.size());

        Geo::Point3<ScalarType> D0=FaceCross[i].D[0];
        Geo::Point3<ScalarType> D1=FaceCross[i].D[1];
        Geo::Point3<ScalarType> D2=FaceCross[i].OrthoDir();
        ScalarType Q=FaceCross[i].Q;

        assert(Tetra[IndexT].size()==4);
        Geo::Point3<ScalarType> P0=TetraVertPos[Tetra[IndexT][0]];
        Geo::Point3<ScalarType> P1=TetraVertPos[Tetra[IndexT][1]];
        Geo::Point3<ScalarType> P2=TetraVertPos[Tetra[IndexT][2]];
        Geo::Point3<ScalarType> P3=TetraVertPos[Tetra[IndexT][3]];

        Tetra::Tetra3<ScalarType> CurrT(P0,P1,P2,P3);
        TetraFrames[IndexT].D[0]=D0;
        TetraFrames[IndexT].D[1]=D1;
        TetraFrames[IndexT].D[2]=D2;

        TetraFrames[IndexT].Pos=CurrT.Barycenter();
        TetraFrames[IndexT].Q=Q;
    }
}

template <class ScalarType>
void InitDefaultFrame(const std::vector<Geo::Point3<ScalarType> > &TetraVertPos,
                      const std::vector<std::vector<int> > &Tetra,
                      std::vector<Field::FrameField<ScalarType> > &TetraFrames)
{
    TetraFrames.clear();
    std::vector<Geo::Point3<ScalarType> > Centers;
    Tetra::FindCenter(TetraVertPos,Tetra,Centers);
    Geo::Point3<ScalarType> D0(1,0,0);
    Geo::Point3<ScalarType> D1(0,1,0);
    Geo::Point3<ScalarType> D2(0,0,1);
    for (size_t i=0;i<Tetra.size();i++)
    {
//        Geo::Point3<ScalarType> P0=TetraVertPos[Tetra[i][0]];
//        Geo::Point3<ScalarType> P1=TetraVertPos[Tetra[i][1]];
//        Geo::Point3<ScalarType> P2=TetraVertPos[Tetra[i][2]];
//        Geo::Point3<ScalarType> P3=TetraVertPos[Tetra[i][3]];
        //Geo::Point3<ScalarType> Bary=(P0+P1+P2+P3)/4.0;
        TetraFrames.push_back(Field::FrameField<ScalarType>(D0,D1,D2,Centers[i]));
    }
}


#endif
